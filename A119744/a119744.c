/*
 * a119744 - explore OEIS A119744 and friends.
 *
 * a(1) = seed; a(n) = 1-based position of the decimal string of a(n-1)
 * in the decimal expansion of Pi written as the digit string
 * "31415926535897932384626433832795..." (leading 3 included).
 * E.g. "13" first occurs at position 111, so 13 -> 111 -> 154 -> ...
 *
 * The pi file is streamed, never loaded into memory, so it may be
 * arbitrarily large (4 TB+). Non-digit bytes (a "3." prefix, newlines,
 * spaces) are ignored. If the file's digits start with "14..." rather
 * than "31...", a leading 3 is supplied automatically so positions
 * match the OEIS convention.
 *
 * Usage: a119744 [-s seed] [-n maxterms] [-b bufMB] pifile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <errno.h>

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#endif

#define MAX_PAT 24  /* digits in a pattern; 24 covers any uint64 */

static size_t buf_size = 8u << 20;  /* 8 MiB read buffer */

/*
 * Stream the digit characters of the pi file through a KMP matcher for
 * `pat` (length m). Returns the 1-based position of the first match,
 * 0 if not found by EOF, or -1 on read error.
 * *digits_seen gets the total count of digit chars scanned.
 * prepend3: feed a synthetic leading '3' before the file contents.
 */
/* Progress display: a "\r"-overwritten status line on stderr, only when
 * stderr is a tty and the scan is long enough to matter. */
#define PROGRESS_MIN (16u << 20)  /* don't bother below 16M digits */

static void progress_clear(int *shown)
{
    if (*shown) {
        fprintf(stderr, "\r%60s\r", "");
        *shown = 0;
    }
}

static int64_t find_first(FILE *f, const char *pat, int m, int prepend3,
                          uint64_t limit, uint64_t prog_total,
                          unsigned char *buf, uint64_t *digits_seen)
{
    int fail[MAX_PAT + 1];
    int k = 0;               /* current KMP state */
    uint64_t pos = 0;        /* digits consumed so far */
    int shown = 0;           /* progress line currently on screen */
    uint32_t last_tenth = (uint32_t)-1;

    /* KMP failure function */
    fail[0] = fail[1] = 0;
    for (int i = 1; i < m; i++) {
        int j = fail[i];
        while (j && pat[i] != pat[j]) j = fail[j];
        fail[i + 1] = (pat[i] == pat[j]) ? j + 1 : 0;
    }

    if (fseeko(f, 0, SEEK_SET) != 0) return -1;

    if (prepend3) {
        pos = 1;
        while (k && pat[k] != '3') k = fail[k];
        if (pat[k] == '3') k++;
        if (k == m) return 1;   /* pattern was "3" itself */
        if (limit && pos >= limit) { *digits_seen = pos; return 0; }
    }

    for (;;) {
        size_t got = fread(buf, 1, buf_size, f);
        if (got == 0) {
            if (ferror(f)) { progress_clear(&shown); return -1; }
            break;  /* EOF */
        }
        for (size_t i = 0; i < got; i++) {
            unsigned char c = buf[i];
            if (c < '0' || c > '9') continue;  /* skip ".", newlines, ... */
            pos++;
            while (k && pat[k] != (char)c) k = fail[k];
            if (pat[k] == (char)c) k++;
            if (k == m) {
                progress_clear(&shown);
                *digits_seen = pos;
                return (int64_t)(pos - (uint64_t)m + 1);
            }
            if (limit && pos >= limit) {
                progress_clear(&shown);
                *digits_seen = pos;
                return 0;
            }
        }
        if (prog_total && pos >= PROGRESS_MIN) {
            uint32_t tenth = (uint32_t)((pos * 1000) / prog_total);
            if (tenth > 1000) tenth = 1000;
            if (tenth != last_tenth) {
                fprintf(stderr, "\rsearching \"%s\": %u.%u%% ",
                        pat, tenth / 10, tenth % 10);
                last_tenth = tenth;
                shown = 1;
            }
        }
    }
    progress_clear(&shown);
    *digits_seen = pos;
    return 0;
}

/* Parse a digit count like "5000000", "500M", "2.5b", "1e10".
 * K/M/B/T = 1e3/1e6/1e9/1e12, case-insensitive. Returns 0 on bad input
 * (0 is not a useful limit, so it doubles as the error value). */
static uint64_t parse_amount(const char *s)
{
    char *end;
    double v = strtod(s, &end);
    if (end == s || v < 0) return 0;
    switch (*end) {
    case 'k': case 'K': v *= 1e3;  end++; break;
    case 'm': case 'M': v *= 1e6;  end++; break;
    case 'b': case 'B': v *= 1e9;  end++; break;
    case 't': case 'T': v *= 1e12; end++; break;
    default: break;  /* plain number or e-notation, handled by strtod */
    }
    if (*end != '\0' || v > 1.8e19) return 0;
    return (uint64_t)v;
}

static void usage(const char *argv0)
{
    fprintf(stderr,
        "usage: %s [-s seed] [-n digits] [--max terms] [-b bufMB] [--short] pifile\n"
        "  -s seed      starting term (default 13, i.e. A119744), or a\n"
        "               range like 1-15 to run each seed in turn; a\n"
        "               sequence that reaches a value from an earlier\n"
        "               seed stops with \"(seq nnn offset mmm)\"\n"
        "  -n digits    search at most this many digits of pi, so a large\n"
        "               file can be probed shallowly (default: whole file);\n"
        "               accepts K/M/B/T suffixes and e-notation, e.g.\n"
        "               500M, 2.5b, 1e10\n"
        "  --max terms  stop after this many terms (default: until loop,\n"
        "               join, or term not found)\n"
        "  -b bufMB     read buffer size in MiB (default 8)\n"
        "  --short      print terms comma-separated on one line; on a loop\n"
        "               just stop; if a term isn't found print \">n\" where\n"
        "               n is the number of pi digits searched\n"
        "  pifile       pi digits as text; may start with \"3.\", \"3\", or\n"
        "               \"14...\"; non-digit bytes are ignored; streamed,\n"
        "               never loaded into memory\n", argv0);
}

int main(int argc, char **argv)
{
    uint64_t seed_lo = 13, seed_hi = 13;
    uint64_t max_terms = 0;    /* 0 = unlimited */
    uint64_t digit_limit = 0;  /* 0 = whole file */
    int short_mode = 0;
    const char *path = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--short")) {
            short_mode = 1;
        } else if (!strcmp(argv[i], "-s") && i + 1 < argc) {
            char *end;
            seed_lo = strtoull(argv[++i], &end, 10);
            if (*end == '-') {
                seed_hi = strtoull(end + 1, &end, 10);
            } else {
                seed_hi = seed_lo;
            }
            if (*end != '\0' || seed_hi < seed_lo) {
                fprintf(stderr, "bad seed or seed range: %s\n", argv[i]);
                usage(argv[0]); return 2;
            }
        } else if (!strcmp(argv[i], "-n") && i + 1 < argc) {
            digit_limit = parse_amount(argv[++i]);
            if (digit_limit == 0) {
                fprintf(stderr, "bad digit count: %s\n", argv[i]);
                usage(argv[0]); return 2;
            }
        } else if (!strcmp(argv[i], "--max") && i + 1 < argc) {
            max_terms = strtoull(argv[++i], NULL, 10);
        } else if (!strcmp(argv[i], "-b") && i + 1 < argc) {
            uint64_t mb = strtoull(argv[++i], NULL, 10);
            if (mb == 0 || mb > 4096) { usage(argv[0]); return 2; }
            buf_size = (size_t)mb << 20;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(argv[0]); return 0;
        } else if (!path) {
            path = argv[i];
        } else {
            usage(argv[0]); return 2;
        }
    }
    if (!path) { usage(argv[0]); return 2; }

    FILE *f = fopen(path, "rb");
    if (!f) { fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno)); return 1; }

    unsigned char *buf = malloc(buf_size);
    if (!buf) { fprintf(stderr, "cannot allocate %zu-byte buffer\n", buf_size); return 1; }
    setvbuf(f, NULL, _IONBF, 0);  /* we do our own buffering */
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fileno(f), 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

    /* Total digits a search can cover, for the stderr progress line.
     * File size in bytes is a close-enough estimate of its digit count;
     * 0 disables progress (stderr not a tty). */
    uint64_t prog_total = 0;
#ifndef _WIN32
    if (isatty(fileno(stderr))) {
        struct stat st;
        if (fstat(fileno(f), &st) == 0 && st.st_size > 0)
            prog_total = (uint64_t)st.st_size;
        if (digit_limit && (prog_total == 0 || digit_limit < prog_total))
            prog_total = digit_limit;
    }
#endif

    /* Detect whether the file's digit stream already includes the
     * leading 3 of pi: peek at the first two digit chars. */
    int prepend3 = 0;
    {
        int d1 = -1, d2 = -1, c;
        while ((c = fgetc(f)) != EOF) {
            if (c >= '0' && c <= '9') {
                if (d1 < 0) d1 = c;
                else { d2 = c; break; }
            }
        }
        if (d1 < 0) { fprintf(stderr, "%s contains no digits\n", path); return 1; }
        if (d1 == '3' && d2 == '1') {
            prepend3 = 0;
        } else if (d1 == '1' && d2 == '4') {
            prepend3 = 1;
            fprintf(stderr, "note: file starts \"14...\"; supplying the "
                            "leading 3 so positions match OEIS.\n");
        } else {
            fprintf(stderr, "warning: file digits start \"%c%c...\", which "
                            "looks like neither \"31\" nor \"14\"; using it "
                            "as-is.\n", d1, d2 < 0 ? '?' : d2);
        }
    }

    /* Seen terms, for loop detection. The published lists have a few
     * dozen terms at most before exceeding any real pi file, so a
     * growable array with linear scan is plenty. */
    uint64_t *seen = NULL;
    size_t nseen = 0, capseen = 0;

    /* Terms from earlier seeds in a range run, so a later sequence can
     * stop as soon as it joins one already computed. */
    struct gterm { uint64_t val, seed, idx; };
    struct gterm *gseen = NULL;
    size_t ngseen = 0, capgseen = 0;

    /* Per-seed outcome, for the end-of-range summary. res_arr resolves
     * joins to the fate of the sequence they merge into. */
    enum { R_LOOP, R_UNBOUND, R_JOIN, R_TRUNC };
    size_t nseeds = (size_t)(seed_hi - seed_lo) + 1;
    unsigned char *raw_arr = malloc(nseeds);
    unsigned char *res_arr = malloc(nseeds);
    if (!raw_arr || !res_arr) { fprintf(stderr, "out of memory\n"); return 1; }
    uint64_t file_digits = 0;  /* learned when a search reaches EOF */

    int exit_code = 0;
    int multi = seed_hi > seed_lo;
    for (uint64_t seed = seed_lo; seed <= seed_hi && exit_code == 0; seed++) {
    nseen = 0;
    int raw = R_TRUNC;       /* outcome if we fall out via max_terms */
    uint64_t jt = seed;      /* join target when raw == R_JOIN */
    uint64_t term = seed;
    if (short_mode) {
        if (multi) printf("s%" PRIu64 ": ", seed);
        printf("%" PRIu64, term);
    } else {
        if (multi) printf("=== seed %" PRIu64 " ===\n", seed);
        printf("a(1) = %" PRIu64 "\n", term);
    }
    fflush(stdout);

    for (uint64_t n = 2; ; n++) {
        if (max_terms && n > max_terms) {
            if (short_mode)
                printf(",...");
            else
                printf("stopped: reached max terms (%" PRIu64 ")\n", max_terms);
            break;
        }

        /* record previous term, checking for a loop */
        int looped = 0;
        for (size_t i = 0; i < nseen; i++) {
            if (seen[i] == term) {
                if (short_mode) {
                    printf(" (loops)");
                } else {
                    size_t cyclen = nseen - i;
                    printf("loop detected: a(%zu) = a(%" PRIu64 ") = %" PRIu64
                           "  (cycle length %zu: ", i + 1, n - 1, term, cyclen);
                    for (size_t j = i; j < nseen; j++)
                        printf("%s%" PRIu64, j > i ? "," : "", seen[j]);
                    printf(")\n");
                }
                raw = R_LOOP;
                looped = 1;
                break;
            }
        }
        if (looped) break;

        /* does this term join an earlier seed's sequence? */
        for (size_t i = 0; i < ngseen; i++) {
            if (gseen[i].val == term) {
                if (short_mode)
                    printf(" (seq %" PRIu64 " offset %" PRIu64 ")",
                           gseen[i].seed, gseen[i].idx);
                else
                    printf("join detected: a(%" PRIu64 ") = %" PRIu64
                           " = term %" PRIu64 " of the seed-%" PRIu64
                           " sequence; it continues identically from there.\n",
                           n - 1, term, gseen[i].idx, gseen[i].seed);
                raw = R_JOIN;
                jt = gseen[i].seed;
                looped = 1;
                break;
            }
        }
        if (looped) break;

        if (nseen == capseen) {
            capseen = capseen ? capseen * 2 : 64;
            seen = realloc(seen, capseen * sizeof *seen);
            if (!seen) { fprintf(stderr, "out of memory\n"); return 1; }
        }
        seen[nseen++] = term;

        char pat[MAX_PAT];
        int m = snprintf(pat, sizeof pat, "%" PRIu64, term);

        uint64_t digits_seen = 0;
        int64_t pos = find_first(f, pat, m, prepend3, digit_limit,
                                 prog_total, buf, &digits_seen);
        if (pos < 0) {
            fprintf(stderr, "read error on %s: %s\n", path, strerror(errno));
            exit_code = 1;
            break;
        }
        if (pos == 0) {
            raw = R_UNBOUND;
            file_digits = digits_seen;
            if (short_mode)
                printf(",>%" PRIu64, digits_seen);
            else
                printf("stopped: \"%s\" not found in the first %" PRIu64
                       " digits of the file; a larger pi file is needed.\n",
                       pat, digits_seen);
            break;
        }
        term = (uint64_t)pos;
        if (short_mode)
            printf(",%" PRIu64, term);
        else
            printf("a(%" PRIu64 ") = %" PRIu64 "   (\"%s\" found at position %"
                   PRIu64 ", scanned %" PRIu64 " digits)\n",
                   n, term, pat, term, digits_seen);
        fflush(stdout);
    }
    if (short_mode)
        printf("\n");

    raw_arr[seed - seed_lo] = (unsigned char)raw;
    res_arr[seed - seed_lo] =
        (raw == R_JOIN) ? res_arr[jt - seed_lo] : (unsigned char)raw;

    /* Register this seed's terms for join detection by later seeds.
     * A sequence stops on its first already-known value, so every
     * recorded term here is globally new — no dup check needed. */
    for (size_t i = 0; i < nseen; i++) {
        if (ngseen == capgseen) {
            capgseen = capgseen ? capgseen * 2 : 64;
            gseen = realloc(gseen, capgseen * sizeof *gseen);
            if (!gseen) { fprintf(stderr, "out of memory\n"); return 1; }
        }
        gseen[ngseen].val = seen[i];
        gseen[ngseen].seed = seed;
        gseen[ngseen].idx = (uint64_t)i + 1;
        ngseen++;
    }
    }

    if (multi && exit_code == 0) {
        int c_loop = 0, c_unb = 0, c_uniq = 0, c_ujoin = 0, c_undet = 0;
        for (size_t i = 0; i < nseeds; i++) {
            if (res_arr[i] == R_LOOP) c_loop++;
            else if (res_arr[i] == R_UNBOUND) {
                c_unb++;
                if (raw_arr[i] == R_UNBOUND) c_uniq++; else c_ujoin++;
            } else c_undet++;
        }

        int show_lists = nseeds <= 25;

#define SEED_LIST(cond) do {                                            \
            if (!show_lists) { printf("\n"); break; }                   \
            int first_ = 1;                                             \
            printf("  [");                                              \
            for (size_t i = 0; i < nseeds; i++)                         \
                if (cond) {                                             \
                    printf("%s%" PRIu64, first_ ? "" : ",",             \
                           seed_lo + (uint64_t)i);                      \
                    first_ = 0;                                         \
                }                                                       \
            printf("]\n");                                              \
        } while (0)

        printf("\n--- summary: seeds %" PRIu64 "-%" PRIu64 " (%zu sequences) ---\n",
               seed_lo, seed_hi, nseeds);
        printf("looping (or joins a looping sequence): %d", c_loop);
        SEED_LIST(res_arr[i] == R_LOOP);
        if (file_digits)
            printf("unbounded (no end within %" PRIu64 " digits): %d",
                   file_digits, c_unb);
        else
            printf("unbounded (no end within the pi file): %d", c_unb);
        SEED_LIST(res_arr[i] == R_UNBOUND);
        printf("  unique unbounded (hit end of file themselves): %d", c_uniq);
        SEED_LIST(res_arr[i] == R_UNBOUND && raw_arr[i] == R_UNBOUND);
        printf("  joined another unbounded sequence: %d", c_ujoin);
        SEED_LIST(res_arr[i] == R_UNBOUND && raw_arr[i] == R_JOIN);
        if (c_undet) {
            printf("undetermined (stopped at --max limit, or joined one): %d",
                   c_undet);
            SEED_LIST(res_arr[i] == R_TRUNC);
        }
#undef SEED_LIST
    }

    free(raw_arr);
    free(res_arr);
    free(gseen);
    free(seen);
    free(buf);
    fclose(f);
    return exit_code;
}
