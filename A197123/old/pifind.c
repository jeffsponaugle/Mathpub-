/*
 * pifind - find every occurrence of one or more digit strings in a pi file.
 *
 *   ./pifind -p pi-10b.txt 8530614 1350168131352524443
 *
 * All numbers are searched in a single pass over the file.  Positions are
 * 1-based from the first digit after the decimal point (same convention as
 * OEIS A197123, pibloom and pisearch).  Each number may be 1..30 digits and
 * may have leading zeros.
 *
 * Build:  cc -O3 -march=native -o pifind pifind.c -lm
 */
#include "pi_common.h"

#define MAX_PAT_LEN 30
#define MAX_PATS    64
#define CTX 3                    /* context digits before/after */

typedef struct {
    char     str[MAX_PAT_LEN + 1];
    uint8_t  dig[MAX_PAT_LEN];   /* digit values 0..9 */
    size_t   len;
    uint64_t found;
} pat_t;

static struct {
    const char *pifile;
    uint64_t ndigits;            /* limit (UINT64_MAX = all) */
    uint64_t skip;               /* digits to skip before searching */
    uint64_t maxprint;           /* per number: stop printing (not counting) */
    size_t chunk;
} opt = { NULL, UINT64_MAX, 0, 1000, 64u << 20 };

static void usage(const char *argv0)
{
    fprintf(stderr,
"Usage: %s -p <pifile> [options] <number> [<number> ...]\n"
"\n"
"Prints every position where each number (1..%d digits, leading zeros allowed,\n"
"up to %d numbers) occurs in the decimal expansion of pi, with %d digits of\n"
"context on each side and the match underlined.  All numbers are searched in\n"
"one pass over the file.  Positions are 1-based after the decimal point.\n"
"\n"
"  -p FILES    pi source file, or comma list of up to 8 continuation files\n"
"  -n COUNT    digits of pi to scan, counted after any skip (500M, 10B, all) [all]\n"
"  -s COUNT    skip this many digits before searching (e.g. -s 1T); reported\n"
"              positions stay absolute (1-based after the decimal point)\n"
"  -m N        print at most N matches per number (still counts) [1000]\n"
"  -c SIZE     read chunk size [64M]\n"
"  -h          this help\n", argv0, MAX_PAT_LEN, MAX_PATS, CTX);
    exit(1);
}

int main(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:n:m:s:c:h")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 's': if (parse_count(optarg, &opt.skip) || opt.skip == UINT64_MAX) die("bad skip count '%s'", optarg); break;
        case 'm': opt.maxprint = strtoull(optarg, NULL, 10); break;
        case 'c': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile) { fprintf(stderr, "missing -p <pifile>\n"); usage(argv[0]); }
    if (optind >= argc) { fprintf(stderr, "at least one number to search for is required\n"); usage(argv[0]); }
    if (argc - optind > MAX_PATS) die("too many numbers (max %d)", MAX_PATS);
    if (opt.ndigits == 0) die("-n must be > 0");

    /* parse the patterns */
    pat_t pats[MAX_PATS];
    int npats = 0;
    size_t Lmax = 0;
    for (int a = optind; a < argc; a++) {
        const char *ps = argv[a];
        size_t L = strlen(ps);
        if (L < 1 || L > MAX_PAT_LEN) die("number '%s' must be 1..%d digits long (got %zu)", ps, MAX_PAT_LEN, L);
        int dup = 0;
        for (int q = 0; q < npats; q++)
            if (pats[q].len == L && !memcmp(pats[q].str, ps, L)) { dup = 1; break; }
        if (dup) { warn("duplicate number '%s' ignored", ps); continue; }
        pat_t *P = &pats[npats++];
        memcpy(P->str, ps, L + 1);
        P->len = L;
        P->found = 0;
        for (size_t i = 0; i < L; i++) {
            if (ps[i] < '0' || ps[i] > '9') die("'%s' is not a decimal number", ps);
            P->dig[i] = (uint8_t)(ps[i] - '0');
        }
        if (L > Lmax) Lmax = L;
    }

    int use_ansi = isatty(1);
    const char *ul_on  = use_ansi ? "\033[4m"  : "[";
    const char *ul_off = use_ansi ? "\033[24m" : "]";
    int show_progress = isatty(2);              /* status line goes to stderr */

    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();

    /* We manage buffering ourselves: matches need CTX digits of context on
     * both sides across chunk boundaries.  Non-final chunks search starts
     * [search_from_p .. ndigits-Lp-CTX] per pattern so every printed match
     * has full after-context; the unsearched tail (sized by the LONGEST
     * pattern) is carried into the next chunk.  At EOF the remaining starts
     * are searched with whatever context exists. */
    /* -p may be a comma-separated list of up to 8 files read as one stream */
    char *pathbuf = strdup(opt.pifile);
    if (!pathbuf) die("out of memory");
    const char *paths[8];
    uint64_t data_bytes[8];              /* digit-bearing bytes per file */
    int nfiles = 0, cur_file = 0;
    for (char *save = pathbuf, *tok; (tok = strsep(&save, ",")) != NULL; ) {
        if (!*tok) continue;
        if (nfiles >= 8) die("too many pi source files (max 8)");
        paths[nfiles++] = tok;
    }
    if (nfiles == 0) die("empty pi source file list");

    int fd = open(paths[0], O_RDONLY);
    if (fd < 0) die("cannot open pi source '%s': %s", paths[0], strerror(errno));
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
    fcntl(fd, F_RDAHEAD, 1);
#endif
    /* skip a leading "3." / "3" (first file only; the rest are continuations) */
    off_t prefix_skip = 0;
    {
        unsigned char head[4];
        ssize_t n = read(fd, head, sizeof head);
        if (n < 0) die("read error: %s", strerror(errno));
        if (n == 0) die("pi source is empty");
        if (n >= 2 && head[0] == '3' && head[1] == '.') prefix_skip = 2;
        else if (n >= 4 && head[0] == '3' && head[1] == '1' && head[2] == '4' && head[3] == '1') prefix_skip = 1;
        else if (!(n >= 3 && head[0] == '1' && head[1] == '4' && head[2] == '1'))
            warn("file does not start with '3.1415' or '1415' - assuming it begins after the decimal point");
        if (lseek(fd, prefix_skip, SEEK_SET) < 0) die("lseek failed: %s", strerror(errno));
    }
    for (int i = 0; i < nfiles; i++) {
        struct stat sti;
        if (stat(paths[i], &sti) != 0) die("cannot stat pi source '%s': %s", paths[i], strerror(errno));
        data_bytes[i] = (uint64_t)sti.st_size;
        if (i == 0) data_bytes[i] = data_bytes[i] > (uint64_t)prefix_skip ? data_bytes[i] - (uint64_t)prefix_skip : 0;
    }

    /* Skip opt.skip digits, backing up CTX digits so the first match still
     * gets its before-context.  Fast path: if a sample of the file shows pure
     * digits, byte offset == digit offset and we can lseek.  Otherwise stream
     * through the prefix counting digits. */
    uint64_t skip_left = 0;              /* digits still to discard (streaming path) */
    uint64_t skip_total = 0;             /* for the skip progress display */
    uint64_t pre_ctx = opt.skip < CTX ? opt.skip : CTX;
    if (opt.skip > 0) {
        uint64_t target = opt.skip - pre_ctx;
        off_t data_start = lseek(fd, 0, SEEK_CUR);
        char *sample = malloc(1 << 22);
        if (!sample) die("out of memory");
        ssize_t sn = read(fd, sample, 1 << 22);
        if (sn < 0) die("read error: %s", strerror(errno));
        int clean = 1;
        for (ssize_t i = 0; i < sn; i++)
            if (sample[i] < '0' || sample[i] > '9') { clean = 0; break; }
        free(sample);
        if (clean) {
            uint64_t rem = target;
            int fi = 0;
            while (fi < nfiles - 1 && rem >= data_bytes[fi]) { rem -= data_bytes[fi]; fi++; }
            if (fi == 0) {
                if (lseek(fd, data_start + (off_t)rem, SEEK_SET) < 0)
                    die("lseek failed: %s", strerror(errno));
            } else {
                close(fd);
                cur_file = fi;
                fd = open(paths[fi], O_RDONLY);
                if (fd < 0) die("cannot open pi source '%s': %s", paths[fi], strerror(errno));
#ifdef POSIX_FADV_SEQUENTIAL
                posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
                fcntl(fd, F_RDAHEAD, 1);
#endif
                if (lseek(fd, (off_t)rem, SEEK_SET) < 0)
                    die("lseek failed on '%s': %s", paths[fi], strerror(errno));
            }
        } else {
            warn("file contains non-digit bytes - skipping by counting digits (slower than a seek)");
            if (lseek(fd, data_start, SEEK_SET) < 0) die("lseek failed: %s", strerror(errno));
            skip_left = target;
            skip_total = target;
        }
    }

    /* total digits for the % display: limit, or file size minus the "3." */
    uint64_t total_digits = opt.ndigits;
    {
        uint64_t all = 0;
        for (int i = 0; i < nfiles; i++) all += data_bytes[i];   /* close enough; newlines are rare */
        uint64_t consumed = 0;                                    /* digits already seeked past */
        if (opt.skip > 0 && skip_left == 0) consumed = opt.skip - pre_ctx;
        uint64_t fd_digits = all > consumed ? all - consumed : 0;
        fd_digits = fd_digits > skip_left ? fd_digits - skip_left : 0;   /* streaming skip not yet consumed */
        if (total_digits == UINT64_MAX || total_digits > fd_digits) total_digits = fd_digits;
    }

    if (opt.ndigits != UINT64_MAX) opt.ndigits += pre_ctx;

    uint8_t *raw = malloc(opt.chunk);
    uint8_t *buf = malloc(opt.chunk + MAX_PAT_LEN + 2 * CTX + 64);
    if (!raw || !buf) die("cannot allocate %zu byte buffers", opt.chunk);

    if (opt.skip) {
        char sb[32];
        printf("Skipping the first %s digits.\n", fmt_u64(opt.skip, sb, sizeof sb));
    }
    if (npats == 1)
        printf("Searching for %s (%zu digits) in %s%s\n", pats[0].str, pats[0].len, opt.pifile,
               opt.ndigits == UINT64_MAX ? "" : " (limited scan)");
    else {
        printf("Searching for %d numbers in %s%s:\n", npats, opt.pifile,
               opt.ndigits == UINT64_MAX ? "" : " (limited scan)");
        for (int q = 0; q < npats; q++) printf("  %s (%zu digits)\n", pats[q].str, pats[q].len);
    }

    uint64_t found_total = 0, delivered = 0, skipped_bytes = 0;
    size_t ndigits = 0;
    uint64_t base = opt.skip - pre_ctx + 1;   /* absolute 1-based position of buf[0] */
    int first_chunk = 1;
    double t0 = now_sec();
    double t_prog = t0;
    int progress_shown = 0;
    int eof = 0, interrupted = 0;

#define CLEAR_PROGRESS() do { if (progress_shown) { fprintf(stderr, "\r\033[K"); progress_shown = 0; } } while (0)

    while (!eof && !interrupted) {
        /* fill buffer with cleaned digits */
        while (ndigits < opt.chunk && !eof) {
            if (delivered >= opt.ndigits) { eof = 1; break; }
            ssize_t n = read(fd, raw, opt.chunk - ndigits < opt.chunk ? opt.chunk - ndigits : opt.chunk);
            if (n < 0) {
                if (errno == EINTR) { if (g_stop) { interrupted = 1; break; } continue; }
                die("read error on pi source: %s", strerror(errno));
            }
            if (n == 0) {
                if (cur_file + 1 < nfiles) {
                    close(fd);
                    cur_file++;
                    fd = open(paths[cur_file], O_RDONLY);
                    if (fd < 0) die("cannot open pi source '%s': %s", paths[cur_file], strerror(errno));
#ifdef POSIX_FADV_SEQUENTIAL
                    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
                    fcntl(fd, F_RDAHEAD, 1);
#endif
                    unsigned char h2[2];
                    ssize_t hn = read(fd, h2, 2);
                    if (hn == 2 && h2[0] == '3' && h2[1] == '.')
                        die("continuation file '%s' starts with \"3.\" - it looks like a from-the-start pi file, not a continuation", paths[cur_file]);
                    if (hn > 0 && lseek(fd, 0, SEEK_SET) < 0) die("lseek failed: %s", strerror(errno));
                    continue;
                }
                eof = 1;
                break;
            }
            size_t out = 0;
            uint64_t room = opt.ndigits - delivered;
            for (ssize_t i = 0; i < n; i++) {
                unsigned d = (unsigned)raw[i] - '0';
                if (d < 10u) {
                    if (skip_left) { skip_left--; continue; }
                    if (out >= room) { eof = 1; break; }
                    buf[ndigits + out++] = (uint8_t)d;
                } else skipped_bytes++;
            }
            ndigits += out;
            delivered += out;
            if (out == 0 && skipped_bytes > (uint64_t)opt.chunk * 4)
                die("no digits found after %" PRIu64 " junk bytes - wrong file?", skipped_bytes);
            if (g_stop) { interrupted = 1; break; }
        }
        if (show_progress && skip_left > 0) {
            double tn = now_sec();
            if (tn - t_prog >= 0.2) {
                char p1[32], p2[32];
                fprintf(stderr, "\r\033[K  skipping: %5.1f%%  %s / %s digits ",
                        100.0 * (double)(skip_total - skip_left) / (double)skip_total,
                        fmt_u64(skip_total - skip_left, p1, sizeof p1), fmt_u64(skip_total, p2, sizeof p2));
                fflush(stderr);
                progress_shown = 1;
                t_prog = tn;
            }
        } else if (show_progress) {
            double tn = now_sec();
            if (tn - t_prog >= 0.2) {
                char p1[32], p2[32], p3[32];
                double pct = total_digits != UINT64_MAX && total_digits > 0
                             ? 100.0 * (double)delivered / (double)total_digits : -1;
                double rate = tn - t0 > 0 ? (double)delivered / (tn - t0) : 0;
                double eta = (pct >= 0 && rate > 0) ? (double)(total_digits - delivered) / rate : -1;
                if (pct >= 0)
                    fprintf(stderr, "\r\033[K  %5.1f%%  %s / %s digits  %.0f M/s  found %" PRIu64 "%s%s ",
                            pct, fmt_u64(delivered, p1, sizeof p1), fmt_u64(total_digits, p2, sizeof p2),
                            rate / 1e6, found_total,
                            eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, p3, sizeof p3) : "");
                else
                    fprintf(stderr, "\r\033[K  %s digits searched  %.0f M/s  found %" PRIu64 " ",
                            fmt_u64(delivered, p1, sizeof p1), rate / 1e6, found_total);
                fflush(stderr);
                progress_shown = 1;
                t_prog = tn;
            }
        }
        if (ndigits == 0) break;

        for (int q = 0; q < npats; q++) {
            pat_t *P = &pats[q];
            const size_t L = P->len;
            if (ndigits < L) continue;

            /* searchable starts for this pattern this round */
            size_t search_from = first_chunk ? (size_t)pre_ctx : CTX + (Lmax - L);
            size_t last_start;
            if (eof || interrupted) {
                last_start = ndigits - L;                       /* final round: take everything */
            } else {
                if (ndigits < L + CTX) continue;
                last_start = ndigits - L - CTX;                 /* keep CTX after-context digits */
            }
            if (last_start + 1 <= search_from) continue;

            const uint8_t *region = buf + search_from;
            size_t region_len = (last_start - search_from) + L;   /* matches must START in range */
            const uint8_t *p;
            while (region_len >= L && (p = memmem(region, region_len, P->dig, L)) != NULL) {
                size_t i = (size_t)(p - buf);
                P->found++;
                found_total++;
                if (P->found <= opt.maxprint) {
                    char b1[32];
                    char pre[CTX + 1], post[CTX + 1], mid[MAX_PAT_LEN + 1];
                    uint64_t abspos = base + i;
                    uint64_t avail = abspos - 1;          /* digits of pi before this position */
                    size_t npre = CTX;
                    if (npre > i) npre = i;               /* stay inside the buffer */
                    if ((uint64_t)npre > avail) npre = (size_t)avail;   /* file start */
                    for (size_t k = 0; k < npre; k++) pre[k] = (char)('0' + buf[i - npre + k]);
                    pre[npre] = 0;
                    size_t npost = ndigits - (i + L); if (npost > CTX) npost = CTX;
                    for (size_t k = 0; k < npost; k++) post[k] = (char)('0' + buf[i + L + k]);
                    post[npost] = 0;
                    for (size_t k = 0; k < L; k++) mid[k] = (char)('0' + buf[i + k]);
                    mid[L] = 0;
                    CLEAR_PROGRESS();
                    char lbl[MAX_PAT_LEN + 8] = "";
                    if (npats > 1) snprintf(lbl, sizeof lbl, " of %.*s", (int)P->len, P->str);
                    printf("Match %" PRIu64 "%s at position %s:  %s%s%s%s%s%s%s\n",
                           P->found, lbl, fmt_u64(abspos, b1, sizeof b1),
                           npre  < CTX ? (abspos - npre == 1 ? "^" : "") : "...", pre,
                           ul_on, mid, ul_off,
                           post, npost < CTX ? (eof ? "$" : "") : "...");
                } else if (P->found == opt.maxprint + 1) {
                    CLEAR_PROGRESS();
                    printf("(more matches of %s - printing stopped after %" PRIu64 ", still counting; raise with -m)\n",
                           P->str, opt.maxprint);
                }
                /* advance one digit so overlapping occurrences are found too */
                size_t adv = (size_t)(p - region) + 1;
                region += adv;
                region_len -= adv;
            }
        }

        if (eof || interrupted) break;

        /* carry tail sized by the longest pattern: its unsearched starts begin
         * at ndigits-Lmax-CTX+1; keep CTX before-context digits ahead of that */
        size_t first_unsearched = ndigits - Lmax - CTX + 1;
        size_t carry_from = first_unsearched >= CTX ? first_unsearched - CTX : 0;
        size_t carry = ndigits - carry_from;
        memmove(buf, buf + carry_from, carry);
        base += carry_from;
        ndigits = carry;
        first_chunk = 0;
    }

    CLEAR_PROGRESS();
    double dt = now_sec() - t0;
    char b1[32], b2[32], d1[32];
    if (interrupted) printf("*** interrupted ***\n");
    int missing = 0;
    for (int q = 0; q < npats; q++) {
        printf("%s occurrence%s of %s", fmt_u64(pats[q].found, b1, sizeof b1),
               pats[q].found == 1 ? "" : "s", pats[q].str);
        if (pats[q].found == 0) { printf("  (NOT FOUND)"); missing++; }
        printf("\n");
    }
    printf("Searched %s digits in %s (%.0f M digits/s)\n",
           fmt_u64(delivered, b2, sizeof b2), fmt_duration(dt, d1, sizeof d1),
           dt > 0 ? delivered / dt / 1e6 : 0);
    if (skipped_bytes) printf("(%s non-digit bytes ignored)\n", fmt_u64(skipped_bytes, b1, sizeof b1));
    free(raw); free(buf); free(pathbuf);
    close(fd);
    return missing == 0 ? 0 : 1;
}
