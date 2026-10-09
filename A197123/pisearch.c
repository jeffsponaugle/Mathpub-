/*
 * pisearch - pass 2 for OEIS A197123: exact check of bloom candidates.
 *
 * Loads the candidate list produced by pibloom into an open-addressing hash
 * table, then slides an L-digit window over pi.  The first time a candidate
 * is seen its position is recorded; the second time it is a confirmed repeat.
 * All repeats are reported (not just the first) and the answer for A197123 is
 * the repeat whose SECOND occurrence is earliest.
 *
 * Build:  cc -O3 -march=native -o pisearch pisearch.c -lm
 */
#include "pi_common.h"          /* digit sources, hashing, window keys, filters, timing, alloc helpers */

typedef struct {
    uint64_t lo;
    uint64_t hi;
    uint64_t first;      /* position of first occurrence, 0 = not seen yet */
    uint64_t second;     /* position of second occurrence, 0 = none yet     */
    uint32_t count;      /* occurrences seen                                */
    uint32_t used;       /* slot occupied                                   */
} entry;                 /* 40 bytes */

typedef struct {
    uint64_t lo, hi, first, second;
} match;

static struct {
    const char *pifile, *outfile;
    const char *candfiles[64];   /* one or more candidate files (repeat -c, or comma lists) */
    int ncand;
    uint64_t ndigits;
    int L;
    int fauto;           /* 1 = auto-detect filter from candidates (default) */
    pfilter F;           /* prefix-range filter (len 0 = none) */
    int stop_first;
    double progress_secs;
    size_t chunk;
} opt = { NULL, NULL, {NULL}, 0, UINT64_MAX, 0, 1, {0,0,0,0,0,""}, 0, 30.0, 64u << 20 };

static entry   *tab;
static uint64_t tab_mask;       /* slots - 1 (power of 2) */
static uint64_t tab_used;
static uint64_t *pre;           /* prefilter bitmap */
static uint64_t pre_mask;       /* bits - 1 (power of 2) */
static winspec  W;

static match   *matches;
static uint64_t nmatches, matches_cap;

/* statistics */
static uint64_t st_windows, st_tested, st_pre_pass, st_hits, st_first_seen, st_extra;
static uint64_t st_max_probe, st_cand_lines, st_cand_dupes, st_cand_bad;
static uint64_t st_last_tested;
static uint64_t g_total_digits;
static int g_tty = -1;                 /* stdout is a terminal: status line overwrites itself */
static int g_status_shown = 0;
/*
 * clear_status - erase the in-place status line on a terminal before
 * printing something that should scroll (a repeat announcement).
 */
static void clear_status(void)
{
    if (g_tty > 0 && g_status_shown) { printf("\r\033[K"); g_status_shown = 0; }
}
static double t_start, t_last;

/*
 * usage - print the command-line help to stderr and exit(1).
 */
static void usage(const char *argv0)
{
    fprintf(stderr,
"pisearch v" PI_TOOLS_VERSION " - A197123 pass 2 (exact repeat search)\n"
"Usage: %s -p <pifile> -l <seqlen> -c <candidates> [options]\n"
"\n"
"  -p SRC      pi digits: text or .ycd file(s)/directories, comma-separated\n"
"  -l N        sequence length (must match the candidate file)\n"
"  -c FILES    candidate file(s) from pibloom/pibloomt (lines '<digits>:<pos>',\n"
"              '#' = comment); repeat -c or give a comma list to verify several\n"
"              partitions' candidates in ONE pass over the digits\n"
"  -n COUNT    digits of pi to scan (e.g. 10B, 1T, all)   [all]\n"
"  -f SPEC     only test windows whose 1-3 digit prefix is SPEC: a value (3, 24),\n"
"              a range (0-4, 10-15) or K/N (2/20); 'auto' (default) uses the\n"
"              tightest range covering every candidate; 'none' disables\n"
"  -o FILE     results file  [results_L<len>.txt]\n"
"  -1          stop at the first repeat found (earliest second occurrence)\n"
"  -P SECS     progress interval in seconds [30]\n"
"  -C SIZE     read chunk size [64M]\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0);
    exit(1);
}

/*
 * parse_args - getopt into the global opt struct; -f accepts a prefix
 * value/range, "auto" (detect from the candidates) or "none".
 */
static void parse_args(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:l:c:n:f:o:1P:C:Vh")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'l': opt.L = atoi(optarg); break;
        case 'c': {
            char *dup = strdup(optarg);
            if (!dup) die("out of memory");
            for (char *save = dup, *tok; (tok = strsep(&save, ",")) != NULL; ) {
                if (!*tok) continue;
                if (opt.ncand >= 64) die("too many candidate files (max 64)");
                opt.candfiles[opt.ncand++] = tok;
            }
        } break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 'f':
            if (!strcmp(optarg, "auto")) { opt.fauto = 1; opt.F.len = 0; }
            else if (!strcmp(optarg, "none")) { opt.fauto = 0; opt.F.len = 0; }
            else {
                if (pfilter_parse(optarg, &opt.F)) die("-f must be a prefix value/range (3, 24, 0-4, 10-15), auto or none");
                opt.fauto = 0;
            }
            break;
        case 'o': opt.outfile = optarg; break;
        case '1': opt.stop_first = 1; break;
        case 'V': print_version("pisearch"); exit(0);
        case 'P': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -P"); break;
        case 'C': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile) { fprintf(stderr, "missing -p\n"); usage(argv[0]); }
    if (opt.ncand == 0) { fprintf(stderr, "missing -c\n"); usage(argv[0]); }
    if (opt.L < 1 || opt.L > MAX_SEQ_LEN) { fprintf(stderr, "missing/invalid -l\n"); usage(argv[0]); }
    if (opt.F.len >= opt.L) die("filter prefix (%d digits) must be shorter than the sequence length (%d)", opt.F.len, opt.L);
}

/* ---------------------------------------------------------------- */
/* Hash table                                                         */
/* ---------------------------------------------------------------- */
/*
 * slot_hash - the 64-bit hash of a window key used for both the hash-table
 * slot and the prefilter bit (same hash as pass 1, so behaviour is
 * consistent).
 */
static inline uint64_t slot_hash(uint64_t hi, uint64_t lo, uint64_t *h2)
{
    uint64_t h1;
    hash_window(lo, hi, &h1, h2);
    return h1;
}

/*
 * tab_insert - insert a candidate key into the open-addressing table
 * (linear probing, <= 50% load) and set its prefilter bit.  Returns 0 if
 * the key was already present (duplicate candidate line), 1 if new.
 * Tracks the longest probe sequence for the load report.
 */
static int tab_insert(uint64_t hi, uint64_t lo)
{
    uint64_t h2;
    uint64_t h = slot_hash(hi, lo, &h2);
    uint64_t i = h & tab_mask;
    uint64_t probes = 1;
    while (tab[i].used) {
        if (tab[i].lo == lo && tab[i].hi == hi) return 0;
        i = (i + 1) & tab_mask;
        probes++;
    }
    tab[i].used = 1; tab[i].lo = lo; tab[i].hi = hi;
    tab[i].first = tab[i].second = 0; tab[i].count = 0;
    tab_used++;
    if (probes > st_max_probe) st_max_probe = probes;
    pre[(h >> 6) & (pre_mask >> 6)] |= 1ULL << (h & 63);
    return 1;
}

/*
 * tab_find - look up a key; NULL if absent.  Linear probing from h.
 */
static inline entry *tab_find(uint64_t hi, uint64_t lo, uint64_t h)
{
    uint64_t i = h & tab_mask;
    while (tab[i].used) {
        if (tab[i].lo == lo && tab[i].hi == hi) return &tab[i];
        i = (i + 1) & tab_mask;
    }
    return NULL;
}

/* ---------------------------------------------------------------- */
/*
 * load_candidates - read the candidate file in two passes: count lines
 * (sizing the table and prefilter), then parse "digits:position" lines,
 * validate them against -l, dedupe and insert.  Reports the tightest
 * prefix range covering every candidate (for -f auto).  Shows progress
 * for both passes; tolerates comments and a small number of bad lines.
 */
static int load_candidates(pfilter *common)   /* tightest prefix range covering all candidates */
{
    /* pass 1: count lines over all files to size the table */
    uint64_t cf_size = 0;
    for (int k = 0; k < opt.ncand; k++) {
        struct stat cst;
        if (stat(opt.candfiles[k], &cst) != 0) die("cannot open candidate file '%s': %s", opt.candfiles[k], strerror(errno));
        if (S_ISREG(cst.st_mode)) cf_size += (uint64_t)cst.st_size;
    }
    {
        char bb[32];
        printf("Reading %d candidate file%s (%s) ...\n", opt.ncand, opt.ncand == 1 ? "" : "s",
               cf_size ? fmt_bytes((double)cf_size, bb, sizeof bb) : "size unknown");
    }
    int tty = isatty(1);
    uint64_t lines = 0;
    {
        char *buf = malloc(1 << 22);
        if (!buf) die("out of memory");
        size_t n;
        uint64_t seen = 0;
        int count_decile = 0;
        double t0 = now_sec(), tlast = t0;
        for (int k = 0; k < opt.ncand; k++) {
        FILE *f = fopen(opt.candfiles[k], "r");
        if (!f) die("cannot open candidate file '%s': %s", opt.candfiles[k], strerror(errno));
        while ((n = fread(buf, 1, 1 << 22, f)) > 0) {
            const char *p = buf, *e = buf + n;
            while ((p = memchr(p, '\n', (size_t)(e - p))) != NULL) { lines++; p++; }
            seen += n;
            double t = now_sec();
            if (t - t0 >= 1.0) {
                char c1[32], c2[32];
                if (tty && t - tlast >= 0.5) {
                    if (cf_size)
                        printf("\r\033[K  counting: %5.1f%%  %s lines so far ", 100.0 * (double)seen / (double)cf_size,
                               fmt_u64(lines, c1, sizeof c1));
                    else
                        printf("\r\033[K  counting: %s read, %s lines so far ", fmt_bytes((double)seen, c2, sizeof c2),
                               fmt_u64(lines, c1, sizeof c1));
                    fflush(stdout);
                    tlast = t;
                } else if (!tty && cf_size) {
                    int dec = (int)(100.0 * (double)seen / (double)cf_size / 10.0);
                    if (dec > count_decile) {
                        printf("  counting: %3d%%  %s lines so far\n", dec * 10, fmt_u64(lines, c1, sizeof c1));
                        count_decile = dec;
                    }
                }
            }
        }
        if (ferror(f)) die("read error on candidate file '%s'", opt.candfiles[k]);
        fclose(f);
        lines++; /* possibly unterminated last line */
        }
        if (tty) { printf("\r\033[K"); fflush(stdout); }
        free(buf);
    }
    uint64_t slots = 1024;
    while (slots < lines * 2 + 16) slots <<= 1;       /* load factor <= 50% */
    tab_mask = slots - 1;
    uint64_t pre_bits = 1ULL << 20;
    while (pre_bits < lines * 64 && pre_bits < (1ULL << 36)) pre_bits <<= 1;
    pre_mask = pre_bits - 1;

    char b1[32], b2[32];
    printf("  candidate lines : ~%s\n", fmt_u64(lines, b1, sizeof b1));
    printf("  hash table      : %s slots x %zu B = %s\n", fmt_u64(slots, b1, sizeof b1), sizeof(entry),
           fmt_bytes((double)slots * sizeof(entry), b2, sizeof b2));
    printf("  prefilter bitmap: %s bits = %s\n", fmt_u64(pre_bits, b1, sizeof b1), fmt_bytes((double)pre_bits / 8, b2, sizeof b2));
    check_fits_in_ram(slots * sizeof(entry) + pre_bits / 8, "candidate hash table + prefilter");
    tab = alloc_huge(slots * sizeof(entry), 1, "candidate hash table");
    pre = alloc_huge(pre_bits / 8, 1, "prefilter bitmap");

    /* pass 2: parse and insert, file by file */
    printf("Loading candidates into hash table ...\n");
    char line[1024];
    uint64_t lineno = 0;
    double lp_t0 = now_sec(), lp_tlast = lp_t0;
    int load_decile = 0;
    const int plen = opt.L >= 4 ? 3 : opt.L >= 3 ? 2 : 1;   /* prefix must stay shorter than L */
    int pmin = 1000, pmax = -1;
    for (int k = 0; k < opt.ncand; k++) {
    FILE *f = fopen(opt.candfiles[k], "r");
    if (!f) die("cannot open candidate file '%s': %s", opt.candfiles[k], strerror(errno));
    uint64_t file_lineno = 0, used_before = tab_used;
    while (fgets(line, sizeof line, f)) {
        lineno++;
        file_lineno++;
        if ((lineno & 0xFFFFF) == 0) {
            double t = now_sec();
            char c1[32], c2[32];
            double pct = 100.0 * (double)lineno / (double)lines;
            if (tty && t - lp_tlast >= 0.5) {
                printf("\r\033[K  loading: %5.1f%%  %s of ~%s lines  %.1f M lines/s ",
                       pct, fmt_u64(lineno, c1, sizeof c1),
                       fmt_u64(lines, c2, sizeof c2), (double)lineno / (t - lp_t0) / 1e6);
                fflush(stdout);
                lp_tlast = t;
            } else if (!tty && t - lp_t0 >= 1.0) {
                int dec = (int)(pct / 10.0);
                if (dec > load_decile) {
                    printf("  loading: %3d%%  %s of ~%s lines  %.1f M lines/s\n",
                           dec * 10, fmt_u64(lineno, c1, sizeof c1), fmt_u64(lines, c2, sizeof c2),
                           (double)lineno / (t - lp_t0) / 1e6);
                    load_decile = dec;
                }
            }
        }
        size_t len = strlen(line);
        if (len && line[len - 1] == '\n') line[--len] = 0;
        if (len && line[len - 1] == '\r') line[--len] = 0;
        if (len == 0 || line[0] == '#') continue;
        if (len == sizeof line - 1) die("'%s' line %" PRIu64 " too long", opt.candfiles[k], file_lineno);
        st_cand_lines++;
        /* find digit run */
        size_t nd = 0;
        while (nd < len && line[nd] >= '0' && line[nd] <= '9') nd++;
        if (nd != (size_t)opt.L || (nd < len && line[nd] != ':')) {
            st_cand_bad++;
            if (st_cand_bad <= 5) warn("'%s' line %" PRIu64 " ignored (expected %d digits then ':'): \"%s\"", opt.candfiles[k], file_lineno, opt.L, line);
            if (st_cand_bad == 6) warn("(further bad-line warnings suppressed)");
            continue;
        }
        uint64_t hi, lo;
        if (win_parse_ascii(&W, line, &hi, &lo) != 0) {   /* cannot happen: validated above */
            st_cand_bad++;
            continue;
        }
        int pv = 0;
        for (int q = 0; q < plen; q++) pv = pv * 10 + (line[q] - '0');
        if (pv < pmin) pmin = pv;
        if (pv > pmax) pmax = pv;
        if (!tab_insert(hi, lo)) st_cand_dupes++;
    }
    if (ferror(f)) die("read error on candidate file '%s'", opt.candfiles[k]);
    fclose(f);
    if (opt.ncand > 1) {
        char c1[32];
        if (tty) printf("\r\033[K");
        printf("  %-40s %s new candidates\n", opt.candfiles[k], fmt_u64(tab_used - used_before, c1, sizeof c1));
    }
    }
    if (tty) { printf("\r\033[K"); fflush(stdout); }
    if (st_cand_bad && st_cand_bad * 100 > st_cand_lines) die("%" PRIu64 " of %" PRIu64 " candidate lines are malformed - wrong -l or wrong file?", st_cand_bad, st_cand_lines);
    if (pmax < 0) { common->len = 0; common->str[0] = 0; }
    else {
        char spec[16];
        snprintf(spec, sizeof spec, "%0*d-%0*d", plen, pmin, plen, pmax);
        pfilter_parse(spec, common);
        /* shrink to the shortest prefix length that expresses the same range */
        while (common->len > 1 && common->lo % 10 == 0 && common->hi % 10 == 9)
            pfilter_set(common, common->len - 1, common->lo / 10, common->hi / 10);
        /* a range covering every possible prefix is no filter at all */
        if (common->lo == 0 && common->hi == (int)pow(10.0, plen) - 1) common->len = 0;
    }
    char b3[32];
    printf("  loaded          : %s unique candidates (%s duplicates, %s bad lines), load %.1f%%, max probe %" PRIu64 "\n",
           fmt_u64(tab_used, b1, sizeof b1), fmt_u64(st_cand_dupes, b2, sizeof b2),
           fmt_u64(st_cand_bad, b3, sizeof b3), 100.0 * (double)tab_used / (double)slots, st_max_probe);
    return tab_used > 0;
}

/* ---------------------------------------------------------------- */
static FILE *outf;

/*
 * record_match - a candidate has just been seen for the second time:
 * store it in the match list, announce it (first 20 on the console), and
 * append "digits:first:second" to the results file.
 */
static void record_match(entry *e)
{
    if (nmatches == matches_cap) {
        matches_cap = matches_cap ? matches_cap * 2 : 1024;
        matches = realloc(matches, matches_cap * sizeof *matches);
        if (!matches) die("out of memory recording matches");
    }
    match *m = &matches[nmatches++];
    m->lo = e->lo; m->hi = e->hi; m->first = e->first; m->second = e->second;
    char s[MAX_SEQ_LEN + 1], b1[32], b2[32];
    win_to_ascii(&W, e->hi, e->lo, s);
    clear_status();
    if (nmatches <= 20)
        printf("*** REPEAT #%" PRIu64 ": %s  at positions %s and %s\n", nmatches, s,
               fmt_u64(e->first, b1, sizeof b1), fmt_u64(e->second, b2, sizeof b2));
    else if (nmatches == 21)
        printf("*** (further repeats are written to %s only; the earliest second occurrence is reported at the end)\n", opt.outfile);
    if (fprintf(outf, "%s:%" PRIu64 ":%" PRIu64 "\n", s, e->first, e->second) < 0)
        die("write error on results file: %s", strerror(errno));
    if (nmatches <= 100) fflush(outf);
}

/*
 * progress - status line with position, % done, tested windows, how many
 * candidates have been seen once, repeats found, rates and ETA.  Overwrites
 * itself on a terminal, logs a line per interval otherwise.
 */
static void progress(int force, uint64_t pos)
{
    double t = now_sec();
    if (!force && t - t_last < opt.progress_secs) return;
    double dt = t - t_last, tot = t - t_start;
    char b1[32], b2[32], b3[32], b4[32], d1[32], d2[32];
    double rate_sec = dt > 0 ? (double)(st_tested - st_last_tested) / dt : 0;
    double win_rate = tot > 0 ? (double)st_windows / tot : 0;
    double eta = -1, pct = -1;
    if (g_total_digits) {
        pct = 100.0 * (double)st_windows / (double)g_total_digits;
        if (pct > 100.0) pct = 100.0;
        if (win_rate > 0 && g_total_digits > st_windows) eta = (double)(g_total_digits - st_windows) / win_rate;
    }
    char pb[24] = "";
    if (pct >= 0) snprintf(pb, sizeof pb, " (%.1f%% done)", pct);
    if (g_tty < 0) g_tty = isatty(1);
    printf("%s[%s] pos %s%s  tested %s  seen-once %s/%s  repeats %" PRIu64 "  %.1f M/s (avg %.1f M win/s)%s%s%s",
           g_tty ? "\r\033[K" : "",
           fmt_duration(tot, d1, sizeof d1), fmt_u64(pos, b1, sizeof b1), pb, fmt_u64(st_tested, b2, sizeof b2),
           fmt_u64(st_first_seen, b3, sizeof b3), fmt_u64(tab_used, b4, sizeof b4), nmatches,
           rate_sec / 1e6, win_rate / 1e6,
           eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, d2, sizeof d2) : "",
           g_tty && !force ? "" : "\n");
    g_status_shown = g_tty && !force;
    fflush(stdout);
    fflush(outf);
    t_last = t; st_last_tested = st_tested;
}

/* ---------------------------------------------------------------- */
/* Scan pipeline: stage 1 hashes a window and prefetches its prefilter
 * word; PPD windows later stage 2 tests the (now cached) bit and, on a
 * hit, prefetches the hash-table slot into a second FIFO probed SPD
 * entries later.  Both queues are strict FIFO, so candidate first/second
 * occurrences are still seen in exact scan order. */
#define PPD 64                     /* primary depth (power of 2) */
#define SPD 16                     /* secondary depth */
typedef struct { uint64_t h, hi, lo, pos; } pent;

static pent  pring[PPD];
static unsigned phead, pcount;
static pent  sring[SPD];
static unsigned shead, stail, scount;
static int   g_stop_scan;          /* set when a -1 stop is requested */
static uint64_t g_pre_wmask;

/*
 * probe_entry - stage 3 of the scan pipeline: look the window up in the
 * table (its slot was prefetched in stage 2) and update the occurrence
 * count: first sighting records the position, second is a confirmed
 * repeat, later ones are counted as extras.  Honours -1 (stop first).
 */
static void probe_entry(const pent *pe)
{
    entry *e = tab_find(pe->hi, pe->lo, pe->h);
    if (e) {
        st_hits++;
        e->count++;
        if (e->count == 1) { e->first = pe->pos; st_first_seen++; }
        else if (e->count == 2) {
            e->second = pe->pos;
            record_match(e);
            if (opt.stop_first) g_stop_scan = 1;
        } else st_extra++;
    }
}

/*
 * stage2 - test the prefilter bit (prefetched PPD windows earlier); on a
 * hit prefetch the table slot and queue the window for probe_entry SPD
 * entries later.  Both queues are strict FIFO so occurrences are still
 * processed in exact scan order.
 */
static inline void stage2(const pent *pe)
{
    if (pre[(pe->h >> 6) & g_pre_wmask] & (1ULL << (pe->h & 63))) {
        st_pre_pass++;
        __builtin_prefetch(&tab[pe->h & tab_mask], 0, 0);
        if (scount == SPD) { probe_entry(&sring[stail]); stail = (stail + 1) % SPD; scount--; }
        sring[shead] = *pe;
        shead = (shead + 1) % SPD;
        scount++;
    }
}

/*
 * drain_pipeline - flush both pipeline queues at the end of a chunk (and
 * at EOF), preserving order.
 */
static void drain_pipeline(void)
{
    unsigned tail = (phead + PPD - pcount) & (PPD - 1);
    while (pcount) { stage2(&pring[tail]); tail = (tail + 1) & (PPD - 1); pcount--; }
    while (scount) { probe_entry(&sring[stail]); stail = (stail + 1) % SPD; scount--; }
}

/*
 * cmp_second - qsort comparator: matches by second occurrence, ascending.
 * The first element after sorting is a(n) by the OEIS definition.
 */
static int cmp_second(const void *a, const void *b)
{
    const match *x = a, *y = b;
    return x->second < y->second ? -1 : x->second > y->second;
}
/*
 * cmp_value - qsort comparator: matches by numeric value, ascending.
 */
static int cmp_value(const void *a, const void *b)
{
    const match *x = a, *y = b;
    if (x->hi != y->hi) return x->hi < y->hi ? -1 : 1;
    return x->lo < y->lo ? -1 : x->lo > y->lo;
}

/*
 * main - pass 2 driver: load candidates, resolve the filter, then scan the
 * digits with the three-stage prefetch pipeline (hash -> prefilter bit ->
 * table probe).  Reports every repeat, names the one with the earliest
 * second occurrence (that is a(n)), and warns if candidates were never
 * seen (usually a shorter -n than pass 1 used).
 */
int main(int argc, char **argv)
{
    parse_args(argc, argv);
    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();
    winspec_init(&W, opt.L);

    char outname[256];
    if (!opt.outfile) { snprintf(outname, sizeof outname, "results_L%d.txt", opt.L); opt.outfile = outname; }

    char b1[64], b2[64], b3[64], d1[32];
    printf("pisearch v%s - A197123 pass 2 (exact repeat search)\n", PI_TOOLS_VERSION);
    printf("  pi source      : %s\n", opt.pifile);
    printf("  sequence length: %d\n", opt.L);
    for (int k = 0; k < opt.ncand; k++) printf("  candidates     : %s\n", opt.candfiles[k]);
    printf("  digits to scan : %s\n", opt.ndigits == UINT64_MAX ? "all" : fmt_u64(opt.ndigits, b1, sizeof b1));
    printf("  results file   : %s\n", opt.outfile);

    pfilter common;
    if (!load_candidates(&common)) die("no usable candidates loaded");
    if (opt.fauto) {
        opt.F = common;
        opt.F.frac_n = 0;
        if (opt.F.len) printf("  filter         : auto - every candidate's %d-digit prefix is in [%s] (%.3g%% of windows), testing only those\n",
                              opt.F.len, opt.F.str, 100.0 * pfilter_fraction(&opt.F));
        else printf("  filter         : none (candidate prefixes span the full range)\n");
    } else if (opt.F.len > 0) {
        /* compare the explicit filter with the candidates' range at the explicit length */
        if (common.len) {
            /* bring both to 3-digit resolution for the comparison */
            int sc = (int)pow(10.0, 3 - common.len), sf = (int)pow(10.0, 3 - opt.F.len);
            int clo = common.lo * sc, chi = common.hi * sc + (sc - 1);
            int flo = opt.F.lo * sf, fhi = opt.F.hi * sf + (sf - 1);
            if (fhi < clo || flo > chi)
                warn("-f %s excludes every candidate (their prefixes are in [%s]) - nothing can match!", opt.F.str, common.str);
            else if (flo > clo || fhi < chi)
                warn("-f %s does not cover all candidates (prefixes in [%s]) - some will be missed", opt.F.str, common.str);
        }
        char fd[32];
        printf("  filter         : %d-digit prefix in %s (%.3g%% of windows)\n", opt.F.len,
               pfilter_describe(&opt.F, fd, sizeof fd), 100.0 * pfilter_fraction(&opt.F));
    } else printf("  filter         : none\n");

    outf = fopen(opt.outfile, "w");
    if (!outf) die("cannot create results file '%s': %s", opt.outfile, strerror(errno));
    fprintf(outf, "# pisearch v%s L=%d pi=%s  format digits:first:second (1-based after decimal point)\n", PI_TOOLS_VERSION, opt.L, opt.pifile);
    for (int k = 0; k < opt.ncand; k++) fprintf(outf, "# candidates=%s\n", opt.candfiles[k]);
    fflush(outf);

    pireader rd;
    pireader_open(&rd, opt.pifile, opt.L, opt.chunk, opt.ndigits);
    {
        char desc[256];
        pistream_describe(&rd.ps, desc, sizeof desc);
        printf("  source format  : %s\n", desc);
        uint64_t fd_digits = rd.digits_est;
        g_total_digits = opt.ndigits;
        if (fd_digits && (g_total_digits == UINT64_MAX || g_total_digits > fd_digits)) g_total_digits = fd_digits;
        if (g_total_digits == UINT64_MAX) g_total_digits = 0;
    }

    printf("Scanning ...\n");
    t_start = t_last = now_sec();
    const pfilter F = opt.F;
    const int flen = F.len;
    int interrupted = 0, stopped = 0;
    uint64_t last_pos = 0;
    const uint64_t pre_wmask = pre_mask >> 6;
    g_pre_wmask = pre_wmask;

    size_t nwin;
    while ((nwin = pireader_fill(&rd)) > 0) {
        const uint8_t *d = rd.buf;
        uint64_t hi, lo;
        win_compute(&W, d, &hi, &lo);
        uint64_t pos = rd.base;
        for (size_t i = 0; i < nwin; i++, pos++) {
            if (flen == 0 || pfilter_match(&F, d + i)) {
                st_tested++;
                uint64_t h2;
                uint64_t h = slot_hash(hi, lo, &h2);
                pent *e = &pring[phead];
                if (pcount == PPD) stage2(e); else pcount++;
                e->h = h; e->hi = hi; e->lo = lo; e->pos = pos;
                __builtin_prefetch(&pre[(h >> 6) & pre_wmask], 0, 0);
                phead = (phead + 1) & (PPD - 1);
            }
            if (i + 1 < nwin)                 /* no roll past the last window */
                win_roll(&W, d + i, &hi, &lo);
        }
        drain_pipeline();
        if (g_stop_scan && opt.stop_first) { stopped = 1; last_pos = pos - 1; goto done; }
        st_windows += nwin;
        last_pos = pos - 1;
        progress(0, last_pos);
        if (g_stop) { interrupted = 1; break; }
    }
done:
    progress(1, last_pos);
    double elapsed = now_sec() - t_start;

    printf("\n%s\n", interrupted ? "*** INTERRUPTED - partial results ***" : stopped ? "Stopped at first repeat" : "Scan complete");
    printf("  windows scanned : %s (positions %s..%s)\n", fmt_u64(st_windows, b1, sizeof b1),
           fmt_u64(rd.first_pos, b3, sizeof b3), fmt_u64(last_pos, b2, sizeof b2));
    printf("  windows tested  : %s\n", fmt_u64(st_tested, b1, sizeof b1));
    printf("  prefilter passed: %s (%.4f%%)\n", fmt_u64(st_pre_pass, b1, sizeof b1), st_tested ? 100.0 * st_pre_pass / st_tested : 0);
    printf("  table hits      : %s (candidates seen once: %s of %s, extra occurrences: %s)\n",
           fmt_u64(st_hits, b1, sizeof b1), fmt_u64(st_first_seen, b2, sizeof b2), fmt_u64(tab_used, d1, sizeof d1), fmt_u64(st_extra, b3, sizeof b3));
    printf("  repeats found   : %" PRIu64 "\n", nmatches);
    printf("  elapsed         : %s  (%.2f M windows/s)\n", fmt_duration(elapsed, d1, sizeof d1), elapsed > 0 ? st_windows / elapsed / 1e6 : 0);
    printf("  peak RSS        : %.0f MiB\n", peak_rss_mb());
    if (!interrupted && !stopped && st_first_seen < tab_used)
        warn("%" PRIu64 " candidates were never seen in the scanned range - was pisearch given fewer digits than pibloom?", tab_used - st_first_seen);

    if (nmatches) {
        char s[MAX_SEQ_LEN + 1];
        qsort(matches, nmatches, sizeof *matches, cmp_second);
        win_to_ascii(&W, matches[0].hi, matches[0].lo, s);
        printf("\nA197123 a(%d) candidate (earliest second occurrence):\n  %s  at positions %s and %s%s\n", opt.L, s,
               fmt_u64(matches[0].first, b1, sizeof b1), fmt_u64(matches[0].second, b2, sizeof b2),
               s[0] == '0' ? "   (note leading zero)" : "");
        if (nmatches > 1) {
            printf("\nAll %" PRIu64 " repeats by second occurrence:\n", nmatches);
            for (uint64_t i = 0; i < nmatches && i < 50; i++) {
                win_to_ascii(&W, matches[i].hi, matches[i].lo, s);
                printf("  %s  %s  %s\n", s, fmt_u64(matches[i].first, b1, sizeof b1), fmt_u64(matches[i].second, b2, sizeof b2));
            }
            if (nmatches > 50) printf("  ... (%" PRIu64 " more in %s)\n", nmatches - 50, opt.outfile);
            qsort(matches, nmatches, sizeof *matches, cmp_value);
            win_to_ascii(&W, matches[0].hi, matches[0].lo, s);
            printf("Numerically smallest repeat: %s (%s, %s)\n", s, fmt_u64(matches[0].first, b1, sizeof b1), fmt_u64(matches[0].second, b2, sizeof b2));
        }
        fprintf(outf, "# summary: %" PRIu64 " repeats; earliest-second-occurrence = ", nmatches);
        qsort(matches, nmatches, sizeof *matches, cmp_second);
        win_to_ascii(&W, matches[0].hi, matches[0].lo, s);
        fprintf(outf, "%s:%" PRIu64 ":%" PRIu64 "%s\n", s, matches[0].first, matches[0].second, interrupted ? " (INTERRUPTED - may not be final)" : "");
    } else {
        printf("\nNo repeats found in the scanned range.\n");
        fprintf(outf, "# no repeats found%s\n", interrupted ? " (INTERRUPTED)" : "");
    }
    if (fclose(outf)) die("close error on results file");
    pireader_close(&rd);
    return interrupted ? 130 : 0;
}
