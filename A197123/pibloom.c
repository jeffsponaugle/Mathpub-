/*
 * pibloom - pass 1 for OEIS A197123 (first n-digit substring to repeat in pi)
 *
 * Reads the digits of pi (text or y-cruncher .ycd, see pi_common.h) and
 * slides an L-digit window over them, inserting every
 * window into a bloom filter (k=4).  A window whose 4 bits are all already
 * set is a *candidate* for a repeat and is written to the candidate file as
 *      <digits>:<position>
 * where position is 1-based from the first digit after the decimal point.
 * The candidate list is a superset of the true repeats (bloom filters never
 * give false negatives); pisearch then does the exact check.
 *
 * Build:  cc -O3 -march=native -o pibloom pibloom.c -lm
 */
#include "pi_common.h"          /* digit sources, hashing, window keys, filters, timing, alloc helpers */

#define K_HASHES 4
#define PIPE_DEPTH 32            /* software prefetch distance (windows) */
#define BLOCK_BITS 512           /* blocked mode: one cache line */

typedef struct {
    uint64_t idx[K_HASHES];      /* bit indices (classic) or block*512+bit (blocked) */
    uint64_t pos;                /* absolute 1-based position */
    uint64_t hi, lo;             /* window value (for output) */
} pipe_ent;

static struct {
    const char *pifile;
    const char *outfile;
    uint64_t ndigits;            /* digits to scan (UINT64_MAX = all) */
    int L;
    pfilter F;                   /* prefix-range filter (len 0 = none) */
    uint64_t bloom_bytes;
    int blocked;
    int prefault;
    double progress_secs;
    size_t chunk;
} opt = { NULL, NULL, UINT64_MAX, 0, {0,0,0,0,0,""}, 0, 0, 1, 30.0, 64u << 20 };

static uint8_t *bloom;
static uint64_t bloom_bits;      /* total bits */
static uint64_t bloom_blocks;    /* blocked mode */

/* statistics */
static uint64_t st_windows;      /* windows scanned (before filter) */
static uint64_t st_tested;       /* windows inserted / tested */
static uint64_t st_candidates;   /* windows that passed (all 4 bits set) */
static uint64_t st_bits_set;     /* exact count of bits flipped 0->1 */
static uint64_t st_last_tested, st_last_cand, st_last_windows, st_last_raw;
static pireader *g_rd;           /* the scan's reader, for MB/s in the status line */
static uint64_t g_total_digits;   /* effective digits to scan, for %% done / ETA */
static double   t_start, t_last;

/*
 * usage - print the command-line help to stderr and exit(1).
 */
static void usage(const char *argv0)
{
    fprintf(stderr,
"pibloom v" PI_TOOLS_VERSION " - A197123 pass 1 (bloom candidate finder)\n"
"Usage: %s -p <pifile> -l <seqlen> -b <bloomsize> [options]\n"
"\n"
"  -p SRC      pi digits: text file(s) (\"3.14159...\") or y-cruncher .ycd\n"
"              file(s)/directories, comma-separated, read as one stream\n"
"  -l N        sequence length to search for (1..%d)\n"
"  -b SIZE     bloom filter size in bytes: 1M .. 16T (e.g. 512M, 64G, 3.5T)\n"
"  -n COUNT    digits of pi to scan (e.g. 500M, 10B, 1T, all)   [all]\n"
"  -f SPEC     only test windows whose 1-3 digit prefix is SPEC: a value (3, 24),\n"
"              an inclusive range (0-4, 10-15, 05-09), or K/N = the K-th of N\n"
"              equal partitions (2/20 = [05-09]); partition work at any granularity\n"
"  -o FILE     candidate output file  [candidates_L<len>_N<count>[_F<d>].txt]\n"
"  -B          blocked bloom: all 4 bits inside one 64-byte cache line\n"
"              (~3-4x faster on huge tables, slightly higher false-positive rate)\n"
"  -P SECS     progress report interval in seconds [30]\n"
"  -c SIZE     read chunk size [64M]\n"
"  -N          do not pre-fault (memset) the bloom memory up front\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0, MAX_SEQ_LEN);
    exit(1);
}

/*
 * parse_args - getopt over the command line into the global opt struct,
 * validating ranges (bloom size 1M..16T, sequence length, filter prefix
 * shorter than the sequence).  Dies with a message on any problem.
 */
static void parse_args(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:l:b:n:f:o:BP:c:NVh")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'l': opt.L = atoi(optarg); break;
        case 'b': if (parse_mem_size(optarg, &opt.bloom_bytes)) die("bad bloom size '%s'", optarg); break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 'f':
            if (pfilter_parse(optarg, &opt.F)) die("bad filter '%s' (use 3, 24, 0-4, 10-15, 2/20 ...)", optarg);
            break;
        case 'o': opt.outfile = optarg; break;
        case 'B': opt.blocked = 1; break;
        case 'P': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -P"); break;
        case 'c': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        case 'N': opt.prefault = 0; break;
        case 'V': print_version("pibloom"); exit(0);
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile) { fprintf(stderr, "missing -p <pifile>\n"); usage(argv[0]); }
    if (opt.L < 1 || opt.L > MAX_SEQ_LEN) { fprintf(stderr, "missing/invalid -l (1..%d)\n", MAX_SEQ_LEN); usage(argv[0]); }
    if (opt.bloom_bytes < (1ULL << 20) || opt.bloom_bytes > (16ULL << 40)) {
        fprintf(stderr, "bloom size must be between 1M and 16T bytes\n"); usage(argv[0]);
    }
    if (opt.ndigits == 0) die("-n must be > 0");
    if (opt.F.len >= opt.L) die("filter prefix (%d digits) must be shorter than the sequence length (%d)", opt.F.len, opt.L);
}

/* ---------------------------------------------------------------- */
/*
 * compute_idx - turn a window key into its K_HASHES bloom bit indices.
 * Classic mode: 4 independent positions over the whole table (double
 * hashing + fastrange).  Blocked mode (-B): one 512-bit block chosen by
 * h1, four bits inside it from h2 - a single cache line per window.
 */
static inline void compute_idx(uint64_t hi, uint64_t lo, uint64_t *idx)
{
    uint64_t h1, h2;
    hash_window(lo, hi, &h1, &h2);
    if (opt.blocked) {
        uint64_t blk = fastrange64(h1, bloom_blocks) * BLOCK_BITS;
        idx[0] = blk + ( h2        & 511);
        idx[1] = blk + ((h2 >> 9)  & 511);
        idx[2] = blk + ((h2 >> 18) & 511);
        idx[3] = blk + ((h2 >> 27) & 511);
    } else {
        idx[0] = fastrange64(h1,           bloom_bits);
        idx[1] = fastrange64(h1 + h2,      bloom_bits);
        idx[2] = fastrange64(h1 + 2 * h2,  bloom_bits);
        idx[3] = fastrange64(h1 + 3 * h2,  bloom_bits);
    }
}

/*
 * prefetch_idx - issue prefetches for the cache lines holding the bits, so
 * the DRAM misses of the next PIPE_DEPTH windows overlap instead of
 * serialising.
 */
static inline void prefetch_idx(const uint64_t *idx)
{
    __builtin_prefetch(&bloom[idx[0] >> 3], 1, 0);
    if (!opt.blocked) {
        __builtin_prefetch(&bloom[idx[1] >> 3], 1, 0);
        __builtin_prefetch(&bloom[idx[2] >> 3], 1, 0);
        __builtin_prefetch(&bloom[idx[3] >> 3], 1, 0);
    }
}

static FILE *outf;
static winspec W;

/*
 * process_ent - test-and-set the window's bits (prefetched PIPE_DEPTH
 * windows ago).  If all were already set the window is a candidate and is
 * written as "digits:position"; otherwise the missing bits are set and
 * the exact count of newly set bits is tracked for the fill statistics.
 */
static inline void process_ent(const pipe_ent *e)
{
    unsigned all = 1;
    for (int i = 0; i < K_HASHES; i++) {
        uint64_t b = e->idx[i];
        uint8_t m = (uint8_t)(1u << (b & 7));
        uint8_t *p = &bloom[b >> 3];
        if (!(*p & m)) { all = 0; *p |= m; st_bits_set++; }
    }
    if (all) {
        char s[MAX_SEQ_LEN + 1];
        win_to_ascii(&W, e->hi, e->lo, s);
        if (fprintf(outf, "%s:%" PRIu64 "\n", s, e->pos) < 0)
            die("write error on candidate file: %s", strerror(errno));
        st_candidates++;
    }
}

/*
 * progress - print the status line (position, % done, tested, candidates,
 * fill, rates, ETA) if the progress interval has elapsed or force is set;
 * also flushes the candidate file so an interrupted run keeps its output.
 */
static void progress(int force, uint64_t pos)
{
    double t = now_sec();
    if (!force && t - t_last < opt.progress_secs) return;
    double dt = t - t_last, tot = t - t_start;
    char b1[32], b2[32], b3[32], b4[32], d1[32], d2[32];
    double fill = (double)st_bits_set / (double)bloom_bits;
    double rate_sec = dt > 0 ? (double)(st_tested - st_last_tested) / dt : 0;
    double rate_all = tot > 0 ? (double)st_tested / tot : 0;
    double win_rate = tot > 0 ? (double)st_windows / tot : 0;
    double eta = -1, pct = -1;
    if (g_total_digits) {
        pct = 100.0 * (double)st_windows / (double)g_total_digits;
        if (pct > 100.0) pct = 100.0;
        if (win_rate > 0 && g_total_digits > st_windows) eta = (double)(g_total_digits - st_windows) / win_rate;
    }
    char pb[24] = "";
    if (pct >= 0) snprintf(pb, sizeof pb, " (%.1f%% done)", pct);
    static int tty = -1;
    if (tty < 0) tty = isatty(1);       /* status line overwrites itself on a terminal */
    double dig_rate = dt > 0 ? (double)(st_windows - st_last_windows) / dt : 0;
    double mb_rate  = (dt > 0 && g_rd) ? (double)(g_rd->ps.total_raw - st_last_raw) / dt / 1e6 : 0;
    printf("%s[%s] pos %s%s  tested %s  cand %s (+%s)  fill %.4f%%  %.1f M/s (avg %.1f M/s)  %.0f M digits/s (%.0f MB/s)%s%s%s",
           tty ? "\r\033[K" : "",
           fmt_duration(tot, d1, sizeof d1), fmt_u64(pos, b1, sizeof b1), pb,
           fmt_u64(st_tested, b2, sizeof b2), fmt_u64(st_candidates, b3, sizeof b3),
           fmt_u64(st_candidates - st_last_cand, b4, sizeof b4),
           fill * 100.0, rate_sec / 1e6, rate_all / 1e6, dig_rate / 1e6, mb_rate,
           eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, d2, sizeof d2) : "",
           tty && !force ? "" : "\n");
    fflush(stdout);
    fflush(outf);
    t_last = t; st_last_tested = st_tested; st_last_cand = st_candidates;
    st_last_windows = st_windows; if (g_rd) st_last_raw = g_rd->ps.total_raw;
}

/*
 * popcount_bloom - count every set bit in the table (64 bits at a time) to
 * cross-check the incrementally tracked bit count at the end of a run.
 */
static uint64_t popcount_bloom(void)
{
    uint64_t sum = 0;
    const uint64_t *w = (const uint64_t *)bloom;
    uint64_t nw = opt.bloom_bytes / 8;
    for (uint64_t i = 0; i < nw; i++) sum += (uint64_t)__builtin_popcountll(w[i]);
    for (uint64_t i = nw * 8; i < opt.bloom_bytes; i++) sum += (uint64_t)__builtin_popcount(bloom[i]);
    return sum;
}

/*
 * main - pass 1 driver: parse options, print the banner with expected
 * fill/false-positive/genuine-repeat estimates, allocate and pre-fault the
 * bloom filter, then stream the digits chunk by chunk.  For every window
 * admitted by the filter the key is rolled in O(1), hashed, prefetched and
 * tested PIPE_DEPTH windows later.  Ends with scan and bloom statistics.
 */
int main(int argc, char **argv)
{
    parse_args(argc, argv);
    setvbuf(stdout, NULL, _IOLBF, 0);
    install_signals();
    winspec_init(&W, opt.L);

    char outname[512];
    if (!opt.outfile) {
        char cnt[32];
        if (opt.ndigits == UINT64_MAX) snprintf(cnt, sizeof cnt, "all");
        else snprintf(cnt, sizeof cnt, "%" PRIu64, opt.ndigits);
        if (opt.F.len) snprintf(outname, sizeof outname, "candidates_L%d_N%s_F%s.txt", opt.L, cnt, opt.F.str);
        else snprintf(outname, sizeof outname, "candidates_L%d_N%s.txt", opt.L, cnt);
        opt.outfile = outname;
    }

    char b1[64], b2[64];
    printf("pibloom v%s - A197123 pass 1 (bloom candidate finder)\n", PI_TOOLS_VERSION);
    printf("  pi source      : %s\n", opt.pifile);
    printf("  sequence length: %d digits%s\n", opt.L, opt.L > LO_DIGITS ? "  (two-word key)" : "");
    printf("  digits to scan : %s\n", opt.ndigits == UINT64_MAX ? "all" : fmt_u64(opt.ndigits, b1, sizeof b1));
    if (opt.F.len) {
        char fd[32];
        printf("  filter         : %d-digit prefix in %s (%.3g%% of windows)\n", opt.F.len,
               pfilter_describe(&opt.F, fd, sizeof fd), 100.0 * pfilter_fraction(&opt.F));
    }
    printf("  bloom size     : %s (%s bits), k=%d, %s\n", fmt_bytes((double)opt.bloom_bytes, b1, sizeof b1),
           fmt_u64(opt.bloom_bytes * 8, b2, sizeof b2), K_HASHES, opt.blocked ? "blocked (1 cache line)" : "classic (4 independent bits)");
    printf("  candidate file : %s\n", opt.outfile);

    /* Expected behaviour, so the user can sanity-check the size before a long run */
    {
        pireader probe;
        pireader_open(&probe, opt.pifile, opt.L, 1 << 16, 1);
        uint64_t file_digits = probe.digits_est;
        char desc[256];
        pistream_describe(&probe.ps, desc, sizeof desc);
        printf("  source format  : %s\n", desc);
        pireader_close(&probe);
        uint64_t n = opt.ndigits == UINT64_MAX ? file_digits : opt.ndigits;
        if (file_digits && n > file_digits) {
            warn("file holds ~%s digits, fewer than requested; scanning what is there", fmt_u64(file_digits, b1, sizeof b1));
            n = file_digits;
        }
        g_total_digits = n;
        double ins = (double)n * pfilter_fraction(&opt.F);
        double m = (double)opt.bloom_bytes * 8.0;
        double fill = 1.0 - exp(-K_HASHES * ins / m);
        double fp = pow(fill, K_HASHES);
        printf("  expected       : ~%.3g insertions, final fill ~%.2f%%, false-positive rate ~%.3g%% (~%.3g spurious candidates)\n",
               ins, fill * 100, fp * 100, fp * ins);
        {
            /* birthday estimate of GENUINE repeats among the tested windows:
             * pairs ~ ins^2 / (2 * space), space = 10^L (10^(L-1) with -f,
             * since both occurrences share the fixed first digit) */
            double space = pow(10.0, opt.L) * pfilter_fraction(&opt.F);
            double real = ins * ins / (2.0 * space);
            if (real > ins) real = ins;   /* saturated: repeats everywhere */
            printf("                   plus ~%.3g genuine repeats expected at this length (birthday estimate)\n", real);
        }
        if (fill > 0.5) warn("bloom filter will be > 50%% full; expect a large candidate file. Increase -b or use -f.");
        if (fill > 0.9) warn("bloom filter far too small for this many digits - results will be nearly useless.");
    }

    bloom_bits = opt.bloom_bytes * 8;
    bloom_blocks = bloom_bits / BLOCK_BITS;
    check_fits_in_ram(opt.bloom_bytes, "bloom filter");
    bloom = alloc_huge(opt.bloom_bytes, opt.prefault, "bloom filter");

    outf = fopen(opt.outfile, "w");
    if (!outf) die("cannot create candidate file '%s': %s", opt.outfile, strerror(errno));
    setvbuf(outf, NULL, _IOFBF, 1 << 20);
    fprintf(outf, "# pibloom v%s L=%d filter=%s bloom_bytes=%" PRIu64 " blocked=%d pi=%s\n",
            PI_TOOLS_VERSION, opt.L, opt.F.len ? opt.F.str : "-", opt.bloom_bytes, opt.blocked, opt.pifile);

    pireader rd;
    pireader_open(&rd, opt.pifile, opt.L, opt.chunk, opt.ndigits);
    g_rd = &rd;

    printf("Scanning ...\n");
    t_start = t_last = now_sec();

    pipe_ent pipe[PIPE_DEPTH];
    const pfilter F = opt.F;
    const int flen = F.len;
    int interrupted = 0;
    uint64_t last_pos = 0;

    size_t nwin;
    while ((nwin = pireader_fill(&rd)) > 0) {
        const uint8_t *d = rd.buf;
        uint64_t hi, lo;
        win_compute(&W, d, &hi, &lo);
        uint64_t pos = rd.base;
        unsigned head = 0, inflight = 0;

        for (size_t i = 0; i < nwin; i++, pos++) {
            if (flen == 0 || pfilter_match(&F, d + i)) {
                pipe_ent *e = &pipe[head];
                if (inflight == PIPE_DEPTH) { process_ent(e); } else inflight++;
                e->hi = hi; e->lo = lo; e->pos = pos;
                compute_idx(hi, lo, e->idx);
                prefetch_idx(e->idx);
                head = (head + 1) & (PIPE_DEPTH - 1);
                st_tested++;
            }
            if (i + 1 < nwin)                 /* no roll past the last window */
                win_roll(&W, d + i, &hi, &lo);
        }
        /* drain pipeline at chunk end */
        unsigned tail = (head + PIPE_DEPTH - inflight) & (PIPE_DEPTH - 1);
        while (inflight--) { process_ent(&pipe[tail]); tail = (tail + 1) & (PIPE_DEPTH - 1); }

        st_windows += nwin;
        last_pos = pos - 1;
        progress(0, last_pos);
        if (g_stop) { interrupted = 1; break; }
    }
    progress(1, last_pos);

    double elapsed = now_sec() - t_start;
    if (fflush(outf) || ferror(outf)) die("write error on candidate file");
    fprintf(outf, "# end windows=%" PRIu64 " tested=%" PRIu64 " candidates=%" PRIu64 " last_pos=%" PRIu64 "%s\n",
            st_windows, st_tested, st_candidates, last_pos, interrupted ? " INTERRUPTED" : "");
    if (fclose(outf)) die("close error on candidate file");

    char d1[32];
    printf("\n%s\n", interrupted ? "*** INTERRUPTED - partial results ***" : "Scan complete");
    printf("  windows scanned   : %s (positions %s..%s)\n", fmt_u64(st_windows, b1, sizeof b1),
           fmt_u64(rd.first_pos, d1, sizeof d1), fmt_u64(last_pos, b2, sizeof b2));
    printf("  windows tested    : %s\n", fmt_u64(st_tested, b1, sizeof b1));
    printf("  candidates        : %s (%.6f%% of tested)\n", fmt_u64(st_candidates, b1, sizeof b1),
           st_tested ? 100.0 * (double)st_candidates / (double)st_tested : 0.0);
    printf("  non-digit bytes   : %s skipped\n", fmt_u64(rd.total_skipped, b1, sizeof b1));
    printf("  elapsed           : %s  (%.2f M windows/s, %.2f M tested/s)\n", fmt_duration(elapsed, d1, sizeof d1),
           elapsed > 0 ? st_windows / elapsed / 1e6 : 0, elapsed > 0 ? st_tested / elapsed / 1e6 : 0);
    printf("  peak RSS          : %.0f MiB\n", peak_rss_mb());

    printf("\nBloom filter statistics\n");
    printf("  bits set (tracked): %s of %s (%.4f%%)\n", fmt_u64(st_bits_set, b1, sizeof b1),
           fmt_u64(bloom_bits, b2, sizeof b2), 100.0 * (double)st_bits_set / (double)bloom_bits);
    {
        double fill_th = 1.0 - exp(-(double)K_HASHES * (double)st_tested / (double)bloom_bits);
        double fill = (double)st_bits_set / (double)bloom_bits;
        printf("  theoretical fill  : %.4f%%   (ratio measured/theory %.4f - ~1.0 means hashes are behaving)\n",
               fill_th * 100, fill_th > 0 ? fill / fill_th : 0);
        printf("  false-pos rate    : %.4g%% per test at final fill  (~%.3g spurious candidates expected over the run)\n",
               pow(fill, K_HASHES) * 100, pow(fill_th, K_HASHES) * (double)st_tested / 2.0);
        printf("  bits per insert   : %.2f\n", st_tested ? (double)bloom_bits / (double)st_tested : 0);
    }
    if (!g_stop || !interrupted) {
        printf("  verifying with full popcount ...");
        fflush(stdout);
        double t0 = now_sec();
        uint64_t pc = popcount_bloom();
        printf(" %s bits (%.1fs)%s\n", fmt_u64(pc, b1, sizeof b1), now_sec() - t0,
               pc == st_bits_set ? " - matches tracked count" : " - MISMATCH with tracked count!");
    }
    printf("\nCandidate file: %s\nNext: pisearch -p %s -l %d -c %s%s\n", opt.outfile, opt.pifile, opt.L, opt.outfile,
           opt.ndigits == UINT64_MAX ? "" : " -n <same count>");
    pireader_close(&rd);
    munmap(bloom, opt.bloom_bytes);
    return interrupted ? 130 : 0;
}
