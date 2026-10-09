/*
 * pibloomt - multithreaded pass 1 for OEIS A197123.
 *
 * 10 worker threads share ONE bloom filter (atomic bit test-and-set):
 *   - without -f: worker t handles windows whose FIRST digit is t
 *   - with -f D : worker t handles windows whose first digit is D and
 *                 whose SECOND digit is t (same total work as pibloom -f D,
 *                 split 10 ways)
 * A separate reader thread streams the pi file into a ring of shared chunks;
 * every worker scans every chunk (memchr for its digit) so the file is read
 * exactly once.
 *
 * Candidate detection is identical to pibloom: a window whose 4 bloom bits
 * are all already set is written to the candidate file.  Bit races between
 * threads are benign: different workers never insert the same window (their
 * digit classes are disjoint), so a race only reorders inserts of different
 * windows - the same ambiguity the sequential scan order already has, and
 * pisearch verifies exactly either way.
 *
 * Build:  cc -O3 -march=native -pthread -o pibloomt pibloomt.c -lm
 */
#include "pi_common.h"
#include <pthread.h>
#include <stdatomic.h>

#define NWORKERS  10
#define K_HASHES  4
#define PIPE_DEPTH 32
#define BLOCK_BITS 512
#define NSLOTS    4

typedef struct {
    uint64_t idx[K_HASHES];
    uint64_t pos;
    uint64_t hi, lo;
} pipe_ent;

static struct {
    const char *pifile;
    const char *outfile;
    uint64_t ndigits;
    int L;
    int flen;                    /* filter prefix length: 0 = none, 1..2 */
    uint8_t fdig[2];
    char fstr[3];
    uint64_t bloom_bytes;
    int blocked;
    int prefault;
    double progress_secs;
    size_t chunk;
} opt = { NULL, NULL, UINT64_MAX, 0, 0, {0,0}, "", 0, 0, 1, 30.0, 64u << 20 };

static uint8_t *bloom;
static uint64_t bloom_bits, bloom_blocks;
static winspec  W;
static FILE    *outf;
static pthread_mutex_t out_mtx = PTHREAD_MUTEX_INITIALIZER;

/* ---------------- chunk ring ---------------- */
typedef struct {
    uint8_t *buf;
    size_t   nwin;               /* windows in this chunk */
    uint64_t base;               /* absolute pos of buf[0] */
    uint64_t seq;
    int      remaining;          /* workers still to process */
} slot_t;

static slot_t slots[NSLOTS];
static pthread_mutex_t ring_mtx = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t  ring_can_fill = PTHREAD_COND_INITIALIZER;
static pthread_cond_t  ring_can_take = PTHREAD_COND_INITIALIZER;
static uint64_t g_total_digits;
static uint64_t chunks_published = 0;
static int producing_done = 0;
static uint64_t g_windows = 0;      /* total windows published */
static uint64_t g_last_pos = 0;

/* ---------------- per-worker stats (cache-line padded) ---------------- */
typedef struct {
    _Atomic uint64_t tested, candidates, bits_set;
    char pad[64 - 3 * sizeof(_Atomic uint64_t) > 0 ? 64 - 3 * sizeof(_Atomic uint64_t) : 8];
} wstat_t;
static wstat_t wstats[NWORKERS];

static void usage(const char *argv0)
{
    fprintf(stderr,
"Usage: %s -p <pifile> -l <seqlen> -b <bloomsize> [options]\n"
"\n"
"Multithreaded pibloom: 1 reader + %d workers over one shared bloom filter.\n"
"Without -f the workers split the search by FIRST digit (full search);\n"
"with -f D (1 digit) they all search first-digit D, split by SECOND digit;\n"
"with -f DD (2 digits) they search that 2-digit prefix, split by THIRD digit.\n"
"\n"
"  -p FILES    pi source file, or comma list of up to 8 continuation files\n"
"  -l N        sequence length (2..%d)\n"
"  -b SIZE     bloom filter size in bytes: 1M .. 1T (e.g. 512M, 64G, 1.5T)\n"
"  -n COUNT    digits of pi to scan (500M, 10B, 1T, all)   [all]\n"
"  -f DIGITS   restrict to windows starting with this 1 or 2 digit prefix\n"
"  -o FILE     candidate output file [candidates_L<len>_N<count>[_F<d>].txt]\n"
"  -B          blocked bloom (1 cache line per window): faster, slightly more\n"
"              spurious candidates\n"
"  -P SECS     progress interval [30]\n"
"  -c SIZE     read chunk size [64M]\n"
"  -N          do not pre-fault the bloom memory\n"
"  -h          this help\n", argv0, NWORKERS, MAX_SEQ_LEN);
    exit(1);
}

static void parse_args(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:l:b:n:f:o:BP:c:Nh")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'l': opt.L = atoi(optarg); break;
        case 'b': if (parse_mem_size(optarg, &opt.bloom_bytes)) die("bad bloom size '%s'", optarg); break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 'f': {
            size_t fl = strlen(optarg);
            if (fl < 1 || fl > 2) die("filter must be 1 or 2 digits");
            for (size_t k = 0; k < fl; k++)
                if (optarg[k] < '0' || optarg[k] > '9') die("filter must be 1 or 2 decimal digits");
            opt.flen = (int)fl;
            opt.fdig[0] = (uint8_t)(optarg[0] - '0');
            opt.fdig[1] = fl == 2 ? (uint8_t)(optarg[1] - '0') : 0;
            snprintf(opt.fstr, sizeof opt.fstr, "%s", optarg);
        } break;
        case 'o': opt.outfile = optarg; break;
        case 'B': opt.blocked = 1; break;
        case 'P': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -P"); break;
        case 'c': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        case 'N': opt.prefault = 0; break;
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile) { fprintf(stderr, "missing -p <pifile>\n"); usage(argv[0]); }
    if (opt.L < 2 || opt.L > MAX_SEQ_LEN) { fprintf(stderr, "missing/invalid -l (2..%d)\n", MAX_SEQ_LEN); usage(argv[0]); }
    if (opt.bloom_bytes < (1ULL << 20) || opt.bloom_bytes > (1ULL << 40)) {
        fprintf(stderr, "bloom size must be between 1M and 1T bytes\n"); usage(argv[0]);
    }
    if (opt.ndigits == 0) die("-n must be > 0");
    /* the split digit sits right after the prefix, so it must be inside the window */
    if (opt.flen + 1 > opt.L) die("sequence length %d too short for a %d-digit filter (need prefix + split digit)", opt.L, opt.flen);
}

/* ---------------- hashing / bloom (atomic) ---------------- */
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
        idx[0] = fastrange64(h1,          bloom_bits);
        idx[1] = fastrange64(h1 + h2,     bloom_bits);
        idx[2] = fastrange64(h1 + 2 * h2, bloom_bits);
        idx[3] = fastrange64(h1 + 3 * h2, bloom_bits);
    }
}

static inline void prefetch_idx(const uint64_t *idx)
{
    __builtin_prefetch(&bloom[idx[0] >> 3], 1, 0);
    if (!opt.blocked) {
        __builtin_prefetch(&bloom[idx[1] >> 3], 1, 0);
        __builtin_prefetch(&bloom[idx[2] >> 3], 1, 0);
        __builtin_prefetch(&bloom[idx[3] >> 3], 1, 0);
    }
}

static inline void process_ent(const pipe_ent *e, wstat_t *st)
{
    unsigned all = 1;
    for (int i = 0; i < K_HASHES; i++) {
        uint64_t b = e->idx[i];
        uint8_t m = (uint8_t)(1u << (b & 7));
        _Atomic uint8_t *p = (_Atomic uint8_t *)&bloom[b >> 3];
        /* cheap read first: avoids dirtying shared cache lines once bits fill in */
        if (atomic_load_explicit(p, memory_order_relaxed) & m) continue;
        uint8_t prev = atomic_fetch_or_explicit(p, m, memory_order_relaxed);
        if (!(prev & m)) {
            all = 0;
            atomic_fetch_add_explicit(&st->bits_set, 1, memory_order_relaxed);
        }
    }
    if (all) {
        char s[MAX_SEQ_LEN + 1];
        win_to_ascii(&W, e->hi, e->lo, s);
        pthread_mutex_lock(&out_mtx);
        if (fprintf(outf, "%s:%" PRIu64 "\n", s, e->pos) < 0)
            die("write error on candidate file: %s", strerror(errno));
        pthread_mutex_unlock(&out_mtx);
        atomic_fetch_add_explicit(&st->candidates, 1, memory_order_relaxed);
    }
}

/* ---------------- reader thread ---------------- */
static void *reader_main(void *arg)
{
    pireader *rd = arg;
    uint64_t seq = 0;
    size_t nwin;
    while (!g_stop && (nwin = pireader_fill(rd)) > 0) {
        slot_t *s = &slots[seq % NSLOTS];
        pthread_mutex_lock(&ring_mtx);
        while (s->remaining != 0)
            pthread_cond_wait(&ring_can_fill, &ring_mtx);
        pthread_mutex_unlock(&ring_mtx);

        memcpy(s->buf, rd->buf, rd->ndigits);
        s->nwin = nwin;
        s->base = rd->base;
        s->seq  = seq;

        pthread_mutex_lock(&ring_mtx);
        s->remaining = NWORKERS;
        chunks_published = ++seq;
        g_windows += nwin;
        g_last_pos = rd->base + nwin - 1;
        pthread_cond_broadcast(&ring_can_take);
        pthread_mutex_unlock(&ring_mtx);
    }
    pthread_mutex_lock(&ring_mtx);
    producing_done = 1;
    pthread_cond_broadcast(&ring_can_take);
    pthread_mutex_unlock(&ring_mtx);
    return NULL;
}

/* ---------------- worker threads ---------------- */
typedef struct { int id; } warg_t;

static void *worker_main(void *argp)
{
    warg_t *wa = argp;
    const int id = wa->id;
    wstat_t *st = &wstats[id];
    const int flen = opt.flen;
    /* the byte we memchr for in the digit stream */
    const uint8_t scan_digit = (uint8_t)(flen ? opt.fdig[0] : id);
    const uint8_t split_digit = (uint8_t)id;    /* the digit this worker owns */
    const uint8_t f1 = opt.fdig[1];             /* only used when flen == 2 */

    pipe_ent pipe[PIPE_DEPTH];
    uint64_t local_tested = 0;

    for (uint64_t seq = 0;; seq++) {
        pthread_mutex_lock(&ring_mtx);
        while (chunks_published <= seq && !producing_done)
            pthread_cond_wait(&ring_can_take, &ring_mtx);
        if (chunks_published <= seq && producing_done) {
            pthread_mutex_unlock(&ring_mtx);
            break;
        }
        pthread_mutex_unlock(&ring_mtx);

        slot_t *s = &slots[seq % NSLOTS];
        const uint8_t *d = s->buf;
        const size_t nwin = s->nwin;
        const uint64_t base = s->base;

        unsigned head = 0, inflight = 0;
        const uint8_t *p = d, *end = d + nwin;
        while ((p = memchr(p, scan_digit, (size_t)(end - p))) != NULL) {
            size_t i = (size_t)(p - d);
            p++;
            if (flen == 1 && d[i + 1] != split_digit) continue;
            if (flen == 2 && (d[i + 1] != f1 || d[i + 2] != split_digit)) continue;
            uint64_t hi, lo;
            win_compute(&W, d + i, &hi, &lo);
            pipe_ent *e = &pipe[head];
            if (inflight == PIPE_DEPTH) process_ent(e, st); else inflight++;
            e->hi = hi; e->lo = lo; e->pos = base + i;
            compute_idx(hi, lo, e->idx);
            prefetch_idx(e->idx);
            head = (head + 1) & (PIPE_DEPTH - 1);
            local_tested++;
            if ((local_tested & 0xFFFFF) == 0)
                atomic_store_explicit(&st->tested, local_tested, memory_order_relaxed);
        }
        unsigned tail = (head + PIPE_DEPTH - inflight) & (PIPE_DEPTH - 1);
        while (inflight--) { process_ent(&pipe[tail], st); tail = (tail + 1) & (PIPE_DEPTH - 1); }

        pthread_mutex_lock(&ring_mtx);
        if (--s->remaining == 0)
            pthread_cond_signal(&ring_can_fill);
        pthread_mutex_unlock(&ring_mtx);
    }
    atomic_store_explicit(&st->tested, local_tested, memory_order_relaxed);
    return NULL;
}

/* ---------------- aggregate stats ---------------- */
static void agg(uint64_t *tested, uint64_t *cand, uint64_t *bits)
{
    uint64_t t = 0, c = 0, b = 0;
    for (int i = 0; i < NWORKERS; i++) {
        t += atomic_load_explicit(&wstats[i].tested, memory_order_relaxed);
        c += atomic_load_explicit(&wstats[i].candidates, memory_order_relaxed);
        b += atomic_load_explicit(&wstats[i].bits_set, memory_order_relaxed);
    }
    *tested = t; *cand = c; *bits = b;
}

/* ---------------- parallel popcount ---------------- */
typedef struct { uint64_t from, to, sum; } pcarg_t;
static void *pc_main(void *a)
{
    pcarg_t *pa = a;
    const uint64_t *w = (const uint64_t *)bloom;
    uint64_t s = 0;
    for (uint64_t i = pa->from; i < pa->to; i++) s += (uint64_t)__builtin_popcountll(w[i]);
    pa->sum = s;
    return NULL;
}
static uint64_t popcount_bloom_parallel(void)
{
    uint64_t nw = opt.bloom_bytes / 8;
    pthread_t th[NWORKERS];
    pcarg_t args[NWORKERS];
    for (int i = 0; i < NWORKERS; i++) {
        args[i].from = nw * i / NWORKERS;
        args[i].to   = nw * (i + 1) / NWORKERS;
        args[i].sum  = 0;
        if (pthread_create(&th[i], NULL, pc_main, &args[i])) die("pthread_create failed");
    }
    uint64_t sum = 0;
    for (int i = 0; i < NWORKERS; i++) { pthread_join(th[i], NULL); sum += args[i].sum; }
    for (uint64_t i = nw * 8; i < opt.bloom_bytes; i++) sum += (uint64_t)__builtin_popcount(bloom[i]);
    return sum;
}

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
        if (opt.flen) snprintf(outname, sizeof outname, "candidates_L%d_N%s_F%s.txt", opt.L, cnt, opt.fstr);
        else snprintf(outname, sizeof outname, "candidates_L%d_N%s.txt", opt.L, cnt);
        opt.outfile = outname;
    }

    char b1[64], b2[64];
    printf("pibloomt - A197123 pass 1, %d worker threads + reader, shared bloom\n", NWORKERS);
    printf("  pi source      : %s\n", opt.pifile);
    printf("  sequence length: %d digits%s\n", opt.L, opt.L > LO_DIGITS ? "  (two-word key)" : "");
    printf("  digits to scan : %s\n", opt.ndigits == UINT64_MAX ? "all" : fmt_u64(opt.ndigits, b1, sizeof b1));
    if (opt.flen == 2)
        printf("  partition      : first two digits == \"%s\" fixed, workers split by THIRD digit\n", opt.fstr);
    else if (opt.flen == 1)
        printf("  partition      : first digit == %s fixed, workers split by SECOND digit\n", opt.fstr);
    else
        printf("  partition      : workers split by FIRST digit (full search)\n");
    printf("  bloom size     : %s (%s bits), k=%d, %s, atomic bit ops\n",
           fmt_bytes((double)opt.bloom_bytes, b1, sizeof b1), fmt_u64(opt.bloom_bytes * 8, b2, sizeof b2),
           K_HASHES, opt.blocked ? "blocked (1 cache line)" : "classic (4 independent bits)");
    printf("  candidate file : %s\n", opt.outfile);

    {
        pireader probe;
        pireader_open(&probe, opt.pifile, opt.L, 1 << 16, 1);
        uint64_t file_digits = probe.file_size ? probe.file_size - probe.total_raw : 0;
        pireader_close(&probe);
        uint64_t n = opt.ndigits == UINT64_MAX ? file_digits : opt.ndigits;
        if (file_digits && n > file_digits) {
            warn("file holds ~%s digits, fewer than requested; scanning what is there", fmt_u64(file_digits, b1, sizeof b1));
            n = file_digits;
        }
        g_total_digits = n;
        double ins = (double)n / pow(10.0, opt.flen);
        double m = (double)opt.bloom_bytes * 8.0;
        double fill = 1.0 - exp(-K_HASHES * ins / m);
        double fp = pow(fill, K_HASHES);
        printf("  expected       : ~%.3g insertions, final fill ~%.2f%%, false-positive rate ~%.3g%% (~%.3g spurious candidates)\n",
               ins, fill * 100, fp * 100, fp * ins);
        {
            /* birthday estimate of GENUINE repeats among the tested windows:
             * pairs ~ ins^2 / (2 * space), space = 10^L (10^(L-1) with -f,
             * since both occurrences share the fixed first digit) */
            double space = pow(10.0, opt.L - opt.flen);
            double real = ins * ins / (2.0 * space);
            if (real > ins) real = ins;   /* saturated: repeats everywhere */
            printf("                   plus ~%.3g genuine repeats expected at this length (birthday estimate)\n", real);
        }
        if (opt.flen == 0)
            printf("                   (no -f: this is a FULL search - 10x the insertions of a filtered run)\n");
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
    fprintf(outf, "# pibloomt L=%d filter=%s bloom_bytes=%" PRIu64 " blocked=%d threads=%d pi=%s\n",
            opt.L, opt.flen ? opt.fstr : "-", opt.bloom_bytes, opt.blocked, NWORKERS, opt.pifile);

    pireader rd;
    pireader_open(&rd, opt.pifile, opt.L, opt.chunk, opt.ndigits);
    for (int i = 0; i < NSLOTS; i++) {
        slots[i].buf = malloc(opt.chunk + (size_t)opt.L + 64);
        if (!slots[i].buf) die("cannot allocate chunk ring");
        slots[i].remaining = 0;
    }

    printf("Scanning ...\n");
    double t_start = now_sec(), t_last = t_start;
    uint64_t last_tested = 0;

    pthread_t rth, wth[NWORKERS];
    warg_t wargs[NWORKERS];
    if (pthread_create(&rth, NULL, reader_main, &rd)) die("pthread_create(reader) failed");
    for (int i = 0; i < NWORKERS; i++) {
        wargs[i].id = i;
        if (pthread_create(&wth[i], NULL, worker_main, &wargs[i])) die("pthread_create(worker) failed");
    }

    /* main thread: progress reporting until reader+workers finish */
    int done = 0;
    while (!done) {
        struct timespec ts = { 0, 200 * 1000 * 1000 };
        nanosleep(&ts, NULL);
        pthread_mutex_lock(&ring_mtx);
        int busy = 0;
        for (int i = 0; i < NSLOTS; i++) busy |= slots[i].remaining != 0;
        done = producing_done && !busy;
        uint64_t pos = g_last_pos, wins = g_windows;
        pthread_mutex_unlock(&ring_mtx);
        double t = now_sec();
        if (!done && t - t_last >= opt.progress_secs) {
            uint64_t tt, cc, bb;
            agg(&tt, &cc, &bb);
            char c1[32], c2[32], c3[32], d1[32], d2[32];
            double win_rate = (t - t_start) > 0 ? (double)wins / (t - t_start) : 0;
            double eta = -1, pct = -1;
            if (g_total_digits) {
                pct = 100.0 * (double)pos / (double)g_total_digits;
                if (pct > 100.0) pct = 100.0;
                if (win_rate > 0 && g_total_digits > wins) eta = (double)(g_total_digits - wins) / win_rate;
            }
            char pb[24] = "";
            if (pct >= 0) snprintf(pb, sizeof pb, " (%.1f%% done)", pct);
            printf("[%s] pos %s%s  tested %s  cand %s  fill %.4f%%  %.1f M tested/s%s%s\n",
                   fmt_duration(t - t_start, d1, sizeof d1), fmt_u64(pos, c1, sizeof c1), pb,
                   fmt_u64(tt, c2, sizeof c2), fmt_u64(cc, c3, sizeof c3),
                   100.0 * (double)bb / (double)bloom_bits,
                   (t - t_last) > 0 ? (double)(tt - last_tested) / (t - t_last) / 1e6 : 0,
                   eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, d2, sizeof d2) : "");
            fflush(outf);
            t_last = t; last_tested = tt;
        }
    }
    pthread_join(rth, NULL);
    for (int i = 0; i < NWORKERS; i++) pthread_join(wth[i], NULL);
    int interrupted = g_stop != 0;

    double elapsed = now_sec() - t_start;
    uint64_t st_tested, st_candidates, st_bits_set;
    agg(&st_tested, &st_candidates, &st_bits_set);

    if (fflush(outf) || ferror(outf)) die("write error on candidate file");
    fprintf(outf, "# end windows=%" PRIu64 " tested=%" PRIu64 " candidates=%" PRIu64 " last_pos=%" PRIu64 "%s\n",
            g_windows, st_tested, st_candidates, g_last_pos, interrupted ? " INTERRUPTED" : "");
    if (fclose(outf)) die("close error on candidate file");

    char d1[32];
    printf("\n%s\n", interrupted ? "*** INTERRUPTED - partial results ***" : "Scan complete");
    printf("  windows scanned   : %s (positions 1..%s)\n", fmt_u64(g_windows, b1, sizeof b1), fmt_u64(g_last_pos, b2, sizeof b2));
    printf("  windows tested    : %s\n", fmt_u64(st_tested, b1, sizeof b1));
    printf("  candidates        : %s (%.6f%% of tested)\n", fmt_u64(st_candidates, b1, sizeof b1),
           st_tested ? 100.0 * (double)st_candidates / (double)st_tested : 0.0);
    printf("  per worker tested : ");
    for (int i = 0; i < NWORKERS; i++)
        printf("%s%.1fM", i ? " " : "", (double)atomic_load(&wstats[i].tested) / 1e6);
    printf("\n");
    printf("  non-digit bytes   : %s skipped\n", fmt_u64(rd.total_skipped, b1, sizeof b1));
    printf("  elapsed           : %s  (%.2f M windows/s through file, %.2f M tested/s)\n",
           fmt_duration(elapsed, d1, sizeof d1),
           elapsed > 0 ? g_windows / elapsed / 1e6 : 0, elapsed > 0 ? st_tested / elapsed / 1e6 : 0);
    printf("  peak RSS          : %.0f MiB\n", peak_rss_mb());

    printf("\nBloom filter statistics\n");
    printf("  bits set (tracked): %s of %s (%.4f%%)\n", fmt_u64(st_bits_set, b1, sizeof b1),
           fmt_u64(bloom_bits, b2, sizeof b2), 100.0 * (double)st_bits_set / (double)bloom_bits);
    {
        double fill_th = 1.0 - exp(-(double)K_HASHES * (double)st_tested / (double)bloom_bits);
        double fill = (double)st_bits_set / (double)bloom_bits;
        printf("  theoretical fill  : %.4f%%   (ratio measured/theory %.4f - ~1.0 means hashes are behaving)\n",
               fill_th * 100, fill_th > 0 ? fill / fill_th : 0);
        printf("  false-pos rate    : %.4g%% per test at final fill\n", pow(fill, K_HASHES) * 100);
        printf("  bits per insert   : %.2f\n", st_tested ? (double)bloom_bits / (double)st_tested : 0);
    }
    if (!interrupted) {
        printf("  verifying with full popcount (%d threads) ...", NWORKERS);
        fflush(stdout);
        double t0 = now_sec();
        uint64_t pc = popcount_bloom_parallel();
        printf(" %s bits (%.1fs)%s\n", fmt_u64(pc, b1, sizeof b1), now_sec() - t0,
               pc == st_bits_set ? " - matches tracked count" : " - MISMATCH with tracked count!");
    }
    printf("\nCandidate file: %s\nNext: pisearch -p %s -l %d -c %s\n", opt.outfile, opt.pifile, opt.L, opt.outfile);
    pireader_close(&rd);
    munmap(bloom, opt.bloom_bytes);
    return interrupted ? 130 : 0;
}
