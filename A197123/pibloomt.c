/*
 * pibloomt - multithreaded pass 1 for OEIS A197123.
 *
 * N worker threads (-t, default 10) share ONE bloom filter (atomic bit
 * test-and-set).  Work is partitioned by digit VALUE, not by file region:
 * the three digits right after the -f prefix form a number v in 0..999 and
 * worker w owns the windows with v in [w*1000/N, (w+1)*1000/N).  With N=10
 * that is exactly "split by the first digit after the prefix", with N=100
 * the first two; any N up to 256 balances to within 0.1%.
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
#include "pi_common.h"          /* digit sources, hashing, window keys, filters, timing, alloc helpers */
#include <pthread.h>            /* reader/worker threads, mutex + condition variables for the chunk ring */
#include <stdatomic.h>          /* atomic bit test-and-set on the shared bloom, lock-free counters */

#define MAX_WORKERS 256               /* -t upper bound */
#define DEFAULT_WORKERS 10
static int nworkers = DEFAULT_WORKERS;
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
    pfilter F;                   /* prefix-range filter (len 0 = none) */
    uint64_t bloom_bytes;
    int blocked;
    int prefault;
    int dryrun;                  /* -T: read, decode, dispatch and hash, but no bloom */
    double progress_secs;
    size_t chunk;
} opt = { NULL, NULL, UINT64_MAX, 0, {0,0,0,0,0,""}, 0, 0, 1, 0, 30.0, 64u << 20 };

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
static wstat_t wstats[MAX_WORKERS];

/*
 * usage - print the command-line help to stderr and exit(1).
 */
static void usage(const char *argv0)
{
    fprintf(stderr,
"pibloomt v" PI_TOOLS_VERSION " - A197123 pass 1, multithreaded\n"
"Usage: %s -p <pifile> -l <seqlen> -b <bloomsize> [options]\n"
"\n"
"Multithreaded pibloom: 1 reader + N workers (-t) over one shared bloom filter.\n"
"Workers partition the windows by the digits that follow the -f prefix (with\n"
"-t 10 that is the next digit; -t 100 the next two), so both occurrences of\n"
"any repeat always land on the same worker in file order.\n"
"\n"
"  -p SRC      pi digits: text or .ycd file(s)/directories, comma-separated\n"
"  -l N        sequence length (2..%d)\n"
"  -b SIZE     bloom filter size in bytes: 1M .. 16T (e.g. 512M, 64G, 3.5T)\n"
"  -n COUNT    digits of pi to scan (500M, 10B, 1T, all)   [all]\n"
"  -f SPEC     restrict to windows whose 1-3 digit prefix is SPEC: a value (3, 24),\n"
"              an inclusive range (0-4, 10-15, 05-09), or K/N = the K-th of N\n"
"              equal partitions (2/20 = [05-09])\n"
"  -o FILE     candidate output file [candidates_L<len>_N<count>[_F<d>].txt]\n"
"  -t N        worker threads, 1..%d [%d]; more threads hide more memory\n"
"              latency on huge tables (needs -l >= prefix length + 3)\n"
"  -B          blocked bloom (1 cache line per window): faster, slightly more\n"
"              spurious candidates\n"
"  -P SECS     progress interval [30]\n"
"  -c SIZE     read chunk size [64M]\n"
"  -N          do not pre-fault the bloom memory\n"
"  -T          dry run: read, decode, dispatch and hash every window but allocate\n"
"              no bloom and insert nothing - measures the source + pipeline rate\n"
"  -V          print version and exit\n"
"  -h          this help\n", argv0, MAX_WORKERS, DEFAULT_WORKERS, MAX_SEQ_LEN);
    exit(1);
}

/*
 * parse_args - getopt into the global opt struct with validation; the
 * sequence must be at least one digit longer than the filter prefix so
 * the worker split digit lies inside the window.
 */
static void parse_args(int argc, char **argv)
{
    int c;
    while ((c = getopt(argc, argv, "p:l:b:n:f:o:t:BP:c:NTVh")) != -1) {
        switch (c) {
        case 'p': opt.pifile = optarg; break;
        case 'l': opt.L = atoi(optarg); break;
        case 'b': if (parse_mem_size(optarg, &opt.bloom_bytes)) die("bad bloom size '%s'", optarg); break;
        case 'n': if (parse_count(optarg, &opt.ndigits)) die("bad digit count '%s'", optarg); break;
        case 'f':
            if (pfilter_parse(optarg, &opt.F)) die("bad filter '%s' (use 3, 24, 0-4, 10-15, 2/20 ...)", optarg);
            break;
        case 'o': opt.outfile = optarg; break;
        case 't': nworkers = atoi(optarg);
                  if (nworkers < 1 || nworkers > MAX_WORKERS) die("-t must be 1..%d", MAX_WORKERS);
                  break;
        case 'B': opt.blocked = 1; break;
        case 'P': opt.progress_secs = atof(optarg); if (opt.progress_secs <= 0) die("bad -P"); break;
        case 'c': { uint64_t v; if (parse_mem_size(optarg, &v) || v < (1u << 20)) die("bad chunk size"); opt.chunk = (size_t)v; } break;
        case 'N': opt.prefault = 0; break;
        case 'T': opt.dryrun = 1; break;
        case 'V': print_version("pibloomt"); exit(0);
        default: usage(argv[0]);
        }
    }
    if (!opt.pifile) { fprintf(stderr, "missing -p <pifile>\n"); usage(argv[0]); }
    if (opt.L < 2 || opt.L > MAX_SEQ_LEN) { fprintf(stderr, "missing/invalid -l (2..%d)\n", MAX_SEQ_LEN); usage(argv[0]); }
    if (opt.dryrun && opt.bloom_bytes == 0) opt.bloom_bytes = 1ULL << 20;   /* size is irrelevant in a dry run */
    if (opt.bloom_bytes < (1ULL << 20) || opt.bloom_bytes > (16ULL << 40)) {
        fprintf(stderr, "bloom size must be between 1M and 16T bytes\n"); usage(argv[0]);
    }
    if (opt.ndigits == 0) die("-n must be > 0");
    /* the three partition digits sit right after the prefix and must be inside the window */
    if (opt.F.len + 3 > opt.L) die("sequence length %d too short: pibloomt needs -l >= prefix length + 3 (use pibloom for shorter sequences)", opt.L);
}

/* ---------------- hashing / bloom (atomic) ---------------- */
/*
 * compute_idx - window key -> K_HASHES bit indices (classic or blocked);
 * identical to pibloom so both tools produce comparable candidates.
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
        idx[0] = fastrange64(h1,          bloom_bits);
        idx[1] = fastrange64(h1 + h2,     bloom_bits);
        idx[2] = fastrange64(h1 + 2 * h2, bloom_bits);
        idx[3] = fastrange64(h1 + 3 * h2, bloom_bits);
    }
}

/*
 * prefetch_idx - prefetch the cache line(s) holding the window's bits.
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

/*
 * process_ent - atomic test-and-set of the window's bits on the SHARED
 * bloom.  A plain load first avoids dirtying a cache line whose bit is
 * already set; otherwise fetch_or reports whether this thread set it.
 * All bits already set => candidate, written under the output mutex.
 * Per-worker counters avoid false sharing.
 */
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
/*
 * reader_main - the single reader thread: pulls chunks from the digit
 * stream (with the L-1 digit overlap already applied), copies each into a
 * free ring slot, and publishes it to all workers.  Blocks only when all
 * NSLOTS chunks are still being processed.  Signals producing_done at EOF.
 */
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
        s->remaining = nworkers;
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

/* ---- SIMD window scanner ----------------------------------------------
 * A worker owns windows whose prefix passes the -f filter AND whose
 * partition value v (3 digits after the prefix) lies in [vlo, vhi].  Rather
 * than memchr-ing for one digit (10% hit density => ~15M call restarts per
 * chunk per worker, which was the pipeline's bottleneck), we test 16
 * positions at a time against per-byte range constraints - a cheap
 * NECESSARY condition covering the prefix digits and the leading digits of
 * v - and run the exact test only on lanes that pass. */
typedef uint8_t v16u __attribute__((vector_size(16)));
typedef int8_t  v16s __attribute__((vector_size(16)));

typedef struct { int off; uint8_t lo, hi; } bytecon;

/* add per-digit range constraints for a number range [lo,hi] of `len`
 * digits placed at byte offset `base`: the leading digit is bounded by the
 * leading digits of lo/hi, and further digits only while they agree */
static int add_range_cons(bytecon *c, int n, int base, int len, int lo, int hi)
{
    int div = 1;
    for (int i = 1; i < len; i++) div *= 10;
    for (int k = 0; k < len; k++) {
        int dlo = (lo / div) % 10, dhi = (hi / div) % 10;
        c[n].off = base + k; c[n].lo = (uint8_t)dlo; c[n].hi = (uint8_t)dhi; n++;
        if (dlo != dhi) break;              /* later digits are unconstrained */
        div /= 10;
    }
    return n;
}

static inline v16u load16(const uint8_t *p) { v16u v; memcpy(&v, p, 16); return v; }

/* bitmask of lanes whose byte is nonzero */
static inline unsigned lanes_set(v16s m)
{
    uint64_t w[2];
    memcpy(w, &m, 16);
    unsigned bits = 0;
    for (int h = 0; h < 2; h++)
        for (int b = 0; b < 8; b++)
            if ((w[h] >> (8 * b)) & 0xFF) bits |= 1u << (h * 8 + b);
    return bits;
}


/*
 * worker_main - one of nworkers threads.  Walks every published chunk in
 * order but owns only the windows whose partition value v (the three
 * digits after the prefix, 0..999) falls in its range [vlo, vhi].  For each
 * leading digit its range covers it memchrs for that digit (SIMD-fast),
 * checks the prefix filter and the v range, computes the key from scratch,
 * then hashes/prefetches/tests through the same PIPE_DEPTH pipeline as
 * pibloom.  Because both occurrences of any repeat share the prefix and v,
 * they always land on the same worker and are seen in file order - no true
 * repeat can be lost to thread interleaving.  Releases the slot to the reader when it is the last done.
 */
static void *worker_main(void *argp)
{
    warg_t *wa = argp;
    const int id = wa->id;
    wstat_t *st = &wstats[id];
    const pfilter F = opt.F;
    const int flen = F.len;
    /* partition value range owned by this worker */
    const int vlo = id * 1000 / nworkers;
    const int vhi = (id + 1) * 1000 / nworkers - 1;
    /* byte-range constraints (necessary conditions) for the SIMD prefilter */
    bytecon cons[8];
    int ncons = 0;
    if (flen) ncons = add_range_cons(cons, ncons, 0, flen, F.lo, F.hi);
    ncons = add_range_cons(cons, ncons, flen, 3, vlo, vhi);
    int maxoff = 0;
    for (int k = 0; k < ncons; k++) if (cons[k].off > maxoff) maxoff = cons[k].off;
    v16u clo[8], chi[8];
    for (int k = 0; k < ncons; k++) { clo[k] = (v16u){0} + cons[k].lo; chi[k] = (v16u){0} + cons[k].hi; }

    pipe_ent pipe[PIPE_DEPTH];
    uint64_t local_tested = 0;
    const int dry = opt.dryrun;
    volatile uint64_t dry_sink = 0;          /* keeps the hash work from being optimised away */

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
        /* positions i in [0, nwin): 16-lane blocks, then a scalar tail.  The
         * bytes read reach d[i + 15 + maxoff]; maxoff <= flen + 2 <= L - 1, so
         * for i + 15 < nwin they stay inside the buffer's overlap. */
        size_t i0 = 0;
        for (; i0 + 16 <= nwin; i0 += 16) {
            v16s m = (v16s){0} - 1;                        /* all lanes on */
            for (int k = 0; k < ncons; k++) {
                v16u x = load16(d + i0 + cons[k].off);
                m &= (v16s)(x >= clo[k]) & (v16s)(x <= chi[k]);
            }
            unsigned bits = lanes_set(m);
            while (bits) {
                int j = __builtin_ctz(bits);
                bits &= bits - 1;
                size_t i = i0 + (size_t)j;
                if (flen && !pfilter_match(&F, d + i)) continue;
                int v = d[i + flen] * 100 + d[i + flen + 1] * 10 + d[i + flen + 2];
                if (v < vlo || v > vhi) continue;
                uint64_t hi, lo;
                win_compute(&W, d + i, &hi, &lo);
                pipe_ent *e = &pipe[head];
                if (dry) {                                     /* -T: hash only */
                    compute_idx(hi, lo, e->idx);
                    dry_sink ^= e->idx[0];
                } else {
                    if (inflight == PIPE_DEPTH) process_ent(e, st); else inflight++;
                    e->hi = hi; e->lo = lo; e->pos = base + i;
                    compute_idx(hi, lo, e->idx);
                    prefetch_idx(e->idx);
                    head = (head + 1) & (PIPE_DEPTH - 1);
                }
                local_tested++;
                if ((local_tested & 0xFFFFF) == 0)
                    atomic_store_explicit(&st->tested, local_tested, memory_order_relaxed);
            }
        }
        for (size_t i = i0; i < nwin; i++) {               /* scalar tail */
            if (flen && !pfilter_match(&F, d + i)) continue;
            int v = d[i + flen] * 100 + d[i + flen + 1] * 10 + d[i + flen + 2];
            if (v < vlo || v > vhi) continue;
            uint64_t hi, lo;
            win_compute(&W, d + i, &hi, &lo);
            pipe_ent *e = &pipe[head];
            if (dry) {
                compute_idx(hi, lo, e->idx);
                dry_sink ^= e->idx[0];
            } else {
                if (inflight == PIPE_DEPTH) process_ent(e, st); else inflight++;
                e->hi = hi; e->lo = lo; e->pos = base + i;
                compute_idx(hi, lo, e->idx);
                prefetch_idx(e->idx);
                head = (head + 1) & (PIPE_DEPTH - 1);
            }
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
/*
 * agg - sum the per-worker tested/candidate/bits-set counters for the
 * progress line and the final report.
 */
static void agg(uint64_t *tested, uint64_t *cand, uint64_t *bits)
{
    uint64_t t = 0, c = 0, b = 0;
    for (int i = 0; i < nworkers; i++) {
        t += atomic_load_explicit(&wstats[i].tested, memory_order_relaxed);
        c += atomic_load_explicit(&wstats[i].candidates, memory_order_relaxed);
        b += atomic_load_explicit(&wstats[i].bits_set, memory_order_relaxed);
    }
    *tested = t; *cand = c; *bits = b;
}

/* ---------------- parallel popcount ---------------- */
typedef struct { uint64_t from, to, sum; } pcarg_t;
/*
 * pc_main - popcount worker: counts set bits over one slice of the table.
 */
static void *pc_main(void *a)
{
    pcarg_t *pa = a;
    const uint64_t *w = (const uint64_t *)bloom;
    uint64_t s = 0;
    for (uint64_t i = pa->from; i < pa->to; i++) s += (uint64_t)__builtin_popcountll(w[i]);
    pa->sum = s;
    return NULL;
}
/*
 * popcount_bloom_parallel - full popcount of the bloom split across
 * the worker threads (a multi-terabyte table takes a while single-threaded);
 * cross-checks the tracked bit count.
 */
static uint64_t popcount_bloom_parallel(void)
{
    uint64_t nw = opt.bloom_bytes / 8;
    pthread_t th[MAX_WORKERS];
    pcarg_t args[MAX_WORKERS];
    for (int i = 0; i < nworkers; i++) {
        args[i].from = nw * i / (uint64_t)nworkers;
        args[i].to   = nw * (i + 1) / (uint64_t)nworkers;
        args[i].sum  = 0;
        if (pthread_create(&th[i], NULL, pc_main, &args[i])) die("pthread_create failed");
    }
    uint64_t sum = 0;
    for (int i = 0; i < nworkers; i++) { pthread_join(th[i], NULL); sum += args[i].sum; }
    for (uint64_t i = nw * 8; i < opt.bloom_bytes; i++) sum += (uint64_t)__builtin_popcount(bloom[i]);
    return sum;
}

/*
 * main - threaded pass 1 driver: options, banner and estimates, bloom
 * allocation, then starts the reader and worker threads and acts as the
 * progress reporter until the ring drains.  Final statistics include the
 * per-worker balance (should be near-equal) and a parallel popcount.
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
    printf("pibloomt v%s - A197123 pass 1, %d worker threads + reader, shared bloom%s\n", PI_TOOLS_VERSION, nworkers,
           opt.dryrun ? "  *** DRY RUN: no bloom, nothing inserted ***" : "");
    printf("  pi source      : %s\n", opt.pifile);
    printf("  sequence length: %d digits%s\n", opt.L, opt.L > LO_DIGITS ? "  (two-word key)" : "");
    printf("  digits to scan : %s\n", opt.ndigits == UINT64_MAX ? "all" : fmt_u64(opt.ndigits, b1, sizeof b1));
    if (opt.F.len) {
        char fd[32];
        printf("  partition      : %d-digit prefix in %s (%.3g%% of windows); %d workers split by the digits after it\n",
               opt.F.len, pfilter_describe(&opt.F, fd, sizeof fd), 100.0 * pfilter_fraction(&opt.F), nworkers);
    }
    else
        printf("  partition      : full search; %d workers split by the leading digits\n", nworkers);
    if (!opt.dryrun) printf("  bloom size     : %s (%s bits), k=%d, %s, atomic bit ops\n",
           fmt_bytes((double)opt.bloom_bytes, b1, sizeof b1), fmt_u64(opt.bloom_bytes * 8, b2, sizeof b2),
           K_HASHES, opt.blocked ? "blocked (1 cache line)" : "classic (4 independent bits)");
    printf("  candidate file : %s\n", opt.outfile);

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
        if (opt.dryrun) {
            printf("  dry run        : ~%.3g windows would be inserted; no bloom is allocated\n", ins);
            goto skip_estimates;
        }
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
        if (opt.F.len == 0)
            printf("                   (no -f: this is a FULL search - 10x the insertions of a filtered run)\n");
        if (fill > 0.5) warn("bloom filter will be > 50%% full; expect a large candidate file. Increase -b or use -f.");
        if (fill > 0.9) warn("bloom filter far too small for this many digits - results will be nearly useless.");
    skip_estimates:;
    }

    bloom_bits = opt.bloom_bytes * 8;
    bloom_blocks = bloom_bits / BLOCK_BITS;
    if (!opt.dryrun) {
        check_fits_in_ram(opt.bloom_bytes, "bloom filter");
        bloom = alloc_huge(opt.bloom_bytes, opt.prefault, "bloom filter");
    }

    outf = fopen(opt.outfile, "w");
    if (!outf) die("cannot create candidate file '%s': %s", opt.outfile, strerror(errno));
    setvbuf(outf, NULL, _IOFBF, 1 << 20);
    fprintf(outf, "# pibloomt v%s L=%d filter=%s bloom_bytes=%" PRIu64 " blocked=%d threads=%d pi=%s\n",
            PI_TOOLS_VERSION, opt.L, opt.F.len ? opt.F.str : "-", opt.bloom_bytes, opt.blocked, nworkers, opt.pifile);

    pireader rd;
    pireader_open(&rd, opt.pifile, opt.L, opt.chunk, opt.ndigits);
    for (int i = 0; i < NSLOTS; i++) {
        slots[i].buf = malloc(opt.chunk + (size_t)opt.L + 64);
        if (!slots[i].buf) die("cannot allocate chunk ring");
        slots[i].remaining = 0;
    }

    printf("Scanning ...\n");
    double t_start = now_sec(), t_last = t_start;
    uint64_t last_tested = 0, last_wins = 0, last_raw = 0;
    const int g_tty = isatty(1);        /* status line overwrites itself on a terminal */

    pthread_t rth, wth[MAX_WORKERS];
    warg_t wargs[MAX_WORKERS];
    if (pthread_create(&rth, NULL, reader_main, &rd)) die("pthread_create(reader) failed");
    for (int i = 0; i < nworkers; i++) {
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
                pct = 100.0 * (double)wins / (double)g_total_digits;
                if (pct > 100.0) pct = 100.0;
                if (win_rate > 0 && g_total_digits > wins) eta = (double)(g_total_digits - wins) / win_rate;
            }
            char pb[24] = "";
            if (pct >= 0) snprintf(pb, sizeof pb, " (%.1f%% done)", pct);
            uint64_t raw_now = rd.ps.total_raw;      /* bytes read so far (racy read, display only) */
            double dig_rate = (t - t_last) > 0 ? (double)(wins - last_wins) / (t - t_last) : 0;
            double mb_rate  = (t - t_last) > 0 ? (double)(raw_now - last_raw) / (t - t_last) / 1e6 : 0;
            printf("%s[%s] pos %s%s  tested %s  cand %s  fill %.4f%%  %.1f M tested/s  %.0f M digits/s (%.0f MB/s)%s%s%s",
                   g_tty ? "\r\033[K" : "",
                   fmt_duration(t - t_start, d1, sizeof d1), fmt_u64(pos, c1, sizeof c1), pb,
                   fmt_u64(tt, c2, sizeof c2), fmt_u64(cc, c3, sizeof c3),
                   100.0 * (double)bb / (double)bloom_bits,
                   (t - t_last) > 0 ? (double)(tt - last_tested) / (t - t_last) / 1e6 : 0,
                   dig_rate / 1e6, mb_rate,
                   eta >= 0 ? "  ETA " : "", eta >= 0 ? fmt_duration(eta, d2, sizeof d2) : "",
                   g_tty ? "" : "\n");
            fflush(stdout);
            fflush(outf);
            t_last = t; last_tested = tt; last_wins = wins; last_raw = raw_now;
        }
    }
    pthread_join(rth, NULL);
    for (int i = 0; i < nworkers; i++) pthread_join(wth[i], NULL);
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
    printf("  windows scanned   : %s (positions %s..%s)\n", fmt_u64(g_windows, b1, sizeof b1),
           fmt_u64(rd.first_pos, d1, sizeof d1), fmt_u64(g_last_pos, b2, sizeof b2));
    printf("  windows tested    : %s\n", fmt_u64(st_tested, b1, sizeof b1));
    printf("  candidates        : %s (%.6f%% of tested)\n", fmt_u64(st_candidates, b1, sizeof b1),
           st_tested ? 100.0 * (double)st_candidates / (double)st_tested : 0.0);
    {
        uint64_t mn = UINT64_MAX, mx = 0, sum = 0;
        for (int i = 0; i < nworkers; i++) {
            uint64_t v = atomic_load(&wstats[i].tested);
            if (v < mn) mn = v;
            if (v > mx) mx = v;
            sum += v;
        }
        if (nworkers <= 20) {
            printf("  per worker tested : ");
            for (int i = 0; i < nworkers; i++)
                printf("%s%.1fM", i ? " " : "", (double)atomic_load(&wstats[i].tested) / 1e6);
            printf("\n");
        } else {
            printf("  per worker tested : %d workers, min %.1fM  avg %.1fM  max %.1fM\n", nworkers,
                   mn / 1e6, (double)sum / nworkers / 1e6, mx / 1e6);
        }
    }
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
    if (opt.dryrun) printf("  (dry run: no bloom was allocated; fill/popcount not applicable)\n");
    if (!interrupted && !opt.dryrun) {
        printf("  verifying with full popcount (%d threads) ...", nworkers);
        fflush(stdout);
        double t0 = now_sec();
        uint64_t pc = popcount_bloom_parallel();
        printf(" %s bits (%.1fs)%s\n", fmt_u64(pc, b1, sizeof b1), now_sec() - t0,
               pc == st_bits_set ? " - matches tracked count" : " - MISMATCH with tracked count!");
    }
    printf("\nCandidate file: %s\nNext: pisearch -p %s -l %d -c %s\n", opt.outfile, opt.pifile, opt.L, opt.outfile);
    pireader_close(&rd);
    if (bloom) munmap(bloom, opt.bloom_bytes);
    return interrupted ? 130 : 0;
}
