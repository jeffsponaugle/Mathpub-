/* dschain.c -- search for chains of primes where p_{k+1} = p_k + digitsum(p_k)
 *
 * Build (with primesieve, recommended):
 *     cc -O3 -march=native -fopenmp -o dschain dschain.c -lprimesieve
 *
 * Build (no primesieve; uses the built-in segmented sieve, ~2-3x slower):
 *     cc -O3 -march=native -fopenmp -DNO_PRIMESIEVE -o dschain dschain.c -lm
 *
 * Usage:
 *     ./dschain --length 7 --start 1e8 --end 1e10 --threads 8
 *     ./dschain -l 8 -s 1e10 -e 3.2e14 -k run8.ckpt      (resumable)
 *     ./dschain -l 9 -s 2 -e 1e12 --first                (find earliest hit)
 *
 * Notes on correctness: candidate chain *starts* come from the sieve
 * over [start,end], but every chain *continuation* is tested with a deterministic
 * 64-bit Miller-Rabin. So chains are never truncated at the end of the range --
 * a chain starting at end-1 is followed as far as it actually goes.
 *
 * Checkpointing: work is divided into fixed chunks. The checkpoint records
 * which chunks are finished (as a watermark plus a short list of out-of-order
 * completions) and, for chunks that were interrupted mid-way, the last prime
 * fully processed in them. Re-running the same command with the same -k file
 * picks up where it stopped. Ctrl-C (or SIGTERM) writes a checkpoint before
 * exiting; a second Ctrl-C exits immediately without writing.
 */

#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>

/* progress line only rewrites itself when stdout is a terminal, so piping
 * or redirecting output gives clean plain text */
static const char *CLR = "";
static const char *NL  = "\n";

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef NO_PRIMESIEVE
#include <primesieve.h>
#endif

/* ------------------------------------------------------------------ */
/* digit sum: 4 digits at a time via lookup table                      */
/* ------------------------------------------------------------------ */

#define DS_TABLE_SIZE 10000
static uint8_t ds_table[DS_TABLE_SIZE];

static void ds_init(void) {
    for (int i = 0; i < DS_TABLE_SIZE; i++) {
        int s = 0, n = i;
        while (n) { s += n % 10; n /= 10; }
        ds_table[i] = (uint8_t)s;
    }
}

static inline unsigned digitsum(uint64_t n) {
    unsigned s = 0;
    while (n) { s += ds_table[n % DS_TABLE_SIZE]; n /= DS_TABLE_SIZE; }
    return s;
}

/* max digit sum of any uint64 is 9*20 = 180 */
#define MAX_DIGITSUM 180

/* ------------------------------------------------------------------ */
/* deterministic Miller-Rabin for all 64-bit n                         */
/* ------------------------------------------------------------------ */

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) {
    return (uint64_t)((__uint128_t)a * b % m);
}

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static const uint32_t SMALL_PRIMES[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97
};
#define N_SMALL (sizeof(SMALL_PRIMES)/sizeof(SMALL_PRIMES[0]))

static int is_prime(uint64_t n) {
    if (n < 2) return 0;
    for (size_t i = 0; i < N_SMALL; i++) {
        uint64_t p = SMALL_PRIMES[i];
        if (n == p) return 1;
        if (n % p == 0) return 0;
    }
    if (n < 97ULL * 97ULL) return 1;

    /* these 12 bases are a deterministic witness set for all n < 3.3e24 */
    static const uint64_t bases[] = {2,3,5,7,11,13,17,19,23,29,31,37};

    uint64_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }

    for (size_t i = 0; i < sizeof(bases)/sizeof(bases[0]); i++) {
        uint64_t x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int witness = 1;
        for (int k = 1; k < r; k++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { witness = 0; break; }
        }
        if (witness) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* segment primality bitmap                                            */
/*                                                                     */
/* A chain step is at most MAX_DIGITSUM, so walks starting in a chunk  */
/* stay close to it. Each chunk sieves a margin past its end and       */
/* answers is-prime for the whole walk with an O(1) bit test instead   */
/* of Miller-Rabin (~100x cheaper). Odds only; anything outside the    */
/* window falls back to Miller-Rabin.                                  */
/* ------------------------------------------------------------------ */

#define SEG_MARGIN 4096         /* covers > 20 links of maximal digit sum */

typedef struct {
    uint64_t base;              /* first odd number covered */
    uint64_t hi;                /* last number covered */
    const uint64_t *bm;         /* 1 bit per odd number, set = prime */
} seg_t;

static inline int seg_prime(const seg_t *s, uint64_t q) {
    if (!(q & 1)) return q == 2;
    if (s && q >= s->base && q <= s->hi) {
        uint64_t idx = (q - s->base) >> 1;
        return (s->bm[idx >> 6] >> (idx & 63)) & 1ULL;
    }
    return is_prime(q);
}

/* ------------------------------------------------------------------ */
/* chain walking                                                       */
/* ------------------------------------------------------------------ */

/* Length of the chain starting at prime p. Fills chain[] if non-NULL.
 * seg may be NULL to test purely with Miller-Rabin. */
static int chain_length(uint64_t p, uint64_t *chain, int cap, const seg_t *seg) {
    int n = 1;
    if (chain && cap > 0) chain[0] = p;
    for (;;) {
        unsigned d = digitsum(p);
        /* p is odd (p>2) so an odd digit sum makes p+d even -> chain ends */
        if (n > 1 || p != 2) {
            if (d & 1) break;
        }
        uint64_t q = p + d;
        if (q < p) break;                 /* overflow guard */
        if (!seg_prime(seg, q)) break;
        p = q;
        if (chain && n < cap) chain[n] = p;
        n++;
    }
    return n;
}

/* Is p the *start* of a maximal chain, i.e. is there no prime r with
 * r + digitsum(r) == p ?  Steps are at most MAX_DIGITSUM. */
static int is_chain_start(uint64_t p) {
    for (unsigned d = 1; d <= MAX_DIGITSUM; d++) {
        if (d >= p) break;
        uint64_t r = p - d;
        if (digitsum(r) == d && is_prime(r)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* candidate prime generation for a chunk                              */
/* ------------------------------------------------------------------ */

#ifdef NO_PRIMESIEVE
/* Built-in segmented sieve. base_primes holds every prime <= sqrt(end);
 * it is filled once in base_init() before any thread starts, then read-only. */
static uint64_t *base_primes = NULL;
static size_t n_base = 0;

static void base_init(uint64_t end) {
    uint64_t lim = (uint64_t)sqrt((double)end) + 2;
    if (lim < 3) lim = 3;
    if (lim > 4000000000ULL) {
        fprintf(stderr,
            "error: --end %" PRIu64 " needs base primes up to %" PRIu64 ",\n"
            "       which the built-in sieve cannot hold. Rebuild against\n"
            "       primesieve (drop -DNO_PRIMESIEVE) or lower --end.\n",
            end, lim);
        exit(1);
    }
    uint8_t *sv = (uint8_t *)calloc(lim + 1, 1);
    if (!sv) { fprintf(stderr, "out of memory building base primes\n"); exit(1); }
    size_t cnt = 0;
    for (uint64_t i = 2; i <= lim; i++) {
        if (sv[i]) continue;
        cnt++;
        for (uint64_t j = i * i; j <= lim; j += i) sv[j] = 1;
    }
    base_primes = (uint64_t *)malloc(cnt * sizeof(uint64_t));
    if (!base_primes) { fprintf(stderr, "out of memory building base primes\n"); exit(1); }
    for (uint64_t i = 2; i <= lim; i++) if (!sv[i]) base_primes[n_base++] = i;
    free(sv);
}
#endif

static uint64_t *gen_primes(uint64_t lo, uint64_t hi, size_t *count) {
#ifndef NO_PRIMESIEVE
    size_t n = 0;
    uint64_t *v = (uint64_t *)primesieve_generate_primes(lo, hi, &n,
                                                         UINT64_PRIMES);
    *count = n;
    return v;
#else
    if (lo < 2) lo = 2;
    if (hi < lo) { *count = 0; return NULL; }

    uint64_t span = hi - lo + 1;
    uint8_t *comp = (uint8_t *)calloc(span, 1);
    if (!comp) { fprintf(stderr, "out of memory sieving segment\n"); exit(1); }

    for (size_t i = 0; i < n_base; i++) {
        uint64_t p = base_primes[i];
        if (p * p > hi) break;
        uint64_t first = (lo + p - 1) / p * p;   /* first multiple of p >= lo */
        if (first < p * p) first = p * p;
        for (uint64_t j = first; j <= hi; j += p) comp[j - lo] = 1;
    }

    size_t n = 0;
    for (uint64_t k = 0; k < span; k++) if (!comp[k]) n++;
    uint64_t *v = (uint64_t *)malloc((n ? n : 1) * sizeof(uint64_t));
    if (!v) { fprintf(stderr, "out of memory collecting segment\n"); exit(1); }
    size_t j = 0;
    for (uint64_t k = 0; k < span; k++) if (!comp[k]) v[j++] = lo + k;
    free(comp);
    *count = n;
    return v;
#endif
}

static void free_primes(uint64_t *v) {
#ifndef NO_PRIMESIEVE
    primesieve_free(v);
#else
    free(v);
#endif
}

/* ------------------------------------------------------------------ */
/* progress reporting                                                  */
/* ------------------------------------------------------------------ */

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void fmt_duration(double s, char *buf, size_t n) {
    if (s < 0 || !isfinite(s)) { snprintf(buf, n, "--:--:--"); return; }
    if (s > 359999) { snprintf(buf, n, ">99h"); return; }
    int t = (int)(s + 0.5);
    snprintf(buf, n, "%02d:%02d:%02d", t / 3600, (t / 60) % 60, t % 60);
}

/* 1234567 -> "1,234,567" */
static void fmt_grouped(uint64_t v, char *buf, size_t n) {
    char tmp[24];
    int len = snprintf(tmp, sizeof tmp, "%" PRIu64, v);
    size_t o = 0;
    for (int i = 0; i < len && o + 2 < n; i++) {
        if (i && (len - i) % 3 == 0) buf[o++] = ',';
        buf[o++] = tmp[i];
    }
    buf[o] = 0;
}

static void fmt_count(double v, char *buf, size_t n) {
    static const char *u[] = {"", "K", "M", "G", "T", "P"};
    int i = 0;
    while (v >= 1000.0 && i < 5) { v /= 1000.0; i++; }
    snprintf(buf, n, "%.2f%s", v, u[i]);
}

/* ------------------------------------------------------------------ */
/* arg parsing                                                         */
/* ------------------------------------------------------------------ */

/* Magnitude suffixes. Case-insensitive; G is accepted as a synonym for B.
 * Note E means exa (1e18) only when it is a trailing suffix -- "3E" is
 * 3e18, while "3E9" is still read as scientific notation for 3000000000. */
static const struct { char c; long double mult; const char *name; } SUFFIXES[] = {
    { 'K', 1e3L,  "thousand"    },
    { 'M', 1e6L,  "million"     },
    { 'B', 1e9L,  "billion"     },
    { 'G', 1e9L,  "billion"     },
    { 'T', 1e12L, "trillion"    },
    { 'P', 1e15L, "quadrillion" },
    { 'E', 1e18L, "quintillion" },
};
#define N_SUFFIX (sizeof(SUFFIXES)/sizeof(SUFFIXES[0]))

#define UINT64_CEIL 18446744073709551615.0L

/* accepts 1234, 1e9, 2.5e8, 3.2e14, 1_000_000, 12M, 4B, 2.5T, 3P, 1E, 900k,
 * "1 000", and an optional trailing "b"/"bytes"-style letter is NOT allowed
 * (typos should be caught, not silently ignored). Exits on garbage. */
static uint64_t parse_num(const char *s, const char *what) {
    char clean[80];
    size_t j = 0;
    for (size_t i = 0; s[i]; i++) {
        if (s[i] == '_' || s[i] == ',' || s[i] == ' ' || s[i] == '\'') continue;
        if (j >= sizeof(clean) - 1) {
            fprintf(stderr, "error: value for %s is too long: %s\n", what, s);
            exit(2);
        }
        clean[j++] = s[i];
    }
    clean[j] = 0;

    if (j == 0) {
        fprintf(stderr, "error: %s needs a value\n", what);
        exit(2);
    }

    errno = 0;
    char *end = NULL;
    long double v = strtold(clean, &end);   /* handles 1e10, 3.2e14, .5e9 */

    if (end == clean) {
        fprintf(stderr, "error: %s is not a number: %s\n", what, s);
        exit(2);
    }
    if (v < 0) {
        fprintf(stderr, "error: %s must not be negative: %s\n", what, s);
        exit(2);
    }
    if (!isfinite((double)v)) {
        fprintf(stderr, "error: %s is not a finite number: %s\n", what, s);
        exit(2);
    }

    if (*end) {
        int matched = 0;
        char c = *end;
        if (c >= 'a' && c <= 'z') c = (char)(c - 'a' + 'A');
        for (size_t i = 0; i < N_SUFFIX; i++) {
            if (SUFFIXES[i].c == c) { v *= SUFFIXES[i].mult; matched = 1; break; }
        }
        if (!matched || end[1]) {
            fprintf(stderr,
                    "error: %s has an unrecognized suffix: %s\n"
                    "       valid suffixes: K M B(G) T P E "
                    "(thousand, million, billion, trillion, quadrillion, quintillion)\n",
                    what, s);
            exit(2);
        }
    }

    if (v > UINT64_CEIL) {
        fprintf(stderr, "error: %s exceeds the 64-bit limit (max ~18.4E): %s\n", what, s);
        exit(2);
    }
    return (uint64_t)v;
}

/* ------------------------------------------------------------------ */
/* shared search state                                                 */
/* ------------------------------------------------------------------ */

#define HIST_MAX 255            /* histogram bucket for lengths >= HIST_MAX */

static uint64_t g_hist[HIST_MAX + 1];
static uint64_t g_primes_tested = 0;
static uint64_t g_found = 0;
static uint64_t g_done_range = 0;
static int      g_longest = 0;
static uint64_t g_longest_at = 0;

/* stop flag, set by a signal */
static volatile sig_atomic_t g_stop = 0;
#define STOP_SIGNAL 1
static volatile sig_atomic_t g_stop_reason = 0;
static volatile sig_atomic_t g_signo = 0;

/* --first: smallest starting prime of any reported hit so far. Threads
 * abandon work above this bound but keep sweeping below it, so when the
 * parallel loop drains every prime <= the bound has been examined and the
 * bound is the confirmed earliest chain in the range. Written only inside
 * critical(output); read lock-free elsewhere with an atomic read. */
static uint64_t g_first_bound = UINT64_MAX;

static void handle_signal(int sig) {
    if (g_stop && g_stop_reason == STOP_SIGNAL) {
        static const char m[] =
            "\nsecond interrupt -- exiting immediately, no checkpoint written\n";
        ssize_t r = write(2, m, sizeof m - 1); (void)r;
        _exit(130);
    }
    g_signo = sig;
    g_stop_reason = STOP_SIGNAL;
    g_stop = 1;
    static const char m[] =
        "\ninterrupt -- stopping workers and writing checkpoint "
        "(Ctrl-C again to abort now)\n";
    ssize_t r = write(2, m, sizeof m - 1); (void)r;
}

/* ------------------------------------------------------------------ */
/* checkpoint bookkeeping                                              */
/*                                                                     */
/* Every read/write of the structures below happens inside the         */
/* critical(output) section, except the lock-free watermark fast path. */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t chunk, last_p; } partial_t;
#define MAX_PARTIALS 8192

static uint64_t *ck_bits    = NULL;     /* one bit per chunk, set when finished */
static uint64_t  ck_nchunks = 0;
static uint64_t  ck_watermark = 0;      /* every chunk < watermark is finished  */
static partial_t ck_part[MAX_PARTIALS]; /* chunks stopped part-way through      */
static int       ck_nparts = 0;
static int       ck_dropped_parts = 0;

static void ck_alloc(uint64_t nchunks) {
    ck_nchunks = nchunks;
    size_t words = (size_t)((nchunks + 63) / 64);
    ck_bits = (uint64_t *)calloc(words ? words : 1, sizeof(uint64_t));
    if (!ck_bits) { fprintf(stderr, "out of memory allocating chunk map\n"); exit(1); }
}

static inline int ck_is_done(uint64_t c) {
    if (c < ck_watermark) return 1;
    return (ck_bits[c >> 6] >> (c & 63)) & 1ULL;
}

static void ck_mark_done(uint64_t c) {
    ck_bits[c >> 6] |= 1ULL << (c & 63);
    while (ck_watermark < ck_nchunks &&
           ((ck_bits[ck_watermark >> 6] >> (ck_watermark & 63)) & 1ULL))
        ck_watermark++;
}

static uint64_t ck_partial_get(uint64_t c) {
    for (int i = 0; i < ck_nparts; i++)
        if (ck_part[i].chunk == c) return ck_part[i].last_p;
    return 0;
}

static void ck_partial_set(uint64_t c, uint64_t last_p) {
    for (int i = 0; i < ck_nparts; i++) {
        if (ck_part[i].chunk == c) {
            if (last_p > ck_part[i].last_p) ck_part[i].last_p = last_p;
            return;
        }
    }
    if (ck_nparts < MAX_PARTIALS) {
        ck_part[ck_nparts].chunk  = c;
        ck_part[ck_nparts].last_p = last_p;
        ck_nparts++;
    } else {
        ck_dropped_parts = 1;   /* only costs redundant work on resume */
    }
}

static void ck_partial_clear(uint64_t c) {
    for (int i = 0; i < ck_nparts; i++) {
        if (ck_part[i].chunk == c) {
            ck_part[i] = ck_part[--ck_nparts];
            return;
        }
    }
}

/* ---- checkpoint file ---- */

#define CK_VERSION 1

typedef struct {
    int      loaded;
    int      version;
    int      complete;
    uint64_t start, end, chunk;
    int      want_len, maximal;
    uint64_t nchunks, watermark;
    double   elapsed;
    uint64_t primes_tested, found;
    int      longest;
    uint64_t longest_at;
    uint64_t hist[HIST_MAX + 1];
    uint64_t *done_above;
    size_t    n_done_above;
    partial_t *partials;
    size_t     n_partials;
} ckpt_t;

typedef struct {
    uint64_t start, end, chunk;
    int      want_len, maximal;
    uint64_t nchunks;
    double   elapsed;
    int      complete;
} ck_meta_t;

/* Written from inside critical(output): the chunk map and partial list are
 * consistent, though chunks still in flight on other threads count as unfinished
 * and will be redone on resume. Temp file + rename keeps the old checkpoint
 * intact if we die mid-write. */
static int ck_write(const char *path, const ck_meta_t *m) {
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);

    FILE *f = fopen(tmp, "w");
    if (!f) {
        fprintf(stderr, "\nwarning: cannot write checkpoint %s: %s\n",
                tmp, strerror(errno));
        return -1;
    }

    fprintf(f, "# dschain checkpoint -- safe to delete, safe to inspect\n");
    fprintf(f, "version %d\n", CK_VERSION);
    fprintf(f, "complete %d\n", m->complete);
    fprintf(f, "start %" PRIu64 "\n", m->start);
    fprintf(f, "end %" PRIu64 "\n", m->end);
    fprintf(f, "chunk %" PRIu64 "\n", m->chunk);
    fprintf(f, "length %d\n", m->want_len);
    fprintf(f, "maximal %d\n", m->maximal);
    fprintf(f, "nchunks %" PRIu64 "\n", m->nchunks);
    fprintf(f, "watermark %" PRIu64 "\n", ck_watermark);
    fprintf(f, "elapsed %.3f\n", m->elapsed);
    fprintf(f, "primes %" PRIu64 "\n", g_primes_tested);
    fprintf(f, "found %" PRIu64 "\n", g_found);
    fprintf(f, "longest %d\n", g_longest);
    fprintf(f, "longest_at %" PRIu64 "\n", g_longest_at);
    for (int i = 1; i <= HIST_MAX; i++)
        if (g_hist[i]) fprintf(f, "hist %d %" PRIu64 "\n", i, g_hist[i]);

    /* finished chunks above the watermark (usually at most one per thread) */
    for (uint64_t c = ck_watermark; c < ck_nchunks; c++)
        if ((ck_bits[c >> 6] >> (c & 63)) & 1ULL)
            fprintf(f, "done %" PRIu64 "\n", c);

    for (int i = 0; i < ck_nparts; i++)
        fprintf(f, "partial %" PRIu64 " %" PRIu64 "\n",
                ck_part[i].chunk, ck_part[i].last_p);

    fflush(f);
    fsync(fileno(f));
    if (fclose(f) != 0) {
        fprintf(stderr, "\nwarning: error closing checkpoint %s: %s\n",
                tmp, strerror(errno));
        return -1;
    }
    if (rename(tmp, path) != 0) {
        fprintf(stderr, "\nwarning: cannot rename %s -> %s: %s\n",
                tmp, path, strerror(errno));
        return -1;
    }
    return 0;
}

static int ck_read(const char *path, ckpt_t *ck) {
    memset(ck, 0, sizeof *ck);
    FILE *f = fopen(path, "r");
    if (!f) return 0;                      /* no checkpoint yet */

    size_t cap_done = 64, cap_part = 64;
    ck->done_above = (uint64_t *)malloc(cap_done * sizeof(uint64_t));
    ck->partials   = (partial_t *)malloc(cap_part * sizeof(partial_t));
    if (!ck->done_above || !ck->partials) {
        fprintf(stderr, "out of memory reading checkpoint\n"); exit(1);
    }

    char line[256];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        uint64_t a, b; int i;
        if      (sscanf(line, "version %d", &ck->version) == 1) ;
        else if (sscanf(line, "complete %d", &ck->complete) == 1) ;
        else if (sscanf(line, "start %" SCNu64, &ck->start) == 1) ;
        else if (sscanf(line, "end %" SCNu64, &ck->end) == 1) ;
        else if (sscanf(line, "chunk %" SCNu64, &ck->chunk) == 1) ;
        else if (sscanf(line, "length %d", &ck->want_len) == 1) ;
        else if (sscanf(line, "maximal %d", &ck->maximal) == 1) ;
        else if (sscanf(line, "nchunks %" SCNu64, &ck->nchunks) == 1) ;
        else if (sscanf(line, "watermark %" SCNu64, &ck->watermark) == 1) ;
        else if (sscanf(line, "elapsed %lf", &ck->elapsed) == 1) ;
        else if (sscanf(line, "primes %" SCNu64, &ck->primes_tested) == 1) ;
        else if (sscanf(line, "found %" SCNu64, &ck->found) == 1) ;
        else if (sscanf(line, "longest %d", &ck->longest) == 1) ;
        else if (sscanf(line, "longest_at %" SCNu64, &ck->longest_at) == 1) ;
        else if (sscanf(line, "hist %d %" SCNu64, &i, &a) == 2) {
            if (i >= 1 && i <= HIST_MAX) ck->hist[i] = a;
        }
        else if (sscanf(line, "done %" SCNu64, &a) == 1) {
            if (ck->n_done_above == cap_done) {
                cap_done *= 2;
                ck->done_above = (uint64_t *)realloc(ck->done_above,
                                                     cap_done * sizeof(uint64_t));
                if (!ck->done_above) { fprintf(stderr, "out of memory\n"); exit(1); }
            }
            ck->done_above[ck->n_done_above++] = a;
        }
        else if (sscanf(line, "partial %" SCNu64 " %" SCNu64, &a, &b) == 2) {
            if (ck->n_partials == cap_part) {
                cap_part *= 2;
                ck->partials = (partial_t *)realloc(ck->partials,
                                                    cap_part * sizeof(partial_t));
                if (!ck->partials) { fprintf(stderr, "out of memory\n"); exit(1); }
            }
            ck->partials[ck->n_partials].chunk  = a;
            ck->partials[ck->n_partials].last_p = b;
            ck->n_partials++;
        }
    }
    fclose(f);
    ck->loaded = 1;
    return 1;
}

/* ------------------------------------------------------------------ */
/* histogram output                                                    */
/* ------------------------------------------------------------------ */

static void print_histogram(const uint64_t *hist, int want_len, int longest) {
    int top = want_len > longest ? want_len : longest;
    if (top > HIST_MAX) top = HIST_MAX;
    if (top < 1) top = 1;

    uint64_t cum[HIST_MAX + 2];
    cum[top + 1] = 0;
    for (int L = top; L >= 1; L--) cum[L] = cum[L + 1] + hist[L];

    if (cum[1] == 0) {
        printf("\nchain-length histogram: nothing examined yet\n");
        return;
    }

    printf("\nchain-length histogram\n");
    printf("  each prime p examined is counted under the length of the chain that\n"
           "  starts at p, so a maximal chain of length 5 also contributes a 4, a 3,\n"
           "  a 2 and a 1 at its later terms\n\n");
    printf("    len              chains           >= len       ratio\n");

    for (int L = 1; L <= top; L++) {
        char cb[32], kb[32];
        fmt_grouped(hist[L], cb, sizeof cb);
        fmt_grouped(cum[L],  kb, sizeof kb);
        char ratio[24];
        if (L == 1 || cum[L - 1] == 0) snprintf(ratio, sizeof ratio, "%10s", "--");
        else snprintf(ratio, sizeof ratio, "%10.3g",
                      (double)cum[L] / (double)cum[L - 1]);
        printf("  %5d  %18s %16s  %s%s\n", L, cb, kb, ratio,
               L == want_len ? "   <-- target" : "");
    }

    /* rough extrapolation: how much more work for one extra link? */
    int hi = 0;
    for (int L = top; L >= 2; L--) if (cum[L] > 0) { hi = L; break; }
    if (hi >= 3 && cum[hi - 1] > 0) {
        double r = (double)cum[hi] / (double)cum[hi - 1];
        if (r > 0 && r < 1)
            printf("\n  ratio at the deep end is about %.4g, so each extra link costs\n"
                   "  roughly %.4g times more primes to test\n", r, 1.0 / r);
    }
}

/* ------------------------------------------------------------------ */

static void usage(const char *argv0) {
    fprintf(stderr,
"Search for chains of primes where each term is the previous term plus its\n"
"own decimal digit sum.\n"
"\n"
"Usage: %s [options]\n"
"  -l, --length N     minimum chain length to report        (default 6)\n"
"  -s, --start N      lowest starting prime to consider     (default 2)\n"
"  -e, --end N        highest starting prime to consider    (default 1e9)\n"
"  -t, --threads N    worker threads                        (default: all cores)\n"
"  -c, --chunk N      range size per work unit              (default 1e7)\n"
"  -a, --all          report mid-chain hits too (default: maximal chains only)\n"
"  -F, --first        find the earliest chain of length >= --length: hits are\n"
"                     printed as found, work above the best hit is dropped, and\n"
"                     the sweep below it finishes so the smallest start is\n"
"                     confirmed, not just the first one stumbled upon\n"
"  -k, --checkpoint F checkpoint to file F, and resume from it if it exists\n"
"  -i, --interval S   seconds between checkpoint writes      (default 60)\n"
"      --restart      ignore any existing checkpoint and start the range over\n"
"  -q, --quiet        no progress line, results only\n"
"  -h, --help         this message\n"
"\n"
"Number formats for --start/--end/--chunk/--threads/--interval:\n"
"  plain      1234        1_000_000     1,000,000\n"
"  scientific 1e9         2.5e8         3.2e14\n"
"  suffixed   900K        12M   4B(=4G)   2.5T   3P   1E\n"
"             K thousand, M million, B or G billion, T trillion,\n"
"             P quadrillion, E quintillion.  Case does not matter.\n"
"             Fractions are fine: 1.5B = 1500000000.\n"
"             Trailing E is exa, but 3E9 still reads as scientific notation.\n"
"             Max is the 64-bit limit, ~18.4E.\n"
"\n"
"Checkpoints:\n"
"  With -k the run saves progress every --interval seconds, on Ctrl-C, and at\n"
"  the end. Re-running the same command with the same -k file continues where\n"
"  it stopped -- --start/--end/--length/--all must match, and the chunk size\n"
"  stored in the file wins over -c. Ctrl-C twice exits without saving.\n"
"  Chunks that were mid-flight when a periodic checkpoint landed are redone,\n"
"  so a resumed run can re-report a hit it already printed.\n"
"\n"
"Statistics:\n"
"  A histogram of chain lengths is printed when the run ends, including when it\n"
"  is interrupted or cut short by --first.\n"
"\n"
"Examples: %s -l 7 -s 1e8 -e 1e11 -t 16\n"
"          %s -l 8 -s 1e10 -e 3.2e14 -k len8.ckpt\n"
"          %s -l 9 -s 2 -e 1e12 --first\n",
        argv0, argv0, argv0, argv0);
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv) {
    int want_len = 6;
    uint64_t start = 2, end = 1000000000ULL;
    uint64_t chunk = 10000000ULL;
    int threads = 0, maximal_only = 1, quiet = 0;
    int stop_on_first = 0, restart = 0, chunk_given = 0;
    const char *ck_path = NULL;
    double ck_interval = 60.0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        int has_next = (i + 1 < argc);
        #define NEXT() (has_next ? argv[++i] : (usage(argv[0]), exit(2), ""))
        if (!strcmp(a, "-l") || !strcmp(a, "--length"))       want_len = (int)parse_num(NEXT(), a);
        else if (!strcmp(a, "-s") || !strcmp(a, "--start"))   start = parse_num(NEXT(), a);
        else if (!strcmp(a, "-e") || !strcmp(a, "--end"))     end = parse_num(NEXT(), a);
        else if (!strcmp(a, "-t") || !strcmp(a, "--threads")) threads = (int)parse_num(NEXT(), a);
        else if (!strcmp(a, "-c") || !strcmp(a, "--chunk"))  { chunk = parse_num(NEXT(), a); chunk_given = 1; }
        else if (!strcmp(a, "-a") || !strcmp(a, "--all"))     maximal_only = 0;
        else if (!strcmp(a, "-F") || !strcmp(a, "--first"))   stop_on_first = 1;
        else if (!strcmp(a, "-k") || !strcmp(a, "--checkpoint")) ck_path = NEXT();
        else if (!strcmp(a, "-i") || !strcmp(a, "--interval")) ck_interval = (double)parse_num(NEXT(), a);
        else if (!strcmp(a, "--restart"))                     restart = 1;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet"))   quiet = 1;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help"))  { usage(argv[0]); return 0; }
        else { fprintf(stderr, "unknown option: %s\n\n", a); usage(argv[0]); return 2; }
        #undef NEXT
    }

    if (want_len < 1) want_len = 1;
    if (start < 2) start = 2;
    if (end < start) { fprintf(stderr, "error: --end must be >= --start\n"); return 2; }
    if (chunk < 1000) chunk = 1000;
    if (ck_interval < 1) ck_interval = 1;
    if (!ck_path && restart)
        fprintf(stderr, "note: --restart has no effect without --checkpoint\n");

    if (isatty(1)) { CLR = "\r\033[K"; NL = ""; }

    /* ---- load a checkpoint before sizing chunks, so the chunk grid matches ---- */
    ckpt_t ck;
    memset(&ck, 0, sizeof ck);
    int resuming = 0;

    if (ck_path && !restart && ck_read(ck_path, &ck)) {
        if (ck.version != CK_VERSION) {
            fprintf(stderr,
                "error: checkpoint %s is version %d, this binary writes version %d.\n"
                "       Delete it or use --restart to start the range over.\n",
                ck_path, ck.version, CK_VERSION);
            return 2;
        }
        if (ck.start != start || ck.end != end ||
            ck.want_len != want_len || ck.maximal != maximal_only) {
            fprintf(stderr,
                "error: checkpoint %s was written for a different search:\n"
                "         checkpoint: --start %" PRIu64 " --end %" PRIu64
                " --length %d %s\n"
                "         this run  : --start %" PRIu64 " --end %" PRIu64
                " --length %d %s\n"
                "       Use a different -k file, or --restart to overwrite it.\n",
                ck_path, ck.start, ck.end, ck.want_len, ck.maximal ? "" : "--all",
                start, end, want_len, maximal_only ? "" : "--all");
            return 2;
        }
        if (chunk_given && chunk != ck.chunk)
            fprintf(stderr, "note: using chunk size %" PRIu64 " from the checkpoint, "
                            "not the %" PRIu64 " you asked for\n", ck.chunk, chunk);
        chunk = ck.chunk;
        resuming = 1;
    }

    ds_init();
#ifdef NO_PRIMESIEVE
    /* chunks sieve SEG_MARGIN past their end for the chain-walk bitmap */
    base_init(end > UINT64_MAX - SEG_MARGIN ? UINT64_MAX : end + SEG_MARGIN);
#endif

    int threads_requested = threads;
    int have_openmp = 0;

#ifdef _OPENMP
    have_openmp = 1;
    (void)threads_requested;
    if (threads > 0) omp_set_num_threads(threads);
    threads = omp_get_max_threads();
#else
    threads = 1;
    if (threads_requested > 1) {
        fprintf(stderr,
"\n"
"*** WARNING: --threads %d ignored -- this binary was built WITHOUT OpenMP. ***\n"
"    The parallel pragmas were compiled out, so the search will run on one\n"
"    core. Rebuild with OpenMP enabled:\n"
"\n"
"      Linux / gcc:   cc -O3 -march=native -fopenmp -o dschain dschain.c -lprimesieve\n"
"\n"
"      macOS / clang: brew install libomp\n"
"                     clang -O3 -Xpreprocessor -fopenmp \\\n"
"                       -I\"$(brew --prefix libomp)/include\" \\\n"
"                       -L\"$(brew --prefix libomp)/lib\" -lomp \\\n"
"                       -o dschain dschain.c -lprimesieve\n"
"\n"
"      macOS / gcc:   brew install gcc && gcc-14 -O3 -fopenmp -o dschain dschain.c -lprimesieve\n"
"\n"
"    Apple's stock clang rejects a bare -fopenmp, which is the usual cause.\n"
"    Or just run `make` -- the Makefile detects all of this.\n"
"\n", threads_requested);
        fflush(stderr);
    }
#endif

    uint64_t total = end - start + 1;
    uint64_t nchunks = (total + chunk - 1) / chunk;

    if (!resuming) {
        /* keep every thread busy: aim for >= 8 chunks per thread */
        if (nchunks < (uint64_t)threads * 8) {
            chunk = total / ((uint64_t)threads * 8) + 1;
            if (chunk < 1000) chunk = 1000;
            nchunks = (total + chunk - 1) / chunk;
        }
        /* and keep the chunk map (1 bit per chunk) to a sane size */
        const uint64_t MAX_CHUNKS = 100000000ULL;   /* 12.5 MB of bitmap */
        if (nchunks > MAX_CHUNKS) {
            chunk = total / MAX_CHUNKS + 1;
            nchunks = (total + chunk - 1) / chunk;
        }
    } else if (nchunks != ck.nchunks) {
        fprintf(stderr, "error: checkpoint chunk grid does not match "
                        "(%" PRIu64 " vs %" PRIu64 " chunks)\n", ck.nchunks, nchunks);
        return 2;
    }

    ck_alloc(nchunks);

    double t0 = now_sec();
    double resumed_elapsed = 0;

    if (resuming) {
        if (ck.complete) {
            printf("checkpoint %s says this range is already finished:\n", ck_path);
            printf("  primes examined  : %" PRIu64 "\n", ck.primes_tested);
            printf("  chains >= %-6d : %" PRIu64 "\n", ck.want_len, ck.found);
            printf("  longest seen     : %d (starting at %" PRIu64 ")\n",
                   ck.longest, ck.longest_at);
            print_histogram(ck.hist, ck.want_len, ck.longest);
            printf("\nUse --restart to run it again.\n");
            return 0;
        }
        ck_watermark = ck.watermark > nchunks ? nchunks : ck.watermark;
        for (size_t i = 0; i < ck.n_done_above; i++) {
            uint64_t c = ck.done_above[i];
            if (c < nchunks) ck_bits[c >> 6] |= 1ULL << (c & 63);
        }
        while (ck_watermark < nchunks &&
               ((ck_bits[ck_watermark >> 6] >> (ck_watermark & 63)) & 1ULL))
            ck_watermark++;
        for (size_t i = 0; i < ck.n_partials && ck_nparts < MAX_PARTIALS; i++) {
            if (ck.partials[i].chunk >= nchunks) continue;
            if (ck_is_done(ck.partials[i].chunk)) continue;
            ck_part[ck_nparts++] = ck.partials[i];
        }
        memcpy(g_hist, ck.hist, sizeof g_hist);
        g_primes_tested = ck.primes_tested;
        g_found         = ck.found;
        g_longest       = ck.longest;
        g_longest_at    = ck.longest_at;
        resumed_elapsed = ck.elapsed;

        uint64_t done_chunks = ck_watermark;
        for (uint64_t c = ck_watermark; c < nchunks; c++)
            if ((ck_bits[c >> 6] >> (c & 63)) & 1ULL) done_chunks++;
        g_done_range = done_chunks >= nchunks ? total : done_chunks * chunk;
        if (g_done_range > total) g_done_range = total;
    }
    free(ck.done_above);
    free(ck.partials);

    printf("dschain: p -> p + digitsum(p)\n");
    {
        char sb2[32], eb2[32], tb2[32];
        fmt_grouped(start, sb2, sizeof sb2);
        fmt_grouped(end,   eb2, sizeof eb2);
        fmt_grouped(total, tb2, sizeof tb2);
        printf("  range     : %s .. %s  (%s integer%s)\n",
               sb2, eb2, tb2, total == 1 ? "" : "s");
    }
    printf("  min length: %d%s\n", want_len,
           stop_on_first ? "   (finding the earliest hit)" : "");
    printf("  threads   : %d%s\n", threads,
           have_openmp ? "" : "   <-- built WITHOUT OpenMP, single-core only");
    printf("  chunks    : %" PRIu64 " x %" PRIu64 "\n", nchunks, chunk);
#ifdef NO_PRIMESIEVE
    printf("  sieve     : built-in segmented (%zu base primes)\n", n_base);
#else
    printf("  sieve     : primesieve\n");
#endif
    printf("  maximal   : %s\n", maximal_only ? "yes" : "no (report all)");
    if (ck_path) {
        printf("  checkpoint: %s every %.0fs\n", ck_path, ck_interval);
        if (resuming) {
            char rb[32];
            fmt_grouped(g_done_range, rb, sizeof rb);
            printf("  resuming  : %.2f%% done (%s integers), %" PRIu64
                   " primes already examined, %.0fs of prior runtime\n",
                   100.0 * (double)g_done_range / (double)total, rb,
                   g_primes_tested, resumed_elapsed);
        }
    }
    printf("\n");
    fflush(stdout);

    /* ---- catch Ctrl-C / SIGTERM ---- */
    {
        struct sigaction sa;
        memset(&sa, 0, sizeof sa);
        sa.sa_handler = handle_signal;
        sigemptyset(&sa.sa_mask);
        sigaction(SIGINT,  &sa, NULL);
        sigaction(SIGTERM, &sa, NULL);
    }

    double last_print = t0, last_ckpt = t0;

    /* monotonic: chunks are handed out in increasing order.  A plain
     * dynamic schedule is nonmonotonic since OpenMP 5.0, and LLVM's libomp
     * then gives each thread its own contiguous 1/T of the range up front
     * (work stealing), so the watermark only advanced with thread 0 and
     * checkpoints had to list millions of finished chunks above it. */
    #pragma omp parallel for schedule(monotonic: dynamic, 1)
    for (uint64_t c = 0; c < nchunks; c++) {
        /* cheap lock-free skips first: after a stop, or on resume, the loop
         * has to race through the remaining indices without taking a lock */
        if (g_stop) continue;
        uint64_t wm;
        #pragma omp atomic read
        wm = ck_watermark;
        if (c < wm) continue;

        uint64_t lo = start + c * chunk;
        if (stop_on_first) {
            uint64_t fb;
            #pragma omp atomic read
            fb = g_first_bound;
            if (lo > fb) continue;      /* whole chunk is above a known hit */
        }

        int skip_this = 0;
        uint64_t skip_upto = 0;
        #pragma omp critical (output)
        {
            if (ck_is_done(c)) skip_this = 1;
            else skip_upto = ck_partial_get(c);
        }
        if (skip_this) continue;

        uint64_t hi = lo + chunk - 1;
        if (hi > end || hi < lo) hi = end;

        /* sieve a margin past the chunk so chain walks are bit tests */
        uint64_t seg_hi = hi + SEG_MARGIN;
        if (seg_hi < hi) seg_hi = UINT64_MAX;     /* overflow clamp */

        size_t np = 0;
        uint64_t *pv = gen_primes(lo, seg_hi, &np);

        seg_t seg;
        seg.base = lo | 1;
        uint64_t nidx  = seg_hi >= seg.base ? ((seg_hi - seg.base) >> 1) + 1 : 1;
        size_t   words = (size_t)((nidx + 63) / 64);
        uint64_t *bm = (uint64_t *)calloc(words, sizeof(uint64_t));
        if (!bm) { fprintf(stderr, "out of memory allocating segment bitmap\n"); exit(1); }
        for (size_t i = 0; i < np; i++) {
            uint64_t q = pv[i];
            if (q & 1) {
                uint64_t idx = (q - seg.base) >> 1;
                bm[idx >> 6] |= 1ULL << (idx & 63);
            }
        }
        seg.hi = seg_hi;
        seg.bm = bm;

        uint64_t local_hist[HIST_MAX + 1];
        memset(local_hist, 0, sizeof local_hist);
        uint64_t local_found = 0, local_tested = 0, last_p = 0;
        int local_longest = 0, aborted = 0;
        uint64_t local_longest_at = 0;

        for (size_t i = 0; i < np && pv[i] <= hi; i++) {   /* margin primes are not starts */
            if ((i & 0x3FF) == 0) {
                if (g_stop) { aborted = 1; break; }
                if (stop_on_first) {
                    uint64_t fb;
                    #pragma omp atomic read
                    fb = g_first_bound;
                    /* primes come in order, so everything from here on is
                     * above an already-reported hit -- drop the rest */
                    if (pv[i] > fb) { aborted = 1; break; }
                }
            }

            uint64_t p = pv[i];
            if (p <= skip_upto) continue;   /* done in an earlier run */

            int n;
            /* fast reject: an odd digit sum makes p+d even, so the chain
             * starting at p has length exactly 1 */
            if (want_len > 1 && (digitsum(p) & 1)) n = 1;
            else n = chain_length(p, NULL, 0, &seg);

            local_tested++;
            local_hist[n > HIST_MAX ? HIST_MAX : n]++;
            if (n > local_longest) { local_longest = n; local_longest_at = p; }

            if (n >= want_len) {
                if (maximal_only && !is_chain_start(p)) { last_p = p; continue; }

                uint64_t ch[128];
                int cap = n < 128 ? n : 128;
                chain_length(p, ch, cap, &seg);

                int report = 1;
                #pragma omp critical (output)
                {
                    if (stop_on_first) {
                        if (p < g_first_bound) {
                            #pragma omp atomic write
                            g_first_bound = p;
                        } else {
                            report = 0;     /* an earlier hit is already known */
                        }
                    }
                    if (report) {
                        printf("%s", CLR);       /* clear progress line */
                        printf("[len %d] ", n);
                        for (int k = 0; k < cap; k++)
                            printf("%s%" PRIu64 "(+%u)", k ? " -> " : "",
                                   ch[k], digitsum(ch[k]));
                        if (cap < n) printf(" -> ...");
                        printf("\n");
                        fflush(stdout);
                    }
                }
                if (report) local_found++;

                if (stop_on_first) {
                    /* the rest of this chunk is above the bound; drop it,
                     * but let other threads keep sweeping the range below */
                    last_p = p;             /* this prime is fully handled */
                    aborted = 1;
                    break;
                }
            }
            last_p = p;
        }

        free(bm);
        if (pv) free_primes(pv);

        #pragma omp critical (output)
        {
            g_primes_tested += local_tested;
            g_found += local_found;
            for (int k = 1; k <= HIST_MAX; k++) g_hist[k] += local_hist[k];
            if (local_longest > g_longest) {
                g_longest = local_longest;
                g_longest_at = local_longest_at;
            }

            if (!aborted) {
                ck_partial_clear(c);
                ck_mark_done(c);
                g_done_range += (hi - lo + 1);
            } else if (last_p) {
                ck_partial_set(c, last_p);
            }

            double t = now_sec();
            if (!quiet && (t - last_print > 0.5 || g_done_range >= total)) {
                last_print = t;
                double frac = (double)g_done_range / (double)total;
                double elapsed = t - t0 + resumed_elapsed;
                double run = t - t0;
                double eta = frac > 1e-9 ? elapsed * (1.0 / frac - 1.0) : -1;
                char eb[32], rb[32], pb[32], sb[32];
                fmt_duration(elapsed, eb, sizeof eb);
                fmt_duration(eta, rb, sizeof rb);
                fmt_count((double)g_primes_tested, pb, sizeof pb);
                fmt_count(run > 0 ? (double)g_done_range / elapsed : 0, sb, sizeof sb);
                printf("%s%6.2f%%  n=%" PRIu64 "  primes=%s  %s/s  "
                       "elapsed %s  eta %s  best=%d  hits=%" PRIu64,
                       CLR, frac * 100.0, start + g_done_range, pb, sb, eb, rb,
                       g_longest, g_found);
                printf("%s", NL);
                fflush(stdout);
            }

            if (ck_path && !g_stop && t - last_ckpt >= ck_interval) {
                ck_meta_t m = { start, end, chunk, want_len, maximal_only,
                                nchunks, t - t0 + resumed_elapsed, 0 };
                ck_write(ck_path, &m);
                last_ckpt = now_sec();
            }
        }
    }

    double t_end = now_sec();
    double elapsed = t_end - t0 + resumed_elapsed;
    int hit_found = stop_on_first && g_first_bound != UINT64_MAX;
    int complete = !g_stop && !hit_found;   /* --first leaves the range unswept above the hit */

    if (ck_path) {
        ck_meta_t m = { start, end, chunk, want_len, maximal_only,
                        nchunks, elapsed, complete };
        ck_write(ck_path, &m);
    }

    if (!quiet) printf("%s", CLR);
    printf("\n");
    if (g_stop_reason == STOP_SIGNAL) {
        printf("stopped by signal %d\n", (int)g_signo);
        if (hit_found)
            printf("earliest hit so far starts at %" PRIu64 " -- the range below it\n"
                   "was not fully swept, so it is NOT confirmed as the smallest\n",
                   g_first_bound);
    } else if (hit_found) {
        printf("earliest chain of length >= %d starts at %" PRIu64 "\n"
               "(confirmed smallest: every prime below it in range was examined)\n",
               want_len, g_first_bound);
    } else {
        printf("done\n");
    }

    printf("  elapsed          : %.2fs%s\n", elapsed,
           resumed_elapsed > 0 ? " (including earlier runs)" : "");
    printf("  primes examined  : %" PRIu64 "\n", g_primes_tested);
    printf("  chains >= %-6d : %" PRIu64 "\n", want_len, g_found);
    printf("  longest seen     : %d (starting at %" PRIu64 ")\n",
           g_longest, g_longest_at);
    printf("  range covered    : %.4f%%\n",
           100.0 * (double)g_done_range / (double)total);

    print_histogram(g_hist, want_len, g_longest);

    if (ck_path) {
        if (complete)
            printf("\ncheckpoint %s marked complete\n", ck_path);
        else
            printf("\ncheckpoint written to %s -- rerun the same command to continue\n",
                   ck_path);
        if (ck_dropped_parts)
            printf("note: more than %d interrupted chunks; some will be redone "
                   "from their start\n", MAX_PARTIALS);
    }
    fflush(stdout);

    free(ck_bits);
    return g_stop_reason == STOP_SIGNAL ? 130 : 0;
}
