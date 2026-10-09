/*
 * wheelchain.c -- constellation search for A016045.
 *
 * a(n) is the smallest prime starting n+1 CONSECUTIVE primes with gaps
 * 2, 4, 6, ... 2n, i.e. members at the oblong offsets
 *
 *      p_i = p + i(i-1),   i = 1..k        (0, 2, 6, 12, ... k(k-1))
 *
 * plus the requirement that no other prime lies inside the span.
 *
 * Unlike a plain segmented sieve, this program never enumerates primes.
 * It treats the problem as a prime k-tuplet search:
 *
 *   1. WHEEL.  Pick M = 2*3*5*...*P and precompute every residue r mod M
 *      for which no member p+off_j is divisible by any prime <= P.  For
 *      k=16 and P=29 only 1,995,840 of 6,469,693,230 residues survive --
 *      a 3200x reduction before any work happens.  Heads are then written
 *      p = r + t*M and we only ever look at those.
 *
 *   2. t-SPACE SIEVE.  For a sieving prime q and offset off_j,
 *
 *          r + off_j + t*M = 0 (mod q)   <=>   t = -(r+off_j)*M^-1 (mod q)
 *
 *      which is a plain stride-q progression in t.  So we run a segmented
 *      sieve over t, one residue class at a time.  The decisive point is
 *      that the bitmap now indexes CANDIDATES, not integers: cost scales
 *      with the 6e-4 of the range that survived the wheel, not with the
 *      range itself.
 *
 *      Since t0 is linear in the block base, we precompute
 *      D[q][j] = -off_j * M^-1 mod q once, and per block need only one
 *      modmul per prime plus an add per offset.
 *
 *   3. Survivors get a Montgomery Miller-Rabin on each member (cheapest
 *      rejection first), and only a full k-tuplet is then checked for
 *      consecutiveness -- that last test is a ~5% filter and worthless
 *      as an early prune.
 *
 * Blocks are processed in ascending t so coverage stays contiguous: you
 * need the SMALLEST occurrence, so a gap in the scanned range is a
 * correctness bug, not just lost work.
 *
 * Build:
 *      cc -O3 -march=native -pthread -o wheelchain wheelchain.c -lm
 *
 * Usage:
 *      ./wheelchain [-w P] [-s B] <chain_length> <start> <end> [threads]
 *
 *      -w P   wheel prime limit (default: chosen from the range)
 *      -s B   sieve candidates with primes up to B (default 200)
 *
 *      -s is worth measuring on your own hardware: the optimum is where
 *      the marginal marking cost equals a Montgomery Miller-Rabin, and
 *      it sits FAR lower than intuition suggests -- around 150-200 here,
 *      where over 10^7 candidates per 10^13 still reach the MR stage.
 *
 * Examples:
 *      ./wheelchain 13 101951000000000 101952000000000     -> a(12)
 *      ./wheelchain 16 221860944705726407 5E18 16          -> hunt a(15)
 */

#define _GNU_SOURCE
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define KMAX        64
#define TBLOCK      (1u << 21)      /* t values per block (256 KiB bitmap) */
#define RCHUNK      64              /* residues per work unit             */
#define MAXWHEEL    41

/* ------------------------------------------------------------------ */
/* 64-bit Montgomery Miller-Rabin                                      */

static inline uint64_t mont_ninv(uint64_t n)   /* n^-1 mod 2^64, n odd */
{
    uint64_t x = n;
    for (int i = 0; i < 6; i++)
        x *= 2 - n * x;
    return x;
}

typedef struct { uint64_t n, ninv, r1, r2; } mont;

static void mont_init(mont *m, uint64_t n)
{
    m->n    = n;
    m->ninv = mont_ninv(n);
    m->r1   = (uint64_t)(~(uint64_t)0 % n) + 1;               /* 2^64 mod n */
    m->r2   = (uint64_t)(((__uint128_t)m->r1 * m->r1) % n);   /* 2^128 mod n */
}

static inline uint64_t mont_redc(const mont *m, __uint128_t t)
{
    uint64_t q = (uint64_t)t * (uint64_t)(0 - m->ninv);
    uint64_t y = (uint64_t)((t + (__uint128_t)q * m->n) >> 64);
    return y >= m->n ? y - m->n : y;
}

static inline uint64_t mont_mul(const mont *m, uint64_t a, uint64_t b)
{
    return mont_redc(m, (__uint128_t)a * b);
}

static inline uint64_t mont_in(const mont *m, uint64_t a)
{
    return mont_mul(m, a % m->n, m->r2);
}

static bool mr_witness(const mont *m, uint64_t a, uint64_t d, int s)
{
    uint64_t one = m->r1, x = one, base = mont_in(m, a);
    if (base == 0)
        return true;                       /* a = 0 mod n: no information */

    for (uint64_t e = d; e; e >>= 1) {
        if (e & 1)
            x = mont_mul(m, x, base);
        base = mont_mul(m, base, base);
    }
    uint64_t nm1 = m->n - one;              /* -1 in Montgomery form */
    if (x == one || x == nm1)
        return true;
    for (int i = 1; i < s; i++) {
        x = mont_mul(m, x, x);
        if (x == nm1)
            return true;
    }
    return false;
}

/* deterministic for n < 3.3e24 */
static const uint64_t MR_BASES[] = {
    2, 325, 9375, 28178, 450775, 9780504, 1795265022
};

static bool is_prime_u64(uint64_t n)
{
    if (n < 2)
        return false;
    for (uint64_t p = 2; p < 40; p++) {
        if (p * p > n)
            return true;
        if (n % p == 0)
            return n == p;
    }
    uint64_t d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }

    mont m;
    mont_init(&m, n);
    for (size_t i = 0; i < sizeof MR_BASES / sizeof *MR_BASES; i++)
        if (!mr_witness(&m, MR_BASES[i], d, s))
            return false;
    return true;
}

/* cheap first filter: strong probable prime base 2 only */
static bool sprp2(uint64_t n)
{
    if (!(n & 1))
        return n == 2;
    uint64_t d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    mont m;
    mont_init(&m, n);
    return mr_witness(&m, 2, d, s);
}

/* ------------------------------------------------------------------ */
/* small helpers                                                       */

/*
 * Barrett reduction: x mod q for q < 2^32, using mag = floor(2^64/q).
 * The sieve does two reductions per prime per residue per block and a
 * hardware 64-bit divide costs ~25ns here -- enough that it dominated
 * everything else at the larger wheels.
 */
static inline uint64_t bmod(uint64_t x, uint64_t q, uint64_t mag)
{
    uint64_t r = x - (uint64_t)(((__uint128_t)x * mag) >> 64) * q;
    while (r >= q)
        r -= q;
    return r;
}

static uint64_t modinv_u64(uint64_t a, uint64_t m)
{
    int64_t t = 0, nt = 1;
    int64_t r = (int64_t)m, nr = (int64_t)(a % m);
    while (nr) {
        int64_t qq = r / nr;
        int64_t tmp = t - qq * nt; t = nt; nt = tmp;
        tmp = r - qq * nr;         r = nr; nr = tmp;
    }
    if (r > 1)
        return 0;
    if (t < 0)
        t += (int64_t)m;
    return (uint64_t)t;
}

static uint32_t *simple_primes(uint32_t limit, size_t *count)
{
    uint8_t *bm = calloc(limit + 1, 1);
    for (uint64_t i = 2; i * i <= limit; i++)
        if (!bm[i])
            for (uint64_t j = i * i; j <= limit; j += i)
                bm[j] = 1;
    size_t n = 0;
    for (uint32_t i = 2; i <= limit; i++)
        n += !bm[i];
    uint32_t *out = malloc(n * sizeof *out);
    n = 0;
    for (uint32_t i = 2; i <= limit; i++)
        if (!bm[i])
            out[n++] = i;
    free(bm);
    *count = n;
    return out;
}

static int parse_u64(const char *arg, uint64_t *out)
{
    char buf[128];
    size_t n = 0;
    for (const char *s = arg; *s; s++) {
        if (*s == ',' || *s == '_' || *s == ' ')
            continue;
        if (n + 1 >= sizeof buf)
            return -1;
        buf[n++] = *s;
    }
    buf[n] = '\0';
    if (!n)
        return -1;
    char *end;
    errno = 0;
    long double v = strtold(buf, &end);
    if (end == buf || v < 0.0L)
        return -1;
    long double mul = 1.0L;
    if (*end) {
        switch (toupper((unsigned char)*end)) {
        case 'K': mul = 1e3L;  break;
        case 'M': mul = 1e6L;  break;
        case 'B':
        case 'G': mul = 1e9L;  break;
        case 'T': mul = 1e12L; break;
        case 'P': mul = 1e15L; break;
        case 'E': mul = 1e18L; break;
        default:  return -1;
        }
        if (*++end)
            return -1;
    }
    long double r = v * mul;
    if (r >= 1.8e19L)
        return -1;
    *out = (uint64_t)(r + 0.5L);
    return 0;
}

/* ------------------------------------------------------------------ */
/* globals                                                             */

static int      g_k, g_nint;
static uint64_t g_off[KMAX];          /* member offsets, off[0] = 0     */
static uint64_t g_int[KMAX * KMAX];   /* interior offsets (must be composite) */
static uint64_t g_span;

static uint64_t  g_M;                 /* wheel modulus                  */
static uint64_t *g_res;               /* admissible residues mod M      */
static size_t    g_nres;

static uint32_t *g_sq;                /* sieving primes                 */
static uint32_t *g_sinv;              /* M^-1 mod q                     */
static uint64_t *g_smagic;            /* floor(2^64/q), for Barrett      */
static uint32_t *g_D;                 /* D[i*k + j] = -off_j*M^-1 mod q */
static size_t    g_nsq;

static uint64_t g_lo, g_hi, g_tlo, g_tcount, g_tblock;
static uint64_t g_nblocks, g_units_per_block, g_nunits;
static uint64_t g_next_unit, g_done_units, g_found, g_tested;
static int      g_progress, g_quiet, g_tty;
static double   g_interval = 0.0;
static double   g_t0, g_last_print;
static long     g_nthreads;
static uint64_t *g_inflight;          /* unit each thread is working on */
static pthread_mutex_t g_out = PTHREAD_MUTEX_INITIALIZER;

/* ------------------------------------------------------------------ */
/* wheel construction                                                  */

static int admissible_w(uint64_t p)
{
    /* number of distinct offsets mod p */
    uint8_t seen[MAXWHEEL > 64 ? MAXWHEEL : 64];
    if (p > 64) {
        int w = 0;
        for (int i = 0; i < g_k; i++) {
            int dup = 0;
            for (int j = 0; j < i; j++)
                if (g_off[i] % p == g_off[j] % p) { dup = 1; break; }
            w += !dup;
        }
        return w;
    }
    memset(seen, 0, sizeof seen);
    int w = 0;
    for (int i = 0; i < g_k; i++) {
        uint64_t rr = g_off[i] % p;
        if (!seen[rr]) { seen[rr] = 1; w++; }
    }
    return w;
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

static void build_wheel(uint32_t wheel_limit)
{
    size_t np;
    uint32_t *wp = simple_primes(wheel_limit, &np);

    g_M = 1;
    g_nres = 1;
    g_res = malloc(sizeof *g_res);
    g_res[0] = 0;

    for (size_t i = 0; i < np; i++) {
        uint64_t p = wp[i];

        /* allowed residues a mod p: no member divisible by p */
        uint8_t bad[MAXWHEEL + 1];
        memset(bad, 0, sizeof bad);
        for (int j = 0; j < g_k; j++)
            bad[(p - g_off[j] % p) % p] = 1;
        uint64_t allow[MAXWHEEL + 1];
        int nallow = 0;
        for (uint64_t a = 0; a < p; a++)
            if (!bad[a])
                allow[nallow++] = a;
        if (!nallow) {
            fprintf(stderr, "wheelchain: pattern inadmissible mod %" PRIu64 "\n", p);
            exit(1);
        }

        uint64_t newM = g_M * p;
        uint64_t Minv = modinv_u64(g_M % p, p);
        size_t newn = g_nres * (size_t)nallow;
        uint64_t *nr = malloc(newn * sizeof *nr);
        if (!nr) {
            fprintf(stderr, "wheelchain: wheel too large (%zu residues)\n", newn);
            exit(1);
        }
        size_t c = 0;
        for (size_t x = 0; x < g_nres; x++) {
            uint64_t r = g_res[x];
            uint64_t rp = r % p;
            for (int a = 0; a < nallow; a++) {
                uint64_t t = ((allow[a] + p - rp) % p) * Minv % p;
                nr[c++] = r + t * g_M;
            }
        }
        free(g_res);
        g_res = nr;
        g_nres = newn;
        g_M = newM;
    }
    free(wp);

    /* ascending, so heads come out roughly ordered within a block */
    qsort(g_res, g_nres, sizeof *g_res, cmp_u64);
}

/* ------------------------------------------------------------------ */
/* full verification of one surviving head                             */

static int verify(uint64_t p, uint64_t *members, int *maxlen)
{
    /* members first: cheap base-2 test, then full MR */
    for (int j = 0; j < g_k; j++)
        if (!sprp2(p + g_off[j]))
            return 0;
    for (int j = 0; j < g_k; j++)
        if (!is_prime_u64(p + g_off[j]))
            return 0;

    /* consecutiveness: every other odd value in the span must be composite */
    for (int i = 0; i < g_nint; i++)
        if (is_prime_u64(p + g_int[i]))
            return 0;

    for (int j = 0; j < g_k; j++)
        members[j] = p + g_off[j];

    /* how much further does it run? */
    int len = g_k;
    uint64_t prev = p + g_off[g_k - 1];
    for (int n = g_k + 1;; n++) {
        uint64_t nxt = p + (uint64_t)n * (n - 1);
        if (nxt < prev)
            break;                              /* overflow guard */
        int ok = is_prime_u64(nxt);
        for (uint64_t d = prev + 2; ok && d < nxt; d += 2)
            if (is_prime_u64(d))
                ok = 0;
        if (!ok)
            break;
        len++;
        prev = nxt;
    }
    *maxlen = len;
    return 1;
}

/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/* progress reporting                                                  */

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* "3d 04h 17m" / "17m 09s" / "9.4s" */
static void fmt_dur(double s, char *buf, size_t n)
{
    if (s < 0 || s != s) { snprintf(buf, n, "--"); return; }
    if (s < 60)            { snprintf(buf, n, "%.1fs", s); return; }
    unsigned long t = (unsigned long)(s + 0.5);
    unsigned d = t / 86400, h = (t % 86400) / 3600,
             m = (t % 3600) / 60, sec = t % 60;
    if (d)      snprintf(buf, n, "%ud %02uh %02um", d, h, m);
    else if (h) snprintf(buf, n, "%uh %02um %02us", h, m, sec);
    else        snprintf(buf, n, "%um %02us", m, sec);
}

/*
 * The frontier is the lowest unit still in flight: everything below it is
 * genuinely finished.  That is the number worth checkpointing -- work
 * units are handed out in ascending order, but they do not RETIRE in
 * order, so "units done" alone would overstate how far the search has
 * actually been verified.
 */
static uint64_t frontier_unit(void)
{
    uint64_t f = __atomic_load_n(&g_next_unit, __ATOMIC_RELAXED);
    if (f > g_nunits)
        f = g_nunits;
    for (long i = 0; i < g_nthreads; i++) {
        uint64_t v = __atomic_load_n(&g_inflight[i], __ATOMIC_RELAXED);
        if (v < f)
            f = v;
    }
    return f;
}

static void report_progress(int final)
{
    if (g_quiet)
        return;

    double t = now_sec();
    if (!final) {
        pthread_mutex_lock(&g_out);
        if (t - g_last_print < g_interval) {
            pthread_mutex_unlock(&g_out);
            return;
        }
        g_last_print = t;
        pthread_mutex_unlock(&g_out);
    }

    uint64_t done = __atomic_load_n(&g_done_units, __ATOMIC_RELAXED);
    double el   = t - g_t0;
    double frac = g_nunits ? (double)done / (double)g_nunits : 1.0;

    /* Exact head value below which the search is genuinely finished --
     * this is the number to record if the run is interrupted. */
    uint64_t fu  = final ? g_nunits : frontier_unit();
    uint64_t blk = fu / g_units_per_block;
    uint64_t tf  = g_tlo + blk * g_tblock;
    uint64_t head = (tf > g_hi / g_M) ? g_hi : tf * g_M;
    if (head < g_lo)  head = g_lo;
    if (head > g_hi)  head = g_hi;

    double rate = el > 0 ? (double)done / el : 0;          /* units/sec  */
    double eta  = rate > 0 ? (g_nunits - done) / rate : -1;
    double nps  = el > 0 ? (double)(g_hi - g_lo) * frac / el : 0;

    char eb[32], elb[32];
    fmt_dur(eta, eb, sizeof eb);
    fmt_dur(el, elb, sizeof elb);

    pthread_mutex_lock(&g_out);
    fprintf(stderr,
            "%s%6.2f%% | done<%" PRIu64 " | %.2e num/s | eta %s | elapsed %s | "
            "tested %" PRIu64 " found %" PRIu64 "  %s",
            g_tty ? "\r" : "", 100.0 * frac, head, nps, eb, elb,
            __atomic_load_n(&g_tested, __ATOMIC_RELAXED),
            __atomic_load_n(&g_found, __ATOMIC_RELAXED),
            g_tty ? "  " : "\n");
    fflush(stderr);
    pthread_mutex_unlock(&g_out);
}

/* ------------------------------------------------------------------ */

static void *worker(void *arg)
{
    int id = (int)(intptr_t)arg;
    uint64_t *bm = malloc((TBLOCK / 64 + 2) * sizeof *bm);
    uint64_t members[KMAX];
    if (!bm)
        return NULL;

    for (;;) {
        uint64_t u = __atomic_fetch_add(&g_next_unit, 1, __ATOMIC_RELAXED);
        if (u >= g_nunits)
            break;
        __atomic_store_n(&g_inflight[id], u, __ATOMIC_RELAXED);

        uint64_t b = u / g_units_per_block;
        uint64_t c = u % g_units_per_block;
        size_t r0 = (size_t)c * RCHUNK;
        size_t r1 = r0 + RCHUNK < g_nres ? r0 + RCHUNK : g_nres;

        uint64_t tstart = g_tlo + b * g_tblock;
        uint64_t len = g_tcount - b * g_tblock;
        if (len > g_tblock)
            len = g_tblock;
        size_t words = (len + 63) / 64;

        for (size_t ri = r0; ri < r1; ri++) {
            uint64_t X = g_res[ri] + tstart * g_M;   /* head at t' = 0 */

            memset(bm, 0, words * sizeof *bm);

            for (size_t i = 0; i < g_nsq; i++) {
                uint64_t q = g_sq[i];
                uint64_t inv = g_sinv[i];
                uint64_t mag = g_smagic[i];
                uint64_t A = bmod(X, q, mag);
                uint64_t u0 = bmod((q - A) % q * inv, q, mag); /* -X*M^-1 mod q */
                const uint32_t *D = g_D + i * g_k;
                for (int j = 0; j < g_k; j++) {
                    uint64_t t = u0 + D[j];
                    if (t >= q)
                        t -= q;
                    /* q divides the member because the member IS q: that
                     * one is prime, so it must not be struck.  Only ever
                     * possible for tiny ranges, but it silently breaks
                     * validation runs against the known terms. */
                    uint64_t base = X + g_off[j];
                    if (base <= q && (q - base) % g_M == 0) {
                        uint64_t tq = (q - base) / g_M;
                        if (t == tq)
                            t += q;
                    }
                    for (; t < len; t += q)
                        bm[t >> 6] |= 1ULL << (t & 63);
                }
            }

            for (size_t wi = 0; wi < words; wi++) {
                uint64_t bits = ~bm[wi];
                uint64_t avail = len - wi * 64;
                if (avail < 64)
                    bits &= (1ULL << avail) - 1;
                while (bits) {
                    int bpos = __builtin_ctzll(bits);
                    bits &= bits - 1;
                    uint64_t p = X + (wi * 64 + bpos) * g_M;
                    if (p < g_lo || p > g_hi)
                        continue;
                    __atomic_fetch_add(&g_tested, 1, __ATOMIC_RELAXED);
                    int maxlen;
                    if (!verify(p, members, &maxlen))
                        continue;
                    __atomic_fetch_add(&g_found, 1, __ATOMIC_RELAXED);
                    pthread_mutex_lock(&g_out);
                    printf("%2d: %" PRIu64, maxlen, p);
                    for (int j = 1; j < g_k; j++)
                        printf(" %" PRIu64, members[j]);
                    if (maxlen > g_k)
                        printf(" ... (extends to %d)", maxlen);
                    putchar('\n');
                    fflush(stdout);
                    pthread_mutex_unlock(&g_out);
                }
            }
        }

        __atomic_fetch_add(&g_done_units, 1, __ATOMIC_RELAXED);
        __atomic_store_n(&g_inflight[id], UINT64_MAX, __ATOMIC_RELAXED);
        report_progress(0);
    }
    __atomic_store_n(&g_inflight[id], UINT64_MAX, __ATOMIC_RELAXED);
    free(bm);
    return NULL;
}

/* ------------------------------------------------------------------ */

static void usage(const char *me)
{
    fprintf(stderr,
        "usage: %s [-w P] [-s B] [-p SECS] [-q] <chain_length> <start> <end> [threads]\n"
        "  -w P     wheel prime limit (default: auto from range)\n"
        "  -s B     sieve candidates with primes up to B (default 200)\n"
        "  -p SECS  progress interval (default 5 on a tty, 60 otherwise)\n"
        "  -q       no progress output\n"
        "  start/end accept K/M/B/T/P/E suffixes\n", me);
}

int main(int argc, char **argv)
{
    uint32_t wheel_limit = 0, sieve_bound = 200;
    int ai = 1;

    while (ai < argc && argv[ai][0] == '-' && argv[ai][1]) {
        if (!strcmp(argv[ai], "-w") && ai + 1 < argc)
            wheel_limit = (uint32_t)strtoul(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "-s") && ai + 1 < argc)
            sieve_bound = (uint32_t)strtoul(argv[++ai], NULL, 10);
        else if (!strcmp(argv[ai], "-p") && ai + 1 < argc)
            g_interval = strtod(argv[++ai], NULL);
        else if (!strcmp(argv[ai], "-q"))
            g_quiet = 1;
        else { usage(argv[0]); return 2; }
        ai++;
    }
    if (argc - ai < 3 || argc - ai > 4) {
        usage(argv[0]);
        return 2;
    }

    long k = strtol(argv[ai], NULL, 10);
    if (k < 2 || k > KMAX) {
        fprintf(stderr, "wheelchain: chain_length must be 2..%d\n", KMAX);
        return 2;
    }
    g_k = (int)k;
    for (int n = 1; n <= g_k; n++)
        g_off[n - 1] = (uint64_t)n * (n - 1);
    g_span = g_off[g_k - 1];

    g_nint = 0;
    for (uint64_t d = 2; d < g_span; d += 2) {
        int is_member = 0;
        for (int j = 0; j < g_k; j++)
            if (g_off[j] == d) { is_member = 1; break; }
        if (!is_member)
            g_int[g_nint++] = d;
    }

    if (parse_u64(argv[ai + 1], &g_lo) || parse_u64(argv[ai + 2], &g_hi)) {
        fprintf(stderr, "wheelchain: bad range\n");
        return 2;
    }
    if (g_lo < 3)
        g_lo = 3;
    if (g_hi < g_lo) {
        fprintf(stderr, "wheelchain: end < start\n");
        return 2;
    }

    long nthreads = (argc - ai == 4) ? strtol(argv[ai + 3], NULL, 10)
                                     : sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1)
        nthreads = 1;

    /* pick a wheel: big enough to matter, small enough that each residue
     * still gets a decent run of t values (and the table fits in RAM) */
    if (!wheel_limit) {
        uint64_t width = g_hi - g_lo;
        uint64_t M = 1;
        size_t nres = 1;
        wheel_limit = 7;
        size_t np;
        uint32_t *wp = simple_primes(MAXWHEEL, &np);
        for (size_t i = 0; i < np; i++) {
            uint64_t p = wp[i];
            int w = admissible_w(p);
            if (w >= (int)p)
                break;
            if (M > UINT64_MAX / p)
                break;
            uint64_t nM = M * p;
            size_t nn = nres * (size_t)(p - w);
            if (nn > 40000000u)
                break;
            if (nM > width / 4096 && i > 3)
                break;
            M = nM; nres = nn; wheel_limit = (uint32_t)p;
        }
        free(wp);
    }

    build_wheel(wheel_limit);

    /* sieving primes: everything above the wheel up to the bound */
    size_t np;
    uint32_t *all = simple_primes(sieve_bound < wheel_limit + 1
                                  ? wheel_limit + 1 : sieve_bound, &np);
    g_nsq = 0;
    for (size_t i = 0; i < np; i++)
        if (all[i] > wheel_limit && all[i] <= sieve_bound)
            g_nsq++;
    g_sq     = malloc(g_nsq * sizeof *g_sq);
    g_sinv   = malloc(g_nsq * sizeof *g_sinv);
    g_smagic = malloc(g_nsq * sizeof *g_smagic);
    g_D    = malloc(g_nsq * (size_t)g_k * sizeof *g_D);
    size_t c = 0;
    for (size_t i = 0; i < np; i++) {
        if (all[i] <= wheel_limit || all[i] > sieve_bound)
            continue;
        uint64_t q = all[i];
        uint64_t inv = modinv_u64(g_M % q, q);
        g_sq[c] = (uint32_t)q;
        g_sinv[c] = (uint32_t)inv;
        g_smagic[c] = (uint64_t)((((__uint128_t)1) << 64) / q); /* floor(2^64/q) */
        for (int j = 0; j < g_k; j++)
            g_D[c * g_k + j] =
                (uint32_t)(((q - g_off[j] % q) % q) * inv % q);
        c++;
    }
    free(all);

    g_tlo    = g_lo / g_M;
    g_tcount = g_hi / g_M - g_tlo + 1;

    /* Full-size blocks unless that would leave too few of them: the
     * verified frontier only advances a whole block at a time, so a run
     * with 3 blocks can report almost no progress.  Never go below 1024
     * t values, or per-block sieve setup starts to dominate. */
    /* Per (residue, block) we redo ~2 reductions per sieving prime, so a
     * block must be long enough that marking dominates that setup.  Aim
     * for ~64 blocks to give the verified frontier some granularity, but
     * never at the cost of more than a few percent throughput. */
    {
        uint64_t floor_t = 8 * (uint64_t)g_nsq;
        if (floor_t < 4096)          /* tiny blocks pay fixed per-block
                                      * overhead (memset, prime walk,
                                      * scan) far too many times */
            floor_t = 4096;
        g_tblock = TBLOCK;
        if (g_tcount / 64 < g_tblock) {
            uint64_t b = (g_tcount + 63) / 64;
            if (b < floor_t)   b = floor_t;
            if (b > g_tcount)  b = g_tcount;
            g_tblock = b;
        }
    }
    g_nblocks = (g_tcount + g_tblock - 1) / g_tblock;
    g_units_per_block = (g_nres + RCHUNK - 1) / RCHUNK;
    g_nunits = g_nblocks * g_units_per_block;

    g_tty = isatty(STDERR_FILENO);
    g_progress = !g_quiet;
    if (g_interval <= 0.0)
        g_interval = g_tty ? 5.0 : 60.0;

    double dens = 1.0;
    {
        size_t n2;
        uint32_t *wp = simple_primes(wheel_limit, &n2);
        for (size_t i = 0; i < n2; i++)
            dens *= 1.0 - (double)admissible_w(wp[i]) / wp[i];
        free(wp);
    }

    fprintf(stderr,
        "k=%d span=%" PRIu64 " interior=%d | wheel<=%u M=%" PRIu64
        " residues=%zu (density %.3e)\n"
        "sieve primes %zu (<=%u) | range %" PRIu64 "..%" PRIu64
        " | t values %" PRIu64 " | %" PRIu64 " work units | %ld thread%s\n",
        g_k, g_span, g_nint, wheel_limit, g_M, g_nres, dens,
        g_nsq, sieve_bound, g_lo, g_hi, g_tcount, g_nunits,
        nthreads, nthreads == 1 ? "" : "s");

    g_nthreads = nthreads;
    g_inflight = malloc((size_t)nthreads * sizeof *g_inflight);
    for (long i = 0; i < nthreads; i++)
        g_inflight[i] = UINT64_MAX;
    g_t0 = now_sec();
    g_last_print = g_t0;

    pthread_t *th = calloc((size_t)nthreads, sizeof *th);
    for (long i = 1; i < nthreads; i++)
        pthread_create(&th[i], NULL, worker, (void *)(intptr_t)i);
    worker((void *)(intptr_t)0);
    for (long i = 1; i < nthreads; i++)
        pthread_join(th[i], NULL);
    free(th);

    double el = now_sec() - g_t0;
    if (g_progress) {
        report_progress(1);
        fprintf(stderr, "\n");
    }
    fflush(stdout);

    char elb[32];
    fmt_dur(el, elb, sizeof elb);
    fprintf(stderr,
            "done in %s | %.3e numbers/s (%.3e candidates/s) | "
            "%" PRIu64 " candidate%s tested, %" PRIu64 " chain%s found\n",
            elb,
            el > 0 ? (double)(g_hi - g_lo) / el : 0.0,
            el > 0 ? (double)(g_hi - g_lo) * dens / el : 0.0,
            g_tested, g_tested == 1 ? "" : "s",
            g_found, g_found == 1 ? "" : "s");
    free(g_inflight);
    return 0;
}
