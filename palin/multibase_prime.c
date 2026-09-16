/*
 * multibase_prime.c
 *
 * Explore the sequence:
 *   a(n) = smallest n-digit number whose decimal digit string is prime
 *          when interpreted in every base 10 through 10+n-1.
 *
 * Known values: 2, 29, 797, 2221, 35771, 175781, 1698311
 *
 * Usage:
 *   multibase_prime [-t threads] [-o] [n]
 *
 *   n        sequence index to search for (1..19). If omitted, starts at
 *            n = 1 and keeps going.
 *   -t N     number of worker threads (default: number of online CPUs).
 *   -o       also append each found value to output.txt.
 *
 * Results go to stdout ("a(n) = value"); a live status line (percent, rate,
 * ETA, elapsed, best-so-far) is repainted on stderr with a carriage return.
 *
 * Method:
 *   For each n, scan the n-digit decimal range [10^(n-1), 10^n - 1] in
 *   fixed-size chunks handed to worker threads via an atomic counter.
 *   Each chunk is sieved (segmented sieve of Eratosthenes) to find numbers
 *   that are prime in base 10; each surviving prime's digit string is then
 *   evaluated in bases 11..10+n-1 and tested with deterministic
 *   Miller-Rabin (first 12 prime witnesses: deterministic below ~3.3e24,
 *   which covers every value this search can practically reach).
 *   The smallest hit is kept with an atomic compare-and-swap min; threads
 *   stop taking chunks that start above the current best, so the final
 *   value after joining all threads is exactly the minimum.
 *
 * The decimal range is limited to n <= 19 so candidates fit in uint64_t;
 * values in higher bases are handled in unsigned __int128.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef unsigned __int128 u128;

#define MAX_N       19               /* 10^19 - 1 still fits in uint64_t   */
#define CHUNK       ((u64)1 << 20)   /* candidates per work unit           */
#define NO_BEST     UINT64_MAX
#define SIEVE_CAP   ((u64)10000000)  /* cap on base-prime sieve limit      */

/* ------------------------------------------------------------------ */
/* Search state for the current n                                      */
/* ------------------------------------------------------------------ */

static int  g_nd;                    /* current digit count n              */
static u64  g_lo, g_hi;              /* decimal range being scanned        */
static bool g_sieve_complete;        /* sieve alone proves base-10 primality */

static u32   *g_base_primes;
static size_t g_nbase;

static _Atomic u64  g_next_chunk;
static _Atomic u64  g_best;
static _Atomic u64  g_done;          /* candidates covered (for progress)  */
static _Atomic bool g_stop_monitor;

static double g_t0;
static FILE  *g_out;                 /* -o: append results here, or NULL   */

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static u64 pow10_u64(int e)
{
    u64 v = 1;
    while (e-- > 0)
        v *= 10;
    return v;
}

/* Format a count as a short human string, e.g. "3.42M". */
static void fmt_count(double v, char *buf, size_t sz)
{
    if (v >= 1e12)      snprintf(buf, sz, "%.2fT", v / 1e12);
    else if (v >= 1e9)  snprintf(buf, sz, "%.2fG", v / 1e9);
    else if (v >= 1e6)  snprintf(buf, sz, "%.2fM", v / 1e6);
    else if (v >= 1e3)  snprintf(buf, sz, "%.1fK", v / 1e3);
    else                snprintf(buf, sz, "%.0f", v);
}

/* Format seconds as h:mm:ss. */
static void fmt_time(double s, char *buf, size_t sz)
{
    if (s < 0 || s > 99.0 * 3600.0 * 100.0) {
        snprintf(buf, sz, "--:--:--");
        return;
    }
    long t = (long)(s + 0.5);
    snprintf(buf, sz, "%ld:%02ld:%02ld", t / 3600, (t / 60) % 60, t % 60);
}

/* ------------------------------------------------------------------ */
/* Primality: deterministic Miller-Rabin on unsigned __int128          */
/* ------------------------------------------------------------------ */

static inline u128 mulmod(u128 a, u128 b, u128 m)
{
    if (m <= (u128)UINT64_MAX) {
        /* a, b < m <= 2^64-1, so the product fits in 128 bits */
        return ((u128)(u64)a * (u64)b) % m;
    }
    /* Rare slow path (values >= 2^64): shift-and-add. */
    u128 r = 0;
    a %= m;
    while (b) {
        if (b & 1) {
            r += a;
            if (r >= m)
                r -= m;
        }
        a <<= 1;
        if (a >= m)
            a -= m;
        b >>= 1;
    }
    return r;
}

static u128 powmod(u128 a, u128 e, u128 m)
{
    u128 r = 1;
    a %= m;
    while (e) {
        if (e & 1)
            r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

/*
 * Deterministic for all inputs < 3,317,044,064,679,887,385,961,981
 * (first 12 prime witnesses) — far beyond anything this search reaches.
 */
static bool is_prime_u128(u128 n)
{
    static const u64 wit[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    size_t i;

    if (n < 2)
        return false;
    for (i = 0; i < sizeof wit / sizeof wit[0]; i++) {
        if (n % wit[i] == 0)
            return n == wit[i];
    }
    if (n < 41 * 41)
        return true;

    u128 d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }

    for (i = 0; i < sizeof wit / sizeof wit[0]; i++) {
        u128 x = powmod(wit[i], d, n);
        if (x == 1 || x == n - 1)
            continue;
        int j;
        for (j = 0; j < r - 1; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1)
                break;
        }
        if (j == r - 1)
            return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Candidate check: digit string prime in bases 11 .. 10+n-1           */
/* (base 10 is established by the sieve / MR before this is called)    */
/* ------------------------------------------------------------------ */

static bool check_candidate(u64 p)
{
    int dig[MAX_N + 1];
    int i, b;
    u64 t = p;

    for (i = 0; i < g_nd; i++) {
        dig[i] = (int)(t % 10);
        t /= 10;
    }

    for (b = 11; b < 10 + g_nd; b++) {
        u128 v = 0;
        for (i = g_nd - 1; i >= 0; i--)
            v = v * (unsigned)b + (unsigned)dig[i];
        if (!is_prime_u128(v))
            return false;
    }
    return true;
}

static void update_best(u64 v)
{
    u64 cur = atomic_load(&g_best);
    while (v < cur &&
           !atomic_compare_exchange_weak(&g_best, &cur, v)) {
        /* cur reloaded by CAS on failure */
    }
}

/* ------------------------------------------------------------------ */
/* Base primes for the segmented sieve                                 */
/* ------------------------------------------------------------------ */

static void build_base_primes(u64 limit)
{
    size_t L = (size_t)limit + 1;
    unsigned char *comp = calloc(L, 1);
    size_t i, j, cnt = 0;

    if (!comp) {
        fprintf(stderr, "out of memory building base primes\n");
        exit(1);
    }
    for (i = 2; i * i < L; i++)
        if (!comp[i])
            for (j = i * i; j < L; j += i)
                comp[j] = 1;
    for (i = 2; i < L; i++)
        if (!comp[i])
            cnt++;

    free(g_base_primes);
    g_base_primes = malloc(cnt * sizeof(u32));
    if (!g_base_primes) {
        fprintf(stderr, "out of memory building base primes\n");
        exit(1);
    }
    g_nbase = 0;
    for (i = 2; i < L; i++)
        if (!comp[i])
            g_base_primes[g_nbase++] = (u32)i;
    free(comp);
}

/* ------------------------------------------------------------------ */
/* Worker threads                                                      */
/* ------------------------------------------------------------------ */

static void *worker(void *arg)
{
    (void)arg;
    unsigned char *comp = malloc(CHUNK / 2 + 2);
    if (!comp)
        return NULL;

    for (;;) {
        u64 c = atomic_fetch_add(&g_next_chunk, 1);
        u128 s128 = (u128)g_lo + (u128)c * CHUNK;
        if (s128 > g_hi)
            break;
        u64 s = (u64)s128;
        if (s > atomic_load(&g_best))
            break;                    /* everything here is above the best */
        u64 e = (g_hi - s >= CHUNK - 1) ? s + CHUNK - 1 : g_hi;

        /* 2 is the only even prime; only reachable when n == 1 */
        if (s <= 2 && 2 <= e && check_candidate(2))
            update_best(2);

        /* segmented sieve of the odd numbers in [s, e] */
        u64 os = s | 1;
        size_t half = (size_t)((e - os) / 2 + 1);
        memset(comp, 0, half);
        for (size_t i = 1; i < g_nbase; i++) {   /* skip prime 2 */
            u64 q = g_base_primes[i];
            u64 qq = q * q;
            if (qq > e)
                break;
            u64 start = qq;
            if (start < s) {
                start = (s + q - 1) / q * q;
                if ((start & 1) == 0)
                    start += q;
            }
            for (u64 m = start; m <= e; m += q << 1)
                comp[(m - os) >> 1] = 1;
        }

        u64 bsnap = atomic_load(&g_best);
        unsigned tick = 0;
        for (u64 p = (os < 3) ? 3 : os; p <= e; p += 2) {
            if (p > bsnap)
                break;
            if (comp[(p - os) >> 1])
                continue;
            /* if the sieve was partial, confirm base-10 primality */
            if (!g_sieve_complete && !is_prime_u128(p))
                continue;
            if (check_candidate(p)) {
                update_best(p);
                break;                /* nothing later in this chunk matters */
            }
            if ((++tick & 1023) == 0)
                bsnap = atomic_load(&g_best);
        }

        atomic_fetch_add(&g_done, e - s + 1);
    }

    free(comp);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Status line                                                         */
/* ------------------------------------------------------------------ */

static void *monitor(void *arg)
{
    (void)arg;
    const struct timespec nap = { 0, 250 * 1000 * 1000 };

    while (!atomic_load(&g_stop_monitor)) {
        nanosleep(&nap, NULL);

        double el    = now_sec() - g_t0;
        u64    done  = atomic_load(&g_done);
        double total = (double)(g_hi - g_lo) + 1.0;
        double pct   = 100.0 * (double)done / total;
        if (pct > 100.0)
            pct = 100.0;
        double rate = (el > 0.001) ? (double)done / el : 0.0;
        double eta  = (rate > 0) ? (total - (double)done) / rate : -1.0;

        char db[24], tb[24], rb[24], eb[24], lb[24], bb[32];
        fmt_count((double)done, db, sizeof db);
        fmt_count(total, tb, sizeof tb);
        fmt_count(rate, rb, sizeof rb);
        fmt_time(eta, eb, sizeof eb);
        fmt_time(el, lb, sizeof lb);
        u64 best = atomic_load(&g_best);
        if (best == NO_BEST)
            snprintf(bb, sizeof bb, "-");
        else
            snprintf(bb, sizeof bb, "%" PRIu64, best);

        fprintf(stderr,
                "\ra(%d): %5.1f%%  %s/%s  %s/s  ETA %s  elapsed %s  best %s        ",
                g_nd, pct, db, tb, rb, eb, lb, bb);
        fflush(stderr);
    }
    return NULL;
}

static void clear_status_line(void)
{
    fprintf(stderr, "\r%*s\r", 100, "");
    fflush(stderr);
}

/* ------------------------------------------------------------------ */
/* Search one sequence index                                           */
/* ------------------------------------------------------------------ */

static void run_search(int n, int nthreads)
{
    g_nd = n;
    g_lo = pow10_u64(n - 1);
    g_hi = pow10_u64(n) - 1;         /* n == 19: computed exactly as
                                        10*pow10(18)-1 without overflow */
    if (n == 19)
        g_hi = 9999999999999999999ULL;

    /* isqrt of g_hi, capped so the base-prime table stays small */
    u64 lim = 1;
    while ((u128)(lim + 1) * (lim + 1) <= g_hi)
        lim++;
    g_sieve_complete = (lim <= SIEVE_CAP);
    if (lim > SIEVE_CAP)
        lim = SIEVE_CAP;
    if (lim < 3)
        lim = 3;
    build_base_primes(lim);

    atomic_store(&g_next_chunk, 0);
    atomic_store(&g_best, NO_BEST);
    atomic_store(&g_done, 0);
    atomic_store(&g_stop_monitor, false);
    g_t0 = now_sec();

    pthread_t mon;
    pthread_t *tid = malloc((size_t)nthreads * sizeof(pthread_t));
    if (!tid) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    pthread_create(&mon, NULL, monitor, NULL);
    for (int i = 0; i < nthreads; i++)
        pthread_create(&tid[i], NULL, worker, NULL);
    for (int i = 0; i < nthreads; i++)
        pthread_join(tid[i], NULL);
    atomic_store(&g_stop_monitor, true);
    pthread_join(mon, NULL);
    free(tid);

    double el  = now_sec() - g_t0;
    u64 best   = atomic_load(&g_best);
    u64 done   = atomic_load(&g_done);
    char db[24], rb[24], lb[24];
    fmt_count((double)done, db, sizeof db);
    fmt_count(el > 0.001 ? (double)done / el : 0.0, rb, sizeof rb);
    fmt_time(el, lb, sizeof lb);

    clear_status_line();
    if (best != NO_BEST) {
        printf("a(%d) = %" PRIu64 "\n", n, best);
        fflush(stdout);
        if (g_out) {
            fprintf(g_out, "a(%d) = %" PRIu64 "\n", n, best);
            fflush(g_out);
        }
        fprintf(stderr, "  [n=%d: scanned %s in %s (%s/s)]\n", n, db, lb, rb);
    } else {
        printf("a(%d): no solution in the %d-digit range\n", n, n);
        fflush(stdout);
        if (g_out) {
            fprintf(g_out, "a(%d): no solution in the %d-digit range\n", n, n);
            fflush(g_out);
        }
    }
}

/* ------------------------------------------------------------------ */

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-t threads] [-o] [n]\n"
        "  a(n) = smallest n-digit number whose decimal digit string is\n"
        "         prime in every base 10 through 10+n-1.\n"
        "  n        sequence index to search (1..%d); omit to start at 1\n"
        "           and keep going\n"
        "  -t N     worker threads (default: online CPU count)\n"
        "  -o       also append each found value to output.txt\n",
        prog, MAX_N);
}

int main(int argc, char **argv)
{
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    int target = 0;
    int opt;

    if (nthreads < 1)
        nthreads = 1;

    while ((opt = getopt(argc, argv, "t:oh")) != -1) {
        switch (opt) {
        case 't':
            nthreads = atoi(optarg);
            if (nthreads < 1 || nthreads > 1024) {
                fprintf(stderr, "invalid thread count '%s'\n", optarg);
                return 1;
            }
            break;
        case 'o':
            g_out = fopen("output.txt", "a");
            if (!g_out) {
                perror("output.txt");
                return 1;
            }
            break;
        case 'h':
            usage(argv[0]);
            return 0;
        default:
            usage(argv[0]);
            return 1;
        }
    }
    if (optind < argc) {
        target = atoi(argv[optind]);
        if (target < 1 || target > MAX_N) {
            fprintf(stderr, "n must be 1..%d\n", MAX_N);
            return 1;
        }
        if (optind + 1 < argc) {
            usage(argv[0]);
            return 1;
        }
    }

    if (target) {
        fprintf(stderr, "searching a(%d) with %d thread(s)\n",
                target, nthreads);
        run_search(target, nthreads);
    } else {
        fprintf(stderr,
                "searching a(1), a(2), ... with %d thread(s); Ctrl-C to stop\n",
                nthreads);
        for (int n = 1; n <= MAX_N; n++)
            run_search(n, nthreads);
        fprintf(stderr,
                "stopping: n > %d exceeds the 64-bit decimal range\n", MAX_N);
    }

    if (g_out)
        fclose(g_out);
    free(g_base_primes);
    return 0;
}
