/*
 * a228768.c -- search for the next term of OEIS A228768
 *
 *   a(n) = smallest number that is an emirp in all bases 2..n
 *          (emirp in base b: p prime, reverse_b(p) prime, reverse_b(p) != p)
 *
 *   known: n:  2   3   4   5    6    7    8    9    10        11          12              13
 *          a: 11  11  53  61  193  193  193  193  93836531  1963209431  14322793967831  14322793967831
 *   wanted: a(14)  (expected ~1e16 per the OEIS heuristic)
 *
 * Strategy
 * --------
 *  1. Leading-digit range filter.  If d is p's leading digit in base b, then d is
 *     the LAST digit of reverse_b(p), so reverse_b(p) = d (mod q) for every prime
 *     q | b.  Hence gcd(d,b) must be 1 or the reversal is composite.  This is
 *     automatic for prime b, but for b in {4,6,8,9,10,12,14} it kills whole
 *     intervals at once, leaving a handful of wide admissible windows.  Run
 *     --windows to list them; they reproduce exactly the four candidate ranges
 *     given in the OEIS comment.
 *
 *  2. Segmented sieve over the surviving intervals to enumerate candidate primes p.
 *
 *  3. Per-candidate reversal tests, cheapest base first (2, then descending, since
 *     a bigger base means fewer digits to reverse), with early exit.  Each reversal
 *     gets trial division by primes < 103 before any modexp; ~90% of candidates die
 *     on the very first base.
 *
 * Build:  cc -O3 -march=native -pthread -o a228768 a228768.c -lm
 *         ...or ~1.7x faster with Kim Walisch's primesieve as the sieve backend:
 *         cc -O3 -march=native -pthread -DUSE_PRIMESIEVE a228768.c -lprimesieve -lm
 * Run:    ./a228768 --self-test
 *         ./a228768 --start 14322793967831 --limit 20000000000000000 --threads 16
 *
 * All arithmetic is 64-bit; reverse_14(p) < 14*p, so p must stay below ~1.3e18.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <signal.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>
#ifdef USE_PRIMESIEVE
#include <primesieve.h>
#endif

typedef uint64_t u64;
typedef unsigned __int128 u128;

/* ------------------------------------------------------------------ */
/* primality                                                          */
/* ------------------------------------------------------------------ */

static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }

/* --- Montgomery arithmetic (odd modulus).  The whole search lives or dies on
   the cost of a 64-bit modexp; hardware divq in the inner loop is ~80 cycles,
   Montgomery REDC is ~10. --------------------------------------------------- */

typedef struct { u64 n, ninv, r1, r2; } mont;

static inline u64 mont_ninv(u64 n)          /* -n^{-1} mod 2^64 */
{
    u64 x = n;                              /* correct mod 2^3 for odd n */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return (u64)0 - x;
}

static inline u64 redc(u128 t, u64 n, u64 ninv)
{
    u64 m = (u64)t * ninv;
    u64 r = (u64)((t + (u128)m * n) >> 64);
    return r >= n ? r - n : r;
}

static inline u64 mmul(u64 a, u64 b, const mont *M) { return redc((u128)a * b, M->n, M->ninv); }

static void mont_init(mont *M, u64 n)
{
    M->n = n;
    M->ninv = mont_ninv(n);
    M->r1 = (u64)((((u128)1 << 64) - n) % n);          /* 2^64 mod n */
    M->r2 = mulmod(M->r1, M->r1, n);                   /* 2^128 mod n, one divq */
}

static inline u64 to_mont(u64 a, const mont *M) { return mmul(a % M->n, M->r2, M); }

static int sprp_mont(const mont *M, u64 a, u64 d, int s)
{
    u64 n = M->n;
    u64 x = to_mont(a, M);
    if (x == 0) return 1;
    u64 r = M->r1, nm1 = n - M->r1;        /* 1 and -1 in Montgomery form */
    u64 y = r;
    for (u64 e = d; e; e >>= 1) {
        if (e & 1) y = mmul(y, x, M);
        x = mmul(x, x, M);
    }
    if (y == r || y == nm1) return 1;
    for (int i = 1; i < s; i++) { y = mmul(y, y, M); if (y == nm1) return 1; }
    return 0;
}

static const unsigned SMALLP[] = {
    3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101
};
#define NSMALLP (sizeof(SMALLP)/sizeof(SMALLP[0]))
#define SMALLP_SQ (101u*101u)

/* deterministic for all n < 2^64 (Jaeschke/Sinclair witness set), odd n > 101 */
static int mr_odd(u64 n)
{
    mont M; mont_init(&M, n);
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const u64 W[] = {2,325,9375,28178,450775,9780504,1795265022};
    for (unsigned i = 0; i < 7; i++) if (!sprp_mont(&M, W[i], d, s)) return 0;
    return 1;
}

static int is_prime_u64(u64 n)
{
    if (n < 2) return 0;
    if (!(n & 1)) return n == 2;
    for (unsigned i = 0; i < NSMALLP; i++) {
        unsigned p = SMALLP[i];
        if (n == p) return 1;
        if (n % p == 0) return 0;
        if ((u64)p * p > n) return 1;
    }
    return mr_odd(n);
}

/* hot path: short-circuiting trial division, then the deterministic test */
static inline int fast_is_prime(u64 q)
{
    if (q < SMALLP_SQ) return is_prime_u64(q);
    if (!(q & 1)) return 0;
    if (q %  3 == 0 || q %  5 == 0 || q %  7 == 0 || q % 11 == 0 ||
        q % 13 == 0 || q % 17 == 0 || q % 19 == 0 || q % 23 == 0 ||
        q % 29 == 0 || q % 31 == 0 || q % 37 == 0 || q % 41 == 0 ||
        q % 43 == 0 || q % 47 == 0 || q % 53 == 0 || q % 59 == 0 ||
        q % 61 == 0 || q % 67 == 0 || q % 71 == 0 || q % 73 == 0 ||
        q % 79 == 0 || q % 83 == 0 || q % 89 == 0 || q % 97 == 0 ||
        q % 101 == 0) return 0;
    return mr_odd(q);
}

/* ------------------------------------------------------------------ */
/* digit reversal                                                     */
/* ------------------------------------------------------------------ */

/* constant bases so the compiler turns the divisions into mul/shift */
#define DEFREV(B) \
    static inline u64 rev##B(u64 n) { u64 r = 0; while (n) { r = r*(B) + n%(B); n /= (B); } return r; }
DEFREV(3)  DEFREV(4)  DEFREV(5)  DEFREV(6)  DEFREV(7)
DEFREV(8)  DEFREV(9)  DEFREV(10) DEFREV(11) DEFREV(12)
DEFREV(13) DEFREV(14)

static inline u64 rev2(u64 n)
{
    u64 v = n;
    v = ((v >> 1) & 0x5555555555555555ull) | ((v & 0x5555555555555555ull) << 1);
    v = ((v >> 2) & 0x3333333333333333ull) | ((v & 0x3333333333333333ull) << 2);
    v = ((v >> 4) & 0x0F0F0F0F0F0F0F0Full) | ((v & 0x0F0F0F0F0F0F0F0Full) << 4);
    v = __builtin_bswap64(v);
    return v >> __builtin_clzll(n);          /* n != 0 guaranteed */
}

static inline u64 revb(u64 n, int b)
{
    switch (b) {
    case  2: return rev2(n);   case  3: return rev3(n);   case  4: return rev4(n);
    case  5: return rev5(n);   case  6: return rev6(n);   case  7: return rev7(n);
    case  8: return rev8(n);   case  9: return rev9(n);   case 10: return rev10(n);
    case 11: return rev11(n);  case 12: return rev12(n);  case 13: return rev13(n);
    case 14: return rev14(n);
    default: { u64 r = 0; while (n) { r = r*b + n%b; n /= b; } return r; }
    }
}

/* ------------------------------------------------------------------ */
/* per-thread statistics (each on its own cache lines, no contention)  */
/* ------------------------------------------------------------------ */

typedef struct __attribute__((aligned(128))) {
    u64 scanned;        /* admissible numbers sieved                    */
    u64 primes;         /* candidate primes handed to the emirp test    */
    u64 depth[16];      /* depth[k] = candidates that cleared exactly k bases */
    u64 best_depth;     /* deepest partial hit seen                     */
    u64 best_p;
    u64 pad[6];
} stats_t;

static stats_t  g_dummy;                    /* used by non-worker paths */
static stats_t *g_stats;
static __thread stats_t *tls = &g_dummy;

/* ------------------------------------------------------------------ */
/* config                                                             */
/* ------------------------------------------------------------------ */

static int NBASE = 14;                 /* test bases 2..NBASE            */
static int BORDER[16], NORDER;         /* test order, cheapest first     */

static void build_order(void)
{
    NORDER = 0;
    BORDER[NORDER++] = 2;              /* bit-reversal, nearly free      */
    for (int b = NBASE; b >= 3; b--)   /* big bases: fewer digits        */
        BORDER[NORDER++] = b;
}

static inline int is_emirp_all(u64 p)
{
    stats_t *st = tls;
    int i;
    for (i = 0; i < NORDER; i++) {
        u64 q = revb(p, BORDER[i]);
        if (q == p || !fast_is_prime(q)) break;
    }
    st->depth[i]++;                             /* cleared exactly i bases */
    if ((u64)i > st->best_depth) { st->best_depth = i; st->best_p = p; }
    return i == NORDER;
}

/* ------------------------------------------------------------------ */
/* leading-digit range filter                                         */
/* ------------------------------------------------------------------ */

static unsigned gcd_u(unsigned a, unsigned b) { while (b) { unsigned t = a % b; a = b; b = t; } return a; }

/* Smallest y >= x whose leading digit in every base 2..NBASE is coprime to that
   base.  *end receives the first point at or above y where some leading digit
   changes, so [y,*end) is uniformly admissible. */
static u64 next_permissible(u64 x, u64 *end)
{
restart:
    if (x < 2) x = 2;
    u64 hi = UINT64_MAX;
    for (int b = 2; b <= NBASE; b++) {
        if (x < (u64)b) { x = b; goto restart; }   /* single digit => palindrome */
        u64 pw = 1;
        while (pw <= x / (u64)b) pw *= (u64)b;     /* pw = b^(digits-1) */
        u64 d = x / pw;
        if (gcd_u((unsigned)d, (unsigned)b) != 1) { x = (d + 1) * pw; goto restart; }
        u64 h = (d + 1) * pw;
        if (h < hi) hi = h;
    }
    *end = hi;
    return x;
}

/* Maximal contiguous admissible windows in [start,limit).  Cheap to build --
   the walker only visits digit boundaries -- and it gives us an exact
   denominator for the progress bar instead of a raw-range percentage. */

typedef struct { u64 lo, hi; } win_t;
static win_t *g_win; static int g_nwin; static u64 g_total_adm;

static void build_windows(u64 start, u64 limit)
{
    int cap = 256;
    g_win = malloc((size_t)cap * sizeof(win_t)); g_nwin = 0; g_total_adm = 0;
    u64 x = start, run_lo = 0, run_hi = 0;
    for (;;) {
        u64 e = 0, s = 0, top = 0; int more = 0;
        if (x < limit) {
            s = next_permissible(x, &e);
            if (s < limit) { top = e < limit ? e : limit; more = 1; }
        }
        if (more && run_hi == s) { run_hi = top; g_total_adm += top - s; x = e; continue; }
        if (run_hi > run_lo) {
            if (g_nwin == cap) { cap *= 2; g_win = realloc(g_win, (size_t)cap * sizeof(win_t)); }
            g_win[g_nwin].lo = run_lo; g_win[g_nwin].hi = run_hi; g_nwin++;
        }
        if (!more) break;
        run_lo = s; run_hi = top; g_total_adm += top - s; x = e;
    }
}

/* which window contains (or next follows) x, 1-based */
static int win_index(u64 x)
{
    for (int i = 0; i < g_nwin; i++) if (x < g_win[i].hi) return i + 1;
    return g_nwin;
}

/* ------------------------------------------------------------------ */
/* base primes for the segmented sieve                                */
/* ------------------------------------------------------------------ */

#ifndef USE_PRIMESIEVE
static uint32_t *bprime;
static uint32_t  nbprime;

static void build_base_primes(u64 limit_sqrt)
{
    uint32_t n = (uint32_t)limit_sqrt + 1;
    uint32_t half = n / 2 + 1;                     /* odds only: i <-> 2i+1 */
    unsigned char *c = calloc(half, 1);
    if (!c) { fprintf(stderr, "oom base sieve\n"); exit(1); }
    for (uint32_t i = 1; (2*i+1)*(2*i+1) <= n; i++)
        if (!c[i]) {
            uint32_t p = 2*i + 1;
            for (uint64_t j = (uint64_t)p*p; j <= n; j += 2*p) c[j/2] = 1;
        }
    uint32_t cnt = 1;                              /* the 2 */
    for (uint32_t i = 1; i < half; i++) if (!c[i] && 2*i+1 <= n) cnt++;
    bprime = malloc((size_t)cnt * sizeof(uint32_t));
    nbprime = 0; bprime[nbprime++] = 2;
    for (uint32_t i = 1; i < half; i++) if (!c[i] && 2*i+1 <= n) bprime[nbprime++] = 2*i + 1;
    free(c);
}
#endif

/* ------------------------------------------------------------------ */
/* search                                                             */
/* ------------------------------------------------------------------ */

#define WINDOW (1u << 25)              /* numbers per sieve window */

static volatile sig_atomic_t stop_flag = 0;
static void on_sigint(int s) { (void)s; stop_flag = 1; }

static u64  g_start, g_limit, g_block;
static u64  g_next_block;              /* atomic cursor            */
static u64 *g_inflight;                /* per-thread block start   */
static int  g_nthreads;
static u64  g_found = 0;               /* smallest hit so far      */
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;

static void report_hit(u64 p)
{
    pthread_mutex_lock(&g_lock);
    if (!g_found || p < g_found) g_found = p;
    printf("\n*** HIT: %" PRIu64 " is an emirp in all bases 2..%d\n", p, NBASE);
    for (int b = 2; b <= NBASE; b++)
        printf("      base %2d: reverse = %" PRIu64 "\n", b, revb(p, b));
    fflush(stdout);
    pthread_mutex_unlock(&g_lock);
}

/* sieve [lo,hi) and emirp-test every prime in it */
#ifdef USE_PRIMESIEVE

/* primesieve does bucketed, pre-sieved, SIMD-assisted segmented sieving; at 1.8e16
   it enumerates ~2.1G numbers/s/core against ~165M/s for the built-in sieve. */
static void scan_interval(u64 lo, u64 hi, primesieve_iterator *it)
{
    if (hi <= lo) return;
    primesieve_jump_to(it, lo, hi);
    u64 p, base = lo, cnt = 0;
    while ((p = primesieve_next_prime(it)) < hi) {
        if (is_emirp_all(p)) report_hit(p);
        if (++cnt == (1u << 20)) {
            tls->primes += cnt; cnt = 0;
            tls->scanned += p - base; base = p;
            if (stop_flag) { return; }
        }
    }
    tls->primes += cnt;
    tls->scanned += hi - base;
}

#else

static void scan_interval(u64 lo, u64 hi, unsigned char *flags)
{
    if (hi <= lo) return;

    if (lo < 128) {                              /* tiny-number path */
        for (u64 n = lo; n < hi && n < 128; n++)
            if (is_prime_u64(n) && is_emirp_all(n)) report_hit(n);
        lo = (hi < 128) ? hi : 128;
        if (lo >= hi) return;
    }

    for (u64 w0 = lo; w0 < hi; w0 += WINDOW) {
        if (stop_flag) return;
        u64 w1 = w0 + WINDOW; if (w1 > hi) w1 = hi;
        u64 o0 = w0 | 1ull;                      /* first odd >= w0 */
        if (o0 >= w1) continue;
        size_t cnt = (size_t)((w1 - o0 + 1) / 2);
        memset(flags, 0, cnt);

        for (uint32_t i = 1; i < nbprime; i++) { /* skip 2 */
            u64 p = bprime[i];
            if (p * p >= w1) break;
            u64 s = (o0 + p - 1) / p * p;        /* first multiple >= o0 */
            if (!(s & 1)) s += p;                /* keep it odd          */
            if (s < p * p) s = p * p;
            if (!(s & 1)) s += p;
            for (u64 m = s; m < w1; m += 2*p) flags[(m - o0) >> 1] = 1;
        }
        u64 np = 0;
        for (size_t i = 0; i < cnt; i++)
            if (!flags[i]) {
                u64 p = o0 + 2*(u64)i;
                np++;
                if (is_emirp_all(p)) report_hit(p);
            }
        tls->primes  += np;
        tls->scanned += w1 - w0;
    }
}

#endif

static void *worker(void *arg)
{
    int id = (int)(intptr_t)arg;
    tls = &g_stats[id];
#ifdef USE_PRIMESIEVE
    primesieve_iterator itv; primesieve_init(&itv);
    primesieve_iterator *flags = &itv;
#else
    unsigned char *flags = malloc(WINDOW / 2 + 2);
    if (!flags) { fprintf(stderr, "oom window\n"); exit(1); }
#endif

    for (;;) {
        if (stop_flag) break;
        u64 blo = __atomic_fetch_add(&g_next_block, g_block, __ATOMIC_RELAXED);
        if (blo >= g_limit) break;
        u64 bhi = blo + g_block; if (bhi > g_limit) bhi = g_limit;
        __atomic_store_n(&g_inflight[id], blo, __ATOMIC_RELAXED);

        u64 x = blo;
        while (x < bhi) {
            u64 e, s = next_permissible(x, &e);
            if (s >= bhi) break;
            u64 top = (e < bhi) ? e : bhi;
            scan_interval(s, top, flags);
            x = e;
            if (stop_flag) break;
        }
        __atomic_store_n(&g_inflight[id], UINT64_MAX, __ATOMIC_RELAXED);
    }
#ifdef USE_PRIMESIEVE
    primesieve_free_iterator(&itv);
#else
    free(flags);
#endif
    return NULL;
}

/* ------------------------------------------------------------------ */
/* live status display                                                */
/* ------------------------------------------------------------------ */

static void fmt_dur(double sec, char *b, size_t n)
{
    if (!(sec >= 0.0) || sec > 3.0e10) { snprintf(b, n, "     --     "); return; }
    unsigned long long s = (unsigned long long)sec;
    unsigned long long d = s / 86400; s %= 86400;
    if (d) snprintf(b, n, "%llud %02lluh %02llum", d, s/3600, (s%3600)/60);
    else   snprintf(b, n, "%02lluh %02llum %02llus", s/3600, (s%3600)/60, s%60);
}

static void fmt_cnt(double v, char *b, size_t n)
{
    static const char suf[] = " KMGTPE";
    int i = 0;
    while (v >= 1000.0 && i < 6) { v /= 1000.0; i++; }
    snprintf(b, n, "%6.2f%c", v, suf[i]);
}

typedef struct { u64 scanned, primes, depth[16], best_depth, best_p; } agg_t;

static void aggregate(agg_t *a)
{
    memset(a, 0, sizeof(*a));
    for (int t = 0; t < g_nthreads; t++) {
        stats_t *st = &g_stats[t];
        a->scanned += __atomic_load_n(&st->scanned, __ATOMIC_RELAXED);
        a->primes  += __atomic_load_n(&st->primes,  __ATOMIC_RELAXED);
        for (int k = 0; k < 16; k++)
            a->depth[k] += __atomic_load_n(&st->depth[k], __ATOMIC_RELAXED);
        u64 bd = __atomic_load_n(&st->best_depth, __ATOMIC_RELAXED);
        if (bd > a->best_depth) { a->best_depth = bd; a->best_p = __atomic_load_n(&st->best_p, __ATOMIC_RELAXED); }
    }
}

static int g_tty = 0;
static int g_status_lines = 0;          /* lines to rewind; 0 = print fresh */

static u64 watermark(void)
{
    u64 w = __atomic_load_n(&g_next_block, __ATOMIC_RELAXED);
    if (w > g_limit) w = g_limit;
    for (int i = 0; i < g_nthreads; i++) {
        u64 v = __atomic_load_n(&g_inflight[i], __ATOMIC_RELAXED);
        if (v < w) w = v;
    }
    return w;
}

/* ------------------------------------------------------------------ */
/* self test                                                          */
/* ------------------------------------------------------------------ */

static u64 brute_term(int n)
{
    int save = NBASE; NBASE = n; build_order();
    u64 x = 2, res = 0;
    while (!res) {
        u64 e, s = next_permissible(x, &e);
        for (u64 p = s; p < e; p++)
            if (is_prime_u64(p) && is_emirp_all(p)) { res = p; break; }
        x = e;
    }
    NBASE = save; build_order();
    return res;
}

static int self_test(void)
{
    static const u64 known[] = {
        11, 11, 53, 61, 193, 193, 193, 193,
        93836531, 1963209431, 14322793967831, 14322793967831
    };
    int fail = 0;

    puts("checking that the published terms satisfy the definition:");
    for (int n = 2; n <= 13; n++) {
        u64 v = known[n-2];
        int save = NBASE; NBASE = n; build_order();
        int ok = is_prime_u64(v) && is_emirp_all(v);
        /* the term must also survive the leading-digit filter, or the filter is wrong */
        u64 e2, s2 = next_permissible(v, &e2);
        int inrange = (s2 <= v && v < e2);
        NBASE = save; build_order();
        printf("  a(%2d) = %-16" PRIu64 " emirp:%s  filter-admissible:%s\n",
               n, v, ok ? "yes" : "NO", inrange ? "yes" : "NO");
        if (!ok || !inrange) fail = 1;
    }

    puts("\nrecomputing small terms from scratch:");
    for (int n = 2; n <= 11; n++) {
        clock_t t0 = clock();
        u64 v = brute_term(n);
        double dt = (double)(clock() - t0) / CLOCKS_PER_SEC;
        int ok = (v == known[n-2]);
        printf("  a(%2d) = %-16" PRIu64 " %s  (%.2fs)\n", n, v, ok ? "ok" : "MISMATCH", dt);
        if (!ok) fail = 1;
    }
    return fail;
}

/* ------------------------------------------------------------------ */

static void usage(const char *me)
{
    fprintf(stderr,
      "usage: %s [options]\n"
      "  --start N      begin search at N            (default 14322793967831)\n"
      "  --limit N      stop below N                 (default 2e16)\n"
      "  --bases N      require emirp in bases 2..N  (default 14)\n"
      "  --threads N    worker threads               (default 4)\n"
      "  --block N      work-unit size               (default 2^30)\n"
      "  --checkpoint F write resume watermark to F\n"
      "  --progress S   status refresh interval in seconds (default 2)\n"
      "  --self-test    verify against the known terms and exit\n"
      "  --windows      list maximal admissible windows and exit\n"
      "  --win-min N    only list windows at least N wide (default 1e9)\n"
      "  --verify N     show the per-base reversals of N and exit\n", me);
    exit(1);
}

int main(int argc, char **argv)
{
    u64 start = 14322793967831ull, limit = 20000000000000000ull, block = 1ull << 30;
    int threads = 4;
    const char *ckpt = NULL;
    double prog_sec = 2.0;
    int mode_test = 0, mode_windows = 0; u64 verify = 0, win_min = 1000000000ull;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--start")   && i+1 < argc) start   = strtoull(argv[++i], 0, 0);
        else if (!strcmp(argv[i], "--limit")   && i+1 < argc) limit   = strtoull(argv[++i], 0, 0);
        else if (!strcmp(argv[i], "--bases")   && i+1 < argc) NBASE   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--threads") && i+1 < argc) threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--block")   && i+1 < argc) block   = strtoull(argv[++i], 0, 0);
        else if (!strcmp(argv[i], "--checkpoint") && i+1 < argc) ckpt = argv[++i];
        else if (!strcmp(argv[i], "--progress") && i+1 < argc) prog_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "--self-test")) mode_test = 1;
        else if (!strcmp(argv[i], "--windows")) mode_windows = 1;
        else if (!strcmp(argv[i], "--win-min") && i+1 < argc) win_min = strtoull(argv[++i], 0, 0);
        else if (!strcmp(argv[i], "--verify")  && i+1 < argc) verify  = strtoull(argv[++i], 0, 0);
        else usage(argv[0]);
    }
    if (NBASE < 2 || NBASE > 14) { fprintf(stderr, "bases must be 2..14\n"); return 1; }
    if (limit > 2000000000000000000ull) { fprintf(stderr, "limit too large for 64-bit reversals\n"); return 1; }
    build_order();

    if (verify) {
        printf("%" PRIu64 " prime: %s\n", verify, is_prime_u64(verify) ? "yes" : "no");
        for (int b = 2; b <= 14; b++) {
            u64 q = revb(verify, b);
            printf("  base %2d  reverse = %-20" PRIu64 " %s\n", b, q,
                   q == verify ? "palindrome" : (is_prime_u64(q) ? "prime" : "composite"));
        }
        return 0;
    }
    if (mode_test) return self_test();

    if (mode_windows) {
        build_windows(start, limit);
        printf("maximal admissible windows in [%" PRIu64 ", %" PRIu64 ") for bases 2..%d\n",
               start, limit, NBASE);
        int n = 0;
        for (int i = 0; i < g_nwin; i++) {
            u64 len = g_win[i].hi - g_win[i].lo;
            if (len < win_min) continue;
            printf(" %3d  [%20" PRIu64 ", %20" PRIu64 ")  %18" PRIu64 "  %9.1f core-hours\n",
                   ++n, g_win[i].lo, g_win[i].hi, len, len / 74.0e6 / 3600.0);
        }
        printf("total admissible: %" PRIu64 "  (%.4f%% of range, ~%.1f core-days at 74M/s)\n",
               g_total_adm, 100.0*(double)g_total_adm/(double)(limit-start),
               g_total_adm/74.0e6/86400.0);
        return 0;
    }

    build_windows(start, limit);
    printf("leading-digit filter keeps %.4f%% of the range: %" PRIu64
           " numbers to sieve in %d window%s\n",
           100.0 * (double)g_total_adm / (double)(limit - start),
           g_total_adm, g_nwin, g_nwin == 1 ? "" : "s");

#ifdef USE_PRIMESIEVE
    printf("sieve backend: primesieve\n");
#else
    build_base_primes((u64)(1.0 + __builtin_sqrt((double)limit)));
    printf("sieve backend: built-in segmented (%u base primes)\n", nbprime);
#endif

    g_start = start; g_limit = limit; g_block = block;
    g_next_block = start; g_nthreads = threads;
    g_inflight = malloc(sizeof(u64) * threads);
    g_stats    = aligned_alloc(128, sizeof(stats_t) * (size_t)threads);
    memset(g_stats, 0, sizeof(stats_t) * (size_t)threads);
    for (int i = 0; i < threads; i++) g_inflight[i] = UINT64_MAX;

    signal(SIGINT, on_sigint);
    printf("searching [%" PRIu64 ", %" PRIu64 ") for emirps in bases 2..%d on %d threads\n",
           start, limit, NBASE, threads);

    pthread_t *th = malloc(sizeof(pthread_t) * threads);
    struct timespec t0; clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < threads; i++) pthread_create(&th[i], 0, worker, (void *)(intptr_t)i);

    /* ---- live status ---- */
    g_tty = isatty(1);
    if (!g_tty && prog_sec < 30.0) prog_sec = 30.0;   /* logs: don't spam */
    double last_t = 0, ema = 0;
    u64 last_scanned = 0;
    agg_t a;

    for (;;) {
        struct timespec ts;
        ts.tv_sec  = (time_t)prog_sec;
        ts.tv_nsec = (long)((prog_sec - (double)ts.tv_sec) * 1e9);
        nanosleep(&ts, 0);

        u64 w = watermark();
        struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
        double now = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
        aggregate(&a);

        double dt   = now - last_t;
        double inst = dt > 0 ? (double)(a.scanned - last_scanned) / dt : 0;
        ema  = (ema == 0) ? inst : 0.7 * ema + 0.3 * inst;
        double avg  = now > 0 ? (double)a.scanned / now : 0;
        last_t = now; last_scanned = a.scanned;

        double pct  = g_total_adm ? 100.0 * (double)a.scanned / (double)g_total_adm : 0;
        double rem  = (double)g_total_adm - (double)a.scanned; if (rem < 0) rem = 0;
        double eta  = ema > 1.0 ? rem / ema : -1;

        char cel[32], ceta[32], cdone[16], ctot[16], cpri[16], crate[16], cavg[16], cprate[16];
        fmt_dur(now, cel, sizeof cel);
        fmt_dur(eta, ceta, sizeof ceta);
        fmt_cnt((double)a.scanned,   cdone,  sizeof cdone);
        fmt_cnt((double)g_total_adm, ctot,   sizeof ctot);
        fmt_cnt((double)a.primes,    cpri,   sizeof cpri);
        fmt_cnt(ema,                 crate,  sizeof crate);
        fmt_cnt(avg,                 cavg,   sizeof cavg);
        fmt_cnt(now > 0 ? (double)a.primes / now : 0, cprate, sizeof cprate);

        /* survivor funnel: how many candidates cleared >= k bases */
        u64 surv[17]; surv[NORDER] = a.depth[NORDER];
        for (int k = NORDER - 1; k >= 0; k--) surv[k] = surv[k+1] + a.depth[k];

        /* projected wall-clock finish */
        char cfin[64] = "";
        if (eta > 0) {
            time_t fin = time(0) + (time_t)eta;
            struct tm tmv; localtime_r(&fin, &tmv);
            strftime(cfin, sizeof cfin, " (%b %d %H:%M)", &tmv);
        }

        const char *EK = g_tty ? "\033[K" : "";
        pthread_mutex_lock(&g_lock);
        if (g_tty && g_status_lines) printf("\033[%dA", g_status_lines);
        int lines = 0;
        printf("%s  elapsed %s   ETA %s%s   window %d/%d   resume at %" PRIu64 "\n",
               EK, cel, ceta, cfin, win_index(w), g_nwin, w); lines++;
        char bar[33];
        { int fill = (int)(pct * 32.0 / 100.0); if (fill > 32) fill = 32; if (fill < 0) fill = 0;
          memset(bar, '.', 32); memset(bar, '#', (size_t)fill); bar[32] = 0; }
        printf("%s  scanned %s / %s admissible  [%s] %6.3f%%\n",
               EK, cdone, ctot, bar, pct); lines++;
        printf("%s  rate    %s/s now, %s/s avg   primes %s tested (%s/s)\n",
               EK, crate, cavg, cpri, cprate); lines++;
        printf("%s  funnel  ", EK);
        for (int k = 1; k <= NORDER && k <= 6; k++)
            printf("%d:%.3g%% ", BORDER[k-1],
                   surv[0] ? 100.0 * (double)surv[k] / (double)surv[0] : 0.0);
        printf("  deepest %llu/%d bases at %" PRIu64 "\n",
               (unsigned long long)a.best_depth, NORDER, a.best_p); lines++;
        if (g_found) { printf("%s  HIT     %" PRIu64 "\n", EK, g_found); lines++; }
        g_status_lines = g_tty ? lines : 0;
        fflush(stdout);
        pthread_mutex_unlock(&g_lock);

        if (ckpt) {
            FILE *f = fopen(ckpt, "w");
            if (f) { fprintf(f, "%" PRIu64 "\n", w); fclose(f); }
        }
        if (w >= limit || stop_flag) break;
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], 0);

    u64 w = watermark();
    struct timespec t1; clock_gettime(CLOCK_MONOTONIC, &t1);
    double now = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);
    agg_t fa; aggregate(&fa);
    char cel[32], cdone[16], cpri[16];
    fmt_dur(now, cel, sizeof cel);
    fmt_cnt((double)fa.scanned, cdone, sizeof cdone);
    fmt_cnt((double)fa.primes,  cpri,  sizeof cpri);
    printf("\n%s\n", stop_flag ? "interrupted" : "range complete");
    printf("  elapsed        %s\n", cel);
    printf("  complete thru  %" PRIu64 "%s\n", w, ckpt ? " (written to checkpoint)" : "");
    printf("  scanned        %s admissible numbers, %s primes tested\n", cdone, cpri);
    printf("  deepest        %llu of %d bases, at %" PRIu64 "\n",
           (unsigned long long)fa.best_depth, NORDER, fa.best_p);
    if (g_found) printf("  SMALLEST HIT   %" PRIu64 "\n", g_found);
    else         printf("  no term found in the scanned range\n");
    return 0;
}
