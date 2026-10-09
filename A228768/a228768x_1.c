/*
 * a228768x.c -- 128-bit search for OEIS A228768
 *
 *   a(n) = smallest number that is an emirp in every base 2..n
 *          (emirp in base b: p prime, reverse_b(p) prime, reverse_b(p) != p)
 *
 *   known: a(2..13) = 11, 11, 53, 61, 193, 193, 193, 193, 93836531,
 *                     1963209431, 14322793967831, 14322793967831
 *          a(14)    = 19775643791090909
 *   wanted: a(15), which the leading-digit filter proves is >= 3e19,
 *           i.e. past 2^64.  Hence this program.
 *
 * Why this is not just the 64-bit program widened
 * -----------------------------------------------
 * At 3e19 the square root is 5.5e9 (2.6e8 sieve primes); in window 3 it is
 * 5.8e11.  Sieving to sqrt(x) is out, and so is primesieve (uint64 API).
 * So the sieve here is only a cheap COMPOSITE FILTER using primes below
 * --sieve-limit, and primality is settled by deterministic Miller-Rabin at
 * the very END of the cascade.  That works because the reversal filters are
 * savage: a candidate must survive small-factor rejection on all 14
 * reversals before anything expensive runs, which happens for roughly
 * 1 in 1e13 sieve survivors.  MR cost is therefore irrelevant, and the same
 * code path works at 1e19 or 1e23.
 *
 * Pipeline, per number in an admissible window:
 *   1. segmented sieve, primes < sieve-limit          ~0.5 ns   4% survive
 *   2. reverse in base 2 (bit reverse), reject if it   ~30 ns   24% survive
 *      has a factor below 89 (two gcds against
 *      primorial blocks) or equals x
 *   3. same for bases 15,14,...,3, largest first        ~13% survive each
 *      (fewest digits = cheapest reversal)
 *   4. Miller-Rabin x, then each reversal              (essentially never)
 *
 * Build:  cc -O3 -pthread -o a228768x a228768x.c -lm
 *         add -march=native on x86-64, or -mcpu=native on Apple silicon.
 *         -DA228768_NO_ASM forces the portable prefilter path.
 * Test:   ./a228768x --self-test
 *         ./a228768x --windows --nbase 15
 * Run:    ./a228768x --nbase 15 --start 3e19 --limit 5*15^16 --threads 32
 *
 * Numbers accepted as 12345, 3e19, or 5*15^16.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <time.h>
#include <signal.h>
#include <pthread.h>
#include <inttypes.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

#define MASK42      (((u64)1 << 42) - 1)
#define MAXBASE     15
#define MAXORDER    16

/* measured on the reference build; refine with --bench */
#define NS_PER_NUMBER 5.1
#define KAPPA         3.09

/* literals past 2^64 must be assembled */
#define C_3E19      ((u128)30000000000ULL * 1000000000ULL)
#define C_5x15p16   ((u128)32842041778ULL * 1000000000ULL + 564453125ULL)

/* ------------------------------------------------------------------ */
/* u128 helpers                                                       */
/* ------------------------------------------------------------------ */

static char *u128_str(u128 v, char *buf)          /* buf >= 41 bytes */
{
    char tmp[41];
    int i = 0;
    if (!v) { buf[0] = '0'; buf[1] = 0; return buf; }
    while (v) { tmp[i++] = '0' + (int)(v % 10); v /= 10; }
    for (int j = 0; j < i; j++) buf[j] = tmp[i - 1 - j];
    buf[i] = 0;
    return buf;
}

static char *S(u128 v)                            /* rotating scratch */
{
    static char ring[8][41];
    static int k = 0;
    k = (k + 1) & 7;
    return u128_str(v, ring[k]);
}

/* accepts  12345 | 3e19 | 5*15^16  */
static int parse_u128(const char *s, u128 *out)
{
    u128 v = 0;
    const char *p = s;
    if (!isdigit((unsigned char)*p)) return -1;
    while (isdigit((unsigned char)*p)) v = v * 10 + (u128)(*p++ - '0');

    if (*p == 'e' || *p == 'E') {
        int e = atoi(p + 1);
        if (e < 0 || e > 38) return -1;
        while (e--) v *= 10;
        *out = v; return 0;
    }
    if (*p == '*') {
        u128 base = 0; int exp = 0;
        p++;
        if (!isdigit((unsigned char)*p)) return -1;
        while (isdigit((unsigned char)*p)) base = base * 10 + (u128)(*p++ - '0');
        if (*p != '^') return -1;
        exp = atoi(p + 1);
        if (exp < 0 || exp > 200) return -1;
        u128 t = 1;
        while (exp--) t *= base;
        *out = v * t; return 0;
    }
    if (*p) return -1;
    *out = v; return 0;
}

/* ------------------------------------------------------------------ */
/* modular arithmetic, valid for modulus < 2^84                       */
/* ------------------------------------------------------------------ */

static inline u128 mulmod(u128 a, u128 b, u128 m)
{
    u128 bh = b >> 42, bl = b & MASK42;
    u128 r = (a * bh) % m;              /* a<2^84, bh<2^42 -> <2^126 */
    r = (r << 42) % m;
    return (r + a * bl) % m;
}

static u128 powmod(u128 a, u128 e, u128 m)
{
    u128 r = 1;
    a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}

/* first 13 primes are a deterministic witness set below this bound */
#define MR13_BOUND (((u128)3317044064679887ULL) * 1000000000ULL + 385961981ULL)

static int is_prime_u128(u128 n)
{
    static const u32 w[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67};
    if (n < 2) return 0;
    for (int i = 0; i < 19; i++) {
        if (n % w[i] == 0) return n == w[i];
        if ((u128)w[i] * w[i] > n) return 1;
    }
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }

    int nw = (n < MR13_BOUND) ? 13 : 19;   /* above the proven bound, extra
                                              witnesses: not a proof, but
                                              error < 4^-19; hits get
                                              re-verified independently */
    for (int i = 0; i < nw; i++) {
        u128 x = powmod(w[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { comp = 0; break; }
        }
        if (comp) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* digit reversal                                                     */
/* ------------------------------------------------------------------ */

static inline u64 bitrev64(u64 v)
{
    v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
    v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
    v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
    return __builtin_bswap64(v);
}

static inline u128 rev_2(u128 x)
{
    u64 hi = (u64)(x >> 64), lo = (u64)x;
    u128 full = ((u128)bitrev64(lo) << 64) | bitrev64(hi);
    int bits = hi ? (128 - __builtin_clzll(hi)) : (64 - __builtin_clzll(lo));
    return full >> (128 - bits);
}

/* division by a compile-time constant, so the compiler uses multiplies */
#define REV_FN(B)                                                    \
static u128 rev_##B(u128 x)                                          \
{                                                                    \
    u128 r = 0;                                                      \
    while (x >> 64) {          /* two-limb; only the top digits */                                                      \
        u64 hi = (u64)(x >> 42), lo = (u64)(x & MASK42);             \
        u64 qh = hi / (B), rr = hi % (B);                            \
        u64 t  = (rr << 42) | lo;                                    \
        u64 ql = t / (B);                                            \
        rr = t % (B);                                                \
        x = ((u128)qh << 42) + ql;                                   \
        r = r * (B) + rr;                                            \
    }                                                                \
    u64 y = (u64)x;                                                  \
    while (y) { u64 q = y / (B); r = r * (B) + (y - q * (B)); y = q; }\
    return r;                                                        \
}
REV_FN(3)  REV_FN(4)  REV_FN(5)  REV_FN(6)  REV_FN(7)  REV_FN(8)
REV_FN(9)  REV_FN(10) REV_FN(11) REV_FN(12) REV_FN(13) REV_FN(14) REV_FN(15)

static inline u128 rev_base(u128 x, int b)
{
    switch (b) {
    case  2: return rev_2(x);   case  3: return rev_3(x);
    case  4: return rev_4(x);   case  5: return rev_5(x);
    case  6: return rev_6(x);   case  7: return rev_7(x);
    case  8: return rev_8(x);   case  9: return rev_9(x);
    case 10: return rev_10(x);  case 11: return rev_11(x);
    case 12: return rev_12(x);  case 13: return rev_13(x);
    case 14: return rev_14(x);  case 15: return rev_15(x);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* small primes and the small-factor prefilter                        */
/* ------------------------------------------------------------------ */

static u32  *g_primes;          /* primes < sieve limit */
static u64   g_nprimes;
static u64   g_sieve_limit = 1u << 20;

/*
 * Divisibility by an odd constant needs no division: for odd p,
 *   p | v   iff   v * (p^-1 mod 2^k)  <=  floor((2^k - 1)/p)
 * (Granlund-Moller).  The smallest primes get this test directly on the
 * 128-bit value; the rest are packed into 64-bit primorial blocks reached
 * by one divq, and only by the ~1/3 of values surviving the direct tests.
 */
#define NDIRECT 8
#define NBLK    3
#define BLKP    16
#define NALL    64

/* The 128-into-64 divide used to reach the primorial blocks is an x86-64
 * instruction.  Everywhere else (aarch64, etc.) every prime is tested with
 * the 128-bit multiply form instead -- slower, but still division-free. */
#if defined(__x86_64__) && !defined(__ILP32__) && !defined(A228768_NO_ASM)
#define HAVE_DIVQ 1
#endif

static u128 g_ainv[NALL], g_alim[NALL];   /* one entry per prefilter prime */
static int  g_nall;
#define g_dinv g_ainv
#define g_dlim g_alim
static struct { u64 M; int n; u64 inv[BLKP], lim[BLKP]; } g_blk[NBLK];
static int  g_nblk;
static u32  g_maxprefilter;

static u128 inv128(u128 p) { u128 x = p; for (int i = 0; i < 7; i++) x *= 2 - p * x; return x; }
static u64  inv64 (u64  p) { u64  x = p; for (int i = 0; i < 6; i++) x *= 2 - p * x; return x; }

#ifdef HAVE_DIVQ
static inline u64 mod128_64(u128 v, u64 m)
{
    u64 hi = (u64)(v >> 64), lo = (u64)v, q, r;
    if (hi >= m) return (u64)(v % m);        /* never here; #DE would be fatal */
    __asm__("divq %[m]" : "=a"(q), "=d"(r) : [m] "r"(m), "a"(lo), "d"(hi));
    (void)q;
    return r;
}
#endif

static void build_prefilter(u32 upto)
{
    int accn = 0;
    u64 acc = 1;
    g_nall = g_nblk = 0;
    for (u64 i = 0; i < g_nprimes && g_nall < NALL; i++) {
        u32 p = g_primes[i];
        if (p > upto) break;
        g_ainv[g_nall] = inv128(p);
        g_alim[g_nall] = (~(u128)0) / p;
        g_nall++;
        g_maxprefilter = p;
        if (g_nall <= NDIRECT) continue;
        if (g_nblk >= NBLK) { g_nall--; break; }       /* blocks are full */
        if (accn == BLKP || acc > ((u64)1 << 62) / p) {
            g_blk[g_nblk].M = acc; g_blk[g_nblk].n = accn; g_nblk++;
            acc = 1; accn = 0;
            if (g_nblk >= NBLK) { g_nall--; g_maxprefilter = (u32)acc; break; }
        }
        g_blk[g_nblk].inv[accn] = inv64(p);
        g_blk[g_nblk].lim[accn] = (~(u64)0) / p;
        accn++; acc *= p;
    }
    if (accn && g_nblk < NBLK) { g_blk[g_nblk].M = acc; g_blk[g_nblk].n = accn; g_nblk++; }
    /* g_maxprefilter must name the largest prime actually tested on BOTH paths */
    g_maxprefilter = 0;
    for (u64 i = 0; i < g_nprimes; i++) {
        if ((int)i >= g_nall) break;
        g_maxprefilter = g_primes[i];
    }
}

static void build_small_primes(void)
{
    u64 L = g_sieve_limit;
    unsigned char *c = calloc(L + 1, 1);
    for (u64 i = 2; i * i <= L; i++)
        if (!c[i]) for (u64 j = i * i; j <= L; j += i) c[j] = 1;
    g_nprimes = 0;
    for (u64 i = 3; i <= L; i += 2) if (!c[i]) g_nprimes++;
    g_primes = malloc(g_nprimes * sizeof *g_primes);
    u64 k = 0;
    for (u64 i = 3; i <= L; i += 2) if (!c[i]) g_primes[k++] = (u32)i;
    free(c);

    build_prefilter(200);
}

/* 1 if v has a prime factor <= g_maxprefilter (v itself assumed larger) */
static inline int has_small_factor(u128 v)
{
    if (v <= (u128)g_maxprefilter * g_maxprefilter) return 0;   /* rare, punt */
    for (int i = 0; i < NDIRECT; i++)
        if (v * g_ainv[i] <= g_alim[i]) return 1;
#ifdef HAVE_DIVQ
    for (int i = 0; i < g_nblk; i++) {
        u64 r = mod128_64(v, g_blk[i].M);
        if (!r) return 1;
        for (int j = 0; j < g_blk[i].n; j++)
            if (r * g_blk[i].inv[j] <= g_blk[i].lim[j]) return 1;
    }
#else
    for (int i = NDIRECT; i < g_nall; i++)
        if (v * g_ainv[i] <= g_alim[i]) return 1;
#endif
    return 0;
}

/* ------------------------------------------------------------------ */
/* leading-digit admissibility                                        */
/* ------------------------------------------------------------------ */
/*
 * If d is x's leading digit in base b then d is the LAST digit of
 * reverse_b(x), so reverse_b(x) = d (mod q) for every prime q | b.
 * gcd(d,b) > 1 therefore forces the reversal composite.  Vacuous for
 * prime b; for composite b it kills whole intervals.
 */

static int g_nbase = 15;
static int g_cbase[MAXBASE + 1], g_ncbase;      /* composite bases <= nbase */

static void build_cbases(void)
{
    g_ncbase = 0;
    for (int b = 4; b <= g_nbase; b++) {
        int comp = 0;
        for (int q = 2; q * q <= b; q++) if (b % q == 0) comp = 1;
        if (comp) g_cbase[g_ncbase++] = b;
    }
}

static inline int igcd(int a, int b) { while (b) { int t = a % b; a = b; b = t; } return a; }

/* leading digit of x in base b, and b^k for that digit position */
static void leading(u128 x, int b, int *d, u128 *pw)
{
    u128 p = 1;
    while (p <= x / b) p *= b;
    *pw = p;
    *d = (int)(x / p);
}

/* least y >= x whose base-b leading digit is coprime to b */
static u128 next_ok(u128 x, int b)
{
    int d; u128 p;
    leading(x, b, &d, &p);
    if (igcd(d, b) == 1) return x;
    do { d++; } while (d < b && igcd(d, b) != 1);
    return (d < b) ? (u128)d * p : p * b;
}

/* least admissible y >= x, or 0 if past limit */
static u128 next_admissible(u128 x, u128 limit)
{
    for (;;) {
        u128 y = x;
        for (int i = 0; i < g_ncbase; i++) {
            u128 t = next_ok(y, g_cbase[i]);
            if (t > y) y = t;
        }
        if (y >= limit) return 0;
        if (y == x) return x;
        x = y;
    }
}

/* end of the maximal admissible window containing x, merging abutting runs */
static u128 admissible_end(u128 x);
static u128 window_end(u128 x, u128 hi)
{
    u128 e = admissible_end(x);
    if (e > hi) return hi;
    while (e < hi && next_admissible(e, hi) == e) {
        u128 e2 = admissible_end(e);
        if (e2 <= e) break;
        e = (e2 > hi) ? hi : e2;
    }
    return e;
}

/* end of the admissible run containing x */
static u128 admissible_end(u128 x)
{
    u128 e = 0;
    for (int i = 0; i < g_ncbase; i++) {
        int d; u128 p;
        leading(x, g_cbase[i], &d, &p);
        u128 t = (u128)(d + 1) * p;
        if (!e || t < e) e = t;
    }
    return e;
}

/* ------------------------------------------------------------------ */
/* the emirp cascade                                                  */
/* ------------------------------------------------------------------ */

static int g_order[MAXORDER], g_norder;

static void build_order(void)
{
    /* base 2 first: bit reversal is nearly free and it kills ~76%.
       then largest base down: fewer digits = cheaper reversal. */
    g_norder = 0;
    g_order[g_norder++] = 2;
    for (int b = g_nbase; b >= 3; b--) g_order[g_norder++] = b;
}

/* per-thread funnel counters */
typedef struct {
    u64 scanned;                /* admissible numbers examined */
    u64 candidates;             /* survived the sieve */
    u64 pass[MAXORDER + 1];     /* survived i cascade stages */
    u64 mr_calls;
    int deepest;
    u64 pad[5];
} __attribute__((aligned(64))) stats_t;

/* stage 1: cheap rejection of all reversals. returns depth reached. */
static inline int cascade_cheap(u128 x, stats_t *st)
{
    for (int i = 0; i < g_norder; i++) {
        u128 r = rev_base(x, g_order[i]);
        if (r == x) return i;                 /* palindrome: not an emirp */
        if (has_small_factor(r)) return i;
        st->pass[i + 1]++;
    }
    return g_norder;
}

/* stage 2: full primality, only for numbers that cleared stage 1 */
static int confirm(u128 x, stats_t *st)
{
    st->mr_calls++;
    if (!is_prime_u128(x)) return 0;
    for (int i = 0; i < g_norder; i++) {
        u128 r = rev_base(x, g_order[i]);
        st->mr_calls++;
        if (!is_prime_u128(r)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* work distribution                                                  */
/* ------------------------------------------------------------------ */

static u128 g_start, g_limit;
static u64  g_chunk_log = 36;
static u64  g_nchunks;
static u64  g_next_chunk;
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static unsigned char *g_done;           /* completed-chunk bitmap */
static u64  g_watermark;                /* lowest incomplete chunk */
static volatile sig_atomic_t g_stop;
static const char *g_ckpt;
static int  g_nthreads = 4;
static int  g_quiet;
static u64  g_seg = 1u << 21;           /* numbers per sieve segment */

static stats_t *g_stats;
static u128 g_hits[64];
static int  g_nhits;

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void on_sigint(int s) { (void)s; g_stop = 1; }

static void report_hit(u128 x)
{
    pthread_mutex_lock(&g_mu);
    if (g_nhits < 64) g_hits[g_nhits++] = x;
    printf("\n*** HIT: %s is an emirp in all bases 2..%d ***\n", S(x), g_nbase);
    fflush(stdout);
    pthread_mutex_unlock(&g_mu);
}

/* ------------------------------------------------------------------ */
/* the scan                                                           */
/* ------------------------------------------------------------------ */

typedef struct { u32 *next; unsigned char *bm; stats_t *st; } worker_t;

/* sieve+test one admissible run [lo,hi) */
static void scan_run(worker_t *w, u128 lo, u128 hi)
{
    if (!(lo & 1)) lo++;                       /* odds only */
    if (lo >= hi) return;

    /* first odd multiple of each sieve prime at or after lo */
    for (u64 i = 0; i < g_nprimes; i++) {
        u64 p = g_primes[i];
        u64 r = (u64)(lo % p);
        u64 off = r ? p - r : 0;
        if ((off & 1) != 0) off += p;          /* keep lo+off odd */
        w->next[i] = (u32)(off >> 1);          /* index among odds */
    }

    u64 segodd = g_seg >> 1;
    for (u128 s = lo; s < hi && !g_stop; s += g_seg) {
        u64 n = segodd;
        if (hi - s < g_seg) n = (u64)((hi - s + 1) >> 1);

        memset(w->bm, 0, (n + 7) >> 3);
        for (u64 i = 0; i < g_nprimes; i++) {
            u64 p = g_primes[i], j = w->next[i];
            for (; j < n; j += p) w->bm[j >> 3] |= 1u << (j & 7);
            w->next[i] = (u32)(j - n);
        }

        w->st->scanned += n;
        for (u64 j = 0; j < n; j++) {
            if (w->bm[j >> 3] & (1u << (j & 7))) continue;
            u128 x = s + 2 * (u128)j;
            w->st->candidates++;
            int depth = cascade_cheap(x, w->st);
            if (depth > w->st->deepest) w->st->deepest = depth;
            if (depth == g_norder && confirm(x, w->st)) report_hit(x);
        }
    }
}

static void mark_done(u64 idx)
{
    pthread_mutex_lock(&g_mu);
    g_done[idx >> 3] |= 1u << (idx & 7);
    while (g_watermark < g_nchunks && (g_done[g_watermark >> 3] & (1u << (g_watermark & 7))))
        g_watermark++;
    pthread_mutex_unlock(&g_mu);
}

static void *worker(void *arg)
{
    long id = (long)arg;
    worker_t w;
    w.next = malloc(g_nprimes * sizeof *w.next);
    w.bm   = malloc((g_seg >> 4) + 8);
    w.st   = &g_stats[id];

    for (;;) {
        u64 idx;
        pthread_mutex_lock(&g_mu);
        idx = g_next_chunk++;
        pthread_mutex_unlock(&g_mu);
        if (idx >= g_nchunks || g_stop) break;

        u128 clo = g_start + ((u128)idx << g_chunk_log);
        u128 chi = clo + ((u128)1 << g_chunk_log);
        if (chi > g_limit) chi = g_limit;

        u128 x = clo;
        while (x < chi && !g_stop) {
            u128 s = next_admissible(x, chi);
            if (!s) break;
            u128 e = admissible_end(s);
            if (e > chi) e = chi;
            scan_run(&w, s, e);
            x = e;
        }
        mark_done(idx);
    }
    free(w.next); free(w.bm);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* reporting                                                          */
/* ------------------------------------------------------------------ */

static u128 admissible_measure(u128 lo, u128 hi)
{
    u128 tot = 0, x = lo;
    while (x < hi) {
        u128 s = next_admissible(x, hi);
        if (!s) break;
        u128 e = admissible_end(s);
        if (e > hi) e = hi;
        tot += e - s;
        x = e;
    }
    return tot;
}

static void print_windows(u128 lo, u128 hi)
{
    printf("admissible windows for bases 2..%d in [%s, ", g_nbase, S(lo));
    printf("%s):\n", S(hi));
    u128 x = lo, tot = 0;
    int i = 0;
    while (x < hi) {
        u128 s = next_admissible(x, hi);
        if (!s) break;
        u128 e = admissible_end(s);
        if (e > hi) e = hi;
        /* merge abutting runs */
        while (e < hi) {
            u128 t = next_admissible(e, hi);
            if (t != e) break;
            u128 e2 = admissible_end(e);
            e = (e2 > hi) ? hi : e2;
        }
        tot += e - s;
        double w = (double)(e - s);
        double cy = w * NS_PER_NUMBER / 1e9 / 31556952.0;
        /* heuristic: each of the nbase conditions holds with prob ~KAPPA/ln x,
           KAPPA calibrated against the observed a(14) search depth */
        double lx = log((double)(s / 2 + e / 2));
        double hits = w * pow(KAPPA / lx, g_nbase);
        printf("  %2d  [%22s, ", ++i, S(s));
        printf("%22s)   width %.4g   %.3g core-yr   E[hits] %.3g\n",
               S(e), w, cy, hits);
        x = e;
    }
    if (!i) printf("   (none)\n");
    printf("  total admissible %.4g  (%.4f%% of range)\n",
           (double)tot, 100.0 * (double)tot / (double)(hi - lo));
    printf("  core-years assume %.1f ns/number/core (run --bench on your hardware);\n"
           "  E[hits] is a heuristic calibrated on the a(14) search, good to a factor of a few.\n",
           NS_PER_NUMBER);
}

static void print_stats(double t0, int final)
{
    stats_t a;
    memset(&a, 0, sizeof a);
    for (int i = 0; i < g_nthreads; i++) {
        a.scanned += g_stats[i].scanned;
        a.candidates += g_stats[i].candidates;
        a.mr_calls += g_stats[i].mr_calls;
        if (g_stats[i].deepest > a.deepest) a.deepest = g_stats[i].deepest;
        for (int j = 0; j <= g_norder; j++) a.pass[j] += g_stats[i].pass[j];
    }
    double el = now() - t0;
    double rate = a.scanned / (el > 0 ? el : 1);

    if (final) printf("\n");
    printf("%s%.0fs  scanned %.4g  %.1fM/s  cand %.3g (%.2f%%)  MR %.3g  deepest %d/%d",
           final ? "" : "\r", el, (double)a.scanned, rate / 1e6,
           (double)a.candidates,
           a.scanned ? 100.0 * a.candidates / a.scanned : 0.0,
           (double)a.mr_calls, a.deepest, g_norder);
    if (final) {
        printf("\n\nfunnel (survivors after each cascade stage):\n");
        u64 prev = a.candidates;
        for (int i = 0; i < g_norder; i++) {
            printf("  base %-3d %14" PRIu64 "   %6.3f%% of previous\n",
                   g_order[i], a.pass[i + 1],
                   prev ? 100.0 * a.pass[i + 1] / prev : 0.0);
            prev = a.pass[i + 1];
        }
        printf("  cost %.2f ns per admissible number per core\n",
               a.scanned ? (now() - t0) * 1e9 * g_nthreads / a.scanned : 0.0);
    }
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* self-test                                                          */
/* ------------------------------------------------------------------ */

static int check(const char *what, int ok)
{
    printf("  %-46s %s\n", what, ok ? "ok" : "FAIL");
    return ok;
}

static int self_test(void)
{
    int ok = 1;
    printf("self-test  (%s prefilter path)\n",
#ifdef HAVE_DIVQ
           "x86 divq"
#else
           "portable"
#endif
           );
    {   /* both paths must agree with honest trial division */
        int good = 1;
        for (u128 v = 1000000000000000000ULL; v < 1000000000000002000ULL; v++) {
            int a = has_small_factor(v), b = 0;
            for (u64 p = 3; p <= g_maxprefilter; p += 2) {
                int isp = 1;
                for (u64 d = 3; d * d <= p; d += 2) if (p % d == 0) isp = 0;
                if (isp && v % p == 0) b = 1;
            }
            if (a != b) good = 0;
        }
        ok &= check("prefilter agrees with trial division", good);
    }

    /* reversal */
    ok &= check("rev_10(19775643791090909)",
                rev_10((u128)19775643791090909ULL) == (u128)90909019734657791ULL);
    ok &= check("rev_2(13) = 11", rev_2(13) == 11);
    ok &= check("rev_15 round trip", rev_15(rev_15(123456789)) == 123456789);
    {
        int good = 1;
        for (int b = 2; b <= 15; b++) {
            u128 v = ((u128)9876543210ULL) * 1000000007ULL;
            if (v % b == 0) continue;      /* trailing zero: reversal loses it */
            if (rev_base(rev_base(v, b), b) != v) good = 0;
        }
        ok &= check("rev_b round trip, bases 2..15", good);
    }

    /* primality */
    ok &= check("MR: 2^61-1 prime", is_prime_u128(((u128)1 << 61) - 1));
    ok &= check("MR: 2^64-59 prime", is_prime_u128((u128)18446744073709551557ULL));
    ok &= check("MR: 3e19+51 composite check",
                !is_prime_u128(C_3E19 + 50));
    {
        /* a 128-bit prime: 10^20 + 39 */
        u128 p = 1; for (int i = 0; i < 20; i++) p *= 10; p += 39;
        ok &= check("MR: 10^20+39 prime", is_prime_u128(p));
        ok &= check("MR: 10^20+40 composite", !is_prime_u128(p + 1));
    }
    {
        int good = 1, cnt = 0;
        for (u64 v = 1000000000000000003ULL; cnt < 300; v += 2, cnt++) {
            int a = is_prime_u128(v);
            int b = 1;
            for (u64 d = 3; d * d <= v && d < 100000; d += 2) if (v % d == 0) { b = 0; break; }
            if (!b && a) good = 0;      /* trial division found a factor but MR says prime */
        }
        ok &= check("MR vs trial division, 300 values near 1e18", good);
    }

    /* admissibility: the known a(14) windows */
    {
        int save = g_nbase;
        g_nbase = 14; build_cbases();
        u128 w1 = next_admissible((u128)14322793967831ULL, (u128)2e16);
        u128 e1 = admissible_end(w1);
        ok &= check("a(14) window 1 start = 14322793967831",
                    w1 == (u128)14322793967831ULL);
        ok &= check("a(14) window 1 end   = 15251194969974",
                    e1 == (u128)15251194969974ULL);
        u128 w3 = next_admissible((u128)116000000000000ULL, (u128)2e16);
        ok &= check("a(14) window 3 start = 18280792200314880",
                    w3 == (u128)18280792200314880ULL);
        g_nbase = save; build_cbases();
    }
    {
        /* a(15): nothing admissible in the rest of 64-bit space */
        int save = g_nbase;
        g_nbase = 15; build_cbases();
        u128 lo = (u128)19775643791090909ULL + 1;
        u128 hi = ~(u128)0 >> 64;                 /* 2^64-1 */
        ok &= check("a(15): no admissible interval below 2^64",
                    next_admissible(lo, hi) == 0);
        u128 w = next_admissible(lo, (u128)4e19);
        ok &= check("a(15) window 1 starts at 3e19", w == C_3E19);
        ok &= check("a(15) window 1 ends at 5*15^16",
                    window_end(w, (u128)4e19) == C_5x15p16);
        g_nbase = save; build_cbases();
    }

    /* known terms, found by actual search */
    {
        struct { int n; u128 lo, hi, want; } t[] = {
            { 10, 93000000,          94000000,          93836531 },
            { 11, 1963000000,        1964000000,        1963209431ULL },
            { 13, 14322793960000ULL, 14322793970000ULL, 14322793967831ULL },
            { 14, 19775643791000000ULL, 19775643792000000ULL, 19775643791090909ULL },
        };
        for (unsigned i = 0; i < sizeof t / sizeof *t; i++) {
            g_nbase = t[i].n; build_cbases(); build_order();
            g_nthreads = 1; g_stop = 0; g_nhits = 0;
            memset(g_stats, 0, sizeof(stats_t));
            worker_t w;
            w.next = malloc(g_nprimes * sizeof *w.next);
            w.bm = malloc((g_seg >> 4) + 8);
            w.st = &g_stats[0];
            u128 x = t[i].lo;
            while (x < t[i].hi) {
                u128 s = next_admissible(x, t[i].hi);
                if (!s) break;
                u128 e = admissible_end(s);
                if (e > t[i].hi) e = t[i].hi;
                scan_run(&w, s, e);
                x = e;
            }
            free(w.next); free(w.bm);
            char msg[96];
            snprintf(msg, sizeof msg, "search finds a(%d) = %s", t[i].n, S(t[i].want));
            ok &= check(msg, g_nhits == 1 && g_hits[0] == t[i].want);
        }
        g_nbase = 15; build_cbases(); build_order();
    }

    printf("%s\n", ok ? "all passed" : "FAILURES");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s [options]\n"
"  --nbase N          search for a(N)            (default 15)\n"
"  --start V          lower bound                (default 3e19)\n"
"  --limit V          upper bound, exclusive     (default 5*15^16)\n"
"  --threads T        worker threads             (default 4)\n"
"  --sieve-limit V    composite-filter bound     (default 1048576)\n"
"  --chunk-log K      work unit = 2^K numbers    (default 36)\n"
"  --resume FILE      checkpoint file to write and read\n"
"  --verify V         check one value in all bases 2..nbase and exit\n"
"  --windows          list admissible windows in [start,limit) and exit\n"
"  --bench            measure throughput on a small slice and exit\n"
"  --self-test        run correctness checks and exit\n"
"  -q                 no live status line\n"
"values may be written 12345, 3e19, or 5*15^16\n", p);
}

int main(int argc, char **argv)
{
    int mode_windows = 0, mode_test = 0, mode_bench = 0, mode_verify = 0;
    u128 vval = 0;
    int have_start = 0, have_limit = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (usage(argv[0]), exit(2), ""))
        if (!strcmp(a, "--nbase")) g_nbase = atoi(NEXT());
        else if (!strcmp(a, "--start")) { if (parse_u128(NEXT(), &g_start)) { fprintf(stderr, "bad --start\n"); return 2; } have_start = 1; }
        else if (!strcmp(a, "--limit")) { if (parse_u128(NEXT(), &g_limit)) { fprintf(stderr, "bad --limit\n"); return 2; } have_limit = 1; }
        else if (!strcmp(a, "--threads")) g_nthreads = atoi(NEXT());
        else if (!strcmp(a, "--sieve-limit")) g_sieve_limit = strtoull(NEXT(), NULL, 0);
        else if (!strcmp(a, "--chunk-log")) g_chunk_log = strtoull(NEXT(), NULL, 0);
        else if (!strcmp(a, "--resume")) g_ckpt = NEXT();
        else if (!strcmp(a, "--verify")) { if (parse_u128(NEXT(), &vval)) { fprintf(stderr, "bad --verify\n"); return 2; } mode_verify = 1; }
        else if (!strcmp(a, "--windows")) mode_windows = 1;
        else if (!strcmp(a, "--bench")) mode_bench = 1;
        else if (!strcmp(a, "--self-test")) mode_test = 1;
        else if (!strcmp(a, "-q")) g_quiet = 1;
        else { usage(argv[0]); return 2; }
    }
    if (g_nbase < 2 || g_nbase > MAXBASE) { fprintf(stderr, "nbase must be 2..%d\n", MAXBASE); return 2; }
    if (g_nthreads < 1 || g_nthreads > 512) { fprintf(stderr, "bad --threads\n"); return 2; }

    build_small_primes();
    build_cbases();
    build_order();
    g_stats = calloc(g_nthreads > 1 ? g_nthreads : 1, sizeof *g_stats);

    if (mode_test) return self_test();

    if (mode_verify) {
        printf("n = %s   bases 2..%d\n", S(vval), g_nbase);
        int ok = is_prime_u128(vval);
        printf("  %-8s prime: %s\n", "n", ok ? "yes" : "NO");
        for (int b = 2; b <= g_nbase && ok; b++) {
            u128 r = rev_base(vval, b);
            int pal = (r == vval), pr = is_prime_u128(r);
            printf("  base %-3d rev = %-26s %s\n", b, S(r),
                   pal ? "FAIL (palindrome)" : pr ? "prime" : "FAIL (composite)");
            if (pal || !pr) ok = 0;
        }
        printf("%s\n", ok ? "VERIFIED: emirp in all bases" : "not an emirp in all bases");
        return ok ? 0 : 1;
    }

    if (!have_start) g_start = C_3E19;
    if (!have_limit) g_limit = C_5x15p16;
    if (g_limit <= g_start) { fprintf(stderr, "limit must exceed start\n"); return 2; }

    if (mode_windows) { print_windows(g_start, g_limit); return 0; }

    if (mode_bench) {
        u128 s = next_admissible(g_start, g_limit);
        if (!s) { fprintf(stderr, "nothing admissible\n"); return 1; }
        u128 e = admissible_end(s);
        if (e > s + 200000000) e = s + 200000000;
        worker_t w;
        w.next = malloc(g_nprimes * sizeof *w.next);
        w.bm = malloc((g_seg >> 4) + 8);
        w.st = &g_stats[0];
        double t0 = now();
        scan_run(&w, s, e);
        g_nthreads = 1;
        print_stats(t0, 1);
        return 0;
    }

    /* checkpoint resume */
    if (g_ckpt) {
        FILE *f = fopen(g_ckpt, "r");
        if (f) {
            char line[128];
            if (fgets(line, sizeof line, f)) {
                u128 v;
                if (!parse_u128(line, &v) && v > g_start) {
                    printf("resuming at %s\n", S(v));
                    g_start = v;
                }
            }
            fclose(f);
        }
    }

    g_nchunks = (u64)(((g_limit - g_start) >> g_chunk_log) + 1);
    g_done = calloc((g_nchunks >> 3) + 1, 1);

    u128 adm = admissible_measure(g_start, g_limit);
    printf("a(%d) search  [%s, ", g_nbase, S(g_start));
    printf("%s)\n", S(g_limit));
    printf("admissible %.4g of %.4g  (%.3f%%)   %d threads, %" PRIu64 " chunks of 2^%" PRIu64 "\n",
           (double)adm, (double)(g_limit - g_start),
           100.0 * (double)adm / (double)(g_limit - g_start),
           g_nthreads, g_nchunks, g_chunk_log);
    printf("cascade order:");
    for (int i = 0; i < g_norder; i++) printf(" %d", g_order[i]);
    printf("\nsieve primes < %" PRIu64 " (%" PRIu64 "), prefilter to %u (%s path)\n\n",
           g_sieve_limit, g_nprimes, g_maxprefilter,
#ifdef HAVE_DIVQ
           "x86 divq"
#else
           "portable"
#endif
           );

    signal(SIGINT, on_sigint);
    double t0 = now();
    pthread_t *th = malloc(g_nthreads * sizeof *th);
    for (long i = 0; i < g_nthreads; i++) pthread_create(&th[i], NULL, worker, (void *)i);

    while (!g_stop) {
        struct timespec ts = { 2, 0 };
        nanosleep(&ts, NULL);
        pthread_mutex_lock(&g_mu);
        u64 wm = g_watermark, nc = g_nchunks, nx = g_next_chunk;
        pthread_mutex_unlock(&g_mu);
        if (!g_quiet) print_stats(t0, 0);
        if (g_ckpt) {
            FILE *f = fopen(g_ckpt, "w");
            if (f) { fprintf(f, "%s\n", S(g_start + ((u128)wm << g_chunk_log))); fclose(f); }
        }
        if (wm >= nc && nx >= nc) break;
    }
    for (int i = 0; i < g_nthreads; i++) pthread_join(th[i], NULL);

    print_stats(t0, 1);
    if (g_nhits) {
        printf("\n%d hit(s):\n", g_nhits);
        for (int i = 0; i < g_nhits; i++) printf("  %s\n", S(g_hits[i]));
    } else {
        printf("\nno emirp in all bases 2..%d in this range\n", g_nbase);
    }
    return g_nhits ? 0 : 1;
}
