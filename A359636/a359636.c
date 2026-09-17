/*
 * a359636.c
 *
 * Compute and extend OEIS A359636:
 *
 *   "a(n) is the least odd prime not in A001359 such that all subsequent
 *    composites in the gap up to the next prime have at least n distinct
 *    prime factors."
 *
 * Known terms (n = 1..8):
 *   7, 19, 643, 51427, 8083633, 1077940147, 75582271489, 34710483181813
 * Upper bounds (David A. Corneth, Jan 2023):
 *   a(9) <= 76340177205657727, a(10) <= 225096507194749219819.
 *
 * Definitions
 * -----------
 * Let p < q be consecutive primes with q - p >= 4 (p is not the lesser member
 * of a twin prime pair, i.e. p is not in A001359).  The gap qualifies for level
 * n when omega(x) >= n for every x in p+1 .. q-1, where omega = A001221 counts
 * distinct prime factors.  a(n) is the smallest p whose gap qualifies.
 *
 * Every known term is a prime gap of 4: p, p+4 prime with p+2 = 0 (mod 3).
 * For example a(8) = 34710483181813:
 *   34710483181814 = 2*11*17*29*47*229*409*727
 *   34710483181815 = 3^2*5*7*19*31*101*263*7043
 *   34710483181816 = 2^3*13*41*43*61*83*139*269
 * and 34710483181817 is prime.
 *
 * Method
 * ------
 * Key reduction.  Every qualifying gap contains an odd multiple of 3, call it
 * m, such that m-1, m, m+1 are all composites of the gap:  a gap of 4 has
 * p = 1 (mod 3), so m = p+2;  a gap of 6 has m = p+2 or p+4;  a gap of 8 or
 * more contains three consecutive odd numbers, one of which is 0 mod 3.
 * Hence it suffices to find every m = 3 (mod 6) with
 *
 *      omega(m-1) >= n,  omega(m) >= n,  omega(m+1) >= n,          (*)
 *
 * and for each such m to look at the prime gap around it.  a(n) is the least
 * p = prevprime(m-1) over all m satisfying (*) whose whole gap qualifies.
 *
 * Sieve.  Index k = (m-3)/6 and sieve the three targets N_t = 6k+3+t
 * (t = -1, 0, +1) simultaneously: for each prime 5 <= p <= T and each target,
 * the k with p | N_t form one residue class mod p, and a byte counter per
 * target (three counters packed in one 32-bit word per k) is incremented along
 * it.  The primes 5..13 and 17..23 are applied as precomputed periodic
 * patterns (periods 5005 and 7429), the rest by ordinary strided marking on
 * L1-sized segments.  The forced divisors 2 | m-1, 2 | m+1, 3 | m are not
 * counted; the thresholds are lowered by one instead.
 *
 * Why a small T suffices.  If N <= Bmax has n distinct prime factors, at most
 * rmax of them can exceed T, where rmax is the largest r with
 * (product of the first n-r primes) * (product of the first r primes above T)
 * <= Bmax.  So N has at least c = n - rmax prime factors <= T, and a k whose
 * three counters are all >= c-1 is a "survivor".  For n = 9, Bmax ~ 7.7e16 and
 * T = 1000 one gets c = 5: every candidate has five prime factors below 1000
 * in each of m-1, m, m+1.  Survivors are rare (order 1e-5), so the sieve with
 * ~140 small primes is the whole cost: about 2 marks per k.
 *
 * Survivors are verified exactly: omega(N_t) >= n by trial division with the
 * same size bound as an early exit (a cofactor needing r more primes but
 * smaller than the r-th power of the next prime fails at once; a cofactor
 * needing exactly one more prime factor succeeds if it exceeds 1; a cofactor
 * needing two is settled by Miller-Rabin plus a perfect-power test).  For the
 * rare k satisfying (*), p = prevprime(m-1) and q = nextprime(m+1) are found
 * with deterministic Miller-Rabin and every composite in (p, q) is checked.
 *
 * The range is cut into chunks handed to worker threads through an atomic
 * counter; an ordered frontier (lowest unfinished chunk) drives the progress
 * line, the checkpoint and the early stop: once the frontier has passed the
 * smallest solution found, that solution is a(n) and the scan ends (unless
 * -a asks for all solutions in the range).
 *
 * Usage
 * -----
 *   a359636 scan N [START] END [-t T] [-c CHUNK] [-T SIEVEMAX] [-S STATE] [-i SECS] [-a] [-q]
 *       Search m = 3 (mod 6) with START <= m <= END+2 for level N; prints every
 *       triple satisfying (*), every qualifying gap, and a(N) if it is in range.
 *       -S FILE keeps a checkpoint (every SECS seconds, default 60, and on
 *       Ctrl-C); rerunning the same command resumes from it.
 *   a359636 hunt N END [-L L] [-Q Q] [-t T]
 *       Structured search for upper bounds: tests every m = 3^e * (N-1 distinct
 *       primes in [5, L]) <= END+2 (default L = 1000); with -Q the largest of
 *       those primes may instead be any prime <= Q.  Cheap, not exhaustive.
 *   a359636 verify N P
 *       Check directly that the gap after the prime P qualifies for level N,
 *       printing the factorization of every composite in the gap.
 *   a359636 omega X ...
 *       Factor the given numbers (trial division + Pollard rho) and print omega.
 *   a359636 selftest [-t T] [-8]
 *       Reproduce a(1..7) (a few seconds), and with -8 also a(8) (~10 min).
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, or X+Y / X-Y of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a359636.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm -o a359636
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <primesieve.h>
#ifdef __APPLE__
#include <pthread/qos.h>
#endif

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

#define NKNOWN 8
static const u64 KNOWN[NKNOWN + 1] = { 0, 7, 19, 643, 51427, 8083633, 1077940147,
                                       75582271489ULL, 34710483181813ULL };

#ifndef SEG
#define SEG      16384          /* k per sieve segment: 64 KB of u32 counters (fits L1) */
#endif
#define PAT_A    5005           /* 5*7*11*13 */
#define PAT_B    7429           /* 17*19*23  */
#define PAT_C    33263          /* 29*31*37  */
#ifndef PATC
#define PATC 1                  /* apply 29,31,37 as a pattern too (else strided) */
#endif
#define FIRST_STRIDED (PATC ? 41 : 29)
#define RING     65536          /* chunk completion ring; bounds how far threads run ahead */
#define MAXSOL   4096

static bool g_quiet = false;
static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;

/* ------------------------------------------------------------------ */
/* Generic helpers                                                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fputs("a359636: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory (%zu bytes)", n);
    return p;
}

/* Parse "123", "2^40", "10^12", "1e12", and X+Y / X-Y of those. */
static u64 parse_num(const char *s)
{
    const char *plus = strchr(s + 1, '+'), *minus = strchr(s + 1, '-');
    if (plus || minus) {
        const char *op = plus ? plus : minus;
        char *a = strndup(s, (size_t)(op - s));
        u64 x = parse_num(a), y = parse_num(op + 1);
        free(a);
        return plus ? x + y : x - y;
    }
    char *end;
    errno = 0;
    unsigned long long base = strtoull(s, &end, 10);
    if (end == s) die("bad number '%s'", s);
    if (*end == '\0') return base;
    if (*end == '^' || *end == 'e' || *end == 'E') {
        bool pw = (*end == '^');
        unsigned long long e = strtoull(end + 1, &end, 10);
        if (*end != '\0') die("bad number '%s'", s);
        u128 r = 1, b = pw ? base : 10;
        for (unsigned long long i = 0; i < e; i++) {
            r *= b;
            if (r >> 64) die("number '%s' does not fit in 64 bits", s);
        }
        if (!pw) r *= base;
        if (r >> 64) die("number '%s' does not fit in 64 bits", s);
        return (u64)r;
    }
    die("bad number '%s'", s);
    return 0;
}

static void fmt_commas(char *out, u64 v)
{
    char tmp[32];
    int n = snprintf(tmp, sizeof tmp, "%" PRIu64, v), o = 0;
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) out[o++] = ',';
        out[o++] = tmp[i];
    }
    out[o] = 0;
}

static void fmt_time(char *out, double s)
{
    if (s < 0 || s > 1e9) { strcpy(out, "?"); return; }
    int t = (int)s;
    if (t >= 86400) sprintf(out, "%dd%02dh%02dm", t / 86400, t % 86400 / 3600, t % 3600 / 60);
    else if (t >= 3600) sprintf(out, "%dh%02dm%02ds", t / 3600, t % 3600 / 60, t % 60);
    else if (t >= 60) sprintf(out, "%dm%02ds", t / 60, t % 60);
    else sprintf(out, "%ds", t);
}

static void on_sigint(int sig) { (void)sig; g_stop = 1; }

/* ------------------------------------------------------------------ */
/* Arithmetic: mulmod, Miller-Rabin, roots, factoring                  */
/* ------------------------------------------------------------------ */

static inline u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static u64 powmod(u64 a, u64 e, u64 m)
{
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

/* Deterministic for all 64-bit n (bases 2..37 suffice below 3.3e24). */
static bool is_prime(u64 n)
{
    static const u64 bases[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    if (n < 2) return false;
    for (int i = 0; i < 12; i++) {
        if (n % bases[i] == 0) return n == bases[i];
    }
    if (n < 41 * 41) return true;
    u64 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 x = powmod(bases[i], d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s && comp; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) comp = false;
        }
        if (comp) return false;
    }
    return true;
}

/* floor(n^(1/r)) for r >= 1 */
static u64 iroot(u64 n, int r)
{
    if (r == 1) return n;
    if (n < 2) return n;
    u64 x = (u64)pow((double)n, 1.0 / r);
    if (x == 0) x = 1;
    for (;;) {                         /* adjust so that x^r <= n < (x+1)^r */
        u128 p = 1;
        bool over = false;
        for (int i = 0; i < r; i++) { p *= x; if (p > n) { over = true; break; } }
        if (over) { x--; continue; }
        u128 q = 1;
        over = false;
        for (int i = 0; i < r; i++) { q *= (x + 1); if (q > n) { over = true; break; } }
        if (!over) { x++; continue; }
        return x;
    }
}

/* Is n = a^j for some j >= 2?  (n > 1) */
static bool is_perfect_power(u64 n)
{
    for (int j = 2; j <= 63; j++) {
        u64 x = iroot(n, j);
        if (x < 2) return false;
        u128 p = 1;
        for (int i = 0; i < j; i++) p *= x;
        if (p == n) return true;
    }
    return false;
}

/* Prime table for trial division and sieving. */
static u32 *g_primes;          /* all primes up to g_plimit */
static size_t g_nprimes;
static u64 g_plimit;

static void build_primes(u64 limit)
{
    size_t n;
    g_primes = (u32 *)primesieve_generate_primes(2, limit, &n, UINT32_PRIMES);
    if (!g_primes) die("primesieve failed");
    g_nprimes = n;
    g_plimit = limit;
}

/*
 * Does N have at least n distinct prime factors?  Trial division with a size
 * bound: after the primes below p have been removed, a cofactor that still
 * needs r more distinct primes must be at least p^r.  Exact for every N whose
 * cube root does not exceed the prime table (the r = 2 case uses Miller-Rabin
 * and a perfect-power test instead of trial division).
 */
static bool omega_at_least(u64 N, int n)
{
    if (n <= 0) return true;
    if (N < 2) return false;
    int cnt = 0;
    u64 lim = iroot(N, n);            /* smallest prime factor must be <= lim */
    for (size_t i = 0; i < g_nprimes; i++) {
        u64 p = g_primes[i];
        int r = n - cnt;               /* distinct primes still needed */
        if (r == 1) return N > 1;      /* any cofactor > 1 supplies one more prime */
        if (p > lim) return false;     /* cofactor too small for r primes above p */
        if (r == 2 && p * p * p > N) { /* cofactor is q, q^2 or q1*q2 territory */
            /* N > 1 here; N has all prime factors >= p and N < p^3, so N is a
               prime, a prime square or a product of two distinct primes. */
            if (is_prime(N)) return false;
            return !is_perfect_power(N);
        }
        if (N % p == 0) {
            cnt++;
            N /= p;
            while (N % p == 0) N /= p;
            r = n - cnt;
            if (r <= 0) return true;
            if (r == 1) return N > 1;
            lim = iroot(N, r);
        }
    }
    /* Prime table exhausted (only possible for huge N): fall back. */
    if (N == 1) return cnt >= n;
    if (is_prime(N)) return cnt + 1 >= n;
    /* composite cofactor with all prime factors > g_plimit: at least 1 prime,
       at least 2 unless a perfect power */
    return cnt + (is_perfect_power(N) ? 1 : 2) >= n;
}

/* Exact omega and factorization (for reporting): trial division + Pollard rho. */
static u64 rho(u64 n)
{
    if (!(n & 1)) return 2;
    for (u64 c = 1;; c++) {
        u64 x = 2, y = 2, d = 1;
        while (d == 1) {
            x = (mulmod(x, x, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            y = (mulmod(y, y, n) + c) % n;
            u64 diff = x > y ? x - y : y - x;
            /* gcd */
            u64 a = diff, b = n;
            while (b) { u64 t = a % b; a = b; b = t; }
            d = a;
        }
        if (d != n) return d;
    }
}

static int factor_full(u64 n, u64 *pr, int *ex)   /* returns number of distinct primes */
{
    int k = 0;
    for (size_t i = 0; i < g_nprimes && (u64)g_primes[i] * g_primes[i] <= n; i++) {
        u64 p = g_primes[i];
        if (n % p == 0) {
            pr[k] = p; ex[k] = 0;
            while (n % p == 0) { n /= p; ex[k]++; }
            k++;
        }
    }
    /* remaining n has no prime factor <= min(g_plimit, sqrt) */
    u64 stack[64];
    int sp = 0;
    if (n > 1) stack[sp++] = n;
    while (sp) {
        u64 x = stack[--sp];
        if (x == 1) continue;
        if (is_prime(x)) {
            int j;
            for (j = 0; j < k; j++) if (pr[j] == x) { ex[j]++; break; }
            if (j == k) { pr[k] = x; ex[k] = 1; k++; }
            continue;
        }
        u64 d = rho(x);
        stack[sp++] = d;
        stack[sp++] = x / d;
    }
    /* sort by prime */
    for (int i = 1; i < k; i++)
        for (int j = i; j > 0 && pr[j - 1] > pr[j]; j--) {
            u64 tp = pr[j]; pr[j] = pr[j - 1]; pr[j - 1] = tp;
            int te = ex[j]; ex[j] = ex[j - 1]; ex[j - 1] = te;
        }
    return k;
}

static void fmt_factorization(char *out, u64 n)
{
    u64 pr[64]; int ex[64];
    int k = factor_full(n, pr, ex);
    int o = 0;
    for (int i = 0; i < k; i++) {
        o += sprintf(out + o, "%s%" PRIu64, i ? "*" : "", pr[i]);
        if (ex[i] > 1) o += sprintf(out + o, "^%d", ex[i]);
    }
    if (k == 0) sprintf(out + o, "%" PRIu64, n);
}

static int omega_exact(u64 n)
{
    u64 pr[64]; int ex[64];
    return factor_full(n, pr, ex);
}

static u64 prev_prime(u64 x)           /* largest prime < x, x > 3 */
{
    x--;
    if (!(x & 1)) x--;
    while (!is_prime(x)) x -= 2;
    return x;
}

static u64 next_prime(u64 x)           /* smallest prime > x */
{
    x++;
    if (!(x & 1)) x++;
    while (!is_prime(x)) x += 2;
    return x;
}

/* ------------------------------------------------------------------ */
/* Gap check                                                           */
/* ------------------------------------------------------------------ */

/* Does the gap (p, q) qualify at level n?  Writes q. */
static bool gap_qualifies(u64 p, int n, u64 *q_out)
{
    u64 q = next_prime(p);
    *q_out = q;
    if (q - p < 4) return false;
    for (u64 x = p + 1; x < q; x++)
        if (!omega_at_least(x, n)) return false;
    return true;
}

/* ------------------------------------------------------------------ */
/* Scan                                                                */
/* ------------------------------------------------------------------ */

typedef struct {
    u64 p, q, m;
} solution_t;

typedef struct {
    int      n;              /* level */
    u64      k_lo, k_hi;     /* k range [k_lo, k_hi), m = 6k+3 */
    u64      chunk;          /* k per chunk */
    u64      nchunks;
    u32      T;              /* sieve prime bound */
    int      c;              /* min # prime factors <= T in each target */
    int      thr;            /* sieved-count threshold c-1 (may be <= 0) */
    /* sieve primes >= 29 up to T */
    u32     *sp;             /* primes */
    u32     *sr;             /* residues, 3 per prime: k = sr[3i+t] (mod p) marks target t */
    size_t   nsp;
    u32      patA[PAT_A], patB[PAT_B];
    u32     *patC;
    /* statistics */
    atomic_ullong survivors, triples, nsol;
    /* solutions */
    solution_t sol[MAXSOL];
    int      nsolrec;
    u64      best_p, best_m;
    pthread_mutex_t sol_mu;
    /* chunk frontier */
    atomic_ullong next_chunk;   /* next chunk index to hand out */
    u64      frontier;          /* all chunks < frontier are complete */
    u8       done[RING];
    pthread_mutex_t fr_mu;
    pthread_cond_t  fr_cv;
    bool     scan_all;
    bool     stop_now;
    double   t0;
    const char *state_file;
    u64      resume_frontier;
} scan_t;

static u64 inv_mod(u64 a, u64 p)       /* a^-1 mod p, p prime */
{
    return powmod(a % p, p - 2, p);
}

/* residue of k for which p | 6k + 3 + t */
static u32 target_residue(u64 p, int t)
{
    u64 inv6 = inv_mod(6, p);
    u64 c = (u64)(3 + t) % p;               /* 6k = -c  (mod p) */
    return (u32)(((p - c) % p) * inv6 % p);
}

static void setup_sieve(scan_t *S)
{
    /* patterns for 5,7,11,13 and 17,19,23 */
    memset(S->patA, 0, sizeof S->patA);
    memset(S->patB, 0, sizeof S->patB);
    static const u32 pa[4] = { 5, 7, 11, 13 }, pb[3] = { 17, 19, 23 };
    for (int i = 0; i < 4; i++) {
        if (pa[i] > S->T) continue;
        for (int t = -1; t <= 1; t++) {
            u32 r = target_residue(pa[i], t), inc = 1u << (8 * (t + 1));
            for (u32 k = r; k < PAT_A; k += pa[i]) S->patA[k] += inc;
        }
    }
    for (int i = 0; i < 3; i++) {
        if (pb[i] > S->T) continue;
        for (int t = -1; t <= 1; t++) {
            u32 r = target_residue(pb[i], t), inc = 1u << (8 * (t + 1));
            for (u32 k = r; k < PAT_B; k += pb[i]) S->patB[k] += inc;
        }
    }
    S->patC = xmalloc(PAT_C * sizeof(u32));
    memset(S->patC, 0, PAT_C * sizeof(u32));
    if (PATC) {
        static const u32 pc[3] = { 29, 31, 37 };
        for (int i = 0; i < 3; i++) {
            if (pc[i] > S->T) continue;
            for (int t = -1; t <= 1; t++) {
                u32 r = target_residue(pc[i], t), inc = 1u << (8 * (t + 1));
                for (u32 k = r; k < PAT_C; k += pc[i]) S->patC[k] += inc;
            }
        }
    }
    /* strided primes FIRST_STRIDED..T */
    size_t cnt = 0;
    for (size_t i = 0; i < g_nprimes && g_primes[i] <= S->T; i++)
        if (g_primes[i] >= FIRST_STRIDED) cnt++;
    S->sp = xmalloc(cnt * sizeof(u32) + 1);
    S->sr = xmalloc(3 * cnt * sizeof(u32) + 1);
    S->nsp = 0;
    for (size_t i = 0; i < g_nprimes && g_primes[i] <= S->T; i++) {
        u32 p = g_primes[i];
        if (p < FIRST_STRIDED) continue;
        S->sp[S->nsp] = p;
        for (int t = -1; t <= 1; t++) S->sr[3 * S->nsp + (t + 1)] = target_residue(p, t);
        S->nsp++;
    }
}

/* c = n - rmax, rmax = largest r with primorial(n-r) * (first r primes > T) <= Bmax */
static int compute_c(int n, u32 T, u64 Bmax)
{
    int rmax = 0;
    for (int r = 0; r <= n; r++) {
        u128 prod = 1;
        bool over = false;
        for (int i = 0; i < n - r; i++) { prod *= g_primes[i]; if (prod > Bmax) { over = true; break; } }
        size_t j = 0;
        while (j < g_nprimes && g_primes[j] <= T) j++;
        for (int i = 0; i < r && !over; i++) {
            if (j + i >= g_nprimes) die("prime table too small in compute_c");
            prod *= g_primes[j + i];
            if (prod > Bmax) over = true;
        }
        if (over) break;
        rmax = r;
    }
    return n - rmax;
}

static void record_solution(scan_t *S, u64 p, u64 q, u64 m)
{
    pthread_mutex_lock(&S->sol_mu);
    bool dup = false;
    for (int i = 0; i < S->nsolrec; i++) if (S->sol[i].p == p) dup = true;
    if (!dup) {
        if (S->nsolrec < MAXSOL) {
            S->sol[S->nsolrec].p = p;
            S->sol[S->nsolrec].q = q;
            S->sol[S->nsolrec].m = m;
            S->nsolrec++;
        }
        atomic_fetch_add(&S->nsol, 1);
        if (!S->best_p || p < S->best_p) { S->best_p = p; S->best_m = m; }
        pthread_mutex_lock(&out_mu);
        if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
        printf("SOLUTION n=%d p=%" PRIu64 " q=%" PRIu64 " gap=%" PRIu64 " (m=%" PRIu64 ")\n",
               S->n, p, q, q - p, m);
        char buf[256];
        for (u64 x = p + 1; x < q; x++) {
            fmt_factorization(buf, x);
            printf("    %" PRIu64 " = %s  omega=%d\n", x, buf, omega_exact(x));
        }
        fflush(stdout);
        pthread_mutex_unlock(&out_mu);
    }
    pthread_mutex_unlock(&S->sol_mu);
}

/* Full check of a survivor k: (*) and then the gap. */
static void check_survivor(scan_t *S, u64 k)
{
    u64 m = 6 * k + 3;
    int n = S->n;
    if (!omega_at_least(m - 1, n)) return;
    if (!omega_at_least(m, n)) return;
    if (!omega_at_least(m + 1, n)) return;
    atomic_fetch_add(&S->triples, 1);
    u64 p = prev_prime(m - 1), q;
    pthread_mutex_lock(&out_mu);
    if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
    printf("triple   n=%d m=%" PRIu64 "  omega(m-1,m,m+1) = %d %d %d  prevprime=%" PRIu64 " (%+" PRId64 ")\n",
           n, m, omega_exact(m - 1), omega_exact(m), omega_exact(m + 1), p, (int64_t)(p - m));
    fflush(stdout);
    pthread_mutex_unlock(&out_mu);
    if (gap_qualifies(p, n, &q)) record_solution(S, p, q, m);
}

static void save_state(scan_t *S);

static void progress(scan_t *S, bool final)
{
    static double last = 0;
    double t = now();
    if (!final && t - last < (stderr_tty ? 1.0 : 60.0)) return;
    last = t;
    u64 fr = S->frontier;
    u64 kdone = fr * S->chunk;
    if (kdone > S->k_hi - S->k_lo) kdone = S->k_hi - S->k_lo;
    u64 mnow = 6 * (S->k_lo + kdone) + 3;
    double frac = (double)kdone / (double)(S->k_hi - S->k_lo);
    double el = t - S->t0;
    double rate = el > 0 ? (double)(kdone - (S->resume_frontier * S->chunk)) / el : 0;
    char c1[32], c2[32], eta[32], elap[32];
    fmt_commas(c1, mnow);
    fmt_time(eta, rate > 0 ? (double)((S->k_hi - S->k_lo) - kdone) / rate : -1);
    fmt_time(elap, el);
    snprintf(c2, sizeof c2, "%.2f", rate * 6 / 1e9);
    char best[48] = "";
    if (S->best_p) snprintf(best, sizeof best, " best=%" PRIu64, S->best_p);
    fprintf(stderr, "%sm=%s (%.2f%%)  %s Gm/s  surv=%llu trip=%llu sol=%llu%s  %s elapsed, ETA %s%s",
            stderr_tty ? "\r\033[K" : "", c1, frac * 100, c2,
            (unsigned long long)atomic_load(&S->survivors),
            (unsigned long long)atomic_load(&S->triples),
            (unsigned long long)atomic_load(&S->nsol),
            best, elap, eta,
            stderr_tty && !final ? "" : "\n");
    fflush(stderr);
}

static void chunk_done(scan_t *S, u64 idx)
{
    pthread_mutex_lock(&S->fr_mu);
    S->done[idx % RING] = 1;
    while (S->frontier < S->nchunks && S->done[S->frontier % RING]) {
        S->done[S->frontier % RING] = 0;
        S->frontier++;
    }
    /* early stop: frontier past the chunk holding the best m */
    if (S->best_p && !S->scan_all) {
        u64 kbest = (S->best_m - 3) / 6;
        u64 cbest = (kbest - S->k_lo) / S->chunk;
        if (S->frontier > cbest) S->stop_now = true;
    }
    pthread_cond_broadcast(&S->fr_cv);
    pthread_mutex_unlock(&S->fr_mu);
}

typedef struct { scan_t *S; int id; } worker_arg;

static void *worker(void *arg)
{
    scan_t *S = ((worker_arg *)arg)->S;
#ifdef __APPLE__
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);   /* prefer the performance cores */
#endif
    u32 *W = xmalloc(SEG * sizeof(u32) + 64);
    u32 *nextj = xmalloc(3 * S->nsp * sizeof(u32) + 4);
    const u32 thrmask = 0x808080u;
    const u32 thradd = (S->thr > 0) ? 0x010101u * (u32)(128 - S->thr) : 0x808080u;
    enum { INC1 = 0x000001u, INC2 = 0x000100u, INC3 = 0x010000u };

    for (;;) {
        u64 idx = atomic_fetch_add(&S->next_chunk, 1);
        if (idx >= S->nchunks || g_stop || S->stop_now) break;
        /* don't run more than RING-1 chunks ahead of the frontier */
        pthread_mutex_lock(&S->fr_mu);
        while (idx >= S->frontier + RING - 1 && !g_stop && !S->stop_now)
            pthread_cond_wait(&S->fr_cv, &S->fr_mu);
        pthread_mutex_unlock(&S->fr_mu);

        u64 k0 = S->k_lo + idx * S->chunk;
        u64 k1 = k0 + S->chunk;
        if (k1 > S->k_hi) k1 = S->k_hi;

        /* first offsets for the strided primes */
        for (size_t i = 0; i < S->nsp; i++) {
            u32 p = S->sp[i];
            u32 k0m = (u32)(k0 % p);
            for (int t = 0; t < 3; t++) {
                u32 r = S->sr[3 * i + t];
                nextj[3 * i + t] = r >= k0m ? r - k0m : r + p - k0m;
            }
        }
        u32 a = (u32)(k0 % PAT_A), b = (u32)(k0 % PAT_B), c = (u32)(k0 % PAT_C);

        for (u64 ks = k0; ks < k1; ks += SEG) {
            u32 len = (u32)(k1 - ks < SEG ? k1 - ks : SEG);
            /* pattern init: W = patA + patB (rotated), one fused pass */
            for (u32 j = 0; j < len;) {
                u32 run = PAT_A - a;
                if (run > PAT_B - b) run = PAT_B - b;
                if (run > len - j) run = len - j;
                const u32 *pa = S->patA + a, *pb = S->patB + b;
                u32 *w = W + j;
                for (u32 i = 0; i < run; i++) w[i] = pa[i] + pb[i];
                j += run; a += run; b += run;
                if (a == PAT_A) a = 0;
                if (b == PAT_B) b = 0;
            }
            if (PATC) {
                for (u32 j = 0; j < len;) {
                    u32 run = PAT_C - c; if (run > len - j) run = len - j;
                    const u32 *pc = S->patC + c;
                    u32 *w = W + j;
                    for (u32 i = 0; i < run; i++) w[i] += pc[i];
                    j += run; c += run; if (c == PAT_C) c = 0;
                }
            }
            /* strided marks: the three targets of a prime share the stride p,
               so after aligning their phases they are marked in one loop */
            for (size_t i = 0; i < S->nsp; i++) {
                const u32 p = S->sp[i];
                u32 j1 = nextj[3 * i], j2 = nextj[3 * i + 1], j3 = nextj[3 * i + 2];
                if (len >= 2 * p) {
                    if (j2 < j1) { W[j2] += INC2; j2 += p; }
                    if (j3 < j1) { W[j3] += INC3; j3 += p; }
                    const u32 d2 = j2 - j1, d3 = j3 - j1;
                    const u32 end = len - (d2 > d3 ? d2 : d3);
                    u32 j = j1;
                    for (; j + 3 * p < end; j += 4 * p) {
                        W[j] += INC1;             W[j + d2] += INC2;             W[j + d3] += INC3;
                        W[j + p] += INC1;         W[j + p + d2] += INC2;         W[j + p + d3] += INC3;
                        W[j + 2 * p] += INC1;     W[j + 2 * p + d2] += INC2;     W[j + 2 * p + d3] += INC3;
                        W[j + 3 * p] += INC1;     W[j + 3 * p + d2] += INC2;     W[j + 3 * p + d3] += INC3;
                    }
                    for (; j < end; j += p) { W[j] += INC1; W[j + d2] += INC2; W[j + d3] += INC3; }
                    j1 = j; j2 = j + d2; j3 = j + d3;
                }
                for (; j1 < len; j1 += p) W[j1] += INC1;
                for (; j2 < len; j2 += p) W[j2] += INC2;
                for (; j3 < len; j3 += p) W[j3] += INC3;
                nextj[3 * i] = j1 - len; nextj[3 * i + 1] = j2 - len; nextj[3 * i + 2] = j3 - len;
            }
            /* scan: a k survives when all three counters are >= thr, i.e. when
               (~(W + thradd)) & 0x808080 == 0; blocks are tested with a min-reduction */
            if (S->thr <= 0) {
                for (u32 j = 0; j < len; j++) { atomic_fetch_add(&S->survivors, 1); check_survivor(S, ks + j); }
            } else {
                u32 j = 0;
                for (; j + 16 <= len; j += 16) {
                    u32 mn = 0xffffffffu;
                    for (u32 i = 0; i < 16; i++) {
                        u32 t = (~(W[j + i] + thradd)) & thrmask;
                        mn = t < mn ? t : mn;
                    }
                    if (mn == 0) {
                        for (u32 i = 0; i < 16; i++)
                            if ((((~(W[j + i] + thradd)) & thrmask)) == 0) {
                                atomic_fetch_add(&S->survivors, 1);
                                check_survivor(S, ks + j + i);
                            }
                    }
                }
                for (; j < len; j++) {
                    if ((((~(W[j] + thradd)) & thrmask)) == 0) {
                        atomic_fetch_add(&S->survivors, 1);
                        check_survivor(S, ks + j);
                    }
                }
            }
        }
        chunk_done(S, idx);
        if (((worker_arg *)arg)->id == 0 && !g_quiet) progress(S, false);
    }
    free(W);
    free(nextj);
    return NULL;
}

/* ---- checkpoint ---- */

static void save_state(scan_t *S)
{
    if (!S->state_file) return;
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", S->state_file);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "a359636: cannot write %s\n", tmp); return; }
    pthread_mutex_lock(&S->fr_mu);
    u64 fr = S->frontier;
    pthread_mutex_unlock(&S->fr_mu);
    fprintf(f, "A359636 state v1\n");
    fprintf(f, "n %d\nk_lo %" PRIu64 "\nk_hi %" PRIu64 "\nchunk %" PRIu64 "\nfrontier %" PRIu64 "\n",
            S->n, S->k_lo, S->k_hi, S->chunk, fr);
    fprintf(f, "survivors %llu\ntriples %llu\n",
            (unsigned long long)atomic_load(&S->survivors), (unsigned long long)atomic_load(&S->triples));
    pthread_mutex_lock(&S->sol_mu);
    for (int i = 0; i < S->nsolrec; i++)
        fprintf(f, "solution %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", S->sol[i].p, S->sol[i].q, S->sol[i].m);
    pthread_mutex_unlock(&S->sol_mu);
    fclose(f);
    rename(tmp, S->state_file);
}

static void load_state(scan_t *S)
{
    FILE *f = fopen(S->state_file, "r");
    if (!f) return;
    char line[512];
    if (!fgets(line, sizeof line, f) || strncmp(line, "A359636 state v1", 16)) die("bad state file");
    int n = -1; u64 klo = 0, khi = 0, chunk = 0, fr = 0;
    unsigned long long surv = 0, trip = 0;
    while (fgets(line, sizeof line, f)) {
        u64 p, q, m;
        if (sscanf(line, "n %d", &n) == 1) continue;
        if (sscanf(line, "k_lo %" SCNu64, &klo) == 1) continue;
        if (sscanf(line, "k_hi %" SCNu64, &khi) == 1) continue;
        if (sscanf(line, "chunk %" SCNu64, &chunk) == 1) continue;
        if (sscanf(line, "frontier %" SCNu64, &fr) == 1) continue;
        if (sscanf(line, "survivors %llu", &surv) == 1) continue;
        if (sscanf(line, "triples %llu", &trip) == 1) continue;
        if (sscanf(line, "solution %" SCNu64 " %" SCNu64 " %" SCNu64, &p, &q, &m) == 3) {
            if (S->nsolrec < MAXSOL) {
                S->sol[S->nsolrec].p = p; S->sol[S->nsolrec].q = q; S->sol[S->nsolrec].m = m; S->nsolrec++;
            }
            atomic_fetch_add(&S->nsol, 1);
            if (!S->best_p || p < S->best_p) { S->best_p = p; S->best_m = m; }
        }
    }
    fclose(f);
    if (n != S->n || klo != S->k_lo || khi != S->k_hi || chunk != S->chunk)
        die("state file %s belongs to a different scan (n=%d k_lo=%" PRIu64 " k_hi=%" PRIu64 " chunk=%" PRIu64 ")",
            S->state_file, n, klo, khi, chunk);
    S->frontier = fr;
    S->resume_frontier = fr;
    atomic_store(&S->next_chunk, fr);
    atomic_store(&S->survivors, surv);
    atomic_store(&S->triples, trip);
    if (!g_quiet) {
        char c1[32];
        fmt_commas(c1, 6 * (klo + fr * chunk) + 3);
        fprintf(stderr, "resuming from %s: m = %s (%.2f%%), %d solution(s) on file\n",
                S->state_file, c1, 100.0 * (double)(fr * chunk) / (double)(khi - klo), S->nsolrec);
    }
}

/* Run the scan; returns a(n) if certain (0 otherwise). */
static u64 run_scan(int n, u64 start, u64 end, int nthreads, u64 chunk, u32 T, const char *state_file,
                    int interval, bool scan_all, bool *certain)
{
    scan_t *S = xmalloc(sizeof *S);
    memset(S, 0, sizeof *S);
    S->n = n;
    S->scan_all = scan_all;
    S->state_file = state_file;
    pthread_mutex_init(&S->sol_mu, NULL);
    pthread_mutex_init(&S->fr_mu, NULL);
    pthread_cond_init(&S->fr_cv, NULL);

    /* m = 3 (mod 6), start <= m <= end + 2 */
    u64 m_lo = start < 9 ? 9 : start;
    while (m_lo % 6 != 3) m_lo++;
    u64 m_hi = end + 2;                   /* inclusive */
    if (m_hi < m_lo) die("empty range");
    S->k_lo = (m_lo - 3) / 6;
    S->k_hi = (m_hi - 3) / 6 + 1;
    u64 Bmax = m_hi + 1;
    if (Bmax > UINT64_MAX - 16) die("range too large for 64-bit arithmetic");

    /* prime table: sieve primes up to T, trial division up to cbrt(Bmax) with margin */
    u64 plim = (u64)iroot(Bmax, 3) + 1000;
    if (plim < 1u << 20) plim = 1u << 20;
    if (plim < 60000) plim = 60000;
    if (!g_primes || g_plimit < plim) { free(g_primes); build_primes(plim); }

    /* T and c: with T = 0 pick the smallest sieve bound that still forces
       min(n,5) prime factors <= T into every target (survivors ~1e-5 then) */
    if (T == 0) {
        static const u32 cand[] = { 100, 150, 200, 250, 300, 400, 500, 600, 700, 800, 1000, 1250, 1500,
                                    2000, 2500, 3000, 4000, 5000, 7000, 10000, 15000, 20000, 30000, 50000 };
        int want = n < 5 ? n : 5;
        T = cand[sizeof cand / sizeof *cand - 1];
        for (size_t i = 0; i < sizeof cand / sizeof *cand; i++)
            if (compute_c(n, cand[i], Bmax) >= want) { T = cand[i]; break; }
    }
    S->T = T;
    S->c = compute_c(n, T, Bmax);
    S->thr = S->c - 1;
    setup_sieve(S);

    if (chunk == 0) chunk = 1u << 26;
    if (chunk % SEG) chunk += SEG - chunk % SEG;
    S->chunk = chunk;
    S->nchunks = (S->k_hi - S->k_lo + chunk - 1) / chunk;

    if (state_file) load_state(S);

    if (!g_quiet) {
        char c1[32], c2[32];
        fmt_commas(c1, m_lo); fmt_commas(c2, m_hi);
        double marks = 0;
        for (size_t i = 0; i < S->nsp; i++) marks += 3.0 / S->sp[i];
        fprintf(stderr, "A359636 level n=%d: m in [%s, %s] (m = 3 mod 6, %" PRIu64 " values), "
                "%d threads, chunk %" PRIu64 "\n", n, c1, c2, S->k_hi - S->k_lo, nthreads, chunk);
        fprintf(stderr, "sieve primes 5..%u (%zu strided), each target needs >= %d prime factors <= %u "
                "(threshold %d after the forced 2/3), %.2f marks per k\n",
                S->T, S->nsp, S->c, S->T, S->thr, marks);
    }

    S->t0 = now();
    pthread_t *th = xmalloc(nthreads * sizeof *th);
    worker_arg *wa = xmalloc(nthreads * sizeof *wa);
    for (int i = 0; i < nthreads; i++) {
        wa[i].S = S; wa[i].id = i;
        if (pthread_create(&th[i], NULL, worker, &wa[i])) die("pthread_create");
    }
    /* main thread: checkpoints */
    double last_save = now();
    for (;;) {
        bool all_done;
        pthread_mutex_lock(&S->fr_mu);
        all_done = S->frontier >= S->nchunks || S->stop_now;
        if (!all_done && !g_stop) {
            struct timespec ts;
            clock_gettime(CLOCK_REALTIME, &ts);
            ts.tv_nsec += 200000000L;
            if (ts.tv_nsec >= 1000000000L) { ts.tv_sec++; ts.tv_nsec -= 1000000000L; }
            pthread_cond_timedwait(&S->fr_cv, &S->fr_mu, &ts);
        }
        pthread_mutex_unlock(&S->fr_mu);
        if (all_done || g_stop) break;
        if (state_file && now() - last_save >= interval) { save_state(S); last_save = now(); }
    }
    if (g_stop) { pthread_mutex_lock(&S->fr_mu); pthread_cond_broadcast(&S->fr_cv); pthread_mutex_unlock(&S->fr_mu); }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    if (state_file) save_state(S);
    if (!g_quiet) progress(S, true);

    u64 covered_k = S->frontier >= S->nchunks ? S->k_hi : S->k_lo + S->frontier * S->chunk;
    u64 covered_m = covered_k ? 6 * covered_k - 3 : 0;   /* all m < 6*covered_k+3 done */
    u64 result = 0;
    *certain = false;
    if (S->best_p && S->best_m <= covered_m) { result = S->best_p; *certain = true; }
    if (!g_quiet) {
        char c1[32];
        fmt_commas(c1, covered_m);
        fprintf(stderr, "scanned m <= %s: %llu survivors, %llu triples, %llu solutions\n", c1,
                (unsigned long long)atomic_load(&S->survivors),
                (unsigned long long)atomic_load(&S->triples),
                (unsigned long long)atomic_load(&S->nsol));
    }
    if (*certain) {
        printf("a(%d) = %" PRIu64 "\n", n, result);
    } else if (S->best_p) {
        printf("a(%d) <= %" PRIu64 " (scan incomplete below it)\n", n, S->best_p);
    } else if (S->frontier >= S->nchunks) {
        printf("a(%d) > %" PRIu64 " (no qualifying gap with p <= %" PRIu64 ")\n", n, end, end);
    } else {
        printf("a(%d) > %" PRIu64 " (scan stopped early; no qualifying gap with p below that)\n", n,
               covered_m > 2 ? covered_m - 2 : 0);
    }
    fflush(stdout);
    free(S->sp); free(S->sr); free(S->patC); free(th); free(wa); free(S);
    return result;
}


/* ------------------------------------------------------------------ */
/* Hunt: structured search for upper bounds                            */
/* ------------------------------------------------------------------ */
/*
 * Enumerate m = 3^e * p_1 * ... * p_(n-1) (e = 1, 2; distinct primes
 * 5 <= p_i <= L) with m <= END+2 and test omega(m-1) >= n, omega(m+1) >= n
 * and then the gap, exactly as in the scan.  This covers only the m whose
 * own prime factors are all small, so it can miss solutions, but it visits
 * about 1e11 candidates for n = 9, L = 1000 below Corneth's bound instead of
 * 1.3e16 sieve positions, and any hit is an upper bound for a(n) and a place
 * where the exhaustive scan may stop.  Work is split over threads by the
 * first two primes chosen.
 */
typedef struct {
    int      n;
    u64      mmax;             /* largest m considered (END + 2) */
    u32     *pr; int npr;      /* primes 5..L */
    u64      Q;                /* if > 0: the last prime is any prime <= Q (not just <= L) */
    u64      ntasks;
    atomic_ullong next_task, tested, triples, nsol;
    u64      best_p, best_m;
    pthread_mutex_t mu;
    double   t0;
} hunt_t;

static void hunt_leaf(hunt_t *H, u64 m, u64 *tested)
{
    (*tested)++;
    if (!omega_at_least(m - 1, H->n)) return;
    if (!omega_at_least(m + 1, H->n)) return;
    atomic_fetch_add(&H->triples, 1);
    u64 p = prev_prime(m - 1), q;
    pthread_mutex_lock(&out_mu);
    if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
    printf("triple   n=%d m=%" PRIu64 "  omega(m-1,m,m+1) = %d %d %d  prevprime=%" PRIu64 " (%+" PRId64 ")\n",
           H->n, m, omega_exact(m - 1), omega_exact(m), omega_exact(m + 1), p, (int64_t)(p - m));
    fflush(stdout);
    pthread_mutex_unlock(&out_mu);
    if (!gap_qualifies(p, H->n, &q)) return;
    pthread_mutex_lock(&H->mu);
    atomic_fetch_add(&H->nsol, 1);
    if (!H->best_p || p < H->best_p) { H->best_p = p; H->best_m = m; }
    pthread_mutex_lock(&out_mu);
    if (stderr_tty && !g_quiet) fputs("\r\033[K", stderr);
    printf("SOLUTION n=%d p=%" PRIu64 " q=%" PRIu64 " gap=%" PRIu64 " (m=%" PRIu64 ")\n", H->n, p, q, q - p, m);
    char buf[256];
    for (u64 x = p + 1; x < q; x++) {
        fmt_factorization(buf, x);
        printf("    %" PRIu64 " = %s  omega=%d\n", x, buf, omega_exact(x));
    }
    fflush(stdout);
    pthread_mutex_unlock(&out_mu);
    pthread_mutex_unlock(&H->mu);
}

/* prod already contains 3^e and `depth` primes; choose the rest from pr[start..] */
static void hunt_dfs(hunt_t *H, int start, int depth, u64 prod, u64 *tested)
{
    int need = H->n - 1 - depth;
    if (need == 0) { hunt_leaf(H, prod, tested); return; }
    if (need == 1 && H->Q) {
        /* last prime free: every prime q with (last chosen prime) < q <= min(Q, mmax/prod) */
        u64 last = start > 0 ? H->pr[start - 1] : 4;
        u64 qmax = H->mmax / prod;
        if (qmax > H->Q) qmax = H->Q;
        size_t lo = 0, hi = g_nprimes;                 /* first prime > last */
        while (lo < hi) { size_t mid = (lo + hi) / 2; if (g_primes[mid] <= last) lo = mid + 1; else hi = mid; }
        for (size_t i = lo; i < g_nprimes && g_primes[i] <= qmax; i++) hunt_leaf(H, prod * g_primes[i], tested);
        return;
    }
    for (int i = start; i < H->npr; i++) {
        u128 q = (u128)prod * H->pr[i];
        u128 mn = q;                                   /* cheapest completion */
        for (int j = 1; j < need; j++) {
            if (i + j >= H->npr) { mn = (u128)H->mmax + 1; break; }
            mn *= H->pr[i + j];
        }
        if (mn > H->mmax) break;
        hunt_dfs(H, i + 1, depth + 1, (u64)q, tested);
    }
}

static void hunt_progress(hunt_t *H, bool final)
{
    static double last = 0;
    double t = now();
    if (!final && t - last < (stderr_tty ? 1.0 : 60.0)) return;
    last = t;
    u64 done = atomic_load(&H->next_task);
    if (done > H->ntasks) done = H->ntasks;
    char el[32], best[48] = "";
    fmt_time(el, t - H->t0);
    if (H->best_p) snprintf(best, sizeof best, " best=%" PRIu64, H->best_p);
    fprintf(stderr, "%stasks %" PRIu64 "/%" PRIu64 " (%.1f%%)  tested=%llu trip=%llu sol=%llu%s  %s elapsed%s",
            stderr_tty ? "\r\033[K" : "", done, H->ntasks, 100.0 * (double)done / (double)H->ntasks,
            (unsigned long long)atomic_load(&H->tested), (unsigned long long)atomic_load(&H->triples),
            (unsigned long long)atomic_load(&H->nsol), best, el, stderr_tty && !final ? "" : "\n");
    fflush(stderr);
}

typedef struct { hunt_t *H; int id; } hunt_arg;

static void *hunt_worker(void *arg)
{
    hunt_t *H = ((hunt_arg *)arg)->H;
#ifdef __APPLE__
    pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
#endif
    u64 npairs = (u64)H->npr * (H->npr - 1) / 2;
    for (;;) {
        u64 t = atomic_fetch_add(&H->next_task, 1);
        if (t >= H->ntasks || g_stop) break;
        int e = (int)(t / npairs);                    /* 0: 3^1, 1: 3^2 */
        u64 pair = t % npairs;
        /* pair -> (i, j), i < j, lexicographic */
        int i = 0;
        u64 acc = 0;
        while (acc + (u64)(H->npr - 1 - i) <= pair) { acc += (u64)(H->npr - 1 - i); i++; }
        int j = i + 1 + (int)(pair - acc);
        u64 base = e ? 9 : 3;
        u128 prod = (u128)base * H->pr[i] * H->pr[j];
        u64 tested = 0;
        if (prod <= H->mmax) hunt_dfs(H, j + 1, 2, (u64)prod, &tested);
        atomic_fetch_add(&H->tested, tested);
        if (((hunt_arg *)arg)->id == 0 && !g_quiet) hunt_progress(H, false);
    }
    return NULL;
}

static int cmd_hunt(int n, u64 end, u32 L, u64 Q, int nthreads)
{
    if (n < 3) die("hunt needs N >= 3");
    hunt_t *H = xmalloc(sizeof *H);
    memset(H, 0, sizeof *H);
    H->n = n;
    H->mmax = end + 2;
    H->Q = Q;
    pthread_mutex_init(&H->mu, NULL);
    u64 plim = (u64)iroot(H->mmax, 3) + 1000;
    if (plim < 1u << 20) plim = 1u << 20;
    if (plim < L + 1000) plim = L + 1000;
    if (plim < Q + 1000) plim = Q + 1000;
    build_primes(plim);
    H->pr = xmalloc((g_nprimes + 1) * sizeof(u32));
    for (size_t i = 0; i < g_nprimes && g_primes[i] <= L; i++)
        if (g_primes[i] >= 5) H->pr[H->npr++] = g_primes[i];
    if (H->npr < n - 1) die("too few primes below L=%u for n=%d", L, n);
    H->ntasks = 2 * ((u64)H->npr * (H->npr - 1) / 2);
    if (!g_quiet) {
        char c1[32];
        fmt_commas(c1, H->mmax);
        fprintf(stderr, "A359636 hunt n=%d: m = 3^e * (%d primes in [5,%u])%s <= %s, %d threads, %" PRIu64 " tasks\n",
                n, Q ? n - 2 : n - 1, L, Q ? " * (one prime <= Q)" : "", c1, nthreads, H->ntasks);
        if (Q) fprintf(stderr, "Q = %" PRIu64 "\n", Q);
    }
    H->t0 = now();
    pthread_t *th = xmalloc(nthreads * sizeof *th);
    hunt_arg *ha = xmalloc(nthreads * sizeof *ha);
    for (int i = 0; i < nthreads; i++) {
        ha[i].H = H; ha[i].id = i;
        if (pthread_create(&th[i], NULL, hunt_worker, &ha[i])) die("pthread_create");
    }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    if (!g_quiet) hunt_progress(H, true);
    if (!g_quiet)
        fprintf(stderr, "hunt tested %llu candidates: %llu triples, %llu solutions\n",
                (unsigned long long)atomic_load(&H->tested), (unsigned long long)atomic_load(&H->triples),
                (unsigned long long)atomic_load(&H->nsol));
    if (H->best_p) printf("a(%d) <= %" PRIu64 " (hunt; m = %" PRIu64 ")\n", n, H->best_p, H->best_m);
    else printf("hunt found no qualifying gap with m <= %" PRIu64 " (n=%d, L=%u)\n", H->mmax, n, L);
    fflush(stdout);
    free(H->pr); free(th); free(ha); free(H);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n < 1 ? 1 : (int)n;
}

static int cmd_verify(int n, u64 p)
{
    build_primes(1u << 21);
    char buf[256];
    if (!is_prime(p)) { printf("%" PRIu64 " is not prime\n", p); return 1; }
    u64 q = next_prime(p);
    printf("p = %" PRIu64 ", next prime q = %" PRIu64 ", gap %" PRIu64 "\n", p, q, q - p);
    bool ok = q - p >= 4;
    if (!ok) printf("p is in A001359 (twin prime), excluded\n");
    for (u64 x = p + 1; x < q; x++) {
        fmt_factorization(buf, x);
        int w = omega_exact(x);
        printf("  %" PRIu64 " = %s  omega=%d%s\n", x, buf, w, w >= n ? "" : "  <-- fails");
        if (w < n) ok = false;
    }
    printf("gap after %" PRIu64 " %s for level %d\n", p, ok ? "QUALIFIES" : "does not qualify", n);
    return ok ? 0 : 1;
}

static int cmd_selftest(int nthreads, bool with8)
{
    int fails = 0;
    int top = with8 ? 8 : 7;
    for (int n = 1; n <= top; n++) {
        bool certain;
        double t = now();
        u64 r = run_scan(n, 0, KNOWN[n], nthreads, n <= 6 ? 1u << 20 : 1u << 26, 0, NULL, 60, false, &certain);
        bool ok = certain && r == KNOWN[n];
        printf("selftest n=%d: got %" PRIu64 "%s expected %" PRIu64 " -> %s  (%.1fs)\n", n, r,
               certain ? "" : " (uncertain)", KNOWN[n], ok ? "ok" : "FAIL", now() - t);
        fflush(stdout);
        if (!ok) fails++;
    }
    /* spot checks of omega_at_least against exact factoring */
    build_primes(1u << 21);
    u64 tests[] = { 8, 9, 10, 20, 21, 22, 644, 645, 646, 51428, 51429, 51430, 8083634, 8083635,
                    1077940148ULL, 75582271490ULL, 34710483181814ULL, 34710483181815ULL, 34710483181816ULL,
                    76340177205657728ULL, 76340177205657729ULL, 76340177205657730ULL,
                    223092870ULL, 6469693230ULL, 200560490130ULL, 1000000007ULL * 998244353ULL,
                    (u64)1000003 * 1000003, 6ULL * 1000003 * 1000033 };
    for (size_t i = 0; i < sizeof tests / sizeof *tests; i++) {
        int w = omega_exact(tests[i]);
        for (int n = 1; n <= 12; n++) {
            bool got = omega_at_least(tests[i], n), exp = w >= n;
            if (got != exp) { printf("omega_at_least(%" PRIu64 ", %d) = %d, expected %d FAIL\n", tests[i], n, got, exp); fails++; }
        }
    }
    printf("selftest %s\n", fails ? "FAILED" : "passed");
    return fails ? 1 : 0;
}

static void usage(void)
{
    fputs("usage:\n"
          "  a359636 scan N [START] END [-t THREADS] [-c CHUNK] [-T SIEVEMAX] [-S STATE] [-i SECS] [-a] [-q]\n"
          "  a359636 hunt N END [-L L] [-Q Q] [-t THREADS] [-q]\n"
          "  a359636 verify N P\n"
          "  a359636 omega X ...\n"
          "  a359636 selftest [-t THREADS] [-8]\n", stderr);
    exit(2);
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    if (argc < 2) usage();
    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
    const char *cmd = argv[1];

    if (!strcmp(cmd, "selftest")) {
        int nthreads = default_threads();
        bool with8 = false;
        for (int i = 2; i < argc; i++) {
            if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-8")) with8 = true;
            else if (!strcmp(argv[i], "-q")) g_quiet = true;
            else usage();
        }
        return cmd_selftest(nthreads, with8);
    }
    if (!strcmp(cmd, "hunt")) {
        if (argc < 4) usage();
        int n = atoi(argv[2]);
        u64 end = 0; u32 L = 1000; u64 Q = 0; int nthreads = default_threads(); int npos = 0;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "-L") && i + 1 < argc) L = (u32)parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-Q") && i + 1 < argc) Q = parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-q")) g_quiet = true;
            else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
            else if (npos == 0) { end = parse_num(argv[i]); npos++; }
            else usage();
        }
        if (!npos) usage();
        if (nthreads < 1) nthreads = 1;
        return cmd_hunt(n, end, L, Q, nthreads);
    }
    if (!strcmp(cmd, "verify")) {
        if (argc != 4) usage();
        return cmd_verify(atoi(argv[2]), parse_num(argv[3]));
    }
    if (!strcmp(cmd, "omega")) {
        build_primes(1u << 21);
        char buf[256];
        for (int i = 2; i < argc; i++) {
            u64 x = parse_num(argv[i]);
            fmt_factorization(buf, x);
            printf("%" PRIu64 " = %s  omega=%d%s\n", x, buf, omega_exact(x), is_prime(x) ? " (prime)" : "");
        }
        return 0;
    }
    if (!strcmp(cmd, "scan")) {
        if (argc < 4) usage();
        int n = atoi(argv[2]);
        if (n < 1 || n > 20) die("N must be 1..20");
        u64 pos[2]; int npos = 0;
        int nthreads = default_threads();
        u64 chunk = 0; u32 T = 0; const char *state = NULL; int interval = 60; bool all = false;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-T") && i + 1 < argc) T = (u32)parse_num(argv[++i]);
            else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
            else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-a")) all = true;
            else if (!strcmp(argv[i], "-q")) g_quiet = true;
            else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
            else if (npos < 2) pos[npos++] = parse_num(argv[i]);
            else usage();
        }
        if (npos == 0) usage();
        u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
        if (nthreads < 1) nthreads = 1;
        if (T && T < 30) T = 30;
        bool certain;
        run_scan(n, start, end, nthreads, chunk, T, state, interval, all, &certain);
        return 0;
    }
    usage();
    return 2;
}
