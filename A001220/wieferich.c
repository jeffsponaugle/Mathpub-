/*
 * wieferich.c
 *
 * Multi-base search for Wieferich primes to base b: primes p with
 *
 *     b^(p-1) == 1 (mod p^2).
 *
 * OEIS has one sequence per base: A001220 (b=2), A014127 (3, Mirimanoff
 * primes), A123692 (5), A212583 (6), A123693 (7), A045616 (10), A111027 (12),
 * A128667 (13), A234810 (14), A242741 (15), A128668 (17), A244260 (18),
 * A090968 (19), A242982 (20), A298951 (22), A128669 (23), A306255 (26),
 * A306256 (30), plus A039951 (the smallest such p for each base n).  Base 2
 * is searched to 2^64 (PrimeGrid, Dec 2022), bases 3, 5, 7 to 1.2e15 and the
 * other bases up to 149 to 2.0e14 (R. Fischer's table, fermatquotient.com,
 * 2025), so the other bases are where new terms can realistically be found.
 *
 * Heuristic
 * ---------
 * For p not dividing b, Fermat gives b^(p-1) = 1 + q p (mod p^2) with the
 * Fermat quotient q = q_p(b) in [0, p); p is a solution iff q = 0.  If q
 * behaves like a random residue, a prime qualifies with probability 1/p, so
 * a search of [x, y] finds ln(ln y / ln x) solutions per base on average:
 * 0.17 from 2e14 to 1e17, but only 0.015 for base 2 from 2^64 to 2^65.
 * Writing q as A in (-p/2, p/2], a "near-Wieferich" prime has |A| <= T; about
 * (2T+1) ln(ln y / ln x) of those appear per base (-n T), a live sanity check.
 *
 * Method
 * ------
 * Primes come from a primesieve iterator per worker (a few percent of the
 * time).  For each prime p < 2^63 the modulus N = p^2 < 2^126 is set up once
 * and shared by all bases:
 *
 *   pinv = p^-1 mod 2^64 (Newton), ninv = -N^-1 = -pinv^2 mod 2^64,
 *   one  = 2^128 mod N (Montgomery form of 1), invN = 1/N as a double.
 *
 * b^(p-1) mod N is computed with Montgomery arithmetic (R = 2^128, two 64-bit
 * limbs) and a fixed window of w exponent bits: w squarings (a dedicated
 * 2-limb squaring with 7 limb products) followed by one multiplication by the
 * plain constant c = b^d, where d is the next window digit and c < 2^32.
 * That small multiplication needs no Montgomery step: t = x*c has at most 158
 * bits, floor(t/N) < 2^32 is estimated with one double multiplication (error
 * below 1) and fixed up by one add or subtract of N.  w is the largest value
 * <= 5 with b^(2^w - 1) < 2^32: 5 for b = 2, 4 for b = 3, 3 for b <= 23, 2 for
 * b <= 1625, 1 above.
 *
 * One exponentiation is a long chain of dependent multiply-adds, so a worker
 * runs LANES (default 8) of them in lock step with the same squaring
 * schedule, which lets the out-of-order core overlap the chains.  Lanes are
 * (p, b) pairs taken from a batch of 64 primes times the bases that share
 * the same w; shorter exponents in a lane group just see leading zero digits.
 *
 * The result r = b^(p-1) mod N is taken out of Montgomery form and
 *
 *   q = (r - 1) * pinv mod 2^64      (exact division: p divides r - 1).
 *
 * If r were wrong, q would be an essentially random 64-bit number, so the
 * test "q < p" re-checks Fermat's little theorem and catches arithmetic or
 * hardware faults with probability 1 - p/2^64 (99.99% at 1e15); failures are
 * reported as FLT errors.  Every q goes into a per-base checksum
 *
 *   cks_b = sum over p of mix64(q ^ p*0x9E3779B97F4A7C15) mod 2^64
 *
 * (mix64 = the splitmix64 finalizer), which does not depend on the order of
 * the primes, so runs with other chunk sizes, thread counts or the CUDA tool
 * give identical values.  verify_wieferich.py computes the same sums with
 * Python integers and its own Miller-Rabin prime list.  Pairs with p | b are
 * skipped; p = 2 is a solution iff b = 1 (mod 4) and is handled directly.
 * Neither enters the checksums.
 *
 * The range is cut into chunks handed out through an atomic counter; chunk
 * results are folded in chunk order, so everything below the frontier is
 * complete and the checkpoint records the frontier.  A solution found in
 * chunk c is CONFIRMED as the smallest one in [START, ...) once the frontier
 * passes c.  With -L FILE one line per chunk (range, prime count, checksums)
 * is appended in frontier order, for comparing independent runs.
 *
 * Usage
 * -----
 *   wieferich scan [START] END [-b BASES] [-t T] [-c CHUNK] [-n T] [-S FILE] [-i SECS] [-L FILE] [-q]
 *       Test every prime in [START, END) (END <= 2^63, rounded up to a chunk
 *       boundary) for every base.  Prints FOUND / NEAR / FLT-ERROR lines as
 *       they occur and a summary with the checksums at the end.  -c CHUNK
 *       (default 10^9), -n T reports near-Wieferich primes with |A| <= T,
 *       -S FILE keeps a checkpoint (every SECS seconds, default 60, and on
 *       Ctrl-C; rerunning the same command resumes, END may change), -L FILE
 *       chunk log, -t T threads (default: all cores), -q no status line.
 *   wieferich check P [-b BASES]
 *       Fermat quotients q_P(b) and A for one prime, from the fast kernel and
 *       from an independent shift-and-add reference.
 *   wieferich bench [N] [-d SPAN] [-b BASES] [-t T]
 *       Time [N, N+SPAN) (defaults 10^15, 10^10) and project longer scans.
 *   wieferich selftest [LIMIT] [-t T]
 *       Arithmetic against the slow reference, a scan of [0, LIMIT) (default
 *       10^9) for all bases 2..30 compared with the known solutions, windows
 *       around the known solutions above LIMIT, and checksums of three windows
 *       (1e15, 1e18, 2^63) against values computed by verify_wieferich.py.
 *
 * BASES is a comma-separated list of bases and ranges (3,5,7 or 2-30); the
 * keyword "oeis" is the 17 bases with OEIS sequences other than 2 (the
 * default), "known" is every base 2..30 that is not a perfect power.  Numbers
 * may be written as decimal, 2^k, 10^k, 1e12, 1.2e15, or X+Y / X-Y of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include wieferich.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm -o wieferich
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
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef unsigned __int128 u128;

#ifndef LANES
#define LANES 8             /* exponentiations run in lock step per worker */
#endif
#define BATCH   64          /* primes set up together; lanes are formed from BATCH x bases */
#define MAXB    64          /* bases per run */
#define MAXW    5           /* largest window */
#define RINGC   4096        /* chunk completion ring: bounds how far workers run ahead */
#define MAXHITS 4096
#define GOLD    0x9E3779B97F4A7C15ULL
#define PMAX    (1ULL << 63) /* p < 2^63 keeps N = p^2 < 2^126 */

static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;

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
    if (stderr_tty) fputs("\r\033[K", stderr);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static const char *fmt_dur(double s, char *buf, size_t n)
{
    if (s < 0 || s != s || s > 1e12) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

static bool mul_ok(u64 a, u64 b, u64 *r)
{
    u128 x = (u128)a * b;
    if (x >> 64) return false;
    *r = (u64)x;
    return true;
}

static bool all_digits(const char *s, size_t n)
{
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

/* decimal, B^E, MeE (M may have a decimal point) */
static bool parse_term(const char *s, size_t n, u64 *out)
{
    char buf[64], *c;
    if (n == 0 || n >= sizeof buf) return false;
    memcpy(buf, s, n);
    buf[n] = 0;
    if ((c = strchr(buf, '^'))) {
        *c = 0;
        if (!all_digits(buf, strlen(buf)) || !all_digits(c + 1, strlen(c + 1))) return false;
        u64 b = strtoull(buf, 0, 10), e = strtoull(c + 1, 0, 10), r = 1;
        while (e--) if (!mul_ok(r, b, &r)) return false;
        *out = r;
        return true;
    }
    if ((c = strpbrk(buf, "eE"))) {
        *c = 0;
        char *dot = strchr(buf, '.');
        size_t ip = dot ? (size_t)(dot - buf) : strlen(buf);
        const char *frac = dot ? dot + 1 : "";
        if ((ip && !all_digits(buf, ip)) || (*frac && !all_digits(frac, strlen(frac)))) return false;
        if (!all_digits(c + 1, strlen(c + 1))) return false;
        long e = strtol(c + 1, 0, 10) - (long)strlen(frac);
        u64 m = 0;
        for (size_t i = 0; i < ip; i++) if (!mul_ok(m, 10, &m) || (m += (u64)(buf[i] - '0')) < (u64)(buf[i] - '0')) return false;
        for (const char *f = frac; *f; f++) if (!mul_ok(m, 10, &m) || (m += (u64)(*f - '0')) < (u64)(*f - '0')) return false;
        for (; e > 0; e--) if (!mul_ok(m, 10, &m)) return false;
        for (; e < 0; e++) { if (m % 10) return false; m /= 10; }
        *out = m;
        return true;
    }
    if (!all_digits(buf, n)) return false;
    errno = 0;
    *out = strtoull(buf, 0, 10);
    return errno == 0;
}

static u64 parse_num(const char *s)
{
    size_t n = strlen(s);
    for (size_t i = n; i-- > 1; ) {           /* split at the last top-level + or - */
        if (s[i] != '+' && s[i] != '-') continue;
        if (s[i - 1] == 'e' || s[i - 1] == 'E') continue;
        u64 a, b;
        if (!parse_term(s, i, &a) || !parse_term(s + i + 1, n - i - 1, &b)) die("bad number '%s'", s);
        if (s[i] == '+') {
            if (a + b < a) die("overflow in '%s'", s);
            return a + b;
        }
        if (b > a) die("negative result in '%s'", s);
        return a - b;
    }
    u64 v;
    if (!parse_term(s, n, &v)) die("bad number '%s'", s);
    return v;
}

static int parse_int(const char *s, int lo, int hi, const char *what)
{
    u64 v = parse_num(s);
    if (v < (u64)lo || v > (u64)hi) die("%s must be between %d and %d", what, lo, hi);
    return (int)v;
}

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n > 0 ? (int)n : 1;
}

static inline u64 mix64(u64 z)
{
    z += GOLD;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* ------------------------------------------------------------------ */
/* Independent primality: deterministic Miller-Rabin below 2^64         */
/* ------------------------------------------------------------------ */

static u64 mulmod64(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }

static u64 powmod64(u64 a, u64 e, u64 m)
{
    u64 r = 1 % m;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod64(r, a, m);
        a = mulmod64(a, a, m);
        e >>= 1;
    }
    return r;
}

static bool is_prime64(u64 n)
{
    static const u64 small[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    if (n < 2) return false;
    for (size_t i = 0; i < sizeof small / sizeof *small; i++) {
        if (n == small[i]) return true;
        if (n % small[i] == 0) return false;
    }
    u64 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (size_t i = 0; i < sizeof small / sizeof *small; i++) {     /* deterministic below 3.3e24 */
        u64 x = powmod64(small[i], d, n);
        if (x == 1 || x == n - 1) continue;
        bool comp = true;
        for (int r = 1; r < s && comp; r++) {
            x = mulmod64(x, x, n);
            if (x == n - 1) comp = false;
        }
        if (comp) return false;
    }
    return true;
}

/* ------------------------------------------------------------------ */
/* Arithmetic mod N = p^2                                               */
/* ------------------------------------------------------------------ */

typedef struct {
    u128 N, one;            /* p^2, 2^128 mod p^2 */
    double invN;            /* 1/N */
    u64 p, pinv, ninv;      /* p, p^-1 mod 2^64, -N^-1 mod 2^64 */
} pctx_t;

static inline void pctx_init(pctx_t *c, u64 p)
{
    u64 x = (3 * p) ^ 2;                        /* 5 correct bits for odd p */
    x *= 2 - p * x; x *= 2 - p * x; x *= 2 - p * x; x *= 2 - p * x;
    c->p = p;
    c->pinv = x;
    c->N = (u128)p * p;
    c->ninv = 0 - x * x;
    c->invN = 1.0 / (double)c->N;
    if (c->N >> 80) {                           /* 2^128 / N < 2^48: estimate the quotient */
        u64 q = (u64)(0x1p128 * c->invN);
        u128 r = (u128)0 - (u128)q * c->N;
        if ((__int128)r < 0) r += c->N;
        else if (r >= c->N) r -= c->N;
        c->one = r;
    } else {
        c->one = ((u128)0 - c->N) % c->N;
    }
}

/* a^2 / 2^128 mod N, a < N < 2^126 */
static inline u128 msqr(u128 a, const pctx_t *m)
{
    const u64 a0 = (u64)a, a1 = (u64)(a >> 64);
    const u64 n0 = (u64)m->N, n1 = (u64)(m->N >> 64);
    const u128 p00 = (u128)a0 * a0, p01 = (u128)a0 * a1, p11 = (u128)a1 * a1;
    const u128 d = p01 << 1;                    /* < 2^127 since a1 < 2^62 */
    const u64 t0 = (u64)p00;
    u128 acc = (u128)(u64)(p00 >> 64) + (u64)d;
    u64 t1 = (u64)acc;
    acc = (acc >> 64) + (u64)(d >> 64) + (u64)p11;
    u64 t2 = (u64)acc;
    u64 t3 = (u64)(acc >> 64) + (u64)(p11 >> 64);
    u64 q = t0 * m->ninv;                       /* two word-by-word REDC steps */
    acc = (u128)q * n0 + t0;
    acc = (u128)q * n1 + t1 + (u64)(acc >> 64);
    t1 = (u64)acc;
    acc = (u128)t2 + (u64)(acc >> 64);
    t2 = (u64)acc;
    t3 += (u64)(acc >> 64);
    q = t1 * m->ninv;
    acc = (u128)q * n0 + t1;
    acc = (u128)q * n1 + t2 + (u64)(acc >> 64);
    const u64 r0 = (u64)acc;
    acc = (u128)t3 + (u64)(acc >> 64);
    const u128 r = ((u128)(u64)acc << 64) | r0;
    return r >= m->N ? r - m->N : r;
}

/* x * c mod N for x < N < 2^126, c < 2^32 (works in or out of Montgomery form) */
static inline u128 msmall(u128 x, u64 c, const pctx_t *m)
{
    const u128 lo = (u128)(u64)x * c;
    const u128 hi = (u128)(u64)(x >> 64) * c;
    const u128 t = lo + (hi << 64);             /* x*c mod 2^128 */
    const u64 top = (u64)(hi >> 64) + (t < lo); /* x*c >> 128 */
    const double td = (double)top * 0x1p128 + (double)(u64)(t >> 64) * 0x1p64 + (double)(u64)t;
    const u64 q = (u64)(td * m->invN);          /* floor(x*c / N) + {-1, 0, 1} */
    u128 r = t - (u128)q * m->N;                /* true value in [-N, 2N): fits in signed 128 */
    if ((__int128)r < 0) r += m->N;
    else if (r >= m->N) r -= m->N;
    return r;
}

/* x / 2^128 mod N (leaves Montgomery form) */
static inline u128 redc(u128 x, const pctx_t *m)
{
    const u64 n0 = (u64)m->N, n1 = (u64)(m->N >> 64);
    u64 t0 = (u64)x, t1 = (u64)(x >> 64);
    u64 q = t0 * m->ninv;
    u128 acc = (u128)q * n0 + t0;
    acc = (u128)q * n1 + t1 + (u64)(acc >> 64);
    t0 = (u64)acc;
    t1 = (u64)(acc >> 64);
    q = t0 * m->ninv;
    acc = (u128)q * n0 + t0;
    acc = (u128)q * n1 + t1 + (u64)(acc >> 64);
    return acc >= m->N ? acc - m->N : acc;
}

/* Reference arithmetic, independent of the Montgomery code: shift and add. */
static u128 mulmod_ref(u128 a, u128 b, u128 m)
{
    u128 r = 0;
    a %= m;
    b %= m;
    while (b) {
        if (b & 1) { r += a; if (r >= m) r -= m; }
        a <<= 1;
        if (a >= m) a -= m;
        b >>= 1;
    }
    return r;
}

static u128 powmod_ref(u64 b, u64 e, u128 m)
{
    u128 r = 1 % m, x = b % m;
    while (e) {
        if (e & 1) r = mulmod_ref(r, x, m);
        x = mulmod_ref(x, x, m);
        e >>= 1;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Bases                                                                */
/* ------------------------------------------------------------------ */

static const u32 OEIS_BASES[] = {3, 5, 6, 7, 10, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 26, 30};
static const char *const OEIS_ANUM[] = {"A014127", "A123692", "A212583", "A123693", "A045616",
    "A111027", "A128667", "A234810", "A242741", "A128668", "A244260", "A090968", "A242982",
    "A298951", "A128669", "A306255", "A306256"};

/* Known solutions (OEIS; R. Fischer's table) for the bases 2..30 that are not perfect powers. */
typedef struct { u32 b; int n; u64 p[8]; } known_t;
static const known_t KNOWN[] = {
    { 2, 2, {1093, 3511}},
    { 3, 2, {11, 1006003}},
    { 5, 7, {2, 20771, 40487, 53471161, 1645333507ULL, 6692367337ULL, 188748146801ULL}},
    { 6, 3, {66161, 534851, 3152573}},
    { 7, 2, {5, 491531}},
    {10, 3, {3, 487, 56598313}},
    {11, 1, {71}},
    {12, 2, {2693, 123653}},
    {13, 3, {2, 863, 1747591}},
    {14, 3, {29, 353, 7596952219ULL}},
    {15, 2, {29131, 119327070011ULL}},
    {17, 5, {2, 3, 46021, 48947, 478225523351ULL}},
    {18, 6, {5, 7, 37, 331, 33923, 1284043}},
    {19, 6, {3, 7, 13, 43, 137, 63061489}},
    {20, 4, {281, 46457, 9377747, 122959073}},
    {21, 1, {2}},
    {22, 5, {13, 673, 1595813, 492366587, 9809862296159ULL}},
    {23, 5, {13, 2481757, 13703077, 15546404183ULL, 2549536629329ULL}},
    {24, 2, {5, 25633}},
    {26, 5, {3, 5, 71, 486999673, 6695256707ULL}},
    {28, 3, {3, 19, 23}},
    {29, 1, {2}},
    {30, 3, {7, 160541, 94727075783ULL}},
};
#define NKNOWN ((int)(sizeof KNOWN / sizeof *KNOWN))

typedef struct {
    int nb;
    u32 b[MAXB];
    int w[MAXB];
    u64 pw[MAXB][1 << MAXW];                    /* b^d for d < 2^w, all < 2^32 */
    int ncls;                                   /* bases grouped by window size */
    int cls_w[MAXW + 1], cls_n[MAXW + 1], cls_j[MAXW + 1][MAXB];
} bases_t;

static bases_t B;

static int window_for(u32 b)
{
    for (int w = MAXW; w > 1; w--) {
        u64 x = 1;
        bool ok = true;
        for (int i = 0; i < (1 << w) - 1 && ok; i++) {
            x *= b;
            if (x >> 32) ok = false;
        }
        if (ok) return w;
    }
    return 1;
}

static int cmp_u32(const void *a, const void *b)
{
    u32 x = *(const u32 *)a, y = *(const u32 *)b;
    return x < y ? -1 : x > y;
}

static void bases_setup(const u32 *list, int n)
{
    u32 tmp[MAXB];
    if (n < 1 || n > MAXB) die("between 1 and %d bases are supported", MAXB);
    memcpy(tmp, list, (size_t)n * sizeof *tmp);
    qsort(tmp, (size_t)n, sizeof *tmp, cmp_u32);
    memset(&B, 0, sizeof B);
    for (int i = 0; i < n; i++) {
        if (tmp[i] < 2) die("bases must be >= 2");
        if (B.nb && B.b[B.nb - 1] == tmp[i]) continue;
        int j = B.nb++;
        B.b[j] = tmp[i];
        B.w[j] = window_for(tmp[i]);
        u64 x = 1;
        for (int d = 0; d < (1 << B.w[j]); d++) { B.pw[j][d] = x; x *= tmp[i]; }
    }
    for (int w = MAXW; w >= 1; w--) {
        int c = B.ncls, k = 0;
        for (int j = 0; j < B.nb; j++) if (B.w[j] == w) B.cls_j[c][k++] = j;
        if (k) { B.cls_w[c] = w; B.cls_n[c] = k; B.ncls++; }
    }
}

static void parse_bases(const char *s)
{
    u32 list[MAXB];
    int n = 0;
    char buf[1024];
    if (strlen(s) >= sizeof buf) die("base list too long");
    strcpy(buf, s);
    for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
        if (!strcmp(tok, "oeis")) {
            for (size_t i = 0; i < sizeof OEIS_BASES / sizeof *OEIS_BASES; i++) {
                if (n >= MAXB) die("too many bases (max %d)", MAXB);
                list[n++] = OEIS_BASES[i];
            }
            continue;
        }
        if (!strcmp(tok, "known")) {
            for (int i = 0; i < NKNOWN; i++) {
                if (n >= MAXB) die("too many bases (max %d)", MAXB);
                list[n++] = KNOWN[i].b;
            }
            continue;
        }
        char *dash = strchr(tok, '-');
        u64 lo, hi;
        if (dash) { *dash = 0; lo = parse_num(tok); hi = parse_num(dash + 1); }
        else lo = hi = parse_num(tok);
        if (lo < 2 || hi < lo || hi >= (1ULL << 32)) die("bad base or range '%s'", tok);
        for (u64 b = lo; b <= hi; b++) {
            if (n >= MAXB) die("too many bases (max %d)", MAXB);
            list[n++] = (u32)b;
        }
    }
    bases_setup(list, n);
}

static const char *anum_for(u32 b)
{
    if (b == 2) return "A001220";
    for (size_t i = 0; i < sizeof OEIS_BASES / sizeof *OEIS_BASES; i++)
        if (OEIS_BASES[i] == b) return OEIS_ANUM[i];
    return NULL;
}

static void print_bases(FILE *f)
{
    for (int j = 0; j < B.nb; j++) fprintf(f, "%s%u", j ? "," : "", B.b[j]);
}

/* ------------------------------------------------------------------ */
/* The kernel                                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    const pctx_t *c;
    const u64 *pw;
    u64 e;                  /* p - 1 */
    int j;                  /* base index */
} lane_t;

/* r[g] = b_g^(p_g - 1) mod p_g^2, normal form */
static inline void run_group(const lane_t *L, int w, u128 *r)
{
    const u64 mask = ((u64)1 << w) - 1;
    u128 x[LANES];
    u64 emax = 0;
    for (int g = 0; g < LANES; g++) emax |= L[g].e;
    int sh = ((64 - __builtin_clzll(emax) + w - 1) / w - 1) * w;
    for (int g = 0; g < LANES; g++) x[g] = msmall(L[g].c->one, L[g].pw[(L[g].e >> sh) & mask], L[g].c);
    for (sh -= w; sh >= 0; sh -= w) {
        for (int s = 0; s < w; s++)
            for (int g = 0; g < LANES; g++) x[g] = msqr(x[g], L[g].c);
        for (int g = 0; g < LANES; g++) x[g] = msmall(x[g], L[g].pw[(L[g].e >> sh) & mask], L[g].c);
    }
    for (int g = 0; g < LANES; g++) r[g] = redc(x[g], L[g].c);
}

/* Fermat quotient via the fast kernel (single lane, padded) */
static u64 fermat_quotient(u64 p, int j, u128 *res)
{
    pctx_t c;
    pctx_init(&c, p);
    lane_t L[LANES];
    for (int g = 0; g < LANES; g++) { L[g].c = &c; L[g].pw = B.pw[j]; L[g].e = p - 1; L[g].j = j; }
    u128 r[LANES];
    run_group(L, B.w[j], r);
    if (res) *res = r[0];
    return ((u64)r[0] - 1) * c.pinv;
}

/* ------------------------------------------------------------------ */
/* Scan state                                                           */
/* ------------------------------------------------------------------ */

typedef struct {
    bool done;
    u64 primes, tests, near, flt;
    u64 cks[MAXB];
} slot_t;

typedef struct { u32 b; u64 p; u64 chunk; bool announced; } hit_t;

static struct {
    u64 start, chunk, end, nchunks;
    atomic_ullong next;
    u64 frontier;                   /* chunks [0, frontier) are complete */
    u64 primes, tests, near, flt;   /* totals over [start, frontier position) */
    u64 cks[MAXB];
    slot_t slot[RINGC];
    int nhits;
    hit_t hits[MAXHITS];
    atomic_ullong tests_live;       /* tests done in all chunks, for the rate display */
} S;

static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static long long g_near = 0;        /* report |A| <= g_near */
static bool g_print = true;         /* print FOUND/NEAR lines (off in selftest) */
static FILE *g_log = NULL;          /* chunk log */

static void record_hit(u32 b, u64 p, u64 idx)
{
    pthread_mutex_lock(&g_mu);
    bool dup = false;
    for (int i = 0; i < S.nhits; i++) if (S.hits[i].b == b && S.hits[i].p == p) dup = true;
    if (!dup && S.nhits < MAXHITS) {
        hit_t *h = &S.hits[S.nhits++];
        h->b = b; h->p = p; h->chunk = idx; h->announced = false;
    }
    if (g_print && !dup) {
        if (stderr_tty) fputs("\r\033[K", stderr);
        const char *a = anum_for(b);
        printf("FOUND base %u p %" PRIu64 "%s%s%s\n", b, p, a ? " (" : "", a ? a : "", a ? ")" : "");
        fflush(stdout);
    }
    pthread_mutex_unlock(&g_mu);
}

static void report_line(const char *fmt, ...)
{
    va_list ap;
    pthread_mutex_lock(&g_mu);
    if (stderr_tty) fputs("\r\033[K", stderr);
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
    pthread_mutex_unlock(&g_mu);
}

typedef struct {
    u64 primes, tests, near, flt;
    u64 cks[MAXB];
    u64 idx;
} acc_t;

static inline void finish_lanes(const lane_t *L, const u128 *r, int n, acc_t *a)
{
    for (int g = 0; g < n; g++) {
        const u64 p = L[g].c->p;
        const u64 q = ((u64)r[g] - 1) * L[g].c->pinv;
        const int j = L[g].j;
        a->tests++;
        if (q >= p) {                                   /* Fermat's little theorem failed */
            a->flt++;
            report_line("FLT-ERROR base %u p %" PRIu64 " r = %016" PRIx64 "%016" PRIx64 "\n",
                        B.b[j], p, (u64)(r[g] >> 64), (u64)r[g]);
            continue;
        }
        a->cks[j] += mix64(q ^ (p * GOLD));
        if (q == 0) {
            record_hit(B.b[j], p, a->idx);
        } else if (g_near > 0) {
            long long A = q <= p / 2 ? (long long)q : -(long long)(p - q);
            if (llabs(A) <= g_near) {
                a->near++;
                if (g_print) report_line("NEAR base %u p %" PRIu64 " A %lld\n", B.b[j], p, A);
            }
        }
    }
}

static void process_batch(const pctx_t *ctx, int np, acc_t *a)
{
    lane_t L[LANES];
    u128 r[LANES];
    for (int c = 0; c < B.ncls; c++) {
        const int w = B.cls_w[c];
        int nl = 0;
        for (int i = 0; i < np; i++) {
            const u64 p = ctx[i].p;
            for (int k = 0; k < B.cls_n[c]; k++) {
                const int j = B.cls_j[c][k];
                if (p <= B.b[j] && B.b[j] % p == 0) continue;
                L[nl].c = &ctx[i]; L[nl].pw = B.pw[j]; L[nl].e = p - 1; L[nl].j = j;
                if (++nl == LANES) {
                    run_group(L, w, r);
                    finish_lanes(L, r, LANES, a);
                    nl = 0;
                }
            }
        }
        if (nl) {
            for (int g = nl; g < LANES; g++) L[g] = L[nl - 1];
            run_group(L, w, r);
            finish_lanes(L, r, nl, a);
        }
    }
    atomic_fetch_add_explicit(&S.tests_live, (u64)np * (u64)B.nb, memory_order_relaxed);
}

static void scan_chunk(primesieve_iterator *it, u64 idx)
{
    const u64 lo = S.start + idx * S.chunk;
    const u64 hi = lo + S.chunk < S.end ? lo + S.chunk : S.end;
    acc_t a;
    memset(&a, 0, sizeof a);
    a.idx = idx;
    pctx_t ctx[BATCH];
    int np = 0;
    primesieve_jump_to(it, lo, hi);
    for (;;) {
        const u64 p = primesieve_next_prime(it);
        if (p >= hi) break;
        a.primes++;
        if (p == 2) {
            for (int j = 0; j < B.nb; j++) if (B.b[j] % 4 == 1) record_hit(B.b[j], 2, idx);
            continue;
        }
        pctx_init(&ctx[np], p);
        if (++np == BATCH) {
            process_batch(ctx, np, &a);
            np = 0;
            if (g_stop) return;                 /* chunk abandoned; redone on resume */
        }
    }
    if (np) process_batch(ctx, np, &a);

    pthread_mutex_lock(&g_mu);
    slot_t *s = &S.slot[idx % RINGC];
    s->done = true;
    s->primes = a.primes; s->tests = a.tests; s->near = a.near; s->flt = a.flt;
    memcpy(s->cks, a.cks, sizeof a.cks);
    while (S.slot[S.frontier % RINGC].done) {   /* fold complete chunks in frontier order */
        slot_t *f = &S.slot[S.frontier % RINGC];
        f->done = false;
        S.primes += f->primes; S.tests += f->tests; S.near += f->near; S.flt += f->flt;
        for (int j = 0; j < B.nb; j++) S.cks[j] += f->cks[j];
        if (g_log) {
            const u64 flo = S.start + S.frontier * S.chunk;
            const u64 fhi = flo + S.chunk < S.end ? flo + S.chunk : S.end;
            fprintf(g_log, "C %" PRIu64 " %" PRIu64 " %" PRIu64, flo, fhi, f->primes);
            for (int j = 0; j < B.nb; j++) fprintf(g_log, " %u:%016" PRIx64, B.b[j], f->cks[j]);
            fputc('\n', g_log);
            fflush(g_log);
        }
        S.frontier++;
    }
    pthread_mutex_unlock(&g_mu);
}

static void *worker(void *arg)
{
    (void)arg;
    primesieve_iterator it;
    primesieve_init(&it);
    while (!g_stop) {
        u64 idx = atomic_fetch_add(&S.next, 1);
        if (idx >= S.nchunks) break;
        for (;;) {                                   /* stay within RINGC of the frontier */
            pthread_mutex_lock(&g_mu);
            u64 f = S.frontier;
            pthread_mutex_unlock(&g_mu);
            if (idx < f + RINGC - 1 || g_stop) break;
            usleep(2000);
        }
        if (g_stop) break;
        scan_chunk(&it, idx);
    }
    primesieve_free_iterator(&it);
    return NULL;
}

/* ---- checkpoint file ---- */

static void save_state(const char *path)
{
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "\nwarning: cannot write %s: %s\n", tmp, strerror(errno)); return; }
    pthread_mutex_lock(&g_mu);
    fprintf(f, "wieferich checkpoint 1\nbases ");
    print_bases(f);
    fprintf(f, "\nstart %" PRIu64 "\nchunk %" PRIu64 "\nend %" PRIu64 "\nfrontier %" PRIu64
               "\nprimes %" PRIu64 "\ntests %" PRIu64 "\nnear %" PRIu64 "\nflt %" PRIu64 "\n",
            S.start, S.chunk, S.end, S.frontier, S.primes, S.tests, S.near, S.flt);
    for (int j = 0; j < B.nb; j++) fprintf(f, "cks %u %016" PRIx64 "\n", B.b[j], S.cks[j]);
    for (int i = 0; i < S.nhits; i++)
        if (S.hits[i].chunk < S.frontier)
            fprintf(f, "hit %u %" PRIu64 " %" PRIu64 "\n", S.hits[i].b, S.hits[i].p, S.hits[i].chunk);
    pthread_mutex_unlock(&g_mu);
    fputs("end\n", f);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        fprintf(stderr, "\nwarning: cannot update %s: %s\n", path, strerror(errno));
}

/* Returns false if there is no checkpoint; dies on a mismatch. */
static bool load_state(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return false;
        die("cannot read %s: %s", path, strerror(errno));
    }
    char line[4096], bases[4096] = "";
    if (!fgets(line, sizeof line, f) || strncmp(line, "wieferich checkpoint 1", 22) != 0)
        die("%s is not a wieferich checkpoint", path);
    u64 start = 0, chunk = 0, end = 0, frontier = 0;
    bool complete = false;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "bases %4000s", bases) == 1) continue;
        if (sscanf(line, "start %" SCNu64, &start) == 1) continue;
        if (sscanf(line, "chunk %" SCNu64, &chunk) == 1) continue;
        if (sscanf(line, "end %" SCNu64, &end) == 1) continue;
        if (sscanf(line, "frontier %" SCNu64, &frontier) == 1) continue;
        if (sscanf(line, "primes %" SCNu64, &S.primes) == 1) continue;
        if (sscanf(line, "tests %" SCNu64, &S.tests) == 1) continue;
        if (sscanf(line, "near %" SCNu64, &S.near) == 1) continue;
        if (sscanf(line, "flt %" SCNu64, &S.flt) == 1) continue;
        if (!strncmp(line, "cks ", 4)) {
            u32 b; u64 v;
            if (sscanf(line, "cks %u %" SCNx64, &b, &v) != 2) die("%s: bad cks line", path);
            int j = -1;
            for (int k = 0; k < B.nb; k++) if (B.b[k] == b) j = k;
            if (j < 0) die("%s: base %u is not in the current base list", path, b);
            S.cks[j] = v;
            continue;
        }
        if (!strncmp(line, "hit ", 4)) {
            hit_t h = {0, 0, 0, false};
            if (sscanf(line, "hit %u %" SCNu64 " %" SCNu64, &h.b, &h.p, &h.chunk) != 3) die("%s: bad hit line", path);
            if (S.nhits < MAXHITS) S.hits[S.nhits++] = h;
            continue;
        }
        if (!strncmp(line, "end", 3)) { complete = true; break; }
    }
    fclose(f);
    if (!complete) die("%s is truncated", path);
    char cur[4096];
    FILE *m = fmemopen(cur, sizeof cur, "w");
    print_bases(m);
    fclose(m);
    if (strcmp(cur, bases) != 0) die("%s was written for bases %s, not %s", path, bases, cur);
    if (start != S.start || chunk != S.chunk)
        die("%s was written for start=%" PRIu64 " chunk=%" PRIu64 ", not start=%" PRIu64 " chunk=%" PRIu64,
            path, start, chunk, S.start, S.chunk);
    if (frontier > S.nchunks) die("%s is past the new END (frontier %" PRIu64 " > %" PRIu64 " chunks)",
                                  path, frontier, S.nchunks);
    (void)end;
    S.frontier = frontier;
    return true;
}

static void on_sigint(int sig)
{
    (void)sig;
    g_stop = 1;
}

/* Announce hits below the frontier (confirmed); returns their number. */
static int announce_locked(void)
{
    int n = 0;
    for (int i = 0; i < S.nhits; i++) {
        hit_t *h = &S.hits[i];
        if (h->chunk >= S.frontier) continue;
        n++;
        if (!h->announced) {
            h->announced = true;
            if (g_print) {
                if (stderr_tty) fputs("\r\033[K", stderr);
                const char *a = anum_for(h->b);
                printf("CONFIRMED base %u p %" PRIu64 "%s%s%s (all primes in [%" PRIu64 ", %" PRIu64 ") tested)\n",
                       h->b, h->p, a ? " (" : "", a ? a : "", a ? ")" : "",
                       S.start, S.start + S.frontier * S.chunk);
                fflush(stdout);
            }
        }
    }
    return n;
}

typedef struct {
    double seconds;
    u64 primes, tests, numbers;
} scan_stats_t;

static scan_stats_t run_scan(u64 start, u64 end, u64 chunk, int threads, const char *state,
                             int interval, bool quiet)
{
    if (end <= start) die("END must be larger than START");
    if (chunk == 0) die("CHUNK must be positive");
    S.start = start;
    S.chunk = chunk;
    S.nchunks = (end - start + chunk - 1) / chunk;
    if (S.nchunks == 0 || (S.nchunks - 1) > (UINT64_MAX - start) / chunk) die("range too large");
    S.end = start + S.nchunks * chunk;
    if (S.end > PMAX || S.end < start) die("END must be <= 2^63 (after rounding up to the chunk size)");
    if (S.end != end && !quiet)
        fprintf(stderr, "note: END rounded up to the chunk boundary %" PRIu64 "\n", S.end);
    S.frontier = 0;
    S.primes = S.tests = S.near = S.flt = 0;
    memset(S.cks, 0, sizeof S.cks);
    memset(S.slot, 0, sizeof S.slot);
    S.nhits = 0;
    atomic_store(&S.tests_live, 0);

    if (state && load_state(state) && !quiet)
        fprintf(stderr, "resuming from %s: %" PRIu64 " of %" PRIu64 " chunks done (position %" PRIu64
                        ", %" PRIu64 " primes tested)\n",
                state, S.frontier, S.nchunks, S.start + S.frontier * S.chunk, S.primes);
    const u64 frontier0 = S.frontier;
    atomic_store(&S.next, S.frontier);
    if (!quiet) {
        fprintf(stderr, "scanning [%" PRIu64 ", %" PRIu64 ") in %" PRIu64 " chunks of %" PRIu64
                        " with %d threads, bases ", S.start, S.end, S.nchunks, S.chunk, threads);
        print_bases(stderr);
        fputc('\n', stderr);
    }

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    pthread_t *th = calloc((size_t)threads, sizeof *th);
    if (!th) die("out of memory");
    for (int i = 0; i < threads; i++)
        if (pthread_create(&th[i], NULL, worker, NULL) != 0) die("pthread_create failed");

    enum { NS = 128 };      /* rate window: one sample per second, ~2 minutes deep */
    double ts[NS];
    u64 ps[NS], qs[NS];
    int ns = 1;
    const double t0 = now();
    const u64 pos0 = S.start + frontier0 * S.chunk;
    double last_status = 0, last_ckpt = t0, last_sample = t0;
    ts[0] = t0; ps[0] = pos0; qs[0] = 0;
    char b1[32], b2[32];

    for (;;) {
        usleep(200000);
        double t = now();
        pthread_mutex_lock(&g_mu);
        u64 f = S.frontier, primes = S.primes, flt = S.flt, near = S.near;
        int nconf = announce_locked(), nhits = S.nhits;
        pthread_mutex_unlock(&g_mu);
        u64 pos = S.start + f * S.chunk;
        if (pos > S.end) pos = S.end;
        if (f >= S.nchunks || g_stop) break;
        u64 tl = atomic_load(&S.tests_live);
        if (t - last_sample >= 1.0) {
            last_sample = t;
            if (ns < NS) { ts[ns] = t; ps[ns] = pos; qs[ns] = tl; ns++; }
            else {
                memmove(ts, ts + 1, (NS - 1) * sizeof *ts);
                memmove(ps, ps + 1, (NS - 1) * sizeof *ps);
                memmove(qs, qs + 1, (NS - 1) * sizeof *qs);
                ts[NS - 1] = t; ps[NS - 1] = pos; qs[NS - 1] = tl;
            }
        }
        double dt = t - ts[0];
        double rate = dt > 0.5 ? (double)(pos - ps[0]) / dt : 0;
        double trate = dt > 0.5 ? (double)(tl - qs[0]) / dt : 0;
        double eta = rate > 0 ? (double)(S.end - pos) / rate : -1;
        bool show = stderr_tty ? (t - last_status >= 1.0) : (t - last_status >= interval);
        if (show && !quiet) {
            last_status = t;
            fprintf(stderr, "%s%s  pos %.6g (%.2f%%)  %.3g/s  %.3g tests/s  %.4g primes  ETA %s  found %d (%d confirmed)",
                    stderr_tty ? "\r\033[K" : "progress: ", fmt_dur(t - t0, b1, sizeof b1), (double)pos,
                    100.0 * (double)(pos - S.start) / (double)(S.end - S.start), rate, trate,
                    (double)primes, fmt_dur(eta, b2, sizeof b2), nhits, nconf);
            if (g_near > 0) fprintf(stderr, "  near %" PRIu64, near);
            if (flt) fprintf(stderr, "  FLT-ERRORS %" PRIu64, flt);
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) {
            last_ckpt = t;
            save_state(state);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);
    pthread_mutex_lock(&g_mu);
    announce_locked();
    pthread_mutex_unlock(&g_mu);
    if (state) save_state(state);
    if (g_stop && S.frontier < S.nchunks && !quiet)
        fprintf(stderr, "stopped at position %" PRIu64 " (%" PRIu64 " of %" PRIu64 " chunks complete)%s\n",
                S.start + S.frontier * S.chunk, S.frontier, S.nchunks, state ? ", checkpoint saved" : "");
    free(th);
    scan_stats_t st;
    st.seconds = t1 - t0;
    st.primes = S.primes;
    st.tests = S.tests;
    st.numbers = (S.frontier - frontier0) * S.chunk;
    return st;
}

static void print_summary(void)
{
    u64 pos = S.start + S.frontier * S.chunk;
    if (pos > S.end) pos = S.end;
    printf("SUMMARY [%" PRIu64 ", %" PRIu64 ") primes %" PRIu64 " tests %" PRIu64 " flt_errors %" PRIu64,
           S.start, pos, S.primes, S.tests, S.flt);
    if (g_near > 0) printf(" near(|A|<=%lld) %" PRIu64, g_near, S.near);
    putchar('\n');
    for (int j = 0; j < B.nb; j++) {
        const char *a = anum_for(B.b[j]);
        printf("base %3u %-8s cks %016" PRIx64 "  solutions:", B.b[j], a ? a : "", S.cks[j]);
        int n = 0;
        for (int i = 0; i < S.nhits; i++)
            if (S.hits[i].b == B.b[j] && S.hits[i].chunk < S.frontier) { printf(" %" PRIu64, S.hits[i].p); n++; }
        if (!n) printf(" none");
        putchar('\n');
    }
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* Commands                                                             */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fputs("usage: wieferich scan [START] END [-b BASES] [-t T] [-c CHUNK] [-n T] [-S FILE] [-i SECS] [-L FILE] [-q]\n"
          "       wieferich check P [-b BASES]\n"
          "       wieferich bench [N] [-d SPAN] [-b BASES] [-t T]\n"
          "       wieferich selftest [LIMIT] [-t T]\n"
          "BASES: list/ranges like 3,5,7 or 2-30; 'oeis' (default) = 3,5,6,7,10,12,13,14,15,17,18,19,20,22,23,26,30;\n"
          "       'known' = all non-power bases 2..30\n", stderr);
    exit(2);
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2];
    int npos = 0, threads = default_threads(), interval = 60;
    u64 chunk = 1000000000ULL;
    const char *state = NULL, *logf = NULL, *bases = "oeis";
    bool quiet = false;
    for (int i = 0; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(a, "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
        else if (!strcmp(a, "-b") && i + 1 < argc) bases = argv[++i];
        else if (!strcmp(a, "-n") && i + 1 < argc) g_near = (long long)parse_num(argv[++i]);
        else if (!strcmp(a, "-S") && i + 1 < argc) state = argv[++i];
        else if (!strcmp(a, "-L") && i + 1 < argc) logf = argv[++i];
        else if (!strcmp(a, "-i") && i + 1 < argc) interval = parse_int(argv[++i], 1, 86400, "interval");
        else if (!strcmp(a, "-q")) quiet = true;
        else if (a[0] != '-' && npos < 2) pos[npos++] = parse_num(a);
        else usage();
    }
    if (npos == 0) usage();
    u64 start = npos == 2 ? pos[0] : 0, end = npos == 2 ? pos[1] : pos[0];
    parse_bases(bases);
    if (logf && !(g_log = fopen(logf, "a"))) die("cannot open %s: %s", logf, strerror(errno));
    scan_stats_t st = run_scan(start, end, chunk, threads, state, interval, quiet);
    print_summary();
    if (!quiet) {
        char b1[32];
        fprintf(stderr, "%s, %.3g primes/s, %.3g tests/s\n", fmt_dur(st.seconds, b1, sizeof b1),
                st.seconds > 0 ? (double)(S.primes) / st.seconds : 0, st.seconds > 0 ? (double)S.tests / st.seconds : 0);
    }
    if (g_log) fclose(g_log);
    return S.flt ? 3 : 0;
}

static int cmd_check(int argc, char **argv)
{
    u64 p = 0;
    const char *bases = "known";
    bool havep = false;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-b") && i + 1 < argc) bases = argv[++i];
        else if (argv[i][0] != '-' && !havep) { p = parse_num(argv[i]); havep = true; }
        else usage();
    }
    if (!havep) usage();
    if (p < 3 || p >= PMAX) die("P must be in [3, 2^63)");
    parse_bases(bases);
    printf("p = %" PRIu64 " (%s)\n", p, is_prime64(p) ? "prime" : "NOT prime: quotients are meaningless");
    int bad = 0;
    for (int j = 0; j < B.nb; j++) {
        u32 b = B.b[j];
        if (b % p == 0) { printf("base %u: p divides b\n", b); continue; }
        u128 rf;
        u64 q = fermat_quotient(p, j, &rf);
        u128 rr = powmod_ref(b, p - 1, (u128)p * p);
        u64 qr = (u64)((rr - 1) / p);
        long long A = q <= p / 2 ? (long long)q : -(long long)(p - q);
        bool ok = rf == rr && q == qr;
        if (!ok) bad++;
        printf("base %3u: q = %" PRIu64 "  A = %lld%s  [reference %s]\n", b, q, A,
               q == 0 ? "  <-- WIEFERICH" : "", ok ? "agrees" : "DISAGREES");
    }
    return bad ? 1 : 0;
}

/* integral of dt / ln t over [x, y] */
static double li_diff(double x, double y)
{
    const int n = 2000;
    double lx = log(x), ly = log(y), h = (ly - lx) / n, s = 0;
    for (int i = 0; i < n; i++) {
        double u = lx + (i + 0.5) * h;
        s += exp(u) / u * h;
    }
    return s;
}

static int cmd_bench(int argc, char **argv)
{
    u64 n0 = 1000000000000000ULL, span = 10000000000ULL;
    int threads = default_threads();
    const char *bases = "oeis";
    bool haven = false;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-d") && i + 1 < argc) span = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(argv[i], "-b") && i + 1 < argc) bases = argv[++i];
        else if (argv[i][0] != '-' && !haven) { n0 = parse_num(argv[i]); haven = true; }
        else usage();
    }
    parse_bases(bases);
    g_print = true;
    u64 chunk = span / (u64)(threads * 8);
    if (chunk < 1000000) chunk = 1000000;
    scan_stats_t st = run_scan(n0, n0 + span, chunk, threads, NULL, 60, true);
    double tps = (double)st.tests / st.seconds, pps = (double)st.primes / st.seconds;
    printf("bench [%" PRIu64 ", +%" PRIu64 "): %" PRIu64 " primes x %d bases in %.2f s with %d threads\n",
           n0, span, st.primes, B.nb, st.seconds, threads);
    printf("  %.4g primes/s, %.4g tests/s (%.4g tests/s per thread), %.4g numbers/s\n",
           pps, tps, tps / threads, (double)st.numbers / st.seconds);
    const double targets[] = {2, 10, 100};
    for (int i = 0; i < 3; i++) {
        double x = (double)n0, y = x * targets[i];
        double tests = li_diff(x, y) * B.nb;
        double exp_new = B.nb * log(log(y) / log(x));
        char b1[32];
        printf("  [N, %gN): %.3g tests, ~%s at this rate; heuristic %.3f new solutions over %d bases\n",
               targets[i], tests, fmt_dur(tests / tps, b1, sizeof b1), exp_new, B.nb);
    }
    return 0;
}

/* ---- selftest ---- */

static u64 rng_state = 0x243F6A8885A308D3ULL;
static u64 rnd(void) { return mix64(rng_state += GOLD); }

static int test_arith(void)
{
    int bad = 0;
    for (int it = 0; it < 200000 && bad < 5; it++) {
        int bits = 2 + (int)(rnd() % 62);                 /* p in [2^(bits-1), 2^bits) */
        u64 p = (rnd() >> (64 - bits)) | (1ULL << (bits - 1)) | 1;
        if (p >= PMAX || p < 3) continue;
        pctx_t c;
        pctx_init(&c, p);
        if (p * c.pinv != 1 || (u64)c.N * c.ninv != UINT64_MAX) { bad++; fprintf(stderr, "inverse p=%" PRIu64 "\n", p); continue; }
        if (c.one != ((u128)0 - c.N) % c.N) { bad++; fprintf(stderr, "one p=%" PRIu64 "\n", p); continue; }
        u128 x = (((u128)rnd() << 64) | rnd()) % c.N;
        u64 k = rnd() >> (rnd() % 64);
        u64 cst = k & 0xFFFFFFFFULL;
        if (msmall(x, cst, &c) != mulmod_ref(x, cst, c.N)) { bad++; fprintf(stderr, "msmall p=%" PRIu64 "\n", p); }
        if (mulmod_ref(msqr(x, &c), c.one, c.N) != mulmod_ref(x, x, c.N)) { bad++; fprintf(stderr, "msqr p=%" PRIu64 "\n", p); }
        if (mulmod_ref(redc(x, &c), c.one, c.N) != x) { bad++; fprintf(stderr, "redc p=%" PRIu64 "\n", p); }
    }
    /* worst cases: x = N-1, largest constants, largest p */
    const u64 edge[] = {3, 5, 4294967291ULL, 4294967311ULL, 18446744073709551557ULL >> 1 | 1, PMAX - 25};
    for (size_t i = 0; i < sizeof edge / sizeof *edge; i++) {
        u64 p = edge[i];
        pctx_t c;
        pctx_init(&c, p);
        u128 x = c.N - 1;
        if (msmall(x, 0xFFFFFFFFULL, &c) != mulmod_ref(x, 0xFFFFFFFFULL, c.N)) { bad++; fprintf(stderr, "edge msmall %" PRIu64 "\n", p); }
        if (mulmod_ref(msqr(x, &c), c.one, c.N) != mulmod_ref(x, x, c.N)) { bad++; fprintf(stderr, "edge msqr %" PRIu64 "\n", p); }
    }
    /* full Fermat quotients for random primes and bases */
    u32 bl[MAXB];
    int nbl = 0;
    for (u32 b = 2; nbl < MAXB; b += 1 + (u32)(rnd() % 97)) bl[nbl++] = b;
    bl[nbl - 1] = 4294967291U;                          /* a base with w = 1 */
    bases_setup(bl, nbl);
    int nq = 0;
    for (int it = 0; it < 3000 && bad < 5; it++) {
        int bits = 3 + (int)(rnd() % 61);
        u64 p = (rnd() >> (64 - bits)) | (1ULL << (bits - 1)) | 1;
        if (p >= PMAX) continue;
        while (!is_prime64(p)) p += 2;
        if (p >= PMAX) continue;
        int j = (int)(rnd() % (u64)B.nb);
        if (B.b[j] % p == 0) continue;
        u128 rf;
        u64 q = fermat_quotient(p, j, &rf);
        u128 rr = powmod_ref(B.b[j], p - 1, (u128)p * p);
        if (rf != rr || q != (u64)((rr - 1) / p)) {
            bad++;
            fprintf(stderr, "quotient mismatch p=%" PRIu64 " b=%u\n", p, B.b[j]);
        }
        nq++;
    }
    fprintf(stderr, "arithmetic: 200000 random reductions + edge cases, %d Fermat quotients vs reference: %s\n",
            nq, bad ? "FAIL" : "ok");
    return bad;
}

static bool known_hit(u32 b, u64 p)
{
    for (int i = 0; i < NKNOWN; i++)
        if (KNOWN[i].b == b)
            for (int k = 0; k < KNOWN[i].n; k++) if (KNOWN[i].p[k] == p) return true;
    return false;
}

/* Checksums from verify_wieferich.py (Python integers, Miller-Rabin primes), bases 'oeis'. */
typedef struct { u64 lo, hi, primes; u64 cks[17]; } ckwin_t;
static const ckwin_t CKWIN[] = {
#include "wieferich_cks.h"
};

static int cmd_selftest(int argc, char **argv)
{
    u64 limit = 1000000000ULL;
    int threads = default_threads();
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (argv[i][0] != '-') limit = parse_num(argv[i]);
        else usage();
    }
    int bad = test_arith();
    g_print = false;

    /* 1. all bases 2..30 over [0, limit) against the known solutions */
    parse_bases("known");
    double t = now();
    u64 chunk = limit / (u64)(threads * 16);
    if (chunk < 1000000) chunk = 1000000;
    run_scan(0, limit, chunk, threads, NULL, 60, true);
    int found = 0, missing = 0, extra = 0;
    for (int i = 0; i < S.nhits; i++) {
        if (known_hit(S.hits[i].b, S.hits[i].p)) found++;
        else { extra++; fprintf(stderr, "  unexpected solution base %u p %" PRIu64 "\n", S.hits[i].b, S.hits[i].p); }
    }
    for (int i = 0; i < NKNOWN; i++)
        for (int k = 0; k < KNOWN[i].n; k++) {
            if (KNOWN[i].p[k] >= S.end) continue;
            bool f = false;
            for (int h = 0; h < S.nhits; h++) if (S.hits[h].b == KNOWN[i].b && S.hits[h].p == KNOWN[i].p[k]) f = true;
            if (!f) { missing++; fprintf(stderr, "  missing base %u p %" PRIu64 "\n", KNOWN[i].b, KNOWN[i].p[k]); }
        }
    bool pi_ok = limit != 1000000000ULL || S.primes == 50847534ULL;
    fprintf(stderr, "scan [0, %" PRIu64 "), %d bases: %" PRIu64 " primes%s, %d known solutions found, %d missing, %d unexpected, %" PRIu64 " FLT errors (%.1f s)\n",
            S.end, B.nb, S.primes, limit == 1000000000ULL ? (pi_ok ? " (= pi(10^9))" : " (WRONG, pi(10^9) = 50847534)") : "",
            found, missing, extra, S.flt, now() - t);
    bad += missing + extra + (int)S.flt + !pi_ok;

    /* 2. windows around the known solutions above the scan */
    t = now();
    const u64 scanned = S.end;
    int wins = 0, wbad = 0;
    for (int i = 0; i < NKNOWN; i++)
        for (int k = 0; k < KNOWN[i].n; k++) {
            u64 p = KNOWN[i].p[k];
            if (p < scanned) continue;
            u32 b = KNOWN[i].b;
            bases_setup(&b, 1);
            run_scan(p > 1000000 ? p - 1000000 : 0, p + 1000000, 200000, threads, NULL, 60, true);
            bool ok = S.nhits == 1 && S.hits[0].p == p && S.flt == 0;
            if (!ok) { wbad++; fprintf(stderr, "  window around base %u p %" PRIu64 ": %d hits\n", b, p, S.nhits); }
            wins++;
        }
    fprintf(stderr, "windows: %d known solutions above %" PRIu64 " re-found through the scan: %s (%.1f s)\n",
            wins, limit, wbad ? "FAIL" : "ok", now() - t);
    bad += wbad;

    /* 3. checksums against verify_wieferich.py */
    t = now();
    parse_bases("oeis");
    int cbad = 0;
    for (size_t w = 0; w < sizeof CKWIN / sizeof *CKWIN; w++) {
        const ckwin_t *c = &CKWIN[w];
        run_scan(c->lo, c->hi, (c->hi - c->lo) / 16, threads, NULL, 60, true);
        bool ok = S.primes == c->primes && S.flt == 0 && S.end == c->hi;
        for (int j = 0; j < B.nb; j++) if (S.cks[j] != c->cks[j]) ok = false;
        if (!ok) {
            cbad++;
            fprintf(stderr, "  checksum window [%" PRIu64 ", %" PRIu64 "): primes %" PRIu64 " (expected %" PRIu64 ")\n",
                    c->lo, c->hi, S.primes, c->primes);
            for (int j = 0; j < B.nb; j++)
                if (S.cks[j] != c->cks[j]) fprintf(stderr, "    base %u: %016" PRIx64 " expected %016" PRIx64 "\n", B.b[j], S.cks[j], c->cks[j]);
        }
    }
    fprintf(stderr, "checksums: %zu windows x %d bases vs verify_wieferich.py: %s (%.1f s)\n",
            sizeof CKWIN / sizeof *CKWIN, B.nb, cbad ? "FAIL" : "ok", now() - t);
    bad += cbad;
    fprintf(stderr, "selftest %s\n", bad ? "FAILED" : "passed");
    return bad ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    if (argc < 2) usage();
    const char *cmd = argv[1];
    if (!strcmp(cmd, "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(cmd, "check")) return cmd_check(argc - 2, argv + 2);
    if (!strcmp(cmd, "bench")) return cmd_bench(argc - 2, argv + 2);
    if (!strcmp(cmd, "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
