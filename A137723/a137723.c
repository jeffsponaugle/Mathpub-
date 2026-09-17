/*
 * a137723.c
 *
 * Compute and extend OEIS A137723:
 *
 *   "First occurrence of a set of n consecutive numbers having at least one
 *    prime gap in their factorization: a(n) = smallest number of this set."
 *
 * Known terms (n = 1..31):
 *   10, 33, 20, 55, 84, 114, 390, 513, 182, 200, 468, 2941, 774, 65522, 1832,
 *   1261, 1130, 1332, 1638, 524289, 1952, 4298, 4524, 69960, 5120, 16385, 2972,
 *   4832, 5352, 10801, 5592;   a(32) > 10^11 (Lucas A. Brown, Oct 2024).
 *
 * Definitions
 * -----------
 * A number has a "prime gap in its factorization" (A073490(m) > 0, called
 * "gapful" below) when its distinct prime factors are not consecutive primes:
 * 10 = 2*5 skips 3, 84 = 2^2*3*7 skips 5.  Numbers with no gap (A073491,
 * "gapless") are 1, the primes, the prime powers, and products of consecutive
 * primes with arbitrary exponents: 12 = 2^2*3, 35 = 5*7, 2940 = 2^2*3*5*7^2.
 *
 * If g < g' are consecutive gapless numbers, then g+1 .. g'-1 is a maximal run
 * of n = g'-g-1 gapful numbers, and a(n) is the smallest g+1 over all runs of
 * length exactly n.  Example: 83 and 89 are consecutive gapless numbers and
 * 84..88 are all gapful, so a(5) = 84.
 *
 * Methods
 * -------
 * scan N      Exhaustive.  Gapless numbers are sparse (about N/ln N primes
 *             plus O(sqrt N) composites), so instead of factoring every
 *             integer we generate them directly: primes with primesieve,
 *             composite gapless numbers by a recursive enumeration of
 *             consecutive-prime products.  Worker threads take chunks of
 *             [1, N] from an atomic counter, merge the two sorted streams and
 *             record every difference.  Every a(n) <= N is found exactly, and
 *             every n not found has a(n) > N.
 *
 * search n    Structural search for EVEN n (odd d = n+1).  The bounding
 *             gapless numbers L and L+d then have different parity, so exactly
 *             one of them is even, and an even gapless number is necessarily
 *             E = 2^e1 3^e2 ... p_k^ek (the first k primes, all exponents >= 1).
 *             There are only ~10^6 such E below 2^64 and ~10^9 below 2^127, so
 *             we enumerate them all, test E+d and E-d for gaplessness, and
 *             verify that the n numbers in between are gapful.  The search is
 *             exhaustive for L <= limit-d, so the smallest hit IS a(n).  All
 *             arithmetic is 128-bit, so the search can go far beyond 2^64.
 *
 *             Odd primes shared by d and E constrain the odd partner O = E+-d:
 *             if p | E and p | d then p | O, and if p | E but p does not divide
 *             d then p does not divide O.  Since O's primes must be
 *             consecutive, the odd primes <= p_k dividing d must form a
 *             contiguous block of primes (otherwise no E with k primes can
 *             work: this bounds k by "kmax"), and O can be prime only when no
 *             odd prime <= p_k divides d.  For d = 33 = 3*11 this forces
 *             E = 2^a whenever O is prime, which is why a(32) is so large:
 *             only 36 powers of two lie below 10^11.  The same happens for
 *             every n == 2 (mod 6).
 *
 * Usage
 * -----
 *   a137723 scan N [-t T] [-m MAXN] [-o FILE]
 *   a137723 search n [n ...] [-l LIMIT] [-t T]      default LIMIT 2^64; above 2^127 uses GMP
 *   a137723 next [-N SCANLIMIT] [-l LIMIT] [-t T]   scan, then search missing even n
 *   a137723 check m                                  factor m, count its prime gaps
 *   a137723 run m                                    show the maximal run around m
 *   a137723 selftest [-t T]
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e15, or 2^54-32.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a137723.c \
 *            -L/opt/homebrew/lib -lprimesieve -o a137723
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
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

#define SMALL_LIMIT 65536u   /* trial-division bound (table of primes below it) */
#define MAXRUN      65536    /* longest run we can record; prime gaps < 2^64 are < 1600 */
#define NKNOWN      31
#define MAXPE       48       /* max distinct primes in a 128-bit number is 26 */

static bool stderr_tty;   /* live \r progress only on a terminal */

static const u64 KNOWN[NKNOWN + 1] = { 0,
    10, 33, 20, 55, 84, 114, 390, 513, 182, 200, 468, 2941, 774, 65522, 1832,
    1261, 1130, 1332, 1638, 524289, 1952, 4298, 4524, 69960, 5120, 16385, 2972,
    4832, 5352, 10801, 5592 };

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
    fputs("a137723: ", stderr);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) die("out of memory");
    return p;
}

static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz);
    if (!p) die("out of memory");
    return p;
}

/* decimal string of a 128-bit value; buf must hold >= 40 chars */
static char *u128_str(u128 v, char *buf)
{
    char tmp[48];
    int i = 0, j = 0;
    if (v == 0) tmp[i++] = '0';
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    while (i) buf[j++] = tmp[--i];
    buf[j] = 0;
    return buf;
}

static int bitlen128(u128 v)
{
    u64 hi = (u64)(v >> 64), lo = (u64)v;
    if (hi) return 128 - __builtin_clzll(hi);
    return lo ? 64 - __builtin_clzll(lo) : 0;
}

static bool parse_dec(const char *s, size_t len, u128 *out)
{
    u128 r = 0;
    if (len == 0) return false;
    for (size_t i = 0; i < len; i++) {
        if (s[i] < '0' || s[i] > '9') return false;
        unsigned dgt = (unsigned)(s[i] - '0');
        if (r > (((u128)0 - 1) - dgt) / 10) return false;
        r = r * 10 + dgt;
    }
    *out = r;
    return true;
}

/* accepts "123", "2^64", "10^12", "1e15", and X+Y / X-Y of those ("2^54-32") */
static bool parse_u128(const char *s, u128 *out)
{
    const char *pm = strpbrk(s, "+-");
    if (pm && pm != s) {
        char head[128];
        size_t hl = (size_t)(pm - s);
        u128 a, b;
        if (hl >= sizeof head) return false;
        memcpy(head, s, hl);
        head[hl] = 0;
        if (!parse_u128(head, &a) || !parse_u128(pm + 1, &b)) return false;
        if (*pm == '+') { if (a > ((u128)0 - 1) - b) return false; *out = a + b; }
        else            { if (b > a) return false; *out = a - b; }
        return true;
    }
    const char *caret = strchr(s, '^');
    const char *ee = strpbrk(s, "eE");
    if (caret) {
        u128 b, e, r = 1;
        if (!parse_dec(s, (size_t)(caret - s), &b)) return false;
        if (!parse_dec(caret + 1, strlen(caret + 1), &e) || e > 400) return false;
        for (unsigned i = 0; i < (unsigned)e; i++) {
            if (b != 0 && r > ((u128)0 - 1) / b) return false;
            r *= b;
        }
        *out = r;
        return true;
    }
    if (ee) {
        u128 m, e;
        if (!parse_dec(s, (size_t)(ee - s), &m)) return false;
        if (!parse_dec(ee + 1, strlen(ee + 1), &e) || e > 40) return false;
        for (unsigned i = 0; i < (unsigned)e; i++) {
            if (m > ((u128)0 - 1) / 10) return false;
            m *= 10;
        }
        *out = m;
        return true;
    }
    return parse_dec(s, strlen(s), out);
}

static u128 arg_u128(const char *s, const char *what)
{
    u128 v;
    if (!parse_u128(s, &v)) die("bad %s: '%s'", what, s);
    return v;
}

static u64 arg_u64(const char *s, const char *what)
{
    u128 v = arg_u128(s, what);
    if (v > UINT64_MAX) die("%s too large (max 2^64-1): '%s'", what, s);
    return (u64)v;
}

static int cmp_u64(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

/* ------------------------------------------------------------------ */
/* Small primes (below SMALL_LIMIT)                                    */
/* ------------------------------------------------------------------ */

static u32     *sp;          /* primes < SMALL_LIMIT                       */
static size_t   nsp;
static uint8_t *small_isp;   /* small_isp[x] = 1 iff x is prime, x < 65536 */
static u32     *r64tab;      /* r64tab[i] = 2^64 mod sp[i]                 */

static void init_small_primes(void)
{
    size_t n = 0;
    u32 *p = (u32 *)primesieve_generate_primes(2, SMALL_LIMIT - 1, &n, UINT32_PRIMES);
    if (!p || n == 0) die("primesieve_generate_primes failed");
    sp = p;
    nsp = n;
    small_isp = xcalloc(SMALL_LIMIT, 1);
    r64tab = xmalloc(n * sizeof(u32));
    for (size_t i = 0; i < n; i++) {
        small_isp[p[i]] = 1;
        r64tab[i] = (u32)((((u128)1) << 64) % p[i]);
    }
}

/* m mod sp[i], using only 64-bit divisions */
static inline u32 mod_small(u128 m, size_t i)
{
    u32 p = sp[i];
    if (m <= UINT64_MAX) return (u32)((u64)m % p);
    u64 hi = (u64)(m >> 64), lo = (u64)m;
    return (u32)(((hi % p) * (u64)r64tab[i] + lo % p) % p);
}

/* ------------------------------------------------------------------ */
/* Modular arithmetic and primality                                    */
/* ------------------------------------------------------------------ */

static inline u64 mulmod64(u64 a, u64 b, u64 m)
{
    if ((a | b) < ((u64)1 << 32)) return (a * b) % m;
    return (u64)(((u128)a * b) % m);
}

static u64 powmod64(u64 a, u64 e, u64 m)
{
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod64(r, a, m);
        a = mulmod64(a, a, m);
        e >>= 1;
    }
    return r;
}

/* strong probable prime test, n odd > 2 */
static bool sprp64(u64 n, u64 a)
{
    a %= n;
    if (a == 0) return true;
    u64 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    u64 x = powmod64(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (int i = 1; i < s; i++) {
        x = mulmod64(x, x, n);
        if (x == n - 1) return true;
        if (x == 1) return false;
    }
    return false;
}

/* deterministic for all n < 2^64 */
static bool is_prime64(u64 n)
{
    static const u32 tiny[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    static const u64 b3[] = { 2, 7, 61 };
    static const u64 b7[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    if (n < 2) return false;
    for (int i = 0; i < 12; i++)
        if (n % tiny[i] == 0) return n == tiny[i];
    if (n < 41u * 41u) return true;
    if (n < 4759123141ULL) {
        for (int i = 0; i < 3; i++) if (!sprp64(n, b3[i])) return false;
        return true;
    }
    for (int i = 0; i < 7; i++) if (!sprp64(n, b7[i])) return false;
    return true;
}

static inline u128 addmod128(u128 a, u128 b, u128 m)   /* a, b < m */
{
    u128 s = a + b;
    if (s < a || s >= m) s -= m;
    return s;
}

static u128 mulmod128(u128 a, u128 b, u128 m)
{
    if (m <= UINT64_MAX) return mulmod64((u64)(a % m), (u64)(b % m), (u64)m);
    if (a >= m) a %= m;
    if (b >= m) b %= m;
    if (a < b) { u128 t = a; a = b; b = t; }
    u128 r = 0;
    while (b) {
        if (b & 1) r = addmod128(r, a, m);
        a = addmod128(a, a, m);
        b >>= 1;
    }
    return r;
}

static u128 powmod128(u128 a, u128 e, u128 m)
{
    u128 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod128(r, a, m);
        a = mulmod128(a, a, m);
        e >>= 1;
    }
    return r;
}

static bool sprp128(u128 n, u128 a)
{
    a %= n;
    if (a == 0) return true;
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    u128 x = powmod128(a, d, n);
    if (x == 1 || x == n - 1) return true;
    for (int i = 1; i < s; i++) {
        x = mulmod128(x, x, n);
        if (x == n - 1) return true;
        if (x == 1) return false;
    }
    return false;
}

/* deterministic below 2^64; above, strong probable prime to the first 24
 * prime bases (a composite passes with probability far below 4^-24)       */
static bool is_prime128(u128 n)
{
    if (n <= UINT64_MAX) return is_prime64((u64)n);
    for (size_t i = 0; i < nsp && sp[i] < 2000; i++)
        if (mod_small(n, i) == 0) return false;
    for (int i = 0; i < 24; i++)
        if (!sprp128(n, sp[i])) return false;
    return true;
}

/* smallest prime > x */
static u128 next_prime128(u128 x)
{
    if (x < 65521) {                 /* 65521 is the largest prime < 2^16 */
        u32 c = (u32)x + 1;
        while (!small_isp[c]) c++;
        return c;
    }
    u128 c = x + 1;
    if (!(c & 1)) c++;
    while (!is_prime128(c)) c += 2;
    return c;
}

/* largest prime < x, 0 if none */
static u128 prev_prime128(u128 x)
{
    if (x <= 2) return 0;
    if (x == 3) return 2;
    u128 c = x - 1;
    if (!(c & 1)) c--;
    while (!is_prime128(c)) c -= 2;
    return c;
}

/* ------------------------------------------------------------------ */
/* Integer roots                                                       */
/* ------------------------------------------------------------------ */

static u128 isqrt128(u128 n)
{
    if (n < 2) return n;
    u128 x = (u128)1 << ((bitlen128(n) + 1) / 2);   /* x >= sqrt(n) */
    for (;;) {
        u128 y = (x + n / x) >> 1;
        if (y >= x) return x;
        x = y;
    }
}

static bool pow_le(u128 x, int t, u128 n)   /* x^t <= n, overflow-safe */
{
    u128 r = 1;
    for (int i = 0; i < t; i++) {
        if (x != 0 && r > n / x) return false;
        r *= x;
    }
    return r <= n;
}

static bool pow_eq(u128 x, int t, u128 n)
{
    u128 r = 1;
    for (int i = 0; i < t; i++) {
        if (x != 0 && r > n / x) return false;
        r *= x;
    }
    return r == n;
}

static u128 iroot128(u128 n, int t)         /* floor(n^(1/t)) */
{
    if (t <= 1 || n < 2) return n;
    if (t == 2) return isqrt128(n);
    int bits = bitlen128(n);
    u128 lo = 1, hi = (u128)1 << ((bits + t - 1) / t);   /* hi^t > n */
    while (hi - lo > 1) {
        u128 mid = lo + (hi - lo) / 2;
        if (pow_le(mid, t, n)) lo = mid; else hi = mid;
    }
    return lo;
}

/* ------------------------------------------------------------------ */
/* Gapless test (A073490(m) == 0 ?)                                    */
/* ------------------------------------------------------------------ */

typedef struct { u128 p; int e; } pe_t;
typedef struct { pe_t f[MAXPE]; int nf; bool gapless; } gl_t;

static void gl_push(gl_t *g, u128 p, int e)
{
    if (g->nf < MAXPE) { g->f[g->nf].p = p; g->f[g->nf].e = e; }
    g->nf++;
}

/* Try to write r as p^e * nextprime(p)^f * ... , a complete factorization
 * into consecutive primes starting at p.  On success fills g.            */
static bool chain_from(u128 r, u128 p, gl_t *g)
{
    gl_t t;
    u128 x = r;
    t.nf = 0;
    for (;;) {
        if (x % p) return false;
        int e = 0;
        do { x /= p; e++; } while (x % p == 0);
        gl_push(&t, p, e);
        if (x == 1) break;
        p = next_prime128(p);
    }
    t.gapless = true;
    *g = t;
    return true;
}

/* r is composite and has no prime factor <= SMALL_LIMIT.  If r is gapless it
 * is p1^e1 ... ps^es with consecutive primes, t = sum(e) total factors,
 * p1 < r^(1/t) < ps, so p1 is one of the (t-1) largest primes below r^(1/t)
 * (or r = p^t).  Since every prime factor exceeds SMALL_LIMIT, r > SMALL_LIMIT^t. */
static bool root_method(u128 r, gl_t *g)
{
    for (int t = 2; t <= 8; t++) {
        u128 c = iroot128(r, t);
        if (c <= SMALL_LIMIT) break;
        bool exact = pow_eq(c, t, r);
        if (exact && is_prime128(c)) {
            g->nf = 0;
            gl_push(g, c, t);
            g->gapless = true;
            return true;
        }
        u128 q = exact ? prev_prime128(c) : (is_prime128(c) ? c : prev_prime128(c));
        for (int j = 0; j < t - 1 && q > SMALL_LIMIT; j++) {
            if (chain_from(r, q, g)) return true;
            q = prev_prime128(q);
        }
    }
    g->gapless = false;
    return false;
}

/* true iff m has no prime gap in its factorization (m is in A073491).
 * If g is given and m is gapless, g receives the full factorization.   */
static bool is_gapless(u128 m, gl_t *g)
{
    gl_t tmp;
    if (!g) g = &tmp;
    g->nf = 0;
    g->gapless = false;
    if (m == 0) return false;
    if (m == 1) { g->gapless = true; return true; }
    if (m <= UINT64_MAX) {
        u64 m64 = (u64)m;
        for (size_t i = 0; i < nsp; i++) {
            u64 p = sp[i];
            if (p * p > m64) {                 /* m is prime */
                gl_push(g, m, 1);
                g->gapless = true;
                return true;
            }
            if (m64 % p == 0) return chain_from(m, p, g);
        }
    } else {
        for (size_t i = 0; i < nsp; i++)
            if (mod_small(m, i) == 0) return chain_from(m, sp[i], g);
    }
    /* no prime factor <= SMALL_LIMIT, and m > SMALL_LIMIT^2 */
    if (is_prime128(m)) {
        gl_push(g, m, 1);
        g->gapless = true;
        return true;
    }
    return root_method(m, g);
}

/* ------------------------------------------------------------------ */
/* Full factorization (for display; Pollard-Brent rho)                 */
/* ------------------------------------------------------------------ */

typedef struct { pe_t f[MAXPE]; int nf; u128 unf; int unf_e; } fac_t;

static void fac_add(fac_t *F, u128 p, int e)
{
    for (int i = 0; i < F->nf; i++)
        if (F->f[i].p == p) { F->f[i].e += e; return; }
    if (F->nf < MAXPE) { F->f[F->nf].p = p; F->f[F->nf].e = e; F->nf++; }
}

static u64 gcd64(u64 a, u64 b)
{
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

static u128 gcd128(u128 a, u128 b)
{
    while (b) { u128 t = a % b; a = b; b = t; }
    return a;
}

static inline u64 addc64(u64 x, u64 c, u64 n)   /* (x + c) mod n, x < n, c < n */
{
    u64 s = x + c;
    if (s < x || s >= n) s -= n;
    return s;
}

/* nontrivial factor of odd composite n (not a prime power) */
static u64 rho64(u64 n)
{
    if (n % 2 == 0) return 2;
    for (u64 c = 1;; c++) {
        u64 y = 2, r = 1, q = 1, g = 1, x = 0, ys = 0;
        do {
            x = y;
            for (u64 i = 0; i < r; i++) y = addc64(mulmod64(y, y, n), c, n);
            u64 k = 0;
            do {
                ys = y;
                u64 lim = (r - k < 128) ? r - k : 128;
                for (u64 i = 0; i < lim; i++) {
                    y = addc64(mulmod64(y, y, n), c, n);
                    q = mulmod64(q, x > y ? x - y : y - x, n);
                }
                g = gcd64(q, n);
                k += lim;
            } while (k < r && g == 1);
            r <<= 1;
        } while (g == 1);
        if (g == n) {
            do {
                ys = addc64(mulmod64(ys, ys, n), c, n);
                g = gcd64(x > ys ? x - ys : ys - x, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
}

/* 128-bit variant with an iteration budget; returns 0 on failure */
static u128 rho128(u128 n, u64 budget)
{
    if (n % 2 == 0) return 2;
    u64 used = 0;
    for (u128 c = 1; used < budget; c++) {
        u128 y = 2, r = 1, q = 1, g = 1, x = 0, ys = 0;
        do {
            x = y;
            for (u128 i = 0; i < r; i++) y = addmod128(mulmod128(y, y, n), c, n);
            u128 k = 0;
            do {
                ys = y;
                u128 lim = (r - k < 128) ? r - k : 128;
                for (u128 i = 0; i < lim; i++) {
                    y = addmod128(mulmod128(y, y, n), c, n);
                    q = mulmod128(q, x > y ? x - y : y - x, n);
                }
                g = gcd128(q, n);
                k += lim;
                used += (u64)lim;
            } while (k < r && g == 1);
            r <<= 1;
        } while (g == 1 && used < budget);
        if (g == 1) return 0;
        if (g == n) {
            do {
                ys = addmod128(mulmod128(ys, ys, n), c, n);
                g = gcd128(x > ys ? x - ys : ys - x, n);
            } while (g == 1);
        }
        if (g != n) return g;
    }
    return 0;
}

/* heavy: also try (slow) 128-bit Pollard rho on cofactors above 2^64 */
static void factor_rec(u128 n, fac_t *F, int mult, bool heavy)
{
    if (n == 1) return;
    if (is_prime128(n)) { fac_add(F, n, mult); return; }
    for (int t = 2; t <= 127; t++) {           /* perfect powers */
        u128 c = iroot128(n, t);
        if (c < 2) break;
        if (pow_eq(c, t, n)) { factor_rec(c, F, mult * t, heavy); return; }
    }
    u128 d = 0;
    if (n <= UINT64_MAX) d = rho64((u64)n);
    else if (heavy) d = rho128(n, (u64)1 << 22);
    if (d == 0 || d == n) { F->unf = n; F->unf_e = mult; return; }
    factor_rec(d, F, mult, heavy);
    factor_rec(n / d, F, mult, heavy);
}

static int cmp_pe(const void *a, const void *b)
{
    const pe_t *x = a, *y = b;
    return x->p < y->p ? -1 : x->p > y->p;
}

static void factor128_ex(u128 m, fac_t *F, bool heavy)
{
    memset(F, 0, sizeof *F);
    if (m < 2) return;
    for (size_t i = 0; i < nsp; i++) {
        if (mod_small(m, i) == 0) {
            int e = 0;
            do { m /= sp[i]; e++; } while (mod_small(m, i) == 0);
            fac_add(F, sp[i], e);
        }
        if (m <= UINT64_MAX && (u64)sp[i] * sp[i] > (u64)m) break;
    }
    factor_rec(m, F, 1, heavy);
    qsort(F->f, (size_t)F->nf, sizeof(pe_t), cmp_pe);
}

static void factor128(u128 m, fac_t *F)
{
    factor128_ex(m, F, true);
}

/* "2^2*3 | 7": consecutive primes joined by '*', a gap shown as " | "
 * (or as '*' too when markers is false).
 * Returns the number of gaps (A073490), or -1 if part is unfactored.    */
static int fmt_factor_ex(u128 m, char *buf, size_t cap, bool markers)
{
    fac_t F;
    size_t pos = 0;
    int gaps = 0;
    if (m < 2) { snprintf(buf, cap, "%s", m ? "1" : "0"); return 0; }
    factor128(m, &F);
    for (int i = 0; i < F.nf && pos < cap; i++) {
        char pb[48];
        if (i) {
            bool gap = next_prime128(F.f[i - 1].p) != F.f[i].p;
            if (gap) gaps++;
            pos += (size_t)snprintf(buf + pos, cap - pos, (gap && markers) ? " | " : "*");
            if (pos >= cap) break;
        }
        pos += (size_t)snprintf(buf + pos, cap - pos, "%s", u128_str(F.f[i].p, pb));
        if (pos >= cap) break;
        if (F.f[i].e > 1) pos += (size_t)snprintf(buf + pos, cap - pos, "^%d", F.f[i].e);
    }
    if (F.unf && pos < cap) {
        char pb[48];
        pos += (size_t)snprintf(buf + pos, cap - pos, "%s[C%s]", F.nf ? "*" : "", u128_str(F.unf, pb));
        if (F.unf_e > 1 && pos < cap) snprintf(buf + pos, cap - pos, "^%d", F.unf_e);
        gaps = -1;
    }
    return gaps;
}

static int fmt_factor(u128 m, char *buf, size_t cap)
{
    return fmt_factor_ex(m, buf, cap, true);
}

/* short description of a gapless number: "prime", "2^16", "5*7", ... */
static void describe_gapless(u128 m, char *buf, size_t cap)
{
    gl_t g;
    if (m == 1) { snprintf(buf, cap, "1"); return; }
    if (!is_gapless(m, &g)) { snprintf(buf, cap, "NOT GAPLESS?!"); return; }
    if (g.nf == 1 && g.f[0].e == 1) { snprintf(buf, cap, "prime"); return; }
    size_t pos = 0;
    for (int i = 0; i < g.nf && i < MAXPE && pos < cap; i++) {
        char pb[48];
        pos += (size_t)snprintf(buf + pos, cap - pos, "%s%s", i ? "*" : "", u128_str(g.f[i].p, pb));
        if (pos < cap && g.f[i].e > 1) pos += (size_t)snprintf(buf + pos, cap - pos, "^%d", g.f[i].e);
    }
}

/* ------------------------------------------------------------------ */
/* scan: exhaustive enumeration of gapless numbers up to N              */
/* ------------------------------------------------------------------ */

typedef struct { u64 *v; size_t n, cap; } vec64_t;

static void vec64_push(vec64_t *v, u64 x)
{
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4096;
        v->v = realloc(v->v, v->cap * sizeof(u64));
        if (!v->v) die("out of memory");
    }
    v->v[v->n++] = x;
}

/* multiply v (whose largest prime is pr[idx-1]) by powers of pr[idx], recurse */
static void comp_extend(vec64_t *out, u128 v, size_t idx, const u64 *pr, size_t np, u64 N)
{
    if (idx >= np) return;    /* pr[idx] > 2 sqrt(N): no product can fit */
    u128 w = v * pr[idx];
    while (w <= N) {
        vec64_push(out, (u64)w);
        comp_extend(out, w, idx + 1, pr, np, N);
        w *= pr[idx];
    }
}

/* all composite gapless numbers <= N (plus the number 1), sorted */
static void gen_composite_gapless(u64 N, vec64_t *out)
{
    u64 P = 2 * (u64)isqrt128(N) + 100000;
    size_t np = 0;
    u64 *pr = (u64 *)primesieve_generate_primes(2, P, &np, UINT64_PRIMES);
    if (!pr || np == 0) die("primesieve_generate_primes failed");
    vec64_push(out, 1);
    for (size_t i = 0; i < np; i++) {
        u128 p = pr[i];
        if (p * p > N) break;
        u128 v = p;
        while (v <= N) {
            if (v != p) vec64_push(out, (u64)v);          /* p^e, e >= 2 */
            comp_extend(out, v, i + 1, pr, np, N);        /* p^e * next primes */
            v *= p;
        }
    }
    primesieve_free(pr);
    qsort(out->v, out->n, sizeof(u64), cmp_u64);
}

static size_t lower_bound64(const u64 *a, size_t n, u64 x)
{
    size_t lo = 0, hi = n;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

typedef struct { u64 first[MAXRUN]; int maxrun; } scanres_t;

static struct {
    u64 N, chunk, nchunks;
    atomic_uint_fast64_t next, done;
    const u64 *comp;
    size_t ncomp;
} SC;

static void *scan_worker(void *arg)
{
    scanres_t *res = arg;
    primesieve_iterator it;
    primesieve_init(&it);
    for (;;) {
        u64 ci = atomic_fetch_add(&SC.next, 1);
        if (ci >= SC.nchunks) break;
        u64 lo = ci * SC.chunk, hi = (ci + 1) * SC.chunk;
        if (lo < 1) lo = 1;
        if (hi > SC.N + 1) hi = SC.N + 1;         /* process pairs (g, g') with lo <= g < hi */
        if (lo >= hi) { atomic_fetch_add(&SC.done, 1); continue; }

        primesieve_jump_to(&it, lo, hi + 1000000);
        u64 p = primesieve_next_prime(&it);
        size_t idx = lower_bound64(SC.comp, SC.ncomp, lo);
        u64 c = idx < SC.ncomp ? SC.comp[idx++] : UINT64_MAX;
        u64 g;
        if (c < p) { g = c; c = idx < SC.ncomp ? SC.comp[idx++] : UINT64_MAX; }
        else       { g = p; p = primesieve_next_prime(&it); }
        while (g < hi) {
            u64 g2;
            if (c < p) { g2 = c; c = idx < SC.ncomp ? SC.comp[idx++] : UINT64_MAX; }
            else       { g2 = p; p = primesieve_next_prime(&it); }
            u64 n = g2 - g - 1;
            if (n) {
                if (n >= MAXRUN) n = MAXRUN - 1;
                if (!res->first[n]) res->first[n] = g + 1;
                if ((int)n > res->maxrun) res->maxrun = (int)n;
            }
            g = g2;
        }
        atomic_fetch_add(&SC.done, 1);
    }
    primesieve_free_iterator(&it);
    return NULL;
}

/* first[n] = a(n) for every n with a(n) <= N (0 otherwise); *maxrun = largest n seen */
static void run_scan(u64 N, int threads, u64 *first, int *maxrun, bool verbose)
{
    double t0 = now();
    vec64_t comp = { 0, 0, 0 };
    gen_composite_gapless(N + 2000000, &comp);
    if (verbose)
        fprintf(stderr, "[scan] %zu composite gapless numbers <= %" PRIu64 " generated in %.2fs\n",
                comp.n - 1, N, now() - t0);

    SC.N = N;
    SC.comp = comp.v;
    SC.ncomp = comp.n;
    u64 target = (u64)threads * 64;
    u64 chunk = (N + target) / target;
    if (chunk < ((u64)1 << 22)) chunk = (u64)1 << 22;
    SC.chunk = chunk;
    SC.nchunks = (N + chunk) / chunk;
    atomic_store(&SC.next, 0);
    atomic_store(&SC.done, 0);

    scanres_t *res = xcalloc((size_t)threads, sizeof *res);
    pthread_t *tid = xmalloc((size_t)threads * sizeof *tid);
    double t1 = now();
    for (int i = 0; i < threads; i++)
        if (pthread_create(&tid[i], NULL, scan_worker, &res[i])) die("pthread_create failed");
    if (verbose) {
        double last = 0;
        for (;;) {
            u64 done = atomic_load(&SC.done);
            double el = now() - t1;
            if (done >= SC.nchunks) break;
            if (done && (stderr_tty || el - last >= 15.0)) {
                double eta = el * (double)(SC.nchunks - done) / (double)done;
                fprintf(stderr, "%s[scan] %5.1f%%  %" PRIu64 "/%" PRIu64 " chunks  %.0fs elapsed  ETA %.0fs   %s",
                        stderr_tty ? "\r" : "", 100.0 * (double)done / (double)SC.nchunks, done, SC.nchunks,
                        el, eta, stderr_tty ? "" : "\n");
                last = el;
            }
            usleep(200000);
        }
        fprintf(stderr, "%s[scan] done: %" PRIu64 " chunks of %" PRIu64 " in %.2fs%40s\n",
                stderr_tty ? "\r" : "", SC.nchunks, chunk, now() - t1, "");
    }
    for (int i = 0; i < threads; i++) pthread_join(tid[i], NULL);

    memset(first, 0, sizeof(u64) * MAXRUN);
    *maxrun = 0;
    for (int i = 0; i < threads; i++) {
        if (res[i].maxrun > *maxrun) *maxrun = res[i].maxrun;
        for (int n = 1; n <= res[i].maxrun; n++)
            if (res[i].first[n] && (!first[n] || res[i].first[n] < first[n]))
                first[n] = res[i].first[n];
    }
    free(res);
    free(tid);
    free(comp.v);
}

/* print the a(n) table; a[] is 128-bit so search results can be merged in */
static void print_table(FILE *fp, const u128 *a, int top, u64 N, const char *note_unknown)
{
    char b1[48], b2[48], b3[48], d1[200], d2[200];
    fprintf(fp, "%4s  %-40s  %-30s  %-30s  %s\n", "n", "a(n)", "L = a(n)-1", "R = a(n)+n", "note");
    for (int n = 1; n <= top; n++) {
        if (!a[n]) {
            fprintf(fp, "%4d  %-40s  %-30s  %-30s  %s\n", n, note_unknown, "", "",
                    (n <= NKNOWN) ? "!! OEIS has a value here" : "");
            continue;
        }
        u128 L = a[n] - 1, R = a[n] + (u128)n;
        describe_gapless(L, d1, sizeof d1);
        describe_gapless(R, d2, sizeof d2);
        char lb[256], rb[256];
        snprintf(lb, sizeof lb, "%s (%s)", u128_str(L, b2), d1);
        snprintf(rb, sizeof rb, "%s (%s)", u128_str(R, b3), d2);
        const char *note;
        char nb[80];
        if (n <= NKNOWN) {
            if (a[n] == KNOWN[n]) note = "ok (OEIS)";
            else { snprintf(nb, sizeof nb, "MISMATCH! OEIS has %" PRIu64, KNOWN[n]); note = nb; }
        } else if (a[n] > N) {
            note = "NEW (structural search)";
        } else {
            note = "new";
        }
        fprintf(fp, "%4d  %-40s  %-30s  %-30s  %s\n", n, u128_str(a[n], b1), lb, rb, note);
    }
}

static void print_summary(FILE *fp, const u128 *a, int top)
{
    int k = 0;
    while (k + 1 <= top && a[k + 1]) k++;
    fprintf(fp, "\nConsecutive terms known: a(1..%d)\n", k);
    if (k) {
        fprintf(fp, "DATA: ");
        for (int n = 1; n <= k; n++) {
            char b[48];
            fprintf(fp, "%s%s", n > 1 ? ", " : "", u128_str(a[n], b));
        }
        fprintf(fp, "\n");
        if (k > NKNOWN) fprintf(fp, "(%d term%s beyond the current OEIS entry)\n", k - NKNOWN, k - NKNOWN > 1 ? "s" : "");
    }
    int miss = 0;
    for (int n = k + 1; n <= top; n++) if (!a[n]) miss++;
    if (miss) {
        fprintf(fp, "Unknown n up to %d:", top);
        int shown = 0;
        for (int n = k + 1; n <= top && shown < 40; n++)
            if (!a[n]) { fprintf(fp, " %d", n); shown++; }
        if (shown < miss) fprintf(fp, " ... (%d total)", miss);
        fprintf(fp, "\n");
    }
}

static void cmd_scan(u64 N, int threads, int maxn, const char *outfile)
{
    static u64 first[MAXRUN];
    static u128 a[MAXRUN];
    int maxrun = 0;
    double t0 = now();
    fprintf(stderr, "[scan] A137723: exhaustive scan of 1..%" PRIu64 " with %d threads\n", N, threads);
    run_scan(N, threads, first, &maxrun, true);
    fprintf(stderr, "[scan] total %.2fs, longest run found: %d\n", now() - t0, maxrun);
    for (int n = 0; n < MAXRUN; n++) a[n] = first[n];
    int top = maxrun;
    if (maxn > 0 && maxn < top) top = maxn;
    if (top < NKNOWN) top = NKNOWN;
    char nb[64];
    snprintf(nb, sizeof nb, "> %" PRIu64, N);
    printf("A137723 scan: every a(n) <= %" PRIu64 " (longest run found: %d)\n", N, maxrun);
    print_table(stdout, a, top, N, nb);
    print_summary(stdout, a, top);
    if (outfile) {
        FILE *fp = fopen(outfile, "w");
        if (!fp) die("cannot write %s", outfile);
        fprintf(fp, "A137723 scan: every a(n) <= %" PRIu64 " (longest run found: %d)\n", N, maxrun);
        print_table(fp, a, top, N, nb);
        print_summary(fp, a, top);
        fclose(fp);
        fprintf(stderr, "[scan] table written to %s\n", outfile);
    }
}

/* ------------------------------------------------------------------ */
/* search: structural search for even n                                */
/* ------------------------------------------------------------------ */

typedef struct { u128 L, R; } pair_t;
typedef struct { pair_t *v; size_t n, cap; u64 tested; } sres_t;

typedef struct { u128 E; int k; bool recurse; } stask_t;

static struct {
    u128 limit;
    u64 d;
    int n, kmax;
    stask_t *tasks;
    size_t ntasks;
    atomic_size_t next, done;
} SQ;

/* Largest k such that the odd primes <= p_k dividing d form a contiguous
 * block of consecutive primes (or none).  INT_MAX if never violated.     */
static int compute_kmax(u64 d)
{
    int first = -1, last = -1;
    for (int k = 2; k <= 60; k++) {
        u32 p = sp[k - 1];
        if (d % p == 0) {
            if (first < 0) first = last = k;
            else if (last == k - 1) last = k;
            else return k - 1;
        }
    }
    return INT_MAX;
}

/* k such that E with more than k primes can never have a PRIME partner */
static int compute_kprime(u64 d)
{
    for (int k = 2; k <= 60; k++)
        if (d % sp[k - 1] == 0) return k - 1;
    return INT_MAX;
}

static void add_pair(sres_t *r, u128 L, u128 R)
{
    if (r->n == r->cap) {
        r->cap = r->cap ? r->cap * 2 : 64;
        r->v = realloc(r->v, r->cap * sizeof(pair_t));
        if (!r->v) die("out of memory");
    }
    r->v[r->n].L = L;
    r->v[r->n].R = R;
    r->n++;
}

static bool interior_gapful(u128 L, u128 R)
{
    for (u128 m = L + 1; m < R; m++)
        if (is_gapless(m, NULL)) return false;
    return true;
}

static void test_E(u128 E, sres_t *res)
{
    u128 O;
    res->tested++;
    O = E + SQ.d;
    if (is_gapless(O, NULL) && interior_gapful(E, O)) add_pair(res, E, O);
    if (E > SQ.d) {
        O = E - SQ.d;
        if (is_gapless(O, NULL) && interior_gapful(O, E)) add_pair(res, O, E);
    }
}

/* E has k-1 primes; append powers of the k-th prime and recurse */
static void dfs(u128 E, int k, sres_t *res)
{
    u128 p = sp[k - 1], w = E;
    for (;;) {
        if (p > SQ.limit / w) break;
        w *= p;
        test_E(w, res);
        if (k + 1 <= SQ.kmax) dfs(w, k + 1, res);
    }
}

static void *search_worker(void *arg)
{
    sres_t *res = arg;
    for (;;) {
        size_t ti = atomic_fetch_add(&SQ.next, 1);
        if (ti >= SQ.ntasks) break;
        stask_t *t = &SQ.tasks[ti];
        test_E(t->E, res);
        if (t->recurse) dfs(t->E, t->k + 1, res);
        atomic_fetch_add(&SQ.done, 1);
    }
    return NULL;
}

static int cmp_pair(const void *a, const void *b)
{
    const pair_t *x = a, *y = b;
    return x->L < y->L ? -1 : x->L > y->L;
}

static void add_task(stask_t **t, size_t *n, size_t *cap, u128 E, int k, bool rec)
{
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 1024;
        *t = realloc(*t, *cap * sizeof(stask_t));
        if (!*t) die("out of memory");
    }
    (*t)[*n].E = E;
    (*t)[*n].k = k;
    (*t)[*n].recurse = rec;
    (*n)++;
}

/* Is m certainly gapful even though factor128 left a composite cofactor?
 * The cofactor's primes all exceed SMALL_LIMIT.  If the small primes found
 * are not consecutive, or the prime following the largest small prime is
 * below SMALL_LIMIT (so it was tested and does not divide m), there is a gap. */
static bool gapful_by_partial(const fac_t *F)
{
    if (F->nf == 0) return false;
    for (int i = 1; i < F->nf; i++)
        if (next_prime128(F->f[i - 1].p) != F->f[i].p) return true;
    return next_prime128(F->f[F->nf - 1].p) < SMALL_LIMIT;
}

/* Independent re-check of one pair using full factorizations (trial division,
 * Pollard rho, Miller-Rabin), i.e. not is_gapless().  Returns 1 = verified,
 * 0 = contradiction found, 2 = could not be certified (unfactored 128-bit part). */
static int verify_pair(u128 L, u128 R, char *why, size_t cap)
{
    char fb[512], b[48];
    fac_t F;
    int g, result = 1;
    g = fmt_factor(L, fb, sizeof fb);
    if (g > 0) { snprintf(why, cap, "L not gapless: %s", fb); return 0; }
    if (g < 0) { snprintf(why, cap, "L not fully factored: %s", fb); result = 2; }
    g = fmt_factor(R, fb, sizeof fb);
    if (g > 0) { snprintf(why, cap, "R not gapless: %s", fb); return 0; }
    if (g < 0) { snprintf(why, cap, "R not fully factored: %s", fb); result = 2; }
    for (u128 m = L + 1; m < R; m++) {
        factor128_ex(m, &F, false);                 /* cheap: no 128-bit rho */
        if (F.unf && gapful_by_partial(&F)) continue;   /* gap already certain */
        g = fmt_factor(m, fb, sizeof fb);           /* full effort */
        if (g == 0) { snprintf(why, cap, "%s is gapless: %s", u128_str(m, b), fb); return 0; }
        if (g < 0) {
            factor128(m, &F);
            if (!gapful_by_partial(&F)) {
                snprintf(why, cap, "%s not fully factored: %s", u128_str(m, b), fb);
                result = 2;
            }
        }
    }
    return result;
}

/* Returns a(n) (as L+1 of the smallest pair) or 0 if none with L <= limit-d.
 * All pairs found are returned in *out (sorted by L) when out != NULL.  */
static u128 run_search(int n, u128 limit, int threads, bool verbose, sres_t *out)
{
    char lb[48];
    u64 d = (u64)n + 1;
    SQ.limit = limit;
    SQ.d = d;
    SQ.n = n;
    SQ.kmax = compute_kmax(d);
    int kprime = compute_kprime(d);

    /* tasks: 2^a (k=1), 2^a 3^b (k=2), 2^a 3^b 5^c (k=3, recursing into 7, 11, ...) */
    stask_t *tasks = NULL;
    size_t ntasks = 0, cap = 0;
    /* every product is guarded against 128-bit overflow before it is formed */
    for (u128 A = 2; A <= limit; A *= 2) {
        add_task(&tasks, &ntasks, &cap, A, 1, false);
        if (SQ.kmax >= 2 && A <= limit / 3) {
            for (u128 B = A * 3; B <= limit; B *= 3) {
                add_task(&tasks, &ntasks, &cap, B, 2, false);
                if (SQ.kmax >= 3 && B <= limit / 5) {
                    for (u128 C = B * 5; C <= limit; C *= 5) {
                        add_task(&tasks, &ntasks, &cap, C, 3, SQ.kmax >= 4);
                        if (C > limit / 5) break;
                    }
                }
                if (B > limit / 3) break;
            }
        }
        if (A > limit / 2) break;
    }
    SQ.tasks = tasks;
    SQ.ntasks = ntasks;
    atomic_store(&SQ.next, 0);
    atomic_store(&SQ.done, 0);

    if (verbose) {
        char db[200];
        fmt_factor_ex(d, db, sizeof db, false);
        fprintf(stderr, "[search] n = %d, d = n+1 = %" PRIu64 " = %s, limit = %s, %d threads\n",
                n, d, db, u128_str(limit, lb), threads);
        fprintf(stderr, "[search] even endpoint E = 2^e1*3^e2*...*p_k^ek: ");
        if (SQ.kmax == INT_MAX) fprintf(stderr, "any k; ");
        else fprintf(stderr, "k <= %d; ", SQ.kmax);
        if (kprime == INT_MAX) fprintf(stderr, "a prime odd partner is possible for every k\n");
        else if (kprime == 1) fprintf(stderr, "a prime odd partner needs k = 1 (E = 2^a)\n");
        else fprintf(stderr, "a prime odd partner needs k <= %d\n", kprime);
        fprintf(stderr, "[search] %zu top-level tasks\n", ntasks);
    }

    sres_t *res = xcalloc((size_t)threads, sizeof *res);
    pthread_t *tid = xmalloc((size_t)threads * sizeof *tid);
    double t0 = now();
    for (int i = 0; i < threads; i++)
        if (pthread_create(&tid[i], NULL, search_worker, &res[i])) die("pthread_create failed");
    if (verbose) {
        double last = 0;
        for (;;) {
            size_t done = atomic_load(&SQ.done);
            if (done >= ntasks) break;
            double el = now() - t0;
            if (el > 2.0 && (stderr_tty || el - last >= 15.0)) {
                fprintf(stderr, "%s[search] %zu/%zu tasks  %.0fs   %s", stderr_tty ? "\r" : "",
                        done, ntasks, el, stderr_tty ? "" : "\n");
                last = el;
            }
            usleep(200000);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(tid[i], NULL);

    sres_t all = { NULL, 0, 0, 0 };
    for (int i = 0; i < threads; i++) {
        all.tested += res[i].tested;
        for (size_t j = 0; j < res[i].n; j++) add_pair(&all, res[i].v[j].L, res[i].v[j].R);
        free(res[i].v);
    }
    free(res);
    free(tid);
    free(tasks);
    qsort(all.v, all.n, sizeof(pair_t), cmp_pair);
    if (verbose)
        fprintf(stderr, "%s[search] enumerated %" PRIu64 " even gapless E in %.2fs, %zu run%s of exactly %d found\n",
                stderr_tty ? "\r" : "", all.tested, now() - t0, all.n, all.n == 1 ? "" : "s", n);

    /* consistency check of every pair (the enumeration never tested E itself) */
    size_t kept = 0;
    for (size_t i = 0; i < all.n; i++) {
        char b1[48], b2[48];
        if (is_gapless(all.v[i].L, NULL) && is_gapless(all.v[i].R, NULL) && interior_gapful(all.v[i].L, all.v[i].R)) {
            all.v[kept++] = all.v[i];
        } else {
            fprintf(stderr, "[search] INTERNAL ERROR: dropped invalid pair L=%s R=%s\n",
                    u128_str(all.v[i].L, b1), u128_str(all.v[i].R, b2));
        }
    }
    all.n = kept;

    u128 result = 0;
    if (all.n) result = all.v[0].L + 1;
    if (out) *out = all; else free(all.v);
    return result;
}

static void cmd_search(int n, u128 limit, int threads)
{
    char b1[48], b2[48], b3[48], d1[200], d2[200];
    if (n < 1) die("n must be >= 1");
    if (n % 2) die("search handles even n only (odd n: the two bounding gapless numbers "
                   "have the same parity, so use 'scan')");
    if (limit > ((u128)1 << 127)) die("limit must be <= 2^127");
    sres_t all;
    u128 a = run_search(n, limit, threads, true, &all);
    u64 d = (u64)n + 1;
    printf("\n=== A137723 search: n = %d (d = %" PRIu64 "), even endpoint <= %s ===\n", n, d, u128_str(limit, b1));
    if (all.n == 0) {
        printf("No run of exactly %d gapful numbers starts at or below %s.\n", n, u128_str(limit - d + 1, b1));
        printf("a(%d) > %s\n", n, u128_str(limit - d + 1, b1));
    } else {
        size_t show = all.n < 40 ? all.n : 40;
        printf("Runs of exactly %d gapful numbers (bounded by gapless L and R = L+%" PRIu64 "):\n", n, d);
        for (size_t i = 0; i < show; i++) {
            describe_gapless(all.v[i].L, d1, sizeof d1);
            describe_gapless(all.v[i].R, d2, sizeof d2);
            printf("  start %-40s L = %s (%s)  R = %s (%s)\n", u128_str(all.v[i].L + 1, b1),
                   u128_str(all.v[i].L, b2), d1, u128_str(all.v[i].R, b3), d2);
        }
        if (show < all.n) printf("  ... %zu more\n", all.n - show);
        char why[600];
        int vr = verify_pair(all.v[0].L, all.v[0].R, why, sizeof why);
        printf("\nIndependent verification of the smallest run (full factorizations): %s%s\n",
               vr == 1 ? "OK" : vr == 0 ? "FAILED: " : "INCOMPLETE: ", vr == 1 ? "" : why);
        printf("a(%d) = %s", n, u128_str(a, b1));
        if (a - 1 <= limit - d) printf("   (exhaustive: every run of exactly %d starting <= %s was enumerated)\n",
                                       n, u128_str(limit - d + 1, b2));
        else printf("   (candidate only; smaller L in (%s, %s] were not all covered)\n",
                    u128_str(limit - d, b2), u128_str(limit, b3));
        if (n <= NKNOWN) printf("OEIS value: %" PRIu64 " -> %s\n", KNOWN[n], a == KNOWN[n] ? "match" : "MISMATCH");
    }
    fflush(stdout);
    free(all.v);
}

/* ------------------------------------------------------------------ */
/* next: scan, then structurally search the missing even n             */
/* ------------------------------------------------------------------ */

static void cmd_next(u64 N, u128 limit, int threads, int maxn, const char *outfile)
{
    static u64 first[MAXRUN];
    static u128 a[MAXRUN];
    int maxrun = 0;
    char lb[48];
    if (limit > ((u128)1 << 127)) die("limit must be <= 2^127");
    fprintf(stderr, "[next] step 1: exhaustive scan to %" PRIu64 "\n", N);
    run_scan(N, threads, first, &maxrun, true);
    for (int n = 0; n < MAXRUN; n++) a[n] = first[n];
    int top = maxrun;
    if (maxn > 0 && maxn < top) top = maxn;
    if (top < NKNOWN + 1) top = NKNOWN + 1;

    fprintf(stderr, "[next] step 2: structural search (even endpoint <= %s) for each missing even n <= %d\n",
            u128_str(limit, lb), top);
    for (int n = 2; n <= top; n += 2) {
        if (a[n]) continue;
        u128 r = run_search(n, limit, threads, true, NULL);
        if (r) {
            a[n] = r;
            fprintf(stderr, "[next]   a(%d) = %s\n", n, u128_str(r, lb));
            fflush(stderr);
        } else {
            fprintf(stderr, "[next]   a(%d) > %s\n", n, u128_str(limit - (u128)n, lb));
        }
    }
    char nb[100];
    snprintf(nb, sizeof nb, "> %" PRIu64 " (odd n) / > limit-n (even n)", N);
    FILE *outs[2] = { stdout, NULL };
    if (outfile && !(outs[1] = fopen(outfile, "w"))) die("cannot write %s", outfile);
    for (int i = 0; i < 2 && outs[i]; i++) {
        fprintf(outs[i], "A137723: exhaustive scan to %" PRIu64 " plus structural search (even endpoint <= %s) for even n\n",
                N, u128_str(limit, lb));
        print_table(outs[i], a, top, N, nb);
        print_summary(outs[i], a, top);
    }
    if (outs[1]) { fclose(outs[1]); fprintf(stderr, "[next] table written to %s\n", outfile); }
}

/* ------------------------------------------------------------------ */
/* check / run                                                         */
/* ------------------------------------------------------------------ */

static void cmd_check(u128 m)
{
    char b[48], fb[600];
    int g = fmt_factor(m, fb, sizeof fb);
    printf("%s = %s\n", u128_str(m, b), fb);
    if (g < 0) printf("A073490(m): unknown (unfactored part)\n");
    else printf("A073490(%s) = %d  ->  %s\n", u128_str(m, b), g, g ? "has a prime gap (gapful)" : "no prime gap (gapless, in A073491)");
    bool gl = is_gapless(m, NULL);
    if (g >= 0 && gl != (g == 0)) printf("!! is_gapless() disagrees with the factorization\n");
    if (m > 2 && is_prime128(m)) printf("%s is prime\n", u128_str(m, b));
}

static void print_run(u128 L, u128 R)
{
    char b1[48], b2[48], d1[200], d2[200], fb[600];
    u128 n = R - L - 1;
    describe_gapless(L, d1, sizeof d1);
    describe_gapless(R, d2, sizeof d2);
    printf("Run of %s gapful numbers: %s .. %s\n", u128_str(n, b1), u128_str(L + 1, b2), u128_str(R - 1, fb));
    printf("  bounded below by %s (%s) and above by %s (%s)\n", u128_str(L, b1), d1, u128_str(R, b2), d2);
    if (n <= 200) {
        for (u128 m = L + 1; m < R; m++) {
            fmt_factor(m, fb, sizeof fb);
            printf("  %s = %s\n", u128_str(m, b1), fb);
        }
    } else {
        printf("  (%s numbers, factorizations omitted)\n", u128_str(n, b1));
    }
    if (n <= NKNOWN && n >= 1) {
        int ni = (int)n;
        if (L + 1 == KNOWN[ni]) printf("  -> this is a(%d) = %" PRIu64 " (OEIS)\n", ni, KNOWN[ni]);
        else printf("  -> a(%d) = %" PRIu64 " (OEIS), so this is a later run of that length\n", ni, KNOWN[ni]);
    } else if (n >= 1) {
        printf("  -> candidate for a(%s) = %s (not in OEIS yet)\n", u128_str(n, b1), u128_str(L + 1, b2));
    }
}

static void cmd_run(u128 m)
{
    char b[48], d[200];
    if (m < 2) die("m must be >= 2");
    if (is_gapless(m, NULL)) {
        describe_gapless(m, d, sizeof d);
        printf("%s is gapless (%s); showing the runs on both sides.\n\n", u128_str(m, b), d);
        u128 L = m - 1;
        while (L >= 1 && !is_gapless(L, NULL)) L--;
        if (L < m - 1) print_run(L, m);
        else printf("(no run below: %s is also gapless)\n", u128_str(L, b));
        u128 R = m + 1;
        while (!is_gapless(R, NULL)) R++;
        printf("\n");
        if (R > m + 1) print_run(m, R);
        else printf("(no run above: %s is also gapless)\n", u128_str(R, b));
        return;
    }
    u128 L = m - 1, R = m + 1;
    while (!is_gapless(L, NULL)) L--;
    while (!is_gapless(R, NULL)) R++;
    print_run(L, R);
}

/* ------------------------------------------------------------------ */
/* GMP mode: structural search for even n above 2^127                  */
/* ------------------------------------------------------------------ */
/*
 * Same theorem as the 128-bit search, applied more selectively so that the
 * enumeration stays small at 1000+ bits.  Let d = n+1 (odd), p_j the smallest
 * odd prime dividing d, and {p_j, ..., p_m} the maximal block of consecutive
 * primes starting at p_j that all divide d.  For an even endpoint
 * E = 2^e1 3^e2 ... p_k^ek the odd partner O = E +- d satisfies:
 *
 *   k <  j : O is coprime to p_1..p_k                 -> may be prime: enumerate E
 *   j<=k<=m: O = p_j^f.. p_k^f.. (times a consecutive continuation beyond p_k)
 *                                                      -> open-ended: enumerate E
 *   m <  k : O = p_j^f * ... * p_m^f EXACTLY            -> enumerate these O instead
 *
 * (and k > kmax is impossible, kmax = index of the next odd prime factor of d
 * minus one).  So E is enumerated only for k <= m ("E side"), and for k > m the
 * partner O ranges over products of the fixed prime set {p_j..p_m} ("O side").
 * For d = 63 = 3^2*7: j = m = 2, kmax = 3: the E side is 2^a and 2^a*3^b, the
 * O side is the powers of 3.  Both are tiny even below 2^100000.
 */
#ifdef HAVE_GMP
#include <gmp.h>
#include <math.h>

#define ZDEPTH 128

typedef struct {
    mpz_t x, p, q, c, tmp, E, O;
    mpz_t stk[ZDEPTH];
} zws_t;

static void zws_init(zws_t *w)
{
    mpz_inits(w->x, w->p, w->q, w->c, w->tmp, w->E, w->O, NULL);
    for (int i = 0; i < ZDEPTH; i++) mpz_init(w->stk[i]);
}

static void zws_clear(zws_t *w)
{
    mpz_clears(w->x, w->p, w->q, w->c, w->tmp, w->E, w->O, NULL);
    for (int i = 0; i < ZDEPTH; i++) mpz_clear(w->stk[i]);
}

/* strong probable prime, 30 rounds (deterministic answer for composites is
 * certain; a "prime" verdict above 2^64 should be certified externally)  */
static bool z_is_prime(const mpz_t m)
{
    return mpz_probab_prime_p(m, 30) > 0;
}

/* largest prime < x (x > 3); out may alias x */
static void z_prevprime(mpz_t out, const mpz_t x)
{
    mpz_sub_ui(out, x, 1);
    if (mpz_even_p(out)) mpz_sub_ui(out, out, 1);
    while (!z_is_prime(out)) mpz_sub_ui(out, out, 2);
}

static u128 z_to_u128(const mpz_t m)
{
    u64 w[2] = { 0, 0 };
    size_t cnt = 0;
    mpz_export(w, &cnt, -1, 8, 0, 0, m);
    return ((u128)w[1] << 64) | w[0];
}

static void zdesc_add(char *d, size_t cap, const mpz_t p, unsigned e)
{
    size_t len = strlen(d);
    if (len + 4 >= cap) return;
    if (len) { d[len++] = '*'; d[len] = 0; }
    size_t need = mpz_sizeinbase(p, 10) + 2;
    if (len + need + 16 >= cap) { snprintf(d + len, cap - len, "..."); return; }
    mpz_get_str(d + len, 10, p);
    len = strlen(d);
    if (e > 1) snprintf(d + len, cap - len, "^%u", e);
}

/* Is r = p^e * nextprime(p)^f * ... exactly?  p is clobbered.  On success the
 * factorization goes to desc and the number of distinct primes to *nprimes. */
static bool z_chain_from(const mpz_t r, mpz_t p, zws_t *w, char *desc, size_t cap, int *nprimes)
{
    int cnt = 0;
    mpz_set(w->x, r);
    if (desc) desc[0] = 0;
    for (;;) {
        if (!mpz_divisible_p(w->x, p)) return false;
        unsigned e = 0;
        do { mpz_divexact(w->x, w->x, p); e++; } while (mpz_divisible_p(w->x, p));
        cnt++;
        if (desc) zdesc_add(desc, cap, p, e);
        if (mpz_cmp_ui(w->x, 1) == 0) break;
        mpz_nextprime(p, p);
    }
    if (nprimes) *nprimes = cnt;
    return true;
}

/* r composite with no prime factor <= SMALL_LIMIT: gap-free only if it is
 * p^t, or a product of consecutive primes whose smallest one is among the t-1
 * largest primes below the t-th root (see root_method() for the argument). */
static bool z_root_method(const mpz_t r, zws_t *w, char *desc, size_t cap)
{
    for (unsigned long t = 2; t < 100000; t++) {
        int exact = mpz_root(w->c, r, t);
        if (mpz_cmp_ui(w->c, SMALL_LIMIT) <= 0) break;
        if (exact && z_is_prime(w->c)) {
            if (desc) { desc[0] = 0; zdesc_add(desc, cap, w->c, (unsigned)t); }
            return true;
        }
        if (exact || !z_is_prime(w->c)) z_prevprime(w->q, w->c);
        else mpz_set(w->q, w->c);
        for (unsigned long jj = 0; jj < t - 1 && mpz_cmp_ui(w->q, SMALL_LIMIT) > 0; jj++) {
            mpz_set(w->p, w->q);
            if (z_chain_from(r, w->p, w, desc, cap, NULL)) return true;
            z_prevprime(w->q, w->q);
        }
    }
    return false;
}

/* A073490(m) == 0 for an arbitrary-size m.  desc (optional) gets the structure. */
static bool z_is_gapfree(const mpz_t m, zws_t *w, char *desc, size_t cap)
{
    if (mpz_sizeinbase(m, 2) <= 127) {
        u128 v = z_to_u128(m);
        bool r = is_gapless(v, NULL);
        if (desc) { if (r) describe_gapless(v, desc, cap); else snprintf(desc, cap, "gapful"); }
        return r;
    }
    for (size_t i = 0; i < nsp; i++)
        if (mpz_divisible_ui_p(m, sp[i])) {
            mpz_set_ui(w->p, sp[i]);
            return z_chain_from(m, w->p, w, desc, cap, NULL);
        }
    if (z_is_prime(m)) { if (desc) snprintf(desc, cap, "probable prime"); return true; }
    return z_root_method(m, w, desc, cap);
}

static bool z_interior_gapful(const mpz_t L, const mpz_t R, zws_t *w)
{
    for (mpz_add_ui(w->tmp, L, 1); mpz_cmp(w->tmp, R) < 0; mpz_add_ui(w->tmp, w->tmp, 1))
        if (z_is_gapfree(w->tmp, w, NULL, 0)) return false;
    return true;
}

/* ---- search state ---- */

typedef struct { int side; unsigned a, b, f1; } ztask_t;
typedef struct { mpz_t L, R; } zpair_t;
typedef struct { zpair_t *v; size_t n, cap; u64 nodes; } zres_t;

static struct {
    mpz_t limit, limitd;          /* limit, limit + d */
    unsigned long d;
    int n, j, m, kmax, kE;
    ztask_t *tasks;
    size_t ntasks;
    atomic_size_t next, done;
} ZQ;

/* j = index of the smallest odd prime factor of d (p_1 = 2, p_2 = 3, ...),
 * m = end of the maximal block p_j, p_j+1, ..., p_m of primes dividing d,
 * kmax = (index of the next odd prime factor of d) - 1, or INT_MAX          */
static void z_analyze(u64 d, int *j, int *m, int *kmax)
{
    *j = 0; *m = 0; *kmax = INT_MAX;
    for (int k = 2; k <= (int)nsp; k++) {
        if (d % sp[k - 1] == 0) {
            if (*j == 0) { *j = k; *m = k; }
            else if (*m == k - 1) *m = k;
            else { *kmax = k - 1; break; }
        }
    }
    if (*j == 0) die("d = %" PRIu64 " has no odd prime factor below 2^16 (n too large for this mode)", d);
}

/* expected number of E = 2^e1..p_k^ek <= X (all e >= 1): volume of a simplex */
static double z_est_E(int k, double lnX)
{
    double s = lnX, denom = 1;
    for (int i = 0; i < k; i++) { s -= log((double)sp[i]); denom *= log((double)sp[i]) * (i + 1); }
    return s <= 0 ? 0 : pow(s, k) / denom;
}

static double z_est_S(int j, int m, double lnX)
{
    double s = lnX, denom = 1;
    int cnt = 0;
    for (int i = j; i <= m; i++) { s -= log((double)sp[i - 1]); cnt++; denom *= log((double)sp[i - 1]) * cnt; }
    return s <= 0 ? 0 : pow(s, cnt) / denom;
}

static void zres_add(zres_t *r, const mpz_t L, const mpz_t R)
{
    if (r->n == r->cap) {
        r->cap = r->cap ? r->cap * 2 : 16;
        r->v = realloc(r->v, r->cap * sizeof(zpair_t));
        if (!r->v) die("out of memory");
    }
    mpz_init_set(r->v[r->n].L, L);
    mpz_init_set(r->v[r->n].R, R);
    r->n++;
}

static void z_test_E(const mpz_t E, zws_t *w, zres_t *res)
{
    res->nodes++;
    mpz_add_ui(w->O, E, ZQ.d);
    if (z_is_gapfree(w->O, w, NULL, 0) && z_interior_gapful(E, w->O, w)) zres_add(res, E, w->O);
    if (mpz_cmp_ui(E, ZQ.d) > 0) {
        mpz_sub_ui(w->O, E, ZQ.d);
        if (z_is_gapfree(w->O, w, NULL, 0) && z_interior_gapful(w->O, E, w)) zres_add(res, w->O, E);
    }
}

/* E has k primes; append powers of p_(k+1) = sp[k] and recurse (E side) */
static void z_dfs(const mpz_t E, int k, zws_t *w, zres_t *res)
{
    if (k >= ZQ.kE || k >= ZDEPTH) return;
    mpz_t *wk = &w->stk[k];
    mpz_set(*wk, E);
    for (;;) {
        mpz_mul_ui(*wk, *wk, sp[k]);
        if (mpz_cmp(*wk, ZQ.limit) > 0) break;
        z_test_E(*wk, w, res);
        z_dfs(*wk, k + 1, w, res);
    }
}

/* O side: O is a product over {p_j..p_m}; its partners E = O +- d must be even
 * gap-free numbers with more than m primes (fewer were covered on the E side) */
static void z_test_O(const mpz_t O, zws_t *w, zres_t *res)
{
    res->nodes++;
    for (int sgn = 0; sgn < 2; sgn++) {
        if (sgn == 0) mpz_add_ui(w->E, O, ZQ.d);
        else { if (mpz_cmp_ui(O, ZQ.d) <= 0) break; mpz_sub_ui(w->E, O, ZQ.d); }
        if (mpz_cmp(w->E, ZQ.limit) > 0 || mpz_cmp_ui(w->E, 2) < 0) continue;
        int k = 0;
        mpz_set_ui(w->p, 2);
        if (!z_chain_from(w->E, w->p, w, NULL, 0, &k)) continue;
        if (k <= ZQ.m) continue;
        if (sgn == 0) { if (z_interior_gapful(O, w->E, w)) zres_add(res, O, w->E); }
        else          { if (z_interior_gapful(w->E, O, w)) zres_add(res, w->E, O); }
    }
}

static void z_dfs_O(const mpz_t O, int i, zws_t *w, zres_t *res)   /* multiply by p_i, i = j+1..m */
{
    if (i > ZQ.m) { z_test_O(O, w, res); return; }
    int depth = i - ZQ.j + 1;
    if (depth >= ZDEPTH) return;
    mpz_t *wk = &w->stk[depth];
    mpz_set(*wk, O);
    for (;;) {
        mpz_mul_ui(*wk, *wk, sp[i - 1]);
        if (mpz_cmp(*wk, ZQ.limitd) > 0) break;
        z_dfs_O(*wk, i + 1, w, res);
    }
}

static void *z_worker(void *arg)
{
    zres_t *res = arg;
    zws_t w;
    zws_init(&w);
    for (;;) {
        size_t ti = atomic_fetch_add(&ZQ.next, 1);
        if (ti >= ZQ.ntasks) break;
        ztask_t *t = &ZQ.tasks[ti];
        if (t->side == 0) {
            mpz_ui_pow_ui(w.stk[0], 2, t->a);
            if (t->b) { mpz_ui_pow_ui(w.tmp, 3, t->b); mpz_mul(w.stk[0], w.stk[0], w.tmp); }
            z_test_E(w.stk[0], &w, res);
            if (t->b) z_dfs(w.stk[0], 2, &w, res);
        } else {
            mpz_ui_pow_ui(w.stk[0], sp[ZQ.j - 1], t->f1);
            z_dfs_O(w.stk[0], ZQ.j + 1, &w, res);
        }
        atomic_fetch_add(&ZQ.done, 1);
    }
    zws_clear(&w);
    return NULL;
}

static int cmp_zpair(const void *a, const void *b)
{
    return mpz_cmp(((const zpair_t *)a)->L, ((const zpair_t *)b)->L);
}

static void ztask_add(ztask_t **t, size_t *n, size_t *cap, int side, unsigned a, unsigned b, unsigned f1)
{
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 1024;
        *t = realloc(*t, *cap * sizeof(ztask_t));
        if (!*t) die("out of memory");
    }
    (*t)[*n].side = side; (*t)[*n].a = a; (*t)[*n].b = b; (*t)[*n].f1 = f1;
    (*n)++;
}

/* Runs the GMP structural search.  Returns the number of pairs; *out receives
 * them sorted by L (caller frees).  a(n) = L+1 of the first pair.           */
static size_t run_search_gmp(int n, const mpz_t limit, int threads, double max_nodes, bool verbose, zres_t *out)
{
    unsigned long d = (unsigned long)n + 1;
    int j, m, kmax;
    z_analyze(d, &j, &m, &kmax);
    mpz_init_set(ZQ.limit, limit);
    mpz_init(ZQ.limitd);
    mpz_add_ui(ZQ.limitd, limit, d);
    ZQ.d = d; ZQ.n = n; ZQ.j = j; ZQ.m = m; ZQ.kmax = kmax;
    ZQ.kE = m < kmax ? m : kmax;
    if (ZQ.kE >= ZDEPTH - 2) die("too many primes in the even endpoint (kE = %d)", ZQ.kE);

    double lnX = mpz_sizeinbase(limit, 2) * log(2.0);
    double estE = 0;
    for (int k = 1; k <= ZQ.kE; k++) estE += z_est_E(k, lnX);
    double estO = z_est_S(j, m, lnX);
    if (verbose) {
        char *ls = mpz_get_str(NULL, 10, limit);
        size_t bits = mpz_sizeinbase(limit, 2);
        fprintf(stderr, "[gmp] n = %d, d = %lu, limit = %s%s (%zu bits), %d threads\n", n, d,
                strlen(ls) > 40 ? "" : ls, strlen(ls) > 40 ? "(long)" : "", bits, threads);
        free(ls);
        fprintf(stderr, "[gmp] smallest odd prime of d: p_%d = %u; block p_%d..p_%d; kmax = ", j, sp[j - 1], j, m);
        if (kmax == INT_MAX) fprintf(stderr, "unbounded\n"); else fprintf(stderr, "%d\n", kmax);
        fprintf(stderr, "[gmp] E side: k <= %d, ~%.3g candidates; O side: products of {", ZQ.kE, estE);
        for (int i = j; i <= m; i++) fprintf(stderr, "%s%u", i > j ? "," : "", sp[i - 1]);
        fprintf(stderr, "} for k > %d, ~%.3g candidates\n", m, estO);
    }
    if (estE + estO > max_nodes)
        die("estimated %.3g candidates exceeds --max-nodes %.3g; lower the limit (kE = %d primes on the E side)",
            estE + estO, max_nodes, ZQ.kE);

    /* tasks */
    ztask_t *tasks = NULL;
    size_t ntasks = 0, cap = 0;
    mpz_t t1, t2;
    mpz_inits(t1, t2, NULL);
    for (unsigned a = 1;; a++) {
        mpz_ui_pow_ui(t1, 2, a);
        if (mpz_cmp(t1, limit) > 0) break;
        ztask_add(&tasks, &ntasks, &cap, 0, a, 0, 0);
        if (ZQ.kE >= 2) {
            for (unsigned b = 1;; b++) {
                mpz_ui_pow_ui(t2, 3, b);
                mpz_mul(t2, t2, t1);
                if (mpz_cmp(t2, limit) > 0) break;
                ztask_add(&tasks, &ntasks, &cap, 0, a, b, 0);
            }
        }
    }
    for (unsigned f = 1;; f++) {
        mpz_ui_pow_ui(t1, sp[j - 1], f);
        if (mpz_cmp(t1, ZQ.limitd) > 0) break;
        ztask_add(&tasks, &ntasks, &cap, 1, 0, 0, f);
    }
    mpz_clears(t1, t2, NULL);
    ZQ.tasks = tasks; ZQ.ntasks = ntasks;
    atomic_store(&ZQ.next, 0);
    atomic_store(&ZQ.done, 0);
    if (verbose) fprintf(stderr, "[gmp] %zu tasks\n", ntasks);

    zres_t *res = xcalloc((size_t)threads, sizeof *res);
    pthread_t *tid = xmalloc((size_t)threads * sizeof *tid);
    double t0 = now();
    for (int i = 0; i < threads; i++)
        if (pthread_create(&tid[i], NULL, z_worker, &res[i])) die("pthread_create failed");
    if (verbose) {
        double last = 0;
        for (;;) {
            size_t done = atomic_load(&ZQ.done);
            if (done >= ntasks) break;
            double el = now() - t0;
            if (el > 2.0 && (stderr_tty || el - last >= 15.0)) {
                fprintf(stderr, "%s[gmp] %zu/%zu tasks  %.0fs   %s", stderr_tty ? "\r" : "", done, ntasks, el,
                        stderr_tty ? "" : "\n");
                last = el;
            }
            usleep(200000);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(tid[i], NULL);

    zres_t all = { NULL, 0, 0, 0 };
    for (int i = 0; i < threads; i++) {
        all.nodes += res[i].nodes;
        for (size_t k = 0; k < res[i].n; k++) {
            zres_add(&all, res[i].v[k].L, res[i].v[k].R);
            mpz_clears(res[i].v[k].L, res[i].v[k].R, NULL);
        }
        free(res[i].v);
    }
    free(res); free(tid); free(tasks);
    qsort(all.v, all.n, sizeof(zpair_t), cmp_zpair);
    if (verbose)
        fprintf(stderr, "%s[gmp] %" PRIu64 " candidates tested in %.2fs, %zu run%s of exactly %d found\n",
                stderr_tty ? "\r" : "", all.nodes, now() - t0, all.n, all.n == 1 ? "" : "s", n);
    mpz_clears(ZQ.limit, ZQ.limitd, NULL);
    *out = all;
    return all.n;
}

static void cmd_search_gmp(int n, const mpz_t limit, int threads, double max_nodes)
{
    if (n < 1) die("n must be >= 1");
    if (n % 2) die("search handles even n only (odd n: the two bounding gapless numbers "
                   "have the same parity, so use 'scan')");
    zres_t all;
    run_search_gmp(n, limit, threads, max_nodes, true, &all);
    unsigned long d = (unsigned long)n + 1;
    zws_t w;
    zws_init(&w);
    mpz_t lim_d;
    mpz_init(lim_d);
    mpz_sub_ui(lim_d, limit, d - 1);          /* limit - d + 1 */
    char *s1 = mpz_get_str(NULL, 10, lim_d);
    printf("\n=== A137723 search (GMP): n = %d (d = %lu), even endpoint <= 2^%zu-ish, %zu bits ===\n",
           n, d, mpz_sizeinbase(limit, 2) - 1, mpz_sizeinbase(limit, 2));
    if (all.n == 0) {
        printf("No run of exactly %d gapful numbers starts at or below %s.\n", n, s1);
        printf("a(%d) > %s\n", n, s1);
    } else {
        char dL[512], dR[512];
        size_t show = all.n < 20 ? all.n : 20;
        printf("Runs of exactly %d gapful numbers (bounded by gapless L and R = L+%lu):\n", n, d);
        for (size_t i = 0; i < show; i++) {
            z_is_gapfree(all.v[i].L, &w, dL, sizeof dL);
            z_is_gapfree(all.v[i].R, &w, dR, sizeof dR);
            char *sL = mpz_get_str(NULL, 10, all.v[i].L), *sR = mpz_get_str(NULL, 10, all.v[i].R);
            mpz_add_ui(w.tmp, all.v[i].L, 1);
            char *sA = mpz_get_str(NULL, 10, w.tmp);
            printf("  start %s\n    L = %s (%s)\n    R = %s (%s)\n", sA, sL, dL, sR, dR);
            free(sL); free(sR); free(sA);
        }
        if (show < all.n) printf("  ... %zu more\n", all.n - show);
        mpz_add_ui(w.tmp, all.v[0].L, 1);
        char *sA = mpz_get_str(NULL, 10, w.tmp);
        printf("\na(%d) = %s\n", n, sA);
        printf("  = %zu digits; exhaustive: every run of exactly %d starting <= %s was enumerated.\n",
               strlen(sA), n, s1);
        char *sL = mpz_get_str(NULL, 10, all.v[0].L), *sR = mpz_get_str(NULL, 10, all.v[0].R);
        printf("  Primality above 2^64 is probabilistic here; certify and cross-check with:\n"
               "    python3 verify_run.py %s %s\n", sL, sR);
        free(sA); free(sL); free(sR);
        for (size_t i = 0; i < all.n; i++) mpz_clears(all.v[i].L, all.v[i].R, NULL);
    }
    free(all.v);
    free(s1);
    mpz_clear(lim_d);
    zws_clear(&w);
    fflush(stdout);
}

/* parse "123", "2^1000", "10^40", "1e50", and X+Y / X-Y of those into an mpz */
static bool parse_mpz(const char *s, mpz_t out)
{
    const char *pm = strpbrk(s, "+-");
    if (pm && pm != s) {
        char head[256];
        size_t hl = (size_t)(pm - s);
        if (hl >= sizeof head) return false;
        memcpy(head, s, hl); head[hl] = 0;
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        bool ok = parse_mpz(head, a) && parse_mpz(pm + 1, b);
        if (ok) {
            if (*pm == '+') mpz_add(out, a, b);
            else { if (mpz_cmp(a, b) < 0) ok = false; else mpz_sub(out, a, b); }
        }
        mpz_clears(a, b, NULL);
        return ok;
    }
    const char *caret = strchr(s, '^');
    const char *ee = strpbrk(s, "eE");
    if (caret) {
        u128 b, e;
        if (!parse_dec(s, (size_t)(caret - s), &b) || b > 1000000) return false;
        if (!parse_dec(caret + 1, strlen(caret + 1), &e) || e > 200000) return false;
        mpz_ui_pow_ui(out, (unsigned long)b, (unsigned long)e);
        return true;
    }
    if (ee) {
        u128 e;
        char head[256];
        size_t hl = (size_t)(ee - s);
        if (hl >= sizeof head) return false;
        memcpy(head, s, hl); head[hl] = 0;
        if (mpz_set_str(out, head, 10) != 0) return false;
        if (!parse_dec(ee + 1, strlen(ee + 1), &e) || e > 100000) return false;
        mpz_t p10;
        mpz_init(p10);
        mpz_ui_pow_ui(p10, 10, (unsigned long)e);
        mpz_mul(out, out, p10);
        mpz_clear(p10);
        return true;
    }
    for (const char *c = s; *c; c++) if (*c < '0' || *c > '9') return false;
    return *s && mpz_set_str(out, s, 10) == 0;
}
/* One display line for a big number: small primes by trial division, then
 * the cofactor (probable prime, gap-free structure, or composite with a gap).
 * Returns the number of prime gaps; *exact tells whether that count is exact. */
static int z_partial_line(const mpz_t m, zws_t *w, char *buf, size_t cap, bool *exact)
{
    size_t pos = 0;
    unsigned long prev = 0;
    int gaps = 0;
    bool any = false;
    *exact = true;
    buf[0] = 0;
    mpz_set(w->x, m);
    for (size_t i = 0; i < nsp && pos < cap; i++) {
        if (!mpz_divisible_ui_p(w->x, sp[i])) continue;
        unsigned e = 0;
        do { mpz_divexact_ui(w->x, w->x, sp[i]); e++; } while (mpz_divisible_ui_p(w->x, sp[i]));
        if (any) {
            bool gap = next_prime128(prev) != sp[i];
            if (gap) gaps++;
            pos += (size_t)snprintf(buf + pos, cap - pos, gap ? " | " : "*");
        }
        if (pos < cap) pos += (size_t)snprintf(buf + pos, cap - pos, "%u", sp[i]);
        if (e > 1 && pos < cap) pos += (size_t)snprintf(buf + pos, cap - pos, "^%u", e);
        prev = sp[i];
        any = true;
        if (mpz_cmp_ui(w->x, 1) == 0) break;
    }
    if (mpz_cmp_ui(w->x, 1) > 0 && pos < cap) {
        /* every prime of the cofactor exceeds SMALL_LIMIT */
        if (any) {
            if (next_prime128(prev) < SMALL_LIMIT) gaps++; else *exact = false;
            pos += (size_t)snprintf(buf + pos, cap - pos, " | ");
        }
        mpz_set(w->E, w->x);                       /* z_is_gapfree clobbers w->x */
        char d2[600];
        char *cs = mpz_get_str(NULL, 10, w->E);
        if (z_is_prime(w->E)) {
            if (pos < cap) snprintf(buf + pos, cap - pos, "%s (probable prime)", cs);
        } else if (z_is_gapfree(w->E, w, d2, sizeof d2)) {
            if (pos < cap) snprintf(buf + pos, cap - pos, "%s", d2);
        } else {
            gaps++;                                /* at least one gap inside the cofactor */
            *exact = false;
            if (pos < cap) snprintf(buf + pos, cap - pos, "[C%s: composite with a prime gap]", cs);
        }
        free(cs);
    }
    return gaps;
}

static void cmd_check_gmp(const mpz_t m)
{
    zws_t w;
    zws_init(&w);
    char line[4096], desc[600];
    bool exact;
    int gaps = z_partial_line(m, &w, line, sizeof line, &exact);
    char *ms = mpz_get_str(NULL, 10, m);
    printf("%s = %s\n", ms, line);
    bool gf = z_is_gapfree(m, &w, desc, sizeof desc);
    printf("A073490 %s %d  ->  %s\n", exact ? "=" : ">=", gaps,
           gf ? "no prime gap (gapless, in A073491)" : "has a prime gap (gapful)");
    if (gf) printf("structure: %s\n", desc);
    if (gf == (gaps > 0)) printf("!! z_is_gapfree() disagrees with the display factorization\n");
    free(ms);
    zws_clear(&w);
}

static void z_print_run(const mpz_t L, const mpz_t R, zws_t *w)
{
    char dL[600], dR[600], line[4096];
    bool exact;
    mpz_t n, m;
    mpz_inits(n, m, NULL);
    mpz_sub(n, R, L);
    mpz_sub_ui(n, n, 1);
    z_is_gapfree(L, w, dL, sizeof dL);
    z_is_gapfree(R, w, dR, sizeof dR);
    char *sn = mpz_get_str(NULL, 10, n), *sL = mpz_get_str(NULL, 10, L), *sR = mpz_get_str(NULL, 10, R);
    mpz_add_ui(m, L, 1);
    char *sA = mpz_get_str(NULL, 10, m);
    printf("Run of %s gapful numbers starting at %s\n", sn, sA);
    printf("  bounded below by %s (%s)\n  and above by %s (%s)\n", sL, dL, sR, dR);
    if (mpz_cmp_ui(n, 300) <= 0) {
        for (mpz_add_ui(m, L, 1); mpz_cmp(m, R) < 0; mpz_add_ui(m, m, 1)) {
            char *s = mpz_get_str(NULL, 10, m);
            z_partial_line(m, w, line, sizeof line, &exact);
            printf("  %s = %s\n", s, line);
            free(s);
        }
    }
    if (mpz_cmp_ui(n, 1) >= 0 && mpz_fits_ulong_p(n))
        printf("  -> candidate for a(%s) = %s\n", sn, sA);
    free(sn); free(sL); free(sR); free(sA);
    mpz_clears(n, m, NULL);
}

static void cmd_run_gmp(const mpz_t m0)
{
    zws_t w;
    zws_init(&w);
    mpz_t L, R;
    mpz_inits(L, R, NULL);
    char desc[600];
    if (z_is_gapfree(m0, &w, desc, sizeof desc)) {
        char *s = mpz_get_str(NULL, 10, m0);
        printf("%s is gapless (%s); showing the runs on both sides.\n\n", s, desc);
        free(s);
        mpz_sub_ui(L, m0, 1);
        while (!z_is_gapfree(L, &w, NULL, 0)) mpz_sub_ui(L, L, 1);
        z_print_run(L, m0, &w);
        printf("\n");
        mpz_add_ui(R, m0, 1);
        while (!z_is_gapfree(R, &w, NULL, 0)) mpz_add_ui(R, R, 1);
        z_print_run(m0, R, &w);
    } else {
        mpz_sub_ui(L, m0, 1);
        while (!z_is_gapfree(L, &w, NULL, 0)) mpz_sub_ui(L, L, 1);
        mpz_add_ui(R, m0, 1);
        while (!z_is_gapfree(R, &w, NULL, 0)) mpz_add_ui(R, R, 1);
        z_print_run(L, R, &w);
    }
    mpz_clears(L, R, NULL);
    zws_clear(&w);
}
#endif /* HAVE_GMP */

/* ------------------------------------------------------------------ */
/* selftest                                                            */
/* ------------------------------------------------------------------ */

static int tfail;

static void expect(bool cond, const char *fmt, ...)
{
    va_list ap;
    if (cond) return;
    tfail++;
    va_start(ap, fmt);
    fputs("  FAIL: ", stdout);
    vprintf(fmt, ap);
    va_end(ap);
    fputc('\n', stdout);
}

/* brute-force gapful indicator for 2..M via smallest-prime-factor sieve */
static uint8_t *brute_gapful(u64 M)
{
    u32 *spf = xcalloc(M + 1, sizeof(u32));
    for (u64 i = 2; i <= M; i++)
        if (!spf[i])
            for (u64 j = i; j <= M; j += i)
                if (!spf[j]) spf[j] = (u32)i;
    uint8_t *g = xcalloc(M + 1, 1);
    for (u64 m = 2; m <= M; m++) {
        u64 x = m;
        u32 prev = 0;
        bool gap = false;
        while (x > 1) {
            u32 p = spf[x];
            if (prev) {
                u64 q = prev + 1;
                while (spf[q] != q) q++;          /* next prime after prev */
                if (q != p) gap = true;
            }
            prev = p;
            while (x % p == 0) x /= p;
        }
        g[m] = gap;
    }
    free(spf);
    return g;
}

static void cmd_selftest(int threads)
{
    double t0 = now();
    printf("A137723 selftest\n");

    /* 1. primality */
    printf("[1] primality tests\n");
    for (u64 x = 0; x < SMALL_LIMIT; x++)
        expect(is_prime64(x) == (x < SMALL_LIMIT && small_isp[x]), "is_prime64(%" PRIu64 ")", x);
    expect(is_prime64(((u64)1 << 61) - 1), "2^61-1 prime");
    expect(is_prime64(18446744073709551557ULL), "2^64-59 prime");
    expect(!is_prime64(3215031751ULL), "3215031751 composite");
    expect(!is_prime64(3825123056546413051ULL), "3825123056546413051 composite (spsp to bases 2..23)");
    expect(!is_prime64(UINT64_MAX), "2^64-1 composite");
    expect(is_prime128(((u128)1 << 89) - 1), "2^89-1 prime");
    expect(is_prime128(((u128)1 << 107) - 1), "2^107-1 prime");
    expect(is_prime128(((u128)1 << 127) - 1), "2^127-1 prime");
    expect(!is_prime128(((u128)1 << 67) - 1), "2^67-1 composite");
    expect(!is_prime128(((u128)1 << 64) + 1), "2^64+1 composite");
    {
        u128 psp;
        parse_u128("318665857834031151167461", &psp);      /* spsp to the first 12 prime bases */
        expect(!is_prime128(psp), "318665857834031151167461 composite");
        parse_u128("3317044064679887385961981", &psp);     /* spsp to the first 13 prime bases */
        expect(!is_prime128(psp), "3317044064679887385961981 composite");
    }
    expect(next_prime128(2) == 3 && next_prime128(65521) == 65537 && prev_prime128(65537) == 65521,
           "next/prev prime around 2^16");
    expect(next_prime128(((u128)1 << 64)) == ((u128)1 << 64) + 13, "nextprime(2^64) = 2^64+13");
    expect(prev_prime128(((u128)1 << 64)) == ((u128)1 << 64) - 59, "prevprime(2^64) = 2^64-59");

    /* 2. gapless test vs brute force */
    u64 M = 2000000;
    printf("[2] is_gapless vs brute force for 1..%" PRIu64 "\n", M);
    uint8_t *bg = brute_gapful(M);
    for (u64 m = 1; m <= M; m++)
        expect(is_gapless(m, NULL) == !bg[m], "is_gapless(%" PRIu64 ")", m);
    {
        gl_t g;
        expect(is_gapless(2940, &g) && g.nf == 4 && g.f[3].p == 7 && g.f[3].e == 2, "2940 = 2^2*3*5*7^2 structure");
        expect(!is_gapless(84, NULL) && !is_gapless(10, NULL) && is_gapless(35, NULL) && is_gapless(1, NULL), "small cases");
    }

    /* 3. root method: consecutive-prime products with large primes */
    printf("[3] gapless test on large consecutive-prime products\n");
    {
        u128 p = prev_prime128((u128)1 << 32), q = next_prime128(p), r = next_prime128(q), s = next_prime128(r);
        expect(is_gapless(p * q, NULL), "p*q (p~2^32)");
        expect(!is_gapless(p * r, NULL), "p*r skips q");
        expect(is_gapless(p * p, NULL), "p^2");
        expect(is_gapless(p * p * q, NULL), "p^2*q");
        expect(is_gapless(p * q * q, NULL), "p*q^2");
        expect(is_gapless(p * q * r, NULL), "p*q*r");
        expect(!is_gapless(p * q * s, NULL), "p*q*s skips r");
        expect(is_gapless(p * p * p, NULL), "p^3");
        expect(!is_gapless(p * p * r, NULL), "p^2*r skips q");
        u128 P = prev_prime128((u128)1 << 62), Q = next_prime128(P), R = next_prime128(Q);
        expect(is_gapless(P * Q, NULL), "P*Q (P~2^62)");
        expect(!is_gapless(P * R, NULL), "P*R skips Q");
        expect(is_gapless(P * P, NULL), "P^2 (~2^124)");
        u128 a = prev_prime128((u128)1 << 41), b = next_prime128(a), c = next_prime128(b);
        expect(is_gapless(a * b * c, NULL), "a*b*c (~2^123)");
        expect(is_gapless(a * a * b, NULL), "a^2*b (~2^123)");
        expect(!is_gapless(a * b * next_prime128(c), NULL), "a*b*d skips c");
        char fb[300];
        expect(fmt_factor(p * q * s, fb, sizeof fb) == 1, "fmt_factor gap count p*q*s: %s", fb);
        expect(fmt_factor(2940, fb, sizeof fb) == 0 && strcmp(fb, "2^2*3*5*7^2") == 0, "fmt_factor(2940) = %s", fb);
        expect(fmt_factor(110, fb, sizeof fb) == 2 && strcmp(fb, "2 | 5 | 11") == 0, "fmt_factor(110) = %s", fb);
    }

    /* 4. scan vs brute force and OEIS */
    u64 S = M - 5000;
    printf("[4] scan to %" PRIu64 " vs brute force and OEIS\n", S);
    {
        static u64 first[MAXRUN];
        static u64 bfirst[MAXRUN];
        int maxrun = 0, bmax = 0;
        u64 g = 1;
        memset(bfirst, 0, sizeof bfirst);
        for (u64 m = 2; m <= M; m++) {
            if (bg[m]) continue;
            if (g <= S) {
                u64 n = m - g - 1;
                if (n) { if (!bfirst[n]) bfirst[n] = g + 1; if ((int)n > bmax) bmax = (int)n; }
            }
            g = m;
        }
        run_scan(S, threads, first, &maxrun, false);
        expect(maxrun == bmax, "scan maxrun %d vs brute %d", maxrun, bmax);
        for (int n = 1; n < MAXRUN; n++)
            expect(first[n] == bfirst[n], "scan first[%d] = %" PRIu64 " vs brute %" PRIu64, n, first[n], bfirst[n]);
        for (int n = 1; n <= NKNOWN; n++)
            expect(first[n] == KNOWN[n], "a(%d): scan %" PRIu64 " vs OEIS %" PRIu64, n, first[n], KNOWN[n]);
        expect(first[32] == 0, "a(32) not below %" PRIu64, S);
    }
    free(bg);

    /* 5. structural search reproduces every even OEIS term */
    printf("[5] structural search for even n <= 30 with limit 2^20\n");
    for (int n = 2; n <= NKNOWN; n += 2) {
        u128 a = run_search(n, (u128)1 << 20, threads, false, NULL);
        char b[48];
        expect(a == KNOWN[n], "search a(%d) = %s vs OEIS %" PRIu64, n, u128_str(a, b), KNOWN[n]);
    }
    {
        expect(compute_kmax(33) == 4 && compute_kprime(33) == 1, "kmax(33)=4, kprime(33)=1");
        expect(compute_kmax(35) == INT_MAX && compute_kprime(35) == 2, "kmax(35)=inf, kprime(35)=2");
        expect(compute_kmax(21) == 3 && compute_kmax(15) == INT_MAX && compute_kmax(55) == 4, "kmax 21/15/55");
        expect(compute_kmax(37) == INT_MAX && compute_kprime(37) == 11, "kmax(37)=inf, kprime(37)=11");
    }

    /* 6. the structural search must agree with the exhaustive scan for every even n */
    {
        static u64 first[MAXRUN];
        int maxrun = 0, checked = 0;
        u64 T = 10000000;
        printf("[6] structural search vs exhaustive scan for every even n up to %" PRIu64 "\n", T);
        run_scan(T, threads, first, &maxrun, false);
        for (int n = 2; n <= maxrun; n += 2) {
            u64 d = (u64)n + 1;
            u128 a = run_search(n, T, threads, false, NULL);
            char b[48];
            if (first[n] && first[n] - 1 <= T - d)
                expect(a == first[n], "n=%d: search %s vs scan %" PRIu64, n, u128_str(a, b), first[n]);
            else
                expect(a == 0 || a == first[n], "n=%d: search %s but scan %" PRIu64, n, u128_str(a, b), first[n]);
            checked++;
        }
        printf("    %d even values of n compared (longest run below %" PRIu64 ": %d)\n", checked, T, maxrun);
    }

#ifdef HAVE_GMP
    /* 7. GMP mode: big-number gap-free test, and the E-side/O-side enumeration
     *    must agree with the 128-bit search and with the exhaustive scan       */
    printf("[7] GMP mode\n");
    {
        zws_t w;
        zws_init(&w);
        mpz_t P, Q, R2, X, Y;
        mpz_inits(P, Q, R2, X, Y, NULL);
        mpz_ui_pow_ui(X, 2, 521); mpz_sub_ui(X, X, 1);
        expect(z_is_gapfree(X, &w, NULL, 0), "2^521-1 (Mersenne prime) gap-free");
        mpz_ui_pow_ui(X, 2, 500); mpz_sub_ui(X, X, 1);
        expect(!z_is_gapfree(X, &w, NULL, 0), "2^500-1 = 3*5*... gapful");
        mpz_ui_pow_ui(X, 3, 100);
        expect(z_is_gapfree(X, &w, NULL, 0), "3^100 gap-free");
        mpz_ui_pow_ui(X, 2, 100); mpz_ui_pow_ui(Y, 3, 50); mpz_mul(X, X, Y); mpz_ui_pow_ui(Y, 5, 20); mpz_mul(X, X, Y);
        expect(z_is_gapfree(X, &w, NULL, 0), "2^100*3^50*5^20 gap-free");
        mpz_mul_ui(Y, X, 11);
        expect(!z_is_gapfree(Y, &w, NULL, 0), "2^100*3^50*5^20*11 gapful");
        mpz_ui_pow_ui(P, 2, 70); mpz_nextprime(P, P); mpz_nextprime(Q, P); mpz_nextprime(R2, Q);
        mpz_mul(X, P, Q);
        expect(z_is_gapfree(X, &w, NULL, 0), "P*Q with P ~ 2^70 gap-free (GMP root method, t=2)");
        mpz_mul(X, P, R2);
        expect(!z_is_gapfree(X, &w, NULL, 0), "P*R skips Q");
        mpz_mul(X, P, Q); mpz_mul(X, X, R2);
        expect(z_is_gapfree(X, &w, NULL, 0), "P*Q*R gap-free (t=3)");
        mpz_mul(X, P, P); mpz_mul(X, X, Q);
        expect(z_is_gapfree(X, &w, NULL, 0), "P^2*Q gap-free");
        mpz_mul(X, P, P); mpz_mul(X, X, R2);
        expect(!z_is_gapfree(X, &w, NULL, 0), "P^2*R skips Q");
        mpz_mul(X, P, P); mpz_mul(X, X, P);
        expect(z_is_gapfree(X, &w, NULL, 0), "P^3 gap-free");
        char desc[256];
        mpz_ui_pow_ui(X, 2, 200); mpz_mul_ui(X, X, 9);
        expect(z_is_gapfree(X, &w, desc, sizeof desc) && strcmp(desc, "2^200*3^2") == 0, "describe 2^200*3^2: %s", desc);
        mpz_clears(P, Q, R2, X, Y, NULL);
        zws_clear(&w);

        mpz_t lim;
        mpz_init(lim);
        for (int n = 2; n <= NKNOWN; n += 2) {
            zres_t all;
            mpz_ui_pow_ui(lim, 2, 20);
            size_t cnt = run_search_gmp(n, lim, threads, 1e12, false, &all);
            bool ok = cnt > 0 && mpz_cmp_ui(all.v[0].L, (unsigned long)KNOWN[n] - 1) == 0;
            expect(ok, "GMP search a(%d) vs OEIS %" PRIu64, n, KNOWN[n]);
            for (size_t i = 0; i < all.n; i++) mpz_clears(all.v[i].L, all.v[i].R, NULL);
            free(all.v);
        }
        {
            zres_t all;
            mpz_ui_pow_ui(lim, 2, 64);
            size_t cnt = run_search_gmp(32, lim, threads, 1e12, false, &all);
            bool ok = cnt == 2 && mpz_cmp_ui(all.v[0].L, 18014398509481951ULL) == 0;
            expect(ok, "GMP search a(32) to 2^64");
            for (size_t i = 0; i < all.n; i++) mpz_clears(all.v[i].L, all.v[i].R, NULL);
            free(all.v);
        }
        {
            static u64 first[MAXRUN];
            int maxrun = 0, checked = 0;
            u64 T = 10000000;
            run_scan(T, threads, first, &maxrun, false);
            mpz_set_ui(lim, T);
            for (int n = 2; n <= maxrun; n += 2) {
                u64 d = (u64)n + 1;
                zres_t all;
                run_search_gmp(n, lim, threads, 1e12, false, &all);
                u64 a = all.n ? mpz_get_ui(all.v[0].L) + 1 : 0;
                if (first[n] && first[n] - 1 <= T - d)
                    expect(a == first[n], "GMP n=%d: search %" PRIu64 " vs scan %" PRIu64, n, a, first[n]);
                else
                    expect(a == 0 || a == first[n], "GMP n=%d: search %" PRIu64 " but scan %" PRIu64, n, a, first[n]);
                for (size_t i = 0; i < all.n; i++) mpz_clears(all.v[i].L, all.v[i].R, NULL);
                free(all.v);
                checked++;
            }
            printf("    %d even values of n compared with the exhaustive scan\n", checked);
        }
        mpz_clear(lim);
    }
#endif

    printf("%s (%d failure%s, %.1fs)\n", tfail ? "SELFTEST FAILED" : "SELFTEST PASSED", tfail, tfail == 1 ? "" : "s", now() - t0);
    exit(tfail ? 1 : 0);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fprintf(stderr,
        "usage:\n"
        "  a137723 scan N [-t T] [-m MAXN] [-o FILE]   exhaustive: every a(n) <= N\n"
        "  a137723 search n [n ...] [-l LIMIT] [-t T]  structural search for even n (default LIMIT 2^64);\n"
        "                                              LIMIT above 2^127 (e.g. 2^1000) uses GMP; --gmp forces it,\n"
        "                                              --max-nodes X caps the estimated enumeration (default 5e9)\n"
        "  a137723 next [-N SCAN] [-l LIMIT] [-t T] [-m MAXN] [-o FILE]\n"
        "                                              scan to SCAN (default 10^10), then search missing even n\n"
        "  a137723 check m                             factor m and count its prime gaps (A073490)\n"
        "  a137723 run m                               show the maximal run of gapful numbers around m\n"
        "                                              (check/run accept numbers above 2^128 when built with GMP)\n"
        "  a137723 selftest [-t T]\n"
        "numbers: decimal, 2^k, 10^k, 1e15, or sums/differences like 2^54-32;  -t defaults to the number of CPUs\n");
    exit(2);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage();
    init_small_primes();
    stderr_tty = isatty(2);
    int threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (threads < 1) threads = 1;
    const char *cmd = argv[1];
    const char *pos[64];
    int npos = 0, maxn = 0;
    const char *outfile = NULL;
    const char *limit_str = NULL;
    u128 limit = (u128)1 << 64;
    u64 scanN = 10000000000ULL;
    bool force_gmp = false;
    double max_nodes = 5e9;

    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) { threads = atoi(argv[++i]); if (threads < 1) threads = 1; }
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) maxn = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "-l") && i + 1 < argc) limit_str = argv[++i];
        else if (!strcmp(argv[i], "--gmp")) force_gmp = true;
        else if (!strcmp(argv[i], "--max-nodes") && i + 1 < argc) max_nodes = atof(argv[++i]);
        else if (!strcmp(argv[i], "-N") && i + 1 < argc) scanN = arg_u64(argv[++i], "scan limit");
        else if (argv[i][0] == '-' && argv[i][1]) usage();
        else if (npos < 64) pos[npos++] = argv[i];
    }

    bool gmp_mode = force_gmp;
#ifdef HAVE_GMP
    mpz_t zlimit;
    mpz_init(zlimit);
    mpz_ui_pow_ui(zlimit, 2, 64);
    if (limit_str) {
        if (!parse_mpz(limit_str, zlimit)) die("bad limit: '%s'", limit_str);
        if (mpz_sizeinbase(zlimit, 2) <= 127) limit = z_to_u128(zlimit);
        else gmp_mode = true;
    }
#else
    if (limit_str) limit = arg_u128(limit_str, "limit");
    if (limit > ((u128)1 << 127)) die("limit above 2^127 needs the GMP build (make with gmp.h available)");
    if (force_gmp) die("this build has no GMP support");
#endif

    if (!strcmp(cmd, "scan")) {
        if (npos != 1) usage();
        u64 N = arg_u64(pos[0], "N");
        if (N < 100) N = 100;
        if (N > ((u64)1 << 63)) die("N must be <= 2^63");
        cmd_scan(N, threads, maxn, outfile);
    } else if (!strcmp(cmd, "search")) {
        if (npos < 1) usage();
        for (int i = 0; i < npos; i++) {
            u128 n = arg_u128(pos[i], "n");
            if (n > 100000) die("n too large");
#ifdef HAVE_GMP
            if (gmp_mode) { cmd_search_gmp((int)n, zlimit, threads, max_nodes); continue; }
#endif
            cmd_search((int)n, limit, threads);
        }
    } else if (!strcmp(cmd, "next")) {
        if (npos) usage();
        if (gmp_mode) die("'next' uses the 128-bit search; use 'search n -l LIMIT' above 2^127");
        cmd_next(scanN, limit, threads, maxn, outfile);
    } else if (!strcmp(cmd, "check") || !strcmp(cmd, "run")) {
        if (npos != 1) usage();
        bool is_check = !strcmp(cmd, "check");
#ifdef HAVE_GMP
        mpz_t zm;
        mpz_init(zm);
        if (!parse_mpz(pos[0], zm)) die("bad m: '%s'", pos[0]);
        if (mpz_sizeinbase(zm, 2) > 127) {
            if (is_check) cmd_check_gmp(zm); else cmd_run_gmp(zm);
            mpz_clear(zm);
            return 0;
        }
        mpz_clear(zm);
#endif
        if (is_check) cmd_check(arg_u128(pos[0], "m"));
        else cmd_run(arg_u128(pos[0], "m"));
    } else if (!strcmp(cmd, "selftest")) {
        cmd_selftest(threads);
    } else {
        usage();
    }
    return 0;
}
