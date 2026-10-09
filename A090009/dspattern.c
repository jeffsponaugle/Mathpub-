/* dspattern.c -- exhaustive search for digit-sum prime chains, by digit pattern
 *
 * Companion to dschain.c.  dschain looks at every prime in a range; dspattern
 * looks only at numbers whose digits allow a chain of the requested length at
 * all.  For long chains that is a vanishing fraction: below 2^64 a chain of 10
 * primes can only start in five narrow digit patterns (176 million numbers in
 * all, searched in under a second), and a chain of 11 or more cannot start
 * anywhere below 8e23.
 *
 * Why patterns exist.  Write p = H*10^4 + x with 0 <= x < 10^4.  Each step
 * adds a digit sum, at most 225 here, so a chain of up to 40 terms moves less
 * than 10^4: it carries into H+1 at most once.  Every term's digit sum is
 * ds(H) plus the digit sum of its low four digits -- or ds(H+1) plus that,
 * after the carry -- and ds(H+1) = ds(H) + 1 - 9t, where t counts the
 * trailing 9s of H.  So whether the first L terms are all odd and prime to 3
 * and 5 (as primes above 5 must be) depends only on the triple
 * (S = ds(H), t, x).  There are few enough triples to classify them all at
 * startup; the admissible ones are the patterns for L.  A chain that never
 * carries does not depend on t, so its pattern is just (S, x).
 *
 * The search.  For each pattern, every H in range with ds(H) = S -- and
 * exactly t trailing 9s, for a chain that carries -- is generated directly,
 * in increasing order, never scanned for.  Each resulting p is sieved against
 * the primes 7..61 on all L terms at once (one mask test per prime),
 * survivors get a base-2 strong probable-prime test per term, and a start
 * whose L terms all pass is proven with deterministic Miller-Rabin and
 * followed to its true length, exactly as dschain does.  Every rejection along
 * the way is a proof of compositeness, so the search is exhaustive: --first
 * reports the smallest start in the range, and an empty result means there is
 * none.
 *
 * Positions are 128-bit, up to 10^30.  Miller-Rabin with the 12 smallest
 * prime bases is deterministic below 3.18e23 (and 64-bit numbers are all
 * below that); with 13 bases it is deterministic below PSI13 =
 * 3317044064679887385961981 (Sorenson and Webster, 2015).  Above PSI13 the
 * search is still exhaustive -- every rejection is still a proof -- but the
 * terms of a hit are only 13-base strong probable primes, and the output
 * says so: certify them with verify_chain.py, whose Pratt certificates work
 * at any size.  Starts below 10^4, where H = 0 and terms could equal a sieve
 * prime, are checked by brute force.
 *
 * The Montgomery arithmetic, primality tests and number parsing are
 * gapsieve.c's.
 *
 * Build: cc -O3 -pthread -o dspattern dspattern.c   (or just: make)
 */

#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE          /* sysconf(_SC_NPROCESSORS_ONLN) on macOS */

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef unsigned __int128 u128;
#define U128_MAX (~(u128)0)

/* Miller-Rabin with the 13 smallest prime bases is proven below this */
#define PSI13 ((u128)179817 << 64 | 0x51adc5b22410a5fdULL)

#define BLOCK       10000u        /* p = H*BLOCK + x */
#define MAX_DS      270           /* bounds the digit sum of anything below 10^30 */
#define MAXLEN      36            /* (MAXLEN-1)*MAX_DS < BLOCK: one carry at most */
#define NDIG        32            /* decimal digits the generator handles */

/* starts are capped here, so that every number a chain check touches
 * (terms, their extension, predecessors) stays below 10^30, well inside the
 * two-word arithmetic and the generator's 32 digits */
#define E30     ((u128)1000000000000000ULL * 1000000000000000ULL)
#define PAT_TOP (E30 - 1 - 4 * (u128)MAXLEN * MAX_DS)

/* ------------------------------------------------------------------ */
/* decimal text                                                        */
/* ------------------------------------------------------------------ */

/* decimal text of v into buf, which needs 40 bytes (gapsieve.c) */
static char *u128_str(char *buf, u128 v)
{
    const uint64_t e19 = 10000000000000000000ULL;
    if (v <= UINT64_MAX) {
        snprintf(buf, 40, "%" PRIu64, (uint64_t)v);
        return buf;
    }
    u128 hi = v / e19;
    uint64_t lo = (uint64_t)(v % e19);
    if (hi <= UINT64_MAX) {
        snprintf(buf, 40, "%" PRIu64 "%019" PRIu64, (uint64_t)hi, lo);
    } else {                                    /* 39 digits, the first 1..3 */
        buf[0] = (char)('0' + (int)(hi / e19));
        snprintf(buf + 1, 39, "%019" PRIu64 "%019" PRIu64,
                 (uint64_t)(hi % e19), lo);
    }
    return buf;
}

/* u128_str into a temporary that lives until the end of the enclosing block */
#define U128S(v) u128_str((char[40]){ 0 }, (v))

/* 1234567 -> "1,234,567" */
static char *grouped(char *buf, u128 v)
{
    char tmp[40];
    u128_str(tmp, v);
    int len = (int)strlen(tmp), o = 0;
    for (int i = 0; i < len; i++) {
        if (i && (len - i) % 3 == 0)
            buf[o++] = ',';
        buf[o++] = tmp[i];
    }
    buf[o] = 0;
    return buf;
}
#define GROUPED(v) grouped((char[56]){ 0 }, (v))

/*
 * gapsieve.c's exact parser: 1e20, 2.5e21, 300T, 1_000_000 or plain digits,
 * converted exactly (no floating point), and "max".  A trailing E alone is
 * exa (1e18), as in dschain; 3E9 is still 3e9.  -1 on garbage.
 */
static int parse_num(const char *arg, u128 *out)
{
    char buf[128];
    size_t n = 0;
    for (const char *s = arg; *s; s++) {
        if (*s == ',' || *s == '_' || *s == ' ' || *s == '\'')
            continue;
        if (n + 1 >= sizeof buf)
            return -1;
        buf[n++] = *s;
    }
    buf[n] = '\0';
    if (!n)
        return -1;
    if (!strcmp(buf, "max")) {
        *out = U128_MAX;
        return 0;
    }

    /* the mantissa's digits, with its decimal point moved into exp */
    u128 m = 0;
    int exp = 0, digits = 0, point = 0;
    const char *s = buf;
    for (; *s; s++) {
        if (*s == '.' && !point) {
            point = 1;
            continue;
        }
        if (!isdigit((unsigned char)*s))
            break;
        unsigned dg = (unsigned)(*s - '0');
        if (m > (U128_MAX - dg) / 10)
            return -1;
        m = m * 10 + dg;
        digits++;
        exp -= point;
    }
    if (!digits)
        return -1;
    if ((*s == 'e' || *s == 'E') && s[1]) {
        int neg = 0, e = 0, ed = 0;
        s++;
        if (*s == '+' || *s == '-')
            neg = *s++ == '-';
        for (; isdigit((unsigned char)*s); s++, ed++)
            if ((e = e * 10 + (*s - '0')) > 1000)
                return -1;
        if (!ed)
            return -1;
        exp += neg ? -e : e;
    }
    if (*s) {
        switch (toupper((unsigned char)*s)) {
        case 'K': exp += 3;  break;
        case 'M': exp += 6;  break;
        case 'B':
        case 'G': exp += 9;  break;
        case 'T': exp += 12; break;
        case 'P': exp += 15; break;
        case 'E': exp += 18; break;
        default:  return -1;
        }
        if (*++s)
            return -1;
    }
    for (; exp > 0; exp--) {
        if (m > U128_MAX / 10)
            return -1;
        m *= 10;
    }
    /* the last digit dropped is the most significant one: round on it */
    int up = 0;
    for (; exp < 0; exp++) {
        up = m % 10 >= 5;
        m /= 10;
    }
    *out = m + (u128)up;
    return 0;
}

static u128 num_arg(const char *arg, const char *what)
{
    u128 v;
    if (parse_num(arg, &v)) {
        fprintf(stderr, "error: %s is not a number: %s\n"
                        "       (1e20, 2.5e21, 300T, 1_000_000, plain digits, or max)\n",
                what, arg);
        exit(2);
    }
    return v;
}

/* ------------------------------------------------------------------ */
/* digit sums                                                          */
/* ------------------------------------------------------------------ */

static uint8_t ds4[BLOCK];
static u128 P10[NDIG + 1];

static void tables_init(void)
{
    for (uint32_t i = 0; i < BLOCK; i++) {
        uint32_t n = i, s = 0;
        while (n) { s += n % 10; n /= 10; }
        ds4[i] = (uint8_t)s;
    }
    P10[0] = 1;
    for (int i = 1; i <= NDIG; i++)
        P10[i] = P10[i - 1] * 10;
}

static inline unsigned digitsum64(uint64_t n)
{
    unsigned s = 0;
    while (n) { s += ds4[n % BLOCK]; n /= BLOCK; }
    return s;
}

static unsigned digitsum(u128 n)
{
    const uint64_t e19 = 10000000000000000000ULL;
    if (n <= UINT64_MAX)
        return digitsum64((uint64_t)n);
    return digitsum(n / e19) + digitsum64((uint64_t)(n % e19));
}

/* ------------------------------------------------------------------ */
/* primality (gapsieve.c)                                              */
/* ------------------------------------------------------------------ */

struct mont { uint64_t n, ninv, r2, one; };

static uint64_t mont_mul(const struct mont *m, uint64_t a, uint64_t b)
{
    /* t + lo*n can need 129 bits when n > 2^63, so sum the halves: the low
     * 64 bits cancel by construction, leaving only their carry. */
    u128 t = (u128)a * b;
    uint64_t lo = (uint64_t)t * m->ninv;
    u128 mn = (u128)lo * m->n;
    u128 r = (t >> 64) + (mn >> 64) + ((uint64_t)t != 0);
    if (r >= m->n)
        r -= m->n;
    return (uint64_t)r;
}

static void mont_init(struct mont *m, uint64_t n)
{
    m->n = n;
    uint64_t inv = n;                 /* Newton: inverse of n mod 2^64 */
    for (int i = 0; i < 5; i++)
        inv *= 2 - n * inv;
    m->ninv = (uint64_t)0 - inv;      /* -n^-1 mod 2^64 */
    m->one = (0 - n) % n;             /* 2^64 mod n */
    u128 r2 = ((u128)m->one << 64) % n;
    m->r2 = (uint64_t)r2;             /* 2^128 mod n */
}

static int sprp(const struct mont *m, uint64_t a)
{
    uint64_t n = m->n;
    if (a >= n) {
        a %= n;
        if (!a)
            return 1;
    }
    uint64_t d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;

    uint64_t x = mont_mul(m, a, m->r2);          /* a in Montgomery form */
    uint64_t r = m->one, nm1 = mont_mul(m, n - 1, m->r2);
    for (uint64_t e = d; e; e >>= 1) {
        if (e & 1)
            r = mont_mul(m, r, x);
        x = mont_mul(m, x, x);
    }
    if (r == m->one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mont_mul(m, r, r);
        if (r == nm1)
            return 1;
    }
    return 0;
}

static const uint32_t small_primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23,
                                         29, 31, 37, 41 };

/* deterministic for all n < 2^64 (the 12-base bound is 3.18e23) */
static int is_prime64(uint64_t n)
{
    if (n < 2)
        return 0;
    for (int i = 0; i < 12; i++) {
        if (n % small_primes[i] == 0)
            return n == small_primes[i];
    }
    if (n < 41 * 41)
        return 1;
    struct mont m;
    mont_init(&m, n);
    for (int i = 0; i < 12; i++)
        if (!sprp(&m, small_primes[i]))
            return 0;
    return 1;
}

static inline uint64_t mmul(uint64_t a, uint64_t b, uint64_t n, uint64_t ninv)
{
    u128 t = (u128)a * b;
    uint64_t lo = (uint64_t)t * ninv;
    u128 mn = (u128)lo * n;
    u128 r = (t >> 64) + (mn >> 64) + ((uint64_t)t != 0);
    if (r >= n)
        r -= n;
    return (uint64_t)r;
}

/* base-2 strong probable-prime test for odd n > 2; "composite" is a proof */
static int sprp2(uint64_t n)
{
    uint64_t inv = n;
    for (int i = 0; i < 5; i++)
        inv *= 2 - n * inv;
    uint64_t ninv = (uint64_t)0 - inv;
    uint64_t one = ((uint64_t)0 - n) % n;        /* 2^64 mod n */
    uint64_t nm1 = n - one;                      /* Montgomery form of -1 */
    uint64_t d = n - 1;
    int s = __builtin_ctzll(d);
    d >>= s;

    int b = 63 - __builtin_clzll(d);
    uint64_t r = one >= n - one ? one - (n - one) : one + one;   /* 2 */
    while (--b >= 0) {
        r = mmul(r, r, n, ninv);
        if ((d >> b) & 1)
            r = r >= n - r ? r - (n - r) : r + r;
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

/* ---------- two-word arithmetic, for numbers from 2^64 up to 10^30 ---------- */

static inline int ctz128(u128 x)                /* x != 0 */
{
    uint64_t lo = (uint64_t)x;
    return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((uint64_t)(x >> 64));
}

static inline int top128(u128 x)                /* top set bit; x != 0 */
{
    uint64_t hi = (uint64_t)(x >> 64);
    return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((uint64_t)x);
}

/* x mod q for q < 2^32, without a 128-bit division */
static inline uint32_t mod_small(u128 x, uint32_t q)
{
    uint64_t r64 = ((uint64_t)0 - q) % q;       /* 2^64 mod q */
    return (uint32_t)(((uint64_t)(x >> 64) % q * r64 + (uint64_t)x % q) % q);
}

/* -n^-1 mod 2^64 for odd n, from its low word */
static inline uint64_t neg_inv64(u128 n)
{
    uint64_t n0 = (uint64_t)n, inv = n0;
    for (int i = 0; i < 5; i++)
        inv *= 2 - n0 * inv;
    return (uint64_t)0 - inv;
}

/*
 * Montgomery product a*b/2^128 mod n, for odd n < 2^127 and a, b < n: CIOS
 * over two 64-bit words, ninv = -n^-1 mod 2^64.  Each pass adds one word of
 * a*b, then a multiple of n that clears the low word, and drops that word.
 * Every sum fits in 128 bits, and the result is below n.
 */
static inline u128 mmul2(u128 a, u128 b, u128 n, uint64_t ninv)
{
    uint64_t a0 = (uint64_t)a, a1 = (uint64_t)(a >> 64);
    uint64_t b0 = (uint64_t)b, b1 = (uint64_t)(b >> 64);
    uint64_t n0 = (uint64_t)n, n1 = (uint64_t)(n >> 64);
    uint64_t t0, t1, t2, m;
    u128 c;

    c = (u128)a0 * b0;                          /* t = a*b0 */
    t0 = (uint64_t)c;
    c = (u128)a1 * b0 + (c >> 64);
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);
    m = t0 * ninv;                              /* t = (t + m*n) / 2^64 */
    c = ((u128)m * n0 + t0) >> 64;
    c += (u128)m * n1 + t1;
    t0 = (uint64_t)c;
    c = (c >> 64) + t2;
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);

    c = (u128)a0 * b1 + t0;                     /* t += a*b1 */
    t0 = (uint64_t)c;
    c = (u128)a1 * b1 + t1 + (c >> 64);
    t1 = (uint64_t)c;
    c = (c >> 64) + t2;
    t2 = (uint64_t)c;
    m = t0 * ninv;                              /* t = (t + m*n) / 2^64 */
    c = ((u128)m * n0 + t0) >> 64;
    c += (u128)m * n1 + t1;
    t0 = (uint64_t)c;
    c = (c >> 64) + t2;
    t1 = (uint64_t)c;
    t2 = (uint64_t)(c >> 64);

    u128 t = (u128)t1 << 64 | t0;               /* t < 2n */
    return t2 || t >= n ? t - n : t;
}

/*
 * 2^128 mod n, the Montgomery form of 1, for odd n < 2^82.  A 128-bit
 * remainder is a slow library loop, so above 2^64 the quotient comes from
 * floating point: 2^128/n in double precision is within 2^13 of the truth,
 * the leftover multiple of n is small enough to estimate the same way, and
 * a final step either way makes the remainder exact.
 */
static inline u128 mont_one2(u128 n)
{
    uint64_t n1 = (uint64_t)(n >> 64);
    if (!n1 || n1 >> 18)                        /* below 2^64, or from 2^82 on, */
        return ((u128)0 - n) % n;               /* where the estimate is too coarse */
    double dn = (double)n1 * 0x1p64 + (double)(uint64_t)n;
    double qd = 0x1p128 / dn;
    uint64_t q = qd < 0x1p64 ? (uint64_t)qd : UINT64_MAX;
    /* 2^128 - q*n, exact as a signed number: |it| < 2^13 n < 2^95 */
    __int128 r = (__int128)((u128)0 - (u128)q * n), sn = (__int128)n;
    double rd = (double)(int64_t)(r >> 64) * 0x1p64 + (double)(uint64_t)r;
    r -= (__int128)(int64_t)(rd / dn) * sn;
    while (r < 0)
        r += sn;
    while (r >= sn)
        r -= sn;
    return (u128)r;
}

/* a + b mod n, for a, b < n < 2^127 */
static inline u128 addmod2(u128 a, u128 b, u128 n)
{
    u128 y = n - b;
    return a >= y ? a - y : a + b;
}

/* sprp2 for odd n > 2 below 2^127, in two-word Montgomery form */
static int sprp2w(u128 n)
{
    uint64_t ninv = neg_inv64(n);
    u128 one = mont_one2(n);                    /* 2^128 mod n */
    u128 nm1 = n - one;
    u128 d = n - 1;
    int s = ctz128(d);
    d >>= s;

    int b = top128(d);
    u128 r = addmod2(one, one, n);              /* 2 */
    while (--b >= 0) {
        r = mmul2(r, r, n, ninv);
        if ((d >> b) & 1)
            r = addmod2(r, r, n);
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul2(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

/* strong probable-prime test to base a < n, n odd, two-word Montgomery */
static int sprpw(u128 n, uint64_t ninv, u128 one, uint32_t a)
{
    u128 nm1 = n - one, d = n - 1;
    int s = ctz128(d);
    d >>= s;

    u128 x = 0;                                 /* a*2^128 mod n */
    for (int b = 31 - __builtin_clz(a); b >= 0; b--) {
        x = addmod2(x, x, n);
        if ((a >> b) & 1)
            x = addmod2(x, one, n);
    }
    u128 r = one;
    for (u128 e = d; e; e >>= 1) {
        if (e & 1)
            r = mmul2(r, x, n, ninv);
        x = mmul2(x, x, n, ninv);
    }
    if (r == one || r == nm1)
        return 1;
    while (--s > 0) {
        r = mmul2(r, r, n, ninv);
        if (r == nm1)
            return 1;
    }
    return 0;
}

/* deterministic for 41 < n < PSI13: the 13 smallest prime bases */
static int is_prime_wide(u128 n)
{
    for (int i = 0; i < 13; i++)
        if (mod_small(n, small_primes[i]) == 0)
            return 0;
    uint64_t ninv = neg_inv64(n);
    u128 one = mont_one2(n);
    for (int i = 0; i < 13; i++)
        if (!sprpw(n, ninv, one, small_primes[i]))
            return 0;
    return 1;
}

/* exact below PSI13; above it a 13-base strong probable-prime test, which
 * has no known counterexample but is not a proof */
static int is_prime(u128 n)
{
    return n <= UINT64_MAX ? is_prime64((uint64_t)n) : is_prime_wide(n);
}

/* the screen: "composite" is a proof, "prime" is only probable */
static inline int sprp2_any(u128 n)
{
    return n <= UINT64_MAX ? sprp2((uint64_t)n) : sprp2w(n);
}

/* ------------------------------------------------------------------ */
/* chains, exactly as dschain follows them                             */
/* ------------------------------------------------------------------ */

/* length of the chain starting at prime p; fills chain[] if non-NULL */
static int chain_length(u128 p, u128 *chain, int cap)
{
    int n = 1;
    if (chain && cap > 0)
        chain[0] = p;
    for (;;) {
        unsigned d = digitsum(p);
        /* p is odd (p > 2), so an odd digit sum makes p + d even */
        if ((n > 1 || p != 2) && (d & 1))
            break;
        u128 q = p + d;
        if (!is_prime(q))
            break;
        p = q;
        if (chain && n < cap)
            chain[n] = p;
        n++;
    }
    return n;
}

/* is there no prime r with r + digitsum(r) == p ? */
static int is_chain_start(u128 p)
{
    for (unsigned d = 1; d <= MAX_DS && d < p; d++) {
        u128 r = p - d;
        if (digitsum(r) == d && is_prime(r))
            return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* numbers with a given digit sum, in increasing order                 */
/* ------------------------------------------------------------------ */

struct dsgen {
    u128 v;
    int8_t d[NDIG];                       /* d[0] is the units digit */
};

/* the smallest y >= lo with digit sum m; 0 if there is none in NDIG digits */
static int dsgen_first(struct dsgen *g, u128 lo, int m)
{
    int d[NDIG], tot = 0;
    for (int i = 0; i < NDIG; i++) {
        d[i] = (int)(lo % 10);
        lo /= 10;
        tot += d[i];
    }
    if (lo)                               /* lo has more than NDIG digits */
        return 0;
    if (tot != m) {
        /* y agrees with lo above position i and is larger at i: the lowest
         * such i that can still reach digit sum m gives the smallest y */
        int below = 0, i;                 /* digit sum of lo below position i */
        for (i = 0; i < NDIG; i++) {
            int above = tot - below - d[i];
            int v = d[i] + 1;
            if (v < m - above - 9 * i)    /* the lower digits hold at most 9i */
                v = m - above - 9 * i;
            int vmax = m - above < 9 ? m - above : 9;
            if (v <= vmax) {
                int r = m - above - v;
                d[i] = v;
                for (int j = 0; j < i; j++) {   /* smallest: ...b999 */
                    d[j] = r < 9 ? r : 9;
                    r -= d[j];
                }
                break;
            }
            below += d[i];
        }
        if (i == NDIG)
            return 0;
    }
    g->v = 0;
    for (int j = NDIG - 1; j >= 0; j--) {
        g->v = g->v * 10 + (u128)d[j];
        g->d[j] = (int8_t)d[j];
    }
    return 1;
}

/* the next larger number with the same digit sum: raise the lowest digit
 * that can rise, and rebuild everything below it as small as possible */
static inline int dsgen_next(struct dsgen *g)
{
    int low = g->d[0];                    /* digit sum below position i */
    u128 lowval = (u128)g->d[0];          /* and the value of those digits */
    for (int i = 1; i < NDIG; i++) {
        if (g->d[i] < 9 && low >= 1) {
            int r = low - 1, a = r / 9, b = r % 9;   /* below i: b, then a 9s */
            for (int j = 0; j < i; j++)
                g->d[j] = (int8_t)(j < a ? 9 : j == a ? b : 0);
            g->d[i]++;
            g->v = g->v - lowval + ((u128)(b + 1) * P10[a] - 1) + P10[i];
            return 1;
        }
        low += g->d[i];
        lowval += (u128)g->d[i] * P10[i];
    }
    return 0;
}

/* how many y in [0, x] have digit sum m */
static u128 ways[NDIG + 1][9 * NDIG + 1];      /* n-digit strings, sum s */

static void ways_init(void)
{
    ways[0][0] = 1;
    for (int n = 1; n <= NDIG; n++)
        for (int s = 0; s <= 9 * NDIG; s++)
            for (int v = 0; v <= 9 && v <= s; v++)
                ways[n][s] += ways[n - 1][s - v];
}

static u128 count_le(u128 x, int m)
{
    if (m < 0 || m > 9 * NDIG)
        return 0;
    int d[NDIG];
    for (int i = 0; i < NDIG; i++) {
        d[i] = (int)(x % 10);
        x /= 10;
    }
    u128 c = 0;
    int rem = m;
    for (int i = NDIG - 1; i >= 0 && rem >= 0; i--) {
        for (int v = 0; v < d[i] && v <= rem; v++)
            c += ways[i][rem - v];
        rem -= d[i];
    }
    return c + (rem == 0);
}

/* ... and of those, how many do not end in 9 (when no9) */
static u128 count_upto(u128 x, int m, int no9)
{
    u128 c = count_le(x, m);
    if (no9 && x >= 9)
        c -= count_le((x - 9) / 10, m - 9);   /* y = 10y' + 9 */
    return c;
}

static u128 count_range(u128 lo, u128 hi, int m, int no9)
{
    if (hi < lo)
        return 0;
    return count_upto(hi, m, no9) - (lo ? count_upto(lo - 1, m, no9) : 0);
}

/* the largest digit sum of any number in [0, x] */
static int max_ds_upto(u128 x)
{
    int best = (int)digitsum(x);
    for (int i = 1; i < NDIG && P10[i] <= x; i++) {
        /* x with the digit at 10^i lowered by one and everything below 9 */
        int s = (int)digitsum((x / P10[i]) * P10[i] - 1);
        if (s > best)
            best = s;
    }
    return best;
}

/* ------------------------------------------------------------------ */
/* patterns                                                            */
/* ------------------------------------------------------------------ */

/*
 * The sieve primes, 7..499, and the products of them the sieve works mod:
 * p is linear in the digits of Q, so p mod M is a constant plus the digits
 * times per-pattern weights.  Every product is below 2^23, so a sum of 32
 * weighted digits stays below 2^32 and the arithmetic stays in 32 bits (on
 * a GPU too).  Stage 1, 7..23, is kept up to date as the digits change and
 * tested for every candidate; stage 2, 29..61, and stage 3, 67..499, are
 * worked out from the digits only for the few that pass the stage before.
 * (Primes above 199 barely matter on the CPU, but they spare the GPU's
 * feeder thread four in five of the starts it would otherwise have to test.)
 */
static const uint32_t SP[] = {
    7, 11, 13, 17, 19, 23,                          /* stage 1 */
    29, 31, 37, 41, 43, 47, 53, 59, 61,             /* stage 2 */
    67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
    139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373,
    379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457,
    461, 463, 467, 479, 487, 491, 499 };
#define NSP ((int)(sizeof SP / sizeof SP[0]))
#define MASKW 8                           /* 64-bit words per mask: q < 512 */

#define NMOD   39
#define STAGE2 1                          /* MODS[1..3] */
#define STAGE3 4                          /* MODS[4..38] */
/* the primes of MODS[k] are SP[MOD_FIRST[k]] .. SP[MOD_FIRST[k+1]-1]:
 * six, then threes to 199, then pairs */
static const uint8_t MOD_FIRST[NMOD + 1] = {
    0, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 43,
    45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77,
    79, 81, 83, 85, 87, 89, 91, 92 };
static uint32_t MODS[NMOD];               /* filled in by mods_init */
#define M0 7436429u                       /* MODS[0] = 7*11*13*17*19*23 */

static void mods_init(void)
{
    for (int k = 0; k < NMOD; k++) {
        MODS[k] = 1;
        for (int i = MOD_FIRST[k]; i < MOD_FIRST[k + 1]; i++)
            MODS[k] *= SP[i];
    }
}

struct pattern {
    int S;                  /* ds(H) */
    int t;                  /* trailing 9s of H, or -1: the chain never carries */
    int m;                  /* digit sum of Q, the number enumerated */
    uint32_t x0;            /* p mod 10^4 */
    u128 p10t;              /* carrying: H = Q*10^t + 10^t - 1 */
    uint32_t off[MAXLEN];   /* term k is p + off[k] */
    int steps[MAXLEN];      /* the digit sum added at each step */
    uint64_t mask[NSP][MASKW]; /* bit r: p = r (mod SP[i]) puts a term on a multiple */
    uint32_t wt[NMOD][NDIG + 1]; /* digit j of Q adds d*wt[k][j] to p mod MODS[k] */
    uint32_t base[NMOD];    /* p mod MODS[k] when Q = 0 */
    uint32_t nine[NDIG + 1];/* 9*(wt[0][0] + ... + wt[0][a-1]) mod M0 */
    u128 qlo, qhi;          /* Q range for the search */
    u128 count;             /* candidates in the range */
};

static struct pattern *g_pat;
static int g_npat, g_nkilled;

/*
 * Walk the first L terms from low digits x0 under a high part with digit
 * sum S and t trailing 9s, in "unbounded" low coordinates (u >= BLOCK means
 * the chain has carried into H+1).  1 if every term is odd and prime to 3
 * and 5, with off[], steps[] and *carries filled in.
 */
static int walk(int S, int t, uint32_t x0, int L, uint32_t *off, int *steps,
                int *carries)
{
    if (!(x0 & 1) || (S + ds4[x0]) % 3 == 0)   /* p = ds(p) (mod 3), and each */
        return 0;                              /* step doubles p mod 3        */
    int S1 = S + 1 - 9 * t;
    uint32_t u = x0;
    *carries = 0;
    for (int k = 0; k < L; k++) {
        uint32_t low = u;
        int hs = S;
        if (u >= BLOCK) {
            low = u - BLOCK;
            hs = S1;
            *carries = 1;
        }
        if (low % 10 == 0 || low % 10 == 5)
            return 0;
        off[k] = u - x0;
        if (k == L - 1)
            break;
        int s = hs + ds4[low];
        if (s & 1)
            return 0;
        steps[k] = s;
        u += (uint32_t)s;
    }
    return 1;
}

/* the Q range and candidate count of a pattern within starts [plo, phi] */
static void pattern_range(struct pattern *P, u128 plo, u128 phi)
{
    P->count = 0;
    P->qlo = 1;
    P->qhi = 0;
    if (phi < plo || phi < P->x0)
        return;
    u128 hlo = plo <= P->x0 ? 0 : (plo - P->x0 + BLOCK - 1) / BLOCK;
    u128 hhi = (phi - P->x0) / BLOCK;
    if (hlo < 1)
        hlo = 1;                          /* H = 0 is the brute-force strip */
    if (hhi < hlo)
        return;
    if (P->t < 0) {
        P->qlo = hlo;
        P->qhi = hhi;
    } else {
        u128 c = P->p10t - 1;
        if (hhi < c)
            return;
        P->qlo = hlo <= c ? 0 : (hlo - c + P->p10t - 1) / P->p10t;
        P->qhi = (hhi - c) / P->p10t;
    }
    P->count = count_range(P->qlo, P->qhi, P->m, P->t >= 0);
}

static void add_pattern(int S, int t, uint32_t x0, int L, const uint32_t *off,
                        const int *steps)
{
    struct pattern P;
    memset(&P, 0, sizeof P);
    P.S = S;
    P.t = t;
    P.x0 = x0;
    P.m = t < 0 ? S : S - 9 * t;
    P.p10t = t < 0 ? 1 : P10[t];
    memcpy(P.off, off, (size_t)L * sizeof off[0]);
    memcpy(P.steps, steps, (size_t)(L - 1) * sizeof steps[0]);
    for (int i = 0; i < NSP; i++) {
        uint32_t q = SP[i], struck = 0;
        for (int k = 0; k < L; k++) {
            uint32_t x = (q - off[k] % q) % q;
            if (!(P.mask[i][x >> 6] >> (x & 63) & 1))
                struck++;
            P.mask[i][x >> 6] |= 1ULL << (x & 63);
        }
        if (struck == q) {                      /* every residue is struck: */
            g_nkilled++;                        /* some term is always a    */
            return;                             /* multiple of q            */
        }
    }
    /* p = Q*10^(t+4) + (10^t - 1)*10^4 + x0 when the chain carries (t >= 0),
     * p = Q*10^4 + x0 when it does not */
    int e = t < 0 ? 4 : t + 4;
    u128 c = t < 0 ? (u128)x0 : (P10[t] - 1) * BLOCK + x0;
    for (int k = 0; k < NMOD; k++) {
        uint64_t w = 1;
        for (int i = 0; i < e; i++)
            w = w * 10 % MODS[k];
        for (int j = 0; j <= NDIG; j++) {
            P.wt[k][j] = (uint32_t)w;
            w = w * 10 % MODS[k];
        }
        P.base[k] = mod_small(c, MODS[k]);
    }
    P.nine[0] = 0;
    for (int a = 1; a <= NDIG; a++)
        P.nine[a] = (uint32_t)((P.nine[a - 1] + 9ULL * P.wt[0][a - 1]) % M0);
    static int cap;
    if (g_npat == cap) {
        cap = cap ? 2 * cap : 256;
        g_pat = realloc(g_pat, (size_t)cap * sizeof *g_pat);
        if (!g_pat) {
            fprintf(stderr, "out of memory for patterns\n");
            exit(1);
        }
    }
    g_pat[g_npat++] = P;
}

/* every pattern for L that a high part H <= hmax can take */
static void classify(int L, u128 hmax)
{
    int smax = max_ds_upto(hmax);
    int tmax = 0;
    for (u128 h = hmax; h; h /= 10)
        tmax++;
    uint32_t off[MAXLEN];
    int steps[MAXLEN], carries;
    for (int S = 0; S <= smax; S++)
        for (int t = 0; t <= tmax && 9 * t <= S; t++)
            for (uint32_t x0 = 1; x0 < BLOCK; x0 += 2) {
                if (!walk(S, t, x0, L, off, steps, &carries))
                    continue;
                if (carries)
                    add_pattern(S, t, x0, L, off, steps);
                else if (t == 0)          /* independent of t: record once */
                    add_pattern(S, -1, x0, L, off, steps);
            }
}

/* ------------------------------------------------------------------ */
/* shared search state                                                 */
/* ------------------------------------------------------------------ */

/* an aligned block of Q: pre*10^r .. pre*10^r + 10^r - 1 */
struct item { int pat, r; u128 qlo, qhi, plo; };

static struct item *g_items;
static uint64_t g_nitems;
static uint64_t g_next;                   /* atomic: next item to claim */
static int g_active;                      /* atomic: workers still running */

static int g_len = 10, g_all = 0, g_first = 0, g_quiet = 0, g_gpu = 0;
static u128 g_start = 2, g_end = U128_MAX;

static uint64_t g_cand;                    /* atomic: candidates so far, live */
static uint64_t g_dcand, g_dsieved, g_dtests; /* atomic: in finished items only */
static uint64_t g_found;                   /* the rest: under g_lock */
static int g_longest;
static u128 g_longest_at;
static u128 g_best = U128_MAX;             /* --first: smallest hit so far */

static volatile sig_atomic_t g_stop = 0;
static volatile sig_atomic_t g_signo = 0;
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
static const char *CLR = "";

/* the item map for checkpoints: one bit per item, and a watermark below
 * which every item is done (items are claimed in order, so few lie above) */
static uint64_t *g_done;
static uint64_t g_wm;

static void handle_signal(int sig)
{
    if (g_stop) {
        static const char m[] = "\nsecond interrupt -- exiting now, "
                                "no checkpoint written\n";
        ssize_t r = write(2, m, sizeof m - 1);
        (void)r;
        _exit(130);
    }
    g_signo = sig;
    g_stop = 1;
    static const char m[] = "\ninterrupt -- stopping workers "
                            "(Ctrl-C again to abort now)\n";
    ssize_t r = write(2, m, sizeof m - 1);
    (void)r;
}

static uint64_t load64(const uint64_t *p) { return __atomic_load_n(p, __ATOMIC_RELAXED); }
static void add64(uint64_t *p, uint64_t v) { __atomic_fetch_add(p, v, __ATOMIC_RELAXED); }

static u128 best_so_far(void)
{
    pthread_mutex_lock(&g_lock);
    u128 b = g_best;
    pthread_mutex_unlock(&g_lock);
    return b;
}

static void mark_done(uint64_t i)
{
    pthread_mutex_lock(&g_lock);
    __atomic_fetch_or(&g_done[i >> 6], 1ULL << (i & 63), __ATOMIC_RELAXED);
    while (g_wm < g_nitems && (g_done[g_wm >> 6] >> (g_wm & 63) & 1))
        g_wm++;
    pthread_mutex_unlock(&g_lock);
}

/* a start p whose first g_len terms are proven prime: follow it and report */
static void report(u128 p)
{
    u128 ch[64];
    int n = chain_length(p, ch, 64);
    if (n < g_len)                        /* cannot happen */
        return;
    if (!g_all && !is_chain_start(p))
        return;
    int cap = n < 64 ? n : 64;
    pthread_mutex_lock(&g_lock);
    int show = 1;
    if (g_first) {
        if (p < g_best)
            g_best = p;
        else
            show = 0;                     /* an earlier hit is already known */
    }
    if (show) {
        printf("%s[len %d] ", CLR, n);
        for (int k = 0; k < cap; k++)
            printf("%s%s(+%u)", k ? " -> " : "", U128S(ch[k]), digitsum(ch[k]));
        if (cap < n)
            printf(" -> ...");
        printf("\n");
        if (ch[cap - 1] >= PSI13)
            printf("  (terms above 3.3e24 are 13-base strong probable primes, which is "
                   "not a proof: certify them with verify_chain.py)\n");
        fflush(stdout);
        g_found++;
        if (n > g_longest || (n == g_longest && p < g_longest_at)) {
            g_longest = n;
            g_longest_at = p;
        }
    }
    pthread_mutex_unlock(&g_lock);
}

static u128 block_count(u128 pre, int r, int m, int no9);

/* a candidate's start: Q = pre*10^r + (the r digits d) */
static u128 start_of(const struct pattern *P, u128 pre, int r, const int8_t *d)
{
    u128 q = 0;
    for (int j = r - 1; j >= 0; j--)
        q = q * 10 + (u128)d[j];
    q += pre * P10[r];
    u128 H = P->t < 0 ? q : q * P->p10t + (P->p10t - 1);
    return H * BLOCK + P->x0;
}

/* stage 1, primes below 64, with constant divisors */
#define STRUCK(P, i, q, x) (((P)->mask[i][0] >> ((x) % q##u)) & 1)

/* stages 2 and 3: MODS[k0..k1-1], worked out from the digits */
static inline int sieve_from(const struct pattern *P, const uint32_t *pref,
                             const int8_t *d, int r, int k0, int k1)
{
    for (int k = k0; k < k1; k++) {
        uint32_t x = pref[k];             /* below 289 * MODS[k] < 2^32 */
        for (int j = 0; j < r; j++)
            x += (uint32_t)d[j] * P->wt[k][j];
        x %= MODS[k];
        for (int i = MOD_FIRST[k]; i < MOD_FIRST[k + 1]; i++) {
            uint32_t y = x % SP[i];
            if (P->mask[i][y >> 6] >> (y & 63) & 1)
                return 0;
        }
    }
    return 1;
}

/* the next r-digit string with the same digit sum (dsgen_next on the free
 * digits), keeping R = p mod M0 up to date; 0 past the block */
static inline int succ(const struct pattern *P, int8_t *d, int r, uint32_t *R)
{
    int low = d[0];
    uint64_t old = (uint64_t)d[0] * P->wt[0][0];  /* the digits below i, weighted */
    for (int i = 1; i < r; i++) {
        if (d[i] < 9 && low >= 1) {
            int rr = low - 1, a = rr / 9, b = rr % 9;   /* below i: b, then a 9s */
            for (int j = 0; j < i; j++)
                d[j] = (int8_t)(j < a ? 9 : j == a ? b : 0);
            d[i]++;
            uint64_t x = (uint64_t)*R + P->wt[0][i] + P->nine[a] +
                         (uint64_t)b * P->wt[0][a] + 9ULL * NDIG * M0 - old;
            *R = (uint32_t)(x % M0);
            return 1;
        }
        low += d[i];
        old += (uint64_t)d[i] * P->wt[0][i];
    }
    return 0;
}

/* test hook: with DSPATTERN_SURVIVORS=file, every start that passes the
 * sieve is appended to file, so CPU and GPU runs can be compared exactly */
static FILE *g_survlog;

static void log_survivor(u128 p)
{
    if (!g_survlog)
        return;
    pthread_mutex_lock(&g_lock);
    fprintf(g_survlog, "%s\n", U128S(p));
    pthread_mutex_unlock(&g_lock);
}

/* the L terms of a start that passed the sieve: screen, then prove */
static int prove_terms(const struct pattern *P, u128 p, uint64_t *tests)
{
    for (int k = 0; k < g_len; k++) {
        (*tests)++;
        if (!sprp2_any(p + P->off[k]))
            return 0;
    }
    for (int k = 0; k < g_len; k++)
        if (!is_prime(p + P->off[k]))
            return 0;
    return 1;
}

/* what the search of an item needs, on the CPU or the GPU */
struct iparams {
    u128 pre;                             /* Q = pre*10^r + y */
    int r, need, no9;                     /* y: r digits, digit sum need, and no
                                             final 9 when the chain carries */
    uint64_t count;                       /* how many such y */
    uint32_t pref[NMOD];                  /* p mod MODS[k] for y = 0 */
};

static void item_params(const struct item *it, struct iparams *ip)
{
    const struct pattern *P = &g_pat[it->pat];
    ip->r = it->r;
    ip->pre = it->qlo / P10[ip->r];
    ip->need = P->m - (int)digitsum(ip->pre);
    ip->no9 = P->t >= 0 && ip->r > 0;
    ip->count = (uint64_t)block_count(ip->pre, ip->r, P->m, P->t >= 0);
    for (int k = 0; k < NMOD; k++)
        ip->pref[k] = (uint32_t)(((uint64_t)mod_small(ip->pre, MODS[k]) * P->wt[k][ip->r] +
                                  P->base[k]) % MODS[k]);
}

#ifdef DSPATTERN_GPU
/* r-digit strings with digit sum s, not ending in 9 when no9 (r >= 1) */
static u128 strings(int r, int s, int no9)
{
    if (s < 0)
        return 0;
    u128 c = ways[r][s];
    if (no9 && s >= 9)
        c -= ways[r - 1][s - 9];
    return c;
}

/* the digits of the k-th y (from 0, in increasing order) of an item */
static void unrank(const struct iparams *ip, uint64_t k, int8_t *d)
{
    int left = ip->need;
    for (int pos = ip->r - 1; pos >= 0; pos--) {
        int dig;
        for (dig = 0; dig <= 9 && dig <= left; dig++) {
            u128 c = pos ? strings(pos, left - dig, ip->no9)
                         : (u128)(left == dig && !(ip->no9 && dig == 9));
            if (k < c)
                break;
            k -= (uint64_t)c;
        }
        d[pos] = (int8_t)dig;
        left -= dig;
    }
}

/* a start that passed the whole sieve, as the k-th candidate of an item */
static void finish_survivor(const struct item *it, const struct iparams *ip,
                            uint64_t k, uint64_t *tests)
{
    const struct pattern *P = &g_pat[it->pat];
    int8_t d[NDIG] = { 0 };
    unrank(ip, k, d);
    u128 p = start_of(P, ip->pre, ip->r, d);
    log_survivor(p);
    if (prove_terms(P, p, tests))
        report(p);
}

#endif

/* an item's work is finished: count it and mark it done */
static void item_finished(uint64_t i, uint64_t cand, uint64_t sieved, uint64_t tests)
{
    add64(&g_dcand, cand);
    add64(&g_dsieved, sieved);
    add64(&g_dtests, tests);
    mark_done(i);
}

/*
 * One item: every Q = pre*10^r + y in the block whose digit sum is P->m (and,
 * for a carrying pattern, whose last digit is not 9), in increasing order.
 * The loop never forms p: it keeps p mod M0 as the digits change and sieves
 * on that, works out the later stages from the digits for the few that
 * pass, and builds p only for the rare starts that survive.  1 if the item
 * ran to its end.
 */
static int run_item(uint64_t idx)
{
    const struct item *it = &g_items[idx];
    const struct pattern *P = &g_pat[it->pat];
    struct iparams ip;
    item_params(it, &ip);
    int r = ip.r;
    uint64_t cand = 0, shown = 0, sieved = 0, tests = 0;
    int stopped = 0;
    if (!ip.count) {
        item_finished(idx, 0, 0, 0);
        return 1;
    }

    int8_t d[NDIG] = { 0 };               /* y, the smallest with digit sum need */
    uint64_t R = ip.pref[0];
    for (int j = 0, left = ip.need; j < r; j++) {
        d[j] = (int8_t)(left < 9 ? left : 9);
        left -= d[j];
        R += (uint64_t)d[j] * P->wt[0][j];
    }
    uint32_t R0 = (uint32_t)(R % M0);

    for (;;) {
        if (!(ip.no9 && d[0] == 9)) {     /* H needs exactly t trailing 9s */
            cand++;
            if (!(STRUCK(P, 0, 7, R0) || STRUCK(P, 1, 11, R0) || STRUCK(P, 2, 13, R0) ||
                  STRUCK(P, 3, 17, R0) || STRUCK(P, 4, 19, R0) || STRUCK(P, 5, 23, R0)) &&
                sieve_from(P, ip.pref, d, r, STAGE2, STAGE3) &&
                sieve_from(P, ip.pref, d, r, STAGE3, NMOD)) {
                sieved++;
                u128 p = start_of(P, ip.pre, r, d);
                log_survivor(p);
                if (prove_terms(P, p, &tests))
                    report(p);
            }
            if (cand == ip.count)
                break;
            if (!(cand & 4095)) {
                if (g_stop) {
                    stopped = 1;
                    break;
                }
                add64(&g_cand, cand - shown);
                shown = cand;
                if (g_first && start_of(P, ip.pre, r, d) > best_so_far())
                    break;                /* the rest is above a known hit */
            }
        }
        if (!succ(P, d, r, &R0))
            break;
    }
    add64(&g_cand, cand - shown);
    if (stopped)
        return 0;                         /* interrupted: redone on resume */
    item_finished(idx, cand, sieved, tests);
    return 1;
}

/* the next item to search, skipping those done before a resume and, under
 * --first, those whose smallest start is above a known hit; 0 at the end */
static int claim_item(uint64_t *out)
{
    for (;;) {
        if (g_stop)
            return 0;
        uint64_t i = __atomic_fetch_add(&g_next, 1, __ATOMIC_RELAXED);
        if (i >= g_nitems)
            return 0;
        if (__atomic_load_n(&g_done[i >> 6], __ATOMIC_RELAXED) >> (i & 63) & 1)
            continue;
        if (g_first && g_items[i].plo > best_so_far())
            continue;                     /* items are sorted by lowest p */
        *out = i;
        return 1;
    }
}

static void *worker(void *arg)
{
    (void)arg;
    uint64_t i;
    while (claim_item(&i))
        run_item(i);
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}

#ifdef DSPATTERN_GPU
/* the Metal GPU search, in dspattern_gpu.m, which includes this file */
static const char *gpu_setup(void);
static const char *gpu_name(void);
static void *gpu_worker(void *arg);
#endif

/* brute force for the strip below 10^4, where H = 0 */
static void direct_scan(uint64_t lo, uint64_t hi)
{
    for (uint64_t p = lo; p <= hi; p++)
        if (is_prime64(p) && chain_length(p, NULL, 0) >= g_len)
            report(p);
}

static int cmp_item(const void *a, const void *b)
{
    const struct item *x = a, *y = b;
    if (x->plo != y->plo)
        return x->plo < y->plo ? -1 : 1;
    if (x->pat != y->pat)
        return x->pat < y->pat ? -1 : 1;
    return x->qlo < y->qlo ? -1 : x->qlo > y->qlo;
}

static uint64_t g_items_cap;
static u128 g_target;                     /* candidates per item, at most */

static void push_item(int pat, u128 qlo, u128 qhi, int r)
{
    const struct pattern *P = &g_pat[pat];
    if (g_nitems == g_items_cap) {
        g_items_cap = g_items_cap ? 2 * g_items_cap : 1024;
        g_items = realloc(g_items, g_items_cap * sizeof *g_items);
        if (!g_items) {
            fprintf(stderr, "out of memory for work items\n");
            exit(1);
        }
    }
    u128 H = P->t < 0 ? qlo : qlo * P->p10t + (P->p10t - 1);
    g_items[g_nitems++] = (struct item){ pat, r, qlo, qhi, H * BLOCK + P->x0 };
}

/* candidates among the aligned block pre*10^r .. pre*10^r + 10^r - 1:
 * numbers with digit sum m (not ending in 9, when no9), straight from ways */
static u128 block_count(u128 pre, int r, int m, int no9)
{
    int need = m - (int)digitsum(pre);
    if (need < 0 || need > 9 * r)
        return 0;
    if (r == 0)
        return !(no9 && pre % 10 == 9);
    u128 c = ways[r][need];
    if (no9 && need >= 9)
        c -= ways[r - 1][need - 9];
    return c;
}

/* An aligned block becomes one item, or is cut by its next digit.  Blocks
 * are also cut until three leading digits are fixed, so that no item spans
 * more than 1% in p: a --first search then sweeps upward almost exactly in
 * order, instead of finishing wide items far above the eventual answer. */
static void split_block(int pat, u128 pre, int r)
{
    const struct pattern *P = &g_pat[pat];
    u128 c = block_count(pre, r, P->m, P->t >= 0);
    if (!c)
        return;
    if ((c <= g_target && pre >= 100) || r == 0) {
        push_item(pat, pre * P10[r], pre * P10[r] + (P10[r] - 1), r);
        return;
    }
    for (int d = 0; d <= 9; d++)
        split_block(pat, pre * 10 + d, r - 1);
}

/*
 * Cut each pattern's Q range into items of at most g_target candidates,
 * sorted by the smallest start each can hold.  The range is covered by
 * aligned decimal blocks (fixed leading digits, any trailing ones), whose
 * counts come straight from the ways table, so no search is needed.
 */
static void make_items(u128 total)
{
    /* items of up to 2^24 candidates (a tenth of a second of one thread),
     * more for huge ranges so the list stays near a million entries */
    g_target = (u128)1 << 24;
    if (total / ((u128)1 << 20) > g_target)
        g_target = total / ((u128)1 << 20);
    for (int i = 0; i < g_npat; i++) {
        const struct pattern *P = &g_pat[i];
        if (!P->count)
            continue;
        u128 lo = P->qlo, hi = P->qhi;
        for (;;) {                        /* the largest aligned block at lo */
            int r = 0;
            while (r + 1 < NDIG && lo % P10[r + 1] == 0 && hi - lo >= P10[r + 1] - 1)
                r++;
            split_block(i, lo / P10[r], r);
            if (hi - lo < P10[r])
                break;
            lo += P10[r];
        }
    }
    if (g_nitems)
        qsort(g_items, g_nitems, sizeof *g_items, cmp_item);
    g_done = calloc(g_nitems / 64 + 1, sizeof *g_done);
    if (!g_done) {
        fprintf(stderr, "out of memory for the item map\n");
        exit(1);
    }
}

/* FNV-1a over the item list: a resumed run must cut the same items */
static uint64_t items_hash(void)
{
    uint64_t h = 1469598103934665603ULL;
    for (uint64_t i = 0; i < g_nitems; i++) {
        const unsigned char *b = (const unsigned char *)&g_items[i].plo;
        for (size_t k = 0; k < sizeof(u128); k++)
            h = (h ^ b[k]) * 1099511628211ULL;
        b = (const unsigned char *)&g_items[i].qhi;
        for (size_t k = 0; k < sizeof(u128); k++)
            h = (h ^ b[k]) * 1099511628211ULL;
    }
    return h;
}

/* ------------------------------------------------------------------ */
/* checkpoints                                                         */
/* ------------------------------------------------------------------ */

#define CK_VERSION 1

struct ckpt {
    int loaded, version, complete, length, all, first;
    u128 start, end, best, longest_at;
    uint64_t nitems, hash, wm, cand, sieved, tests, found;
    int longest;
    double elapsed;
    uint64_t *done;
    size_t ndone;
};

/* Written under g_lock; items in flight count as unfinished and are redone.
 * Temp file + fsync + rename, so a crash never leaves a damaged file. */
static int ck_write(const char *path, int complete, double elapsed, uint64_t hash)
{
    char tmp[4096];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) {
        fprintf(stderr, "\nwarning: cannot write checkpoint %s: %s\n", tmp,
                strerror(errno));
        return -1;
    }
    fprintf(f, "# dspattern checkpoint -- safe to inspect\n");
    fprintf(f, "version %d\n", CK_VERSION);
    fprintf(f, "complete %d\n", complete);
    fprintf(f, "length %d\n", g_len);
    fprintf(f, "all %d\n", g_all);
    fprintf(f, "first %d\n", g_first);
    fprintf(f, "start %s\n", U128S(g_start));
    fprintf(f, "end %s\n", U128S(g_end));
    fprintf(f, "items %" PRIu64 " %016" PRIx64 "\n", g_nitems, hash);
    fprintf(f, "elapsed %.3f\n", elapsed);
    fprintf(f, "candidates %" PRIu64 "\n", load64(&g_dcand));
    fprintf(f, "sieved %" PRIu64 "\n", load64(&g_dsieved));
    fprintf(f, "tests %" PRIu64 "\n", load64(&g_dtests));
    fprintf(f, "found %" PRIu64 "\n", g_found);
    if (g_longest)
        fprintf(f, "longest %d %s\n", g_longest, U128S(g_longest_at));
    if (g_best != U128_MAX)
        fprintf(f, "best %s\n", U128S(g_best));
    fprintf(f, "watermark %" PRIu64 "\n", g_wm);
    for (uint64_t i = g_wm; i < g_nitems; i++)
        if (g_done[i >> 6] >> (i & 63) & 1)
            fprintf(f, "done %" PRIu64 "\n", i);
    fflush(f);
    fsync(fileno(f));
    if (fclose(f) != 0 || rename(tmp, path) != 0) {
        fprintf(stderr, "\nwarning: cannot finish checkpoint %s: %s\n", path,
                strerror(errno));
        return -1;
    }
    return 0;
}

static int ck_read(const char *path, struct ckpt *ck)
{
    memset(ck, 0, sizeof *ck);
    ck->best = U128_MAX;
    FILE *f = fopen(path, "r");
    if (!f)
        return 0;
    size_t cap = 0;
    char line[512], w[128];
    while (fgets(line, sizeof line, f)) {
        uint64_t a;
        int i;
        if (line[0] == '#' || line[0] == '\n')
            continue;
        if (sscanf(line, "version %d", &ck->version) == 1) ;
        else if (sscanf(line, "complete %d", &ck->complete) == 1) ;
        else if (sscanf(line, "length %d", &ck->length) == 1) ;
        else if (sscanf(line, "all %d", &ck->all) == 1) ;
        else if (sscanf(line, "first %d", &ck->first) == 1) ;
        else if (sscanf(line, "start %127s", w) == 1) parse_num(w, &ck->start);
        else if (sscanf(line, "end %127s", w) == 1) parse_num(w, &ck->end);
        else if (sscanf(line, "items %" SCNu64 " %" SCNx64, &ck->nitems, &ck->hash) == 2) ;
        else if (sscanf(line, "elapsed %lf", &ck->elapsed) == 1) ;
        else if (sscanf(line, "candidates %" SCNu64, &ck->cand) == 1) ;
        else if (sscanf(line, "sieved %" SCNu64, &ck->sieved) == 1) ;
        else if (sscanf(line, "tests %" SCNu64, &ck->tests) == 1) ;
        else if (sscanf(line, "found %" SCNu64, &ck->found) == 1) ;
        else if (sscanf(line, "longest %d %127s", &i, w) == 2) {
            ck->longest = i;
            parse_num(w, &ck->longest_at);
        } else if (sscanf(line, "best %127s", w) == 1) parse_num(w, &ck->best);
        else if (sscanf(line, "watermark %" SCNu64, &ck->wm) == 1) ;
        else if (sscanf(line, "done %" SCNu64, &a) == 1) {
            if (ck->ndone == cap) {
                cap = cap ? 2 * cap : 64;
                ck->done = realloc(ck->done, cap * sizeof *ck->done);
                if (!ck->done) {
                    fprintf(stderr, "out of memory reading checkpoint\n");
                    exit(1);
                }
            }
            ck->done[ck->ndone++] = a;
        }
    }
    fclose(f);
    ck->loaded = 1;
    return 1;
}

/* ------------------------------------------------------------------ */
/* output                                                              */
/* ------------------------------------------------------------------ */

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void fmt_duration(double s, char *buf, size_t n)
{
    if (s < 0 || !isfinite(s)) {
        snprintf(buf, n, "--:--:--");
        return;
    }
    if (s > 359999) {
        snprintf(buf, n, ">99h");
        return;
    }
    int t = (int)(s + 0.5);
    snprintf(buf, n, "%02d:%02d:%02d", t / 3600, (t / 60) % 60, t % 60);
}

static void usage(const char *argv0)
{
    fprintf(stderr,
"Exhaustive search for chains of primes where each term is the previous term\n"
"plus its own decimal digit sum, testing only numbers whose digits allow a\n"
"chain of the requested length (dschain scans every prime instead).\n"
"\n"
"Usage: %s [options]\n"
"  -l, --length N     minimum chain length to find          (default 10)\n"
"  -s, --start N      lowest starting prime to consider     (default 2)\n"
"  -e, --end N        highest starting prime, or max        (default max: 1e30,\n"
"                     less a margin; above 3.3e24, where Miller-Rabin is not a\n"
"                     proof, hits are reported as probable primes to certify)\n"
"  -t, --threads N    CPU worker threads   (default: all cores; with --gpu, all\n"
"                     but two, and 0 runs the GPU alone)\n"
"      --gpu          also search on the GPU (macOS builds, Metal)\n"
"  -a, --all          report mid-chain starts too (default: maximal chains only)\n"
"  -F, --first        report only the earliest chain; the search is exhaustive,\n"
"                     so it is the confirmed smallest start in the range\n"
"  -k, --checkpoint F checkpoint to file F, and resume from it if it exists\n"
"  -i, --interval S   seconds between checkpoint writes      (default 60)\n"
"      --restart      ignore an existing checkpoint and start over\n"
"  -p, --patterns     list the digit patterns for --length and exit\n"
"  -q, --quiet        no progress output\n"
"  -h, --help         this message\n"
"\n"
"Numbers are exact: 1e20, 2.5e21, 300T, 18_446_744_073_709_551_616, or max.\n"
"\n"
"Examples: %s -l 10 -s 1e19 -e 1e21 --first -k a10.ckpt\n"
"          %s -l 9 -s 2 -e 1e11 --first        (reproduces a(9) in a blink)\n"
"          %s -l 10 --patterns\n",
            argv0, argv0, argv0, argv0);
}

static void list_patterns(void)
{
    for (int i = 0; i < g_npat; i++) {
        const struct pattern *P = &g_pat[i];
        printf("  low digits %04u  ds(H)=%-3d ", P->x0, P->S);
        if (P->t < 0)
            printf("t=any ");
        else
            printf("t=%-3d ", P->t);
        printf(" ds(p)=%-3u steps", P->S + ds4[P->x0]);
        for (int k = 0; k < g_len - 1; k++)
            printf(" +%d", P->steps[k]);
        printf("\n      %s candidates in range", GROUPED(P->count));
        struct dsgen g;
        if (P->count && dsgen_first(&g, P->qlo, P->m)) {
            /* the first that really has t trailing 9s */
            while (P->t >= 0 && g.d[0] == 9 && dsgen_next(&g))
                ;
            u128 H = P->t < 0 ? g.v : g.v * P->p10t + (P->p10t - 1);
            printf(", the smallest %s", U128S(H * BLOCK + P->x0));
        }
        printf("\n");
    }
}

static void print_status(double t0, double prior, u128 total, int live)
{
    double el = now_sec() - t0;
    uint64_t c = load64(&g_cand);
    double frac = total ? (double)c / (double)total : 1;
    double rate = el + prior > 0 ? (double)c / (el + prior) : 0;
    char ebuf[32], rbuf[32];
    fmt_duration(el + prior, ebuf, sizeof ebuf);
    fmt_duration(frac > 1e-9 && frac < 1 ? (el + prior) * (1 / frac - 1) : -1,
                 rbuf, sizeof rbuf);
    pthread_mutex_lock(&g_lock);
    u128 below = g_wm < g_nitems ? g_items[g_wm].plo : g_end;
    char hb[64];
    if (g_first)
        snprintf(hb, sizeof hb, "%s", g_best == U128_MAX ? "no hit yet" : "hit!");
    else
        snprintf(hb, sizeof hb, "hits=%" PRIu64, g_found);
    printf("%s%6.2f%%  %.3gM cand/s  complete below %.4Le  elapsed %s  eta %s  %s%s",
           CLR, 100 * frac, rate / 1e6, (long double)below, ebuf, rbuf, hb,
           live ? "" : "\n");
    fflush(stdout);
    pthread_mutex_unlock(&g_lock);
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    int threads = -1, list = 0, restart = 0;
    const char *ck_path = NULL;
    double ck_interval = 60;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        int has_next = i + 1 < argc;
#define NEXT() (has_next ? argv[++i] : (usage(argv[0]), exit(2), ""))
        if (!strcmp(a, "-l") || !strcmp(a, "--length"))
            g_len = (int)num_arg(NEXT(), a);
        else if (!strcmp(a, "-s") || !strcmp(a, "--start"))
            g_start = num_arg(NEXT(), a);
        else if (!strcmp(a, "-e") || !strcmp(a, "--end"))
            g_end = num_arg(NEXT(), a);
        else if (!strcmp(a, "-t") || !strcmp(a, "--threads"))
            threads = (int)num_arg(NEXT(), a);
        else if (!strcmp(a, "-a") || !strcmp(a, "--all"))
            g_all = 1;
        else if (!strcmp(a, "-F") || !strcmp(a, "--first"))
            g_first = 1;
        else if (!strcmp(a, "-k") || !strcmp(a, "--checkpoint"))
            ck_path = NEXT();
        else if (!strcmp(a, "-i") || !strcmp(a, "--interval"))
            ck_interval = (double)num_arg(NEXT(), a);
        else if (!strcmp(a, "--restart"))
            restart = 1;
        else if (!strcmp(a, "-p") || !strcmp(a, "--patterns"))
            list = 1;
        else if (!strcmp(a, "--gpu"))
            g_gpu = 1;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet"))
            g_quiet = 1;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "unknown option: %s\n\n", a);
            usage(argv[0]);
            return 2;
        }
#undef NEXT
    }
    if (g_len < 4 || g_len > MAXLEN) {
        fprintf(stderr, "error: --length must be between 4 and %d "
                        "(dschain is the tool for shorter chains)\n", MAXLEN);
        return 2;
    }
    if (g_start < 2)
        g_start = 2;
    if (g_end > PAT_TOP) {
        if (g_end != U128_MAX)
            fprintf(stderr, "note: --end lowered to %s, the top of the search\n",
                    U128S(PAT_TOP));
        g_end = PAT_TOP;
    }
    if (g_end < g_start) {
        fprintf(stderr, "error: --end must be >= --start\n");
        return 2;
    }
    if (threads == 0 && !g_gpu) {
        fprintf(stderr, "error: --threads 0 needs --gpu\n");
        return 2;
    }
    if (threads < 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = n > 0 ? (int)n : 1;
        if (g_gpu)                        /* leave room for the GPU's feeder */
            threads -= 2;
        if (threads < 1)
            threads = 1;
    }
#ifndef DSPATTERN_GPU
    if (g_gpu) {
        fprintf(stderr, "error: this build has no GPU support; on a Mac, build it "
                        "with make\n");
        return 2;
    }
#endif
    if (ck_interval < 1)
        ck_interval = 1;
    if (isatty(1))
        CLR = "\r\033[K";

    tables_init();
    ways_init();
    mods_init();
    if (getenv("DSPATTERN_SURVIVORS") && !(g_survlog = fopen(getenv("DSPATTERN_SURVIVORS"), "w"))) {
        fprintf(stderr, "error: cannot write %s\n", getenv("DSPATTERN_SURVIVORS"));
        return 2;
    }

    /* patterns cover starts in [plo, g_end]; brute force does below 10^4 */
    u128 plo = g_start > BLOCK ? g_start : BLOCK;
    int have_pat = plo <= g_end;
    classify(g_len, have_pat ? g_end / BLOCK : 0);
    u128 total = 0;
    int ncarry = 0;
    for (int i = 0; i < g_npat; i++) {
        if (have_pat)
            pattern_range(&g_pat[i], plo, g_end);
        total += g_pat[i].count;
        ncarry += g_pat[i].t >= 0;
    }

    printf("dspattern: p -> p + digitsum(p), chains of length >= %d\n", g_len);
    printf("  range     : %s .. %s\n", GROUPED(g_start), GROUPED(g_end));
    printf("  patterns  : %d (%d of them carry past the low four digits)", g_npat,
           ncarry);
    if (g_nkilled)
        printf("; %d more ruled out by a prime <= 61", g_nkilled);
    printf("\n");
    printf("  candidates: %s numbers in range fit a pattern", GROUPED(total));
    if (total)
        printf(" (1 in %.3Lg)", (long double)(g_end - g_start + 1) / (long double)total);
    printf("\n");
    if (list) {
        printf("\n");
        list_patterns();
        return 0;
    }

    make_items(total);
    uint64_t hash = items_hash();

    /* ---- resume ---- */
    struct ckpt ck;
    double prior = 0;
    memset(&ck, 0, sizeof ck);
    if (ck_path && !restart && ck_read(ck_path, &ck)) {
        if (ck.version != CK_VERSION) {
            fprintf(stderr, "error: checkpoint %s is version %d; this binary writes "
                            "version %d.  Use --restart to start over.\n",
                    ck_path, ck.version, CK_VERSION);
            return 2;
        }
        if (ck.length != g_len || ck.all != g_all || ck.first != g_first ||
            ck.start != g_start ||
            ck.end != g_end || ck.nitems != g_nitems || ck.hash != hash) {
            fprintf(stderr, "error: checkpoint %s was written for a different search:\n"
                            "         checkpoint: --length %d --start %s --end %s%s%s\n"
                            "         this run  : --length %d --start %s --end %s%s%s\n"
                            "       (or by a build that cuts the work differently).\n"
                            "       Use a different -k file, or --restart to overwrite it.\n",
                    ck_path, ck.length, U128S(ck.start), U128S(ck.end),
                    ck.all ? " --all" : "", ck.first ? " --first" : "", g_len,
                    U128S(g_start), U128S(g_end), g_all ? " --all" : "",
                    g_first ? " --first" : "");
            return 2;
        }
        for (uint64_t i = 0; i < ck.wm && i < g_nitems; i++)
            g_done[i >> 6] |= 1ULL << (i & 63);
        for (size_t i = 0; i < ck.ndone; i++)
            if (ck.done[i] < g_nitems)
                g_done[ck.done[i] >> 6] |= 1ULL << (ck.done[i] & 63);
        while (g_wm < g_nitems && (g_done[g_wm >> 6] >> (g_wm & 63) & 1))
            g_wm++;
        g_cand = g_dcand = ck.cand;
        g_dsieved = ck.sieved;
        g_dtests = ck.tests;
        g_found = ck.found;
        g_longest = ck.longest;
        g_longest_at = ck.longest_at;
        if (g_first)
            g_best = ck.best;
        prior = ck.elapsed;
        free(ck.done);
        if (ck.complete) {
            printf("\ncheckpoint %s says this search is already complete:\n", ck_path);
            if (g_first && g_best != U128_MAX)
                printf("  earliest chain of length >= %d starts at %s (confirmed)\n",
                       g_len, U128S(g_best));
            else
                printf("  chains >= %d found: %" PRIu64 "\n", g_len, g_found);
            printf("Use --restart to run it again.\n");
            return 0;
        }
    }

#ifdef DSPATTERN_GPU
    if (g_gpu) {
        const char *why = gpu_setup();
        if (why) {
            fprintf(stderr, "error: cannot use the GPU: %s\n", why);
            return 1;
        }
        printf("  threads   : GPU (%s) plus %d CPU thread%s\n", gpu_name(), threads,
               threads == 1 ? "" : "s");
    } else
#endif
    printf("  threads   : %d\n", threads);
    printf("  mode      : %s%s\n", g_all ? "all starts" : "maximal chains only",
           g_first ? ", earliest only" : "");
    if (ck_path) {
        printf("  checkpoint: %s every %.0fs\n", ck_path, ck_interval);
        if (ck.loaded)
            printf("  resuming  : %" PRIu64 " of %" PRIu64 " work items done, "
                   "%.0fs of prior runtime%s\n", g_wm, g_nitems, prior,
                   g_best != U128_MAX ? ", with a hit recorded" : "");
    }
    printf("\n");
    fflush(stdout);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = handle_signal;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    double t0 = now_sec();
    if (g_start < BLOCK)                  /* the strip where H = 0 */
        direct_scan((uint64_t)g_start, g_end < BLOCK - 1 ? (uint64_t)g_end : BLOCK - 1);

    int nt = threads + (g_gpu ? 1 : 0);   /* the GPU is fed by one more */
    pthread_t *tid = malloc((size_t)(nt ? nt : 1) * sizeof *tid);
    if (!tid) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    g_active = nt;
    for (int i = 0; i < nt; i++) {
        void *(*fn)(void *) = worker;
#ifdef DSPATTERN_GPU
        if (g_gpu && i == 0)
            fn = gpu_worker;
#endif
        pthread_create(&tid[i], NULL, fn, NULL);
    }

    /* the main thread reports progress and writes checkpoints until the
     * workers are done, or have stopped for a signal */
    int live = !g_quiet && isatty(1);
    double last_line = now_sec(), last_ck = last_line;
    for (;;) {
        struct timespec ts = { 0, 250000000 };
        nanosleep(&ts, NULL);
        if (!__atomic_load_n(&g_active, __ATOMIC_ACQUIRE))
            break;
        double t = now_sec();
        if (live)
            print_status(t0, prior, total, 1);
        else if (!g_quiet && t - last_line >= 60) {
            print_status(t0, prior, total, 0);
            last_line = t;
        }
        if (ck_path && !g_stop && t - last_ck >= ck_interval) {
            pthread_mutex_lock(&g_lock);
            ck_write(ck_path, 0, t - t0 + prior, hash);
            pthread_mutex_unlock(&g_lock);
            last_ck = now_sec();
        }
    }
    for (int i = 0; i < nt; i++)
        pthread_join(tid[i], NULL);
    free(tid);

    double el = now_sec() - t0 + prior;
    int hit = g_first && g_best != U128_MAX;
    int complete = !g_stop;               /* --first: everything below the hit */
    if (ck_path)
        ck_write(ck_path, complete, el, hash);
    if (live)
        printf("%s", CLR);

    printf("\n");
    if (g_stop) {
        printf("stopped by signal %d -- the search is incomplete\n", (int)g_signo);
        if (hit)
            printf("earliest hit so far starts at %s -- NOT confirmed as the smallest\n",
                   U128S(g_best));
    } else if (hit) {
        printf("earliest chain of length >= %d starts at %s\n"
               "(confirmed smallest: every candidate below it in range was tested)\n",
               g_len, U128S(g_best));
        if (g_best + (u128)MAXLEN * MAX_DS >= PSI13)
            printf("(its terms reach 3.3e24, beyond proven Miller-Rabin: certify them "
                   "with verify_chain.py)\n");
    } else if (!g_found) {
        printf("no chain of length >= %d starts in this range (exhaustive)\n", g_len);
    } else {
        printf("done\n");
    }
    printf("  elapsed          : %.2fs%s\n", el, prior > 0 ? " (including earlier runs)" : "");
    printf("  candidates       : %s%s\n", GROUPED(load64(&g_dcand)),
           g_stop ? " (in finished work items)" : "");
    printf("  after sieve      : %s\n", GROUPED(load64(&g_dsieved)));
    printf("  base-2 tests     : %s\n", GROUPED(load64(&g_dtests)));
    printf("  chains >= %-6d : %" PRIu64 "\n", g_len, g_found);
    if (g_found)
        printf("  longest found    : %d (starting at %s)\n", g_longest,
               U128S(g_longest_at));
    if (ck_path)
        printf("\ncheckpoint %s %s\n", ck_path,
               complete ? "marked complete" : "written -- rerun the same command to continue");
    return g_stop ? 130 : 0;
}
