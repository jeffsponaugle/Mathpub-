/* u128.h -- 128-bit modular arithmetic and deterministic primality.
 *
 * Montgomery multiplication with R = 2^128. Valid for odd moduli n < 2^126.
 * Miller-Rabin with the first 13 primes as bases is DETERMINISTIC for
 * n < 3,317,044,064,679,887,385,961,981 (~3.3e24), which covers the whole
 * CPAP-7 target range (< 1e23) with room to spare.
 */
#ifndef U128_H
#define U128_H

#include <stdint.h>
#include <stdio.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

/* 128 x 128 -> 256 */
static inline void mul256(u128 a, u128 b, u128 *hi, u128 *lo)
{
    u64 a0 = (u64)a, a1 = (u64)(a >> 64);
    u64 b0 = (u64)b, b1 = (u64)(b >> 64);
    u128 p00 = (u128)a0 * b0;
    u128 p01 = (u128)a0 * b1;
    u128 p10 = (u128)a1 * b0;
    u128 p11 = (u128)a1 * b1;
    u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
    *lo = ((u128)(u64)mid << 64) | (u64)p00;
    *hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
}

typedef struct { u128 n, ninv, r1, r2; } mont_t;

/* -n^-1 mod 2^128 by Newton iteration */
static inline u128 mont_ninv(u128 n)
{
    u128 x = 1;
    for (int i = 0; i < 7; i++) x *= 2 - n * x;   /* x -> n^-1 mod 2^(2^i) */
    return ~x + 1;                                /* negate */
}

static inline u128 redc(u128 hi, u128 lo, u128 n, u128 ninv)
{
    u128 m = lo * ninv;
    u128 mh, ml;
    mul256(m, n, &mh, &ml);
    u128 carry = (lo + ml) < lo;      /* lo + ml is 0 mod 2^128 */
    u128 t = hi + mh + carry;
    if (t >= n) t -= n;
    return t;
}

static inline void mont_init(mont_t *m, u128 n)
{
    m->n = n;
    m->ninv = mont_ninv(n);
    m->r1 = ((u128)0 - n) % n;                 /* 2^128 mod n */
    u128 r = m->r1;                            /* r2 = 2^256 mod n */
    for (int i = 0; i < 128; i++) { r <<= 1; if (r >= n || r < m->r1) r %= n; r %= n; }
    m->r2 = r;
}

static inline u128 mont_mul(const mont_t *m, u128 a, u128 b)
{
    u128 hi, lo;
    mul256(a, b, &hi, &lo);
    return redc(hi, lo, m->n, m->ninv);
}

static inline u128 to_mont(const mont_t *m, u128 a)   { return mont_mul(m, a % m->n, m->r2); }
static inline u128 from_mont(const mont_t *m, u128 a) { return redc(0, a, m->n, m->ninv); }

static inline u128 mont_pow(const mont_t *m, u128 a, u128 e)
{
    u128 r = m->r1;                            /* 1 in Montgomery form */
    while (e) {
        if (e & 1) r = mont_mul(m, r, a);
        a = mont_mul(m, a, a);
        e >>= 1;
    }
    return r;
}

static const uint32_t MR_BASES[13] =
    {2,3,5,7,11,13,17,19,23,29,31,37,41};

/* small primes for trial division; also used to reject fast */
static const uint32_t SMALLP[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
    193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
    293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401
};
#define NSMALLP (sizeof SMALLP / sizeof SMALLP[0])

static inline int is_prime_u128(u128 n)
{
    if (n < 2) return 0;
    for (size_t i = 0; i < NSMALLP; i++) {
        u128 p = SMALLP[i];
        if (p * p > n) return 1;
        if (n % p == 0) return n == p;
    }
    mont_t m;
    mont_init(&m, n);
    u128 d = n - 1;
    int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 13; i++) {
        u128 a = to_mont(&m, MR_BASES[i]);
        u128 x = mont_pow(&m, a, d);
        if (x == m.r1 || x == m.n - m.r1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) {
            x = mont_mul(&m, x, x);
            if (x == m.n - m.r1) { ok = 1; break; }
        }
        if (!ok) return 0;
    }
    return 1;
}

/* decimal I/O */
static inline int u128_str(u128 v, char *out)
{
    char t[48]; int n = 0;
    if (!v) t[n++] = '0';
    while (v) { t[n++] = '0' + (int)(v % 10); v /= 10; }
    for (int i = 0; i < n; i++) out[i] = t[n - 1 - i];
    out[n] = 0;
    return n;
}

static inline void u128_print(FILE *f, u128 v)
{
    char b[48]; u128_str(v, b); fputs(b, f);
}

#endif
