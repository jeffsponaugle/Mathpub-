/* emirpx.c -- verify that a number is prime and an emirp in bases 2..N
 *
 * An emirp in base b: the number is prime, and the value obtained by
 * reversing its base-b digit string is a *different* prime.
 *
 * 128-bit arithmetic throughout (Montgomery multiplication, no 256-bit
 * types, no inline asm) -- builds clean on x86-64 and arm64.
 *
 *   cc -O2 -Wall -Wextra -o emirpx emirpx.c
 *
 * Usage:  emirpx <number> <max_base> [options]
 *         emirpx - <max_base> [options]     (numbers on stdin, one per line)
 *
 * Exit:   0 = prime and emirp in every base 2..max_base
 *         1 = failed somewhere
 *         2 = usage / parse / range error
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

typedef unsigned __int128 u128;
typedef unsigned long long u64;

/* Everything must stay below 2^127 so that modular doubling and the
 * reversal accumulator cannot wrap. */
#define U128_TOP  ((u128)1 << 127)
#define U128_MAXV (U128_TOP - 1)

/* Miller-Rabin with the first 13 primes is deterministic below this. */
#define MR_DET_LIMIT (((u128)3317044ULL * (u128)1000000000000000000ULL) \
                      + (u128)64679887385961981ULL)

/* ------------------------------------------------------------------ */
/* 128x128 -> 256 multiply, built from 64-bit limbs                    */
/* ------------------------------------------------------------------ */

static inline void mul128(u128 a, u128 b, u128 *hi, u128 *lo)
{
    u64 a0 = (u64)a, a1 = (u64)(a >> 64);
    u64 b0 = (u64)b, b1 = (u64)(b >> 64);

    u128 p00 = (u128)a0 * b0;
    u128 p01 = (u128)a0 * b1;
    u128 p10 = (u128)a1 * b0;
    u128 p11 = (u128)a1 * b1;

    u128 l = p00;
    u128 h = p11;
    u128 t;

    t = (u128)(u64)p01 << 64;
    l += t;  h += (l < t);
    h += (p01 >> 64);

    t = (u128)(u64)p10 << 64;
    l += t;  h += (l < t);
    h += (p10 >> 64);

    *hi = h;
    *lo = l;
}

/* ------------------------------------------------------------------ */
/* Montgomery arithmetic mod odd n < 2^127                             */
/* ------------------------------------------------------------------ */

typedef struct {
    u128 n;      /* modulus (odd)            */
    u128 ninv;   /* -n^-1 mod 2^128          */
    u128 r1;     /* 2^128 mod n              */
    u128 r2;     /* 2^256 mod n              */
} mont_t;

static void mont_init(mont_t *m, u128 n)
{
    m->n = n;

    u128 inv = n;                       /* correct mod 2^3 for odd n */
    for (int i = 0; i < 6; i++)         /* 3 -> 6 -> ... -> 192 bits */
        inv *= (u128)2 - n * inv;
    m->ninv = (u128)0 - inv;

    m->r1 = ((u128)0 - n) % n;          /* (2^128 - n) mod n == 2^128 mod n */

    u128 x = m->r1;                     /* double 128 times => 2^256 mod n */
    for (int i = 0; i < 128; i++) {
        x <<= 1;
        if (x >= n) x -= n;
    }
    m->r2 = x;
}

static inline u128 mont_mul(const mont_t *m, u128 a, u128 b)
{
    u128 hi, lo, mnhi, mnlo;

    mul128(a, b, &hi, &lo);
    mul128(lo * m->ninv, m->n, &mnhi, &mnlo);

    u128 t = lo + mnlo;
    unsigned carry = (t < lo);

    u128 s = hi + mnhi;
    unsigned ovf = (s < hi);
    u128 s2 = s + carry;
    ovf |= (s2 < s);

    if (ovf || s2 >= m->n) s2 -= m->n;
    return s2;
}

static inline u128 to_mont(const mont_t *m, u128 x)   { return mont_mul(m, x, m->r2); }

static u128 mont_pow(const mont_t *m, u128 base, u128 e)
{
    u128 r = m->r1;             /* 1 in Montgomery form */
    while (e) {
        if (e & 1) r = mont_mul(m, r, base);
        base = mont_mul(m, base, base);
        e >>= 1;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Primality                                                           */
/* ------------------------------------------------------------------ */

static const unsigned SMALL_P[] = {
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,
     43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101,
    103, 107, 109, 113
};
#define NSMALL (sizeof(SMALL_P)/sizeof(SMALL_P[0]))
#define NDET   13              /* bases 2..41: deterministic < 3.317e24 */

static int mr_witness(const mont_t *m, u128 a, u128 d, int s)
{
    u128 n  = m->n;
    u128 n1 = n - 1;

    a %= n;
    if (a == 0) return 0;

    u128 x = mont_pow(m, to_mont(m, a), d);
    u128 one = m->r1, minus1 = to_mont(m, n1);

    if (x == one || x == minus1) return 0;
    for (int i = 1; i < s; i++) {
        x = mont_mul(m, x, x);
        if (x == minus1) return 0;
        if (x == one)    return 1;
    }
    return 1;                    /* composite */
}

/* returns 1 = prime, 0 = composite; *proven = 0 when only strong-probable */
static int is_prime(u128 n, int *proven)
{
    if (proven) *proven = 1;
    if (n < 2) return 0;

    for (size_t i = 0; i < NSMALL; i++) {
        if (n % SMALL_P[i] == 0) return n == SMALL_P[i];
        if ((u128)SMALL_P[i] * SMALL_P[i] > n) return 1;
    }

    mont_t m;
    mont_init(&m, n);

    u128 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }

    size_t nb = NDET;
    if (n >= MR_DET_LIMIT) {            /* past the proven range */
        nb = NSMALL;
        if (proven) *proven = 0;
    }

    for (size_t i = 0; i < nb; i++)
        if (mr_witness(&m, SMALL_P[i], d, s)) return 0;

    return 1;
}

/* ------------------------------------------------------------------ */
/* Base conversion                                                     */
/* ------------------------------------------------------------------ */

#define MAXDIG 128

/* digits[0] is least significant; returns count */
static int to_digits(u128 n, unsigned b, unsigned *digits)
{
    int k = 0;
    if (n == 0) { digits[k++] = 0; return k; }
    while (n) {
        digits[k++] = (unsigned)(n % b);
        n /= b;
    }
    return k;
}

/* reverse the base-b digit string; 0 on overflow past 2^127-1 */
static int reverse_base(u128 n, unsigned b, u128 *out)
{
    u128 rev = 0;
    while (n) {
        unsigned d = (unsigned)(n % b);
        n /= b;
        if (rev > (U128_MAXV - d) / b) return 0;
        rev = rev * b + d;
    }
    *out = rev;
    return 1;
}

static const char DSYM[] = "0123456789abcdefghijklmnopqrstuvwxyz";

/* render digits (most significant first) into buf */
static void render(const unsigned *digits, int k, unsigned b, char *buf, size_t sz)
{
    size_t o = 0;
    if (b <= 36) {
        for (int i = k - 1; i >= 0 && o + 2 < sz; i--)
            buf[o++] = DSYM[digits[i]];
    } else {
        for (int i = k - 1; i >= 0 && o + 8 < sz; i--)
            o += (size_t)snprintf(buf + o, sz - o, "%s%u",
                                  (i == k - 1) ? "" : ".", digits[i]);
    }
    buf[o] = '\0';
}

static char *u128_dec(u128 v, char *buf, size_t sz)
{
    char tmp[48];
    int i = 0;
    if (v == 0) tmp[i++] = '0';
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    size_t o = 0;
    while (i > 0 && o + 1 < sz) buf[o++] = tmp[--i];
    buf[o] = '\0';
    return buf;
}

/* ------------------------------------------------------------------ */
/* Input parsing                                                       */
/* ------------------------------------------------------------------ */

/* accepts decimal, or 0x-prefixed hex; ',' and '_' ignored */
static int parse_u128(const char *s, u128 *out)
{
    while (isspace((unsigned char)*s)) s++;
    unsigned base = 10;
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) { base = 16; s += 2; }

    u128 v = 0;
    int any = 0;
    for (; *s; s++) {
        if (*s == ',' || *s == '_' || *s == '\'') continue;
        if (isspace((unsigned char)*s)) break;
        unsigned d;
        if (isdigit((unsigned char)*s))       d = (unsigned)(*s - '0');
        else if (base == 16 && isxdigit((unsigned char)*s))
                                              d = (unsigned)(tolower(*s) - 'a' + 10);
        else return 0;
        if (d >= base) return 0;
        if (v > (U128_MAXV - d) / base) return 0;   /* overflow / too big */
        v = v * base + d;
        any = 1;
    }
    if (!any) return 0;
    *out = v;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Checking                                                            */
/* ------------------------------------------------------------------ */

typedef struct {
    int quiet;          /* -q: verdict only                         */
    int nodigits;       /* -d: values only, no digit strings        */
    int stop;           /* -s: stop at first failing base           */
    int allow_pal;      /* --allow-palindrome                       */
} opts_t;

static int check_number(u128 n, unsigned maxbase, const opts_t *o)
{
    char db[512], rb[512], v1[48], v2[48];
    unsigned digits[MAXDIG], rdig[MAXDIG];
    int proven;

    if (n < 2 || n >= U128_TOP) {
        printf("%s: out of range (need 2 <= n < 2^127)\n", u128_dec(n, v1, sizeof v1));
        return 0;
    }

    int prime = is_prime(n, &proven);
    printf("n = %s%s\n", u128_dec(n, v1, sizeof v1),
           prime ? (proven ? "  [prime]" : "  [probable prime]") : "  [COMPOSITE]");
    if (!prime) return 0;

    int ok = 1;
    if (!o->quiet) {
        if (o->nodigits) printf("\n%5s  %-40s %s\n", "base", "reversal (decimal)", "result");
        else             printf("\n%5s  %-30s %-30s %s\n", "base", "digits", "reversed", "result");
    }

    for (unsigned b = 2; b <= maxbase; b++) {
        u128 rev;
        if (!reverse_base(n, b, &rev)) {
            printf("%5u  reversal overflows 128 bits -- cannot check\n", b);
            ok = 0;
            if (o->stop) break;
            continue;
        }

        int pal   = (rev == n);
        int rprime = 0, rproven = 1;
        if (!pal || o->allow_pal) rprime = is_prime(rev, &rproven);

        int good = rprime && (o->allow_pal || !pal);
        if (!good) ok = 0;

        if (!o->quiet) {
            const char *verdict = good ? (rproven ? "emirp" : "emirp (prp)")
                                       : (pal ? "FAIL palindrome" : "FAIL reversal composite");
            if (o->nodigits) {
                printf("%5u  %-40s %s\n", b, u128_dec(rev, v2, sizeof v2), verdict);
            } else {
                int k  = to_digits(n,   b, digits);
                int rk = to_digits(rev, b, rdig);
                render(digits, k, b, db, sizeof db);
                render(rdig, rk, b, rb, sizeof rb);
                printf("%5u  %-30s %-30s %s\n", b, db, rb, verdict);
            }
        }
        if (!good && o->stop) break;
    }

    printf("\n%s is %san emirp in all bases 2..%u\n",
           u128_dec(n, v1, sizeof v1), ok ? "" : "NOT ", maxbase);
    return ok;
}

static void usage(const char *p)
{
    fprintf(stderr,
        "usage: %s <number|-> <max_base> [options]\n"
        "  <number>   decimal or 0xhex, up to 2^127-1 ('-' reads numbers from stdin)\n"
        "  <max_base> check bases 2 through this value\n"
        "options:\n"
        "  -q                   verdict only, no per-base table\n"
        "  -d                   omit digit strings, print reversal values only\n"
        "  -s                   stop at the first failing base\n"
        "  --allow-palindrome   accept palindromic representations as emirps\n"
        "exit: 0 = emirp in all bases, 1 = not, 2 = usage error\n", p);
}

int main(int argc, char **argv)
{
    if (argc < 3) { usage(argv[0]); return 2; }

    opts_t o = {0, 0, 0, 0};
    const char *numarg = argv[1];

    char *end;
    unsigned long mb = strtoul(argv[2], &end, 10);
    if (*end || mb < 2 || mb > 100000) {
        fprintf(stderr, "bad max_base '%s' (need 2..100000)\n", argv[2]);
        return 2;
    }
    unsigned maxbase = (unsigned)mb;

    for (int i = 3; i < argc; i++) {
        if      (!strcmp(argv[i], "-q")) o.quiet = 1;
        else if (!strcmp(argv[i], "-d")) o.nodigits = 1;
        else if (!strcmp(argv[i], "-s")) o.stop = 1;
        else if (!strcmp(argv[i], "--allow-palindrome")) o.allow_pal = 1;
        else { fprintf(stderr, "unknown option '%s'\n", argv[i]); usage(argv[0]); return 2; }
    }

    int allok = 1;

    if (!strcmp(numarg, "-")) {
        char line[600];
        while (fgets(line, sizeof line, stdin)) {
            char *p = line;
            while (isspace((unsigned char)*p)) p++;
            if (*p == '\0' || *p == '#') continue;
            u128 n;
            if (!parse_u128(p, &n)) {
                fprintf(stderr, "skipping unparseable line: %s", line);
                allok = 0;
                continue;
            }
            if (!check_number(n, maxbase, &o)) allok = 0;
            printf("\n");
        }
    } else {
        u128 n;
        if (!parse_u128(numarg, &n)) {
            fprintf(stderr, "bad number '%s' (must be < 2^127)\n", numarg);
            return 2;
        }
        allok = check_number(n, maxbase, &o);
    }

    return allok ? 0 : 1;
}
