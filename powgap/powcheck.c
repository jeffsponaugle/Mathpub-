/*
 * powcheck.c -- verify a powgap chain by hand.
 *
 * Given a starting number p, walk the pattern p, p+2, p+6, p+14, p+30, ...
 * (offset k is 2^(k+1)-2, so the gaps are 2, 4, 8, 16, ...), printing each
 * member, whether it is prime, and whether any other prime interrupts the
 * gap, until the chain breaks.  Companion checker for powgap / OEIS A090807.
 *
 * Numbers accept the same formats as powgap: commas and underscores are
 * ignored (so powgap output can be pasted directly), and K M B/G T P E
 * suffixes and scientific notation work: 1997, 3,637,803,390,827, 1.5T.
 *
 * Build:  gcc -O3 -march=native -o powcheck powcheck.c
 *         clang -O3 -mcpu=native -o powcheck powcheck.c   (Apple silicon)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#define MAXLEN 24

/* ------------------------------------------------------------------ */
/* number parsing and formatting (same as powgap)                       */
/* ------------------------------------------------------------------ */

static int parse_qty(const char *s, uint64_t *out)
{
    char buf[160];
    size_t n = 0;
    for (const char *p = s; *p && n < sizeof buf - 1; ++p)
        if (*p != ',' && *p != '_' && !isspace((unsigned char)*p)) buf[n++] = *p;
    buf[n] = 0;
    if (!n) return -1;

    unsigned long long mult = 1;
    size_t len = strlen(buf);
    switch (toupper((unsigned char)buf[len - 1])) {
        case 'K': mult = 1000ULL;                   break;
        case 'M': mult = 1000000ULL;                break;
        case 'G': case 'B': mult = 1000000000ULL;   break;
        case 'T': mult = 1000000000000ULL;          break;
        case 'P': mult = 1000000000000000ULL;       break;
        case 'E': mult = 1000000000000000000ULL;    break;
        default:  mult = 0;
    }
    if (mult) buf[--len] = 0; else mult = 1;
    if (!len) return -1;

    char *endp;
    if (strpbrk(buf, ".eE")) {                 /* fractional / scientific */
        errno = 0;
        long double v = strtold(buf, &endp);
        if (*endp || errno || v < 0) return -1;
        long double r = v * (long double)mult;
        if (r >= 1.8446744073709551615e19L) return -1;
        *out = (uint64_t)r;
    } else {
        errno = 0;
        unsigned long long v = strtoull(buf, &endp, 10);
        if (*endp || errno) return -1;
        unsigned __int128 r = (unsigned __int128)v * mult;
        if (r > UINT64_MAX) return -1;
        *out = (uint64_t)r;
    }
    return 0;
}

static void fmt_u64(uint64_t v, char *out)      /* 1234567 -> 1,234,567 */
{
    char t[24]; int n = snprintf(t, sizeof t, "%llu", (unsigned long long)v);
    int o = 0;
    for (int i = 0; i < n; i++) {
        if (i && (n - i) % 3 == 0) out[o++] = ',';
        out[o++] = t[i];
    }
    out[o] = 0;
}

/* ------------------------------------------------------------------ */
/* primality: deterministic Miller-Rabin for the full 64-bit range      */
/* ------------------------------------------------------------------ */

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m)
{ return (uint64_t)((__uint128_t)a * b % m); }

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}

static const uint32_t TRIALP[] = {
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
    101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
    193,197,199,211,223,227,229,233,239,241,251
};
#define NTRIALP (sizeof TRIALP / sizeof TRIALP[0])

static int is_prime_u64(uint64_t n)
{
    if (n < 2) return 0;
    for (size_t i = 0; i < NTRIALP; i++) {
        uint32_t p = TRIALP[i];
        if ((uint64_t)p * p > n) return 1;
        if (n % p == 0) return n == p;
    }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const uint64_t base[] = {2,325,9375,28178,450775,9780504,1795265022};
    for (int i = 0; i < 7; i++) {
        uint64_t a = base[i] % n;
        if (!a) continue;
        uint64_t x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { ok = 1; break; }
        }
        if (!ok) return 0;
    }
    return 1;
}

/* small factor for a friendlier "composite" message; 0 if none found */
static uint64_t small_factor(uint64_t n)
{
    for (size_t i = 0; i < NTRIALP; i++)
        if (n % TRIALP[i] == 0 && n != TRIALP[i]) return TRIALP[i];
    for (uint64_t d = 257; d < 1000000 && d * d <= n; d += 2)
        if (n % d == 0) return d;
    return 0;
}

/* ------------------------------------------------------------------ */

static void check(uint64_t p)
{
    char a[32], b[32];
    fmt_u64(p, a);
    printf("p = %s\n", a);

    int len = 0, consec = 0, consec_broken = 0;
    uint64_t prev = 0;
    for (int k = 0; k <= MAXLEN; k++) {
        uint64_t off = (1ULL << (k + 1)) - 2;
        if (off > UINT64_MAX - p) {
            printf("  (next offset exceeds 2^64 -- stopping)\n");
            break;
        }
        uint64_t m = p + off;
        fmt_u64(m, a);
        printf("  %2d:  %*s%s", k + 1, (int)(21 - strlen(a)), "", a);
        if (k) printf("   gap %3llu = 2^%d", (unsigned long long)(1ULL << k), k);
        else   printf("                ");

        if (!is_prime_u64(m)) {
            uint64_t f = m < 2 ? 0 : small_factor(m);
            if (m < 2)      printf("   not prime -- chain ends\n");
            else if (f) {
                fmt_u64(m / f, b);
                printf("   composite (%llu x %s) -- chain ends\n",
                       (unsigned long long)f, b);
            } else
                printf("   composite -- chain ends\n");
            break;
        }
        printf("   prime");
        len++;

        if (k) {
            int interposed = 0;
            for (uint64_t x = prev + 2; x < m; x += 2) {
                if (is_prime_u64(x)) {
                    fmt_u64(x, b);
                    printf("%s%s", interposed ? ", " : "   [intervening prime: ",
                           b);
                    interposed = 1;
                }
            }
            if (interposed) { printf("]"); consec_broken = 1; }
        }
        if (!consec_broken) consec = len;
        printf("\n");
        prev = m;
        if (k == MAXLEN) printf("  (pattern limit of %d members reached)\n", MAXLEN + 1);
    }

    if (len == 0)
        printf("  -> %s is not prime: not a chain\n", a);
    else
        printf("  -> chain length %d (consecutive-prime prefix %d)"
               "  [powgap: len %d (consec %d)]\n", len, consec, len, consec);
}

static void usage(const char *p)
{
    fprintf(stderr,
"usage: %s <p> [<p> ...]\n"
"\n"
"  For each starting number p, walks the powgap pattern p, p+2, p+6, p+14,\n"
"  ... (gaps 2, 4, 8, 16, ...) printing every member, its primality, and any\n"
"  intervening primes, until the chain breaks.\n"
"\n"
"  Numbers may contain commas/underscores and powgap-style suffixes:\n"
"  1997   3,637,803,390,827   1.5T   6_824_897\n", p);
}

int main(int argc, char **argv)
{
    if (argc < 2) { usage(argv[0]); return 1; }
    int bad = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(argv[0]); return 0;
        }
        uint64_t p;
        if (parse_qty(argv[i], &p)) {
            fprintf(stderr, "cannot parse number: %s\n", argv[i]);
            bad = 1;
            continue;
        }
        if (i > 1) printf("\n");
        check(p);
    }
    return bad;
}
