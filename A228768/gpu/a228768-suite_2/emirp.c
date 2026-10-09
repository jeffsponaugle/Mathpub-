/*
 * emirp.c -- verify that a number is an emirp in every base from 2 to N.
 *
 * An emirp in base b is a prime whose base-b digit reversal is also prime.
 * By the strict definition the reversal must differ from the original, i.e.
 * base-b palindromic primes do not count; pass -p to accept them.
 *
 * Build:  cc -O2 -o emirp emirp.c
 * Usage:  ./emirp [-q] [-p] <number> [max_base]
 *
 * Primality is deterministic Miller-Rabin over the first 12 prime bases,
 * which is proven correct for all n < 3.3e24, so it covers all of uint64.
 * (primesieve is great for generating primes in bulk but has no single-value
 *  primality test, so it isn't needed here.)
 *
 * Exit status: 0 = emirp in all requested bases, 1 = not, 2 = usage/range error.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <inttypes.h>

typedef unsigned __int128 u128;

#define MAX_DIGITS 64           /* base 2 is the worst case for uint64_t */
#define MAX_BASE   65536

/* ---------------------------------------------------------------- primality */

static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m)
{
    return (uint64_t)((u128)a * b % m);
}

static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m)
{
    uint64_t r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static int miller_witness(uint64_t a, uint64_t n, uint64_t d, int s)
{
    uint64_t x = powmod(a, d, n);
    if (x == 1 || x == n - 1) return 0;          /* not a witness */
    for (int i = 1; i < s; i++) {
        x = mulmod(x, x, n);
        if (x == n - 1) return 0;
    }
    return 1;                                     /* composite */
}

static int is_prime(uint64_t n)
{
    static const uint64_t small[] = {2,3,5,7,11,13,17,19,23,29,31,37};

    if (n < 2) return 0;
    for (size_t i = 0; i < sizeof small / sizeof *small; i++) {
        if (n % small[i] == 0) return n == small[i];
    }

    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }

    for (size_t i = 0; i < sizeof small / sizeof *small; i++) {
        if (miller_witness(small[i], n, d, s)) return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------- digits */

/* Fills d[] least-significant first; returns the digit count. */
static unsigned to_digits(uint64_t n, unsigned b, uint32_t *d)
{
    unsigned len = 0;
    if (n == 0) { d[len++] = 0; return len; }
    while (n) { d[len++] = (uint32_t)(n % b); n /= b; }
    return len;
}

/* Interprets d[] (least-significant first) as most-significant first, i.e.
 * the digit reversal.  Returns 0 on uint64_t overflow. */
static int reversed_value(const uint32_t *d, unsigned len, unsigned b, uint64_t *out)
{
    u128 v = 0;
    for (unsigned i = 0; i < len; i++) {
        v = v * b + d[i];
        if (v > UINT64_MAX) return 0;
    }
    *out = (uint64_t)v;
    return 1;
}

/* Renders digits into buf. rev=0 prints most-significant first (the normal
 * representation); rev=1 prints the reversal.  Bases above 36 use bracketed
 * decimal digits, e.g. [12][0][7].  Returns the string length. */
static size_t render(const uint32_t *d, unsigned len, unsigned b, int rev,
                     char *buf, size_t cap)
{
    static const char sym[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    size_t p = 0;
    for (unsigned i = 0; i < len; i++) {
        uint32_t dig = rev ? d[i] : d[len - 1 - i];
        if (b <= 36) {
            if (p + 1 < cap) buf[p] = sym[dig];
            p++;
        } else {
            int k = snprintf(buf + (p < cap ? p : cap - 1),
                             p < cap ? cap - p : 1, "[%" PRIu32 "]", dig);
            p += (size_t)k;
        }
    }
    if (p < cap) buf[p] = '\0'; else buf[cap - 1] = '\0';
    return p;
}

/* --------------------------------------------------------------------- main */

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-q] [-p] <number> [max_base]\n"
        "  <number>    the number to test (2 .. 2^64-1)\n"
        "  [max_base]  test bases 2..max_base (default 36)\n"
        "  -q          summary only, no per-base table\n"
        "  -p          accept base-b palindromes as emirps\n", prog);
}

int main(int argc, char **argv)
{
    int quiet = 0, allow_pal = 0;
    int argi = 1;

    while (argi < argc && argv[argi][0] == '-' && argv[argi][1]) {
        for (const char *c = argv[argi] + 1; *c; c++) {
            if (*c == 'q') quiet = 1;
            else if (*c == 'p') allow_pal = 1;
            else { usage(argv[0]); return 2; }
        }
        argi++;
    }
    if (argi >= argc) { usage(argv[0]); return 2; }

    errno = 0;
    char *end;
    uint64_t n = strtoull(argv[argi], &end, 0);
    if (errno || *end || end == argv[argi] || n < 2) {
        fprintf(stderr, "error: bad number '%s'\n", argv[argi]);
        return 2;
    }
    argi++;

    unsigned long max_base = 36;
    if (argi < argc) {
        errno = 0;
        max_base = strtoul(argv[argi], &end, 0);
        if (errno || *end || max_base < 2 || max_base > MAX_BASE) {
            fprintf(stderr, "error: max_base must be 2..%d\n", MAX_BASE);
            return 2;
        }
    }

    printf("n = %" PRIu64 "  (bases 2..%lu, palindromes %s)\n",
           n, max_base, allow_pal ? "allowed" : "rejected");

    if (!is_prime(n)) {
        printf("\n%" PRIu64 " is NOT prime -- it cannot be an emirp in any base.\n", n);
        return 1;
    }
    printf("%" PRIu64 " is prime.\n", n);

    /* Column width: base 2 always produces the longest string. */
    uint32_t d[MAX_DIGITS];
    unsigned wlen = to_digits(n, 2, d);
    int w = (int)wlen;
    if (w < 12) w = 12;

    if (!quiet) {
        printf("\n%5s  %-*s  %-*s  %-20s  %s\n",
               "base", w, "n in base", w, "reversed", "reversed value", "result");
    }

    char buf_fwd[MAX_DIGITS * 8], buf_rev[MAX_DIGITS * 8];
    int all_ok = 1;
    unsigned long first_fail = 0;

    for (unsigned long b = 2; b <= max_base; b++) {
        unsigned len = to_digits(n, (unsigned)b, d);
        uint64_t r;
        int ok, pal = 0;
        const char *verdict;

        render(d, len, (unsigned)b, 0, buf_fwd, sizeof buf_fwd);
        render(d, len, (unsigned)b, 1, buf_rev, sizeof buf_rev);

        if (!reversed_value(d, len, (unsigned)b, &r)) {
            ok = 0;
            verdict = "OVERFLOW (reversal exceeds 2^64)";
            r = 0;
        } else {
            pal = (r == n);
            if (pal && !allow_pal) {
                ok = 0;
                verdict = "FAIL (palindrome)";
            } else if (!is_prime(r)) {
                ok = 0;
                verdict = "FAIL (reversal composite)";
            } else {
                ok = 1;
                verdict = pal ? "ok (prime palindrome)" : "ok (both prime)";
            }
        }

        if (!ok && all_ok) { all_ok = 0; first_fail = b; }

        if (!quiet) {
            printf("%5lu  %-*s  %-*s  %-20" PRIu64 "  %s\n",
                   b, w, buf_fwd, w, buf_rev, r, verdict);
        }
    }

    if (n <= max_base) {
        printf("\nnote: for bases > %" PRIu64 " the number is a single digit,\n"
               "      so its reversal is itself (a palindrome).\n", n);
    }

    if (all_ok) {
        printf("\nRESULT: %" PRIu64 " is an emirp in every base 2..%lu.\n", n, max_base);
        return 0;
    }
    printf("\nRESULT: %" PRIu64 " is NOT an emirp in all bases 2..%lu"
           " (first failure at base %lu).\n", n, max_base, first_fail);
    return 1;
}
