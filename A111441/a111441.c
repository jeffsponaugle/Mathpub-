/* A111441: numbers k such that the sum of the squares of the first k primes
 * is divisible by k.
 *
 * Single-threaded reference version using primesieve.
 *
 * The running sum S = sum_{i<=n} p_i^2 is kept as a 192-bit integer
 * (lo: 128 bits, hi: 64 bits).  S exceeds 2^128 only around n ~ 3e12,
 * but the wide accumulator makes the code valid to n ~ 1e17.
 *
 * No term is divisible by 2 or 3 (see OEIS comments), so only
 * n == 1 or 5 (mod 6) is tested.
 *
 * Usage: ./a111441 [max_n]
 *
 * Build: cc -O3 -o a111441 a111441.c -lprimesieve
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <primesieve.h>

typedef unsigned __int128 u128;

typedef struct {
    u128     lo;
    uint64_t hi;   /* S = hi * 2^128 + lo */
} u192;

static inline void u192_add_u128(u192 *s, u128 x)
{
    s->lo += x;
    if (s->lo < x)
        s->hi++;
}

/* S mod n for 64-bit n */
static inline uint64_t u192_mod(const u192 *s, uint64_t n)
{
    uint64_t rlo = (uint64_t)(s->lo % n);
    if (s->hi == 0)
        return rlo;
    /* 2^64 mod n, then 2^128 mod n */
    uint64_t t64  = (uint64_t)((((u128)(UINT64_MAX % n)) + 1) % n);
    uint64_t t128 = (uint64_t)(((u128)t64 * t64) % n);
    u128 r = ((u128)(s->hi % n) * t128) % n;
    return (uint64_t)((r + rlo) % n);
}

int main(int argc, char **argv)
{
    uint64_t max_n = (argc > 1) ? strtoull(argv[1], NULL, 10) : 1000000000ULL;

    primesieve_iterator it;
    primesieve_init(&it);

    u192 sum = {0, 0};
    uint64_t n = 0;
    time_t t0 = time(NULL);

    while (n < max_n) {
        uint64_t p = primesieve_next_prime(&it);
        n++;
        u192_add_u128(&sum, (u128)p * p);

        uint64_t m = n % 6;
        if ((m == 1 || m == 5) && u192_mod(&sum, n) == 0)
            printf("HIT  n = %llu  (p_n = %llu, %lds elapsed)\n",
                   (unsigned long long)n, (unsigned long long)p,
                   (long)(time(NULL) - t0));
    }

    /* n = 1 is a term but 1 % 6 == 1, so it is caught above. */
    fprintf(stderr, "done: checked n <= %llu in %lds (last prime %llu)\n",
            (unsigned long long)max_n, (long)(time(NULL) - t0),
            (unsigned long long)primesieve_next_prime(&it));
    primesieve_free_iterator(&it);
    return 0;
}
