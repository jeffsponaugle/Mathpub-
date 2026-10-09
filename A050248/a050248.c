/*
 * a050248.c — explore OEIS A050248: integer averages of the first k primes.
 *
 * A term appears whenever k divides S(k) = p_1 + p_2 + ... + p_k; the term
 * is the average S(k)/k.  Related sequences: k values are A045345, the
 * sums S(k) are A050247.
 *
 * Known terms (16): 2, 38, 110, 3066, 60020, 740282, 2340038, 29380602,
 * 957565746, 31043311588, 569424748566, 7207204117608, 10871205353578,
 * 196523412770096, 2665506690112870, 122498079071529726
 * a(17) > 125237452139872271 (Paul W. Dyson, Sep 2022), i.e. the search
 * frontier is beyond p ~ 2.5e17 — don't expect new terms on a laptop,
 * but everything up to ~1e12 reproduces in minutes.
 *
 * The running sum S(k) overflows 64 bits near p ~ 6e9, so the accumulator
 * is unsigned __int128 (good to ~3.4e38 — sum of all primes below ~1e19).
 *
 * Build:   make            (uses primesieve)
 * Usage:   ./a050248 [limit]        sieve primes p <= limit (default 1e10)
 *          ./a050248 -q [limit]     quiet: suppress progress lines
 */

#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <primesieve.h>

typedef unsigned __int128 u128;

/* Format a u128 as decimal into buf (must hold >= 40 bytes). */
static const char *u128_str(u128 v, char *buf)
{
    char *p = buf + 40;
    *--p = '\0';
    do {
        *--p = '0' + (char)(v % 10);
        v /= 10;
    } while (v);
    return p;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    uint64_t limit = 10000000000ULL; /* 1e10 */
    bool quiet = false;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-q") == 0) {
            quiet = true;
        } else {
            char *end;
            double d = strtod(argv[i], &end); /* accept 1e12 style */
            if (*end != '\0' || d < 2) {
                fprintf(stderr, "usage: %s [-q] [prime_limit]\n", argv[0]);
                return 1;
            }
            limit = (uint64_t)d;
        }
    }

    primesieve_iterator it;
    primesieve_init(&it);

    u128 sum = 0;              /* S(k) = sum of first k primes  */
    uint64_t k = 0;            /* how many primes so far        */
    uint64_t p;                /* current prime                 */
    unsigned nterms = 0;
    char b1[40], b2[40];

    const uint64_t progress_every = 250000000ULL; /* primes between reports */
    uint64_t next_report = progress_every;
    double t0 = now_sec();

    printf("A050248: integer averages of the first k primes, p <= %" PRIu64 "\n\n", limit);
    printf("%4s %20s %20s %26s %26s\n", "n", "k (A045345)", "p_k", "S(k) (A050247)", "a(n) = S(k)/k");

    while ((p = primesieve_next_prime(&it)) <= limit) {
        sum += p;
        k++;
        if (sum % k == 0) {
            nterms++;
            printf("%4u %20" PRIu64 " %20" PRIu64 " %26s %26s\n",
                   nterms, k, p, u128_str(sum, b1), u128_str(sum / k, b2));
            fflush(stdout);
        }
        if (!quiet && k == next_report) {
            double dt = now_sec() - t0;
            fprintf(stderr, "  ... k = %" PRIu64 ", p = %" PRIu64
                    ", avg = %s  [%.1fs, %.1fM primes/s]\n",
                    k, p, u128_str(sum / k, b1), dt, (double)k / dt / 1e6);
            next_report += progress_every;
        }
    }

    primesieve_free_iterator(&it);

    double dt = now_sec() - t0;
    printf("\nDone: %u term(s) among the first %" PRIu64 " primes (p <= %" PRIu64 ") in %.1fs.\n",
           nterms, k, limit, dt);
    printf("Final S(k) = %s, average = %s.\n",
           u128_str(sum, b1), u128_str(sum / k, b2));
    return 0;
}
