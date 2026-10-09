/* primechain.c
 *
 * Search for chains of primes under the recurrence
 *
 *     p_{n+1} = p_n * digitsum(p_n) + 1
 *
 * where every term must itself be prime. Prints every chain whose length
 * is >= the requested target length.
 *
 * Usage:
 *     ./primechain <target_len> [lo] [hi]
 *
 *     target_len : minimum chain length to report (>= 1)
 *     lo, hi     : range of starting primes to scan (default 2 .. 100000000)
 *
 * Example:
 *     ./primechain 6 2 20000000
 *
 * Notes:
 *   - Primality uses deterministic Miller-Rabin (the 12 base witnesses
 *     {2..37} are provably exact for all n < 3.3e24, covering all 64 bits).
 *   - Modular multiply uses unsigned __int128 to avoid overflow.
 *   - Terms can exceed 2^64 for long chains (roughly length >= 8). The
 *     program detects this and warns rather than reporting wrong results;
 *     to go further you need a bignum library such as GMP.
 *
 * Build:  cc -O2 -o primechain primechain.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

/* (a * b) mod m, overflow-safe via 128-bit intermediate */
static inline u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((u128)a * b % m);
}

/* (a ^ e) mod m */
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

/* Deterministic Miller-Rabin, exact for all 64-bit n. */
static int is_prime(u64 n) {
    if (n < 2) return 0;
    static const u64 sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n % sp[i] == 0) return n == sp[i];
    }
    u64 d = n - 1;
    int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }
    for (int i = 0; i < 12; i++) {
        u64 a = sp[i];
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int j = 0; j < r - 1; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

static inline int digit_sum(u64 n) {
    int s = 0;
    while (n) { s += (int)(n % 10); n /= 10; }
    return s;
}

/* Walk the chain from a (assumed prime). Fills chain[] with the terms and
 * returns the length. *overflowed is set if we stopped only because the next
 * term would exceed 2^64 (so the true length may be longer). cap is the
 * size of chain[]. */
static int chain_walk(u64 a, u64 *chain, int cap, int *overflowed) {
    *overflowed = 0;
    int len = 0;
    u64 p = a;
    for (;;) {
        if (len < cap) chain[len] = p;
        len++;
        int ds = digit_sum(p);
        /* would p*ds + 1 overflow u64? */
        if (p > (UINT64_MAX - 1) / (u64)ds) { *overflowed = 1; return len; }
        u64 nxt = p * (u64)ds + 1;
        if (!is_prime(nxt)) return len;
        p = nxt;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <target_len> [lo] [hi]\n", argv[0]);
        return 1;
    }
    int target = atoi(argv[1]);
    if (target < 1) target = 1;
    u64 lo = (argc > 2) ? strtoull(argv[2], NULL, 10) : 2;
    u64 hi = (argc > 3) ? strtoull(argv[3], NULL, 10) : 100000000ULL;
    if (lo < 2) lo = 2;

    enum { CAP = 64 };
    u64 chain[CAP];
    long found = 0, overflow_hits = 0;

    for (u64 p = lo; p <= hi; p++) {
        /* Parity prune: an odd prime with an odd digit sum produces an even
         * (hence composite) successor, so it can never start a chain of
         * length >= 2. Skip it before the costly primality test. */
        if (target >= 2 && p != 2 && (digit_sum(p) & 1)) continue;
        if (!is_prime(p)) continue;

        int of;
        int len = chain_walk(p, chain, CAP, &of);
        if (of) overflow_hits++;
        if (len >= target) {
            found++;
            printf("length %d  start %llu:\n  ", len, (unsigned long long)p);
            for (int i = 0; i < len && i < CAP; i++)
                printf("%llu%s", (unsigned long long)chain[i],
                       (i + 1 < len && i + 1 < CAP) ? " -> " : "");
            if (of) printf("  [WARNING: next term exceeds 2^64 -- length may be longer]");
            printf("\n");
        }
    }

    fprintf(stderr, "\nscanned starts [%llu, %llu]; %ld chain(s) of length >= %d found.\n",
            (unsigned long long)lo, (unsigned long long)hi, found, target);
    if (overflow_hits)
        fprintf(stderr, "note: %ld chain(s) hit the 2^64 ceiling; use GMP to go further.\n",
                overflow_hits);
    return 0;
}
