// brute.c -- O(n log n) reference for A061955: V(n) mod n by direct summation.
//
// V(n) = sum_{k=1..n} rev(k) * 2^{S(k-1)}, S(m) = total bit length of 1..m,
// rev(k) = k's binary digits reversed (kept at full width bitlen(k)).
// Valid for n < 2^32 (all products fit in 64 bits).
//
//   brute n1 n2 ...        print "n residue" for each n
//   brute -r LO HI         print "n residue" for every odd n in [LO, HI)
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t brute(uint64_t n)
{
    uint64_t v = 0, p = 1 % n;            // p = 2^{S(k-1)} mod n
    for (uint64_t k = 1; k <= n; k++) {
        int d = 0;
        uint64_t r = 0;
        for (uint64_t kk = k; kk; kk >>= 1, d++)
            r = (r << 1) | (kk & 1);       // reversed binary of k
        v = (v + (r % n) * p) % n;
        p = (p << d) % n;
    }
    return v;
}

int main(int argc, char **argv)
{
    if (argc == 4 && !strcmp(argv[1], "-r")) {
        uint64_t lo = strtoull(argv[2], 0, 0), hi = strtoull(argv[3], 0, 0);
        for (uint64_t n = lo | 1; n < hi; n += 2)
            printf("%llu %llu\n", (unsigned long long)n, (unsigned long long)brute(n));
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        uint64_t n = strtoull(argv[i], 0, 0);
        printf("%llu %llu\n", (unsigned long long)n, (unsigned long long)brute(n));
    }
    return 0;
}
