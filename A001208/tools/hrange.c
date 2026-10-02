/* hrange.c -- compute the h-range of a postage stamp basis.
 * usage: hrange h a1 a2 ... ak      (prints largest n with 1..n all sums of <= h elements)
 *        hrange -f file             (file lines: h a1 a2 ... ak ; prints "h k n a1..ak" per line)
 * Method: min-stamp-count DP (unbounded knapsack) over [0, h*ak]. uint16 counts (h up to 65000). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static long hrange(int h, int k, const long *a) {
    long lim = (long)h * a[k - 1] + 1;
    uint16_t *c = malloc((lim + 1) * sizeof(uint16_t));
    for (long i = 0; i <= lim; i++) c[i] = 0xFFFF;
    c[0] = 0;
    for (int i = 0; i < k; i++) {
        long ai = a[i];
        for (long x = ai; x <= lim; x++) {
            unsigned v = c[x - ai] + 1u;
            if (v < c[x]) c[x] = (uint16_t)v;
        }
    }
    long n = 0;
    while (n + 1 <= lim && c[n + 1] <= h) n++;
    free(c);
    return n;
}

int main(int argc, char **argv) {
    if (argc >= 3 && !strcmp(argv[1], "-f")) {
        FILE *f = fopen(argv[2], "r");
        if (!f) { perror("open"); return 1; }
        char line[4096];
        while (fgets(line, sizeof line, f)) {
            long a[64]; int k = 0; int h;
            char *p = line, *e;
            h = (int)strtol(p, &e, 10); if (e == p) continue; p = e;
            for (;;) { long v = strtol(p, &e, 10); if (e == p) break; a[k++] = v; p = e; }
            if (k == 0) continue;
            long n = hrange(h, k, a);
            printf("%d %d %ld", h, k, n);
            for (int i = 0; i < k; i++) printf(" %ld", a[i]);
            printf("\n");
            fflush(stdout);
        }
        return 0;
    }
    if (argc < 3) { fprintf(stderr, "usage: hrange h a1 ... ak | hrange -f file\n"); return 1; }
    int h = atoi(argv[1]); long a[64]; int k = 0;
    for (int i = 2; i < argc; i++) a[k++] = atol(argv[i]);
    printf("%ld\n", hrange(h, k, a));
    return 0;
}
