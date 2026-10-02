/* hsearch2.c -- coordinate-scan local search for h-bases with large h-range (lower bounds / targets).
 * usage: hsearch2 -h H -k K -max MAXELEM -lim LIM -sec S [-j T] [-seed "a1 ... ak"]... [-kick K]
 * Starts from each seed (and scaled/perturbed copies), repeatedly scans every single coordinate over its
 * whole feasible interval (exact 1-D optimisation), then applies random kicks and re-descends.
 * Prints every improvement "best <range>: a1 ... ak". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

static int H, K, NT = 8; static long MAXE = 0, LIM = 0; static double SECS = 60; static int KICK = 2;
static long best_global = 0; static long best_basis[64]; static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;
static long seeds[64][64]; static int nseeds = 0;
static inline uint64_t rng(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x; }
static double tnow(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static long hrange(const long *a, int k, uint16_t *cnt) {
    for (long x = 0; x <= LIM; x++) cnt[x] = 0xFFFF;
    cnt[0] = 0;
    for (int i = 0; i < k; i++) {
        long ai = a[i]; if (ai > LIM) continue;
        for (long x = ai; x <= LIM; x++) { unsigned v = cnt[x - ai] + 1u; if (v < cnt[x]) cnt[x] = (uint16_t)v; }
    }
    long r = 0; while (r + 1 <= LIM && cnt[r + 1] <= H) r++;
    return r;
}
static void report(long r, const long *a, int k) {
    pthread_mutex_lock(&lk);
    if (r > best_global) { best_global = r; memcpy(best_basis, a, sizeof(long) * k);
        printf("best %ld:", r); for (int i = 0; i < k; i++) printf(" %ld", a[i]); printf("%s\n", r >= LIM ? "   [HIT LIM - raise -lim]" : ""); fflush(stdout); }
    pthread_mutex_unlock(&lk);
}
/* exact coordinate descent; returns range */
static long descend(long *a, uint16_t *cnt, uint64_t *seed, double deadline) {
    long cur = hrange(a, K, cnt);
    int improved = 1;
    while (improved && tnow() < deadline) {
        improved = 0;
        int order[64]; for (int i = 0; i < K; i++) order[i] = i;
        for (int i = K - 1; i > 1; i--) { int j = 1 + rng(seed) % i; int t = order[i]; order[i] = order[j]; order[j] = t; }
        for (int oi = 1; oi < K; oi++) {
            int i = order[oi]; if (i == 0) continue;
            long lo = a[i - 1] + 1, hi = (i < K - 1) ? a[i + 1] - 1 : MAXE;
            if (hi < lo) continue;
            long span = hi - lo + 1, stride = 1; if (span > 6000) stride = span / 6000 + 1;   /* coarse pass for huge spans */
            long bestx = a[i], bestr = cur, save = a[i];
            for (long x = lo; x <= hi; x += stride) { a[i] = x; long r = hrange(a, K, cnt); if (r > bestr) { bestr = r; bestx = x; } }
            if (stride > 1) { long c = bestx; for (long x = c - stride; x <= c + stride; x++) { if (x < lo || x > hi) continue; a[i] = x; long r = hrange(a, K, cnt); if (r > bestr) { bestr = r; bestx = x; } } }
            a[i] = bestx;
            if (bestr > cur) { cur = bestr; improved = 1; report(cur, a, K); } else a[i] = save;
            if (tnow() > deadline) break;
        }
    }
    return cur;
}
static void *worker(void *arg) {
    long id = (long)(uintptr_t)arg;
    uint64_t seed = (uint64_t)id * 0x9E3779B97F4A7C15ULL ^ (uint64_t)time(NULL);
    long a[64], b[64]; uint16_t *cnt = malloc(sizeof(uint16_t) * (LIM + 2));
    double deadline = tnow() + SECS; int round = 0;
    while (tnow() < deadline) {
        if (nseeds > 0 && (round == 0 || rng(&seed) % 3)) {
            memcpy(a, seeds[(id + round) % nseeds], sizeof(long) * K);
            if (round > 0) { double f = 0.85 + 0.3 * (rng(&seed) % 1000) / 1000.0; for (int i = 1; i < K; i++) { a[i] = (long)(a[i] * pow(f, (double)i / (K - 1))); } }
        } else {
            a[0] = 1; double ratio = pow((double)MAXE, 1.0 / (K - 1));
            for (int i = 1; i < K; i++) a[i] = (long)(pow(ratio, i) * (0.6 + 0.8 * (rng(&seed) % 1000) / 1000.0));
        }
        for (int i = 1; i < K; i++) { if (a[i] <= a[i - 1]) a[i] = a[i - 1] + 1; if (a[i] > MAXE) a[i] = MAXE - (K - 1 - i); }
        long cur = descend(a, cnt, &seed, deadline);
        /* kick loop */
        int stale = 0;
        while (stale < 6 && tnow() < deadline) {
            memcpy(b, a, sizeof(long) * K);
            int nk = 1 + rng(&seed) % KICK;
            for (int t = 0; t < nk; t++) { int i = 1 + rng(&seed) % (K - 1); long span = 1 + b[i] / (4 + rng(&seed) % 12); b[i] += (long)(rng(&seed) % (2 * span + 1)) - span; }
            for (int i = 1; i < K; i++) { if (b[i] < 2) b[i] = 2; if (b[i] > MAXE) b[i] = MAXE; }
            for (int i = 1; i < K; i++) for (int j = i; j > 0 && b[j] < b[j - 1]; j--) { long t = b[j]; b[j] = b[j - 1]; b[j - 1] = t; }
            for (int i = 1; i < K; i++) if (b[i] <= b[i - 1]) b[i] = b[i - 1] + 1;
            if (b[K - 1] > MAXE) continue;
            long r = descend(b, cnt, &seed, deadline);
            if (r > cur) { cur = r; memcpy(a, b, sizeof(long) * K); stale = 0; } else stale++;
        }
        round++;
    }
    free(cnt); return NULL;
}
int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h")) H = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k")) K = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-max")) MAXE = atol(argv[++i]);
        else if (!strcmp(argv[i], "-lim")) LIM = atol(argv[++i]);
        else if (!strcmp(argv[i], "-sec")) SECS = atof(argv[++i]);
        else if (!strcmp(argv[i], "-j")) NT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-kick")) KICK = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-seed")) { char *s = strdup(argv[++i]); int n = 0; for (char *t = strtok(s, " ,"); t; t = strtok(NULL, " ,")) seeds[nseeds][n++] = atol(t); if (n == K) nseeds++; else fprintf(stderr, "seed has %d elements, expected %d\n", n, K); }
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    if (H < 1 || K < 2 || MAXE < 2 || LIM < 2) { fprintf(stderr, "usage: hsearch2 -h H -k K -max MAXELEM -lim LIM -sec S [-j T] [-seed \"...\"]...\n"); return 1; }
    pthread_t th[64]; for (long i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, (void *)(i + 1));
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    printf("done h=%d k=%d best range %ld:", H, K, best_global); for (int i = 0; i < K; i++) printf(" %ld", best_basis[i]); printf("\n"); return 0;
}
