/* hsearch.c -- heuristic search (simulated annealing + seeded local search) for h-bases
 * with large h-range.  Produces LOWER BOUNDS for n_h(k) and a target for the exact search.
 *
 * usage: hsearch -h H -k K -max MAXELEM -lim LIM -sec SECONDS [-j threads] [-seed "a1 a2 ... ak"]...
 *        [-scan "a1 ... a_{k-1}"]   (scan all a_k in (a_{k-1}, MAXELEM] for the given (k-1)-prefix)
 * LIM = upper limit for the range evaluation (table size); choose ~1.5x the expected range.
 * Prints every improvement: "best <range>: a1 ... ak". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

static int H, K, NT = 8; static long MAXE = 0, LIM = 0; static double SECS = 60;
static long best_global = 0; static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;
static long seeds[64][64]; static int nseeds = 0;

static inline uint64_t rng(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x; }
static double tnow(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

/* h-range of sorted basis a[0..k-1] (a[0]=1) computed up to LIM; also coverage count */
static long hrange(const long *a, int k, uint16_t *cnt, long *cov) {
    for (long x = 0; x <= LIM; x++) cnt[x] = 0xFFFF;
    cnt[0] = 0;
    for (int i = 0; i < k; i++) {
        long ai = a[i]; if (ai > LIM) continue;
        for (long x = ai; x <= LIM; x++) { unsigned v = cnt[x - ai] + 1u; if (v < cnt[x]) cnt[x] = (uint16_t)v; }
    }
    long r = 0; while (r + 1 <= LIM && cnt[r + 1] <= H) r++;
    long c = 0; if (cov) { for (long x = 1; x <= LIM; x++) c += (cnt[x] <= H); *cov = c; }
    return r;
}
static void sort_basis(long *a, int k) { for (int i = 1; i < k; i++) { long v = a[i]; int j = i - 1; while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; } a[j + 1] = v; } }
static int valid(const long *a, int k) { if (a[0] != 1) return 0; for (int i = 1; i < k; i++) if (a[i] <= a[i - 1] || a[i] > MAXE) return 0; return 1; }

static void report(long r, const long *a, int k) {
    pthread_mutex_lock(&lk);
    if (r > best_global) {
        best_global = r;
        printf("best %ld:", r); for (int i = 0; i < k; i++) printf(" %ld", a[i]); printf("   [%s]\n", r >= LIM ? "HIT LIM - raise -lim" : "ok"); fflush(stdout);
    }
    pthread_mutex_unlock(&lk);
}

static void *worker(void *arg) {
    long id = (long)(uintptr_t)arg;
    uint64_t seed = (uint64_t)id * 0x9E3779B97F4A7C15ULL ^ (uint64_t)time(NULL) ^ ((uint64_t)clock() << 20);
    long a[64], b[64]; uint16_t *cnt = malloc(sizeof(uint16_t) * (LIM + 2));
    double start = tnow();
    int restart = 0;
    while (tnow() - start < SECS) {
        /* init: seeds (perturbed) or random geometric-ish */
        if (nseeds > 0 && (restart % 2 == 0 || (int)(rng(&seed) % 3))) {
            memcpy(a, seeds[rng(&seed) % nseeds], sizeof(long) * K);
            if (restart > 0) for (int i = 1; i < K; i++) { long d = (long)(rng(&seed) % 7) - 3; a[i] += d * (1 + a[i] / 200); }
        } else {
            a[0] = 1; double ratio = pow((double)MAXE, 1.0 / (K - 1));
            for (int i = 1; i < K; i++) { double v = pow(ratio, i) * (0.6 + 0.8 * (rng(&seed) % 1000) / 1000.0); a[i] = (long)v; if (a[i] < 2) a[i] = 2; if (a[i] > MAXE) a[i] = MAXE; }
        }
        sort_basis(a, K); for (int i = 1; i < K; i++) if (a[i] <= a[i - 1]) a[i] = a[i - 1] + 1; if (a[K - 1] > MAXE) a[K - 1] = MAXE;
        if (!valid(a, K)) { restart++; continue; }
        long cov; long r = hrange(a, K, cnt, &cov); double cur = r + 0.001 * cov;
        report(r, a, K);
        long iters = 20000 + 2000000 / (LIM / 1000 + 1); if (iters > 400000) iters = 400000;
        for (long it = 0; it < iters && tnow() - start < SECS; it++) {
            memcpy(b, a, sizeof(long) * K);
            int i = 1 + rng(&seed) % (K - 1);
            long span = 1 + b[i] / (8 + (rng(&seed) % 40));            /* relative step */
            long d = (long)(rng(&seed) % (2 * span + 1)) - span; if (d == 0) d = 1;
            if (rng(&seed) % 10 == 0) b[i] = 2 + rng(&seed) % (MAXE - 1); else b[i] += d;
            if (b[i] < 2) b[i] = 2; if (b[i] > MAXE) b[i] = MAXE;
            sort_basis(b, K);
            int dup = 0; for (int j = 1; j < K; j++) if (b[j] == b[j - 1]) dup = 1; if (dup) continue;
            long cov2; long r2 = hrange(b, K, cnt, &cov2); double s2 = r2 + 0.001 * cov2;
            double temp = (0.02 * r) * (1.0 - (double)it / iters) + 0.5;
            if (s2 >= cur || exp((s2 - cur) / temp) > (rng(&seed) % 1000000) / 1e6) { memcpy(a, b, sizeof(long) * K); cur = s2; }
            if (r2 > best_global) report(r2, b, K);
        }
        restart++;
    }
    free(cnt); return NULL;
}

int main(int argc, char **argv) {
    const char *scan = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h")) H = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k")) K = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-max")) MAXE = atol(argv[++i]);
        else if (!strcmp(argv[i], "-lim")) LIM = atol(argv[++i]);
        else if (!strcmp(argv[i], "-sec")) SECS = atof(argv[++i]);
        else if (!strcmp(argv[i], "-j")) NT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-seed")) { char *s = strdup(argv[++i]); int n = 0; for (char *t = strtok(s, " ,"); t; t = strtok(NULL, " ,")) seeds[nseeds][n++] = atol(t); if (n == K) nseeds++; else fprintf(stderr, "seed has %d elements, expected %d\n", n, K); }
        else if (!strcmp(argv[i], "-scan")) scan = argv[++i];
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    if (H < 1 || K < 2 || MAXE < 2 || LIM < 2) { fprintf(stderr, "usage: hsearch -h H -k K -max MAXELEM -lim LIM -sec S [-j T] [-seed \"...\"] [-scan \"...\"]\n"); return 1; }
    if (scan) {
        long a[64]; int n = 0; char *s = strdup(scan); for (char *t = strtok(s, " ,"); t; t = strtok(NULL, " ,")) a[n++] = atol(t);
        if (n != K - 1) { fprintf(stderr, "scan prefix must have k-1 elements\n"); return 1; }
        uint16_t *cnt = malloc(sizeof(uint16_t) * (LIM + 2)); long bestr = 0, besta = 0;
        for (long x = a[K - 2] + 1; x <= MAXE; x++) { a[K - 1] = x; long r = hrange(a, K, cnt, NULL); if (r > bestr) { bestr = r; besta = x; printf("scan best %ld with a_k=%ld\n", r, x); fflush(stdout); } }
        a[K - 1] = besta; memcpy(seeds[nseeds++], a, sizeof(long) * K);
        printf("scan done: best %ld:", bestr); for (int i = 0; i < K; i++) printf(" %ld", a[i]); printf("\n"); fflush(stdout);
        best_global = bestr;
    }
    pthread_t th[64]; for (long i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, (void *)(i + 1));
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    printf("done h=%d k=%d best range %ld\n", H, K, best_global); return 0;
}
