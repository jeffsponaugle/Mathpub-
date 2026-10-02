/* sah.c -- simulated annealing for h-bases: look for k-element sets with h-range >= target.
 * usage: sah h k target seconds [maxelem] [threads]   (default maxelem = target, threads = 8) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
static int H, K, TGT, MAXE, NT = 8; static double SECS;
static int best_global = 0; static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;
static inline uint64_t rng(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x; }
static double tnow(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
/* h-range via min-stamp DP up to TGT+1 */
static void eval(const int *a, int k, uint8_t *cnt, int *range, int *cov) {
    int lim = TGT + 1;
    memset(cnt, 255, lim + 1); cnt[0] = 0;
    for (int i = 0; i < k; i++) { int ai = a[i]; if (ai > lim) continue; for (int x = ai; x <= lim; x++) { unsigned v = cnt[x - ai] + 1u; if (v < cnt[x]) cnt[x] = (uint8_t)v; } }
    int r = 0; while (r + 1 <= lim && cnt[r + 1] <= H) r++;
    int c = 0; for (int x = 1; x <= TGT; x++) c += (cnt[x] <= H);
    *range = r; *cov = c;
}
static void *worker(void *arg) {
    uint64_t seed = (uint64_t)(uintptr_t)arg * 0x9E3779B97F4A7C15ULL ^ (uint64_t)time(NULL);
    int a[64], b[64]; uint8_t *cnt = malloc(TGT + 2);
    double start = tnow();
    while (tnow() - start < SECS) {
        a[0] = 1;
        for (int i = 1; i < K; i++) { int v; again: v = 2 + rng(&seed) % (MAXE - 1); for (int j = 0; j < i; j++) if (a[j] == v) goto again; a[i] = v; }
        int r, c; eval(a, K, cnt, &r, &c); double cur = r + 0.01 * c;
        long iters = 300000; double temp;
        for (long it = 0; it < iters && tnow() - start < SECS; it++) {
            memcpy(b, a, sizeof(int) * K);
            int i = 1 + rng(&seed) % (K - 1); int v, tries = 0;
            do { if (rng(&seed) % 3) v = b[i] + (int)(rng(&seed) % 21) - 10; else v = 2 + rng(&seed) % (MAXE - 1); tries++; int dup = 0; for (int j = 0; j < K; j++) if (j != i && b[j] == v) dup = 1; if (v >= 2 && v <= MAXE && !dup) break; } while (tries < 50);
            if (tries >= 50) continue;
            b[i] = v;
            int r2, c2; eval(b, K, cnt, &r2, &c2); double s2 = r2 + 0.01 * c2;
            temp = 3.0 * (1.0 - (double)it / iters) + 0.05;
            if (s2 >= cur || exp((s2 - cur) / temp) > (rng(&seed) % 1000000) / 1e6) { memcpy(a, b, sizeof(int) * K); cur = s2; }
            if (r2 > best_global) {
                pthread_mutex_lock(&lk);
                if (r2 > best_global) { best_global = r2; int s[64]; memcpy(s, b, sizeof(int) * K);
                    for (int x = 0; x < K; x++) for (int y = x + 1; y < K; y++) if (s[y] < s[x]) { int t = s[x]; s[x] = s[y]; s[y] = t; }
                    printf("%s range %d:", r2 >= TGT ? "FOUND" : "best", r2); for (int x = 0; x < K; x++) printf(" %d", s[x]); printf("\n"); fflush(stdout); }
                pthread_mutex_unlock(&lk);
            }
        }
    }
    free(cnt); return NULL;
}
int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: sah h k target seconds [maxelem] [threads]\n"); return 1; }
    H = atoi(argv[1]); K = atoi(argv[2]); TGT = atoi(argv[3]); SECS = atof(argv[4]); MAXE = argc > 5 ? atoi(argv[5]) : TGT; if (argc > 6) NT = atoi(argv[6]);
    pthread_t th[64]; for (long i = 0; i < NT; i++) pthread_create(&th[i], NULL, worker, (void *)(i + 1));
    for (int i = 0; i < NT; i++) pthread_join(th[i], NULL);
    printf("done h=%d k=%d best range %d (target %d)\n", H, K, best_global, TGT); return 0;
}
