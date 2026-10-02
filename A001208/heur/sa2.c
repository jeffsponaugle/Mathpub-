/* sa2.c -- simulated annealing for 2-bases: look for k-element sets with 2-range >= target.
 * usage: sa2 k target seconds [maxelem]
 * Objective: lexicographic (range, #covered in [1,target]). Moves: replace one element. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

static int K, TGT, MAXE;
static double SECS;
static int best_global = 0; static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;

static inline uint64_t rng(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x; }

static void eval(const int *a, int k, int *range, int *cov) {
    static __thread uint8_t m[1024];
    memset(m, 0, 2 * MAXE + 2);
    m[0] = 1;
    for (int i = 0; i < k; i++) { m[a[i]] = 1; for (int j = i; j < k; j++) m[a[i] + a[j]] = 1; }
    int r = 0; while (r + 1 <= 2 * MAXE && m[r + 1]) r++;
    int c = 0; for (int x = 1; x <= TGT; x++) c += m[x];
    *range = r; *cov = c;
}
static double score(int r, int c) { return r * 1.0 + c * 0.02; }

static void *worker(void *arg) {
    uint64_t seed = (uint64_t)(uintptr_t)arg * 0x9E3779B97F4A7C15ULL + time(NULL);
    int a[64], b[64];
    double t0 = (double)clock() / CLOCKS_PER_SEC; (void)t0;
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); double start = ts.tv_sec + ts.tv_nsec * 1e-9;
    for (;;) {
        clock_gettime(CLOCK_MONOTONIC, &ts); if (ts.tv_sec + ts.tv_nsec * 1e-9 - start > SECS) break;
        /* random admissible-ish init: 1 plus random distinct values */
        a[0] = 1;
        for (int i = 1; i < K; i++) { int v; again: v = 2 + rng(&seed) % (MAXE - 1); for (int j = 0; j < i; j++) if (a[j] == v) goto again; a[i] = v; }
        int r, c; eval(a, K, &r, &c); double cur = score(r, c), bestl = cur;
        double temp = 2.0;
        for (long it = 0; it < 400000; it++) {
            memcpy(b, a, sizeof(int) * K);
            int i = 1 + rng(&seed) % (K - 1);
            int v; int tries = 0;
            do { if (rng(&seed) & 1) v = b[i] + (int)(rng(&seed) % 11) - 5; else v = 2 + rng(&seed) % (MAXE - 1); tries++; int dup = 0; for (int j = 0; j < K; j++) if (j != i && b[j] == v) dup = 1; if (v >= 2 && v <= MAXE && !dup) break; } while (tries < 50);
            if (tries >= 50) continue;
            b[i] = v;
            int r2, c2; eval(b, K, &r2, &c2); double s2 = score(r2, c2);
            if (s2 >= cur || exp((s2 - cur) / temp) > (rng(&seed) % 1000000) / 1e6) { memcpy(a, b, sizeof(int) * K); cur = s2; if (cur > bestl) bestl = cur; }
            temp = 2.0 * (1.0 - (double)it / 400000) + 0.05;
            if (r2 >= TGT) {
                pthread_mutex_lock(&lk);
                printf("FOUND range %d:", r2); int s[64]; memcpy(s, b, sizeof(int) * K);
                for (int x = 0; x < K; x++) for (int y = x + 1; y < K; y++) if (s[y] < s[x]) { int t = s[x]; s[x] = s[y]; s[y] = t; }
                for (int x = 0; x < K; x++) printf(" %d", s[x]); printf("\n"); fflush(stdout);
                pthread_mutex_unlock(&lk);
            }
            if (r2 > best_global) { pthread_mutex_lock(&lk); if (r2 > best_global) { best_global = r2; printf("best so far range %d\n", r2); fflush(stdout); } pthread_mutex_unlock(&lk); }
        }
    }
    return NULL;
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: sa2 k target seconds [maxelem]\n"); return 1; }
    K = atoi(argv[1]); TGT = atoi(argv[2]); SECS = atof(argv[3]); MAXE = argc > 4 ? atoi(argv[4]) : TGT;
    if (MAXE > 500) MAXE = 500;
    pthread_t th[8]; for (long i = 0; i < 8; i++) pthread_create(&th[i], NULL, worker, (void *)(i + 1));
    for (int i = 0; i < 8; i++) pthread_join(th[i], NULL);
    printf("done, best range %d (target %d)\n", best_global, TGT);
    return 0;
}
