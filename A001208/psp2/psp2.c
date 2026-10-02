/*
 * psp2.c -- exhaustive search for extremal additive 2-bases
 *           (global postage stamp problem with h = 2 stamps, OEIS A001212).
 *
 * Find every set A = {a_0=0 < a_1=1 < a_2 < ... < a_k} of k positive
 * integers such that every integer in [1,T] is a sum of at most two elements
 * of A.  Running with T = best known + 1 and finding nothing proves
 * n_2(k) = best known; running with T = best known lists all extremal bases.
 *
 * Challis-style depth-first search ("K-program") with exact prunes only:
 *   admissibility        a_{j+1} <= n_2(A_j) + 1  (prefix must cover [1, a_{j+1}-1])
 *   element-wise bound   a_j <= n_2(j-1) + 1       (known A001212 values)
 *   gaps test            #uncovered values in [1,T] after A_j <= sum_{i=j+1..k} (i+1)
 *                        (a new element a_i creates at most i+1 new sums)
 *   last element         the smallest remaining gap g must be a_k + a_p (or 2 a_k),
 *                        so a_k is one of at most k+1 candidates g - a_p, g/2.
 * Coverage of [0,T] is a 256-bit mask (T <= 255); per-thread level stacks, no copying.
 *
 * Build:  clang -O3 -mcpu=apple-m1 -std=c11 -pthread -o psp2 psp2.c
 * Usage:  ./psp2 -k K -t T [-j threads] [-d splitdepth] [-q]
 *         ./psp2 -selftest
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>

#define NW 4
#define MAXK 40
typedef struct { uint64_t w[NW]; } bits;

static const int N2[] = {0, 2, 4, 8, 12, 16, 20, 26, 32, 40, 46, 54, 64, 72, 80,
                         92, 104, 116, 128, 140, 152, 164, 180, 196, 212};
#define N2_KNOWN 24

static int K, T;
static int UB[MAXK + 2];
static int CAP[MAXK + 2];
static int SPLIT = 8, NTHREADS = 10, QUIET = 0;
static bits TMASK;

static inline void bits_set(bits *b, int i) { b->w[i >> 6] |= 1ULL << (i & 63); }
static inline void bits_shl_or(bits *dst, const bits *src, int s) {
    int ws = s >> 6, bs = s & 63;
    for (int i = NW - 1; i >= ws; i--) {
        uint64_t v = src->w[i - ws] << bs;
        if (bs && i - ws - 1 >= 0) v |= src->w[i - ws - 1] >> (64 - bs);
        dst->w[i] |= v;
    }
}
static inline int bits_popcount_masked(const bits *b) {
    int c = 0; for (int i = 0; i < NW; i++) c += __builtin_popcountll(b->w[i] & TMASK.w[i]); return c;
}
static inline int bits_first_zero(const bits *b) {
    for (int i = 0; i < NW; i++) { uint64_t v = ~b->w[i]; if (v) return i * 64 + __builtin_ctzll(v); }
    return 64 * NW;
}

typedef struct { bits E, S; } lvl;
typedef struct { uint64_t nodes[MAXK + 2], pruned[MAXK + 2], found; } stats;
typedef struct { int a[MAXK + 2]; lvl L; } prefix_t;

static prefix_t *prefixes; static size_t nprefixes, prefix_cap;
static _Atomic size_t next_prefix;
static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;

static void report(const int *a, int j, stats *st) {
    st->found++;
    if (QUIET) return;
    pthread_mutex_lock(&out_lock);
    printf("FOUND k=%d range>=%d:", K, T);
    for (int i = 1; i <= j; i++) printf(" %d", a[i]);
    printf("\n"); fflush(stdout);
    pthread_mutex_unlock(&out_lock);
}

/* last level: choose a_K to cover every remaining gap */
static inline void leaf(int j, int *a, const lvl *L, stats *st) {
    bits G; int any = 0;
    for (int i = 0; i < NW; i++) { G.w[i] = TMASK.w[i] & ~L->S.w[i]; any |= (G.w[i] != 0); }
    int last = a[j];
    if (!any) { /* already covered with K-1 elements */ a[j + 1] = last + 1; report(a, j + 1, st); return; }
    int g = 0; for (int i = 0; i < NW; i++) if (G.w[i]) { g = i * 64 + __builtin_ctzll(G.w[i]); break; }
    int hi = g; if (hi > UB[j + 1]) hi = UB[j + 1];
    /* candidates z = g - a_p (p = 0..j) and z = g/2 */
    for (int p = 0; p <= j + 1; p++) {
        int z;
        if (p <= j) z = g - a[p]; else { if (g & 1) break; z = g / 2; }
        if (z <= last || z > hi) continue;
        st->nodes[j + 1]++;
        bits M = {{0}}; bits_shl_or(&M, &L->E, z); if (2 * z <= T) bits_set(&M, 2 * z);
        int ok = 1; for (int i = 0; i < NW; i++) if (G.w[i] & ~M.w[i]) { ok = 0; break; }
        if (ok) { a[j + 1] = z; report(a, j + 1, st); }
    }
}

static void dfs(int j, int *a, lvl *L, stats *st) {
    if (j + 1 == K) { leaf(j, a, &L[j], st); return; }
    int last = a[j];
    int hi = bits_first_zero(&L[j].S);           /* = range + 1 */
    if (hi > UB[j + 1]) hi = UB[j + 1];
    for (int x = last + 1; x <= hi; x++) {
        L[j + 1].S = L[j].S;
        bits_shl_or(&L[j + 1].S, &L[j].E, x);
        if (2 * x <= T) bits_set(&L[j + 1].S, 2 * x);
        st->nodes[j + 1]++;
        int gaps = (T + 1) - bits_popcount_masked(&L[j + 1].S);
        if (gaps > CAP[j + 1]) { st->pruned[j + 1]++; continue; }
        L[j + 1].E = L[j].E; bits_set(&L[j + 1].E, x);
        a[j + 1] = x;
        dfs(j + 1, a, L, st);
    }
}

static void gen_prefixes(int j, int *a, lvl *L) {
    if (j == SPLIT) {
        if (nprefixes == prefix_cap) { prefix_cap = prefix_cap ? 2 * prefix_cap : 4096; prefixes = realloc(prefixes, prefix_cap * sizeof(prefix_t)); }
        memcpy(prefixes[nprefixes].a, a, sizeof(int) * (MAXK + 2));
        prefixes[nprefixes].L = L[j];
        nprefixes++; return;
    }
    int last = a[j];
    int hi = bits_first_zero(&L[j].S); if (hi > UB[j + 1]) hi = UB[j + 1];
    for (int x = last + 1; x <= hi; x++) {
        L[j + 1].S = L[j].S; bits_shl_or(&L[j + 1].S, &L[j].E, x); if (2 * x <= T) bits_set(&L[j + 1].S, 2 * x);
        int gaps = (T + 1) - bits_popcount_masked(&L[j + 1].S);
        if (gaps > CAP[j + 1]) continue;
        L[j + 1].E = L[j].E; bits_set(&L[j + 1].E, x);
        a[j + 1] = x;
        gen_prefixes(j + 1, a, L);
    }
}

static stats *tstats;
static void *worker(void *arg) {
    stats *st = arg;
    int a[MAXK + 2]; lvl L[MAXK + 2];
    for (;;) {
        size_t i = atomic_fetch_add(&next_prefix, 1);
        if (i >= nprefixes) break;
        memcpy(a, prefixes[i].a, sizeof a);
        L[SPLIT] = prefixes[i].L;
        dfs(SPLIT, a, L, st);
    }
    return NULL;
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static uint64_t run(int k, int t, int verbose) {
    K = k; T = t;
    if (T > 64 * NW - 1) { fprintf(stderr, "T too large for NW=%d\n", NW); exit(1); }
    memset(&TMASK, 0, sizeof TMASK); for (int i = 0; i <= T; i++) bits_set(&TMASK, i);
    for (int j = 0; j <= K + 1; j++) {
        UB[j] = (j >= 1 && j - 1 <= N2_KNOWN) ? N2[j - 1] + 1 : T; if (UB[j] > T) UB[j] = T;
        int c = 0; for (int i = j + 1; i <= K; i++) c += i + 1; CAP[j] = c;
    }
    int save = SPLIT; if (SPLIT >= K - 1) SPLIT = K - 2; if (SPLIT < 1) SPLIT = 1;
    double t0 = now();
    nprefixes = 0;
    int a[MAXK + 2] = {0}; lvl L[MAXK + 2]; memset(L, 0, sizeof L);
    a[1] = 1; bits_set(&L[1].E, 0); bits_set(&L[1].E, 1); bits_set(&L[1].S, 0); bits_set(&L[1].S, 1); bits_set(&L[1].S, 2);
    gen_prefixes(1, a, L);
    double t1 = now();
    if (verbose) { printf("k=%d T=%d: %zu prefixes at depth %d (%.2fs)\n", K, T, nprefixes, SPLIT, t1 - t0); fflush(stdout); }
    atomic_store(&next_prefix, 0);
    tstats = calloc(NTHREADS, sizeof(stats));
    pthread_t th[256];
    for (int i = 0; i < NTHREADS; i++) pthread_create(&th[i], NULL, worker, &tstats[i]);
    for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
    double t2 = now();
    stats tot; memset(&tot, 0, sizeof tot);
    for (int i = 0; i < NTHREADS; i++) { tot.found += tstats[i].found; for (int j = 0; j <= K + 1; j++) { tot.nodes[j] += tstats[i].nodes[j]; tot.pruned[j] += tstats[i].pruned[j]; } }
    uint64_t total = 0; for (int j = 0; j <= K + 1; j++) total += tot.nodes[j];
    if (verbose) {
        printf("nodes per level:"); for (int j = 2; j <= K; j++) if (tot.nodes[j]) printf(" L%d=%llu", j, (unsigned long long)tot.nodes[j]); printf("\n");
        printf("gap-pruned per level:"); for (int j = 2; j <= K; j++) if (tot.pruned[j]) printf(" L%d=%llu", j, (unsigned long long)tot.pruned[j]); printf("\n");
        printf("RESULT k=%d T=%d found=%llu total_nodes=%llu time=%.2fs rate=%.2fM nodes/s threads=%d\n",
               K, T, (unsigned long long)tot.found, (unsigned long long)total, t2 - t0, total / (t2 - t0) / 1e6, NTHREADS);
        fflush(stdout);
    }
    free(tstats); SPLIT = save;
    return tot.found;
}

static void selftest(void) {
    struct { int k, n, count; } cases[] = {
        {3, 8, 1}, {4, 12, 1}, {5, 16, 1}, {6, 20, 5}, {7, 26, 3}, {8, 32, 2}, {9, 40, 1},
        {10, 46, 2}, {11, 54, 4}, {12, 64, 1}, {13, 72, 1}, {14, 80, 3}, {15, 92, 1}, {16, 104, 1}, {17, 116, 1}
    };
    int ok = 1, sq = QUIET; QUIET = 1;
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        double t0 = now();
        uint64_t f = run(cases[i].k, cases[i].n, 0);
        uint64_t g = run(cases[i].k, cases[i].n + 1, 0);
        int good = (f == (uint64_t)cases[i].count && g == 0);
        printf("selftest k=%2d n=%3d: found %llu (expect %d), at n+1: %llu (expect 0) %.2fs %s\n", cases[i].k, cases[i].n,
               (unsigned long long)f, cases[i].count, (unsigned long long)g, now() - t0, good ? "OK" : "FAIL");
        fflush(stdout);
        if (!good) ok = 0;
    }
    QUIET = sq;
    printf(ok ? "SELFTEST PASSED\n" : "SELFTEST FAILED\n");
    if (!ok) exit(1);
}

int main(int argc, char **argv) {
    int k = 0, t = 0, st = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-k")) k = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t")) t = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-j")) NTHREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-d")) SPLIT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-q")) QUIET = 1;
        else if (!strcmp(argv[i], "-selftest")) st = 1;
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    if (NTHREADS < 1) NTHREADS = 1; if (NTHREADS > 256) NTHREADS = 256;
    if (st) { selftest(); return 0; }
    if (k < 3 || t < 1) { fprintf(stderr, "usage: psp2 -k K -t T [-j threads] [-d splitdepth] [-q] | -selftest\n"); return 1; }
    run(k, t, 1);
    return 0;
}
