/*
 * psphs.c -- exhaustive search for extremal h-bases with SMALL h (2..8) and
 *            moderate k: the "rows" of the postage stamp table
 *            (A001212 h=2, A001213 h=3, A001214 h=4, A001215 h=5, A001216 h=6,
 *             A005342 h=7, A005343 h=8).
 *
 * Find every A = {a_1=1 < a_2 < ... < a_k} such that every integer in [1,T]
 * is a sum of at most h elements of A (repetition allowed).
 *
 * Bitset DFS (Challis "K-program" style):
 *   S[m] = set of sums of exactly m elements of A u {0}  (= sums of <= m elements of A),
 *   kept as bit masks over [0, 64*NW).  Adding x:  S[m] |= S_old[m-c] << (c*x), c=1..m.
 *   Prunes (all exact, never discard a candidate that could reach T):
 *     admissibility      a_{j+1} <= n_h(A_j) + 1     (= first gap of S[h])
 *     element-wise bound a_j <= n_h(j-1) + 1          (known extremal values)
 *     gaps test          #uncovered in [1,T] <= sum_{i=j+1..k} C(i+h-1, h-1)
 *                        (each new a_i creates at most C(i+h-1,h-1) new sums)
 *
 * Build: clang -O3 -mcpu=apple-m1 -std=c11 -pthread -o psphs psphs.c
 * Usage: ./psphs -h H -k K -t T [-j threads] [-d splitdepth] [-q]
 *        ./psphs -selftest
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>

#define MAXNW 16                /* up to 1024-bit masks: T <= 1023 */
#define MAXH 8
#define MAXK 48

static int H, K, T, NW;
static int UB[MAXK + 2];
static long CAP[MAXK + 2];
static int SPLIT = 6, NTHREADS = 10, QUIET = 0;
static uint64_t TMASK[MAXNW];

/* known extremal ranges n_h(k) (OEIS, proven values only), index [h][k] */
static int known(int h, int k) {
    static const int r2[] = {0,2,4,8,12,16,20,26,32,40,46,54,64,72,80,92,104,116,128,140,152,164,180,196,212};
    static const int r3[] = {0,3,7,15,24,36,52,70,93,121,154,186,225,271,323,385};
    static const int r4[] = {0,4,10,26,44,70,108,162,228,310,422,550,700};
    static const int r5[] = {0,5,14,35,71,126,211,336,524,726,1016};
    static const int r6[] = {0,6,18,52,114,216,388,638,1007,1545};
    static const int r7[] = {0,7,23,69,165,345,664,1137,1911};
    static const int r8[] = {0,8,28,89,234,512,1045,2001,3485};
    const int *r = NULL; int n = 0;
    switch (h) {
        case 2: r = r2; n = 24; break; case 3: r = r3; n = 15; break; case 4: r = r4; n = 12; break;
        case 5: r = r5; n = 10; break; case 6: r = r6; n = 9; break; case 7: r = r7; n = 8; break; case 8: r = r8; n = 8; break;
    }
    if (k == 0) return 0;
    if (r && k <= n) return r[k];
    return -1;
}

typedef struct {
    int j;
    int a[MAXK + 2];
    uint64_t S[MAXH + 1][MAXNW];
} node;

typedef struct { uint64_t nodes[MAXK + 2], pruned[MAXK + 2], found; } stats;

static node *prefixes; static size_t nprefixes, prefix_cap;
static _Atomic size_t next_prefix;
static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;

static inline void shl_or(uint64_t *dst, const uint64_t *src, int s) {
    if (s >= 64 * NW) return;
    int ws = s >> 6, bs = s & 63;
    for (int i = NW - 1; i >= ws; i--) {
        uint64_t v = src[i - ws] << bs;
        if (bs && i - ws - 1 >= 0) v |= src[i - ws - 1] >> (64 - bs);
        dst[i] |= v;
    }
}
static inline int first_zero(const uint64_t *b) {
    for (int i = 0; i < NW; i++) { uint64_t v = ~b[i]; if (v) return i * 64 + __builtin_ctzll(v); }
    return 64 * NW;
}
static inline int covered_count(const uint64_t *b) {
    int c = 0; for (int i = 0; i < NW; i++) c += __builtin_popcountll(b[i] & TMASK[i]); return c;
}
/* add element x to node (in place) */
static inline void add_elem(node *nd, int x) {
    for (int m = H; m >= 1; m--)
        for (int c = 1; c <= m; c++)
            shl_or(nd->S[m], nd->S[m - c], c * x);
}

static void report(const node *nd, stats *st) {
    st->found++;
    if (QUIET) return;
    pthread_mutex_lock(&out_lock);
    printf("FOUND h=%d k=%d range>=%d:", H, K, T);
    for (int i = 1; i <= nd->j; i++) printf(" %d", nd->a[i]);
    printf("\n"); fflush(stdout);
    pthread_mutex_unlock(&out_lock);
}

static void dfs(const node *nd, stats *st) {
    int j = nd->j, last = nd->a[j];
    int hi = first_zero(nd->S[H]);            /* = range + 1 */
    if (hi > UB[j + 1]) hi = UB[j + 1];
    if (j + 1 == K) {
        for (int x = last + 1; x <= hi; x++) {
            node c = *nd; add_elem(&c, x); st->nodes[j + 1]++;
            if (covered_count(c.S[H]) == T + 1) { c.a[j + 1] = x; c.j = j + 1; report(&c, st); }
        }
        return;
    }
    for (int x = last + 1; x <= hi; x++) {
        node c = *nd; add_elem(&c, x); st->nodes[j + 1]++;
        long gaps = (T + 1) - covered_count(c.S[H]);
        if (gaps > CAP[j + 1]) { st->pruned[j + 1]++; continue; }
        c.a[j + 1] = x; c.j = j + 1;
        dfs(&c, st);
    }
}

static void gen_prefixes(const node *nd) {
    if (nd->j == SPLIT) {
        if (nprefixes == prefix_cap) { prefix_cap = prefix_cap ? 2 * prefix_cap : 4096; prefixes = realloc(prefixes, prefix_cap * sizeof(node)); }
        prefixes[nprefixes++] = *nd; return;
    }
    int j = nd->j, last = nd->a[j];
    int hi = first_zero(nd->S[H]); if (hi > UB[j + 1]) hi = UB[j + 1];
    for (int x = last + 1; x <= hi; x++) {
        node c = *nd; add_elem(&c, x);
        long gaps = (T + 1) - covered_count(c.S[H]);
        if (gaps > CAP[j + 1]) continue;
        c.a[j + 1] = x; c.j = j + 1;
        gen_prefixes(&c);
    }
}

static stats *tstats;
static void *worker(void *arg) {
    stats *st = arg;
    for (;;) { size_t i = atomic_fetch_add(&next_prefix, 1); if (i >= nprefixes) break; dfs(&prefixes[i], st); }
    return NULL;
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
static long binom(int n, int r) { if (r < 0 || r > n) return 0; long v = 1; for (int i = 1; i <= r; i++) v = v * (n - r + i) / i; return v; }

static uint64_t run(int h, int k, int t, int verbose) {
    H = h; K = k; T = t;
    NW = (T + 64) / 64; if (NW > MAXNW) { fprintf(stderr, "T too large\n"); exit(1); }
    memset(TMASK, 0, sizeof TMASK);
    for (int i = 0; i <= T; i++) TMASK[i >> 6] |= 1ULL << (i & 63);
    for (int j = 0; j <= K + 1; j++) {
        int kn = (j >= 1) ? known(H, j - 1) : -1;
        UB[j] = (kn >= 0) ? kn + 1 : T;
        if (UB[j] > T) UB[j] = T;
        long c = 0; for (int i = j + 1; i <= K; i++) c += binom(i + H - 1, H - 1);
        CAP[j] = c;
    }
    int split = SPLIT; if (split >= K) split = K - 1; if (split < 1) split = 1;
    int save_split = SPLIT; SPLIT = split;
    double t0 = now();
    nprefixes = 0;
    node root; memset(&root, 0, sizeof root);
    root.j = 1; root.a[0] = 0; root.a[1] = 1;
    root.S[0][0] = 1;                       /* {0} */
    for (int m = 1; m <= H; m++) root.S[m][0] = (m + 1 < 64) ? ((1ULL << (m + 1)) - 1) : ~0ULL;  /* sums of <= m ones: [0,m] */
    gen_prefixes(&root);
    double t1 = now();
    if (verbose) { printf("h=%d k=%d T=%d: %zu prefixes at depth %d (%.2fs)\n", H, K, T, nprefixes, SPLIT, t1 - t0); fflush(stdout); }
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
        printf("RESULT h=%d k=%d T=%d found=%llu total_nodes=%llu time=%.2fs rate=%.2fM nodes/s threads=%d\n",
               H, K, T, (unsigned long long)tot.found, (unsigned long long)total, t2 - t0, total / (t2 - t0) / 1e6, NTHREADS);
        fflush(stdout);
    }
    free(tstats); SPLIT = save_split;
    return tot.found;
}

static void selftest(void) {
    /* h, k, n(h,k), expected number of extremal bases per Challis-Robinson tables (where listed) */
    struct { int h, k, n, cnt; } cases[] = {
        {2,6,20,5},{2,7,26,3},{2,10,46,2},{2,11,54,4},{2,12,64,1},
        {3,3,15,1},{3,4,24,1},{3,5,36,1},{3,6,52,2},{3,7,70,1},{3,8,93,1},{3,9,121,1},{3,10,154,1},
        {4,3,26,1},{4,4,44,1},{4,5,70,1},{4,6,108,2},{4,7,162,3},{4,8,228,1},
        {5,3,35,1},{5,4,71,2},{5,5,126,1},{5,6,211,2},{5,7,336,1},
        {6,3,52,1},{6,4,114,1},{6,5,216,1},{6,6,388,1},
        {7,3,69,1},{7,4,165,1},{7,5,345,1},
        {8,3,89,1},{8,4,234,1},{8,5,512,2},
    };
    int ok = 1, sq = QUIET; QUIET = 1;
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        double t0 = now();
        uint64_t f = run(cases[i].h, cases[i].k, cases[i].n, 0);
        uint64_t g = run(cases[i].h, cases[i].k, cases[i].n + 1, 0);
        int good = (f == (uint64_t)cases[i].cnt && g == 0);
        printf("selftest h=%d k=%2d n=%4d: found %llu (expect %d), at n+1: %llu (expect 0) %.2fs %s\n", cases[i].h, cases[i].k, cases[i].n,
               (unsigned long long)f, cases[i].cnt, (unsigned long long)g, now() - t0, good ? "OK" : "FAIL");
        fflush(stdout);
        if (!good) ok = 0;
    }
    QUIET = sq;
    printf(ok ? "SELFTEST PASSED\n" : "SELFTEST FAILED\n");
}

int main(int argc, char **argv) {
    int h = 0, k = 0, t = 0, st = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h")) h = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k")) k = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t")) t = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-j")) NTHREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-d")) SPLIT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-q")) QUIET = 1;
        else if (!strcmp(argv[i], "-selftest")) st = 1;
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    if (NTHREADS < 1) NTHREADS = 1; if (NTHREADS > 256) NTHREADS = 256;
    if (st) { selftest(); return 0; }
    if (h < 2 || h > MAXH || k < 3 || t < 1) { fprintf(stderr, "usage: psphs -h H -k K -t T [-j threads] [-d depth] [-q] | -selftest\n"); return 1; }
    run(h, k, t, 1);
    return 0;
}
