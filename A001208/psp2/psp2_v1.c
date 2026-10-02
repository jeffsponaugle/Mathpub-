/*
 * psp2.c -- exhaustive search for extremal additive 2-bases
 *           (global postage stamp problem with h = 2 stamps, OEIS A001212).
 *
 * Find every set A = {a_0=0 < a_1=1 < a_2 < ... < a_k} of k positive
 * integers such that every integer in [1,T] is a sum of at most two elements
 * of A (a_i + a_j with 0 <= i <= j <= k).  Running with T = best known + 1
 * and finding nothing proves n_2(k) = best known.
 *
 * Algorithm: Challis-style depth-first search (the "K-program"):
 *   - admissibility: a_{j+1} <= n_2(A_j) + 1   (prefix must cover [1,a_{j+1}-1])
 *   - element-wise upper bound: a_j <= n_2(j-1) + 1  (known A001212 values)
 *   - gaps test: the number of uncovered values in [1,T] after A_j must not
 *     exceed sum_{i=j+1..k} (i+1), the most new sums the remaining elements
 *     can contribute (a_i + a_p for p = 0..i-1, plus 2 a_i).
 *   - optional largest-element split (-M): fix a_k = M first, so that the
 *     sums M + a_p are placed immediately and the gaps test bites earlier.
 *
 * Coverage of [0,T] is kept as a 256-bit mask (T <= 255), elements as a
 * second 256-bit mask; adding x costs one 256-bit shift-or.
 *
 * Build:  clang -O3 -mcpu=apple-m1 -std=c11 -pthread -o psp2 psp2.c
 * Usage:  ./psp2 -k K -t T [-j threads] [-d splitdepth] [-M] [-q]
 *         ./psp2 -selftest
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>

#define NW 4                    /* 256-bit masks: T <= 255 */
#define MAXK 40

typedef struct { uint64_t w[NW]; } bits;

/* Known extremal 2-ranges n_2(k), k = 0..24 (A001212, n_2(0)=0). */
static const int N2[] = {0, 2, 4, 8, 12, 16, 20, 26, 32, 40, 46, 54, 64, 72, 80,
                         92, 104, 116, 128, 140, 152, 164, 180, 196, 212};
#define N2_KNOWN 24

static int K, T;                /* length and target range */
static int UB[MAXK + 2];        /* UB[j]: max allowed a_j */
static int CAP[MAXK + 2];       /* CAP[j]: max new values the elements a_{j+1..k} can add */
static int SPLIT = 9;           /* prefix depth handed to threads */
static int NTHREADS = 10;
static int QUIET = 0;
static int MSPLIT = 0;          /* fix largest element first */
static bits TMASK;              /* bits 0..T set */

static inline void bits_clear(bits *b) { memset(b, 0, sizeof *b); }
static inline void bits_set(bits *b, int i) { b->w[i >> 6] |= 1ULL << (i & 63); }
static inline int bits_test(const bits *b, int i) { return (b->w[i >> 6] >> (i & 63)) & 1; }
/* dst |= (src << s), bits shifted beyond 64*NW are dropped */
static inline void bits_shl_or(bits *dst, const bits *src, int s) {
    int ws = s >> 6, bs = s & 63;
    for (int i = NW - 1; i >= ws; i--) {
        uint64_t v = src->w[i - ws] << bs;
        if (bs && i - ws - 1 >= 0) v |= src->w[i - ws - 1] >> (64 - bs);
        dst->w[i] |= v;
    }
}
static inline int bits_popcount_masked(const bits *b) {
    int c = 0;
    for (int i = 0; i < NW; i++) c += __builtin_popcountll(b->w[i] & TMASK.w[i]);
    return c;
}
/* first zero bit (there is always one below 64*NW unless all set) */
static inline int bits_first_zero(const bits *b) {
    for (int i = 0; i < NW; i++) {
        uint64_t v = ~b->w[i];
        if (v) return i * 64 + __builtin_ctzll(v);
    }
    return 64 * NW;
}

/* ---------- search state ---------- */
typedef struct {
    int j;                      /* number of elements placed (a_1..a_j) */
    int a[MAXK + 2];
    bits E;                     /* element mask (bits at 0 and a_1..a_j, and M if MSPLIT) */
    bits S;                     /* sum mask (sums <= 255) */
} node;

typedef struct {
    uint64_t nodes[MAXK + 2];   /* nodes visited per level */
    uint64_t pruned_gap[MAXK + 2];
    uint64_t found;
} stats;

static node *prefixes;
static size_t nprefixes, prefix_cap;
static _Atomic size_t next_prefix;
static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;

static void report(const node *nd, stats *st) {
    /* exact range: first zero of full sum mask */
    int r = bits_first_zero(&nd->S) - 1;
    st->found++;
    if (QUIET) return;
    pthread_mutex_lock(&out_lock);
    printf("FOUND k=%d range>=%d (range in [0,255] window = %d):", K, T, r);
    for (int i = 1; i <= nd->j; i++) printf(" %d", nd->a[i]);
    if (MSPLIT) printf(" %d", nd->a[K]);
    printf("\n");
    fflush(stdout);
    pthread_mutex_unlock(&out_lock);
}

/* depth-first search continuing from node nd; last level handled inline */
static void dfs(node *nd, stats *st) {
    int j = nd->j;
    int last = nd->a[j];
    int r = bits_first_zero(&nd->S) - 1;          /* range of current set */
    int hi = r + 1;                                /* admissibility */
    if (hi > UB[j + 1]) hi = UB[j + 1];
    int klast = MSPLIT ? K - 1 : K;                /* last free level */
    if (MSPLIT && hi > nd->a[K] - (K - 1 - j)) hi = nd->a[K] - (K - 1 - j); /* room for the rest below M */
    if (j + 1 == klast) {
        /* final free element: just test coverage */
        for (int x = last + 1; x <= hi; x++) {
            node c = *nd;
            bits_shl_or(&c.S, &nd->E, x);
            if (2 * x <= T) bits_set(&c.S, 2 * x);
            st->nodes[j + 1]++;
            if (bits_popcount_masked(&c.S) == T + 1) {
                c.a[j + 1] = x; c.j = j + 1;
                report(&c, st);
            }
        }
        return;
    }
    for (int x = last + 1; x <= hi; x++) {
        node c = *nd;
        bits_shl_or(&c.S, &nd->E, x);
        if (2 * x <= T) bits_set(&c.S, 2 * x);
        st->nodes[j + 1]++;
        int gaps = (T + 1) - bits_popcount_masked(&c.S);
        if (gaps > CAP[j + 1]) { st->pruned_gap[j + 1]++; continue; }
        bits_set(&c.E, x);
        c.a[j + 1] = x; c.j = j + 1;
        dfs(&c, st);
    }
}

/* enumerate admissible prefixes to depth SPLIT (no gaps test needed; apply anyway) */
static void gen_prefixes(node *nd) {
    if (nd->j == SPLIT) {
        if (nprefixes == prefix_cap) {
            prefix_cap = prefix_cap ? prefix_cap * 2 : 1024;
            prefixes = realloc(prefixes, prefix_cap * sizeof(node));
        }
        prefixes[nprefixes++] = *nd;
        return;
    }
    int j = nd->j, last = nd->a[j];
    int r = bits_first_zero(&nd->S) - 1;
    int hi = r + 1;
    if (hi > UB[j + 1]) hi = UB[j + 1];
    if (MSPLIT && hi > nd->a[K] - (K - 1 - j)) hi = nd->a[K] - (K - 1 - j);
    for (int x = last + 1; x <= hi; x++) {
        node c = *nd;
        bits_shl_or(&c.S, &nd->E, x);
        if (2 * x <= T) bits_set(&c.S, 2 * x);
        int gaps = (T + 1) - bits_popcount_masked(&c.S);
        if (gaps > CAP[j + 1]) continue;
        bits_set(&c.E, x);
        c.a[j + 1] = x; c.j = j + 1;
        gen_prefixes(&c);
    }
}

static stats *tstats;
static void *worker(void *arg) {
    stats *st = (stats *)arg;
    for (;;) {
        size_t i = atomic_fetch_add(&next_prefix, 1);
        if (i >= nprefixes) break;
        node nd = prefixes[i];
        dfs(&nd, st);
    }
    return NULL;
}

static double now(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void setup_tables(void) {
    bits_clear(&TMASK);
    for (int i = 0; i <= T; i++) bits_set(&TMASK, i);
    for (int j = 0; j <= K + 1; j++) {
        /* a_j <= n_2(j-1) + 1 if known, else <= T */
        UB[j] = (j - 1 <= N2_KNOWN && j >= 1) ? N2[j - 1] + 1 : T;
        if (UB[j] > T) UB[j] = T;
    }
    for (int j = 0; j <= K + 1; j++) {
        int c = 0;
        for (int i = j + 1; i <= (MSPLIT ? K - 1 : K); i++) c += i + 1;
        if (MSPLIT) c += 0; /* M's sums with later elements are counted in their (i+1) */
        CAP[j] = c;
    }
}

/* run one search; returns number of bases found with range >= T */
static uint64_t run(int k, int t, int msplit, int verbose) {
    K = k; T = t; MSPLIT = msplit;
    if (T > 64 * NW - 1) { fprintf(stderr, "T too large for NW=%d\n", NW); exit(1); }
    setup_tables();
    double t0 = now();
    nprefixes = 0;
    int split = SPLIT; if (split >= (MSPLIT ? K - 1 : K)) split = (MSPLIT ? K - 1 : K) - 1;
    if (split < 1) split = 1;
    SPLIT = split;
    node root; memset(&root, 0, sizeof root);
    root.j = 1; root.a[0] = 0; root.a[1] = 1;
    bits_set(&root.E, 0); bits_set(&root.E, 1);
    bits_set(&root.S, 0); bits_set(&root.S, 1); bits_set(&root.S, 2);
    if (!MSPLIT) {
        gen_prefixes(&root);
    } else {
        /* a_K = M ranges over [ceil(T/2), UB[K]]; M must be <= n_2(A_{K-1})+1 which the DFS enforces via hi */
        for (int M = (T + 1) / 2; M <= UB[K]; M++) {
            node rm = root;
            rm.a[K] = M;
            bits_set(&rm.E, M);
            bits_set(&rm.S, M); if (M + 1 <= T) bits_set(&rm.S, M + 1); if (2 * M <= T) bits_set(&rm.S, 2 * M);
            gen_prefixes(&rm);
        }
    }
    double t1 = now();
    if (verbose) { printf("k=%d T=%d msplit=%d: %zu prefixes at depth %d (%.2fs)\n", K, T, MSPLIT, nprefixes, SPLIT, t1 - t0); fflush(stdout); }
    atomic_store(&next_prefix, 0);
    tstats = calloc(NTHREADS, sizeof(stats));
    pthread_t th[256];
    for (int i = 0; i < NTHREADS; i++) pthread_create(&th[i], NULL, worker, &tstats[i]);
    for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
    double t2 = now();
    stats tot; memset(&tot, 0, sizeof tot);
    for (int i = 0; i < NTHREADS; i++) {
        tot.found += tstats[i].found;
        for (int j = 0; j <= K + 1; j++) { tot.nodes[j] += tstats[i].nodes[j]; tot.pruned_gap[j] += tstats[i].pruned_gap[j]; }
    }
    uint64_t total = 0; for (int j = 0; j <= K + 1; j++) total += tot.nodes[j];
    if (verbose) {
        printf("nodes per level:");
        for (int j = 2; j <= K; j++) if (tot.nodes[j]) printf(" L%d=%llu", j, (unsigned long long)tot.nodes[j]);
        printf("\n");
        printf("gap-pruned per level:");
        for (int j = 2; j <= K; j++) if (tot.pruned_gap[j]) printf(" L%d=%llu", j, (unsigned long long)tot.pruned_gap[j]);
        printf("\n");
        printf("RESULT k=%d T=%d found=%llu total_nodes=%llu time=%.2fs rate=%.2fM nodes/s threads=%d\n",
               K, T, (unsigned long long)tot.found, (unsigned long long)total, t2 - t0, total / (t2 - t0) / 1e6, NTHREADS);
        fflush(stdout);
    }
    free(tstats);
    return tot.found;
}

static void selftest(void) {
    /* expected number of extremal bases, from Challis-Robinson / Kohonen tables */
    struct { int k, n, count; } cases[] = {
        {3, 8, 1}, {4, 12, 1}, {5, 16, 1}, {6, 20, 5}, {7, 26, 3}, {8, 32, 2}, {9, 40, 1},
        {10, 46, 2}, {11, 54, 4}, {12, 64, 1}, {13, 72, 1}, {14, 80, 3}, {15, 92, 1}, {16, 104, 1}
    };
    int ok = 1;
    int saveq = QUIET; QUIET = 1;
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++) {
        for (int ms = 0; ms <= 1; ms++) {
            uint64_t f = run(cases[i].k, cases[i].n, ms, 0);
            uint64_t g = run(cases[i].k, cases[i].n + 1, ms, 0);
            int good = (f == (uint64_t)cases[i].count && g == 0);
            printf("selftest k=%2d n=%3d msplit=%d: found %llu (expect %d), found at n+1: %llu (expect 0)  %s\n",
                   cases[i].k, cases[i].n, ms, (unsigned long long)f, cases[i].count, (unsigned long long)g, good ? "OK" : "FAIL");
            if (!good) ok = 0;
        }
    }
    QUIET = saveq;
    printf(ok ? "SELFTEST PASSED\n" : "SELFTEST FAILED\n");
    if (!ok) exit(1);
}

int main(int argc, char **argv) {
    int k = 0, t = 0, st = 0, ms = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-k")) k = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t")) t = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-j")) NTHREADS = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-d")) SPLIT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-M")) ms = 1;
        else if (!strcmp(argv[i], "-q")) QUIET = 1;
        else if (!strcmp(argv[i], "-selftest")) st = 1;
        else { fprintf(stderr, "bad arg %s\n", argv[i]); return 1; }
    }
    if (NTHREADS < 1) NTHREADS = 1; if (NTHREADS > 256) NTHREADS = 256;
    if (st) { selftest(); return 0; }
    if (k < 3 || t < 1) { fprintf(stderr, "usage: psp2 -k K -t T [-j threads] [-d splitdepth] [-M] [-q] | -selftest\n"); return 1; }
    run(k, t, ms, 1);
    return 0;
}
