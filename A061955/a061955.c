// a061955.c -- fast search for OEIS A061955: n such that n | V(n), where V(n) is
// the base-2 concatenation rev(n) rev(n-1) ... rev(1), rev(k) = binary of k
// written least-significant digit first (full width, trailing zeros kept).
//
// V(n) = sum_{k=1..n} rev(k) * 2^{S(k-1)},  S(m) = total bit length of 1..m.
// V(n) is odd (its last bit is rev(1) = 1), so only odd n can be terms.
//
// Group k by bit length dp: k = 2^{dp-1} + m, 0 <= m < M.  Then
//   S(k-1) = A_dp + dp*m,   A_dp = (dp-2) 2^{dp-1} + 1,
//   rev(k) = 1 + 2 rev_{dp-1}(m)            (rev_j = j-bit reversal).
// With x = 2^dp define, for m < 2^j,
//   G_j = sum x^m,   R_j = sum rev_j(m) x^m,   y_j = x^{2^j}.
// Splitting m on its top bit gives
//   G_{j+1} = G_j (1 + y_j),  R_{j+1} = 2 R_j (1 + y_j) + y_j G_j,  y_{j+1} = y_j^2,
// so a full block (M = 2^{dp-1}) contributes 2^{A_dp} (G_{dp-1} + 2 R_{dp-1}) after
// dp-1 steps of 3 modular multiplications, and 2^{A_{dp+1}} = 2^{A_dp} * y_{dp-1}.
// The last, partial block is split into aligned segments along the set bits j
// of M: segment [P, P + 2^j) contributes x^P (c_j G_j + 2^{dp-j} R_j) with
// c_j = rev_dp(2^{dp-1} + P).  Total cost per n: about 1.5 log2(n)^2 Montgomery
// multiplications, versus O(n) for direct evaluation.
//
//   a061955 res  LO HI           print "n residue" (V(n) mod n), odd n in [LO,HI)
//   a061955 resl LO HI           same, via the interleaved multi-lane path
//   a061955 bench LO COUNT       time COUNT odd n starting at LO
//   a061955 search LO HI [THREADS] [CHUNK]
//                                multithreaded search, hits -> stdout + hits.txt
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

#ifndef LANES
#define LANES 4
#endif

static inline u64 inv64(u64 n)          // n^-1 mod 2^64, n odd
{
    u64 x = n;                          // correct to 3 bits
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}
static inline u64 addm(u64 a, u64 b, u64 n) { u64 s = a + b; return s >= n ? s - n : s; }
static inline u64 subm(u64 a, u64 b, u64 n) { return a >= b ? a - b : a - b + n; }

// Montgomery product a*b/2^64 mod n  (n odd < 2^63, a,b < n, ni = n^-1 mod 2^64)
static inline u64 mmul(u64 a, u64 b, u64 n, u64 ni)
{
    u128 t = (u128)a * b;
    u64 m = (u64)t * ni;                // m*n == t (mod 2^64)
    u64 hi = (u64)(t >> 64), mh = (u64)(((u128)m * n) >> 64);
    u64 r = hi - mh;                    // (t - m*n) / 2^64, in (-n, n)
    return hi < mh ? r + n : r;
}

static inline int bitlen(u64 n) { return 64 - __builtin_clzll(n); }

// Add the last block (k = 2^{d-1} .. n) to total.  All values in Montgomery form.
static inline u64 last_block(u64 n, u64 ni, u64 one, int d, u64 powA, u64 total)
{
    u64 p2[66];                         // p2[i] = 2^i
    p2[0] = one;
    for (int i = 1; i <= d; i++) p2[i] = addm(p2[i - 1], p2[i - 1], n);
    u64 M = n - (1ull << (d - 1)) + 1;  // 1 <= M <= 2^{d-1}
    int J = 63 - __builtin_clzll(M);
    u64 s = one;                        // 1 + sum over set bits b of M of 2^{d-1-b}
    for (u64 b = M; b; b &= b - 1)
        s = addm(s, p2[d - 1 - __builtin_ctzll(b)], n);
    // Walk j upward; at each set bit j: acc = acc * y_j + x^0-relative segment sum.
    u64 y = p2[d], G = one, R = 0, acc = 0;
    for (int j = 0;; j++) {
        if (M >> j & 1) {
            s = subm(s, p2[d - 1 - j], n);          // now s = c_j
            u64 T = addm(mmul(s, G, n, ni), mmul(p2[d - j], R, n, ni), n);
            acc = addm(mmul(acc, y, n, ni), T, n);
        }
        if (j == J) break;
        u64 t = mmul(y, G, n, ni), u = mmul(y, R, n, ni);
        G = addm(G, t, n);
        R = addm(R, u, n); R = addm(R, R, n); R = addm(R, t, n);
        y = mmul(y, y, n, ni);
    }
    return addm(total, mmul(powA, acc, n, ni), n);
}

// V(n) mod n in Montgomery form (0 iff n | V(n)), n odd >= 3.
static u64 total_mont(u64 n, u64 ni, u64 one)
{
    int d = bitlen(n);
    u64 total = 0, powA = one, x = one;
    for (int dp = 1; dp < d; dp++) {
        x = addm(x, x, n);              // x = 2^dp
        u64 y = x, G = one, R = 0;
        for (int j = 0; j < dp - 1; j++) {
            u64 t = mmul(y, G, n, ni), u = mmul(y, R, n, ni);
            G = addm(G, t, n);
            R = addm(R, u, n); R = addm(R, R, n); R = addm(R, t, n);
            y = mmul(y, y, n, ni);
        }
        total = addm(total, mmul(powA, addm(G, addm(R, R, n), n), n, ni), n);
        powA = mmul(powA, y, n, ni);    // 2^{A_{dp+1}} = 2^{A_dp} * x^{2^{dp-1}}
    }
    return last_block(n, ni, one, d, powA, total);
}

static u64 residue(u64 n)               // plain V(n) mod n, n odd
{
    if (n == 1) return 0;
    u64 ni = inv64(n), one = (0 - n) % n;
    return mmul(total_mont(n, ni, one), 1, n, ni);
}

// LANES odd n's with common bit length d >= 2, full blocks interleaved for ILP.
// Writes Montgomery-form totals (0 iff divisible).
static void totals_lanes(const u64 *nn, int d, u64 *out)
{
    u64 n[LANES], ni[LANES], one[LANES], total[LANES], powA[LANES], x[LANES];
    for (int l = 0; l < LANES; l++) {
        n[l] = nn[l]; ni[l] = inv64(n[l]); one[l] = (0 - n[l]) % n[l];
        total[l] = 0; powA[l] = one[l]; x[l] = one[l];
    }
    for (int dp = 1; dp < d; dp++) {
        u64 y[LANES], G[LANES], R[LANES];
        for (int l = 0; l < LANES; l++) {
            x[l] = addm(x[l], x[l], n[l]);
            y[l] = x[l]; G[l] = one[l]; R[l] = 0;
        }
        for (int j = 0; j < dp - 1; j++)
            for (int l = 0; l < LANES; l++) {
                u64 t = mmul(y[l], G[l], n[l], ni[l]), u = mmul(y[l], R[l], n[l], ni[l]);
                G[l] = addm(G[l], t, n[l]);
                u64 r = addm(R[l], u, n[l]);
                r = addm(r, r, n[l]);
                R[l] = addm(r, t, n[l]);
                y[l] = mmul(y[l], y[l], n[l], ni[l]);
            }
        for (int l = 0; l < LANES; l++) {
            u64 bs = addm(G[l], addm(R[l], R[l], n[l]), n[l]);
            total[l] = addm(total[l], mmul(powA[l], bs, n[l], ni[l]), n[l]);
            powA[l] = mmul(powA[l], y[l], n[l], ni[l]);
        }
    }
    for (int l = 0; l < LANES; l++)
        out[l] = last_block(n[l], ni[l], one[l], d, powA[l], total[l]);
}

// ---------------------------------------------------------------- search driver

static pthread_mutex_t g_mx = PTHREAD_MUTEX_INITIALIZER;
static FILE *g_hits;
static u64 g_lo, g_hi, g_chunk, g_nchunks;
static atomic_ullong g_next;
static atomic_uchar *g_done;
static atomic_ullong g_tested;

static void report(u64 n)
{
    pthread_mutex_lock(&g_mx);
    printf("HIT %llu\n", (unsigned long long)n);
    fflush(stdout);
    if (g_hits) { fprintf(g_hits, "%llu\n", (unsigned long long)n); fflush(g_hits); }
    pthread_mutex_unlock(&g_mx);
}

// Test every odd n in [a, b).
static void scan(u64 a, u64 b)
{
    u64 n = a | 1, buf[LANES], out[LANES];
    while (n < b) {
        int d = bitlen(n);
        u64 last = n + 2 * (LANES - 1);
        if (d >= 2 && last < b && bitlen(last) == d) {
            for (int l = 0; l < LANES; l++) buf[l] = n + 2 * l;
            totals_lanes(buf, d, out);
            for (int l = 0; l < LANES; l++)
                if (out[l] == 0) report(buf[l]);
            n += 2 * LANES;
        } else {
            if (residue(n) == 0) report(n);
            n += 2;
        }
    }
}

static void *worker(void *arg)
{
    (void)arg;
    for (;;) {
        u64 c = atomic_fetch_add(&g_next, 1);
        if (c >= g_nchunks) break;
        u64 a = g_lo + c * g_chunk, b = a + g_chunk;
        if (b > g_hi || b < a) b = g_hi;
        scan(a, b);
        atomic_fetch_add(&g_tested, (b - a) / 2);
        atomic_store(&g_done[c], 1);
    }
    return 0;
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int search(u64 lo, u64 hi, int nthreads, u64 chunk)
{
    g_lo = lo; g_hi = hi; g_chunk = chunk;
    g_nchunks = (hi - lo + chunk - 1) / chunk;
    g_done = calloc(g_nchunks, 1);
    g_hits = fopen("hits.txt", "a");
    fprintf(stderr, "search [%llu, %llu)  threads=%d  chunk=%llu  chunks=%llu\n",
            (unsigned long long)lo, (unsigned long long)hi, nthreads,
            (unsigned long long)chunk, (unsigned long long)g_nchunks);
    pthread_t th[256];
    for (int i = 0; i < nthreads; i++) pthread_create(&th[i], 0, worker, 0);
    double t0 = now(), tlast = t0;
    u64 contig = 0;                     // chunks [0, contig) all done
    for (;;) {
        sleep(1);
        while (contig < g_nchunks && atomic_load(&g_done[contig])) contig++;
        double t = now();
        if (contig == g_nchunks || t - tlast >= 30) {
            u64 upto = lo + contig * chunk;
            if (upto > hi || contig == g_nchunks) upto = hi;
            double rate = atomic_load(&g_tested) / (t - t0);
            fprintf(stderr, "[%8.0fs] complete below %llu  (%.3g odd n/s)\n", t - t0,
                    (unsigned long long)upto, rate);
            FILE *f = fopen("progress.txt", "w");
            if (f) { fprintf(f, "%llu\n", (unsigned long long)upto); fclose(f); }
            tlast = t;
        }
        if (contig == g_nchunks) break;
    }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], 0);
    fprintf(stderr, "done [%llu, %llu) in %.1fs\n", (unsigned long long)lo,
            (unsigned long long)hi, now() - t0);
    if (g_hits) fclose(g_hits);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc >= 4 && !strcmp(argv[1], "res")) {
        u64 lo = strtoull(argv[2], 0, 0), hi = strtoull(argv[3], 0, 0);
        for (u64 n = lo | 1; n < hi; n += 2)
            printf("%llu %llu\n", (unsigned long long)n, (unsigned long long)residue(n));
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "resl")) {
        u64 lo = strtoull(argv[2], 0, 0), hi = strtoull(argv[3], 0, 0), buf[LANES], out[LANES];
        u64 n = lo | 1;
        while (n < hi) {
            u64 last = n + 2 * (LANES - 1);
            if (bitlen(n) >= 2 && last < hi && bitlen(last) == bitlen(n)) {
                for (int l = 0; l < LANES; l++) buf[l] = n + 2 * l;
                totals_lanes(buf, bitlen(n), out);
                for (int l = 0; l < LANES; l++)
                    printf("%llu %llu\n", (unsigned long long)buf[l],
                           (unsigned long long)mmul(out[l], 1, buf[l], inv64(buf[l])));
                n += 2 * LANES;
            } else {
                printf("%llu %llu\n", (unsigned long long)n, (unsigned long long)residue(n));
                n += 2;
            }
        }
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "bench")) {
        u64 lo = strtoull(argv[2], 0, 0) | 1, cnt = strtoull(argv[3], 0, 0);
        double t0 = now();
        u64 acc = 0;
        for (u64 i = 0; i < cnt; i++) acc += residue(lo + 2 * i) == 0;
        double t1 = now();
        u64 buf[LANES], out[LANES];
        for (u64 i = 0; i + LANES <= cnt; i += LANES) {
            for (int l = 0; l < LANES; l++) buf[l] = lo + 2 * (i + l);
            totals_lanes(buf, bitlen(buf[0]), out);
            for (int l = 0; l < LANES; l++) acc += out[l] == 0;
        }
        double t2 = now();
        printf("bitlen %d: scalar %.3f us/n, %d-lane %.3f us/n  (hits %llu)\n", bitlen(lo),
               1e6 * (t1 - t0) / cnt, LANES, 1e6 * (t2 - t1) / cnt, (unsigned long long)acc);
        return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "search")) {
        u64 lo = strtoull(argv[2], 0, 0), hi = strtoull(argv[3], 0, 0);
        int th = argc > 4 ? atoi(argv[4]) : (int)sysconf(_SC_NPROCESSORS_ONLN);
        u64 chunk = argc > 5 ? strtoull(argv[5], 0, 0) : (1ull << 22);
        if (lo < 1) lo = 1;
        return search(lo, hi, th, chunk);
    }
    fprintf(stderr, "usage: %s res|resl LO HI | bench LO COUNT | search LO HI [THREADS] [CHUNK]\n",
            argv[0]);
    return 1;
}
