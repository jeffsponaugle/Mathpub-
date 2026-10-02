/*
 * psph.c -- exhaustive search for extremal h-bases A_k
 *           (the "global postage stamp problem"; OEIS A001208..A001216, A053346, A053348, ...)
 *
 * Problem
 * -------
 * A basis is A_k = {a_1 = 1 < a_2 < ... < a_k}.  Its h-range n_h(A_k) is the largest n such
 * that every integer 1..n is a sum of AT MOST h elements of A_k (repetition allowed).
 * n_h(k) = max over all A_k of n_h(A_k) is the extremal h-range.
 *
 * This program finds ALL bases A_k with n_h(A_k) >= TGT.
 *   TGT = best known value      -> lists every extremal basis (verification of the tables)
 *   TGT = best known value + 1  -> an empty result proves n_h(k) = best known value
 *
 * Algorithm (Challis 1993 "H-program" + Challis & Robinson 2010 "difficult target")
 * ----------------------------------------------------------------------------------
 * 1. Depth-first enumeration of the prefix a_2 < a_3 < ... < a_{k-1}.  A prefix A_j is
 *    *admissible* iff it generates 1..a_{j+1}-1, i.e. a_{j+1} <= n_h(A_j) + 1 =: g_j (g_j is the
 *    first gap of A_j).  If a_{j+1} > g_j then n_h(A_k) = n_h(A_j) <= n_h(j) < TGT, so this
 *    prune is exact whenever TGT > n_h(k-1) (always the case for TGT near n_h(k)).
 *    n_h(A_j) is computed with the min-stamps table  cnt_j[x] = min(cnt_{j-1}[x], cnt_j[x-a_j]+1)
 *    (uint8, values > h are "infinite").  The table of level j is a pure function of A_j, so it
 *    is shared by all siblings a_{j+1}; tables are extended LAZILY, in chunks, only as far as a
 *    descendant actually needs them, and the extension of level j stops at its first gap.
 *    The first gap of A_j is <= min(h*a_j, n_h(j)) + 1 where n_h(j) is the known extremal value
 *    for j denominations (built-in table, overridable with -b); if no gap is found below that
 *    limit the program aborts, because a known bound would have been violated.
 * 2. Lower bounds (exact necessary conditions):  n_h(A_k) <= h*a_k, so a_k >= L_k = ceil(TGT/h).
 *    Since a_k <= g_{k-1} <= h*a_{k-1} + 1 we need a_{k-1} >= L_{k-1} = ceil((L_k - 1)/h), and
 *    recursively a_j >= L_j = ceil((L_{j+1} - 1)/h).  Equivalently every prefix must satisfy
 *    g_j >= L_{j+1}; in particular the (k-1)-prefix must have n_h(A_{k-1}) >= ceil(TGT/h) - 1.
 * 3. Leaf: for a surviving (k-1)-prefix the table cnt_{k-1} is extended to TGT and every
 *    a_k in [max(a_{k-1}+1, L_k), g_{k-1}] is tested with the Challis difficult target
 *       X = (C_k - 1) a_k + (C_{k-1} - 1) a_{k-1} + ... + (C_2 - 1) a_2 + (a_2 - 1),
 *       C_k = floor(TGT/a_k),  C_{i-1} = floor(a_i / a_{i-1}).
 *    X < TGT, and X is representable iff  min_c (cnt_{k-1}[X - c a_k] + c) <= h  over
 *    0 <= c <= C_k - 1 (c copies of a_k).  Only c >= cmin = ceil((X - h a_{k-1})/(a_k - a_{k-1}))
 *    can work because cnt_{k-1}[y] >= y / a_{k-1}.  If X is not representable the candidate is
 *    rejected (this kills almost everything).  Survivors go through a second stage of further
 *    hard targets (same test, other remainders) and finally a FULL check of 1..TGT, done as a
 *    streaming block DP  u[x] = min(cnt_{k-1}[x], u[x - a_k] + 1)  with early exit at the
 *    first gap.  A basis that passes the full check is reported together with its exact
 *    n_h(A_k) (independent DP up to h*a_k + 1).
 * 4. Threads: all admissible prefixes down to a split depth d are enumerated into a work list;
 *    worker threads take items with an atomic counter and continue the DFS below.
 *
 * Correctness notes: every prune above is a necessary condition for n_h(A_k) >= TGT, so the
 * search is exhaustive.  All ranges are 64-bit; table entries are 8-bit with saturation.
 *
 * Usage
 * -----
 *   psph -h H -k K -t TGT [-j threads] [-d splitdepth] [-b j:n_h(j)]... [-p progress_sec] [-v]
 *   psph -r H a1 a2 ... ak        compute n_H(A) exactly
 *   psph -selftest [quick|full]   verify against the Challis / Challis-Robinson tables
 *   psph -table H K               print the built-in bound n_H(j) for j < K
 *
 * Build:  clang -O3 -mcpu=apple-m1 -std=c11 -pthread -o psph psph.c        (8-bit cells, h <= 254)
 *         clang -O3 -mcpu=apple-m1 -std=c11 -pthread -DWIDE -o psph16 psph.c (16-bit cells, h <= 65534)
 */
#define _DARWIN_C_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdatomic.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <sys/resource.h>
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#define HAVE_NEON 1
#endif

typedef uint8_t  u8;
typedef int64_t  i64;
typedef uint64_t u64;
/* stamp-count cell: 8-bit by default (h <= 254); compile with -DWIDE for 16-bit cells (h <= 65534, 2x memory) */
#ifdef WIDE
typedef uint16_t cnt_t;
#define CNT_MAX 65535u
#define HMAX 65534
#ifdef HAVE_NEON
typedef uint16x8_t vcnt_t;
#define VL 8
#define vld_c  vld1q_u16
#define vst_c  vst1q_u16
#define vdup_c vdupq_n_u16
#define vqadd_c vqaddq_u16
#define vmin_c vminq_u16
#define vmax_c vmaxq_u16
#define vmaxv_c vmaxvq_u16
#define vceq_c vceqq_u16
#endif
#else
typedef uint8_t cnt_t;
#define CNT_MAX 255u
#define HMAX 254
#ifdef HAVE_NEON
typedef uint8x16_t vcnt_t;
#define VL 16
#define vld_c  vld1q_u8
#define vst_c  vst1q_u8
#define vdup_c vdupq_n_u8
#define vqadd_c vqaddq_u8
#define vmin_c vminq_u8
#define vmax_c vmaxq_u8
#define vmaxv_c vmaxvq_u8
#define vceq_c vceqq_u8
#endif
#endif

#define MAXK    24
#define PAD     64
#define CHUNK   2048          /* max bytes per DP step while hunting for the first gap */
#define INF_I64 ((i64)1 << 60)

static inline i64 imin(i64 a, i64 b) { return a < b ? a : b; }
static inline i64 imax(i64 a, i64 b) { return a > b ? a : b; }
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }
static double cpu_s(void) { struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }
static double proc_cpu_s(void) { struct rusage ru; getrusage(RUSAGE_SELF, &ru); return (double)ru.ru_utime.tv_sec + 1e-6 * (double)ru.ru_utime.tv_usec + (double)ru.ru_stime.tv_sec + 1e-6 * (double)ru.ru_stime.tv_usec; }

/* ===================================================================================== */
/*  Known extremal values n_h(k)  (Challis 1993; Challis & Robinson, JIS 13 (2010) 10.2.3;   */
/*  Challis addendum 2013).  Used ONLY as upper bounds a_{j+1} <= n_h(j)+1 and as table caps. */
/* ===================================================================================== */
static const i64 K3_TAB[] = {0, 3,8,15,26,35,52,69,89,112,146,172,212,259,302,354,418,476,548,633,
    714,805,902,1012,1127,1254,1382,1524,1678,1841,2010,2188,2382,2584,2801,3020,3256,3508,3772,
    4043,4326,4628,4941,5272,5606,5960,6334,6723,7120};                              /* h=1..48 */
static const i64 K4_TAB[] = {0, 4,12,24,44,71,114,165,234,326,427,547,708,873,1094,1383,1650,1935,
    2304,2782,3324,3812,4368,5130,5892,6745,7880,8913,9919,11081,12376,13932,15657,17242,18892,
    21061,23445,25553,27978,31347,33981,36806,39914,43592,47536,51218,54900,59702,63891,69362,
    74348,81303,86751,92199,97836};                                                  /* h=1..54 */
static const i64 K5_TAB[] = {0, 5,16,36,70,126,216,345,512,797,1055,1475,2047,2659,3403,4422,5629,
    6865,8669,10835,12903,15785,18801,22456,26469,31108,36949,42744,49436,57033,66771,75558,86303,
    96852,110253,123954,140688,158389,178811,197293,223580,247194,273443,300747,331461,368894,
    401350,443231,490325,536399,586322,634430,699698,754166,823136,892139,968914,1052562,1150377,
    1236682,1325927,1420882,1547688,1678695,1782370,1888725,2036874,2165553,
    /* addendum h=68..90 */
    2330896,2496702,2653201,2846834,3047485,3250580,3429203,3629795,3864527,4103963,4416370,
    4643287,4975426,5223883,5519971,5796515,6139689,6513282,6912409,7258582,7677138,8029729,
    8525267};                                                                        /* h=1..90 */
/* h=20: the paper gives 45754 with basis {1,17,93,436,2898,6897} (verified by DP); OEIS A001211
   prints 45745, which is smaller than what that basis achieves, so 45754 is used here. */
static const i64 K6_TAB[] = {0, 6,20,52,108,211,388,664,1045,1617,2510,3607,5118,7066,9748,12793,
    17061,22342,28874,36560,45754,57814,72997,87555,106888,129783,156744};           /* h=1..26 */
static const i64 K7_TAB[] = {0, 7,26,70,162,336,638,1137,2001,3191,5047,7820,11568,17178,24466};/*1..14*/
static const i64 K8_TAB[] = {0, 8,32,93,228,524,1007,1911,3485};                    /* h=1..8  */
/* rows (fixed h), k = 1.. */
static const i64 H2_TAB[]  = {0, 2,4,8,12,16,20,26,32,40,46,54,64,72,80,92,104,116,128,140,152,164,180,196,212};
static const i64 H3_TAB[]  = {0, 3,7,15,24,36,52,70,93,121,154,186,225,271,323,385};
static const i64 H4_TAB[]  = {0, 4,10,26,44,70,108,162,228,310,422,550,700};
static const i64 H5_TAB[]  = {0, 5,14,35,71,126,211,336,524,726,1016};
static const i64 H6_TAB[]  = {0, 6,18,52,114,216,388,638,1007,1545};
static const i64 H7_TAB[]  = {0, 7,23,69,165,345,664,1137,1911};
#define NEL(t) ((int)(sizeof(t)/sizeof((t)[0])) - 1)

/* k=3, h >= 23 (Challis-Robinson Appendix; k=3 is completely solved) */
static i64 k3_formula(int h, i64 *A)
{
    static const int c[9][6] = {{3,1,1,0,0,0},{3,1,1,0,0,1},{5,2,1,1,0,1},{5,2,1,1,0,2},{7,3,1,2,0,2},
                                {6,2,2,2,1,2},{8,3,2,3,1,2},{8,3,2,3,1,3},{10,4,2,4,1,3}};
    int t = h / 9, r = h % 9;
    i64 a2 = 6*t + c[r][0];
    i64 a3 = (2*t + c[r][1]) + (2*t + c[r][2]) * a2;
    i64 n  = (4*t + c[r][3]) + (2*t + c[r][4]) * a2 + (3*t + c[r][5]) * a3;
    if (A) { A[0] = 1; A[1] = a2; A[2] = a3; }
    return n;
}

/* k=4, 55 <= h <= 302 (Challis addendum 2013, types A/B/C), h = 12t + r */
static i64 k4_formula(int h, i64 *A)
{
    struct row { int r; char ty; int c[10]; int t0, t1; };
    static const struct row R[] = {
        {0,'A',{2,1,0,1,0,1,-3,0,4,-1},4,5},   {0,'A',{1,0,0,0,0,0,-2,0,1,1},6,11},
        {0,'B',{2,2,-1,3,-1,0,-1,-2,-1,4},12,25},
        {1,'A',{1,0,2,1,1,0,0,0,1,0},5,25},
        {2,'A',{2,1,1,1,1,1,-3,1,4,0},5,6},    {2,'A',{1,0,2,1,1,0,0,0,1,1},7,20},
        {2,'B',{5,3,-1,6,-1,0,0,-2,-1,5},21,25},
        {3,'A',{3,1,2,2,1,1,-1,0,4,0},1,24},
        {4,'A',{3,1,2,2,1,1,-1,0,4,1},2,24},
        {5,'A',{3,1,2,2,1,1,-1,0,4,2},4,24},
        {6,'A',{3,1,2,2,1,1,-1,0,4,3},5,24},
        {7,'A',{7,3,2,5,1,2,-1,0,7,1},2,11},   {7,'A',{8,4,1,7,1,0,0,1,1,5},12,24},
        {8,'A',{7,3,3,5,2,2,-1,1,7,1},1,16},   {8,'A',{8,4,1,7,1,0,0,1,1,6},17,24},
        {9,'A',{7,3,3,5,2,2,-1,1,7,2},1,21},   {9,'A',{8,4,1,7,1,0,0,1,1,7},22,24},
        {10,'A',{7,3,3,5,2,2,-1,1,7,3},4,19},  {10,'C',{11,6,1,10,1,0,0,3,0,7},20,24},
        {11,'A',{10,4,3,7,2,2,0,1,7,3},2,7},   {11,'B',{11,4,2,10,1,2,3,1,1,6},8,22},
        {11,'A',{11,5,2,9,2,0,1,2,1,7},23,23}, {11,'B',{12,5,1,12,1,0,3,0,-1,10},24,24},
    };
    int t = h / 12, r = h % 12;
    for (size_t i = 0; i < sizeof(R)/sizeof(R[0]); i++) {
        const struct row *w = &R[i];
        if (w->r != r || t < w->t0 || t > w->t1) continue;
        const int *c = w->c;
        i64 a2 = 9*t + c[0], a3, a4, n;
        if (w->ty == 'B') a3 = (2*t + c[1]) + (3*t + c[2]) * a2;
        else              a3 = (4*t + c[1]) + (3*t + c[2]) * a2;
        a4 = (7*t + c[3]) + (2*t + c[4]) * a2 + (2*t + c[5]) * a3;
        if (w->ty == 'A')      n = (2*t + c[6]) + (t + c[7]) * a2 + (6*t + c[8]) * a3 + (3*t + c[9]) * a4;
        else if (w->ty == 'B') n = (4*t + c[6]) + (3*t + c[7]) * a2 + (2*t + c[8]) * a3 + (3*t + c[9]) * a4;
        else                   n = (t + c[6]) + (4*t + c[7]) * a2 + (6*t + c[8]) * a3 + (3*t + c[9]) * a4;
        if (A) { A[0] = 1; A[1] = a2; A[2] = a3; A[3] = a4; }
        return n;
    }
    return -1;
}

/* n_h(k) if known (proven), else -1 */
static i64 known_n(int h, int k)
{
    if (h < 1 || k < 1) return -1;
    if (k == 1) return h;
    if (k == 2) return ((i64)h*h + 6*h + 1) / 4;
    if (k == 3) return h <= NEL(K3_TAB) ? K3_TAB[h] : k3_formula(h, NULL);
    if (k == 4) return h <= NEL(K4_TAB) ? K4_TAB[h] : (h <= 302 ? k4_formula(h, NULL) : -1);
    if (k == 5 && h <= NEL(K5_TAB)) return K5_TAB[h];
    if (k == 6 && h <= NEL(K6_TAB)) return K6_TAB[h];
    if (k == 7 && h <= NEL(K7_TAB)) return K7_TAB[h];
    if (k == 8 && h <= NEL(K8_TAB)) return K8_TAB[h];
    if (h == 2 && k <= NEL(H2_TAB)) return H2_TAB[k];
    if (h == 3 && k <= NEL(H3_TAB)) return H3_TAB[k];
    if (h == 4 && k <= NEL(H4_TAB)) return H4_TAB[k];
    if (h == 5 && k <= NEL(H5_TAB)) return H5_TAB[k];
    if (h == 6 && k <= NEL(H6_TAB)) return H6_TAB[k];
    if (h == 7 && k <= NEL(H7_TAB)) return H7_TAB[k];
    return -1;
}

/* ===================================================================================== */
/*  Exact h-range of a basis (independent, simple in-place DP; used for -r and for solutions) */
/* ===================================================================================== */
static void dp_inplace(u8 *T, i64 a, i64 lim)
{   /* T[x] = min(T[x], T[x-a]+1) for x = a..lim, ascending (unbounded knapsack) */
    i64 x = a;
    while (x <= lim) {
        i64 len = imin(a, lim - x + 1), i = 0;
        u8 *out = T + x; const u8 *back = T + x - a;
#ifdef HAVE_NEON
        uint8x16_t one = vdupq_n_u8(1);
        for (; i + 16 <= len; i += 16)
            vst1q_u8(out + i, vminq_u8(vld1q_u8(out + i), vqaddq_u8(vld1q_u8(back + i), one)));
#endif
        for (; i < len; i++) { unsigned v = (unsigned)back[i] + 1u; if (v > 255u) v = 255u; if (v < out[i]) out[i] = (u8)v; }
        x += len;
    }
}

static i64 hrange_exact(const i64 *a, int k, int h)
{
    i64 lim = (i64)h * a[k-1] + 1;
    if (h >= 255) {                      /* 16-bit counts (the 8-bit tables cannot express "> h" here) */
        uint16_t *T = (uint16_t *)malloc(sizeof(uint16_t) * ((size_t)lim + 1));
        if (!T) { fprintf(stderr, "out of memory\n"); exit(2); }
        for (i64 x = 0; x <= lim; x++) T[x] = (uint16_t)(x < 65535 ? x : 65535);
        for (int j = 1; j < k; j++) { i64 aj = a[j]; for (i64 x = aj; x <= lim; x++) { unsigned v = T[x-aj] + 1u; if (v < T[x]) T[x] = (uint16_t)v; } }
        i64 n = 0;
        while (n + 1 <= lim && T[n+1] <= h) n++;
        free(T);
        return n;
    }
    u8 *T = (u8 *)malloc((size_t)lim + 1 + PAD);
    if (!T) { fprintf(stderr, "out of memory\n"); exit(2); }
    for (i64 x = 0; x <= lim; x++) T[x] = (u8)(x < 255 ? x : 255);
    for (int j = 1; j < k; j++) dp_inplace(T, a[j], lim);
    i64 n = 0;
    while (n + 1 <= lim && T[n+1] <= h) n++;
    free(T);
    return n;
}

/* ===================================================================================== */
/*  Search state (one search at a time)                                                    */
/* ===================================================================================== */
static int H, K, D, VERB = 0, PROG = 30, S2N = -1, NOFULL = 0, S2CAP = 64;
static i64 TGT;
static i64 LOW[MAXK + 2];           /* LOW[j]: a_j >= LOW[j]                                  */
static i64 NB[MAXK + 2];            /* NB[j]: n_h(j) if known (bound a_{j+1} <= NB[j]+1), else -1 */
static i64 NB_USER[MAXK + 2];       /* -b overrides                                            */
static cnt_t *T1;                   /* shared table for A_1 = {1}: T1[x] = min(x,CNT_MAX)      */
static int QUIET_SOL = 0;           /* selftest: do not print solutions                        */

typedef struct { i64 a[MAXK + 1]; } item_t;
static item_t *ITEMS; static i64 NITEMS, ITEMCAP;
static atomic_llong NEXT_ITEM, DONE_ITEMS;

typedef struct { int k; i64 n; i64 a[MAXK + 1]; } sol_t;
static sol_t *SOLS; static int NSOLS, SOLCAP;
static pthread_mutex_t SOL_MX = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    int id;
    cnt_t *T[MAXK + 1];             /* T[j]: min-stamps table of prefix A_j (T[1] shared)     */
    i64  P[MAXK + 1];               /* T[j] is valid for indices 0..P[j]                      */
    i64  a[MAXK + 1];
    i64  gap[MAXK + 1];             /* first gap of A_j = n_h(A_j)+1                          */
    int  valid[MAXK + 1];
    cnt_t *U0, *U1;                 /* block buffers for the full check                       */
    i64 *hard; int *hardT; int nhard; /* stage-2 deep holes: y < a with the largest T[y]       */
    i64 *reach;                     /* r[m] = n_m(A_{k-1}), m = 0..H                            */
    /* statistics */
    u64 visited[MAXK + 1], nodes[MAXK + 1], pruned[MAXK + 1], degen[MAXK + 1];
    u64 cand, rpass, xpass, s2pass, fullfail, sols, direct, probes, fc_bytes, tab_bytes;
    u64 fhist[18];
    double t_tab, t_cand, t_busy; u64 items_done;
} ctx_t;

static ctx_t *CTX;

/* ----- DP kernel: out[i] = min(prev[i], sat(back[i]+1)), i < len; if check, return first i with out[i] > h */
static inline i64 dp_chunk(cnt_t *restrict out, const cnt_t *restrict prev, const cnt_t *restrict back, i64 len, int h, int check)
{
    i64 i = 0; int any = 0;
#ifdef HAVE_NEON
    {
        vcnt_t one = vdup_c(1), mx = vdup_c(0);
        for (; i + VL <= len; i += VL) {
            vcnt_t p = vld_c(prev + i);
            vcnt_t b = vqadd_c(vld_c(back + i), one);
            vcnt_t r = vmin_c(p, b);
            vst_c(out + i, r);
            mx = vmax_c(mx, r);
        }
        if (check && vmaxv_c(mx) > h) any = 1;
    }
#endif
    for (; i < len; i++) {
        unsigned v = (unsigned)back[i] + 1u; if (v > CNT_MAX) v = CNT_MAX;
        unsigned p = prev[i];
        cnt_t r = (cnt_t)(p < v ? p : v);
        out[i] = r;
        if (check && r > h) any = 1;
    }
    if (!any) return -1;
    for (i = 0; i < len; i++) if (out[i] > h) return i;
    return -1;
}

/* extend T_j (no gap checking) so that indices 0..upto are valid */
static void extend(ctx_t *c, int j, i64 upto)
{
    if (upto > TGT) upto = TGT;
    if (c->P[j] >= upto) return;
    if (j == 1) { c->P[1] = upto; return; }          /* T1 is precomputed up to TGT */
    extend(c, j - 1, upto);
    const i64 a = c->a[j];
    cnt_t *T = c->T[j]; const cnt_t *Tp = c->T[j - 1];
    i64 x = c->P[j] + 1;
    if (x < a) {
        i64 end = imin(a - 1, upto);
        memcpy(T + x, Tp + x, sizeof(cnt_t) * (size_t)(end - x + 1));
        x = end + 1;
    }
    while (x <= upto) {
        i64 len = imin(a, upto - x + 1);
        dp_chunk(T + x, Tp + x, T + x - a, len, H, 0);
        x += len;
    }
    c->P[j] = upto;
}

/* extend T_j chunk by chunk up to `limit`; return the first x in (P_old, limit] with T_j[x] > H, or 0 */
static i64 find_gap(ctx_t *c, int j, i64 limit)
{
    const i64 a = c->a[j];
    cnt_t *T = c->T[j]; const cnt_t *Tp = c->T[j - 1];
    i64 x = c->P[j] + 1;
    if (limit > TGT) limit = TGT;
    if (x < a && x <= limit) {
        i64 end = imin(a - 1, limit);
        extend(c, j - 1, end);
        memcpy(T + x, Tp + x, sizeof(cnt_t) * (size_t)(end - x + 1));
        c->P[j] = end;
        for (i64 y = x; y <= end; y++) if (T[y] > H) return y;
        x = end + 1;
    }
    while (x <= limit) {
        i64 len = imin(imin(a, CHUNK), limit - x + 1);
        extend(c, j - 1, x + len - 1);
        i64 r = dp_chunk(T + x, Tp + x, T + x - a, len, H, 1);
        c->P[j] = x + len - 1;
        if (r >= 0) return x + r;
        x += len;
    }
    return 0;
}

static void node_enter(ctx_t *c, int j, i64 a)
{
    c->a[j] = a; c->P[j] = 0; c->T[j][0] = 0; c->gap[j] = 0; c->valid[j] = 1;
    for (int t = j + 1; t <= K; t++) c->valid[t] = 0;
}

/* limit below which the first gap of A_j must lie */
static inline i64 gap_limit(int j, i64 a)
{
    i64 lim = (i64)H * a;
    if (NB[j] >= 0 && NB[j] < lim) lim = NB[j];
    lim += 1;
    if (lim > TGT) lim = TGT;
    return lim;
}

static void fatal_bound(ctx_t *c, int j)
{
    fprintf(stderr, "\nFATAL: prefix A_%d = {", j);
    for (int i = 1; i <= j; i++) fprintf(stderr, "%s%lld", i > 1 ? "," : "", (long long)c->a[i]);
    fprintf(stderr, "} has no gap up to %lld, but the bound table says n_%d(%d) = %lld.\n"
            "The bound is wrong (or TGT <= n_h(k-1)); rerun with -b %d:<larger value>.\n",
            (long long)gap_limit(j, c->a[j]), H, j, (long long)NB[j], j);
    exit(3);
}

static void report_solution(ctx_t *c, i64 ak)
{
    i64 A[MAXK + 1];
    for (int i = 1; i < K; i++) A[i-1] = c->a[i];
    A[K-1] = ak;
    i64 n = hrange_exact(A, K, H);
    pthread_mutex_lock(&SOL_MX);
    if (NSOLS == SOLCAP) { SOLCAP = SOLCAP ? 2*SOLCAP : 64; SOLS = (sol_t *)realloc(SOLS, (size_t)SOLCAP * sizeof(sol_t)); }
    SOLS[NSOLS].k = K; SOLS[NSOLS].n = n;
    for (int i = 0; i < K; i++) SOLS[NSOLS].a[i] = A[i];
    NSOLS++;
    if (!QUIET_SOL) {
        printf("SOLUTION h=%d k=%d n_h=%lld :", H, K, (long long)n);
        for (int i = 0; i < K; i++) printf(" %lld", (long long)A[i]);
        printf("\n"); fflush(stdout);
    }
    pthread_mutex_unlock(&SOL_MX);
    c->sols++;
    if (n < TGT) { fprintf(stderr, "INTERNAL ERROR: reported basis has n_h=%lld < TGT=%lld\n", (long long)n, (long long)TGT); exit(4); }
}

/* max of n cells; first index of value v (v must occur) */
static inline int vec_max_c(const cnt_t *p, i64 n)
{
    i64 i = 0; int m = 0;
#ifdef HAVE_NEON
    vcnt_t mx = vdup_c(0);
    for (; i + VL <= n; i += VL) mx = vmax_c(mx, vld_c(p + i));
    m = vmaxv_c(mx);
#endif
    for (; i < n; i++) if (p[i] > m) m = p[i];
    return m;
}
static inline i64 vec_find_c(const cnt_t *p, i64 n, cnt_t v)
{
    i64 i = 0;
#ifdef HAVE_NEON
    vcnt_t vv = vdup_c(v);
    for (; i + VL <= n; i += VL) if (vmaxv_c(vceq_c(vld_c(p + i), vv))) break;
#endif
    for (; i < n; i++) if (p[i] == v) return i;
    return n;
}

/* Full check of 1..TGT for a_k = a.  Streaming block DP  u[x] = min(T[x], u[x-a]+1)  over the range
   [0, E], E = min(TGT, h*a_{k-1}), with early exit at the first gap.  Beyond E every prefix value exceeds
   h, so u[x] = u[x-a] + 1 there, i.e. u[x] = u[z] + j for the unique z in (E-a, E] with x = z + j*a; the
   first gap above E is therefore  min_{z in (E-a, E]} ( z + (h+1-u[z])*a ),  computed in O(a).
   Returns the first gap in [1, TGT], or 0 if there is none. */
static i64 full_check(ctx_t *c, i64 a, i64 E)
{
    const cnt_t *T = c->T[K-1];
    const cnt_t *last = T; i64 start_last = 0, len_last = imin(a, E + 1);   /* block 0 = T[0..a) */
    const cnt_t *prev = NULL;
    cnt_t *cur = c->U0, *other = c->U1;
    for (i64 start = a; start <= E; start += a) {
        i64 len = imin(a, E + 1 - start);
        i64 r = dp_chunk(cur, T + start, last, len, H, 1);
        c->fc_bytes += (u64)len;
        if (r >= 0) return start + r;
        prev = last; last = cur; start_last = start; len_last = len;
        cnt_t *t = cur; cur = other; other = t;
        if (cur == last) cur = other;              /* never overwrite the block we just produced */
    }
    (void)len_last;
    /* Window (E-a, E].  z + (h+1-u[z])*a is minimised at the smallest z carrying the maximum u in the
       window (two z differ by less than a, one unit of u is worth a).  The window is the tail of `prev`
       (positions below start_last) followed by `last`. */
    i64 z0 = imax(E - a + 1, 1);
    int umax = 0; i64 zbest = -1;
    if (z0 < start_last) {                                  /* part 1: prev[z - (start_last - a)] */
        const cnt_t *p = prev + (z0 - (start_last - a)); i64 n = start_last - z0;
        int m = vec_max_c(p, n);
        if (m > umax) { umax = m; zbest = z0 + vec_find_c(p, n, (cnt_t)m); }
    }
    {                                                       /* part 2: last[z - start_last], z in [max(z0,start_last), E] */
        i64 zs = imax(z0, start_last);
        const cnt_t *p = last + (zs - start_last); i64 n = E - zs + 1;
        int m = vec_max_c(p, n);
        if (m > umax || zbest < 0) { umax = m; zbest = zs + vec_find_c(p, n, (cnt_t)m); }
    }
    c->fc_bytes += (u64)(E - z0 + 1);
    i64 best = zbest + ((i64)H + 1 - umax) * a;
    return best <= TGT ? best : 0;
}

/* representable(x) for the full basis with a_k = a, using the prefix table */
static inline int leaf_rep(const cnt_t *T, i64 x, i64 a, i64 ak1, u64 *probes)
{
    i64 Cx = x / a;                                   /* max copies of a */
    i64 num = x - (i64)H * ak1, den = a - ak1;
    i64 cmin = num <= 0 ? 0 : (num + den - 1) / den;
    i64 cmax = imin(Cx, (i64)H);
    for (i64 cc = cmax; cc >= cmin; cc--) {
        (*probes)++;
        if ((i64)T[x - cc * a] + cc <= H) return 1;
    }
    return 0;
}

/* stage-2 deep-hole list: keep the NDH positions y (< current a) with the largest T[y], at most one per
   residue class mod a_{k-1} (positions in the same class share their rescue structure) */
static inline void deep_insert(i64 *dy, int *dT, int *nd, int *dmin, int NDH, i64 y, int t, i64 ak1)
{
    i64 res = y % ak1;
    for (int i = 0; i < *nd; i++)
        if (dy[i] % ak1 == res) {
            if (t > dT[i]) { dy[i] = y; dT[i] = t; if (i == *dmin) { *dmin = 0; for (int j = 1; j < *nd; j++) if (dT[j] < dT[*dmin]) *dmin = j; } }
            return;
        }
    if (*nd < NDH) { dy[*nd] = y; dT[*nd] = t; (*nd)++; }
    else if (t > dT[*dmin]) { dy[*dmin] = y; dT[*dmin] = t; }
    else return;
    *dmin = 0; for (int j = 1; j < *nd; j++) if (dT[j] < dT[*dmin]) *dmin = j;
}

static void leaf(ctx_t *c)
{
    const int km1 = K - 1;
    const i64 ak1 = c->a[km1];
    i64 lo = imax(ak1 + 1, LOW[K]);
    i64 hi = imin(c->gap[km1], TGT);
    if (lo > hi) return;
    double t0 = cpu_s();
    /* probes only ever touch y <= (h-c)*a_{k-1} <= h*a_{k-1}: the table is needed only up to E */
    const i64 E = imin(TGT, (i64)H * ak1);
    i64 p_old = c->P[km1];
    extend(c, km1, E);
    c->tab_bytes += (u64)(c->P[km1] - p_old);
    const cnt_t *T = c->T[km1];
    /* stage-2 "deep holes": the S2N positions y < a with the largest prefix cost T[y]; maintained
       incrementally as a grows (y = a-1 becomes eligible).  Start with all y < lo. */
    /* number of deep holes: a full check costs ~E/16 vector ops, a hole test a few probes */
    const int NDH = (S2N >= 0) ? S2N : (int)imax(8, imin(S2CAP, E / 1024));
    int nd = 0, dmin = 0;
    i64 *dy = c->hard; int *dT = c->hardT;
    if (NDH > 0 && lo > 1) {
        /* only positions within a few levels of the maximum can enter the list: scan level by level */
        int tmax = vec_max_c(T + 1, lo - 1);
        for (int lvl = tmax; lvl >= 0 && lvl >= tmax - 3; lvl--) {
            for (i64 y = 1; y < lo; y++) if (T[y] == lvl) deep_insert(dy, dT, &nd, &dmin, NDH, y, lvl, ak1);
            if (nd >= NDH) break;
        }
    }
    /* reach sequence r[m] = n_m(A_{k-1}) = (first y with T[y] > m) - 1, m = 0..H */
    i64 *r = c->reach;
    const i64 g = c->gap[km1];
    {
        int m = 0;
        for (i64 y = 1; y < g && m <= H; y++) while (m <= H && (int)T[y] > m) { r[m] = y - 1; m++; }
        while (m <= H) { r[m] = g - 1; m++; }
    }
    double t1 = cpu_s(); c->t_tab += t1 - t0;
    /* prefix part of the Challis difficult target */
    i64 Z = (K >= 3) ? c->a[2] - 1 : 0;
    for (int i = 3; i <= km1; i++) Z += (c->a[i] / c->a[i-1] - 1) * c->a[i-1];
    u64 probes = 0;
    int ma = 0;                                       /* ma = min m with r[m] >= a-1 (nondecreasing in a) */
    for (i64 a = lo; a <= hi; a++) {
        c->cand++;
        if (NDH > 0 && a > lo) {                      /* y = a-1 becomes an eligible remainder */
            int t = T[a - 1];
            if (nd < NDH || t > dT[dmin]) deep_insert(dy, dT, &nd, &dmin, NDH, a - 1, t, ak1);
        }
        i64 Ck  = TGT / a;
        /* Stage 1: Challis difficult target X (c copies of a, c <= Ck-1 since X < Ck*a) */
        i64 Ck1 = a / ak1;
        i64 X = (Ck - 1) * a + (Ck1 - 1) * ak1 + Z;
        {
            i64 num = X - (i64)H * ak1, den = a - ak1;
            i64 cmin = num <= 0 ? 0 : (num + den - 1) / den;
            int rep = 0;
            for (i64 cc = Ck - 1; cc >= cmin; cc--) {
                probes++;
                if ((i64)T[X - cc * a] + cc <= H) { rep = 1; break; }
            }
            if (!rep) continue;
        }
        c->xpass++;
        /* Stage 2a: deep holes in the two blocks just below the top block */
        {
            int ok = 1;
            for (int i = 0; i < nd && ok; i++) if (!leaf_rep(T, (Ck - 1) * a + dy[i], a, ak1, &probes)) ok = 0;
            if (ok && Ck >= 3) for (int i = 0; i < nd && ok; i++) if (!leaf_rep(T, (Ck - 2) * a + dy[i], a, ak1, &probes)) ok = 0;
            if (!ok) continue;
        }
        c->rpass++;
        /* Stage 2b: reach test.  With c copies of a the block [c a, (c+1) a) is fully covered iff H-c >= ma
           (ma = min m with r[m] >= a-1); cstar = H-ma+1 is the first block that is not, and
           xs = cstar*a + r[ma-1] + 1 its first number not covered that way.  If cstar > Ck (or xs > TGT)
           all of 1..TGT is covered (sufficient condition; still verified by the full check). */
        while (r[ma] < a - 1) ma++;                   /* ma <= H because r[H] = g-1 >= a-1 */
        {
            i64 cstar = (i64)H - ma + 1;
            if (cstar <= Ck) {
                i64 xs = cstar * a + r[ma - 1] + 1;
                if (xs <= TGT) { if (!leaf_rep(T, xs, a, ak1, &probes)) continue; }
                else c->direct++;
            } else c->direct++;
        }
        /* Stage 2c: deep holes in the top block, and TGT itself */
        {
            int ok = 1;
            for (int i = 0; i < nd && ok; i++) { i64 x = Ck * a + dy[i]; if (x <= TGT && !leaf_rep(T, x, a, ak1, &probes)) ok = 0; }
            if (ok && !leaf_rep(T, TGT, a, ak1, &probes)) ok = 0;
            if (!ok) continue;
        }
        c->s2pass++;
        if (NOFULL) continue;                         /* profiling only: never report */
        i64 g = full_check(c, a, E);
        if (g) {
            c->fullfail++;
            i64 dist = Ck - g / a; if (dist > 17) dist = 17; if (dist < 0) dist = 0;
            c->fhist[dist]++;
            continue;
        }
        report_solution(c, a);
    }
    c->probes += probes;
    c->t_cand += cpu_s() - t1;
}

/* enumerate a_j (j <= K-1) given valid A_{j-1} with gap[j-1]; recurse / leaf */
static void dfs(ctx_t *c, int j)
{
    i64 lo = imax(c->a[j-1] + 1, LOW[j]);
    i64 hi = c->gap[j-1];
    for (i64 a = lo; a <= hi; a++) {
        c->visited[j]++;
        node_enter(c, j, a);
        i64 lim = gap_limit(j, a);
        i64 g = find_gap(c, j, lim);
        if (g == 0) {
            if (lim >= TGT) {                         /* prefix alone already reaches TGT */
                c->degen[j]++;
                if (VERB) { fprintf(stderr, "degenerate: A_%d = {", j); for (int i = 1; i <= j; i++) fprintf(stderr, "%s%lld", i>1?",":"", (long long)c->a[i]); fprintf(stderr, "} has n_h >= TGT\n"); }
                continue;
            }
            fatal_bound(c, j);
        }
        c->gap[j] = g;
        if (g < LOW[j+1]) { c->pruned[j]++; continue; }
        c->nodes[j]++;
        if (j == K - 1) leaf(c); else dfs(c, j + 1);
    }
}

/* ----- work items ----- */
static void push_item(ctx_t *c)
{
    if (NITEMS == ITEMCAP) { ITEMCAP = ITEMCAP ? 2*ITEMCAP : 1024; ITEMS = (item_t *)realloc(ITEMS, (size_t)ITEMCAP * sizeof(item_t)); if (!ITEMS) { fprintf(stderr, "out of memory\n"); exit(2); } }
    for (int i = 0; i <= MAXK; i++) ITEMS[NITEMS].a[i] = (i <= D) ? c->a[i] : 0;
    NITEMS++;
}

static void gen_items(ctx_t *c, int j)
{
    i64 lo = imax(c->a[j-1] + 1, LOW[j]);
    i64 hi = c->gap[j-1];
    for (i64 a = lo; a <= hi; a++) {
        c->visited[j]++;
        node_enter(c, j, a);
        i64 lim = gap_limit(j, a);
        i64 g = find_gap(c, j, lim);
        if (g == 0) { if (lim >= TGT) { c->degen[j]++; continue; } fatal_bound(c, j); }
        c->gap[j] = g;
        if (g < LOW[j+1]) { c->pruned[j]++; continue; }
        c->nodes[j]++;
        if (j == D) push_item(c); else gen_items(c, j + 1);
    }
}

static void setup_item(ctx_t *c, const item_t *it)
{
    for (int j = 2; j <= D; j++) {
        if (c->valid[j] && c->a[j] == it->a[j]) continue;
        node_enter(c, j, it->a[j]);
        i64 g = find_gap(c, j, gap_limit(j, it->a[j]));
        if (g == 0 || g < LOW[j+1]) { fprintf(stderr, "INTERNAL ERROR: work item does not reproduce\n"); exit(5); }
        c->gap[j] = g;
    }
}

static void *worker(void *arg)
{
    ctx_t *c = (ctx_t *)arg;
    double t0 = cpu_s();
    for (;;) {
        i64 i = atomic_fetch_add(&NEXT_ITEM, 1);
        if (i >= NITEMS) break;
        setup_item(c, &ITEMS[i]);
        if (D == K - 1) leaf(c); else dfs(c, D + 1);
        c->items_done++;
        atomic_fetch_add(&DONE_ITEMS, 1);
    }
    c->t_busy += cpu_s() - t0;
    return NULL;
}

static void ctx_init(ctx_t *c, int id)
{
    memset(c, 0, sizeof(*c));
    c->id = id;
    c->T[1] = T1; c->P[1] = TGT; c->a[1] = 1; c->gap[1] = H + 1; c->valid[1] = 1;
    for (int j = 2; j <= K - 1; j++) {
        c->T[j] = (cnt_t *)malloc(sizeof(cnt_t) * ((size_t)TGT + 1 + PAD));
        if (!c->T[j]) { fprintf(stderr, "out of memory (tables)\n"); exit(2); }
        c->T[j][0] = 0;
    }
    i64 maxA = (NB[K-1] >= 0) ? NB[K-1] + 1 : TGT;
    if (maxA > TGT) maxA = TGT;
    c->U0 = (cnt_t *)malloc(sizeof(cnt_t) * ((size_t)maxA + 1 + PAD));
    c->U1 = (cnt_t *)malloc(sizeof(cnt_t) * ((size_t)maxA + 1 + PAD));
    c->hard = (i64 *)malloc(sizeof(i64) * 512);
    c->hardT = (int *)malloc(sizeof(int) * 512);
    c->reach = (i64 *)malloc(sizeof(i64) * (size_t)(H + 2));
    if (!c->U0 || !c->U1 || !c->hard || !c->hardT || !c->reach) { fprintf(stderr, "out of memory (buffers)\n"); exit(2); }
}

static void ctx_free(ctx_t *c)
{
    for (int j = 2; j <= K - 1; j++) free(c->T[j]);
    free(c->U0); free(c->U1); free(c->hard); free(c->hardT); free(c->reach);
}

static void fmt_time(double s, char *buf, size_t n)
{
    if (s < 120) snprintf(buf, n, "%.1fs", s);
    else if (s < 7200) snprintf(buf, n, "%.1fm", s / 60);
    else snprintf(buf, n, "%.2fh", s / 3600);
}

/* ===================================================================================== */
/*  Driver                                                                                  */
/* ===================================================================================== */
typedef struct {
    double wall; u64 visited[MAXK + 1], nodes[MAXK + 1], pruned[MAXK + 1], degen[MAXK + 1];
    u64 cand, rpass, xpass, s2pass, fullfail, direct, probes, fc_bytes, tab_bytes; i64 nitems; int split;
    double t_tab, t_cand, t_busy, cpu, t_gen; u64 fhist[18];
} stats_t;
static stats_t ST;

static void run_search(int h, int k, i64 tgt, int nthr, int split)
{
    if (k < 2 || k > MAXK) { fprintf(stderr, "k must be in 2..%d\n", MAXK); exit(1); }
    if (h < 1 || h > HMAX) { fprintf(stderr, "h must be in 1..%d (compile with -DWIDE for 16-bit stamp counts)\n", HMAX); exit(1); }
    H = h; K = k; TGT = tgt;
    if (TGT < 1) { fprintf(stderr, "TGT must be >= 1\n"); exit(1); }
    double t_start = now_s(), cpu_start = proc_cpu_s();
    memset(&ST, 0, sizeof ST);
    NSOLS = 0; NITEMS = 0;

    for (int j = 1; j <= K; j++) NB[j] = (NB_USER[j] >= 0) ? NB_USER[j] : known_n(H, j);
    LOW[K] = (TGT + H - 1) / H;
    for (int j = K - 1; j >= 2; j--) LOW[j] = imax(2, (LOW[j+1] - 1 + H - 1) / H);
    LOW[1] = 1;
    if (K >= 3 && NB[K-1] >= 0 && TGT <= NB[K-1])
        fprintf(stderr, "WARNING: TGT=%lld <= n_%d(%d)=%lld: the admissibility prune is not exact in this degenerate regime\n",
                (long long)TGT, H, K-1, (long long)NB[K-1]);
    if (VERB || 1) {
        fprintf(stderr, "[psph] h=%d k=%d TGT=%lld threads=%d", H, K, (long long)TGT, nthr);
        i64 kn = known_n(H, K);
        if (kn >= 0) fprintf(stderr, "  (known n_%d(%d)=%lld%s)", H, K, (long long)kn, TGT == kn ? ", listing extremal bases" : TGT == kn + 1 ? ", proof run" : "");
        fprintf(stderr, "\n[psph] bounds a_{j+1} <= n_h(j)+1:");
        for (int j = 1; j <= K - 1; j++) { if (NB[j] >= 0) fprintf(stderr, " n(%d)=%lld", j, (long long)NB[j]); else fprintf(stderr, " n(%d)=?", j); }
        fprintf(stderr, "\n[psph] lower bounds a_j >= LOW[j]:");
        for (int j = 2; j <= K; j++) fprintf(stderr, " L%d=%lld", j, (long long)LOW[j]);
        fprintf(stderr, "\n");
    }

    T1 = (cnt_t *)malloc(sizeof(cnt_t) * ((size_t)TGT + 1 + PAD));
    if (!T1) { fprintf(stderr, "out of memory (T1)\n"); exit(2); }
    for (i64 x = 0; x <= TGT; x++) T1[x] = (cnt_t)(x < (i64)CNT_MAX ? x : (i64)CNT_MAX);

    CTX = (ctx_t *)calloc((size_t)nthr, sizeof(ctx_t));
    for (int i = 0; i < nthr; i++) ctx_init(&CTX[i], i);

    /* work items */
    if (K == 2) { D = 1; NITEMS = 0; ctx_t *c = &CTX[0]; push_item(c); }
    else {
        int dmin = 2, dmax = K - 1;
        if (split > 0) { D = split; if (D < dmin) D = dmin; if (D > dmax) D = dmax; gen_items(&CTX[0], 2); }
        else {
            for (D = dmin; ; D++) {
                NITEMS = 0;
                ctx_t *c = &CTX[0];
                for (int j = 2; j <= MAXK; j++) { c->visited[j] = c->nodes[j] = c->pruned[j] = c->degen[j] = 0; }
                gen_items(c, 2);
                if (NITEMS >= (i64)64 * nthr || D >= dmax) break;
            }
        }
        for (int j = 2; j <= K; j++) CTX[0].valid[j] = 0;
    }
    double t_gen = now_s() - t_start; ST.t_gen = t_gen;
    fprintf(stderr, "[psph] split depth %d: %lld work items (generated in %.2fs)\n", D, (long long)NITEMS, t_gen);

    atomic_store(&NEXT_ITEM, 0); atomic_store(&DONE_ITEMS, 0);
    pthread_t *th = (pthread_t *)malloc(sizeof(pthread_t) * (size_t)nthr);
    for (int i = 0; i < nthr; i++) pthread_create(&th[i], NULL, worker, &CTX[i]);

    /* progress */
    double last = now_s();
    for (;;) {
        long long done = atomic_load(&DONE_ITEMS);
        if (done >= NITEMS) break;
        usleep(20000);
        double t = now_s();
        if (PROG > 0 && t - last >= PROG) {
            last = t;
            u64 cand = 0, xp = 0, sols = 0;
            for (int i = 0; i < nthr; i++) { cand += CTX[i].cand; xp += CTX[i].s2pass; sols += CTX[i].sols; }
            done = atomic_load(&DONE_ITEMS);
            double el = t - t_start, frac = NITEMS ? (double)done / (double)NITEMS : 1.0;
            char b1[32], b2[32]; fmt_time(el, b1, sizeof b1); fmt_time(frac > 0 ? el / frac - el : 0, b2, sizeof b2);
            fprintf(stderr, "[psph] %lld/%lld items (%.1f%%) elapsed %s, ETA %s, leaf cand %.3g, full checks %llu, sols %llu\n",
                    done, (long long)NITEMS, 100.0 * frac, b1, b2, (double)cand, (unsigned long long)xp, (unsigned long long)sols);
        }
    }
    for (int i = 0; i < nthr; i++) pthread_join(th[i], NULL);
    free(th);
    if (VERB) for (int i = 0; i < nthr; i++) fprintf(stderr, "[psph] thread %d: %llu items, busy %.2fs (tables %.2fs, candidates %.2fs)\n", i, (unsigned long long)CTX[i].items_done, CTX[i].t_busy, CTX[i].t_tab, CTX[i].t_cand);
    ST.wall = now_s() - t_start; ST.cpu = proc_cpu_s() - cpu_start; ST.nitems = NITEMS; ST.split = D;
    for (int i = 0; i < nthr; i++) {
        ctx_t *c = &CTX[i];
        for (int j = 0; j <= MAXK; j++) { ST.visited[j] += c->visited[j]; ST.nodes[j] += c->nodes[j]; ST.pruned[j] += c->pruned[j]; ST.degen[j] += c->degen[j]; }
        ST.cand += c->cand; ST.rpass += c->rpass; ST.direct += c->direct; ST.xpass += c->xpass; ST.s2pass += c->s2pass; ST.fullfail += c->fullfail; ST.probes += c->probes; ST.fc_bytes += c->fc_bytes; ST.tab_bytes += c->tab_bytes;
        ST.t_tab += c->t_tab; ST.t_cand += c->t_cand; ST.t_busy += c->t_busy;
        for (int j = 0; j < 18; j++) ST.fhist[j] += c->fhist[j];
    }
    for (int i = 0; i < nthr; i++) ctx_free(&CTX[i]);
    free(CTX); CTX = NULL;
    free(T1); T1 = NULL;
}

static void print_stats(void)
{
    fprintf(stderr, "[psph] ---- statistics h=%d k=%d TGT=%lld ----\n", H, K, (long long)TGT);
    fprintf(stderr, "[psph] %-6s %16s %16s %16s %10s\n", "level", "visited", "admissible", "pruned(LOW)", "degenerate");
    for (int j = 2; j <= K - 1; j++)
        fprintf(stderr, "[psph] a_%-4d %16llu %16llu %16llu %10llu\n", j, (unsigned long long)ST.visited[j],
                (unsigned long long)ST.nodes[j], (unsigned long long)ST.pruned[j], (unsigned long long)ST.degen[j]);
    fprintf(stderr, "[psph] leaf candidates a_%d tested: %llu   X-pass: %llu (%.3g%%)   deephole-pass: %llu (%.3g%%)   stage2-pass: %llu   full-check fail: %llu   direct: %llu   solutions: %d\n",
            K, (unsigned long long)ST.cand, (unsigned long long)ST.xpass, ST.cand ? 100.0 * (double)ST.xpass / (double)ST.cand : 0.0,
            (unsigned long long)ST.rpass, ST.cand ? 100.0 * (double)ST.rpass / (double)ST.cand : 0.0,
            (unsigned long long)ST.s2pass, (unsigned long long)ST.fullfail, (unsigned long long)ST.direct, NSOLS);
    fprintf(stderr, "[psph] probes: %llu (%.2f per candidate)\n", (unsigned long long)ST.probes, ST.cand ? (double)ST.probes / (double)ST.cand : 0.0);
    if (VERB && ST.fullfail) {
        fprintf(stderr, "[psph] full-check failure position (blocks of a_k below the top block):");
        for (int j = 0; j < 18; j++) if (ST.fhist[j]) fprintf(stderr, " [%d]=%llu", j, (unsigned long long)ST.fhist[j]);
        fprintf(stderr, "\n");
    }
    double busy = ST.t_busy > 0 ? ST.t_busy : 1e-9;
    fprintf(stderr, "[psph] leaf table bytes %.3g, full-check bytes %.3g\n", (double)ST.tab_bytes, (double)ST.fc_bytes);
    fprintf(stderr, "[psph] wall %.2fs, CPU %.2fs (%.1f cores), items %lld (split %d, gen %.2fs); %.3g candidates per CPU-second\n",
            ST.wall, ST.cpu, ST.cpu / ST.wall, (long long)ST.nitems, ST.split, ST.t_gen, (double)ST.cand / (ST.cpu > 0 ? ST.cpu : 1e-9));
    fprintf(stderr, "[psph] worker CPU split: leaf tables %.1f%%, leaf candidates (X test, stage 2, full checks) %.1f%%, prefix DP+overhead %.1f%%\n",
            100.0 * ST.t_tab / busy, 100.0 * ST.t_cand / busy, 100.0 * (busy - ST.t_tab - ST.t_cand) / busy);
}

/* ===================================================================================== */
/*  Selftest                                                                                */
/* ===================================================================================== */
static const char *SELFTEST_QUICK[] = {
    /* "h k n : basis ; basis ..."  -- all extremal bases from Challis-Robinson 2010 App. A */
    "2 3 8 : 1 3 4", "2 4 12 : 1 3 5 6", "2 5 16 : 1 3 5 7 8",
    "2 6 20 : 1 2 5 8 9 10 ; 1 3 4 8 9 11 ; 1 3 4 9 11 16 ; 1 3 5 6 13 14 ; 1 3 5 7 9 10",
    "2 7 26 : 1 2 5 8 11 12 13 ; 1 3 4 9 10 12 13 ; 1 3 5 7 8 17 18",
    "2 8 32 : 1 2 5 8 11 14 15 16 ; 1 3 5 7 9 10 21 22",
    "3 3 15 : 1 4 5", "3 4 24 : 1 4 7 8", "3 5 36 : 1 4 6 14 15",
    "3 6 52 : 1 3 7 9 19 24 ; 1 4 6 14 17 29", "3 7 70 : 1 4 5 15 18 27 34",
    "3 8 93 : 1 3 6 10 24 26 39 41", "3 9 121 : 1 3 8 9 14 32 36 51 53",
    "4 3 26 : 1 5 8", "4 4 44 : 1 3 11 18", "4 5 70 : 1 3 11 15 32",
    "4 6 108 : 1 4 9 16 38 49 ; 1 5 8 27 29 44",
    "4 7 162 : 1 4 9 24 35 49 51 ; 1 4 10 15 37 50 71 ; 1 5 8 25 31 52 71",
    "4 8 228 : 1 3 8 19 33 39 92 102",
    "5 3 35 : 1 6 7", "5 4 71 : 1 4 12 21 ; 1 5 12 28", "5 5 126 : 1 4 9 31 51",
    "5 6 211 : 1 4 13 24 56 61 ; 1 5 8 33 54 67", "5 7 336 : 1 4 13 24 30 87 106",
    "5 8 524 : 1 6 8 33 48 77 183 236",
    "6 3 52 : 1 7 12", "6 4 114 : 1 4 19 33", "6 5 216 : 1 7 12 43 52", "6 6 388 : 1 7 11 48 83 115",
    "6 7 638 : 1 4 18 31 104 145 170",
    /* k=3 */
    "7 3 69 : 1 8 13", "8 3 89 : 1 9 14", "9 3 112 : 1 9 20", "10 3 146 : 1 10 26",
    "11 3 172 : 1 9 30 ; 1 10 26", "12 3 212 : 1 11 37", "13 3 259 : 1 13 34", "14 3 302 : 1 12 52",
    "15 3 354 : 1 12 52", "16 3 418 : 1 15 54", "17 3 476 : 1 14 61", "18 3 548 : 1 15 80",
    "19 3 633 : 1 18 65", "20 3 714 : 1 17 91", "21 3 805 : 1 17 91", "22 3 902 : 1 19 102 ; 1 20 92",
    /* k=4 */
    "7 4 165 : 1 5 24 37", "8 4 234 : 1 6 25 65", "9 4 326 : 1 5 34 60", "10 4 427 : 1 6 41 67",
    "11 4 547 : 1 7 48 85", "12 4 708 : 1 7 48 126", "13 4 873 : 1 9 56 155", "14 4 1094 : 1 8 61 164",
    "15 4 1383 : 1 12 65 240", "16 4 1650 : 1 11 78 216", "17 4 1935 : 1 11 90 252",
    "18 4 2304 : 1 16 73 338", "19 4 2782 : 1 10 99 360", "20 4 3324 : 1 16 103 488",
    /* k=5 */
    "7 5 345 : 1 8 11 64 102", "8 5 512 : 1 9 15 78 115 ; 1 9 15 80 118", "9 5 797 : 1 9 23 108 181",
    "10 5 1055 : 1 8 27 119 194", "11 5 1475 : 1 10 34 165 270", "12 5 2047 : 1 10 26 195 320",
    "13 5 2659 : 1 13 34 242 409", "14 5 3403 : 1 11 48 278 720", "15 5 4422 : 1 14 50 325 782",
    "16 5 5629 : 1 14 61 381 984", "17 5 6865 : 1 13 67 326 1191", "18 5 8669 : 1 14 75 500 1306",
    "19 5 10835 : 1 14 89 523 1892", "20 5 12903 : 1 14 102 589 1912",
    /* k=6 */
    "7 6 664 : 1 7 12 64 113 193", "8 6 1045 : 1 9 14 65 170 297", "9 6 1617 : 1 6 31 48 256 373",
    "10 6 2510 : 1 9 31 96 366 411", "11 6 3607 : 1 7 41 105 490 815", "12 6 5118 : 1 6 47 120 565 946",
    /* k=7 */
    "7 7 1137 : 1 7 18 62 104 244 259 ; 1 8 13 66 115 254 415",
};
static const char *SELFTEST_FULL[] = {   /* slower cases (minutes each) */
    "6 8 1007 : 1 5 18 29 97 170 219 308",
    "8 7 2001 : 1 6 28 47 127 412 602",
    "9 7 3191 : 1 7 30 86 189 607 920",
    "7 8 1911 : 1 4 17 31 117 209 513 550 ; 1 6 20 41 109 228 509 580",
};

static int parse_case(const char *s, int *h, int *k, i64 *n, i64 bases[][MAXK], int *nb)
{
    const char *p = s; char *e;
    *h = (int)strtol(p, &e, 10); p = e; *k = (int)strtol(p, &e, 10); p = e; *n = strtoll(p, &e, 10); p = e;
    while (*p == ' ') p++;
    if (*p != ':') return 0;
    p++; *nb = 0;
    for (;;) {
        for (int i = 0; i < *k; i++) { bases[*nb][i] = strtoll(p, &e, 10); if (e == p) return 0; p = e; }
        (*nb)++;
        while (*p == ' ') p++;
        if (*p == ';') { p++; continue; }
        return 1;
    }
}

static int cmp_sol(const void *x, const void *y)
{
    const sol_t *a = (const sol_t *)x, *b = (const sol_t *)y;
    for (int i = 0; i < a->k; i++) if (a->a[i] != b->a[i]) return a->a[i] < b->a[i] ? -1 : 1;
    return 0;
}
static int cmp_basis(const void *x, const void *y)
{
    const i64 *a = (const i64 *)x, *b = (const i64 *)y;
    for (int i = 0; i < MAXK; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

static int selftest(int full, int nthr)
{
    int fails = 0, total = 0;
    int saved_quiet = QUIET_SOL, saved_prog = PROG;
    QUIET_SOL = 1; PROG = 0;
    /* 1. formula checks */
    printf("== formula checks ==\n");
    for (int h = 23; h <= 60; h++) {
        i64 A[3]; i64 n = k3_formula(h, A); i64 r = hrange_exact(A, 3, h);
        i64 tab = h <= NEL(K3_TAB) ? K3_TAB[h] : -1;
        int ok = (r == n) && (tab < 0 || tab == n);
        total++; if (!ok) { fails++; printf("FAIL k=3 formula h=%d: formula n=%lld, DP=%lld, table=%lld\n", h, (long long)n, (long long)r, (long long)tab); }
    }
    for (int h = 15; h <= 302; h++) {
        i64 A[4]; i64 n = k4_formula(h, A);
        if (n < 0) { if (h >= 55) { fails++; total++; printf("FAIL k=4 formula has no row for h=%d\n", h); } continue; }
        i64 r = hrange_exact(A, 4, h);
        i64 tab = h <= NEL(K4_TAB) ? K4_TAB[h] : -1;
        int ok = (r == n) && (tab < 0 || tab == n);
        total++; if (!ok) { fails++; printf("FAIL k=4 formula h=%d: formula n=%lld (A=%lld,%lld,%lld), DP=%lld, table=%lld\n", h, (long long)n, (long long)A[1], (long long)A[2], (long long)A[3], (long long)r, (long long)tab); }
    }
    /* table consistency: columns vs rows where they overlap */
    for (int h = 2; h <= 7; h++) for (int k = 3; k <= 8; k++) {
        i64 col = -1;
        if (k == 3) col = K3_TAB[h]; else if (k == 4) col = K4_TAB[h]; else if (k == 5) col = K5_TAB[h];
        else if (k == 6) col = K6_TAB[h]; else if (k == 7) col = K7_TAB[h]; else if (k == 8 && h <= 8) col = K8_TAB[h];
        const i64 *row = h == 2 ? H2_TAB : h == 3 ? H3_TAB : h == 4 ? H4_TAB : h == 5 ? H5_TAB : h == 6 ? H6_TAB : H7_TAB;
        int rn = h == 2 ? NEL(H2_TAB) : h == 3 ? NEL(H3_TAB) : h == 4 ? NEL(H4_TAB) : h == 5 ? NEL(H5_TAB) : h == 6 ? NEL(H6_TAB) : NEL(H7_TAB);
        if (k <= rn && col >= 0 && row[k] != col) { fails++; printf("FAIL table inconsistency h=%d k=%d: column %lld row %lld\n", h, k, (long long)col, (long long)row[k]); }
        total++;
    }
    printf("formula/table checks: %d tests, %d failures\n", total, fails);
    /* 2. search cases */
    const char **cases = SELFTEST_QUICK; int ncases = (int)(sizeof(SELFTEST_QUICK)/sizeof(SELFTEST_QUICK[0]));
    for (int pass = 0; pass < (full ? 2 : 1); pass++) {
        if (pass == 1) { cases = SELFTEST_FULL; ncases = (int)(sizeof(SELFTEST_FULL)/sizeof(SELFTEST_FULL[0])); }
        for (int ci = 0; ci < ncases; ci++) {
            int h, k, nb; i64 n; i64 bases[16][MAXK]; memset(bases, 0, sizeof bases);
            if (!parse_case(cases[ci], &h, &k, &n, bases, &nb)) { printf("bad case string: %s\n", cases[ci]); fails++; continue; }
            /* verify the listed bases first */
            int ok = 1;
            for (int b = 0; b < nb; b++) { i64 r = hrange_exact(bases[b], k, h); if (r != n) { ok = 0; printf("FAIL listed basis h=%d k=%d #%d has range %lld != %lld\n", h, k, b, (long long)r, (long long)n); } }
            double t0 = now_s();
            run_search(h, k, n, nthr, 0);
            double t1 = now_s();
            /* compare sets */
            qsort(SOLS, (size_t)NSOLS, sizeof(sol_t), cmp_sol);
            qsort(bases, (size_t)nb, sizeof(bases[0]), cmp_basis);
            int same = (NSOLS == nb);
            for (int b = 0; same && b < nb; b++) { for (int i = 0; i < k; i++) if (SOLS[b].a[i] != bases[b][i]) same = 0; if (SOLS[b].n != n) same = 0; }
            if (!same) {
                ok = 0;
                printf("FAIL h=%d k=%d TGT=%lld: expected %d bases, found %d:", h, k, (long long)n, nb, NSOLS);
                for (int b = 0; b < NSOLS; b++) { printf(" {"); for (int i = 0; i < k; i++) printf("%s%lld", i?",":"", (long long)SOLS[b].a[i]); printf("}(n=%lld)", (long long)SOLS[b].n); }
                printf("\n");
            }
            u64 cand1 = ST.cand, xp1 = ST.xpass; double cpu1 = ST.cpu;
            run_search(h, k, n + 1, nthr, 0);
            double t2 = now_s();
            if (NSOLS != 0) { ok = 0; printf("FAIL h=%d k=%d TGT=%lld (known+1): found %d bases!\n", h, k, (long long)n + 1, NSOLS); for (int b = 0; b < NSOLS; b++) { printf("   {"); for (int i = 0; i < k; i++) printf("%s%lld", i?",":"", (long long)SOLS[b].a[i]); printf("} n=%lld\n", (long long)SOLS[b].n); } }
            total++; if (!ok) fails++;
            printf("%s h=%2d k=%d n=%-8lld bases=%d  [TGT=n: %.2fs wall %.2fs cpu, %llu cand, %llu X-pass]  [TGT=n+1: %.2fs wall %.2fs cpu]\n",
                   ok ? "PASS" : "FAIL", h, k, (long long)n, nb, t1 - t0, cpu1, (unsigned long long)cand1, (unsigned long long)xp1, t2 - t1, ST.cpu);
            fflush(stdout);
        }
    }
    printf("== selftest: %d tests, %d failures ==\n", total, fails);
    QUIET_SOL = saved_quiet; PROG = saved_prog;
    return fails;
}

/* ===================================================================================== */
static void usage(void)
{
    fprintf(stderr,
        "psph -- exhaustive search for extremal postage stamp bases (global postage stamp problem)\n"
        "  psph -h H -k K -t TGT [-j threads] [-d splitdepth] [-b j:value]... [-p sec] [-s2 N] [-v]\n"
        "        find all bases A_K with n_H(A_K) >= TGT\n"
        "        -j   worker threads (default 10)\n"
        "        -d   prefix depth of work items (default: auto, smallest depth giving >= 64*threads items)\n"
        "        -b   override/supply the known extremal value n_H(j) used as bound a_{j+1} <= n_H(j)+1\n"
        "        -p   progress interval in seconds (default 30, 0 = off)\n"
        "        -s2  number of stage-2 deep-hole remainders (default: auto = clamp(h*a_{k-1}/1024, 8, s2cap), 0 = off)\n"
        "        -s2cap  upper limit of the automatic -s2 choice (default 64, max 512)\n"
        "  psph -r H a1 a2 ... ak        exact n_H(A)\n"
        "  psph -selftest [quick|full]   verify against the Challis / Challis-Robinson tables\n"
        "  psph -table H K               print built-in bounds n_H(j), j < K\n");
}

int main(int argc, char **argv)
{
    int h = 0, k = 0, nthr = 10, split = 0;
    i64 tgt = 0;
    for (int j = 0; j <= MAXK + 1; j++) NB_USER[j] = -1;
    if (argc < 2) { usage(); return 1; }
    if (!strcmp(argv[1], "-r")) {
        if (argc < 4) { usage(); return 1; }
        int hh = atoi(argv[2]); int kk = argc - 3; i64 A[MAXK + 1];
        if (kk > MAXK) { fprintf(stderr, "too many denominations\n"); return 1; }
        for (int i = 0; i < kk; i++) A[i] = atoll(argv[3 + i]);
        for (int i = 1; i < kk; i++) if (A[i] <= A[i-1]) { fprintf(stderr, "basis must be strictly increasing\n"); return 1; }
        if (A[0] != 1) { fprintf(stderr, "a_1 must be 1\n"); return 1; }
        i64 n = hrange_exact(A, kk, hh);
        printf("n_%d(A) = %lld   (A = {", hh, (long long)n);
        for (int i = 0; i < kk; i++) printf("%s%lld", i?",":"", (long long)A[i]);
        printf("}, h*a_k = %lld)\n", (long long)hh * A[kk-1]);
        i64 kn = known_n(hh, kk);
        if (kn >= 0) printf("known extremal n_%d(%d) = %lld%s\n", hh, kk, (long long)kn, n > kn ? "  ** EXCEEDED **" : n == kn ? "  (basis is extremal)" : "");
        return 0;
    }
    if (!strcmp(argv[1], "-selftest")) {
        int full = (argc > 2 && !strcmp(argv[2], "full"));
        for (int i = 2; i + 1 < argc; i++) if (!strcmp(argv[i], "-j")) nthr = atoi(argv[i+1]);
        return selftest(full, nthr) ? 1 : 0;
    }
    if (!strcmp(argv[1], "-table")) {
        if (argc < 4) { usage(); return 1; }
        int hh = atoi(argv[2]), kk = atoi(argv[3]);
        for (int j = 1; j < kk; j++) { i64 v = known_n(hh, j); printf("n_%d(%d) = %lld%s\n", hh, j, (long long)v, v < 0 ? " (unknown)" : ""); }
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-h") && i + 1 < argc) h = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) k = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-t") && i + 1 < argc) tgt = atoll(argv[++i]);
        else if (!strcmp(argv[i], "-j") && i + 1 < argc) nthr = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) split = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) PROG = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s2") && i + 1 < argc) S2N = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-v")) VERB++;
        else if (!strcmp(argv[i], "-nofull")) NOFULL = 1;   /* profiling: skip the full check (results invalid) */
        else if (!strcmp(argv[i], "-s2cap") && i + 1 < argc) { S2CAP = atoi(argv[++i]); if (S2CAP < 8) S2CAP = 8; if (S2CAP > 512) S2CAP = 512; }
        else if (!strcmp(argv[i], "-b") && i + 1 < argc) {
            int j = 0; long long v = 0;
            if (sscanf(argv[++i], "%d:%lld", &j, &v) != 2 || j < 1 || j > MAXK) { fprintf(stderr, "bad -b argument (want j:value)\n"); return 1; }
            NB_USER[j] = v;
        }
        else if (!strcmp(argv[i], "-help") || !strcmp(argv[i], "--help")) { usage(); return 0; }
        else { fprintf(stderr, "unknown argument %s\n", argv[i]); usage(); return 1; }
    }
    if (h <= 0 || k <= 0) { usage(); return 1; }
    if (nthr < 1) nthr = 1; if (nthr > 256) nthr = 256;
    if (S2N > 512) S2N = 512;
    if (tgt <= 0) { i64 kn = known_n(h, k); if (kn < 0) { fprintf(stderr, "-t TGT required (no known value for h=%d k=%d)\n", h, k); return 1; } tgt = kn; fprintf(stderr, "[psph] using TGT = known value %lld\n", (long long)tgt); }
    run_search(h, k, tgt, nthr, split);
    print_stats();
    fprintf(stderr, "[psph] %d bases with n_%d(A_%d) >= %lld\n", NSOLS, h, k, (long long)tgt);
    return 0;
}
