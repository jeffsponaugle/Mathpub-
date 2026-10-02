/*
 * gpu_leaf.cu -- CUDA leaf for psph.c (exhaustive search for extremal postage stamp bases)
 *
 * The CPU program (psph.c) enumerates admissible prefixes A_{k-1} = {1 = a_1 < ... < a_{k-1}}; for every
 * surviving prefix the "leaf" tests all candidates a = a_k in [lo, hi].  This file is that leaf on the GPU.
 *
 * Per CPU thread, prefixes are accumulated in a batch (pinned host memory): the min-stamps table T[0..E]
 * of the prefix (E = min(TGT, h*a_{k-1}); T[x] = exact min number of prefix stamps for x, values > h are
 * "infinite"), the candidate range, a_{k-1}, the Challis prefix constant Z, the reach sequence r[0..h]
 * (r_m = n_m(A_{k-1})) and a list of "deep hole" remainders y (positions with large T[y]).  When the batch is
 * full it is uploaded and run on the thread's own stream:
 *
 *   kernel 1 (filter, one thread per candidate): the Challis difficult target X, the deep-hole targets
 *            (C_k-1)a+y and (C_k-2)a+y, the reach target c*a + r_{m-1} + 1, the top-block deep holes
 *            C_k a + y <= TGT and TGT itself; each is an EXACT representability test
 *               representable(t)  <=>  exists c in [0, t/a] with T[t - c a] + c <= h,
 *            exact because T is the exact table of the prefix and a representation of t with c copies of a
 *            leaves t - c a <= (h-c) a_{k-1} <= E (indices above E are treated as unrepresentable, which is
 *            exact: T[y] >= y/a_{k-1} > h there).  Candidates failing any test are rejected; only an
 *            unrepresentable number <= TGT ever causes a rejection.  Survivors are appended with atomicAdd.
 *   kernel 2 (full check, one thread -- or one warp -- per survivor): the exact streaming block DP
 *            u[x] = min(T[x], u[x-a]+1) over [0, E] (in place in a buffer of a cells, 4 bytes = 4 (or 2)
 *            cells per step with the CUDA SIMD-in-a-word intrinsics), early exit at the first gap, then the
 *            O(a) window formula for the first gap above E -- a line-by-line port of full_check() in psph.c.
 *
 * Survivors of kernel 2 are handed to the report callback; psph.c re-verifies each with an independent DP.
 *
 * Exactness is independent of the CHOICE of targets (any set of numbers <= TGT may be tested), so the GPU
 * chooses the deep holes per prefix (CPU: incrementally per candidate); the X test, the reach test and the
 * full check are replicated exactly.  The candidate range and all bounds come from psph.c unchanged.
 *
 * Statistics mirror psph.c (cand, X-pass, deephole-pass, stage2-pass, full-check fail, direct, probes).
 *
 * Environment knobs (all optional):
 *   PSPH_GPU_BATCH_MB      table bytes per batch (default 64)
 *   PSPH_GPU_BATCH_PREFIXES prefixes per batch (default 8192)
 *   PSPH_GPU_MAX_CAND      candidates per batch (default 1<<26)
 *   PSPH_GPU_SCRATCH_MB    full-check scratch per thread (default 64; grows if one survivor needs more)
 *   PSPH_GPU_FULL          thread | warp | auto   (full-check kernel flavour, default auto)
 *   PSPH_GPU_INT64=1       force the 64-bit arithmetic filter kernel (default: 32-bit when provably safe)
 *   PSPH_GPU_PARANOID=1    re-run every batch on the host with the same (shared) code and compare
 *   PSPH_GPU_VERBOSE=1     per-batch log lines
 */
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <algorithm>
#include <vector>
#include <sys/time.h>
#include "gpu_leaf.h"

typedef gpu_cell_t cell_t;
typedef long long i64;
typedef unsigned long long u64;

#ifdef WIDE
#define CS 2                              /* bytes per cell */
#define CPW 2                             /* cells per 32-bit word */
#define CMASK 0xFFFFu
#define REPL(x) ((unsigned)(x) * 0x00010001u)
#else
#define CS 1
#define CPW 4
#define CMASK 0xFFu
#define REPL(x) ((unsigned)(x) * 0x01010101u)
#endif
#define MAXA 32                           /* >= MAXK+1 of psph.c */
#define TABPAD 32                         /* bytes of zero padding after every table (straddling word loads) */

#define CUDA_CHECK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "[gpu] FATAL %s:%d: %s -> %s\n", __FILE__, __LINE__, #x, cudaGetErrorString(e_)); exit(7); } } while (0)

static double wall_s(void) { struct timeval tv; gettimeofday(&tv, NULL); return (double)tv.tv_sec + 1e-6 * (double)tv.tv_usec; }

/* ===================================================================================== */
/*  SIMD-in-a-word helpers, identical semantics on host and device                           */
/* ===================================================================================== */
__host__ __device__ __forceinline__ unsigned v_addus(unsigned a, unsigned b)
{
#ifdef __CUDA_ARCH__
#ifdef WIDE
    return __vaddus2(a, b);
#else
    return __vaddus4(a, b);
#endif
#else
    unsigned r = 0;
    for (int k = 0; k < CPW; k++) { unsigned s = ((a >> (k*CS*8)) & CMASK) + ((b >> (k*CS*8)) & CMASK); if (s > CMASK) s = CMASK; r |= s << (k*CS*8); }
    return r;
#endif
}
__host__ __device__ __forceinline__ unsigned v_minu(unsigned a, unsigned b)
{
#ifdef __CUDA_ARCH__
#ifdef WIDE
    return __vminu2(a, b);
#else
    return __vminu4(a, b);
#endif
#else
    unsigned r = 0;
    for (int k = 0; k < CPW; k++) { unsigned x = (a >> (k*CS*8)) & CMASK, y = (b >> (k*CS*8)) & CMASK; r |= (x < y ? x : y) << (k*CS*8); }
    return r;
#endif
}
/* per cell: CMASK where a > b, else 0 */
__host__ __device__ __forceinline__ unsigned v_cmpgtu(unsigned a, unsigned b)
{
#ifdef __CUDA_ARCH__
#ifdef WIDE
    return __vcmpgtu2(a, b);
#else
    return __vcmpgtu4(a, b);
#endif
#else
    unsigned r = 0;
    for (int k = 0; k < CPW; k++) { unsigned x = (a >> (k*CS*8)) & CMASK, y = (b >> (k*CS*8)) & CMASK; if (x > y) r |= CMASK << (k*CS*8); }
    return r;
#endif
}
__host__ __device__ __forceinline__ unsigned fshift_r(unsigned lo, unsigned hi, unsigned sh)
{   /* low 32 bits of ((hi:lo) >> sh), 0 <= sh < 32 */
#ifdef __CUDA_ARCH__
    return __funnelshift_r(lo, hi, sh);
#else
    return (unsigned)(((((u64)hi) << 32) | lo) >> (sh & 31));
#endif
}
__host__ __device__ __forceinline__ int ffs_u(unsigned x)
{
#ifdef __CUDA_ARCH__
    return __ffs(x);
#else
    return __builtin_ffs((int)x);
#endif
}
template<typename P> __host__ __device__ __forceinline__ unsigned ld_ro(const P *p)
{
#ifdef __CUDA_ARCH__
    return (unsigned)__ldg(p);
#else
    return (unsigned)*p;
#endif
}
__host__ __device__ __forceinline__ unsigned lowmask(unsigned ncells)   /* mask of the ncells low cells, 1 <= ncells <= CPW */
{
    return ncells >= CPW ? 0xFFFFFFFFu : ((1u << (ncells * CS * 8)) - 1u);
}
__host__ __device__ __forceinline__ i64 imin64(i64 a, i64 b) { return a < b ? a : b; }
__host__ __device__ __forceinline__ i64 imax64(i64 a, i64 b) { return a > b ? a : b; }

/* ===================================================================================== */
/*  Batch data layout (shared by host and device)                                           */
/* ===================================================================================== */
struct gpu_prefix {                       /* 64 bytes */
    u64 tab_off;                          /* offset (cells) of T[0] in the batch table buffer; 16-byte aligned */
    i64 E, lo, hi, ak1, Z;
    unsigned deep_off, ndeep;             /* deep holes deep[deep_off .. +ndeep), ascending */
    unsigned reach_off, pad0;             /* reach[reach_off .. +h+1) */
};

struct gpu_dev_stats { u64 xpass, rpass, s2pass, direct, probes, fc_cells; };

/* ===================================================================================== */
/*  Filter (kernel 1) core: exact port of the candidate loop of leaf() / leaf_rep() in psph.c  */
/*  Returns 0 = X fail, 1 = deep-hole (2a) fail, 2 = reach fail, 3 = top-block/TGT fail, 4 = survivor */
/* ===================================================================================== */
template<typename I>
__host__ __device__ __forceinline__ int rep_probe(const cell_t *__restrict__ T, I E, I x, I m /* = x / a */, I a, I ak1, int H, unsigned &probes)
{
    I num = x - (I)H * ak1, den = a - ak1;
    I cmin = num <= 0 ? (I)0 : (num + den - 1) / den;
    I cmax = m < (I)H ? m : (I)H;
    for (I cc = cmax; cc >= cmin; cc--) {
        probes++;
        I idx = x - cc * a;
        if (idx <= E && (I)ld_ro(T + idx) + cc <= (I)H) return 1;
    }
    return 0;
}

template<typename I>
__host__ __device__ __forceinline__ int filter_one(const gpu_prefix &pm, const cell_t *__restrict__ T, const unsigned *__restrict__ dy, const i64 *__restrict__ r,
                                                   I a, int H, I TGT, unsigned &probes, int &direct)
{
    const I E = (I)pm.E, ak1 = (I)pm.ak1, Z = (I)pm.Z;
    const int nd = (int)pm.ndeep;
    I Ck  = TGT / a;
    I Ck1 = a / ak1;
    I X = (Ck - 1) * a + (Ck1 - 1) * ak1 + Z;
    /* Stage 1: Challis difficult target (c copies of a, c <= Ck-1 since X < Ck a) */
    {
        I num = X - (I)H * ak1, den = a - ak1;
        I cmin = num <= 0 ? (I)0 : (num + den - 1) / den;
        int rep = 0;
        for (I cc = Ck - 1; cc >= cmin; cc--) {
            probes++;
            I idx = X - cc * a;
            if (idx <= E && (I)ld_ro(T + idx) + cc <= (I)H) { rep = 1; break; }
        }
        if (!rep) return 0;
    }
    /* Stage 2a: deep holes in the two blocks below the top block (only remainders y < a; list ascending) */
    for (int i = 0; i < nd; i++) { I y = (I)dy[i]; if (y >= a) break; if (!rep_probe<I>(T, E, (Ck - 1) * a + y, Ck - 1, a, ak1, H, probes)) return 1; }
    if (Ck >= 3) for (int i = 0; i < nd; i++) { I y = (I)dy[i]; if (y >= a) break; if (!rep_probe<I>(T, E, (Ck - 2) * a + y, Ck - 2, a, ak1, H, probes)) return 1; }
    /* Stage 2b: reach test.  ma = min m with r[m] >= a-1 (r[H] = g-1 >= a-1, so ma <= H; the bound is a safety net:
       whatever ma comes out, xs is some number <= TGT and its test is exact) */
    {
        int ma = 0;
        while (ma < H && r[ma] < (i64)a - 1) ma++;
        I cstar = (I)H - ma + 1;
        if (cstar <= Ck) {
            I xs = cstar * a + (I)r[ma > 0 ? ma - 1 : 0] + 1;
            if (xs <= TGT) { if (!rep_probe<I>(T, E, xs, xs / a, a, ak1, H, probes)) return 2; }
            else direct = 1;
        } else direct = 1;
    }
    /* Stage 2c: deep holes in the top block, and TGT itself */
    for (int i = 0; i < nd; i++) { I y = (I)dy[i]; if (y >= a) break; I x = Ck * a + y; if (x <= TGT && !rep_probe<I>(T, E, x, Ck, a, ak1, H, probes)) return 3; }
    if (!rep_probe<I>(T, E, TGT, Ck, a, ak1, H, probes)) return 3;
    return 4;
}

template<typename I>
__global__ void __launch_bounds__(256) k_filter(const gpu_prefix *__restrict__ P, int np, const u64 *__restrict__ cst, u64 total,
                                                const cell_t *__restrict__ tab, const unsigned *__restrict__ deep, const i64 *__restrict__ reach,
                                                int H, i64 TGT, unsigned *__restrict__ surv, unsigned *__restrict__ nsurv, unsigned surv_cap,
                                                gpu_dev_stats *__restrict__ st)
{
    u64 t = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    unsigned probes = 0; int xp = 0, rp = 0, sp = 0, dr = 0;
    if (t < total) {
        int lo = 0, hi = np - 1;                       /* largest p with cst[p] <= t */
        while (lo < hi) { int mid = (lo + hi + 1) >> 1; if (cst[mid] <= t) lo = mid; else hi = mid - 1; }
        const int p = lo;
        const gpu_prefix pm = P[p];
        I a = (I)pm.lo + (I)(t - cst[p]);
        int res = filter_one<I>(pm, tab + pm.tab_off, deep + pm.deep_off, reach + pm.reach_off, a, H, (I)TGT, probes, dr);
        xp = res >= 1; rp = res >= 2; sp = res >= 4;
        if (sp) { unsigned s = atomicAdd(nsurv, 1u); if (s < surv_cap) surv[s] = (unsigned)t; }
    }
    /* warp-reduce the counters (all 32 lanes reach this point: blockDim is a multiple of 32, no early return) */
    u64 pr = probes; unsigned c1 = (unsigned)xp, c2 = (unsigned)rp, c3 = (unsigned)sp, c4 = (unsigned)dr;
    for (int o = 16; o > 0; o >>= 1) {
        pr += __shfl_down_sync(0xffffffffu, pr, o);
        c1 += __shfl_down_sync(0xffffffffu, c1, o); c2 += __shfl_down_sync(0xffffffffu, c2, o);
        c3 += __shfl_down_sync(0xffffffffu, c3, o); c4 += __shfl_down_sync(0xffffffffu, c4, o);
    }
    if ((threadIdx.x & 31) == 0) {
        if (pr) atomicAdd(&st->probes, pr);
        if (c1) atomicAdd(&st->xpass, (u64)c1);
        if (c2) atomicAdd(&st->rpass, (u64)c2);
        if (c3) atomicAdd(&st->s2pass, (u64)c3);
        if (c4) atomicAdd(&st->direct, (u64)c4);
    }
}

/* ===================================================================================== */
/*  Full check (kernel 2) core: exact port of full_check() in psph.c                           */
/*  NL = 1: one thread per survivor;  NL = 32: one warp per survivor (lanes stride over the words of a block).  */
/*  u: word buffer of >= roundup(a*CS, 4) + 16 bytes (16-byte aligned) for this survivor.                 */
/*  Returns the first gap in [1, TGT], or 0 if there is none.                                            */
/* ===================================================================================== */
template<int NL>
__host__ __device__ i64 fullcheck_core(const cell_t *__restrict__ Tc, i64 a, i64 E, int H, i64 TGT, unsigned *__restrict__ u, int lane, u64 *cells)
{
    const unsigned *Tw = reinterpret_cast<const unsigned *>(Tc);
    const unsigned Hrep = REPL(H), ONES = REPL(1);
    int nblocks = 0;
    i64 s_last = 0, len_last = imin64(a, E + 1);
    for (i64 start = a; start <= E; start += a) {
        const i64 len = imin64(a, E + 1 - start);
        const unsigned nw = (unsigned)((len + CPW - 1) / CPW);
        const u64 boff = (u64)start * CS;
        const unsigned q0 = (unsigned)(boff >> 2), sh = (unsigned)(boff & 3) * 8u;
        i64 mygap = -1;
        for (unsigned w = (unsigned)lane; w < nw; w += NL) {
            unsigned tlo = ld_ro(Tw + q0 + w);
            unsigned tv = sh ? fshift_r(tlo, ld_ro(Tw + q0 + w + 1), sh) : tlo;
            unsigned prev = nblocks ? u[w] : ld_ro(Tw + w);
            unsigned v = v_minu(tv, v_addus(prev, ONES));
            unsigned gt = v_cmpgtu(v, Hrep);
            i64 nvalid = len - (i64)w * CPW;
            if (nvalid < CPW) { unsigned lm = lowmask((unsigned)nvalid); gt &= lm; v = (v & lm) | (prev & ~lm); }
            u[w] = v;
            if (gt) { int b = ffs_u(gt) - 1; i64 pos = start + (i64)w * CPW + b / (8 * CS); if (mygap < 0 || pos < mygap) mygap = pos; }
        }
        if (lane == 0) *cells += (u64)len;
        nblocks++; s_last = start; len_last = len;
#ifdef __CUDA_ARCH__
        if (NL > 1) {
            i64 g = mygap < 0 ? (i64)0x7fffffffffffffffLL : mygap;
            for (int o = 16; o > 0; o >>= 1) { i64 og = __shfl_xor_sync(0xffffffffu, g, o); if (og < g) g = og; }
            if (g != (i64)0x7fffffffffffffffLL) return g;
            __syncwarp();
        } else
#endif
        if (mygap >= 0) return mygap;
    }
    /* Window (E-a, E]: the first gap above E is min_z z + (H+1-u[z])*a, attained at the smallest z carrying the
       maximum u in the window.  The buffer holds, for offset i, the last block value at i (i < len_last) and the
       value of the block before (i >= len_last) -- or T itself when only one block / no block was computed. */
    int umax = -1; i64 zbest = -1;
    if (nblocks == 0) {                                     /* a > E: window is [1, E] of T */
        for (i64 z = 1 + lane; z <= E; z += NL) { int v = (int)ld_ro(Tc + z); if (v > umax) { umax = v; zbest = z; } }
    } else {
        const cell_t *uc = reinterpret_cast<const cell_t *>(u);
        /* part 1: z in [E-a+1, s_last-1]  <->  i in [len_last, a)  (empty if the last block is full) */
        for (i64 i = len_last + lane; i < a; i += NL) {
            int v = (nblocks >= 2) ? (int)uc[i] : (int)ld_ro(Tc + i);
            if (v > umax) { umax = v; zbest = s_last - a + i; }
        }
        /* part 2: z in [s_last, E]  <->  i in [0, len_last) */
        for (i64 i = lane; i < len_last; i += NL) { int v = (int)uc[i]; if (v > umax) { umax = v; zbest = s_last + i; } }
    }
    if (lane == 0) *cells += (u64)(E - imax64(E - a + 1, 1) + 1);
#ifdef __CUDA_ARCH__
    if (NL > 1) {
        int m = umax;
        for (int o = 16; o > 0; o >>= 1) { int om = __shfl_xor_sync(0xffffffffu, m, o); if (om > m) m = om; }
        i64 z = (umax == m && zbest >= 0) ? zbest : (i64)0x7fffffffffffffffLL;
        for (int o = 16; o > 0; o >>= 1) { i64 oz = __shfl_xor_sync(0xffffffffu, z, o); if (oz < z) z = oz; }
        umax = m; zbest = z;
    }
#endif
    i64 best = zbest + ((i64)H + 1 - umax) * a;
    return best <= TGT ? best : 0;
}

template<int NL>
__global__ void __launch_bounds__(128) k_full(const gpu_prefix *__restrict__ P, const cell_t *__restrict__ tab,
                                              const unsigned *__restrict__ sp, const unsigned *__restrict__ sj, const u64 *__restrict__ off, int n,
                                              unsigned char *__restrict__ scratch, int H, i64 TGT, i64 *__restrict__ res, gpu_dev_stats *__restrict__ st)
{
    int idx, lane;
    if (NL == 1) { idx = blockIdx.x * blockDim.x + threadIdx.x; lane = 0; }
    else { idx = blockIdx.x * (blockDim.x >> 5) + (threadIdx.x >> 5); lane = threadIdx.x & 31; }
    if (idx >= n) return;
    const gpu_prefix pm = P[sp[idx]];
    const i64 a = pm.lo + (i64)sj[idx];
    unsigned *u = reinterpret_cast<unsigned *>(scratch + off[idx]);
    u64 cells = 0;
    i64 g = fullcheck_core<NL>(tab + pm.tab_off, a, pm.E, H, TGT, u, lane, &cells);
    if (lane == 0) { res[idx] = g; atomicAdd(&st->fc_cells, cells); }
}

/* ===================================================================================== */
/*  Host side                                                                               */
/* ===================================================================================== */
struct tbatch {
    int tid, inited;
    cudaStream_t stream; cudaEvent_t ev[4];
    cell_t *h_tab, *d_tab; size_t tab_cap, tab_used;             /* cells */
    gpu_prefix *h_meta, *d_meta; int np, np_cap;
    u64 *h_cst, *d_cst;
    unsigned *h_deep, *d_deep; size_t deep_cap, deep_used;
    i64 *h_reach, *d_reach; size_t reach_cap, reach_used;
    int64_t *h_a;                                                /* np_cap * MAXA prefix copies for reporting */
    u64 total; i64 max_ak1, max_hi;
    unsigned *d_surv, *h_surv; unsigned surv_cap; unsigned *d_nsurv, *h_nsurv;
    unsigned *h_sp, *d_sp, *h_sj, *d_sj; u64 *h_off, *d_off; i64 *h_res, *d_res; size_t chunk_cap;
    unsigned char *d_scratch; size_t scratch_cap;
    gpu_dev_stats *d_st, *h_st;
    unsigned char *h_scratch; size_t h_scratch_cap;              /* paranoid mode */
    gpu_stats_t st;
};

static int G_H, G_K, G_NDEEP_CAP, G_INITED = 0;
static i64 G_TGT;
static gpu_report_fn G_REPORT; static void *G_USER;
static size_t CFG_BATCH_BYTES = (size_t)64 << 20, CFG_SCRATCH = (size_t)64 << 20;
static int CFG_BATCH_PREFIXES = 8192; static u64 CFG_MAX_CAND = (u64)1 << 26;
static int CFG_FULL = 0 /* 0 auto, 1 thread, 2 warp */, CFG_INT64 = 0, CFG_PARANOID = 0, CFG_VERBOSE = 0;
static tbatch *TB = NULL; static int NTB = 0;

static size_t env_size(const char *name, size_t def) { const char *s = getenv(name); if (!s || !*s) return def; double v = atof(s); return v <= 0 ? def : (size_t)v; }

template<typename T> static void host_alloc(T **p, size_t n) { void *v = NULL; CUDA_CHECK(cudaMallocHost(&v, n * sizeof(T) + 64)); *p = (T *)v; }
template<typename T> static void dev_alloc(T **p, size_t n) { void *v = NULL; CUDA_CHECK(cudaMalloc(&v, n * sizeof(T) + 64)); *p = (T *)v; }
template<typename T> static void host_free(T *p) { if (p) CUDA_CHECK(cudaFreeHost(p)); }
template<typename T> static void dev_free(T *p) { if (p) CUDA_CHECK(cudaFree(p)); }

static void tb_init(tbatch *b, int tid)
{
    memset(b, 0, sizeof(*b));
    b->tid = tid; b->inited = 1;
    CUDA_CHECK(cudaStreamCreateWithFlags(&b->stream, cudaStreamNonBlocking));
    for (int i = 0; i < 4; i++) CUDA_CHECK(cudaEventCreate(&b->ev[i]));
    b->tab_cap = CFG_BATCH_BYTES / CS;
    host_alloc(&b->h_tab, b->tab_cap); dev_alloc(&b->d_tab, b->tab_cap);
    b->np_cap = CFG_BATCH_PREFIXES;
    host_alloc(&b->h_meta, (size_t)b->np_cap); dev_alloc(&b->d_meta, (size_t)b->np_cap);
    host_alloc(&b->h_cst, (size_t)b->np_cap + 1); dev_alloc(&b->d_cst, (size_t)b->np_cap + 1);
    b->deep_cap = (size_t)b->np_cap * 16; host_alloc(&b->h_deep, b->deep_cap); dev_alloc(&b->d_deep, b->deep_cap);
    b->reach_cap = (size_t)b->np_cap * 16; host_alloc(&b->h_reach, b->reach_cap); dev_alloc(&b->d_reach, b->reach_cap);
    b->h_a = (int64_t *)malloc(sizeof(int64_t) * (size_t)b->np_cap * MAXA);
    b->surv_cap = 1u << 16; host_alloc(&b->h_surv, b->surv_cap); dev_alloc(&b->d_surv, b->surv_cap);
    host_alloc(&b->h_nsurv, 1); dev_alloc(&b->d_nsurv, 1);
    b->chunk_cap = 1u << 14;
    host_alloc(&b->h_sp, b->chunk_cap); dev_alloc(&b->d_sp, b->chunk_cap);
    host_alloc(&b->h_sj, b->chunk_cap); dev_alloc(&b->d_sj, b->chunk_cap);
    host_alloc(&b->h_off, b->chunk_cap); dev_alloc(&b->d_off, b->chunk_cap);
    host_alloc(&b->h_res, b->chunk_cap); dev_alloc(&b->d_res, b->chunk_cap);
    b->scratch_cap = CFG_SCRATCH; CUDA_CHECK(cudaMalloc((void **)&b->d_scratch, b->scratch_cap));
    host_alloc(&b->h_st, 1); dev_alloc(&b->d_st, 1);
    if (!b->h_a) { fprintf(stderr, "[gpu] out of host memory\n"); exit(7); }
}

static void tb_free(tbatch *b)
{
    if (!b->inited) return;
    if (cudaStreamSynchronize(b->stream) != cudaSuccess) { b->inited = 0; return; }   /* runtime already torn down */
    host_free(b->h_tab); dev_free(b->d_tab); host_free(b->h_meta); dev_free(b->d_meta); host_free(b->h_cst); dev_free(b->d_cst);
    host_free(b->h_deep); dev_free(b->d_deep); host_free(b->h_reach); dev_free(b->d_reach); free(b->h_a);
    host_free(b->h_surv); dev_free(b->d_surv); host_free(b->h_nsurv); dev_free(b->d_nsurv);
    host_free(b->h_sp); dev_free(b->d_sp); host_free(b->h_sj); dev_free(b->d_sj); host_free(b->h_off); dev_free(b->d_off); host_free(b->h_res); dev_free(b->d_res);
    dev_free(b->d_scratch); host_free(b->h_st); dev_free(b->d_st); free(b->h_scratch);
    for (int i = 0; i < 4; i++) CUDA_CHECK(cudaEventDestroy(b->ev[i]));
    CUDA_CHECK(cudaStreamDestroy(b->stream));
    b->inited = 0;
}

static void tb_grow_chunk(tbatch *b, size_t need)
{
    if (need <= b->chunk_cap) return;
    CUDA_CHECK(cudaStreamSynchronize(b->stream));
    host_free(b->h_sp); dev_free(b->d_sp); host_free(b->h_sj); dev_free(b->d_sj); host_free(b->h_off); dev_free(b->d_off); host_free(b->h_res); dev_free(b->d_res);
    b->chunk_cap = need + need / 4;
    host_alloc(&b->h_sp, b->chunk_cap); dev_alloc(&b->d_sp, b->chunk_cap);
    host_alloc(&b->h_sj, b->chunk_cap); dev_alloc(&b->d_sj, b->chunk_cap);
    host_alloc(&b->h_off, b->chunk_cap); dev_alloc(&b->d_off, b->chunk_cap);
    host_alloc(&b->h_res, b->chunk_cap); dev_alloc(&b->d_res, b->chunk_cap);
}

static inline size_t scratch_need(i64 a) { size_t r4 = ((size_t)a * CS + 3) & ~(size_t)3; return (r4 + 16 + 15) & ~(size_t)15; }   /* a cells rounded to a word, + straddle word + slack, 16-byte multiple */

/* ---- host reference (paranoid mode): same code as the kernels, run sequentially on the CPU ---- */
static void paranoid_check(tbatch *b, int use64, const std::vector<unsigned> &gsurv, const std::vector<i64> &gres, const gpu_dev_stats &dst)
{
    gpu_dev_stats hs; memset(&hs, 0, sizeof hs);
    std::vector<unsigned> hsurv;
    for (int p = 0; p < b->np; p++) {
        const gpu_prefix &pm = b->h_meta[p];
        const cell_t *T = b->h_tab + pm.tab_off;
        for (i64 a = pm.lo; a <= pm.hi; a++) {
            unsigned probes = 0; int dr = 0, res;
            if (use64) res = filter_one<i64>(pm, T, b->h_deep + pm.deep_off, b->h_reach + pm.reach_off, a, G_H, G_TGT, probes, dr);
            else        res = filter_one<int>(pm, T, b->h_deep + pm.deep_off, b->h_reach + pm.reach_off, (int)a, G_H, (int)G_TGT, probes, dr);
            hs.probes += probes; hs.xpass += res >= 1; hs.rpass += res >= 2; hs.s2pass += res >= 4; hs.direct += dr != 0;
            if (res >= 4) hsurv.push_back((unsigned)(b->h_cst[p] + (u64)(a - pm.lo)));
        }
    }
    int bad = 0;
    if (hs.probes != dst.probes || hs.xpass != dst.xpass || hs.rpass != dst.rpass || hs.s2pass != dst.s2pass || hs.direct != dst.direct) {
        fprintf(stderr, "[gpu] PARANOID: counter mismatch host/device: probes %llu/%llu xpass %llu/%llu rpass %llu/%llu s2pass %llu/%llu direct %llu/%llu\n",
                hs.probes, dst.probes, hs.xpass, dst.xpass, hs.rpass, dst.rpass, hs.s2pass, dst.s2pass, hs.direct, dst.direct); bad = 1;
    }
    if (hsurv != gsurv) { fprintf(stderr, "[gpu] PARANOID: survivor list mismatch (host %zu, device %zu)\n", hsurv.size(), gsurv.size()); bad = 1; }
    /* full checks of the device's survivors with the host copy of the same core */
    for (size_t s = 0; s < gsurv.size() && !bad; s++) {
        unsigned t = gsurv[s];
        int p = (int)(std::upper_bound(b->h_cst, b->h_cst + b->np, (u64)t) - b->h_cst) - 1;
        const gpu_prefix &pm = b->h_meta[p];
        i64 a = pm.lo + (i64)(t - b->h_cst[p]);
        size_t need = scratch_need(a);
        if (need > b->h_scratch_cap) { free(b->h_scratch); b->h_scratch_cap = need * 2; void *v = NULL; if (posix_memalign(&v, 64, b->h_scratch_cap)) { fprintf(stderr, "[gpu] out of host memory\n"); exit(7); } b->h_scratch = (unsigned char *)v; }
        u64 cells = 0;
        i64 g = fullcheck_core<1>(b->h_tab + pm.tab_off, a, pm.E, G_H, G_TGT, (unsigned *)b->h_scratch, 0, &cells);
        if (g != gres[s]) { fprintf(stderr, "[gpu] PARANOID: full check mismatch for prefix %d a=%lld: host gap %lld, device gap %lld\n", p, a, g, gres[s]); bad = 1; }
    }
    if (bad) { fprintf(stderr, "[gpu] PARANOID: FATAL mismatch in batch %llu of thread %d\n", b->st.batches, b->tid); exit(8); }
}

/* ---- the pipeline for one batch ---- */
static void run_batch(tbatch *b)
{
    if (b->np == 0) return;
    double t0 = wall_s();
    cudaStream_t s = b->stream;
    b->h_cst[b->np] = b->total;
    /* 32-bit arithmetic is exact when every intermediate fits: values are bounded by TGT + max_hi (top-block
       targets) and H*max_ak1 + TGT (|num| in cmin) */
    int use64 = CFG_INT64 || (G_TGT + b->max_hi >= ((i64)1 << 31) - 16) || ((i64)G_H * b->max_ak1 + G_TGT >= ((i64)1 << 31) - 16);

    CUDA_CHECK(cudaEventRecord(b->ev[0], s));
    CUDA_CHECK(cudaMemcpyAsync(b->d_tab, b->h_tab, b->tab_used * CS, cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(b->d_meta, b->h_meta, (size_t)b->np * sizeof(gpu_prefix), cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(b->d_cst, b->h_cst, ((size_t)b->np + 1) * sizeof(u64), cudaMemcpyHostToDevice, s));
    if (b->deep_used) CUDA_CHECK(cudaMemcpyAsync(b->d_deep, b->h_deep, b->deep_used * sizeof(unsigned), cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaMemcpyAsync(b->d_reach, b->h_reach, b->reach_used * sizeof(i64), cudaMemcpyHostToDevice, s));
    CUDA_CHECK(cudaEventRecord(b->ev[1], s));

    /* kernel 1, re-run with a larger survivor list if it overflowed (never lose a survivor) */
    unsigned nsurv = 0;
    for (int attempt = 0; ; attempt++) {
        CUDA_CHECK(cudaMemsetAsync(b->d_nsurv, 0, sizeof(unsigned), s));
        CUDA_CHECK(cudaMemsetAsync(b->d_st, 0, sizeof(gpu_dev_stats), s));
        unsigned blocks = (unsigned)((b->total + 255) / 256);
        if (use64) k_filter<i64><<<blocks, 256, 0, s>>>(b->d_meta, b->np, b->d_cst, b->total, b->d_tab, b->d_deep, b->d_reach, G_H, G_TGT, b->d_surv, b->d_nsurv, b->surv_cap, b->d_st);
        else       k_filter<int><<<blocks, 256, 0, s>>>(b->d_meta, b->np, b->d_cst, b->total, b->d_tab, b->d_deep, b->d_reach, G_H, G_TGT, b->d_surv, b->d_nsurv, b->surv_cap, b->d_st);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaMemcpyAsync(b->h_nsurv, b->d_nsurv, sizeof(unsigned), cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaEventRecord(b->ev[2], s));
        CUDA_CHECK(cudaStreamSynchronize(s));
        b->st.filter_runs++;
        nsurv = *b->h_nsurv;
        if (nsurv <= b->surv_cap) break;
        if (attempt > 4) { fprintf(stderr, "[gpu] FATAL: survivor list keeps overflowing\n"); exit(7); }
        host_free(b->h_surv); dev_free(b->d_surv);
        b->surv_cap = nsurv + nsurv / 4 + 1024;
        host_alloc(&b->h_surv, b->surv_cap); dev_alloc(&b->d_surv, b->surv_cap);
        if (CFG_VERBOSE) fprintf(stderr, "[gpu] thread %d: survivor list overflow (%u), re-running the filter with cap %u\n", b->tid, nsurv, b->surv_cap);
    }
    { float ms; CUDA_CHECK(cudaEventElapsedTime(&ms, b->ev[0], b->ev[1])); b->st.t_h2d += ms * 1e-3; CUDA_CHECK(cudaEventElapsedTime(&ms, b->ev[1], b->ev[2])); b->st.t_filter += ms * 1e-3; }
    CUDA_CHECK(cudaMemcpyAsync(b->h_st, b->d_st, sizeof(gpu_dev_stats), cudaMemcpyDeviceToHost, s));
    if (nsurv) CUDA_CHECK(cudaMemcpyAsync(b->h_surv, b->d_surv, nsurv * sizeof(unsigned), cudaMemcpyDeviceToHost, s));
    CUDA_CHECK(cudaStreamSynchronize(s));
    gpu_dev_stats dst = *b->h_st;
    std::vector<unsigned> surv(b->h_surv, b->h_surv + nsurv);
    std::sort(surv.begin(), surv.end());

    /* kernel 2 over the (sorted) survivors, in chunks that fit the scratch */
    std::vector<i64> gres(nsurv, -1);
    tb_grow_chunk(b, surv.size());                                 /* a chunk never exceeds the survivor count */
    size_t pos = 0; int p = 0;
    while (pos < surv.size()) {
        size_t n = 0, used = 0; i64 amax = 0;
        while (pos + n < surv.size()) {
            unsigned t = surv[pos + n];
            while (p + 1 < b->np && b->h_cst[p + 1] <= t) p++;
            i64 a = b->h_meta[p].lo + (i64)(t - b->h_cst[p]);
            size_t need = scratch_need(a);
            if (used + need > b->scratch_cap) {
                if (n > 0) break;
                CUDA_CHECK(cudaStreamSynchronize(s)); dev_free(b->d_scratch); b->scratch_cap = need * 2;
                CUDA_CHECK(cudaMalloc((void **)&b->d_scratch, b->scratch_cap));
                if (CFG_VERBOSE) fprintf(stderr, "[gpu] thread %d: scratch grown to %zu MB\n", b->tid, b->scratch_cap >> 20);
            }
            b->h_sp[n] = (unsigned)p; b->h_sj[n] = (unsigned)(a - b->h_meta[p].lo); b->h_off[n] = used;
            used += need; if (a > amax) amax = a; n++;
        }
        int warp = (CFG_FULL == 2) || (CFG_FULL == 0 && (n < 2048 || amax > 65536));
        CUDA_CHECK(cudaMemcpyAsync(b->d_sp, b->h_sp, n * sizeof(unsigned), cudaMemcpyHostToDevice, s));
        CUDA_CHECK(cudaMemcpyAsync(b->d_sj, b->h_sj, n * sizeof(unsigned), cudaMemcpyHostToDevice, s));
        CUDA_CHECK(cudaMemcpyAsync(b->d_off, b->h_off, n * sizeof(u64), cudaMemcpyHostToDevice, s));
        CUDA_CHECK(cudaMemsetAsync(b->d_st, 0, sizeof(gpu_dev_stats), s));
        CUDA_CHECK(cudaEventRecord(b->ev[0], s));
        if (warp) k_full<32><<<(unsigned)((n + 3) / 4), 128, 0, s>>>(b->d_meta, b->d_tab, b->d_sp, b->d_sj, b->d_off, (int)n, b->d_scratch, G_H, G_TGT, b->d_res, b->d_st);
        else      k_full<1><<<(unsigned)((n + 127) / 128), 128, 0, s>>>(b->d_meta, b->d_tab, b->d_sp, b->d_sj, b->d_off, (int)n, b->d_scratch, G_H, G_TGT, b->d_res, b->d_st);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaEventRecord(b->ev[1], s));
        CUDA_CHECK(cudaMemcpyAsync(b->h_res, b->d_res, n * sizeof(i64), cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaMemcpyAsync(b->h_st, b->d_st, sizeof(gpu_dev_stats), cudaMemcpyDeviceToHost, s));
        CUDA_CHECK(cudaStreamSynchronize(s));
        { float ms; CUDA_CHECK(cudaEventElapsedTime(&ms, b->ev[0], b->ev[1])); b->st.t_full += ms * 1e-3; }
        b->st.fc_cells += b->h_st->fc_cells;
        for (size_t i = 0; i < n; i++) gres[pos + i] = b->h_res[i];
        pos += n;
    }
    if (CFG_PARANOID) paranoid_check(b, use64, surv, gres, dst);

    /* results: solutions to the callback, failures to the histogram */
    p = 0;
    for (size_t i = 0; i < surv.size(); i++) {
        unsigned t = surv[i];
        while (p + 1 < b->np && b->h_cst[p + 1] <= t) p++;
        const gpu_prefix &pm = b->h_meta[p];
        i64 a = pm.lo + (i64)(t - b->h_cst[p]);
        i64 g = gres[i];
        if (g < 0) { fprintf(stderr, "[gpu] FATAL: missing full-check result\n"); exit(7); }
        if (g == 0) { b->st.sols++; G_REPORT(b->tid, b->h_a + (size_t)p * MAXA, (int64_t)a, G_USER); }
        else {
            b->st.fullfail++;
            i64 Ck = G_TGT / a, dist = Ck - g / a; if (dist > 17) dist = 17; if (dist < 0) dist = 0;
            b->st.fhist[dist]++;
        }
    }
    b->st.cand += b->total; b->st.xpass += dst.xpass; b->st.rpass += dst.rpass; b->st.s2pass += dst.s2pass; b->st.direct += dst.direct; b->st.probes += dst.probes;
    b->st.batches++; b->st.prefixes += (u64)b->np; b->st.tab_bytes += b->tab_used * CS;
    b->st.t_host += wall_s() - t0;
    if (CFG_VERBOSE) fprintf(stderr, "[gpu] thread %d batch %llu: %d prefixes, %.3g MB tables, %llu candidates, X-pass %llu, stage2-pass %llu (%u survivors), sols so far %llu, %s, %.3fs\n",
                             b->tid, b->st.batches, b->np, (double)(b->tab_used * CS) / 1048576.0, b->total, dst.xpass, dst.s2pass, nsurv, b->st.sols, use64 ? "int64" : "int32", wall_s() - t0);
    b->np = 0; b->tab_used = 0; b->deep_used = 0; b->reach_used = 0; b->total = 0; b->max_ak1 = 0; b->max_hi = 0;
}

/* ===================================================================================== */
/*  API                                                                                     */
/* ===================================================================================== */
extern "C" int gpu_init(int h, int k, int64_t tgt, int nthreads, int ndeep_cap, gpu_report_fn report, void *user)
{
    if (!G_INITED) {
        CFG_BATCH_BYTES = env_size("PSPH_GPU_BATCH_MB", 64) << 20;
        CFG_SCRATCH = env_size("PSPH_GPU_SCRATCH_MB", 64) << 20;
        CFG_BATCH_PREFIXES = (int)env_size("PSPH_GPU_BATCH_PREFIXES", 8192);
        CFG_MAX_CAND = (u64)env_size("PSPH_GPU_MAX_CAND", (size_t)1 << 26);
        if (CFG_MAX_CAND > ((u64)1 << 31) - 4096) CFG_MAX_CAND = ((u64)1 << 31) - 4096;
        const char *f = getenv("PSPH_GPU_FULL"); CFG_FULL = (f && !strcmp(f, "thread")) ? 1 : (f && !strcmp(f, "warp")) ? 2 : 0;
        CFG_INT64 = getenv("PSPH_GPU_INT64") && atoi(getenv("PSPH_GPU_INT64"));
        CFG_PARANOID = getenv("PSPH_GPU_PARANOID") && atoi(getenv("PSPH_GPU_PARANOID"));
        CFG_VERBOSE = getenv("PSPH_GPU_VERBOSE") && atoi(getenv("PSPH_GPU_VERBOSE"));
        cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);       /* waiting CPU threads sleep instead of spinning */
        int dev = 0; cudaDeviceProp pr; CUDA_CHECK(cudaGetDevice(&dev)); CUDA_CHECK(cudaGetDeviceProperties(&pr, dev));
        fprintf(stderr, "[gpu] %s (sm_%d%d, %d SMs, %.1f GB), cells %d bytes, batch %zu MB / %d prefixes / %llu candidates, scratch %zu MB, full-check %s%s%s\n",
                pr.name, pr.major, pr.minor, pr.multiProcessorCount, (double)pr.totalGlobalMem / 1073741824.0, CS,
                CFG_BATCH_BYTES >> 20, CFG_BATCH_PREFIXES, CFG_MAX_CAND, CFG_SCRATCH >> 20,
                CFG_FULL == 1 ? "thread" : CFG_FULL == 2 ? "warp" : "auto", CFG_INT64 ? ", int64 forced" : "", CFG_PARANOID ? ", PARANOID host re-check" : "");
        G_INITED = 1;
    }
    if (tgt >= ((int64_t)1 << 31) - 16) { fprintf(stderr, "[gpu] TGT must be < 2^31 in GPU mode\n"); return 1; }
    if (nthreads > NTB) {
        tbatch *nb = (tbatch *)calloc((size_t)nthreads, sizeof(tbatch));
        if (!nb) { fprintf(stderr, "[gpu] out of memory\n"); exit(7); }
        if (TB) { memcpy(nb, TB, sizeof(tbatch) * (size_t)NTB); free(TB); }
        TB = nb; NTB = nthreads;
    }
    G_H = h; G_K = k; G_TGT = (i64)tgt; G_NDEEP_CAP = ndeep_cap; G_REPORT = report; G_USER = user;
    for (int i = 0; i < nthreads; i++) { if (!TB[i].inited) tb_init(&TB[i], i); memset(&TB[i].st, 0, sizeof(gpu_stats_t)); TB[i].np = 0; TB[i].tab_used = 0; TB[i].deep_used = 0; TB[i].reach_used = 0; TB[i].total = 0; TB[i].max_ak1 = 0; TB[i].max_hi = 0; }
    return 0;
}

extern "C" void gpu_submit_prefix(int tid, const gpu_cell_t *T, int64_t E, int64_t lo, int64_t hi, const int64_t *a, int64_t Z,
                                  const int64_t *reach, const int64_t *deep_y, int ndeep)
{
    tbatch *b = &TB[tid];
    if (lo > hi) return;
    if (ndeep > G_NDEEP_CAP) ndeep = G_NDEEP_CAP;
    const size_t cells = (size_t)E + 1;
    const size_t need = ((cells * CS + TABPAD + 15) & ~(size_t)15) / CS;    /* cells incl. padding, 16-byte multiple */
    const u64 ncand = (u64)(hi - lo + 1);
    const int HP1 = G_H + 1;
    if (b->np > 0 && (b->np >= b->np_cap || b->tab_used + need > b->tab_cap || b->total + ncand > CFG_MAX_CAND)) run_batch(b);
    if (need > b->tab_cap) {                                            /* a single table larger than the batch buffer */
        CUDA_CHECK(cudaStreamSynchronize(b->stream));
        host_free(b->h_tab); dev_free(b->d_tab);
        b->tab_cap = need + need / 8;
        host_alloc(&b->h_tab, b->tab_cap); dev_alloc(&b->d_tab, b->tab_cap);
        if (CFG_VERBOSE) fprintf(stderr, "[gpu] thread %d: table buffer grown to %zu MB\n", tid, (b->tab_cap * CS) >> 20);
    }
    if (ncand > CFG_MAX_CAND) {                                        /* candidates are independent: split the range */
        int64_t mid = lo + (int64_t)(CFG_MAX_CAND - 1);
        gpu_submit_prefix(tid, T, E, lo, mid, a, Z, reach, deep_y, ndeep);
        gpu_submit_prefix(tid, T, E, mid + 1, hi, a, Z, reach, deep_y, ndeep);
        return;
    }
    if (b->deep_used + (size_t)ndeep > b->deep_cap) {
        CUDA_CHECK(cudaStreamSynchronize(b->stream));
        size_t nc = (b->deep_used + (size_t)ndeep) * 2; unsigned *nh; host_alloc(&nh, nc); memcpy(nh, b->h_deep, b->deep_used * sizeof(unsigned));
        host_free(b->h_deep); dev_free(b->d_deep); b->h_deep = nh; dev_alloc(&b->d_deep, nc); b->deep_cap = nc;
    }
    if (b->reach_used + (size_t)HP1 > b->reach_cap) {
        CUDA_CHECK(cudaStreamSynchronize(b->stream));
        size_t nc = (b->reach_used + (size_t)HP1) * 2; i64 *nh; host_alloc(&nh, nc); memcpy(nh, b->h_reach, b->reach_used * sizeof(i64));
        host_free(b->h_reach); dev_free(b->d_reach); b->h_reach = nh; dev_alloc(&b->d_reach, nc); b->reach_cap = nc;
    }
    gpu_prefix &pm = b->h_meta[b->np];
    pm.tab_off = b->tab_used; pm.E = (i64)E; pm.lo = (i64)lo; pm.hi = (i64)hi; pm.ak1 = (i64)a[G_K - 1]; pm.Z = (i64)Z;
    pm.deep_off = (unsigned)b->deep_used; pm.ndeep = (unsigned)ndeep; pm.reach_off = (unsigned)b->reach_used; pm.pad0 = 0;
    memcpy(b->h_tab + b->tab_used, T, cells * CS);
    memset(b->h_tab + b->tab_used + cells, 0, (need - cells) * CS);
    b->tab_used += need;
    for (int i = 0; i < ndeep; i++) b->h_deep[b->deep_used + (size_t)i] = (unsigned)deep_y[i];
    b->deep_used += (size_t)ndeep;
    for (int i = 0; i < HP1; i++) b->h_reach[b->reach_used + (size_t)i] = (i64)reach[i];
    b->reach_used += (size_t)HP1;
    int64_t *ac = b->h_a + (size_t)b->np * MAXA;
    for (int i = 0; i < G_K && i < MAXA; i++) ac[i] = a[i];
    b->h_cst[b->np] = b->total;
    b->total += ncand;
    if (pm.ak1 > b->max_ak1) b->max_ak1 = pm.ak1;
    if ((i64)hi > b->max_hi) b->max_hi = (i64)hi;
    b->np++;
}

extern "C" void gpu_flush(int tid) { if (tid >= 0 && tid < NTB && TB[tid].inited) run_batch(&TB[tid]); }

extern "C" void gpu_thread_stats(int tid, gpu_stats_t *out)
{
    if (tid >= 0 && tid < NTB && TB[tid].inited) *out = TB[tid].st; else memset(out, 0, sizeof(*out));
}

extern "C" void gpu_finish(void)
{
    for (int i = 0; i < NTB; i++) if (TB[i].inited) run_batch(&TB[i]);
    CUDA_CHECK(cudaDeviceSynchronize());
}

extern "C" void gpu_print_stats(FILE *f)
{
    gpu_stats_t s; memset(&s, 0, sizeof s);
    for (int i = 0; i < NTB; i++) if (TB[i].inited) {
        const gpu_stats_t &t = TB[i].st;
        s.cand += t.cand; s.xpass += t.xpass; s.rpass += t.rpass; s.s2pass += t.s2pass; s.fullfail += t.fullfail; s.direct += t.direct; s.probes += t.probes; s.sols += t.sols;
        s.batches += t.batches; s.prefixes += t.prefixes; s.tab_bytes += t.tab_bytes; s.fc_cells += t.fc_cells; s.filter_runs += t.filter_runs;
        s.t_h2d += t.t_h2d; s.t_filter += t.t_filter; s.t_full += t.t_full; s.t_host += t.t_host;
    }
    double tg = s.t_h2d + s.t_filter + s.t_full;
    fprintf(f, "[gpu] batches %llu (%llu prefixes, %.3g candidates, %.3g table bytes uploaded, %llu filter launches), full-check cells %.3g\n",
            s.batches, s.prefixes, (double)s.cand, (double)s.tab_bytes, s.filter_runs, (double)s.fc_cells);
    fprintf(f, "[gpu] GPU time (sum over streams) %.2fs: H2D %.2fs, filter %.2fs, full check %.2fs; host pipeline wall (sum over threads) %.2fs; %.3g candidates per GPU-second (filter), %.3g per pipeline-second\n",
            tg, s.t_h2d, s.t_filter, s.t_full, s.t_host, s.t_filter > 0 ? (double)s.cand / s.t_filter : 0.0, s.t_host > 0 ? (double)s.cand / s.t_host : 0.0);
}

extern "C" void gpu_shutdown(void)
{
    for (int i = 0; i < NTB; i++) tb_free(&TB[i]);
    free(TB); TB = NULL; NTB = 0;
}
