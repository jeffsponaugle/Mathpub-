/*
 * a007508_cuda.cu - GPU twin prime counter for OEIS A007508 (DGX Spark / GB10)
 *
 * Same problem, k-space layout, chunking and log format as the CPU tool
 * (../cpu/a007508.c): a twin pair (p, p+2), p > 5, has p = 30k + r with
 * r in {11, 17, 29}; candidates are indexed by k and class j, bit 3k + j.
 * Here a set bit means "composite" (p or p+2 has a factor), so marking is an
 * atomicOr and survivors are the zero bits.
 *
 * A chunk [K0, K1) is processed in super-chunks of NSUB sub-chunks of NSEG
 * segments of SEG_K k-units.  Sieving primes are tiered by size:
 *   q <= QSPLIT       sieved inside kernel_segment (shared-memory atomics).
 *   QSPLIT < q <= QL2 kernel_large (per sub-chunk): one thread per prime marks
 *                     its hits in a bitmap that stays in the 24 MB L2; the
 *                     prime keeps its next hit as state.  ~40 ps per hit.
 *   q > QL2           record lists, ~16 bytes of memory traffic per hit:
 *   kernel_generate   (per super-chunk) one thread per sieving prime
 *                     walks the prime's hits in the super-chunk
 *                     and appends 4-byte records to a list per sub-chunk,
 *                     through shared-memory bins flushed in bulk (one global
 *                     atomic per bin flush).  Stores the next hit as state.
 *   kernel_distribute (per sub-chunk) splits the sub-chunk's list by segment.
 *   kernel_segment    (per sub-chunk, one block per segment) shared-memory
 *                     bitmap = presieve patterns 7..97 | L2 bitmap | list,
 *                     then primes 101..QSPLIT with shared atomics (warp per
 *                     prime below PWARP, thread per prime above); then the
 *                     zero bits are counted.
 * The tier bounds are compile-time (-DQSPLIT=, -DQL2=).
 *
 * Build:  nvcc -O3 -std=c++17 -arch=sm_121 -o a007508_cuda a007508_cuda.cu \
 *              -I$HOME/local/include -L$HOME/local/lib -lprimesieve
 * Usage:  same options as the CPU tool (no -t/-k/-S); see usage().
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cinttypes>
#include <ctime>
#include <vector>
#include <algorithm>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cuda_runtime.h>
#include <primesieve.h>

typedef uint64_t u64; typedef uint32_t u32; typedef uint8_t u8; typedef int64_t i64; typedef int32_t i32;
typedef unsigned __int128 u128;

#ifndef SEG_K
#define SEG_K     (1u << 18)              /* k-units per segment: 786,432 bits = 96 KB shared */
#endif
#ifndef NSEG
#define NSEG      64                      /* segments per sub-chunk (power of two; SUB_K = 2^24 k-units) */
#endif
#ifndef NSUB
#define NSUB      32                      /* sub-chunks per super-chunk (large primes visited once per super-chunk) */
#endif
#ifndef BLOCKDIM
#define BLOCKDIM  1024
#endif
#ifndef QSPLIT
#define QSPLIT    (1u << 21)              /* primes <= QSPLIT: sieved per segment in shared memory */
#endif
#ifndef QL2
#define QL2       (1u << 26)              /* QSPLIT < q <= QL2: per sub-chunk into an L2-resident bitmap; above: record lists */
#endif
#ifndef PWARP
#define PWARP     4096                    /* primes < PWARP: one warp per prime */
#endif
#define SEG_BITS  (3u * SEG_K)
#define SEG_WORDS (SEG_BITS / 32)
#define SEG_BYTES (SEG_BITS / 8)
#define SUB_K     ((u64)SEG_K * NSEG)
#define SUB_WORDS ((u64)SEG_WORDS * NSEG)
static_assert((NSEG & (NSEG - 1)) == 0 && (SEG_K & (SEG_K - 1)) == 0, "NSEG and SEG_K must be powers of two");
#ifndef APPLY_L2
#define APPLY_L2  1                       /* 1: large-prime records -> L2 bitmap atomics; 0: per-segment lists */
#endif
#define SUB_LOG   (__builtin_ctz(SEG_K) + __builtin_ctz(NSEG))
#define SUPER_K   (SUB_K * NSUB)
#define GEN_THREADS 256
#define GEN_B     8                       /* hits per thread per generator phase */
#define BIN_CAP   256                     /* entries per shared-memory bin */

#define CUDA_CHECK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { fprintf(stderr, "CUDA error %s at line %d: %s\n", #x, __LINE__, cudaGetErrorString(e_)); exit(1); } } while (0)

/* ======================================================================= */
/*  Wheel tables (host), copied to constant memory                         */
/* ======================================================================= */

static const int R3[3] = {11, 17, 29};
static const int W8[8] = {1, 7, 11, 13, 17, 19, 23, 29};
static int8_t WIDX[30];
static i32 WM[8][6], WC[8][6], DKM[8][6], DKC[8][6];
static u8  WJ[8][6];
static i32 MAXDM;

__constant__ i32 c_WM[8][6], c_WC[8][6], c_DKM[8][6], c_DKC[8][6];
__constant__ u8  c_WJ[8][6];

typedef struct { int m, c, j, t; } wentry_t;

static int modinv30(int w) { for (int x = 1; x < 30; x++) if ((w * x) % 30 == 1) return x; return -1; }

static void wheel_init(void)
{
    memset(WIDX, -1, sizeof WIDX);
    for (int i = 0; i < 8; i++) WIDX[W8[i]] = (int8_t)i;
    MAXDM = 0;
    for (int wi = 0; wi < 8; wi++) {
        int w = W8[wi], winv = modinv30(w);
        wentry_t e[6]; int cnt = 0;
        for (int j = 0; j < 3; j++)
            for (int ti = 0; ti < 2; ti++) {
                int t = ti ? -2 : 0, d = R3[j] - t, m = (d * winv) % 30, num = m * w - d;
                if (num % 30 != 0) { fprintf(stderr, "wheel bug\n"); exit(1); }
                e[cnt].m = m; e[cnt].c = num / 30; e[cnt].j = j; e[cnt].t = t; cnt++;
            }
        for (int a = 1; a < 6; a++)
            for (int b = a; b > 0 && e[b - 1].m > e[b].m; b--) { wentry_t tmp = e[b]; e[b] = e[b - 1]; e[b - 1] = tmp; }
        for (int u = 0; u < 6; u++) { WM[wi][u] = e[u].m; WC[wi][u] = e[u].c; WJ[wi][u] = (u8)e[u].j; }
        for (int u = 0; u < 6; u++) {
            int v = (u + 1) % 6;
            DKM[wi][u] = e[v].m - e[u].m + (v == 0 ? 30 : 0);
            DKC[wi][u] = e[v].c - e[u].c + (v == 0 ? w : 0);
            if (DKM[wi][u] > MAXDM) MAXDM = DKM[wi][u];
        }
    }
    CUDA_CHECK(cudaMemcpyToSymbol(c_WM, WM, sizeof WM));
    CUDA_CHECK(cudaMemcpyToSymbol(c_WC, WC, sizeof WC));
    CUDA_CHECK(cudaMemcpyToSymbol(c_DKM, DKM, sizeof DKM));
    CUDA_CHECK(cudaMemcpyToSymbol(c_DKC, DKC, sizeof DKC));
    CUDA_CHECK(cudaMemcpyToSymbol(c_WJ, WJ, sizeof WJ));
}

/* ======================================================================= */
/*  Presieve patterns 7..97, composite polarity (1 = killed)               */
/* ======================================================================= */

#define NPAT 7
static const int PATP[NPAT][5] = {{7,11,13,17,19},{23,29,31,0,0},{37,41,43,0,0},{47,53,59,0,0},{61,67,71,0,0},{73,79,83,0,0},{89,97,0,0,0}};
static u64 PATLEN[NPAT];          /* byte period = 3 * product */
static u8 *PATH[NPAT];            /* host pattern, absolute k, bit 3k+j */
static u8 *PATD[NPAT];            /* device pattern, re-aligned to the chunk start */
static u64 PATLEN_MAX;

static void presieve_init(void)
{
    PATLEN_MAX = 0;
    for (int g = 0; g < NPAT; g++) {
        u64 per = 1;
        for (int i = 0; i < 5; i++) if (PATP[g][i]) per *= (u64)PATP[g][i];
        PATLEN[g] = 3 * per;
        if (PATLEN[g] > PATLEN_MAX) PATLEN_MAX = PATLEN[g];
        u8 *b = (u8 *)calloc(PATLEN[g], 1);
        for (int i = 0; i < 5; i++) {
            int q = PATP[g][i]; if (!q) continue;
            for (int j = 0; j < 3; j++)
                for (int k0 = 0; k0 < q; k0++) {
                    int n = (30 * k0 + R3[j]) % q;
                    if (n != 0 && n != q - 2) continue;
                    for (u64 k = (u64)k0; k < 8 * per; k += (u64)q) { u64 bb = 3 * k + (u64)j; b[bb >> 3] |= (u8)(1u << (bb & 7)); }
                }
        }
        PATH[g] = b;
        CUDA_CHECK(cudaMalloc(&PATD[g], PATLEN[g] + 8));
    }
}

/* re-align pattern g so that chunk-relative bit b is absolute bit 3*K0 + b; upload */
static void presieve_align_upload(u64 K0)
{
    static u8 *tmp = NULL;
    if (!tmp) tmp = (u8 *)malloc(PATLEN_MAX + 8);
    u64 B = 3 * K0; unsigned ph = (unsigned)(B & 7);
    for (int g = 0; g < NPAT; g++) {
        u64 len = PATLEN[g], off = (B >> 3) % len;
        const u8 *src = PATH[g];
        for (u64 i = 0; i < len + 8; i++) {           /* 8 bytes of wrap-around slack */
            u64 i0 = (off + i) % len, i1 = (i0 + 1) % len;
            tmp[i] = ph ? (u8)((src[i0] >> ph) | (src[i1] << (8 - ph))) : src[i0];
        }
        CUDA_CHECK(cudaMemcpy(PATD[g], tmp, len + 8, cudaMemcpyHostToDevice));
    }
}

/* ======================================================================= */
/*  Device kernels                                                         */
/* ======================================================================= */

/* First hit >= K0 of the wheel for prime (qd, wi), fast path (K0 >= q). Returns k (absolute) and u. */
__device__ __forceinline__ void first_hit(u32 qd, u32 wi, u64 K0, u64 *kout, u32 *uout)
{
    u64 q = 30ull * qd + (u64)(wi == 0 ? 1 : wi == 1 ? 7 : wi == 2 ? 11 : wi == 3 ? 13 : wi == 4 ? 17 : wi == 5 ? 19 : wi == 6 ? 23 : 29);
    u64 base = K0 % q;
    u64 best = ~0ull; u32 bu = 0;
    #pragma unroll
    for (int u = 0; u < 6; u++) {
        i64 v = (i64)c_WM[wi][u] * (i64)qd + c_WC[wi][u];
        if (v >= (i64)q) v %= (i64)q;
        u64 rho = (u64)v;
        u64 off = rho >= base ? rho - base : rho + q - base;
        if (off < best) { best = off; bu = (u32)u; }
    }
    *kout = K0 + best; *uout = bu;
}

/* state = (k_rel << 3) | u, k_rel = next hit relative to chunk start */
__global__ void kernel_init_state(const u32 *__restrict__ qd, const u8 *__restrict__ wi, u64 *__restrict__ state, u32 n, u64 K0)
{
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    u64 k; u32 u;
    first_hit(qd[i], wi[i], K0, &k, &u);
    state[i] = ((k - K0) << 3) | u;
}

/* Medium primes: mark hits in [sub_lo, sub_hi) (chunk-relative k) into the sub-chunk bitmap. */
__global__ void __launch_bounds__(256)
kernel_large(const u32 *__restrict__ qd, const u8 *__restrict__ wi, u64 *__restrict__ state, u32 n,
             u64 sub_lo, u64 sub_hi, u32 *__restrict__ bitmap)
{
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    u64 s = __ldcs(&state[i]);
    u64 k = s >> 3;
    if (k >= sub_hi) return;
    u32 u = (u32)(s & 7), w = wi[i];
    const i64 q_d = qd[i];
    do {
        u64 b = 3 * (k - sub_lo) + c_WJ[w][u];
        atomicOr(&bitmap[b >> 5], 1u << (b & 31));
        k = (u64)((i64)k + (i64)c_DKM[w][u] * q_d + c_DKC[w][u]);
        u = (u == 5) ? 0 : u + 1;
    } while (k < sub_hi);
    __stcs(&state[i], (k << 3) | u);
}

/* sub_lo mod q for the small primes, once per sub-chunk */
__global__ void kernel_submod(const u32 *__restrict__ sq, u32 *__restrict__ rsub, u32 nsmall, u64 sub_lo_abs)
{
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nsmall) rsub[i] = (u32)(sub_lo_abs % sq[i]);
}


/* Large primes: walk the hits of prime i in [super_lo, super_hi) (chunk-relative
 * k) and append records (koff_in_subchunk << 2 | j) to glist[t] for sub-chunk t.
 * Phases of at most GEN_B hits per thread, block-wide flush of the shared bins. */
__global__ void __launch_bounds__(GEN_THREADS)
kernel_generate(const u32 *__restrict__ qd, const u8 *__restrict__ wi, u64 *__restrict__ state, u32 n,
                u64 super_lo, u64 super_hi, u32 *__restrict__ glist, u32 *__restrict__ gcount, u32 gcap,
                u32 *__restrict__ overflow)
{
    extern __shared__ u32 binmem[];                 /* NSUB * BIN_CAP */
    u32 (*bins)[BIN_CAP] = (u32 (*)[BIN_CAP])binmem;
    __shared__ u32 bcnt[NSUB];
    __shared__ int s_active;
    const u32 tid = threadIdx.x, lane = tid & 31, warp = tid >> 5;
    const u32 i = blockIdx.x * GEN_THREADS + tid;
    bool has = false, touched = false;
    u64 k = 0; u32 u = 0, w = 0; i64 q_d = 0;
    if (i < n) {
        u64 s = state[i];
        k = s >> 3; u = (u32)(s & 7);
        if (k < super_hi) { has = true; touched = true; w = wi[i]; q_d = qd[i]; }
    }
    if (tid < NSUB) bcnt[tid] = 0;
    __syncthreads();
    for (;;) {
        if (tid == 0) s_active = 0;
        __syncthreads();
        if (has) {
            for (int b = 0; b < GEN_B; b++) {
                u64 rel = k - super_lo;
                u32 t = (u32)(rel >> SUB_LOG);
                u32 rec = ((u32)(rel & (SUB_K - 1)) << 2) | c_WJ[w][u];
                u32 pos = atomicAdd(&bcnt[t], 1u);
                if (pos >= BIN_CAP) { atomicSub(&bcnt[t], 1u); break; }      /* bin full: retry next phase */
                bins[t][pos] = rec;
                k = (u64)((i64)k + (i64)c_DKM[w][u] * q_d + c_DKC[w][u]);
                u = (u == 5) ? 0 : u + 1;
                if (k >= super_hi) { has = false; break; }
            }
            if (has) s_active = 1;
        }
        __syncthreads();
        for (u32 t = warp; t < NSUB; t += GEN_THREADS / 32) {
            u32 c = bcnt[t];
            if (c > BIN_CAP) c = BIN_CAP;
            if (c) {
                u32 base = 0;
                if (lane == 0) base = atomicAdd(&gcount[t], c);
                base = __shfl_sync(0xffffffffu, base, 0);
                if (base + c > gcap) { if (lane == 0) atomicMax(overflow, base + c); }
                else for (u32 e = lane; e < c; e += 32) glist[(u64)t * gcap + base + e] = bins[t][e];
            }
        }
        __syncthreads();
        if (tid < NSUB) bcnt[tid] = 0;
        int active = s_active;
        __syncthreads();
        if (!active) break;
    }
    if (touched) state[i] = (k << 3) | u;
}

/* Split the records of one sub-chunk by segment; output = bit offset in segment. */
__global__ void __launch_bounds__(256)
kernel_distribute(const u32 *__restrict__ glist_t, const u32 *__restrict__ pn, u32 gcap, u32 *__restrict__ seglist,
                  u32 *__restrict__ segcount, u32 segcap, u32 *__restrict__ overflow)
{
    __shared__ u32 hcnt[NSEG], hbase[NSEG];
    const u32 tid = threadIdx.x;
    u32 n = *pn; if (n > gcap) n = gcap;
    const u32 PER = 8, TILE = 256 * PER;
    for (u32 tile = blockIdx.x * TILE; tile < n; tile += gridDim.x * TILE) {
        if (tid < NSEG) hcnt[tid] = 0;
        __syncthreads();
        u32 recs[PER], rank[PER], segs[PER];
        #pragma unroll
        for (u32 p = 0; p < PER; p++) {
            u32 idx = tile + p * 256 + tid;
            segs[p] = NSEG;
            if (idx < n) {
                u32 r = glist_t[idx];
                u32 seg = r >> (2 + __builtin_ctz(SEG_K));
                segs[p] = seg;
                rank[p] = atomicAdd(&hcnt[seg], 1u);
                recs[p] = 3 * ((r >> 2) & (SEG_K - 1)) + (r & 3);
            }
        }
        __syncthreads();
        if (tid < NSEG) hbase[tid] = hcnt[tid] ? atomicAdd(&segcount[tid], hcnt[tid]) : 0;
        __syncthreads();
        #pragma unroll
        for (u32 p = 0; p < PER; p++)
            if (segs[p] < NSEG) {
                u32 o = hbase[segs[p]] + rank[p];
                if (o < segcap) seglist[(u64)segs[p] * segcap + o] = recs[p];
                else atomicMax(overflow, o + 1);
            }
        __syncthreads();
    }
}

/* Apply one sub-chunk's records to the L2-resident bitmap. */
__global__ void __launch_bounds__(256)
kernel_apply(const u32 *__restrict__ glist_t, const u32 *__restrict__ pn, u32 gcap, u32 *__restrict__ bitmap)
{
    u32 n = *pn; if (n > gcap) n = gcap;
    for (u32 i = blockIdx.x * 256 + threadIdx.x; i < n; i += gridDim.x * 256) {
        u32 r = __ldcs(&glist_t[i]);
        u32 b = 3 * (r >> 2) + (r & 3);
        atomicOr(&bitmap[b >> 5], 1u << (b & 31));
    }
}

struct PatPtrs { const u8 *p[NPAT]; u32 len[NPAT]; };

/* One block per segment: shared bitmap = large-prime bits | patterns, then the
 * small primes (101..QSPLIT), then count the zero bits of the valid range. */
__global__ void __launch_bounds__(BLOCKDIM)
kernel_segment(const u32 *__restrict__ bitmap, const u32 *__restrict__ seglist, const u32 *__restrict__ segcount, u32 segcap, PatPtrs pats,
               const u32 *__restrict__ sq, const u32 *__restrict__ sqd, const u8 *__restrict__ swi,
               const u32 *__restrict__ srho, const u32 *__restrict__ rsub, const u32 *__restrict__ segmod,
               const u32 *__restrict__ sinv, u32 nsmall, u32 nwarp,
               u64 sub_lo, u64 chunk_len, unsigned long long *__restrict__ count)
{
    extern __shared__ u32 sh[];
    __shared__ u32 s_scratch[32];
    const u32 seg = blockIdx.x, tid = threadIdx.x;
    const u64 seg_k0 = sub_lo + (u64)seg * SEG_K;          /* chunk-relative k */
    if (seg_k0 >= chunk_len) return;

    /* --- 1. large-prime bits | presieve patterns --- */
    {
        const u64 byte0 = (seg_k0 * 3) >> 3;                /* chunk-relative byte of word 0 */
        u32 p0[NPAT];
        #pragma unroll
        for (int g = 0; g < NPAT; g++) p0[g] = (u32)(byte0 % pats.len[g]) + 4 * tid;
        #pragma unroll
        for (int g = 0; g < NPAT; g++) while (p0[g] >= pats.len[g]) p0[g] -= pats.len[g];
        const u32 *bm = bitmap + (u64)seg * SEG_WORDS;
        for (u32 w = tid; w < SEG_WORDS; w += BLOCKDIM) {
            u32 v = __ldcs(&bm[w]);
            #pragma unroll
            for (int g = 0; g < NPAT; g++) {
                const u8 *p = pats.p[g] + p0[g];
                v |= (u32)p[0] | ((u32)p[1] << 8) | ((u32)p[2] << 16) | ((u32)p[3] << 24);
                p0[g] += 4 * BLOCKDIM;
                if (p0[g] >= pats.len[g]) p0[g] -= pats.len[g];
            }
            sh[w] = v;
        }
    }
    __syncthreads();
    if (segcount) {
        const u32 *lst = seglist + (u64)seg * segcap;
        const u32 cnt = segcount[seg] < segcap ? segcount[seg] : segcap;
        for (u32 i = tid; i < cnt; i += BLOCKDIM) { u32 b = __ldcs(&lst[i]); atomicOr(&sh[b >> 5], 1u << (b & 31)); }
        __syncthreads();
    }

    /* --- 2. small primes --- */
    {
        const u32 warp = tid >> 5, lane = tid & 31, nwarps = BLOCKDIM / 32;
        /* warp per prime: lane -> (progression u = lane % 6, sub-stride lane / 6), 30 lanes busy */
        for (u32 i = warp; i < nwarp; i += nwarps) {
            const u32 q = sq[i];
            u32 x = rsub[i] + seg * segmod[i];
            u32 r = x - __umulhi(x, sinv[i]) * q;
            if (r >= q) r -= q;                                /* seg_k0 mod q */
            if (lane < 30) {
                const u32 u = lane % 6, sub = lane / 6, w = swi[i];
                u32 rho = srho[6 * i + u];
                u32 off = rho >= r ? rho - r : rho + q - r;   /* first hit, k-offset in segment */
                const u32 stride = 15 * q;                     /* 5 sub-strides of 3q bits */
                for (u32 b = 3 * off + c_WJ[w][u] + sub * 3 * q; b < SEG_BITS; b += stride)
                    atomicOr(&sh[b >> 5], 1u << (b & 31));
            }
        }
        /* thread per prime */
        for (u32 i = nwarp + tid; i < nsmall; i += BLOCKDIM) {
            const u32 q = sq[i];
            u32 x = rsub[i] + seg * segmod[i];
            u32 r = x - __umulhi(x, sinv[i]) * q;
            if (r >= q) r -= q;
            const u32 w = swi[i], stride = 3 * q;
            #pragma unroll
            for (u32 u = 0; u < 6; u++) {
                u32 rho = srho[6 * i + u];
                u32 off = rho >= r ? rho - r : rho + q - r;
                for (u32 b = 3 * off + c_WJ[w][u]; b < SEG_BITS; b += stride)
                    atomicOr(&sh[b >> 5], 1u << (b & 31));
            }
        }
    }
    __syncthreads();

    /* --- 3. count survivors (zero bits) below the chunk end --- */
    {
        u64 valid = chunk_len - seg_k0;
        u32 vbits = valid >= (u64)SEG_K ? SEG_BITS : (u32)(3 * valid);
        u32 c = 0;
        for (u32 w = tid; w < SEG_WORDS; w += BLOCKDIM) {
            u32 x = ~sh[w];
            u32 lo = 32 * w;
            if (lo + 32 <= vbits) c += __popc(x);
            else if (lo < vbits) c += __popc(x & ((1u << (vbits - lo)) - 1));
        }
        for (int o = 16; o > 0; o >>= 1) c += __shfl_down_sync(0xffffffffu, c, o);
        const u32 warp = tid >> 5, lane = tid & 31;
        if (lane == 0) s_scratch[warp] = c;
        __syncthreads();
        if (warp == 0) {
            u32 s = (lane < BLOCKDIM / 32) ? s_scratch[lane] : 0;
            for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(0xffffffffu, s, o);
            if (lane == 0) atomicAdd(count, (unsigned long long)s);
        }
    }
}

/* ======================================================================= */
/*  Host driver for one chunk                                              */
/* ======================================================================= */

static u64 isqrt128(u128 x)
{
    u64 r = (u64)sqrtl((long double)x);
    while ((u128)r * r > x) r--;
    while ((u128)(r + 1) * (r + 1) <= x) r++;
    return r;
}

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static int g_verbose = 0;

/* Count twin pairs with K0 <= k < K1 on the GPU.  Requires K0 >= sqrt(30*K1). */
static u64 gpu_count(u64 K0, u64 K1)
{
    if (K1 <= K0) return 0;
    const u64 len = K1 - K0;
    const u64 qmax = isqrt128((u128)30 * K1 + 1);
    if (K0 < qmax) { fprintf(stderr, "chunk start %" PRIu64 " is below sqrt(end): use the CPU tool for this range\n", K0); exit(1); }
    double t0 = now();

    /* sieving primes 101..qmax */
    size_t np = 0;
    u64 *primes = (u64 *)primesieve_generate_primes(101, qmax, &np, UINT64_PRIMES);
    size_t nsmall = 0;
    while (nsmall < np && primes[nsmall] <= QSPLIT) nsmall++;
    size_t nwarp = 0;
    while (nwarp < nsmall && primes[nwarp] < PWARP) nwarp++;
    size_t nmed_end = nsmall;
    while (nmed_end < np && primes[nmed_end] <= QL2) nmed_end++;
    size_t nmed = nmed_end - nsmall, nlarge = np - nmed_end;

    std::vector<u32> sq(nsmall), sqd(nsmall), srho(6 * nsmall), segmod(nsmall), sinv(nsmall);
    std::vector<u8> swi(nsmall);
    for (size_t i = 0; i < nsmall; i++) {
        u64 q = primes[i];
        sq[i] = (u32)q; sqd[i] = (u32)(q / 30); swi[i] = (u8)WIDX[q % 30];
        for (int u = 0; u < 6; u++) {
            i64 v = (i64)WM[swi[i]][u] * (i64)sqd[i] + WC[swi[i]][u];
            v %= (i64)q; if (v < 0) v += (i64)q;
            srho[6 * i + u] = (u32)v;
        }
        segmod[i] = (u32)((u64)SEG_K % q);
        sinv[i] = (u32)((1ull << 32) / q);
    }
    std::vector<u32> lqd(nlarge + nmed); std::vector<u8> lwi(nlarge + nmed);   /* medium first, then large */
    for (size_t i = 0; i < nmed + nlarge; i++) { u64 q = primes[nsmall + i]; lqd[i] = (u32)(q / 30); lwi[i] = (u8)WIDX[q % 30]; }
    primesieve_free(primes);

    u32 *d_sq, *d_sqd, *d_srho, *d_segmod, *d_sinv, *d_lqd; u8 *d_swi, *d_lwi; u64 *d_state;
    unsigned long long *d_count;
    CUDA_CHECK(cudaMalloc(&d_sq, nsmall * 4 + 4));      CUDA_CHECK(cudaMalloc(&d_sqd, nsmall * 4 + 4));
    CUDA_CHECK(cudaMalloc(&d_srho, nsmall * 24 + 4));   CUDA_CHECK(cudaMalloc(&d_segmod, nsmall * 4 + 4));
    CUDA_CHECK(cudaMalloc(&d_sinv, nsmall * 4 + 4));
    CUDA_CHECK(cudaMalloc(&d_swi, nsmall + 4));
    CUDA_CHECK(cudaMalloc(&d_lqd, (nmed + nlarge) * 4 + 4));     CUDA_CHECK(cudaMalloc(&d_lwi, nmed + nlarge + 4));
    CUDA_CHECK(cudaMalloc(&d_state, (nmed + nlarge) * 8 + 8));
    CUDA_CHECK(cudaMalloc(&d_count, 8));
    u32 *d_bitmap[2], *d_rsub2[2];
    for (int p = 0; p < 2; p++) { CUDA_CHECK(cudaMalloc(&d_bitmap[p], SUB_WORDS * 4)); CUDA_CHECK(cudaMalloc(&d_rsub2[p], nsmall * 4 + 4)); }
    /* record lists: expected hits per sub-chunk from primes in (QSPLIT, qmax] */
    double hits_per_k = nlarge ? 6.0 * (log(log((double)qmax)) - log(log((double)QL2))) : 0.0;
    u32 gcap = (u32)(hits_per_k * (double)SUB_K * 1.25) + (1u << 20);
    u32 segcap = (u32)(hits_per_k * (double)SEG_K * 1.25) + (1u << 16);
    u32 *d_glist, *d_gcount, *d_seglist[2], *d_segcount[2], *d_overflow;
    CUDA_CHECK(cudaMalloc(&d_glist, (u64)NSUB * gcap * 4));
    CUDA_CHECK(cudaMalloc(&d_gcount, NSUB * 4));
    for (int p = 0; p < 2; p++) { CUDA_CHECK(cudaMalloc(&d_seglist[p], APPLY_L2 ? 4 : (u64)NSEG * segcap * 4)); CUDA_CHECK(cudaMalloc(&d_segcount[p], NSEG * 4)); }
    CUDA_CHECK(cudaMalloc(&d_overflow, 4));
    CUDA_CHECK(cudaMemset(d_overflow, 0, 4));
    CUDA_CHECK(cudaMemcpy(d_sq, sq.data(), nsmall * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_sqd, sqd.data(), nsmall * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_srho, srho.data(), nsmall * 24, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_segmod, segmod.data(), nsmall * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_sinv, sinv.data(), nsmall * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_swi, swi.data(), nsmall, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_lqd, lqd.data(), (nmed + nlarge) * 4, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_lwi, lwi.data(), nmed + nlarge, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(d_count, 0, 8));
    presieve_align_upload(K0);
    PatPtrs pp;
    for (int g = 0; g < NPAT; g++) { pp.p[g] = PATD[g]; pp.len[g] = (u32)PATLEN[g]; }

    if (nmed + nlarge) kernel_init_state<<<(unsigned)((nmed + nlarge + 255) / 256), 256>>>(d_lqd, d_lwi, d_state, (u32)(nmed + nlarge), K0);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaFuncSetAttribute(kernel_segment, cudaFuncAttributeMaxDynamicSharedMemorySize, SEG_BYTES));
    CUDA_CHECK(cudaFuncSetAttribute(kernel_generate, cudaFuncAttributeMaxDynamicSharedMemorySize, NSUB * BIN_CAP * 4));
    double t1 = now();

    u64 nsub = (len + SUB_K - 1) / SUB_K;
    static cudaStream_t sA = NULL, sB = NULL;
    static cudaEvent_t evMark[2], evSeg[2];
    if (!sA) {
        CUDA_CHECK(cudaStreamCreate(&sA)); CUDA_CHECK(cudaStreamCreate(&sB));
        for (int p = 0; p < 2; p++) { CUDA_CHECK(cudaEventCreateWithFlags(&evMark[p], cudaEventDisableTiming)); CUDA_CHECK(cudaEventCreateWithFlags(&evSeg[p], cudaEventDisableTiming)); }
    }
    std::vector<cudaEvent_t> ev;
    if (g_verbose) { ev.resize(4 * nsub); for (auto &e : ev) CUDA_CHECK(cudaEventCreate(&e)); }
    for (u64 sub = 0; sub < nsub; sub++) {
        const int p = (int)(sub & 1);
        u64 sub_lo = sub * SUB_K, sub_hi = sub_lo + SUB_K;
        if (sub_hi > len) sub_hi = len;
        if (sub >= 2) CUDA_CHECK(cudaStreamWaitEvent(sA, evSeg[p], 0));        /* buffer p free again */
        if (sub % NSUB == 0) {                                                  /* new super-chunk: generate */
            u64 super_hi = sub_lo + SUPER_K; if (super_hi > len) super_hi = len;
            CUDA_CHECK(cudaMemsetAsync(d_gcount, 0, NSUB * 4, sA));
            if (g_verbose) cudaEventRecord(ev[4 * sub], sA);
            if (nlarge) kernel_generate<<<(unsigned)((nlarge + GEN_THREADS - 1) / GEN_THREADS), GEN_THREADS, NSUB * BIN_CAP * 4, sA>>>(
                d_lqd + nmed, d_lwi + nmed, d_state + nmed, (u32)nlarge, sub_lo, super_hi, d_glist, d_gcount, gcap, d_overflow);
            if (g_verbose) cudaEventRecord(ev[4 * sub + 1], sA);
        } else if (g_verbose) { cudaEventRecord(ev[4 * sub], sA); cudaEventRecord(ev[4 * sub + 1], sA); }
        u32 t = (u32)(sub % NSUB);
        CUDA_CHECK(cudaMemsetAsync(d_bitmap[p], 0, SUB_WORDS * 4, sA));
        if (nmed) kernel_large<<<(unsigned)((nmed + 255) / 256), 256, 0, sA>>>(d_lqd, d_lwi, d_state, (u32)nmed, sub_lo, sub_hi, d_bitmap[p]);
#if APPLY_L2
        if (nlarge) kernel_apply<<<1024, 256, 0, sA>>>(d_glist + (u64)t * gcap, d_gcount + t, gcap, d_bitmap[p]);
#else
        CUDA_CHECK(cudaMemsetAsync(d_segcount[p], 0, NSEG * 4, sA));
        if (nlarge) kernel_distribute<<<512, 256, 0, sA>>>(d_glist + (u64)t * gcap, d_gcount + t, gcap, d_seglist[p], d_segcount[p], segcap, d_overflow);
#endif
        if (nsmall) kernel_submod<<<(unsigned)((nsmall + 255) / 256), 256, 0, sA>>>(d_sq, d_rsub2[p], (u32)nsmall, K0 + sub_lo);
        if (g_verbose) cudaEventRecord(ev[4 * sub + 2], sA);
        CUDA_CHECK(cudaEventRecord(evMark[p], sA));
        CUDA_CHECK(cudaStreamWaitEvent(sB, evMark[p], 0));
        unsigned nseg = (unsigned)((sub_hi - sub_lo + SEG_K - 1) / SEG_K);
        kernel_segment<<<nseg, BLOCKDIM, SEG_BYTES, sB>>>(d_bitmap[p], d_seglist[p], APPLY_L2 ? NULL : d_segcount[p], segcap, pp, d_sq, d_sqd, d_swi, d_srho, d_rsub2[p], d_segmod, d_sinv,
                                                         (u32)nsmall, (u32)nwarp, sub_lo, len, d_count);
        if (g_verbose) cudaEventRecord(ev[4 * sub + 3], sB);
        CUDA_CHECK(cudaEventRecord(evSeg[p], sB));
        CUDA_CHECK(cudaGetLastError());
    }
    CUDA_CHECK(cudaStreamSynchronize(sA));
    CUDA_CHECK(cudaStreamSynchronize(sB));
    unsigned long long cnt = 0;
    CUDA_CHECK(cudaMemcpy(&cnt, d_count, 8, cudaMemcpyDeviceToHost));
    double t2 = now();
    if (g_verbose) {
        float tg = 0, td = 0, ts = 0, ms;
        for (u64 sub = 0; sub < nsub; sub++) {
            cudaEventElapsedTime(&ms, ev[4 * sub], ev[4 * sub + 1]); tg += ms;
            cudaEventElapsedTime(&ms, ev[4 * sub + 1], ev[4 * sub + 2]); td += ms;
            cudaEventElapsedTime(&ms, ev[4 * sub + 2], ev[4 * sub + 3]); ts += ms;
        }
        for (auto &e : ev) cudaEventDestroy(e);
        fprintf(stderr, "# chunk [%" PRIu64 ", %" PRIu64 "): %.1e numbers, %zu small + %zu medium + %zu large primes, setup %.2f s, sieve %.2f s (stream A: generate %.2f, L2+distribute %.2f; stream B: segment %.2f), lists %u/%u, %.3e numbers/s\n",
                K0, K1, (double)len * 30, nsmall, nmed, nlarge, t1 - t0, t2 - t1, tg / 1000, td / 1000, ts / 1000, gcap, segcap, (double)len * 30 / (t2 - t1));
    }

    cudaFree(d_sq); cudaFree(d_sqd); cudaFree(d_srho); cudaFree(d_segmod); cudaFree(d_sinv);
    cudaFree(d_swi); cudaFree(d_lqd); cudaFree(d_lwi); cudaFree(d_state); cudaFree(d_count);
    cudaFree(d_glist); cudaFree(d_gcount); cudaFree(d_overflow);
    for (int p = 0; p < 2; p++) { cudaFree(d_seglist[p]); cudaFree(d_segcount[p]); cudaFree(d_bitmap[p]); cudaFree(d_rsub2[p]); }
    return (u64)cnt;
}

/* ======================================================================= */
/*  Decades, logging, main (same conventions as the CPU tool)              */
/* ======================================================================= */

static u128 pow10_128(int n) { u128 r = 1; while (n-- > 0) r *= 10; return r; }
static u64 Kbound(int n) { return n < 1 ? 0 : (u64)((pow10_128(n) - 10) / 30); }

static u64 default_chunk_k(int n)
{
    int e = n - 4; if (e < 10) e = 10; if (e > 14) e = 14;
    return (u64)(pow10_128(e) / 30);
}

/* twin pairs (p, p+2) with p + 2 < 10^n for n <= 4, by trial division */
static u64 small_decade_count(int n)
{
    u64 lim = (u64)pow10_128(n), c = 0, lo = n > 1 ? (u64)pow10_128(n - 1) : 0;
    auto isp = [](u64 x) { if (x < 2) return false; for (u64 d = 2; d * d <= x; d++) if (x % d == 0) return false; return true; };
    for (u64 p = lo; p + 2 < lim; p++) if (isp(p) && isp(p + 2)) c++;
    return c;
}

static inline int  bit_get(const u8 *b, u64 i) { return (b[i >> 3] >> (i & 7)) & 1; }
static inline void bit_set(u8 *b, u64 i)       { b[i >> 3] |= (u8)(1u << (i & 7)); }

static void trim_partial_line(const char *path)
{
    int fd = open(path, O_RDWR);
    if (fd < 0) return;
    struct stat st;
    if (fstat(fd, &st) != 0 || st.st_size == 0) { close(fd); return; }
    char c;
    if (pread(fd, &c, 1, st.st_size - 1) == 1 && c == '\n') { close(fd); return; }
    off_t keep = 0, pos = st.st_size; char buf[4096];
    while (pos > 0 && keep == 0) {
        size_t n = (size_t)(pos < (off_t)sizeof buf ? pos : (off_t)sizeof buf);
        pos -= (off_t)n;
        if (pread(fd, buf, n, pos) != (ssize_t)n) break;
        for (size_t i = n; i-- > 0; ) if (buf[i] == '\n') { keep = pos + (off_t)i + 1; break; }
    }
    fprintf(stderr, "# log %s: dropping %lld bytes of partial last line\n", path, (long long)(st.st_size - keep));
    if (ftruncate(fd, keep) != 0) perror("ftruncate");
    close(fd);
}

static u64 load_chunks(const char *path, int n, u64 chunk_k, u64 nchunks, u8 *done, u64 *total)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    u64 resumed = 0; char line[256];
    while (fgets(line, sizeof line, f)) {
        int ln; u64 idx, ck, cnt;
        if (!strchr(line, '\n')) continue;
        if (sscanf(line, "k %d %" SCNu64 " %" SCNu64 " %" SCNu64, &ln, &idx, &ck, &cnt) != 4 || ln != n) continue;
        if (ck != chunk_k) { fprintf(stderr, "log %s: decade %d was run with chunk %" PRIu64 " k-units, now %" PRIu64 "\n", path, n, ck, chunk_k); exit(1); }
        if (idx >= nchunks || bit_get(done, idx)) continue;
        bit_set(done, idx); *total += cnt; resumed++;
    }
    fclose(f);
    return resumed;
}

static int load_decades(const char *path, u64 *a_out)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    int best = 0; u64 a = 0; char line[256];
    while (fgets(line, sizeof line, f)) {
        int n; u64 v;
        if (!strchr(line, '\n')) continue;
        if (sscanf(line, "d %d %" SCNu64, &n, &v) == 2 && n > best) { best = n; a = v; }
    }
    fclose(f); *a_out = a; return best;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [options] N\n"
        "  N            compute a(1)..a(N), twin prime pairs below 10^n (N <= 20)\n"
        "  -s chunk     chunk size in numbers (default 10^max(10,n-4), capped at 10^14)\n"
        "  -f k -a A    start at decade k with a(k-1) = A\n"
        "  -l file      append-only chunk log (same format as the CPU tool); rerun to resume\n"
        "  -r A:B       only chunks A <= idx < B of a single decade\n"
        "  -v           per-chunk timing on stderr\n", prog);
    exit(2);
}

int main(int argc, char **argv)
{
    int from = 0, have_a = 0, have_r = 0, opt;
    u64 chunk_numbers = 0, a_known = 0, r_first = 0, r_last = UINT64_MAX;
    const char *logpath = NULL;
    while ((opt = getopt(argc, argv, "s:f:a:l:r:vh")) != -1) {
        switch (opt) {
        case 's': chunk_numbers = strtoull(optarg, NULL, 10); break;
        case 'f': from = atoi(optarg); break;
        case 'a': a_known = strtoull(optarg, NULL, 10); have_a = 1; break;
        case 'l': logpath = optarg; break;
        case 'r': if (sscanf(optarg, "%" SCNu64 ":%" SCNu64, &r_first, &r_last) != 2 || r_last <= r_first) usage(argv[0]); have_r = 1; break;
        case 'v': g_verbose = 1; break;
        default: usage(argv[0]);
        }
    }
    if (optind != argc - 1) usage(argv[0]);
    int N = atoi(argv[optind]);
    if (N < 1 || N > 20) { fprintf(stderr, "N must be 1..20\n"); return 2; }

    wheel_init();
    presieve_init();

    FILE *logf = NULL;
    if (logpath) { trim_partial_line(logpath); logf = fopen(logpath, "a"); if (!logf) { perror(logpath); return 1; } }

    u64 cumulative = 0; int start = 1;
    if (from > 0) {
        if (from > 1 && !have_a) {
            u64 v; int k = logpath ? load_decades(logpath, &v) : 0;
            if (k == from - 1) { a_known = v; have_a = 1; }
            else { fprintf(stderr, "-f %d needs -a a(%d)\n", from, from - 1); return 2; }
        }
        start = from; cumulative = from > 1 ? a_known : 0;
    } else if (logpath) {
        u64 v; int k = load_decades(logpath, &v);
        if (k > 0) { start = k + 1; cumulative = v; }
    }
    if (start > N) { fprintf(stderr, "nothing to do: a(%d) = %" PRIu64 " already known\n", start - 1, cumulative); return 0; }
    if (have_r && start != N) { fprintf(stderr, "-r needs a single decade: use -f %d (with -a) and N = %d\n", N, N); return 2; }

    printf("# A007508: twin prime pairs below 10^n  (GPU, seg=%u k, nseg=%d, qsplit=%u)\n", SEG_K, NSEG, QSPLIT);
    printf("# n  a(n)  [decade count]  [seconds]\n");
    if (start > 1) printf("# resuming at decade %d with a(%d) = %" PRIu64 "\n", start, start - 1, cumulative);
    fflush(stdout);
    double t_start = now();

    for (int n = start; n <= N; n++) {
        double t0 = now();
        u64 dec = 0, resumed = 0, ndone = 0, nchunks, chunk_k = 0, first = 0, last = 0;
        int complete = 1;
        if (n <= 4) {
            dec = small_decade_count(n); nchunks = 0;
        } else {
            u64 K0 = Kbound(n - 1), K1 = Kbound(n);
            chunk_k = chunk_numbers ? chunk_numbers / 30 : default_chunk_k(n);
            if (chunk_k < 64) chunk_k = 64;
            nchunks = (K1 - K0 + chunk_k - 1) / chunk_k;
            first = 0; last = nchunks;
            if (have_r) { first = r_first < nchunks ? r_first : nchunks; last = r_last < nchunks ? r_last : nchunks; }
            u8 *done = NULL;
            if (logpath) {
                done = (u8 *)calloc((nchunks + 7) / 8 + 1, 1);
                resumed = load_chunks(logpath, n, chunk_k, nchunks, done, &dec);
                if (resumed) printf("# decade %d: %" PRIu64 " of %" PRIu64 " chunks resumed from %s\n", n, resumed, nchunks, logpath);
            }
            for (u64 i = first; i < last; i++) {
                if (done && bit_get(done, i)) continue;
                u64 k_lo = K0 + i * chunk_k, k_hi = k_lo + chunk_k;
                if (k_hi > K1 || k_hi < k_lo) k_hi = K1;
                u64 c = gpu_count(k_lo, k_hi);
                dec += c; ndone++;
                if (logf) { fprintf(logf, "k %d %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", n, i, chunk_k, c); fflush(logf); }
            }
            complete = (resumed + ndone) == nchunks;
            free(done);
        }
        if (complete) {
            cumulative += dec;
            printf("%2d %" PRIu64 "  [%" PRIu64 "]  [%.2f]\n", n, cumulative, dec, now() - t0);
            if (logf) { fprintf(logf, "d %d %" PRIu64 "\n", n, cumulative); fflush(logf); }
        } else {
            printf("%2d partial  [%" PRIu64 " in chunks %" PRIu64 "..%" PRIu64 " of %" PRIu64 "; %" PRIu64 " chunks logged]  [%.2f]\n",
                   n, dec, first, last, nchunks, resumed + ndone, now() - t0);
        }
        fflush(stdout);
    }
    fprintf(stderr, "# total time %.2f s\n", now() - t_start);
    if (logf) fclose(logf);
    return 0;
}
