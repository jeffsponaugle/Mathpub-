/*
 * gapsieve_cuda.cu -- the CUDA half of gapsieve, for NVIDIA GPUs (written for
 * the DGX Spark's GB10).  The kernels are gapsieve.metal's, ported: for each
 * unit of consecutive segments and each shape, roffs_kernel works out where
 * each class starts in each pattern group, sieve_kernel ANDs the rotated
 * patterns in two phases and lists the candidates that survive, then a few
 * rounds of pass_kernel test one member each with the base-2 strong test,
 * and tail_kernel tests the members left.  The host (gapsieve_cuda_host.c)
 * proves the rare survivors on the CPU.
 *
 * Every unit goes on one stream, so units run one after another with the
 * next already queued: on Apple GPUs, units running side by side slowed each
 * other badly.  The screening kernels read their candidate count from device
 * memory and loop over it, so a unit needs no round trip to the host.
 *
 * Numbers stay below 2^82 (PSI13 in gapsieve.c), in three 32-bit limbs.
 */
#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gapsieve_cuda.h"

#define KBLOCKS (1u << 17)      /* candidates per class per segment: gapsieve.c */
#define WPT 4                   /* 32-bit words per thread, first sieve phase */
#define TPG 256                 /* threads per sieve block: 1024 words */
#define HLIST 4096              /* survivors copied back with every unit */

/* ---------- arithmetic: gapsieve.metal's, in CUDA ---------- */

struct u96 { uint32_t l0, l1, l2; };

__device__ static inline u96 to96(uint64_t lo, uint64_t hi)
{
    u96 r;
    r.l0 = (uint32_t)lo;
    r.l1 = (uint32_t)(lo >> 32);
    r.l2 = (uint32_t)hi;
    return r;
}

/* base + off + moff, as (lo, hi) */
__device__ static inline void add3(uint64_t blo, uint64_t bhi, uint64_t off,
                                   uint32_t moff, uint64_t *lo, uint64_t *hi)
{
    uint64_t l = blo + off, h = bhi + (l < off);
    uint64_t l2 = l + moff;
    h += l2 < l;
    *lo = l2;
    *hi = h;
}

__device__ static inline u96 add96(u96 a, u96 b)
{
    u96 r;
    uint64_t c = (uint64_t)a.l0 + b.l0;
    r.l0 = (uint32_t)c;
    c = (uint64_t)a.l1 + b.l1 + (c >> 32);
    r.l1 = (uint32_t)c;
    r.l2 = a.l2 + b.l2 + (uint32_t)(c >> 32);
    return r;
}

__device__ static inline u96 sub96(u96 a, u96 b)
{
    u96 r;
    uint64_t c = (uint64_t)a.l0 - b.l0;
    r.l0 = (uint32_t)c;
    c = (uint64_t)a.l1 - b.l1 - (c >> 63);
    r.l1 = (uint32_t)c;
    r.l2 = a.l2 - b.l2 - (uint32_t)(c >> 63);
    return r;
}

__device__ static inline bool ge96(u96 a, u96 b)
{
    if (a.l2 != b.l2)
        return a.l2 > b.l2;
    if (a.l1 != b.l1)
        return a.l1 > b.l1;
    return a.l0 >= b.l0;
}

__device__ static inline bool eq96(u96 a, u96 b)
{
    return a.l0 == b.l0 && a.l1 == b.l1 && a.l2 == b.l2;
}

/* a*b/2^96 mod n, less than 2n, for a, b < 2n (CIOS, three limbs) */
__device__ static inline u96 mmul3(u96 a, u96 b, u96 n, uint32_t ninv)
{
    uint32_t t0 = 0, t1 = 0, t2 = 0, t3 = 0;
    uint32_t bl[3] = { b.l0, b.l1, b.l2 };
#pragma unroll
    for (int i = 0; i < 3; i++) {
        uint32_t bi = bl[i];
        uint64_t x = (uint64_t)a.l0 * bi + t0;          /* t += a*b[i] */
        t0 = (uint32_t)x;
        x = (uint64_t)a.l1 * bi + t1 + (x >> 32);
        t1 = (uint32_t)x;
        x = (uint64_t)a.l2 * bi + t2 + (x >> 32);
        t2 = (uint32_t)x;
        x = (uint64_t)t3 + (x >> 32);
        t3 = (uint32_t)x;
        uint32_t t4 = (uint32_t)(x >> 32);
        uint32_t m = t0 * ninv;                         /* t = (t + m*n) / 2^32 */
        x = (uint64_t)m * n.l0 + t0;
        x = (uint64_t)m * n.l1 + t1 + (x >> 32);
        t0 = (uint32_t)x;
        x = (uint64_t)m * n.l2 + t2 + (x >> 32);
        t1 = (uint32_t)x;
        x = (uint64_t)t3 + (x >> 32);
        t2 = (uint32_t)x;
        t3 = t4 + (uint32_t)(x >> 32);
    }
    u96 r;                                              /* t3 is 0 here */
    r.l0 = t0;
    r.l1 = t1;
    r.l2 = t2;
    return r;
}

/* base-2 strong probable-prime test for odd n, 2 < n < 2^94 */
__device__ static bool sprp3(u96 n)
{
    uint32_t inv = n.l0;                                /* n^-1 mod 2^32 */
    for (int i = 0; i < 4; i++)
        inv *= 2 - n.l0 * inv;
    uint32_t ninv = 0 - inv;
    int top = n.l2 ? 95 - __clz(n.l2) : n.l1 ? 63 - __clz(n.l1) : 31 - __clz(n.l0);
    u96 one = { 0, 0, 0 };                              /* 2^96 mod n, by doubling */
    if (top >= 64)
        one.l2 = 1u << (top - 64);
    else if (top >= 32)
        one.l1 = 1u << (top - 32);
    else
        one.l0 = 1u << top;
    for (int i = top; i < 96; i++) {
        one = add96(one, one);
        if (ge96(one, n))
            one = sub96(one, n);
    }
    u96 nm1 = sub96(n, one);                            /* Montgomery form of -1 */
    u96 n2 = add96(n, n);

    u96 d = n;                                          /* n - 1 = d * 2^s */
    d.l0 &= ~1u;
    int s = d.l0 ? __ffs(d.l0) - 1 : d.l1 ? 31 + __ffs(d.l1) : 63 + __ffs(d.l2);
    uint64_t dlo = (uint64_t)d.l1 << 32 | d.l0, dhi = d.l2;
    dlo = s >= 64 ? dhi >> (s - 64) : (dlo >> s) | (s ? dhi << (64 - s) : 0);
    dhi = s >= 64 ? 0 : dhi >> s;
    int b = dhi ? 127 - __clzll(dhi) : 63 - __clzll(dlo);

    u96 r = add96(one, one);                            /* 2, below 2n */
    while (--b >= 0) {
        r = mmul3(r, r, n, ninv);
        uint64_t bit = b >= 64 ? (dhi >> (b - 64)) & 1 : (dlo >> b) & 1;
        u96 dbl = add96(r, r);
        if (ge96(dbl, n2))
            dbl = sub96(dbl, n2);
        r = bit ? dbl : r;
    }
    if (ge96(r, n))
        r = sub96(r, n);
    if (eq96(r, one) || eq96(r, nm1))
        return true;
    while (--s > 0) {
        r = mmul3(r, r, n, ninv);
        if (ge96(r, n))
            r = sub96(r, n);
        if (eq96(r, nm1))
            return true;
    }
    return false;
}

/* ---------- kernels ---------- */

__constant__ uint32_t c_moffs[2][CU_MAXK];              /* member offsets, per shape */

/* one thread per (group, class, segment): class_bitmap's R, where the class
   starts; bq holds each segment's base mod each prime */
__global__ void roffs_kernel(cu_params sp, const uint32_t *res, const uint32_t *bq,
                             const uint32_t *q, const uint32_t *qinv,
                             const uint32_t *crt, const uint4 *ginfo, uint32_t *roff)
{
    uint64_t id = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t total = (uint64_t)sp.nseg * sp.A * sp.ngroups;
    if (id >= total)
        return;
    uint32_t gi = (uint32_t)(id % sp.ngroups);
    uint64_t rest = id / sp.ngroups;
    uint32_t j = (uint32_t)(rest % sp.A), sg = (uint32_t)(rest / sp.A);
    uint32_t r = res[j];
    uint32_t s0 = r >= sp.bmodW ? r - sp.bmodW : r + sp.W - sp.bmodW;
    uint4 g = ginfo[gi];
    const uint32_t *b = bq + (size_t)sg * sp.nps;
    uint64_t R = 0;
    for (uint32_t i = g.y; i < g.y + g.z; i++) {
        uint32_t qq = q[i];
        uint32_t t = b[i] + s0 % qq;
        if (t >= qq)
            t -= qq;
        R += (uint64_t)((uint64_t)qinv[i] * t % qq) * crt[i];
    }
    roff[((size_t)sg * sp.A + j) * sp.ngroups + gi] = (uint32_t)(R % g.x);
}

/*
 * One block per 1024 words of one class of one segment.  Phase one: each
 * thread ANDs the first sp.dense rotated patterns into its WPT words.  The
 * first patterns hold the smallest primes and strike the most, so the words
 * still standing are then gathered in shared memory (one atomic per warp),
 * and phase two gives each of them a thread for the remaining patterns.
 * Surviving candidates go on the list as offsets from the unit's base; the
 * count keeps rising past cap so the host can see an overflow.
 */
__global__ void __launch_bounds__(TPG)
sieve_kernel(cu_params sp, const uint32_t *res, const uint32_t *pat32,
             const uint4 *ginfo, const uint32_t *roff, uint64_t *cand, uint32_t *count)
{
    __shared__ uint32_t tw[TPG * WPT], tv[TPG * WPT];
    __shared__ uint32_t tn;
    uint32_t lid = threadIdx.x, lane = lid & 31;
    uint32_t w0 = (blockIdx.x * TPG + lid) * WPT, j = blockIdx.y, sg = blockIdx.z;
    if (lid == 0)
        tn = 0;
    const uint32_t *ro = roff + ((size_t)sg * sp.A + j) * sp.ngroups;
    uint32_t acc[WPT];
#pragma unroll
    for (int t = 0; t < WPT; t++)
        acc[t] = ~0u;
    for (uint32_t gi = 0; gi < sp.dense; gi++) {
        uint32_t R = ro[gi];
        const uint32_t *src = pat32 + 2 * (size_t)ginfo[gi].w + (R >> 5) + w0;
        uint32_t sh = R & 31, cur = src[0];
#pragma unroll
        for (int t = 0; t < WPT; t++) {
            uint32_t nxt = src[t + 1];
            acc[t] &= (cur >> sh) | ((nxt << 1) << (31 - sh));
            cur = nxt;
        }
    }
    __syncthreads();
    uint32_t nz = 0;                        /* gather: one atomic per warp */
#pragma unroll
    for (int t = 0; t < WPT; t++)
        nz += acc[t] != 0;
    uint32_t x = nz;
#pragma unroll
    for (int o = 1; o < 32; o <<= 1) {
        uint32_t y = __shfl_up_sync(0xffffffffu, x, o);
        if (lane >= (uint32_t)o)
            x += y;
    }
    uint32_t tot = __shfl_sync(0xffffffffu, x, 31), at = 0;
    if (lane == 0 && tot)
        at = atomicAdd(&tn, tot);
    uint32_t k = __shfl_sync(0xffffffffu, at, 0) + x - nz;
#pragma unroll
    for (int t = 0; t < WPT; t++)
        if (acc[t]) {
            tw[k] = w0 + t;
            tv[k++] = acc[t];
        }
    __syncthreads();

    uint32_t n = tn, r = res[j];
    uint64_t s0 = (uint64_t)sg * sp.segsize +
                  (r >= sp.bmodW ? r - sp.bmodW : r + sp.W - sp.bmodW);
    for (uint32_t i = lid; i < n; i += TPG) {
        uint32_t w = tw[i], a = tv[i];
        for (uint32_t gi = sp.dense; gi < sp.ngroups && a; gi++) {
            uint32_t R = ro[gi];
            const uint32_t *src = pat32 + 2 * (size_t)ginfo[gi].w + (R >> 5) + w;
            uint32_t sh = R & 31;
            a &= (src[0] >> sh) | ((src[1] << 1) << (31 - sh));
        }
        for (; a; a &= a - 1) {
            uint64_t off = s0 + (uint64_t)(w * 32 + __ffs(a) - 1) * sp.W;
            if (off <= sp.lim) {
                uint32_t idx = atomicAdd(count, 1);
                if (idx < sp.cap)
                    cand[idx] = off;
            }
        }
    }
}

/* member m for every remaining candidate; keep those whose member passes */
__global__ void pass_kernel(cu_params sp, const uint64_t *in, const uint32_t *cnt_in,
                            uint64_t *out, uint32_t *cnt_out, int d, uint32_t m)
{
    uint32_t n = min(*cnt_in, sp.cap), moff = c_moffs[d][m];
    for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += gridDim.x * blockDim.x) {
        uint64_t off = in[i], lo, hi;
        add3(sp.base_lo, sp.base_hi, off, moff, &lo, &hi);
        if (sprp3(to96(lo, hi)))
            out[atomicAdd(cnt_out, 1)] = off;
    }
}

/* every member left, for each remaining candidate, stopping at the first
   that fails: by now few candidates remain */
__global__ void tail_kernel(cu_params sp, const uint64_t *in, const uint32_t *cnt_in,
                            uint64_t *out, uint32_t *cnt_out, int d, uint32_t m0,
                            uint32_t nm)
{
    uint32_t n = min(*cnt_in, sp.cap);
    for (uint32_t i = blockIdx.x * blockDim.x + threadIdx.x; i < n;
         i += gridDim.x * blockDim.x) {
        uint64_t off = in[i], lo, hi;
        bool ok = true;
        for (uint32_t m = 0; m < nm && ok; m++) {
            add3(sp.base_lo, sp.base_hi, off, c_moffs[d][m0 + m], &lo, &hi);
            ok = sprp3(to96(lo, hi));
        }
        if (ok)
            out[atomicAdd(cnt_out, 1)] = off;
    }
}

__global__ void sprp_kernel(const uint64_t *n, size_t count, uint8_t *ok)
{
    size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < count)
        ok[i] = sprp3(to96(n[2 * i], n[2 * i + 1])) ? 1 : 0;
}

/* ---------- device memory and the C interface ---------- */

static cudaStream_t g_stream;
static int g_blocks;                    /* grid for the screening kernels */
static uint32_t g_ngroups, g_nps, g_A[2];
static uint4 *d_ginfo;
static uint32_t *d_q, *d_qinv, *d_crt, *d_pat[2], *d_res[2];

static struct slot {
    uint32_t unit, cap;
    uint32_t *bq, *roff[2], *cnt;       /* cnt: per shape, c0 c1 stat */
    uint64_t *l0[2], *l1[2];
    uint32_t *h_bq, *h_cnt;             /* pinned: stat[2], nfin[2] */
    uint64_t *h_list[2];                /* pinned: HLIST survivors per shape */
    int fin[2];                         /* the survivors are in list 1 */
    cudaEvent_t done;
    const uint64_t *dfin[2];
} S[CU_MAXSLOT];

#define CK(x) do { if ((x) != cudaSuccess) return -1; } while (0)

int cu_init(char *name, size_t n)
{
    cudaDeviceProp p;
    int dev = 0;
    CK(cudaSetDevice(dev));
    CK(cudaGetDeviceProperties(&p, dev));
    snprintf(name, n, "%s", p.name);
    g_blocks = p.multiProcessorCount * (p.maxThreadsPerMultiProcessor / 256);
    CK(cudaStreamCreateWithFlags(&g_stream, cudaStreamNonBlocking));
    return 0;
}

int cu_tables(uint32_t ngroups, const uint32_t *ginfo, uint32_t nps,
              const uint32_t *q, const uint32_t *qinv, const uint32_t *crt)
{
    g_ngroups = ngroups;
    g_nps = nps;
    CK(cudaMalloc(&d_ginfo, ngroups * sizeof(uint4)));
    CK(cudaMemcpy(d_ginfo, ginfo, ngroups * sizeof(uint4), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_q, nps * 4));
    CK(cudaMalloc(&d_qinv, nps * 4));
    CK(cudaMalloc(&d_crt, nps * 4));
    CK(cudaMemcpy(d_q, q, nps * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_qinv, qinv, nps * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_crt, crt, nps * 4, cudaMemcpyHostToDevice));
    return 0;
}

int cu_shape(int d, const uint32_t *pat, size_t npat, const uint32_t *res,
             uint32_t A, const uint32_t *offs, uint32_t k)
{
    if (k > CU_MAXK)
        return -1;
    g_A[d] = A;
    CK(cudaMalloc(&d_pat[d], npat * 4 + 64));
    CK(cudaMemcpy(d_pat[d], pat, npat * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_res[d], A * 4));
    CK(cudaMemcpy(d_res[d], res, A * 4, cudaMemcpyHostToDevice));
    CK(cudaMemcpyToSymbol(c_moffs, offs, k * 4, (size_t)d * CU_MAXK * 4));
    return 0;
}

int cu_grow(int k, uint32_t cap)
{
    struct slot *X = &S[k];
    for (int d = 0; d < 2; d++) {
        cudaFree(X->l0[d]);
        cudaFree(X->l1[d]);
        X->l0[d] = X->l1[d] = NULL;
        if (!d_pat[d])
            continue;
        CK(cudaMalloc(&X->l0[d], (size_t)cap * 8));
        CK(cudaMalloc(&X->l1[d], (size_t)cap * 8));
    }
    X->cap = cap;
    return 0;
}

int cu_slot(int k, uint32_t unit, uint32_t cap)
{
    struct slot *X = &S[k];
    X->unit = unit;
    CK(cudaMalloc(&X->bq, (size_t)unit * g_nps * 4));
    CK(cudaMallocHost(&X->h_bq, (size_t)unit * g_nps * 4));
    CK(cudaMalloc(&X->cnt, 2 * 4 * 4));
    CK(cudaMallocHost(&X->h_cnt, 4 * 4));
    for (int d = 0; d < 2; d++) {
        if (!d_pat[d])
            continue;
        CK(cudaMalloc(&X->roff[d], (size_t)unit * g_A[d] * g_ngroups * 4));
        CK(cudaMallocHost(&X->h_list[d], HLIST * 8));
    }
    CK(cudaEventCreateWithFlags(&X->done, cudaEventDisableTiming));
    return cu_grow(k, cap);
}

int cu_submit(int k, const uint32_t *bq, uint32_t nseg, uint32_t nps,
              const struct cu_params *p[2], int split)
{
    struct slot *X = &S[k];
    memcpy(X->h_bq, bq, (size_t)nseg * nps * 4);
    CK(cudaMemcpyAsync(X->bq, X->h_bq, (size_t)nseg * nps * 4,
                       cudaMemcpyHostToDevice, g_stream));
    for (int d = 0; d < 2; d++) {
        if (!p[d])
            continue;
        cu_params sp = *p[d];
        uint32_t *c[2] = { X->cnt + 4 * d, X->cnt + 4 * d + 1 }, *stat = X->cnt + 4 * d + 2;
        uint64_t *l[2] = { X->l0[d], X->l1[d] };
        uint64_t total = (uint64_t)sp.nseg * sp.A * sp.ngroups;
        roffs_kernel<<<(unsigned)((total + 255) / 256), 256, 0, g_stream>>>(
            sp, d_res[d], X->bq, d_q, d_qinv, d_crt, d_ginfo, X->roff[d]);
        CK(cudaMemsetAsync(c[0], 0, 4, g_stream));
        sieve_kernel<<<dim3(KBLOCKS / 32 / WPT / TPG, sp.A, sp.nseg), TPG, 0, g_stream>>>(
            sp, d_res[d], d_pat[d], d_ginfo, X->roff[d], l[0], c[0]);
        CK(cudaMemcpyAsync(stat, c[0], 4, cudaMemcpyDeviceToDevice, g_stream));
        int in = 0;
        for (uint32_t m = 0; m < sp.k; in ^= 1) {
            CK(cudaMemsetAsync(c[in ^ 1], 0, 4, g_stream));
            if ((int)m < split) {
                pass_kernel<<<g_blocks, 256, 0, g_stream>>>(sp, l[in], c[in], l[in ^ 1],
                                                              c[in ^ 1], d, m);
                m++;
            } else {
                tail_kernel<<<g_blocks, 256, 0, g_stream>>>(sp, l[in], c[in], l[in ^ 1],
                                                              c[in ^ 1], d, m, sp.k - m);
                m = sp.k;
            }
        }
        X->fin[d] = in;                 /* the last round's output */
        CK(cudaMemcpyAsync(X->h_cnt + d, stat, 4, cudaMemcpyDeviceToHost, g_stream));
        CK(cudaMemcpyAsync(X->h_cnt + 2 + d, c[in], 4, cudaMemcpyDeviceToHost, g_stream));
        CK(cudaMemcpyAsync(X->h_list[d], l[in], HLIST * 8 < (size_t)sp.cap * 8 ? HLIST * 8
                           : (size_t)sp.cap * 8, cudaMemcpyDeviceToHost, g_stream));
    }
    CK(cudaGetLastError());
    CK(cudaEventRecord(X->done, g_stream));
    return 0;
}

int cu_wait(int k, uint32_t stat[2], uint32_t nfin[2])
{
    struct slot *X = &S[k];
    CK(cudaEventSynchronize(X->done));
    CK(cudaGetLastError());
    for (int d = 0; d < 2; d++) {
        stat[d] = X->h_cnt[d];
        nfin[d] = X->h_cnt[2 + d];
    }
    return 0;
}

const uint64_t *cu_list(int k, int d, uint32_t n)
{
    struct slot *X = &S[k];
    if (n > HLIST) {                    /* rare: many survivors (short chains) */
        static uint64_t *big[2];
        static size_t nbig[2];
        if (n > nbig[d]) {
            free(big[d]);
            big[d] = (uint64_t *)malloc((size_t)n * 8);
            nbig[d] = big[d] ? n : 0;
            if (!big[d])
                return NULL;
        }
        uint64_t *src = X->fin[d] ? X->l1[d] : X->l0[d];
        if (cudaMemcpy(big[d], src, (size_t)n * 8, cudaMemcpyDeviceToHost) != cudaSuccess)
            return NULL;
        return big[d];
    }
    return X->h_list[d];
}

int cu_sprp_batch(const uint64_t *n, size_t count, uint8_t *ok)
{
    uint64_t *dn;
    uint8_t *dok;
    CK(cudaMalloc(&dn, count * 16));
    CK(cudaMalloc(&dok, count));
    CK(cudaMemcpy(dn, n, count * 16, cudaMemcpyHostToDevice));
    sprp_kernel<<<(unsigned)((count + 255) / 256), 256>>>(dn, count, dok);
    CK(cudaGetLastError());
    CK(cudaMemcpy(ok, dok, count, cudaMemcpyDeviceToHost));
    cudaFree(dn);
    cudaFree(dok);
    return 0;
}
