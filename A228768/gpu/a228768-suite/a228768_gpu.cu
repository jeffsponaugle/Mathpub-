/*
 * a228768_gpu.cu -- GPU sieves for the A228768 a(15) search.
 *
 * Split of work:
 *   GPU   both segmented sieves (x-side and base-2 reversal side) plus
 *         stream compaction of the survivors.  All 32/64-bit index-space
 *         arithmetic, no 128-bit values ever reach the device.
 *   CPU   the cascade over bases 3..15 and the deterministic Miller-Rabin
 *         confirmation, via a228768_check() from a228768x.c -- one source
 *         of truth for the number theory.
 *
 * That split is roughly balanced on a DGX Spark: the sieves emit about 1.1%
 * of odd numbers as survivors, and the 20 Arm cores can chew through those
 * at about the rate the GPU produces them.
 *
 * LIBRARIES: nothing beyond the CUDA Toolkit.  No Thrust, no cuB, no
 * third-party dependency -- compaction is a warp-aggregated atomic, which
 * is a dozen lines and avoids version-coupling the build.
 *
 * Build on DGX Spark (GB10 is compute capability 12.1 = sm_121, aarch64):
 *     cc  -O3 -c -DA228768_LIB a228768x.c -o a228768lib.o
 *     nvcc -O3 -arch=sm_121 -Xcompiler -fopenmp \
 *          a228768_gpu.cu a228768lib.o -o a228768gpu -lm
 *
 * Build WITHOUT a GPU, to check the algorithm (kernels run serially on the
 * host; results are bit-for-bit identical):
 *     cc  -O3 -c -DA228768_LIB a228768x.c -o a228768lib.o
 *     c++ -O3 -DEMULATE -fopenmp a228768_gpu.cu a228768lib.o -o a228768emu -lm
 *
 * Always run --verify first on new hardware.  It recomputes one block's
 * survivors by honest trial division and compares against the GPU output.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

static double wall(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

typedef unsigned long long u64;
typedef unsigned int       u32;
typedef unsigned __int128  u128;          /* host side only */

/* ------------------------------------------------------------------ */
/* host <-> CPU-library interface                                     */
/* ------------------------------------------------------------------ */

extern "C" {
void a228768_setup(int nbase, unsigned long sieve_limit);
int  a228768_check(u128 x);
unsigned long a228768_nprimes(void);
const u32 *a228768_prime_list(void);
int  a228768_next_admissible(u128 x, u128 limit, u128 *out);
void a228768_admissible_end(u128 x, u128 *out);
}

/* ------------------------------------------------------------------ */
/* CUDA / emulation shims                                             */
/* ------------------------------------------------------------------ */

#ifdef EMULATE
  #define DEVFN  static inline
  #define HOSTDEV static inline
  static inline u64 dev_umul64hi(u64 a, u64 b) { return (u64)(((u128)a * b) >> 64); }
  static inline u32 dev_atomicOr(u32 *p, u32 v) { u32 o = *p; *p |= v; return o; }
  static inline u32 dev_atomicAdd(u32 *p, u32 v) { u32 o = *p; *p += v; return o; }
  #define DEV_CTZ(x) __builtin_ctz(x)
  static inline u64 dev_brevll(u64 v) {
      v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
      v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
      v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
      return __builtin_bswap64(v);
  }
  #define CUDA_OK(x) do { (void)(x); } while (0)
  typedef int cuda_stream_t;
  static inline void stream_create(cuda_stream_t *s) { *s = 0; }
  static inline void stream_sync(cuda_stream_t s) { (void)s; }
#else
  #include <cuda_runtime.h>
  #define DEVFN  __device__ inline
  #define HOSTDEV __host__ __device__ inline
  #define dev_umul64hi __umul64hi
  #define dev_atomicOr atomicOr
  #define dev_atomicAdd atomicAdd
  #define dev_brevll __brevll
  #define DEV_CTZ(x) (__ffs(x) - 1)
  #define CUDA_OK(call) do {                                              \
      cudaError_t e_ = (call);                                            \
      if (e_ != cudaSuccess) {                                            \
          fprintf(stderr, "CUDA error %s at %s:%d\n",                     \
                  cudaGetErrorString(e_), __FILE__, __LINE__);            \
          exit(1);                                                        \
      }                                                                   \
  } while (0)
  typedef cudaStream_t cuda_stream_t;
  static inline void stream_create(cuda_stream_t *s) { CUDA_OK(cudaStreamCreate(s)); }
  static inline void stream_sync(cuda_stream_t s) { CUDA_OK(cudaStreamSynchronize(s)); }
#endif

static unsigned g_gputhreads = 512;   /* threads per CUDA block */

/* odds per CUDA block; 2^18 bits = 32 KB of shared memory */
#define SUBLOG 18
#define SUB    (1u << SUBLOG)
#define SUBW   (SUB / 32)

/* ------------------------------------------------------------------ */
/* device helpers                                                     */
/* ------------------------------------------------------------------ */

/* a mod p, using M = floor(2^64/p).  The quotient estimate is short by at
 * most two, so the correction loop runs at most twice. */
DEVFN u64 fmod64(u64 a, u32 p, u64 M)
{
    u64 q = dev_umul64hi(a, M);
    u64 r = a - q * (u64)p;
    while (r >= p) r -= p;
    return r;
}

DEVFN u64 mulmod_small(u64 a, u64 b, u32 p, u64 M)
{
    return fmod64(a * b, p, M);           /* a,b < p < 2^20, product < 2^40 */
}

DEVFN u64 powmod_small(u64 a, u64 e, u32 p, u64 M)
{
    u64 r = 1;
    a = fmod64(a, p, M);
    while (e) {
        if (e & 1) r = mulmod_small(r, a, p, M);
        a = mulmod_small(a, a, p, M);
        e >>= 1;
    }
    return r;
}

#define REVPAD_CASE(B) \
    case B: { u64 r = 0; for (u32 i = 0; i < k; i++) { u64 d = L % (B); L /= (B); r = r * (B) + d; } return r; }

/* reverse the low k base-b digits of L (zero padded) */
DEVFN u64 revpad(u64 L, u32 b, u32 k)
{
    switch (b) {
    REVPAD_CASE(3)  REVPAD_CASE(4)  REVPAD_CASE(5)  REVPAD_CASE(6)
    REVPAD_CASE(7)  REVPAD_CASE(8)  REVPAD_CASE(9)  REVPAD_CASE(10)
    REVPAD_CASE(11) REVPAD_CASE(12) REVPAD_CASE(13) REVPAD_CASE(14)
    REVPAD_CASE(15)
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* kernel bodies, written as (block, thread) functions so the same code  */
/* runs on the device and, serially, on the host                        */
/* ------------------------------------------------------------------ */

/*
 * GPU-side cascade stages.
 *
 * Over a block, x = H*b^k + L with H fixed, so
 *     rev_b(x) = revpad_k(L) * b^ndigits(H) + rev_b(H)  =  A*P + C
 * where A < b^k fits in 64 bits and P, C are per-block 128-bit constants.
 * For a small prime q,  q | rev_b(x)  <=>  A == T (mod q)  with
 *     T = -C * P^-1  (mod q)
 * precomputed on the host.  So the device never touches a 128-bit value:
 * it reverses the low digits into A and compares residues.
 *
 * Rejecting only when a small q divides the reversal is always safe -- a
 * genuine emirp has that reversal prime, hence not divisible by q.
 */
#define MAXQ  16
#define MAXGB 10

#define TMODE_CMP    0     /* reject when A == T (mod q) */
#define TMODE_ALWAYS 1     /* q | P and q | C: every A is divisible */
#define TMODE_NEVER  2     /* q | P, q does not divide C: never divisible */

struct BaseCtx {
    u32 b, k;
    u64 bk, Bm;            /* b^k, and base mod b^k */
    u32 T[2][MAXQ];        /* [wrapped] -- L may roll into H+1 mid-block */
    u32 mode[2][MAXQ];
};

struct SieveArgs {
    const u32 *p;          /* sieve primes */
    const u64 *mag;        /* floor(2^64/p) */
    const u32 *r64;        /* 2^64 mod p */
    u32 *cinv;             /* 2^-t mod p  (per window) */
    u32 *xoff;             /* first odd-index of a multiple, per prime */
    u32 *roff;             /* first R-index (relative to half), per prime */
    u32 *rbm;              /* R-side bitmap, half bits */
    u32 *out;              /* sieve survivors: odd-indices within the block */
    u32 *outn;
    u32 *out2, *outn2;     /* survivors of the device cascade */
    u64 nprimes;
    u64 half;              /* 2^(s-1): odd positions in the block */
    u64 baseh, basel;      /* block base as two 64-bit halves */
    u64 K;                 /* bitrev_t(base >> s) */
    int t, s;
    u32 outcap;
    /* pre-sieve: primes 3..19 are ~40% of all marks and are supplied as a
     * repeating pattern instead of being marked one bit at a time.  The
     * x-side pattern is fixed (odd n divisible by p sits at index (p-1)/2
     * mod p); the R-side residues depend on K, so that one is regenerated
     * per block by k_rpat. */
    const u32 *xpat;
    u32 *rpat;
    u64 patbits, xphase;
    int npre;
    struct BaseCtx *bctx;   /* cascade stages run on the device */
    int ngb, nq;
    const u32 *q;
    const u64 *qinv, *qlim;
};

/* one thread per prime: 2^-t mod p */
DEVFN void body_cinv(u64 i, const SieveArgs a)
{
    if (i >= a.nprimes) return;
    u32 p = a.p[i]; u64 M = a.mag[i];
    u64 c = powmod_small(2, (u64)a.t % (p - 1), p, M);
    a.cinv[i] = (u32)powmod_small(c, p - 2, p, M);   /* inverse of 2^t */
}

/* one thread per prime: starting offsets for both sieves in this block */
DEVFN void body_prep(u64 i, const SieveArgs a)
{
    if (i >= a.nprimes) return;
    u32 p = a.p[i]; u64 M = a.mag[i];

    /* x-side.  Odd position j holds x = base + 2j + 1, so
       p | x  <=>  j == -(base+1)/2 (mod p). */
    u64 bh = fmod64(a.baseh, p, M);
    u64 bl = fmod64(a.basel + 1, p, M);
    u64 bp = fmod64(bh * (u64)a.r64[i] + bl, p, M);
    u64 inv2 = (p + 1) / 2;
    a.xoff[i] = (u32)mulmod_small(bp ? p - bp : 0, inv2, p, M);

    /* R-side.  rev = R*2^t + K, so p | rev  <=>  R == -K * 2^-t (mod p). */
    u64 kp = fmod64(a.K, p, M);
    u64 R0 = mulmod_small(kp ? p - kp : 0, a.cinv[i], p, M);
    u64 hp = fmod64(a.half, p, M);
    u32 ro = (u32)(R0 + p - hp);
    a.roff[i] = ro >= p ? ro - p : ro;
}

DEVFN u32 getbits32(const u32 *pat, u64 bp)
{
    u64 w = bp >> 5;
    u32 o = (u32)(bp & 31);
    return o ? ((pat[w] >> o) | (pat[w + 1] << (32 - o))) : pat[w];
}

/* seed the sub-segment: pattern if we have one, zeros otherwise */
DEVFN void body_fill(u32 tid, u32 nthr, u32 *sm, const u32 *pat,
                     u64 patbits, u64 phase0, u64 first)
{
    if (!pat) { for (u32 w = tid; w < SUBW; w += nthr) sm[w] = 0; return; }
    u64 ph = (phase0 + first) % patbits;
    for (u32 w = tid; w < SUBW; w += nthr) {
        u64 bp = ph + (u64)w * 32;
        if (bp >= patbits) bp %= patbits;
        sm[w] = getbits32(pat, bp);
    }
}

/* rebuild the R-side pattern for this block: three phases, so three
   launches -- the cost is trivial against a 2^27-odd block */
DEVFN void body_rpat_clear(u64 gid, u64 stride, const SieveArgs a)
{
    u64 words = (a.patbits + 31) / 32 + 2;
    for (u64 w = gid; w < words; w += stride) a.rpat[w] = 0;
}
DEVFN void body_rpat_mark(u64 gid, u64 stride, const SieveArgs a)
{
    for (int i = 0; i < a.npre; i++) {
        u32 p = a.p[i];
        for (u64 j = (u64)a.roff[i] + gid * p; j < a.patbits; j += stride * p)
            dev_atomicOr(&a.rpat[j >> 5], 1u << (j & 31));
    }
}
DEVFN void body_rpat_pad(u64 gid, const SieveArgs a)
{
    if (gid) return;                      /* replicate 32 bits past the end */
    for (u32 b = 0; b < 32; b++)
        if ((a.rpat[b >> 5] >> (b & 31)) & 1) {
            u64 d = a.patbits + b;
            a.rpat[d >> 5] |= 1u << (d & 31);
        }
}

/* mark composites into a shared-memory sub-segment.
 * off[] holds each prime's first position relative to the block start;
 * `first` is this sub-segment's position within the block. */
DEVFN void body_mark(u32 tid, u32 nthr, u32 *sm, const SieveArgs a, const u32 *off, u64 first)
{
    /* primes below index npre came in with the pattern */
    /* Work per prime spans four orders of magnitude: prime 3 has SUB/3
     * multiples in this sub-segment, a prime above SUB has at most one.
     * Parallelising over primes therefore leaves one lane doing ~36x the
     * average.  So: small primes are marked cooperatively by the whole
     * block (each lane takes every nthr'th multiple), and only once a
     * prime is sparse enough do we switch to one-prime-per-lane.
     * The split point is uniform across the block, so it never diverges. */
    u64 i = a.npre;
    for (; i < a.nprimes; i++) {
        u32 p = a.p[i];
        if ((u64)p * nthr >= SUB) break;
        u32 g = (u32)fmod64(first, p, a.mag[i]);
        u32 j = off[i] + p - g;
        if (j >= p) j -= p;                       /* both < p: one subtract */
        for (u32 k = j + tid * p; k < SUB; k += nthr * p)
            dev_atomicOr(&sm[k >> 5], 1u << (k & 31));
    }
    for (i += tid; i < a.nprimes; i += nthr) {
        u32 p = a.p[i];
        u32 g = (u32)fmod64(first, p, a.mag[i]);
        u32 j = off[i] + p - g;
        if (j >= p) j -= p;
        for (; j < SUB; j += p) dev_atomicOr(&sm[j >> 5], 1u << (j & 31));
    }
}

DEVFN void body_store(u32 tid, u32 nthr, const u32 *sm, u32 *dst)
{
    for (u32 w = tid; w < SUBW; w += nthr) dst[w] = sm[w];
}

/* combine: a survivor is an odd x the x-sieve left clear whose base-2
 * reversal the R-sieve also left clear.  Emits the odd-index within block. */
/* 1 if every device-side stage clears (or there are none) */
DEVFN int gpu_cascade(u64 off, const SieveArgs a)
{
    for (int s = 0; s < a.ngb; s++) {
        const struct BaseCtx *c = &a.bctx[s];
        u64 L = c->Bm + off;
        int w = 0;
        if (L >= c->bk) { L -= c->bk; w = 1; }
        u64 A = revpad(L, c->b, c->k);
        for (int j = 0; j < a.nq; j++) {
            u32 mode = c->mode[w][j];
            if (mode == TMODE_NEVER) continue;
            if (mode == TMODE_ALWAYS) return 0;
            u64 v = A + a.q[j] - c->T[w][j];
            if (v * a.qinv[j] <= a.qlim[j]) return 0;   /* q divides it */
        }
    }
    return 1;
}

DEVFN void body_emit(u32 tid, u32 nthr, const u32 *sm, const SieveArgs a, u64 first)
{
    u64 half = a.half;
    int s = a.s;
    for (u32 w = tid; w < SUBW; w += nthr) {
        u32 bits = ~sm[w];
        while (bits) {
            u32 b = (u32)DEV_CTZ(bits);
            bits &= bits - 1;
            u64 idx = first + (u64)w * 32 + b;         /* odd index in block */
            u64 L = 2 * idx + 1;                       /* offset from base */
            u64 R = dev_brevll(L) >> (64 - s);
            u64 ridx = R - half;
            if (a.rbm[ridx >> 5] & (1u << (ridx & 31))) continue;
            u32 pos = dev_atomicAdd(a.outn, 1u);
            if (pos < a.outcap) a.out[pos] = (u32)idx;
        }
    }
}

/*
 * Cascade over the compacted survivor list.
 *
 * Running this inside body_emit was a mistake: that walks the sieve bitmap,
 * where ~1% of bits are survivors, so about one lane per warp had a
 * candidate and the other 31 idled through the whole reversal.  Over a
 * dense list every lane has work.
 */
DEVFN void body_cascade(u64 gid, u64 stride, const SieveArgs a)
{
    u32 n = *a.outn;
    for (u64 i = gid; i < n; i += stride) {
        u32 idx = a.out[i];
        if (!gpu_cascade(2 * (u64)idx + 1, a)) continue;
        u32 pos = dev_atomicAdd(a.outn2, 1u);
        if (pos < a.outcap) a.out2[pos] = idx;
    }
}

/* ------------------------------------------------------------------ */
/* launchers                                                          */
/* ------------------------------------------------------------------ */

#ifdef EMULATE

static u32 g_sm[SUBW];

static void run_cinv(SieveArgs a, cuda_stream_t st)  { (void)st; for (u64 i = 0; i < a.nprimes; i++) body_cinv(i, a); }
static void run_prep(SieveArgs a, cuda_stream_t st)  { (void)st; for (u64 i = 0; i < a.nprimes; i++) body_prep(i, a); }

/* emulate with the real block width so the same code path is exercised;
 * marking is idempotent OR, so serial lanes give an identical bitmap */
static void run_rsieve(SieveArgs a, cuda_stream_t st)
{
    (void)st;
    u64 nsub = a.half / SUB;
    for (u64 sub = 0; sub < nsub; sub++) {
        for (u32 t = 0; t < g_gputhreads; t++) body_fill(t, g_gputhreads, g_sm, a.rpat, a.patbits, 0, sub * SUB);
        for (u32 t = 0; t < g_gputhreads; t++) body_mark(t, g_gputhreads, g_sm, a, a.roff, sub * SUB);
        for (u32 t = 0; t < g_gputhreads; t++) body_store(t, g_gputhreads, g_sm, a.rbm + sub * SUBW);
    }
}

static void run_cascade_gpu(SieveArgs a, cuda_stream_t st)
{
    (void)st;
    body_cascade(0, 1, a);
}

static void run_rpat(SieveArgs a, cuda_stream_t st)
{
    (void)st;
    if (!a.npre) return;
    body_rpat_clear(0, 1, a);
    body_rpat_mark(0, 1, a);
    body_rpat_pad(0, a);
}

static void run_xsieve(SieveArgs a, cuda_stream_t st)
{
    (void)st;
    u64 nsub = a.half / SUB;
    for (u64 sub = 0; sub < nsub; sub++) {
        for (u32 t = 0; t < g_gputhreads; t++) body_fill(t, g_gputhreads, g_sm, a.xpat, a.patbits, a.xphase, sub * SUB);
        for (u32 t = 0; t < g_gputhreads; t++) body_mark(t, g_gputhreads, g_sm, a, a.xoff, sub * SUB);
        for (u32 t = 0; t < g_gputhreads; t++) body_emit(t, g_gputhreads, g_sm, a, sub * SUB);
    }
}
static void dev_sync(void) {}

#else

__global__ void k_cinv(SieveArgs a)
{
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    body_cinv(i, a);
}
__global__ void k_prep(SieveArgs a)
{
    u64 i = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    body_prep(i, a);
}
__global__ void k_rsieve(SieveArgs a)
{
    __shared__ u32 sm[SUBW];
    u64 first = (u64)blockIdx.x * SUB;
    body_fill(threadIdx.x, blockDim.x, sm, a.rpat, a.patbits, 0, first);
    __syncthreads();
    body_mark(threadIdx.x, blockDim.x, sm, a, a.roff, first);
    __syncthreads();
    body_store(threadIdx.x, blockDim.x, sm, a.rbm + (u64)blockIdx.x * SUBW);
}
__global__ void k_xsieve(SieveArgs a)
{
    __shared__ u32 sm[SUBW];
    u64 first = (u64)blockIdx.x * SUB;
    body_fill(threadIdx.x, blockDim.x, sm, a.xpat, a.patbits, a.xphase, first);
    __syncthreads();
    body_mark(threadIdx.x, blockDim.x, sm, a, a.xoff, first);
    __syncthreads();
    body_emit(threadIdx.x, blockDim.x, sm, a, first);
}

static void run_cinv(SieveArgs a, cuda_stream_t st)
{
    k_cinv<<<(unsigned)((a.nprimes + 255) / 256), 256, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}
static void run_prep(SieveArgs a, cuda_stream_t st)
{
    k_prep<<<(unsigned)((a.nprimes + 255) / 256), 256, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}
static void run_rsieve(SieveArgs a, cuda_stream_t st)
{
    k_rsieve<<<(unsigned)(a.half / SUB), g_gputhreads, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}
static void run_xsieve(SieveArgs a, cuda_stream_t st)
{
    k_xsieve<<<(unsigned)(a.half / SUB), g_gputhreads, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}

__global__ void k_rpat_clear(SieveArgs a)
{
    u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    body_rpat_clear(gid, (u64)gridDim.x * blockDim.x, a);
}
__global__ void k_rpat_mark(SieveArgs a)
{
    u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    body_rpat_mark(gid, (u64)gridDim.x * blockDim.x, a);
}
__global__ void k_rpat_pad(SieveArgs a)
{
    body_rpat_pad((u64)blockIdx.x * blockDim.x + threadIdx.x, a);
}
static void run_rpat(SieveArgs a, cuda_stream_t st)
{
    if (!a.npre) return;
    k_rpat_clear<<<256, 256, 0, st>>>(a);
    k_rpat_mark<<<256, 256, 0, st>>>(a);
    k_rpat_pad<<<1, 1, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}

__global__ void k_cascade(SieveArgs a)
{
    u64 gid = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    u64 stride = (u64)gridDim.x * blockDim.x;
    body_cascade(gid, stride, a);
}

/* grid-stride, so the launch does not need the survivor count on the host */
static void run_cascade_gpu(SieveArgs a, cuda_stream_t st)
{
    k_cascade<<<2048, 256, 0, st>>>(a);
    CUDA_OK(cudaGetLastError());
}
static void dev_sync(void) { CUDA_OK(cudaDeviceSynchronize()); }

#endif

/* ------------------------------------------------------------------ */
/* host driver                                                        */
/* ------------------------------------------------------------------ */

static int g_nbase = 15;
static int g_s = 23;                       /* log2 of the block, in numbers */
static u64 g_sieve_limit = 1u << 19;
static int g_verify = 0;
static u64 g_maxblocks = 0;
static int g_sieve_only = 0;
static const char *g_ckpt = NULL;
static int g_presieve = 1;
static int g_ngb = 6;      /* cascade stages run on the device */
static int g_nq = 10;      /* small primes used by those stages */
static double g_t_gpu = 0, g_t_cpu = 0;

static char *u128str(u128 v, char *b)
{
    char t[42]; int i = 0;
    if (!v) { b[0] = '0'; b[1] = 0; return b; }
    while (v) { t[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    for (int j = 0; j < i; j++) b[j] = t[i - 1 - j];
    b[i] = 0; return b;
}
static char *S(u128 v) { static char r[8][42]; static int k; k = (k + 1) & 7; return u128str(v, r[k]); }

static int parse_u128(const char *s, u128 *out)
{
    u128 v = 0; const char *p = s; int frac = 0;
    if (*p < '0' || *p > '9') return -1;
    while (*p >= '0' && *p <= '9') v = v * 10 + (u128)(*p++ - '0');
    if (*p == '.') { p++; while (*p >= '0' && *p <= '9') { v = v * 10 + (u128)(*p++ - '0'); frac++; } }
    if (*p == 'e' || *p == 'E') {
        int e = atoi(p + 1) - frac;
        if (e < 0 || e > 38) return -1;
        while (e--) v *= 10;
        *out = v; return 0;
    }
    if (*p == '*') {
        u128 b = 0; p++;
        while (*p >= '0' && *p <= '9') b = b * 10 + (u128)(*p++ - '0');
        if (*p != '^') return -1;
        int e = atoi(p + 1); u128 t = 1;
        while (e-- > 0) t *= b;
        *out = v * t; return 0;
    }
    while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') p++;
    return *p || frac ? -1 : (*out = v, 0);
}

static int bitlen128(u128 x)
{
    u64 hi = (u64)(x >> 64);
    return hi ? 128 - __builtin_clzll(hi) : 64 - __builtin_clzll((u64)x);
}
static u64 hbitrev(u64 v, int n)
{
    v = ((v >> 1) & 0x5555555555555555ULL) | ((v & 0x5555555555555555ULL) << 1);
    v = ((v >> 2) & 0x3333333333333333ULL) | ((v & 0x3333333333333333ULL) << 2);
    v = ((v >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((v & 0x0F0F0F0F0F0F0F0FULL) << 4);
    return __builtin_bswap64(v) >> (64 - n);
}
static u128 hrev2(u128 x)
{
    int m = bitlen128(x);
    u64 lo = (u64)x, hi = (u64)(x >> 64);
    u128 full = ((u128)hbitrev(lo, 64) << 64) | hbitrev(hi, 64);
    return full >> (128 - m);
}

/* ---- host: build the per-block cascade constants ---- */

static u128 hrev_base(u128 x, u32 b, int *ndig)
{
    u128 r = 0; int n = 0;
    while (x) { r = r * b + (u32)(x % b); x /= b; n++; }
    if (ndig) *ndig = n;
    return r;
}
static u32 modinv_u32(u32 a, u32 q)          /* a^-1 mod prime q */
{
    u32 r = 1, e = q - 2, base = a % q;
    while (e) { if (e & 1) r = (u32)((u64)r * base % q); base = (u32)((u64)base * base % q); e >>= 1; }
    return r;
}

static void build_bctx(struct BaseCtx *c, u32 b, u128 base, int slog,
                       const u32 *qs, int nq)
{
    u64 bk = 1; u32 k = 0;
    while (bk < ((u64)1 << slog)) { bk *= b; k++; }
    c->b = b; c->k = k; c->bk = bk;
    c->Bm = (u64)(base % bk);
    u128 H = base / bk;
    for (int w = 0; w < 2; w++) {
        int nd = 0;
        u128 C = hrev_base(H + w, b, &nd);
        u128 P = 1;
        for (int i = 0; i < nd; i++) P *= b;
        for (int j = 0; j < nq; j++) {
            u32 q = qs[j];
            u32 Pq = (u32)(P % q), Cq = (u32)(C % q);
            if (Pq == 0) {
                c->mode[w][j] = Cq == 0 ? TMODE_ALWAYS : TMODE_NEVER;
                c->T[w][j] = 0;
            } else {
                c->mode[w][j] = TMODE_CMP;
                c->T[w][j] = (u32)((u64)((q - Cq) % q) * modinv_u32(Pq, q) % q);
            }
        }
    }
}

/* ---- CPU stage: cascade + confirmation over the survivor list ---- */

struct CascadeJob {
    u128 base, lo, hi;
    const u32 *out;
    u32 n;
    volatile u32 next;
    int nbase;
};

static void *cascade_worker(void *arg)
{
    CascadeJob *j = (CascadeJob *)arg;
    for (;;) {
        u32 i = __sync_fetch_and_add(&j->next, 1u);
        if (i >= j->n) break;
        u128 xv = j->base + 2 * (u128)j->out[i] + 1;
        if (xv < j->lo || xv >= j->hi) continue;
        if (a228768_check(xv)) {
            printf("\n*** HIT: %s is an emirp in all bases 2..%d ***\n", S(xv), j->nbase);
            fflush(stdout);
        }
    }
    return NULL;
}

static void run_cascade(u128 base, u128 lo, u128 hi, const u32 *out, u32 n, int nthreads)
{
    CascadeJob j;
    j.base = base; j.lo = lo; j.hi = hi;
    j.out = out; j.n = n; j.next = 0; j.nbase = g_nbase;
    if (nthreads < 1) nthreads = 1;
    if ((u32)nthreads > n) nthreads = (int)(n ? n : 1);
    pthread_t *th = (pthread_t *)malloc((size_t)nthreads * sizeof *th);
    for (int i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, cascade_worker, &j);
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    free(th);
}

/* device-or-host allocation */
static void *xalloc(size_t n)
{
#ifdef EMULATE
    void *p = calloc(1, n);
    if (!p) { fprintf(stderr, "out of memory\n"); exit(1); }
    return p;
#else
    void *p = NULL;
    CUDA_OK(cudaMallocManaged(&p, n));       /* GB10 is cache-coherent: no copies */
    memset(p, 0, n);
    return p;
#endif
}

int main(int argc, char **argv)
{
    u128 start = 0, limit = 0;
    int threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    start = (u128)30000000000ULL * 1000000000ULL;                       /* 3e19 */
    limit = (u128)32842041778ULL * 1000000000ULL + 564453125ULL;        /* 5*15^16 */

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        #define NEXT() (i + 1 < argc ? argv[++i] : (fprintf(stderr, "missing value\n"), exit(2), ""))
        if (!strcmp(a, "--nbase")) g_nbase = atoi(NEXT());
        else if (!strcmp(a, "--start")) { if (parse_u128(NEXT(), &start)) { fprintf(stderr, "bad --start\n"); return 2; } }
        else if (!strcmp(a, "--limit")) { if (parse_u128(NEXT(), &limit)) { fprintf(stderr, "bad --limit\n"); return 2; } }
        else if (!strcmp(a, "--block-log")) g_s = atoi(NEXT());
        else if (!strcmp(a, "--sieve-limit")) g_sieve_limit = strtoull(NEXT(), NULL, 0);
        else if (!strcmp(a, "--threads")) threads = atoi(NEXT());
        else if (!strcmp(a, "--gpu-threads")) g_gputhreads = (unsigned)atoi(NEXT());
        else if (!strcmp(a, "--blocks")) g_maxblocks = strtoull(NEXT(), NULL, 0);
        else if (!strcmp(a, "--verify")) g_verify = 1;
        else if (!strcmp(a, "--sieve-only")) g_sieve_only = 1;
        else if (!strcmp(a, "--gpu-bases")) g_ngb = atoi(NEXT());
        else if (!strcmp(a, "--resume")) g_ckpt = NEXT();
        else if (!strcmp(a, "--no-presieve")) g_presieve = 0;
        else if (!strcmp(a, "--gpu-primes")) g_nq = atoi(NEXT());
        else {
            fprintf(stderr,
"a228768gpu (build %s %s)\n"
"usage: %s [options]\n"
"  --nbase N          search for a(N)                      (default 15)\n"
"  --start V          lower bound; V may be 3e19 or 5*15^16\n"
"  --limit V          upper bound, exclusive\n"
"  --block-log K      GPU block = 2^K numbers              (default 28)\n"
"  --sieve-limit V    composite-filter bound               (default 524288)\n"
"  --threads T        CPU cascade threads                  (default: all cores)\n"
"  --gpu-threads T    threads per CUDA block               (default 512)\n"
"  --gpu-bases N      cascade stages run on the GPU        (default 6)\n"
"  --gpu-primes N     small primes used by those stages    (default 10)\n"
"  --resume FILE      checkpoint file; written as blocks complete, read at start\n"
"  --no-presieve      disable the small-prime pattern      (for A/B timing)\n"
"  --blocks N         stop after N blocks\n"
"  --sieve-only       skip the cascade; measures raw sieve throughput\n"
"  --verify           check one block against trial division, then continue\n",
                __DATE__, __TIME__, argv[0]);
            return 2;
        }
    }
    if (g_s < SUBLOG + 1 || g_s > 30) { fprintf(stderr, "--block-log must be %d..30\n", SUBLOG + 1); return 2; }

    if (g_ckpt) {                       /* resume where a previous run stopped */
        FILE *f = fopen(g_ckpt, "r");
        if (f) {
            char line[64];
            if (fgets(line, sizeof line, f)) {
                u128 v;
                if (!parse_u128(line, &v) && v > start && v < limit) {
                    printf("resuming at %s\n", S(v));
                    start = v;
                }
            }
            fclose(f);
        }
    }

    a228768_setup(g_nbase, (unsigned long)g_sieve_limit);
    u64 np = a228768_nprimes();
    const u32 *plist = a228768_prime_list();

    SieveArgs a;
    memset(&a, 0, sizeof a);
    a.nprimes = np;
    u32 *dp   = (u32 *)xalloc(np * sizeof(u32));
    u64 *dmag = (u64 *)xalloc(np * sizeof(u64));
    u32 *dr64 = (u32 *)xalloc(np * sizeof(u32));
    for (u64 i = 0; i < np; i++) {
        u32 p = plist[i];
        dp[i] = p;
        dmag[i] = (u64)(~(u128)0 / p >> 64);            /* floor(2^64/p) */
        dr64[i] = (u32)((~(u64)0 % p + 1) % p);         /* 2^64 mod p */
    }
    a.p = dp; a.mag = dmag; a.r64 = dr64;
    a.cinv = (u32 *)xalloc(np * sizeof(u32));

    /* pre-sieve pattern over primes 3..19, in odd-index space:
       2J+1 is divisible by p exactly when J == (p-1)/2 (mod p) */
    u64 patbits = 1;
    int npre = 0;
    if (g_presieve)
        for (u64 i = 0; i < np && plist[i] <= 19; i++) { patbits *= plist[i]; npre++; }
    else patbits = 1;
    u64 patwords = (patbits + 31) / 32 + 2;
    u32 *xpat = (u32 *)xalloc(patwords * sizeof(u32));
    if (npre) {
        for (int i = 0; i < npre; i++) {
            u64 p = plist[i];
            for (u64 j = (p - 1) / 2; j < patbits; j += p) xpat[j >> 5] |= 1u << (j & 31);
        }
        for (u32 b = 0; b < 32; b++)
            if ((xpat[b >> 5] >> (b & 31)) & 1) {
                u64 d = patbits + b;
                xpat[d >> 5] |= 1u << (d & 31);
            }
    }
    a.xpat = npre ? xpat : NULL;
    a.patbits = patbits;
    a.npre = npre;

    u64 half = (u64)1 << (g_s - 1);
    a.half = half; a.s = g_s;
    a.outcap = (u32)(half / 16);

    /* two slots, so the CPU can cascade block N's survivors while the GPU
       sieves block N+1 */
    struct Slot { u32 *rbm, *out, *outn, *out2, *outn2, *xoff, *roff, *rpat; cuda_stream_t st; } slot[2];
    for (int k = 0; k < 2; k++) {
        slot[k].rbm  = (u32 *)xalloc(half / 8);
        slot[k].out  = (u32 *)xalloc((size_t)a.outcap * sizeof(u32));
        slot[k].outn = (u32 *)xalloc(sizeof(u32));
        slot[k].out2 = (u32 *)xalloc((size_t)a.outcap * sizeof(u32));
        slot[k].outn2 = (u32 *)xalloc(sizeof(u32));
        slot[k].rpat = (u32 *)xalloc(patwords * sizeof(u32));
        slot[k].xoff = (u32 *)xalloc(np * sizeof(u32));
        slot[k].roff = (u32 *)xalloc(np * sizeof(u32));
        stream_create(&slot[k].st);
    }

    /* small primes used by the device cascade */
    static const u32 QS[MAXQ] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59};
    if (g_nq > MAXQ) g_nq = MAXQ;
    u32 *dq    = (u32 *)xalloc(MAXQ * sizeof(u32));
    u64 *dqinv = (u64 *)xalloc(MAXQ * sizeof(u64));
    u64 *dqlim = (u64 *)xalloc(MAXQ * sizeof(u64));
    for (int j = 0; j < g_nq; j++) {
        u64 q = QS[j], inv = q;
        for (int it = 0; it < 6; it++) inv *= 2 - q * inv;   /* q^-1 mod 2^64 */
        dq[j] = (u32)q; dqinv[j] = inv; dqlim[j] = ~(u64)0 / q;
    }
    a.q = dq; a.qinv = dqinv; a.qlim = dqlim; a.nq = g_nq;
    if (g_ngb > MAXGB) g_ngb = MAXGB;
    if (g_ngb > g_nbase - 2) g_ngb = g_nbase - 2;
    a.ngb = g_ngb;
    struct BaseCtx *bctx[2];
    for (int k = 0; k < 2; k++) bctx[k] = (struct BaseCtx *)xalloc(MAXGB * sizeof(struct BaseCtx));

    printf("device cascade: %d base(s) %d..%d against %d small primes\n",
           a.ngb, g_nbase, g_nbase - a.ngb + 1, a.nq);
    printf("a(%d) GPU sieve   block 2^%d   sieve primes < %llu (%llu)   %s build\n",
           g_nbase, g_s, (unsigned long long)g_sieve_limit, (unsigned long long)np,
#ifdef EMULATE
           "EMULATED (host, serial)"
#else
           "CUDA"
#endif
           );
    printf("range [%s, ", S(start));
    printf("%s)\n", S(limit));

    u128 done_to = 0;
    double last_ckpt = 0;
    int last_t = -1, cur = 0, pend = -1;
    u128 pend_base = 0, pend_lo = 0, pend_hi = 0;
    u32 pend_n = 0;
    const u32 *pend_list = NULL;
    u64 nblocks = 0, tot_cand = 0;
    double t0 = wall();

    u128 x = start;
    while (x < limit) {
        u128 run;
        if (!a228768_next_admissible(x, limit, &run)) break;
        u128 rend;
        a228768_admissible_end(run, &rend);
        if (rend > limit) rend = limit;

        u128 pos = run;
        while (pos < rend) {
            int m = bitlen128(pos);
            int s = g_s;
            if (s > m - 1) { fprintf(stderr, "block too large for %d-bit values\n", m); return 2; }
            if (m - s > 60) { fprintf(stderr, "raise --block-log: t=%d exceeds 60\n", m - s); return 2; }
            u128 base = pos & ~(((u128)1 << s) - 1);
            u128 bend = base + ((u128)1 << s);
            if (bend > rend) bend = rend;

            a.baseh = (u64)(base >> 64);
            a.basel = (u64)base;
            a.t = m - s;
            a.K = hbitrev((u64)(base >> s), a.t);

            a.bctx = bctx[cur];
            a.ngb = g_verify ? 0 : g_ngb;   /* --verify compares sieves only */
            for (int gi = 0; gi < a.ngb; gi++)
                build_bctx(&a.bctx[gi], (u32)(g_nbase - gi), base, s, QS, a.nq);

            a.rbm = slot[cur].rbm; a.out = slot[cur].out;
            a.outn = slot[cur].outn;
            a.out2 = slot[cur].out2; a.outn2 = slot[cur].outn2;
            a.xoff = slot[cur].xoff; a.roff = slot[cur].roff;
            a.rpat = slot[cur].rpat;
            a.xphase = npre ? (u64)((base >> 1) % patbits) : 0;

            if (a.t != last_t) { dev_sync(); run_cinv(a, slot[cur].st); stream_sync(slot[cur].st); last_t = a.t; }
            *a.outn = 0;
            *a.outn2 = 0;
            run_prep(a, slot[cur].st);
            run_rpat(a, slot[cur].st);
            run_rsieve(a, slot[cur].st);
            run_xsieve(a, slot[cur].st);
            if (a.ngb) run_cascade_gpu(a, slot[cur].st);

            /* GPU is busy with this block: cascade the previous one now */
            if (pend >= 0) {
                double tc = wall();
                if (!g_sieve_only)
                    run_cascade(pend_base, pend_lo, pend_hi, pend_list, pend_n, threads);
                g_t_cpu += wall() - tc;
                pend = -1;
                /* everything below pend_hi is fully processed */
                done_to = pend_hi;
                if (g_ckpt && wall() - last_ckpt > 10.0) {
                    FILE *f = fopen(g_ckpt, "w");
                    if (f) { fprintf(f, "%s\n", S(done_to)); fclose(f); }
                    last_ckpt = wall();
                }
            }
            double tg = wall();
            stream_sync(slot[cur].st);
            g_t_gpu += wall() - tg;      /* GPU time NOT hidden by the cascade */

            u32 n = a.ngb ? *a.outn2 : *a.outn;
            if (n > a.outcap || *a.outn > a.outcap) {
                fprintf(stderr, "survivor buffer overflow (%u)\n", n); return 1;
            }
            tot_cand += n;

            if (g_verify) {
                /* honest reference: a survivor is an odd x with no prime
                   factor below the sieve limit, whose base-2 reversal also
                   has none.  Compare as sets. */
                u32 *ref = (u32 *)malloc((size_t)a.outcap * sizeof(u32));
                u32 rn = 0, bad = 0;
                for (u64 j = 0; j < half; j++) {
                    u128 xv = base + 2 * (u128)j + 1;
                    int ok = 1;
                    for (u64 i = 0; i < np && ok; i++) if (xv % dp[i] == 0) ok = 0;
                    if (ok) {
                        u128 r = hrev2(xv);
                        for (u64 i = 0; i < np && ok; i++) if (r % dp[i] == 0) ok = 0;
                    }
                    if (ok && rn < a.outcap) ref[rn++] = (u32)j;
                }
                if (rn != n) bad = 1;
                else {
                    char *seen = (char *)calloc(half, 1);
                    for (u32 i = 0; i < n; i++) seen[slot[cur].out[i]] = 1;   /* ngb=0 in verify mode */
                    for (u32 i = 0; i < rn; i++) if (!seen[ref[i]]) bad = 1;
                    free(seen);
                }
                printf("verify block at %s: gpu %u, reference %u -- %s\n",
                       S(base), n, rn, bad ? "MISMATCH" : "match");
                free(ref);
                if (bad) return 1;
            }

            pend = cur; pend_base = base; pend_lo = pos; pend_hi = bend; pend_n = n;
            pend_list = a.ngb ? slot[cur].out2 : slot[cur].out;
            cur ^= 1;

            pos = bend;
            if (++nblocks % 8 == 0) {
                double el = wall() - t0;
                double scanned = (double)nblocks * (double)half;
                printf("\r%llu blocks  %.4g odds  %.1f M odds/s  survivors %.3f%%   ",
                       (unsigned long long)nblocks, scanned,
                       el > 0 ? scanned / el / 1e6 : 0.0,
                       100.0 * (double)tot_cand / scanned);
                fflush(stdout);
            }
            if (g_maxblocks && nblocks >= g_maxblocks) { x = limit; break; }
        }
        x = rend;
        if (g_maxblocks && nblocks >= g_maxblocks) break;
    }

    if (pend >= 0) {
        if (!g_sieve_only) run_cascade(pend_base, pend_lo, pend_hi, pend_list, pend_n, threads);
        done_to = pend_hi;
    }
    if (g_ckpt && done_to) {
        FILE *f = fopen(g_ckpt, "w");
        if (f) { fprintf(f, "%s\n", S(done_to)); fclose(f); }
        printf("checkpoint: complete through %s\n", S(done_to));
    }

    double el = wall() - t0;
    double scanned = (double)nblocks * (double)half;
    printf("\n%llu blocks, %.4g odd numbers, survivors %.4f%%, %.1f s\n",
           (unsigned long long)nblocks, scanned,
           scanned > 0 ? 100.0 * (double)tot_cand / scanned : 0.0, el);
    printf("%.1f M odds/s%s\n", el > 0 ? scanned / el / 1e6 : 0.0,
           g_sieve_only ? "   (SIEVE ONLY -- no cascade, results not checked)" : "");
    /* The two stages overlap, so wall = max(GPU, CPU), never the sum.
       g_t_cpu is cascade time that ran while the GPU was busy; g_t_gpu is
       the GPU time left over after that overlap. */
    printf("time: wall %.1f s;  CPU cascade %.1f s (overlapped);  residual GPU wait %.1f s\n",
           el, g_t_cpu, g_t_gpu);
    if (g_t_gpu < 0.05 * el)
        printf("  --> CPU-BOUND: device work hides entirely inside the cascade.\n"
               "      Raise --gpu-bases to move more of it onto the GPU.\n");
    else if (g_t_cpu < 0.05 * el)
        printf("  --> GPU-BOUND: device busy ~%.1f s, cascade almost free.\n"
               "      Lower --gpu-bases; the CPU has idle capacity.\n", el);
    else
        printf("  --> GPU-BOUND: device busy ~%.1f s, with %.1f s of cascade hidden inside it.\n"
               "      Near the crossover; try --gpu-bases one lower.\n", el, g_t_cpu);
    return 0;
}
