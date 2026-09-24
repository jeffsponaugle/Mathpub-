/*
 * a158939_cuda.cu  --  GPU (CUDA) version of the A158939 scan (see ../a158939.c)
 *
 * Same problem and output as the CPU tool: for every prime p find
 * L(p) = largest n with g_1 < g_2 < ... < g_n (the gaps after p strictly
 * increasing), and report the first prime of every run length, the exact
 * prime count and the exact run-length histogram of the range scanned.
 *
 * Sieve layout
 * ------------
 * Only numbers coprime to 210 are represented: 48 bits per period of 210
 * numbers, bit j of period P standing for 210*P + R[j] with R the 48 residues
 * coprime to 210.  A bit that is still 0 after sieving is a prime.  For a
 * sieving prime q >= 11 the multiples q*m with gcd(m,210) = 1 fall into 48
 * arithmetic progressions of period index (stride q) with a fixed residue
 * index each, so a mark is "atomicOr(word[i >> 5], 1 << (i & 31))" with the
 * bit index i advancing by 48*q.  The primes 2, 3, 5, 7 are added by the host.
 *
 * The range is processed in chunks of CHUNK_SEGS segments of SEG_WORDS words
 * (3.01e6 numbers per segment, 2.89e8 per chunk).  Per chunk:
 *
 *   1. kernel_large / kernel_bucket:  the sieving primes q > QSPLIT (up to
 *      sqrt of the chunk end) mark their few multiples in the chunk into an
 *      8 MB bitmap in global memory with atomics; 8 MB stays inside the
 *      GB10's 24 MB L2 (pinned with a persisting access window), where the
 *      atomics are cheap.  Each prime keeps a 4-byte state: chunk (mod 64)
 *      and byte offset of its next multiple, wheel index.  Primes below
 *      CHUNK_NUMS/6 hit every chunk and live in a flat array (one thread
 *      each); larger ones live as 8-byte entries in 64 bucket lists, one per
 *      upcoming chunk, so a chunk only touches the primes that hit it and
 *      re-appends each to the list of its next chunk (block-staged so the
 *      appends are coalesced).
 *   2. kernel_segment: one block per segment, the segment's words in shared
 *      memory.  It ORs in a periodic pattern for the primes 7..19 and the
 *      large-prime bits from the chunk bitmap, marks the primes 23..QSPLIT
 *      (warp per prime below PWARP, thread per prime above), then walks the
 *      finished bitmap: each thread takes a slice of words, extracts primes
 *      and gaps, and runs the same increasing-run bookkeeping as the CPU tool
 *      (a run that ends at q_k settles the run lengths of the L primes before
 *      it; a thread continues past its slice until the run in progress ends).
 *      Every block also sieves OVERLAP words past its own segment so that the
 *      last slices can settle their runs; a run still open at the end of the
 *      overlap is flagged for the host to settle on the CPU (never happens
 *      in practice).  Per-segment results (first prime of each run length,
 *      histogram, prime count) go to global memory; the host merges them.
 *
 * 2, 3, 5 are not in the wheel; the host adds them (L(2)=2, L(3)=1, L(5)=2)
 * when the scan starts at 0.  Chunks complete in order, so every term found
 * is confirmed immediately.  The checkpoint (-S) records the chunk frontier,
 * the table, the histogram and the prime count; the large-prime state is
 * recomputed on resume.
 *
 * Usage
 *   a158939_cuda scan [START] END [-n N] [-r NMIN] [-S FILE] [-i SECS] [-Q QSPLIT] [-q]
 *   a158939_cuda bench [N] [-d SPAN] [-Q QSPLIT]
 *   a158939_cuda selftest
 *
 * Build (DGX Spark / GB10, CUDA 13):
 *   nvcc -O3 -std=c++17 -arch=sm_121 a158939_cuda.cu -o a158939_cuda
 */

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cstdarg>
#include <cctype>
#include <cmath>
#include <csignal>
#include <cerrno>
#include <ctime>
#include <cinttypes>
#include <climits>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <cuda_runtime.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t  u8;
typedef unsigned __int128 u128;

/* ---- geometry (compile-time, overridable with -D) ---- */
#ifndef SEG_WORDS
#define SEG_WORDS   21504           /* words owned by a segment: 84 KB = 14336 periods of 210 = 3,010,560 numbers */
#endif
#ifndef OVERLAP
#define OVERLAP     528             /* extra words sieved past a segment: 352 periods = 73,920 numbers */
#endif
#ifndef CHUNK_SEGS
#define CHUNK_SEGS  96              /* segments per chunk: 7.875 MB bitmap (must stay in L2) */
#endif
#ifndef BLOCKDIM
#define BLOCKDIM    1024
#endif
#ifndef PWARP
#define PWARP       2048            /* primes below this: one warp per prime */
#endif
#define WHEEL       210
#define NRES        48              /* residues coprime to 210 */
#define SHARED_WORDS (SEG_WORDS + OVERLAP)
#define SEG_PERIODS (SEG_WORDS * 32 / NRES)              /* 14336 */
#define OVL_PERIODS (OVERLAP * 32 / NRES)                /* 352 */
#define SEG_NUMS    ((u64)SEG_PERIODS * WHEEL)           /* 3,010,560 */
#define CHUNK_WORDS ((u64)SEG_WORDS * CHUNK_SEGS)
#define CHUNK_NUMS  ((u64)SEG_PERIODS * CHUNK_SEGS * WHEEL)  /* 289,013,760 */
#define PAT_PERIOD  46189           /* 11*13*17*19: period of the presieve pattern in wheel periods */
#define PAT_BYTES   (PAT_PERIOD * 6)
#define NMAX        64              /* run lengths 1..NMAX-1 tracked */
#define NKNOWN      15

static const u64 KNOWN[NKNOWN + 1] = {
    0, 3, 2, 17, 347, 2903, 15373, 128981, 1319407, 17797517, 94097537,
    6927837557ULL, 48486712783ULL, 968068681511ULL, 1472840004017ULL,
    129001208165717ULL
};
/* A133697(n) = pi(A158939(n+2)) for n = 0..13 */
static const u64 A133697[14] = {
    1, 7, 69, 420, 1796, 12073, 101397, 1139211, 5440508, 320620306ULL,
    2058187481ULL, 36451609409ULL, 54594153615ULL, 4100904808215ULL
};
/* pi(10^k), k = 6..16 */
static const u64 PI10[17] = {
    0, 0, 0, 0, 0, 0, 78498ULL, 664579ULL, 5761455ULL, 50847534ULL, 455052511ULL,
    4118054813ULL, 37607912018ULL, 346065536839ULL, 3204941750802ULL,
    29844570422669ULL, 279238341033925ULL
};

static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;
static int g_report = 14;
static bool g_l2_persist = true;       /* keep the bitmap resident in L2 with a persisting access window (-N disables) */
static bool g_quiet_init = false;

#define CUDA_CHECK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) die("CUDA error %s at line %d: %s", #x, __LINE__, cudaGetErrorString(e_)); } while (0)

/* ------------------------------------------------------------------ */
/* Generic helpers (same conventions as a158939.c)                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void die(const char *fmt, ...)
{
    va_list ap;
    if (stderr_tty) fputs("\r\033[K", stderr);
    fputs("a158939_cuda: ", stderr);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static const char *fmt_dur(double s, char *buf, size_t n)
{
    if (s < 0 || s != s || s > 1e9) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

static bool mul_ok(u64 a, u64 b, u64 *r)
{
    u128 x = (u128)a * b;
    if (x >> 64) return false;
    *r = (u64)x;
    return true;
}

static bool all_digits(const char *s, size_t n)
{
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

static bool parse_term(const char *s, size_t n, u64 *out)
{
    char buf[64], *c;
    if (n == 0 || n >= sizeof buf) return false;
    memcpy(buf, s, n);
    buf[n] = 0;
    if ((c = strchr(buf, '^'))) {
        *c = 0;
        if (!all_digits(buf, strlen(buf)) || !all_digits(c + 1, strlen(c + 1))) return false;
        u64 b = strtoull(buf, 0, 10), e = strtoull(c + 1, 0, 10), r = 1;
        while (e--) if (!mul_ok(r, b, &r)) return false;
        *out = r;
        return true;
    }
    if ((c = strpbrk(buf, "eE"))) {
        *c = 0;
        if (!all_digits(c + 1, strlen(c + 1))) return false;
        long e = strtol(c + 1, 0, 10);
        char mant[64];
        int md = 0, frac = -1;
        for (const char *q = buf; *q; q++) {
            if (*q == '.') { if (frac >= 0) return false; frac = 0; continue; }
            if (!isdigit((unsigned char)*q)) return false;
            mant[md++] = *q;
            if (frac >= 0) frac++;
        }
        if (md == 0 || md > 19) return false;
        mant[md] = 0;
        if (frac > 0) e -= frac;
        if (e < 0) return false;
        u64 m = strtoull(mant, 0, 10);
        while (e-- > 0) if (!mul_ok(m, 10, &m)) return false;
        *out = m;
        return true;
    }
    if (!all_digits(buf, n) || n > 20) return false;
    errno = 0;
    u64 v = strtoull(buf, 0, 10);
    if (errno) return false;
    *out = v;
    return true;
}

static u64 parse_num(const char *s)
{
    size_t n = strlen(s);
    for (size_t i = n; i-- > 1; ) {
        if (s[i] != '+' && s[i] != '-') continue;
        u64 a, b;
        if (!parse_term(s, i, &a) || !parse_term(s + i + 1, n - i - 1, &b)) die("bad number '%s'", s);
        if (s[i] == '+') { if (a + b < a) die("overflow in '%s'", s); return a + b; }
        if (b > a) die("negative result in '%s'", s);
        return a - b;
    }
    u64 v;
    if (!parse_term(s, n, &v)) die("bad number '%s'", s);
    return v;
}

static int parse_int(const char *s, int lo, int hi, const char *what)
{
    u64 v = parse_num(s);
    if (v < (u64)lo || v > (u64)hi) die("%s must be between %d and %d", what, lo, hi);
    return (int)v;
}

/* deterministic Miller-Rabin below 2^64, for the host-side odds and ends */
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }
static u64 powmod(u64 a, u64 e, u64 m)
{
    u64 r = 1;
    a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static bool is_prime64(u64 n)
{
    static const u64 small[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    static const u64 bases[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    if (n < 2) return false;
    for (size_t i = 0; i < sizeof small / sizeof *small; i++)
        if (n % small[i] == 0) return n == small[i];
    u64 d = n - 1;
    int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    for (size_t i = 0; i < sizeof bases / sizeof *bases; i++) {
        u64 a = bases[i] % n;
        if (a == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int j = 1; j < r && composite; j++) { x = mulmod(x, x, n); if (x == n - 1) composite = false; }
        if (composite) return false;
    }
    return true;
}
static u64 next_prime64(u64 n)
{
    if (n < 2) return 2;
    u64 m = n + 1;
    if (!(m & 1)) m++;
    while (!is_prime64(m)) m += 2;
    return m;
}
/* run length of p by next-prime search; gaps[0..L] filled (L+1 values) */
static int run_length64(u64 p, u16 *gaps, int maxgaps)
{
    u64 prev = p, prevgap = 0;
    int L = 0;
    for (;;) {
        u64 q = next_prime64(prev), g = q - prev;
        if (L >= maxgaps) return -1;
        gaps[L] = (u16)g;
        if (g <= prevgap) return L;
        L++;
        prevgap = g;
        prev = q;
    }
}

/* ------------------------------------------------------------------ */
/* Device side                                                         */
/* ------------------------------------------------------------------ */

__constant__ u8  c_R[NRES]    = { 1,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,
                                  107,109,113,121,127,131,137,139,143,149,151,157,163,167,169,173,179,181,187,191,193,197,199,209 };
__constant__ u8  c_GAP[NRES];       /* R[j+1] - R[j] (2 for the wrap 209 -> 211) */
__constant__ u8  c_RIDX[WHEEL];     /* residue -> index, 255 if not coprime to 210 */
__constant__ u8  c_INV[NRES];       /* inverse of R[j] modulo 210 */
__constant__ u8  c_T[NRES][NRES];   /* floor(R[b]*R[j] / 210) */
__constant__ u8  c_J[NRES][NRES];   /* index of (R[b]*R[j]) mod 210 */

/* Per-segment results written by kernel_segment, merged on the host. */
struct SegResult {
    u32 first[NMAX];        /* offset (from seg_lo) of the first prime with run length n, ~0u = none */
    u32 firstidx[NMAX];     /* its index among the segment's primes inside the bounds (1-based) */
    u32 runs[NMAX];         /* run ends of length L whose L primes are all inside the slice bounds */
    u32 extra[NMAX];        /* run lengths from final run ends, per prime inside the bounds */
    u32 count;              /* primes in [seg_lo, seg_hi) ∩ [start, end) */
    u32 unsettled;          /* offset of the first prime of a run still open at the bitmap end, ~0u = none */
    u32 overflow;           /* a run reached NMAX-1 gaps */
    u32 pad;
};

/* First multiplier m >= max(q, ceil(lo/q)) with gcd(m, 210) = 1. */
__device__ __forceinline__ u64 first_multiplier(u32 q, u64 lo)
{
    u64 m = (lo + q - 1) / q;
    if (m < q) m = q;
    u32 r = (u32)(m % WHEEL);
    while (c_RIDX[r] == 255) { m++; r = (r == WHEEL - 1) ? 0 : r + 1; }
    return m;
}

/* bit index of number n inside the chunk starting at base (n coprime to 210, n >= base, base multiple of 210) */
__device__ __forceinline__ u32 bit_index(u64 n, u64 base)
{
    u32 d = (u32)(n - base);
    return NRES * (d / WHEEL) + c_RIDX[d % WHEEL];
}

/* Large-prime state, one u32 per prime: bits 31..26 = index (mod 64) of the
 * chunk holding the prime's next multiple, bits 25..0 = that multiple's bit
 * index inside the chunk (< 2^26).  The multiplier's residue class is
 * recovered from the bit and q.  Consecutive multiples are at most 10q
 * apart, so the chunk distance is below 64 for every q below 1.8e9.
 *
 * Primes below CHUNK_NUMS/10 ("dense") hit every chunk and keep their state
 * in a flat array read and rewritten by kernel_large.  Larger primes
 * ("sleepers") live as 8-byte entries (state << 32 | q) in a ring of SLOTS
 * bucket lists, one per upcoming chunk, so a chunk touches only the primes
 * that hit it (kernel_bucket); each entry is re-appended to the list of the
 * chunk holding its next multiple, block-staged so the appends are coalesced. */
#define SLOTS        64
#define OVERLAP_NUMS ((u64)OVL_PERIODS * WHEEL)
#define BLOCKDIM_B   1024

__device__ __forceinline__ u32 pack_state(u64 n)
{
    u64 kc = n / CHUNK_NUMS;
    return ((u32)(kc & (SLOTS - 1)) << 26) | bit_index(n, kc * CHUNK_NUMS);
}

/* the chunk whose extended range [chunk_lo, hi_ext) first contains n */
__device__ __forceinline__ u64 slot_of(u64 n)
{
    return n < OVERLAP_NUMS ? 0 : (n - OVERLAP_NUMS) / CHUNK_NUMS;
}

__device__ __forceinline__ void append_entry(u64 *entries, u32 *counts, u32 cap, u32 *err, u64 slot, u64 entry)
{
    u32 s = (u32)(slot % SLOTS);
    u32 pos = atomicAdd(&counts[s], 1u);
    if (pos < cap) entries[(u64)s * cap + pos] = entry;
    else atomicOr(err, 1u);
}

/* State of the active large primes [0, nactive) for a scan (re)starting at chunk k0. */
__global__ void kernel_init_state(const u32 *__restrict__ primes, u32 *__restrict__ state, u32 nd, u32 nactive, u64 k0,
                                  u64 *__restrict__ entries, u32 *__restrict__ counts, u32 cap, u32 *__restrict__ err)
{
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nactive) return;
    const u64 chunk_lo = k0 * CHUNK_NUMS;
    u32 q = primes[i];
    u64 n = (u64)q * first_multiplier(q, chunk_lo);
    u32 st = pack_state(n);
    if (i < nd) { state[i] = st; return; }
    u64 sl = slot_of(n);
    if (sl < k0) sl = k0;                             /* the first chunk's overlap belongs to no earlier chunk here */
    append_entry(entries, counts, cap, err, sl, ((u64)st << 32) | q);
}

/* Newly active primes [nprev, nactive): their first multiple q*q is in the current chunk or its overlap. */
__global__ void kernel_activate(const u32 *__restrict__ primes, u32 *__restrict__ state, u32 nd, u32 nprev, u32 nactive, u64 k,
                                u64 *__restrict__ entries, u32 *__restrict__ counts, u32 cap, u32 *__restrict__ err)
{
    u32 i = nprev + blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nactive) return;
    u32 q = primes[i];
    u64 n = (u64)q * q;
    u32 st = pack_state(n);
    if (i < nd) { state[i] = st; return; }
    u64 sl = slot_of(n);
    if (sl < k) sl = k;
    append_entry(entries, counts, cap, err, sl, ((u64)st << 32) | q);
}

/* (chunk_lo/210) mod q for every medium sieving prime (one 64-bit division each, once per chunk) */
__global__ void kernel_chunkmod(const u32 *__restrict__ sprimes, u32 *__restrict__ rchunk, u32 nsmall, u64 chunk_lo)
{
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nsmall) rchunk[i] = (u32)((chunk_lo / WHEEL) % sprimes[i]);
}

/* x mod q for x < 2^32, with inv = floor(2^32 / q) */
__device__ __forceinline__ u32 barrett32(u32 x, u32 q, u32 inv)
{
    u32 r = x - __umulhi(x, inv) * q;
    return r >= q ? r - q : r;
}

/* Mark the multiples of q from n (multiplier residue index jm) below hi_ext into the
 * chunk bitmap; returns the packed state of the first multiple >= chunk_hi. */
__device__ __forceinline__ u32 mark_prime(u32 q, u64 n, u32 jm, u64 chunk_lo, u64 chunk_hi, u64 hi_ext,
                                          u32 *__restrict__ bitmap, const u8 *s_gap, const u8 *s_ridx)
{
    u64 nsave = 0;
    bool saved = false;
    do {
        if (!saved && n >= chunk_hi) { nsave = n; saved = true; }
        u32 d = (u32)(n - chunk_lo);
        u32 i = NRES * (d / WHEEL) + s_ridx[d % WHEEL];
        atomicOr(&bitmap[i >> 5], 1u << (i & 31));
        n += (u64)q * s_gap[jm];
        jm = (jm == NRES - 1) ? 0 : jm + 1;
    } while (n < hi_ext);
    if (!saved) nsave = n;
    return pack_state(nsave);
}

/* Multiple encoded by state s for chunk k (in the chunk, or in its overlap when the target is k+1),
 * and the residue index of its multiplier. */
__device__ __forceinline__ bool state_position(u32 s, u32 q, u64 k, u64 chunk_lo, u64 chunk_hi, u64 *n, u32 *jm,
                                               const u8 *s_ridx, const u8 *s_inv, const u8 *s_R)
{
    const u32 tgt = s >> 26, kk = (u32)(k & (SLOTS - 1)), i = s & 0x3FFFFFFu;
    u64 base;
    if (tgt == kk) base = chunk_lo;
    else if (tgt == ((kk + 1) & (SLOTS - 1)) && i < OVERLAP * 32) base = chunk_hi;
    else return false;
    const u32 rn = s_R[i % NRES];
    *n = base + (u64)WHEEL * (i / NRES) + rn;
    const u32 qi = s_inv[s_ridx[q % WHEEL]];                    /* q^-1 mod 210 */
    *jm = s_ridx[(rn * qi) % WHEEL];                            /* m = n/q  =>  m mod 210 = rn * q^-1 */
    return true;
}

/* Dense large primes: one thread per prime, state always rewritten (coalesced). */
__global__ void __launch_bounds__(256)
kernel_large(const u32 *__restrict__ primes, u32 *__restrict__ state, u32 n_dense,
             u64 k, u64 chunk_lo, u64 chunk_hi, u64 hi_ext, u32 *__restrict__ bitmap)
{
    __shared__ u8 s_ridx[WHEEL], s_gap[NRES], s_inv[NRES], s_R[NRES];
    for (u32 t = threadIdx.x; t < WHEEL; t += blockDim.x) s_ridx[t] = c_RIDX[t];
    if (threadIdx.x < NRES) { s_gap[threadIdx.x] = c_GAP[threadIdx.x]; s_inv[threadIdx.x] = c_INV[threadIdx.x]; s_R[threadIdx.x] = c_R[threadIdx.x]; }
    __syncthreads();
    u32 i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_dense) return;
    const u32 s = __ldcs(&state[i]);
    const u32 q = __ldcs(&primes[i]);
    u64 n;
    u32 jm;
    if (!state_position(s, q, k, chunk_lo, chunk_hi, &n, &jm, s_ridx, s_inv, s_R)) { __stcs(&state[i], s); return; }
    __stcs(&state[i], mark_prime(q, n, jm, chunk_lo, chunk_hi, hi_ext, bitmap, s_gap, s_ridx));
}

/* Sleeping large primes due in chunk k: process the n_in entries of its bucket list,
 * re-append each to the list of the chunk holding its next multiple. */
__global__ void __launch_bounds__(BLOCKDIM_B)
kernel_bucket(u64 *__restrict__ entries, u32 *__restrict__ counts, u32 cap, u32 *__restrict__ err, u32 n_in,
              u64 k, u64 chunk_lo, u64 chunk_hi, u64 hi_ext, u32 *__restrict__ bitmap)
{
    __shared__ u32 s_cnt[SLOTS], s_off[SLOTS], s_base[SLOTS];
    __shared__ u64 s_stage[BLOCKDIM_B];
    __shared__ u8 s_rel[BLOCKDIM_B];
    __shared__ u8 s_ridx[WHEEL], s_gap[NRES], s_inv[NRES], s_R[NRES];
    const u32 tid = threadIdx.x;
    if (tid < SLOTS) s_cnt[tid] = 0;
    if (tid < WHEEL) s_ridx[tid] = c_RIDX[tid];
    if (tid < NRES) { s_gap[tid] = c_GAP[tid]; s_inv[tid] = c_INV[tid]; s_R[tid] = c_R[tid]; }
    __syncthreads();

    const u32 slot = (u32)(k % SLOTS);
    const u32 i = blockIdx.x * BLOCKDIM_B + tid;
    const bool have = i < n_in;
    u32 rel = 0, rank = 0;
    u64 out = 0;
    if (have) {
        const u64 e = entries[(u64)slot * cap + i];
        const u32 q = (u32)e, s = (u32)(e >> 32);
        u64 n;
        u32 jm;
        if (state_position(s, q, k, chunk_lo, chunk_hi, &n, &jm, s_ridx, s_inv, s_R)) {
            const u32 st = mark_prime(q, n, jm, chunk_lo, chunk_hi, hi_ext, bitmap, s_gap, s_ridx);
            const u32 tgt = st >> 26, bi = st & 0x3FFFFFFu;
            u64 kc = k + (((tgt - (u32)(k % SLOTS)) & (SLOTS - 1)));        /* absolute chunk of the next multiple */
            u64 nn = kc * CHUNK_NUMS + (u64)WHEEL * (bi / NRES);               /* (residue irrelevant for the slot) */
            u64 sl = slot_of(nn);
            if (sl <= k) sl = k + 1;                                            /* in this chunk's overlap: again next chunk */
            rel = (u32)(sl - k);
            out = ((u64)st << 32) | q;
            rank = atomicAdd(&s_cnt[rel], 1u);
        } else {
            atomicOr(err, 2u);                                                  /* entry in the wrong list */
        }
    }
    __syncthreads();
    if (tid == 0) {
        u32 acc = 0;
        for (u32 r = 0; r < SLOTS; r++) { s_off[r] = acc; acc += s_cnt[r]; }
    }
    if (tid < SLOTS) {
        u32 c = s_cnt[tid];
        u32 base = c ? atomicAdd(&counts[(u32)((k + tid) % SLOTS)], c) : 0;
        s_base[tid] = base;
        if (c && base + c > cap) atomicOr(err, 1u);
    }
    __syncthreads();
    if (have && out) {
        u32 p = s_off[rel] + rank;
        s_stage[p] = out;
        s_rel[p] = (u8)rel;
    }
    __syncthreads();
    const u32 total = s_off[SLOTS - 1] + s_cnt[SLOTS - 1];
    if (tid < total) {
        const u32 r = s_rel[tid];
        const u32 pos = s_base[r] + (tid - s_off[r]);
        if (pos < cap) entries[(u64)((k + r) % SLOTS) * cap + pos] = s_stage[tid];
    }
}

/* Block-wide exclusive prefix sum (blockDim.x <= 1024, multiple of 32); scratch: 32 words.
 * Returns the exclusive prefix; *tot receives the block total (all threads). */
__device__ __forceinline__ u32 block_exscan(u32 v, u32 *scratch, u32 *tot)
{
    const int lane = threadIdx.x & 31, warp = threadIdx.x >> 5, nwarps = blockDim.x >> 5;
    u32 x = v;
    for (int o = 1; o < 32; o <<= 1) { u32 y = __shfl_up_sync(0xffffffffu, x, o); if (lane >= o) x += y; }
    if (lane == 31) scratch[warp] = x;
    __syncthreads();
    if (warp == 0) {
        u32 w = (lane < nwarps) ? scratch[lane] : 0;
        for (int o = 1; o < 32; o <<= 1) { u32 y = __shfl_up_sync(0xffffffffu, w, o); if (lane >= o) w += y; }
        if (lane < nwarps) scratch[lane] = w;
    }
    __syncthreads();
    u32 prefix = (warp ? scratch[warp - 1] : 0) + x - v;
    *tot = scratch[nwarps - 1];
    __syncthreads();
    return prefix;
}

/* number offset (from the segment start) of bit i of the segment bitmap */
__device__ __forceinline__ u32 numoff(u32 i, const u8 *s_R) { return WHEEL * (i / NRES) + s_R[i % NRES]; }

/* Record the first `count` primes of the run that starts at offset o (a prime)
 * as the first occurrences of run lengths run, run-1, ...: prime t of the walk
 * has run length run - t and index idx0 + t.  Stops at slice_hi. */
__device__ __forceinline__ void record_from(const u32 *sh, u32 o, u32 idx0, u32 count, u32 run, u32 slice_hi,
                                            unsigned long long *s_first, const u8 *s_R, const u8 *s_ridx)
{
    u32 ib = NRES * (o / WHEEL) + s_ridx[o % WHEEL];
    u32 w = ib >> 5, bi = ib & 31;
    u32 x = (~sh[w]) & (~0u << bi);
    for (u32 t = 0; t < count; t++) {
        while (!x) { if (++w >= SHARED_WORDS) return; x = ~sh[w]; }
        u32 bit = __ffs(x) - 1;
        x &= x - 1;
        u32 off = numoff(32 * w + bit, s_R);
        if (off >= slice_hi) return;
        atomicMin(&s_first[run - t], ((unsigned long long)off << 32) | (idx0 + t));
    }
}

/*
 * One block per segment.  sprimes[0..nsmall) are the sieving primes 23..QSPLIT
 * (ascending, nwarpprimes of them below PWARP).  start/end bound the primes to
 * be accounted for.
 */
__global__ void __launch_bounds__(BLOCKDIM)
kernel_segment(const u32 *__restrict__ bitmap, const u8 *__restrict__ pattern,
               const u32 *__restrict__ sprimes, const u32 *__restrict__ sinv, const u32 *__restrict__ ssegmod,
               const u32 *__restrict__ rchunk, u32 nsmall, u32 nwarpprimes,
               u64 chunk_lo, u64 start, u64 end, SegResult *__restrict__ results)
{
    extern __shared__ u32 sh[];                         /* SHARED_WORDS words of bitmap */
    __shared__ unsigned long long s_first[NMAX];       /* (offset << 32) | prime index */
    __shared__ u32 s_runs[NMAX], s_extra[NMAX];
    __shared__ u32 s_unsettled, s_overflow, s_scratch[32];
    __shared__ u8 s_R[NRES], s_ridx[WHEEL];
    __shared__ u8 s_T[NRES][NRES], s_J[NRES][NRES];

    const u32 seg = blockIdx.x;
    const u64 seg_lo = chunk_lo + (u64)seg * SEG_NUMS;
    const u32 word0 = seg * SEG_WORDS;
    const u32 tid = threadIdx.x;

    /* --- 0. tables --- */
    for (u32 t = tid; t < NRES * NRES; t += BLOCKDIM) { (&s_T[0][0])[t] = (&c_T[0][0])[t]; (&s_J[0][0])[t] = (&c_J[0][0])[t]; }
    for (u32 t = tid; t < WHEEL; t += BLOCKDIM) s_ridx[t] = c_RIDX[t];
    if (tid < NRES) s_R[tid] = c_R[tid];

    /* --- 1. presieve pattern (11, 13, 17, 19) and large-prime bits --- */
    {
        u32 pbase = (u32)(((seg_lo / WHEEL) % PAT_PERIOD) * 6);   /* byte offset of this segment's first period */
        for (u32 w = tid; w < SHARED_WORDS; w += BLOCKDIM) {
            u32 kk = pbase + 4 * w;
            if (kk >= PAT_BYTES) kk -= PAT_BYTES;
            u32 v = 0;
            #pragma unroll
            for (int t = 0; t < 4; t++) {
                v |= (u32)pattern[kk] << (8 * t);
                kk++;
                if (kk == PAT_BYTES) kk = 0;
            }
            sh[w] = v | __ldcs(&bitmap[word0 + w]);
        }
        if (tid < NMAX) { s_first[tid] = ~0ull; s_runs[tid] = 0; s_extra[tid] = 0; }
        if (tid == 0) { s_unsettled = ~0u; s_overflow = 0; }
    }
    __syncthreads();

    /* --- 2. small primes: warp per prime (q < PWARP), then thread per prime ---
     * For q = 210a + r (r coprime to 210, class b) and the multiplier class j,
     * the multiples q*(210t + R[j]) have period index a*R[j] + T[b][j] + q*t
     * and residue index J[b][j].  With P0 = seg_lo/210 (mod q, from rchunk and
     * ssegmod by Barrett), the first period offset in the segment is
     * (a*R[j] + T[b][j] - P0) mod q, the bit index 48*that + J[b][j], and
     * the stride 48*q bits.  Only q itself (multiplier 1) must not be marked. */
    const u32 SHBITS = SHARED_WORDS * 32;
    {
        const u32 warp = tid >> 5, lane = tid & 31, nwarps = BLOCKDIM / 32;
        for (u32 i = warp; i < nwarpprimes; i += nwarps) {
            const u32 q = sprimes[i], inv = sinv[i];
            const u32 a = q / WHEEL, b = s_ridx[q % WHEEL];
            const u32 p0 = barrett32(rchunk[i] + seg * ssegmod[i], q, inv);   /* P0 mod q */
            const bool self = (u64)q >= seg_lo && (u64)q < seg_lo + SEG_NUMS + OVERLAP_NUMS;
            for (u32 j = lane; j < NRES; j += 32) {
                u32 c = a * s_R[j] + s_T[b][j];
                if (c >= q) c -= q;
                u32 d = c >= p0 ? c - p0 : c - p0 + q;
                u32 bi = NRES * d + s_J[b][j];
                if (self && j == 0) bi += NRES * q;                          /* skip the multiple q*1 = q */
                for (; bi < SHBITS; bi += NRES * q)
                    atomicOr(&sh[bi >> 5], 1u << (bi & 31));
            }
        }
        for (u32 i = nwarpprimes + tid; i < nsmall; i += BLOCKDIM) {
            const u32 q = sprimes[i], inv = sinv[i];
            const u32 a = q / WHEEL, b = s_ridx[q % WHEEL];
            const u32 p0 = barrett32(rchunk[i] + seg * ssegmod[i], q, inv);
            const bool self = (u64)q >= seg_lo && (u64)q < seg_lo + SEG_NUMS + OVERLAP_NUMS;
            for (u32 j = 0; j < NRES; j++) {
                u32 c = a * s_R[j] + s_T[b][j];
                if (c >= q) c -= q;
                u32 d = c >= p0 ? c - p0 : c - p0 + q;
                u32 bi = NRES * d + s_J[b][j];
                if (self && j == 0) bi += NRES * q;
                for (; bi < SHBITS; bi += NRES * q)
                    atomicOr(&sh[bi >> 5], 1u << (bi & 31));
            }
        }
    }
    __syncthreads();
    if (seg_lo == 0 && tid == 0) sh[0] = (sh[0] & ~0x1Eu) | 1u;   /* 1 is not prime; 11, 13, 17, 19 are */
    __syncthreads();

    /* --- 3. prime count in the owned region ∩ [start, end) --- */
    u32 lo_off = (start > seg_lo) ? (u32)(start - seg_lo) : 0;              /* first number accounted */
    u32 hi_off = (end < seg_lo + SEG_NUMS) ? (u32)(end - seg_lo) : (u32)SEG_NUMS;   /* exclusive */
    if (lo_off > hi_off) lo_off = hi_off;
    const u32 WPT = SEG_WORDS / BLOCKDIM;                  /* words per thread slice */
    const u32 w0 = tid * WPT;
    u32 prefix;
    {
        u32 c = 0;
        for (u32 w = w0; w < w0 + WPT; w++) {
            u32 x = ~sh[w];
            if (x) {
                u32 wlo = numoff(32 * w, s_R), whi = numoff(32 * w + 32, s_R);   /* word w covers offsets [wlo, whi) */
                if (wlo >= lo_off && whi <= hi_off) c += __popc(x);
                else if (whi > lo_off && wlo < hi_off) {
                    while (x) {
                        u32 bit = __ffs(x) - 1;
                        x &= x - 1;
                        u32 off = numoff(32 * w + bit, s_R);
                        c += (off >= lo_off && off < hi_off);
                    }
                }
            }
        }
        u32 total;
        prefix = block_exscan(c, s_scratch, &total);
        if (tid == 0) results[seg].count = total;
    }

    /* --- 4. run lengths: thread per slice of words, continuing until settled --- */
    u32 c1 = 0, c2 = 0, c3 = 0, c4 = 0, c5 = 0, c6 = 0, c7 = 0;
    {
        u32 slice_lo = numoff(32 * w0, s_R), slice_hi = numoff(32 * (w0 + WPT), s_R);   /* offsets */
        if (slice_lo < lo_off) slice_lo = lo_off;
        if (slice_hi > hi_off) slice_hi = hi_off;
#ifdef SKIP_EXTRACT
        slice_hi = slice_lo;
#endif
        if (slice_lo < slice_hi) {
            u32 run = 0, M = 0, prev = 0, prevgap = 0, cnt = 0, run_out = 0, run_start = 0, run_start_cnt = 0;
            bool have_prev = false, done = false;
            for (u32 w = w0; w < SHARED_WORDS && !done; w++) {
                u32 x = ~sh[w];
                while (x) {
                    const u32 bit = __ffs(x) - 1;
                    x &= x - 1;
                    const u32 off = numoff(32 * w + bit, s_R);
                    if (!have_prev) {
                        if (off < slice_lo) continue;
                        if (off >= slice_hi) { done = true; break; }   /* no prime in the slice */
                        have_prev = true;
                        prev = off;
                        cnt = 1;
                        run_start = off;
                        run_start_cnt = 1;
                        continue;
                    }
                    const u32 g = off - prev;
                    if (g > prevgap) {
                        run++;                                   /* prev joins the run's primes */
                        run_out += (prev >= slice_hi);
                    } else {
                        if (run > NMAX - 1) { s_overflow = 1; run = NMAX - 1; }
                        if (prev < slice_hi) {
                            switch (run) {
                            case 1: c1++; break; case 2: c2++; break; case 3: c3++; break; case 4: c4++; break;
                            case 5: c5++; break; case 6: c6++; break; case 7: c7++; break;
                            default: atomicAdd(&s_runs[run], 1u);
                            }
                        } else {
                            for (u32 n = run_out + 1; n <= run; n++) atomicAdd(&s_extra[n], 1u);
                        }
                        if (run > M) {
                            record_from(sh, run_start, prefix + run_start_cnt, run - M, run, slice_hi, s_first, s_R, s_ridx);
                            M = run;
                        }
                        run = 1;
                        run_out = (prev >= slice_hi);
                        run_start = prev;
                        run_start_cnt = cnt;
                        if (prev >= slice_hi) { done = true; break; }   /* every prime < slice_hi is settled */
                    }
                    prevgap = g;
                    prev = off;
                    if (off < slice_hi) cnt++;
                }
            }
            if (!done && have_prev) {
                u32 startp = run ? run_start : prev;
                if (startp < slice_hi) atomicMin(&s_unsettled, startp);
            }
        }
    }
    {
        u32 v[7] = { c1, c2, c3, c4, c5, c6, c7 };
        #pragma unroll
        for (int n = 0; n < 7; n++) {
            u32 s = v[n];
            for (int o = 16; o > 0; o >>= 1) s += __shfl_down_sync(0xffffffffu, s, o);
            if ((tid & 31) == 0 && s) atomicAdd(&s_runs[n + 1], s);
        }
    }
    __syncthreads();
    if (tid < NMAX) {
        results[seg].first[tid] = (u32)(s_first[tid] >> 32);
        results[seg].firstidx[tid] = (u32)s_first[tid];
        results[seg].runs[tid] = s_runs[tid];
        results[seg].extra[tid] = s_extra[tid];
    }
    if (tid == 0) { results[seg].unsettled = s_unsettled; results[seg].overflow = s_overflow; }
}

/* ------------------------------------------------------------------ */
/* Host side: sieving primes, pattern, GPU context                     */
/* ------------------------------------------------------------------ */

static u8 H_R[NRES], H_RIDX[WHEEL];       /* residues coprime to 210 and the inverse map (255 = not coprime) */
static void init_wheel_tables(void)
{
    int j = 0;
    memset(H_RIDX, 255, sizeof H_RIDX);
    for (int r = 0; r < WHEEL; r++)
        if (r % 2 && r % 3 && r % 5 && r % 7) { H_R[j] = (u8)r; H_RIDX[r] = (u8)j; j++; }
}

/* all primes <= n (simple odd sieve; n <= ~2e8 takes well under a second) */
static std::vector<u32> primes_upto(u32 n)
{
    std::vector<u32> out;
    if (n < 2) return out;
    std::vector<u8> comp((n >> 1) + 1, 0);           /* comp[i] <-> 2i+1 */
    for (u64 i = 1; (2 * i + 1) * (2 * i + 1) <= n; i++)
        if (!comp[i])
            for (u64 j = (2 * i + 1) * (2 * i + 1) >> 1; j <= (n >> 1); j += 2 * i + 1) comp[j] = 1;
    out.push_back(2);
    for (u64 i = 1; 2 * i + 1 <= n; i++) if (!comp[i]) out.push_back((u32)(2 * i + 1));
    return out;
}

static std::vector<u8> build_pattern(void)
{
    std::vector<u8> pat(PAT_BYTES, 0);
    static const u32 P[4] = { 11, 13, 17, 19 };
    for (u32 p = 0; p < PAT_PERIOD; p++)
        for (int j = 0; j < NRES; j++) {
            u32 n = WHEEL * p + H_R[j];
            for (int k = 0; k < 4; k++) if (n % P[k] == 0) { pat[6 * p + (j >> 3)] |= (u8)(1 << (j & 7)); break; }
        }
    return pat;
}

struct Gpu {
    u32 *d_bitmap, *d_sprimes, *d_lprimes, *d_mstate;
    u32 *d_sinv, *d_ssegmod, *d_rchunk;     /* per medium prime: floor(2^32/q), SEG_PERIODS mod q, (chunk_lo/210) mod q */
    u8  *d_pattern;
    u64 *d_entries;                 /* SLOTS bucket lists of cap entries for the sleeping large primes */
    u32 *d_counts, *d_err, *h_cnt;  /* per-list counts, error flags, pinned copy of {count, err} */
    u32 cap, nd;                    /* list capacity; number of dense large primes (q < CHUNK_NUMS/6) */
    SegResult *d_results, *h_results;
    std::vector<u32> lprimes;       /* large sieving primes (> qsplit), ascending */
    u32 nsmall, nwarpprimes, nlarge;
    u32 qsplit;
    cudaStream_t sA;
    cudaEvent_t ev_t0, ev_t1, ev_t2, ev_t3;
};

static void gpu_init(Gpu &g, u64 max_number, u32 qsplit)
{
    init_wheel_tables();
    {
        u8 gap[NRES], inv[NRES], T[NRES][NRES], J[NRES][NRES];
        for (int j = 0; j < NRES; j++) {
            gap[j] = (u8)((j + 1 < NRES ? H_R[j + 1] : H_R[0] + WHEEL) - H_R[j]);
            inv[j] = 0;
            for (int x = 1; x < WHEEL; x++) if ((H_R[j] * x) % WHEEL == 1) { inv[j] = (u8)x; break; }
            for (int b = 0; b < NRES; b++) {
                u32 prod = (u32)H_R[b] * H_R[j];
                T[b][j] = (u8)(prod / WHEEL);
                J[b][j] = H_RIDX[prod % WHEEL];
            }
        }
        CUDA_CHECK(cudaMemcpyToSymbol(c_GAP, gap, sizeof gap));
        CUDA_CHECK(cudaMemcpyToSymbol(c_RIDX, H_RIDX, sizeof H_RIDX));
        CUDA_CHECK(cudaMemcpyToSymbol(c_INV, inv, sizeof inv));
        CUDA_CHECK(cudaMemcpyToSymbol(c_T, T, sizeof T));
        CUDA_CHECK(cudaMemcpyToSymbol(c_J, J, sizeof J));
    }

    u32 sq = (u32)sqrt((double)max_number);
    while ((u64)sq * sq > max_number) sq--;
    while ((u64)(sq + 1) * (sq + 1) <= max_number) sq++;
    std::vector<u32> all = primes_upto(sq);
    std::vector<u32> sp;
    for (u32 q : all) {
        if (q < 23) continue;
        if (q <= qsplit) sp.push_back(q); else g.lprimes.push_back(q);
    }
    g.nsmall = (u32)sp.size();
    g.nwarpprimes = 0;
    while (g.nwarpprimes < g.nsmall && sp[g.nwarpprimes] < PWARP) g.nwarpprimes++;
    g.nlarge = (u32)g.lprimes.size();
    g.qsplit = qsplit;

    std::vector<u8> pat = build_pattern();
    CUDA_CHECK(cudaStreamCreate(&g.sA));
    CUDA_CHECK(cudaMalloc(&g.d_bitmap, (CHUNK_WORDS + OVERLAP) * sizeof(u32)));
    CUDA_CHECK(cudaMalloc(&g.d_results, CHUNK_SEGS * sizeof(SegResult)));
    CUDA_CHECK(cudaMallocHost(&g.h_results, CHUNK_SEGS * sizeof(SegResult)));
    CUDA_CHECK(cudaEventCreate(&g.ev_t0)); CUDA_CHECK(cudaEventCreate(&g.ev_t1));
    CUDA_CHECK(cudaEventCreate(&g.ev_t2)); CUDA_CHECK(cudaEventCreate(&g.ev_t3));
    /* sleeping large primes: q >= CHUNK_NUMS/6 may skip chunks; capacity from the expected hits per chunk at END */
    {
        const u32 qdense = (u32)(CHUNK_NUMS / 10);        /* consecutive multiples are at most 10q apart */
        g.nd = (u32)(std::upper_bound(g.lprimes.begin(), g.lprimes.end(), qdense - 1) - g.lprimes.begin());
        double expect = 0;
        for (u32 i = g.nd; i < g.nlarge; i++) expect += (double)CHUNK_NUMS * NRES / (double)WHEEL / g.lprimes[i];
        u64 cap = (u64)(2.0 * expect) + (1u << 20);
        if (cap > 0xFFFFFFF0ULL) die("bucket capacity too large");
        g.cap = (u32)cap;
        CUDA_CHECK(cudaMalloc(&g.d_entries, (u64)SLOTS * g.cap * sizeof(u64)));
        CUDA_CHECK(cudaMalloc(&g.d_counts, SLOTS * sizeof(u32)));
        CUDA_CHECK(cudaMalloc(&g.d_err, sizeof(u32)));
        CUDA_CHECK(cudaMemset(g.d_counts, 0, SLOTS * sizeof(u32)));
        CUDA_CHECK(cudaMemset(g.d_err, 0, sizeof(u32)));
        CUDA_CHECK(cudaMallocHost(&g.h_cnt, 2 * sizeof(u32)));
        if (!g_quiet_init)
            fprintf(stderr, "large primes: %u dense (< %u), %u sleepers in %d bucket lists of %u entries (%.1f GB)\n",
                    g.nd, qdense, g.nlarge - g.nd, SLOTS, g.cap, (double)SLOTS * g.cap * 8 / 1e9);
    }
    CUDA_CHECK(cudaMalloc(&g.d_pattern, PAT_BYTES));
    CUDA_CHECK(cudaMemcpy(g.d_pattern, pat.data(), PAT_BYTES, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&g.d_sprimes, (g.nsmall + 1) * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(g.d_sprimes, sp.data(), g.nsmall * sizeof(u32), cudaMemcpyHostToDevice));
    {
        std::vector<u32> inv(g.nsmall), segmod(g.nsmall);
        for (u32 i = 0; i < g.nsmall; i++) {
            const u32 q = sp[i];
            inv[i] = (u32)((1ULL << 32) / q);
            segmod[i] = (u32)(SEG_PERIODS % q);
        }
        CUDA_CHECK(cudaMalloc(&g.d_sinv, (g.nsmall + 1) * sizeof(u32)));
        CUDA_CHECK(cudaMalloc(&g.d_ssegmod, (g.nsmall + 1) * sizeof(u32)));
        CUDA_CHECK(cudaMalloc(&g.d_rchunk, (g.nsmall + 1) * sizeof(u32)));
        CUDA_CHECK(cudaMemcpy(g.d_sinv, inv.data(), g.nsmall * sizeof(u32), cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(g.d_ssegmod, segmod.data(), g.nsmall * sizeof(u32), cudaMemcpyHostToDevice));
    }
    CUDA_CHECK(cudaMalloc(&g.d_lprimes, (g.nlarge + 1) * sizeof(u32)));
    if (g.nlarge) CUDA_CHECK(cudaMemcpy(g.d_lprimes, g.lprimes.data(), g.nlarge * sizeof(u32), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&g.d_mstate, (g.nlarge + 1) * sizeof(u32)));
    CUDA_CHECK(cudaFuncSetAttribute(kernel_segment, cudaFuncAttributeMaxDynamicSharedMemorySize, SHARED_WORDS * sizeof(u32)));
    if (g_l2_persist) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
        size_t bytes = (CHUNK_WORDS + OVERLAP) * sizeof(u32);
        size_t persist = bytes < (size_t)prop.persistingL2CacheMaxSize ? bytes : (size_t)prop.persistingL2CacheMaxSize;
        CUDA_CHECK(cudaDeviceSetLimit(cudaLimitPersistingL2CacheSize, persist));
        size_t win = bytes;
        if (win > (size_t)prop.accessPolicyMaxWindowSize) win = prop.accessPolicyMaxWindowSize;
        cudaStreamAttrValue attr;
        memset(&attr, 0, sizeof attr);
        attr.accessPolicyWindow.base_ptr = g.d_bitmap;
        attr.accessPolicyWindow.num_bytes = win;
        attr.accessPolicyWindow.hitRatio = 1.0f;
        attr.accessPolicyWindow.hitProp = cudaAccessPropertyPersisting;
        attr.accessPolicyWindow.missProp = cudaAccessPropertyStreaming;
        CUDA_CHECK(cudaStreamSetAttribute(g.sA, cudaStreamAttributeAccessPolicyWindow, &attr));
    }
}

/* number of large primes q with q*q < limit */
static u32 active_count(const Gpu &g, u64 limit)
{
    u32 sq = (u32)sqrt((double)limit);
    while ((u64)sq * sq >= limit && sq) sq--;
    while ((u64)(sq + 1) * (sq + 1) < limit) sq++;
    return (u32)(std::upper_bound(g.lprimes.begin(), g.lprimes.end(), sq) - g.lprimes.begin());
}

/* ------------------------------------------------------------------ */
/* Scan state and merge                                                */
/* ------------------------------------------------------------------ */

struct Cand {
    u64 p;                  /* 0 = none */
    u64 pi;                 /* primes in [start, p] (exact: segments merge in order) */
    u16 gaps[NMAX + 1];
};

static struct {
    u64 start, end;         /* requested bounds */
    u64 next_chunk;         /* next chunk index to process (chunks [k_lo, next) done) */
    u64 k_lo, k_hi;         /* chunk index range [k_lo, k_hi) */
    u64 primes;             /* primes in [start, pos) */
    u64 resume_pos;         /* numbers below this were scanned by an earlier run (= start for a fresh run) */
    u64 hist[NMAX];
    Cand best[NMAX];
    bool overflow;
    u64 cpu_settled;        /* primes settled by the CPU fallback */
} S;

static void print_cand(int n, const Cand *c, const char *tag)
{
    if (stderr_tty) fputs("\r\033[K", stderr);
    printf("%s a(%d) = %" PRIu64 "  gaps:", tag, n, c->p);
    for (int i = 0; i < n; i++) printf(" %u", c->gaps[i]);
    printf(" | %u  (pi(p) = %" PRIu64 ")\n", c->gaps[n], c->pi);
    fflush(stdout);
}

/* record prime p (with prime index pi) as the first of run length n if unset */
static void record_first(int n, u64 p, u64 pi)
{
    if (n < 1 || n >= NMAX || S.best[n].p) return;
    Cand *c = &S.best[n];
    c->p = p;
    c->pi = pi;
    int L = run_length64(p, c->gaps, NMAX + 1);
    if (L != n) die("GPU/CPU disagreement: GPU says L(%" PRIu64 ") = %d, next-prime search says %d", p, n, L);
    if (n >= g_report) print_cand(n, c, "CONFIRMED");
}

/* CPU fallback: settle every prime in [lo, hi) directly; pi_before = primes below lo */
static void settle_cpu(u64 lo, u64 hi, u64 pi_before)
{
    u64 p = is_prime64(lo) ? lo : next_prime64(lo);
    u64 pi = pi_before;
    for (; p < hi; p = next_prime64(p)) {
        u16 gaps[NMAX + 1];
        pi++;
        int L = run_length64(p, gaps, NMAX + 1);
        if (L < 0) die("run longer than %d at %" PRIu64, NMAX, p);
        S.hist[L]++;
        S.cpu_settled++;
        record_first(L, p, pi);
    }
}

/* merge one chunk's segment results (in order) */
static void merge_chunk(const SegResult *results, u64 chunk_lo)
{
    for (u32 s = 0; s < CHUNK_SEGS; s++) {
        const SegResult *r = &results[s];
        const u64 seg_lo = chunk_lo + (u64)s * SEG_NUMS;
        const u64 seg_hi = seg_lo + SEG_NUMS;
        if (seg_hi <= S.resume_pos || seg_lo >= S.end) continue;
        if (r->overflow) S.overflow = true;
        /* first occurrences (segments and chunks arrive in order, so first seen = smallest) */
        for (int n = 1; n < NMAX; n++)
            if (r->first[n] != ~0u) record_first(n, seg_lo + r->first[n], S.primes + r->firstidx[n]);
        /* histogram: primes with run length n = run ends of length >= n, plus final-run extras */
        u64 acc = 0;
        for (int n = NMAX - 1; n >= 1; n--) { acc += r->runs[n]; S.hist[n] += acc + r->extra[n]; }
        if (r->unsettled != ~0u) {
            /* a run was still open at the end of the block's bitmap: settle [p0, hi) on the CPU */
            u64 p0 = seg_lo + r->unsettled, hi = seg_hi < S.end ? seg_hi : S.end;
            u64 before = S.primes;                  /* primes below seg_lo (∩ [start,·)) */
            u64 lo_n = seg_lo > S.resume_pos ? seg_lo : S.resume_pos;
            for (u64 p = is_prime64(lo_n) ? lo_n : next_prime64(lo_n); p < p0; p = next_prime64(p)) before++;
            settle_cpu(p0, hi, before);
        }
        S.primes += r->count;
    }
}

/* ------------------------------------------------------------------ */
/* Checkpoint                                                          */
/* ------------------------------------------------------------------ */

static void save_state(const char *path)
{
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "\nwarning: cannot write %s: %s\n", tmp, strerror(errno)); return; }
    u64 pos = S.next_chunk * CHUNK_NUMS;
    if (pos > S.end) pos = S.end;                /* a partial last chunk is done only up to END */
    fprintf(f, "A158939 cuda checkpoint 2\nstart %" PRIu64 "\nend %" PRIu64 "\npos %" PRIu64
               "\nprimes %" PRIu64 "\ncpu_settled %" PRIu64 "\n", S.start, S.end, pos, S.primes, S.cpu_settled);
    for (int n = 1; n < NMAX; n++)
        if (S.hist[n]) fprintf(f, "hist %d %" PRIu64 "\n", n, S.hist[n]);
    for (int n = 1; n < NMAX; n++) {
        if (!S.best[n].p) continue;
        fprintf(f, "cand %d %" PRIu64 " %" PRIu64, n, S.best[n].p, S.best[n].pi);
        for (int i = 0; i <= n; i++) fprintf(f, " %u", S.best[n].gaps[i]);
        fputc('\n', f);
    }
    fputs("end\n", f);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        fprintf(stderr, "\nwarning: cannot update %s: %s\n", path, strerror(errno));
}

static bool load_state(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return false;
        die("cannot read %s: %s", path, strerror(errno));
    }
    char line[2048];
    if (!fgets(line, sizeof line, f) || strncmp(line, "A158939 cuda checkpoint 2", 25) != 0)
        die("%s is not an A158939 cuda checkpoint", path);
    u64 start = 0, end = 0, pos = 0, primes = 0, cpu = 0;
    bool complete = false;
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "start %" SCNu64, &start) == 1) continue;
        if (sscanf(line, "end %" SCNu64, &end) == 1) continue;
        if (sscanf(line, "pos %" SCNu64, &pos) == 1) continue;
        if (sscanf(line, "primes %" SCNu64, &primes) == 1) continue;
        if (sscanf(line, "cpu_settled %" SCNu64, &cpu) == 1) continue;
        if (strncmp(line, "hist ", 5) == 0) {
            int n = 0; u64 c = 0;
            if (sscanf(line, "hist %d %" SCNu64, &n, &c) != 2 || n < 1 || n >= NMAX) die("%s: bad hist line", path);
            S.hist[n] = c;
            continue;
        }
        if (strncmp(line, "cand ", 5) == 0) {
            char *e;
            long n = strtol(line + 5, &e, 10);
            if (n < 1 || n >= NMAX) die("%s: bad cand line", path);
            S.best[n].p = strtoull(e, &e, 10);
            S.best[n].pi = strtoull(e, &e, 10);
            for (int i = 0; i <= n; i++) S.best[n].gaps[i] = (u16)strtoul(e, &e, 10);
            continue;
        }
        if (strncmp(line, "end", 3) == 0) { complete = true; break; }
    }
    fclose(f);
    if (!complete) die("%s is truncated", path);
    if (start != S.start) die("%s was written for START %" PRIu64 ", not %" PRIu64, path, start, S.start);
    if (end != S.end) fprintf(stderr, "note: checkpoint END was %" PRIu64 ", continuing to %" PRIu64 "\n", end, S.end);
    if (pos < S.start) die("%s: position %" PRIu64 " is below START", path, pos);
    S.next_chunk = pos / CHUNK_NUMS;             /* the chunk containing pos is rescanned from pos */
    S.resume_pos = pos;
    S.primes = primes;
    S.cpu_settled = cpu;
    return true;
}

static void on_sigint(int sig) { (void)sig; g_stop = 1; }

/* ------------------------------------------------------------------ */
/* Scan driver                                                         */
/* ------------------------------------------------------------------ */

struct ScanStats { double seconds; u64 numbers; u64 primes; };

static Gpu G;
static bool G_ready = false;
static bool g_profile = false;
static double prof_large, prof_segment, prof_memset, prof_host;   /* seconds */

static void gpu_free(Gpu &g)
{
    cudaDeviceSynchronize();
    cudaFree(g.d_bitmap); cudaFree(g.d_results); cudaFreeHost(g.h_results);
    cudaEventDestroy(g.ev_t0); cudaEventDestroy(g.ev_t1); cudaEventDestroy(g.ev_t2); cudaEventDestroy(g.ev_t3);
    cudaFree(g.d_pattern); cudaFree(g.d_sprimes); cudaFree(g.d_lprimes); cudaFree(g.d_mstate);
    cudaFree(g.d_sinv); cudaFree(g.d_ssegmod); cudaFree(g.d_rchunk);
    cudaFree(g.d_entries); cudaFree(g.d_counts); cudaFree(g.d_err); cudaFreeHost(g.h_cnt);
    cudaStreamDestroy(g.sA);
    g.lprimes.clear();
}

static ScanStats run_scan(u64 start, u64 end, u32 qsplit, const char *state, int interval, int stop_n, bool quiet)
{
    if (end <= start) die("END must be larger than START");
    memset(&S, 0, sizeof S);
    S.start = start;
    S.end = end;
    S.k_lo = start / CHUNK_NUMS;
    S.k_hi = (end + CHUNK_NUMS - 1) / CHUNK_NUMS;
    S.next_chunk = S.k_lo;
    S.resume_pos = start;
    const bool resumed = state && load_state(state);
    if (resumed)
        fprintf(stderr, "resuming from %s at position %" PRIu64 " (%" PRIu64 " primes so far)\n",
                state, S.resume_pos, S.primes);
    const u64 max_number = S.k_hi * CHUNK_NUMS + OVERLAP_NUMS;
    if (G_ready) gpu_free(G);
    g_quiet_init = quiet;
    gpu_init(G, max_number, qsplit);
    G_ready = true;

    if (!quiet)
        fprintf(stderr, "scanning [%" PRIu64 ", %" PRIu64 ") in %" PRIu64 " chunks of %" PRIu64
                        " (%u small + %u large sieving primes, QSPLIT %u)%s\n",
                start, end, S.k_hi - S.k_lo, (u64)CHUNK_NUMS, G.nsmall, G.nlarge, qsplit,
                stop_n ? ", stopping when a(N) is found" : "");

    /* the primes 2, 3, 5, 7 are outside the wheel: L(2)=2 (1,2,2), L(3)=1, L(5)=2 (2,4,2), L(7)=1 (4,2) */
    if (!resumed) {
        static const u64 sp[4] = { 2, 3, 5, 7 };
        static const int sL[4] = { 2, 1, 2, 1 };
        for (int i = 0; i < 4; i++)
            if (sp[i] >= start && sp[i] < end) { S.primes++; S.hist[sL[i]]++; record_first(sL[i], sp[i], S.primes); }
    }

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    /* state of the large primes for the first chunk: dense array + bucket lists */
    u64 k = S.next_chunk;
    u32 nprev = active_count(G, k * CHUNK_NUMS + CHUNK_NUMS + OVERLAP_NUMS);
    CUDA_CHECK(cudaMemsetAsync(G.d_counts, 0, SLOTS * sizeof(u32), G.sA));
    CUDA_CHECK(cudaMemsetAsync(G.d_err, 0, sizeof(u32), G.sA));
    if (nprev) {
        kernel_init_state<<<(nprev + 255) / 256, 256, 0, G.sA>>>(G.d_lprimes, G.d_mstate, G.nd, nprev, k,
                                                                 G.d_entries, G.d_counts, G.cap, G.d_err);
        CUDA_CHECK(cudaGetLastError());
    }

    const double t0 = now();
    const u64 k0 = k;
    double last_status = 0, last_ckpt = t0;
    enum { NS = 128 };
    double ts[NS]; u64 ks[NS]; int ns = 1;
    ts[0] = t0; ks[0] = k0;
    char b1[32], b2[32];

    for (; k < S.k_hi && !g_stop; k++) {
        const u64 chunk_lo = k * CHUNK_NUMS, chunk_hi = chunk_lo + CHUNK_NUMS;
        const u64 hi_ext = chunk_hi + OVERLAP_NUMS;
        const u32 nact = active_count(G, hi_ext);
        if (g_profile) CUDA_CHECK(cudaEventRecord(G.ev_t0, G.sA));
        if (nact > nprev)
            kernel_activate<<<(nact - nprev + 255) / 256, 256, 0, G.sA>>>(G.d_lprimes, G.d_mstate, G.nd, nprev, nact, k,
                                                                          G.d_entries, G.d_counts, G.cap, G.d_err);
        CUDA_CHECK(cudaMemsetAsync(G.d_bitmap, 0, (CHUNK_WORDS + OVERLAP) * sizeof(u32), G.sA));
        if (g_profile) CUDA_CHECK(cudaEventRecord(G.ev_t1, G.sA));
        const u32 n_dense = nact < G.nd ? nact : G.nd;
        if (n_dense)
            kernel_large<<<(n_dense + 255) / 256, 256, 0, G.sA>>>(G.d_lprimes, G.d_mstate, n_dense,
                                                                   k, chunk_lo, chunk_hi, hi_ext, G.d_bitmap);
        /* sleepers due in this chunk: read the list length, process, clear the list */
        CUDA_CHECK(cudaMemcpyAsync(G.h_cnt, G.d_counts + (k % SLOTS), sizeof(u32), cudaMemcpyDeviceToHost, G.sA));
        CUDA_CHECK(cudaMemcpyAsync(G.h_cnt + 1, G.d_err, sizeof(u32), cudaMemcpyDeviceToHost, G.sA));
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaStreamSynchronize(G.sA));
        if (G.h_cnt[1])
            die("bucket list error %u at chunk %" PRIu64 " (capacity %u entries)", G.h_cnt[1], k, G.cap);
        const u32 n_in = G.h_cnt[0];
        if (n_in)
            kernel_bucket<<<(n_in + BLOCKDIM_B - 1) / BLOCKDIM_B, BLOCKDIM_B, 0, G.sA>>>(
                G.d_entries, G.d_counts, G.cap, G.d_err, n_in, k, chunk_lo, chunk_hi, hi_ext, G.d_bitmap);
        CUDA_CHECK(cudaMemsetAsync(G.d_counts + (k % SLOTS), 0, sizeof(u32), G.sA));
        if (g_profile) CUDA_CHECK(cudaEventRecord(G.ev_t2, G.sA));
        kernel_chunkmod<<<(G.nsmall + 255) / 256, 256, 0, G.sA>>>(G.d_sprimes, G.d_rchunk, G.nsmall, chunk_lo);
        const u64 lo_bound = (chunk_lo < S.resume_pos) ? S.resume_pos : start;   /* rescanned chunk: count from resume_pos */
        kernel_segment<<<CHUNK_SEGS, BLOCKDIM, SHARED_WORDS * sizeof(u32), G.sA>>>(
            G.d_bitmap, G.d_pattern, G.d_sprimes, G.d_sinv, G.d_ssegmod, G.d_rchunk, G.nsmall, G.nwarpprimes,
            chunk_lo, lo_bound, end, G.d_results);
        if (g_profile) CUDA_CHECK(cudaEventRecord(G.ev_t3, G.sA));
        CUDA_CHECK(cudaMemcpyAsync(G.h_results, G.d_results, CHUNK_SEGS * sizeof(SegResult), cudaMemcpyDeviceToHost, G.sA));
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaStreamSynchronize(G.sA));
        if (g_profile) {
            float ms;
            CUDA_CHECK(cudaEventElapsedTime(&ms, G.ev_t0, G.ev_t1)); prof_memset += ms * 1e-3;
            CUDA_CHECK(cudaEventElapsedTime(&ms, G.ev_t1, G.ev_t2)); prof_large += ms * 1e-3;
            CUDA_CHECK(cudaEventElapsedTime(&ms, G.ev_t2, G.ev_t3)); prof_segment += ms * 1e-3;
        }
        double th0 = now();
        merge_chunk(G.h_results, chunk_lo);
        prof_host += now() - th0;
        nprev = nact;
        S.next_chunk = k + 1;
        if (stop_n > 0 && S.best[stop_n].p) g_stop = 1;

        double t = now();
        if (t - ts[ns - 1] >= 1.0) {
            if (ns < NS) { ts[ns] = t; ks[ns] = S.next_chunk; ns++; }
            else { memmove(ts, ts + 1, (NS - 1) * sizeof *ts); memmove(ks, ks + 1, (NS - 1) * sizeof *ks); ts[NS - 1] = t; ks[NS - 1] = S.next_chunk; }
        }
        bool show = stderr_tty ? (t - last_status >= 1.0) : (t - last_status >= interval);
        if (show && !quiet) {
            last_status = t;
            double rate = (t - ts[0] > 0.5) ? (double)(S.next_chunk - ks[0]) * CHUNK_NUMS / (t - ts[0]) : 0;
            double eta = rate > 0 ? (double)(S.k_hi - S.next_chunk) * CHUNK_NUMS / rate : -1;
            int maxc = 0;
            for (int n = 1; n < NMAX; n++) if (S.best[n].p) maxc = n;
            fprintf(stderr, "%s%s  pos %.6g (%.2f%%)  %.3g/s  %.4g primes  ETA %s  found to n=%d",
                    stderr_tty ? "\r\033[K" : "progress: ", fmt_dur(t - t0, b1, sizeof b1),
                    (double)(S.next_chunk * CHUNK_NUMS), 100.0 * (double)(S.next_chunk - S.k_lo) / (double)(S.k_hi - S.k_lo),
                    rate, (double)S.primes, fmt_dur(eta, b2, sizeof b2), maxc);
            if (stop_n > 0 && !S.best[stop_n].p) fprintf(stderr, "  a(%d) > %.4g", stop_n, (double)(S.next_chunk * CHUNK_NUMS));
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) { last_ckpt = t; save_state(state); }
    }
    double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);
    if (state) save_state(state);
    if (S.overflow) fprintf(stderr, "warning: a run exceeded %d gaps; longer run lengths were clamped\n", NMAX - 1);
    if (S.cpu_settled) fprintf(stderr, "note: %" PRIu64 " primes were settled by the CPU fallback\n", S.cpu_settled);
    if (g_stop && S.next_chunk < S.k_hi && !quiet)
        fprintf(stderr, "stopped at position %" PRIu64 "%s\n", S.next_chunk * CHUNK_NUMS, state ? ", checkpoint saved" : "");
    ScanStats st;
    st.seconds = t1 - t0;
    st.numbers = (S.next_chunk - k0) * CHUNK_NUMS;
    st.primes = S.primes;
    return st;
}

/* position up to which the scan is complete (exclusive) */
static u64 scan_pos(void) { u64 p = S.next_chunk * CHUNK_NUMS; return p < S.end ? p : S.end; }

static void print_hist(void)
{
    const u64 pos = scan_pos();
    printf("\n%" PRIu64 " primes in [%" PRIu64 ", %" PRIu64 ")", S.primes, S.start, pos);
    if (S.start == 0)
        for (int k = 6; k <= 16; k++)
            if (pos == (u64)pow(10, k)) printf(S.primes == PI10[k] ? "  = pi(10^%d)" : "  BUT pi(10^%d) = %" PRIu64, k, PI10[k]);
    printf("\nrun-length distribution (model: gaps i.i.d. exponential, P(L = n) = n/(n+1)!):\n");
    printf("  n  primes with L(p) = n     model        ratio\n");
    double fact = 1;
    for (int n = 1; n < NMAX; n++) {
        fact *= (n + 1);
        double model = (double)S.primes * n / fact;
        if (S.hist[n] == 0 && model < 0.01) break;
        printf(" %2d  %-24" PRIu64 " %-12.4g %.3f\n", n, S.hist[n], model, model > 0 ? (double)S.hist[n] / model : 0);
    }
    fflush(stdout);
}

static void print_table(void)
{
    int maxn = 0;
    for (int n = 1; n < NMAX; n++) if (S.best[n].p) maxn = n;
    printf("\n %2s  %-24s %-21s %s\n", "n", S.start ? "a(n) (>= START)" : "a(n)",
           S.start ? "pi(a(n)) (from START)" : "pi(a(n))", "gaps | gap that ends the run");
    for (int n = 1; n <= maxn; n++) {
        const Cand *c = &S.best[n];
        if (!c->p) { printf(" %2d  (none found)\n", n); continue; }
        printf(" %2d  %-24" PRIu64 " %-21" PRIu64, n, c->p, c->pi);
        for (int i = 0; i < n; i++) printf(" %u", c->gaps[i]);
        printf(" | %u", c->gaps[n]);
        if (S.start == 0 && n <= NKNOWN)
            printf(c->p == KNOWN[n] ? "  = OEIS" : "  DIFFERS from OEIS %" PRIu64, KNOWN[n]);
        else if (S.start == 0)
            printf("  NEW");
        if (S.start == 0 && n >= 2 && n - 2 < 14)
            printf(c->pi == A133697[n - 2] ? ", pi = A133697(%d)" : ", pi DIFFERS from A133697(%d) = %" PRIu64, n - 2, A133697[n - 2]);
        else if (S.start == 0 && n >= 2)
            printf(", pi = new A133697(%d)", n - 2);
        putchar('\n');
    }
    printf(" %2d  > %" PRIu64 "  (no%s prime below this position has run length %d)\n",
           maxn + 1, scan_pos(), S.start ? " such" : "", maxn + 1);
    print_hist();
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    fprintf(stderr,
        "usage: a158939_cuda scan [START] END [-n N] [-r NMIN] [-S FILE] [-i SECS] [-Q QSPLIT] [-N] [-q]\n"
        "       a158939_cuda bench [N] [-d SPAN] [-Q QSPLIT] [-N] [-p]\n"
        "       a158939_cuda selftest [-Q QSPLIT]\n");
    exit(2);
}

static u32 parse_qsplit(const char *s)
{
    u64 v = parse_num(s);
    if (v < PWARP || v > 100000000ULL) die("QSPLIT must be between %d and 10^8", PWARP);
    return (u32)v;
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2];
    int npos = 0, interval = 60, stop_n = 0;
    u32 qsplit = 4000000;
    const char *state = NULL;
    bool quiet = false;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-n") && i + 1 < argc) stop_n = parse_int(argv[++i], 1, NMAX - 1, "N");
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) g_report = parse_int(argv[++i], 1, 1000, "NMIN");
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = parse_int(argv[++i], 1, 86400, "SECS");
        else if (!strcmp(argv[i], "-Q") && i + 1 < argc) qsplit = parse_qsplit(argv[++i]);
        else if (!strcmp(argv[i], "-1") || !strcmp(argv[i], "-L")) ;   /* accepted for compatibility */
        else if (!strcmp(argv[i], "-N")) g_l2_persist = false;
        else if (!strcmp(argv[i], "-q")) quiet = true;
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos < 2) pos[npos++] = parse_num(argv[i]);
        else usage();
    }
    if (npos == 0) usage();
    u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
    ScanStats st = run_scan(start, end, qsplit, state, interval, stop_n, quiet);
    print_table();
    char b[32];
    fprintf(stderr, "%" PRIu64 " numbers, %" PRIu64 " primes in %s (%.3g numbers/s)\n",
            st.numbers, st.primes, fmt_dur(st.seconds, b, sizeof b), st.seconds > 0 ? (double)st.numbers / st.seconds : 0);
    return 0;
}

static int cmd_bench(int argc, char **argv)
{
    u64 N = 1000000000000000ULL, span = 2000000000000ULL;
    u32 qsplit = 4000000;
    int npos = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-d") && i + 1 < argc) span = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-p")) g_profile = true;
        else if (!strcmp(argv[i], "-1") || !strcmp(argv[i], "-L")) ;   /* accepted for compatibility */
        else if (!strcmp(argv[i], "-N")) g_l2_persist = false;
        else if (!strcmp(argv[i], "-Q") && i + 1 < argc) qsplit = parse_qsplit(argv[++i]);
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos++ == 0) N = parse_num(argv[i]);
        else usage();
    }
    g_report = NMAX;
    ScanStats st = run_scan(N, N + span, qsplit, NULL, 60, 0, true);
    double rate = (double)st.numbers / st.seconds;
    printf("bench: [%" PRIu64 ", %" PRIu64 ") QSPLIT %u: %.2f s, %.4g numbers/s, %.4g primes/s\n",
           N, N + span, qsplit, st.seconds, rate, (double)st.primes / st.seconds);
    if (g_profile)
        printf("       GPU time: memset %.2f s, large %.2f s, segment %.2f s; host merge %.2f s (%.0f chunks, %u large primes active)\n",
               prof_memset, prof_large, prof_segment, prof_host, (double)st.numbers / CHUNK_NUMS,
               active_count(G, N + span));
    printf("       time for a scan [0, X) at this rate:");
    static const double X[] = { 1e15, 2e15, 5e15, 1e16 };
    char b[32];
    for (size_t i = 0; i < sizeof X / sizeof *X; i++) printf("  %.0e: %s", X[i], fmt_dur(X[i] / rate, b, sizeof b));
    putchar('\n');
    print_hist();
    return 0;
}

static int cmd_selftest(int argc, char **argv)
{
    u32 qsplit = 4000000;
    int fails = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-Q") && i + 1 < argc) qsplit = parse_qsplit(argv[++i]);
        else usage();
    }
#define CHECK(cond, ...) do { if (cond) printf("ok    " __VA_ARGS__); else { printf("FAIL  " __VA_ARGS__); fails++; } putchar('\n'); fflush(stdout); } while (0)
    CHECK(parse_num("2^54-32") == 18014398509481952ULL && parse_num("1e12") == 1000000000000ULL &&
          parse_num("1.29e14") == 129000000000000ULL && parse_num("10^15+7") == 1000000000000007ULL, "number parsing");
    CHECK(is_prime64(2) && !is_prime64(1) && !is_prime64(9) && is_prime64(1000000007ULL) && !is_prime64(3215031751ULL) &&
          is_prime64(18446744073709551557ULL) && is_prime64(129001208165717ULL), "Miller-Rabin spot checks");
    {
        bool ok = true;
        for (int n = 1; n <= NKNOWN; n++) { u16 g[NMAX + 1]; if (run_length64(KNOWN[n], g, NMAX + 1) != n) ok = false; }
        CHECK(ok, "run lengths of a(1..15) by next-prime search");
    }
    g_report = NMAX;
    char b[32];

    /* (a) [0, 1e11): exact values from the CPU tool (a158939.c) */
    {
        ScanStats st = run_scan(0, 100000000000ULL, qsplit, NULL, 60, 0, true);
        printf("      [0, 1e11): %s, %.3g numbers/s\n", fmt_dur(st.seconds, b, sizeof b), st.numbers / st.seconds);
        CHECK(S.primes == 4118054813ULL, "pi(10^11) = %" PRIu64 " (expected 4118054813)", S.primes);
        bool ok = true;
        for (int n = 1; n <= 12; n++) {
            if (S.best[n].p != KNOWN[n]) { ok = false; printf("      a(%d): GPU %" PRIu64 " expected %" PRIu64 "\n", n, S.best[n].p, KNOWN[n]); }
            if (n >= 2 && S.best[n].pi != A133697[n - 2]) { ok = false; printf("      pi(a(%d)): GPU %" PRIu64 " expected %" PRIu64 "\n", n, S.best[n].pi, A133697[n - 2]); }
        }
        CHECK(ok && S.best[13].p == 0, "a(1..12) and pi(a(n)) = A133697 below 10^11");
        static const u64 H[13] = { 0, 2103848608ULL, 1377895718ULL, 489399710ULL, 120501561ULL, 22556442ULL, 3381243ULL,
                                   421291ULL, 45508ULL, 4337ULL, 374ULL, 20ULL, 1ULL };
        ok = true;
        for (int n = 1; n <= 12; n++) if (S.hist[n] != H[n]) { ok = false; printf("      hist[%d]: GPU %" PRIu64 " expected %" PRIu64 "\n", n, S.hist[n], H[n]); }
        CHECK(ok && S.hist[13] == 0, "run-length histogram of [0, 1e11) matches the CPU tool exactly");
    }
    /* (b) odd bounds inside chunks: [12345000000, 13345000000) */
    {
        run_scan(12345000000ULL, 13345000000ULL, qsplit, NULL, 60, 0, true);
        CHECK(S.primes == 42962974ULL, "primes in [12345000000, 13345000000) = %" PRIu64 " (expected 42962974)", S.primes);
        static const u64 F[11] = { 0, 12345000011ULL, 12345000061ULL, 12345000403ULL, 12345000713ULL, 12345000709ULL,
                                   12345006277ULL, 12345463129ULL, 12345463127ULL, 12358263017ULL, 12525319013ULL };
        static const u64 FI[11] = { 0, 1, 2, 19, 35, 34, 278, 20093, 20092, 570629, 7758320 };
        static const u64 H[11] = { 0, 21967314ULL, 14376273ULL, 5095117ULL, 1250868ULL, 233575ULL, 34929ULL, 4355ULL, 494ULL, 48ULL, 1ULL };
        bool ok = true, okh = true;
        for (int n = 1; n <= 10; n++) {
            if (S.best[n].p != F[n] || S.best[n].pi != FI[n]) { ok = false; printf("      n=%d: GPU %" PRIu64 " (pi %" PRIu64 ") expected %" PRIu64 " (pi %" PRIu64 ")\n", n, S.best[n].p, S.best[n].pi, F[n], FI[n]); }
            if (S.hist[n] != H[n]) { okh = false; printf("      hist[%d]: GPU %" PRIu64 " expected %" PRIu64 "\n", n, S.hist[n], H[n]); }
        }
        CHECK(ok && S.best[11].p == 0, "first primes >= START and their indices match the CPU tool");
        CHECK(okh && S.hist[11] == 0, "run-length histogram matches the CPU tool exactly");
    }
    /* (c) [1e15, 1e15 + 2e12): the CPU tool's bench window */
    {
        ScanStats st = run_scan(1000000000000000ULL, 1002000000000000ULL, qsplit, NULL, 60, 0, true);
        printf("      [1e15, 1e15+2e12): %s, %.3g numbers/s\n", fmt_dur(st.seconds, b, sizeof b), st.numbers / st.seconds);
        CHECK(S.primes == 57904107084ULL, "primes in [1e15, 1e15+2e12) = %" PRIu64 " (expected 57904107084)", S.primes);
        static const u64 H[14] = { 0, 29414157759ULL, 19357439615ULL, 6979639719ULL, 1754686141ULL, 337953666ULL, 52549953ULL,
                                   6834766ULL, 764019ULL, 74339ULL, 6523ULL, 539ULL, 41ULL, 4ULL };
        bool ok = true;
        for (int n = 1; n <= 13; n++) if (S.hist[n] != H[n]) { ok = false; printf("      hist[%d]: GPU %" PRIu64 " expected %" PRIu64 "\n", n, S.hist[n], H[n]); }
        CHECK(ok && S.hist[14] == 0, "run-length histogram of [1e15, 1e15+2e12) matches the CPU tool exactly");
    }
    /* (d) a run that ends inside a chunk, then resumed with a larger END, must equal a straight run */
    {
        const char *st = "selftest_resume.state";
        remove(st);
        run_scan(0, 12345678901ULL, qsplit, st, 60, 0, true);
        run_scan(0, 23456789012ULL, qsplit, st, 60, 0, true);
        u64 primes_r = S.primes, hist_r[NMAX];
        Cand best_r[NMAX];
        memcpy(hist_r, S.hist, sizeof hist_r);
        memcpy(best_r, S.best, sizeof best_r);
        run_scan(0, 23456789012ULL, qsplit, NULL, 60, 0, true);
        bool ok = primes_r == S.primes;
        for (int n = 1; n < NMAX; n++)
            if (hist_r[n] != S.hist[n] || best_r[n].p != S.best[n].p || best_r[n].pi != S.best[n].pi) ok = false;
        remove(st);
        CHECK(ok, "resume across END=12345678901 (inside a chunk) to 23456789012 equals a straight run (%" PRIu64 " primes)", S.primes);
    }
#undef CHECK
    printf(fails ? "SELFTEST FAILED (%d)\n" : "selftest passed\n", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) usage();
    if (!strcmp(argv[1], "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(argv[1], "bench")) return cmd_bench(argc - 2, argv + 2);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
