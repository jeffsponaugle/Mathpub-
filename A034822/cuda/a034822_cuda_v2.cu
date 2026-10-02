/*
 * a034822_cuda.cu -- GPU exhaustive search for palindromic squares of even length L
 * (OEIS A034822 / A263618).  Same method as ../a034822.c (see the header there):
 *
 *   n = t*10^h + B, B < 10^h.  B fixes the low h digits of S = n^2, hence (for a
 *   palindrome) the high h digits P, hence t to a window of width ~10^(D-2h).
 *   B is enumerated by a digit DFS; one window serves the 4 roots B, B + 5*10^(h-1),
 *   10^h - B, 5*10^(h-1) - B; 11 | n fixes t mod 11, so each root has <= 2 candidates.
 *
 * GPU layout: one thread per "base" = the low IB digits of B (unit digit 1..5,
 * unit 5 with tens digit 0..4).  The thread runs the remaining K = h - IB digit
 * levels (top digit 0..4) as template-unrolled loops.  sqrt(P*10^e) over the
 * thread's contiguous P range [Pb, Pb + 10^K) is a quadratic Taylor polynomial in
 * 2^-36 fixed point (one exact 128-bit isqrt per thread), evaluated once per
 * level-(h-1) parent and linearly for its 5 children.
 *
 * Filters per candidate (all 64-bit): digits h,h+1 of S vs L-1-h,L-2-h through
 * D = t^2 + floor(2tB/10^h) - P*10^e mod 2^64; then 4 digit pairs; then an exact
 * 128-bit test of S digits h..2h-1 vs their mirrors.  Survivors (essentially only
 * true palindromes) are returned and fully verified on the host.
 *
 * Only even L with h = L/4 rounded down (window <= 22) and h <= 17 are supported;
 * instantiated for L = 40..70.
 *
 * Usage:
 *   a034822_cuda search L [-b bases_per_launch] [-S state] [-r a:b] [-B blockdim]
 *   a034822_cuda selftest            (L = 40..52 against A263618)
 * Output: "FOUND L n n^2" lines, "DONE ..." at the end; progress on stderr.
 */
#include <cuda_runtime.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef uint64_t u64;
typedef int64_t i64;
typedef uint32_t u32;
typedef unsigned __int128 u128;
typedef __int128 i128;

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); exit(1); } } while (0)

static const int KNOWN[71] = {0,
    4,0,3,0,7,1,5,0,11,0,5,1,19,0,13,1,25,0,18,0,48,1,31,0,70,1,44,2,105,
    0,70,1,153,1,98,3,209,0,132,0,291,1,181,1,384,0,234,2,496,1,301,1,
    636,0,383,0,798,1,474,1,981,0,578,0,1199,2,701,1,1443,-1};

__host__ __device__ constexpr u64 p10(int k) { return k == 0 ? 1ull : 10ull * p10(k - 1); }

#define SLACK 64ull
#ifndef LEAF_UNROLL
#define LEAF_UNROLL 1
#endif
constexpr int kLeafUnroll = LEAF_UNROLL;
#ifndef MINB
#define MINB 1
#endif
#define MAXOUT 4096

struct Out {
    unsigned long long nodes, filt, passB, pass3;
    unsigned int nout;
    u64 dbg[8];
    u64 t[MAXOUT], b[MAXOUT];
};

/* ---------------- device helpers ---------------- */

__device__ __forceinline__ double d_u128_to_double(u128 x)
{
    return (double)(u64)(x >> 64) * 18446744073709551616.0 + (double)(u64)x;
}

/* floor(x/d), quotient < 2^64 */
__device__ __forceinline__ u64 d_udiv128_64(u128 x, u64 d)
{
    double qd = d_u128_to_double(x) / (double)d;
    u64 q = qd >= 18446744073709549568.0 ? ~0ull : (u64)qd;
    i128 r = (i128)(x - (u128)q * d);
    double rd = r < 0 ? -d_u128_to_double((u128)(-r)) : d_u128_to_double((u128)r);
    i64 adj = (i64)(rd / (double)d);
    q += (u64)adj;
    r -= (i128)adj * (i128)d;
    while (r < 0) { q--; r += d; }
    while (r >= (i128)d) { q++; r -= d; }
    return q;
}

__device__ __forceinline__ u64 d_isqrt128(u128 x)
{
    if (x == 0) return 0;
    u64 r = (u64)sqrt(d_u128_to_double(x));
    if (r == 0) r = 1;
    u128 s = ((u128)r + d_udiv128_64(x, r)) >> 1;
    r = (u64)s;
    while ((u128)r * r > x) r--;
    while ((u128)(r + 1) * (r + 1) <= x) r++;
    return r;
}

__constant__ unsigned char c_rev2[100];
__constant__ u64 c_f64[20]; /* floor(2^64/10^k) */

template <int L, int H>
struct Cfg {
    static constexpr int e = L - 3 * H;
    static constexpr u64 PE = p10(e);          /* 10^e */
    static constexpr u64 PEX = p10(e - 2);     /* 10^(e-2) */
    static constexpr u64 P4 = (e >= 9) ? p10(e - 4) : 0;
    static constexpr u64 PH = p10(H);
    static constexpr u64 H5 = 5 * p10(H - 1);
    static constexpr u64 Q25 = 25 * p10(H - 2);
    static constexpr u64 H5M11 = H5 % 11;
    static constexpr u64 PHM11 = PH % 11;
};

struct Base {
    u64 Pb, Rb, S1, FRAC, S2; /* fixed point 2^-36; S2 = s2 * 2^(36+SQ) */
    u32 rb11, rb100;
};

struct Thr {
    unsigned long long nodes, filt, passB, pass3;
};

__device__ __forceinline__ u64 addm(u64 a, u64 b, u64 m) { u64 c = a + b; return c >= m ? c - m : c; }
__device__ __forceinline__ u64 subm(u64 a, u64 b, u64 m) { return a >= b ? a - b : a + m - b; }

/* exact stage 3; returns 1 if S digits h..min(2h-1, L-1-2h) mirror correctly */
template <int L, int H>
__device__ __forceinline__ int exact_check(u64 t, u64 Bp, u64 Bq, u64 P, u64 y)
{
    typedef Cfg<L, H> C;
    u128 Cc = (u128)2 * t * Bp + Bq;
    u64 X = d_udiv128_64(Cc, C::PH);
    u64 c0 = (u64)(Cc - (u128)X * C::PH);
    u128 Y = (u128)t * t + X;
#ifdef DBG_B
    if (Bp == DBG_B + C::H5) printf("DBGX X=%llu c0=%llu Yhi=%llu Ylo=%llu q=%llu want=%llu\n", (unsigned long long)X,
        (unsigned long long)c0, (unsigned long long)(u64)(Y >> 64), (unsigned long long)(u64)Y,
        (unsigned long long)d_udiv128_64(Y, C::PEX), (unsigned long long)(P * 100 + y));
#endif
    if (d_udiv128_64(Y, C::PEX) != P * 100 + y) return 0;
    u64 yq = d_udiv128_64(Y, 10000000000000000000ull);
    u64 Ym = (u64)(Y - (u128)yq * 10000000000000000000ull);
    int pmax = 2 * H - 1;
    if (pmax > L - 1 - 2 * H) pmax = L - 1 - 2 * H;
    u64 pw = 100; /* 10^(p-h) */
    for (int p = H + 2; p <= pmax; p++, pw *= 10) {
        int q = L - 1 - p - 2 * H;
        if (q > 18) continue;
        u64 pq = 1;
        for (int j = 0; j < q; j++) pq *= 10;
#ifdef DBG_B
        if (Bp == DBG_B + C::H5) printf("DBGL p=%d q=%d pw=%llu pq=%llu a=%llu b=%llu Ym=%llu yq=%llu\n", p, q,
            (unsigned long long)pw, (unsigned long long)pq, (unsigned long long)((c0 / pw) % 10),
            (unsigned long long)((Ym / pq) % 10), (unsigned long long)Ym, (unsigned long long)yq);
#endif
        if ((c0 / pw) % 10 != (Ym / pq) % 10) return 0;
    }
    return 1;
}

template <int L, int H>
__device__ __forceinline__ void rare(u64 t, u64 Bp, u64 Bq, u64 P, u64 D, u32 y, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    st.filt++;
#ifdef DBG_T
    if (t == DBG_T) { out->dbg[0] = D; out->dbg[1] = Bp; out->dbg[2] = Bq; out->dbg[3] = P; out->dbg[4] = y; out->dbg[5] = 1; }
#endif
    if (C::P4) {
        u64 low4 = (Bq % 10000 + 2 * (Bp % 10000) * (t % 10000)) % 10000;
        u64 y4 = D / C::P4, rem = D - y4 * C::P4;
        if (!(rem < SLACK || rem > C::P4 - SLACK) &&
            y4 != 100u * c_rev2[low4 % 100] + c_rev2[low4 / 100])
            return;
    }
    st.passB++;
#ifdef DBG_B
    if (Bp == DBG_B + C::H5) printf("DBG passB t=%llu\n", (unsigned long long)t);
#endif
    if (!exact_check<L, H>(t, Bp, Bq, P, y)) return;
    st.pass3++;
    unsigned int k = atomicAdd(&out->nout, 1u);
    if (k < MAXOUT) { out->t[k] = t; out->b[k] = Bp; }
}

template <int L, int H>
__device__ __forceinline__ void root(u64 Bp, u64 Bq, u64 f, u32 b100, u32 q100, u64 P, u64 W0lo,
                                     u64 tlo, u64 T2, u32 k, u32 tm100, u64 thi, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    u64 t0 = tlo + k, t1 = t0 + 11;
    u64 sq0 = T2 + 2 * tlo * k + (u64)(k * k);
    u64 sq1 = sq0 + 22 * t0 + 121;
    u64 D0 = sq0 + __umul64hi(2 * t0, f) - W0lo;
    u64 D1 = sq1 + __umul64hi(2 * t1, f) - W0lo;
    u32 m = (2 * b100) % 100;
    u32 tt = tm100 + k;
    u32 l0 = (q100 + m * tt) % 100, l1 = (q100 + m * (tt + 11)) % 100;
    u32 y0 = c_rev2[l0], y1 = c_rev2[l1];
    bool h0 = (t0 <= thi) & (D0 + SLACK - (u64)y0 * C::PEX < C::PEX + 2 * SLACK);
    bool h1 = (t1 <= thi) & (D1 + SLACK - (u64)y1 * C::PEX < C::PEX + 2 * SLACK);
#ifdef DBG_B
    if (Bp == DBG_B + C::H5) printf("DBGR t0=%llu D0=%llu y0=%u l0=%u q100=%u m=%u tt=%u f=%llu h0=%d thi=%llu\n",
        (unsigned long long)t0, (unsigned long long)D0, y0, l0, q100, m, tt, (unsigned long long)f, (int)h0, (unsigned long long)thi);
#endif
    if (h0) rare<L, H>(t0, Bp, Bq, P, D0, y0, st, out);
    if (h1) rare<L, H>(t1, Bp, Bq, P, D1, y1, st, out);
}

/* leaf: B complete (h digits), P, Q = floor(B^2/10^h), f ~ B 2^64/10^h, b11 = B mod 11,
 * V = sqrt(P*10^e) - Rb and V1 = the same for P+1, fixed point 2^-36 */
template <int L, int H>
__device__ __forceinline__ void leaf(u64 B, u64 P, u64 Q, u64 f, u32 b11, u64 V, u64 V1,
                                     const Base &bs, u32 b100, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    st.nodes++;
    const i64 EPSF = 1ll << 26; /* 2^-10 */
    i64 a = ((i64)V - EPSF) >> 36;
    i64 bb = ((i64)V1 + EPSF) >> 36;
    u64 tlo = bs.Rb + (u64)a, thi = bs.Rb + (u64)bb;
    u64 W0lo = P * C::PE, T2 = tlo * tlo;
    u32 nb100 = (100 - b100) % 100;
    u32 q100 = (u32)(Q % 100);
    u32 qB = (q100 + b100) % 100;
    u32 qC = (q100 + 200 - 2 * b100) % 100;
    u32 qD = (q100 + 100 - b100) % 100;
    u32 tm11 = (u32)((bs.rb11 + 11 + a) % 11), tm100 = (u32)((bs.rb100 + 100 + a) % 100);
    u32 rA = b11, rB = (b11 + C::H5M11) % 11, rC = (C::PHM11 + 11 - b11) % 11, rD = (C::H5M11 + 11 - b11) % 11;
    if (!(H & 1)) { rA = (11 - rA) % 11; rB = (11 - rB) % 11; rC = (11 - rC) % 11; rD = (11 - rD) % 11; }
    u32 kA = (rA + 11 - tm11) % 11, kB = (rB + 11 - tm11) % 11, kC = (rC + 11 - tm11) % 11, kD = (rD + 11 - tm11) % 11;
#ifdef DBG_B
    if (B == DBG_B) printf("DBG B=%llu P=%llu Q=%llu tlo=%llu thi=%llu a=%lld V=%llu V1=%llu k=%u %u %u %u tm11=%u b11=%u Rb=%llu\n",
        (unsigned long long)B, (unsigned long long)P, (unsigned long long)Q, (unsigned long long)tlo, (unsigned long long)thi,
        (long long)a, (unsigned long long)V, (unsigned long long)V1, kA, kB, kC, kD, tm11, b11, (unsigned long long)bs.Rb);
#endif
    root<L, H>(B, Q, f, b100, q100, P, W0lo, tlo, T2, kA, tm100, thi, st, out);
    root<L, H>(B + C::H5, Q + B + C::Q25, f + (1ull << 63), b100, qB, P, W0lo, tlo, T2, kB, tm100, thi, st, out);
    root<L, H>(C::PH - B, C::PH - 2 * B + Q, 0 - f, nb100, qC, P, W0lo, tlo, T2, kC, tm100, thi, st, out);
    root<L, H>(C::H5 - B, C::Q25 - B + Q, (1ull << 63) - f, nb100, qD, P, W0lo, tlo, T2, kD, tm100, thi, st, out);
}

/* DFS level I (digits 0..I-1 of B fixed) */
template <int L, int H, int I>
struct Level {
    __device__ __forceinline__ static void run(u64 Bi, u64 Qi, u64 P, u64 f, u32 b11, const Base &bs, int SQ,
                                               u32 b100, Thr &st, Out *out)
    {
        typedef Cfg<L, H> C;
        if constexpr (I == H - 1) {
            /* parent of 5 leaves: P's last digit s is added by the leaf's top digit;
             * V(u0 + s) = V0 + s*D0 with u0 = P - Pb */
            u64 u0 = P - bs.Pb;
            u128 q2 = (u128)(u0 * u0) * bs.S2;
            u64 V0 = bs.FRAC + u0 * bs.S1 - (u64)(q2 >> SQ);
            u128 q3 = (u128)(2 * u0 + 1) * bs.S2;
            u64 D0 = bs.S1 - (u64)(q3 >> SQ);
#pragma unroll (kLeafUnroll)
            for (u32 d = 0; d <= 4; d++) {
                u64 v = Qi + 2 * d * Bi + (u64)(d * d) * p10(I);
                u32 s = (u32)(v % 10);
                u64 Vs = V0 + s * D0;
                leaf<L, H>(Bi + d * p10(I), P + s, v / 10, f + d * c_f64[H - I], (b11 + d * (u32)(p10(I) % 11)) % 11,
                           Vs, Vs + D0, bs, b100, st, out);
            }
        } else {
#pragma unroll 1
            for (u32 d = 0; d <= 9; d++) {
                u64 v = Qi + 2 * d * Bi + (u64)(d * d) * p10(I);
                u64 s = v % 10;
                Level<L, H, I + 1>::run(Bi + d * p10(I), v / 10, P + s * p10(H - 1 - I), f + d * c_f64[H - I],
                                        (b11 + d * (u32)(p10(I) % 11)) % 11, bs, SQ, b100, st, out);
            }
        }
    }
};

template <int L, int H, int IB>
__global__ void __launch_bounds__(256, MINB) search_kernel(u64 base0, u64 nbases, int SQ, Out *out)
{
    typedef Cfg<L, H> C;
    u64 idx = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    Thr st = {0, 0, 0, 0};
    if (idx < nbases) {
        u64 k = base0 + idx;
        u64 per = p10(IB - 1);
        u64 u = 1 + k / per, rest = k % per;
        bool valid = !(u == 5 && (rest % 10) > 4);
        if (valid) {
            u64 Bi = 0, Qi = 0, P = 0, f = 0;
            u32 b11 = 0;
            u64 rr = rest;
            for (int i = 0; i < IB; i++) {
                u64 d = (i == 0) ? u : rr % 10;
                if (i > 0) rr /= 10;
                u64 v = Qi + 2 * d * Bi + d * d * p10(i);
                P += (v % 10) * p10(H - 1 - i);
                Bi += d * p10(i);
                f += d * c_f64[H - i];
                b11 = (u32)((b11 + d * (p10(i) % 11)) % 11);
                Qi = v / 10;
            }
            /* base: sqrt(P*10^e) = Rb + frac, s1 = sqrt'/1, s2 = -sqrt''/2 per unit P */
            Base bs;
            u128 X = (u128)P * C::PE;
            u64 r = d_isqrt128(X);
            double R = (double)r;
            double frac = d_u128_to_double(X - (u128)r * r) / (2.0 * R);
            double Rreal = R + frac, Pd = (double)P;
            double s1 = Rreal / (2.0 * Pd), s2 = Rreal / (8.0 * Pd * Pd);
            bs.Pb = P;
            bs.Rb = r;
            bs.FRAC = (u64)(frac * 68719476736.0);
            bs.S1 = (u64)(s1 * 68719476736.0);
            bs.S2 = (u64)(s2 * ldexp(1.0, 36 + SQ));
            bs.rb11 = (u32)(r % 11);
            bs.rb100 = (u32)(r % 100);
            u32 b100 = (u32)(Bi % 100);
            Level<L, H, IB>::run(Bi, Qi, P, f, b11, bs, SQ, b100, st, out);
        }
    }
    /* warp-reduce stats */
    for (int o = 16; o > 0; o >>= 1) {
        st.nodes += __shfl_down_sync(0xffffffff, st.nodes, o);
        st.filt += __shfl_down_sync(0xffffffff, st.filt, o);
        st.passB += __shfl_down_sync(0xffffffff, st.passB, o);
        st.pass3 += __shfl_down_sync(0xffffffff, st.pass3, o);
    }
    if ((threadIdx.x & 31) == 0) {
        atomicAdd(&out->nodes, st.nodes);
        atomicAdd(&out->filt, st.filt);
        atomicAdd(&out->passB, st.passB);
        atomicAdd(&out->pass3, st.pass3);
    }
}

/* ---------------- host ---------------- */

static int square_digits(u128 n, char *out)
{
    u64 a[6] = {0}, s[12] = {0};
    int na = 0;
    while (n) { a[na++] = (u64)(n % 1000000000u); n /= 1000000000u; }
    for (int i = 0; i < na; i++) {
        u64 carry = 0;
        for (int j = 0; j < na; j++) {
            u64 v = s[i + j] + a[i] * a[j] + carry;
            s[i + j] = v % 1000000000u;
            carry = v / 1000000000u;
        }
        int k = i + na;
        while (carry) { u64 v = s[k] + carry; s[k] = v % 1000000000u; carry = v / 1000000000u; k++; }
    }
    int ns = 2 * na;
    while (ns > 0 && s[ns - 1] == 0) ns--;
    char buf[128];
    int len = 0;
    for (int i = 0; i < ns; i++) {
        u64 v = s[i];
        for (int k = 0; k < 9; k++) { buf[len++] = (char)('0' + v % 10); v /= 10; }
    }
    while (len > 1 && buf[len - 1] == '0') len--;
    for (int i = 0; i < len; i++) out[i] = buf[len - 1 - i];
    out[len] = 0;
    return len;
}

static void u128_str(u128 n, char *out)
{
    char buf[64];
    int len = 0;
    do { buf[len++] = (char)('0' + (int)(n % 10)); n /= 10; } while (n);
    for (int i = 0; i < len; i++) out[i] = buf[len - 1 - i];
    out[len] = 0;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static volatile sig_atomic_t stop_flag = 0;
static void on_sig(int s) { (void)s; stop_flag = 1; }

typedef void (*launch_fn)(u64 base0, u64 nb, int SQ, Out *out, int bd);

template <int L, int H, int IB>
static void launch(u64 base0, u64 nb, int SQ, Out *out, int bd)
{
    u64 blocks = (nb + bd - 1) / bd;
    search_kernel<L, H, IB><<<(unsigned)blocks, bd>>>(base0, nb, SQ, out);
}

struct Variant { int L, H, IB; launch_fn fn; };

/* IB from the cubic Taylor bound: (h+e)/2 + 3(1-IB) - log10(16) < -8, K = h - IB <= 7 */
static constexpr int ib_for(int L, int H)
{
    int e = L - 3 * H;
    for (int ib = 2; ib <= H - 1; ib++) {
        double lg = (H + e) / 2.0 + 3.0 * (1 - ib) - 1.2041199826559248;
        if (lg < -8 && H - ib <= 7) return ib;
    }
    return H - 1;
}

#define VAR(L, H) { L, H, ib_for(L, H), launch<L, H, ib_for(L, H)> }
static const Variant VARIANTS[] = {
    VAR(40, 10), VAR(42, 10), VAR(44, 11), VAR(46, 11), VAR(48, 12), VAR(50, 12), VAR(52, 13), VAR(54, 13),
    VAR(56, 14), VAR(58, 14), VAR(60, 15), VAR(62, 15), VAR(64, 16), VAR(66, 16), VAR(68, 17), VAR(70, 17),
};

static u64 P10H[20];

static u64 run(int L, u64 bpl, const char *state, u64 ra, u64 rb, int bd, int quiet)
{
    const Variant *vv = NULL;
    for (size_t i = 0; i < sizeof VARIANTS / sizeof VARIANTS[0]; i++)
        if (VARIANTS[i].L == L) vv = &VARIANTS[i];
    if (!vv) { fprintf(stderr, "L=%d not compiled in (even 40..70)\n", L); exit(1); }
    int H = vv->H, IB = vv->IB, e = L - 3 * H;
    u64 nbases = 5 * P10H[IB - 1];
    if (rb == 0 || rb > nbases) rb = nbases;
    /* SQ: S2 = s2 * 2^(36+SQ) < 2^63 for the largest s2 = s1max/(4*10^(h-1)) */
    double s1max = sqrt(pow(10.0, e - H + 1)) / 2.0;
    double s2max = s1max / (4.0 * pow(10.0, H - 1));
    int SQ = (int)floor(log2(ldexp(1.0, 62) / (s2max * ldexp(1.0, 36))));
    if (SQ > 90) SQ = 90;
    /* u0^2 * S2 must fit in 128 bits: u0 < 10^K */
    int K = H - IB;
    if (2 * K * log2(10.0) + 63 > 127) { fprintf(stderr, "K too large\n"); exit(1); }
    if (SQ >= 64 + 63) { fprintf(stderr, "bad SQ\n"); exit(1); }

    u64 pos = ra, found = 0;
    unsigned long long nodes = 0, filt = 0, passB = 0, pass3 = 0;
    if (state) {
        FILE *f = fopen(state, "r");
        if (f) {
            int sL; u64 sra, srb, spos, sfound;
            if (fscanf(f, "A034822cuda L=%d range=%" SCNu64 ":%" SCNu64 " pos=%" SCNu64 " found=%" SCNu64
                          " nodes=%llu filt=%llu passB=%llu pass3=%llu",
                       &sL, &sra, &srb, &spos, &sfound, &nodes, &filt, &passB, &pass3) == 9 &&
                sL == L && sra == ra && srb == rb) {
                pos = spos;
                found = sfound;
                if (!quiet) fprintf(stderr, "resumed at base %" PRIu64 "\n", pos);
            } else { fprintf(stderr, "state file mismatch\n"); exit(1); }
            fclose(f);
        }
    }
    if (!quiet)
        fprintf(stderr, "L=%d h=%d e=%d IB=%d K=%d SQ=%d bases [%" PRIu64 ",%" PRIu64 ") of %" PRIu64
                        ", %" PRIu64 " per launch\n", L, H, e, IB, K, SQ, ra, rb, nbases, bpl);
    Out *d_out, *h_out;
    CK(cudaMalloc(&d_out, sizeof(Out)));
    CK(cudaMallocHost(&h_out, sizeof(Out)));
    double t0 = now_sec(), tlast = t0;
    u64 pos0 = pos;
    while (pos < rb && !stop_flag) {
        u64 nb = rb - pos < bpl ? rb - pos : bpl;
        CK(cudaMemset(d_out, 0, sizeof(Out)));
        vv->fn(pos, nb, SQ, d_out, bd);
        CK(cudaGetLastError());
        CK(cudaMemcpy(h_out, d_out, sizeof(Out), cudaMemcpyDeviceToHost));
        if (h_out->nout > MAXOUT) { fprintf(stderr, "output overflow\n"); exit(1); }
        for (unsigned i = 0; i < h_out->nout; i++) {
            u128 n = (u128)h_out->t[i] * P10H[H] + h_out->b[i];
            char ss[128], ns[64];
            int len = square_digits(n, ss);
            int pal = 1;
            for (int a = 0, b = len - 1; a < b; a++, b--) if (ss[a] != ss[b]) pal = 0;
            if (len == L && pal) {
                found++;
                u128_str(n, ns);
                if (!quiet) { printf("FOUND %d %s %s\n", L, ns, ss); fflush(stdout); }
            }
        }
#ifdef DBG_T
        if (h_out->dbg[5]) fprintf(stderr, "DBGT D=%llu Bp=%llu Bq=%llu P=%llu y=%llu\n", (unsigned long long)h_out->dbg[0],
            (unsigned long long)h_out->dbg[1], (unsigned long long)h_out->dbg[2], (unsigned long long)h_out->dbg[3], (unsigned long long)h_out->dbg[4]);
#endif
        nodes += h_out->nodes; filt += h_out->filt; passB += h_out->passB; pass3 += h_out->pass3;
        pos += nb;
        double t = now_sec();
        if (state) {
            char tmp[1024];
            snprintf(tmp, sizeof tmp, "%s.tmp", state);
            FILE *f = fopen(tmp, "w");
            if (f) {
                fprintf(f, "A034822cuda L=%d range=%" PRIu64 ":%" PRIu64 " pos=%" PRIu64 " found=%" PRIu64
                           " nodes=%llu filt=%llu passB=%llu pass3=%llu\n",
                        L, ra, rb, pos, found, nodes, filt, passB, pass3);
                fclose(f);
                rename(tmp, state);
            }
        }
        if (!quiet && (t - tlast > 10 || pos >= rb)) {
            double frac = (double)(pos - ra) / (double)(rb - ra);
            double rate = (double)(pos - pos0) / (t - t0);
            fprintf(stderr, "\r%8.4f%%  base %" PRIu64 "  leaves %.3e  %.3e leaves/s  found %" PRIu64
                            "  eta %.0fs    ", 100 * frac, pos, (double)nodes,
                    (double)nodes / (t - t0 + 1e-9),
                    found, rate > 0 ? (rb - pos) / rate : 0.0);
            tlast = t;
        }
    }
    double el = now_sec() - t0;
    if (!quiet) fprintf(stderr, "\n");
    printf("%s L=%d count=%" PRIu64 " bases=[%" PRIu64 ",%" PRIu64 ") pos=%" PRIu64
           " leaves=%llu filt=%llu passB=%llu pass3=%llu time=%.2fs rate=%.3e leaves/s\n",
           pos >= rb ? "DONE" : "PARTIAL", L, found, ra, rb, pos, nodes, filt, passB, pass3, el,
           (double)nodes / el);
    fflush(stdout);
    cudaFree(d_out);
    cudaFreeHost(h_out);
    return found;
}

int main(int argc, char **argv)
{
    P10H[0] = 1;
    for (int i = 1; i < 20; i++) P10H[i] = P10H[i - 1] * 10;
    unsigned char rev2[100];
    for (int k = 0; k < 100; k++) rev2[k] = (unsigned char)(10 * (k % 10) + k / 10);
    u64 f64[20];
    for (int k = 0; k < 20; k++) f64[k] = k == 0 ? ~0ull : (u64)(((u128)1 << 64) / P10H[k]);
    CK(cudaMemcpyToSymbol(c_rev2, rev2, sizeof rev2));
    CK(cudaMemcpyToSymbol(c_f64, f64, sizeof f64));
    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);
    int bd = 256;
    if (argc >= 2 && !strcmp(argv[1], "selftest")) {
        int maxL = argc >= 3 ? atoi(argv[2]) : 52, bad = 0;
        for (int L = 40; L <= maxL; L += 2) {
            u64 c = run(L, 1 << 20, NULL, 0, 0, bd, 1);
            fprintf(stderr, "L=%d count=%" PRIu64 " expected %d %s\n", L, c, KNOWN[L], (int)c == KNOWN[L] ? "ok" : "MISMATCH");
            if ((int)c != KNOWN[L]) bad++;
        }
        fprintf(stderr, bad ? "SELFTEST FAILED\n" : "SELFTEST PASSED\n");
        return bad != 0;
    }
    if (argc >= 3 && !strcmp(argv[1], "search")) {
        int L = atoi(argv[2]);
        u64 bpl = 1 << 20, ra = 0, rb = 0;
        const char *state = NULL;
        for (int i = 3; i < argc; i++) {
            if (!strcmp(argv[i], "-b") && i + 1 < argc) bpl = strtoull(argv[++i], NULL, 10);
            else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
            else if (!strcmp(argv[i], "-B") && i + 1 < argc) bd = atoi(argv[++i]);
            else if (!strcmp(argv[i], "-r") && i + 1 < argc) {
                if (sscanf(argv[++i], "%" SCNu64 ":%" SCNu64, &ra, &rb) != 2) goto usage;
            } else goto usage;
        }
        run(L, bpl, state, ra, rb, bd, 0);
        return 0;
    }
usage:
    fprintf(stderr, "usage: %s search L [-b bases_per_launch] [-S state] [-r a:b] [-B blockdim]\n"
                    "       %s selftest [maxL]\n", argv[0], argv[0]);
    return 2;
}
