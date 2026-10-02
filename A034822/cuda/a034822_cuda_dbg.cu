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
    u64 dbg[40];
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
__constant__ u64 c_p10[20];

template <int L, int H>
struct Cfg {
    static constexpr int e = L - 3 * H;
    static constexpr u64 PE = p10(e);          /* 10^e */
    static constexpr u64 PEX = p10(e - 2);     /* 10^(e-2) */
    static constexpr int NB = (e >= 12) ? 6 : (e >= 9) ? 4 : 0;
    static constexpr u64 PN = p10(e - NB);
    static constexpr u64 PH = p10(H);
    static constexpr u64 H5 = 5 * p10(H - 1);
    static constexpr u64 Q25 = 25 * p10(H - 2);
    static constexpr u64 H5M11 = H5 % 11;
    static constexpr u64 PHM11 = PH % 11;
};

struct Base {
    u64 Pb, Rb, S1, FRAC, S2; /* fixed point 2^-36; S2 = s2 * 2^(36+SQ) */
    u32 rb11, rb100, rb1000;
};

#define HITBUF 16
struct Hit { u64 t, Bp, Bq, D; u32 y; };
struct Thr {
    unsigned long long nodes, filt, passB, pass3;
    unsigned long long sc[32];
    int slot;
    int nh;
    Hit hb[HITBUF];
};

__device__ __forceinline__ u64 addm(u64 a, u64 b, u64 m) { u64 c = a + b; return c >= m ? c - m : c; }
__device__ __forceinline__ u64 subm(u64 a, u64 b, u64 m) { return a >= b ? a - b : a + m - b; }

/* exact stage: S digits h..min(2h-1, L-1-2h) must mirror; integer-only (u128 divides are rare) */
template <int L, int H>
__device__ __forceinline__ int exact_check(u64 t, u64 Bp, u64 Bq, u64 P, u64 y)
{
    typedef Cfg<L, H> C;
    u128 Cc = (u128)2 * t * Bp + Bq;
    u64 X = (u64)(Cc / C::PH);
    u64 c0 = (u64)(Cc - (u128)X * C::PH);
    u128 Y = (u128)t * t + X;
    if ((u64)(Y / C::PEX) != P * 100 + y) return 0;
    u64 Ym = (u64)(Y % 10000000000000000000ull);
    int pmax = 2 * H - 1;
    if (pmax > L - 1 - 2 * H) pmax = L - 1 - 2 * H;
    u64 pw = 100;
    for (int p = H + 2; p <= pmax; p++, pw *= 10) {
        int q = L - 1 - p - 2 * H;
        if (q > 18) continue;
        u64 pq = c_p10[q];
        if ((c0 / pw) % 10 != (Ym / pq) % 10) return 0;
    }
    return 1;
}

__device__ __forceinline__ u32 rev_digits(u32 x, int n)
{
    u32 r = 0;
    for (int i = 0; i < n; i++) { r = 10 * r + x % 10; x /= 10; }
    return r;
}

/* filter hit: check S digits h..h+NB-1 against their mirrors in 64-bit, then exact */
template <int L, int H>
__device__ __forceinline__ void rare(u64 t, u64 Bp, u64 Bq, u64 P, u64 D, u32 y, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    st.filt++;
    if (C::NB) {
        constexpr u64 M = p10(C::NB);
        u64 low = (Bq % M + 2 * (Bp % M) * (t % M)) % M;
        u64 yN = D / C::PN, rem = D - yN * C::PN;
        if (!(rem < SLACK || rem > C::PN - SLACK) && yN != rev_digits((u32)low, C::NB)) return;
    }
    st.passB++;
    if (!exact_check<L, H>(t, Bp, Bq, P, y)) return;
    st.pass3++;
    unsigned int k = atomicAdd(&out->nout, 1u);
    if (k < MAXOUT) { out->t[k] = t; out->b[k] = Bp; }
}

__constant__ u64 c_ce[20];  /* 10^i * E' mod 10^h, E' = -1 (mod 2^h), 1 (mod 5^h) */

/* main filter: 3 digit pairs after P.  low3 = (Bq + 2 t Bp) mod 1000 = S digits h..h+2;
 * the top digits after P are D / 10^(e-3) and must equal rev3(low3); s_y3p[l] = rev3(l) * 10^(e-3). */
__device__ __forceinline__ u32 rev3(u32 l) { u32 a = l / 100, r = l - 100 * a, b = r / 10, c = r - 10 * b; return 100 * c + 10 * b + a; }
__device__ __forceinline__ u32 red(u32 x, u32 m) { return x >= m ? x - m : x; } /* x < 2m */

template <int L, int H, int NE>
__device__ __forceinline__ void root(u64 Bp, u64 Bq, u64 f, u64 G, u32 m, u32 q1000, u32 k, u64 E0,
                                     u64 A2, u64 tlo, u32 tm1000, u64 thi, const u64 *y3p, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    constexpr u64 PX3 = p10(C::e - 3);
    const u64 fh = f >> 32;
    u64 t0 = tlo + k;
    int myslot = st.slot++; (void)myslot;
    /* D = t^2 + floor(2 t f/2^64) - P*10^e (mod 2^64), E0 = tlo^2 - P*10^e, A2 = 2 tlo */
    u64 D0 = E0 + (A2 + k) * k + G + ((2ull * k * fh) >> 32);
    u32 tt = tm1000 + k;
    u32 l0 = (q1000 + m * tt) % 1000;
    bool h0 = (t0 <= thi) & (D0 + SLACK - y3p[l0] < PX3 + 2 * SLACK);
#ifdef FILTER_CHECK
    {
        /* exact D mod 2^64 (exact in the accept region): t^2 + floor((2 t Bp + Bq)/10^h) - P 10^e,
         * with P 10^e = tlo^2 - E0 (mod 2^64) */
        u128 Cc = (u128)2 * t0 * Bp + Bq;
        u64 X = (u64)(Cc / C::PH);
        u64 Dex = t0 * t0 + X - (tlo * tlo - E0);
        u64 w = y3p[l0];
        if (t0 <= thi && Dex >= w && Dex < w + PX3 && !h0) atomicAdd((unsigned long long *)&out->dbg[0], 1ull);
    }
#endif
    st.sc[(myslot & 7) * 2] += h0; if (h0) { int i = st.nh++; st.hb[i].t = t0; st.hb[i].Bp = Bp; st.hb[i].Bq = Bq; st.hb[i].D = D0; st.hb[i].y = rev3(l0) / 10; }
    if (NE == 2) {
        u64 t1 = t0 + 11;
        u32 k1 = k + 11;
        u64 D1 = E0 + (A2 + k1) * k1 + G + ((2ull * k1 * fh) >> 32);
        u32 l1 = (q1000 + m * (tt + 11)) % 1000;
        bool h1 = (t1 <= thi) & (D1 + SLACK - y3p[l1] < PX3 + 2 * SLACK);
        st.sc[(myslot & 7) * 2 + 1] += h1; if (h1) { int i = st.nh++; st.hb[i].t = t1; st.hb[i].Bp = Bp; st.hb[i].Bq = Bq; st.hb[i].D = D1; st.hb[i].y = rev3(l1) / 10; }
    }
}

/* the four roots x, x + H5, 10^h - x, H5 - x (x < H5) sharing one window */
template <int L, int H, int NE>
__device__ __forceinline__ void group(u64 x, u64 Qx, u64 fx, u32 bx11, u32 bx1000, u32 qx1000, u64 E0,
                                      u64 tlo, u32 tm11, u32 tm1000, u64 thi, const u64 *y3p, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    u64 A2 = 2 * tlo;
    u64 G = __umul64hi(A2, fx);
    /* residues mod 11: root r needs t = -(-1)^h r (mod 11); k = (that - tlo) mod 11 */
    u32 rA = bx11, rB = red(bx11 + (u32)C::H5M11, 11), rC = red((u32)C::PHM11 + 11 - bx11, 11),
        rD = red((u32)C::H5M11 + 11 - bx11, 11);
    if (!(H & 1)) { rA = red(11 - rA, 11); rB = red(11 - rB, 11); rC = red(11 - rC, 11); rD = red(11 - rD, 11); }
    u32 kA = red(rA + 11 - tm11, 11), kB = red(rB + 11 - tm11, 11), kC = red(rC + 11 - tm11, 11), kD = red(rD + 11 - tm11, 11);
    u32 m = red(2 * bx1000, 1000), nm = red(1000 - m, 1000); /* 2 b mod 1000 for x, x+H5 and for the negatives */
    root<L, H, NE>(x, Qx, fx, G, m, qx1000, kA, E0, A2, tlo, tm1000, thi, y3p, st, out);
    root<L, H, NE>(x + C::H5, Qx + x + C::Q25, fx + (1ull << 63), G + tlo, m, red(qx1000 + bx1000, 1000), kB,
                   E0, A2, tlo, tm1000, thi, y3p, st, out);
    root<L, H, NE>(C::PH - x, C::PH - 2 * x + Qx, 0 - fx, 2 * tlo - G - 1, nm, red(qx1000 + 1000 - m, 1000), kC,
                   E0, A2, tlo, tm1000, thi, y3p, st, out);
    root<L, H, NE>(C::H5 - x, C::Q25 - x + Qx, (1ull << 63) - fx, tlo - G - 1, nm, red(qx1000 + 1000 - bx1000, 1000), kD,
                   E0, A2, tlo, tm1000, thi, y3p, st, out);
}

/* digit state for B (low I digits) and B2 = B E' mod 10^h */
struct DS {
    u64 B, Q, P, f, B2f, B2, Q2, f2;
    u32 b11, b211;
};

__device__ __forceinline__ bool s2_level(int i, u64 B, bool unit5, int H)
{
    return !unit5 && i >= 1 && i <= H - 2 && B != 0 && (__ffsll((long long)B) - 1) == i - 1;
}

template <int L, int H, int NE>
__device__ __forceinline__ void leafs(const DS &a, bool s2, const Base &bs, int SQ, u32 b10k, u32 b210k,
                                      const u64 *y3p, Thr &st, Out *out)
{
    typedef Cfg<L, H> C;
    constexpr int I = H - 1;
    constexpr u32 TOP10K = (u32)(p10(I) % 10000);
    /* V(u0 + s) = V0 + s*D0 with u0 = P - Pb */
    u64 u0 = a.P - bs.Pb;
    u128 q2 = (u128)(u0 * u0) * bs.S2;
    u64 V0 = bs.FRAC + u0 * bs.S1 - (u64)(q2 >> SQ);
    u128 q3 = (u128)(2 * u0 + 1) * bs.S2;
    u64 D0 = bs.S1 - (u64)(q3 >> SQ);
    u32 q10k = (u32)(a.Q % 10000), q210k = (u32)(a.Q2 % 10000);
    u32 b2top = (u32)(a.B2f / p10(I));
    u64 W0base = a.P * C::PE;
    const u32 b1000 = b10k % 1000, b21000 = b210k % 1000;
    const i64 EPSF = 1ll << 26; /* 2^-10 */
#pragma unroll (kLeafUnroll)
    for (u32 d = 0; d <= 4; d++) {
        u64 v = a.Q + 2 * d * a.B + (u64)(d * d) * p10(I);
        u32 vm = (q10k + 2 * d * b10k + d * d * TOP10K) % 10000; /* v mod 10^4 */
        u32 s = vm % 10;
        u64 B = a.B + d * p10(I), Q = v / 10, P = a.P + s;
        u64 f = a.f + d * c_f64[1];
        u32 b11 = (a.b11 + d * (u32)(p10(I) % 11)) % 11;
        st.nodes++;
        st.slot = 0;
        u64 Vs = V0 + s * D0;
        i64 wa = ((i64)Vs - EPSF) >> 36;
        i64 wb = ((i64)(Vs + D0) + EPSF) >> 36;
        u64 tlo = bs.Rb + (u64)wa, thi = bs.Rb + (u64)wb;
        u64 E0 = tlo * tlo - (W0base + s * C::PE);
        u32 tm11 = (u32)((bs.rb11 + 11 + wa) % 11), tm1000 = (u32)((bs.rb1000 + 1000 + wa) % 1000);
#ifdef NO_ROOTS
        st.filt += (B ^ Q ^ f ^ b11 ^ P ^ E0 ^ tlo ^ tm11 ^ tm1000 ^ thi ^ vm) & 1;
        if (0) {
#else
        group<L, H, NE>(B, Q, f, b11, b1000, vm / 10, E0, tlo, tm11, tm1000, thi, y3p, st, out);
        if (s2) {
#endif
            u32 d2 = (b2top + d) % 10; /* E' = 1 (mod 10): top digit of B2 shifts with d */
            u64 B2 = a.B2 + d2 * p10(I);
            u64 v2 = a.Q2 + 2 * d2 * a.B2 + (u64)(d2 * d2) * p10(I);
            u32 vm2 = (q210k + 2 * d2 * b210k + d2 * d2 * TOP10K) % 10000;
            u64 Q2 = v2 / 10, f2 = a.f2 + d2 * c_f64[1];
            u32 b211 = (a.b211 + d2 * (u32)(p10(I) % 11)) % 11;
            u32 q21000 = vm2 / 10;
            if (B2 >= C::H5) {
                Q2 = Q2 - B2 + C::Q25; B2 -= C::H5; f2 -= 1ull << 63;
                b211 = (b211 + 11 - (u32)C::H5M11) % 11; q21000 = (q21000 + 1000 - b21000) % 1000;
            }
            group<L, H, NE>(B2, Q2, f2, b211, b21000, q21000, E0, tlo, tm11, tm1000, thi, y3p, st, out);
        }
        if (st.nh) {
#pragma unroll 1
            for (int i = 0; i < st.nh; i++) rare<L, H>(st.hb[i].t, st.hb[i].Bp, st.hb[i].Bq, P, st.hb[i].D, st.hb[i].y, st, out);
            st.nh = 0;
    for (int q = 0; q < 32; q++) st.sc[q] = 0;
    st.slot = 0;
        }
    }
}

template <int H, int I>
__device__ __forceinline__ void dstep(const DS &a, DS &n, u32 d)
{
    u64 v = a.Q + 2ull * d * a.B + (u64)(d * d) * p10(I);
    u64 s = v % 10;
    n.P = a.P + s * p10(H - 1 - I);
    n.B = a.B + d * p10(I);
    n.Q = v / 10;
    n.f = a.f + d * c_f64[H - I];
    n.b11 = (a.b11 + d * (u32)(p10(I) % 11)) % 11;
    u64 b2f = (a.B2f + d * c_ce[I]) % p10(H);
    n.B2f = b2f;
    u32 d2 = (u32)((b2f / p10(I)) % 10);
    u64 v2 = a.Q2 + 2ull * d2 * a.B2 + (u64)(d2 * d2) * p10(I);
    n.B2 = a.B2 + d2 * p10(I);
    n.Q2 = v2 / 10;
    n.f2 = a.f2 + d2 * c_f64[H - I];
    n.b211 = (a.b211 + d2 * (u32)(p10(I) % 11)) % 11;
}

/* runtime-level version for the base replay */
__device__ __forceinline__ void dstep_rt(const DS &a, DS &n, u32 d, int I, int H)
{
    u64 pI = c_p10[I], ph = c_p10[H];
    u64 v = a.Q + 2ull * d * a.B + (u64)(d * d) * pI;
    n.P = a.P + (v % 10) * c_p10[H - 1 - I];
    n.B = a.B + d * pI;
    n.Q = v / 10;
    n.f = a.f + d * c_f64[H - I];
    n.b11 = (u32)((a.b11 + d * (pI % 11)) % 11);
    u64 b2f = (a.B2f + d * c_ce[I]) % ph;
    n.B2f = b2f;
    u32 d2 = (u32)((b2f / pI) % 10);
    u64 v2 = a.Q2 + 2ull * d2 * a.B2 + (u64)(d2 * d2) * pI;
    n.B2 = a.B2 + d2 * pI;
    n.Q2 = v2 / 10;
    n.f2 = a.f2 + d2 * c_f64[H - I];
    n.b211 = (u32)((a.b211 + d2 * (pI % 11)) % 11);
}

/* dstep with a runtime level I in [I0, H-2], dispatched to compile-time constants */
template <int H, int I0>
__device__ __forceinline__ void dstep_ct(const DS &a, DS &n, u32 d, int I)
{
    if constexpr (I0 >= H - 1) {
        dstep_rt(a, n, d, I, H);
    } else {
        if (I == I0) dstep<H, I0>(a, n, d);
        else dstep_ct<H, I0 + 1>(a, n, d, I);
    }
}

/* iterative DFS over levels IB..H-2 (stack in local memory, touched once per 5 leaves);
 * the leaf parent level H-1 runs in registers (leafs) */
template <int L, int H, int IB, int NE>
__device__ __forceinline__ void dfs_iter(const DS &a0, bool s2_0, const Base &bs, int SQ, u32 b1000, u32 b21000,
                                         const u64 *y3p, Thr &st, Out *out)
{
    constexpr int K = H - 1 - IB; /* internal levels IB .. H-2 */
    if constexpr (K <= 0) {
        leafs<L, H, NE>(a0, s2_0, bs, SQ, b1000, b21000, y3p, st, out);
    } else {
        DS stk[K];
        unsigned char dig[K], dmx[K];
        bool s2s[K];
        stk[0] = a0;
        {
            bool r = s2_level(IB, a0.B, false, H);
            dmx[0] = r ? 4 : 9;
            s2s[0] = s2_0 | r;
            dig[0] = 0;
        }
        int lv = 0;
#pragma unroll 1
        while (lv >= 0) {
            if (dig[lv] > dmx[lv]) { lv--; continue; }
            u32 d = dig[lv]++;
            DS n;
            dstep_ct<H, IB>(stk[lv], n, d, IB + lv);
            if (lv == K - 1) {
                leafs<L, H, NE>(n, s2s[lv], bs, SQ, b1000, b21000, y3p, st, out);
            } else {
                bool r = s2_level(IB + lv + 1, n.B, false, H);
                bool s2n = s2s[lv] | r;
                lv++;
                stk[lv] = n;
                dmx[lv] = r ? 4 : 9;
                s2s[lv] = s2n;
                dig[lv] = 0;
            }
        }
    }
}

template <int L, int H, int IB>
__global__ void __launch_bounds__(256, MINB) search_kernel(u64 base0, u64 nbases, int SQ, int use_s2, Out *out)
{
    typedef Cfg<L, H> C;
    __shared__ u64 s_y3p[1000];
    for (int i = threadIdx.x; i < 1000; i += blockDim.x) s_y3p[i] = (u64)rev3(i) * p10(C::e - 3);
    __syncthreads();
    u64 idx = (u64)blockIdx.x * blockDim.x + threadIdx.x;
    Thr st;
    st.nodes = st.filt = st.passB = st.pass3 = 0;
    st.nh = 0;
    if (idx < nbases) {
        u64 k = base0 + idx;
        u64 per = p10(IB - 1);
        u64 u = 1 + k / per, rest = k % per;
        bool unit5 = (u == 5);
        bool valid = !(unit5 && (rest % 10) > 4);
        DS a = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, n;
        bool s2 = false;
        u64 rr = rest;
        for (int i = 0; i < IB && valid; i++) {
            u32 d = (i == 0) ? (u32)u : (u32)(rr % 10);
            if (i > 0) rr /= 10;
            if (use_s2 && s2_level(i, a.B, unit5, H)) { if (d > 4) valid = false; s2 = true; }
            dstep_rt(a, n, d, i, H);
            a = n;
        }
        if (valid) {
            Base bs;
            u128 X = (u128)a.P * C::PE;
            u64 r = d_isqrt128(X);
            double R = (double)r;
            double frac = d_u128_to_double(X - (u128)r * r) / (2.0 * R);
            double Rreal = R + frac, Pd = (double)a.P;
            double s1 = Rreal / (2.0 * Pd), s2d = Rreal / (8.0 * Pd * Pd);
            bs.Pb = a.P;
            bs.Rb = r;
            bs.FRAC = (u64)(frac * 68719476736.0);
            bs.S1 = (u64)(s1 * 68719476736.0);
            bs.S2 = (u64)(s2d * ldexp(1.0, 36 + SQ));
            bs.rb11 = (u32)(r % 11);
            bs.rb100 = (u32)(r % 100);
            bs.rb1000 = (u32)(r % 1000);
            u32 b1000 = (u32)(a.B % 10000), b21000 = (u32)(a.B2 % 10000); /* mod 10^4 */
            /* window width <= 11 everywhere in this thread: one candidate per root */
            bool one = s1 < 9.9;
            if (one) dfs_iter<L, H, IB, 1>(a, s2, bs, SQ, b1000, b21000, s_y3p, st, out);
            else dfs_iter<L, H, IB, 2>(a, s2, bs, SQ, b1000, b21000, s_y3p, st, out);
        }
    }
    for (int o = 16; o > 0; o >>= 1) {
        st.nodes += __shfl_down_sync(0xffffffff, st.nodes, o);
        st.filt += __shfl_down_sync(0xffffffff, st.filt, o);
        st.passB += __shfl_down_sync(0xffffffff, st.passB, o);
        st.pass3 += __shfl_down_sync(0xffffffff, st.pass3, o);
    }
    for (int q = 0; q < 16; q++) if (st.sc[q]) atomicAdd((unsigned long long *)&out->dbg[8 + q], st.sc[q]);
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
static int g_use_s2 = 1;

template <int L, int H, int IB>
static void launch(u64 base0, u64 nb, int SQ, Out *out, int bd)
{
    u64 blocks = (nb + bd - 1) / bd;
    search_kernel<L, H, IB><<<(unsigned)blocks, bd>>>(base0, nb, SQ, g_use_s2, out);
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
    {
        /* E' = -1 (mod 2^h), 1 (mod 5^h); c_ce[k] = 10^k E' mod 10^h */
        u64 m2 = 1ull << H, m5 = P10H[H] >> H, M = P10H[H], e1 = 1;
        while (e1 % m2 != m2 - 1) e1 += m5;
        u64 ES = e1 % M, ce[20] = {0};
        if ((u64)(((u128)ES * ES) % M) != 1) { fprintf(stderr, "bad E'\n"); exit(1); }
        for (int k = 0; k < H; k++) ce[k] = (u64)(((u128)P10H[k] * ES) % M);
        CK(cudaMemcpyToSymbol(c_ce, ce, sizeof ce));
    }
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
    static unsigned long long slotacc[16];
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
#ifdef FILTER_CHECK
        if (h_out->dbg[0]) { fprintf(stderr, "FILTER_CHECK: %llu exact hits rejected by the filter!\n", (unsigned long long)h_out->dbg[0]); exit(3); }
#endif
        for (int q = 0; q < 16; q++) slotacc[q] += h_out->dbg[8 + q];
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
    for (int q = 0; q < 16; q++) fprintf(stderr, "slot%d/%d %llu\n", q / 2, q % 2, slotacc[q]);
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
    CK(cudaMemcpyToSymbol(c_p10, P10H, sizeof P10H));
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
