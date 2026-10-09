/*
 * dspattern.metal -- the GPU half of dspattern (macOS, Metal).
 *
 * The GPU runs the enumeration and the sieve: one thread takes K consecutive
 * candidates of one work item (an aligned block of Q with a fixed prefix
 * and r free digits whose digit sum is `need`).  It finds the digits of its
 * first candidate from the counting table ("unranking"), then does exactly
 * what run_item does on the CPU: step to the next digit string with the same
 * digit sum, keep p mod M0 = 7*11*13*17*19*23 up to date from per-pattern
 * digit weights, test 7..23 on every candidate, and work out 29..61 and
 * 67..199 from the digits for the few that pass.  Everything is 32-bit
 * arithmetic, which is what Apple GPUs do natively.  Survivors go out as
 * (item, rank) pairs; the CPU rebuilds p from the rank, and screens and
 * proves it.  The layouts below must match dspattern_gpu.m.
 */
#include <metal_stdlib>
using namespace metal;

#define NDIG 32
#ifndef NSP                     /* the host sets these three from dspattern.c */
#define NSP 92
#define NMOD 39
#define MASKW 8
#endif
#define WS 289                  /* ways rows: digit sums 0 .. 9*NDIG */
#define M0 7436429u             /* 7*11*13*17*19*23 */
#ifndef VS                      /* list_kernel: the low digits, V; the host */
#define VS 6                    /* sets both of these */
#define VN 1000000u             /* 10^VS */
#endif
#define SMAX (9 * VS)
#define VLMAX 73                /* list directory entries: 9*VS + 1 for VS <= 8 */

struct gpat {                   /* one pattern */
    ulong mask[NSP][MASKW];     /* bit x: p = x (mod SP[i]) strikes a term */
    uint wt[NMOD][NDIG + 1];    /* digit j of Q adds d*wt[k][j] mod MODS[k] */
    uint nine[NDIG + 1];        /* 9*(wt[0][0] + ... + wt[0][a-1]) mod M0 */
    uint nineu[4][NDIG + 1];    /* 9*(wt[k][VS] + ... + wt[k][VS+a-1]) mod MODS[k] */
    uint pad;
};

struct gitem {                  /* one work item of a unit */
    uint pat, r, need, no9;
    uint tv, pad;               /* list_kernel: the item's V table, as an offset */
    uint pref[(NMOD + 1) & ~1]; /* p mod MODS[k] when the free digits are 0 */
    ulong count;                /* candidates in the item */
    ulong chunk0;               /* its first thread (chunk) in the unit */
};

struct gparams {
    uint nitems, K, cap, pad;
    ulong nchunks;
};

struct gconst {                 /* the sieve primes and their products */
    uint sp[NSP];
    uint mods[NMOD];
    uint first[NMOD + 1];
};

struct gsurv { uint item, pad; ulong rank; };

struct gvlist { uint start[VLMAX], len[VLMAX]; };   /* V lists by digit sum */

/* digit i of the 32 packed in four words as nibbles (selects, not an
   array, so the digits stay in registers) */
static inline uint nib(uint w0, uint w1, uint w2, uint w3, int i)
{
    uint w = i < 8 ? w0 : i < 16 ? w1 : i < 24 ? w2 : w3;
    return (w >> ((uint)(i & 7) * 4)) & 15u;
}

/* the nibbles [0, n) that fall in word k */
static inline uint lowmask(int n, int k)
{
    int lo = 8 * k;
    return n <= lo ? 0u : n >= lo + 8 ? 0xffffffffu : (1u << ((n - lo) * 4)) - 1;
}

/* in word k: nibbles [0, a) set to 9, and nibble a to b */
static inline uint lowpat(int a, uint b, int k)
{
    uint p = 0x99999999u & lowmask(a, k);
    if (a >= 8 * k && a < 8 * k + 8)
        p |= b << ((uint)(a - 8 * k) * 4);
    return p;
}

/* bit 4j set where nibble j of w is nonzero */
static inline uint nonzero(uint w)
{
    return (w | (w >> 1) | (w >> 2) | (w >> 3)) & 0x11111111u;
}

/* the lowest set bit of a 128-bit value in four words, as a bit index */
static inline int lowest(uint a0, uint a1, uint a2, uint a3)
{
    return a0 ? (int)ctz(a0) : a1 ? 32 + (int)ctz(a1) : a2 ? 64 + (int)ctz(a2)
                                                        : 96 + (int)ctz(a3);
}

/* the bits of word k at positions above bit p */
static inline uint above(uint a, int k, int p)
{
    int lo = 32 * k;
    return p < lo ? a : p >= lo + 31 ? 0u : a & ~((2u << (uint)(p - lo)) - 1u);
}

kernel void enum_kernel(constant gparams &gp [[buffer(0)]],
                        device const gpat *pats [[buffer(1)]],
                        device const gitem *items [[buffer(2)]],
                        device const ulong *W0 [[buffer(3)]],
                        device const ulong *W1 [[buffer(4)]],
                        constant gconst &gc [[buffer(5)]],
                        device gsurv *out [[buffer(6)]],
                        device atomic_uint *nout [[buffer(7)]],
                        uint tid [[thread_position_in_grid]])
{
    if (tid >= gp.nchunks)
        return;
    uint lo = 0, hi = gp.nitems - 1;            /* the item holding chunk tid */
    while (lo < hi) {
        uint mid = (lo + hi + 1) >> 1;
        if (items[mid].chunk0 <= tid)
            lo = mid;
        else
            hi = mid - 1;
    }
    device const gitem &it = items[lo];
    device const gpat &P = pats[it.pat];
    ulong k0 = (ulong)(tid - (uint)it.chunk0) * gp.K;
    uint n = (uint)min((ulong)gp.K, it.count - k0);
    int r = (int)it.r;
    bool no9 = it.no9 != 0;
    device const ulong *W = no9 ? W1 : W0;

    /* the digits of candidate k0: at each position, from the top, skip the
       blocks of candidates with each smaller digit */
    uint w0 = 0, w1 = 0, w2 = 0, w3 = 0;
    int left = (int)it.need;
    ulong k = k0;
    for (int pos = r - 1; pos >= 0; pos--) {
        int dig = 0;
        for (; dig <= 9 && dig <= left; dig++) {
            ulong c = pos ? W[pos * WS + (left - dig)]
                          : (ulong)(left == dig && !(no9 && dig == 9));
            if (k < c)
                break;
            k -= c;
        }
        uint v = (uint)dig << ((uint)(pos & 7) * 4);
        if (pos < 8) w0 |= v; else if (pos < 16) w1 |= v; else if (pos < 24) w2 |= v; else w3 |= v;
        left -= dig;
    }
    uint R = it.pref[0];                        /* p mod M0 */
    for (int j = 0; j < r; j++)
        R += nib(w0, w1, w2, w3, j) * P.wt[0][j];
    R %= M0;
    uint c01 = r >= 2 ? (P.wt[0][1] + M0 - P.wt[0][0]) % M0 : 0;
    uint m7 = (uint)P.mask[0][0], m11 = (uint)P.mask[1][0], m13 = (uint)P.mask[2][0];
    uint m17 = (uint)P.mask[3][0], m19 = (uint)P.mask[4][0], m23 = (uint)P.mask[5][0];

    uint done = 0;
    for (;;) {
        /* stage 1, without branches: every lane does the same work */
        bool valid = !(no9 && (w0 & 15u) == 9);         /* H needs exactly t trailing 9s */
        uint struck = ((m7 >> (R % 7u)) | (m11 >> (R % 11u)) | (m13 >> (R % 13u)) |
                       (m17 >> (R % 17u)) | (m19 >> (R % 19u)) | (m23 >> (R % 23u))) & 1u;
        if (valid && !struck) {
            bool ok = true;                     /* stages 2 and 3: rare */
            for (uint mk = 1; mk < NMOD && ok; mk++) {
                uint x = it.pref[mk];           /* stays below 2^32 */
                for (int j = 0; j < r; j++)
                    x += nib(w0, w1, w2, w3, j) * P.wt[mk][j];
                x %= gc.mods[mk];
                for (uint i = gc.first[mk]; i < gc.first[mk + 1]; i++) {
                    uint y = x % gc.sp[i];
                    if ((P.mask[i][y >> 6] >> (y & 63)) & 1) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                uint idx = atomic_fetch_add_explicit(nout, 1, memory_order_relaxed);
                if (idx < gp.cap) {
                    out[idx].item = lo;
                    out[idx].rank = k0 + done;
                }
            }
        }
        if (valid && ++done == n)
            break;

        /* the next digit string with the same digit sum, without loops.
           Usually the units digit can fall by one and the tens rise by one;
           when every lane can, do just that */
        uint d0 = w0 & 15u, d1 = (w0 >> 4) & 15u;
        if (simd_all(d0 >= 1 && d1 < 9)) {
            w0 += 15u;
            R += c01;
            if (R >= M0)
                R -= M0;
            continue;
        }
        /* In general: z is the lowest nonzero digit and i the lowest digit
           above z that is not 9, so the digits between are all 9 and the
           digit sum below i is d_z + 9(i-1-z).  Raise digit i and rebuild
           the digits below it as small as possible (b, then a 9s). */
        int z = lowest(nonzero(w0), nonzero(w1), nonzero(w2), nonzero(w3)) >> 2;
        uint n0 = nonzero(w0 ^ 0x99999999u), n1 = nonzero(w1 ^ 0x99999999u);
        uint n2 = nonzero(w2 ^ 0x99999999u), n3 = nonzero(w3 ^ 0x99999999u);
        int zb = 4 * z + 3;
        int i = lowest(above(n0, 0, zb), above(n1, 1, zb), above(n2, 2, zb), above(n3, 3, zb)) >> 2;
        uint dz = nib(w0, w1, w2, w3, z);
        int rr = (int)dz + 9 * (i - 1 - z) - 1, a = rr / 9;
        uint b = (uint)(rr - 9 * a);
        uint m;
        m = lowmask(i, 0); w0 = (w0 & ~m) | (lowpat(a, b, 0) & m);
        m = lowmask(i, 1); w1 = (w1 & ~m) | (lowpat(a, b, 1) & m);
        m = lowmask(i, 2); w2 = (w2 & ~m) | (lowpat(a, b, 2) & m);
        m = lowmask(i, 3); w3 = (w3 & ~m) | (lowpat(a, b, 3) & m);
        uint inc = 1u << ((uint)(i & 7) * 4);
        if (i < 8) w0 += inc; else if (i < 16) w1 += inc; else if (i < 24) w2 += inc; else w3 += inc;
        /* the weighted digits below i were d_z*wt[z] plus 9s from z+1 to i-1 */
        uint old = (dz * P.wt[0][z] + P.nine[i] + M0 - P.nine[z + 1]) % M0;
        R = (R + P.wt[0][i] + P.nine[a] + b * P.wt[0][a] + (M0 - old)) % M0;
    }
}

/*
 * list_kernel: the same search, with the lanes of a SIMD group in lockstep.
 * Write the free digits as y = U*10^VS + V.  For a given U the candidates are
 * the V with digit sum need - ds(U) (not ending in 9 when the chain
 * carries): one of the precomputed sorted lists, the same for every U with
 * that digit sum.  A SIMD group takes K consecutive candidates, shares one U
 * at a time, and its lanes stride through U's list together.  Stage 1 needs
 * no division: tvp holds V's share of p mod 7, 11, 13, 17 19 and 23, packed
 * in 26 bits (one table per trailing-9 count t, which fixes the digit
 * weights), and each prime's mask is rotated by U's share once per U, so a
 * candidate costs one coalesced load and six shifts.  Moving to the next U
 * is the same on every lane, so the lanes never diverge except for the rare
 * candidates that pass stage 1.  Ranks and survivors are exactly
 * enum_kernel's.
 */

/* bit x of the result is bit (x + k) mod q of mask: the mask for a residue
   of V, given U's share k */
static inline uint rot(uint mask, uint k, uint q)
{
    return ((mask >> k) | (mask << (q - k))) & ((1u << q) - 1u);
}

kernel void list_kernel(constant gparams &gp [[buffer(0)]],
                        device const gpat *pats [[buffer(1)]],
                        device const gitem *items [[buffer(2)]],
                        device const ulong *W0 [[buffer(3)]],
                        device const ulong *W1 [[buffer(4)]],
                        constant gconst &gc [[buffer(5)]],
                        device gsurv *out [[buffer(6)]],
                        device atomic_uint *nout [[buffer(7)]],
                        device const uint *listv [[buffer(8)]],
                        constant gvlist *vl [[buffer(9)]],
                        device const uint *tvp [[buffer(10)]],
                        uint tg [[threadgroup_position_in_grid]],
                        uint sg [[simdgroup_index_in_threadgroup]],
                        uint nsg [[simdgroups_per_threadgroup]],
                        uint lane [[thread_index_in_simdgroup]],
                        uint sw [[threads_per_simdgroup]])
{
    uint chunk = tg * nsg + sg;                 /* one chunk per SIMD group */
    if (chunk >= gp.nchunks)
        return;
    uint lo = 0, hi = gp.nitems - 1;            /* the item holding it */
    while (lo < hi) {
        uint mid = (lo + hi + 1) >> 1;
        if (items[mid].chunk0 <= chunk)
            lo = mid;
        else
            hi = mid - 1;
    }
    device const gitem &it = items[lo];
    device const gpat &P = pats[it.pat];
    ulong k0 = (ulong)(chunk - (uint)it.chunk0) * gp.K;
    uint left = (uint)min((ulong)gp.K, it.count - k0);
    int r = (int)it.r, m = r - VS, need = (int)it.need;
    uint no9 = it.no9;
    device const ulong *W = no9 ? W1 : W0;

    /* the digits of candidate k0, as in enum_kernel (every lane alike) */
    uint w0 = 0, w1 = 0, w2 = 0, w3 = 0;
    int rest = need;
    ulong k = k0;
    for (int pos = r - 1; pos >= 0; pos--) {
        int dig = 0;
        for (; dig <= 9 && dig <= rest; dig++) {
            ulong c = pos ? W[pos * WS + (rest - dig)]
                          : (ulong)(rest == dig && !(no9 && dig == 9));
            if (k < c)
                break;
            k -= c;
        }
        uint v = (uint)dig << ((uint)(pos & 7) * 4);
        if (pos < 8) w0 |= v; else if (pos < 16) w1 |= v; else if (pos < 24) w2 |= v; else w3 |= v;
        rest -= dig;
    }
    /* y = U*10^VS + V.  V's digit sum picks the list; V's place in it is
       its rank among the strings with that digit sum, from the counts */
    int sigma = 0;
    for (int j = 0; j < VS; j++)
        sigma += (int)nib(w0, w1, w2, w3, j);
    uint idx = 0;
    for (int pos = VS - 1, sr = sigma; pos >= 0; pos--) {
        int dv = (int)nib(w0, w1, w2, w3, pos);
        for (int dig = 0; dig < dv && dig <= sr; dig++)
            idx += pos ? (uint)W[pos * WS + (sr - dig)]
                       : (uint)(sr == dig && !(no9 && dig == 9));
        sr -= dv;
    }
    const uint sh = 4 * VS;
    uint u0 = (w0 >> sh) | (w1 << (32 - sh)), u1 = (w1 >> sh) | (w2 << (32 - sh));
    uint u2 = (w2 >> sh) | (w3 << (32 - sh)), u3 = w3 >> sh;
    int su = need - sigma;
    ulong base = k0 - idx;                      /* the rank of U's first candidate */

    /* U's digits' share of p mod M0 and of the stage-2 moduli */
    uint RU0 = it.pref[0], RU1 = it.pref[1], RU2 = it.pref[2], RU3 = it.pref[3];
    for (int j = 0; j < m; j++) {
        uint d = nib(u0, u1, u2, u3, j);
        RU0 += d * P.wt[0][VS + j];
        RU1 += d * P.wt[1][VS + j];
        RU2 += d * P.wt[2][VS + j];
        RU3 += d * P.wt[3][VS + j];
    }
    RU0 %= M0;
    RU1 %= gc.mods[1];
    RU2 %= gc.mods[2];
    RU3 %= gc.mods[3];

    device const uint *tv = tvp + it.tv;
    device const uint *lv = listv + no9 * VN;
    int shi = need, slo = max(0, need - (no9 ? SMAX - 1 : SMAX));
    uint m7 = (uint)P.mask[0][0], m11 = (uint)P.mask[1][0], m13 = (uint)P.mask[2][0];
    uint m17 = (uint)P.mask[3][0], m19 = (uint)P.mask[4][0], m23 = (uint)P.mask[5][0];

    for (;;) {
        uint r7 = rot(m7, RU0 % 7u, 7u), r11 = rot(m11, RU0 % 11u, 11u);
        uint r13 = rot(m13, RU0 % 13u, 13u), r17 = rot(m17, RU0 % 17u, 17u);
        uint r19 = rot(m19, RU0 % 19u, 19u), r23 = rot(m23, RU0 % 23u, 23u);
        uint len = vl[no9].len[sigma], st = vl[no9].start[sigma];
        uint cnt = min(len - idx, left);
        for (uint e = idx + lane; e < idx + cnt; e += sw) {
            uint x = tv[st + e];
            uint struck = (r7 >> (x & 7u)) | (r11 >> ((x >> 3) & 15u)) |
                          (r13 >> ((x >> 7) & 15u)) | (r17 >> ((x >> 11) & 31u)) |
                          (r19 >> ((x >> 16) & 31u)) | (r23 >> ((x >> 21) & 31u));
            if (struck & 1u)
                continue;
            /* stage 2 from U's residues plus V's digits */
            uint vv = lv[st + e], x1 = RU1, x2 = RU2, x3 = RU3;
            for (int j = 0; j < VS; j++, vv /= 10u) {
                uint d = vv % 10u;
                x1 += d * P.wt[1][j];
                x2 += d * P.wt[2][j];
                x3 += d * P.wt[3][j];
            }
            x1 %= gc.mods[1];
            x2 %= gc.mods[2];
            x3 %= gc.mods[3];
            bool ok = true;
            for (uint i = gc.first[1]; i < gc.first[4] && ok; i++) {
                uint xx = i < gc.first[2] ? x1 : i < gc.first[3] ? x2 : x3;
                uint y = xx % gc.sp[i];
                ok = !((P.mask[i][y >> 6] >> (y & 63)) & 1);
            }
            /* stage 3 from all the digits: rarer still */
            for (uint mk = 4; mk < NMOD && ok; mk++) {
                uint xx = it.pref[mk];
                vv = lv[st + e];
                for (int j = 0; j < VS; j++, vv /= 10u)
                    xx += (vv % 10u) * P.wt[mk][j];
                for (int j = 0; j < m; j++)
                    xx += nib(u0, u1, u2, u3, j) * P.wt[mk][VS + j];
                xx %= gc.mods[mk];
                for (uint i = gc.first[mk]; i < gc.first[mk + 1]; i++) {
                    uint y = xx % gc.sp[i];
                    if ((P.mask[i][y >> 6] >> (y & 63)) & 1) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                uint at = atomic_fetch_add_explicit(nout, 1, memory_order_relaxed);
                if (at < gp.cap) {
                    out[at].item = lo;
                    out[at].rank = base + e;
                }
            }
        }
        left -= cnt;
        if (!left)
            break;
        base += len;

        /* The next U whose digit sum su leaves a list (slo <= su <= shi):
           raise the lowest digit i that can rise, with the digits below it
           as small as possible -- digit sum t, as b then a 9s.  Usually
           that is just U + 1.  The residues move by the difference. */
        int i = 0, below = 0;
        uint o0 = 0, o1 = 0, o2 = 0, o3 = 0;    /* the digits below i, weighted */
        for (; i < m; i++) {
            uint di = nib(u0, u1, u2, u3, i);
            int ps = su - below + 1;            /* digits from i up, after raising i */
            if (di < 9 && ps <= shi && slo - ps <= 9 * i)
                break;
            below += (int)di;
            o0 += di * P.wt[0][VS + i];
            o1 += di * P.wt[1][VS + i];
            o2 += di * P.wt[2][VS + i];
            o3 += di * P.wt[3][VS + i];
        }
        if (i >= m)
            break;                              /* cannot happen: the count is exact */
        int ps = su - below + 1, t = max(0, slo - ps), a = t / 9;
        uint b = (uint)(t - 9 * a), msk;
        msk = lowmask(i, 0); u0 = (u0 & ~msk) | (lowpat(a, b, 0) & msk);
        msk = lowmask(i, 1); u1 = (u1 & ~msk) | (lowpat(a, b, 1) & msk);
        msk = lowmask(i, 2); u2 = (u2 & ~msk) | (lowpat(a, b, 2) & msk);
        msk = lowmask(i, 3); u3 = (u3 & ~msk) | (lowpat(a, b, 3) & msk);
        uint inc = 1u << ((uint)(i & 7) * 4);
        if (i < 8) u0 += inc; else if (i < 16) u1 += inc; else if (i < 24) u2 += inc; else u3 += inc;
        su = ps + t;
        RU0 = (RU0 + P.wt[0][VS + i] + P.nineu[0][a] + b * P.wt[0][VS + a] + (M0 - o0 % M0)) % M0;
        uint M1 = gc.mods[1], M2 = gc.mods[2], M3 = gc.mods[3];
        RU1 = (RU1 + P.wt[1][VS + i] + P.nineu[1][a] + b * P.wt[1][VS + a] + (M1 - o1 % M1)) % M1;
        RU2 = (RU2 + P.wt[2][VS + i] + P.nineu[2][a] + b * P.wt[2][VS + a] + (M2 - o2 % M2)) % M2;
        RU3 = (RU3 + P.wt[3][VS + i] + P.nineu[3][a] + b * P.wt[3][VS + a] + (M3 - o3 % M3)) % M3;
        sigma = need - su;
        idx = 0;
    }
}
