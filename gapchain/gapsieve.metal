/*
 * gapsieve.metal -- the GPU half of gapsieve (macOS, Metal).  The GPU works
 * on units of consecutive segments of one shape (a segment is a whole number
 * of wheel turns, so a unit's segments share their classes): roffs_kernel
 * (where each class starts in each pattern group, per segment), sieve_kernel
 * (AND the rotated patterns, 32-bit words, in two phases, and list the
 * candidates that survive), then a few rounds of prep_kernel and pass_kernel
 * (the base-2 strong test of one member for every remaining candidate,
 * keeping those that pass), and tail_kernel for the members left, once few
 * candidates remain.  The host proves the rare survivors on the CPU.
 *
 * Numbers stay below 2^82.  MSL has no 128-bit integers, and Apple GPUs have
 * no FP64 for the CPU's shortcut to 2^128 mod n.  sprp2w is gapsieve.c's
 * two-word arithmetic, ported, on (lo, hi) pairs of 64-bit words; the
 * kernels use sprp3, the same test in three 32-bit limbs, which is faster
 * (sprp2w stays as the reference it was checked against).
 */
#include <metal_stdlib>
using namespace metal;

struct u128 { ulong lo, hi; };

static inline u128 mk(ulong lo, ulong hi) { u128 r; r.lo = lo; r.hi = hi; return r; }
static inline u128 mul64(ulong a, ulong b) { return mk(a * b, mulhi(a, b)); }

static inline u128 add64(u128 a, ulong b)
{
    u128 r;
    r.lo = a.lo + b;
    r.hi = a.hi + (r.lo < b ? 1ul : 0ul);
    return r;
}

static inline u128 sub128(u128 a, u128 b)
{
    u128 r;
    r.lo = a.lo - b.lo;
    r.hi = a.hi - b.hi - (a.lo < b.lo ? 1ul : 0ul);
    return r;
}

static inline bool ge128(u128 a, u128 b) { return a.hi > b.hi || (a.hi == b.hi && a.lo >= b.lo); }
static inline bool eq128(u128 a, u128 b) { return a.hi == b.hi && a.lo == b.lo; }
static inline int top128(u128 x) { return x.hi ? 127 - (int)clz(x.hi) : 63 - (int)clz(x.lo); }

/* a + b mod n, for a, b < n < 2^127 */
static inline u128 addmod2(u128 a, u128 b, u128 n)
{
    u128 y = sub128(n, b);
    if (ge128(a, y))
        return sub128(a, y);
    u128 r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo ? 1ul : 0ul);
    return r;
}

/* Montgomery product a*b/2^128 mod n: gapsieve.c's mmul2 (CIOS, two words) */
static inline u128 mmul2(u128 a, u128 b, u128 n, ulong ninv)
{
    u128 c = mul64(a.lo, b.lo);                     /* t = a*b0 */
    ulong t0 = c.lo;
    c = add64(mul64(a.hi, b.lo), c.hi);
    ulong t1 = c.lo, t2 = c.hi;
    ulong m = t0 * ninv;                            /* t = (t + m*n) / 2^64 */
    c = add64(mul64(m, n.lo), t0);
    c = add64(add64(mul64(m, n.hi), t1), c.hi);
    t0 = c.lo;
    c = add64(mk(c.hi, 0), t2);
    t1 = c.lo;
    t2 = c.hi;

    c = add64(mul64(a.lo, b.hi), t0);               /* t += a*b1 */
    t0 = c.lo;
    c = add64(add64(mul64(a.hi, b.hi), t1), c.hi);
    t1 = c.lo;
    c = add64(mk(c.hi, 0), t2);
    t2 = c.lo;
    m = t0 * ninv;                                  /* t = (t + m*n) / 2^64 */
    c = add64(mul64(m, n.lo), t0);
    c = add64(add64(mul64(m, n.hi), t1), c.hi);
    t0 = c.lo;
    c = add64(mk(c.hi, 0), t2);
    t1 = c.lo;
    t2 = c.hi;

    u128 t = mk(t0, t1);
    return (t2 != 0 || ge128(t, n)) ? sub128(t, n) : t;
}

/* 2^128 mod n by doubling up from n's top bit: no FP64 here */
static inline u128 r128mod(u128 n)
{
    int b = top128(n);
    u128 x = b >= 64 ? mk(0, 1ul << (b - 64)) : mk(1ul << b, 0);
    for (int i = b; i < 128; i++)
        x = addmod2(x, x, n);
    return x;
}

/* gapsieve.c's sprp2w: base-2 strong probable-prime test, odd n > 2 */
static inline bool sprp2w(u128 n)
{
    ulong n0 = n.lo, inv = n0;
    for (int i = 0; i < 5; i++)
        inv *= 2 - n0 * inv;
    ulong ninv = 0 - inv;
    u128 one = r128mod(n);
    u128 nm1 = sub128(n, one);
    u128 d = sub128(n, mk(1, 0));
    int s = d.lo ? (int)ctz(d.lo) : 64 + (int)ctz(d.hi);
    d = s >= 64 ? mk(d.hi >> (s - 64), 0)
                : mk((d.lo >> s) | (s ? d.hi << (64 - s) : 0), d.hi >> s);
    int b = top128(d);
    u128 r = addmod2(one, one, n);                  /* 2 */
    while (--b >= 0) {
        r = mmul2(r, r, n, ninv);
        ulong bit = b >= 64 ? (d.hi >> (b - 64)) & 1 : (d.lo >> b) & 1;
        u128 dbl = addmod2(r, r, n);
        r = bit ? dbl : r;                          /* a select, not a branch */
    }
    if (eq128(r, one) || eq128(r, nm1))
        return true;
    while (--s > 0) {
        r = mmul2(r, r, n, ninv);
        if (eq128(r, nm1))
            return true;
    }
    return false;
}

/*
 * The same test in three 32-bit limbs, which the kernels use.  Apple GPUs
 * multiply 32 bits at a time, so each 64-bit product above costs several
 * multiplies; limbs of 32 bits, with 64-bit sums that cannot overflow,
 * need about half as many instructions.  With R = 2^96 and n < 2^94,
 * Montgomery products of values below 2n stay below 2n, so they are
 * reduced only for the final comparisons.  Every number tested is below
 * 2^82 (PSI13 in gapsieve.c).
 */
struct u96 { uint l0, l1, l2; };

static inline u96 to96(u128 x)
{
    u96 r;
    r.l0 = (uint)x.lo;
    r.l1 = (uint)(x.lo >> 32);
    r.l2 = (uint)x.hi;
    return r;
}

static inline u96 add96(u96 a, u96 b)
{
    u96 r;
    ulong c = (ulong)a.l0 + b.l0;
    r.l0 = (uint)c;
    c = (ulong)a.l1 + b.l1 + (c >> 32);
    r.l1 = (uint)c;
    r.l2 = a.l2 + b.l2 + (uint)(c >> 32);
    return r;
}

static inline u96 sub96(u96 a, u96 b)
{
    u96 r;
    ulong c = (ulong)a.l0 - b.l0;
    r.l0 = (uint)c;
    c = (ulong)a.l1 - b.l1 - (c >> 63);
    r.l1 = (uint)c;
    r.l2 = a.l2 - b.l2 - (uint)(c >> 63);
    return r;
}

static inline bool ge96(u96 a, u96 b)
{
    if (a.l2 != b.l2)
        return a.l2 > b.l2;
    if (a.l1 != b.l1)
        return a.l1 > b.l1;
    return a.l0 >= b.l0;
}

static inline bool eq96(u96 a, u96 b) { return a.l0 == b.l0 && a.l1 == b.l1 && a.l2 == b.l2; }

/* a*b/2^96 mod n, less than 2n, for a, b < 2n (CIOS, three limbs) */
static inline u96 mmul3(u96 a, u96 b, u96 n, uint ninv)
{
    uint t0 = 0, t1 = 0, t2 = 0, t3 = 0;
    uint bl[3] = { b.l0, b.l1, b.l2 };
    for (int i = 0; i < 3; i++) {
        uint bi = bl[i];
        ulong x = (ulong)a.l0 * bi + t0;            /* t += a*b[i] */
        t0 = (uint)x;
        x = (ulong)a.l1 * bi + t1 + (x >> 32);
        t1 = (uint)x;
        x = (ulong)a.l2 * bi + t2 + (x >> 32);
        t2 = (uint)x;
        x = (ulong)t3 + (x >> 32);
        t3 = (uint)x;
        uint t4 = (uint)(x >> 32);
        uint m = t0 * ninv;                         /* t = (t + m*n) / 2^32 */
        x = (ulong)m * n.l0 + t0;
        x = (ulong)m * n.l1 + t1 + (x >> 32);
        t0 = (uint)x;
        x = (ulong)m * n.l2 + t2 + (x >> 32);
        t1 = (uint)x;
        x = (ulong)t3 + (x >> 32);
        t2 = (uint)x;
        t3 = t4 + (uint)(x >> 32);
    }
    u96 r;                                          /* t3 is 0 here */
    r.l0 = t0;
    r.l1 = t1;
    r.l2 = t2;
    return r;
}

/* base-2 strong probable-prime test for odd n, 2 < n < 2^94 */
static inline bool sprp3(u96 n)
{
    uint inv = n.l0;                                /* n^-1 mod 2^32 */
    for (int i = 0; i < 4; i++)
        inv *= 2 - n.l0 * inv;
    uint ninv = 0 - inv;
    int top = n.l2 ? 95 - (int)clz(n.l2) : n.l1 ? 63 - (int)clz(n.l1)
                                                : 31 - (int)clz(n.l0);
    u96 one = { 0, 0, 0 };                          /* 2^96 mod n, by doubling */
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
    u96 nm1 = sub96(n, one);                        /* Montgomery form of -1 */
    u96 n2 = add96(n, n);

    u96 d = n;                                      /* n - 1 = d * 2^s */
    d.l0 &= ~1u;
    int s = d.l0 ? (int)ctz(d.l0) : d.l1 ? 32 + (int)ctz(d.l1) : 64 + (int)ctz(d.l2);
    ulong dlo = ((ulong)d.l1 << 32 | d.l0), dhi = d.l2;
    dlo = s >= 64 ? dhi >> (s - 64) : (dlo >> s) | (s ? dhi << (64 - s) : 0ul);
    dhi = s >= 64 ? 0 : dhi >> s;
    int b = dhi ? 127 - (int)clz(dhi) : 63 - (int)clz(dlo);

    u96 r = add96(one, one);                        /* 2, below 2n */
    while (--b >= 0) {
        r = mmul3(r, r, n, ninv);
        ulong bit = b >= 64 ? (dhi >> (b - 64)) & 1 : (dlo >> b) & 1;
        u96 dbl = add96(r, r);
        if (ge96(dbl, n2))
            dbl = sub96(dbl, n2);
        r = bit ? dbl : r;                          /* a select, not a branch */
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

/* one unit of one shape: gpu_params in gapsieve_gpu.m */
struct segparams {
    ulong base_lo, base_hi;     /* the unit's start */
    ulong lim;                  /* the largest candidate offset in range */
    ulong segsize;              /* a unit is nseg segments of this length */
    uint W, bmodW, A, ngroups, nps, k;
    uint cap;                   /* capacity of each candidate list */
    uint nseg;
    uint dense, pad;            /* pattern groups sieved before compacting */
};

/* one thread per (group, class, segment): class_bitmap's R, where the class
   starts; bq holds each segment's base mod each prime */
kernel void roffs_kernel(constant segparams &sp [[buffer(0)]],
                         device const uint *res [[buffer(1)]],
                         device const uint *bq [[buffer(2)]],
                         device const uint *q [[buffer(3)]],
                         device const uint *qinv [[buffer(4)]],
                         device const uint *crt [[buffer(5)]],
                         device const uint4 *ginfo [[buffer(6)]],   /* P, first, n, word offset */
                         device uint *roff [[buffer(7)]],
                         device uint *order [[buffer(8)]],
                         uint3 id [[thread_position_in_grid]])
{
    uint gi = id.x, j = id.y, sg = id.z;
    if (gi == 0 && j == 0 && sg == 0)
        order[0] = 1;           /* see gpu_encode */
    if (gi >= sp.ngroups || j >= sp.A)
        return;
    uint r = res[j];
    uint s0 = r >= sp.bmodW ? r - sp.bmodW : r + sp.W - sp.bmodW;
    uint4 g = ginfo[gi];
    device const uint *b = bq + sg * sp.nps;
    ulong R = 0;
    for (uint i = g.y; i < g.y + g.z; i++) {
        uint qq = q[i];
        uint t = b[i] + s0 % qq;
        if (t >= qq)
            t -= qq;
        R += (ulong)((ulong)qinv[i] * t % qq) * crt[i];
    }
    roff[(sg * sp.A + j) * sp.ngroups + gi] = (uint)(R % g.x);
}

#define WPT 4                   /* 32-bit words per thread, first phase */
#define TPG 256                 /* threads per threadgroup: 1024 words */

/*
 * One threadgroup per 1024 words of one class of one segment.  Phase one:
 * each thread ANDs the first sp.dense rotated patterns into its WPT words,
 * 32 bits at a time (Apple GPUs are 32-bit machines).  The first patterns
 * hold the smallest primes and strike the most: after three, about 90% of
 * words are already zero.  So the words still standing are gathered in
 * threadgroup memory, and phase two gives each of them one thread for the
 * remaining patterns, with no lane left idle on a zero word.  Surviving
 * candidates go on the list as offsets from the unit's base; the count
 * keeps rising past cap so the host can see an overflow.
 */
kernel void sieve_kernel(constant segparams &sp [[buffer(0)]],
                         device const uint *res [[buffer(1)]],
                         device const uint *pat32 [[buffer(2)]],
                         device const uint4 *ginfo [[buffer(3)]],
                         device const uint *roff [[buffer(4)]],
                         device ulong *cand [[buffer(5)]],
                         device atomic_uint *count [[buffer(6)]],
                         uint3 id [[thread_position_in_grid]],
                         uint lid [[thread_index_in_threadgroup]])
{
    threadgroup uint tw[TPG * WPT], tv[TPG * WPT];
    threadgroup atomic_uint tn;
    uint w0 = id.x * WPT, j = id.y, sg = id.z;  /* the grid is exact */
    if (lid == 0)
        atomic_store_explicit(&tn, 0, memory_order_relaxed);
    device const uint *ro = roff + (sg * sp.A + j) * sp.ngroups;
    uint acc[WPT];
    for (int t = 0; t < WPT; t++)
        acc[t] = ~0u;
    for (uint gi = 0; gi < sp.dense; gi++) {
        uint R = ro[gi];
        device const uint *src = pat32 + 2 * ginfo[gi].w + (R >> 5) + w0;
        uint sh = R & 31, cur = src[0];
        for (int t = 0; t < WPT; t++) {
            uint nxt = src[t + 1];
            acc[t] &= (cur >> sh) | ((nxt << 1) << (31 - sh));
            cur = nxt;
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    uint nz = 0;                    /* gather: one atomic per SIMD group */
    for (int t = 0; t < WPT; t++)
        nz += acc[t] != 0;
    uint k = simd_prefix_exclusive_sum(nz), tot = simd_sum(nz), at = 0;
    if (simd_is_first() && tot)
        at = atomic_fetch_add_explicit(&tn, tot, memory_order_relaxed);
    k += simd_broadcast_first(at);
    for (int t = 0; t < WPT; t++)
        if (acc[t]) {
            tw[k] = w0 + t;
            tv[k++] = acc[t];
        }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint n = atomic_load_explicit(&tn, memory_order_relaxed);
    uint r = res[j];
    ulong s0 = (ulong)sg * sp.segsize +
               (r >= sp.bmodW ? r - sp.bmodW : r + sp.W - sp.bmodW);
    for (uint i = lid; i < n; i += TPG) {
        uint w = tw[i], a = tv[i];
        for (uint gi = sp.dense; gi < sp.ngroups && a; gi++) {
            uint R = ro[gi];
            device const uint *src = pat32 + 2 * ginfo[gi].w + (R >> 5) + w;
            uint sh = R & 31;
            a &= (src[0] >> sh) | ((src[1] << 1) << (31 - sh));
        }
        for (; a; a &= a - 1) {
            ulong off = s0 + (ulong)(w * 32 + ctz(a)) * sp.W;
            if (off <= sp.lim) {
                uint idx = atomic_fetch_add_explicit(count, 1, memory_order_relaxed);
                if (idx < sp.cap)
                    cand[idx] = off;
            }
        }
    }
}

/* size the next pass from the last one's count; the first pass also
   records the sieve's raw count, for the overflow check */
kernel void prep_kernel(device const uint *cnt_in [[buffer(0)]],
                        device atomic_uint *cnt_out [[buffer(1)]],
                        device uint *args [[buffer(2)]],
                        constant uint &cap [[buffer(3)]],
                        device uint *stat [[buffer(4)]],
                        constant uint &first [[buffer(5)]],
                        device uint *order [[buffer(6)]],
                        uint i [[thread_position_in_grid]])
{
    order[0] = 1;
    uint n = min(cnt_in[0], cap);
    if (first)
        stat[0] = cnt_in[0];
    args[0] = max((n + 255) / 256, 1u);
    args[1] = 1;
    args[2] = 1;
    atomic_store_explicit(cnt_out, 0, memory_order_relaxed);
}

/* one member for every remaining candidate; keep those whose member passes */
kernel void pass_kernel(constant segparams &sp [[buffer(0)]],
                        device const ulong *in [[buffer(1)]],
                        device const uint *cnt_in [[buffer(2)]],
                        device ulong *out [[buffer(3)]],
                        device atomic_uint *cnt_out [[buffer(4)]],
                        constant uint &moff [[buffer(5)]],
                        uint i [[thread_position_in_grid]])
{
    if (i >= min(cnt_in[0], sp.cap))
        return;
    ulong off = in[i];
    if (sprp3(to96(add64(add64(mk(sp.base_lo, sp.base_hi), off), moff))))
        out[atomic_fetch_add_explicit(cnt_out, 1, memory_order_relaxed)] = off;
}

/* every member left, for each remaining candidate, stopping at the first
   that fails: by now few candidates remain, and a round each would cost
   more in dispatches than the idle lanes cost here */
kernel void tail_kernel(constant segparams &sp [[buffer(0)]],
                        device const ulong *in [[buffer(1)]],
                        device const uint *cnt_in [[buffer(2)]],
                        device ulong *out [[buffer(3)]],
                        device atomic_uint *cnt_out [[buffer(4)]],
                        constant uint *moffs [[buffer(5)]],
                        constant uint &nm [[buffer(6)]],
                        uint i [[thread_position_in_grid]])
{
    if (i >= min(cnt_in[0], sp.cap))
        return;
    ulong off = in[i];
    u128 p = add64(mk(sp.base_lo, sp.base_hi), off);
    for (uint m = 0; m < nm; m++)
        if (!sprp3(to96(add64(p, moffs[m]))))
            return;
    out[atomic_fetch_add_explicit(cnt_out, 1, memory_order_relaxed)] = off;
}
