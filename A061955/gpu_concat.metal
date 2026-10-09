// gpu_concat.metal -- GPU kernel for the base-2 reversed-digit families (see concat.c).
// One thread per odd n: computes V(n) mod n with 64-bit Montgomery arithmetic, using
// the same recurrences as concat.c's full_block2 / part_block2.
//   fam: 0 = Lk (A061955), 1 = Ld (A029519), 2 = Rk (A029495), 3 = Rd (A061931)
//   mode 0: store the plain residue for every thread; mode 1: append n with V(n) == 0.
#include <metal_stdlib>
using namespace metal;

struct Params {
    ulong n0;       // first n (odd); thread i tests n0 + 2 i
    uint count;     // number of threads with work
    uint fam;
    uint mode;
    uint maxhits;
};

static inline ulong addm(ulong a, ulong b, ulong q) { ulong s = a + b; return s >= q ? s - q : s; }
static inline ulong subm(ulong a, ulong b, ulong q) { return a >= b ? a - b : a - b + q; }
static inline ulong halfm(ulong a, ulong q) { return (a & 1) ? (a >> 1) + (q >> 1) + 1 : a >> 1; }

static inline ulong mmul(ulong a, ulong b, ulong q, ulong qi)   // a b / 2^64 mod q
{
    ulong lo = a * b, hi = mulhi(a, b);
    ulong m = lo * qi;
    ulong mh = mulhi(m, q);
    ulong r = hi - mh;
    return hi < mh ? r + q : r;
}

static inline int bitlen64(ulong n)
{
    uint h = uint(n >> 32), l = uint(n);
    return h ? 64 - int(clz(h)) : (l ? 32 - int(clz(l)) : 0);
}

static inline int ctz64(ulong n)
{
    uint l = uint(n);
    return l ? int(ctz(l)) : 32 + int(ctz(uint(n >> 32)));
}

// Complete block of Dp-digit words; x = 2^Dp, xh = 2^(Dp-1) (Montgomery form).
static inline void full_block(bool left, bool drop, int Dp, ulong x, ulong xh, ulong one,
                              ulong q, ulong qi, thread ulong &bV, thread ulong &bX)
{
    if (!drop) {
        ulong y = x, G = one, R = 0;
        for (int j = 0; j < Dp - 1; j++) {
            ulong t = mmul(y, G, q, qi), u = mmul(y, R, q, qi);
            ulong r = addm(R, u, q);
            r = addm(r, r, q);
            R = addm(r, left ? t : G, q);
            G = addm(G, t, q);
            y = mmul(y, y, q, qi);
        }
        bV = addm(G, addm(R, R, q), q);
        bX = y;
        return;
    }
    ulong y = xh, G = one, X = one, B = 0;
    for (int j = 0; j < Dp - 1; j++) {
        ulong mu = addm(y, y, q);
        y = mmul(y, y, q, qi);
        ulong Gm = mmul(G, mu, q, qi), Xn = mmul(mu, X, q, qi);
        if (left) {
            ulong U = addm(addm(B, B, q), X, q);
            B = subm(addm(addm(U, mmul(U, mu, q, qi), q), Gm, q), addm(Xn, Xn, q), q);
        } else {
            ulong Bs = addm(B, mmul(B, mu, q, qi), q);
            B = addm(addm(Bs, Bs, q), G, q);
        }
        G = addm(G, Gm, q);
        X = Xn;
    }
    if (left) {
        ulong a = addm(subm(G, X, q), addm(B, B, q), q);
        bV = addm(addm(a, a, q), one, q);
    } else {
        bV = addm(G, addm(B, B, q), q);
    }
    bX = addm(X, X, q);
}

// Last, incomplete block: k = 2^(D-1) .. n, N = n+1 < 2^D.  x = 2^D, xh = 2^(D-1).
static inline void part_block(bool left, bool drop, int D, ulong N, ulong x, ulong xh, ulong one,
                              ulong q, ulong qi, thread ulong &bV, thread ulong &bX)
{
    ulong cb = 0, p = xh;                     // cb = rev_D(N): bit i of N has weight 2^(D-1-i)
    for (int i = 0; i < D; i++) {
        if ((N >> i) & 1) cb = addm(cb, p, q);
        p = halfm(p, q);
    }
    ulong G = one, B = 0, X = drop ? one : x, y = xh, p1 = x;   // p1 = 2^(D-j)
    bV = 0; bX = one;
    for (int j = 0; j < D - 1; j++) {
        ulong p0 = halfm(p1, q);              // 2^(D-1-j)
        if ((N >> j) & 1) {
            cb = subm(cb, p0, q);             // rev_D(Q_j)
            ulong sB = mmul(p1, B, q, qi), V, Xs;
            if (!drop) {
                V = addm(mmul(cb, G, q, qi), sB, q);
                Xs = X;
            } else {
                int vz = j + 1 + ctz64(N >> (j + 1));
                ulong L0 = p1;
                for (int k = j; k < vz; k++) L0 = halfm(L0, q);   // 2^(D - vz)
                Xs = mmul(X, L0, q, qi);
                if (left)
                    V = addm(mmul(addm(mmul(subm(G, X, q), cb, q, qi), sB, q), L0, q, qi), cb, q);
                else
                    V = addm(mmul(cb, G, q, qi), sB, q);
            }
            if (left) bV = addm(mmul(bV, Xs, q, qi), V, q);
            else bV = addm(mmul(V, bX, q, qi), bV, q);
            bX = mmul(bX, Xs, q, qi);
        }
        if (j == D - 2) break;
        if (!drop) {
            ulong yy = X, t = mmul(yy, G, q, qi), u = mmul(yy, B, q, qi);
            ulong r = addm(B, u, q);
            r = addm(r, r, q);
            B = addm(r, left ? t : G, q);
            G = addm(G, t, q);
            X = mmul(yy, yy, q, qi);
        } else {
            ulong mu = addm(y, y, q);
            y = mmul(y, y, q, qi);
            ulong Gm = mmul(G, mu, q, qi), Xn = mmul(mu, X, q, qi);
            if (left) {
                ulong U = addm(addm(B, B, q), X, q);
                B = subm(addm(addm(U, mmul(U, mu, q, qi), q), Gm, q), addm(Xn, Xn, q), q);
            } else {
                ulong Bs = addm(B, mmul(B, mu, q, qi), q);
                B = addm(addm(Bs, Bs, q), G, q);
            }
            G = addm(G, Gm, q);
            X = Xn;
        }
        p1 = p0;
    }
}

kernel void concat2(constant Params &P [[buffer(0)]],
                    device ulong *res [[buffer(1)]],
                    device atomic_uint *nhits [[buffer(2)]],
                    device ulong *hits [[buffer(3)]],
                    uint gid [[thread_position_in_grid]])
{
    if (gid >= P.count) return;
    const ulong n = P.n0 + 2 * ulong(gid);
    const bool left = P.fam < 2, drop = (P.fam & 1) != 0;
    ulong V = 0;
    if (n > 1) {
        const ulong q = n;
        ulong qi = q;
        for (int i = 0; i < 5; i++) qi *= 2 - q * qi;
        const ulong one = (0 - q) % q;
        const int D = bitlen64(n);
        ulong accV = 0, accX = one, x = one, xh = one, bV, bX;
        for (int Dp = 1; Dp <= D; Dp++) {
            xh = x;
            x = addm(x, x, q);                // x = 2^Dp, xh = 2^(Dp-1)
            if (Dp < D || n + 1 == (ulong(1) << D))
                full_block(left, drop, Dp, x, xh, one, q, qi, bV, bX);
            else
                part_block(left, drop, D, n + 1, x, xh, one, q, qi, bV, bX);
            if (left) { accV = addm(mmul(bV, accX, q, qi), accV, q); accX = mmul(bX, accX, q, qi); }
            else      { accV = addm(mmul(accV, bX, q, qi), bV, q); accX = mmul(accX, bX, q, qi); }
        }
        V = accV == 0 ? 0 : mmul(accV, 1, q, qi);
    }
    if (P.mode == 0) {
        res[gid] = V;
    } else if (V == 0) {
        uint k = atomic_fetch_add_explicit(nhits, 1, memory_order_relaxed);
        if (k < P.maxhits) hits[k] = n;
    }
}
