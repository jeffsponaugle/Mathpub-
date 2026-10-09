// concat.c -- fast search for the OEIS families "n divides the concatenation of 1..n".
//
// Family (tag, first A-number; base b = 2..25 is A-number first + b - 2):
//   Rn A029447  right concatenation 1 2 ... n, normal digits
//   Ln A029471  left  concatenation n ... 2 1, normal digits
//   Rk A029495  right, digits reversed, all zeros kept
//   Ld A029519  left,  digits reversed, least significant zeros of each k dropped
//   Rd A061931  right, digits reversed, least significant zeros dropped
//   Lk A061955  left,  digits reversed, all zeros kept
//
// Method: O(b log_b(n)^2) multiplications mod n per n, instead of O(n).
// A digit string is kept mod q as (V, X) = (value, b^length); concatenation is
// (V1,X1)(V2,X2) = (V1 X2 + V2, X1 X2).  For D-digit words, an aligned range
// k = P + m, 0 <= m < b^j, P == 0 mod b^j, has values val(P) + s*v_j(m), where v_j(m)
// is m's j-digit reversal (s = b^(D-j)) or m itself (s = 1, normal digits), and word
// lengths that depend only on m -- except for m = 0, whose trailing zeros come from P.
// The "interior" I_j (words m = 1 .. b^j-1 in concatenation order) is summarized by
//   A = sum b^position,  B = sum v_j(m) b^position,  X = b^length,  G = A + X.
// Splitting m on its top digit t, with L = b^(length of word t*b^j), mu = X L,
// P = sum_{t<b} mu^t, Q = sum_t t mu^t, S = sum_{t>=1} t mu^(b-1-t):
//   X' = mu^(b-1) X,   G' = P G,
//   reversed, left : B' = P (b B + X) + Q G - b X'
//   reversed, right: B' = b P B + S G
//   normal,   left : B' = P B + b^j (Q G + P X - b X')
//   normal,   right: B' = P B + b^j S G
// At each digit position the ranges with t >= 1 are summed in closed form from the same
// powers of mu, and the range with t = 0 (starting at P itself) is added separately.
// Odd moduli use Montgomery arithmetic; even n = 2^s q is tested mod 2^64 first, then mod q.
//
// Usage:
//   concat res    SEQ LO HI          print "n V(n)mod n" for n in [LO,HI)  (scalar path)
//   concat resl   SEQ LO HI          same, odd n through the 4-lane path
//   concat bench  SEQ LO COUNT       time the search loop on COUNT consecutive n
//   concat resv   SEQ LO HI          same, through the floating-point vector lanes
//   concat search SEQ LO HI [-t THREADS] [-c CHUNK] [-o DIR] [-x 0|1]
// SEQ is an A-number (e.g. A029519) or TAG:BASE (e.g. Ld:2).  -x 0 forces the 64-bit
// Montgomery engine; -x 1 (default where vectors exist) uses floating-point lanes.
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef unsigned long long ull;

#define AI static inline __attribute__((always_inline))
#define MAXL 4    // moduli evaluated together for instruction-level parallelism
#define MAXD 64   // max digits of n
#define MAXB 37   // max base + 1

typedef struct { int b, left, rev, drop; char tag[3]; int anum; } fam_t;

static const struct { int first; const char *tag; } FAMS[] = {
    {29447, "Rn"}, {29471, "Ln"}, {29495, "Rk"}, {29519, "Ld"}, {61931, "Rd"}, {61955, "Lk"},
};

// ------------------------------------------------------------------ arithmetic

enum { MONT, POW2 };                          // odd q (Montgomery, R = 2^64) or mod 2^64
typedef struct { u64 q, qi, one, bm; } ring_t;   // bm = b in ring form

AI u64 inv64(u64 n)                           // n^-1 mod 2^64, n odd
{
    u64 x = n;
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}

AI u64 mmul(u64 a, u64 b, u64 q, u64 qi)      // a b / 2^64 mod q, q odd < 2^63
{
    u128 t = (u128)a * b;
    u64 m = (u64)t * qi;
    u64 hi = (u64)(t >> 64), mh = (u64)(((u128)m * q) >> 64);
    u64 r = hi - mh;
    return hi < mh ? r + q : r;
}

AI u64 rmul(const int mode, const ring_t *R, u64 a, u64 b)
{
    return mode == POW2 ? a * b : mmul(a, b, R->q, R->qi);
}
AI u64 radd(const int mode, const ring_t *R, u64 a, u64 b)
{
    if (mode == POW2) return a + b;
    u64 s = a + b;
    return s >= R->q ? s - R->q : s;
}
AI u64 rsub(const int mode, const ring_t *R, u64 a, u64 b)
{
    if (mode == POW2) return a - b;
    return a >= b ? a - b : a - b + R->q;
}
AI u64 rsmall(const int mode, const ring_t *R, u64 a, unsigned c)   // c*a for small c
{
    u64 r = 0;
    for (; c; c >>= 1) {
        if (c & 1) r = radd(mode, R, r, a);
        a = radd(mode, R, a, a);
    }
    return r;
}
AI u64 rmulb(const int mode, const ring_t *R, int b, u64 a)        // b*a
{
    return b == 2 ? radd(mode, R, a, a) : rmul(mode, R, a, R->bm);
}

static void ring_init(const int mode, ring_t *R, u64 q, int b)
{
    if (mode == POW2) { R->q = R->qi = 0; R->one = 1; R->bm = (u64)b; return; }
    R->q = q; R->qi = inv64(q); R->one = (0 - q) % q;
    R->bm = rsmall(MONT, R, R->one, (unsigned)b);
}

static inline int ndig(u64 n, int b)
{
    if (b == 2) return n ? 64 - __builtin_clzll(n) : 0;
    int d = 0;
    for (; n; n /= b) d++;
    return d;
}

// ------------------------------------------------------------------ core evaluation

// Complete block of all Dp-digit words in base 2, reversed digits: (bV, bX) per lane.
// Zeros kept: all words have length Dp, so with y_j = 2^(Dp 2^j) the full sets F_j
// (m = 0 .. 2^j-1) satisfy G' = G (1+y), R' = 2 R (1+y) + (left ? y G : G), y' = y^2,
// and the block is F_{Dp-1} with offset 1 and scale 2.
// Zeros dropped: the interior recurrence of the header with b = 2 (P = 1+mu, Q = mu, S = 1);
// the block is I_{Dp-1} plus the one-digit word "1" for k = 2^(Dp-1).
AI void full_block2(const fam_t *F, const int mode, const int nl, const ring_t *Rv,
                    u64 (*pw)[MAXD + 1], const int Dp, u64 *bV, u64 *bX)
{
    const int left = F->left;
    if (!F->drop) {
        u64 y[MAXL], G[MAXL], R[MAXL];
        for (int l = 0; l < nl; l++) { y[l] = pw[l][Dp]; G[l] = Rv[l].one; R[l] = 0; }
        for (int j = 0; j < Dp - 1; j++)
            for (int l = 0; l < nl; l++) {
                const ring_t *Q = &Rv[l];
                u64 t = rmul(mode, Q, y[l], G[l]), u = rmul(mode, Q, y[l], R[l]);
                u64 r = radd(mode, Q, R[l], u);
                r = radd(mode, Q, r, r);
                R[l] = radd(mode, Q, r, left ? t : G[l]);
                G[l] = radd(mode, Q, G[l], t);
                y[l] = rmul(mode, Q, y[l], y[l]);
            }
        for (int l = 0; l < nl; l++) {
            bV[l] = radd(mode, &Rv[l], G[l], radd(mode, &Rv[l], R[l], R[l]));
            bX[l] = y[l];
        }
        return;
    }
    u64 G[MAXL], B[MAXL], X[MAXL], y[MAXL];
    for (int l = 0; l < nl; l++) { G[l] = X[l] = Rv[l].one; B[l] = 0; y[l] = pw[l][Dp - 1]; }
    for (int j = 0; j < Dp - 1; j++)
        for (int l = 0; l < nl; l++) {
            const ring_t *Q = &Rv[l];
            u64 mu = radd(mode, Q, y[l], y[l]);       // = X 2^(Dp-j) = 2^((Dp-1) 2^j + 1)
            y[l] = rmul(mode, Q, y[l], y[l]);
            u64 Gm = rmul(mode, Q, G[l], mu);
            u64 Xn = rmul(mode, Q, mu, X[l]);
            if (left) {
                u64 U = radd(mode, Q, radd(mode, Q, B[l], B[l]), X[l]);
                B[l] = rsub(mode, Q, radd(mode, Q, radd(mode, Q, U, rmul(mode, Q, U, mu)), Gm),
                            radd(mode, Q, Xn, Xn));
            } else {
                u64 Bs = radd(mode, Q, B[l], rmul(mode, Q, B[l], mu));
                B[l] = radd(mode, Q, radd(mode, Q, Bs, Bs), G[l]);
            }
            G[l] = radd(mode, Q, G[l], Gm);
            X[l] = Xn;
        }
    for (int l = 0; l < nl; l++) {
        const ring_t *Q = &Rv[l];
        if (left) {                           // (A + 2B) * 2 + 1, A = G - X
            u64 a = radd(mode, Q, rsub(mode, Q, G[l], X[l]), radd(mode, Q, B[l], B[l]));
            bV[l] = radd(mode, Q, radd(mode, Q, a, a), Q->one);
        } else {                              // X + A + 2B = G + 2B
            bV[l] = radd(mode, Q, G[l], radd(mode, Q, B[l], B[l]));
        }
        bX[l] = radd(mode, Q, X[l], X[l]);
    }
}

// Last, incomplete block (base 2, reversed digits): k = 2^(D-1) .. n with N = n+1 < 2^D.
// Each set bit j < D-1 of N contributes the range [Q_j, Q_j + 2^j), Q_j = N with bits
// <= j cleared; its words have values rev_D(Q_j) + 2^(D-j) rev_j(m).
AI void part_block2(const fam_t *F, const int mode, const int nl, const ring_t *Rv,
                    u64 (*pw)[MAXD + 1], const int D, const u64 *nv, u64 *bV, u64 *bX)
{
    const int left = F->left;
    u64 N[MAXL], cb[MAXL], G[MAXL], B[MAXL], X[MAXL], y[MAXL];
    for (int l = 0; l < nl; l++) {
        const ring_t *Q = &Rv[l];
        N[l] = nv[l] + 1;
        bV[l] = cb[l] = B[l] = 0;
        bX[l] = G[l] = Q->one;
        X[l] = F->drop ? Q->one : pw[l][D];  // keep: X holds y_j = 2^(D 2^j)
        y[l] = pw[l][D - 1];                 // drop: mu_j = 2 y_j, y_j = 2^((D-1) 2^j)
        for (u64 t = N[l]; t; t &= t - 1) cb[l] = radd(mode, Q, cb[l], pw[l][D - 1 - __builtin_ctzll(t)]);
    }
    for (int j = 0; j < D - 1; j++) {
        for (int l = 0; l < nl; l++) {
            if (!(N[l] >> j & 1)) continue;
            const ring_t *Q = &Rv[l];
            u64 V, Xs;
            cb[l] = rsub(mode, Q, cb[l], pw[l][D - 1 - j]);          // rev_D(Q_j)
            const u64 sB = rmul(mode, Q, pw[l][D - j], B[l]);
            if (!F->drop) {                                          // full set F_j
                V = radd(mode, Q, rmul(mode, Q, cb[l], G[l]), sB);
                Xs = X[l];
            } else {                                                 // I_j and the word Q_j
                const u64 L0 = pw[l][D - (j + 1 + __builtin_ctzll(N[l] >> (j + 1)))];
                Xs = rmul(mode, Q, X[l], L0);
                if (left)
                    V = radd(mode, Q, rmul(mode, Q, radd(mode, Q, rmul(mode, Q, rsub(mode, Q, G[l], X[l]), cb[l]), sB), L0),
                             cb[l]);
                else
                    V = radd(mode, Q, rmul(mode, Q, cb[l], G[l]), sB);
            }
            if (left) bV[l] = radd(mode, Q, rmul(mode, Q, bV[l], Xs), V);
            else bV[l] = radd(mode, Q, rmul(mode, Q, V, bX[l]), bV[l]);
            bX[l] = rmul(mode, Q, bX[l], Xs);
        }
        if (j == D - 2) break;
        for (int l = 0; l < nl; l++) {
            const ring_t *Q = &Rv[l];
            if (!F->drop) {
                u64 y = X[l], t = rmul(mode, Q, y, G[l]), u = rmul(mode, Q, y, B[l]);
                u64 r = radd(mode, Q, B[l], u);
                r = radd(mode, Q, r, r);
                B[l] = radd(mode, Q, r, left ? t : G[l]);
                G[l] = radd(mode, Q, G[l], t);
                X[l] = rmul(mode, Q, y, y);
            } else {
                u64 mu = radd(mode, Q, y[l], y[l]);
                y[l] = rmul(mode, Q, y[l], y[l]);
                u64 Gm = rmul(mode, Q, G[l], mu);
                u64 Xn = rmul(mode, Q, mu, X[l]);
                if (left) {
                    u64 U = radd(mode, Q, radd(mode, Q, B[l], B[l]), X[l]);
                    B[l] = rsub(mode, Q, radd(mode, Q, radd(mode, Q, U, rmul(mode, Q, U, mu)), Gm),
                                radd(mode, Q, Xn, Xn));
                } else {
                    u64 Bs = radd(mode, Q, B[l], rmul(mode, Q, B[l], mu));
                    B[l] = radd(mode, Q, radd(mode, Q, Bs, Bs), G[l]);
                }
                G[l] = radd(mode, Q, G[l], Gm);
                X[l] = Xn;
            }
        }
    }
}

// pm[l][0..top] = powers of mu[l], by a shallow product tree (depth about 2 log2 top).
AI void powers(const int mode, const int nl, const ring_t *Rv, const u64 *mu, const int top,
               u64 (*pm)[MAXB + 1])
{
    for (int l = 0; l < nl; l++) { pm[l][0] = Rv[l].one; pm[l][1] = mu[l]; }
    for (int t = 2; t <= top; t++)
        for (int l = 0; l < nl; l++)
            pm[l][t] = (t & 1) ? rmul(mode, &Rv[l], pm[l][t - 1], mu[l])
                               : rmul(mode, &Rv[l], pm[l][t >> 1], pm[l][t >> 1]);
}

// Sums over powers pm[0..T-1]:  *S0 = sum pm[u];  *S1 = sum (u+1) pm[u] (up = 1) or
// sum (T-u) pm[u] (up = 0).  Additions only.
AI void psums(const int mode, const ring_t *R, const u64 *pm, const int T, const int up, u64 *S0, u64 *S1)
{
    u64 run = 0, s = 0;
    if (up) for (int u = T - 1; u >= 0; u--) { run = radd(mode, R, run, pm[u]); s = radd(mode, R, s, run); }
    else    for (int u = 0; u < T; u++)      { run = radd(mode, R, run, pm[u]); s = radd(mode, R, s, run); }
    *S0 = run; *S1 = s;
}

// One step of the recurrences for general b, given pm = powers of y (keep) or mu (drop).
// P = sum_{t<b} pm[t], S = sum_{v<=b-2} (pm[0] + ... + pm[v]); then Q = (b-1) P - S
// for left concatenation and S itself plays the role of Q for right concatenation.
//   keep (full sets F_j): G' = P G,  R' = b P R + Q G
//   drop (interiors I_j): G' = P G,  X' = pm[b-1] X,
//        left B' = b (P B + P G - X') + P X - P G - S G,  right B' = b P B + S G
AI void stepb(const int mode, const ring_t *R, const int b, const int left, const int drop,
              const u64 *pm, u64 *G, u64 *B, u64 *X)
{
    u64 P = 0, S = 0;
    for (int t = 0; t < b; t++) {
        P = radd(mode, R, P, pm[t]);
        if (t < b - 1) S = radd(mode, R, S, P);
    }
    const u64 PG = rmul(mode, R, P, *G), PB = rmul(mode, R, P, *B), SG = rmul(mode, R, S, *G);
    if (!drop) {
        *B = left ? rsub(mode, R, rsub(mode, R, rmulb(mode, R, b, radd(mode, R, PB, PG)), PG), SG)
                  : radd(mode, R, rmulb(mode, R, b, PB), SG);
    } else {
        const u64 Xn = rmul(mode, R, pm[b - 1], *X);
        if (left) {
            const u64 PX = rmul(mode, R, P, *X);
            *B = rsub(mode, R, rsub(mode, R, radd(mode, R, rmulb(mode, R, b, rsub(mode, R, radd(mode, R, PB, PG), Xn)), PX), PG), SG);
        } else {
            *B = radd(mode, R, rmulb(mode, R, b, PB), SG);
        }
        *X = Xn;
    }
    *G = PG;
}

// Complete block of all Dp-digit words, reversed digits, any base b >= 2: (bV, bX).
// Keep: full sets F_j over y_j = b^(Dp b^j); the block is F_{Dp-1} at offsets t = 1..b-1
// (value t + b rev(m)), so V = G S1 + b R S0 over powers of Y = y_{Dp-1}.
// Drop: interiors I_j with mu_j = X_j b^(Dp-j); the block adds the words t*b^(Dp-1) = "t".
AI void full_blockb(const fam_t *F, const int mode, const int nl, const ring_t *Rv,
                    u64 (*pw)[MAXD + 1], const int Dp, u64 *bV, u64 *bX)
{
    const int b = F->b, left = F->left, drop = F->drop;
    u64 G[MAXL], B[MAXL], X[MAXL], mu[MAXL], pm[MAXL][MAXB + 1];
    for (int l = 0; l < nl; l++) { G[l] = X[l] = Rv[l].one; B[l] = 0; mu[l] = pw[l][Dp]; }
    for (int j = 0; j < Dp - 1; j++) {
        if (drop)
            for (int l = 0; l < nl; l++) mu[l] = rmul(mode, &Rv[l], X[l], pw[l][Dp - j]);
        powers(mode, nl, Rv, mu, drop ? b - 1 : b, pm);
        for (int l = 0; l < nl; l++) {
            stepb(mode, &Rv[l], b, left, drop, pm[l], &G[l], &B[l], &X[l]);
            if (!drop) mu[l] = pm[l][b];          // y_{j+1} = y_j^b
        }
    }
    if (drop)
        for (int l = 0; l < nl; l++) mu[l] = rmul(mode, &Rv[l], X[l], pw[l][1]);
    powers(mode, nl, Rv, mu, b - 1, pm);
    for (int l = 0; l < nl; l++) {
        const ring_t *R = &Rv[l];
        u64 S0, S1;
        psums(mode, R, pm[l], b - 1, left, &S0, &S1);
        if (!drop) {
            bV[l] = radd(mode, R, rmul(mode, R, G[l], S1), rmulb(mode, R, b, rmul(mode, R, B[l], S0)));
        } else if (left) {                    // alpha S1 + b (B L) S0, alpha = (G - X) L + 1, L = b
            u64 al = radd(mode, R, rmulb(mode, R, b, rsub(mode, R, G[l], X[l])), R->one);
            bV[l] = radd(mode, R, rmul(mode, R, al, S1), rmul(mode, R, rmulb(mode, R, b, rmulb(mode, R, b, B[l])), S0));
        } else {                              // G S1 + b B S0
            bV[l] = radd(mode, R, rmul(mode, R, G[l], S1), rmulb(mode, R, b, rmul(mode, R, B[l], S0)));
        }
        bX[l] = pm[l][b - 1];
    }
}

// Last, incomplete block (reversed digits, any b): k = b^(D-1) .. n, N = n+1 < b^D with
// digits e[i].  At position j the ranges t = (j == D-1 ? 1 : 0) .. e_j - 1 start at
// Q_j + t b^j, Q_j = N with digits <= j cleared; word values rev_D(Q_j) + t b^(D-1-j) + s rev_j(m)
// with s = b^(D-j).  Keep: each range is a full set F_j.  Drop: ranges t >= 1 are I_j plus
// a word of length D-j; the range t = 0 is I_j plus the word Q_j of length D - v_b(Q_j).
AI void part_blockb(const fam_t *F, const int mode, const int nl, const ring_t *Rv,
                    u64 (*pw)[MAXD + 1], const int D, int (*e)[MAXD], int (*vz)[MAXD],
                    u64 *bV, u64 *bX)
{
    const int b = F->b, left = F->left, drop = F->drop;
    u64 G[MAXL], B[MAXL], X[MAXL], mu[MAXL], cb[MAXL], pm[MAXL][MAXB + 1];
    for (int l = 0; l < nl; l++) {
        const ring_t *R = &Rv[l];
        G[l] = X[l] = bX[l] = R->one; B[l] = bV[l] = cb[l] = 0;
        mu[l] = pw[l][D];
        for (int i = 0; i < D; i++)
            if (e[l][i]) cb[l] = radd(mode, R, cb[l], rsmall(mode, R, pw[l][D - 1 - i], e[l][i]));
    }
    for (int j = 0; j < D; j++) {
        const int top = j == D - 1;
        if (drop)
            for (int l = 0; l < nl; l++) mu[l] = rmul(mode, &Rv[l], X[l], pw[l][D - j]);
        powers(mode, nl, Rv, mu, (drop || top) ? b - 1 : b, pm);
        for (int l = 0; l < nl; l++) {
            const int ej = e[l][j];
            if (!ej) continue;
            const ring_t *R = &Rv[l];
            const u64 om = pw[l][D - 1 - j], sB = rmul(mode, R, pw[l][D - j], B[l]);
            cb[l] = rsub(mode, R, cb[l], rsmall(mode, R, om, ej));       // rev_D(Q_j)
            if (!drop) {                      // ranges u = 0..T-1, t = t0 + u, each F_j, X = y
                const int t0 = top, T = ej - t0;
                if (T <= 0) continue;
                u64 S0, S1;                   // S1 - S0 = sum u y^u (left), sum u y^(T-1-u) (right)
                psums(mode, R, pm[l], T, left, &S0, &S1);
                const u64 S1u = rsub(mode, R, S1, S0);
                const u64 c = radd(mode, R, cb[l], t0 ? om : 0);
                const u64 V = radd(mode, R, rmul(mode, R, radd(mode, R, rmul(mode, R, c, G[l]), sB), S0),
                                   rmul(mode, R, om, rmul(mode, R, G[l], S1u)));
                const u64 XT = pm[l][T];
                if (left) bV[l] = radd(mode, R, rmul(mode, R, bV[l], XT), V);
                else bV[l] = radd(mode, R, rmul(mode, R, V, bX[l]), bV[l]);
                bX[l] = rmul(mode, R, bX[l], XT);
                continue;
            }
            const u64 L = pw[l][D - j], A = rsub(mode, R, G[l], X[l]);
            const int T = ej - 1;             // ranges t = 1 .. ej-1
            if (T > 0) {
                u64 S0, S1;
                psums(mode, R, pm[l], T, left, &S0, &S1);
                const u64 XT = pm[l][T];
                if (left) {
                    u64 al = radd(mode, R, rmul(mode, R, A, L), R->one);
                    u64 c0 = radd(mode, R, rmul(mode, R, cb[l], al), rmul(mode, R, sB, L));
                    u64 VT = radd(mode, R, rmul(mode, R, c0, S0), rmul(mode, R, om, rmul(mode, R, al, S1)));
                    bV[l] = radd(mode, R, rmul(mode, R, bV[l], XT), VT);
                    bX[l] = rmul(mode, R, bX[l], XT);
                } else {
                    u64 c0 = radd(mode, R, rmul(mode, R, cb[l], G[l]), sB);
                    u64 VT = radd(mode, R, rmul(mode, R, c0, S0), rmul(mode, R, om, rmul(mode, R, G[l], S1)));
                    bV[l] = radd(mode, R, rmul(mode, R, VT, bX[l]), bV[l]);
                    bX[l] = rmul(mode, R, XT, bX[l]);
                }
            }
            if (!top) {                       // t = 0: the range starts with the word Q_j
                const u64 L0 = pw[l][D - vz[l][j]], X0 = rmul(mode, R, X[l], L0);
                if (left) {
                    u64 V0 = radd(mode, R, rmul(mode, R, radd(mode, R, rmul(mode, R, A, cb[l]), sB), L0), cb[l]);
                    bV[l] = radd(mode, R, rmul(mode, R, bV[l], X0), V0);
                    bX[l] = rmul(mode, R, bX[l], X0);
                } else {
                    u64 V0 = radd(mode, R, rmul(mode, R, cb[l], G[l]), sB);
                    bV[l] = radd(mode, R, rmul(mode, R, V0, bX[l]), bV[l]);
                    bX[l] = rmul(mode, R, X0, bX[l]);
                }
            }
        }
        if (top) break;
        for (int l = 0; l < nl; l++) {
            stepb(mode, &Rv[l], b, left, drop, pm[l], &G[l], &B[l], &X[l]);
            if (!drop) mu[l] = pm[l][b];
        }
    }
}

// out[l] = V(nv[l]) in ring Rv[l] (Montgomery form for MONT).  All nv[l] must have
// the same number of base-b digits.  nv[l] < 2^58.
AI void core(const fam_t *F, const int mode, const int nl, const u64 *nv, const ring_t *Rv,
             u64 *out)
{
    const int b = F->b, left = F->left, rev = F->rev, drop = F->drop;
    const int D = ndig(nv[0], b);
    int e[MAXL][MAXD], vz[MAXL][MAXD], full[MAXL];
    u64 pw[MAXL][MAXD + 1], accV[MAXL], accX[MAXL];

    for (int l = 0; l < nl; l++) {
        const ring_t *R = &Rv[l];
        u64 t = nv[l] + 1;                    // digits of N = n+1 bound the last block
        if (b == 2) {
            for (int i = 0; i < D; i++) e[l][i] = (int)(t >> i & 1);
            t >>= D;
        } else
            for (int i = 0; i < D; i++) { e[l][i] = (int)(t % b); t /= b; }
        full[l] = t != 0;                     // N == b^D: last block is complete
        for (int i = D - 1, z = D; i >= 0; i--) { vz[l][i] = z; if (e[l][i]) z = i; }
        pw[l][0] = R->one;                    // pw[i] = b^i
        for (int i = 1; i <= D; i++) pw[l][i] = rmulb(mode, R, b, pw[l][i - 1]);
        accV[l] = 0; accX[l] = R->one;
    }

    for (int Dp = 1; Dp <= D; Dp++) {         // block of Dp-digit words
        u64 G[MAXL], B[MAXL], X[MAXL], bV[MAXL], bX[MAXL], cb[MAXL];
        int part[MAXL];
        if (b == 2 && rev && Dp < D) {
            full_block2(F, mode, nl, Rv, pw, Dp, bV, bX);
            goto join;
        }
        if (b != 2 && rev && Dp < D) {
            full_blockb(F, mode, nl, Rv, pw, Dp, bV, bX);
            goto join;
        }
        if (rev && D >= 2) {
            int anyfull = 0;
            for (int l = 0; l < nl; l++) anyfull |= full[l];
            if (!anyfull) {
                if (b == 2) part_block2(F, mode, nl, Rv, pw, D, nv, bV, bX);
                else part_blockb(F, mode, nl, Rv, pw, D, e, vz, bV, bX);
                goto join;
            }
        }
        for (int l = 0; l < nl; l++) {
            const ring_t *R = &Rv[l];
            part[l] = Dp == D && !full[l];
            G[l] = X[l] = bX[l] = R->one;     // I_0 is empty
            B[l] = bV[l] = cb[l] = 0;
            if (part[l])                      // cb = val(N), reduced digit by digit below
                for (int i = 0; i < D; i++)
                    if (e[l][i])
                        cb[l] = radd(mode, R, cb[l],
                                     rsmall(mode, R, rev ? pw[l][D - 1 - i] : pw[l][i], e[l][i]));
        }
        for (int j = 0; j < Dp; j++) {
            u64 pm[MAXL][MAXB];               // powers of mu
            const int top = j == Dp - 1;
            for (int l = 0; l < nl; l++) {
                const ring_t *R = &Rv[l];
                const u64 L = drop ? pw[l][Dp - j] : pw[l][Dp];
                pm[l][0] = R->one;
                pm[l][1] = rmul(mode, R, X[l], L);
                for (int t = 2; t < b; t++) pm[l][t] = rmul(mode, R, pm[l][t - 1], pm[l][1]);
            }
            // ranges [Q + t b^j, Q + (t+1) b^j) at digit position j
            for (int l = 0; l < nl; l++) {
                if (!top && !part[l]) continue;
                const ring_t *R = &Rv[l];
                const int ej = (top && !part[l]) ? b : e[l][j];
                if (ej == 0) continue;
                const u64 om = rev ? pw[l][Dp - 1 - j] : pw[l][j];   // weight of digit j
                if (part[l]) cb[l] = rsub(mode, R, cb[l], rsmall(mode, R, om, ej));  // val(Q)
                const u64 L = drop ? pw[l][Dp - j] : pw[l][Dp];
                const u64 A = rsub(mode, R, G[l], X[l]);
                const u64 sB = rev ? rmul(mode, R, pw[l][Dp - j], B[l]) : B[l];
                const int T = ej - 1;         // ranges with t = 1 .. ej-1
                if (T > 0) {
                    u64 run = 0, S1 = 0;
                    if (left)
                        for (int u = T - 1; u >= 0; u--) {
                            run = radd(mode, R, run, pm[l][u]);
                            S1 = radd(mode, R, S1, run);   // sum (u+1) mu^u
                        }
                    else
                        for (int u = 0; u < T; u++) {
                            run = radd(mode, R, run, pm[l][u]);
                            S1 = radd(mode, R, S1, run);   // sum (T-u) mu^u
                        }
                    const u64 S0 = run, XT = pm[l][T];
                    if (left) {
                        u64 al = radd(mode, R, rmul(mode, R, A, L), R->one);
                        u64 c0 = radd(mode, R, rmul(mode, R, cb[l], al), rmul(mode, R, sB, L));
                        u64 VT = radd(mode, R, rmul(mode, R, c0, S0),
                                      rmul(mode, R, om, rmul(mode, R, al, S1)));
                        bV[l] = radd(mode, R, rmul(mode, R, bV[l], XT), VT);
                        bX[l] = rmul(mode, R, bX[l], XT);
                    } else {
                        u64 c0 = radd(mode, R, rmul(mode, R, cb[l], G[l]), sB);
                        u64 VT = radd(mode, R, rmul(mode, R, c0, S0),
                                      rmul(mode, R, om, rmul(mode, R, G[l], S1)));
                        bV[l] = radd(mode, R, rmul(mode, R, VT, bX[l]), bV[l]);
                        bX[l] = rmul(mode, R, XT, bX[l]);
                    }
                }
                if (!top) {                   // t = 0: range starts at Q itself
                    const u64 L0 = drop ? pw[l][Dp - vz[l][j]] : pw[l][Dp];
                    const u64 X0 = rmul(mode, R, X[l], L0);
                    if (left) {
                        u64 V0 = radd(mode, R, rmul(mode, R, radd(mode, R, rmul(mode, R, A, cb[l]), sB), L0),
                                      cb[l]);
                        bV[l] = radd(mode, R, rmul(mode, R, bV[l], X0), V0);
                        bX[l] = rmul(mode, R, bX[l], X0);
                    } else {
                        u64 V0 = radd(mode, R, rmul(mode, R, cb[l], G[l]), sB);
                        bV[l] = radd(mode, R, rmul(mode, R, V0, bX[l]), bV[l]);
                        bX[l] = rmul(mode, R, X0, bX[l]);
                    }
                }
            }
            if (top) break;
            for (int l = 0; l < nl; l++) {    // I_j -> I_{j+1}
                const ring_t *R = &Rv[l];
                const u64 mu = pm[l][1];
                u64 Xn, Gn, Bn;
                if (rev && b == 2) {
                    const u64 Gm = rmul(mode, R, G[l], mu);
                    Xn = rmul(mode, R, mu, X[l]);
                    Gn = radd(mode, R, G[l], Gm);
                    if (left) {
                        u64 U = radd(mode, R, radd(mode, R, B[l], B[l]), X[l]);
                        Bn = rsub(mode, R, radd(mode, R, radd(mode, R, U, rmul(mode, R, U, mu)), Gm),
                                  radd(mode, R, Xn, Xn));
                    } else {
                        u64 Bs = radd(mode, R, B[l], rmul(mode, R, B[l], mu));
                        Bn = radd(mode, R, radd(mode, R, Bs, Bs), G[l]);
                    }
                } else {
                    u64 P = 0, W = 0, run = 0;    // W = Q (left) or S (right)
                    if (left)
                        for (int t = b - 1; t >= 1; t--) {
                            run = radd(mode, R, run, pm[l][t]);
                            W = radd(mode, R, W, run);
                        }
                    else
                        for (int t = 0; t <= b - 2; t++) {
                            run = radd(mode, R, run, pm[l][t]);
                            W = radd(mode, R, W, run);
                        }
                    for (int t = 0; t < b; t++) P = radd(mode, R, P, pm[l][t]);
                    Xn = rmul(mode, R, pm[l][b - 1], X[l]);
                    Gn = rmul(mode, R, P, G[l]);
                    const u64 WG = rmul(mode, R, W, G[l]);
                    if (rev && left)
                        Bn = rsub(mode, R,
                                  radd(mode, R, rmul(mode, R, P, radd(mode, R, rmulb(mode, R, b, B[l]), X[l])), WG),
                                  rmulb(mode, R, b, Xn));
                    else if (rev)
                        Bn = radd(mode, R, rmulb(mode, R, b, rmul(mode, R, P, B[l])), WG);
                    else if (left)
                        Bn = radd(mode, R, rmul(mode, R, P, B[l]),
                                  rmul(mode, R, pw[l][j],
                                       rsub(mode, R, radd(mode, R, WG, rmul(mode, R, P, X[l])),
                                            rmulb(mode, R, b, Xn))));
                    else
                        Bn = radd(mode, R, rmul(mode, R, P, B[l]), rmul(mode, R, pw[l][j], WG));
                }
                X[l] = Xn; G[l] = Gn; B[l] = Bn;
            }
        }
    join:
        for (int l = 0; l < nl; l++) {        // join the block to the words so far
            const ring_t *R = &Rv[l];
            if (left) {
                accV[l] = radd(mode, R, rmul(mode, R, bV[l], accX[l]), accV[l]);
                accX[l] = rmul(mode, R, bX[l], accX[l]);
            } else {
                accV[l] = radd(mode, R, rmul(mode, R, accV[l], bX[l]), bV[l]);
                accX[l] = rmul(mode, R, accX[l], bX[l]);
            }
        }
    }
    for (int l = 0; l < nl; l++) out[l] = accV[l];
}

static void core_mont4(const fam_t *F, const u64 *n, const ring_t *R, u64 *o) { core(F, MONT, 4, n, R, o); }
static void core_mont1(const fam_t *F, const u64 *n, const ring_t *R, u64 *o) { core(F, MONT, 1, n, R, o); }
static void core_pow2(const fam_t *F, const u64 *n, const ring_t *R, u64 *o) { core(F, POW2, 1, n, R, o); }

// V(n) mod n, any n >= 1 (verification path: CRT of the 2-part and the odd part)
static u64 full_residue(const fam_t *F, u64 n)
{
    if (n == 1) return 0;
    int s = __builtin_ctzll(n);
    u64 q = n >> s, r2 = 0, rq = 0, v;
    ring_t R;
    if (s) { ring_init(POW2, &R, 0, F->b); core_pow2(F, &n, &R, &v); r2 = v & ((1ull << s) - 1); }
    if (q > 1) { ring_init(MONT, &R, q, F->b); core_mont1(F, &n, &R, &v); rq = mmul(v, 1, q, R.qi); }
    if (!s) return rq;
    u64 k = ((r2 - rq) * inv64(q)) & ((1ull << s) - 1);
    return rq + q * k;
}

// ------------------------------------------------------------------ floating-point lanes
// For n < 2^43 (reversed-digit families) residues are kept as exact integers in doubles,
// signed, |x| <~ n/2 after reduction.  a*b mod n:  h = a*b, l = fma(a, b, -h) is the exact
// low part, q = round(h/n), r = fma(-q, n, h) + l.  For |a|, |b| <= 4n the quotient
// estimate is within 1/2 + 16 n 2^-52 of h/n, so |r| <= n/2 + 24 n^2 2^-52 < 0.52 n and every
// intermediate is an exact integer below 2^53.  Sums are re-reduced before they can exceed
// 4n as a multiplicand.  A vector holds VW different moduli: AVX-512 8, NEON 2, otherwise 1.
// round(x) for |x| < 2^51 is (x + 1.5*2^52) - 1.5*2^52 in round-to-nearest mode.
#if defined(__AVX512F__) && !defined(NO_FPV)
#include <immintrin.h>
#define VW 8
typedef __m512d vd;
#define VFMS(a, b, c) _mm512_fmsub_pd(a, b, c)            // a*b - c
#define VFNM(a, b, c) _mm512_fnmadd_pd(a, b, c)           // c - a*b
#define VBC(x) _mm512_set1_pd(x)
#define VSEL(m, a, b) _mm512_mask_blend_pd((__mmask8)(m), b, a)   // bit l of m ? a : b
#elif defined(__aarch64__) && !defined(NO_FPV)
#include <arm_neon.h>
#define VW 2
typedef float64x2_t vd;
#define VFMS(a, b, c) vfmaq_f64(vnegq_f64(c), a, b)
#define VFNM(a, b, c) vfmsq_f64(c, a, b)
#define VBC(x) vdupq_n_f64(x)
static inline vd vsel2(unsigned m, vd a, vd b)
{
    uint64x2_t k = {(m & 1) ? ~0ull : 0, (m & 2) ? ~0ull : 0};
    return vbslq_f64(k, a, b);
}
#define VSEL(m, a, b) vsel2(m, a, b)
#else
#include <math.h>
#define VW 1
typedef double vd;
#define VFMS(a, b, c) fma(a, b, -(c))
#define VFNM(a, b, c) fma(-(a), b, c)
#define VBC(x) ((double)(x))
#define VSEL(m, a, b) (((m) & 1) ? (a) : (b))
#endif

#define FPV_MAXN (1ull << 43)
typedef union { vd v; double a[VW]; } vu;
typedef struct { vd n, ni; } vring_t;

AI vd vround(vd x) { return (x + VBC(6755399441055744.0)) - VBC(6755399441055744.0); }
AI vd vmm(const vring_t *R, vd a, vd b)       // a b mod n, |a|,|b| <= 4n -> |r| < 0.52 n
{
    vd h = a * b, l = VFMS(a, b, h);
    return VFNM(vround(h * R->ni), R->n, h) + l;
}
AI vd vnr(const vring_t *R, vd x) { return VFNM(vround(x * R->ni), R->n, x); }   // |x| < 2^52

typedef struct {
    vring_t R;
    int b, D;
    vu pw[MAXD + 1];                           // b^i mod n
    int e[VW][MAXD], vz[VW][MAXD], full[VW];   // digits of N = n+1 (full: N = b^D)
} vctx_t;

// powers pm[0..top] of mu
AI void vpowers(const vring_t *R, vd mu, int top, vu *pm)
{
    pm[0].v = VBC(1.0);
    pm[1].v = mu;
    for (int t = 2; t <= top; t++) pm[t].v = (t & 1) ? vmm(R, pm[t - 1].v, mu) : vmm(R, pm[t >> 1].v, pm[t >> 1].v);
}

// One step I_j -> I_{j+1} (drop) or F_j -> F_{j+1} (keep), as stepb() in the integer code.
AI void vstep(const vring_t *R, int b, int left, int drop, const vu *pm, vd *G, vd *B, vd *X)
{
    vd P = VBC(0.0), S = VBC(0.0);
    for (int t = 0; t < b; t++) {
        P = P + pm[t].v;
        if (t < b - 1) S = S + P;
    }
    P = vnr(R, P); S = vnr(R, S);
    const vd PG = vmm(R, P, *G), PB = vmm(R, P, *B), SG = vmm(R, S, *G), bb = VBC(b);
    if (!drop) {
        *B = left ? vnr(R, bb * (PB + PG) - PG - SG) : vnr(R, bb * PB + SG);
    } else {
        const vd Xn = vmm(R, pm[b - 1].v, *X);
        if (left) *B = vnr(R, bb * (PB + PG - Xn) + vmm(R, P, *X) - PG - SG);
        else *B = vnr(R, bb * PB + SG);
        *X = Xn;
    }
    *G = PG;
}

// Complete block of Dp-digit words (all lanes), any base: as full_blockb / full_block2.
static void vfull(const fam_t *F, const vctx_t *C, int Dp, vd *bV, vd *bX)
{
    const vring_t *R = &C->R;
    const int b = F->b, left = F->left, drop = F->drop;
    const vd one = VBC(1.0);
    if (b == 2 && !drop) {                     // full sets, 3 products per step
        vd y = C->pw[Dp].v, G = one, Rr = VBC(0.0);
        for (int j = 0; j < Dp - 1; j++) {
            vd t = vmm(R, y, G), u = vmm(R, y, Rr);
            Rr = vnr(R, (Rr + u) * VBC(2.0) + (left ? t : G));
            G = vnr(R, G + t);
            y = vmm(R, y, y);
        }
        *bV = vnr(R, G + Rr + Rr);
        *bX = y;
        return;
    }
    if (b == 2) {                              // interiors, mu = 2 y, y = 2^((Dp-1) 2^j)
        vd y = C->pw[Dp - 1].v, G = one, X = one, B = VBC(0.0);
        for (int j = 0; j < Dp - 1; j++) {
            vd mu = y + y;
            y = vmm(R, y, y);
            vd Gm = vmm(R, G, mu), Xn = vmm(R, mu, X);
            if (left) {
                vd U = B + B + X;
                B = vnr(R, U + vmm(R, U, mu) + Gm - Xn - Xn);
            } else {
                B = vnr(R, (B + vmm(R, B, mu)) * VBC(2.0) + G);
            }
            G = vnr(R, G + Gm);
            X = Xn;
        }
        *bV = left ? vnr(R, (G - X + B + B) * VBC(2.0) + one) : vnr(R, G + B + B);
        *bX = vnr(R, X + X);
        return;
    }
    vu pm[MAXB + 1];
    vd G = one, B = VBC(0.0), X = one, mu = C->pw[Dp].v;
    for (int j = 0; j < Dp - 1; j++) {
        if (drop) mu = vmm(R, X, C->pw[Dp - j].v);
        vpowers(R, mu, drop ? b - 1 : b, pm);
        vstep(R, b, left, drop, pm, &G, &B, &X);
        if (!drop) mu = pm[b].v;
    }
    if (drop) mu = vmm(R, X, C->pw[1].v);
    vpowers(R, mu, b - 1, pm);
    vd S0 = VBC(0.0), S1 = VBC(0.0);           // u = 0 .. b-2
    for (int u = 0; u < b - 1; u++) {
        S0 = S0 + pm[u].v;
        S1 = S1 + VBC(left ? u + 1 : b - 1 - u) * pm[u].v;
    }
    S0 = vnr(R, S0); S1 = vnr(R, S1);
    const vd bb = VBC(b);
    if (!drop) *bV = vnr(R, vmm(R, G, S1) + bb * vmm(R, B, S0));
    else if (left) *bV = vnr(R, vmm(R, vnr(R, bb * (G - X) + one), S1) + vmm(R, vnr(R, VBC(b * b) * B), S0));
    else *bV = vnr(R, vmm(R, G, S1) + bb * vmm(R, B, S0));
    *bX = pm[b - 1].v;
}

// Last block k = b^(D-1) .. n with per-lane digits of N (lanes with N = b^D carry digit b
// at position D-1), as part_blockb / part_block2 with per-lane masks.
static void vpart(const fam_t *F, const vctx_t *C, vd *bVo, vd *bXo)
{
    const vring_t *R = &C->R;
    const int b = F->b, left = F->left, drop = F->drop, D = C->D;
    const vd one = VBC(1.0);
    vu pm[MAXB + 1];
    vd cb = VBC(0.0);
    for (int i = 0; i < D; i++) {               // cb = rev_D(N) mod n
        double ed[VW];
        for (int l = 0; l < VW; l++) ed[l] = C->e[l][i];
        vu t; memcpy(t.a, ed, sizeof ed);
        cb = vnr(R, cb + t.v * C->pw[D - 1 - i].v);
    }
    vd G = one, B = VBC(0.0), X = one, mu = C->pw[D].v, bV = VBC(0.0), bX = one;
    for (int j = 0; j < D; j++) {
        const int top = j == D - 1;
        if (drop) mu = vmm(R, X, C->pw[D - j].v);
        vpowers(R, mu, (drop || top) ? b - 1 : b, pm);
        unsigned mact = 0;
        vu ev;
        for (int l = 0; l < VW; l++) { ev.a[l] = C->e[l][j]; if (C->e[l][j]) mact |= 1u << l; }
        if (mact) {
            const vd om = C->pw[D - 1 - j].v, sB = vmm(R, C->pw[D - j].v, B);
            cb = VSEL(mact, vnr(R, cb - ev.v * om), cb);           // rev_D(Q_j)
            const int t0 = drop ? 1 : top;     // ranges t = t0 .. e_j-1 summed in closed form
            int T[VW];
            unsigned mt = 0;
            for (int l = 0; l < VW; l++) { T[l] = C->e[l][j] - t0; if (T[l] > 0) mt |= 1u << l; }
            if (mt) {
                // S0 = sum_{u<T} mu^u;  S1 = keep: left sum u mu^u, right sum (T-1-u) mu^u
                //                             drop: left sum (u+1) mu^u, right sum (T-u) mu^u
                vd S0 = VBC(0.0), S1 = VBC(0.0);
                for (int u = 0; u < b; u++) {
                    unsigned m = 0;
                    vu cu;
                    for (int l = 0; l < VW; l++) {
                        if (u < T[l]) m |= 1u << l;
                        cu.a[l] = left ? (drop ? u + 1 : u) : (drop ? T[l] - u : T[l] - 1 - u);
                    }
                    if (!m) break;
                    S0 = VSEL(m, S0 + pm[u].v, S0);
                    S1 = VSEL(m, S1 + cu.v * pm[u].v, S1);
                }
                S0 = vnr(R, S0); S1 = vnr(R, S1);
                vu XT;
                for (int l = 0; l < VW; l++) XT.a[l] = T[l] > 0 ? pm[T[l]].a[l] : 1.0;
                vd V;
                if (!drop) {
                    vd c = cb;
                    if (top) c = cb + om;
                    V = vnr(R, vmm(R, vnr(R, vmm(R, c, G) + sB), S0) + vmm(R, om, vmm(R, G, S1)));
                } else if (left) {
                    const vd L = C->pw[D - j].v;
                    const vd al = vnr(R, vmm(R, G - X, L) + one);
                    const vd c0 = vnr(R, vmm(R, cb, al) + vmm(R, sB, L));
                    V = vnr(R, vmm(R, c0, S0) + vmm(R, om, vmm(R, al, S1)));
                } else {
                    const vd c0 = vnr(R, vmm(R, cb, G) + sB);
                    V = vnr(R, vmm(R, c0, S0) + vmm(R, om, vmm(R, G, S1)));
                }
                const vd nV = left ? vnr(R, vmm(R, bV, XT.v) + V) : vnr(R, vmm(R, V, bX) + bV);
                bX = VSEL(mt, vmm(R, bX, XT.v), bX);
                bV = VSEL(mt, nV, bV);
            }
            if (drop && !top) {                // t = 0: the range starts with the word Q_j
                vu L0;
                for (int l = 0; l < VW; l++) L0.a[l] = C->pw[D - C->vz[l][j]].a[l];
                const vd X0 = vmm(R, X, L0.v);
                vd nV;
                if (left) {
                    const vd V0 = vnr(R, vmm(R, vnr(R, vmm(R, G - X, cb) + sB), L0.v) + cb);
                    nV = vnr(R, vmm(R, bV, X0) + V0);
                } else {
                    const vd V0 = vnr(R, vmm(R, cb, G) + sB);
                    nV = vnr(R, vmm(R, V0, bX) + bV);
                }
                bX = VSEL(mact, vmm(R, bX, X0), bX);
                bV = VSEL(mact, nV, bV);
            }
        }
        if (top) break;
        if (b == 2 && !drop) {                 // vstep() specialised: P = 1 + y, Q = y, S = 1
            vd t = vmm(R, mu, G), u = vmm(R, mu, B);
            B = vnr(R, (B + u) * VBC(2.0) + (left ? t : G));
            G = vnr(R, G + t);
            mu = pm[2].v;
        } else {
            vstep(R, b, left, drop, pm, &G, &B, &X);
            if (!drop) mu = pm[b].v;
        }
    }
    *bVo = bV; *bXo = bX;
}

// V(n) mod n for VW moduli with equal digit counts (reversed-digit families, n < 2^43);
// returns residues in [0, n).
// dg[l] (optional) = base-b digits of n_l + 1, least significant first, positions 0..D
// (position D is nonzero exactly when n_l + 1 = b^D); without it they are computed here.
static void core_fpv(const fam_t *F, const u64 *nv, int D, const unsigned char *const *dg, u64 *res)
{
    vctx_t C;
    const int b = F->b, left = F->left;
    vu nn, ni;
    for (int l = 0; l < VW; l++) { nn.a[l] = (double)nv[l]; ni.a[l] = 1.0 / (double)nv[l]; }
    C.R.n = nn.v; C.R.ni = ni.v;
    const vring_t *R = &C.R;
    C.b = b;
    C.D = D;
    C.pw[0].v = VBC(1.0);
    for (int i = 1; i <= D; i++) C.pw[i].v = vnr(R, C.pw[i - 1].v * VBC(b));
    int anyfull = 0;
    for (int l = 0; l < VW; l++) {
        if (dg) {
            for (int i = 0; i < D; i++) C.e[l][i] = dg[l][i];
            C.full[l] = dg[l][D] != 0;
        } else {
            u64 t = nv[l] + 1;
            for (int i = 0; i < D; i++) { C.e[l][i] = (int)(t % b); t /= b; }
            C.full[l] = t != 0;
        }
        if (C.full[l]) { C.e[l][D - 1] = b; anyfull = 1; }     // whole top digit range
        for (int i = D - 1, z = D; i >= 0; i--) { C.vz[l][i] = z; if (C.e[l][i]) z = i; }
    }
    vd accV = VBC(0.0), accX = VBC(1.0), bV, bX;
    for (int Dp = 1; Dp <= D; Dp++) {
        if (Dp < D) vfull(F, &C, Dp, &bV, &bX);
        else if (anyfull && b == 2) {
            int all = 1;
            for (int l = 0; l < VW; l++) all &= C.full[l];
            if (all) vfull(F, &C, Dp, &bV, &bX);
            else vpart(F, &C, &bV, &bX);
        } else vpart(F, &C, &bV, &bX);
        if (left) { accV = vnr(R, vmm(R, bV, accX) + accV); accX = vmm(R, bX, accX); }
        else      { accV = vnr(R, vmm(R, accV, bX) + bV); accX = vmm(R, accX, bX); }
    }
    vu r; r.v = accV;
    for (int l = 0; l < VW; l++) {
        double x = r.a[l];
        if (x < 0) x += nn.a[l];
        res[l] = nv[l] == 1 ? 0 : (u64)x;
    }
}

// ------------------------------------------------------------------ search

static u64 gcd64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

// Necessary conditions read off the last word.  Left: V == 1 (mod b).  Right, reversed:
// V == rev(n) (mod b^len(n)), so gcd(n, b^len(n)) must divide rev(n).
static int prefilter(const fam_t *F, u64 n)
{
    const u64 b = (u64)F->b;
    if (gcd64(n, b) == 1) return 1;
    if (F->left) return 0;
    if (!F->rev) return 1;
    int d = 0, v = 0;
    u64 r = 0, m = n, g = 1;
    for (u64 t = n; t; t /= b) { r = r * b + t % b; d++; }
    for (u64 t = n; t % b == 0; t /= b) v++;
    for (int i = F->drop ? d - v : d; i > 0; i--) {   // g = gcd(n, b^len)
        u64 h = gcd64(m, b);
        if (h == 1) break;
        m /= h; g *= h;
    }
    return r % g == 0;
}

static fam_t g_F;
static pthread_mutex_t g_mx = PTHREAD_MUTEX_INITIALIZER;
static FILE *g_hits;
static u64 g_lo, g_hi, g_chunk, g_nchunks;
static atomic_ullong g_next, g_tested;
static atomic_uchar *g_done;

static void report(u64 n)
{
    pthread_mutex_lock(&g_mx);
    printf("HIT %s %llu\n", g_F.tag, (ull)n);
    fflush(stdout);
    if (g_hits) { fprintf(g_hits, "%llu\n", (ull)n); fflush(g_hits); }
    pthread_mutex_unlock(&g_mx);
}

typedef struct { u64 n, q; } cand_t;

static void run_group(const fam_t *F, const cand_t *c, int k)
{
    u64 nv[MAXL], out[MAXL];
    ring_t R[MAXL];
    for (int l = 0; l < k; l++) { nv[l] = c[l].n; ring_init(MONT, &R[l], c[l].q, F->b); }
    if (k == MAXL) core_mont4(F, nv, R, out);
    else for (int l = 0; l < k; l++) core_mont1(F, &nv[l], &R[l], &out[l]);
    for (int l = 0; l < k; l++)
        if (out[l] == 0) report(nv[l]);
}

static int g_fpv = VW > 1;                     // use the floating-point lanes when allowed
static int g_pow2 = 0;                         // 2-adic prefilter for even n before the lanes

static int g_noeval;                           // diagnostics: -z 1 skips the evaluation

// v_p(x) for x > 0 and a prime p dividing some base <= 36 (constant divisors, no div op)
static int vp(u64 x, int p)
{
    int c = 0;
#define VPC(P) case P: while (x % P == 0) { x /= P; c++; } return c;
    switch (p) {
    case 2: return __builtin_ctzll(x);
    VPC(3) VPC(5) VPC(7) VPC(11) VPC(13) VPC(17) VPC(19) VPC(23) VPC(29) VPC(31)
    }
#undef VPC
    while (x % p == 0) { x /= p; c++; }
    return c;
}

// prefilter() for n with gcd(n mod b, b) > 1, given the digits of n+1: same condition,
// gcd(n, b^len(n)) | rev(n), checked prime by prime from p-adic valuations.
static int prefilter_slow(const fam_t *F, u64 n, const unsigned char *dN, int Dn)
{
    if (F->left) return 0;
    if (!F->rev) return 1;
    const int b = F->b;
    int d[MAXD + 2], v = 0;
    for (int i = 0; i <= Dn; i++) d[i] = dN[i];
    int i = 0;                                                 // n = N - 1 (borrow)
    while (d[i] == 0) d[i++] = b - 1;
    d[i]--;
    while (v < Dn && d[v] == 0) v++;
    u64 r = 0;
    for (i = 0; i < Dn; i++) r = r * b + (u64)d[i];          // rev(n), d[0] most significant
    const int len = F->drop ? Dn - v : Dn;
    for (int p = 2, bb = b; bb > 1; p++) {
        if (bb % p) continue;
        int eb = 0;
        while (bb % p == 0) { bb /= p; eb++; }
        int en = vp(n, p), cap = len * eb;
        if ((en < cap ? en : cap) > vp(r, p)) return 0;
    }
    return 1;
}

typedef struct { u64 n; unsigned char dg[MAXD + 1]; } vcand_t;

static void run_fpv(const fam_t *F, const vcand_t *c, int k, int D)
{
    u64 nv[VW], res[VW];
    const unsigned char *dg[VW];
    if (k <= 0 || g_noeval) return;
    for (int l = 0; l < VW; l++) { int i = l < k ? l : k - 1; nv[l] = c[i].n; dg[l] = c[i].dg; }
    core_fpv(F, nv, D, dg, res);
    for (int l = 0; l < k; l++)
        if (res[l] == 0) report(nv[l]);
}

// Floating-point lanes: walk n with the digits of n+1 kept incrementally, so the hot path
// has no divisions (n mod b is a table lookup on the last digit).
static u64 scan_fpv(const fam_t *F, u64 a, u64 z)
{
    const int b = F->b;
    unsigned char dN[MAXD + 2] = {0}, cop[MAXB];
    int lenN = 0, Dn = ndig(a, b), k = 0, kd = 0;
    for (int r = 0; r < b; r++) cop[r] = gcd64((u64)r, (u64)b) == 1;
    for (u64 t = a + 1; t; t /= b) dN[lenN++] = (unsigned char)(t % b);
    vcand_t buf[VW];
    u64 cnt = 0;
    for (u64 n = a; n < z; n++) {
        const int r = dN[0] ? dN[0] - 1 : b - 1;              // n mod b
        if (cop[r] || prefilter_slow(F, n, dN, Dn)) {
            if (k && Dn != kd) { run_fpv(F, buf, k, kd); k = 0; }
            buf[k].n = n;
            memcpy(buf[k].dg, dN, (size_t)Dn + 1);
            kd = Dn; k++; cnt++;
            if (k == VW) { run_fpv(F, buf, k, kd); k = 0; }
        }
        Dn = lenN;                                            // digit count of n+1
        int i = 0;                                            // N -> N+1
        while (dN[i] == b - 1) dN[i++] = 0;
        dN[i]++;
        if (i >= lenN) lenN = i + 1;
    }
    if (k) run_fpv(F, buf, k, kd);
    return cnt;
}

// Test every n in [a, z); returns the number of full evaluations.
static u64 scan(const fam_t *F, u64 a, u64 z)
{
    if (g_fpv && F->rev && z <= FPV_MAXN && !g_pow2) return scan_fpv(F, a, z);
    cand_t buf[MAXL];
    int k = 0, kd = 0;
    u64 cnt = 0;
    ring_t R2;
    ring_init(POW2, &R2, 0, F->b);
    for (u64 n = a; n < z; n++) {
        if (!prefilter(F, n)) continue;
        u64 q = n;
        if (!(n & 1)) {                       // cheap 2-adic test first
            int s = __builtin_ctzll(n);
            u64 v;
            q = n >> s;
            core_pow2(F, &n, &R2, &v);
            if (v & ((1ull << s) - 1)) continue;
            if (q == 1) { report(n); continue; }
        }
        int d = ndig(n, F->b);
        cnt++;
        if (k && d != kd) { run_group(F, buf, k); k = 0; }
        buf[k].n = n; buf[k].q = q; k++; kd = d;
        if (k == MAXL) { run_group(F, buf, k); k = 0; }
    }
    if (k) run_group(F, buf, k);
    return cnt;
}

static void *worker(void *arg)
{
    (void)arg;
    for (;;) {
        u64 c = atomic_fetch_add(&g_next, 1);
        if (c >= g_nchunks) break;
        u64 a = g_lo + c * g_chunk, z = a + g_chunk;
        if (z > g_hi || z < a) z = g_hi;
        atomic_fetch_add(&g_tested, scan(&g_F, a, z));
        atomic_store(&g_done[c], 1);
    }
    return 0;
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static int search(u64 lo, u64 hi, int nth, u64 chunk, const char *dir)
{
    char path[1024], prog[1100];
    if (dir) mkdir(dir, 0777);
    snprintf(path, sizeof path, "%s/hits_%s.txt", dir ? dir : ".", g_F.tag[0] ? g_F.tag : "x");
    if (g_F.anum) snprintf(path, sizeof path, "%s/hits_A%06d.txt", dir ? dir : ".", g_F.anum);
    snprintf(prog, sizeof prog, "%s.progress", path);
    g_lo = lo; g_hi = hi; g_chunk = chunk;
    g_nchunks = (hi - lo + chunk - 1) / chunk;
    g_done = calloc(g_nchunks, 1);
    g_hits = fopen(path, "a");
    if (g_hits) { fprintf(g_hits, "# search [%llu, %llu)\n", (ull)lo, (ull)hi); fflush(g_hits); }
    fprintf(stderr, "search %s base %d [%llu, %llu) threads=%d chunk=%llu -> %s\n", g_F.tag, g_F.b,
            (ull)lo, (ull)hi, nth, (ull)chunk, path);
    pthread_t th[1024];
    for (int i = 0; i < nth; i++) pthread_create(&th[i], 0, worker, 0);
    double t0 = now(), tlast = t0;
    u64 contig = 0;
    for (;;) {
        usleep(100000);
        while (contig < g_nchunks && atomic_load(&g_done[contig])) contig++;
        double t = now();
        if (contig == g_nchunks || t - tlast >= 60) {
            u64 upto = contig == g_nchunks ? hi : lo + contig * chunk;
            fprintf(stderr, "[%8.0fs] complete below %llu  (%.3g n/s)\n", t - t0, (ull)upto,
                    (double)(upto - lo) / (t - t0));
            FILE *f = fopen(prog, "w");
            if (f) { fprintf(f, "%s base %d complete [%llu, %llu)\n", g_F.tag, g_F.b, (ull)lo, (ull)upto); fclose(f); }
            tlast = t;
        }
        if (contig == g_nchunks) break;
    }
    for (int i = 0; i < nth; i++) pthread_join(th[i], 0);
    if (g_hits) { fprintf(g_hits, "# done [%llu, %llu) in %.0fs\n", (ull)lo, (ull)hi, now() - t0); fclose(g_hits); }
    fprintf(stderr, "done [%llu, %llu) in %.1fs\n", (ull)lo, (ull)hi, now() - t0);
    return 0;
}

// ------------------------------------------------------------------ main

static int parse_seq(const char *s, fam_t *F)
{
    memset(F, 0, sizeof *F);
    if (s[0] == 'A' || s[0] == 'a') {
        int a = atoi(s + 1);
        for (int i = 0; i < 6; i++)
            if (a >= FAMS[i].first && a < FAMS[i].first + 24) {
                strcpy(F->tag, FAMS[i].tag);
                F->b = a - FAMS[i].first + 2;
                F->anum = a;
            }
        if (!F->b) return 0;
    } else {
        if (strlen(s) < 4 || s[2] != ':') return 0;
        memcpy(F->tag, s, 2);
        F->b = atoi(s + 3);
        for (int i = 0; i < 6; i++)
            if (!strcmp(F->tag, FAMS[i].tag) && F->b >= 2 && F->b <= 25) F->anum = FAMS[i].first + F->b - 2;
    }
    if (!strchr("LR", F->tag[0]) || !strchr("nkd", F->tag[1]) || F->b < 2 || F->b >= MAXB) return 0;
    F->left = F->tag[0] == 'L';
    F->rev = F->tag[1] != 'n';
    F->drop = F->tag[1] == 'd';
    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 5 || !parse_seq(argv[2], &g_F)) {
        fprintf(stderr,
                "usage: %s res|resl|resv|bench|search SEQ LO HI|COUNT [-t THREADS] [-c CHUNK] [-o DIR] [-x 0|1]\n"
                "  SEQ = A-number (A029447..A029542, A061931..A061978) or TAG:BASE, TAG in\n"
                "        Rn Ln (normal digits), Rk Lk (reversed, zeros kept), Rd Ld (reversed, zeros dropped)\n",
                argv[0]);
        return 1;
    }
    const fam_t *F = &g_F;
    u64 lo = strtoull(argv[3], 0, 0), hi = strtoull(argv[4], 0, 0);
    if (lo < 1) lo = 1;
    if (hi > (1ull << 58)) { fprintf(stderr, "n must stay below 2^58\n"); return 1; }

    if (!strcmp(argv[1], "res")) {
        for (u64 n = lo; n < hi; n++) printf("%llu %llu\n", (ull)n, (ull)full_residue(F, n));
        return 0;
    }
    if (!strcmp(argv[1], "resl")) {
        u64 n = lo;
        while (n < hi) {                      // 4 consecutive odd n with equal digit count
            u64 nv[MAXL], out[MAXL];
            ring_t R[MAXL];
            int ok = (n & 1) && n + 2 * (MAXL - 1) < hi && n > 1;
            for (int l = 0; ok && l < MAXL; l++) ok = ndig(n + 2 * l, F->b) == ndig(n, F->b);
            if (!ok) { printf("%llu %llu\n", (ull)n, (ull)full_residue(F, n)); n++; continue; }
            for (int l = 0; l < MAXL; l++) { nv[l] = n + 2 * l; ring_init(MONT, &R[l], nv[l], F->b); }
            core_mont4(F, nv, R, out);
            for (int l = 0; l < MAXL; l++) {
                printf("%llu %llu\n", (ull)nv[l], (ull)mmul(out[l], 1, nv[l], R[l].qi));
                if (l < MAXL - 1) printf("%llu %llu\n", (ull)nv[l] + 1, (ull)full_residue(F, nv[l] + 1));
            }
            n += 2 * MAXL - 1;
        }
        return 0;
    }
    for (int i = 5; i + 1 < argc; i += 2)
        if (!strcmp(argv[i], "-x")) g_fpv = atoi(argv[i + 1]) && VW > 0;
        else if (!strcmp(argv[i], "-p")) g_pow2 = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-z")) g_noeval = atoi(argv[i + 1]);
    if (!strcmp(argv[1], "pftest")) {         // prefilter_slow (digits) vs prefilter (gcd)
        u64 bad = 0;
        for (u64 n = lo; n < hi; n++) {
            unsigned char dN[MAXD + 2] = {0};
            int len = 0;
            for (u64 t = n + 1; t; t /= F->b) dN[len++] = (unsigned char)(t % F->b);
            int r = (int)(n % F->b), fast = gcd64((u64)r, (u64)F->b) == 1 || prefilter_slow(F, n, dN, ndig(n, F->b));
            if (fast != prefilter(F, n)) { if (bad++ < 5) printf("prefilter mismatch n=%llu\n", (ull)n); }
        }
        printf("pftest %s base %d [%llu, %llu): %llu mismatches\n", F->tag, F->b, (ull)lo, (ull)hi, (ull)bad);
        return bad != 0;
    }
    if (!strcmp(argv[1], "resv")) {           // residues through the floating-point lanes
        if (!F->rev || hi > FPV_MAXN) { fprintf(stderr, "resv: reversed families, n < 2^43\n"); return 1; }
        u64 n = lo;
        while (n < hi) {
            u64 nv[VW], res[VW];
            int k = 0, d = ndig(n, F->b);
            while (k < VW && n + k < hi && ndig(n + k, F->b) == d) { nv[k] = n + k; k++; }
            for (int l = k; l < VW; l++) nv[l] = nv[k - 1];
            core_fpv(F, nv, d, NULL, res);
            for (int l = 0; l < k; l++) printf("%llu %llu\n", (ull)nv[l], (ull)res[l]);
            n += k;
        }
        return 0;
    }
    if (!strcmp(argv[1], "bench")) {
        double t0 = now();
        u64 cnt = scan(F, lo, lo + hi);
        double t = now() - t0;
        printf("%s base %d near %.3g: %.3f us per n, %.3f us per evaluated n (%llu of %llu)\n", F->tag,
               F->b, (double)lo, 1e6 * t / hi, cnt ? 1e6 * t / cnt : 0.0, (ull)cnt, (ull)hi);
        return 0;
    }
    if (!strcmp(argv[1], "search")) {
        int nth = (int)sysconf(_SC_NPROCESSORS_ONLN);
        u64 chunk = 1ull << 22;
        const char *dir = 0;
        for (int i = 5; i + 1 < argc; i += 2) {
            if (!strcmp(argv[i], "-t")) nth = atoi(argv[i + 1]);
            else if (!strcmp(argv[i], "-c")) chunk = strtoull(argv[i + 1], 0, 0);
            else if (!strcmp(argv[i], "-o")) dir = argv[i + 1];
        }
        if (nth < 1) nth = 1;
        if (nth > 1024) nth = 1024;
        return search(lo, hi, nth, chunk, dir);
    }
    fprintf(stderr, "unknown mode %s\n", argv[1]);
    return 1;
}
