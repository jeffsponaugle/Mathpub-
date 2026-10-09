// core.h - shared (C++ host / Metal device) core of the double-base palindrome search.
//
// Numbers are held in "base-Q limbs": N = sum l[i] * rho^i with rho = Q^c <= 2^63,
// so the base-Q digits of N can be read straight out of the limbs (no big division).
//
// Search model (see README.md for the full derivation):
//   N is a base-P palindrome of length L.  Fix the k outer digit pairs (a "node"):
//       N = C + M * P^k,   C = sum_i d_i (P^(L-1-i) + P^i),   M = middle palindrome of length m = L-2k.
//   The node interval [C, C + (P^m-1) P^k] fixes the top s base-Q digits T of N; if N is also a
//   base-Q palindrome its low s digits are rev(T), i.e.  M * P^k == rev(T) - C  (mod Q^s).
//   All middles M are pre-bucketed by V = M*P^k mod Q^c, so each node costs one bucket probe.
#pragma once

#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef ulong u64;
typedef uint u32;
#define PCONST constant
#define PDEV device
#define TH thread
#define FN inline
#define CLZ64(x) ((u32)clz(x))
#define MULHI64(a, b) mulhi(a, b)
#define RBIT64(x) reverse_bits(x)
#define UMIN(a, b) min(a, b)
#elif defined(__CUDACC_RTC__) || defined(__CUDACC__)
typedef unsigned long long u64;
typedef unsigned int u32;
#define PCONST
#define PDEV
#define TH
#define FN __device__ __forceinline__
#define CLZ64(x) ((u32)__clzll((long long)(x)))
#define MULHI64(a, b) __umul64hi((a), (b))
#define RBIT64(x) __brevll(x)
#define UMIN(a, b) ((a) < (b) ? (a) : (b))
#else
#include <cstdint>
typedef uint64_t u64;
typedef uint32_t u32;
#define PCONST
#define PDEV
#define TH
#define FN static inline
#define CLZ64(x) ((u32)__builtin_clzll(x))
#define MULHI64(a, b) ((u64)(((unsigned __int128)(a) * (u64)(b)) >> 64))
#if defined(__clang__)
#define RBIT64(x) __builtin_bitreverse64(x)
#else  // gcc has no bit-reversal builtin
static inline u64 rbit64_portable(u64 x) {
  x = ((x >> 1) & 0x5555555555555555ull) | ((x & 0x5555555555555555ull) << 1);
  x = ((x >> 2) & 0x3333333333333333ull) | ((x & 0x3333333333333333ull) << 2);
  x = ((x >> 4) & 0x0f0f0f0f0f0f0f0full) | ((x & 0x0f0f0f0f0f0f0f0full) << 4);
  return __builtin_bswap64(x);
}
#define RBIT64(x) rbit64_portable(x)
#endif
#define UMIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef NLIMB
#define NLIMB 3  // limbs per number: N < rho^3 (>= ~10^54 for every supported base); build with -DNLIMB=4 for more
#endif
#define MAXQP 66   // size of the Q^j tables (j = 0..c, c <= 63)
#define MAXSEG 6   // max number of base-Q lengths spanned by one base-P length
#define MAXCEN 64  // max number of allowed centre digits

// Division by an invariant d of numerators n < 2^63 (Granlund-Montgomery, N = 63):
//   kind 0: d = 2^sh            -> n >> sh
//   kind 1: m = ceil(2^(63+l)/d) -> mulhi(n, m) >> (l-1),  l = ceil(log2 d)
struct FastDiv {
  u64 m;
  u32 sh;
  u32 kind;
};

struct NumQ {
  u64 l[NLIMB];
};

// One base-Q length LQ inside the current base-P length.
struct SegInfo {
  NumQ lo;      // Q^(LQ-1)
  NumQ hi;      // Q^LQ - 1
  u32 LQ;       // number of base-Q digits
  u32 allowed;  // 0 => skip (parity / prime rules)
  u32 S;        // number of top digits extracted = min(c, LQ)
  u32 ia;       // top-S extraction: a = LQ - S, ia = a / c
  u32 oa;       //                   oa = a % c
  u32 iz;       // zeroing the low (LQ - t) digits: iz = (LQ-t)/c
  u32 oz;       //                                  oz = (LQ-t)%c
  u32 pad;
};

struct Params {
  u64 rho;  // Q^c
  u64 Dt;   // number of table buckets
  u64 gt;   // gcd(P^k, Q^t)
  FastDiv gtdiv;
  u32 Q, c, P, k;
  u32 m, hm, hasCenter, NC;  // middle length, middle pairs, centre?, #centre digits
  u32 t, ueff, gmode, nseg;  // key digits, check digits, gcd(P,Q)>1?, #segments
  u32 qbits, revR, recLoN, pad1;  // Q = 2^qbits (0 if not), reversal chunk digits, #low reconstruction rows
  FastDiv recDiv;                  // division by recLoN (index -> (hi, lo) reconstruction split)
  u64 qpow[MAXQP];
  FastDiv qdiv[MAXQP];
  u64 gs[MAXQP];  // gcd(P^k, Q^s)
  FastDiv gsdiv[MAXQP];
  u32 ndlo[68];  // ndlo[b] = #base-Q digits of 2^(b-1)
  u32 centerDig[MAXCEN];
  NumQ Wmax;  // (P^m - 1) * P^k
  SegInfo seg[MAXSEG];
  // ---- v2 (residue-class partitioned) search; t = r + tp there --------------------------------
  u64 Qr;       // Q^r: number of classes (N mod Q^r)
  u32 r, tp;    // class digits, key digits inside a class
  u32 khi, klo; // outer digits split: A = (A_hi, A_lo)
  u32 h1, nB1;  // middle split: B = (B1: first h1 middle pairs, B2: the rest + centre)
  u32 nB2, w;   // |B2|, classes per batch
  u32 v2u, pad3;  // check digits stored in an entry (low 20 bits)
  FastDiv nB2div;  // idxB -> (b1, b2)
  NumQ Zb2;        // bound on the B2 part of M * P^k: N in [C + B1part, C + B1part + Zb2)
  u32 r0, pad4;    // low base-Q digits skipped by the class (1 when Q = 2 and P is odd: N, V odd)
};

// Parameter access.  The GPU build defines SPEC_* (values of the current search) so that the
// Metal compiler sees constants; the CPU build reads the Params struct.
#ifdef SPEC_Q
#define PRM_R0 ((u32)SPEC_R0)
#define PRM_R ((u32)SPEC_R)
#define PRM_TP ((u32)SPEC_TP)
#define PRM_V2U ((u32)SPEC_V2U)
#define PRM_Q ((u32)SPEC_Q)
#define PRM_P ((u32)SPEC_P)
#define PRM_C ((u32)SPEC_C)
#define PRM_RHO ((u64)SPEC_RHO)
#define PRM_QBITS ((u32)SPEC_QBITS)
#define PRM_GMODE ((u32)SPEC_GMODE)
#define PRM_T ((u32)SPEC_T)
#define PRM_UEFF ((u32)SPEC_UEFF)
#define PRM_NSEG ((u32)SPEC_NSEG)
#define PRM_REVR ((u32)SPEC_REVR)
#else
#define PRM_R0 (p.r0)
#define PRM_R (p.r)
#define PRM_TP (p.tp)
#define PRM_V2U (p.v2u)
#define PRM_Q (p.Q)
#define PRM_P (p.P)
#define PRM_C (p.c)
#define PRM_RHO (p.rho)
#define PRM_QBITS (p.qbits)
#define PRM_GMODE (p.gmode)
#define PRM_T (p.t)
#define PRM_UEFF (p.ueff)
#define PRM_NSEG (p.nseg)
#define PRM_REVR (p.revR)
#endif

FN u64 fdiv(u64 n, FastDiv d) { return d.kind ? (MULHI64(n, d.m) >> d.sh) : (n >> d.sh); }

// x / Q^j, x mod Q^j, x * Q^j  (shifts and masks when Q is a power of two)
#define QDIV(x, j) (PRM_QBITS ? ((x) >> ((j) * PRM_QBITS)) : fdiv((x), p.qdiv[(j)]))
#define QMUL(x, j) (PRM_QBITS ? ((x) << ((j) * PRM_QBITS)) : (x) * p.qpow[(j)])
#define QMOD(x, j) (PRM_QBITS ? ((x) & ((((u64)1) << ((j) * PRM_QBITS)) - 1)) : (x) - fdiv((x), p.qdiv[(j)]) * p.qpow[(j)])

#if NLIMB == 3
FN u64 limb(NumQ x, u32 i) { return i == 0 ? x.l[0] : (i == 1 ? x.l[1] : x.l[2]); }
#else
FN u64 limb(NumQ x, u32 i) { return i == 0 ? x.l[0] : (i == 1 ? x.l[1] : (i == 2 ? x.l[2] : x.l[3])); }
#endif

FN NumQ nq_add(NumQ a, NumQ b, u64 rho) {
  NumQ r;
  u64 cy = 0;
  for (int i = 0; i < NLIMB; i++) {
    u64 s = a.l[i] + b.l[i] + cy;
    cy = (s >= rho) ? 1 : 0;
    r.l[i] = cy ? s - rho : s;
  }
  return r;
}

FN NumQ nq_sub1(NumQ a, u64 rho) {  // a - 1, a > 0
  bool borrow = true;
  for (int i = 0; i < NLIMB; i++) {
    if (borrow) {
      if (a.l[i] > 0) {
        a.l[i] -= 1;
        borrow = false;
      } else {
        a.l[i] = rho - 1;
      }
    }
  }
  return a;
}

FN int nq_cmp(NumQ a, NumQ b) {
  for (int i = NLIMB - 1; i >= 0; i--)
    if (a.l[i] != b.l[i]) return a.l[i] < b.l[i] ? -1 : 1;
  return 0;
}

// x with its low z digits cleared (z = iz*c + oz)
FN NumQ nq_zero_low(NumQ x, u32 iz, u32 oz, PCONST const Params& p) {
  NumQ r;
  for (int i = 0; i < NLIMB; i++) {
    u64 v = x.l[i];
    if ((u32)i < iz)
      v = 0;
    else if ((u32)i == iz)
      v = v - QMOD(v, oz);
    r.l[i] = v;
  }
  return r;
}

// Value of the top S digits of x (positions a..a+S-1, a = ia*c + oa), S <= c, x < Q^(a+S).
FN u64 nq_top(NumQ x, u32 ia, u32 oa, PCONST const Params& p) {
  u64 v = QDIV(limb(x, ia), oa);
  if (ia + 1 < NLIMB) v += QMUL(limb(x, ia + 1), PRM_C - oa);
  return v;
}

// Digits [pos, pos+cnt) of x as a number (cnt <= c).
FN u64 nq_get(NumQ x, u32 pos, u32 cnt, PCONST const Params& p) {
  u32 ia = pos / PRM_C, oa = pos - ia * PRM_C;
  u64 lo = QDIV(limb(x, ia), oa);
  u32 avail = PRM_C - oa;
  if (cnt <= avail) return QMOD(lo, cnt);
  u32 need = cnt - avail;
  u64 hi = (ia + 1 < NLIMB) ? limb(x, ia + 1) : 0;
  hi = QMOD(hi, need);
  return lo + QMUL(hi, avail);
}

// Number of base-Q digits of d >= 1 (d < rho).
FN u32 ndigits(u64 d, PCONST const Params& p) {
  u32 b = 64 - CLZ64(d);
  u32 j = p.ndlo[b];
  while (j < 64 && p.qpow[j] <= d) j++;
  return j;
}

// Shared leading digits of the S-digit numbers x <= y; T receives their common prefix.
FN u32 shared_prefix(u64 x, u64 y, u32 S, PCONST const Params& p, TH u64& T) {
  if (x == y) {
    T = x;
    return S;
  }
  if (PRM_QBITS) {
    u32 lb = 64 - CLZ64(x ^ y);
    u32 j = (lb + PRM_QBITS - 1) / PRM_QBITS;
    T = x >> (j * PRM_QBITS);
    return S - j;
  }
  u32 j = ndigits(y - x, p);
  for (;;) {
    u64 qx = fdiv(x, p.qdiv[j]), qy = fdiv(y, p.qdiv[j]);
    if (qx == qy) {
      T = qx;
      return S - j;
    }
    j++;
  }
}

// rev_s(T): the s-digit reversal of T (T < Q^s).
FN u64 rev_digits(u64 T, u32 s, PCONST const Params& p, PDEV const u32* revtab) {
  if (s == 0) return 0;
  if (PRM_Q == 2) return RBIT64(T) >> (64 - s);
  u64 R = 0;
  u32 r = PRM_REVR;
  while (s >= r) {
    u64 q = QDIV(T, r);
    u64 ch = T - QMUL(q, r);
    R = QMUL(R, r) + revtab[ch];
    T = q;
    s -= r;
  }
  if (s) R = QMUL(R, s) + QDIV((u64)revtab[T], r - s);
  return R;
}

FN bool nq_is_pal(NumQ N, PCONST const Params& p, PDEV const u32* revtab) {
  int top = NLIMB - 1;
  while (top > 0 && N.l[top] == 0) top--;
  if (N.l[top] == 0) return true;
  u32 L = (u32)top * PRM_C + ndigits(N.l[top], p);
  u32 hlen = L / 2, pos = 0;
  while (pos < hlen) {
    u32 cnt = UMIN(PRM_C, hlen - pos);
    u64 lowv = nq_get(N, pos, cnt, p);
    u64 highv = nq_get(N, L - pos - cnt, cnt, p);
    if (rev_digits(highv, cnt, p, revtab) != lowv) return false;
    pos += cnt;
  }
  return true;
}

// N = C + M * P^k for the middle with table index idx = hi * recLoN + lo
// (recHi[hi] / recLo[lo] hold the precomputed contributions of the two halves of the middle).
FN NumQ reconstruct(NumQ C, u32 idx, PCONST const Params& p, PDEV const NumQ* recHi, PDEV const NumQ* recLo) {
  u32 hi = (u32)fdiv(idx, p.recDiv);
  u32 lo = idx - hi * p.recLoN;
  return nq_add(nq_add(C, recHi[hi], PRM_RHO), recLo[lo], PRM_RHO);
}

// A bucket probe prepared for one piece of a node interval.
struct ProbeReq {
  u64 key;  // bucket
  u64 dh;   // required check digits (mod Q^uu)
  u32 uu;   // number of usable check digits
  u32 flag; // 1 = probe, 2 = overflow (s < t; never produced by the planner's t)
};

// Piece [a, b] of a node interval, entirely inside one level-t block of segment sg.
FN ProbeReq prep_piece(NumQ C, NumQ a, NumQ b, PCONST const SegInfo& sg, PCONST const Params& p,
                       PDEV const u32* revtab) {
  ProbeReq rq;
  rq.key = 0;
  rq.dh = 0;
  rq.uu = 0;
  rq.flag = 0;
  u64 x = nq_top(a, sg.ia, sg.oa, p), y = nq_top(b, sg.ia, sg.oa, p);
  u64 T;
  u32 s = shared_prefix(x, y, sg.S, p, T);
  if (s < PRM_T) {
    rq.flag = 2;
    return rq;
  }
  u64 R = rev_digits(T, s, p, revtab);
  u64 c0 = C.l[0];
  u64 diff = (R >= c0) ? (R - c0) : (R + (PRM_RHO - c0));  // (rev(T) - C) mod Q^c
  if (PRM_GMODE) {
    u64 q = fdiv(diff, p.gsdiv[s]);
    if (diff != q * p.gs[s]) return rq;  // N == C (mod P^k) contradicts N == rev(T) (mod Q^s)
  }
  u64 qt = QDIV(diff, PRM_T);
  u64 lowt = diff - QMUL(qt, PRM_T);
  rq.key = PRM_GMODE ? fdiv(lowt, p.gtdiv) : lowt;
  rq.uu = UMIN(s - PRM_T, PRM_UEFF);
  rq.dh = QMOD(qt, rq.uu);
  rq.flag = 1;
  return rq;
}

// Full test of one candidate (entry index idx of node C).
template <typename Sink>
FN void check_candidate(NumQ C, u32 idx, PCONST const Params& p, PDEV const u32* revtab, PDEV const NumQ* recHi,
                        PDEV const NumQ* recLo, TH Sink& sink) {
  NumQ N = reconstruct(C, idx, p, recHi, recLo);
  if (nq_is_pal(N, p, revtab)) sink.emit(N);
}

// Scan bucket entries [beg, end) of a prepared probe; entries whose check digits agree are
// handed to sink.candidate(C, idx) (immediate full test on the CPU, deferred on the GPU).
template <typename Sink>
FN void scan_bucket(ProbeReq rq, u32 beg, u32 end, NumQ C, PCONST const Params& p, PDEV const u64* ents,
                    TH Sink& sink) {
  sink.probes += 1;
  for (u32 i = beg; i < end; i++) {
    u64 e = ents[i];
    u64 chk = e & 0xffffffffull;
    if (QMOD(chk, rq.uu) != rq.dh) continue;
    sink.cands += 1;
    sink.candidate(C, (u32)(e >> 32));
  }
}

template <typename Sink>
FN void run_probe(ProbeReq rq, NumQ C, PCONST const Params& p, PDEV const u32* offs, PDEV const u64* ents,
                  TH Sink& sink) {
  scan_bucket(rq, offs[rq.key], offs[rq.key + 1], C, p, ents, sink);
}

// Enumerate the pieces of node C (split at base-Q length and level-t block boundaries)
// and hand each prepared probe to out.add(rq).
template <typename Out>
FN void node_requests(NumQ C, PCONST const Params& p, PDEV const u32* revtab, TH Out& out) {
  NumQ lo = C, hi = nq_add(C, p.Wmax, PRM_RHO);
  for (u32 si = 0; si < PRM_NSEG; si++) {
    PCONST const SegInfo& sg = p.seg[si];
    if (nq_cmp(hi, sg.lo) < 0) break;
    if (!sg.allowed || nq_cmp(lo, sg.hi) > 0) continue;
    NumQ a = nq_cmp(lo, sg.lo) < 0 ? sg.lo : lo;
    NumQ b = nq_cmp(hi, sg.hi) > 0 ? sg.hi : hi;
    for (;;) {  // top piece first
      NumQ beta = nq_zero_low(b, sg.iz, sg.oz, p);
      bool last = nq_cmp(beta, a) <= 0;
      ProbeReq rq = prep_piece(C, last ? a : beta, b, sg, p, revtab);
      if (rq.flag) out.add(rq);
      if (last) break;
      b = nq_sub1(beta, PRM_RHO);
    }
  }
}

// =============================================================================================
// v2: residue-class partitioned meet in the middle (coprime P, Q).
//
//   N = C_hi + C_lo + M * P^k.  Every N lies in one "A_hi piece" (A_hi interval cut at base-Q
//   length and level-r boundaries), whose top r digits fix N mod Q^r = Rr.  The class of N is
//   cls = V(B) mod Q^r = (Rr - C_hi - C_lo) mod Q^r.  Classes are processed in batches
//   [c0, c0 + w): only middles B and outer parts A_lo of those classes are generated, so the
//   bucket table of a batch is small (cache resident) while the total work stays ~ |A| + |B|.
// =============================================================================================
struct AhiPiece {
  NumQ C;   // C_hi
  NumQ a;   // piece range (one base-Q length, one level-(r0+r) block)
  NumQ b;
  u64 Rr;   // class base: ((N - C_hi) mod Q^(r0+r)) div Q^r0 for every N in the piece
  u32 seg;  // segment index
  u32 pad;
};

// v2 entry: (index of B << 20) | check digits (index < 2^44)
#define V2_IDX(e) ((e) >> 20)
#define V2_CHK(e) ((e)&0xfffffull)

// Stage-1 test of a candidate: N lies in [X, X + Zb2) with X = C + B1-part, so the shared top
// digits of X and X + Zb2 - 1 are N's top digits; their reversal must equal N's low digits,
// which are known from C and the entry ((C + V) mod Q^nlow).  Returns false if they differ.
FN bool v2_prefilter(NumQ C, u64 idxB, u64 Vlow, u32 nlow, PCONST const SegInfo& sg, PCONST const Params& p,
                     PDEV const u32* revtab, PDEV const NumQ* recB1) {
  u64 b1 = fdiv(idxB, p.nB2div);
  NumQ X = nq_add(C, recB1[b1], PRM_RHO);
  NumQ Y = nq_sub1(nq_add(X, p.Zb2, PRM_RHO), PRM_RHO);
  if (nq_cmp(X, sg.lo) < 0 || nq_cmp(Y, sg.hi) > 0) return true;  // not inside one base-Q length
  u64 x = nq_top(X, sg.ia, sg.oa, p), y = nq_top(Y, sg.ia, sg.oa, p);
  u64 T;
  u32 s = shared_prefix(x, y, sg.S, p, T);
  u64 R = rev_digits(T, s, p, revtab);
  u32 n = UMIN(s, nlow);
  u64 c0 = QMOD(C.l[0], nlow);
  u64 Nlow = c0 + QMOD(Vlow, nlow);
  return QMOD(Nlow, n) == QMOD(R, n);
}

// N = C + M * P^k for middle index idxB = b1 * nB2 + b2 (B2 part rebuilt from its digits).
FN NumQ v2_rebuild(NumQ C, u64 idxB, PCONST const Params& p, PDEV const NumQ* recB1, PDEV const NumQ* midW,
                   PDEV const NumQ* cenW) {
  u64 b1 = fdiv(idxB, p.nB2div);
  u64 b2 = idxB - b1 * p.nB2;
  NumQ N = nq_add(C, recB1[b1], PRM_RHO);
  if (p.hasCenter) {
    u64 q = b2 / p.NC;
    u32 cd = p.centerDig[(u32)(b2 - q * p.NC)];
    b2 = q;
    if (cd) N = nq_add(N, cenW[cd], PRM_RHO);
  }
  for (int j = (int)p.hm - 1; j >= (int)p.h1; j--) {
    u64 q = b2 / PRM_P;
    u32 e = (u32)(b2 - q * PRM_P);
    b2 = q;
    if (e) N = nq_add(N, midW[(u32)j * PRM_P + e], PRM_RHO);
  }
  return N;
}

// One node piece [a, b] (inside one level-(r+tp) block): probe the batch table.
template <typename Sink>
FN void v2_piece(NumQ C, NumQ a, NumQ b, PCONST const SegInfo& sg, u64 c0, PCONST const Params& p,
                 PDEV const u32* revtab, PDEV const u32* boffs, PDEV const u64* bents, TH Sink& sink) {
  u64 x = nq_top(a, sg.ia, sg.oa, p), y = nq_top(b, sg.ia, sg.oa, p);
  u64 T;
  u32 s = shared_prefix(x, y, sg.S, p, T);
  u32 rt = PRM_R0 + PRM_R + PRM_TP;
  if (s < rt) {
    sink.overflow();
    return;
  }
  u64 R = rev_digits(T, s, p, revtab);
  u64 c0l = C.l[0];
  u64 diff = (R >= c0l) ? (R - c0l) : (R + (PRM_RHO - c0l));
  u64 cls = QMOD(QDIV(diff, PRM_R0), PRM_R) - c0;  // in [0, w) by construction of the batch
  if (cls >= (u64)p.w) {               // cannot happen; never index outside the batch table
    sink.overflow();
    return;
  }
  u64 key = QMUL(cls, PRM_TP) + QMOD(QDIV(diff, PRM_R0 + PRM_R), PRM_TP);
  u32 uu = UMIN(s - rt, PRM_V2U);
  u64 dh = QMOD(QDIV(diff, rt), uu);
  u64 dlow = QMOD(diff, rt);  // V mod Q^(r+tp) for every entry of this bucket
  u32 beg = boffs[key], end = boffs[key + 1];
  sink.probes += 1;
  for (u32 i = beg; i < end; i++) {
    u64 e = bents[i];
    u64 chk = V2_CHK(e);
    if (QMOD(chk, uu) != dh) continue;
    sink.cands += 1;
    // stage 1: the leading middle digits (B1) fix more top digits of N than the node did;
    // compare them with the low digits of N = C + V known from the entry
    if (!v2_prefilter(C, V2_IDX(e), dlow + QMUL(chk, rt), rt + PRM_V2U, sg, p, revtab, sink.recB1)) continue;
    sink.pre += 1;
    sink.candidate(C, V2_IDX(e));
  }
}

// Node C = C_hi + C_lo, restricted to A_hi piece pc.
template <typename Sink>
FN void v2_node(NumQ C, NumQ pa, NumQ pb, u32 segi, u64 c0, PCONST const Params& p, PDEV const u32* revtab,
                PDEV const u32* boffs, PDEV const u64* bents, TH Sink& sink) {
  NumQ lo = C, hi = nq_add(C, p.Wmax, PRM_RHO);
  NumQ a = nq_cmp(lo, pa) < 0 ? pa : lo;
  NumQ b = nq_cmp(hi, pb) > 0 ? pb : hi;
  if (nq_cmp(a, b) > 0) return;
  PCONST const SegInfo& sg = p.seg[segi];
  for (;;) {  // split at level-(r+tp) blocks (iz/oz are set for t = r + tp)
    NumQ beta = nq_zero_low(b, sg.iz, sg.oz, p);
    bool last = nq_cmp(beta, a) <= 0;
    v2_piece(C, last ? a : beta, b, sg, c0, p, revtab, boffs, bents, sink);
    if (last) break;
    b = nq_sub1(beta, PRM_RHO);
  }
}

// key / check of a middle with V = M * P^k mod Q^c inside batch [c0, c0 + w)
FN u64 v2_key(u64 V, u64 c0, PCONST const Params& p) {
  return QMUL(QMOD(QDIV(V, PRM_R0), PRM_R) - c0, PRM_TP) + QMOD(QDIV(V, PRM_R0 + PRM_R), PRM_TP);
}
FN u64 v2_entry(u64 V, u64 idxB, PCONST const Params& p) {
  return (idxB << 20) | QMOD(QDIV(V, PRM_R0 + PRM_R + PRM_TP), PRM_V2U);
}
