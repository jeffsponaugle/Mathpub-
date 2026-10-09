// plan2.cpp - parameter choice and precomputation for the v2 (class-partitioned) search.
#include "plan2.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <thread>

static double dpow(double b, double e) { return std::pow(b, e); }

// largest u with Q^u < 2^20 (check digits stored in a v2 entry)
static u32 check_digits20(u32 Q) {
  u64 v = 1;
  u32 u = 0;
  while (v * Q < (1ull << 20)) {
    v *= Q;
    u++;
  }
  return u;
}

Plan2 plan2_choose(const SearchOpts& o, const V2Opts& vo, u32 L) {
  Plan2 p2;
  const u32 P = o.P, Q = o.Q;
  if (std::gcd(P, Q) != 1) {
    p2.why = "v2 needs coprime bases";
    return p2;
  }
  // basics from a v1 plan (parity rules, centres, leading digits, segments)
  SearchOpts ob = o;
  ob.noTableLimit = true;
  if (vo.k >= 0) ob.forceK = vo.k;  // (forcing k also disables the brute-force shortcut)
  Plan probe = make_plan(ob, L);
  if (probe.skip || probe.brute) {
    p2.why = probe.skip ? probe.skipReason : "small length (brute force)";
    p2.base = probe;
    return p2;
  }
  const Params& pp = probe.prm;
  u32 LQlo = pp.seg[0].LQ, LQhi = pp.seg[pp.nseg - 1].LQ, LQmin = 0;
  for (u32 i = 0; i < pp.nseg; i++)
    if (pp.seg[i].allowed) {
      LQmin = pp.seg[i].LQ;
      break;
    }
  (void)LQlo;
  (void)LQhi;
  u32 NCfull = 0;
  {
    // centre digits allowed (same rule as v1)
    bool oddOnly = (Q == 2) && (P % 2 == 1) && (L % 2 == 1);
    for (u32 d = 0; d < P; d++)
      if (!oddOnly || d % 2 == 1) NCfull++;
  }
  const double nTop = (double)probe.topDigits.size();
  const u32 c = pp.c;
  // Q = 2 and P odd: N, C and the B1 part are fixed mod 2 (N odd, C even), so classes skip bit 0
  const u32 r0 = (Q == 2 && P % 2 == 1 && L % 2 == 1) ? 1 : 0;
  double mem = (double)physical_memory_bytes();
  double memBudget = o.memBudget > 0 ? o.memBudget : 0.35 * mem;

  struct Best {
    double cost = 1e300;
    u32 k, khi, h1, r, tp;
    u64 w;
    double nA, nB, nb, batchE;
  } best;
  for (u32 k = 1; 2 * k < L; k++) {
    if (vo.k >= 0 && (int)k != vo.k) continue;
    u32 m = L - 2 * k, hm = m / 2;
    double NC = (m % 2) ? (double)NCfull : 1.0;
    double nB = dpow(P, hm) * NC;
    double nA = nTop * dpow(P, k - 1);
    if (nB > 1.7e13) continue;  // 44-bit entry index
    mpz_class Wk = (mpow(P, m) - 1) * mpow(P, k);
    int rtMax = std::min<int>((int)LQmin - (int)ndigits_base(Wk, Q), (int)c - 1) - (int)r0;
    if (rtMax < 2) continue;
    for (u32 khi = 1; khi < k; khi++) {
      if (vo.khi >= 0 && (int)khi != vo.khi) continue;
      double nAhi = nTop * dpow(P, khi - 1), nAlo = dpow(P, k - khi);
      mpz_class Whi = (mpow(P, L - 2 * khi) - 1) * mpow(P, khi);
      int rMax = std::min<int>((int)LQmin - (int)ndigits_base(Whi, Q) - (int)r0, rtMax - 1);
      if (rMax < 1) continue;
      for (u32 h1 = 0; h1 <= hm; h1++) {
        if (vo.h1 >= 0 && (int)h1 != vo.h1) continue;
        double nB1 = dpow(P, h1), nB2 = dpow(P, hm - h1) * NC;
        if (nB2 > 4.0e9 || nB1 > 4.0e9 || nAlo > 4.0e9 || nAhi > 2.0e9) continue;
        for (int r = 1; r <= rMax; r++) {
          if (vo.r >= 0 && r != vo.r) continue;
          double Qr = dpow(Q, r);
          if (Qr > (double)(1u << 30)) break;
          // classes per batch so that a batch table holds ~batchEntries middles
          double perClass = nB / Qr;
          u64 w = (u64)std::max(1.0, std::floor(vo.batchEntries / std::max(perClass, 1e-9)));
          w = std::min<u64>(w, (u64)Qr);
          double nb = std::ceil(Qr / (double)w);
          double batchE = perClass * (double)w;
          // key digits inside a class: ~2 entries per bucket
          int tp = vo.tp;
          if (tp < 0) {
            double want = std::max(1.0, batchE / 2.0 / (double)w);
            tp = (int)std::floor(std::log(want) / std::log((double)Q));
            tp = std::max(0, std::min(tp, rtMax - r));
          }
          if (r + tp > rtMax) continue;
          double buckets = (double)w * dpow(Q, tp);
          if (buckets > 3.9e9) continue;
          double nPieces = nAhi * 1.5;
          // batch tables: 3 in flight on the GPU, one per worker thread on the CPU
          double tables = o.gpu ? 3.0 : (double)std::max(1, o.threads);
          double memUse = nAlo * 32 + Qr * 8 + nB2 * 44 + nPieces * 112 + nB1 * 40 + tables * (batchE * 8 * 1.3 + buckets * 4);
          if (memUse > memBudget) continue;
          // candidates (full tests): entries agreeing with a node on all s known digits
          double sAvg = 0.5 * (double)(LQmin + pp.seg[pp.nseg - 1].LQ) - std::log(Wk.get_d()) / std::log((double)Q) - 1.0;
          double bucketSize = batchE / std::max(1.0, buckets);
          double cands = nA * std::min(bucketSize, nB * std::pow((double)Q, -sAvg));
          // cost model (seconds); constants calibrated on an M4 Max (GPU) / 14 threads (CPU)
          double batchBytes = batchE * 8 + buckets * 4;
          double cost;
          if (o.gpu) {
            // measured (M4 Max, L=41..43): ~1 ns per node, ~0.45 ns per candidate (stage-1 prefilter
            // inline + ~5% deferred full tests), ~0.6 ns per generated middle
            double tNode = batchBytes <= 48e6 ? 0.8e-9 : batchBytes <= 512e6 ? 1.0e-9 : 1.4e-9;
            tNode += 0.02e-9 * bucketSize;
            // the stage-1 prefilter passes ~ a few * P^-h1 of the candidates (h1 leading middle digits)
            double pass = std::min(1.0, 4.0 * std::pow((double)P, -(double)h1));
            double tCand = 0.35e-9 + pass * 1.5e-9;
            cost = nA * tNode + cands * tCand + nB * 0.6e-9 + nb * ((nPieces + nB1) * 0.03e-9 + 60e-6 + buckets * 0.1e-9);
          } else {
            double thr = std::max(1, o.threads);
            double eff = std::min(1.0, nb / thr);  // batches are the unit of parallelism
            double tNode = (batchBytes <= 16e6 ? 1.5e-9 : 2.8e-9) * 14.0 / thr;
            double pass = std::min(1.0, 4.0 * std::pow((double)P, -(double)h1));
            cost = (nA * tNode + cands * (0.6e-9 + pass * 1.8e-9) * 14.0 / thr + nB * 3e-9 * 14.0 / thr +
                    nb * ((nPieces + nB1) * 0.4e-9 * 14.0 / thr + buckets * 0.5e-9 * 14.0 / thr)) / eff;
          }
          if (cost < best.cost) best = {cost, k, khi, h1, (u32)r, (u32)tp, w, nA, nB, nb, batchE};
        }
      }
    }
  }
  if (best.cost >= 1e300) {
    p2.why = "no feasible v2 parameters";
    p2.base = probe;
    return p2;
  }
  // full plan for the chosen k with t = r + tp
  ob.forceK = (int)best.k;
  ob.forceT = (int)(r0 + best.r + best.tp);
  p2.base = make_plan(ob, L);
  Params& p = p2.base.prm;
  if (p.t != r0 + best.r + best.tp) {
    p2.why = "key digits clamped by the planner";
    return p2;
  }
  p.r0 = r0;
  p.r = best.r;
  p.tp = best.tp;
  p.Qr = (u64)dpow(Q, best.r);
  p.khi = best.khi;
  p.klo = best.k - best.khi;
  p.h1 = best.h1;
  p.nB1 = (u32)dpow(P, best.h1);
  p.nB2 = (u32)(dpow(P, p.hm - best.h1) * p.NC);
  p.w = (u32)best.w;
  p.v2u = std::min(check_digits20(Q), c - r0 - best.r - best.tp);
  p.nB2div = make_fastdiv(p.nB2);
  // B2 part of M * P^k is below P^(k + m - h1)
  p.Zb2 = to_numq(mpow(P, best.k + p.m - best.h1), p);
  p2.nA = (u64)best.nA;
  p2.nB = (u64)best.nB;
  p2.nAlo = (u64)dpow(P, p.klo);
  p2.nBatches = (u64)best.nb;
  p2.estSeconds = best.cost;
  p2.batchEntriesEst = best.batchE;
  p2.ok = true;
  return p2;
}

// ---- builders ------------------------------------------------------------------------------
template <class F>
static void enum_digits(u32 n, u32 P, F&& f) {  // all n-digit tuples, first digit most significant
  std::vector<u32> d(n, 0);
  u64 idx = 0;
  for (;;) {
    f(idx, d);
    idx++;
    int i = (int)n - 1;
    while (i >= 0 && ++d[i] == P) d[i--] = 0;
    if (i < 0) break;
  }
}

void plan2_build(Plan2& p2, int threads) {
  (void)threads;
  Plan& pl = p2.base;
  Params& p = pl.prm;
  const u32 P = p.P;
  const u64 Qr = p.Qr;
  // ---- A_hi pieces --------------------------------------------------------------------------
  p2.ahi.clear();
  const NumQ Whi = pl.WmaxAt[p.khi];
  std::vector<u32> segZi(p.nseg), segZo(p.nseg);
  const u32 rr = p.r0 + p.r;  // digits fixed per A_hi piece
  const u64 Qrr = p.qpow[rr], Q0 = p.qpow[p.r0];
  for (u32 si = 0; si < p.nseg; si++) {
    u32 z = p.seg[si].LQ - rr;
    segZi[si] = z / p.c;
    segZo[si] = z % p.c;
  }
  auto add_hi = [&](const NumQ& Chi) {
    NumQ lo = Chi, hi = nq_add(Chi, Whi, p.rho);
    for (u32 si = 0; si < p.nseg; si++) {
      const SegInfo& sg = p.seg[si];
      if (nq_cmp(hi, sg.lo) < 0) break;
      if (!sg.allowed || nq_cmp(lo, sg.hi) > 0) continue;
      NumQ a = nq_cmp(lo, sg.lo) < 0 ? sg.lo : lo;
      NumQ b = nq_cmp(hi, sg.hi) > 0 ? sg.hi : hi;
      for (;;) {
        NumQ beta = nq_zero_low(b, segZi[si], segZo[si], p);
        bool last = nq_cmp(beta, a) <= 0;
        NumQ pa = last ? a : beta;
        AhiPiece pc{};
        pc.C = Chi;
        pc.a = pa;
        pc.b = b;
        pc.seg = si;
        u64 top = nq_top(pa, sg.ia, sg.oa, p);   // top S digits
        u64 Tr = top / p.qpow[sg.S - rr];         // top r0 + r digits
        u64 Nlow = rev_digits(Tr, rr, p, pl.revtab.data());  // N mod Q^(r0+r)
        u64 x = (Nlow + Qrr - Chi.l[0] % Qrr) % Qrr;
        if (x % Q0 != (p.r0 ? 1 : 0)) throw std::runtime_error("v2: unexpected low digits");
        pc.Rr = (x / Q0) % Qr;
        p2.ahi.push_back(pc);
        if (last) break;
        b = nq_sub1(beta, p.rho);
      }
    }
  };
  for (u32 d0 : pl.topDigits) {
    NumQ C0 = nq_add(NumQ{}, pl.outerW[d0], p.rho);
    if (p.khi == 1) {
      add_hi(C0);
      continue;
    }
    enum_digits(p.khi - 1, P, [&](u64, const std::vector<u32>& d) {
      NumQ C = C0;
      for (u32 i = 0; i + 1 < p.khi; i++)
        if (d[i]) C = nq_add(C, pl.outerW[(size_t)(i + 1) * P + d[i]], p.rho);
      add_hi(C);
    });
  }
  // ---- A_lo CSR by C_lo mod Q^r -------------------------------------------------------------
  {
    std::vector<NumQ> all;
    all.reserve(p2.nAlo);
    enum_digits(p.klo, P, [&](u64, const std::vector<u32>& d) {
      NumQ C{};
      for (u32 i = 0; i < p.klo; i++)
        if (d[i]) C = nq_add(C, pl.outerW[(size_t)(p.khi + i) * P + d[i]], p.rho);
      all.push_back(C);
    });
    p2.aloOffs.assign(Qr + 1, 0);
    auto key = [&](const NumQ& C) {
      if (C.l[0] % p.qpow[p.r0] != 0) throw std::runtime_error("v2: C_lo low digits not zero");
      return (C.l[0] / p.qpow[p.r0]) % Qr;
    };
    for (auto& C : all) p2.aloOffs[key(C) + 1]++;
    for (u64 i = 0; i < Qr; i++) p2.aloOffs[i + 1] += p2.aloOffs[i];
    std::vector<u32> cur(p2.aloOffs.begin(), p2.aloOffs.end() - 1);
    p2.aloC.assign(all.size(), NumQ{});
    for (auto& C : all) p2.aloC[cur[key(C)]++] = C;
  }
  // ---- B1 / B2 ------------------------------------------------------------------------------
  {
    const u32 h1 = p.h1, hm = p.hm;
    p2.b1V.assign(p.nB1, 0);
    p2.recB1.assign(p.nB1, NumQ{});
    enum_digits(h1, P, [&](u64 idx, const std::vector<u32>& d) {
      u64 V = 0;
      NumQ R{};
      for (u32 j = 0; j < h1; j++) {
        V += pl.midV[(size_t)j * P + d[j]];
        if (V >= p.rho) V -= p.rho;
        if (d[j]) R = nq_add(R, pl.midW[(size_t)j * P + d[j]], p.rho);
      }
      if (V % p.qpow[p.r0] != 0) throw std::runtime_error("v2: B1 low digits not zero");
      p2.b1V[idx] = V;
      p2.recB1[idx] = R;
    });
    std::vector<u64> V2(p.nB2);
    u32 h2 = hm - h1;
    enum_digits(h2, P, [&](u64 idx, const std::vector<u32>& d) {
      u64 V = 0;
      for (u32 j = 0; j < h2; j++) {
        V += pl.midV[(size_t)(h1 + j) * P + d[j]];
        if (V >= p.rho) V -= p.rho;
      }
      for (u32 ci = 0; ci < p.NC; ci++) {
        u64 v = V;
        if (p.hasCenter) {
          v += pl.cenV[p.centerDig[ci]];
          if (v >= p.rho) v -= p.rho;
        }
        V2[idx * p.NC + ci] = v;
      }
    });
    auto key2 = [&](u64 v) { return (v / p.qpow[p.r0]) % Qr; };
    p2.b2Offs.assign(Qr + 1, 0);
    for (u64 v : V2) p2.b2Offs[key2(v) + 1]++;
    for (u64 i = 0; i < Qr; i++) p2.b2Offs[i + 1] += p2.b2Offs[i];
    std::vector<u32> cur(p2.b2Offs.begin(), p2.b2Offs.end() - 1);
    p2.b2V.assign(p.nB2, 0);
    p2.b2Idx.assign(p.nB2, 0);
    for (u64 b2 = 0; b2 < p.nB2; b2++) {
      u32 pos = cur[key2(V2[b2])]++;
      p2.b2V[pos] = V2[b2];
      p2.b2Idx[pos] = (u32)b2;
    }
  }
}

void print_plan2(const Plan2& p2, FILE* f) {
  if (!p2.ok) {
    fprintf(f, "  L=%u: v2 not applicable (%s)\n", p2.base.L, p2.why.c_str());
    return;
  }
  const Params& p = p2.base.prm;
  fprintf(f,
          "  L=%u v2: k=%u (hi %u + lo %u) m=%u middle=%u+%u pairs (+centre x%u)  |A|=%.3g |B|=%.3g\n"
          "        classes Q^r=%llu (r=%u, skip %u) w=%u/batch -> %llu batches of ~%.3g middles, key digits %u (+%u check)  est %.3gs\n",
          p2.base.L, p.k, p.khi, p.klo, p.m, p.h1, p.hm - p.h1, p.hasCenter ? p.NC : 0, (double)p2.nA,
          (double)p2.nB, (unsigned long long)p.Qr, p.r, p.r0, p.w, (unsigned long long)p2.nBatches, p2.batchEntriesEst,
          p.tp, p.v2u, p2.estSeconds);
}
