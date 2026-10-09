// plan.cpp - parameter selection and precomputation (GMP on the host).
#include "plan.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>

static u64 gcd64(u64 a, u64 b) { return std::gcd(a, b); }

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif
#include <unistd.h>
u64 physical_memory_bytes() {
#ifdef __APPLE__
  u64 mem = 0;
  size_t len = sizeof(mem);
  if (sysctlbyname("hw.memsize", &mem, &len, nullptr, 0) == 0) return mem;
#endif
  long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGE_SIZE);
  return (pages > 0 && psz > 0) ? (u64)pages * (u64)psz : (u64)8 << 30;
}

u32 digits_per_limb(u32 Q) {
  // largest c with Q^c <= 2^63
  unsigned __int128 v = 1;
  u32 c = 0;
  while (v * Q <= ((unsigned __int128)1 << 63)) {
    v *= Q;
    c++;
  }
  return c;
}

FastDiv make_fastdiv(u64 d) {
  FastDiv f{};
  if (d == 0) throw std::runtime_error("division by zero");
  if ((d & (d - 1)) == 0) {
    f.kind = 0;
    f.sh = (u32)__builtin_ctzll(d);
    f.m = 0;
    return f;
  }
  u32 l = 64 - (u32)__builtin_clzll(d - 1);  // ceil(log2 d)
  // m = ceil(2^(63+l) / d)  (< 2^64 for non powers of two)
  unsigned __int128 num = (unsigned __int128)1 << (63 + l);
  unsigned __int128 m = (num + d - 1) / d;
  f.kind = 1;
  f.m = (u64)m;
  f.sh = l - 1;
  return f;
}

NumQ to_numq(const mpz_class& x, const Params& p) {
  NumQ r{};
  mpz_class t = x;
  if (t < 0) throw std::runtime_error("to_numq: negative");
  for (int i = 0; i < NLIMB; i++) r.l[i] = mpz_fdiv_q_ui(t.get_mpz_t(), t.get_mpz_t(), p.rho);
  if (t != 0)
    throw std::runtime_error("number too large for " + std::to_string(NLIMB) +
                             " limbs (search beyond ~10^54: rebuild with make NLIMB=4)");
  return r;
}

mpz_class from_numq(const NumQ& x, const Params& p) {
  mpz_class r = 0, rho;
  mpz_set_ui(rho.get_mpz_t(), p.rho);
  for (int i = NLIMB - 1; i >= 0; i--) {
    r *= rho;
    mpz_class li;
    mpz_set_ui(li.get_mpz_t(), x.l[i]);
    r += li;
  }
  return r;
}

std::string to_base(const mpz_class& x, u32 b) {
  if (b <= 62) {
    std::vector<char> buf(mpz_sizeinbase(x.get_mpz_t(), b) + 2);
    mpz_get_str(buf.data(), (int)b, x.get_mpz_t());
    return std::string(buf.data());
  }
  return "?";
}

u32 ndigits_base(const mpz_class& x, u32 b) {
  if (x == 0) return 1;
  // mpz_sizeinbase may overestimate by one for non powers of 2
  u32 n = (u32)mpz_sizeinbase(x.get_mpz_t(), b);
  mpz_class pw;
  mpz_ui_pow_ui(pw.get_mpz_t(), b, n - 1);
  if (pw > x) n--;
  return n;
}

bool is_pal_base(const mpz_class& x, u32 b) {
  if (x < 0) return false;
  std::vector<u32> d;
  mpz_class t = x;
  if (t == 0) return true;
  while (t > 0) d.push_back((u32)mpz_fdiv_q_ui(t.get_mpz_t(), t.get_mpz_t(), b));
  for (size_t i = 0, j = d.size() - 1; i < j; i++, j--)
    if (d[i] != d[j]) return false;
  return true;
}

mpz_class mpow(u32 b, u32 e) {
  mpz_class r;
  mpz_ui_pow_ui(r.get_mpz_t(), b, e);
  return r;
}

void init_radix(Params& p, u32 Q, std::vector<u32>& revtab) {
  p.Q = Q;
  p.c = digits_per_limb(Q);
  if (const char* e = std::getenv("DBPAL_LIMB_DIGITS"))  // testing only: small limbs use every limb at small N
    p.c = std::max(1u, std::min(p.c, (u32)std::atoi(e)));
  p.qbits = ((Q & (Q - 1)) == 0) ? (u32)__builtin_ctz(Q) : 0;
  u64 v = 1;
  for (u32 j = 0; j < MAXQP; j++) {
    if (j <= p.c) {
      p.qpow[j] = v;
      p.qdiv[j] = make_fastdiv(v);
      if (j < p.c) v *= Q;
    } else {
      p.qpow[j] = ~0ull;  // sentinel: larger than any value we compare against
      p.qdiv[j] = make_fastdiv(1ull << 63);
    }
  }
  p.rho = p.qpow[p.c];
  // ndlo[b] = number of base-Q digits of 2^(b-1)
  for (u32 b = 1; b <= 64; b++) {
    u64 x = 1ull << (b - 1);
    u32 nd = 0;
    unsigned __int128 w = 1;
    while (w <= x) {
      w *= Q;
      nd++;
    }
    p.ndlo[b] = nd;
  }
  p.ndlo[0] = 0;
  p.ndlo[65] = p.ndlo[66] = p.ndlo[67] = 0;
  // reversal chunk: largest r with Q^r <= 65536
  u32 r = 1;
  while (r < p.c && (u64)p.qpow[r + 1] <= 65536) r++;
  p.revR = r;
  u32 n = (u32)p.qpow[r];
  revtab.assign(n, 0);
  for (u32 x = 0; x < n; x++) {
    u32 y = x, rv = 0;
    for (u32 i = 0; i < r; i++) {
      rv = rv * Q + y % Q;
      y /= Q;
    }
    revtab[x] = rv;
  }
}

// Number of digits per check field: largest u with Q^u <= 2^32 - 1.
u32 check_digits(u32 Q) {
  u64 v = 1;
  u32 u = 0;
  while (v * Q <= 0xffffffffull) {
    v *= Q;
    u++;
  }
  return u;
}

void fill_seginfo(SegInfo& sg, u32 LQ, const Params& p) {
  mpz_class lo = mpow(p.Q, LQ - 1), hi = mpow(p.Q, LQ) - 1;
  sg.lo = to_numq(lo, p);
  sg.hi = to_numq(hi, p);
  sg.LQ = LQ;
  sg.S = std::min(p.c, LQ);
  u32 a = LQ - sg.S;
  sg.ia = a / p.c;
  sg.oa = a % p.c;
  u32 z = (LQ >= p.t) ? LQ - p.t : 0;
  sg.iz = z / p.c;
  sg.oz = z % p.c;
}

Plan make_plan(const SearchOpts& o, u32 L) {
  Plan pl;
  pl.L = L;
  pl.partI = o.partI;
  pl.partN = std::max<u32>(1, o.partN);
  const u32 P = o.P, Q = o.Q;
  Params& p = pl.prm;
  init_radix(p, Q, pl.revtab);
  p.P = P;
  p.gmode = (gcd64(P, Q) > 1) ? 1 : 0;
  pl.lowN = mpow(P, L - 1);
  pl.highN = mpow(P, L) - 1;

  // ---- parity rules -------------------------------------------------------------------
  // even length palindromes in base b are divisible by b+1; the lowest digit equals the
  // (nonzero) leading digit, so N is never divisible by the base.
  if (L % 2 == 0 && L >= 2 && (P + 1) % Q == 0) {
    pl.skip = true;
    pl.skipReason = "even base-P length => N divisible by P+1, a multiple of Q";
    return pl;
  }
  if (o.primesOnly && L % 2 == 0 && L > 2) {
    pl.skip = true;
    pl.skipReason = "even base-P length => N divisible by P+1 (not prime)";
    return pl;
  }

  u32 LQlo = ndigits_base(pl.lowN, Q), LQhi = ndigits_base(pl.highN, Q);
  if (LQhi - LQlo + 1 > MAXSEG) throw std::runtime_error("too many base-Q lengths per base-P length");

  // allowed centre digits: Q = 2 => N odd (lowest bit = leading bit); P odd => digit sum odd => centre odd.
  // (Only for Q = 2: in base 4, 8, 16 the lowest digit equals the leading digit but may be even.)
  std::vector<u32> centers;
  bool oddOnly = (Q == 2) && (P % 2 == 1) && (L % 2 == 1);
  for (u32 d = 0; d < P; d++)
    if (!oddOnly || (d % 2 == 1)) centers.push_back(d);
  // allowed leading digits: Q = 2 and P even => N odd => lowest P-digit (= leading digit) odd
  pl.topDigits.clear();
  for (u32 d = 1; d < P; d++)
    if (!(Q == 2 && P % 2 == 0 && d % 2 == 0)) pl.topDigits.push_back(d);

  // ---- small lengths: brute force ----------------------------------------------------
  double halfCount = std::pow((double)P, (L + 1) / 2);
  if (L < 5 || (halfCount < 4e6 && o.forceK < 0)) {
    pl.brute = true;
    return pl;
  }

  // ---- choose k and t ------------------------------------------------------------------
  const u32 umax = check_digits(Q);
  SearchOpts o2 = o;
  if (o2.tableMax == 0) {
    // ~40% of physical memory at ~9 bytes per entry, capped by 32-bit offsets
    double mem = (double)physical_memory_bytes();
    o2.tableMax = (u64)std::min(0.40 * mem / 9.0, 3.9e9);
  }
  const SearchOpts& oo = o2;
  struct Choice {
    u32 k, t;
    double cost, nodes, cands;
    u64 table, Dt;
  };
  std::vector<Choice> choices;
  for (u32 k = 1; 2 * k < L; k++) {
    if (o.forceK >= 0 && (int)k != o.forceK) continue;
    u32 m = L - 2 * k;
    u32 hm = m / 2;
    bool hasC = m % 2 == 1;
    double tab = std::pow((double)P, hm) * (hasC ? (double)centers.size() : 1.0);
    if (!o.noTableLimit && (tab > (double)oo.tableMax || tab > 4.0e9)) continue;
    if (k > 40 || hm > 40) continue;
    mpz_class W = (mpow(P, m) - 1) * mpow(P, k);
    u32 dW = ndigits_base(W, Q);
    // smallest base-Q length that is actually searched
    u32 LQmin = 0;
    for (u32 LQ = LQlo; LQ <= LQhi; LQ++) {
      bool allowed = true;
      if (LQ % 2 == 0 && LQ >= 2 && (Q + 1) % P == 0) allowed = false;
      if (o.primesOnly && LQ % 2 == 0 && LQ > 2) allowed = false;
      if (allowed) {
        LQmin = LQ;
        break;
      }
    }
    if (LQmin == 0) {
      pl.skip = true;
      pl.skipReason = "no admissible base-Q length";
      return pl;
    }
    int tmax = (int)LQmin - (int)dW;  // => node interval width < Q^(LQ - t): at most 2 pieces
    tmax = std::min(tmax, (int)p.c - 1);
    if (tmax < 0) tmax = 0;
    mpz_class Pk = mpow(P, k);
    auto Dt_of = [&](u32 t) -> double {
      mpz_class Qt = mpow(Q, t), g;
      mpz_gcd(g.get_mpz_t(), Pk.get_mpz_t(), Qt.get_mpz_t());
      mpz_class D = Qt / g;
      return D.get_d();
    };
    u32 t = 0;
    if (o.forceT >= 0) {
      t = (u32)std::min(o.forceT, tmax);
    } else {
      // largest t <= tmax with Dt <= max(tab/2, 1) (buckets of ~2+ entries keep the offset array small)
      for (int tt = tmax; tt >= 0; tt--) {
        if (Dt_of(tt) <= std::max(tab / 2.0, 1.0) && Dt_of(tt) < 4.0e9) {
          t = (u32)tt;
          break;
        }
      }
    }
    double Dt = Dt_of(t);
    double B = tab / Dt;
    double nodes = (double)pl.topDigits.size() * std::pow((double)P, k - 1);
    if (p.gmode) {
      // consistency pruning keeps roughly P^k / gcd(P^k, Q^s) nodes
      double s = (double)k * std::log((double)P) / std::log((double)Q);
      mpz_class Qs = mpow(Q, (u32)s), g;
      mpz_gcd(g.get_mpz_t(), Pk.get_mpz_t(), Qs.get_mpz_t());
      nodes /= std::max(1.0, g.get_d());
    }
    // candidates: entries agreeing with a node on all s known digits ~ nodes * table / Q^s
    double sAvg = 0.5 * (double)(LQmin + LQhi) - std::log(W.get_d()) / std::log((double)Q) - 1.0;
    double cands = nodes * std::min(B, tab * std::pow((double)Q, -sAvg));
    // estimated seconds (calibrated on an M4 Max: 32-core GPU / 14 CPU threads)
    double bytes = 8.0 * tab + 4.0 * Dt;
    double tNode, tCand;
    if (o.gpu) {
      tNode = bytes <= 48e6 ? 0.18 : bytes <= 512e6 ? 0.45 : bytes <= 4e9 ? 0.95 : 1.2;
      tCand = 0.25;
    } else {
      tNode = 2.65 * 14.0 / std::max(1, o.threads);
      tCand = 1.83 * 14.0 / std::max(1, o.threads);
    }
    tNode += 0.03 * B;
    double cost = 1e-9 * (nodes * tNode + cands * tCand + 5.3 * tab + 1.0 * Dt);
    choices.push_back({k, t, cost, nodes, cands, (u64)tab, (u64)Dt});
  }
  if (choices.empty()) throw std::runtime_error("no feasible depth k (increase --table-max)");
  Choice best = *std::min_element(choices.begin(), choices.end(),
                                  [](const Choice& a, const Choice& b) { return a.cost < b.cost; });

  // ---- fill Params ------------------------------------------------------------------------
  const u32 k = best.k, m = L - 2 * k, hm = m / 2;
  p.k = k;
  p.m = m;
  p.hm = hm;
  p.hasCenter = m % 2;
  p.NC = p.hasCenter ? (u32)centers.size() : 1;
  if (p.NC > MAXCEN) throw std::runtime_error("too many centre digits");
  for (u32 i = 0; i < MAXCEN; i++) p.centerDig[i] = (i < centers.size()) ? centers[i] : 0;
  p.t = best.t;
  p.ueff = std::min(umax, p.c - p.t);
  p.Dt = best.Dt;
  pl.tableEntries = best.table;
  pl.nodesEst = best.nodes;
  pl.candsEst = best.cands;
  pl.costEst = best.cost;

  mpz_class Pk = mpow(P, k);
  {
    mpz_class Qt = mpow(Q, p.t), g;
    mpz_gcd(g.get_mpz_t(), Pk.get_mpz_t(), Qt.get_mpz_t());
    p.gt = g.get_ui();
    p.gtdiv = make_fastdiv(p.gt);
    mpz_class D = Qt / g;
    if (D.get_ui() != p.Dt) throw std::runtime_error("Dt mismatch");
  }
  for (u32 s = 0; s < MAXQP; s++) {
    if (s <= p.c) {
      mpz_class Qs = mpow(Q, s), g;
      mpz_gcd(g.get_mpz_t(), Pk.get_mpz_t(), Qs.get_mpz_t());
      p.gs[s] = g.get_ui();
    } else {
      p.gs[s] = 1;
    }
    p.gsdiv[s] = make_fastdiv(p.gs[s]);
  }

  // segments
  p.nseg = 0;
  for (u32 LQ = LQlo; LQ <= LQhi; LQ++) {
    SegInfo& sg = p.seg[p.nseg++];
    memset(&sg, 0, sizeof(sg));
    fill_seginfo(sg, LQ, p);
    bool allowed = true;
    if (LQ % 2 == 0 && LQ >= 2 && (Q + 1) % P == 0) allowed = false;
    if (o.primesOnly && LQ % 2 == 0 && LQ > 2) allowed = false;
    sg.allowed = allowed ? 1 : 0;
    if (p.t > sg.S) throw std::runtime_error("t larger than extracted digits");
  }

  // weights
  mpz_class rhoz;
  mpz_set_ui(rhoz.get_mpz_t(), p.rho);
  pl.outerW.assign((size_t)k * P, NumQ{});
  for (u32 i = 0; i < k; i++) {
    mpz_class w = mpow(P, L - 1 - i) + mpow(P, i);
    for (u32 d = 0; d < P; d++) pl.outerW[(size_t)i * P + d] = to_numq(w * d, p);
  }
  pl.midW.assign((size_t)std::max(hm, 1u) * P, NumQ{});
  pl.midV.assign((size_t)std::max(hm, 1u) * P, 0);
  for (u32 j = 0; j < hm; j++) {
    mpz_class u = mpow(P, j) + mpow(P, m - 1 - j);
    mpz_class us = u * Pk;
    mpz_class uv = us % rhoz;
    for (u32 e = 0; e < P; e++) {
      pl.midW[(size_t)j * P + e] = to_numq(us * e, p);
      mpz_class v = (uv * e) % rhoz;
      pl.midV[(size_t)j * P + e] = v.get_ui();
    }
  }
  pl.cenW.assign(P, NumQ{});
  pl.cenV.assign(P, 0);
  if (p.hasCenter) {
    mpz_class u = mpow(P, hm) * Pk;
    for (u32 d = 0; d < P; d++) {
      pl.cenW[d] = to_numq(u * d, p);
      mpz_class v = (u * d) % rhoz;
      pl.cenV[d] = v.get_ui();
    }
  }
  // reconstruction tables: index = hi * recLoN + lo, lo = (last g2 middle digits, centre)
  {
    u32 g2best = 0;
    double best_sz = 1e300;
    for (u32 g2 = 0; g2 <= hm; g2++) {
      double sz = std::pow((double)P, hm - g2) + std::pow((double)P, g2) * p.NC;
      if (sz < best_sz) {
        best_sz = sz;
        g2best = g2;
      }
    }
    u32 g1 = hm - g2best, g2 = g2best;
    u64 nHi = 1, nLo = p.NC;
    for (u32 i = 0; i < g1; i++) nHi *= P;
    for (u32 i = 0; i < g2; i++) nLo *= P;
    p.recLoN = (u32)nLo;
    p.recDiv = make_fastdiv(nLo);
    pl.recHi.assign(nHi, NumQ{});
    pl.recLo.assign(nLo, NumQ{});
    for (u64 x = 0; x < nHi; x++) {  // digits e_0..e_{g1-1}, e_0 most significant
      NumQ acc{};
      u64 r = x;
      for (int j = (int)g1 - 1; j >= 0; j--) {
        u32 e = (u32)(r % P);
        r /= P;
        if (e) acc = nq_add(acc, pl.midW[(size_t)j * P + e], p.rho);
      }
      pl.recHi[x] = acc;
    }
    for (u64 x = 0; x < nLo; x++) {  // (e_g1..e_{hm-1}) * NC + ci
      NumQ acc{};
      u64 r = x;
      if (p.hasCenter) {
        u32 ci = (u32)(r % p.NC);
        r /= p.NC;
        u32 cd = p.centerDig[ci];
        if (cd) acc = nq_add(acc, pl.cenW[cd], p.rho);
      }
      for (int j = (int)hm - 1; j >= (int)g1; j--) {
        u32 e = (u32)(r % P);
        r /= P;
        if (e) acc = nq_add(acc, pl.midW[(size_t)j * P + e], p.rho);
      }
      pl.recLo[x] = acc;
    }
  }
  pl.WmaxAt.assign(k + 1, NumQ{});
  for (u32 i = 0; i <= k; i++) pl.WmaxAt[i] = to_numq((mpow(P, L - 2 * i) - 1) * mpow(P, i), p);
  p.Wmax = pl.WmaxAt[k];

  if (p.gmode) {
    pl.gtab.assign(k + 1, std::vector<u64>(p.c + 1, 1));
    for (u32 i = 0; i <= k; i++)
      for (u32 s = 0; s <= p.c; s++) {
        mpz_class a = mpow(P, i), b = mpow(Q, s), g;
        mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        pl.gtab[i][s] = g.get_ui();
      }
  }
  return pl;
}

void print_plan(const Plan& pl, FILE* f) {
  const Params& p = pl.prm;
  if (pl.skip) {
    fprintf(f, "  L=%u: skipped (%s)\n", pl.L, pl.skipReason.c_str());
    return;
  }
  if (pl.brute) {
    fprintf(f, "  L=%u: brute force over %.3g palindromes\n", pl.L, std::pow((double)p.P, (pl.L + 1) / 2));
    return;
  }
  fprintf(f, "  L=%u: k=%u m=%u (hm=%u, centre=%u x%u) table=%llu buckets=%llu (t=%u, B=%.2f) check=%u digits\n",
          pl.L, p.k, p.m, p.hm, p.hasCenter, p.NC, (unsigned long long)pl.tableEntries,
          (unsigned long long)p.Dt, p.t, (double)pl.tableEntries / (double)p.Dt, p.ueff);
  fprintf(f, "        nodes~%.3g cands~%.3g est %.3gs  gmode=%u  base-%u lengths:", pl.nodesEst, pl.candsEst, pl.costEst,
          p.gmode, p.Q);
  for (u32 i = 0; i < p.nseg; i++) fprintf(f, " %u%s", p.seg[i].LQ, p.seg[i].allowed ? "" : "(skip)");
  fprintf(f, "\n");
}
