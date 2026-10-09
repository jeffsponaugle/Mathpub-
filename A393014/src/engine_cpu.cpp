// engine_cpu.cpp - multi-threaded CPU enumeration of the outer digits.
#include <chrono>
#include <memory>
#include <cmath>
#include <mutex>
#include <thread>

#include "engine.h"

bool node_viable(const Plan& pl, const NumQ& C, u32 depth) {
  const Params& p = pl.prm;
  NumQ lo = C, hi = nq_add(C, pl.WmaxAt[depth], p.rho);
  bool any = false;
  int segFound = -1, nsegHit = 0;
  for (u32 si = 0; si < p.nseg; si++) {
    const SegInfo& sg = p.seg[si];
    if (nq_cmp(hi, sg.lo) < 0) break;
    if (nq_cmp(lo, sg.hi) > 0) continue;
    nsegHit++;
    segFound = (int)si;
    if (sg.allowed) any = true;
  }
  if (!any) return false;
  if (!p.gmode || nsegHit != 1 || depth == 0) return true;
  // whole interval inside one base-Q length: compare residues modulo gcd(P^depth, Q^s)
  const SegInfo& sg = p.seg[segFound];
  u64 x = nq_top(lo, sg.ia, sg.oa, p), y = nq_top(hi, sg.ia, sg.oa, p);
  u64 T;
  u32 s = shared_prefix(x, y, sg.S, p, T);
  u64 g = pl.gtab[depth][s];
  if (g <= 1) return true;
  u64 R = rev_digits(T, s, p, pl.revtab.data());
  u64 c0 = C.l[0];
  u64 diff = (R >= c0) ? (R - c0) : (R + (p.rho - c0));
  return diff % g == 0;
}

static void dfs_prefix(const Plan& pl, u32 depth, u32 k0, const NumQ& C, u64& pruned,
                       const std::function<void(const NumQ&)>& emit) {
  if (depth == k0) {
    emit(C);
    return;
  }
  const Params& p = pl.prm;
  if (depth == 0) {
    for (u32 d : pl.topDigits) {
      NumQ C2 = nq_add(C, pl.outerW[d], p.rho);
      if (node_viable(pl, C2, 1))
        dfs_prefix(pl, 1, k0, C2, pruned, emit);
      else
        pruned++;
    }
    return;
  }
  for (u32 d = 0; d < p.P; d++) {
    NumQ C2 = d ? nq_add(C, pl.outerW[(size_t)depth * p.P + d], p.rho) : C;
    if (node_viable(pl, C2, depth + 1))
      dfs_prefix(pl, depth + 1, k0, C2, pruned, emit);
    else
      pruned++;
  }
}

u64 gen_prefixes(const Plan& pl, u32 k0, const std::function<void(const NumQ&)>& emit) {
  u64 pruned = 0;
  NumQ zero{};
  dfs_prefix(pl, 0, k0, zero, pruned, emit);
  return pruned;
}

// Depth-first walk below a prefix.  Leaves are processed in batches so that the bucket loads
// of many nodes are in flight at once (prefetch), instead of one cache miss at a time.
struct ReqBuf {
  static constexpr int MAXN = 64;           // nodes per batch
  static constexpr int MAXR = 4 * MAXN + 8; // pieces per batch
  ProbeReq rq[MAXR];
  u32 node[MAXR];
  u32 beg[MAXR], end[MAXR];
  u32 n = 0, cur = 0;
  u64 over = 0;
  void add(ProbeReq r) {
    if (r.flag == 2) {
      over++;
      return;
    }
    rq[n] = r;
    node[n] = cur;
    n++;
  }
};

struct Walker {
  const Plan& pl;
  const TableView& tv;
  CpuSink sink;
  u64 nodes = 0, pruned = 0;
  NumQ batch[ReqBuf::MAXN];
  u32 nb = 0;
  ReqBuf rb;

  void flush() {
    const Params& p = pl.prm;
    const u32* revtab = pl.revtab.data();
    rb.n = 0;
    for (u32 i = 0; i < nb; i++) {
      rb.cur = i;
      node_requests(batch[i], p, revtab, rb);
    }
    for (u32 r = 0; r < rb.n; r++) __builtin_prefetch(&tv.offs[rb.rq[r].key]);
    for (u32 r = 0; r < rb.n; r++) {
      rb.beg[r] = tv.offs[rb.rq[r].key];
      rb.end[r] = tv.offs[rb.rq[r].key + 1];
      __builtin_prefetch(&tv.ents[rb.beg[r]]);
    }
    for (u32 r = 0; r < rb.n; r++) scan_bucket(rb.rq[r], rb.beg[r], rb.end[r], batch[rb.node[r]], p, tv.ents, sink);
    sink.over += rb.over;
    rb.over = 0;
    nb = 0;
  }

  void walk(u32 depth, const NumQ& C) {
    const Params& p = pl.prm;
    if (depth == p.k) {
      nodes++;
      batch[nb++] = C;
      if (nb == ReqBuf::MAXN) flush();
      return;
    }
    const NumQ* w = &pl.outerW[(size_t)depth * p.P];
    for (u32 d = 0; d < p.P; d++) {
      NumQ C2 = d ? nq_add(C, w[d], p.rho) : C;
      if (depth + 1 < p.k && !node_viable(pl, C2, depth + 1)) {
        pruned++;
        continue;
      }
      walk(depth + 1, C2);
    }
  }
};

RunStats run_cpu(const Plan& pl, const TableView& tv, int threads, std::vector<NumQ>& found) {
  auto t0 = std::chrono::steady_clock::now();
  const Params& p = pl.prm;
  if (threads <= 0) threads = (int)std::thread::hardware_concurrency();
  RunStats st;
  // prefixes: deep enough for load balance
  u32 k0 = 1;
  double cnt = (double)pl.topDigits.size();
  while (k0 < p.k && cnt < 2000.0 * threads) {
    cnt *= p.P;
    k0++;
  }
  std::vector<NumQ> prefixes;
  u64 seen = 0;
  st.pruned += gen_prefixes(pl, k0, [&](const NumQ& C) {
    if (seen++ % pl.partN == pl.partI) prefixes.push_back(C);  // --part i/n: interleaved prefixes
  });

  std::atomic<size_t> next{0};
  std::mutex mu;
  std::vector<std::thread> th;
  for (int w = 0; w < threads; w++)
    th.emplace_back([&] {
      std::vector<NumQ> local;
      auto wk = std::make_unique<Walker>(Walker{pl, tv, CpuSink{0, 0, 0, &local, &pl}, 0, 0, {}, 0, {}});
      for (;;) {
        size_t i = next.fetch_add(1);
        if (i >= prefixes.size()) break;
        wk->walk(k0, prefixes[i]);
      }
      wk->flush();
      Walker& wkr = *wk;
      std::lock_guard<std::mutex> lk(mu);
      st.nodes += wkr.nodes;
      st.pruned += wkr.pruned;
      st.probes += wkr.sink.probes;
      st.cands += wkr.sink.cands;
      st.overflow += wkr.sink.over;
      found.insert(found.end(), local.begin(), local.end());
    });
  for (auto& t : th) t.join();
  st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return st;
}

// ---------------------------------------------------------------------------------------------
// Brute force: all base-P palindromes of length L (leading digit restricted as in the plan).
struct BruteCtx {
  const Plan& pl;
  std::vector<NumQ> w;  // [i*P + d]: d * (P^(L-1-i) + P^i), last position = centre P^(L/2) if odd
  u32 npos;             // number of free positions
  std::vector<NumQ>* found;
  u64 count = 0;
  void walk(u32 i, const NumQ& N) {
    const Params& p = pl.prm;
    if (i == npos) {
      count++;
      if (nq_is_pal(N, p, pl.revtab.data())) found->push_back(N);
      return;
    }
    u32 d0 = (i == 0) ? 1 : 0;
    for (u32 d = d0; d < p.P; d++)  // no digit filters: this is the independent reference
      walk(i + 1, d ? nq_add(N, w[(size_t)i * p.P + d], p.rho) : N);
  }
};

RunStats run_brute(const Plan& pl, int threads, std::vector<NumQ>& found) {
  auto t0 = std::chrono::steady_clock::now();
  const Params& p = pl.prm;
  if (threads <= 0) threads = (int)std::thread::hardware_concurrency();
  RunStats st;
  u32 L = pl.L, P = p.P;
  u32 npos = (L + 1) / 2;
  std::vector<NumQ> w((size_t)npos * P);
  for (u32 i = 0; i < npos; i++) {
    mpz_class a, b;
    mpz_ui_pow_ui(a.get_mpz_t(), P, L - 1 - i);
    mpz_ui_pow_ui(b.get_mpz_t(), P, i);
    mpz_class wi = (L - 1 - i == i) ? a : a + b;
    for (u32 d = 0; d < P; d++) w[(size_t)i * P + d] = to_numq(wi * d, p);
  }
  // parallelise over the first two positions
  std::vector<std::pair<u32, u32>> jobs;
  for (u32 d0 = (L == 1 ? 0 : 1); d0 < P; d0++) {  // every leading digit (no plan filters)
    if (npos >= 2)
      for (u32 d1 = 0; d1 < P; d1++) jobs.push_back({d0, d1});
    else
      jobs.push_back({d0, ~0u});
  }
  std::atomic<size_t> next{0};
  std::mutex mu;
  std::vector<std::thread> th;
  for (int t = 0; t < threads; t++)
    th.emplace_back([&] {
      std::vector<NumQ> local;
      BruteCtx bc{pl, w, npos, &local};
      for (;;) {
        size_t j = next.fetch_add(1);
        if (j >= jobs.size()) break;
        auto [d0, d1] = jobs[j];
        NumQ N{};
        if (d0) N = nq_add(N, w[d0], p.rho);
        if (d1 == ~0u) {
          bc.walk(1, N);
        } else {
          if (d1) N = nq_add(N, w[(size_t)P + d1], p.rho);
          bc.walk(2, N);
        }
      }
      std::lock_guard<std::mutex> lk(mu);
      st.nodes += bc.count;
      found.insert(found.end(), local.begin(), local.end());
    });
  for (auto& t : th) t.join();
  st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return st;
}
