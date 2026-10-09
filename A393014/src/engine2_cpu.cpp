// engine2_cpu.cpp - v2 (class-partitioned) search on the CPU.  Each worker takes whole batches of
// classes [c0, c0 + w): generate the middles of those classes, bucket them, then walk the outer
// parts of the same classes and probe.
#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>

#include "engine2.h"

struct CpuSink2 {
  u64 probes = 0, cands = 0, pre = 0, over = 0;
  const Plan2* p2;
  std::vector<NumQ>* found;
  const NumQ* recB1;
  void overflow() { over++; }
  void candidate(NumQ C, u64 idxB) {
    const Params& p = p2->base.prm;
    NumQ N = v2_rebuild(C, idxB, p, recB1, p2->base.midW.data(), p2->base.cenW.data());
    if (nq_is_pal(N, p, p2->base.revtab.data())) found->push_back(N);
  }
};

// cyclic residue range [lo, lo + w) of a CSR keyed by x mod Qr -> up to two index ranges
static inline void csr_range(const std::vector<u32>& offs, u64 Qr, u64 lo, u64 w, u32 r[4]) {
  if (lo + w <= Qr) {
    r[0] = offs[lo];
    r[1] = offs[lo + w];
    r[2] = r[3] = 0;
  } else {
    r[0] = offs[lo];
    r[1] = offs[Qr];
    r[2] = offs[0];
    r[3] = offs[lo + w - Qr];
  }
}

RunStats run_cpu2(const Plan2& p2, int threads, std::vector<NumQ>& found) {
  auto t0 = std::chrono::steady_clock::now();
  const Params& p = p2.base.prm;
  const u64 Qr = p.Qr, w = p.w;
  const u64 nb = (Qr + w - 1) / w;
  std::vector<u64> mine;  // --part i/n: interleaved batches
  for (u64 bi = p2.base.partI; bi < nb; bi += p2.base.partN) mine.push_back(bi);
  const u64 nBuckets = w * p.qpow[p.tp];
  if (threads <= 0) threads = (int)std::thread::hardware_concurrency();
  RunStats st;
  std::atomic<u64> next{0};
  std::mutex mu;
  std::vector<std::thread> th;
  for (int t = 0; t < threads; t++)
    th.emplace_back([&] {
      std::vector<u32> offs(nBuckets + 1);
      std::vector<u64> ents;
      std::vector<NumQ> local;
      CpuSink2 sink{0, 0, 0, 0, &p2, &local, p2.recB1.data()};
      u64 nodes = 0;
      for (;;) {
        u64 j = next.fetch_add(1);
        if (j >= mine.size()) break;
        u64 bi = mine[j];
        u64 c0 = bi * w, ww = std::min(w, Qr - c0);
        // ---- middles of classes [c0, c0 + ww): count, scan, scatter
        std::fill(offs.begin(), offs.end(), 0);
        auto for_b = [&](auto&& f) {
          for (u32 b1 = 0; b1 < p.nB1; b1++) {
            u64 V1 = p2.b1V[b1];
            u64 lo = (c0 + Qr - (V1 / p.qpow[p.r0]) % Qr) % Qr;
            u32 rg[4];
            csr_range(p2.b2Offs, Qr, lo, ww, rg);
            for (int part = 0; part < 2; part++)
              for (u32 i = rg[2 * part]; i < rg[2 * part + 1]; i++) {
                u64 V = V1 + p2.b2V[i];
                if (V >= p.rho) V -= p.rho;
                f(V, (u64)b1 * p.nB2 + p2.b2Idx[i]);
              }
          }
        };
        for_b([&](u64 V, u64) { offs[v2_key(V, c0, p) + 1]++; });
        for (u64 i = 0; i < nBuckets; i++) offs[i + 1] += offs[i];
        ents.resize(offs[nBuckets]);
        for_b([&](u64 V, u64 idx) { ents[offs[v2_key(V, c0, p)]++] = v2_entry(V, idx, p); });
        for (u64 i = nBuckets; i > 0; i--) offs[i] = offs[i - 1];
        offs[0] = 0;
        // ---- outer parts of the same classes
        for (const AhiPiece& pc : p2.ahi) {
          u64 base = (pc.Rr + Qr - c0) % Qr;  // C_lo key for class c0
          u64 lo = (base + Qr - (ww - 1)) % Qr;                    // ... for class c0 + ww - 1
          u32 rg[4];
          csr_range(p2.aloOffs, Qr, lo, ww, rg);
          for (int part = 0; part < 2; part++)
            for (u32 i = rg[2 * part]; i < rg[2 * part + 1]; i++) {
              NumQ C = nq_add(pc.C, p2.aloC[i], p.rho);
              nodes++;
              v2_node(C, pc.a, pc.b, pc.seg, c0, p, p2.base.revtab.data(), offs.data(), ents.data(), sink);
            }
        }
      }
      std::lock_guard<std::mutex> lk(mu);
      st.nodes += nodes;
      st.probes += sink.probes;
      st.cands += sink.cands;
      st.overflow += sink.over;
      found.insert(found.end(), local.begin(), local.end());
    });
  for (auto& t : th) t.join();
  st.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return st;
}
