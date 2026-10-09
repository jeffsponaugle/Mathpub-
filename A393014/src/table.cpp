// table.cpp - multi-threaded counting-sort construction of the middle table.
#include "table.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <thread>

size_t table_offs_bytes(const Plan& pl) { return (size_t)(pl.prm.Dt + 1) * sizeof(u32); }
size_t table_ents_bytes(const Plan& pl) { return (size_t)pl.tableEntries * sizeof(u64); }

static inline void key_check(const Params& p, u64 V, u64& key, u32& chk) {
  u64 qt = fdiv(V, p.qdiv[p.t]);
  u64 lowt = V - qt * p.qpow[p.t];
  key = p.gmode ? fdiv(lowt, p.gtdiv) : lowt;
  chk = (u32)(qt - fdiv(qt, p.qdiv[p.ueff]) * p.qpow[p.ueff]);
}

template <class Body>
static void parallel_chunks(u64 nchunks, int threads, Body body) {
  std::atomic<u64> next{0};
  std::vector<std::thread> th;
  for (int w = 0; w < threads; w++)
    th.emplace_back([&] {
      for (;;) {
        u64 c = next.fetch_add(1);
        if (c >= nchunks) break;
        body(c);
      }
    });
  for (auto& t : th) t.join();
}

double build_table(const Plan& pl, TableView tv, int threads) {
  auto t0 = std::chrono::steady_clock::now();
  const Params& p = pl.prm;
  if (threads <= 0) threads = (int)std::thread::hardware_concurrency();
  const u64 Dt = p.Dt;
  if (pl.tableEntries >= 0xffffffffull) throw std::runtime_error("table too large for 32-bit offsets");
  memset(tv.offs, 0, (size_t)(Dt + 1) * sizeof(u32));

  // split the middle digits into chunks over the leading pd digits
  u32 pd = 0;
  u64 nchunks = 1;
  while (pd < p.hm && nchunks < 4096) {
    nchunks *= p.P;
    pd++;
  }

  // pass 1: counts in offs[key + 1]
  parallel_chunks(nchunks, threads, [&](u64 c) {
    enum_middles_range(pl, c, pd, [&](u64, u64 V) {
      u64 key;
      u32 chk;
      key_check(p, V, key, chk);
      __atomic_fetch_add(&tv.offs[key + 1], 1u, __ATOMIC_RELAXED);
    });
  });
  // exclusive prefix sum: offs[key] = start of bucket key
  {
    u64 nb = Dt + 1;
    int T = threads;
    std::vector<u64> part(T + 1, 0);
    u64 per = (nb + T - 1) / T;
    std::vector<std::thread> th;
    for (int w = 0; w < T; w++)
      th.emplace_back([&, w] {
        u64 a = std::min(nb, (u64)w * per), b = std::min(nb, (u64)(w + 1) * per);
        u64 s = 0;
        for (u64 i = a; i < b; i++) s += tv.offs[i];
        part[w + 1] = s;
      });
    for (auto& t : th) t.join();
    th.clear();
    for (int w = 0; w < T; w++) part[w + 1] += part[w];
    if (part[T] != pl.tableEntries) throw std::runtime_error("table count mismatch");
    for (int w = 0; w < T; w++)
      th.emplace_back([&, w] {
        u64 a = std::min(nb, (u64)w * per), b = std::min(nb, (u64)(w + 1) * per);
        u64 s = part[w];
        for (u64 i = a; i < b; i++) {
          s += tv.offs[i];
          tv.offs[i] = (u32)s;
        }
      });
    for (auto& t : th) t.join();
    // offs[i] now holds sum_{j<=i} count[j] = start of bucket i (counts were stored at key+1)
  }
  // pass 2: scatter; offs[key] advances to the end of bucket key
  parallel_chunks(nchunks, threads, [&](u64 c) {
    enum_middles_range(pl, c, pd, [&](u64 idx, u64 V) {
      u64 key;
      u32 chk;
      key_check(p, V, key, chk);
      u32 pos = __atomic_fetch_add(&tv.offs[key], 1u, __ATOMIC_RELAXED);
      tv.ents[pos] = (idx << 32) | chk;
    });
  });
  // offs[key] == start of bucket key+1: shift right to restore starts
  memmove(tv.offs + 1, tv.offs, (size_t)Dt * sizeof(u32));
  tv.offs[0] = 0;
  if (tv.offs[Dt] != pl.tableEntries) throw std::runtime_error("table fill mismatch");
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}
