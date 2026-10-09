// table.h - the bucketed table of middle palindromes.
#pragma once
#include "plan.h"

// Bucket layout (CSR): bucket `key` holds ents[offs[key] .. offs[key+1]).
// Entry = (index of the middle << 32) | check, where
//   V     = M * P^k mod Q^c
//   key   = (V mod Q^t) / gcd(P^k, Q^t)
//   check = floor(V / Q^t) mod Q^ueff
struct TableView {
  u32* offs = nullptr;  // Dt + 1 entries
  u64* ents = nullptr;  // tableEntries entries
};

size_t table_offs_bytes(const Plan& pl);
size_t table_ents_bytes(const Plan& pl);

// Fill a table in caller-provided storage. Returns the build time in seconds.
double build_table(const Plan& pl, TableView tv, int threads);

// Enumerate every middle (index, V) for a plan; used by the builder and by tests.
template <class F>
void enum_middles_range(const Plan& pl, u64 prefix, u32 prefixDigits, F&& f);

#include <vector>
template <class F>
static void enum_rec(const Plan& pl, u32 j, u64 idx, u64 V, F& f) {
  const Params& p = pl.prm;
  if (j == p.hm) {
    if (p.hasCenter) {
      for (u32 ci = 0; ci < p.NC; ci++) {
        u64 v = V + pl.cenV[p.centerDig[ci]];
        if (v >= p.rho) v -= p.rho;
        f(idx * p.NC + ci, v);
      }
    } else {
      f(idx, V);
    }
    return;
  }
  for (u32 e = 0; e < p.P; e++) {
    u64 v = V + pl.midV[(size_t)j * p.P + e];
    if (v >= p.rho) v -= p.rho;
    enum_rec(pl, j + 1, idx * p.P + e, v, f);
  }
}

template <class F>
void enum_middles_range(const Plan& pl, u64 prefix, u32 prefixDigits, F&& f) {
  const Params& p = pl.prm;
  std::vector<u32> e(prefixDigits);
  u64 r = prefix;
  for (int j = (int)prefixDigits - 1; j >= 0; j--) {
    e[j] = (u32)(r % p.P);
    r /= p.P;
  }
  u64 V = 0;
  for (u32 j = 0; j < prefixDigits; j++) {
    V += pl.midV[(size_t)j * p.P + e[j]];
    if (V >= p.rho) V -= p.rho;
  }
  enum_rec(pl, prefixDigits, prefix, V, f);
}
