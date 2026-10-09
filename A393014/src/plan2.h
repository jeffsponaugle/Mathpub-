// plan2.h - v2: residue-class partitioned meet in the middle (coprime bases).
//
// Time ~ |A| + |B| ~ N^(1/4) with memory ~ N^(1/8): instead of one table holding every middle,
// the middles are generated class by class (V(B) mod Q^r) and each batch of classes gets a small,
// cache-resident bucket table; the outer parts A = (A_hi, A_lo) of the same classes are found by a
// lookup on C_lo mod Q^r.  See README.md.
#pragma once
#include "plan.h"

struct Plan2 {
  Plan base;  // radix, segments (t = r + tp), weights, Wmax = node width at depth k, v2 fields in prm
  std::vector<AhiPiece> ahi;  // A_hi pieces
  std::vector<NumQ> aloC;     // C_lo values in CSR order (key C_lo mod Q^r)
  std::vector<u32> aloOffs;   // Q^r + 1
  std::vector<u64> b1V;       // V1 = M1 * P^k mod Q^c per b1
  std::vector<NumQ> recB1;    // M1 * P^k per b1
  std::vector<u64> b2V;       // V2 in CSR order (key V2 mod Q^r)
  std::vector<u32> b2Idx;     // b2 index in CSR order
  std::vector<u32> b2Offs;    // Q^r + 1
  u64 nA = 0, nB = 0, nAlo = 0, nBatches = 0;
  double estSeconds = 0, batchEntriesEst = 0;
  bool ok = false;
  std::string why;
};

struct V2Opts {
  int k = -1, khi = -1, h1 = -1, r = -1, tp = -1;
  double batchEntries = 4e6;  // target entries per batch table
};

// Choose parameters (or honour overrides) and estimate the cost; no big allocations.
Plan2 plan2_choose(const SearchOpts& o, const V2Opts& vo, u32 L);
// Build the A_hi / A_lo / B1 / B2 structures for a chosen plan.
void plan2_build(Plan2& p2, int threads);
void print_plan2(const Plan2& p2, FILE* f);
