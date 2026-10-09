// engine.h - node enumeration (CPU) and the common result type.
#pragma once
#include <atomic>
#include <functional>
#include <vector>

#include "plan.h"
#include "table.h"

struct RunStats {
  u64 nodes = 0;      // leaves (depth-k nodes) handed to process_node
  u64 probes = 0;     // bucket probes
  u64 cands = 0;      // entries passing the check digits (full palindrome test)
  u64 pruned = 0;     // interior nodes cut by the admissibility / gcd checks
  u64 overflow = 0;   // pieces with s < t (should stay 0)
  double seconds = 0;
};

struct CpuSink {
  u64 probes = 0, cands = 0, over = 0;
  std::vector<NumQ>* found;
  const Plan* pl;
  void emit(NumQ n) { found->push_back(n); }
  void candidate(NumQ C, u32 idx) {
    check_candidate(C, idx, pl->prm, pl->revtab.data(), pl->recHi.data(), pl->recLo.data(), *this);
  }
};

// Can a depth-i node with outer contribution C still contain a solution?
// (some admissible base-Q length must meet its interval and, when gcd(P,Q) > 1,
//  the low digits fixed by C must agree with the reversed top digits)
bool node_viable(const Plan& pl, const NumQ& C, u32 depth);

// Generate the viable nodes at depth k0 (prefixes for the workers / the GPU).
// Calls emit(C) for each, in order.  Returns the number of pruned interior nodes.
u64 gen_prefixes(const Plan& pl, u32 k0, const std::function<void(const NumQ&)>& emit);

// Run the whole length on the CPU.  `found` receives raw hits (unverified).
RunStats run_cpu(const Plan& pl, const TableView& tv, int threads, std::vector<NumQ>& found);

// Exhaustive reference: every base-P palindrome of length L, checked in base Q.
RunStats run_brute(const Plan& pl, int threads, std::vector<NumQ>& found);
