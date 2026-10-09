// plan.h - host-side planning: choose depth k / key size t for one base-P length L and
// precompute every constant the CPU and GPU engines need.
#pragma once
#include <gmpxx.h>

#include <string>
#include <vector>

#include "core.h"

static_assert(sizeof(FastDiv) == 16, "FastDiv layout");
static_assert(sizeof(NumQ) == 8 * NLIMB, "NumQ layout");
static_assert(sizeof(SegInfo) == 2 * sizeof(NumQ) + 32, "SegInfo layout");

struct SearchOpts {
  u32 P = 9, Q = 2;          // P: palindrome base enumerated + tabulated, Q: base checked via top digits
  bool primesOnly = false;   // only report primes (and skip lengths that cannot hold primes)
  u64 tableMax = 0;          // max table entries (8 bytes each + bucket offsets); 0 = from RAM size
  bool gpu = false;          // cost model: GPU (Metal) or CPU
  int forceK = -1;           // override depth
  int forceT = -1;           // override key digits
  bool noTableLimit = false; // (v2) plan without the v1 table-size limit
  u32 partI = 0, partN = 1;  // search only part partI of partN (interleaved prefixes / batches)
  double memBudget = 0;      // (v2) planning memory budget in bytes; 0 = 35% of RAM
  int threads = 0;           // CPU threads (0 = all)
  bool verbose = false;
};

// Everything that depends on (P, Q, L, k).
struct Plan {
  u32 L = 0;              // base-P length
  bool brute = false;     // too small for the table method
  bool skip = false;      // provably no solutions (parity rules)
  std::string skipReason;
  Params prm{};
  std::vector<NumQ> outerW;    // [i*P + d] = d * (P^(L-1-i) + P^i), i < k
  std::vector<NumQ> midW;      // [j*P + e] = e * (P^(k+j) + P^(k+m-1-j)), j < hm
  std::vector<NumQ> cenW;      // [d] = d * P^(k+hm)
  std::vector<NumQ> recHi;     // M*P^k contribution of the leading middle digits (index / recLoN)
  std::vector<NumQ> recLo;     // ... of the trailing middle digits and the centre (index % recLoN)
  std::vector<u64> midV;       // [j*P + e] = e * (P^j + P^(m-1-j)) * P^k mod rho
  std::vector<u64> cenV;       // [d] = d * P^(hm) * P^k mod rho
  std::vector<NumQ> WmaxAt;    // [i] = (P^(L-2i) - 1) * P^i, i = 0..k
  std::vector<u32> revtab;     // reversal of Q^revR-digit chunks
  std::vector<u32> topDigits;  // allowed values of the leading base-P digit
  std::vector<std::vector<u64>> gtab;  // gtab[i][s] = gcd(P^i, Q^s) (gmode pruning)
  u64 tableEntries = 0;
  double nodesEst = 0, candsEst = 0, costEst = 0;
  mpz_class lowN, highN;  // this length covers [P^(L-1), P^L)
  u32 partI = 0, partN = 1;  // see SearchOpts
};

u64 physical_memory_bytes();

// Radix helpers -------------------------------------------------------------
u32 digits_per_limb(u32 Q);
FastDiv make_fastdiv(u64 d);
NumQ to_numq(const mpz_class& x, const Params& p);
mpz_class from_numq(const NumQ& x, const Params& p);
std::string to_base(const mpz_class& x, u32 b);
bool is_pal_base(const mpz_class& x, u32 b);
u32 ndigits_base(const mpz_class& x, u32 b);

mpz_class mpow(u32 b, u32 e);
u32 check_digits(u32 Q);                                 // largest u with Q^u < 2^32
void fill_seginfo(SegInfo& sg, u32 LQ, const Params& p);  // uses p.t for the zero-low level

// Fill the Q-dependent parts of Params (rho, qpow, qdiv, ndlo, revR); used by the planner and tests.
void init_radix(Params& p, u32 Q, std::vector<u32>& revtab);

// Build a plan for base-P length L.  Returns plan.brute / plan.skip for small / impossible lengths.
Plan make_plan(const SearchOpts& o, u32 L);
void print_plan(const Plan& pl, FILE* f);
