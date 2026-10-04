# A005234 — primes p such that p#+1 is prime — feasibility assessment (2026-10-02)

## Current state of the sequence
- OEIS data already runs to a(27) = 9562633 (9562633#+1, 4,151,498 digits, PrimeGrid, Aug 2025).
  a(23)..a(27) were all added 2024-2025 by Jeppe Stig Nielsen from PrimeGrid's BOINC
  Primorial Prime Search (PRS). b-file is complete to n=27; nothing stale to fix.
- PrimeGrid PRS leading edge on 2026-10-02: n = 10,233,653 (server_status_subprojects.php).
  PRS is a contiguous, double-checked search of both p#+1 and p#-1, using PRST (gwnum, CPU only).
  Frontier moved 3.6M (Jul 2024) -> 10.23M (Oct 2025).

## What a(28) would cost
- a(28) must be the *next* p past the frontier, so a valid new term needs a contiguous search
  from 10.23M upward (a lucky hit far out would not be accepted as a(28) until the gap closes).
- Heuristic density: P(p#+1 prime) ~ e^gamma * ln p / p, so expected count of p#+1 primes in
  [a,b] ~ e^gamma * ln(b/a). Covering 10.2M -> 20.5M gives ~1.2 expected primes.
- That range holds ~617k primes p; after sieving to ~2^70, ~1/3 survive -> ~200k PRP tests
  of 4.4M to 8.9M digits each. PrimeGrid reports 1.4-12.5 h wall per task (multi-core desktop).
- On the 96-core box: ~8 concurrent tests at ~3 h each -> ~60 tests/day -> ~9 years.
  Cost per test grows ~n^2 log n, so the back half of the range is far worse than the front.
- GPU: no usable path. PRST/LLR2 are CPU-only; gpuowl/prpll are Mersenne-only. A generic-modulus
  CUDA NTT for 15-30 Mbit numbers is a major project and GB10 (273 GB/s) would not beat a good CPU.

## Verdict
Not feasible as an independent effort; PrimeGrid will reach any a(28) years before we could.
Only realistic way to "extend A005234" with Jeff's name on it: attach the 96-core box (and/or
Sparks' Arm cores) to PrimeGrid's PRS subproject. Discoverer credit goes to the user whose
machine returns the first-found prime; expected share = share of PRS throughput.

Possible tiny OEIS contribution: a comment "a(28) > 10233653 (PrimeGrid PRS frontier, Oct 2026)",
but it is a moving target and of little value.
