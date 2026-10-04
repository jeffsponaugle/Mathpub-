# A007540 — Wilson primes: feasibility study (2026-10-02)

Wilson primes are primes p with (p-1)! ≡ -1 (mod p²). Known terms: 5, 13, 563.
The search bound is p < 2·10¹³ (Costa, Gerbicz, Harvey, *A search for Wilson primes*,
Math. Comp. 83 (2014); arXiv:1209.3436). Nobody has extended it since; OEIS, Wikipedia and
Zimmermann's records page all still quote 2·10¹³.

**Verdict:** with the hardware on hand this bound cannot be moved by a meaningful amount.
Doubling it to 4·10¹³ is 2–4 *years* of work on both DGX Sparks even with a GPU
multiplier (4 years CPU-only), and about a year on the 96-core/2 TB box. A token 5 %
extension is 1–3 months on the Sparks, a +10 % extension about a month on the big box. Details and measurements below. Nothing has been run beyond benchmarks and
nothing is being submitted.

## Why it is so expensive

The only practical algorithm is the one in the paper: for a block of primes
M < p ≤ N compute M! modulo ∏p² (Stage 1), then descend a remainder tree to get
(p-1)! mod p² for every p in the block (Stage 2), with a cyclotomic trick
(Stage 3) that lets each prime be done using ((p-1)/e)! for the largest allowed
e | p-1 (e ∈ {2,4,…,42,…,90}). Stage 2 is cheap per prime (≈ log⁴ p). Stage 1 costs
≈ (M/e)·log² per block *regardless of block size*, so its cost per prime is
∝ p / (block size). The block size is limited by RAM (≈ 270 bytes/prime measured),
and at p ≈ 2·10¹³ the search is deep in the regime where Stage 1 dominates. The
2014 run already consumed 1.1 million core-hours on 2008–2010 clusters
(464k Stage 1, 655k Stage 2, 44k Stage 3), and every doubling of the bound costs
a bit more than the whole previous search.

## What was done

1. **Recovered the original search code.** Robert Gerbicz published it:
   <https://github.com/gerbicz/Wilson-prime> (`pw13.c`, LGPL, 5.9k lines, OpenMP, uses
   GMP internals) plus David Harvey's `ntt` 0.1.2 multiplication library and two
   independent verification programs (`wilsontest.c`, p^{1/2} algorithm). It is vendored
   under `third_party/Wilson-prime/`.
2. **Made it build on arm64 (Apple Silicon and the Sparks).** Two issues:
   * `ntt-internal.h` used `#ifdef AVOID_128_BIT` in two places where the rest of the
     header uses `#if` (the macro is defined to 0), so the x86-only `mulq/divq` asm macros
     were referenced → build failure. Patched (`ntt-avoid128.patch`).
   * The shipped `tunetab.c` has `tune_tab[1] = SIZE_MAX`, which **disables the NTT
     multiplier entirely** (everything falls back to single-threaded GMP). The tuner must
     be run per machine. `build.sh` does this automatically.
   `pw13.c` includes `gmp-impl.h`/`longlong.h`, so `build.sh` builds GMP 6.3.0 from source
   and links statically.
3. **Validated the output.**
   * Primes in [10⁶, 1.2·10⁶] (all 14,440 of them): 42 sampled residues match a Python
     brute-force computation of (p-1)! mod p².
   * Primes near 4·10¹¹ (class e = 6): all 400 residues compared agree exactly with
     Gerbicz's independent `wilsontest` program (different algorithm). Data in
     `runs/2026-10-02/`. `tools/wilson_brute.c` (exact 128-bit arithmetic, pthreads,
     p < 2³⁹; 1e10 in 30 s on 8 M1 cores) reproduces the header test values of `pw13.c`
     and gives exactly the same residues as pw13 and wilsontest for 400000000063 and
     400000000171 (`runs/2026-10-02/brute_4e11.txt`, ~20 min each on 8 M1 cores).
4. **Benchmarked a real workunit on the Sparks** (class e = 6, 2·10⁶ primes starting at
   4·10¹¹, 20 threads):

   | stage | atom1 | atom2 (instrumented) | note |
   |---|---|---|---|
   | Stage 1 ((q-1)/6)! mod ∏p², q ≈ 4·10¹¹ | 3,470 s | 4,041 s | 78 % of the time is 8,863 multiplications of 2²⁰–2²¹ limbs (the Barrett reductions mod the 154 Mbit modulus) |
   | Stage 2 (remainder tree, 2·10⁶ primes) | – | 125 s | ≈ 62 µs/prime |
   | Stage 3 (cyclotomic step) | – | 2 s | |
   | peak RSS | – | 0.54 GB | ≈ 270 B/prime |

   A second workunit on atom1 with 2·10⁷ primes (10× larger block, same start) took
   3,996 s for Stage 1 — the same as the 2·10⁶ block, confirming that Stage 1 cost is
   independent of block size — and 1,468 s for Stage 2 (73 µs/prime, i.e. linear in block
   size with a small log factor), 13 s for Stage 3, peak RSS 2.45 GB (122 B/prime).
   Log: `runs/2026-10-02/pw13_atom1_e6_W2e7.log`.

   On the M1 Pro the 8-thread run deadlocked inside libgomp after a few minutes; the
   Mac is not usable for this code.
5. **Measured the multiplication kernel**, which is where the time goes
   (`tools/nttbench.c`, `runs/2026-10-02/nttbench_atom1.txt`). 1-Gbit × 1-Gbit on a Spark:
   GMP 1 thread 4.8 s; ntt 1 thread 3.4 s; ntt 8–20 threads 1.2 s. The library's parallel
   speedup saturates at ≈ 4× on the GB10's 10 Cortex-X925 + 10 Cortex-A725 cores
   (pinning to the big cores does not help). For comparison the authors reported 6× on
   8 cores in 2012.
6. **Measured what the GB10 GPU could do** (`tools/gpu/nttgpu_bench.cu`,
   `runs/2026-10-02/gpu_nttbench_atom1.txt`): 215 GB/s copy bandwidth, 250·10⁹ 62-bit
   Montgomery mulmods/s, and a validated radix-2 Stockham NTT of length 2²⁶ over a 62-bit
   prime in 116 ms (4.5 ms/pass, bandwidth-bound). A radix-8/16 shared-memory kernel would
   be ~3× faster. Rough conclusion: a CUDA drop-in for `ntt_mpn_mul_bonus` could make the
   big multiplications 3–6× faster than the 20-thread CPU library, i.e. Stage 1 perhaps
   3× faster overall. It does not change the picture by an order of magnitude.

## Cost model

Using the measured numbers (Stage 1 = 4,000 s at M/e = 6.7·10¹⁰ scaling as
(M/e)·log², independent of block size; Stage 2 = 73 µs/prime scaling as log³ p;
memory 122 B/prime at the 2·10⁷ block; class shares from a sieve of
97k primes near 2.5·10¹³ with the paper's e thresholds — e = 2 alone is 11 % of primes
and has M/e = 1.25·10¹³):

| configuration | ms/prime @2.5·10¹³ | +5 % (→2.1·10¹³) | +10 % | +25 % | ×2 (→4·10¹³) |
|---|---|---|---|---|---|
| 1 Spark, CPU, blocks of 2·10⁸ primes (≈25 GB) | 1.35 | 419 d | 858 d | 2,289 d | 12,062 d |
| 1 Spark, CPU, blocks of 8·10⁸ (≈100 GB) | 0.42 | 136 d | 277 d | 729 d | 3,654 d |
| 1 Spark + GPU multiply (Stage 1 ×3, Stage 2 ×1.5), 8·10⁸ blocks | 0.18 | 59 d | 120 d | 313 d | 1,502 d |
| 1 Spark + GPU, optimistic (×6 / ×2), 8·10⁸ blocks | 0.11 | 37 d | 74 d | 191 d | 893 d |
| 96-core/2 TB box, 4 processes × 24 threads, 3·10⁹-prime blocks (assumes x86 ≈ Spark per process) | 0.05 | 17 d | 34 d | 87 d | 400 d |

Divide the Spark columns by 2 for both machines. Stage 1 is 65–90 % of the cost in every
Spark row, so the only real levers are RAM per block (the 2 TB box) and raw big-multiply
speed. The 96-core row is a guess: its CPU model, memory channels and the ntt library's
scaling on it are unknown; the row assumes it equals one Spark per process.

## If you still want to run something

* The cheapest honest contribution is a +5–10 % extension on the 2 TB box (weeks) or the
  two Sparks (months). It would update the OEIS comment "next term > 2·10¹³" and nothing
  else. I would not do it on the Sparks without the GPU multiplier.
* GPU route: write a CUDA `ntt_mpn_mul_bonus` (split into ≤ 93-bit coefficients, 3–6
  62-bit NTT primes, 4-step transforms, CRT), keep the `preserve_op`/overlap semantics by
  copying, and swap it in above ~2¹⁸ limbs. The GPU microbenchmark here is the starting
  point. Expect a few days of work and validation against GMP at every size.
* Related sequences with much lower bounds that the same machinery would move cheaply:
  A157250 Wilson numbers (complete only to 5·10⁸), A197632 Lerch primes (a(6) > 2·10⁸).

## Running the code

```
./build.sh [threads]      # builds GMP 6.3.0, tunes ntt for this machine, builds pw13, wilsontest, tools
mkdir run && cd run
printf '20\n6\n400000000000 400220000000\n2000000\n1\n0\n36000\n' | ../pw13
#        ^cores ^e  ^start ^end            ^primes/block ^print all residues ^no savefile ^wall-time limit (s)
```
Outputs in the run directory: `wilsonres.txt` (all residues if requested), `wilson.txt`
(near misses, |k| ≤ max(100, p/5000)), `wilsondone.txt`. The program checkpoints to
`parallelwilsonwork.txt` and resumes interactively. Keep `interval` ≥ 10⁶: it is also the
Stage 1 sieve chunk size and tiny values make Stage 1 crawl. For e > 2 the program skips
p below a built-in per-e threshold (`treshold_for_e[]` in `pw13.c`; e = 6 starts at
3.8·10¹¹), because smaller p are handled by a smaller e.

Independent check of a range: `printf '400000000000 400000030000\n2000\n1\n' | ../wilsontest`
(single-threaded, ~1 h for 30,000 numbers at 4·10¹¹).
