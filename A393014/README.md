# dbpal — numbers that are palindromic in two bases

`dbpal` searches exhaustively, length by length, for numbers that are palindromes in two bases
at once (optionally only primes). It was built to extend
[OEIS A393014](https://oeis.org/A393014) (*primes that are palindromic in both base 2 and base 9*)
and the sequences around it:

| sequence | what |
|---|---|
| A393014 | primes palindromic in bases 2 and 9 |
| A259385 | numbers palindromic in bases 2 and 9 |
| A046472 / A007632 | primes / numbers palindromic in bases 2 and 10 |
| A046473 / A007633 | … bases 3 and 10 |
| A046477 / A029804 | … bases 8 and 10 |
| A046478 / A029965 | … bases 9 and 10 |
| A046474..A046484, A029961..A029970, A029731, A060792, A259374..A259390, … | the rest of the family (any pair of bases 2..36 works) |

It runs on the CPU (all cores), on Apple GPUs through Metal and on NVIDIA GPUs through CUDA, and
can split a search across several machines. Results and the search log are in [`results/`](results/); see
[RESULTS.md](RESULTS.md) for the full list.

### Results so far

| sequence | before | now (exhaustive) | new |
|---|---|---|---|
| A393014 primes, bases 2 & 9 | a(9) = 2.07·10³³ | **no a(10) below 9⁴⁹ ≈ 5.7·10⁴⁶** | bound |
| A259385 all, bases 2 & 9 | b-file to 2.7·10³¹ (69 terms) | complete below 9⁴⁹ (97 terms) | 28 terms |
| A046478 primes, bases 9 & 10 | 8 terms, to 1.4·10¹² | **a(9) = 314419338131817706717607718131833914413**; a(10) > 10⁴⁴ | 1 term |
| A029965 all, bases 9 & 10 | to 5.5·10¹⁷ (66 terms) | complete below 10⁴⁴ (139 terms) | 73 terms |
| A046479 primes, bases 10 & 11 | 9 terms, a(10) > 10³⁶ | **a(10) = 1926434676726841632361486276764346291, a(11) = 333884947479382241444142283974749488333**; a(12) > 11³⁸ ≈ 3.7·10³⁹ | 2 terms |
| A046473 primes, bases 3 & 10 | a(8) > 10³² | a(8) > 10⁴⁰ | bound |
| A007633 all, bases 3 & 10 | to 3.4·10³³ (99 terms) | complete below 10⁴⁰ (117 terms) | 18 terms |
| A046477 primes, bases 8 & 10 | a(11) > 10³⁴ | a(11) > 10⁴⁴ (being re-run after the base-8 filter fix) | bound |
| A046472 primes, bases 2 & 10 | 6 terms | no new prime among the 183 known A007632 terms (to 10⁵⁵) | — |
| A060792 all, bases 2 & 3 | 17 terms (since 2014), a(18) > 3⁹³ | **a(18) = 93814833782752683486286194707262031368542962943489**; complete below 3¹⁰⁵ | 1 term |
| A046476 primes, bases 7 & 10 | a(9) = 9.7·10²⁸ | a(10) > 10⁴² | bound |
| A029964 all, bases 7 & 10 | to 9.1·10¹⁸ (65 terms) | complete below 10⁴² (140 terms) | 75 terms |

(Longer searches are still running; RESULTS.md is regenerated from the logs by `tools/report.py`.)

## Quick start

```bash
make                      # dbpal (CPU + GPU: Metal on macOS, CUDA on Linux); `make cpu`: CPU only
./dbpal selftest  -b 2,9 --Lmin 5 --Lmax 15      # table method vs brute force, every depth
./dbpal selftest2 -b 2,9 --Lmin 7 --Lmax 15      # class-partitioned method vs brute force
./dbpal search -b 2,9 --below 9^41 --out results/b2_9.tsv --log results/search_log.jsonl
./dbpal plan   -b 2,9 --Lmin 45 --Lmax 51        # parameters and time estimates, no search
python3 tools/oeis.py compare results/b2_9.tsv   # compare with every matching OEIS sequence
python3 tools/verify.py results/b2_9.tsv         # independent check (pure Python, BPSW)
python3 tools/campaign.py run --bases 3,10 --below 10^42   # resumable, length by length
python3 tools/report.py                          # RESULTS.md from results/ and the log
python3 tools/bfile.py results/b2_9.tsv          # OEIS b-files (complete up to the certified bound)
python3 tools/stats.py results/b2_9.tsv          # found vs. heuristically expected, per length
python3 tools/site.py                            # results/site/index.html (shareable page)
```

**Several machines.** `--part i/n` searches only every n-th prefix (v1) or batch (v2) of each
length; a length counts as complete when all n parts are logged. With a copy of the repository on
another Mac:

```bash
python3 tools/campaign.py run --bases 2,9 --below 9^51 --part 0/2          # here
ssh other 'cd dbpal && python3 tools/campaign.py run --bases 2,9 --below 9^51 --part 1/2'
bash tools/push_log.sh user@other         # let the other machine see what is done here
python3 tools/merge_remote.py user@other  # pull its results and log into results/
```

Requirements: a C++17 compiler (clang++ or g++) and GMP (`brew install gmp`, or `libgmp-dev`).
GPU engine: macOS with Metal, or Linux with an NVIDIA GPU and the CUDA toolkit in `/usr/local/cuda`
(override with `make CUDA=...`; on Linux build with `make CXX=g++` if clang is absent). The kernels
are compiled at run time from source embedded in the binary (Metal: `newLibraryWithSource`; CUDA:
NVRTC for the device's own architecture), with the constants of each search baked in, so no offline
GPU toolchain is needed. `DBPAL_CUDA_DEVICE=n` picks the CUDA device.

Output (`--out`, tab separated): `n  isprime  b1  b2  n_in_base_b1  n_in_base_b2`.
The log (`--log`) gets one JSON line per finished length; a length that is logged as `done` has
been searched completely, so the log certifies search bounds ("no further terms below X").

## How it works

Write `P` for the base whose palindromes are enumerated and `Q` for the other one (by default
`Q` is the base that is a power of two, e.g. A393014 uses `P = 9`, `Q = 2`).
A base-`P` palindrome of length `L` is split into `k` outer digit pairs and a middle palindrome
of length `m = L - 2k`:

    N = C + M · P^k,     C = Σ_{i<k} d_i (P^(L-1-i) + P^i).

**Top digits give bottom digits.** Fixing the outer digits (a *node*) confines `N` to an
interval of width `≈ P^(L-k)`; all numbers in it share their top `s ≈ k·log_Q P` base-`Q` digits
`T`. If `N` is also a base-`Q` palindrome, its bottom `s` digits are the reversal of `T`:

    M · P^k ≡ rev(T) − C   (mod Q^s).

**Meet in the middle.** All middles `M` are bucketed in advance by `V = M·P^k mod Q^c`, so each
node costs one bucket probe instead of `P^(m/2)` palindrome tests. Entries whose stored check
digits agree with all `s` known digits get a full test. With the outer and middle halves
balanced this takes about `N^(1/4)` work, against `N^(1/2)` for enumerating the palindromes of
one base. Large `N` would need a table of size `N^(1/4)`, so dbpal has two variants:

* **v1** — one table of all middles (as large as memory allows). Time `≈ N^(1/2) / table size`.
* **v2** — residue-class partitioning (in the spirit of Schroeppel–Shamir). Split the outer
  digits into `A_hi` and `A_lo` and the middle into `B1` and `B2`. The class `c = V mod Q^r` of a
  solution can be computed from both sides: `c = V(B1) + V(B2)` and
  `c = rev_r(top r digits fixed by A_hi) − C_hi − C_lo (mod Q^r)`. Both sides are therefore
  enumerated class by class through small lookup tables. Each batch of classes gets its own
  small, cache-resident bucket table, so time stays `≈ N^(1/4)` and memory drops to `≈ N^(1/8)`.
  For bases 2,9 at `N ≈ 10^46` this is ~40 min instead of ~65 h.

**Bases with a common factor** (2 & 10, 8 & 10, 10 & 15, …): `N mod P^k` is fixed by `C` and
`N mod Q^s` by `rev(T)`, and the two must agree modulo `gcd(P^k, Q^s)`. That check prunes the
tree at every depth (for 2 & 10 half of all children die at every level), which brings the
exponent down to about `N^0.21`.

**Cheap filters.** Even-length palindromes in base `b` are divisible by `b+1`. So (for example)
every number palindromic in bases 2 and 9 has odd base-9 length and an odd centre digit, and
primes need odd lengths in both bases. Lengths, centre digits and leading digits that cannot
work are skipped.

**Arithmetic.** Numbers are stored in base-`Q^c` limbs (`Q^c ≤ 2^63`), so base-`Q` digits are
read straight out of the limbs and no multi-precision division is ever needed. All divisions
are by invariant integers, done with Granlund–Montgomery multipliers. The same `core.h` is
compiled for the CPU and, via run-time compilation with the search parameters baked in as
constants, for Metal. On the GPU, candidates go through a cheap inline pre-filter: the leading
middle digits fix more top digits of `N`, which must mirror its known low digits. The few
survivors are tested fully in a separate pass, which avoids divergence.

## Verification

`make test` runs all of the following:

* `selftest` / `selftest2`: the table methods against an unfiltered brute force for every depth /
  split, on many base pairs, CPU and GPU (thousands of configurations).
* CPU and GPU engines give identical results (`tests/cmp_engines.sh`).
* `tests/oeis_families.py`: 14 sequences of the family are searched over their whole b-file range
  and must match the OEIS exactly (dense ones such as A097856 included). This check is not
  optional: it caught an over-eager parity filter for bases 4, 8 and 16 that the brute-force
  self-tests shared with the engines.
* Every reported number is re-checked with GMP and by `tools/verify.py` (independent pure-Python
  palindrome test + Baillie–PSW); `tools/bfile.py` refuses to write a b-file unless the OEIS
  terms are an exact prefix of the computed list.

## Files

    src/core.h          shared CPU/GPU arithmetic and per-node logic (v1 + v2)
    src/plan*.cpp       parameter choice, cost model, precomputation (GMP)
    src/table.cpp       v1 table (multi-threaded counting sort)
    src/engine*_cpu.cpp CPU engines (v1 batched/prefetching walker, v2 batch loop)
    src/gpu_metal.mm    Metal host code;  src/kernels.metal  Metal kernels
    src/gpu_cuda.cpp    CUDA host code (driver API + NVRTC);  src/kernels.cu  CUDA kernels
    tools/oeis.py       OEIS fetch / family catalog / comparison
    tools/verify.py     independent verifier;  tools/campaign.py  resumable searches
    tools/report.py, bfile.py, stats.py, site.py   RESULTS.md, b-files, statistics, results page
    tools/merge_remote.py, push_log.sh             multi-machine campaigns
    tests/              run_all.sh (make test), cmp_engines.sh, oeis_families.py, limbs.py
    results/            b<B1>_<B2>.tsv results, search_log.jsonl (certified search bounds), bfiles/
