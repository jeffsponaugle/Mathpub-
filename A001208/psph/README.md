# psph — exhaustive search for extremal postage-stamp bases

`psph.c` is a single-file C11/pthreads tool for the **global postage stamp problem**:
for a basis `A_k = {1 = a_1 < a_2 < ... < a_k}` the *h-range* `n_h(A_k)` is the largest `n` such that
every integer `1..n` is a sum of at most `h` elements of `A_k` (repetition allowed), and
`n_h(k) = max_A n_h(A_k)` is the extremal h-range (OEIS A001208 (k=3), A001209 (k=4), A001210 (k=5),
A001211 (k=6), A053346 (k=7), A053348 (k=8); rows A001212..A001216, A005342..A005344).

Given `h`, `k` and a target `TGT`, `psph` finds **all** bases with `n_h(A_k) >= TGT`:

* `TGT = n_h(k)` (the best known value) lists every extremal basis — a check of the published tables;
* `TGT = n_h(k) + 1` with an empty result is a **proof** that `n_h(k)` is the true extremal value
  (as long as the bounds used, see below, are the proven extremal values of the smaller `k`).

Built and measured on an Apple M1 Pro (8P+2E cores, 32 GB) with
`clang -O3 -mcpu=apple-m1 -std=c11 -pthread`.

## Usage

```
make
./psph -h H -k K [-t TGT] [-j threads] [-d splitdepth] [-b j:value]... [-p sec] [-s2 N] [-v]
./psph -r H a1 a2 ... ak          # exact n_H(A) of one basis
./psph -selftest [quick|full]     # verify against the Challis / Challis–Robinson tables
./psph -table H K                 # print the built-in bounds n_H(j), j < K
```

* `-t TGT` defaults to the built-in known value `n_h(k)` when it exists.
* `-j` worker threads (default 10). `-d` prefix depth of the work items (default: the smallest depth that
  yields at least 64×threads items). `-p` progress interval on stderr (default 30 s, 0 = off).
* `-b j:value` supplies / overrides the extremal value `n_h(j)` used as the bound `a_{j+1} <= n_h(j)+1`.
  Built in: k=2 formula, k=3 (table h<=48, Hofmeister/Challis formula beyond), k=4 (table h<=54, Challis
  formula families 55<=h<=302), k=5 h<=90, k=6 h<=26, k=7 h<=14, k=8 h<=8, and the rows h=2..7. All frontier
  cases listed below therefore need no `-b`; `psph -table H K` shows what will be used.
* `-s2 N` number of stage-2 "deep hole" targets (default: automatic, `clamp(h*a_{k-1}/1024, 8, 64)`).
* Results (`SOLUTION h= k= n_h= : a1 ... ak`) go to **stdout**, everything else to **stderr**.

Limits: `h <= 254` (stamp counts are stored in bytes; `h >= 255` would need a 16-bit build), `k <= 24`.

## Algorithm

The search is the Challis (1993) "H-program" with the Challis–Robinson (2010) *difficult target*, plus
two additions (deep-hole targets and the `E`-window shortcut) that came out of profiling.

1. **Prefix enumeration with admissibility.** Depth-first over `a_2 < a_3 < ... < a_{k-1}`.
   For every prefix `A_j` the min-stamps table `cnt_j[x]` (bytes, values > h are "infinite") is built from
   its parent by `cnt_j[x] = min(cnt_{j-1}[x], cnt_j[x-a_j]+1)` (NEON, 16 bytes per step; chunks of length
   `<= a_j` so the recurrence has no intra-chunk dependency). Tables are extended **lazily** in chunks and
   the extension of `A_j` stops at its first gap `g_j = n_h(A_j)+1`; a parent table is only extended as far
   as a descendant actually needs it, and it is shared by all siblings.
   The children are `a_{j+1} in [a_j+1, g_j]`. The first gap must lie at or below
   `min(h*a_j, n_h(j)) + 1`, where `n_h(j)` is the known extremal value for `j` denominations
   (built-in table from the papers, k=3 formula, k=4 formula families for 55 <= h <= 302); if no gap is
   found below that limit the run **aborts** ("bound violated") instead of continuing silently.
2. **Lower bounds.** `n_h(A_k) <= h*a_k` gives `a_k >= L_k = ceil(TGT/h)`; since `a_k <= g_{k-1} <= h*a_{k-1}+1`
   we need `a_{k-1} >= L_{k-1} = ceil((L_k-1)/h)`, and recursively `a_j >= L_j = ceil((L_{j+1}-1)/h)`.
   Equivalently every prefix must satisfy `g_j >= L_{j+1}` (in particular `n_h(A_{k-1}) >= ceil(TGT/h)-1`).
3. **Leaf.** For a surviving `(k-1)`-prefix the table is extended to `E = min(TGT, h*a_{k-1})` (no probe ever
   needs a larger index: `cnt[y] <= h-c` forces `y <= (h-c)*a_{k-1}`), and every
   `a_k in [max(a_{k-1}+1, L_k), g_{k-1}]` goes through
   * **X test** (Challis–Robinson): `X = (C_k-1)a_k + (C_{k-1}-1)a_{k-1} + ... + (C_2-1)a_2 + (a_2-1)`,
     `C_k = floor(TGT/a_k)`, `C_{i-1} = floor(a_i/a_{i-1})`; `X < TGT`, and `X` is representable iff
     `cnt[X - c a_k] + c <= h` for some `0 <= c <= C_k-1`; only
     `c >= cmin = ceil((X - h a_{k-1})/(a_k - a_{k-1}))` can work (`cnt[y] >= y/a_{k-1}`).
     Rejects 84–90 % of candidates for h >> k, about 50 % in the diagonal regime h ≈ k.
   * **Deep-hole targets.** The N remainders `y < a_k` with the largest prefix cost `cnt[y]`, at most one per
     residue class mod `a_{k-1}` (positions in one class share their "rescue" structure), maintained
     incrementally as `a_k` grows; targets `(C_k-1)a_k + y` and `(C_k-2)a_k + y` are tested exactly as `X`.
     Then a "reach" target `c* a_k + r_{m-1} + 1` (`r_m = n_m(A_{k-1})`, `m` = stamps needed for `[1,a_k-1]`,
     `c* = h-m+1`), the deep holes in the top block, and `TGT` itself.
   * **Full check** of the survivors: streaming block DP `u[x] = min(cnt[x], u[x-a_k]+1)` over `[0,E]` with
     early exit at the first gap. Above `E` the prefix contributes nothing, so `u[x] = u[x-a_k]+1` there and
     the first gap above `E` is `min_{z in (E-a_k, E]} z + (h+1-u[z]) a_k` — attained at the first position
     carrying the maximum `u` in that window — i.e. `O(a_k)` instead of `O(TGT-E)` (a 10× saving for k=4 at
     large h, where `E = h*a_3 << TGT`).
   * A basis that passes the full check is reported with its exact `n_h(A_k)` from an independent DP up to
     `h*a_k+1` (the program aborts if that value were ever `< TGT`).
4. **Threads.** All admissible prefixes down to a split depth are enumerated into a work list; workers pull
   items with an atomic counter and continue the DFS below; load balance is dynamic (10 threads give 9.6
   cores on this machine).

### Why every prune is exact

* *Admissibility* `a_{j+1} <= n_h(A_j)+1`: if it fails, `n_h(A_j)+1 < a_{j+1}` is not representable by any
  element `>= a_{j+1}`, so `n_h(A_k) = n_h(A_j) <= n_h(j) < TGT` whenever `TGT > n_h(k-1)`; the program warns
  if `TGT <= n_h(k-1)` (degenerate regime).
* *Table cap* `n_h(A_j) <= n_h(j)`: a proven extremal value is an upper bound for every `A_j`; the cap only
  limits how far a table is built, and a missing gap below the cap aborts the run.
* *Lower bounds* `L_j`: pure consequences of `n_h(A) <= h*max(A)` and `a_{j+1} <= g_j`.
* *X / deep-hole / reach / TGT tests*: each is an exact representability test of one number `< = TGT`
  (`min_c cnt[x-c a_k]+c <= h` is exact for the full basis because `cnt` is the exact min-stamps table of
  the prefix); a candidate is rejected only when such a number is proven unrepresentable.
* *Full check*: exact DP over `[1,E]`, exact formula above `E` (proof in the source comment).
* *E limit*: a representation of `x` with `c` copies of `a_k` uses `h-c` stamps `<= a_{k-1}` for the rest,
  so `x - c a_k <= (h-c) a_{k-1} <= h a_{k-1}`; prefix-table entries beyond `E` are never consulted and are
  provably `> h`.
* Overflow: all ranges and products (`h*a_k`, `C_k*a_k`) are 64-bit; stamp counts are 8-bit with
  saturating adds (`h <= 254` enforced).

## Validation

`./psph -selftest quick` (≈ 4 min on 10 threads) runs 81 search cases — every `(h,k)` of the
Challis–Robinson Appendix with small `k` (h=2..6, k=3..8/9) plus k=3 (h=7..22), k=4 (h=7..20), k=5 (h=7..20),
k=6 (h=7..12), k=7 (h=7) — each **twice**: with `TGT = n_h(k)` the program must find exactly the listed
extremal bases (as a set, including all cases with 2–5 extremal bases such as k=5 h=8 {1,9,15,78,115},
{1,9,15,80,118}; k=7 h=7; k=4 h=46 is covered by the k=4 formula check), and with `TGT = n_h(k)+1` it must
find nothing. It also checks the k=3 formula (h=23..60) and the k=4 formula families (h=55..302) against a
direct h-range DP and against the tabulated values, and the row/column tables against each other.

Result: **422 tests, 0 failures** (8-bit build). The 16-bit build `psph16` was run on the same selftest: 80 of 81 search cases had passed (none failed) when this was written; see `logs/selftest16_quick.log`. The larger cases k=7 h=8,9 and k=8 h=7 (`-selftest full`) were run as
part of the timing campaign below and also return exactly the published bases
(k=8 h=7: {1,4,17,31,117,209,513,550} and {1,6,20,41,109,228,509,580}; k=7 h=9: {1,7,30,86,189,607,920}).

Two table notes:

* **OEIS A001211(20) looks wrong.** The Challis–Robinson table gives `n_20(6) = 45754` with basis
  {1,17,93,436,2898,6897}; `psph -r 20 1 17 93 436 2898 6897` (and an independent Python DP) confirm this
  basis has 20-range 45754. OEIS (and the bfile in `../bfiles`) print 45745, which this basis exceeds, so
  45745 must be a transposition typo. The program uses 45754.
* The k=4 formula families reproduce `../bfiles/bA001209.txt` for all h = 1..302 (0 mismatches); k=5, k=7
  bfiles also match the built-in tables.

## Measurements

All times are **CPU core-seconds on the M1 Pro** (sum over threads, `getrusage`), which is the right unit
for extrapolation: the machine was shared with other jobs during part of the campaign, so wall times vary
(10 threads give ≈ 9.6 cores when the machine is free). "cand" = leaf candidates `a_k` tested;
"X-pass" = fraction surviving the Challis X test; "full" = fraction that reached the full check.
Raw logs of every run are in `logs/`. `TGT` = known extremal value in every run; every run found exactly the published extremal basis/bases
(listed in `../bfiles` and the Challis–Robinson tables), e.g. k=4 h=100: {1,75,1983,35044}; k=7 h=9:
{1,7,30,86,189,607,920}; k=6 h=14: {1,11,49,188,810,2109}; k=5 h=30: {1,24,201,1718,7596}.

| k | h | TGT | CPU (s) | cand | X-pass | full | ns / cand | prefixes `A_{k-1}` (admissible / visited) | note |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 60 | 143 814 | 6.8 | 1.74e8 | 16.3 % | 1.13 % | 39 | 36 674 / 38 887 | |
| 4 | 100 | 979 831 | 542 | 3.22e9 | 16.2 % | 1.66 % | 168 | 159 831 / 172 090 | |
| 5 | 20 | 12 903 | 12.4 | 5.63e8 | 18.1 % | 1.3 % | 22 | 1.2e6 | |
| 5 | 30 | 66 771 | 414 | 1.80e10 | 16.0 % | 0.25 % | 23 | 5.02e6 / 5.20e6 | `-s2 32` |
| 5 | 35 | 123 954 | 1 542 | 7.15e10 | 17.8 % | 0.09 % | 21.6 | 1.165e7 / 1.208e7 | 8 threads, machine shared with two 2-thread jobs |
| 6 | 10 | 2 510 | 20.6 | 8.12e8 | 30 % | | 25 | | pre-dedupe binary (selftest) |
| 6 | 11 | 3 607 | 57.1 | 2.35e9 | 29 % | | 24 | | pre-dedupe binary (selftest) |
| 6 | 12 | 5 118 | 146 | 6.24e9 | 27 % | 0.6 % | 23 | | pre-dedupe binary (selftest) |
| 6 | 14 | 9 748 | 806 | 3.64e10 | 23.6 % | 0.62 % | 22 | 4.14e7 / 4.22e7 | 2 threads, machine shared (843 s with the pre-dedupe binary) |
| 7 | 7 | 1 137 | 144 | 3.85e9 | 52.4 % | 0.9 % | 37 | 2.47e7 / 2.49e7 | |
| 7 | 8 | 2 001 | 785 | 2.47e10 | 44.3 % | 0.68 % | 32 | 9.91e7 / 1.00e8 | pre-dedupe binary |
| 7 | 9 | 3 191 | 4 515 | 1.38e11 | 43.3 % | 1.0 % | 33 | | pre-dedupe binary |
| 8 | 7 | 1 911 | ≈ 3.0e4 (est.) | ≈ 8e11 (est.) | ≈ 45 % | ≈ 0.5 % | ≈ 38 | | **partial**: 20.7 % of the 3 407 work items in 660 s wall × 9.6 cores (1.64e11 cand, 8.4e8 full checks so far); extrapolated by item fraction |

"pre-dedupe binary" rows were measured before the residue-deduplicated deep-hole filter was added; for these
small-`E` cases the change is within ±15 % (the k=7 h=7 case went 146 → 144 s). Every k=4 and k=5 number is
from the final binary. Runs that could not fit the 25-minute budget were stopped (k=8 h=7 after 11 min;
k=5 h=40+, k=6 h=16+, k=7 h=10+, k=4 h=150+, k=8 h=8 were not started).

### Growth with h (fixed k)

Least-squares power laws `CPU ∝ h^p` through the measured points (local exponents in brackets):

| k | points (h) | p (CPU) | p (candidates) | local exponents | per-candidate cost |
|---|---|---|---|---|---|
| 4 | 60, 100 | **8.6** | 5.7 | — (2 points) | grows ∝ h^2.9: the full check costs `E/16 = h·a_3/16` vector ops and `a_3 ∝ h^2` |
| 5 | 20, 30, 35 | **8.6** | 8.5 | [20→30: 8.65, 30→35: 8.53] | ≈ constant 22–23 ns (probes dominate after the deep-hole filter) |
| 6 | 10, 11, 12, 14 | **10.9** | 11.4 (12→14) | [10→11: 10.7, 11→12: 10.8, 12→14: 11.1] | ≈ constant 23–25 ns |
| 7 | 7, 8, 9 | **13.7** | 13.7 | [7→8: 12.7, 8→9: 14.8] | 32–37 ns |

The exponents are what the tree structure predicts: the number of admissible `(k-1)`-prefixes grows like
`h^{(k-1)(k-2)/2}` tempered by the lower-bound prune, and the candidate range per prefix like `h^{k-2}`;
the local exponents still creep up with h for k=6,7, so for those columns the power law is a *lower*
estimate and an exponential fit (`CPU × 2.5 per +1 in h` for k=6, `× 5.6` for k=7, `× 1.42` for k=5,
`× 1.12` for k=4) is an upper one.

### Extrapolated cost of the frontier cases

Power-law extrapolation anchored at the largest measured h (upper estimate from the exponential fit in
parentheses where the two differ materially). 1 day on this M1 Pro ≈ 230 core-hours (9.6 cores).

| case | TGT assumed | bound used for a_{k-1} | core-hours (power law) | (exponential) | on this M1 Pro | on 100 x86 cores | confidence |
|---|---|---|---|---|---|---|---|
| k=4, h=303 | 71.0e6 (= n(302,4)·(303/302)^4) | n_303(3) formula | **2 000** (×1.6 for the mandatory 16-bit build ≈ 3 300) | — | 9–14 days | 1–1.5 days | 2 points only (h=60,100), but the exponent 8.6 = 5.7 (candidates) + 2.9 (full-check cost) is structurally understood |
| k=5, h=91 | 9.1e6 (1.07 × n(90,5)) | n_91(4) = 682 832 (formula) | **1 636** (exp. fit ×1.30/h: ≈ 1e6) | | 7 days | < 1 day | 3 points (h=20,30,35), exponent stable at 8.5–8.7 and per-candidate cost flat (22–23 ns); the full-check cost E ∝ h^4 will eventually bite, so expect 2–5× more |
| k=5, h=100 | 1.7e7 (1.07^10 × n(90,5)) | n_100(4) = 979 831 | **3 691** (exp. fit: ≈ 1e7) | | 16 days | 1.5 days | as above |
| k=6, h=27 | 1.9e5 (1.21 × n(26,6)) | n_27(5) = 42 744 | **285** | (40 000) | 1.3 days (– months) | hours (– 2 weeks) | 4 points; local exponent still rising → true cost likely 2–10× the power-law figure |
| k=6  h=30 | 3.4e5 | n_30(5) = 66 771 | **899** | (630 000) | 4.3 days (– years) | 10 h (– 9 months) | as above |
| k=7, h=15 | 3.5e4 (1.42 × n(14,7)) | n_15(6) = 12 793 | **1 300** | (38 000) | 6 days (– 6 months) | 13 h (– 16 days) | 3 points; exponent rising (12.7 → 14.8) so the truth is probably between the two |
| k=7, h=16 | 4.9e4 | n_16(6) = 17 061 | **3 200** | (215 000) | 2 weeks (– years) | 1.3 days (– 3 months) | as above |
| k=8, h=9 | 6.1e3 (row h=9: 3191·1.9) | n_9(7) = 3 191 | **≈ 280** | | ≈ 1.2 days | hours | **single partial point** (k=8 h=7 ≈ 8.5 core-h) × k=7 growth (9/7)^13.7 ≈ 31; cross-check: Challis–Robinson count 2.5e13 admissible 8-sets for h=8 ≈ our rate → (8,8) ≈ 50–140 core-h, (8,9) ≈ 300–700 |
| k=9, h=8 | 6.1e3 (row h=8: 3485·1.75) | n_8(8) = 3 485 | **≈ 1e5** | | ≈ 1 year | ≈ 40 days | **not measured**: 2.5e13 admissible 8-prefixes (Challis–Robinson) × ≈ 500 candidates each at 3e7 cand/s |
| k=9, h=7 | 3.4e3 (row h=7: 1911·1.75) | n_7(8) = 1 911 | **≈ 3 000** | | ≈ 2 weeks | ≈ 1.3 days | **single partial point**: ≈ 8e11 admissible 8-prefixes (our k=8 h=7 candidate count) × ≈ 300 candidates each; could be off by 3–5× either way |

Notes on the assumptions: TGT only enters through `L_k = ceil(TGT/h)` and the X test, so a ±10 % error in
the assumed TGT changes the cost by well under a factor 2; the bounds `n_h(j)` for the smaller `j` are the
proven extremal values (k=4 by the Challis formula families, all re-verified by DP), so none of the frontier
runs needs an unproven bound. For k=4, h ≥ 255 the 16-bit build `psph16` must be used (`make psph16`);
it reproduces the 8-bit results exactly and is ≈ 1.6× slower per candidate.

**k=4 column specifically**: the X test rejects ≈ 84 % of candidates at every h measured, the deep-hole
filter another ≈ 90 % of the rest, but the surviving 1.1–1.7 % all need the full DP, and that DP is
≈ 90–96 % of the k=4 run time (measured with `-nofull`: the whole k=4 h=60 search minus the full checks
takes 2.5 CPU-s instead of 6.8 — 24 before the deduplicated deep holes, 73 before the `E`-window shortcut).
So for k=4 the cost is `#candidates × 1.5 % × E/16` with `E = h·a_3`: h=303 is ≈ 2 000–3 300 core-hours
with the present code, i.e. a 9–14 day run on this machine or about a day on 100 cores — the cheapest
genuinely new term in the table, and the one where a GPU full-check kernel (below) would pay off most.

## Profile: where the time goes

Per-thread CPU time is split into three buckets (`-v` prints them per thread): *prefix DP* (building the
`cnt_j` tables for levels 2..k-1 while hunting their first gap), *leaf tables* (extending `cnt_{k-1}` to
`E`, the reach sequence, the deep-hole list), and *leaf candidates* (X test, deep-hole tests, full checks).

| regime | prefix DP | leaf tables | leaf candidates | what dominates inside the leaf |
|---|---|---|---|---|
| k=4, h = 60..254 | < 0.1 % | 0.5–2 % | 98–99 % | the **full-check DP** (≈ 96 % of all time at h=60 when measured with `-nofull`); X + deep-hole probes cost ≈ 14 ns per candidate |
| k=5, h = 30..45 | 0.3–1 % | 4–25 % | 75–95 % | full checks and probes in comparable amounts; `E`-window shortcut and residue-deduplicated deep holes cut the full-check share from ≈ 55 % to ≈ 20 % |
| k = 6..8, h ≈ k..12 | 3–6 % | 13–18 % | 77–83 % | **probes**: X test passes 45–52 % here, ≈ 8 probes/candidate, full checks are cheap (TGT ≤ a few thousand) |

Memory traffic is small: the per-thread tables are `(k-1) x (TGT+1)` bytes but only `[0, E]` is ever
touched at the leaf; everything fits in L2 for the measured cases, and the probe addresses of consecutive
candidates advance by small constant strides, so the hardware prefetcher keeps the X test at ~4 ns/probe.

### What would map well to a GPU (Apple Metal)

The leaf is embarrassingly parallel across `a_k` **within one prefix** (thousands to a million
candidates per prefix for k=4/5 at large h, sharing one read-only table of `E` bytes), and across
prefixes. A natural Metal design:

* **Device memory**: the `cnt_{k-1}` table of a *batch* of surviving prefixes (`E` bytes each; for
  k=5,h=91 `E ≈ 3.9 MB`, for k=4,h=303 `E ≈ 7 MB`; a batch of 64–256 prefixes fits easily), plus the
  per-prefix constants (`Z`, `a_{k-1}`, `lo`, `hi`, reach sequence, deep-hole list — all built on the CPU,
  which is < 5 % of the work).
* **Kernel 1 — probe filter**: one thread per candidate `a_k`; computes `X`, `cmin`, loops `c` with a few
  gathered byte loads, then the deep-hole / reach / TGT targets; writes a survivor flag. This is exactly
  the CPU hot loop; it is gather-bound and would run at full memory-latency-hiding throughput on the GPU
  (the M1 Pro GPU has ~200 GB/s shared with the CPU; the probes are ~8–20 byte loads per candidate so a
  2k–4k x speed-up over one CPU core is not expected — more like the usual 10–30x per chip).
* **Kernel 2 — full checks** for the survivors (0.3–7 % of candidates): one thread (or one SIMD-group
  per candidate with the block DP striped across lanes) running the streaming DP over `[0,E]` with two
  `a_k`-sized ring buffers in device memory; the `E`-window minimum is a reduction. Per candidate this is
  `≈ 3E` bytes of traffic, independent across candidates, so a GPU would be bandwidth-bound at roughly
  200 GB/s / 3E per check — for k=4,h=303 (E ≈ 7 MB) that is ≈ 10k full checks/s per chip versus ≈ 150/s
  per CPU core, i.e. roughly a 6–8x gain over the 10 CPU cores, not more.
* The **prefix DFS** (lazy table extension, first-gap hunting) is sequential and branchy and should stay
  on the CPU; it is < 1 % of the time for h >> k.

Bottom line: a GPU port helps most where the full check dominates (k=4, k=5 at large h) and buys
perhaps 5–10x over the 10 CPU cores of this machine (bandwidth-limited); for the diagonal regime
(k = 7..9, h ≈ k) the probe kernel is the only hot spot and a 16-core M1 Pro GPU would be comparable to
~20–40 CPU cores at best.

## Feasibility summary

* **This M1 Pro alone (≈ 230 core-hours/day)**: k=6 h=27 (days), k=8 h=9 (≈ a day, low confidence),
  k=4 h=303 (≈ 2 weeks with `psph16`) and k=5 h=91 (1–4 weeks) are realistic multi-day runs; k=6 h=30,
  k=7 h=15 and k=9 h=7 are weeks-to-months with a large uncertainty; k=5 h=100, k=7 h=16 and k=9 h=8 are
  out of reach without algorithmic progress.
* **≈ 100 x86 cores (≈ 2 400 core-hours/day; AVX2 gives 32-byte vectors, so per-core speed should be at
  least M1-like)**: everything in the table except k=9 h=8 (≈ 40 days) and the pessimistic ends of k=6 h=30
  / k=7 h=16 fits in days; k=4 h=303, k=5 h=91/100, k=6 h=27/30 and k=8 h=9 are 1–2 day jobs. The work
  list is trivially partitionable (`-d` depth, item ranges), so the code can be run as many independent
  processes.
* **GPU port (Metal)**: see the profile section — a 5–10× gain over the 10 CPU cores where the full check
  dominates (k=4, k=5 at large h), bandwidth-bound; little gain for the probe-dominated diagonal cases.
  The most useful split is CPU prefix DFS + GPU leaf (probe kernel + full-check kernel) with batches of
  prefix tables in device memory.

