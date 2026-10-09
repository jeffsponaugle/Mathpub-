# A253316 — counting 2n × 2n Takuzu grids

[OEIS A253316](https://oeis.org/A253316): number of 2n × 2n *Takuzu* (Binairo)
grids — 0/1 matrices in which

1. every row and column has n zeros and n ones,
2. no row or column has three consecutive equal entries,
3. all rows are distinct and all columns are distinct.

| n | grid | a(n) | source |
|---|------|------|--------|
| 0 | 0×0 | 1 | OEIS |
| 1 | 2×2 | 2 | OEIS, reproduced |
| 2 | 4×4 | 72 | OEIS, reproduced |
| 3 | 6×6 | 4140 | OEIS, reproduced |
| 4 | 8×8 | 4111116 | OEIS, reproduced |
| 5 | 10×10 | 48183195384 | OEIS, reproduced |
| **6** | **12×12** | **4972955852202492** | **computed with this tool**; matches the value added to OEIS by Kirill Khoruzhii (Sep 2026) by a different method |

`a253316.c` is a multi-threaded C program that computes these terms exactly.
It checkpoints its progress, can be stopped and resumed, and prints progress
with an ETA.

## Build

```
make            # cc -O3 -mcpu=native (or -march=native) a253316.c -lpthread -lm
make test       # self-test (n = 1..5) + symmetry consistency check on n = 6
```

It needs a C11 compiler with `unsigned __int128` (clang or gcc) and pthreads.
It is tested on macOS/arm64 (the inner loop uses NEON there; other platforms use a
portable scalar fallback).

## Running

```
./a253316                 # compute a(6) with all CPUs, checkpoint a253316_n6.ckpt
./a253316 -n 5            # any 1 <= n <= 6
./a253316 -t 10 -p 60     # 10 threads, progress line every 60 s
```

| option | meaning |
|---|---|
| `-n N` | half size (grid is 2N × 2N), 1..6, default 6 |
| `-t T` | worker threads (default: all CPUs) |
| `-c FILE` | checkpoint file (default `a253316_n<N>.ckpt`) |
| `-p SEC` | progress interval (default 10 s) |
| `-m GB` | memory budget for half-grid buffers (default 50 % of RAM) |
| `--no-checkpoint` | neither read nor write a checkpoint |
| `--no-symmetry` | process all seams instead of symmetry classes (8× slower, cross-check) |
| `--selftest` | brute force (n ≤ 4) + known terms (n ≤ 5), symmetric and full modes |
| `--check-orbits K` | recount every seam of K random symmetry orbits, all must agree |
| `--join M` | pair-counting method: `auto` (default), `bitset`, `ie` (inclusion–exclusion), `brute` (one pair at a time; slow reference) |
| `--bench K`, `--unit X Y` | single-thread timing of sample work units |
| `-q` | print only the result |

Progress lines go to stderr, the result to stdout:

```
[10:12:40] units 812/5548  work  41.37%  active 14  elapsed 14m30s  ETA 20m35s (~Mon 10:33)  partial 206433...
a(6) = 4972955852202492
```

`work` is the fraction of the *estimated* total work (a per-unit cost model
calibrated on benchmarks). It includes partial progress of the units in flight.
The ETA is extrapolated from this session's rate. `partial` is the sum over
finished units, so it is a lower bound until the end.

### Checkpoint / restart

The work is split into independent units (one per symmetry class of "seams",
see below). Each finished unit is appended to the checkpoint file as one
line and `fsync`ed:

```
A253316-checkpoint v1 n=6 lines=208 units=5548 sym=1
# U <x> <y> <orbit weight> <count for seam (x,y)> <seconds>
U 38 169 4 3822300903925 128.825
...
# RESULT a(6) = 4972955852202492
```

* Ctrl-C (or SIGTERM) stops the run. Units in flight are discarded, and
  everything already written is kept. A second Ctrl-C exits immediately.
* Run the same command again to resume. Finished units are skipped.
* The file is locked (`flock`), so two instances cannot share it by accident.
* The header must match (same n and settings), and a torn last line from a
  crash is ignored.
* The answer is `sum(weight × count)` over all `U` lines, so it can be
  recomputed from the file alone:
  `awk '/^U/{s+=$4*$5} END{printf "%.0f\n", s}'` (awk uses doubles; fine as a
  sanity check, but not for all 16 digits).

## Method

**Seams.** Split the grid between rows n and n+1. Fix the two middle rows
`(x, y)` (distinct valid lines; 208 × 207 ordered pairs for n = 6). All
remaining constraints separate into a *top half* (rows 1..n ending in x) and a
*bottom half* (rows n+1..2n starting with y):

* Each half is enumerated by DFS. The half's rows must be distinct and valid,
  its columns must have no three equal consecutive entries, the triples crossing
  the seam must be valid, and the top half must not use row y (nor the bottom
  half row x).
* **Column balance** holds exactly when the per-column one-counts of the two
  halves add up to n. Halves are grouped by the count vector (a dense
  mixed-radix index), so only halves in matching groups are combined.
* **Distinct columns:** columns j,k of the whole grid coincide exactly when they
  coincide in both halves. Each half records a bit mask of its equal column
  pairs. A clash can only occur on pairs where both seam rows agree, so the
  masks are reduced to those pairs (after that about 1.3 pairs per half remain
  on average).
* **Distinct rows:** besides the within-half checks, no row of the top half may
  equal a row of the bottom half. Each half carries its n−1 non-seam rows.

There are two independent ways to count the valid pairs in a group of A top
halves and B bottom halves:

* **bitset:** for one side, build a bitset per feature (equal column pair or
  row value) over that side's halves. For each half on the other side, OR the
  bitsets of its ~6 features and popcount; that gives the number of conflicting
  partners. The cost is about A·B·6/64 bit operations (NEON, tiled to stay in
  cache).
* **ie** (inclusion–exclusion over shared rows): the count is
  Σ_R (−1)^|R| · E(T_R, B_R). Here R runs over sets of row values shared by
  both sides, T_R are the halves containing all of R, and E counts pairs with
  disjoint equal-pair masks. E is computed from histograms of the (few hundred)
  distinct masks and a precomputed disjointness table. The cost is linear in
  A+B instead of proportional to A·B.

`auto` (the default) uses `ie` for groups with A·B > 6400·(A+B), i.e. both
sides above roughly 13 000 halves, and `bitset` for the rest. The threshold is
tunable with `--ie-min`.

**Symmetry.** Complementing all entries, mirroring left↔right and mirroring
top↔bottom map seams to seams. These generate a group of order 8. Only one seam
per orbit is computed and weighted by the orbit size: 5548 units instead of
43056 for n = 6. `--check-orbits` verifies that all seams of an orbit give the
same count. `--no-symmetry` recomputes everything without the reduction.

**Scale for n = 6.** Up to 1.0·10⁸ halves per side for the largest seams.
There are 7.40·10¹⁵ candidate (top, bottom) pairs over all 43056 seams, i.e.
pairs whose column counts fit, and 67 % of them survive the distinctness
tests. Memory per thread is bounded
by `-m`: a unit that would not fit is processed in several passes over a hash
partition of the count vectors.

## Verification

* `--selftest`: n = 1..4 against an independent brute-force backtracking
  counter, and n = 1..5 against the OEIS terms. Each is run both with the
  symmetry reduction and over all seams.
* `--check-orbits`: on n = 6, every seam of random orbits (8 seams each) gives
  the identical count.
* `--join brute`: on n = 6 the pair-by-pair reference join gives the same count
  as the bitset join, e.g. seam (29,141): 4,542,667,299 both ways. This covers
  the 66-bit column-pair masks, which only occur for n = 6.
* All three join methods pass `--selftest`.
* The full n = 6 run was repeated with `--join ie`, which shares no pair-counting
  code with the bitset join. All 5548 unit counts are identical, and both
  totals are 4972955852202492:

  ```
  $ ./compare_ckpt.sh a253316_n6.ckpt a253316_n6_ie.ckpt
  a253316_n6.ckpt: 5548 units, total = 4972955852202492
  a253316_n6_ie.ckpt: 5548 units, total = 4972955852202492
  5548 units in both files, 0 mismatches
  ```
* Interrupt/resume was exercised during the production run. A unit computed in
  a benchmark and again in the run gave the same count.

## Performance

Apple M4 Max (10 performance + 4 efficiency cores, 36 GB):

| n | units | wall time (14 threads) |
|---|------|-----------|
| 5 | 931 | 0.3 s |
| 6 | 5548 | 31 min with `--join bitset` (production run: 4.5 min + 26.5 min across one interrupt/resume) |
| 6 | 5548 | 38.7 min with `--join ie` (verification run) |

Peak resident memory was about 12 GB. Single-thread time for the pair counting
of individual n = 6 seams:

| seam (rank by size) | bitset | ie | auto (default) |
|---|---|---|---|
| (54,153) largest | 77.1 s | 41.3 s | 35.7 s |
| (54,143) #73 | 13.1 s | 14.8 s | 11.0 s |
| (42,115) #1317 | 1.8 s | 4.3 s | 1.9 s |
| (29,141) smallest | 0.25 s | 0.90 s | ≈0.3 s |

The default `auto` mode was not timed on a full run. From these numbers it
should beat the 31-minute bitset run, since the largest seams dominate the
running time.

## Files

* `a253316.c`, `Makefile`: the program.
* `compare_ckpt.sh`: compares two checkpoints unit by unit and recomputes both
  totals exactly with `bc`.
* `a253316_n6.ckpt`, `run_n6.log`: the production run (bitset join).
* `a253316_n6_ie.ckpt`, `run_n6_ie.log`: the verification run (ie join).
* `tools/variants_bruteforce.c`, `tools/symmetry_classes.c`: brute-force counts
  for n ≤ 4 of related sequences (no distinctness rule; up to symmetry).
* `submissions.md`: what can be submitted to OEIS, cross-checked against the
  current entries.

Running `./a253316` in this directory finds the completed checkpoint and prints
a(6) immediately. To recompute from scratch, pass a new checkpoint name, e.g.
`./a253316 -c fresh.ckpt`.

n = 7 (14 × 14) is far beyond this machine. There are 33,880 work units with
2.3·10¹⁰ half-grids each (median), up to about 5 TB per unit, and roughly 9–12
CPU-years in total. It is plausible only on a multi-terabyte server over weeks,
after porting the program to n = 7.
