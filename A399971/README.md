# A399971 — compact Numberlink solutions

a(n) = number of compact Numberlink solutions on an n×n grid, up to the 8
symmetries of the square (rules: https://oeis.org/A399971/a399971.txt).

| n | a(n) | not reduced for symmetry | starting positions (reduced) |
|---|------|--------------------------|------------------------------|
| 1 | 0 | 0 | 0 |
| 2 | 1 | 2 | 1 |
| 3 | 13 | 86 | 13 |
| 4 | 2139 | 16418 | 2139 |
| 5 | 2004574 | 16029978 | 2004408 |
| 6 | 11263684837 | 90107368610 | 11262875246 |
| 7 | ≈ 3.51·10^14 (sampling estimate, ±0.5%) | | |

All path covers of the grid (every cell used, every path ≥ 2 cells, before the
compactness test), from `dp -P`: 6, 242, 99636, 232227462, 3952987691824,
448633115022881616 for n = 2..7.

By number of links (see `results/links_by_count.txt`): the minimum is 2, 2, 3,
3, 3, 3 for n = 2..7, so the draft's guess 1 + floor(n/2) fails from n = 6 on
(38 and 44 three-link solutions, all re-checked exhaustively); the maximum is
floor(n²/2) (A007590), not ceiling(n²/2) (A000982) as the draft says.

## Programs

* `a399971.c` — the solver. Canonical (symmetry-reduced) enumeration of path
  covers whose paths are induced, plus an exact compactness test per cover
  (search for a routing of the same endpoint pairs that leaves a cell empty).
  Multithreaded; n = 6 takes about 3½ minutes on a 14-core laptop.
* `verify.c` — independent brute force, no shortcuts: every partition into
  paths, every routing, Burnside for the symmetry reduction. Practical to n = 5.
* `dp.c` — transfer-matrix experiment that tracks every alternative routing
  along the sweep. Exact (agrees for n ≤ 5) but its state grows faster than the
  enumeration, so it is not the way to n = 7. `dp -P n` counts all path covers.
* `run_ranges.sh`, `merge_ranges.py` — split a run into resumable chunks,
  possibly on several machines, and add the pieces up.

```
cc -O3 -march=native -pthread -o a399971 a399971.c -lm
cc -O2 -o verify verify.c
cc -O3 -march=native -o dp dp.c
./a399971 1 6
```

Useful options: `-e SECS` estimates a(n) and the full running time by random
sampling; `-x K` re-decides a sample of covers with an exhaustive router;
`-L K` keeps only solutions with at most K links (fast for small K);
`-P` / `-p` print starting positions / drawings; building with `-DXCHECK` gives
a variant with different canonical representatives and the original witness
search, for cross-checking.

## Verification

* n ≤ 5: `verify.c` reproduces every count above (and the full-solution counts).
* n = 6: computed by four differently configured builds (row order; row order
  reversed with the original search; ring order; ring order reversed with the
  original search), all identical, with ~3.7 million covers re-decided by
  exhaustive routing and no disagreement; a chunked run merged by
  `merge_ranges.py` gives the same totals.
* n = 7: 54176 covers from the middle of the search re-decided by the
  exhaustive router (about half of them non-compact), no disagreement.
* The solver also checks internally that every compactness search finds the
  cover it started from.

## Computing a(7)

About 1.0·10^8 thread-seconds on an M-series laptop (≈ 83 days on 14 threads,
≈ 56 days on a 24-thread M2 Ultra, ≈ 33 days on both); the first chunk on the
M2 Ultra (units 0..495, 165941115102 solutions) took 39 minutes. The work is 992055
independent units:

```
./run_ranges.sh 7 2000 0 999       # machine A: chunks 0..999
./run_ranges.sh 7 2000 1000 1999   # machine B: chunks 1000..1999
./merge_ranges.py results/n7/chunk_*.txt   # after copying all chunk files together
```

Every machine must use the same chunk count (2000 here), since it determines
which units a chunk number means.

Interrupted runs resume where they stopped (finished chunks are skipped).
Chunks computed by different builds can be mixed as long as n, the split (16)
and the cell order (2) are the same; `merge_ranges.py` refuses to add ranges of
different runs and reports gaps.

A run of `./run_ranges.sh 7 2000` on the Mac Studio (2026-10-04/05) was stopped
after 10 chunks (units 0..4959, results in `~/A399971/results/n7` there);
running the same command resumes it. On machines with many threads use fewer,
bigger chunks (e.g. 100 in total) so threads do not sit idle at the end of each
chunk; chunk numbers only mean the same units for the same chunk count.

`results/` holds the n = 6 outputs of the different builds and the merged
chunked run.

## Approaches that did not pay off

* `dp.c`: tracking all alternative routings across the sweep is exact, but the
  sets of routings are large (≈170 on average already at n = 5) and rarely
  shared between covers, so its state count grows much faster than the number
  of covers.
* Pre-routing paths that are forced as soon as their neighbourhood is decided
  (sharing that part of the compactness test between covers): correct, saves
  29% of the search steps at n = 6, but costs more than it saves.
