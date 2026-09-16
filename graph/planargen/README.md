# planargen — counting connected planar graphs (OEIS A003094) with a pruned geng

Target: a(14) of [A003094](https://oeis.org/A003094) (unlabeled connected
planar simple graphs on n nodes), and as by-products row 14 of
[A049334](https://oeis.org/A049334) (by nodes and edges), and rows 13–14 of
[A005470](https://oeis.org/A005470) / [A039735](https://oeis.org/A039735)
(all planar graphs) via the Euler transform.

Status of the OEIS entries (checked 2026-09-09): a(0..13) of A003094 are
known; a(13) = 5942258308 was added by Georg Grasegger (Jul 2023) with the
naive pipeline `geng -c 13 | planarg -q | countg -q`, i.e. by filtering all
~5.05e13 connected graphs on 13 vertices.  For n = 14 that pipeline would
have to sift ~2.9e16 graphs, which is why a(14) is open.

## Result (run of 2026-09-09/10, Apple M1 Pro, 44.4 core-hours, 1000 pieces)

**a(14) = 117243822184** connected planar graphs on 14 unlabeled vertices.
Row 14 of A049334 (k = 13..36 edges): 3159, 39260, 300748, 1799700, 9123403,
40216577, 153876251, 505912342, 1416544333, 3360456524, 6737246119,
11399818779, 16256246886, 19492037645, 19578032113, 16377361650, 11313238350,
6375615561, 2880692140, 1017576903, 270473444, 50848432, 6022143, 339722.
By Euler transform: A005470(13) = 6295835195, A005470(14) = 123897775473,
and rows 13–14 of A039735 (see `results_n14/RESULT.txt`, which also has the
b-file lines).  Checks that passed on the result itself: the ten known
entries k = 13..22 of row 14 (Grasegger's b-file) and the triangulation
count A000109(14) = 339722 at k = 36; 12 of the 1000 production pieces
recomputed with nauty's planarity tester were identical; a(14)/a(13) = 19.73
continues the ratio sequence 16.57, 17.96, 18.96 smoothly.

## One command per term

Single process, from this directory; the answer is the `>Z` line on stderr.
The edge range is C(n,2) minus the planar range n-1..3n-6, i.e. the edge
counts of the *complements* geng generates (see "Run" below); times are for
one Apple M1 Pro core.

| n | command | a(n) = A003094(n) | time |
|---|---------|-------------------|------|
| 1 | `./geng_coplanar_lr -u 1 0:0` | 1 | instant |
| 2 | `./geng_coplanar_lr -u 2 0:0` | 1 | instant |
| 3 | `./geng_coplanar_lr -u 3 0:1` | 2 | instant |
| 4 | `./geng_coplanar_lr -u 4 0:3` | 6 | instant |
| 5 | `./geng_coplanar_lr -u 5 1:6` | 20 | instant |
| 6 | `./geng_coplanar_lr -u 6 3:10` | 99 | instant |
| 7 | `./geng_coplanar_lr -u 7 6:15` | 646 | instant |
| 8 | `./geng_coplanar_lr -u 8 10:21` | 5,974 | 0.01 s |
| 9 | `./geng_coplanar_lr -u 9 15:28` | 71,885 | 0.05 s |
| 10 | `./geng_coplanar_lr -u 10 21:36` | 1,052,805 | 0.8 s |
| 11 | `./geng_coplanar_lr -u 11 28:45` | 17,449,299 | 15 s |
| 12 | `./geng_coplanar_lr -u 12 36:55` | 313,372,298 | 5.5 min |
| 13 | `./geng_coplanar_lr -u 13 45:66` | 5,942,258,308 | 2.0 h |

Add `-v` for the per-edge counts and check them with
`python3 verify.py <n> coplanar <stderr file>` (rows of A049334 through n=13
are in `a049334_rows.json`).

For n = 14 and 15 use the detached parallel launcher, which splits the run
into MOD independent, restartable pieces, runs P at a time at low priority,
and writes `../runs/n<n>/RESULT.txt` (total, per-edge verification, Euler
transform for A005470/A039735, b-file lines) when the last piece finishes:

```
./launch.sh 14 1000 8                 # a(14): ~44 core-hours -> ~6 h on 8 cores (done 2026-09-10: 117243822184)
./launch.sh 15 4000 24                # a(15): 1134 core-hours (done 2026-09-12 on a Mac Studio: 2385519362107)
```
Two hosts share a run by piece range (same n and MOD on both):
```
./launch.sh 15 4000 20 0 1999         # host A
./launch.sh 15 4000 20 2000 3999      # host B
# then copy host B's runs/n15/part_* into host A's runs/n15 and run
python3 aggregate.py 15 4000 ../runs/n15 && python3 euler.py ../runs/n15/rows_with_n15.json 15 && python3 bfile.py ../runs/n15/rows_with_n15.json 15
```
For n = 15 the only automatic reference is the triangulation count
A000109(15) = 2406841 at 39 edges; `launch_n14.sh`/`finish_n14.sh` are the
n=14-specific versions used for the recorded run.

## Result for n = 15 (run of 2026-09-10/12, Apple Mac Studio, 1134 core-hours)

**a(15) = 2385519362107** connected planar graphs on 15 unlabeled vertices,
computed with the byte-identical binary (sha256 020c8ead...) in 4000 pieces,
24 processes in parallel, 47.3 h wall clock (Sep 10 11:53 to Sep 12 11:10
PDT), mean piece 17 min.  A005470(15) = 2516143483576 by Euler transform;
rows 15 of A049334 and A039735 and the four b-files (n <= 15) are in
`results_n15/`.  Checks on row 15 itself, all exact: k = 14, 15, 16, 17
equal the counts of trees (A000055), connected unicyclic graphs (A001429)
and connected graphs with n+1 and n+2 edges (A001435, A001436), which are
all planar; k = 39 equals the triangulation count A000109(15) = 2406841; the
column sum for k = 14 over all rows equals A046091(14) = 444855; the ratios
a(n)/a(n-1) = 18.96, 19.73, 20.35 and a(n)*n!/A096332(n) = 1.1990, 1.1996,
1.2034 continue smoothly; k = 38 was reproduced with plantri (see below);
and 9 of the 4000 pieces recomputed with nauty's planarity tester as the
backend were identical to the production run.  Row 15 (k = 14..39): 7741, 110381, 959374,
6499706, 37304702, 186382186, 810410809, 3040436758, 9769142263,
26778167163, 62540416563, 124447818228, 211016092810, 304772256777,
374432406051, 390255221658, 343573902943, 253865707989, 156007244550,
78725191219, 32043431201, 10253244873, 2481451306, 426720176, 46427839,
2406841.

## Independent check of the dense end with plantri

`check_plantri.sh n [2]` uses plantri (Brinkmann & McKay's triangulation
generator, a different program and method) plus nauty's `deledgeg`/`shortg`:
every planar graph is a spanning subgraph of a simple triangulation, and a
triangulation stays connected after deleting up to two edges, so the
connected planar graphs with 3n-7 (3n-8) edges are exactly the distinct
graphs obtained by deleting one (two) edges from the triangulations.
Results, all equal to the generator's values:

| n  | T(n,3n-7) one edge deleted | T(n,3n-8) two edges deleted |
|----|----------------------------|-----------------------------|
| 12 | 108597 (known)             |                             |
| 13 | 797583 (known)             | 6135196 (known)             |
| 14 | **6022143**                | **50848432**                |
| 15 | **46427839** (93.9M graphs deduplicated) | not run (1.8e9 graphs) |

## The idea

1. **Prune inside geng instead of filtering its output.**  Planarity is
   hereditary (closed under deleting vertices), and geng's canonical
   augmentation builds every graph by adding one vertex at a time to an
   induced subgraph, so a `PREPRUNE`/`PRUNE` hook that rejects non-planar
   intermediate graphs makes geng generate *only* planar graphs, each exactly
   once.  Search-tree size drops from ~3e16 to ~1e11 nodes at n = 14.

2. **Generate complements.**  geng always adds a vertex of *maximum* degree
   (`geng.c`: `xlb >= dmax`, ties broken by canonical labelling).  For sparse
   classes that is the expensive direction: thousands of neighbourhood
   subsets per parent.  Complements of planar graphs are also a hereditary
   class, and "maximum degree in the complement" is "minimum degree in the
   planar graph"; every planar graph has a vertex of degree <= 5, so the new
   vertex needs at most 6 neighbours and only ~10–500 subsets per parent
   have to be tried.  `coplanar.c` complements each candidate and tests the
   complement.  Measured: 2.8× fewer candidate tests than the direct prune
   at n = 10, growing with n.

3. **Cheap tests first.**  Edge bound (<= 3n-6), connectivity of the final
   graph, and two sufficient conditions for planarity of parent + new vertex
   (degree <= 1; degree 2 with adjacent neighbours) skip ~45% of tests.

4. **Test planarity after canonicity at the final level.**  geng calls
   `PREPRUNE` before and `PRUNE` after its canonicity test; at the last level
   the planarity test is deferred to `PRUNE`, so non-canonical candidates
   never get one.

5. **Faster planarity test.**  nauty's tester (`planarity.c`) does ~100
   malloc/free pairs and 185 `ASSERT`s per call.  `planarity_arena.c` is a
   copy routed through a bump allocator (`arena.c`) and compiled with
   `-DNDEBUG`; `lrplanar.c` is a from-scratch left-right planarity test
   (Brandes 2009) on bitmask graphs with no allocation, about 2× faster
   again.  `lrtest` validates it against nauty's tester: 0 mismatches on all
   graphs with n <= 9, all 11,716,571 connected graphs with n = 10, and
   3,000,000 random graphs with n = 11..16 near the planarity threshold.

Single-core time to count all connected planar graphs on 10 / 11 vertices
(Apple M1 Pro):

| variant                                              | n=10    | n=11   |
|------------------------------------------------------|---------|--------|
| naive `geng -c n \| planarg -u`                      | 103 s   | (~3 h) |
| `geng_direct` (planarity PREPRUNE, geng -c)          | 13.6 s  |        |
| `geng_coplanar`, nauty tester with malloc            | 13.7 s  | 65 s   |
| + arena allocator + trivial-planarity shortcuts      | 3.1 s   |        |
| + `-DNDEBUG`                                         | 2.5 s   | 51 s   |
| + planarity deferred to PRUNE at final level         | 1.41 s  | 27.6 s |
| `geng_coplanar_lr` (left-right test)                 | 0.78 s  | 15.1 s |

n = 12 with `geng_coplanar_lr`: 323 core-seconds.  Cost grows ~22× per
vertex (the count grows ~19× and the per-graph cost slowly), so n = 13 is
~2 core-hours, **n = 14 about 45 core-hours (~2 core-days)** and n = 15
roughly 45 core-days.

## Validation

`verify.py` compares the per-edge-count output (`-v`) with A049334; all of
the following matched exactly (every edge count, every total):

| n  | A003094      | run                                   |
|----|--------------|---------------------------------------|
| 7–11 | 646 … 17449299 | both backends, single process       |
| 12 | 313372298    | 32 and 16 pieces (`run_split.sh`), both backends |
| 13 | 5942258308   | 128 pieces, both backends: all 34 edge counts of Grasegger's row 13 reproduced; `geng_coplanar_lr` 2.01 core-hours (pieces 47–67 s), `geng_coplanar` with planarity before canonicity 9.5 core-hours; per-edge outputs identical |

Backend agreement at n = 14: pieces 0, 1, 2 of 1000 were run with both
`geng_coplanar` (nauty tester) and `geng_coplanar_lr`; totals (100948133,
101643519, 128927287) and all per-edge counts were identical.  The LR pieces
took 144–184 cpu-seconds each, so all 1000 pieces are ~45 core-hours; the
expected a(14) is about 1.1e11.

Independent checks done at n = 14 itself: the subclasses with maximum
degree <= 3 (`geng -c -D3 14 | planarg` vs `geng_coplanar_lr -d10 14 55:78`,
both 227508 graphs) and maximum degree <= 4 (`check_D4.sh`: `geng -c -D4 14 |
planarg | countg --e` vs `geng_coplanar_lr -d9 14 55:78`, both 186092397
graphs), with identical edge distributions in both cases; the max-degree-4
subclass was also rerun split 1000 ways (the production partition
parameters) and summed to the same 186092397.  Still to come
when the full run finishes: the partial row 14 of
A049334 for k <= 22 edges in the OEIS b-file; and the count at 3n-6 = 36
edges, which must equal the number of triangulations A000109(14) = 339722.

Code hygiene: `make geng_coplanar_lr_asan` builds the generator with
AddressSanitizer and UBSan; it reproduced the n = 10 and n = 11 counts and
production piece 0/1000 of n = 14 exactly with no sanitizer reports.

## Build

```
cd ../nauty2_9_3 && ./configure && make      # nauty itself (needs geng's W1 objects)
cd ../planargen && make                      # geng_coplanar, geng_coplanar_lr, geng_direct, lrtest
```
Works unchanged on aarch64 Linux (DGX Spark): gcc, `-march=native`.

## Run

Vertex count n, complement edge range C(n,2)-(3n-6) : C(n,2)-(n-1):

```
./geng_coplanar_lr -u -v 12 36:55            # a(12), per-edge counts on stderr
./run_split.sh 13 128 8 out13                # 128 pieces, 8 at a time
python3 aggregate.py 13 128 out13            # sums pieces, verifies, writes rows JSON
python3 euler.py out13/rows_with_n13.json 13 # A005470 and A039735 rows by Euler transform
```
`-v` edge counts are for the complements: k complement edges = C(n,2)-k
planar edges (`verify.py`/`aggregate.py` translate).  Pieces are independent
and restartable; a piece is complete when its file has a `>Z` line.

### Setting up on another machine (Mac or Linux)

```
rsync -a --exclude runs ~/src/math/graph/ othermachine:src/math/graph/
ssh othermachine
cd src/math/graph/nauty2_9_3 && make clean >/dev/null 2>&1; ./configure && make -j8
cd ../planargen && make clean && make
./geng_coplanar_lr -u -v 11 28:45          # must report 17449299 graphs (~15 s on an M1 core)
./launch.sh 14 1000 12                      # P = number of performance cores; detached
# ... later:  cat ../runs/n14/RESULT.txt    (written automatically when all pieces finish)
```
Two machines that run the same n, edge range and MOD produce identical
pieces, so `grep -h '>Z' runs/n14/part_*` on both can be compared piece by
piece as an independent check.

### n = 14 on two DGX Sparks (20 cores each)

```
./launch.sh 14 1000 20 0 499          # host A
./launch.sh 14 1000 20 500 999        # host B
# then copy host B's runs/n14/part_* into host A's runs/n14 and
python3 aggregate.py 14 1000 ../runs/n14
python3 euler.py ../runs/n14/rows_with_n14.json 14
```
Expected ~45 core-hours total, i.e. a few hours wall clock.  Use the same
n, edge range and MOD everywhere (geng's res/MOD partition depends on them;
`-x`/`-X` must also be identical if you change them).  A single M1 Pro does
it overnight (~6 h on 8 cores).

## What did not help / notes

- GPUs: the work is dominated by nauty's partition refinement and small
  DFS-based planarity tests on 14-vertex graphs, branchy pointer code with
  microsecond granularity.  The Sparks' value here is 2×20 ARM cores, not
  the GB10 GPU.  A batched GPU planarity filter would be possible but the
  filter is no longer the bottleneck after items 4–5.
- Restricting the new vertex's degree further than min-degree+1 is not
  possible without changing geng's canonical deletion rule.
- Files: `coplanar.c` (hooks), `direct.c` (baseline), `lrplanar.[ch]`,
  `lrtest.c`, `arena.c`, `planarity_arena.c`/`planarity.h` (patched copies of
  nauty's), `run_split.sh`, `aggregate.py`, `verify.py`, `euler.py`,
  `a049334_rows.json` (OEIS reference rows 1..13).
