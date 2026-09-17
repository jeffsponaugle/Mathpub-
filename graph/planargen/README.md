# planargen — counting connected planar graphs (OEIS A003094) with a pruned geng

Counts of unlabeled connected planar simple graphs on n nodes,
[A003094](https://oeis.org/A003094), and as by-products the tables by nodes
and edges ([A049334](https://oeis.org/A049334) connected,
[A039735](https://oeis.org/A039735) all) and the counts of all planar graphs
([A005470](https://oeis.org/A005470)) via the Euler transform.  The tool is
nauty's `geng` with compile-time pruning hooks, run on complements, with a
custom planarity test.

New terms computed with it (both later verified in several independent ways,
see "Verification"):

| n  | A003094(n) connected planar | A005470(n) all planar | run |
|----|-----------------------------|-----------------------|-----|
| 13 | 5942258308 (known, Grasegger 2023; reproduced) | 6295835195 (new, Euler transform) | 2.0 core-hours |
| 14 | **117243822184**            | **123897775473**      | 44.4 core-hours, Apple M1 Pro, 2026-09-09/10 |
| 15 | **2385519362107**           | **2516143483576**     | 1134 core-hours, Apple Mac Studio, 2026-09-10/12 |

Before this, a(13) was the last known term; it had been computed in 2023
with the naive pipeline `geng -c 13 | planarg -q | countg -q`, which filters
all ~5e13 connected graphs on 13 vertices.  That pipeline would need about
350 core-years for n = 14 (see "Cost").  a(14) and the related rows were
submitted to OEIS on 2026-09-11.

## Results

### n = 14 (2026-09-09 20:44 to 2026-09-10 02:46 PDT, Apple M1 Pro)

1000 pieces, 8 processes, 44.4 core-hours (pieces 121–247 s, median 158 s).
Binary `geng_coplanar_lr`, sha256 `020c8ead…` (kept in `runs/n14/`).

Row 14 of A049334, k = 13..36 edges: 3159, 39260, 300748, 1799700, 9123403,
40216577, 153876251, 505912342, 1416544333, 3360456524, 6737246119,
11399818779, 16256246886, 19492037645, 19578032113, 16377361650, 11313238350,
6375615561, 2880692140, 1017576903, 270473444, 50848432, 6022143, 339722.

### n = 15 (2026-09-10 11:53 to 2026-09-12 11:10 PDT, Apple Mac Studio)

4000 pieces, 24 processes, 47.3 h wall clock, 1134.07 core-hours (mean piece
17 min).  Byte-identical binary (same sha256 `020c8ead…`), piece header
`-X0x200000d0D14 n=15 e=66-91`.

Row 15 of A049334, k = 14..39 edges: 7741, 110381, 959374, 6499706,
37304702, 186382186, 810410809, 3040436758, 9769142263, 26778167163,
62540416563, 124447818228, 211016092810, 304772256777, 374432406051,
390255221658, 343573902943, 253865707989, 156007244550, 78725191219,
32043431201, 10253244873, 2481451306, 426720176, 46427839, 2406841.

### Files

- `results_n15/`: `RESULT.txt` (full provenance and checks), the four
  b-files covering n <= 15 (`b003094.txt`, `b005470.txt`, `b049334.txt`
  terms 1–289, `b039735.txt` terms 1–289), `rows_with_n15.json`.
- `results_n14/`: the same for n <= 14 (`oeis/` holds the b-files), plus the
  max-degree-4 cross-check output.
- Rows 13–15 of A039735 (all planar graphs by nodes and edges) are in the
  b-file; the existing OEIS b-file stopped inside row 9.

## Verification

Nothing about these counts follows from theory; the only meaning of
"correct" is agreement between independent computations.  What was checked:

**The method reproduces everything known.**  The production binary
reproduces A003094(n) for n = 1..13 and every entry of A049334 rows 1–13,
including all 34 edge counts of Grasegger's row 13 (128 pieces, 2.0
core-hours).  Both planarity backends (`geng_coplanar` with nauty's tester,
`geng_coplanar_lr` with the left-right test) give identical per-edge output
for all n <= 13.

**Independent values inside the new rows, all exact:**

| row | edges k | value | independent source |
|-----|---------|-------|--------------------|
| 14 | 13..22 | 3159 … 3360456524 | Grasegger's partial row 14 (naive pipeline) |
| 14 | 34, 35 | 50848432, 6022143 | plantri triangulations minus two / one edge (`check_plantri.sh`) |
| 14 | 36 | 339722 | triangulations A000109(14) |
| 15 | 14 | 7741 | trees A000055(15) |
| 15 | 15 | 110381 | connected unicyclic graphs A001429(15) |
| 15 | 16, 17 | 959374, 6499706 | connected graphs with n+1, n+2 edges A001435/A001436 (all planar) |
| 15 | 38 | 46427839 | plantri triangulations minus one edge (93.9M graphs deduplicated) |
| 15 | 39 | 2406841 | triangulations A000109(15) |
| all | column 14 | 444855 | connected planar graphs with 14 edges, A046091(14) |

**Subclasses at n = 14 against nauty's own pipeline**: maximum degree <= 3
(`geng -c -D3 14 | planarg`, 227508 graphs) and <= 4 (`check_D4.sh`,
186092397 graphs) have identical edge distributions to
`geng_coplanar_lr -d10` / `-d9`; the latter also under the production
1000-way partition.

**Backend agreement on production pieces**: 15 of the 1000 pieces of n = 14
and 9 of the 4000 pieces of n = 15 (0, 500, …, 3999) recomputed with nauty's
planarity tester instead of the left-right test; per-edge counts identical in
every case (`check_backend.sh`).

**Planarity test**: `lrtest` compares `lrplanar.c` with nauty's tester: 0
disagreements on all 288248 graphs with 5–9 vertices, all 12005168 graphs
(connected or not) with 10 vertices, and ~14 million random graphs with
11–16 vertices, including disconnected and non-planar disconnected ones.

**Smoothness**: a(n)/a(n-1) = 16.57, 17.96, 18.96, 19.73, 20.35 for
n = 11..15; a(n)·n!/A096332(n) (unlabeled vs labeled) = 1.2048, 1.1990,
1.1996, 1.2034 for n = 12..15.

**Code hygiene**: `make geng_coplanar_lr_asan` (AddressSanitizer + UBSan)
reproduced n = 10, 11 and production piece 0/1000 of n = 14 with no reports.
The comments in `lrplanar.c` were verified to be comment-only (identical
object code).

Not independently confirmed: the middle of each row (k = 23..33 of row 14,
k = 18..37 of row 15), i.e. most of the count, rests on this method alone,
as every single-method OEIS term does.

## One command per term

Single process, from this directory; the answer is the `>Z` line on stderr.
The edge range is C(n,2) minus the planar range n-1..3n-6, i.e. the edge
counts of the *complements* geng generates; times are for one Apple M1 Pro
core.

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
`python3 verify.py <n> coplanar <stderr file>` (rows of A049334 through
n = 14 are in `a049334_rows.json`).

For n >= 14 use the detached parallel launcher: it splits the run into MOD
independent, restartable pieces, runs P at a time at low priority, survives
closing the terminal, and writes `../runs/n<n>/RESULT.txt` (total, per-edge
verification, Euler transform for A005470/A039735, b-file lines) when the
last piece finishes:

```
./launch.sh 14 1000 8       # a(14): 44.4 core-hours; ~6 h on 8 cores
./launch.sh 15 4000 24      # a(15): 1134 core-hours; 47 h on 24 processes
```
Two hosts share a run by piece range (same n and MOD on both), then the part
files are copied into one directory and aggregated:
```
./launch.sh 15 4000 20 0 1999         # host A
./launch.sh 15 4000 20 2000 3999      # host B
python3 aggregate.py 15 4000 ../runs/n15 && python3 euler.py ../runs/n15/rows_with_n15.json 15 && python3 bfile.py ../runs/n15/rows_with_n15.json 15
```
Progress: `ls ../runs/n15/part_*_of_4000.txt | wc -l`.  Stop with
`pkill -f 'geng_coplanar_lr -u -v 15 '`; relaunching skips completed pieces.
MOD only sets the granularity (any value gives the same total): pieces of
2–20 minutes keep the load balanced and limit the loss on interruption; the
fixed cost per piece is ~0.2 s.  `launch_n14.sh`/`finish_n14.sh` are the
n=14-specific versions used for the recorded run.

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
   at n = 10, growing with n.  Because geng works on the complements, its
   `-c` option is not used: connectivity of the planar graph is tested in the
   hook on completed graphs only (it is not hereditary), and disconnected
   intermediate graphs are allowed, as in `geng -c` itself.

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
   (de Fraysseix–Rosenstiehl, Brandes 2009) on bitmask graphs with no
   allocation, about 2× faster again, and heavily commented.  Every
   intermediate and final test is a complete planarity test; the only
   bypasses are exact (edge bound; fewer than 9 edges; the two sufficient
   conditions above).

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

## Cost

Measured with `geng_coplanar_lr` (M1 Pro core-equivalents), against the
naive pipeline's cost estimated from the exact numbers of connected graphs
with <= 3n-6 edges times planarg's measured per-graph cost:

| n  | graphs the naive pipeline filters | naive pipeline | this tool | speedup |
|----|-----------------------------------|----------------|-----------|---------|
| 12 | 4.4e10  | 5 core-days     | 5.4 min          | 1,500× |
| 13 | 5.3e12  | 2 core-years    | 2.0 core-hours   | 8,800× |
| 14 | 8.5e14  | 350 core-years  | 44.4 core-hours  | 70,000× |
| 15 | —       | —               | 1134 core-hours  |        |

Cost grows ~22–26× per vertex (the count grows ~20× and the per-graph cost
slowly), so n = 16 would be roughly 30,000 core-hours (3.4 core-years):
about a month on two 20-core DGX Sparks.

## Build

```
cd ../nauty2_9_3 && ./configure && make      # nauty itself (needs geng's W1 objects)
cd ../planargen && make                      # geng_coplanar, geng_coplanar_lr, geng_direct, lrtest
```
Works unchanged on aarch64 Linux (DGX Spark): gcc, `-march=native`.  For
`check_plantri.sh`, build plantri in `../plantri/` (`gcc -O3 -o plantri
plantri58/plantri.c`).

## Run

Vertex count n, complement edge range C(n,2)-(3n-6) : C(n,2)-(n-1):

```
./geng_coplanar_lr -u -v 12 36:55            # a(12), per-edge counts on stderr
./run_split.sh 13 128 8 out13                # 128 pieces, 8 at a time
python3 aggregate.py 13 128 out13            # sums pieces, verifies, writes rows JSON
python3 euler.py out13/rows_with_n13.json 13 # A005470 and A039735 rows by Euler transform
python3 bfile.py out13/rows_with_n13.json 13 # OEIS b-file lines
```
`-v` edge counts are for the complements: k complement edges = C(n,2)-k
planar edges (`verify.py`/`aggregate.py` translate).  Pieces are independent
and restartable; a piece is complete when its file has a `>Z` line.
`euler.py`/`bfile.py` need every row 1..n in the rows JSON; `aggregate.py`
builds it from `a049334_rows.json` plus the new row.

### Setting up on another machine (Mac or Linux)

```
rsync -a --exclude runs ~/src/math/graph/ othermachine:src/math/graph/
ssh othermachine
cd src/math/graph/nauty2_9_3 && make clean >/dev/null 2>&1; ./configure && make -j8
cd ../planargen && make clean && make
./geng_coplanar_lr -u -v 11 28:45          # must report 17449299 graphs (~15 s on an M1 core)
./launch.sh 15 4000 20                      # P = number of performance cores; detached
```
Two machines that run the same n, edge range and MOD produce identical
pieces, so their part files can be compared piece by piece.  Copying the
binary instead of rebuilding (as done for the n = 15 run) gives a
byte-identical executable on any Apple-silicon Mac.

### Cross-check scripts

- `check_backend.sh`: recompute selected production pieces with nauty's
  planarity tester and diff against the run's files.
- `check_D4.sh` + `compare_naive.py`: max-degree-4 subclass at n = 14 against
  `geng -c -D4 14 | planarg | countg --e`.
- `check_plantri.sh n [2]`: T(n,3n-7) and T(n,3n-8) from plantri
  triangulations by deleting one/two edges and removing isomorphs with
  `shortg` (every planar graph is a spanning subgraph of a triangulation,
  and a triangulation stays connected after deleting two edges).
- `lrtest`: `geng ... | ./lrtest` compares the two planarity testers.

## lrplanar_sg: the planarity test for graphs of any size

At Brendan McKay's request (he wants to include the left-right test in nauty,
without the n <= WORDSIZE limit) `lrplanar_sg.c`/`.h` is a second
implementation of the same algorithm for nauty's `sparsegraph`
representation.  `lrplanar.c` is untouched and remains what the a(14) and
a(15) runs used.  Differences:

- `boolean lrplanar_sg(sparsegraph *sg)`: any size (vertex and edge counts
  must fit in an `int`); loops and parallel edges are ignored, so multigraphs
  are answered for their underlying simple graph.  `lrplanar_dense(g,m,n)`
  wraps it for dense graphs; `lrplanar_freedyn()` frees the work space.
- Work space via nauty's `DYNALLSTAT`/`DYNALLOC1` (thread-local with TLS),
  grown on demand and kept between calls; about 10 ints per vertex and 11
  per edge (phase-1 arrays are reused in phase 2).
- Both depth-first searches are iterative (explicit stacks): a path on 10^7
  vertices is fine.  Outgoing edges are sorted by nesting depth with one
  global counting sort, so high-degree vertices cost linear time.
- Edges are identified without a v x v table: each vertex's neighbour list
  is copied once without loops/duplicates into its slice of a work array,
  and the DFS distinguishes tree, back and already-oriented edges by heights.

Measured on the M1 Pro (input reading included in wall time and memory):

| graph | vertices | edges | answer | wall | peak RSS |
|-------|----------|-------|--------|------|----------|
| random Apollonian network | 10^7 | 3·10^7 | planar | 8.7 s | 1.65 GB |
| same + K5 on 5 random vertices | 10^7 | 3·10^7+10 | non-planar | 4.5 s | 1.46 GB |
| path | 10^7 | 10^7 | planar | 0.8 s | 0.92 GB |
| 2000 x 2000 grid | 4·10^6 | 8·10^6 | planar | 0.4 s | 0.54 GB |
| star K_{1,10^6} | 10^6 | 10^6 | planar | 0.06 s | 83 MB |

On small graphs it runs at the speed of the bitmask version (0.7 us per
graph on the 10-vertex graphs, versus 3.3 us for nauty's tester).  A side
finding: nauty's `planarity.c` is quadratic in vertex degree (11 s for a
star with 10^5 leaves, 46 s for 2·10^5), so the reference tester is skipped
for the high-degree families in the test suite.

### Test suite: `make tests` (or `./run_tests.sh [quick|full]`)

`lrtest_sg` reads graph6/sparse6 of any size (loops and parallel edges
allowed) and compares three testers on every graph: nauty's `planarity.c`
(the original malloc version, linked from the nauty tree), `lrplanar_sg`,
and `lrplanar` (for simple graphs with n <= 32); `-e p|n` additionally
checks against a known answer.  `biggraphs` generates families with known
planarity of any size (grids, cylinders, tori, ladders, Moebius ladders,
paths, cycles, stars, wheels, K_{2,n}, K_{3,n}, random trees and forests,
random Apollonian networks, those with random edges deleted, with K5 or
K_{3,3} planted, subdivided K5, disconnected graphs with a far-away K5, and
planar and non-planar multigraphs with loops and parallel edges).
`run_tests.sh full` (4.7 minutes) runs 86 checks in eight tiers:

1. every graph, connected or not, on 1..10 vertices (12,005,168 at n = 10), three testers;
2. 720,000 random graphs near the planarity threshold, n = 11..64, including sparse ones with many components;
3. 350,000 random regular multigraphs with loops and parallel edges;
4. random sparse graphs with 100..10,000 vertices;
5. 21 known-answer families at 10^4..10^5 vertices;
6. 20 huge known-answer graphs at 10^6..10^7 vertices, with time and memory;
7. regression: `geng_coplanar_lrsg` (the enumeration with the sparse tester
   as backend) reproduces A003094 and the A049334 rows for n = 9, 10, 11;
8. the same harness built with AddressSanitizer/UBSan on a sample of tiers 1-5.

Result on 2026-09-17: 86 passed, 0 failed, no sanitizer reports.

## Notes

- GPUs: the work is dominated by nauty's partition refinement and small
  DFS-based planarity tests on 14–15-vertex graphs, branchy pointer code
  with microsecond granularity.  The Sparks' value here is 2×20 ARM cores,
  not the GB10 GPU.
- Restricting the new vertex's degree further than min-degree+1 is not
  possible without changing geng's canonical deletion rule.
- Possible follow-ups with the same data: A126201 (rooted connected planar
  graphs, needs vertex-orbit counts per graph, ~100 core-hours at n = 14);
  the diagonal T(n,3n-7) and the disconnected counts A005470 - A003094 are
  not in OEIS.
- Files: `coplanar.c` (hooks), `direct.c` (baseline), `lrplanar.[ch]`
  (bitmask test used for the runs), `lrplanar_sg.[ch]` (any-size sparsegraph
  test), `lrtest.c`, `lrtest_sg.c`, `biggraphs.c`, `run_tests.sh`,
  `arena.c`, `planarity_arena.c`/`planarity.h` (patched copies of nauty's),
  `launch.sh`/`finish.sh` (detached runs), `run_split.sh`, `aggregate.py`,
  `verify.py`, `euler.py`, `bfile.py`, `check_*.sh`, `compare_naive.py`,
  `a049334_rows.json` (OEIS reference rows 1..14), `results_n14/`,
  `results_n15/`.
