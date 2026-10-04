# A001208 and the postage-stamp family: status and extension analysis

Date: 2026-09-30. Hardware used for all measurements: Apple M1 Pro (10 CPU cores, 16-core GPU via Metal), 32 GB.

## 1. The problem and the sequences

n(h,k) = largest N such that every integer 1..N is a sum of at most h stamps drawn from k
denominations (repetition allowed), maximised over the choice of denominations.
OEIS convention: a(n) = N (one less than the "smallest unobtainable value" used by Lunnon).

| OEIS | Meaning | Terms in OEIS | Proven frontier (literature) |
|---|---|---|---|
| A014616 | k=2, h=n | closed form | all |
| A001208 | k=3, h=n | 1000 (b-file) | all (Challis recursion, h>=23) |
| A001209 | k=4, h=n | 54 | h<=302 (Challis & Robinson 2013 addendum formulas) |
| A001210 | k=5, h=n | 67 | h<=90 (2013 addendum) |
| A001211 | k=6, h=n | 25 | h<=26 (2013 addendum); **h=27 proven here (2026-10-03: a(27) = 186942)** |
| A053346 | k=7, h=n | 13 | h<=14 (2013 addendum) |
| A053348 | k=8, h=n | 8 | h<=8 |
| A001212 | h=2, k=n | 24 | k<=24 (Kohonen & Corander 2014, 606 CPU-days) |
| A001213 | h=3, k=n | 25 | **k<=15 only** |
| A001214 | h=4, k=n | 25 | **k<=12 only** |
| A001215 | h=5, k=n | 20 | **k<=10 only** |
| A001216 | h=6, k=n | 10 | **k<=9 only** |
| A005342 | h=7, k=n | 8 | k<=8 |
| A005343 | h=8, k=n | 8 | k<=8 |
| A005344 | h=9, k=n | 7 | k<=7 |
| A075060 | h=10, k=n | 7 | k<=7 |
| A084192 / A084193 | the array by antidiagonals | 71 / 71 | extendable to 109 / 113 from proven values |
| A196416 | Lunnon's array (N+1) | 66 | extendable to 145 |
| A195618, A196069, A196094, A234941 | A001208+1, A001209+1, A001212+1, A001212+2 | | follow the main entries |

Sources: Challis, "Two new techniques for computing extremal h-bases A_k", Comp. J. 36 (1993);
Challis & Robinson, "Some extremal postage stamp bases", J. Integer Seq. 13 (2010) Art. 10.2.3,
with the July 2013 addendum (challis6-new.pdf); Kohonen & Corander, J. Integer Seq. 17 (2014)
Art. 14.3.4 (k=24 for h=2); Kohonen, J. Integer Seq. 17 (2014) Art. 14.6.8 and arXiv:1503.03416
(restricted 2-bases to k=47); Al Zimmermann's "Son of Darts" contest (2009-2010), final report at
azspcs.com; Friedman's Math Magic page 0403 (last updated 2010-06-21).

## 2. Finding: OEIS carries unproven contest values as if exact

Friedman's table was updated the day after the Son of Darts contest ended; its entries for
3 stamps with 16..25 denominations, 4 stamps with 13..25, 5 stamps with 11..20 and 6 stamps with 10
are the contest's best scores minus one. Robert Price copied them into A001213, A001214, A001215
and A001216 in July 2013 ("from Friedman"). Challis & Robinson's exhaustive results stop at
(3,15), (4,12), (5,10), (6,9), and nothing published since proves any further term. The contest
history itself shows the values are heuristic: A001213(17) was *corrected upward* in Nov 2009
when a contestant beat the previous OEIS value.

So the following OEIS terms are lower bounds, not proven values:
A001213(16..25), A001214(13..25), A001215(11..20), A001216(10). They are very likely correct
(a year-long contest with strong entrants), but the entries should say so. Suggested OEIS action:
add a comment to each entry and keep the terms, or move them to a comment as "a(16) >= 450" etc.

The contest's best-known values beyond the OEIS data (useful as lower-bound comments) are in
`bfiles/sonofdarts_best.txt` (3 stamps to k=40, 4 stamps to k=30, 5 stamps to k=20).

## 2b. Finding: OEIS A001211(20) is a typo

OEIS A001211 (6 denominations) has a(20) = 45745. Challis & Robinson's table gives n(20,6) = 45754 with the
extremal basis {1, 17, 93, 436, 2898, 6897}; two independent h-range computations here confirm that this
basis realises every value 1..45754 with at most 20 stamps, so the OEIS term is too small by a digit
transposition (Friedman's page has the same 45745, and "3560" for 36560 in the same row, so the typo was
copied from there). The b-file in `bfiles/bA001211.txt` carries 45754. A084192/A084193/A196416 are not
affected (that entry lies beyond their extended range). Credit for spotting it: the psph selftest built
by the search tool (section 4), which cross-checks rows against columns.

## 3. Free extensions from the literature (verified here, not yet submitted)

Every basis below was re-verified by direct h-range computation (`tools/hrange`, 16-bit counts).

| Sequence | Current OEIS | New terms available | Source | File |
|---|---|---|---|---|
| A001210 (k=5) | 67 terms | a(68)..a(90) | 2013 addendum table | `bfiles/bA001210.txt` |
| A001211 (k=6) | 25 terms | a(26) = 156744; **a(27) = 186942 (new, this work)** | 2013 addendum; exhaustive search Oct 2-3 2026 (`logs/n27k6_final_audit.txt`) | `bfiles/bA001211.txt` |
| A053346 (k=7) | 13 terms | a(14) = 24466 | 2013 addendum | `bfiles/bA053346.txt` |
| A001209 (k=4) | 54 terms | a(55)..a(302) | formula families A/B/C + 2013 coefficient table; all 248 formula bases verified, and the formulas reproduce the 1993 table wherever they overlap | `bfiles/bA001209.txt` |
| A084192 | 71 terms | to 109 | proven array entries only (stops at (h=5,k=11)) | `bfiles/bA084192.txt` |
| A084193 | 71 terms | to 113 | stops at (h=7,k=9) | `bfiles/bA084193.txt` |
| A196416 | 66 terms | to 145 | stops at (h=7,k=9) | `bfiles/bA196416.txt` |

Note on A001209: the extremal bases for h in 55..302 are Challis's computational results expressed
as formulas (types A, B, C after Mossige and Selmer); they are proven by his exhaustive runs, not by
theorem, which is the same status as the other terms.

## 4. Tools built

* `psp2/psp2.c` — exhaustive search for 2-bases (A001212): Challis-style DFS with admissibility,
  element-wise upper bounds, gaps test, first-gap candidate rule for the last element; 256-bit masks;
  pthreads. Selftest reproduces every published extremal basis and count for k <= 17.
* `psp2gpu/psp2gpu.m` — Metal port (CPU enumerates prefixes, one GPU thread per prefix runs the
  DFS below it). Same selftest passes.
* `psph_small/psphs.c` — same DFS generalised to h = 2..8 (the rows A001213..A005343), masks
  S_1..S_h of sums of at most m elements; selftest passes for every published row value.
* `psph/` — general fixed-k, large-h searcher (Challis difficult-target method) for the columns
  A001209/A001210/A001211/A053346/A053348 and the near-diagonal entries. See `psph/README.md`.
* `tools/hrange.c`, `tools/make_bfiles.py` — verification and b-file generation.
* `heur/sa2.c` — simulated annealing for 2-bases (lower bounds / sanity).

## 5. Measurements and feasibility

(see section 6 for the summary table; raw logs in the scratch directory were copied to `logs/`)

### 5.1 h = 2 (A001212, next term k = 25, lower bound 228)

Proof-mode node counts (target = known value + 1) of `psp2`:

| k | known n(2,k) | proof-mode nodes (T = n+1) | growth | CPU time, 10 threads (clean run) | GPU (Metal) time |
|---|---|---|---|---|---|
| 16 | 104 | 1.72e9 | | 2.1 s (840M nodes/s) | 3.6 s both runs |
| 17 | 116 | 9.16e9 | x5.3 | 12.3 s (750M nodes/s) | 13.2 s (694M nodes/s) |
| 18 | 128 | 5.78e10 | x6.3 | 123 s (470M/s, machine shared) | 83 s (694M nodes/s) |
| 19 | 140 | 4.24e11 | x7.3 | 784 s (540M/s, machine shared) | not run |
| 20 | 152 | 3.57e12 | x8.4 | not run | 6708 s (532M/s, Metal, overnight) |

Checks: for every k <= 17 the searcher finds exactly the published extremal bases (1,1,1,5,3,2,1,2,4,1,1,3,1,1,1
bases for k = 3..17) and nothing at n+1. Node counts are about 2%, 1.1%, 0.7% of the admissible-prefix
counts A167809(k) for k = 16, 17, 18: the gaps test is pulling its weight, but only in the last two or
three levels (for k = 18, levels 16-17 are 85% of all nodes).

Extrapolation: the growth factor rises by about 1.0 per k (5.3, 6.3, 7.3, 8.4 measured through k = 20).
If it keeps rising (9.5, 10.6, ... 13.9) the counts are 5.4e16 for k = 24 and 7.5e17 for k = 25; if it
saturates at x10 (the growth rate of the admissible-prefix counts A167809), 3.6e16 and 3.6e17. Either way
A001212(25) is a 4e17-8e17-node computation with this pruning.

Calibration against the literature: Kohonen & Corander needed 606 CPU-days (2.6 GHz Opterons, 2013) for
k = 24, i.e. 5e7 CPU-seconds, with a C++ implementation of Challis's K-program. Kohonen (2015, Table 1)
reports 1.5e12 prefixes visited to enumerate all 19-prefixes of range >= 123 with the same program; our
k = 19 run at the higher target 141 visits 4.2e11, so the two implementations are in the same class.
Our extrapolated 1e16-5e16 candidate evaluations for k = 24 would need 2e8-1e9 evaluations per second per
2013 core, which is too fast, so either Kohonen's count excludes rejected candidates (ours includes them,
about 10x more), or the growth factor saturates sooner than assumed. Reading the 1993 paper (paywalled;
Challis's address is in the JIS paper) before any k = 25 campaign is worthwhile.

Cost of A001212(25) (proof that n(2,25) = 228, assuming the heuristic does not find 229):

| resource | rate | time for 3e17 nodes | time for 6e17 nodes |
|---|---|---|---|
| this M1 Pro, CPU only | 8e8/s | 12 years | 24 years |
| this M1 Pro, CPU + Metal GPU | 1.5e9/s | 6 years | 13 years |
| Jeff's fleet CPUs (M1 Pro + Mac Studio + 96-core x86 + 2 Spark CPUs ~ 160 cores) | ~1.3e10/s | 9 months | 18 months |
| fleet CPUs + 2 Spark GPUs with a CUDA port (est. 5-10x the M1 Pro GPU each) | ~2.5e10/s | 5 months | 9 months |
| the same, if the growth factor saturates and the count is nearer 1e17 | | ~2 months | |

The Metal port runs at the same speed as the 10 CPU cores (~690M nodes/s), not 10x faster: the DFS is
branchy, keeps a per-thread mask stack in private memory, and diverges within SIMD groups. A level-synchronous
formulation or a CUDA port on the Sparks (far more registers, better divergence handling) should do
better; a 2-3x gain on the CPU side is also available from vectorising the 256-bit shift-or and
closing the last two levels in candidate form.

Known structure that could cut the search: every extremal basis for k <= 24 except k = 10 is restricted
(range = 2 * largest element), and the restricted case is already solved for k <= 47 by Kohonen's
mirror-image meet-in-the-middle; the expensive part is excluding NON-restricted bases, for which no
analogue of the mirror theorem exists (Kohonen posed this as an open problem). We did not find an
elementary prune that beats the gaps test (element-wise lower bounds derived from counting turn out to
be implied by it).


### 5.2 Rows h = 3..8 (A001213..A005343)

Proof-mode node counts of `psphs` (target = known value + 1, masks up to 1024 bits), M1 Pro, 10 threads
(timings were taken with other jobs running, so rates are lower bounds):

| h | k | known n | nodes | time | nodes/s | growth per k |
|---|---|---|---|---|---|---|
| 3 | 10 | 154 | 6.87e8 | 7.6 s | 9.1e7 | |
| 3 | 11 | 186 | 1.82e10 | 184 s | 9.9e7 | x26.5 |
| 4 | 8 | 228 | 1.39e8 | 2.2 s | 6.4e7 | |
| 4 | 9 | 310 | 7.09e9 | 113 s | 6.3e7 | x51 |
| 5 | 7 | 336 | 5.52e7 | 1.6 s | 3.4e7 | |
| 5 | 8 | 524 | 4.78e9 | 142 s | 3.4e7 | x87 |
| 6 | 6 | 388 | 5.31e6 | 0.1 s | 4.5e7 | |
| 6 | 7 | 638 | 5.99e8 | 21 s | 2.8e7 | x113 |

Extrapolating with constant growth factors (optimistic: for h=2 the factor itself grows with k):

| target | what it would settle | est. nodes | M1 Pro CPU time | verdict |
|---|---|---|---|---|
| (3,16): A001213(16) = 450? | first unproven h=3 term | ~2e17 | ~70 years | not feasible with this method; Challis's own (3,15) run is consistent with ~1e16 nodes over several years |
| (4,13): A001214(13) = 878? | first unproven h=4 term | ~5e16 | ~25 years | not feasible |
| (5,11): A001215(11) = 1393? | first unproven h=5 term | ~3e15 | ~3 years | fleet + GPU: weeks; marginal |
| (6,10): A001216(10) = 2287? | the only unproven h=6 term | ~9e14 | ~1 year | fleet + GPU: days to weeks; the best "prove a contest value" candidate |

The admissible-prefix counts for h >= 3 grow much faster than for h = 2 because every new element
creates C(i+h-1, h-1) sums instead of i+1, so the gaps test only bites in the last two levels.
Beating these numbers needs a genuinely better pruning idea for h >= 3 (none in the literature that
we found; Challis's difficult-target technique is designed for h >> k).


### 5.3 Columns k = 4..8 and near-diagonal entries (psph)

`psph` (built and measured by a subagent; full details, selftest results and profile in `psph/README.md`)
implements Challis's H-program approach: DFS over admissible prefixes with lazily extended min-stamp tables,
the Challis-Robinson "difficult target" X test at the leaf, extra deep-hole targets, and an exact full check
for survivors. It reproduces every published extremal basis in 422 selftests (and found the A001211(20) typo).
CPU core-seconds on the M1 Pro with TGT = known value:

| k | h | CPU core-s | leaf candidates | X-test pass | growth fit |
|---|---|---|---|---|---|
| 4 | 60 / 100 | 6.8 / 542 | 1.7e8 / 3.2e9 | 16 % | CPU ~ h^8.6 (full checks dominate: 90-96 % of time) |
| 5 | 20 / 30 | 12 / 414 | 5.6e8 / 1.8e10 | 16-18 % | CPU ~ h^8.6 (x1.42 per +1 in h) |
| 6 | 10 / 11 / 12 / 14 | 21 / 57 / 146 / (run) | 8e8 .. 3.6e10 | 24-30 % | CPU ~ h^11 and still steepening (x2.5 per +1) |
| 7 | 7 / 8 / 9 | 144 / 785 / 4515 | 3.9e9 .. 1.4e11 | 43-52 % | CPU ~ h^13.7 (x5.6 per +1) |
| 8 | 7 | ~3.0e4 (20 % done, extrapolated) | ~8e11 | ~45 % | one partial point |

Extrapolated cost of the frontier terms (core-hours; this M1 Pro delivers ~230 core-hours/day; power-law
estimate first, exponential-fit upper estimate in parentheses; see psph/README.md for assumptions):

| target | new OEIS term(s) | est. core-hours | this M1 Pro | ~100 x86 cores | confidence |
|---|---|---|---|---|---|
| k=4, h=303 | A001209(303) | 2000-3300 | 9-14 days | ~1 day | 2 points, structurally understood; needs the 16-bit build |
| k=5, h=91 | A001210(91) | 1700 (x2-5 likely) | 1-4 weeks | 1-3 days | 2-3 points |
| k=5, h=100 | A001210(100) | 3800 (x2-5) | 3-8 weeks | ~1 week | as above |
| k=6, h=27 | A001211(27) | **DONE 2026-10-03: n(27,6) = 186942**, basis {1,19,194,1095,7370,27669}; actual 3326 core-h (1.1e14 leaf candidates, 7.4e10 full checks), ~41 h wall on mathd+mathb+mathg (see DEPLOY.md) | | | |
| k=6, h=30 | A001211(30) | 1000 (630 000) | 4.5 days (years) | 10 h (9 months) | as above |
| k=7, h=15 | A053346(15) | 1300 (38 000) | 6 days (6 months) | 13 h (16 days) | 3 points |
| k=8, h=9 | A053348(9) = A005344(8) | **DONE 2026-10-02: n(9,8) = 6082**, basis {1,8,27,88,197,521,1226,1461}; ~1.5e14 leaf candidates, ~1 day on 4 x86 boxes + M1 Pro + Orin (see DEPLOY.md) | | | |
| k=9, h=7 | A005342(9), unlocks A084193/A196416 antidiagonal 16 | ~3000 | ~2 weeks | ~1.3 days | single partial point, x3-5 either way |
| k=9, h=8 | A005343(9) | ~1e5 | ~1 year | ~40 days | not measured; from Challis-Robinson's count of 2.5e13 admissible 8-sets |

For a genuinely new term TGT is unknown: run a heuristic (or the column's formula pattern) to get a good
basis first, then prove with TGT = that + 1 and iterate upward if something better turns up; a +-10 %
error in TGT changes the cost by less than 2x. The GPU helps the k=4/5 cases most (full-check DP is
bandwidth-bound; 5-10x over the 10 CPU cores expected), and little in the diagonal regime.


## 6. Summary of extension potential

Ranked by value per unit of effort:

**Tier 0 — ready now, no computation (all data verified here):**
1. Correct **A001211(20)** from 45745 to 45754 (typo; basis exhibited and verified).
2. Add **A001210** a(68)-a(90), **A001211** a(26), **A053346** a(14) (Challis & Robinson 2013 addendum, never entered).
3. Extend **A001209** to 302 terms (formula families, each basis verified) and the arrays
   **A084192/A084193/A196416** to 109/113/145 terms from proven values.
4. Add "lower bound, not proven" comments to **A001213(16+), A001214(13+), A001215(11+), A001216(10)**.

**Tier 1 — new terms feasible on this M1 Pro in days to a few weeks (psph, CPU):**
5. **A053348(9) = A005344(8)** via (k=8, h=9): **DONE — 6082** (2026-10-02; b-files and submission draft in bfiles/).
6. **A001209(303)** via k=4: ~2 weeks here, ~1 day on 100 cores; cheapest genuinely new column term.
7. **A001211(27)** via (k=6, h=27): **DONE — 186942** (2026-10-03; 3326 core-hours on three Xeon boxes; b-file and submission draft in bfiles/).
8. **A001210(91)** via (k=5, h=91): 1-4 weeks here, days on 100 cores; then h=92, 93, ... at ~x1.4 each.
9. **A005342(9)** via (k=9, h=7): ~2 weeks here; also the first missing entry for A084193/A196416.

**Tier 2 — fleet-scale campaigns (bring the Mac Studio, the 96-core box and the Sparks online):**
10. **A001212(25)** (expected answer 228): 3e17-6e17 node visits with the present pruning = months on the
    whole fleet even with a CUDA port; the Metal GPU here is only as fast as the 10 CPU cores. Worth doing
    only after recovering Challis's 1993 pruning details or finding a non-restricted analogue of Kohonen's
    mirror bound. A001212(26) (expected 244) is ~10x more again.
11. **A005343(9)** via (k=9, h=8): ~1e5 core-hours (~40 days on 100 cores).
12. **A001216(10)** proof of the contest value 2287: ~9e14 nodes with the small-h searcher, ~1 year here,
    days-weeks on the fleet; the only contest row value whose proof is within reach.
13. **A053346(15)** via (k=7, h=15): days to months (wide uncertainty).

**Not feasible with the methods at hand:** proving A001213(16) (~2e17 nodes), A001214(13) (~5e16),
A001215(11) (~3e15, marginal) by the direct DFS; these rows need a new pruning idea for h >= 3.

Recommended next step if you want to start computing: run `psph` for (k=8, h=9) and k=4 h=303 on this
machine now (both fit in a week), and submit the Tier-0 items to OEIS in parallel.


## 7. Suggested OEIS actions (drafts; nothing submitted yet)

**Consolidated, paste-ready package: `submission.md` (2026-10-03; includes an OEIS status check of every entry).**

1. **A001210** — extend with a(68)..a(90) from Challis & Robinson's July 2013 addendum
   (b-file `bfiles/bA001210.txt`, all h-ranges independently verified). Edit line:
   "a(68)-a(90) from the July 2013 addendum to Challis and Robinson, added by Jeff Sponaugle".
2. **A001211** — correct a(20) from 45745 to 45754 (see 2b), add a(26) = 156744 (addendum) and **a(27) = 186942 (new; exhaustive search, draft `bfiles/A001211_27_submission.md`)**. **A053346** — add a(14) = 24466 (same source).
3. **A001209** — replace the 54-term b-file by `bfiles/bA001209.txt` (h = 1..302) and add the comment:
   "a(55)-a(302) are given by three formula families (types A, B, C) with the coefficient table in the
   July 2013 addendum to Challis and Robinson; each formula basis was verified to have the stated
   h-range." (Robinson's existing comment mentions only a(29)..a(254).)
4. **A001213, A001214, A001215, A001216** — add the comment: "Terms a(16) onward [resp. a(13), a(11), a(10)]
   are the best values found in Al Zimmermann's 'Son of Darts' programming contest (2010) as tabulated by
   Friedman; they have not been proved optimal, so they are lower bounds. Exhaustive results are known only
   for k <= 15 [12, 10, 9] (Challis and Robinson 2010, addendum 2013)." Editors may prefer to truncate the
   data and move the contest values to comments as "a(16) >= 450, ...". Contest solutions achieving each
   value are listed in `bfiles/sonofdarts_best.txt` (all 111 re-verified).
5. **A084192, A084193, A196416** — extend the b-files from proven values only
   (`bfiles/bA084192.txt` etc.; 109, 113 and 145 terms).
6. **A001212** — comment: "a(25) >= 228 (the extremal restricted basis, Kohonen 2014); simulated annealing
   (this work) finds no 25-element basis with range 229." (Only if we decide it is worth saying.)

## 8. Heuristic probes

`heur/sa2` (h = 2 simulated annealing, 8 threads) recovers the known optimum for k = 10, 14 and 24 within
minutes, which makes its k = 25 result informative: in 400 s and again in 900 s it reaches 228 (the
restricted extremal basis) and nothing higher; for k = 26 (600 s) it reaches 244, again exactly the
restricted optimum A006638(26). This supports the conjecture that the extremal unrestricted bases for
k = 25, 26 are the restricted ones (so A001212(25) = 228 and A001212(26) = 244 are the expected answers,
and the exhaustive runs would be confirmations). `heur/sah` (general h) is much weaker (337 of 385 for (3,15); 1339 of 1545 for
(6,9)), so it says nothing about the contest rows; the contest entrants clearly used far better search.

