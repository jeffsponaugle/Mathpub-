# OEIS follow-on work: A171740's cross-references

A171740 ("first number palindromic with n digits in more than one base") has
23 cross-referenced sequences from the same author (James G. Merickel,
2009–2012), in two families. This directory holds the tools, validation and
results for extending and checking them with the machinery built for
A171740 (see `../README.md`).

| Family | Question | Sequences |
|---|---|---|
| Base-gap | smallest L-digit palindrome in two bases differing by d | A216841 (L=2), A216840 (3), A216843 (4), A216899–A216909 (5–15), A216910 (17); **no sequence for L=16** |
| Multi-base | least number that is a k-digit palindrome in at least n bases | A171701 (k=2), A171702 (3), A171703 (4), A171704 (5), A171705 (6), A171706 (7); A171741 / A171742 are the 3-fold / 4-fold cross-sections |

## Layout

```
A171701/   a171701.c — divisor sieve for k = 2; b171701.txt (1152 terms)
A171702/   b171702.txt (223 terms; a(101)-a(223) new, from multipal k = 3)
A171703/   b171703.txt (22 terms; a(10)-a(15) proven exact, a(16)-a(22) new)
A171740/   b171740.txt (n = 1..35, consecutive; replacement for the published file)
A216/      gappal.c  — the A171740 engine restricted to one base gap (-D d)
           run_gap.py, batch.sh — parallel, resumable drivers
           known/    — OEIS data and b-files used for validation
           results/  — L<L>.tsv: "d value b1 b2 seconds", one line per term
multipal/  multipal.c  — CPU sieve (dense counting or radix-sort mode)
           multipalg.m — Metal GPU sieve (bucket scatter + threadgroup sort)
           runs/       — outputs of the proof runs
```

## Tools

### gappal (base-gap family)

`gappal.c` is `a171740.c` with one option added: `-D d` keeps only the partner
base `b1 = b2 − d`, so each stage is a single base pair. Everything else is
the validated A171740 engine: stages ascend in the larger base, the first
match is kept only after every stage whose floor lies below it is cleared,
MITM for L ≥ 8, 192-bit exact values. Even L with d = 1 is rejected up front
(an even-length palindrome in base b is divisible by b+1, so bases b and b+1
cannot share one).

One engine change was needed: near the threshold where a pair first becomes
possible the value window is very narrow, and the A171740 rule "largest MITM
table the memory cap allows" built millions of table entries per stage to
probe a handful of prefixes. `gappal` picks the table depth minimizing
(table entries + ~4 × prefixes). The search is exact at any depth. This took
L=8, d=1000 from minutes to 8 seconds.

    ./gappal -n 8 -D 2 -M -q          # -> RESULT L=8 d=2 v=2727014580 b1=17 b2=19
    ./run_gap.py -L 8 --dmax 1001 -j 22

### multipal / multipalg (multi-base family)

A k-digit palindrome in base b is v = Σ dᵢ·Pᵢ(b) with Pᵢ = b^(k−1−i) + bⁱ.
Because v increases lexicographically in its digits, the palindromes of each
base inside a value chunk [S, E) are found level by level with one division
per bound, and the innermost digit sweeps an arithmetic progression. Every
k-digit palindrome of every base is generated chunk by chunk, and the number
of bases producing each value is counted exactly:

- **CPU (`multipal`)**: a u16 counting array where palindromes are dense
  (k = 2, 3: around 10¹⁰ a number is on average a 3-digit palindrome in ~4
  bases), LSD radix sort where they are sparse (k ≥ 4).
- **GPU (`multipalg`)**: the host lists the chunk's arithmetic runs; kernel
  `gen` scatters values into 65,536 buckets by high bits; kernel `bsort`
  bitonic-sorts each bucket in threadgroup memory and counts runs. Every
  value reported with ≥ 3 representations is re-checked on the host in exact
  arithmetic before it is recorded.

Output: for each multiplicity c the least v with exactly c representations
(`EXACT`), and the least v with at least n (`ATLEAST`, which is a(n) when the
search starts at 1). The covered range is contiguous from the start, so an
interrupted GPU run resumes from the "done below" value in its status line.

Throughput: GPU 1.6×10⁹ palindromes/s (M4 Max) and 2.6×10⁹/s (M2 Ultra) on
k=5; CPU ~3×10⁸/s on 10 laptop cores.

## Validation

- **gappal**: every known term of all 15 sequences checked, 1,184 in total:
  all b-file terms for L = 2–7 up to d = 200 (Chai Wah Wu's b-files) and all
  data terms for L = 8–17 (Merickel, 2012). All match, before and after the
  adaptive-depth change.
- **multipal (CPU)**: all 100 terms of A171702's b-file (Resta), all terms of
  A171703 (9), A171704 (4), A171705 (3), A171706 (3), A171701 (via k=2), and
  the A171741 / A171742 cross-sections. All match, except the two data errors
  below.
- **multipalg (GPU) vs multipal (CPU)**: full multiplicity histograms (count
  and least value for every c) identical on k=4 (10⁸ to 3.46×10¹²), k=5
  (10⁶ to 4.84×10⁹) and k=6 (10¹⁰ to 2.85×10¹²).

## Findings

**Data errors in OEIS entries** (both now corrected on OEIS)
- **A171701 a(41)** was listed as 30160; the correct value is **20160**. The
  printed value also exceeded a(42) = 20160, impossible for an "at least n"
  sequence. Corrected in rev 19 (Sep 30).
- **A171741 a(3)** was listed as 154; the correct value is **121**
  (232 in base 7, 171 in base 8, 121 in base 10). 154 (bases 6, 8, 9) is the
  second such number. A171702's own a(3) = 121 agrees. Corrected in rev 28
  (Oct 2).

For what is live, pending review, or not yet submitted on OEIS, see
`SUBMISSIONS.md`.

**New upper bounds** (verified; the entries currently give lower bounds only).
Each comes from v = L^(k−1), which in base L/d − 1 has digits
d^(k−1) × (a row of Pascal's triangle), a palindrome for every divisor d of L
small enough that no digit reaches the base:
- A171705(4) = A171742(6) ≤ 40968⁵ ≈ 1.154×10²³
- A171706(4) ≤ 327696⁶ ≈ 1.238×10³³
- A171741(8) ≤ 229644⁷ ≈ 3.368×10³⁷

**New terms / proofs**
- **A171704**: two new terms, the first since 2009. See the A171704 section below.
- **A171702**: from 100 terms to 223 (exhaustive to 10¹²); b-file in `A171702/`.
- **A171701**: b-file extended from 44 to 1,152 terms (exhaustive to 10¹⁰).
- **A171703**: from 9 proven terms to 22. See the A171703 section below.
- **Base-gap family**:
  - New b-files for A216902–A216905 (L = 8–11, about 1,000 terms each; OEIS
    had 10–19) and A216906–A216909 (L = 12–15, about 200 terms each; OEIS had
    7, 11, 4 and 5).
  - **The new length-16 sequence is complete for n = 2..60** (59 terms).
  - A216910 (L = 17): n = 1..60 (60 terms; OEIS had 3).
  - All of Merickel's 2012 values are confirmed.
  - **The family is complete**: every length from 8 to 17 now has a b-file.

## A171703: the first number that is a 4-digit palindrome in at least n bases

The entry had 9 terms (to 3456649728000 = 15120³) and a 2014 comment by
Hiroaki Yamanouchi giving upper bounds for a(10)–a(15). All six bounds turn
out to be exact, and the search has added seven more terms:

| n | a(n) | Bases | From the cube construction | Sporadic | Status |
|---|---|---|---|---|---|
| 10, 11 | 27653197824000 = 30240³ | 11 | 10 | 11906 | **proven** (was a bound) |
| 12 | 575100098496000 = 83160³ | 12 | 12 | — | **proven** (was a bound) |
| 13, 14 | 1733189185728000 = 120120³ | 14 | 13 | 18815 | **proven** (was a bound) |
| 15 | 8400090327552000 = 203280³ | 15 | 14 | 12671 | **proven** (was a bound) |
| 16 | 35158608576000000 = 327600³ | 16 | 16 | — | **new, proven** |
| 17 | 40254862491648000 = 342720³ | 17 | 16 | 43349 | **new, proven** |
| 18 | 205045005215232000 = 589680³ | 18 | 18 | — | **new, proven** |
| 19, 20 | 374368864117248000 = 720720³ | 20 | 20 | — | **new, proven** |
| 21 | 1168773800173056000 = 1053360³ | 21 | 21 | — | **new, proven** |
| 22 | 2994950912937984000 = 1441440³ | 22 | 22 | — | **new, proven** |

**The construction.** (b+1)³ is 1 3 3 1 in base b, and d³(b+1)³ is
d³ × (1 3 3 1) whenever 3d³ < b. So L³ is a 4-digit palindrome in base
L/d − 1 for every divisor d of L with 3d⁴ + d < L, and cubes of highly
composite numbers collect many bases. Every term from a(5) on is such a
cube, but coincidental "sporadic" bases decide which cube comes first. For
example, 342720³ has only 16 construction bases; its sporadic 17th puts
a(17) three times below the best pure-construction candidate for 17.

**Proof.**
- a(10)–a(15): exhaustive CPU sieve below 8.4×10¹⁵ (1.1×10¹² palindromes,
  615 s on the Studio), confirmed independently by the GPU sieve.
- a(16)–a(20): the GPU sieve continued contiguously from 8.4×10¹⁵, in two
  ranges split at 1.5×10¹⁷ (Studio, then laptop, up to 3.744×10¹⁷; 1.8×10¹³
  palindromes). No value below 1.5×10¹⁷ has 18 bases, and in the laptop's
  range the first with 18 is 589680³. Nothing has 19 until 720720³, which
  has exactly 20. Every reported value was re-checked on the host and
  independently in Python. Independent CPU sieves reproduced a(16)–a(18)
  (8.4×10¹⁵ to a(18)), a(19)–a(21) (a(18) to a(21)) and a(22) (a(21) to
  a(22), 4.7×10¹³ palindromes) exactly, so **every term through a(22) is
  confirmed on two architectures**.

- a(21): the GPU sieve continued to 1053360³ (Studio 3.744×10¹⁷ to
  8.7×10¹⁷, laptop above; 2.6×10¹³ palindromes). Nothing below 1053360³ has
  21 bases (the most is 20, at 720720³ and 837175343380992000), and 1053360³
  has exactly 21.

- a(22): the GPU sieve continued to 1441440³ (laptop 1.169×10¹⁸ to 2.1×10¹⁸,
  Studio above; 4.6×10¹³ palindromes). Nothing below 1441440³ has 22 bases
  (the most is 21), and 1441440³ has exactly 22.

**Paused (2026-09-30 ~07:45, by request): the a(23) search.** The
construction's bound is 1965600³ = 7594259452416000000, and two stretches
above a(22) are already covered:
- 2994950912937984001 to 3248644529349623809: at most 20 bases;
- 5200000000000000000 to 5347099262633639936: at most 21 bases (e.g.
  5257132155072000000).

To resume without redoing work, sweep [3248644529349623809,
5200000000000000000) and [5347099262633639936, 7594259452416000001) with
`multipalg2 -k 4 -l <lo> -u <hi> -r 21 -H 8`.

## A171704: the first number that is a 5-digit palindrome in at least n bases

The entry had four terms from 2009 and a 2015 comment that a(5) > 10¹⁴.

| n | a(n) | Bases, with digits | Status |
|---|---|---|---|
| 1 | 17 | 2: 1 0 0 0 1 | OEIS |
| 2 | 2293 | 5: 3 3 1 3 3 · 6: 1 4 3 4 1 | OEIS |
| 3 | 267140 | 14: 6 13 4 13 6 · 15: 5 4 2 4 5 · 18: 2 9 14 9 2 | OEIS |
| 4 | 4838419019 | 91, 125, 144, 161 | OEIS |
| 5 | **4922057407205376** = 8376⁴ | 2093, 2791, 4187, 5906, 8375 | **new, proven** |
| 6 | **34777153514704896** = 13656⁴ | 3218, 3413, 4551, 4828, 6827, 13655 | **new, proven** |

**Why these numbers.** (b+1)⁴ is 1 4 6 4 1 in base b, and d⁴(b+1)⁴ is
d⁴ × (1 4 6 4 1) as long as 6d⁴ < b. So L⁴ is a 5-digit palindrome in base
L/d − 1 for every small divisor d of L:

- **a(5) = 8376⁴.** 8376 = 2³·3·349, so d = 1, 2, 3, 4 give four bases:
  - 8375: 1 4 6 4 1
  - 4187: 16 64 96 64 16
  - 2791: 81 324 486 324 81
  - 2093: 256 1024 1536 1024 256

  The fifth base, 5906 (digits 4 268 4497 268 4), is a coincidence outside
  the pattern.
- **a(6) = 13656⁴.** 13656 = 2³·3·569, and d = 1, 2, 3, 4 give four bases:
  - 13655: 1 4 6 4 1
  - 6827: 16 64 96 64 16
  - 4551: 81 324 486 324 81
  - 3413: 256 1024 1536 1024 256

  Base 3218 comes from a different palindromic polynomial: the value equals
  81·(2b² + 3b + 2)² at b = 3218, with digits 324 972 1377 972 324.
  Base 4828 (digits 64 32 132 32 64) is sporadic.
- Doubling every digit keeps a palindrome as long as no digit reaches the
  base, so 2·13656⁴ = 69554307029409792 has the same six bases. The Studio's
  GPU range found it first, before the laptop's lower range reached 13656⁴.

The pure construction would only have reached a(5) ≤ 18780⁴ ≈ 1.24×10¹⁷
(five divisor bases, no coincidences). The exhaustive search found the true
a(5) 25 times lower, because 8376⁴ gets its fifth base by coincidence.

**Proof.**
- **a(5):** `multipal` (CPU) searched every 5-digit palindrome in every base
  below 10¹² and found at most 4 bases. `multipalg` (GPU) then swept upward
  from 10¹² contiguously. Its first value with 5 bases is 4922057407205376,
  re-checked on the host, and independently in Python.
- **a(6):** the same contiguous sweep (10¹² to 3.5×10¹⁶, 1.8×10¹³
  palindromes, 2.8 hours on the laptop GPU) found exactly one value with 6
  bases, 34777153514704896, and it lies below every other 6-base value found.
  The Studio GPU's range above 3.5×10¹⁶ found only its double.
- An independent CPU re-run from 10¹² up to a(6) (1.8×10¹³ palindromes in
  two parts, 2246 s and 7886 s on the Studio) reproduced both a(5) and
  a(6) exactly, so both terms meet the two-architecture standard.

**a(7) > 1.244×10¹⁷.** The full search covered every value below
124389107494560001 (5.0×10¹³ palindromes in all), and nothing there has 7
bases. Below that bound:
- exactly two values have 6 bases: a(6) = 13656⁴ and 2·13656⁴;
- 15 values have exactly 5 bases:
  - a(5) = 4922057407205376
  - 8451077831905536, 14471373983256576, 16677267707126016
  - 16889601600000000, 23520349636263936, 34899554288722176
  - 57600990724260096, 57779667567968256, 65490367914680831
  - 94307980945784831, 115559335135936512, 116997070082015231
  - 120617634521846016, and 124389107494560000 = 18780⁴ (the old
    construction bound).

## Status

As of 2026-09-30, 06:50 PDT. Results land in `A216/results/` and
`multipal/runs/`; b-files are rebuilt with `A216/make_bfiles.py`.

**Base-gap family**

| Run | Machine | State |
|---|---|---|
| L = 8–11, d ≤ 1000 | Studio CPU | **done** (1000 / 1000 / 999 / 1000 terms) |
| L = 12, 13, d ≤ 200 | laptop CPU | **done** (199 / 200 terms) |
| L = 14, d ≤ 200 | laptop CPU | **done** (199 terms; OEIS had 4). The hardest term, n = 151, needed 12.6 CPU-hours. |
| L = 15, d ≤ 200 | laptop CPU | **done** (200 terms; OEIS had 5) |
| L = 16 (new), d ≤ 60 | Studio CPU | **done** (59 terms) |
| L = 17, d ≤ 60 | Studio CPU | **done** (60 terms; OEIS had 3) |

**A171703 (4-digit palindromes in ≥ n bases)**

| Range | Machine | State |
|---|---|---|
| 1 → 8.4×10¹⁵ | Studio CPU + laptop GPU | **done: a(10)–a(15) proven on both** |
| 8.4×10¹⁵ → 1.169×10¹⁸ | both GPUs | **done: a(16)–a(21) proven** |
| 2.05×10¹⁷ → 1.169×10¹⁸, CPU re-check of a(19)–a(21) | Studio CPU | **done** (identical results; 3.34×10¹³ palindromes, 7.1 h) |
| 1.169×10¹⁸ → 2.995×10¹⁸, CPU re-check of a(22) | Studio CPU | **done** 2026-10-06 (identical results; 4.7×10¹³ palindromes) |
| 1.169×10¹⁸ → 2.995×10¹⁸ | both GPUs | **done: a(22) proven** |
| 2.995×10¹⁸ → 7.594×10¹⁸ (a(23) search) | both GPUs | **paused** by request; resume points in the A171703 section |

CPU re-checks already done: a(16)–a(18) (identical results).

**A171704 (5-digit palindromes in ≥ n bases)**: **done.** a(5) and a(6) are
proven and confirmed on both CPU and GPU; a(7) > 1.244×10¹⁷.

**A171702 (3-digit palindromes in ≥ n bases)**: **done: extended from 100
to 223 terms.** The dense-mode CPU sieve covered every value from Resta's
a(100) = 13967553601 to 10¹² (4.4×10¹² palindromes, 53 min on the laptop).
- a(101) = 16908091201
- a(221) = a(222) = a(223) = 963761198401, a 3-digit palindrome in 223 bases
- The new terms take 21 distinct values.
- Spot checks by independent base counting in Python: a(101), a(150),
  a(200), a(223).
- New b-file: `A171702/b171702.txt` (n = 1..223, Resta's terms kept).
