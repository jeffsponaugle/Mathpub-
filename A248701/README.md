# A248701 — primes with n increasing gaps before them and n decreasing gaps after

Tool for computing and extending [OEIS A248701](https://oeis.org/A248701):

> Smallest prime such that the preceding n prime gaps are increasing and the
> following n prime gaps are decreasing.

The same scan extends the three companion sequences that the entry
cross-references: [A248702](https://oeis.org/A248702) (a valley instead of a
peak), [A248703](https://oeis.org/A248703) (strict peak) and
[A248704](https://oeis.org/A248704) (strict valley).

## Definition

Write the primes as p₁ = 2 < p₂ < … and the gaps as gᵢ = pᵢ₊₁ − pᵢ. The prime
pₖ sits between the gap gₖ₋₁ *into* it and the gap gₖ *out of* it. Let

    L = length of the longest non-decreasing run of gaps that ends with g_(k-1)
    M = length of the longest non-increasing run of gaps that starts with g_k

The **peak depth** of pₖ is min(L, M), and A248701(n) is the smallest prime of
peak depth ≥ n. "Increasing" in the OEIS name is non-strict (a(2) = 7 has gaps
2, 2 | 4, 2) and nothing is required of gₖ₋₁ against gₖ; this is exactly what
the PARI and Maple programs in the entry test.

Example from the entry: 359 lies among 337, 347, 349, 353, **359**, 367, 373,
379, 383 with gaps 10, 2, 4, 6 | 8, 6, 6, 4. The run into 359 is 2 ≤ 4 ≤ 6
(the 10 before it breaks it), so L = 3; out of it 8 ≥ 6 ≥ 6 ≥ 4, so M = 4; the
depth is 3 and a(3) = 359.

The companions use the same gap stream:

| sequence | shape | runs | OEIS index | a(n) = smallest prime with |
|---|---|---|---|---|
| A248701 | peak | non-decreasing in, non-increasing out | offset 1 | peak depth ≥ n |
| A248702 | valley | non-increasing in, non-decreasing out | offset 0, a(0) = 2 | valley depth ≥ n |
| A248703 | strict peak | strictly increasing in, strictly decreasing out | offset 1 | strict peak depth ≥ n+1 |
| A248704 | strict valley | strictly decreasing in, strictly increasing out | offset 1 | strict valley depth ≥ n |

A248703's index counts strict *steps*: its a(1) = 23 has gaps 2, 4 | 6, 2,
two gaps and one strict step on each side, so its index is one less than the
depth. Each a(n) is non-decreasing in n, since depth ≥ n+1 implies depth ≥ n.

## Results

Everything below was computed with this tool on Sep 16–17 2026: the primes
below 1.357 × 10¹³ on an Apple M1 Pro (10 cores, sharing the machine with
another job), the rest on a Mac Studio resuming from the checkpoint. All
1,962,383,843,874 primes below 60233937494016 (6.02 × 10¹³) have been
examined, in 2335 s of scanning altogether, and every one of the 26 terms in
the four OEIS entries is reproduced.

**Three new terms of A248701** (the OEIS entry, keyword `more`, ends at a(7)):

| n | a(n) | gaps into a(n) | gaps out of a(n) |
|---|------|----------------|------------------|
| 8 | **12579905251** | 4 6 12 12 18 18 20 24 | 46 14 12 10 8 6 6 4 |
| 9 | **5108217950351** | 18 18 20 22 24 24 24 26 34 | 62 58 56 46 36 30 24 14 6 |
| 10 | **59852066157421** | 8 10 24 26 34 38 40 50 90 102 | 88 68 52 30 24 20 18 18 10 8 |

and **a(11) > 3.86 × 10¹⁵** (see below). The sequence is now

```
3, 7, 359, 7853, 96401, 2812099, 294276293, 12579905251, 5108217950351, 59852066157421
```

**New terms of the companions** from the same scan:

| sequence | new terms | gaps into it | gaps out of it |
|---|---|---|---|
| A248702 | a(7) = **938665577** | 94 80 24 12 12 12 4 | 2 10 12 12 18 36 36 |
|         | a(8) = **2400369437** | 76 30 26 24 12 12 10 6 | 6 6 8 10 30 32 34 36 |
|         | a(9) = **299917793009** | 80 70 44 36 36 24 10 6 6 | 8 12 12 12 16 24 24 24 26 |
|         | a(10) = **384671458489889** | 64 60 60 44 36 34 18 14 10 6 | 2 6 10 12 12 20 36 40 60 104 |
|         | a(11) = **3860678242735729** | 114 94 78 50 40 30 30 18 14 12 6 | 10 12 12 18 18 18 26 40 66 108 114 |
| A248703 | a(7) = **149822520893** | 6 10 12 30 36 54 60 72 | 78 56 40 20 16 14 10 2 |
|         | a(8) = **13193280477899** | 6 10 12 14 16 18 30 102 150 | 90 68 36 24 22 14 10 8 6 |
|         | a(9) = **746882911420231** | 6 8 12 16 20 22 32 34 36 38 | 130 48 38 36 18 16 14 12 10 8 |
| A248704 | a(7) = **3531448007** | 44 40 30 26 16 14 4 | 6 8 10 12 26 30 58 |
|         | a(8) = **17190066197** | 66 50 24 18 16 14 10 6 | 6 8 10 12 14 16 18 36 |
|         | a(9) = **37148264596189** | 62 58 48 32 28 26 18 16 2 | 4 6 12 18 20 28 30 60 122 |
|         | a(10) = **1958854030679863** | 70 56 42 40 38 30 12 10 8 6 | 16 20 22 26 30 54 58 74 88 168 |

The last four come from the continuing hunt for a(11) on the Mac Studio (24
threads, resumed from 6.02 × 10¹³ on Sep 17 2026, `scan_1e17_records.txt`).
Since the A248702(11) record was folded, every prime below 3860678242735729
has been examined, so **a(11) > 3.86 × 10¹⁵** for A248701 and likewise
A248702(12), A248703(10) and A248704(11) are all above that bound. The
A248702 entry only had the comment "a(7) >= 8960453, if it exists"; it exists.

```
A248702: 2, 3, 19, 43, 2687, 179819, 1107791, 938665577, 2400369437, 299917793009, 384671458489889, 3860678242735729
A248703: 23, 1439, 21433, 1130863, 19881311, 331542583, 149822520893, 13193280477899, 746882911420231
A248704: 3, 19, 1429, 25243, 340577, 1107791, 3531448007, 17190066197, 37148264596189, 1958854030679863
```

The scan output with every record below 10¹³ and its surrounding gaps is in
`scan_1e13.txt` and the Mac Studio's final summary in `scan_1e15_summary.txt`
(the record lines for A248704(9) and a(10) are in the Studio's
`scan_1e15.txt`). The terms are collected in `DATA.txt` and the b-files
`b248701.txt` … `b248704.txt`; `verify_py_windows.txt` lists all fifteen new
terms with their windows of consecutive primes.

### Depth statistics below 6.02 × 10¹³

Number of primes p < 60233937494016 of depth ≥ d (every prime except 2 has
depth ≥ 1; depths below 4 are only examined below 2³⁰, see `-d`):

| d | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| peak (A248701) | 4000374598 | 177063528 | 5565658 | 131963 | 2492 | 45 | 1 |
| valley (A248702) | 3914675711 | 169843380 | 5245379 | 122209 | 2204 | 31 | 0 |
| strict peak (A248703) | 2604151963 | 87029116 | 1936140 | 30587 | 362 | 3 | 0 |
| strict valley (A248704) | 2563256345 | 84140091 | 1836529 | 28707 | 300 | 3 | 0 |

For independent continuous gaps P(depth ≥ d) would be 1/(d!)², so the counts
would fall by a factor d² per step. The observed ratios
count(≥ d)/count(≥ d+1) for the peak shape are 22.6, 31.8, 42.2 and 53.0 for
d = 4..7, i.e. 0.83–0.90 d²: ties help a non-strict run. The strict shapes
go the other way, 29.9, 45.0, 63.3 and 84.5, i.e. 1.2–1.3 d², because a tie
breaks a strict run. The peak ratios 55 and 45 for d = 8, 9 rest on only 45
and 1 hits.

### Where a(11) should be

With 45 primes of peak depth ≥ 9 among the 1.96 × 10¹² primes below
6.02 × 10¹³, and ratios of about 80 and 100 for the next two steps, a prime
has depth ≥ 11 with probability about 3 × 10⁻¹⁵. That gives an expected 0.09
such primes below 10¹⁵, 0.8 below 10¹⁶, 2.4 below 3 × 10¹⁶ and 8 below 10¹⁷,
so a(11) is probably between 10¹⁵ and 3 × 10¹⁶ (median about 8 × 10¹⁵). The
Mac Studio scanned 3.5 × 10¹⁰ numbers/s around 5 × 10¹³, and the sieve slows
by about 1.7× toward 10¹⁶, so reaching 10¹⁶ takes about a week and 10¹⁷ about
two months. The checkpointed command, which also picks up the companions'
next terms on the way and stops by itself at a(11), is

    ./a248701 scan 1e17 -n 11 -S a248701.state > scan_1e17.txt 2> scan_1e17.log

It has been running on the Studio (24 threads) since Sep 17 2026. As of Sep 19
it had found A248702(10), A248703(9), A248704(10) and A248702(11) but no prime
of peak depth 11 below 3.86 × 10¹⁵ (expected number 0.33). Given that, the
remaining odds are about 39 % of a hit before 10¹⁶, 87 % before 3 × 10¹⁶ and
99.9 % before 10¹⁷.

### Performance

primesieve's iterator, single thread, plain prime enumeration:

| around | numbers/s | ns per number |
|---|---|---|
| 10⁹–10¹⁰ | 5 × 10⁹ | 0.20 |
| 10¹² | 4.2 × 10⁹ | 0.24 |
| 10¹⁴ | 3.3 × 10⁹ | 0.30 |
| 10¹⁶ | 2.5 × 10⁹ | 0.40 |
| 10¹⁸ | 1.0 × 10⁹ | 0.96 |

The depth bookkeeping costs a few percent on top of that with the default
`-d 4`. The 10¹³ scan ran at 1.6 × 10¹⁰ numbers/s with 10 threads while
another 9-thread job had the rest of the machine (the process got about four
cores); on an idle machine expect roughly twice that. On the Mac Studio the
continuation from 1.357 × 10¹³ to 6.02 × 10¹³ ran at 3.5 × 10¹⁰ numbers/s
with the default thread count.

## Verification

* `a248701 selftest` (1.5 s) checks the `show` depths of the OEIS examples,
  compares the chunk scanner with a straightforward array-based reference
  implementation on four 3 × 10⁶ ranges (at 0, 10⁹, 10¹² and 1.2 × 10¹⁴),
  checks that a chunk equals the sum of its halves, that a scan of [0, 3 × 10⁸)
  gives identical records and counts with one chunk, chunk 777777 and chunk
  2²⁰ on 10 threads, that a scan starting at 3 × 10⁸ adds up with one ending
  there, that `-d 4` and `-d 1` agree above 2³⁰ on every depth ≥ 4, that the
  `-n` early stop works, that the checkpoint file round-trips, and that the
  scan below 4 × 10⁸ reproduces all 26 known terms and the depth histogram
  computed independently with numpy.
* `verify_a248701.py` is an independent implementation (segmented numpy sieve,
  vectorised run lengths, pure-Python Miller–Rabin for single windows) sharing
  nothing with the C code but the definitions. `verify_a248701.py scan 1.72e10`
  reproduces a(8) = 12579905251 and the four companion terms below 1.72 × 10¹⁰,
  and its depth counts agree with the C tool's count for count (32 numbers,
  `cross_1.72e10.txt` vs `verify_py_1.72e10.txt`).
  `verify_a248701.py window P N` confirmed for every new term that the 2N+1
  listed numbers are consecutive primes (all numbers between them composite)
  with the stated gap pattern (`verify_py_windows.txt`).
* `verify_a248701.py scan 5.2e12 268435456 --procs=8` (79 minutes on a Mac
  Studio, `verify_py_5.2e12.txt`) independently reproduces a(8) = 12579905251
  and **a(9) = 5108217950351**, together with A248702(7..9), A248703(7) and
  A248704(7..8), and gives a(10) > 5.2 × 10¹². The C tool over the same range
  (`cross_5.2e12.txt`, 288 s) gives the same 184,126,901,158 primes, the same
  largest gap 652, and the same depth histogram for every shape and every depth
  from 4 to 9 (24 numbers, e.g. 230 primes of peak depth ≥ 8 and 1 of depth
  ≥ 9). The seven terms above the Python-scanned range (A248703(8), A248704(9),
  a(10), A248702(10), A248703(9), A248704(10) and A248702(11)) have verified
  windows, but their minimality rests on the C scan alone (the scan that
  agreed with Python on every count below 5.2 × 10¹²). Repeating the
  Python scan to 6.03 × 10¹³ would take about 15 hours with eight processes on
  the Studio: `verify_a248701.py scan 6.03e13 268435456 --procs=8`.

## Usage

```
make                        # needs primesieve: brew install primesieve
./a248701 selftest
./a248701 scan 1e13 -S a248701.state          # all four sequences below 10^13, checkpointed
./a248701 scan 1e17 -n 11 -S a248701.state    # continue, stop once A248701(11) is known
./a248701 show 5108217950351                  # depths of a prime and the gaps around it
python3 verify_a248701.py scan 1.72e10 --procs=4
python3 verify_a248701.py window 12579905251 8
```

`scan [START] END` prints every new record of depth ≥ 6 (`-r`) as it is found,
then the depth-count table and the four a(n) tables with lower bounds for the
first unknown term. `-S FILE` writes a checkpoint every 60 s (`-i`) and on
Ctrl-C; rerunning the same command (with the same or a larger END) resumes
from it. `-d DMIN` (default 4) is the smallest depth examined above 2³⁰:
examining every depth costs about as much as the sieve itself, and the small
depths are all known; chunks starting below 2³⁰, which hold every known term,
always examine everything. `-t` sets the thread count (default: all cores),
`-c` the chunk size (default: range/(64·threads), clamped to [2²⁰, 2³⁴]).

## How it works

Worker threads take chunks from an atomic counter. A thread sieves its chunk
with primesieve, starting 256 primes below it and running 256 primes past it so
the run lengths at the edges are exact, and feeds the gaps to an O(1)-per-gap
state machine: the lengths of the non-decreasing, non-increasing, strictly
increasing and strictly decreasing runs ending at the current gap, plus rings
of the last 256 run lengths and primes. When the run ending at gap j has
length r, every n ≤ r names the centre k = j − n + 1 whose following n gaps
are that run's tail; the centre has depth ≥ n exactly when the opposite run
ending at gap k − 1 is ≥ n. Each (centre, n) pair is examined once, so a chunk
yields, for each shape and depth, the number of centres of depth ≥ d and the
smallest one. Chunk results are folded into the global totals in increasing
order behind a completion frontier, so the first record folded is the
smallest, the totals always describe an initial segment of the range, and a
checkpoint is just the frontier plus the totals.

## Files

* `a248701.c`, `Makefile` — the tool
* `verify_a248701.py` — independent Python/numpy check
* `scan_1e13.txt`, `scan_1e13.log`, `scan_1e13.state` — the complete scan below 10¹³
* `scan_1e15.txt`, `scan_1e15.log`, `a248701.state` — the M1 Pro leg of the continuation (to 1.357 × 10¹³); the complete files live on the Mac Studio
* `scan_1e15_summary.txt` — the Mac Studio's summary of the continuation, complete below 6.02 × 10¹³
* `scan_1e17_records.txt` — records found so far by the a(11) hunt on the Studio (complete below 3.86 × 10¹⁵ as of Sep 19 2026)
* `cross_1.72e10.txt`, `verify_py_1.72e10.txt`, `cross_5.2e12.txt`, `verify_py_5.2e12.txt`, `verify_py_windows.txt` — cross-checks
* `selftest.txt` — selftest output
* `DATA.txt`, `b248701.txt`, `b248702.txt`, `b248703.txt`, `b248704.txt` — terms

## Notes for the OEIS submissions

* A248701: add a(8) = 12579905251, a(9) = 5108217950351, a(10) = 59852066157421;
  b-file n = 1..10; comment "a(11) > 3.86*10^15" (raise it to wherever the hunt
  has got to, or add a(11)); the `more` keyword stays.
  A comment that the monotonicity is weak and that the two gaps adjacent to
  a(n) are not compared (a(7) has 60 | 24 around it, a(8) has 24 | 46) would
  remove the ambiguity in the name; A248703 is the strict version.
* A248702: add a(7) = 938665577, a(8) = 2400369437, a(9) = 299917793009,
  a(10) = 384671458489889, a(11) = 3860678242735729; b-file n = 0..11; replace
  the comment "a(7) >= 8960453, if it exists" by "a(12) > 3.86*10^15".
* A248703: add a(7) = 149822520893, a(8) = 13193280477899, a(9) = 746882911420231;
  "a(10) > 3.86*10^15".
  Worth a comment that the index counts strict steps, so a(n) has n+1 gaps on
  each side.
* A248704: add a(7) = 3531448007, a(8) = 17190066197, a(9) = 37148264596189,
  a(10) = 1958854030679863; "a(11) > 3.86*10^15".
