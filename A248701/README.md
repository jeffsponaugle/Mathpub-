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

Everything below was computed with this tool on an Apple M1 Pro (10 cores,
Sep 16 2026, sharing the machine with another job). The scan of all primes
below 10¹³ (346,065,536,839 primes) took 630 s and reproduces every one of the
26 terms in the four OEIS entries.

**Two new terms of A248701** (the OEIS entry, keyword `more`, ends at a(7)):

| n | a(n) | gaps into a(n) | gaps out of a(n) |
|---|------|----------------|------------------|
| 8 | **12579905251** | 4 6 12 12 18 18 20 24 | 46 14 12 10 8 6 6 4 |
| 9 | **5108217950351** | 18 18 20 22 24 24 24 26 34 | 62 58 56 46 36 30 24 14 6 |

and **a(10) > 10¹³**. The sequence is now

```
3, 7, 359, 7853, 96401, 2812099, 294276293, 12579905251, 5108217950351
```

**New terms of the companions** from the same scan:

| sequence | new terms | gaps into it | gaps out of it |
|---|---|---|---|
| A248702 | a(7) = **938665577** | 94 80 24 12 12 12 4 | 2 10 12 12 18 36 36 |
|         | a(8) = **2400369437** | 76 30 26 24 12 12 10 6 | 6 6 8 10 30 32 34 36 |
|         | a(9) = **299917793009** | 80 70 44 36 36 24 10 6 6 | 8 12 12 12 16 24 24 24 26 |
| A248703 | a(7) = **149822520893** | 6 10 12 30 36 54 60 72 | 78 56 40 20 16 14 10 2 |
| A248704 | a(7) = **3531448007** | 44 40 30 26 16 14 4 | 6 8 10 12 26 30 58 |
|         | a(8) = **17190066197** | 66 50 24 18 16 14 10 6 | 6 8 10 12 14 16 18 36 |

with A248702(10) > 10¹³, A248703(8) > 10¹³ and A248704(9) > 10¹³. The
A248702 entry only had the comment "a(7) >= 8960453, if it exists"; it exists.

```
A248702: 2, 3, 19, 43, 2687, 179819, 1107791, 938665577, 2400369437, 299917793009
A248703: 23, 1439, 21433, 1130863, 19881311, 331542583, 149822520893
A248704: 3, 19, 1429, 25243, 340577, 1107791, 3531448007, 17190066197
```

The full scan output with every record and its surrounding gaps is in
`scan_1e13.txt`; the terms are collected in `DATA.txt` and the b-files
`b248701.txt` … `b248704.txt`.

### Depth statistics below 10¹³

Number of primes p < 10¹³ of depth ≥ d (each prime except 2 has depth ≥ 1;
depths below 4 are only examined below 2³⁰, see `-d`):

| d | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|
| peak (A248701) | 711387410 | 31635394 | 999414 | 23913 | 442 | 6 |
| valley (A248702) | 695705113 | 30316986 | 941602 | 21986 | 405 | 7 |
| strict peak (A248703) | 452401589 | 14964865 | 328619 | 5162 | 68 | 0 |
| strict valley (A248704) | 445357981 | 14469023 | 312373 | 4756 | 62 | 0 |

For independent continuous gaps P(depth ≥ d) would be 1/(d!)², i.e. the count
would fall by a factor d² per step. The observed ratios are about 0.8 d² for
the non-strict shapes (ties help a non-strict run) and about d² for the strict
ones: 22 → 32 → 42 → 54 → 74 for the peak counts above.

### Where a(10) should be

Six primes of peak depth ≥ 9 below 10¹³ (one per 5.8 × 10¹⁰ primes), and a
ratio of about 80 for the next step, give about 0.6 primes of depth ≥ 10 in
[10¹³, 10¹⁴] and about 6 in [10¹⁴, 10¹⁵]. So a(10) is below 10¹⁴ with
probability about 45 % and below 10¹⁵ with probability above 99.8 %; the
expected position is around 10¹⁴. Each further term costs another factor of
roughly 100 in the search bound: a(11) should be near 10¹⁶ (weeks of sieving at
the rates below) and a(12) beyond 10¹⁸.

A checkpointed continuation of the scan toward 10¹⁵ with `-n 10` (stop as
soon as a(10) is known) was started on Sep 16 2026; see `scan_1e15.txt`.

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
cores); on an idle machine expect roughly twice that.

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
* A full independent Python scan to 5.2 × 10¹² (for the minimality of a(9))
  runs at low priority with four processes; its output goes to
  `verify_py_5.2e12.txt`.

## Usage

```
make                        # needs primesieve: brew install primesieve
./a248701 selftest
./a248701 scan 1e13 -S a248701.state          # all four sequences below 10^13, checkpointed
./a248701 scan 1e15 -n 10 -S a248701.state    # continue, stop once A248701(10) is known
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
* `scan_1e15.txt`, `scan_1e15.log`, `a248701.state` — the continuation toward a(10)
* `cross_1.72e10.txt`, `verify_py_1.72e10.txt`, `verify_py_windows.txt`, `verify_py_5.2e12.txt` — cross-checks
* `selftest.txt` — selftest output
* `DATA.txt`, `b248701.txt`, `b248702.txt`, `b248703.txt`, `b248704.txt` — terms

## Notes for the OEIS submissions

* A248701: add a(8) = 12579905251, a(9) = 5108217950351; b-file n = 1..9;
  comment "a(10) > 10^13" (to be raised when the continuation finishes); the
  `more` keyword stays.
* A248702: add a(7) = 938665577, a(8) = 2400369437, a(9) = 299917793009;
  replace the comment "a(7) >= 8960453, if it exists" by "a(10) > 10^13".
* A248703: add a(7) = 149822520893; "a(8) > 10^13".
* A248704: add a(7) = 3531448007, a(8) = 17190066197; "a(9) > 10^13".
