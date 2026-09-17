# A335406 — first run of n equal prime gaps

Tool for computing [OEIS A335406](https://oeis.org/A335406):

> First position of n in the sequence of run-lengths of the sequence of prime gaps.

Take the primes p₁ = 2 < p₂ = 3 < … and their gaps gᵢ = pᵢ₊₁ − pᵢ (A001223), and cut
the gap sequence into maximal runs of equal values:

    (1), (2,2), (4), (2), (4), (2), (4), (6), (2), (6), (4), (2), (4), (6,6), (2), …

The run lengths are A333254 = 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, …, and
a(n) is the index of the first run of length exactly n. The OEIS entry (Sep 2026,
keyword `hard,more`) has a(1..5) = 1, 2, 49, 633353, 6706139.

A run of n equal gaps d beginning at gap gᵢ is the same thing as n+1 consecutive
primes pᵢ, pᵢ+d, …, pᵢ+nd in arithmetic progression (a "CPAP-(n+1)"). So the first
run of length n begins at the smallest CPAP-(n+1), A006560(n+1); its gap index i is
A089180(n−1); and because a new run begins exactly at the gaps that differ from
their predecessor,

    a(n) = i − E(i−1),        E(m) = #{ j ≤ m : gⱼ = gⱼ₊₁ }.

Example: 251, 257, 263, 269 are the first four consecutive primes in arithmetic
progression. 251 = p₅₄, and among g₁..g₅₃ exactly five gaps equal their successor
(j = 2, 15, 36, 39, 46), so a(3) = 54 − 5 = 49.

## Results (Apple M1 Pro, 10 cores, Sep 16 2026)

**Known terms reproduced**, with the gap index, start prime and common gap of each run:

| n | a(n) | gap index i (A089180(n−1)) | first prime (A006560(n+1)) | gap |
|---|------|------|------|-----|
| 1 | 1 | 1 | 2 | 1 |
| 2 | 2 | 2 | 3 | 2 |
| 3 | 49 | 54 | 251 | 6 |
| 4 | 633353 | 654926 | 9843019 | 30 |
| 5 | 6706139 | 6904737 | 121174811 | 30 |

(A089180's DATA gives 654926 for the third term; the "659426" in that entry's
EXAMPLE is a transposition typo.)

**No new term is computable, but a rigorous lower bound is.** Every prime below
10¹³ was enumerated (346,065,536,839 primes, 396 s); the longest run of equal gaps
is 5, so the first run of length 6 begins beyond 10¹³ and its index exceeds the
number of runs below 10¹³:

    a(6) > 339,595,430,914.

Cumulative statistics at the powers of ten (all exact; π(10ᵏ) agrees with the known
values; runs of length 3, 4, 5 are maximal CPAP-4, CPAP-5, CPAP-6):

| x | π(x) | E(x) | runs | length 3 | length 4 | length 5 |
|---|------|------|------|------|------|------|
| 10⁷ | 664,579 | 21,837 | 642,742 | 1,567 | 1 | 0 |
| 10⁸ | 5,761,455 | 167,032 | 5,594,423 | 10,211 | 8 | 0 |
| 10⁹ | 50,847,534 | 1,328,401 | 49,519,133 | 71,962 | 83 | 1 |
| 10¹⁰ | 455,052,511 | 10,801,518 | 444,250,993 | 521,186 | 940 | 19 |
| 10¹¹ | 4,118,054,813 | 89,648,445 | 4,028,406,368 | 3,920,209 | 10,369 | 239 |
| 10¹² | 37,607,912,018 | 756,279,950 | 36,851,632,068 | 30,263,946 | 101,738 | 2,527 |
| 10¹³ | 346,065,536,839 | 6,470,105,925 | 339,595,430,914 | 239,112,966 | 970,975 | 24,383 |

The fraction E/π of primes whose gap repeats the previous gap falls slowly
(2.90% at 10⁸, 1.87% at 10¹³; roughly 0.6/ln x). The full (length, gap) tables are in
`scan_1e12.txt` and `scan_1e13.txt` (lines `ld POS LEN GAP COUNT`); for instance all
24,383 runs of five equal gaps below 10¹³ have gap 30 (24,315), 60 (67) or 90 (1).

## Why a(6) is out of reach

a(6) is the run index of the smallest CPAP-7. Two independent obstacles:

1. **The smallest CPAP-7 is unknown.** Its common difference is a multiple of 210,
   so six consecutive prime-free stretches of 209 integers are needed. The smallest
   known CPAP-7 starts at 71137654873189893604531 ≈ 7.1·10²² (P. Zimmermann, 2022;
   OEIS A006560 comment) and was found by a targeted search, not an exhaustive one.
2. **Even with the prime known, a(6) needs E at that point**, i.e. the number of
   repeated gaps among all primes below it. Unlike π(x), which combinatorial
   algorithms (primecount) reach at 10²² in hours, no method for E(x) is known other
   than enumerating every prime gap.

`model_a335406.py` estimates where the first CPAP-7 is. It integrates a
Hardy–Littlewood/Cramér random model with a wheel modulo 30030 (exact treatment of
the residue classes mod 2·3·5·7·11·13 for the seven members, the 1254 interior
integers and the two neighbouring gaps that make the run maximal; primes > 13 through
the singular series and per-position corrections). Checked against the exact
(length, gap) counts from the scan, the raw model is accurate to 1% for small gaps
but over-predicts runs of large gaps; the deficit is regular — for single gaps
ln(observed/model) = −F(u)/ln x with u = d/ln x, F ≈ 0.29·u^1.78, and each further
link of a chain costs another factor exp(−G(u)/ln x), G ≈ 0.46·u^1.6 — and after
fitting F and G the model reproduces every chain with ≥ 100 occurrences to within
about 10% (runs of 4 and 5 gaps of 30 and 60; runs of 2–3 gaps of 120–210; the
CPAP-7 regime, d = 210 near 10²⁰, has u ≈ 4.5, inside the calibrated range).

| x | expected maximal CPAP-7 below x (calibrated) | P(first CPAP-7 < x) |
|---|------|------|
| 10¹⁸ | 4.7·10⁻⁵ | 0.00005 |
| 10¹⁹ | 0.0024 | 0.002 |
| 10²⁰ | 0.098 | 0.09 |
| 10²¹ | 3.4 | 0.97 |
| 10²² | 102 | 1 |
| 7.1·10²² (Zimmermann) | ≈ 1700 | 1 |

Quantiles for the start P of the first CPAP-7: 10% 1.0·10²⁰, 50% 3.5·10²⁰, 90%
7.7·10²⁰ (the raw, uncalibrated model gives 0.5–3.8·10²⁰). About 1700 CPAP-7 are
expected below the known one, so it is almost certainly not the smallest. The
implied a(6) ≈ li(P)·(1 − E/π) is **2·10¹⁸ to 2·10¹⁹** (median 7.5·10¹⁸).

Cost of computing it: the scan runs at 2.5·10¹⁰ numbers/s on this machine (10
threads; primesieve slows further beyond 10¹⁸, so these are lower bounds):

| sieve to | time at 2.5·10¹⁰/s | what it buys |
|---|---|---|
| 10¹⁴ | 1.1 h | a(6) > ≈3.2·10¹² |
| 10¹⁵ | 11 h | a(6) > ≈3·10¹³ |
| 10¹⁶ | 4.6 days | a(6) > ≈2.8·10¹⁴ |
| 10¹⁸ | 1.3 years | a(6) > ≈2.4·10¹⁶ |
| 1.0·10²⁰ (10% quantile) | 130 years | maybe a(6) |
| 3.5·10²⁰ (median) | 440 years | probably a(6) |
| 7.1·10²² (known CPAP-7) | 90,000 years | a(6) for certain |

A thousand-core cluster would bring the median case to about five months, so a(6)
is a large distributed-computing project, not a workstation computation, and a(7)
(first CPAP-8, expected far beyond 10²⁵) is hopeless. Longer runs here only raise
the lower bound; I did not start any.

## Tool

`a335406.c` — single file, C11, pthreads, primesieve; exact integer counting, no
primality tests. The range is cut into chunks; workers walk the primes of a chunk with
a primesieve iterator, starting from the last prime below the chunk so the gap into
the first prime is known, and count primes, repeated gaps, runs by (length, gap) and
the first run of each length that begins in the chunk (as offsets). A run that begins
in a chunk is followed past its end; one that begins earlier belongs to the previous
chunk. Finished chunks are folded into the totals in order, which turns offsets into
absolute indices via a(n) = i − E(i−1) and prints each a(n) as soon as it is known.

    make                       # cc -O2 ... -lprimesieve
    ./a335406 selftest         # a(1..5), A089180, A006560, first 60 terms of A333254,
                               # chunked scan vs reference pass for several chunk sizes,
                               # thread counts, marks and split ranges (1.5 s)
    ./a335406 scan 1e13 -m pow10 -H    # exhaustive scan, totals at every power of ten,
                                       # (length, gap) table; -t T threads, -c CHUNK, -q
    ./a335406 runs 60          # first 60 terms of A333254
    ./a335406 naive 2e8 -H     # single-pass reference implementation

Numbers may be written as `1e13`, `10^13`, `2^40`, `X+Y`, `X−Y`.

`verify_a335406.py [N]` — independent check: numpy sieve, its own run logic, compares
every reported quantity (π, E, runs, histogram, (length, gap) table, first
occurrences) with `a335406 naive` and `a335406 scan`, and with the OEIS terms. Passes
at N = 2·10⁸.

`model_a335406.py SCANFILE [RATE]` — the heuristic model and its calibration; the
output for the 10¹³ scan is `model_1e13.txt`.

Files: `scan_1e12.txt`, `scan_1e13.txt` (+ `.log`) — scan outputs; `model_1e13.txt`;
`DATA.txt` — the known terms.

## Possible OEIS additions

- A335406: "a(6) > 339595430914: no run of six equal gaps begins at a prime below
  10^13 (exhaustive). a(6) is the index of the run starting at the smallest CPAP-7,
  A006560(7), and equals A089180(5) − E where E counts the j < A089180(5) with
  A001223(j) = A001223(j+1)." (a(n) = A089180(n−1) − #{j < A089180(n−1): A001223(j) = A001223(j+1)}.)
- A089180: a(5) > π(10^13) = 346065536839; EXAMPLE typo 659426 → 654926.
- The counts of runs by length at the powers of ten (maximal CPAP-k counts) may
  match or extend existing sequences; not checked.
