# A053686 — record prime gaps that repeat before the next record

Tool to reproduce and extend [OEIS A053686](https://oeis.org/A053686):

> Record gaps between consecutive primes that repeat at least once before a new record occurs.

OEIS terms (12, keyword `more`): 2, 4, 6, 14, 34, 36, 52, 86, 132, 154, 250, 336.
**New here: a(13) = 906**, and a(14) ≥ 1132 (see Results).

Let g₁ = 1 < g₂ = 2 < g₃ = 4 < … be the record ("maximal") prime gaps
([A005250](https://oeis.org/A005250)) and Pₙ the prime after which gₙ first occurs
([A002386](https://oeis.org/A002386)). The record gₙ is a term iff it occurs a second
time before the next record, i.e. iff some prime p with Pₙ < p < Pₙ₊₁ is followed by
the gap gₙ. Equivalently (David Wilson's comment) the terms are the values that occur
at least twice in [A085237](https://oeis.org/A085237), the nondecreasing gaps.

| k | a(k) | record n | first occurrence Pₙ (A133788) | second occurrence after | repeats before the next record |
|--:|--:|--:|--:|--:|--:|
| 1 | 2 | 2 | 3 | 5 | 1 |
| 2 | 4 | 3 | 7 | 13 | 2 |
| 3 | 6 | 4 | 23 | 31 | 6 |
| 4 | 14 | 6 | 113 | 293 | 2 |
| 5 | 34 | 10 | 1327 | 8467 | 1 |
| 6 | 36 | 11 | 9551 | 12853 | 2 |
| 7 | 52 | 13 | 19609 | 25471 | 1 |
| 8 | 86 | 15 | 155921 | 338033 | 1 |
| 9 | 132 | 20 | 1357201 | 1561919 | 1 |
| 10 | 154 | 22 | 4652353 | 11113933 | 2 |
| 11 | 250 | 29 | 387096133 | 428045491 | 1 |
| 12 | 336 | 34 | 3842610773 | 4275912661 | 1 |
| **13** | **906** | 61 | 218209405436543 | 543684371469023 | 1 |

(The "second occurrence" and "repeats" columns are computed here and are not in the
OEIS entry; a(13) is new.) The other 49 records below 924 — 1, 8, 18, 20, 22, 44, 72,
96, 112, 114, 118, 148, 180, 210, 220, 222, 234, 248, 282, 288, 292, 320, 354, 382, 384,
394, 456, 464, 468, 474, 486, 490, 500, 514, 516, 532, 534, 540, 582, 588, 602, 652,
674, 716, 766, 778, 804, 806, 916, 924 — do not repeat.

## Could it be extended?

Yes — a one-day scan decided nine more records and found one new term. The situation in
the OEIS was a little odd:

* All maximal gaps below 2⁶⁴ are known (80 records, the last being 1550 after
  1.836·10¹⁹; the table is embedded in the tool), so each record n is a yes/no
  question — does gₙ occur again in (Pₙ, Pₙ₊₁)? — that a sieve over that interval
  settles.
* A085237 was computed to its term 778 (Donovan Johnson 2008, Charles Greathouse 2011)
  and lists 354, 382, …, 766 exactly once each. That means records 35..57 do **not**
  repeat and **a(13) ≥ 778**, but this was never carried over into A053686, whose last
  extension (a(12), 2008) stops at 4.3·10⁹.
* Deciding the next six records — 778, 804, 806, 906, 916, 924 — needs every gap up to
  P₆₄ = 1693182318746371 ≈ 1.7·10¹⁵: about **22 hours** on this M1 Pro, six on a Mac
  Studio. Record 64 (1132) would need 4.4·10¹⁶ (24 days here, about a week on the
  Studio); the search for maximal gaps itself took a distributed effort to reach 2⁶⁴.

The run (below) settled all six: 906 repeats, the other five do not, so a(13) = 906 and
a(14) ≥ 1132; on the way it independently re-verified the part of A085237 that A053686
silently depends on.

## Method

Deciding record n means looking at every prime gap in [Pₙ, Pₙ₊₁), so the cost is one
pass of a prime sieve over the interval. Nothing cleverer applies: the target is one
exact gap size in the extreme tail, and a partial sieve (small primes only) leaves ~5%
survivors against ~3% primes, so it would not cut the primality work enough to pay for
itself. primesieve sieves at ~2·10⁹ numbers/s per core near 10¹⁵ on an M1 Pro and the
per-prime bookkeeping here is a subtraction and a compare.

The range is cut into chunks (default 10¹⁰) handed to pthreads through an atomic counter;
a worker sieves its chunk with a primesieve iterator, reads primes straight from the
iterator's buffer and looks at every gap (p, q) with p in the chunk (the gap crossing the
chunk end belongs to the chunk that owns p, so nothing is missed or double counted).
Gaps within 198 of the record of their interval are rare and get a closer look:

* the gap after Pₙ₊₁ must be exactly gₙ₊₁ — every record in the range has to turn up
  where A002386 says, which is a strong check on the sieve and on the table;
* a gap equal to gₙ after a prime p ≠ Pₙ is a repeat, i.e. gₙ is a term;
* a gap larger than gₙ would be an unlisted maximal gap or a sieve error (`ANOMALY`);
* all of them go into a per-record histogram of gaps gₙ − 2k, k < 100, which shows how
  close the near misses come and feeds the model below.

Completed chunks are folded into the totals in order (the frontier); that is when their
`REPEAT` lines are printed and when the record interval that ends in the chunk is
announced `DECIDED`. With `-S FILE` the position scanned to, the prime count, the
histogram and the repeat positions are checkpointed (every `-i` seconds and on Ctrl-C);
rerunning with the same START resumes from that position, and END and the chunk size
may be changed between runs, so a run can be stopped at 10¹³ and continued to 10¹⁵ later.
The prime count over [START, 10ᵏ) is compared with π(10ᵏ) − π(START−1) at the end.

## Results (Sep 16–17 2026)

**a(13) = 906.** The record gap 906 first occurs after P₆₁ = 218209405436543 and occurs
again after **543684371469023** (543684371469023 + 906 = 543684371469929 is the next
prime), before the record 916 after P₆₂ = 1189459969825483. Records 58, 59, 60, 62, 63
(778, 804, 806, 916, 924) do not repeat before their next records, so **a(14) ≥ 1132**;
record 64 (1132, after 1693182318746371) would need the scan carried to 4.4·10¹⁶.

Checks on the new term, all independent of the scan:

* `a053686 verify 543684371469023 906`: deterministic Miller–Rabin next-prime search and
  GMP `mpz_nextprime` both give 543684371469929 as the next prime;
* `verify_a053686.py 543684371469023 906`: pure-Python Miller–Rabin, same;
* `primesieve 543684371469023 543684371469929 -p` lists exactly those two primes;
* `openssl prime` confirms both primes.

Checks on the scan itself: all 30 record gaps in [4302407359, 1693182318746372) were found
exactly where A002386 lists them, no gap anywhere exceeded the record of its interval,
and the number of primes examined, **49,749,425,527,899**, equals
π(1693182318746371) − π(4302407358) computed with primecount (49749629143526 −
203615627). Every number in the range was therefore sieved and every gap seen once.

The run: `scan r35 1e13` on this M1 Pro (5 m 20 s, 3.1·10¹⁰ numbers/s, `scan_1e13.txt`),
then the checkpoint was moved to a Mac Studio, which finished
`scan r35 1693182318746372 -S a053686.state -H -i 60` at 7.7–8.2·10¹⁰ numbers/s in
about six hours. Its final table is in `scan_1.69e15_summary.txt`; the full log and
state are on that machine. The new part:

| n | gap | first after | next record after | repeats | runner-up gap (count) | gaps within 100 of the record | |
|--:|--:|--:|--:|--:|--:|--:|:--|
| 55 | 674 | 7177162611713 | 13829048559701 | 0 | 662 (1) | 49 | agrees with A085237 |
| 56 | 716 | 13829048559701 | 19581334192423 | 0 | 700 (1) | 17 | agrees with A085237 |
| 57 | 766 | 19581334192423 | 42842283925351 | 0 | 744 (1) | 26 | agrees with A085237 |
| 58 | 778 | 42842283925351 | 90874329411493 | 0 | 758 (1) | 29 | new |
| 59 | 804 | 90874329411493 | 171231342420521 | 0 | 788 (1) | 34 | new |
| 60 | 806 | 171231342420521 | 218209405436543 | 0 | 776 (1) | 45 | new |
| **61** | **906** | 218209405436543 | 1189459969825483 | **1** | 880 (1) | 52 | **new term, repeat after 543684371469023** |
| 62 | 916 | 1189459969825483 | 1686994940955803 | 0 | 902 (1) | 56 | new |
| 63 | 924 | 1686994940955803 | 1693182318746371 | 0 | 790 (1) | 0 | new |

Records 35..54 (354 … 652), scanned in the first leg, do not repeat either (table in
`scan_1e13.txt`), so A085237's implication for records 35..57 is confirmed
independently. Note how close the near misses come: the runner-up in 906's interval was
880, in 916's 902, in 778's 758.

### Was a new term likely?

`tail_a053686.py` fits the histogram of near misses: in each interval the number of gaps
of size g below the record is modelled as c′(g)·A·exp(b·(gₙ − g)), exponential in the
deficit below the record times the Hardy–Littlewood factor
c′(g) = ∏_{p | g, p > 2} (p−1)/(p−2) (gaps divisible by 3 are twice as common as their
neighbours), and the expected number of repeats is c′(gₙ)·A. Fitted on records 35..54
the slope is b ≈ 0.049 (counts double every 14 below the record; 1/b ≈ 20, well short
of ln x ≈ 25–30, which is why the naive Cramér model badly overpredicts large gaps).
Calibration: for records 20..54 the model expects 3.4 repeats and 5 occurred; for
records 35..54 it expects 1.7 and 0 occurred (Poisson probability 18%). Before the long
run it gave a 40–50% chance of at least one new term among 778..924 and named 804, 906
and 924 (divisible by 3) as the likely ones; 906 was the one. Run it on the Studio's
`a053686.log` (which has the `hist` lines) to see the fitted expectation for 906's
interval itself.

### Cost of going further

At 8·10¹⁰ numbers/s (Mac Studio), from P₆₄:

| record n | gap | needs the scan to | time |
|--:|--:|--:|--:|
| 64 | 1132 | 43841547845541059 | 6 days |
| 65 | 1184 | 55350776431903243 | 8 days |
| 66 | 1198 | 80873624627234849 | 12 days |
| 67 | 1220 | 203986478517455989 | 29 days |

so a(14) (whether 1132 repeats) is a week on the Studio, but the sieve slows beyond 10¹⁶
and the projections are optimistic. The command, continuing from the same checkpoint:

    caffeinate -i ./a053686 scan r35 43841547845541060 -S a053686.state -H -i 60 | tee -a a053686.log

## Tool

`a053686.c` — single file, C11, pthreads, primesieve; GMP optional (adds a second
next-prime check to `verify`).

    make && make test                 # selftest, ~10 s
    ./a053686 scan r35 1e13 -S a053686.state -H     # the run above (START may be rN = A002386(N))
    ./a053686 scan 0 1e11 -q -H > scan_0_1e11.txt   # the known part, with histogram lines
    ./a053686 verify 4275912661 336   # is (P, P+G) a prime gap, and a repeat of its record?
    ./a053686 bench 1e15              # throughput and the time to each decision
    ./tail_a053686.py scan_1e13.txt   # near-tail model from the hist lines

Options for `scan`: `-t T` threads (default all cores), `-c CHUNK` (default 10¹⁰),
`-n NEAR` also print gaps within NEAR of the record, `-S FILE` checkpoint/resume
(same START; END and CHUNK may change), `-i SECS` checkpoint and progress interval
(default 60), `-H` histogram lines, `-q` no status line, `-s KIB` primesieve sieve size. Numbers may be written as
`1e13`, `10^13`, `2^40`, `X+Y`, `X-Y`, `rN`.

`verify_a053686.py [P G | sieve N | table]` — independent check in Python/numpy
(default: table + sieve to 2·10⁸, comparing with `a053686 scan 0 N -c 1e6 -q -H`).

Files: `scan_1e13.txt`, `scan_1e13.log` — the 10¹³ leg on the M1 Pro;
`scan_1.69e15_summary.txt` — final table of the full run (Mac Studio); `scan_0_1e11.txt`
— the known part with histogram lines; `a053686.state` — checkpoint (on the Studio it
covers the full range and is the one to continue from); `DATA.txt` — the terms.

## Possible OEIS additions

- A053686: **a(13) = 906**. Example/comment: "906, the record gap after 218209405436543,
  occurs again after 543684371469023 (the next prime is 543684371469929), before the
  next record 916 after 1189459969825483." And: "a(14) ≥ 1132: the records 778, 804,
  806, 916 and 924 do not repeat before the next record (exhaustive scan to
  A002386(64) = 1693182318746371)."
- A133788: a(13) = 218209405436543.
- A085237: a(80)..a(86) = 804, 806, 906, 906, 916, 924, 1132 (778 is not repeated, 906
  is; a(87) is 1132 again or 1184).
- The second occurrences 5, 13, 31, 293, 8467, 12853, 25471, 338033, 1561919, 11113933,
  428045491, 4275912661, 543684371469023 (the prime after which the repeat of a(k)
  occurs) do not seem to be in the OEIS and could accompany A133788.
