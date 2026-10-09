# gapchain

A threaded search for **runs of consecutive primes whose gaps increase by two
each step**: 2, 4, 6, 8, ... ([A016045](https://oeis.org/A016045)), and for
the mirror image, gaps that **decrease** by two down to 2: ..., 8, 6, 4, 2
([A263049](https://oeis.org/A263049)). The increasing search is the default;
`--gaps dec` or `--gaps both` selects the other.

A prime `p` heads a chain of length `L` when the next `L-1` primes — with no
primes skipped in between — sit at exactly these offsets:

```
p, p+2, p+6, p+12, p+20, p+30, ...
        the m-th prime after p is at p + m(m+1)
```

Equivalently, the gap between the m-th and (m+1)-th member of the chain is
`2m`. The smallest interesting example is a chain of five:

```
347, 349, 353, 359, 367
   +2   +4   +6   +8
```

The search reports every such chain in a range, and finishes with a histogram
of how many chains of each length it found.

## Relation to OEIS A016045

[**A016045**](https://oeis.org/A016045) records the *first* prime that heads
each chain length:

> a(n) is the smallest prime prime(k) such that the gaps between the primes
> prime(k), prime(k+1), prime(k+2), ..., prime(k+n) are 2, 4, 6, ... 2n.

Mind the indexing difference. OEIS counts **gaps**; this program counts
**primes**. A chain of `L` primes has `L-1` gaps, so:

```
gapchain --length  L  =  n + 1
OEIS index         n  =  L - 1
```

| OEIS n | `--length` | a(n) — first prime heading the chain |
| -----: | ---------: | -----------------------------------: |
|      1 |          2 |                                    3 |
|      2 |          3 |                                    5 |
|      3 |          4 |                                   17 |
|      4 |          5 |                                  347 |
|      5 |          6 |                               13,901 |
|      6 |          7 |                              128,981 |
|      7 |          8 |                              128,981 |
|      8 |          9 |                          113,575,727 |
|      9 |         10 |                        2,426,256,797 |
|     10 |         11 |                      137,168,442,221 |
|     11 |         12 |                    4,656,625,081,181 |
|     12 |         13 |                  101,951,758,179,851 |
|     13 |         14 |                  484,511,389,338,941 |
|     14 |         15 |              221,860,944,705,726,407 |
|     15 |         16 | **1,397,398,433,200,922,807** (new) |
|     16 |         17 | **254,345,176,718,302,423,991** (new) |
|     17 |         18 | **307,223,174,679,659,356,577** (new) |

Until September 2026 the OEIS entry listed terms through a(14). **a(15) =
1397398433200922807**, the
first chain of 16 consecutive primes (gaps 2 through 30), was found with
`gapsieve`; see [*New results*](#new-results). A chain meeting gaps 2..2n also
meets 2..2(n-1), so the sequence is non-decreasing. **a(16) =
254345176718302423991**, the first chain of 17 primes, and **a(17) =
307223174679659356577**, the first chain of 18, were both found and confirmed
on 2026-09-27; see [*New results*](#new-results),
[*The a(16) searches*](#the-a16-searches) and
[*Search status*](#search-status). All three are now in the entry; see
[*In the OEIS*](#in-the-oeis).

`--first` computes a single term directly — it is exactly the "smallest prime
such that..." question the sequence asks. For `a(8)`, which needs a chain of
`8 + 1 = 9` primes:

```bash
./gapchain_ps --length 9 --end 2e8 --first
```

```
 9: 113575727 113575729 113575733 113575739 113575747 113575757 113575769 113575783 113575799
```

To get the whole head of the sequence from one pass, sort the output first:
threads emit chains out of order, and a chain reported at length 8 is also the
answer for length 7, so each chain has to fill in every shorter length it
covers.

```bash
./gapchain_ps -l 2 -e 3e9 2>/dev/null | sort -k2 -n |
  awk '{for (L = 2; L <= $1 + 0; L++) if (!(L in a)) {a[L] = $2; print L, $2}}' |
  sort -n
```

```
2 3
3 5
4 17
5 347
6 13901
7 128981
8 128981
9 113575727
10 2426256797
```

The sequence is due to Robert G. Wilson v, corrected and extended by Jud
McCranie, with `a(11)`–`a(14)` contributed by Dmitry Petukhov and
`a(15)`–`a(17)`, from these searches, by Jeff Sponaugle. The entry cites
Steven Kahan, "2-4-6-8...Prime Gaps That Appreciate", *Journal of Recreational
Mathematics* **25**(1), 44–46, 1993.

## Decreasing gaps: OEIS A263049

[**A263049**](https://oeis.org/A263049) is the mirrored problem:

> a(n) = smallest prime prime(k) such that the gaps between the primes
> prime(k), prime(k+1), prime(k+2), ..., prime(k+n) are 2n, 2n-2, ... 6, 4, 2.

The example given in the entry is 1979, 1987, 1993, 1997, 1999, with gaps 8,
6, 4, 2, which makes a(4) = 1979. A decreasing chain of `L` primes starting
at `p` has members `p + m(2L-1-m)`. The same indexing rule applies:
`--length L` finds `a(L-1)`.

| OEIS n | `--length` | a(n) |
| -----: | ---------: | ---: |
|      1 |          2 |                           3 |
|      2 |          3 |                           7 |
|      3 |          4 |                          31 |
|      4 |          5 |                       1,979 |
|      5 |          6 |                      41,203 |
|      6 |          7 |                     752,251 |
|      7 |          8 |                   5,647,457 |
|      8 |          9 |                  32,465,047 |
|      9 |         10 |                 245,333,233 |
|     10 |         11 |                 245,333,213 |
|     11 |         12 |              27,797,667,517 |
|     12 |         13 |         196,559,847,120,517 |
|     13 |         14 |       3,040,284,075,731,561 |
|     14 |         15 |     253,253,149,671,986,983 |
|     15 |         16 | **253,253,149,671,986,953** (new) |
|     16 |         17 | **134,491,669,675,212,201,371** (new) |

The sequence is due to Vasily Danilov, with a(12)–a(14) from Dmitry
Petukhov, and a(15) and a(16), from these searches, from Jeff Sponaugle.

**a(15) = 253253149671986953** was found with `gapsieve`; see
[*New results*](#new-results). a(16), the first 17-prime decreasing chain, is
at least a(15) − 32, and the search from there found it: **a(16) =
134491669675212201371**, on 2026-09-27; see [*New results*](#new-results).
Both are now in the entry; see [*In the OEIS*](#in-the-oeis). The next open
term is a(17), the first 18-prime decreasing chain.

A decreasing chain differs from an increasing one in two ways.

- **It grows at the front.** Its last gap is always 2, so a longer chain
  contains the shorter pattern as its tail, not its head. That is why a(10)
  is *smaller* than a(9). The 11-prime chain from 245,333,213 has the 10-prime
  pattern starting at its second prime, 245,333,233, and that is a(9). This
  sequence is not monotone. It can dip by at most one gap: a(n) ≥ a(n-1) -
  2n.
- **The term is the head of the exact pattern.** So `--first --gaps dec`
  ranks candidates by where the exact L-prime pattern starts, which is what
  the OEIS entry records. When that pattern ends a longer chain, the output
  says so.

`gapsieve --first --gaps dec` reproduces every known term, a(1) through a(14).
Details are in [README_SIEVE.md](README_SIEVE.md).

## New results

Two new terms were found with these tools on 2026-09-23 and three more on
2026-09-27. Before then, both entries listed terms only through a(14);
A016045 had last been revised in 2021 and A263049 in 2016. All five are now
in the OEIS (see [*In the OEIS*](#in-the-oeis), below).

| Sequence | Term | Value |
| -------- | ---- | ----- |
| A016045 (increasing gaps) | a(15) | 1397398433200922807 |
| A263049 (decreasing gaps) | a(15) | 253253149671986953 |
| A263049 (decreasing gaps) | a(16) | 134491669675212201371 (2026-09-27) |
| A016045 (increasing gaps) | a(16) | 254345176718302423991 (2026-09-27) |
| A016045 (increasing gaps) | a(17) | 307223174679659356577 (2026-09-27) |

### In the OEIS

Jeff Sponaugle submitted all five terms, and the OEIS editors approved them.
Each entry's extensions line credits the new terms to Jeff Sponaugle. The revision
numbers and dates below come from the entries' histories
([A016045](https://oeis.org/history?seq=A016045),
[A263049](https://oeis.org/history?seq=A263049)), in US Eastern time.

| Entry | Added | Submitted | Approved |
| --- | --- | --- | --- |
| A016045 | a(15) | revision 32, 2026-09-24 | revision 39, 2026-09-24 |
| A263049 | a(15) | revision 37, 2026-09-24 | revision 39, 2026-09-24 |
| A263049 | a(16), a comment on the gaps around its chain, and the keyword `hard` | revision 40, 2026-09-27 | revision 50, 2026-09-28 |
| A016045 | a(16) and a(17) | revision 40, 2026-09-27 | revision 43, 2026-09-27 |

Along the way, Michel Marcus renamed both entries to write prime(k) for the
kth prime instead of p(k). A263049's old b-file, which covered only a(1) to
a(14), was removed, so both entries' b-files are now generated from their
data. The A263049 comment gives the gap before a(16) as 10 and the gap after
its chain as 48. Miller-Rabin and Baillie-PSW, as in `verify_chain.py`, agree
on both.

Further OEIS submissions these searches support are drafted, ready to paste,
in [submission.md](submission.md). They include new terms for A349121 and
A090870, two new sequences, and comments, programs and cross-references.

### A016045 a(15) = 1397398433200922807

Sixteen consecutive primes with gaps 2, 4, 6, ..., 30:

```
1397398433200922807 1397398433200922809 1397398433200922813 1397398433200922819
1397398433200922827 1397398433200922837 1397398433200922849 1397398433200922863
1397398433200922879 1397398433200922897 1397398433200922917 1397398433200922939
1397398433200922963 1397398433200922989 1397398433200923017 1397398433200923047
```

The next prime is 1397398433200923107, 60 further on rather than 32, so the
chain is exactly 16 primes long and a(16) is a different, larger chain. a(15)
is about 6.3 times a(14), far below the roughly 8×10¹⁸ that the sequence's
earlier growth suggested.

**How it was found.** The search was a single `--first` run from a(14):

```bash
./gapsieve -l 16 -s 221860944705726407 -e 1.8e19 --first -c a15.ck -r --log a15.log
```

Starting at a(14) is enough because the sequence is non-decreasing. The run
began on a 14-core Apple Silicon laptop. After about 0.1% of the range it moved
through its checkpoint to house1, a 32-thread x86-64 Linux machine, built with
GCC and `-march=native` and running at about 112 T/s. The hit came after 6.61%
of the range up to 1.8×10¹⁹; house1's final session took 2h49m.

**How it was checked.**

- **The chain.** primesieve 12.14's exact sieve lists these sixteen numbers as
  consecutive primes with the gaps above. A separate Python Miller-Rabin
  (deterministic below 3.3×10²⁴) finds exactly these sixteen primes between
  a(15) and a(15) + 240, and none in between. `gapchain_ps`, enumerating the
  real primes with primesieve, reports the same chain.
- **Minimality.** The whole range from a(14) up to the candidate was searched
  again on the Apple Silicon laptop, which differs from house1 in CPU
  architecture, compiler (clang) and build. It took 2h09m of search time, and
  the only chain found was this one. Nothing below a(14) can qualify, because
  the sequence is non-decreasing.

### A016045 a(16) = 254345176718302423991

Seventeen consecutive primes with gaps 2, 4, 6, ..., 32:

```
254345176718302423991 254345176718302423993 254345176718302423997
254345176718302424003 254345176718302424011 254345176718302424021
254345176718302424033 254345176718302424047 254345176718302424063
254345176718302424081 254345176718302424101 254345176718302424123
254345176718302424147 254345176718302424173 254345176718302424201
254345176718302424231 254345176718302424263
```

The next prime is 254345176718302424299, 36 further on rather than 34, so the
chain is exactly 17 primes long and a(17) is a different, larger chain. a(16)
is about 182 times a(15).

**How it was found.** The laptop found it at 16:54 on 2026-09-27, 18 minutes
into its part of the last open range below the 18-prime chain (next section),
on its GPU plus 12 CPU threads (a16w80):

```bash
./gapsieve -l 17 -s 252971263259530362880 -e 255239263259530362880 --first -c a16w80.ck -r --log a16w80.log --gpu
```

**How it was checked.**

- **The chain.** `gapsieve` proves each member prime with Miller-Rabin on the
  13 smallest prime bases, which is exact below 3.3×10²⁴, and checks that
  every number between the members is composite. `verify_chain.py`, an
  independent Python check, tested all 273 numbers from the head to the last
  member with the same Miller-Rabin and separately with Baillie-PSW. Both
  found exactly these seventeen primes, and the next prime 36 past the chain.
- **Minimality.** Every start below 2×10²⁰ had already been searched with no
  chain. From 2×10²⁰ to the find, four runs cover every start, each checked
  against its checkpoint: the Studio's a16w18 (2×10²⁰ to its split), atom1's
  a16w66 (from there to 2.329×10²⁰), the laptop's a16w19 (to
  252971263259530362880, where it stopped at 15:38) and the laptop's a16w80.
  In `--first` mode a run finishes every segment below the one holding its
  chain before it stops. a16w80's checkpoint shows 1,080,692 segments done,
  and the chain is in segment 1,080,666. A CPU-only rerun of that segment and
  the ten before it, through a separate code path from the GPU kernels,
  reports the same chain first.

### A016045 a(17) = 307223174679659356577

Eighteen consecutive primes with gaps 2, 4, 6, ..., 34:

```
307223174679659356577 307223174679659356579 307223174679659356583
307223174679659356589 307223174679659356597 307223174679659356607
307223174679659356619 307223174679659356633 307223174679659356649
307223174679659356667 307223174679659356687 307223174679659356709
307223174679659356733 307223174679659356759 307223174679659356787
307223174679659356817 307223174679659356849 307223174679659356883
```

All eighteen are the chain a(17) asks for (gaps 2 to 34). Its first seventeen
are a 17-prime chain too, so when it was found it was the candidate for both
a(16) and a(17). At 16:54 the laptop found the smaller 17-prime chain above,
which is a(16). This one remains a(17)'s candidate: only an 18-prime chain
below it could replace it.

**How it was found.** atom1, one of the two DGX Sparks, found it at 14:39 on
2026-09-27 with the new CUDA build, 83 minutes into its first piece, 3×10²⁰
to 3.594×10²⁰ (a16w40):

```bash
./gapsieve -l 17 -s 3e20 -e 3.594e20 --first -c a16w40.ck -r --log a16w40.log --gpu -t 0
```

**How it was checked.** `gapsieve` proves each member prime with Miller-Rabin
on the 13 smallest prime bases, which is exact below 3.3×10²⁴, and checks
that every number between the members is composite. An independent Python
check tested all 579 numbers from the prime before the head to six primes
past the chain, with the same Miller-Rabin and separately with Baillie-PSW;
both found exactly these primes. The gap after the chain is 30, not 36, so
it has exactly 18 primes. The prime before the head is 24 below it, but an
increasing chain can only grow at its end.

**Why it is the smallest.** Every start below 2×10²⁰ had been searched, and
atom1 searched from 3×10²⁰ up to the chain. All of 2×10²⁰ to 3×10²⁰ was then
searched in 27 pieces on thirteen machines, the last finishing at 17:08 on
2026-09-27. Their checkpoints show that each piece starts at or below where
the one before it reached, so no start was skipped. The only chain in that
range is a(16)'s, and it has exactly 17 primes. Every 18-prime chain begins
with a 17-prime one, so no 18-prime chain starts below this one: it is a(17),
confirmed at 17:08.

### A263049 a(15) = 253253149671986953

Sixteen consecutive primes with gaps 30, 28, 26, ..., 2:

```
253253149671986953 253253149671986983 253253149671987011 253253149671987037
253253149671987061 253253149671987083 253253149671987103 253253149671987121
253253149671987137 253253149671987151 253253149671987163 253253149671987173
253253149671987181 253253149671987187 253253149671987191 253253149671987193
```

**How it was found.** The new decreasing-gap search was being validated by
re-deriving the known a(14) from scratch on the laptop:

```bash
./gapsieve -l 15 -s 3040284075731531 -e 2.6e17 --first --gaps dec
```

In 32.5 minutes at about 128 T/s it returned exactly a(14) =
253253149671986983. It also reported that this 15-prime pattern is the tail of
a 16-prime chain starting 30 lower. The known a(14) chain had been part of a
longer one all along.

**Why it is the smallest.** A 16-prime decreasing chain contains a 15-prime
pattern starting at its second prime, 30 above its head. So
a(15) ≥ a(14) − 30, which is exactly where this chain starts. No further search
is needed.

**How it was checked.** primesieve's exact sieve, the Python Miller-Rabin,
`gapchain_ps` and `gapsieve` all agree. The prime before the head is
253253149671986951, only 2 below, so the chain cannot extend further back, and
a(16) is a different chain.

### A263049 a(16) = 134491669675212201371

Seventeen consecutive primes with gaps 32, 30, 28, ..., 2:

```
134491669675212201371 134491669675212201403 134491669675212201433
134491669675212201461 134491669675212201487 134491669675212201511
134491669675212201533 134491669675212201553 134491669675212201571
134491669675212201587 134491669675212201601 134491669675212201613
134491669675212201623 134491669675212201631 134491669675212201637
134491669675212201641 134491669675212201643
```

**How it was found.** The laptop's GPU found it at 00:49 on 2026-09-27, in
its piece of the search from 1.07×10²⁰ to 1.85×10²⁰ (d16w14), after 35.25%
of the piece:

```bash
nice -n 19 ./gapsieve -l 17 -s 1.07e20 -e 1.85e20 --first --gaps dec -c d16w14.ck -r --log d16w14.log --gpu
```

**How it was checked.** `gapsieve` proves each member prime with Miller-Rabin
on the 13 smallest prime bases, which is exact below 3.3×10²⁴, and checks
that every number between the members is composite. An independent Python
check (`verify_chain.py`) tested all 273 numbers from the head to the last
member, with the same Miller-Rabin and separately with Baillie-PSW. Both
found exactly these 17 primes. The prime before the head is
134491669675212201361, 10 below it, so the chain cannot extend further back.

**Why it is the smallest.** Every chain start below it has been searched,
with no other chain found:
- below 2⁶⁴, on the Studio and the laptop (see
  [*The a(16) searches*](#the-a16-searches)), with the last million checked
  in Python;
- 2⁶⁴ to 6.6×10¹⁹, in pieces on nine machines, finished at 18:59 on
  2026-09-26;
- 6.6×10¹⁹ to 1.07×10²⁰, by the M1 Max, the Studio and the laptop, finished at
  09:47 on 2026-09-27, when this became a(16);
- 1.07×10²⁰ up to the chain, by the laptop run that found it, which searched
  every segment below the chain's before stopping.

The pieces meet at their boundaries, each checked against its checkpoint.

### The a(16) searches

Both next terms need chains of 17 primes. Neither search found one starting
anywhere below 2⁶⁴, so both a(16) terms are above 2⁶⁴ ≈ 1.8447×10¹⁹. Both
searches then continued above 2⁶⁴ (below), and by the evening of 2026-09-26
**A016045's a(16) was above 10²⁰ and A263049's above 6.6×10¹⁹**. At 00:49
the next morning the laptop found A263049's a(16), at
134,491,669,675,212,201,371, confirmed at 09:47 once every start below it had
been searched (see above). At 16:54 that day the laptop found A016045's
a(16), at 254,345,176,718,302,423,991, confirmed the same minute, since every
start below it had already been searched.
[*Search status*](#search-status) has the latest bounds. The search up to 2⁶⁴
was done in pieces:

| Term | Chain starts searched | Machine | Time |
| ---- | --------------------- | ------- | ---- |
| A263049 a(16) | a(15) − 32 to 1.8×10¹⁹ | Mac Studio, 24 threads, ~205 T/s | 24h04m, finished 2026-09-24 13:55 |
| | 1.8×10¹⁹ to 1.844×10¹⁹ | Mac Studio | 36m |
| | 1.844×10¹⁹ to 2⁶⁴ − 273 | the 14-core Apple Silicon laptop | 44s, 2026-09-25 08:17 |
| | the last million below 2⁶⁴ | a separate Python check (below) | seconds |
| A016045 a(16) | a(15) to 1.35×10¹⁹ | the laptop, ~164 T/s | 26h16m, finished 2026-09-25 00:38 |
| | 1.35×10¹⁹ to 1.844×10¹⁹ | Mac Studio | 6h42m |
| | 1.844×10¹⁹ to 2⁶⁴ − 9,507 | the laptop | 43s, 2026-09-25 08:16 |
| | the last million below 2⁶⁴ | the Python check | seconds |

Running past its share, the laptop also re-searched Studio's A016045 part from
1.35×10¹⁹ to 1.732×10¹⁹, and found nothing there either.

The commands, in table order:

```bash
./gapsieve -l 17 -s 253253149671986921 -e 1.8e19 --first --gaps dec -c d16.ck -r --log d16.log
```

```bash
./gapsieve -l 17 -s 18000000000000000000 -e 1.844e19 --first --gaps dec -c d16b.ck -r --log d16.log
```

```bash
./gapsieve -l 17 -s 18440000000000000000 -e 18446744073709551343 --first --gaps dec -c d16top.ck -r --log d16top.log
```

```bash
./gapsieve -l 17 -s 1397398433200922807 -e 1.8e19 --first -c a16.ck -r --log a16.log
```

```bash
./gapsieve -l 17 -s 1.35e19 -e 1.844e19 --first -c a16hi.ck -r --log a16hi.log
```

```bash
./gapsieve -l 17 -s 18440000000000000000 -e 18446744073709542109 --first -c a16top.ck -r --log a16top.log
```

**Reaching 2⁶⁴.** `gapsieve` used to stop at 1.844×10¹⁹. For the last
stretch it was extended to search as high as 64-bit arithmetic allows, which
is where the two ends in the last commands above come from (they were given
as `-e max` to that 64-bit build). Every number it checked for a chain start
had to stay below 2⁶⁴, so it stopped 273 short for decreasing chains and
9,507 short for increasing ones, which it follows up to 80 primes past the 17
to report their full length. That build was checked in the last 10⁹ below 2⁶⁴ against a
brute force over primesieve's exact list of primes. All 311 chains of length
5 of both shapes and all 2,038 increasing chains of length 4 matched, and a
decreasing-only run matched right up to 2⁶⁴ − 21. Below the old ceiling it
gives exactly the same output as before.

The starts it cannot reach, whose chains would cross 2⁶⁴, were covered by a
short Python program, [`tail_check.py`](tail_check.py). It lists every prime
from 2⁶⁴ − 1,001,000 to 2⁶⁴ + 1,000 with a Miller–Rabin test that is
deterministic far past 2⁶⁴. It agrees with primesieve on all 22,495 of those
primes below 2⁶⁴, and it finds the known largest 64-bit prime, 2⁶⁴ − 59.
Among starts in the last million below 2⁶⁴, the longest increasing chain has 3
primes, and the longest start of the 17-prime decreasing pattern has 4. To
rerun it, which takes a few seconds and needs the `primesieve` command:

```bash
python3 tail_check.py
```

A search that covers a range without a hit shows that the term exceeds the
bound, not that it doesn't exist.

### Continuing above 2⁶⁴

`gapsieve` now works with 128-bit numbers, up to 3.3×10²⁴, so both searches
can continue. Numbers above 2⁶⁴ cost more to test: at `-l 17` the laptop
runs about 115 T/s there instead of 160 ([README_SIEVE.md](README_SIEVE.md)
has the details and the checks behind the 128-bit version). Starting each
search at 1.8446744×10¹⁹ overlaps the finished range by 7×10¹⁰, which costs
under a second:

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 1e20 --first -c a16w.ck -r --log a16w.log
```

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 1e20 --first --gaps dec -c d16w.ck -r --log d16w.log
```

For a rough idea of the odds, take the chain counts measured at 1.2×10¹⁹
(length 12 was the longest with any) and extend them to length 17. That puts
each a(16) below 10²⁰ with about a 25% chance, below 2×10²⁰ with 45%, below
5×10²⁰ with 75% and below 10²¹ with 90%. Treat those as good to within a
factor of 2–3 in range.

**Checking a new machine.** Every machine is checked before it joins. This
run lists all 234 chains of length 7 in a window above 2⁶⁴, and the digest of
its output must be `71fe408bdec8f5c6` (on Linux, use `sha256sum`):

```bash
./gapsieve -l 7 -g both -s 2e20 -e 200000000200000000000 -q | LC_ALL=C sort | shasum -a 256 | cut -c1-16
```

The A016045 search above 2⁶⁴ runs in pieces, one per machine. On
2026-09-26, once the Mac Studio searched on its GPU at about 1,280 T/s, the
pieces were rebalanced so that all of them finish together, at about 19:40
that evening instead of about 03:00 the next morning. They did: every start
below 10²⁰ had been searched by 19:39, with no chain found. A watcher on each of
mathb, matha and mathd (`split_a16w1.sh`, `split_a16w2.sh`, `split_a16w5.sh`)
stops its run at a split point, a segment boundary. The Studio searches the
rest of each piece after its own, in turn (`studio_queue.sh`).

| Range | Machine | Status |
| ----- | ------- | ------ |
| 1.8446744×10¹⁹ to 3.99614×10¹⁹ | the laptop for the first 1.18%, then mathb (2 × 24-core Xeon Platinum 8268, 96 threads, ~194 T/s) | started 2026-09-25 12:10 as 1.8446744×10¹⁹ to 4.5×10¹⁹; moved to mathb at 13:14. Its watcher stopped it there at 19:38, as planned |
| 4.5×10¹⁹ to 5.53934×10¹⁹ | matha (2 × 12-core Xeon Platinum 8158, 48 threads, ~101 T/s) | started 2026-09-25 15:05 as 4.5×10¹⁹ to 6.5×10¹⁹. Its watcher stopped it at 55,393,432,389,038,899,200 (8,175,065 segments) at 19:36, as planned; the first split was at 5.75×10¹⁹ |
| 5.75×10¹⁹ to 6.06276×10¹⁹ | mathd (~186 T/s) | started 2026-09-26 14:56 as 5.75×10¹⁹ to 6.5×10¹⁹, by `queue.sh` when mathd's first range finished. Its watcher stopped it at 60,627,646,902,545,285,120 (2,460,084 segments) at 19:39, as planned |
| 6.5×10¹⁹ to 8.45×10¹⁹ | Mac Studio (M2 Ultra) | started 2026-09-25 12:15 as 6.5×10¹⁹ to 10²⁰ on the CPU (~142 T/s); a watcher stops it once it passes 8.45×10¹⁹. Paused 2026-09-26 11:34 at 34.09% to test and tune the GPU search; resumed 15:17 on the GPU (`--gpu`, ~914 T/s) with a new watcher (`handoff2.sh`), and restarted 15:56 at 40.14% with the faster GPU build (~1,280 T/s). Passed 8.45×10¹⁹ at 17:08 |
| 8.45×10¹⁹ to 10²⁰ | mathd (2 × 24-core Xeon Platinum 8168, 96 threads) | started 2026-09-25 15:23; **finished 2026-09-26 14:55, no chain** |
| 3.99614×10¹⁹ to 4.5×10¹⁹ | Mac Studio (a16w6) | finished 2026-09-26, no chain |
| 5.53934×10¹⁹ to 5.75×10¹⁹ | Mac Studio (a16w7) | finished 2026-09-26, no chain |
| 6.06276×10¹⁹ to 6.5×10¹⁹ | Mac Studio (a16w8) | finished 2026-09-26 19:37, no chain |

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 4.5e19 --first -c a16w1.ck -r --log a16w1.log
```

```bash
./gapsieve -l 17 -s 4.5e19 -e 6.5e19 --first -c a16w2.ck -r --log a16w2.log
```

```bash
./gapsieve -l 17 -s 5.75e19 -e 6.5e19 --first -c a16w5.ck -r --log a16w5.log
```

```bash
./gapsieve -l 17 -s 6.5e19 -e 1e20 --first -c a16w3.ck -r --log a16w3.log --gpu
```

```bash
./gapsieve -l 17 -s 8.45e19 -e 1e20 --first -c a16w4.ck -r --log a16w4.log
```

The Studio's queued pieces, run one after another:

```bash
./gapsieve -l 17 -s 39961410245763891200 -e 4.5e19 --first -c a16w6.ck -r --log a16w6.log --gpu
```

```bash
./gapsieve -l 17 -s 55393432389038899200 -e 5.75e19 --first -c a16w7.ck -r --log a16w7.log --gpu
```

```bash
./gapsieve -l 17 -s 60627646902545285120 -e 6.5e19 --first -c a16w8.ck -r --log a16w8.log --gpu
```

**Above 10²⁰.** When each machine's work to 10²⁰ ends, it starts a piece of
10²⁰ to 2×10²⁰ (`next_a16w1.sh`, `next_a16w2.sh`, `next_a16w5.sh` and
`next_studio.sh`). The pieces are sized to each machine's speed, so that all
should finish at about 11:00 on 2026-09-27. The fastest machine takes the
lowest range, so the proven bound rises fastest, and a chain found there is
confirmed soonest.

| Range | Machine | Status |
| ----- | ------- | ------ |
| 10²⁰ to 1.73×10²⁰ | Mac Studio, GPU (a16w9) | started 2026-09-26 19:37 |
| 1.73×10²⁰ to 1.84×10²⁰ | mathb (a16w10) | started 2026-09-26 19:38 |
| 1.84×10²⁰ to 1.945×10²⁰ | mathd (a16w11) | started 2026-09-26 19:39 |
| 1.945×10²⁰ to 2×10²⁰ | matha (a16w12) | started 2026-09-26 19:36 |

```bash
./gapsieve -l 17 -s 1e20 -e 1.73e20 --first -c a16w9.ck -r --log a16w9.log --gpu
```

```bash
./gapsieve -l 17 -s 1.73e20 -e 1.84e20 --first -c a16w10.ck -r --log a16w10.log
```

```bash
./gapsieve -l 17 -s 1.84e20 -e 1.945e20 --first -c a16w11.ck -r --log a16w11.log
```

```bash
./gapsieve -l 17 -s 1.945e20 -e 2e20 --first -c a16w12.ck -r --log a16w12.log
```

**2×10²⁰ to 3×10²⁰, with the whole fleet.** When A263049's machines were
freed on 2026-09-27, A016045 was re-planned around all eleven. First, the
Studio resumes its piece a16w9 after its help with A263049 (`studio_help.sh`),
and a watcher (`split_a16w9.sh`) stops it at 166,246,697,682,959,073,280
(52,107,046 segments). The laptop searches the rest of that piece (a16w17)
once the range below A263049's find is done (`laptop_queue3.sh`). Then 2×10²⁰
to 3×10²⁰ is divided among all eleven by speed and by when each is free, the
fastest at the bottom, all to finish at about 18:20 on 2026-09-27. (The
earlier four-machine plan for this range, with pieces a16w13 to a16w16,
never started.) The scripts are `next3_studio.sh`, `laptop_queue3.sh`,
`next3_m1max.sh`, `next3_a16w10.sh`, `next3_a16w11.sh`, `next3_a16w12.sh`, and
`switch_to_a16.sh` for the five machines that moved over at once.

| Range | Machine | Status |
| ----- | ------- | ------ |
| 1.66247×10²⁰ to 1.73×10²⁰ | the laptop, GPU (a16w17) | started 2026-09-27 09:46; finished 11:14, no chain |
| 2×10²⁰ to 2.329×10²⁰ | Mac Studio, GPU (a16w18) | started 2026-09-27 11:17, when a16w9 reached its split |
| 2.329×10²⁰ to 2.656×10²⁰ | the laptop, GPU (a16w19) | started 2026-09-27 11:14 |
| 2.656×10²⁰ to 2.791×10²⁰ | the M1 Max, GPU (a16w20) | started 2026-09-27 09:47 |
| 2.791×10²⁰ to 2.838×10²⁰ | mathb (a16w21) | queued, after a16w10 |
| 2.838×10²⁰ to 2.882×10²⁰ | mathd (a16w22) | queued, after a16w11 |
| 2.882×10²⁰ to 2.92×10²⁰ | mathc (a16w23) | started 2026-09-27 07:57 |
| 2.92×10²⁰ to 2.948×10²⁰ | matha (a16w24) | started 2026-09-27 10:46 |
| 2.948×10²⁰ to 2.969×10²⁰ | the M1 Pro (a16w25) | started 2026-09-27 07:57 |
| 2.969×10²⁰ to 2.981×10²⁰ | the WSL laptop (a16w26) | started 2026-09-27 07:57. Its power adapter came unplugged and it hibernated on battery, from about 10:55 to 11:20 and again from 12:05 to 12:48; the run carried on each time it woke |
| 2.981×10²⁰ to 2.9925×10²⁰ | math4 (a16w27) | started 2026-09-27 07:57 |
| 2.9925×10²⁰ to 3×10²⁰ | the Intel Mac (a16w28) | started 2026-09-27 07:59 |

```bash
nice -n 19 ./gapsieve -l 17 -s 166246697682959073280 -e 1.73e20 --first -c a16w17.ck -r --log a16w17.log --gpu
```

```bash
./gapsieve -l 17 -s 2e20 -e 2.329e20 --first -c a16w18.ck -r --log a16w18.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 2.329e20 -e 2.656e20 --first -c a16w19.ck -r --log a16w19.log --gpu
```

```bash
./gapsieve -l 17 -s 2.656e20 -e 2.791e20 --first -c a16w20.ck -r --log a16w20.log --gpu
```

```bash
./gapsieve -l 17 -s 2.791e20 -e 2.838e20 --first -c a16w21.ck -r --log a16w21.log
```

```bash
./gapsieve -l 17 -s 2.838e20 -e 2.882e20 --first -c a16w22.ck -r --log a16w22.log
```

```bash
./gapsieve -l 17 -s 2.882e20 -e 2.92e20 --first -c a16w23.ck -r --log a16w23.log
```

```bash
./gapsieve -l 17 -s 2.92e20 -e 2.948e20 --first -c a16w24.ck -r --log a16w24.log
```

```bash
nice -n 19 ./gapsieve -l 17 -s 2.948e20 -e 2.969e20 --first -c a16w25.ck -r --log a16w25.log -t 8
```

```bash
./gapsieve -l 17 -s 2.969e20 -e 2.981e20 --first -c a16w26.ck -r --log a16w26.log
```

```bash
./gapsieve -l 17 -s 2.981e20 -e 2.9925e20 --first -c a16w27.ck -r --log a16w27.log
```

```bash
./gapsieve -l 17 -s 2.9925e20 -e 3e20 --first -c a16w28.ck -r --log a16w28.log
```

**3×10²⁰ to 5×10²⁰, and 5×10²⁰ to 7×10²⁰, with thirteen machines.** On
2026-09-27 at 13:16 two DGX Sparks joined: atom1 and atom2, each with an
NVIDIA GB10 running the CUDA build (`make CUDA=1`; see
[README_SIEVE.md](README_SIEVE.md)) at about 1,400 T/s. So 3×10²⁰ to 5×10²⁰
(planned at 11:25 for eleven machines, to finish at 09:25 on 2026-09-28) was
re-planned for thirteen. The Sparks start at once, at the bottom, and the
other eleven start their pieces as their pieces of the range below finish,
between about 18:10 and 19:45. All should finish at about 00:45 on
2026-09-28. Each machine then goes straight on to a piece of 5×10²⁰ to
7×10²⁰, split by speed, all to finish at about 09:20. One script per machine
(`queue5_*.sh`) runs its pieces in turn. It starts nothing after a piece that
didn't finish cleanly or reported a candidate chain. On a Mac the whole queue
runs under `caffeinate`.

| Range | Machine | Status |
| ----- | ------- | ------ |
| 3×10²⁰ to 3.594×10²⁰ | atom1, DGX Spark, GPU (a16w40) | started 2026-09-27 13:16; **found an 18-prime chain at 307,223,174,679,659,356,577 at 14:39**, after 12.16% of the piece |
| 3.594×10²⁰ to 4.15×10²⁰ | atom2, DGX Spark, GPU (a16w41) | started 2026-09-27 13:16; stopped at 14:50, at 13.94%, above the find |
| 4.15×10²⁰ to 4.443×10²⁰ | Mac Studio, GPU (a16w42) | queued, after a16w18 |
| 4.443×10²⁰ to 4.733×10²⁰ | the laptop, GPU (a16w43) | queued, after a16w19 |
| 4.733×10²⁰ to 4.835×10²⁰ | the M1 Max, GPU (a16w44) | queued, after a16w20 |
| 4.835×10²⁰ to 4.88×10²⁰ | mathb (a16w45) | queued, after a16w21 |
| 4.88×10²⁰ to 4.922×10²⁰ | mathd (a16w46) | queued, after a16w22 |
| 4.922×10²⁰ to 4.947×10²⁰ | mathc (a16w47) | queued, after a16w23 |
| 4.947×10²⁰ to 4.97×10²⁰ | matha (a16w48) | queued, after a16w24 |
| 4.97×10²⁰ to 4.983×10²⁰ | the M1 Pro (a16w49) | queued, after a16w25 |
| 4.983×10²⁰ to 4.989×10²⁰ | the WSL laptop (a16w50) | queued, after a16w26 |
| 4.989×10²⁰ to 4.996×10²⁰ | math4 (a16w51) | queued, after a16w27 |
| 4.996×10²⁰ to 5×10²⁰ | the Intel Mac (a16w52) | queued, after a16w28 |

| Range | Machine | Status |
| ----- | ------- | ------ |
| 5×10²⁰ to 5.442×10²⁰ | atom1, DGX Spark, GPU (a16w53) | queued, after a16w40 |
| 5.442×10²⁰ to 5.856×10²⁰ | atom2, DGX Spark, GPU (a16w54) | queued, after a16w41 |
| 5.856×10²⁰ to 6.252×10²⁰ | Mac Studio, GPU (a16w55) | queued, after a16w42 |
| 6.252×10²⁰ to 6.643×10²⁰ | the laptop, GPU (a16w56) | queued, after a16w43 |
| 6.643×10²⁰ to 6.779×10²⁰ | the M1 Max, GPU (a16w57) | queued, after a16w44 |
| 6.779×10²⁰ to 6.839×10²⁰ | mathb (a16w58) | queued, after a16w45 |
| 6.839×10²⁰ to 6.895×10²⁰ | mathd (a16w59) | queued, after a16w46 |
| 6.895×10²⁰ to 6.927×10²⁰ | mathc (a16w60) | queued, after a16w47 |
| 6.927×10²⁰ to 6.958×10²⁰ | matha (a16w61) | queued, after a16w48 |
| 6.958×10²⁰ to 6.975×10²⁰ | the M1 Pro (a16w62) | queued, after a16w49 |
| 6.975×10²⁰ to 6.985×10²⁰ | the WSL laptop (a16w63) | queued, after a16w50 |
| 6.985×10²⁰ to 6.994×10²⁰ | math4 (a16w64) | queued, after a16w51 |
| 6.994×10²⁰ to 7×10²⁰ | the Intel Mac (a16w65) | queued, after a16w52 |

```bash
./gapsieve -l 17 -s 3e20 -e 3.594e20 --first -c a16w40.ck -r --log a16w40.log --gpu -t 0
```

```bash
./gapsieve -l 17 -s 5e20 -e 5.442e20 --first -c a16w53.ck -r --log a16w53.log --gpu -t 0
```

```bash
./gapsieve -l 17 -s 3.594e20 -e 4.15e20 --first -c a16w41.ck -r --log a16w41.log --gpu -t 0
```

```bash
./gapsieve -l 17 -s 5.442e20 -e 5.856e20 --first -c a16w54.ck -r --log a16w54.log --gpu -t 0
```

```bash
./gapsieve -l 17 -s 4.15e20 -e 4.443e20 --first -c a16w42.ck -r --log a16w42.log --gpu
```

```bash
./gapsieve -l 17 -s 5.856e20 -e 6.252e20 --first -c a16w55.ck -r --log a16w55.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 4.443e20 -e 4.733e20 --first -c a16w43.ck -r --log a16w43.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 6.252e20 -e 6.643e20 --first -c a16w56.ck -r --log a16w56.log --gpu
```

```bash
./gapsieve -l 17 -s 4.733e20 -e 4.835e20 --first -c a16w44.ck -r --log a16w44.log --gpu
```

```bash
./gapsieve -l 17 -s 6.643e20 -e 6.779e20 --first -c a16w57.ck -r --log a16w57.log --gpu
```

```bash
./gapsieve -l 17 -s 4.835e20 -e 4.88e20 --first -c a16w45.ck -r --log a16w45.log
```

```bash
./gapsieve -l 17 -s 6.779e20 -e 6.839e20 --first -c a16w58.ck -r --log a16w58.log
```

```bash
./gapsieve -l 17 -s 4.88e20 -e 4.922e20 --first -c a16w46.ck -r --log a16w46.log
```

```bash
./gapsieve -l 17 -s 6.839e20 -e 6.895e20 --first -c a16w59.ck -r --log a16w59.log
```

```bash
./gapsieve -l 17 -s 4.922e20 -e 4.947e20 --first -c a16w47.ck -r --log a16w47.log
```

```bash
./gapsieve -l 17 -s 6.895e20 -e 6.927e20 --first -c a16w60.ck -r --log a16w60.log
```

```bash
./gapsieve -l 17 -s 4.947e20 -e 4.97e20 --first -c a16w48.ck -r --log a16w48.log
```

```bash
./gapsieve -l 17 -s 6.927e20 -e 6.958e20 --first -c a16w61.ck -r --log a16w61.log
```

```bash
nice -n 19 ./gapsieve -l 17 -s 4.97e20 -e 4.983e20 --first -c a16w49.ck -r --log a16w49.log -t 8
```

```bash
nice -n 19 ./gapsieve -l 17 -s 6.958e20 -e 6.975e20 --first -c a16w62.ck -r --log a16w62.log -t 8
```

```bash
./gapsieve -l 17 -s 4.983e20 -e 4.989e20 --first -c a16w50.ck -r --log a16w50.log
```

```bash
./gapsieve -l 17 -s 6.975e20 -e 6.985e20 --first -c a16w63.ck -r --log a16w63.log
```

```bash
./gapsieve -l 17 -s 4.989e20 -e 4.996e20 --first -c a16w51.ck -r --log a16w51.log
```

```bash
./gapsieve -l 17 -s 6.985e20 -e 6.994e20 --first -c a16w64.ck -r --log a16w64.log
```

```bash
./gapsieve -l 17 -s 4.996e20 -e 5e20 --first -c a16w52.ck -r --log a16w52.log
```

```bash
./gapsieve -l 17 -s 6.994e20 -e 7e20 --first -c a16w65.ck -r --log a16w65.log
```

**After the find.** At 14:39 on 2026-09-27 atom1 found an 18-prime chain at
3.072×10²⁰ (see [*New results*](#new-results)). Nothing above it can change
a(16) or a(17), so at 14:50 every queued piece above it was stopped, along
with atom2's. The one range still open below the find, 2×10²⁰ to 3×10²⁰, was
then shared among all thirteen machines, to finish at about 16:50 instead of
19:45. A watcher on each of the eleven machines searching it
(`csplit_*.sh`) stops its run at a split point. The two Sparks search the
eleven tops above those points, one after another (`cqueue_atom1.sh`,
`cqueue_atom2.sh`), and stop at once if a top reports a chain, since that
one would be smaller.

| Range | Machine | Status |
| ----- | ------- | ------ |
| 225753294219027415040 to 232900000000000000000 | atom1, top of Studio's a16w18 (a16w66) | started 2026-09-27 14:52, one after another on its Spark |
| 258536912857616220160 to 265600000000000000000 | atom2, top of laptop's a16w19 (a16w67) | started 2026-09-27 14:52, one after another on its Spark |
| 276734927651240673280 to 279100000000000000000 | atom1, top of M1 Max's a16w20 (a16w68) | started 2026-09-27 14:52, one after another on its Spark |
| 282814879427290071040 to 283800000000000000000 | atom2, top of mathb's a16w21 (a16w69) | started 2026-09-27 14:52, one after another on its Spark |
| 287231052619728814080 to 288200000000000000000 | atom2, top of mathd's a16w22 (a16w70) | started 2026-09-27 14:52, one after another on its Spark |
| 291517862366314496000 to 292000000000000000000 | atom2, top of mathc's a16w23 (a16w71) | started 2026-09-27 14:52, one after another on its Spark |
| 294217455178150051840 to 294800000000000000000 | atom1, top of matha's a16w24 (a16w72) | started 2026-09-27 14:52, one after another on its Spark |
| 296597663014124257280 to 296900000000000000000 | atom1, top of M1 Pro's a16w25 (a16w73) | started 2026-09-27 14:52, one after another on its Spark |
| 297783853005522206720 to 298100000000000000000 | atom2, top of WSL's a16w26 (a16w74) | started 2026-09-27 14:52, one after another on its Spark |
| 299065474174207262720 to 299250000000000000000 | atom1, top of math4's a16w27 (a16w75) | started 2026-09-27 14:52, one after another on its Spark |
| 299882968349079306240 to 300000000000000000000 | atom2, top of Intel's a16w28 (a16w76) | started 2026-09-27 14:52, one after another on its Spark |

At 15:38 the laptop had to leave, so its run (a16w19) was stopped there,
searched to 252,971,263,259,530,362,880. The rest of its share, up to its split
point 258,536,912,857,616,220,160, was first divided among atom1, atom2 and
the Studio (a16w77, a16w78, a16w79), which moved the confirmation to about
17:15. The laptop came back at 16:36, before any of those three had started,
so they were never run. The same range was split four ways instead, sized so
that all four finish together at about 17:07. The laptop started its part at
once (`lrun_a16w80.sh`). The other three start theirs when their own
confirmation work ends (`lhelp2_*.sh`).

| Range | Machine | Starts |
| --- | --- | --- |
| 252971263259530362880 to 255239263259530362880 | laptop (a16w80) | 2026-09-27 16:36 |
| 255239263259530362880 to 256365263259530362880 | Studio (a16w81) | when a16w18 stops at its split, about 16:53 |
| 256365263259530362880 to 257486263259530362880 | atom2 (a16w82) | when its queue of tops ends, about 16:54 |
| 257486263259530362880 to 258536912857616220160 | atom1 (a16w83) | when its queue of tops ends, about 16:55 |

```bash
./gapsieve -l 17 -s 252971263259530362880 -e 255239263259530362880 --first -c a16w80.ck -r --log a16w80.log --gpu
```

```bash
./gapsieve -l 17 -s 255239263259530362880 -e 256365263259530362880 --first -c a16w81.ck -r --log a16w81.log --gpu
```

```bash
./gapsieve -l 17 -s 256365263259530362880 -e 257486263259530362880 --first -c a16w82.ck -r --log a16w82.log --gpu -t 0
```

```bash
./gapsieve -l 17 -s 257486263259530362880 -e 258536912857616220160 --first -c a16w83.ck -r --log a16w83.log --gpu -t 0
```

At 16:54 the laptop's part (a16w80) stopped at a 17-prime chain at
254,345,176,718,302,423,991, 60.58% of the way through. That is a(16); see
[*New results*](#new-results). a(17) still needs every start below its
candidate searched for 18-prime chains, so the laptop went on from just above
the find (a16w84, `lrun_a16w84.sh`), and the other three parts continue as
they were. They also run with `--first`, so if one stops at a 17-prime chain,
the rest of its part will be searched from just above that chain.

```bash
./gapsieve -l 17 -s 254345176718302423992 -e 255239263259530362880 --first -c a16w84.ck -r --log a16w84.log --gpu
```

None stopped early. All four finished between 17:06 and 17:08 with no chain,
which confirmed a(17).

The A263049 search above 2⁶⁴ started on one machine, with a first piece sized
to end when the A016045 machines free up. On 2026-09-26 the pieces were
rebalanced around the laptop's GPU in the same way, so that all of them finish
at about 18:50 that evening instead of about 04:40 the next morning. Every
start below 6.6×10¹⁹ had been searched by 18:59, with no chain found. Watchers
(`split_d16w1.sh`, `split_d16w2.sh`, `split_d16w5.sh`, `split_d16w6.sh`,
`split_d16w7.sh`) stop five of the runs at split points, and the laptop
searches the rest of each piece after its own (`laptop_queue.sh`).

| Range | Machine | Status |
| ----- | ------- | ------ |
| 1.8446744×10¹⁹ to 2.86376×10¹⁹ | mathc (2 × 12-core Xeon Platinum 8158, 48 threads, ~105 T/s) | started 2026-09-25 15:47 as 1.8446744×10¹⁹ to 3.1×10¹⁹. Its watcher stopped it there at 18:51, as planned |
| 3.1×10¹⁹ to 3.59575×10¹⁹ | a second MacBook (M1 Pro, 10 cores, ~58 T/s) | started 2026-09-25 19:03 as 3.1×10¹⁹ to 3.8×10¹⁹; restarted 21:31 at 7.59% with `nice -n 19` and 8 threads, after running flat out starved the Mac's own services and made it unusable. Its watcher stopped it there at 18:52, as planned |
| 3.8×10¹⁹ to 4.45×10¹⁹ | an M1 Max MacBook Pro (10 cores, 32-core GPU) | started 2026-09-25 20:11 on the CPU (~59 T/s); restarted 2026-09-26 16:23 at 66.41% on its GPU (`--gpu`, ~614 T/s); finished 17:23, no chain |
| 4.45×10¹⁹ to 5.75×10¹⁹ | the 14-core M4 Max laptop | started 2026-09-25 20:22 on the CPU (~118 T/s); paused 07:21 to 09:49 on 2026-09-26 while the laptop was in use, then resumed at 35.87% with `nice -n 19`. Paused 09:56 at 36.27% to develop the GPU search; resumed 15:17 on the GPU (`--gpu`, ~640 T/s), and restarted 15:56 at 47.08% with the faster GPU build (~1,280 T/s); finished 17:27, no chain |
| 5.75×10¹⁹ to 5.99374×10¹⁹ | math4 (2 × 4-core Xeon E5-2637 v3, 16 threads, ~31 T/s) | started 2026-09-25 20:45 as 5.75×10¹⁹ to 6.05×10¹⁹. Its watcher stopped it there at 18:51, as planned |
| 6.05×10¹⁹ to 6.30844×10¹⁹ | a Windows laptop under WSL2 (Core i7-1360P, 16 threads, ~32 T/s) | started 2026-09-25 21:08 as 6.05×10¹⁹ to 6.4×10¹⁹. Its watcher stopped it there at 18:59, as planned |
| 6.4×10¹⁹ to 6.55317×10¹⁹ | a 2016 Intel MacBook Pro (Core i7-6700HQ, 8 threads, ~20 T/s; binary cross-compiled for x86-64 macOS 12) | started 2026-09-25 21:17 as 6.4×10¹⁹ to 6.6×10¹⁹. Its watcher stopped it there at 18:50, as planned |
| 2.86376×10¹⁹ to 3.1×10¹⁹ | the laptop (d16w8) | finished 2026-09-26, no chain |
| 3.59575×10¹⁹ to 3.8×10¹⁹ | the laptop (d16w9) | finished 2026-09-26, no chain |
| 5.99374×10¹⁹ to 6.05×10¹⁹ | the laptop (d16w10) | finished 2026-09-26, no chain |
| 6.30844×10¹⁹ to 6.4×10¹⁹ | the laptop (d16w11) | finished 2026-09-26, no chain |
| 6.55317×10¹⁹ to 6.6×10¹⁹ | the laptop (d16w12) | finished 2026-09-26 18:49, no chain |

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 3.1e19 --first --gaps dec -c d16w1.ck -r --log d16w1.log
```

```bash
nice -n 19 ./gapsieve -l 17 -s 3.1e19 -e 3.8e19 --first --gaps dec -c d16w2.ck -r --log d16w2.log -t 8
```

On a machine someone uses, or one running management and security software,
run gapsieve at the lowest priority and leave a couple of cores free, as on
the M1 Pro. At normal priority on every core it can starve the machine's own
services. On that Mac the backlog reached a load average of 350, the prompt
and SSH stopped responding, and it took about two minutes to drain once the
search stopped.

On 2026-09-27 the WSL laptop's power adapter came unplugged. It ran on its
battery until it hibernated, twice: from about 10:55 to 11:20, and from 12:05
to 12:48, until the adapter was reconnected. Its run was frozen in hibernation and
carried on when the laptop woke, losing only the time: a checkpoint was not
even needed. A laptop that searches overnight should be checked for mains
power, not just for its sleep settings.

On a Mac, keep it awake across a handoff. Each run holds off idle sleep with
`caffeinate -i`, but only while it runs. When the Intel Mac's watcher stopped
its run at 18:50:22, the next run's `caffeinate` started eight seconds later,
and by then the Mac, idle for hours, was already going to sleep. It slept for
an hour and a half. It woke briefly every 20 minutes or so for maintenance,
in a low-power state that limited its CPU to 40% speed, and so searched
almost nothing. Waking it (`caffeinate -u -t 2`) and holding a system-sleep
assertion for the life of the run (`caffeinate -s -w PID`) put it back to
full speed at 20:26. A chain of runs should hold one `caffeinate -s` across
all of them.

```bash
./gapsieve -l 17 -s 3.8e19 -e 4.45e19 --first --gaps dec -c d16w3.ck -r --log d16w3.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 4.45e19 -e 5.75e19 --first --gaps dec -c d16w4.ck -r --log d16w4.log --gpu
```

```bash
./gapsieve -l 17 -s 5.75e19 -e 6.05e19 --first --gaps dec -c d16w5.ck -r --log d16w5.log
```

```bash
./gapsieve -l 17 -s 6.05e19 -e 6.4e19 --first --gaps dec -c d16w6.ck -r --log d16w6.log
```

```bash
./gapsieve -l 17 -s 6.4e19 -e 6.6e19 --first --gaps dec -c d16w7.ck -r --log d16w7.log
```

The laptop's queued pieces, run one after another:

```bash
nice -n 19 ./gapsieve -l 17 -s 28637618584061050880 -e 3.1e19 --first --gaps dec -c d16w8.ck -r --log d16w8.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 35957542650153533440 -e 3.8e19 --first --gaps dec -c d16w9.ck -r --log d16w9.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 59937403886031994880 -e 6.05e19 --first --gaps dec -c d16w10.ck -r --log d16w10.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 63084449125441863680 -e 6.4e19 --first --gaps dec -c d16w11.ck -r --log d16w11.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 65531703868629975040 -e 6.6e19 --first --gaps dec -c d16w12.ck -r --log d16w12.log --gpu
```

**Above 6.6×10¹⁹.** As each machine's current work ends, it starts a piece of
6.6×10¹⁹ to 2×10²⁰ (`next_m1max.sh`, `next_laptop.sh` and `next_d16w1.sh`,
`next_d16w2.sh`, `next_d16w5.sh`, `next_d16w6.sh`, `next_d16w7.sh`). The
pieces are sized to each machine's speed and start time, so that all should
finish around noon on 2026-09-27. The M1 Max frees up first, so it takes the
lowest range. Once warm, though, the M1 Max ran at about 438 T/s instead of
614, and the Intel Mac lost 90 minutes asleep, which would have left their
pieces finishing hours after the rest. So that evening the laptop took on the
top of both, after its own piece (`laptop_queue2.sh`), and all three should
finish together at about 13:40 on 2026-09-27.

| Range | Machine | Status |
| ----- | ------- | ------ |
| 6.6×10¹⁹ to 9.81916×10¹⁹ | the M1 Max, GPU (d16w13) | started 2026-09-26 17:23 as 6.6×10¹⁹ to 1.07×10²⁰; about 438 T/s once warm, down from 614 as the MacBook heats up (2 CPU threads instead of 8 gave 428, so it went back to 8). Its watcher (`split2_d16w13.sh`) stopped it at 98,191,564,376,758,026,240, until the find; then (`split3_d16w13.sh`) at 92,036,734,533,827,297,280 (20,479,471 segments), so that the range below the find is finished sooner. Stopped there 09:47 on 2026-09-27, no chain: **a(16) confirmed** |
| 1.07×10²⁰ to 1.85×10²⁰ | the laptop, GPU (d16w14) | started 2026-09-26 18:49; **found a 17-prime chain at 134,491,669,675,212,201,371 at 00:49 on 2026-09-27**, after 35.25% of the piece |
| 1.85×10²⁰ to 1.914×10²⁰ | mathc (d16w15) | started 2026-09-26 18:51; stopped 2026-09-27 07:57 at 76.03%, above the find; mathc moved to A016045 |
| 1.914×10²⁰ to 1.949×10²⁰ | the M1 Pro (d16w16) | started 2026-09-26 18:52; stopped 2026-09-27 07:57 at 76.43%, above the find; the M1 Pro moved to A016045 |
| 1.949×10²⁰ to 1.9685×10²⁰ | the WSL laptop (d16w17) | started 2026-09-26 19:00; stopped 2026-09-27 07:57 at 78.81%, above the find; the WSL laptop moved to A016045 |
| 1.9685×10²⁰ to 1.9875×10²⁰ | math4 (d16w18) | started 2026-09-26 18:51; stopped 2026-09-27 07:57 at 75.21%, above the find; math4 moved to A016045 |
| 1.9875×10²⁰ to 1.99943×10²⁰ | the Intel Mac (d16w19) | started 2026-09-26 18:50 as 1.9875×10²⁰ to 2×10²⁰, but the Mac fell asleep seconds later (below) and ran only in brief maintenance wakes until 20:26. Stopped 2026-09-27 07:59 at 66.42%, above the find; the Intel Mac moved to A016045 |
| 9.20367×10¹⁹ to 9.81916×10¹⁹ | the Studio, GPU (d16w29) | started 2026-09-27 07:54, pausing the Studio's A016045 piece, to finish the range below the find sooner; finished 09:13, no chain |
| 9.81916×10¹⁹ to 1.07×10²⁰ | the laptop, GPU (d16w20) | started 2026-09-27 07:50. Its queue had stopped at the find, as designed, and was restarted by hand. Finished 09:46, no chain |
| 1.99943×10²⁰ to 2×10²⁰ | the laptop (d16w21) | not run: above the find |

```bash
./gapsieve -l 17 -s 6.6e19 -e 1.07e20 --first --gaps dec -c d16w13.ck -r --log d16w13.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 1.07e20 -e 1.85e20 --first --gaps dec -c d16w14.ck -r --log d16w14.log --gpu
```

```bash
./gapsieve -l 17 -s 1.85e20 -e 1.914e20 --first --gaps dec -c d16w15.ck -r --log d16w15.log
```

```bash
nice -n 19 ./gapsieve -l 17 -s 1.914e20 -e 1.949e20 --first --gaps dec -c d16w16.ck -r --log d16w16.log -t 8
```

```bash
./gapsieve -l 17 -s 1.949e20 -e 1.9685e20 --first --gaps dec -c d16w17.ck -r --log d16w17.log
```

```bash
./gapsieve -l 17 -s 1.9685e20 -e 1.9875e20 --first --gaps dec -c d16w18.ck -r --log d16w18.log
```

```bash
./gapsieve -l 17 -s 1.9875e20 -e 2e20 --first --gaps dec -c d16w19.ck -r --log d16w19.log
```

The laptop's queued pieces after d16w14:

```bash
nice -n 19 ./gapsieve -l 17 -s 98191564376758026240 -e 1.07e20 --first --gaps dec -c d16w20.ck -r --log d16w20.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 199943325641973104640 -e 2e20 --first --gaps dec -c d16w21.ck -r --log d16w21.log --gpu
```

**After the find.** A263049's search was to go on to 3×10²⁰, but at
00:49 on 2026-09-27 the laptop found a chain at 1.345×10²⁰, and nothing above
it can change a(16). So at 07:50 the search narrowed to the one range still
open below the find, 6.6×10¹⁹ to 1.07×10²⁰. Studio joined in (d16w29), and the
M1 Max's split moved down to 9.20367×10¹⁹, so that all of it should be
searched by about 09:50. It was, at 09:47, with no other chain: a(16) is
confirmed. The five machines searching above the find (mathc,
the M1 Pro, the WSL laptop, math4 and the Intel Mac) moved to A016045 at
07:57, with `switch_to_a16.sh`. The laptop and the M1 Max follow once the
range below the find is done.

The queue design caught the find as intended. The laptop's queue
(`laptop_queue2.sh`) started nothing after d16w14 reported a candidate. That
left the next piece down, d16w20, unstarted until it was restarted by hand
seven hours later. So a queue should stop only the work above a find, not
the work below it.

```bash
./gapsieve -l 17 -s 92036734533827297280 -e 98191564376758026240 --first --gaps dec -c d16w29.ck -r --log d16w29.log --gpu
```

```bash
nice -n 19 ./gapsieve -l 17 -s 98191564376758026240 -e 1.07e20 --first --gaps dec -c d16w20.ck -r --log d16w20.log --gpu
```

A Mac without Apple's command-line tools can still run gapsieve: build it on
another Mac for that one's architecture and macOS version, copy the binary
over, and check it with the fingerprint command. For a 2016 Intel MacBook Pro
on macOS 12:

```bash
cc -O3 -arch x86_64 -march=skylake -mmacosx-version-min=12.0 -pthread -o gapsieve gapsieve.c
```

### Search status

As of 2026-09-28 08:55 PDT.

**A263049**: done. **a(16) = 134,491,669,675,212,201,371**, confirmed at
09:47 on 2026-09-27.

**A016045**: done through a(17). **a(16) = 254,345,176,718,302,423,991**,
found by the laptop and confirmed at 16:54 on 2026-09-27. **a(17) =
307,223,174,679,659,356,577**, found by atom1 at 14:39 and confirmed at 17:08,
once every start below it had been searched. The next term, a(18), needs a
chain of 19 primes. No search for it is running.

**OEIS**: both entries now include all five new terms: A016045 through a(17)
(revision 43) and A263049 through a(16) (revision 50). See
[*In the OEIS*](#in-the-oeis).

**The fleet.** The WSL laptop and the Intel Mac were taken out of the fleet at
17:15 on 2026-09-27, after their last pieces had finished. Future searches use
the other eleven machines.

## How the search works

`primesieve` streams primes in order and the program slides a window of
`W = --length + 8` consecutive primes across that stream, testing at each
position whether the forward gaps are 2, 4, 6, ...

The useful property is that a prime sitting *inside* a chain fails on its very
first test — its next gap is already 4 or more — so every chain is reported
exactly once from its head, and no backward context is ever needed. That makes
the range splittable into independent segments, which is what makes the whole
thing both trivially threadable and cheap to checkpoint.

## Two tools, one search

The directory holds two programs with the same CLI and the same output:

- **`gapchain_ps`** enumerates every prime with primesieve and slides a gap
  window across the stream. Simple, exact, and the reference implementation —
  but it pays to look at every prime, which caps it at ~26 billion
  integers/second near 2×10¹⁷.
- **`gapsieve`** never enumerates primes. It sieves for the chain *pattern*
  directly: a wheel over the primes 2..19 discards ~99.87% of integers as
  impossible chain heads outright (p ≡ 2 mod 3 is forced, and so on). It
  then sieves every member of every candidate against the primes up to 199
  using precomputed bit patterns, a whole 64-bit word at a time, and screens
  the rare survivors with batched base-2 strong tests. A full hit is then
  verified *exactly*: every member deterministically prime (the 12-base
  test, valid far beyond 2⁶⁴), no interloper prime between members, and the
  chain extended to its true length. Runs at ~145 **trillion**
  integers/second on the same machine, ~5,000× faster, in ~17 MB of memory.
  See [README_SIEVE.md](README_SIEVE.md) for the details.

The sieve can only ever strike a candidate whose member is a genuine multiple
of a sieve prime, so no true chain is lost; the one theoretical exception
(a member *equal* to a sieve prime) only exists below the sieve depth, and
that strip is searched exactly by a direct scan. The `--presieve` and
`--depth` knobs trade sieve work against screening work. Results are
identical at any values.

Both tools have been verified to produce byte-identical output on every range
tested, from `3..10⁶` through windows at 2.2×10¹⁷, for increasing and
decreasing chains alike, and `gapsieve` reproduces the known terms of both
sequences. Use `gapsieve` for real searches; keep `gapchain_ps` as an
independent cross-check. To confirm a hit, rerun a small window around it
with the same `--gaps`.

## Building

```bash
make              # both programs
make gapsieve     # just the sieve, which needs nothing else
make clean
```

On a Mac, `make` builds `gapsieve` with the GPU search (`--gpu`). It compiles
`gapsieve_gpu.m`, which includes `gapsieve.c`, links Apple's Metal framework,
and embeds the kernels from `gapsieve.metal` with `xxd`, so the binary stands
alone. `make GPU=0` builds the CPU-only version. Linux builds are unchanged.

The Makefile finds primesieve (needed only by `gapchain_ps`) through
`pkg-config`, then Homebrew, then the compiler's default paths. It also
compiles for the CPU doing the build (`-march=native`, or `-mcpu=native` where
a compiler only knows that spelling). On a 32-thread x86 Linux machine that
doubled `gapsieve`'s speed, from 55 to 112 T/s, because compilers otherwise
target a generic x86-64 without AVX. On Apple Silicon it makes no difference.

A binary built this way may not run on an older CPU, so run `make` on each
machine rather than copying binaries around. `make NATIVE=` builds a portable
binary. After changing flags, run `make clean` first, since make only notices
source changes. The manual commands below build the portable way; add
`-march=native` to match the Makefile.

`gapsieve` has no dependencies:

```bash
cc -O3 -pthread gapsieve.c -o gapsieve
```

With the GPU search, on a Mac:

```bash
xxd -i gapsieve.metal > gapsieve_metal.h
```

```bash
cc -O3 -fobjc-arc -pthread -o gapsieve gapsieve_gpu.m -framework Metal -framework Foundation
```

For `gapchain_ps`, the only dependency is
[primesieve](https://github.com/kimwalisch/primesieve)
(tested against 12.14). If `pkg-config` can see it, one line works on both
platforms:

```bash
cc -O3 -pthread gapchain_ps.c $(pkg-config --cflags --libs primesieve) -lm -o gapchain_ps
```

### macOS

```bash
brew install primesieve
```

Homebrew on Apple Silicon installs under `/opt/homebrew`, on Intel under
`/usr/local`. Apple's `cc` does not search either by default, so pass the
paths explicitly if you are not using `pkg-config`:

```bash
cc -O3 -pthread -I/opt/homebrew/include gapchain_ps.c \
   -L/opt/homebrew/lib -lprimesieve -o gapchain_ps
```

If your editor reports `primesieve.h file not found` while the command line
build succeeds, it is clangd missing the include path. Drop a
`compile_flags.txt` next to the source containing `-I/opt/homebrew/include`.

### Linux

Install the development package:

```bash
sudo apt install libprimesieve-dev      # Debian, Ubuntu
sudo dnf install primesieve-devel       # Fedora, RHEL
sudo pacman -S primesieve               # Arch
```

Then the headers and library are on the default search path:

```bash
cc -O3 -pthread gapchain_ps.c -lprimesieve -lm -o gapchain_ps
```

If your distribution's package is old, or you want to tune the build for your
CPU, compile primesieve from source and point at it:

```bash
cc -O3 -pthread -I/path/to/primesieve/include gapchain_ps.c \
   -L/path/to/primesieve/lib -lprimesieve -lstdc++ -lm -o gapchain_ps
```

`-lstdc++` is only needed when linking primesieve statically — it is a C++
library behind a C API.

### Linux with an NVIDIA GPU (CUDA)

On a Linux machine with an NVIDIA GPU and the CUDA toolkit, `make CUDA=1`
builds `gapsieve` with the CUDA GPU search (`--gpu`), the same search the Macs
run on Metal. It was written for the DGX Spark (GB10, compute capability
12.1, CUDA 13). The Makefile compiles for the GPU in the machine doing the
build (`-arch=native`, which needs CUDA 11.5 or later). To build for another
GPU, set `CUDA_ARCH`, for example `make CUDA=1 CUDA_ARCH=-arch=sm_87` for a
Jetson Orin.

```bash
make CUDA=1 gapsieve
```

By hand:

```bash
cc -O3 -mcpu=native -pthread -c -o gapsieve_cuda_host.o gapsieve_cuda_host.c
nvcc -O3 -arch=sm_121 -o gapsieve gapsieve_cuda.cu gapsieve_cuda_host.o -lpthread
```

Check a new GPU with the fingerprint command in *Continuing above 2⁶⁴*, run
with `--gpu -t 0`.

## Usage

```
gapchain_ps -l LEN -e END [options]
```

Every parameter is a flag; nothing is positional, and flags may appear in any
order. Each has a short and a long form.

### Required

| Flag | Meaning |
| ---- | ------- |
| `-l N`, `--length N` | Chain length: the minimum number of **primes** in a run before it is reported. Must be 2–1000. This is `n+1` in A016045 terms. |
| `-e N`, `--end N`    | Highest prime to consider as a chain head, inclusive. `gapsieve` also takes `max`; see below. |

### Optional

| Flag | Meaning |
| ---- | ------- |
| `-s N`, `--start N`      | Lowest prime to consider as a chain head. Default 3; values below 3 are clamped to 3. |
| `-t N`, `--threads N`    | Worker threads. Default: number of online CPUs. |
| `-g W`, `--gaps W`       | Which chains: `inc` (2, 4, 6, ...; A016045; the default), `dec` (..., 6, 4, 2; A263049) or `both`. `both` does two searches in one pass and takes as long as two separate runs. |
| `-f`, `--first`          | Report only the smallest qualifying chain (for `both`, one of each shape), then stop. See below. `gapchain_ps` supports this for `inc` only. |
| `-c F`, `--checkpoint F` | Checkpoint to `F` periodically and on SIGINT/SIGTERM. |
| `-i S`, `--interval S`   | Checkpoint interval in seconds. Default 60. `0` writes only on exit. |
| `-r`, `--resume`         | Resume from `F`. Required if `F` already exists. |
| `--gpu`                  | `gapsieve` on a Mac only: also search on the GPU (Metal on a Mac, CUDA with `make CUDA=1`), nine to twelve times as fast as the CPU alone. `-t` then sets the CPU threads beside it: all cores but two by default, `0` for the GPU alone. See [README_SIEVE.md](README_SIEVE.md#on-the-gpu-macos). |
| `-q`, `--quiet`          | Suppress progress output. |
| `-h`, `--help`           | Print usage. |

`--start` and `--end` bound the chain **head**, not its members. A chain headed
at or just below `--end` is still reported in full even though its later primes
lie past it, so splitting a search into adjacent ranges neither loses a chain
nor cuts one in half.

`--length` is a floor, not an exact match: asking for 9 also reports every
chain of 10, 11, and longer. It is *not* a runtime lever — scanning speed is
effectively identical for 8 or 800, because almost every prime fails the test
on its first gap either way. What it really controls is how much output you
get, and the range width is what controls how long the search takes.

#### Number formats for `--start` and `--end`

Both accept, in any mix:

| Form                  | Example               | Value |
| --------------------- | --------------------- | ----- |
| plain digits          | `113575727`           | 113,575,727 |
| digit separators      | `1_000_000`, `1,000,000` | 1,000,000 |
| scientific / e notation | `1e14`, `2.5e15`, `3.0e15` | 10¹⁴, 2.5×10¹⁵, 3×10¹⁵ |
| magnitude suffix      | `500K`, `20M`, `3G` or `3B`, `300T`, `2P` | 10³, 10⁶, 10⁹, 10¹², 10¹⁵ |

`gapsieve` converts every form exactly, rounding a fraction to the nearest
integer, so `1e20` is exactly 10²⁰ and `2.5e15` exactly 2.5×10¹⁵.
`gapchain_ps` converts digit-only strings exactly, so a literal value above 2⁵³
such as `12345678901234567` survives intact, but its other forms go through
floating point and are rounded to the nearest integer. That's fine for round
magnitudes like `1e14`; use plain digits there when you need an exact odd
starting point.

`gapchain_ps` accepts values up to 1.844 × 10¹⁹. `gapsieve` accepts numbers
of any size and searches as high as its proven primality test allows, just
below ψ₁₃ ≈ 3.317 × 10²⁴. `--end max` asks for exactly that: `gapsieve`
lowers any end above the limit to the last start whose chain it can check
(ψ₁₃ − 9,507 at `-l 17`, ψ₁₃ − 273 for `--gaps dec` alone) and says so.

### Finding only the first chain: `--first`

`--first` answers "what is the **smallest** prime heading a chain of at least
this length?" — the A016045 question — and stops as soon as that is settled
rather than scanning to `--end`.

```bash
./gapchain_ps -l 11 -e 2e11 --first
```

```
11: 137168442221 137168442223 137168442227 ... 137168442331

first chain of length >= 11 is at 137168442221, after searching 68.58% of the range in 00:00:07
```

Only the winning chain is printed to stdout — not every chain encountered
along the way — so the output is the answer and nothing else.

**Stopping correctly is the subtle part.** Threads work on different segments
at once, so the first hit *found* is very often not the first hit *in the
range*: a thread scanning segment 7 can report while segments 3 and 5 are
still being scanned. Stopping there would return the wrong prime.

Because segments are handed out in ascending order, the moment any segment `b`
produces a hit the search can conclude:

- every **unclaimed** segment is above `b`, so no new segment is started;
- every **in-flight segment above `b`** is irrelevant and is abandoned
  immediately, mid-scan;
- every **in-flight segment below `b`** is allowed to run to completion,
  because it may still hold an earlier chain.

Whenever one of those earlier segments does produce a better hit, it replaces
the candidate and the same rules apply again from the new, lower segment. The
run ends when no thread is left that could beat the current best, so the
reported prime is the true minimum, not merely the first one noticed. You can
watch this happen — each improvement is announced on stderr:

```
candidate: 40079796917 (length 7, segment 2); finishing earlier segments
candidate: 40066317017 (length 8, segment 1); finishing earlier segments
candidate: 40055475827 (length 7, segment 0)
```

That trace also shows that length does not break ties: the length-8 chain in
segment 1 loses to a length-7 chain in segment 0, because `--first` ranks by
position, and both satisfy the `-l 7` floor.

A worker also stops scanning its own segment at its first hit, since it scans
in increasing order and anything further on is larger.

If the range holds no qualifying chain, `--first` says so and exits 0 — a
complete search that found nothing is a legitimate result, not an error.

`--first` combines with checkpointing, and the mode is recorded in the
checkpoint so a `--first` run cannot be resumed as a counting run or vice
versa. Note that a `--first` checkpoint carries no histogram: only segments
that were scanned to the end count toward the frontier, and the segment
holding the winner was deliberately cut short.

### Exit codes

| Code | Meaning |
| ---- | ------- |
| 0    | Search finished. Under `--first` this covers both "found it" and "no such chain in the range" — check stdout, which is empty in the latter case. |
| 1    | Fatal error: unreadable or mismatched checkpoint, out of memory. |
| 2    | Bad command line. |
| 130  | Interrupted by SIGINT or SIGTERM; a checkpoint was written if `-c` was given. |

## Output

Chains go to **stdout**, one per line: the length, then the primes.

```
 9: 113575727 113575729 113575733 113575739 113575747 113575757 113575769 113575783 113575799
```

With `--gaps dec` or `both`, every line starts with its shape. The increasing
format above stays unchanged when only `inc` is searched:

```
inc 11: 137168442221 137168442223 137168442227 ... 137168442331
dec 11: 245333213 245333233 245333251 245333267 ... 245333323
```

A decreasing chain is printed whole, from its first prime, including the
primes before the pattern if it extends backwards.

A trailing `+` means the chain filled the entire `--length + 8` window, so its
true length may be greater than the number shown. To resolve one, rerun that
neighbourhood with a larger `--length` — the window grows with it.

Under `--first`, stdout carries exactly one such line: the winning chain.

Progress, the header, and the final statistics go to **stderr**, so
`> chains.txt` captures exactly the chains and nothing else.

### Progress and ETA

When stderr is a terminal, a single line repaints four times a second:

```
  3.14%  at 1.00941e+14  15.65G/s  eta 00:30:56  elapsed 00:01:00  64 chains
```

The fields are percent of the range completed, the position the completed
prefix has reached, numbers scanned per second, estimated time remaining,
elapsed time, and chains found so far. The rate and ETA are measured from the
start of *this* run, so after a resume they describe the remaining work, not
the whole range.

When stderr is redirected to a file, the same line is written once a minute
with a newline instead — useful for a `nohup`'d multi-day run. `-q` silences
both.

### Final statistics

Printed on normal completion and on Ctrl-C:

```
completed 100.00% of 3..200000000000 (11921/11921 segments), 00:00:09 elapsed
chains by length:
       7            6551   87.69%
       8             823   11.02%
       9              80    1.07%
      10              15    0.20%
      11               2    0.03%
   total            7471
```

Only non-empty lengths are listed. The highest possible bucket is
`--length + 8`, and if it is non-empty it is labelled with a `+` to match
the per-chain marker — those chains filled the window and may be longer. The
run above found nothing that long, so its buckets stop at 11.

## Checkpointing

Long searches are meant to be interruptible. Add `--checkpoint FILE`:

```bash
./gapchain_ps -l 16 -s 2.2e17 -e 5e17 --checkpoint run.ck > chains.txt
```

The file is rewritten every 60 seconds (`--interval` to change), and once more
on exit — whether that exit is normal completion or Ctrl-C.

To resume, repeat the *same* search and add `--resume`:

```bash
./gapchain_ps -l 16 -s 2.2e17 -e 5e17 --checkpoint run.ck --resume >> chains.txt
```

Note the `>>`. The program never touches your output file, so append unless
you want to discard the chains found before the interruption.

### How it works, and what it guarantees

The range is cut into fixed-size segments handed to threads in order. The
segment size is chosen at startup — 2²⁴ integers for low searches, growing
with √`end` up to 2³² for high ones, because each segment pays an O(√`end`)
iterator setup that must be amortized (at 2.2×10¹⁷, 2²⁴ segments spend 80% of
their time on setup). The chosen size is recorded in the checkpoint, and
resume adopts whatever the file says, so checkpoints written by older builds
still work. A checkpoint stores the number of *leading* segments that are
completely finished, plus the chain histogram for exactly those segments.
Segments that finish out of order wait to be counted until the frontier
reaches them, so the file always describes a consistent prefix of the range,
never a patchwork.

Ordering is deliberate: the frontier is sampled first, then stdout is flushed,
then the file is written to `FILE.tmp`, fsynced, and renamed into place. So
every chain the checkpoint claims to have found is already on disk, and a
crash mid-write cannot corrupt an existing checkpoint.

The guarantee on resume is **no gaps, but possibly duplicates**. Work between
the last checkpoint and the interruption is redone, so those chains are
printed a second time. Sort and `uniq` if you need exact counts from a
resumed output file:

```bash
sort -u chains.txt | wc -l
```

A clean Ctrl-C redoes very little. Each CPU thread abandons the segment it
is on, the GPU finishes its units in flight, and the final checkpoint records
everything finished below the lowest abandoned segment. Segments that finished
beyond it are searched again on resume, so their chains are printed twice.
With CPU threads alone that's rare, since they work through the segments
almost in step. With `--gpu` and CPU threads it is normal, because the GPU runs
far ahead of the CPU threads: stopped after 3 seconds, a test run at `-l 8`
printed 565 of its 2,258 chains twice. `--first` runs print nothing per
segment, so there it only means a fraction of a second of repeated work.
`SIGKILL`, a crash or a power loss repeats up to the last `-i` seconds of
work. A hard kill can also leave one
truncated line at the end of the output file, from a partially flushed buffer;
the chain in it is re-emitted in full on resume.

### Interrupting

The first Ctrl-C asks the workers to stop at their next segment boundary,
which takes under a second per thread. It then writes the checkpoint and
prints statistics. A second Ctrl-C exits immediately, leaving the previous
checkpoint intact.

### Safety rails

- If `FILE` exists and you did not pass `--resume`, the program refuses to
  start rather than overwrite a long run's progress.
- If `FILE` describes a different `--length`, range, or `--first` setting,
  resuming is a fatal error and both searches are printed for comparison.
- Resuming a finished checkpoint reprints its statistics and exits 0.
- The thread count is *not* part of the checkpoint; resuming with a different
  number of threads is fine.

### Checkpoint file format

Plain text, safe to read and to hand-edit if you know what you are doing:

```
gapchain_ps-checkpoint 1
k 9
lo 100000000000000
hi 130000000000000
segsize 16777216
nseg 1788140
first 0
done 6104
found 4
len 9 3
len 10 1
len 11 0
...
```

`k` is `--length`, `lo`/`hi` are the range, `segsize` is the segment size this
run was cut into (resume adopts it from the file), and `first` records whether
the run used `--first`. `done` is the completed-segment frontier, so progress
in absolute terms is `lo + done × segsize`. Lowering `done` by hand
re-searches from further back, which is a legitimate way to recover from a
truncated output file — the histogram will then over-count, but the chains
will be complete.

## Performance notes

Throughput depends strongly on how high you are searching, because sieving
itself slows down as the primes grow. Measured on Apple Silicon (14 threads
unless noted):

| Around   | gapchain_ps                 | gapsieve (l=16) |
| -------- | --------------------------- | --------------- |
| 10¹²     | ~33 billion/s (8 threads) |               |
| 10¹⁴     | ~16 billion/s             |               |
| 2.2×10¹⁷ | ~26 billion/s, 2.8 GB     | ~145 trillion/s, 17 MB |
| 8×10¹⁸   | ~9 billion/s, 5.3 GB      | ~145 trillion/s, 17 MB |

The 2.2×10¹⁷ figure depends on the adaptive segment sizing described above;
with the old fixed 2²⁴ segments the same machine managed only ~1.4 billion/s
there. Memory grows with √`end` because each thread's iterator holds the
sieving primes up to √`end`.

At the 10¹⁴ rate a 10¹⁵-wide window takes roughly 18 hours. Plan long runs
from a short measured sample near the range you actually care about rather
than from a figure taken lower down — and for the same reason, expect the ETA
to drift upward across a very wide range, since it extrapolates from the rate
observed so far.

Scaling with cores is close to linear. The only shared state is a mutex around
stdout, and it does not become a bottleneck even when chains are common: a run
printing 80,000 chains took the same wall time as one printing none over the
same range.

Memory is negligible — a window of `--length + 8` primes plus one
primesieve iterator per thread.

## Files

| File | |
| ---- | --- |
| `gapsieve.c`    | The constellation sieve — use this for real searches. |
| `gapsieve_gpu.m`, `gapsieve.metal` | The Mac GPU search (`--gpu`): the host side, which includes `gapsieve.c`, and the Metal kernels. |
| `gapsieve_cuda_host.c`, `gapsieve_cuda.cu`, `gapsieve_cuda.h` | The CUDA GPU search (`--gpu`, `make CUDA=1`): the host side, which includes `gapsieve.c`; the CUDA kernels; and the interface between them. |
| `tail_check.py` | The check of the last million starts below 2⁶⁴, for the a(16) > 2⁶⁴ bounds. |
| `verify_chain.py` | An independent check of a chain, with decreasing gaps by default or increasing gaps with `inc`: every number from its head to its last member, by Miller-Rabin and by Baillie-PSW. `python3 verify_chain.py 134491669675212201371 17`, `python3 verify_chain.py 254345176718302423991 17 inc` |
| `submission.md` | Every OEIS submission from these searches: the ones already approved, and ready-to-paste text for the rest, with the data and checks behind each. |
| `gapchain_ps.c` | The enumeration search, built on primesieve; the independent cross-check. |
| `gapsieve`, `gapchain_ps` | The compiled binaries. |
| `Makefile`      | `make` builds both binaries; `make clean` removes them. |
| `README.md`     | This file. |
