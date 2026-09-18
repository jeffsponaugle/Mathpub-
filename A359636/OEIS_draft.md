# A359636 — draft OEIS submission text

STATUS: FINAL (Sep 18 2026). a(9) = 31610535900218923 was found by an exhaustive
scan of every m == 3 (mod 6) up to 3.161e16 (GPU version of the tool, 10 h 42 min
on a DGX Spark) and verified three independent ways (see README.md).

## DATA

    7, 19, 643, 51427, 8083633, 1077940147, 75582271489, 34710483181813, 31610535900218923

(OEIS terms field: the same list, comma separated; a(9) has 17 digits, so the
DATA line stays within the usual limit.)

## b-file  (b359636.txt)

    1 7
    2 19
    3 643
    4 51427
    5 8083633
    6 1077940147
    7 75582271489
    8 34710483181813
    9 31610535900218923

## EXAMPLE (addition)

    a(9) = 31610535900218923: 31610535900218924 = 2^2*11*17*19*53*79*149*359*9931,
    31610535900218925 = 3*5^2*13*23*31*37*157*1559*5021 and
    31610535900218926 = 2*7*29*61*103*197*223*311*907 each have 9 distinct prime
    factors, and 31610535900218927 is prime.

## COMMENTS (suggested)

1. Replace the existing bound comment

       a(9) <= 76340177205657727, a(10) <= 225096507194749219819. - David A. Corneth, Jan 12 2023

   by

       a(10) <= 225096507194749219819. - David A. Corneth, Jan 12 2023

   (or keep it and add "a(9) = 31610535900218923 confirmed by exhaustive search" after it).

2. New comment:

       Every qualifying gap contains an odd multiple of 3, m, with m-1, m, m+1 all
       composites of the gap (m = p+2 for a gap of 4, since then p == 1 (mod 3)), so
       a(n) can be found by searching m == 3 (mod 6) with omega(m-1), omega(m),
       omega(m+1) all >= n and then checking the surrounding prime gap. - Jeff Sponaugle, Sep 2026

3. New comment (observation):

       The gaps at a(1)-a(9) all have length 4.

4. New comment (lower bound for a(10), a by-product of the a(9) search):

       a(10) > 31610555571634177. - Jeff Sponaugle, Sep 18 2026

   Justification (for the submitter, not for the entry): the exhaustive level-9
   scan tested every m == 3 (mod 6) with m <= 31610555571634179 and recorded all
   98 values with omega(m-1), omega(m), omega(m+1) >= 9.  A level-10 gap would
   have to contain such an m with all three omegas >= 10; among the 98 only four
   have a single member with omega = 10 (m = 5791430882243721, 8382386489666475,
   21469520444799855, 24048717372672915) and none has all three, so no gap with
   p <= 31610555571634177 qualifies at level 10.  Corneth's a(10) <= 225096507194749219819
   stands as the upper bound.

## EXTENSIONS

    a(9) from Jeff Sponaugle, Sep 18 2026

## PROG / LINKS (optional)

    (C) See the a-file / GitHub link: a359636.c, a three-target sieve over
    m == 3 (mod 6) with exact verification of survivors; also reproduces a(1)-a(8).

## Notes for the submitter

* The exhaustive scan output (`gpu_n9.txt`), progress log (`gpu_n9.log`) and
  state file (`gpu_n9.state`, final frontier: all m <= 31610555571634179 covered)
  document the search; 98 triples m with omega(m-1), omega(m), omega(m+1) >= 9
  occur below the term, only one of them bounded by primes.
* Independent verification of the gap: `python3 verify_a359636.py verify 9 31610535900218923`.
* Corneth's a(10) bound is beyond 64-bit arithmetic; this tool does not address a(10).
