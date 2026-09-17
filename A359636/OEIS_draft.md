# A359636 — draft OEIS submission text

STATUS: DRAFT. The value of a(9) below is the smallest known qualifying gap
(found by the structured `hunt` search with a free largest prime on Sep 17 2026
and verified twice). It becomes a term only when the exhaustive scan
(`scan_n9.txt`) has covered every m <= 33955545649252305; if that scan finds a
smaller qualifying gap first, all occurrences of 33955545649252303 below must be
replaced by the smaller prime.

## DATA

    7, 19, 643, 51427, 8083633, 1077940147, 75582271489, 34710483181813, 33955545649252303

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
    9 33955545649252303

## EXAMPLE (addition)

    a(9) = 33955545649252303: 33955545649252304 = 2^4*7^2*11*17*29*37*163*1021*1297,
    33955545649252305 = 3*5*13*19*61*71*127*251*66383 and
    33955545649252306 = 2*23*41*47*53*113*181*313*1129 each have 9 distinct prime
    factors, and 33955545649252307 is prime.

## COMMENTS (suggested)

1. Replace the existing bound comment

       a(9) <= 76340177205657727, a(10) <= 225096507194749219819. - David A. Corneth, Jan 12 2023

   by

       a(10) <= 225096507194749219819. - David A. Corneth, Jan 12 2023

   (or keep it and add "a(9) = 33955545649252303 confirmed by exhaustive search" after it).

2. New comment:

       Every qualifying gap contains an odd multiple of 3, m, with m-1, m, m+1 all
       composites of the gap (m = p+2 for a gap of 4, since then p == 1 (mod 3)), so
       a(n) can be found by searching m == 3 (mod 6) with omega(m-1), omega(m),
       omega(m+1) all >= n and then checking the surrounding prime gap. - Jeff Sponaugle, Sep 2026

3. New comment (observation):

       The gaps at a(1)-a(9) all have length 4.

## EXTENSIONS

    a(9) from Jeff Sponaugle, Sep 2026   [use the exact date of the confirmed scan, OEIS style "Sep 23 2026"]

## PROG / LINKS (optional)

    (C) See the a-file / GitHub link: a359636.c, a three-target sieve over
    m == 3 (mod 6) with exact verification of survivors; also reproduces a(1)-a(8).

## Notes for the submitter

* The exhaustive scan log (`scan_n9.log`) and output (`scan_n9.txt`) document the
  search: every m == 3 (mod 6) below the term was covered; the state file records
  the final frontier.  Quote the number of triples found below a(9) in the
  submission comments if useful.
* Independent verification of the gap: `python3 verify_a359636.py verify 9 33955545649252303`.
* Corneth's a(10) bound is beyond 64-bit arithmetic; this tool does not address a(10).
