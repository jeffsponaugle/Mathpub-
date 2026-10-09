# OEIS submissions

Every OEIS submission this project's data supports, checked against the live
OEIS on **October 3, 2026**. For each entry the current page was read and its
b-file downloaded and compared term by term with the copy taken on
September 29. The data files are in [`oeis/`](oeis/), and
[`oeis/README.md`](oeis/README.md) has the evidence behind each.

## Status

| # | sequence | what it counts | on the OEIS now | proposed | status |
|---|---|---|---|---|---|
| 1 | [A090009](https://oeis.org/A090009) | the earliest chain of ≥ n primes | a(1)–a(11); a(10), a(11) and comment from you | — | **done, approved** |
| 2 | [A320879](https://oeis.org/A320879) | chains of ≥ 8 primes | your b-file, 4,379 terms, with the 16 missing terms | (link wording, optional) | **done, approved** |
| 3 | [A320880](https://oeis.org/A320880) | chains of ≥ 9 primes | your b-file, 886 terms < 10¹⁸ | — | **done, approved** |
| 4 | [A320878](https://oeis.org/A320878) | chains of ≥ 7 primes | your b-file, 10,000 terms | — | **done, approved** |
| 5 | new | chains of ≥ 10 primes | does not exist | new entry, 665 terms < 10²² | **ready** (re-check running) |
| 6 | [A048527](https://oeis.org/A048527) | chains of exactly 6 primes | b-file, 3,000 terms (Blomberg) | 10,000 terms | **ready** |
| 7 | [A048523](https://oeis.org/A048523) | chains of exactly 2 primes | b-file, 1,000 terms (Dale) | 10,000 terms | **ready** |
| 8 | [A048519](https://oeis.org/A048519) | chains of ≥ 2 primes | b-file, 1,000 terms (Dale) | 10,000 terms | **ready** |
| 9 | [A320881](https://oeis.org/A320881) | numbers p + digitsum(p) | 63 terms, no uploaded b-file | 10,000 terms | **ready** |
| 10 | [A048524](https://oeis.org/A048524) | chains of exactly 3 primes | b-file, 10,000 terms (Blomberg) | 20,000 terms | optional |
| 11 | [A048525](https://oeis.org/A048525) | chains of exactly 4 primes | b-file, 10,000 terms (Blomberg) | 20,000 terms | optional |
| 12 | [A048526](https://oeis.org/A048526) | chains of exactly 5 primes | b-file, 10,000 terms (Blomberg) | 20,000 terms | optional |
| 13 | A090009, A320880 | — | — | cross-references to the new entry; extend Hasler's comment | after #5 |
| 14 | Puzzle 163 (primepuzzles.net) | the earliest chain of 9..13 primes | 9 (Resta), 10 (Wroblewski, 1.37×10²⁷) | 10, 11, and none of 12, 13 below 10³⁰ | **ready** (email) |
| — | A048520, A047791, A062028, A230093 | — | complete enough | nothing | not needed |

"Chains" in A048519..A048527 and A320878..A320880 include starts inside
longer chains, as those entries define them. A320878..A320880 are named by the
primes the iteration yields ("yields 6 primes in a row" is a chain of 7), and
A048523..A048527 by iterations ("only k iterations" is a chain of exactly
k + 1).

## Done

### 1. A090009

On the OEIS now, approved: DATA `2, 11, 11, 277, 37783, 516493, 286330897,
286330897, 56676324799, 63782998969989799429, 8968999974999983898865998461`;
your comment (parity, mod 3, mod 5, "a(12) > 10^30"); EXTENSIONS "a(10)-a(11)
from Jeff Sponaugle, Sep 29 2026"; keywords `base,more,nonn,hard`. Its b-file
is the one OEIS synthesizes from the DATA line. Nothing more to do until #5
has an A-number (see #13).

### 2. A320879 — correction

On the OEIS now: "Jeff Sponaugle, Table of n, a(n) for n = 1..4379 (terms
1..445 from Lars Blomberg, except 197..212 corrected)". The b-file is
identical to `oeis/b320879.txt`, including the 16 terms Blomberg's file
lacked.

Optional wording fix: Blomberg's terms were not wrong at 197..212. His file
skipped 16 terms there, so his terms 197..445 are now 213..461. A more exact
note:

> (terms 1..196 and 213..461 from Lars Blomberg; 197..212, missing from his
> file, added)

### 3. A320880

On the OEIS now: "Jeff Sponaugle, Table of n, a(n) for n = 1..886 (Terms <
10^18)"; EXTENSIONS "a(16) onward from Jeff Sponaugle, Sep 29 2026". The
b-file is identical to `oeis/b320880.txt`.

### 4. A320878

On the OEIS now: "Jeff Sponaugle, Table of n, a(n) for n = 1..10000 (first 200
from M. F. Hasler, next 7426 from Lars Blomberg)". The b-file is identical to
`oeis/b320878.txt` (the 10,000 terms are all terms up to 136510061910829).

## Ready

### 5. New sequence: chains of at least 10 primes

Not on the OEIS: a search for its second term, 96989986598717899429, finds
nothing, and a(10) appears only in A090009. It is the next member of
A320878..A320880, and Hasler's comment in A090009 ("a(k) = A[k](1) for A[3..9]
= A048524 .. A048527, A320878 .. A320880") would continue with it as A[10].

File: `oeis/b_new10.txt`, all **665 terms below 10²²**, found on the Mac
Studio (`oeis/len10all.log`: 100,476,518,473,002 candidates, every one
tested). Every term is proven to start a chain of at least 10 primes with
Pratt certificates (`verify_chain.py`). An independent re-check with the other
GPU kernel is running on the laptop (`oeis/new10-recheck.log`, about three
hours); submit once it reports the same 665. A heuristic count predicted
about 654.

Draft entry:

```
%N Primes such that iteration of A062028 (n + its digit sum) yields 9 primes in a row.
%S 63782998969989799429,96989986598717899429,268997958095988899329,296888698475698799329,
%T 478776289862999899429,498077995966789799329,586289088899957899429,679957758586388899429,
%U 688177977888837899329,736794879935979899329,789867584696559799429
%O 1,1
%C a(1) = A090009(10) = start of the first chain of 10 primes under iteration of A062028.
%C Every term below 10^22 starts a chain of exactly 10 primes: the first chain of 11 primes starts at A090009(11) = 8968999974999983898865998461.
%C All terms below 10^22 end in 9329 or 9429. A chain of 10 primes needs digit sums that stay even through a carry into the higher digits; below 2^64 only five patterns of the last four digits (9319, 9329, 9419, 9429, 9743, each with a prescribed digit sum and number of trailing 9s above them) allow that, and none of the 176288301 numbers below 2^64 that fit one of them starts a chain of 10 primes.
%F a(n) = A320880(k) for those k with A062028(A320880(k)) also in A320880.
%e a(1) = 63782998969989799429: adding digit sums 142, 140, 136, 146, 148, 116, 124, 122 and 118 gives the primes 63782998969989799571, ...799711, ...799847, ...799993, ...800141, ...800257, ...800381, ...800503 and 63782998969989800621, whose digit sum, 119, is odd, so the chain ends there with 10 primes.
%o (PARI) is(n, k=10) = for(i=1, k, if(!isprime(n), return(0)); n += sumdigits(n)); 1
%H Jeff Sponaugle, <a href="/A3xxxxx/b3xxxxx.txt">Table of n, a(n) for n = 1..665</a> (terms < 10^22)
%H Carlos Rivera, <a href="https://www.primepuzzles.net/puzzles/puzz_163.htm">Puzzle 163. P+SOD(P)</a>, The Prime Puzzles & Problems Connection.
%Y Cf. A062028, A090009, A048519, A048523, A048524, A048525, A048526, A048527, A320878, A320879, A320880.
%K nonn,base
%A Jeff Sponaugle, Oct 03 2026
```

(The `%F` formula matches the style of A320880's "A320880 = { n in A320879 |
A062028(n) in A320879 }".)

### 6. A048527 — b-file to 10,000 terms

On the OEIS now: "Lars Blomberg, Table of n, a(n) for n = 1..3000", ending at
a(3000) = 265337253901. File: `oeis/b048527.txt`, 10,000 terms, ending at
a(10000) = 1332424832891; it begins with Blomberg's 3,000 unchanged. Checked:
`dschain` against `dspattern`. Suggested link: "Jeff Sponaugle, Table of n,
a(n) for n = 1..10000 (first 3000 from Lars Blomberg)".

### 7. A048523 — b-file to 10,000 terms

On the OEIS now: "Harvey P. Dale, Table of n, a(n) for n = 1..1000", ending at
a(1000) = 83449. File: `oeis/b048523.txt`, 10,000 terms, ending at a(10000) =
1310999; begins with Dale's 1,000 unchanged. Checked: a brute-force program
(byte sieve, direct digit sums). Suggested link: "Jeff Sponaugle, Table of n,
a(n) for n = 1..10000 (first 1000 from Harvey P. Dale)".

### 8. A048519 — b-file to 10,000 terms

On the OEIS now: "Harvey P. Dale, Table of n, a(n) for n = 1..1000", ending at
a(1000) = 70111. File: `oeis/b048519.txt`, 10,000 terms, ending at a(10000) =
1101937; begins with Dale's 1,000 unchanged. Checked: the definition computed
directly (a sieve and p + digitsum(p)). Suggested link: "Jeff Sponaugle,
Table of n, a(n) for n = 1..10000 (first 1000 from Harvey P. Dale)".

### 9. A320881 — first b-file, 10,000 terms

On the OEIS now: only the 63 DATA terms, with the b-file synthesized from them.
File: `oeis/b320881.txt`, 10,000 terms, ending at a(10000) = 105121; begins
with the 63 published terms. Checked: its primes agree with A048520's
independent b-file (M. F. Hasler) for all 4,259 terms in range. Suggested
link: "Jeff Sponaugle, Table of n, a(n) for n = 1..10000".

## Optional

### 10–12. A048524, A048525, A048526 — b-files to 20,000 terms

Each has "Lars Blomberg, Table of n, a(n) for n = 1..10000" now, the customary
size. Files `oeis/b048524.txt` (a(20000) = 36179369), `oeis/b048525.txt`
(a(20000) = 666370813) and `oeis/b048526.txt` (a(20000) = 34502052131), each
beginning with Blomberg's 10,000 unchanged and well under the 1 MB limit.
Checked: brute force (A048524), `dschain` against `dspattern` (A048525,
A048526). Suggested link: "Jeff Sponaugle, Table of n, a(n) for n = 1..20000
(first 10000 from Lars Blomberg)". Worth sending only if the editors want
longer tables.

## After the new entry is approved

### 13. Cross-references

- **A090009:** add the new A-number to CROSSREFS, and extend Hasler's comment
  so that A[10] is the new entry.
- **A320880:** add the new A-number to CROSSREFS (and to A320878, A320879,
  A048519, A048523..A048527 if the editors want the family fully linked).

## Not the OEIS

### 14. Rivera's Puzzle 163

The page (checked October 3) lists L = 9 (G. Resta, 2011) and L = 10 (J.
Wroblewski, October 2025, 1367618706414097919365699429, from a restricted
search), and nothing for 11, 12 or 13. It does not mention a(10) or a(11).
A note to Carlos Rivera could say:

- The earliest chain of 10 primes starts at 63782998969989799429 (exhaustive
  search; Wroblewski's chain is a later one).
- The earliest chain of 11 primes starts at 8968999974999983898865998461.
- No chain of 12 or more primes starts below 10^30, which also settles 13 up
  to that bound.
- Both are A090009(10) and A090009(11), with the method in its comments.

## Not needed

A048520 (10,000 terms, Hasler) and A047791 (10,000 terms, Zumkeller) already
have full b-files; A062028 and A230093 are the plain map n + digitsum(n) and
its preimage count.
