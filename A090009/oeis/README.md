# OEIS contributions

Data computed with `dschain` and `dspattern` for A090009 and the sequences it
cross-references. a(10) and a(11) were submitted to A090009 in September 2026;
the rest has not been submitted yet. The table of submissions and their status
is in the main README, under [*OEIS submissions*](../README.md#oeis-submissions);
each section below says what to send and what backs it up. Every b-file is in the OEIS `n a(n)` format with
offset 1, and every one of them was checked twice, by two different programs
(see *Validation* below).

A note on counting, because the family's names differ: A048523..A048527 are
"primes for which only k iterations of p -> p + digitsum(p) are possible", i.e.
chains of exactly k+1 primes; A320878..A320880 are primes for which the
iteration "yields 6, 7, 8 primes in a row" beyond p, i.e. chains of at least 7,
8, 9 primes. Both count every prime, including those inside a longer chain.
A090009 counts the primes in the chain, from its true start.

## A090009 — new terms a(10), a(11) (submitted, September 2026)

```
a(10) = 63782998969989799429
a(11) = 8968999974999983898865998461
```

Neither number appears anywhere in the OEIS (searched September 2026).

- **Primality.** `python3 verify_chain.py <p>` proves every term of both
  chains prime with Pratt certificates (and 13-base Miller–Rabin where that
  is a proof), confirms the digit sums, that the chains end (the last digit
  sum is odd), and that no prime r has r + digitsum(r) = p.
- **Minimality.** a(10): two independent exhaustive searches of all
  3,577,235,998 candidates above a(9) (`dspattern`, and a separate GMP-based
  program) find exactly one chain of ten, a(10). a(11): the Mac Studio's
  search (`../a11.log`) and an independent re-search with a different GPU
  kernel (`../a11-recheck.log`) test all 18.4 trillion candidates below it and
  find exactly one chain of eleven or more, a(11). `indep_count.py`, which
  shares no code with `dspattern`, counts exactly the same 18,411,673,898,715
  candidates.
- **Method.** See `../README.md`, *Why long chains are rare*: a chain of L
  primes can only start at p = H·10^4 + x with (digitsum(H), trailing 9s of H,
  x) in a small, exhaustively computed set of patterns; below 2^64 there are
  five for L = 10 and none for L = 11.

Suggested additions to the entry:

- DATA: append `63782998969989799429, 8968999974999983898865998461`.
- COMMENT: "No chain of 12 or more such primes starts below 10^30." (One digit
  pattern allows 12 terms below 10^30; all 378,958,525 numbers fitting it have
  a term with a factor below 500 or failing a base-2 strong test.)
- EXTENSIONS: "a(10)-a(11) from <name>, Sep 29 2026".

Carlos Rivera's Puzzle 163 (primepuzzles.net), linked from A320878..A320880,
asks for the earliest start of chains of 9..13 primes. It lists a chain of 10
found in October 2025 by J. Wroblewski, 1367618706414097919365699429, from a
restricted search; a(10) above is far smaller, a(11) answers 11, and the
bound above answers 12 and 13 up to 10^30. (`dspattern` finds Wroblewski's
chain as the first hit in a narrow window around it.)

## New sequence: chains of at least 10 primes

The natural next member of A320878..A320880: "Primes such that iteration of
A062028 (n + its digit sum) yields 9 primes in a row", i.e. primes starting a
chain of at least 10 primes, counting starts inside longer chains. Its first
term is A090009's a(10) = 63782998969989799429. (There are none below 2^64,
and the same patterns that rule them out make them rare above: a heuristic
count, which matched observed counts closely, expects about 20 below 10^21
and 650 below 10^22.)

`b_new10.txt` has all **665 terms below 10^22**, found on the Mac Studio
(`len10all.log`: all 100,476,518,473,002 candidates tested, in 53 minutes),
against about 654 expected. Every term starts a chain of exactly 10 primes,
proven with Pratt certificates by `verify_chain.py`, and every one ends in
9329 or 9429. A re-check with the other GPU kernel is in `new10-recheck.log`.
The draft entry is in [`../submission.md`](../submission.md).

A "yields 10 primes in a row" sequence would start with a(11), but only one
term is known, well short of what a new entry needs.

## A320879 — correction and extension (chains of at least 8 primes)

The current b-file (445 terms, "Terms < 10^14") **omits 16 terms** between
20335560943853 and 22016355916853:

```
20717323653863  20860820504861  21014922200887  21074359310861
21146004494809  21219833610851  21351253491851  21414751343861
21430195082809  21450616717853  21521592413809  21530162924861
21550256446853  21601353248809  21900547065863  21917284010861
```

Each starts a chain of exactly 8 primes (all proven by `verify_chain.py`),
and all of them are in A320878's b-file, which is correct. `b320879.txt`
here has all **4,379 terms below 10^16**, including those 16.

## A320880 — extension (chains of at least 9 primes)

15 terms known (all below 10^14). `b320880.txt` has all **886 terms below
10^18**.

## A320878 — extension (chains of at least 7 primes)

7,626 terms (below 10^14). `b320878.txt` has the first **10,000 terms**, up
to a(10000) = 136510061910829.

## A048523..A048527 — b-file extensions (chains of exactly 2..6 primes)

| sequence | published b-file | here |
|---|---|---|
| A048523 (exactly 2) | 1,000 terms | `b048523.txt`, 10,000 terms, a(10000) = 1310999 |
| A048524 (exactly 3) | 10,000 terms | `b048524.txt`, 20,000 terms, a(20000) = 36179369 |
| A048525 (exactly 4) | 10,000 terms | `b048525.txt`, 20,000 terms, a(20000) = 666370813 |
| A048526 (exactly 5) | 10,000 terms | `b048526.txt`, 20,000 terms, a(20000) = 34502052131 |
| A048527 (exactly 6) | 3,000 terms | `b048527.txt`, 10,000 terms, a(10000) = 1332424832891 |

Each begins with the published b-file's terms, unchanged.

## A048519 — b-file extension (chains of at least 2 primes)

Primes p with p + digitsum(p) prime. The published b-file has 1,000 terms;
`b048519.txt` has **10,000**, up to a(10000) = 1101937, beginning with the
published ones unchanged.

## A320881 — b-file extension (numbers equal to a prime plus its digit sum)

Not a chain sequence, but in the same family and trivial to extend: its
b-file has only 63 terms. `b320881.txt` has **10,000**, up to a(10000) =
105121. (Every value up to N comes from a prime below N, so the list is
complete.) Its primes are A048520, whose independent 10,000-term b-file
agrees with ours for all 4,259 terms in range.

## Validation

- **The published data, reproduced.** Regenerated over their published
  ranges, A048523..A048527, A320878 and A320880 match their b-files term for
  term (44,641 terms); A320879 matches except for the 16 missing terms above.
- **The extensions, twice each.** A048519: `dschain` against the definition
  computed directly (a sieve and p + digitsum(p)). A048523 and A048524: a
  brute-force program
  (byte sieve, direct digit sums) sharing no code with either tool.
  A048525..A048527: `dschain` (every prime, primesieve) against `dspattern`
  (digit patterns). A320878..A320880: `dspattern`'s GPU list kernel against
  its CPU enumeration (A320878, A320879) and its first GPU kernel (A320880).
  Every term of A320879 and A320880 was also proven individually with
  `verify_chain.py`.
