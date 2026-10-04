# OEIS submissions from the multi-base Wieferich search

Status as of **Oct 03 2026**. The search is paused (atom2 lent to another sequence) with the 16
Phase-1 bases complete up to 3.0848e14. Everything below has been checked against the live OEIS
entries, their edit histories and R. Fischer's site on Oct 03 2026.

## 1. Already done

**A306256 (base 30): a(4) = 303632117562967 is in the OEIS.** Your edit (revisions #9–#14, Oct 02
2026) was approved the same day: Amiram Eldar adjusted the keywords (removed `bref`), Hugo
Pfoertner reviewed it and Michael De Vlieger approved it (revision #18). The live entry now has the
new term, your comment "a(4) was found by an exhaustive search of all primes below 3.04*10^14.",
the EXAMPLE and "a(4) from _Jeff Sponaugle_, Oct 02 2026". It still also carries the old comment
"No more terms up to 9.8*10^13." (see 3.1).

## 2. State of the OEIS entries (checked Oct 03 2026)

* **No changes since Sep 30** to A001220 and the 16 other base entries (A014127, A123692, A212583,
  A123693, A045616, A111027, A128667, A234810, A242741, A128668, A244260, A090968, A242982, A298951,
  A128669, A306255), and none of them has a pending edit.
* **A039951 has a pending edit by Richard Fischer** (revisions #128–#130, Sep 30 / Oct 02 2026):
  a new a-file link "Table of n = 1..10007, updated Sep 30 2026". It does not touch a(47) or a(72),
  but the OEIS allows only one pending edit per entry, so an edit of ours has to wait until his is
  approved.
* **Fischer's site:** the table for bases up to 10125 was updated Sep 30 2026; the table for bases up to
  1052 and his news page are still from Sep 2025 (base 30 still lists only 7, 160541, 94727075783;
  bases up to 149 searched to 2.00E+14).

## 3. What the search has established (verified)

* **Base 30:** every prime below 3.0848e14 has been tested. The independent scan `v30` covered [0, 3.0364e14)
  with 9,396,616,950,003 primes, equal to π(303640000000000) by primecount. Phase 1 covered [2e14, 3.0848e14).
  The only solutions are 7, 160541, 94727075783, 303632117562967, so **a(5) > 3.0848e14**.
* **Bases 6, 10, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 26, 47, 72:** every prime in [2e14, 3.0848e14)
  has been tested: 3,271,298,909,197 primes, which equals π(3.0848e14) − π(2e14) by primecount, in 10,848
  contiguous chunks, with 0 FLT errors. There are no solutions, and 119 chunks recomputed by the
  independent CPU code all match. Together with Fischer's search to 2e14, this gives **no further
  terms below 3.08*10^14**.
* Bases 3, 5, 7: not searched yet (Phase 2). A partial range [7e14, 7.38e14) for the 16 bases is
  done but is not contiguous with the rest yet, so it supports no claim.

### 3.1 Ready now (optional): A306256 bound

Replace `No more terms up to 9.8*10^13.` by

    No further terms below 3.08*10^14. - _Jeff Sponaugle_, Oct 03 2026

(or append "There are no further terms below 3.08*10^14." to your existing a(4) comment and delete
the old line). Recommendation: make this one edit later together with the bound updates in 4, when
the search reaches 7e14 or 1.2e15, rather than an edit for a 1% gain now.

## 4. After the search resumes: bound updates (recommended at 1.2e15 = end of Phase 1)

One edit per entry: replace the stale bound comment by

    No further terms below X: all primes below 2*10^14 were searched by R. Fischer (see link) and
    those from 2*10^14 to X by an exhaustive search. - _Jeff Sponaugle_, <date>

with X = 3.08*10^14 today, 7*10^14 when the lower Phase-1 range finishes, 1.2*10^15 at the end of
Phase 1, and 3.2*10^15 at the end of Phase 2. For base 30 (A306256) the search below 2e14 is our own,
so the sentence about Fischer can be dropped there.

| entry | base | current bound comment | proposed |
|:--|--:|:--|:--|
| A212583 | 6 | "Next term > 4.119*10^13. [See Fischer link]" | text above |
| A045616 | 10 | "No further terms below 1.172*10^14 (as of Feb 2020, cf. Fischer's table)." | text above |
| A111027 | 12 | "Richard Fischer has carried this search to 4.8 * 10^13 (as of January 2014)." | text above |
| A128667 | 13 | "No further terms up to 3.127*10^13." | text above |
| A234810 | 14 | (none) | add text above |
| A242741 | 15 | "According to Richard Fischer there is no other term up to approximately 5*10^13." | text above |
| A128668 | 17 | "Mossinghoff showed that there are no further terms up to 10^14." | text above |
| A244260 | 18 | "According to Richard Fischer there is no other term up to approximately 5*10^13." | text above |
| A090968 | 19 | "No further terms up to 3.127*10^13." | text above |
| A242982 | 20 | "According to Richard Fischer, there is no other term up to approximately 5*10^13." | text above |
| A298951 | 22 | "Next term, if it exists, is larger than 8.72*10^13." | text above |
| A128669 | 23 | "No further terms up to 3.127*10^13." | text above |
| A306255 | 26 | "No more terms up to 9.8*10^13." | text above |
| A306256 | 30 | "No more terms up to 9.8*10^13." | "No further terms below X." (see 3.1) |
| A039951 | 47, 72 | "a(47) > 1.4*10^14, a(72) > 1.4*10^14 (see Fischer's tables)." | "a(47) > X, a(72) > X." (after Fischer's pending edit is approved) |

A014127 (3), A123692 (5) and A123693 (7) get the same treatment only after Phase 2 has searched
past 1.2e15 for them (their entries quote 9.7e14, 9.7e14 and 1.2e15).

Suggested text for the "Discussion" box of these edits:

    Bound from an exhaustive GPU search of all primes in [2*10^14, X) for bases 6, 10, 12, 13, 14,
    15, 17, 18, 19, 20, 22, 23, 26, 30, 47, 72 (DGX Spark, Montgomery arithmetic mod p^2). The prime
    count of the searched range equals pi(X) - pi(2*10^14) from primecount, every test re-checks
    Fermat's little theorem, and 1% of the range was recomputed by an independent CPU program with
    identical checksums. Below 2*10^14 the bound rests on R. Fischer's table (linked in the entry).

## 5. Not OEIS: tell Richard Fischer (optional)

His tables list base 30 with three solutions, and 303632117562967 is larger than every solution in
his record list for bases up to 1052 (the last is base 770, p = 194092950328757, Nov 2024). His address is
on fermatquotient.com. Draft for you to send:

    Dear Mr. Fischer,

    while searching for Wieferich primes to several bases with a GPU program, I found a fourth
    base-30 solution: p = 303632117562967 (30^(p-1) == 1 mod p^2, not mod p^3). It was confirmed by
    an independent exhaustive base-30 search of every prime below 3.0364e14 (which also confirms your
    result that there is none between 94727075783 and 2e14) and by independent CPU and Python checks.
    It is now a(4) of OEIS A306256. The search has covered all primes in [2e14, 3.08e14) for bases
    6, 10, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 26, 30, 47 and 72 so far, with no other solutions,
    and will continue to 1.2e15 and beyond. Thank you for your tables, which made it clear where to
    look.

    Best regards,
    Jeff Sponaugle
