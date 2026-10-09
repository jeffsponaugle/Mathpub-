# OEIS submissions: digit-content counts of the first 10^n primes

Status as of **2026-10-03**, cross-checked against the live OEIS entries
(fetched via `https://oeis.org/search?q=id:<A-number>&fmt=text`) and against
this project's computation (`digitprimes_results.txt`, independently verified
by the 278 checks in `verification.txt`).

## 1. Already submitted and published ✓

The term extensions computed by this project's a(13) run (first 10^13 primes,
all p ≤ p(10^13) = 323,780,508,946,331) **are already live on OEIS**, credited
in each entry's Extensions section:

- "no digit d" family: `a(13) from _Jeff Sponaugle_, Sep 20 2026`
- "at least one digit d" family: `a(12)-a(13) from _Jeff Sponaugle_, Sep 20-22 2026`

Nothing from the a(13) run remains unsubmitted. (b-files are not needed:
with ≤ 14 terms, OEIS synthesizes the "Table of n, a(n)" file automatically.)

## 2. Cross-check of every published term vs. this computation

Every term currently published in all 20 sequences was compared against this
project's run. **All 270 published terms match; no sequence has terms we lack.**

| Sequence | Definition | Offset | Published terms | Last published term | Cross-check |
|---|---|---|---|---|---|
| A231412 | no digit 0 | 0 | 14 (a(0)–a(13)) | 2,572,084,882,736 | all terms match |
| A231726 | ≥1 digit 0 | 1 | 13 (a(1)–a(13)) | 7,427,915,117,264 | all terms match |
| A228413 | no digit 1 | 0 | 14 (a(0)–a(13)) | 1,266,309,423,393 | all terms match |
| A231787 | ≥1 digit 1 | 1 | 13 (a(1)–a(13)) | 8,733,690,576,607 | all terms match |
| A228414 | no digit 2 | 0 | 14 (a(0)–a(13)) | 1,761,252,701,672 | all terms match |
| A231788 | ≥1 digit 2 | 1 | 13 (a(1)–a(13)) | 8,238,747,298,328 | all terms match |
| A228415 | no digit 3 | 0 | 14 (a(0)–a(13)) | 1,770,118,421,758 | all terms match |
| A231789 | ≥1 digit 3 | 1 | 13 (a(1)–a(13)) | 8,229,881,578,242 | all terms match |
| A228416 | no digit 4 | 0 | 14 (a(0)–a(13)) | 2,565,954,604,082 | all terms match |
| A231790 | ≥1 digit 4 | 1 | 13 (a(1)–a(13)) | 7,434,045,395,918 | all terms match |
| A228417 | no digit 5 | 0 | 14 (a(0)–a(13)) | 2,567,010,423,857 | all terms match |
| A231792 | ≥1 digit 5 | 1 | 13 (a(1)–a(13)) | 7,432,989,576,143 | all terms match |
| A228418 | no digit 6 | 0 | 14 (a(0)–a(13)) | 2,567,938,055,999 | all terms match |
| A231793 | ≥1 digit 6 | 1 | 13 (a(1)–a(13)) | 7,432,061,944,001 | all terms match |
| A228419 | no digit 7 | 0 | 14 (a(0)–a(13)) | 1,926,738,759,318 | all terms match |
| A231794 | ≥1 digit 7 | 1 | 13 (a(1)–a(13)) | 8,073,261,240,682 | all terms match |
| A228420 | no digit 8 | **1** | 14 (a(0)–a(13)) | 2,570,658,222,631 | all terms match |
| A231795 | ≥1 digit 8 | 1 | 13 (a(1)–a(13)) | 7,429,341,777,369 | all terms match |
| A228421 | no digit 9 | 0 | 14 (a(0)–a(13)) | 1,928,516,761,196 | all terms match |
| A231796 | ≥1 digit 9 | 1 | 13 (a(1)–a(13)) | 8,071,483,238,804 | all terms match |

Endpoint anchor: the run's crossing primes p(10^0)–p(10^13) all match
OEIS A006988 (the 10^n-th prime), independently published.

## 3. Remaining submission opportunities

### 3.1 New terms: a(14) for all 20 sequences  — the main one

Each "no digit d" entry still carries keywords `more,hard`: more terms are
wanted. a(14) means digit statistics over the first 10^14 primes, i.e. all
primes up to p(10^14) = 3,475,385,758,524,527 (A006988(14)) — about 10.7× the
a(13) sieving volume.

- Command: `./digitprimes --target 14` (resumes from the existing
  `digitprimes.ckpt`, which already sits at the 10^13 mark, so only the new
  decade is computed).
- Estimated runtime on this M4 Max: roughly 30–40 hours (throughput declines
  slowly with height; a(13) averaged ~1.15 G primes/s).
- The run self-validates: it must reproduce every a(n ≤ 13) term and hit
  p(10^14) exactly, and `gen_verification.py` can be extended with A006988(14)
  for the endpoint check.
- Submits as: one term appended to each of the 20 sequences
  (a(14) on the "no d" family, a(14) on the "≥1 d" family), with an
  Extensions line `a(14) from _Jeff Sponaugle_, <date>`.

### 3.2 Exact complement formula on the "at least one d" entries

The ten "≥1 digit d" entries currently have only the asymptotic
`a(n) ~ 10^n` (Greathouse, 2014). The exact identity is submittable as a
%F Formula line on each:

| Entry | Formula to add |
|---|---|
| A231726 | a(n) = 10^n − A231412(n) |
| A231787 | a(n) = 10^n − A228413(n) |
| A231788 | a(n) = 10^n − A228414(n) |
| A231789 | a(n) = 10^n − A228415(n) |
| A231790 | a(n) = 10^n − A228416(n) |
| A231792 | a(n) = 10^n − A228417(n) |
| A231793 | a(n) = 10^n − A228418(n) |
| A231794 | a(n) = 10^n − A228419(n) |
| A231795 | a(n) = 10^n − A228420(n) |
| A231796 | a(n) = 10^n − A228421(n) |

This is how the a(12)–a(13) values were in fact derived, and it makes future
extensions of the "≥1" family automatic whenever the "no d" family is
extended. A matching named crossref (e.g. "Cf. A228413 (complement count)")
would tighten the currently generic %Y lists.

### 3.3 Offset correction: A228420

A228420's offset reads `1,2`, but its terms plainly begin at n = 0
(a(0) = 1: the single first prime, 2, contains no 8), and all nine sibling
"no digit d" entries use offset `0,2`. Submittable as an offset edit to
`0,2` with that one-line justification. (The second component, 2 = position
of the first term exceeding 1, is already correct.)

### 3.4 Keyword `hard` on the "≥1 digit d" family

The "no digit d" entries carry `hard`; the ten "≥1 digit d" entries do not,
yet each of their new terms requires exactly the same 10^n-prime enumeration
(or an extension of the complement sequence). Proposing `hard` on
A231726, A231787–A231790, A231792–A231796 is a defensible editorial edit —
reviewers may instead view them as trivially derived from the "no d" family,
so expect discussion.

### 3.5 Program contributions (optional)

All 20 entries already have a Mathematica program (%t); only A228421 has an
additional %o. The C++ tool here is too long for an inline program section,
but if the repo is published (e.g. GitHub), a %H Links line such as
`Jeff Sponaugle, <a href="...">digitprimes</a> (C++ program using primesieve)`
on each entry would document how the large terms were obtained — the same
pattern Lucas A. Brown used for the a(12) batch.

### Not worth submitting

- **b-files**: auto-synthesized for sequences this short; only useful when a
  sequence has far more terms than fit in the %S/%T/%U lines.
- **a(14) for A006988 itself**: p(10^14) and beyond are already published.

## 4. Submission mechanics

1. Log in (or register) at https://oeis.org — contributions are tied to a
   registered user name (the existing extensions are credited to
   _Jeff Sponaugle_).
2. On each sequence page choose "edit", add the new material in the draft
   (terms to the Data section, %F/%E/%H lines as appropriate), and give a
   short justification in the pink-box comment, e.g.: "a(14) computed by
   exhaustive enumeration of the first 10^14 primes with a primesieve-based
   C++ program; run reproduces all existing terms a(0)-a(13) and the
   endpoint p(10^14) = A006988(14)."
3. Each draft goes through editor review; batch the 20 edits over a few days
   (as was done for a(13)) to keep the review load reasonable.
4. Keep `verification.txt` regenerated (`python3 gen_verification.py`) after
   any new run — editors appreciate the cross-checks being stated.
