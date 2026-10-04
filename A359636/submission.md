# A359636 — what can still be submitted (as of Oct 3 2026)

## Already in the entry (OEIS revision #31, Sep 18 2026)

* DATA now ends `..., 34710483181813, 31610535900218923` and EXTENSIONS has
  `a(9) from Jeff Sponaugle, Sep 18 2026`.  **a(9) is done.**
* COMMENTS has `a(10) > 31610555571634177. - Jeff Sponaugle, Sep 18 2026`.
* David A. Corneth's comment `a(9) <= 76340177205657727, a(10) <= 225096507194749219819`
  is still there unchanged (the a(9) half is now obsolete).
* No b-file was uploaded (the site synthesizes one from the 9 terms; fine), no
  EXAMPLE for a(9), no program link.

## 1. Improved lower bound for a(10)  — the main item

Replace the existing comment by

    a(10) > 408697968217030657. - Jeff Sponaugle, Oct 03 2026

Evidence: the GPU scan `gpu_n10.*` on atom1 (`/home/jbs/A359636/cuda`, state
frontier 31719137) tested every m == 3 (mod 6) with m <= 408697968217030659 for
omega(m-1), omega(m), omega(m+1) >= 10 and found none.  Every qualifying
level-10 gap contains such an m (the gap holds an odd multiple of 3 with both
neighbours inside it; for p <= M-2 that m is <= M since M == 3 (mod 6)), so no
prime p <= 408697968217030657 starts a qualifying gap.  Independently, the
level-9 scan had already shown the same below 3.16e16 (none of its 98 triples
has all three omegas >= 10), and the second machine's partial stretch
664309293225476097 <= m <= 740968712792834043 also contains no such m (it does
not extend the contiguous bound and is not claimed).

Option: finishing both halves of the paused scan (about 3.5 days on each Spark)
would give `a(10) > 999999999999999997`, i.e. a(10) > 10^18 - 3; resume commands
are in README.md (Status section).  Submit the bound above now, or wait for the
finished scan and submit the stronger one — editors prefer one edit over two.

## 2. EXAMPLE line for a(9)

    a(9) = 31610535900218923: 31610535900218924 = 2^2*11*17*19*53*79*149*359*9931, 31610535900218925 = 3*5^2*13*23*31*37*157*1559*5021 and 31610535900218926 = 2*7*29*61*103*197*223*311*907 each have 9 distinct prime factors, and 31610535900218927 is prime.

(Verified by `a359636 verify 9 31610535900218923`, by `verify_a359636.py`, and
by re-detection with the CPU sieve; see README.md.)

## 3. Optional comments

a) The reduction that makes the search tractable (useful to anyone extending
the sequence; keep or drop at the editors' taste):

    Every qualifying gap contains an odd multiple of 3, m, with m-1, m, m+1 all composites of the gap (for a gap of 4, m = p+2 since then p == 1 (mod 3)), so a(n) can be found by searching m == 3 (mod 6) with omega(m-1), omega(m), omega(m+1) >= n and then checking the surrounding prime gap. - Jeff Sponaugle, Oct 03 2026

b) Observation:

    The gaps at a(1)-a(9) all have length 4.

c) Housekeeping of the obsolete bound: suggest trimming Corneth's comment to

    a(10) <= 225096507194749219819. - David A. Corneth, Jan 12 2023

(or leaving it; editors sometimes keep superseded bounds for history).

## 4. Optional program link

If the tools are published (e.g. GitHub), a LINKS entry such as

    Jeff Sponaugle, <a href="...">C and CUDA programs</a> (three-target sieve over m == 3 (mod 6) with exact verification; reproduces a(1)-a(9)).

Nothing to submit for the cross-referenced sequences (A001359, A075590,
A185032); this work does not touch them.

## Files backing the submission

* `gpu_n9.txt/.log/.state` — the level-9 scan that found a(9) (98 triples, 1 solution).
* `gpu_n10.txt/.log/.state` (atom1) — level-10 coverage m <= 408697968217030659, 0 triples.
* `gpu_n10b.txt/.log/.state` (atom2) — level-10 coverage 664309293225476097..740968712792834043, 0 triples.
* `OEIS_draft.md` — the earlier a(9) draft (now submitted); `b359636.txt` — b-file for n = 1..9.
