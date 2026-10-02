# A306256 — draft OEIS submission text (Wieferich primes to base 30)

STATUS: FINAL (Oct 01 2026). a(4) = 303632117562967 was found by the multi-base GPU search
(Phase 1: 16 bases, all primes from 2e14 upward, DGX Spark) early on Oct 01 2026 and confirmed
by an independent exhaustive base-30 scan of every prime below 3.0364e14:

* that scan tested 9,396,616,950,003 primes, which equals pi(303640000000000) from primecount,
  in 30,364 contiguous chunks with no FLT errors; its only solutions are 7, 160541,
  94727075783 and 303632117562967 (so it also confirms R. Fischer's "none up to 2e14");
* 316 of its chunks (1%) were recomputed by the independent CPU code with identical prime
  counts and checksums;
* the solution itself checks out with the CPU tool (fast kernel and a shift-and-add
  reference), Python big integers with their own Miller–Rabin test, and OpenSSL primality, and a
  CPU rescan of its GPU chunk [303630000000000, 303640000000000) gives the same prime count
  (299887113), checksum (1829f95106870656) and solution.

## DATA

    7, 160541, 94727075783, 303632117562967

## COMMENTS (suggested)

1. Replace

       No more terms up to 9.8*10^13.

   by

       a(4) = 303632117562967 was found by an exhaustive search of all primes below
       3.04*10^14. No further terms below X. - Jeff Sponaugle, Oct 01 2026

   with X = the frontier of the running search at submission time (3.05*10^14 on Oct 01
   2026; the queued runs reach 7*10^14, 1.2*10^15 and then 3.2*10^15).

## EXAMPLE (optional)

    a(4) = 303632117562967: 30^(p-1) == 1 (mod p^2) but not (mod p^3); the Fermat quotient
    (30^(p-1) - 1)/p is divisible by p, and (30^(p-1) - 1)/p^2 == 180437496968356 (mod p).

## EXTENSIONS

    a(4) from Jeff Sponaugle, Oct 01 2026

## Notes for the submitter

* Records on atom2 in `/home/jbs/A001220/runs`: `v30.txt` / `v30.chunks` / `v30.dcheck`
  (independent base-30 scan of [0, 3.0364e14)); `p1a.txt` (the original FOUND / CONFIRMED
  lines), `p1a.chunks` (contiguous chunk log from 2e14 with per-base checksums; chunks before
  2.3982e14 were computed on atom1 before the checkpoint moved to atom2), `p1a.dcheck`.
* Re-verify: `./wieferich check 303632117562967 -b 30` and
  `python3 verify_wieferich.py check 303632117562967 30`.
* Context: under the 1/p heuristic the 16 Phase-1 bases had 0.20 expected solutions in
  [2e14, 3.05e14), so a find this early was an 18% event.
* It is also larger than every solution in R. Fischer's table for bases up to 1052 (his
  record list ends with base 770, p = 194092950328757, Nov 2024); he may want to know
  (address on fermatquotient.com).
* A039951(30) = 7 is unchanged.
