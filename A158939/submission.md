# A158939 — OEIS submission notes

Checked against the live entries on Oct 3 2026 (A158939 #38, A229832 #35, A133697 #27).

## Where the entries stand

Everything computed here is already in the OEIS. Jeff Sponaugle proposed the edits on
Sep 19 and Sep 22 2026 and each was approved the same day:

| entry | what went in | revisions |
|---|---|---|
| [A158939](https://oeis.org/A158939) | a(16) = 17293451238695141, a(17) = 52461866207504471, keyword `hard`, comment "a(18) > 5.30*10^16" | #27–#30 Sep 19 (approved by Alois P. Heinz); #31–#38 Sep 22 (edited by Michel Marcus, reviewed by Pontus von Brömssen, approved by Alois P. Heinz) |
| [A229832](https://oeis.org/A229832) | a(15) = 17293451238695143, a(16) = 52461866207504473 | #29–#31 Sep 19; #32–#35 Sep 22 (approved by Michael De Vlieger) |
| [A133697](https://oeis.org/A133697) | a(14) = 475618519121221, a(15) = 1400080864310974 | #21–#23 Sep 19; #24–#27 Sep 22 (approved by Michael De Vlieger) |

No draft is pending on any of the three. Their b-files are synthesized from the data
lines (17, 16 and 16 terms), so none was uploaded and none is needed at this length.

Correction to the local files: A158939 has offset 1 and starts at a(1) = 3, and always
has (archived copies from 2015, 2021 and 2025 agree). The "a(0) = 7" in the earlier
README table, OEIS_notes.md and b-file here was a misreading; b158939.txt is now n = 1..17.

Related entries, untouched by this work: A158940 (decreasing gaps, #7 of 2022,
"a(15) > 1.3*10^14" since 2016), A068843 (nondecreasing gaps, #26 of Dec 2024),
A348927 and A349121 (gaps in arithmetic progression, Nov 2021).

## What can still be submitted

Items 1 and 2 are small edits worth doing now; 3, 4 and 5 are optional; 6 should
wait; 7 and 8 need new runs (not started — your call, as usual).

### 1. Program link in A158939 (recommended)

The C program, the CUDA program, the Python verifier, the b-file and these notes are
public in the repository (checked: the directory URL resolves):

```
%H A158939 Jeff Sponaugle, <a href="https://github.com/jeffsponaugle/Mathpub-/tree/main/A158939">C and CUDA programs, verifier and notes</a>
```

Alternatively upload `a158939.c` as an a-file ("C program"); the editors assign the
file name. The PARI program in the entry (Charles R Greathouse IV) stays — it is the
definition in executable form.

### 2. Example for a(17) in A158939 (recommended)

The entry gives examples for a(8) and a(14) in this exact form. The gaps below were
re-derived from the definition by three independent programs (`a158939 verify` with
deterministic Miller-Rabin and with GMP, and `verify_run.py` in pure Python):

```
%e A158939 a(17)=52461866207504471 is the first prime to be followed by n=17 monotonic increasing prime gaps: 2,4,6,8,12,18,22,30,36,38,40,44,46,50,54,58,84.
```

If a second one is wanted:

```
%e A158939 a(16)=17293451238695141 is the first prime to be followed by n=16 monotonic increasing prime gaps: 2,4,6,8,22,26,30,34,36,44,46,54,56,64,86,108.
```

### 3. Density comment in A158939 (optional)

Exact counts from the scan of [0, 2*10^16) (GPU plus the CPU-scanned sliver; the total
prime count equals pi(2*10^16) = 547863431950008):

| exactly n increasing gaps | 13 | 14 | 15 | 16 |
|--:|--:|--:|--:|--:|
| primes below 2*10^16 | 23856 | 1467 | 84 | 1 |

```
%C A158939 Below 2*10^16 there are 23856 primes followed by exactly 13 monotonic increasing gaps, 1467 by exactly 14, 84 by exactly 15 and one (a(16)) by exactly 16. The counts fall by a factor of about 16 to 18 per step in n, which puts a(18) near 2*10^18. - _Jeff Sponaugle_, Oct 03 2026
```

This is the kind of comment that was accepted in A252768; it tells readers why the
entry is `hard` and where the next term is expected. Drop the last sentence if an
editor objects to the heuristic.

### 4. A229832: existence is now a theorem (optional)

A229832 still carries Sondow's 2013 comment "I make the stronger conjecture that the
sequence a(n) is infinite." Since a(n) is the prime after A158939(n+1) (Chris Boyd's
comment) and Banks, Freiberg and Turnage-Butterbaugh prove that A158939(n) exists for
every n (comment in A158939), A229832(n) exists for every n:

```
%C A229832 a(n) exists for every n, because A158939(n+1) does (Banks, Freiberg, & Turnage-Butterbaugh; see A158939). - _Jeff Sponaugle_, Oct 03 2026
%H A229832 William D. Banks, Tristan Freiberg, and Caroline L. Turnage-Butterbaugh, <a href="https://arxiv.org/abs/1311.7003">Consecutive primes in tuples</a>, arXiv:1311.7003 [math.NT], 2013-2014.
```

Adding keyword `hard` to A229832 and A133697 would match A158939 (which got `hard`
on Sep 22 at your suggestion); editors may or may not care.

### 5. A133697 cleanup (optional)

The 2008 comment "a(9) > 120000000. - Robert G. Wilson v, Mar 01 2008" is obsolete
(a(9) = 320620306 has been in the data since 2021) and can be deleted in the same edit
as the `hard` keyword.

### 6. a(18) bound — wait

The entry says a(18) > 5.30*10^16 (the Sep 22 stop at 53027239660093440). The paused
checkpoint on atom2 is at 54557804097699840 (5.46*10^16), only 3% further: not worth
an edit. Update the comment when the scan passes 10^17 (about 3.5 days of one Spark
from the checkpoint at 1.5e11 numbers/s) or finds a(18). Resume command:

```
cd /home/jbs/A158939/cuda && setsid nohup ./a158939_cuda scan 2e18 -Q 1e6 -S gpu_a18.state -i 60 >> gpu_a18.txt 2>> gpu_a18.log < /dev/null &
```

When a(18) appears: `./a158939 verify P`, `./verify_run.py P 18`, `primecount P` for
pi(a(18)) = A133697(16), and the next prime after a(18) is A229832(17); submit all three
entries as on Sep 22. Expected location: median about 2*10^18, 90% point about 8*10^18
(about five months of one Spark; splits across machines with `scan START END`).

### 7. A158940 — decreasing gaps (needs a run)

"First primes followed by sequences of exactly n monotonic decreasing prime gaps",
known to a(14) = 23840790158827 with "a(15) > 1.3*10^14" (Giovanni Resta, 2016) and no
progress since. Our tools track increasing runs only; mirroring the bookkeeping for
decreasing runs is a small change in both `a158939.c` and the CUDA extraction kernel
(same ring-buffer logic with the comparison reversed; the register budget of the
extraction kernel is the only thing to watch). If decreasing runs are about as dense as
increasing ones, a(15) sits in [1.3*10^14, 10^15] and a(16) near 10^16 to 10^17: a
15-hour GPU pass to 10^16 would very likely settle a(15) and has a fair chance at
a(16); another week reaches 10^17. Verification would be the same chain as here.

### 8. A068843 — nondecreasing gaps (needs a run, cheap)

"Smallest prime in the first occurrence of a nondecreasing difference for a set of
exactly n successive primes", known to a(15) = 2552790756469 (Resta, 2017), keyword
`hard`, no bound for a(16) in the entry. Note the different "exactly": the run must be
maximal on both sides (the gap before the run must be larger, see the entry's
examples). Ties make these runs more common, so the terms grow by a factor of only
5 to 20 per step; a(16) should be near 10^13 to 10^14 and a(17) near 10^15: an hour on
the Mac or minutes on the Spark once the tracker is added (one extra comparison and
the preceding-gap check).

## Verification chain behind the submitted terms

* Run lengths 16 at 17293451238695141 and 17 at 52461866207504471 re-derived from the
  definition by three independent implementations: `a158939 verify` (deterministic
  Miller-Rabin next-prime search), GMP `mpz_nextprime`, `verify_run.py` (pure Python).
* Minimality: GPU scan of [0, 10^16) (prime count = pi(10^16) = 279238341033925
  exactly), CPU scan of the sliver [10^16, 10^16 + 41615360) (1129839 primes, no run
  longer than 9), GPU scan of [10^16 + 41615360, 2*10^16) (total count = pi(2*10^16)
  = 547863431950008 exactly), continued GPU scan of [2*10^16, 54557804097699840).
* pi(a(16)) = 475618519121221 and pi(a(17)) = 1400080864310974 from primecount 8.7;
  the GPU's running prime index plus the sliver's 1129839 primes gives the same values.
* Both tools reproduce a(1)..a(15) with pi(a(n)) = A133697(n-2) for n = 2..15, and
  agree exactly (prime counts, first occurrences, run-length histograms) on [0, 10^13),
  on a window around a(15), and on windows at 3*10^16 and 10^18.
* The 15 consecutive weak primes of A229832(15): 17293451238695143, ...147, ...153,
  ...161, ...183, ...209, ...239, ...273, ...309, ...353, ...399, ...453, ...509,
  ...573, ...659; A229832(16) = 52461866207504473 starts the corresponding 16.
