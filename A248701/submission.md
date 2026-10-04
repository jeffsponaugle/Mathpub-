# A248701–A248704 — OEIS submission notes

Checked against the live entries on Oct 3 2026.

## Where the entries stand

| entry | revision | in the OEIS now | still missing |
|---|---|---|---|
| [A248701](https://oeis.org/A248701) | #47, Sep 18 2026 | a(8)–a(10), clarifying comment, b-file n = 1..10, `%E a(8)-a(10) from Jeff Sponaugle, Sep 17 2026` | the bound in the comment is stale ("a(11) > 6.02*10^13"; proven ≥ 3.86*10^15) |
| [A248702](https://oeis.org/A248702) | #25, Sep 20 2026 | a(7)–a(11), b-file n = 0..11, two `%E` lines | the obsolete comment "a(7) >= 8960453, if it exists" (R. J. Mathar, 2014) is still there; no bound for a(12) |
| [A248703](https://oeis.org/A248703) | #25, Sep 20 2026 | **offset changed to 2 by the editors** (a(n) has n gaps on each side, a(3) = 1439), a(8)–a(10) in the new indexing, b-file n = 2..10 | the `%E` lines still say "a(7)-a(8)" and "a(9)" (old indices); no bound for a(11) |
| [A248704](https://oeis.org/A248704) | #19, Sep 17 2026 | a(7)–a(9), b-file n = 1..9, `%E a(7)-a(9)` | **a(10) = 1958854030679863 (found Sep 19) is not submitted**; no bound for a(11) |

Nobody else has touched the four entries since your edits, and the A248703
offset change is the only editorial change of substance. The tool, the Python
verifier, `README.md` and `b248703.txt` in this directory now use the new
A248703 indexing; the older output files (`scan_1e13.txt`,
`scan_1e15_summary.txt`, `scan_1e17_records.txt`) label A248703 terms with the
old index, one less.

## Before submitting: pin down the current bound

Every bound below is written as 3.86*10^15, the frontier the a(11) hunt had
reached on Sep 19 (all primes below 3860678242735729 examined). The run on the
Mac Studio has moved on since; use the real frontier, rounded *down* to three
significant digits:

* if the run finished (it stops by itself at a(11) or at 10^17):
  `grep '^# scanned' scan_1e17.txt` gives `scanned [START, X]`; the bound is X;
* if it is still running or was interrupted: `grep '^lo ' a248701.state` gives
  the first number not yet covered, so the bound is lo − 1 (chunks below `lo`
  are all folded; the file is rewritten every 60 s and on Ctrl-C);
* `tail -1 scan_1e17.log` shows the same position in the status line.

If `scan_1e17.txt` contains a line `A248701(11) = P`, skip to section 5 first.

## 1. A248704: a(10) = 1958854030679863  (recommended — the only unsubmitted term)

Found Sep 19 2026 by the exhaustive scan; every prime below it was examined.
Its window of 21 consecutive primes was re-verified with the pure-Python
Miller–Rabin check (`verify_py_windows.txt`): gaps 70, 56, 42, 40, 38, 30, 12,
10, 8, 6 before it (strictly decreasing) and 16, 20, 22, 26, 30, 54, 58, 74,
88, 168 after it (strictly increasing).

```
%S A248704 3,19,1429,25243,340577,1107791,3531448007,17190066197,37148264596189,1958854030679863
%C A248704 a(11) > 3.86*10^15. - _Jeff Sponaugle_, Oct 03 2026
%E A248704 a(10) from _Jeff Sponaugle_, Oct 03 2026
```

Upload `b248704.txt` from this directory (n = 1..10) to replace the current
b-file. Optional example line:

```
%e A248704 a(10) = 1958854030679863: the ten gaps before it are 70, 56, 42, 40, 38, 30, 12, 10, 8, 6 and the ten after it are 16, 20, 22, 26, 30, 54, 58, 74, 88, 168.
```

## 2. A248702: retire the obsolete comment  (recommended)

The comment "a(7) >= 8960453, if it exists. - _R. J. Mathar_, Dec 04 2014" is
superseded by the data (a(7) = 938665577 is in the entry). Replace it with the
current bound and say in the edit summary that the 2014 bound is now in the
data:

```
%C A248702 a(12) > 3.86*10^15. - _Jeff Sponaugle_, Oct 03 2026
```

Optional, the same clarification that A248701 already carries, adapted to the
valley (the name says "monotonically", which readers may take as strict):

```
%C A248702 "Monotonically decreasing" and "monotonically increasing" are meant in the weak sense: the n gaps before a(n) are nonincreasing and the n gaps after it are nondecreasing (a(9) = 299917793009 has gaps 80, 70, 44, 36, 36, 24, 10, 6, 6 before it and 8, 12, 12, 12, 16, 24, 24, 24, 26 after it). The two gaps adjacent to a(n) are not compared with each other: a(6) = 1107791 has gaps 4 and 2 around it, a(9) has 6 and 8. A248704 is the version with strict inequalities.
```

## 3. A248703: fix the extension indices, add the bound  (recommended)

With the new offset 2 the terms you added are a(8), a(9) (Sep 17) and a(10)
(Sep 20), but the extension lines still read "a(7)-a(8)" and "a(9)". Edit them
to

```
%E A248703 a(8)-a(9) from _Jeff Sponaugle_, Sep 17 2026
%E A248703 a(10) from _Jeff Sponaugle_, Sep 20 2026
```

and add the bound

```
%C A248703 a(11) > 3.86*10^15. - _Jeff Sponaugle_, Oct 03 2026
```

Optional comment spelling out the shape, since the name only says
"decreasing on either side":

```
%C A248703 a(n) is the smallest prime p such that the n prime gaps before p are strictly increasing and the n prime gaps after p are strictly decreasing; the two gaps adjacent to p are not compared with each other. A248701 is the version with weak inequalities.
```

## 4. A248701: raise the bound  (recommended), extras optional

Edit the last sentence of the existing comment rather than adding a new one
(it also lacks its full stop):

```
old:  a(11) > 6.02*10^13 - _Jeff Sponaugle_, Sep 17 2026
new:  a(11) > 3.86*10^15. - _Jeff Sponaugle_, Oct 03 2026
```

Optional example for the largest term:

```
%e A248701 a(10) = 59852066157421: the ten gaps before it are 8, 10, 24, 26, 34, 38, 40, 50, 90, 102 and the ten after it are 88, 68, 52, 30, 24, 20, 18, 18, 10, 8.
```

Optional program link, if this directory is pushed to the public repository
used for A252768 (check the URL first):

```
%H A248701 Jeff Sponaugle, <a href="https://github.com/jeffsponaugle/Mathpub-/tree/main/A248701">C program, Python verifier and notes</a>
```

The same link can go on the other three entries. Alternatively upload
`a248701.c` as an a-file ("C program for A248701-A248704"); the editors
assign the file name.

## 5. If a(11) turns up

Before submitting a new term:

1. `./a248701 show P` must report peak depth ≥ 11;
2. `python3 verify_a248701.py window P 11` must print the window and
   "checked with Miller-Rabin" (independent primality test of the 23 primes
   and of every number between them);
3. the record line is only printed once every earlier chunk is folded, so the
   term is the smallest; `# scanned [.., X]` in `scan_1e17.txt` gives the bound
   for a(12).

Then for A248701: append P to `%S`, add `%E A248701 a(11) from _Jeff
Sponaugle_, <date>`, upload a b-file n = 1..11 (`./a248701` output plus the
header of `b248701.txt`), and turn the bound sentence into "a(12) > X". Keep
`more`: a(12) needs another factor of roughly 120 in the search bound (around
10^18), i.e. months of sieving.

## 6. Not worth submitting

* The depth histograms and the heuristic growth model (count(≥ d)/count(≥ d+1)
  ≈ 0.85 d² for the weak shapes, ≈ 1.25 d² for the strict ones): README only.
* Sequences of all primes of a given depth (the full sets whose first terms
  are a(n)): plain intersections of gap conditions, low interest.
