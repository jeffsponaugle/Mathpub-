# OEIS submission package — postage stamp problem family (A001208 and cross-references)

Prepared 2026-10-03 by Jeff Sponaugle (with Claude). Everything below is verified and ready to enter; nothing has been
submitted yet. Each block is written in OEIS field style so it can be pasted into the edit form. One sequence per
submission; sign new lines "_Jeff Sponaugle_, Oct 03 2026". Supporting files are in `bfiles/`, logs in `logs/`.

## 0. OEIS status check (2026-10-03, via oeis.org JSON + revision histories)

No data or comment changes affect this package since the analysis of 2026-09-30:

| entry | revision | last edit | terms | note |
|---|---|---|---|---|
| A001211 | 42 | 2026-05-30 (S. A. Irvine, batch) | 25 | a(20) still 45745 (typo); a(26) absent |
| A001210 | 44 | 2026-05-30 (batch) | 67 (b-file R. Price) | a(68..90) absent |
| A001209 | 50 | 2026-05-30 (batch) | 54 (b-file R. Price) | comment still says a(29)..a(254) computable |
| A053346 / A053348 | 25 / 25 | 2025-12-01 (A. Zabolotskii) | 13 / 8 | a(14) absent / a(9) absent |
| A005344 / A005342 / A005343 / A075060 | 29 / 30 / 30 / 20 | 2026-05-30, 2025-11-05 | 7 / 8 / 8 / 7 | unchanged |
| A001213 / A001214 / A001215 / A001216 | 38 / 46 / 46 / 38 | 2026-05-30 (batch) | 25 / 25 / 20 / 10 | contest values still shown as exact; no caveat comment |
| A001212 | 112 | 2026-08-14 (M. De Vlieger: archived Friedman link) | 24 | unchanged |
| A084192 / A084193 / A196416 | 21 / 27 / 14 | 2019 / 2021 / 2017 | 71 / 71 / 66 | unchanged |

The 2026-05-30 edits are a same-second batch across the whole family (link housekeeping); the 2025-12-01 edits on
A053346/A053348 predate our analysis and are already reflected in README.md.

---

## 1. NEW TERM: A053348(9) = 6082  (8 denominations, 9 stamps)   — highest priority

Exhaustive search completed 2026-10-02 (audit `logs/n9k8_final_audit.txt`: 11003/11003 depth-4 work items, 479 reported
bases). Details and verification notes: `bfiles/A053348_A005344_submission.md`.

```
%S A053348 8,32,93,228,524,1007,1911,3485,6082
%H A053348 Jeff Sponaugle, <a href="/A053348/b053348.txt">Table of n, a(n) for n = 1..9</a>
%e A053348 a(9) = 6082: the unique extremal 9-basis with 8 denominations is {1, 8, 27, 88, 197, 521, 1226, 1461}: every
           integer 1..6082 is a sum of at most 9 of these denominations and 6083 is not. - _Jeff Sponaugle_, Oct 03 2026
%C A053348 a(9) was found by an exhaustive search over all admissible bases (Challis's H-program method with the
           Challis-Robinson "difficult target" test and an exact final check), using the proven values of A005344(1..7)
           as bounds; about 1.5*10^14 leaf candidates were examined. - _Jeff Sponaugle_, Oct 03 2026
%E A053348 a(9) from _Jeff Sponaugle_, Oct 03 2026
```
b-file: `bfiles/bA053348.txt` (n = 1..9).
Uniqueness of the basis: established (2026-10-03). The ~18 % of items first searched only with target 6083 were re-run at
target 6082 on the Jetson (8 chunks, 0 bases); all other items were searched with targets <= 6082 and reported exactly one
basis of range 6082. The word "unique" in the example line is therefore justified.

## 2. NEW TERM: A005344(8) = 6082  (n denominations, 9 stamps) — same result, transposed sequence

```
%S A005344 9,34,112,326,797,1617,3191,6082
%H A005344 Jeff Sponaugle, <a href="/A005344/b005344.txt">Table of n, a(n) for n = 1..8</a>
%e A005344 a(8) = 6082 = A053348(9), attained only by {1, 8, 27, 88, 197, 521, 1226, 1461}. - _Jeff Sponaugle_, Oct 03 2026
%E A005344 a(8) from _Jeff Sponaugle_, Oct 03 2026
```
b-file: `bfiles/bA005344.txt`. (The entry's only link is Friedman's page; the Challis & Robinson JIS 2010 link used by
the sibling entries could be added at the same time.)

## 3. NEW TERM + CORRECTION + LITERATURE TERM: A001211 (6 denominations, n stamps)

Three changes, best made in one edit:
* **a(20): 45745 -> 45754** (OEIS typo). Challis & Robinson's table gives n(20,6) = 45754; the basis {1, 17, 93, 436, 2898, 6897}
  has 20-range 45754 (verified by direct computation), so 45745 cannot be the maximum.
* **a(26) = 156744** from the Challis & Robinson July 2013 addendum (never entered; h-range re-verified).
* **a(27) = 186942**, new: exhaustive search completed 2026-10-03 (audit `logs/n27k6_final_audit.txt`: 30/30 chunks of 100000
  depth-4 work items, boundary windows checked, 1147 reported bases; 1.1*10^14 leaf candidates, 7.4*10^10 full checks,
  3326 core-hours on three Xeon servers). Details: `bfiles/A001211_27_submission.md`.

```
%S A001211 6,20,52,108,211,388,664,1045,1617,2510,3607,5118,7066,9748,12793,17061,22342,28874,36560,45754,
%T A001211 57814,72997,87555,106888,129783,156744,186942
%H A001211 Jeff Sponaugle, <a href="/A001211/b001211.txt">Table of n, a(n) for n = 1..27</a>
%H A001211 M. F. Challis and J. P. Robinson, <a href="https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6-new.pdf">
           Addendum (July 2013) to "Some Extremal Postage Stamp Bases"</a>
%e A001211 a(20) = 45754 is attained by {1, 17, 93, 436, 2898, 6897}. a(27) = 186942: the extremal 27-basis with 6 denominations
           is {1, 19, 194, 1095, 7370, 27669}; every integer 1..186942 is a sum of at most 27 of these denominations and 186943
           is not. - _Jeff Sponaugle_, Oct 03 2026
%C A001211 a(27) was found by an exhaustive search (Challis's H-program method with the Challis-Robinson "difficult target"
           test and an exact final check), using the proven values n(27,k), k <= 5 (A014616, A001208, A001209, A001210) as
           bounds. - _Jeff Sponaugle_, Oct 03 2026
%E A001211 a(20) corrected (45745 was a typo for Challis and Robinson's 45754), a(26) from the Challis-Robinson 2013 addendum,
           and a(27) from _Jeff Sponaugle_, Oct 03 2026
```
b-file: `bfiles/bA001211.txt` (n = 1..27). Uniqueness of the a(27) basis: certain for chunks searched with target <= 186942
(0-13, 15-21); chunks 14, 22-29 ran with target 186943 and could hide a second basis of range exactly 186942 (~15 h on a
96-thread box to settle; not needed for the entry).

## 4. LITERATURE TERMS: A001210 a(68)..a(90)  (5 denominations)

From the Challis & Robinson July 2013 addendum; every basis re-verified by direct h-range computation (`tools/hrange`).
Values: 2330896, 2496702, 2653201, 2846834, 3047485, 3250580, 3429203, 3629795, 3864527, 4103963, 4416370, 4643287,
4975426, 5223883, 5519971, 5796515, 6139689, 6513282, 6912409, 7258582, 7677138, 8029729, 8525267.

```
%H A001210 Jeff Sponaugle, <a href="/A001210/b001210.txt">Table of n, a(n) for n = 1..90</a>   (replaces R. Price's 1..67)
%H A001210 M. F. Challis and J. P. Robinson, <a href="https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6-new.pdf">
           Addendum (July 2013) to "Some Extremal Postage Stamp Bases"</a>
%E A001210 a(68)-a(90) from the July 2013 addendum to Challis and Robinson, added by _Jeff Sponaugle_, Oct 03 2026
```
(The data line itself stays at its current length; the terms go in the b-file.) b-file: `bfiles/bA001210.txt`.

## 5. LITERATURE TERM: A053346(14) = 24466  (7 denominations)

```
%S A053346 7,26,70,162,336,638,1137,2001,3191,5047,7820,11568,17178,24466
%H A053346 Jeff Sponaugle, <a href="/A053346/b053346.txt">Table of n, a(n) for n = 1..14</a>
%H A053346 M. F. Challis and J. P. Robinson, <a href="https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6-new.pdf">
           Addendum (July 2013) to "Some Extremal Postage Stamp Bases"</a>
%E A053346 a(14) from the July 2013 addendum to Challis and Robinson, added by _Jeff Sponaugle_, Oct 03 2026
```
b-file: `bfiles/bA053346.txt`.

## 6. B-FILE EXTENSION: A001209 to n = 302  (4 denominations)

a(55)..a(302) follow from the three formula families (types A, B, C) with the coefficient table of the 2013 addendum;
all 302 formula bases were verified to have the stated h-range (`bfiles/SUMMARY.txt`: 0 mismatches), and the formulas
reproduce the 1993/2010 tables wherever they overlap. Robinson's existing comment mentions only a(29)..a(254).

```
%H A001209 Jeff Sponaugle, <a href="/A001209/b001209.txt">Table of n, a(n) for n = 1..302</a>   (replaces R. Price's 1..54)
%H A001209 M. F. Challis and J. P. Robinson, <a href="https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6-new.pdf">
           Addendum (July 2013) to "Some Extremal Postage Stamp Bases"</a>
%C A001209 a(55)-a(302) are given by the three formula families (types A, B, C) with the coefficient table in the July 2013
           addendum to Challis and Robinson; each formula basis was verified by direct computation to have the stated
           range. - _Jeff Sponaugle_, Oct 03 2026
%E A001209 b-file extended to n = 302 using the Challis-Robinson formulas by _Jeff Sponaugle_, Oct 03 2026
```
b-file: `bfiles/bA001209.txt`. Optional lower-bound comment (not proven, so only as a comment): "a(303) >= 71148327,
a(304) >= 72060639, a(305) >= 72972951, a(306) >= 73885263 (the formula basis {1, 228, 17657, 912312} continued; ranges
verified)". Our exhaustive proof of a(303) was started but not finished (resumable, see DEPLOY.md).

## 7. ARRAY B-FILES: A084192 (to 109 terms), A084193 (113), A196416 (145)

Extended from proven values only (stop at the first unknown antidiagonal entry: (h,k) = (5,11) for A084192, (7,9) for
A084193/A196416); the existing OEIS data is reproduced exactly as a prefix (`bfiles/SUMMARY.txt`). The new entries n(9,8)
and n(27,6) lie beyond these antidiagonals and do not enter the b-files.

```
%H A084192 Jeff Sponaugle, <a href="/A084192/b084192.txt">Table of n, a(n) for n = 0..108</a>
%H A084193 Jeff Sponaugle, <a href="/A084193/b084193.txt">Table of n, a(n) for n = 0..112</a>
%H A196416 Jeff Sponaugle, <a href="/A196416/b196416.txt">Table of n, a(n) for n = 0..144</a>
```
(Check each file's index range against the entry's offset before uploading: A084192/A084193 have offset 0, A196416 offset 0.)
Files: `bfiles/bA084192.txt`, `bfiles/bA084193.txt`, `bfiles/bA196416.txt`.

## 8. CAVEAT COMMENTS: A001213, A001214, A001215, A001216 carry unproven contest values

Friedman's table (source of these terms, entered Jul 2013) was updated the day after Al Zimmermann's "Son of Darts"
contest ended (2010-06-20); the entries for 3 stamps with 16..25 denominations, 4 stamps with 13..25, 5 stamps with
11..20 and 6 stamps with 10 are the contest's best scores minus one. Challis & Robinson's exhaustive results stop at
(3,15), (4,12), (5,10), (6,9) and nothing published since proves a further term. All 111 contest bases were re-verified
to attain their values (`bfiles/sonofdarts_best.txt`), so the terms are correct as lower bounds.

Proposed comment (adapt the indices per entry):
```
%C A001213 Terms a(16) onward are the best values found in Al Zimmermann's "Son of Darts" programming contest (2010), as
           tabulated by Friedman; they have not been proved optimal and should be regarded as lower bounds. Exhaustive
           results are known only for n <= 15 (Challis and Robinson 2010, addendum 2013). - _Jeff Sponaugle_, Oct 03 2026
%C A001214 ... a(13) onward ... exhaustive results only for n <= 12 ...
%C A001215 ... a(11) onward ... exhaustive results only for n <= 10 ...
%C A001216 ... a(10) ... exhaustive results only for n <= 9 ...
```
Editors may instead truncate the data and keep the contest values as "a(16) >= 450, ..." comments; `bfiles/sonofdarts_best.txt`
also lists the contest's best-known values beyond the OEIS data (3 stamps to k = 40, 4 stamps to k = 30, 5 stamps to k = 20)
if lower-bound comments are wanted. Add the link
`Al Zimmermann's Programming Contests, <a href="http://azspcs.com/Contest/SonOfDarts">Son of Darts</a> (2010)` where used.

## 9. OPTIONAL lower-bound comments (only if an editor wants them)

* A001212: a(25) >= 228 and a(26) >= 244 (Kohonen's extremal restricted 2-bases, A006638); a simulated-annealing search
  over unrestricted bases found nothing higher (this work). No proof.
* A001210: a(91) >= 8897042 (the h = 90 extremal basis {1, 70, 1412, 32085, 371775} has 91-range 8897042). Weak bound; skip.

## 10. Suggested submission order

1. A053348 (new term) and A005344 (same term) — one draft each, b-files attached.
2. A001211 — correction + a(26) + a(27) + b-file, one draft (explain the a(20) typo in the edit summary with the basis).
3. A001210, A053346, A001209 — literature b-files (cite the addendum; mention the independent verification).
4. A084192 / A084193 / A196416 — array b-files.
5. A001213-A001216 caveat comments (expect discussion; offer the truncation alternative).
Keep the drafts small; OEIS editors prefer one change set per sequence. The addendum URL
https://cs.uwaterloo.ca/journals/JIS/VOL13/Challis/challis6-new.pdf was verified reachable (HTTP 200) on 2026-10-03.
B-file index ranges verified: A084192 n = 0..108, A084193 n = 0..112, A196416 n = 0..144, A001210 n = 1..90, A001209 n = 1..302.
