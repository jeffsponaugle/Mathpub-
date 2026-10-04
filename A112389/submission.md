# OEIS submission draft — LEGO 2×4 entries (A123829, A123831, A112390; links in A112389)

Draft of 2026-10-03. Everything below was computed or checked with the tools in this
directory: `lego.c` (modes `-S`, `-r`, `-b`), `rows_a112390.py`, `sym_vs_lasse.py`. See
README.md for the methods and checks, and `results/` for the run outputs.

## Status on the OEIS (checked 2026-10-03, evening)

* **A112389** (all buildings). #69, May 11 2026: a(1..10), keywords `nonn,hard,more`.
  No pending edit and no draft. a(11) is still unknown, so there is nothing to add to DATA;
  section 4 has optional link fixes.
* **A112390** (triangle by height). #26, Jun 01 2024: rows 2..6. No pending edit.
  **T(5,2) is wrong** (section 1).
* **A123829** (180°-symmetric buildings). #9, Feb 15 2025: a(1..9). No pending edit.
  Section 2.
* **A123831** (180°-symmetric, on a fixed base brick). #9, Nov 29 2017: a(1..8), keyword
  `more`. No pending edit. Section 3.
* **A123830** (on a fixed base brick). #20, Aug 17 2018: a(1..9). a(10) needs F for the
  refinements <1541> and <13421>, which nobody has yet, so nothing to submit.
* A272690 (towers) has a closed form, so nothing to do.
* Deleuran's GitHub (github.com/LasseD/A112389): last commit 2026-09-18, and his table
  (`sum-for-size.py`) is unchanged.

## Summary

| item | status |
|---|---|
| A112390: T(5,2) = 248688 → **298688** | **to submit** (the current row 5 does not sum to A112389(5)) |
| A112390: rows 7–8 | **to submit** (new; recomputed here independently) |
| A112390: rows 9–10 | optional: taken from Deleuran's refinement table (consistent, see §1) |
| A123829: a(10) = 145036229, **a(11) = 2556043967** | **to submit** (a(10) is also in Deleuran's GitHub README; a(11) is new) |
| A123831: **a(9) = 34006742, a(10) = 471305645** | **to submit** (new) |
| A112389: link to Deleuran's project; third author on the arXiv link | optional |
| A112389 a(11), A123830 a(10) | not possible yet: 7 refinements of size 11 are still missing |

Next step: submit sections 1–3 as three separate edits; section 4 is optional. A123829(10)
and A112390 rows 9–10 are Deleuran's numbers, so consider a short note to him first
(section 5).

Before citing the program link, push the current `lego.c` (memo fix and `-S` mode),
README.md, the helper scripts and `results/`. `e6f20f4` has an older version of the
directory.

---

## 1. A112390 — correct T(5,2), add rows 7–8

Current entry: #26, rows 2..6, keywords `nonn,tabl,more`, offset `2,1`.

**The error.** Row 5 reads 248688, 3203175, 4425804, 2238736, which sums to 10116403. But
A112389(5) = 10166403. The height-2 buildings with 5 bricks have layer profiles <41>, <32>,
<23>, <14>, with 550, 148794, 148794 and 550 buildings, so T(5,2) = **298688**. The other
three entries of row 5 and rows 2, 3, 4, 6 are correct; they were recomputed here and agree
with Deleuran's table. The wrong value most likely came from a one-digit transcription
error (9 → 4).

### DATA (replace: rows 2..7, 21 terms, 196 characters)

```
24, 500, 1060, 11707, 59201, 48672, 298688, 3203175, 4425804, 2238736, 7946227, 162216127, 359949655, 282010252, 102981504, 217401554, 8179857079, 26945799212, 29086956358, 16580215072, 4737148480
```

### b-file (upload)

`b112390.txt`: rows 2..8 (28 terms, indices 2..29, following %O 2,1). Row 8:

```
6071726520, 413098084569, 1936986606505, 2730620544737, 2040743061784, 928646763600, 217908828672
```

Optional, instead: `results/b112390_rows2-10.txt` adds rows 9 and 10 (45 terms) from
Deleuran's table:

```
row 9:  172218479774, 20947504620559, 136237276258906, 243794643167537, 223089242301790, 132014359862752, 50351768747072, 10023806116096
row 10: 4945620426880, 1069899937649401, 9458503109455193, 21085020016366854, 22932624809115449, 16278809832673900, 8094551995005200, 2667346176625088, 461095081334784
```

Both rows sum to A112389(9) and A112389(10), and their last entries equal A272690(9) and
A272690(10). They have not been recomputed here. Row 9 would take an estimated hour with
`python3 rows_a112390.py results/sym11.txt 9 results/deleuran-sum-for-size.py`. Row 10 contains
refinements such as <55>, <541> and <532> that take this method many hours each. If rows 9–10
are included, add Deleuran's link (section 4) and the extra EXTENSIONS line below.

### EXAMPLE (replace)

```
Triangle begins:
         24;
        500,        1060;
      11707,       59201,       48672;
     298688,     3203175,     4425804,     2238736;
    7946227,   162216127,   359949655,   282010252,   102981504;
  217401554,  8179857079, 26945799212, 29086956358, 16580215072, 4737148480;
```

### COMMENTS (add)

```
Row sums give A112389. T(n,n) = A272690(n), the buildings of maximal height. - _Jeff Sponaugle_, Oct 03 2026
```

### EXTENSIONS (add)

```
T(5,2) corrected and rows 7-8 added by _Jeff Sponaugle_, Oct 03 2026
```

and, only if rows 9–10 are included:

```
Rows 9-10 from the refinement counts of Lasse Deleuran, Oct 03 2026
```

### Notes to the editors

```
T(5,2) was 248688, which makes row 5 sum to 10116403 instead of A112389(5) = 10166403. The height-2 buildings with 5 bricks have layer profiles <41>, <32>, <23>, <14> with 550, 148794, 148794, 550 buildings, so T(5,2) = 298688. Rows 2-8 were recomputed layer profile by layer profile with a C program: Redelmeier enumeration of the bottleneck-free profiles (with some layers counted in closed form by Moebius inversion over the partitions of the remaining components), Burnside over the rotations, and composition at layers holding a single brick. Each row sums to A112389(n), T(n,n) = A272690(n), and all profile counts agree with the independent refinement table of Lasse Deleuran (github.com/LasseD/A112389).
```

---

## 2. A123829 — a(10), a(11)

Current entry: #9, a(1..9), keyword `nonn`.

### DATA (append two terms)

```
1, 2, 44, 185, 3276, 15682, 282377, 1480410, 26264942, 145036229, 2556043967
```

### EXTENSIONS (add)

Deleuran published a(10) on GitHub (2025) but did not submit it, so credit him:

```
a(10) from Lasse Deleuran, confirmed by _Jeff Sponaugle_; a(11) from _Jeff Sponaugle_, Oct 03 2026
```

### LINKS (add)

```
Lasse Deleuran, <a href="https://github.com/LasseD/A112389">Counting models of 2 X 4 LEGO bricks</a>, GitHub, 2025-2026.
Jeff Sponaugle, <a href="https://github.com/jeffsponaugle/Mathpub-/tree/main/A112389">C program and notes</a>
```

### COMMENTS (add, optional)

```
a(n) = (S(n) + R(n))/2, where S(n) and R(n) count the buildings, up to translation only, that are invariant under a rotation by 180 and 90 degrees, respectively. R(n) = 0 unless n == 0 (mod 4); R(8) = 122. - _Jeff Sponaugle_, Oct 03 2026
```

`R(n) = 0` unless 4 | n: a brick is never mapped to itself by a 90° rotation, and orbits of
size 2 would overlap, so every orbit has 4 bricks. `R(8) = 122` comes from
`./lego -S 8` (S90total).

### Keywords

Optional: add `more`. a(12) is computable with the same program in roughly 9–10 hours on 10
threads (the time grew 14.7× from n = 9 to 10 and 13.1× from 10 to 11).

### Notes to the editors

```
Computed by enumerating rotation-invariant buildings directly: for each of the 4 classes of rotation centres, a Redelmeier enumeration over brick orbits under the rotation, keeping only connected unions, then Burnside over the rotation group. The program reproduces a(1)-a(9). a(10) = 145036229 agrees with the value published by Lasse Deleuran (github.com/LasseD/A112389). For n <= 11, the counts per layer profile agree with all 192 symmetric counts in Deleuran's table, and a(11) equals his partial size-11 total 2101349123 plus 454694844, the symmetric counts of the 9 layer profiles (with reversals) that his table does not yet contain.
```

---

## 3. A123831 — a(9), a(10)

Current entry: #9, a(1..8), keyword `nonn,more`.

### DATA (append two terms)

```
2, 24, 88, 1079, 5563, 74597, 412309, 5633649, 34006742, 471305645
```

### EXTENSIONS (add)

```
a(9)-a(10) from _Jeff Sponaugle_, Oct 03 2026
```

`more` stays.

### LINKS (add, optional)

Same program link as in section 2.

### Notes to the editors

```
The base brick is the only brick of the bottom layer, so the rotation must fix it and its centre is the rotation centre. a(n) is half the number of 180-degree-invariant buildings (up to translation) with n+1 bricks whose bottom layer is a single brick; the base can lie in either orientation. These were enumerated directly (Redelmeier over brick orbits). The program reproduces a(1)-a(8). a(9) = 34006742 also follows from the refinement table of Lasse Deleuran (github.com/LasseD/A112389), and a(10) = 471305645 follows from his table plus our counts for the layer profiles it does not yet contain.
```

---

## 4. A112389 — optional link fixes (no data change)

* Add a link to the ongoing computation of a(11):

  ```
  Lasse Deleuran, <a href="https://github.com/LasseD/A112389">Counting models of 2 X 4 LEGO bricks</a>, GitHub, 2025-2026.
  ```

* The arXiv link names two authors, but version 2 of arXiv:2605.07380 (30 Jul 2026) lists
  three. Change the link text to:

  ```
  Alexander Gunning, Anthony J. Guttmann and Rasmus M. Nilsson, <a href="https://arxiv.org/abs/2605.07380">Counting LEGO configurations</a>, arXiv:2605.07380 [math.CO], 2026.
  ```

Do not propose a(11) or a partial value. a(11) is still unknown, and the extrapolation
8.3649*10^18 is in the linked paper.

---

## 5. Not submittable yet, and a note to Deleuran

* **A112389 a(11).** As of 2026-09-18 Deleuran's table lacks 9 bottleneck-free refinements
  of size 11. We computed two exactly: <32321> = 16226885812568655 (861639) and
  <353> = 12004326444628855 (2459880). We also have the symmetric counts of all 9 (see
  README.md). The others are estimated at hours to days each, and <2531> and <1541> are out of
  reach with the current method.
* **A123830 a(10)** needs the same missing refinements (<1541> and <13421> start with a single
  brick).
* Suggested short note to Deleuran (not sent):

  ```
  Hi Lasse, two of the size-11 refinements missing from sum-for-size.py:
  <32321> 16226885812568655 (861639), <353> 12004326444628855 (2459880),
  computed with an independent program (layer aggregation by Moebius inversion over component
  partitions; it reproduces your table for all refinements of size <= 8 and every symmetric
  count you list up to size 11). Symmetric counts of all 9 missing refinements:
  <353> 2459880, <4421> 29863601, <3431> 3571587, <2531> 3560118, <2441> 64268008,
  <2432> 121440336, <1541> 993142, <32321> 861639, <13421> 1559051.
  I would like to submit A123829 a(10) (your 145036229, confirmed) and a(11) = 2556043967 to the
  OEIS, crediting you for a(10). Unless you prefer to do that yourself? Also, A112390 has
  T(5,2) = 248688; it should be 298688 (your table gives it).
  ```

---

## Files backing the submission

* `results/sym11.txt`: output of `./lego -SV 11`. It has the A123829 and A123831 totals for
  each n (n = 11 took 42 min on 10 threads) and S180/S90 for every layer profile.
* `results/rows_a112390.txt`: A112390 rows 2..8, recomputed profile by profile and compared
  with Deleuran's table (no mismatches).
* `b112390.txt`: rows 2..8 (to upload). `results/b112390_rows2-10.txt` also has rows 9–10
  from Deleuran's table.
* `results/deleuran-sum-for-size.py`: copy of Deleuran's table (public domain), unchanged
  since 2026-09-18.
* `results/missing11.txt`: the full runs for <32321> and <353>.
