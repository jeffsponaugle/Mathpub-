# OEIS submission candidates

Cross-checked against the live OEIS on **2026-10-03**. Entries were fetched in raw
format (`oeis.org/search?q=id:…&fmt=text`), and every candidate sequence was
searched by its terms. OEIS changes over time, so re-check right before
submitting.

## Summary

| # | Candidate | Status in OEIS (2026-10-03) | What we can submit | Ready? |
|---|-----------|------------------------------|--------------------|--------|
| 1 | **A253316 a(6)** = 4972955852202492 | **Already present.** Added by Kirill Khoruzhii (entry revision #30, Sep 28 2026) with his C++ program. Our value matches exactly. | Edit: "a(6) confirmed" line plus a link to our much faster C program | Yes |
| 2 | Takuzu grids **without** the distinct rows/columns rule | Not in OEIS | New sequence: 1, 2, 90, 11222, 12413918 | Yes (5 terms); better with n = 5, 6 |
| 3 | Takuzu grids **up to symmetry** (rotations, reflections, 0↔1 swap) | Not in OEIS | New sequence: 1, 1, 10, 321, 259400 | Yes (5 terms); better with n = 5, 6 |
| 4 | Takuzu grids up to rotations and reflections only | Not in OEIS | New sequence: 1, 1, 13, 589, 516478 | Same as #3 |
| 5 | Grids with distinct rows only (= distinct columns only) | Not in OEIS | New sequence: 1, 2, 72, 5868, 6417612 | Possible, but likely seen as derivative |
| 6 | Cross-references in A177790 | A177790 does not cite A253316 or A398298 | Small edit: "Cf. A253316, A398298." | Optional |
| 7 | A253316 a(7) | Unknown; entry still has keyword `more` | Nothing yet | No: weeks of compute on the 4 TB machine |
| 8 | Rectangular 2m × 2n Takuzu arrays (table) | Not in OEIS | Nothing yet | No: program handles square grids only |

Recommended order: #1 now. Then #2 and #3, ideally after computing n = 5 and
n = 6 (see "Extending the new sequences").

---

## 1. A253316 — a(6): independent confirmation and program

### What OEIS has now

```
%I A253316 #30 Sep 28 2026 14:46:38
%S A253316 1,2,72,4140,4111116,48183195384,4972955852202492
%N A253316 Number of 2n X 2n Takuzu grids.
%H A253316 Kirill Khoruzhii, <a href="/A253316/a253316_2.cpp.txt">C++ program for computing a(n), n = 0..6</a>
%H A253316 Wikipedia, <a href="https://en.wikipedia.org/wiki/Takuzu">Takuzu</a>
%Y A253316 Cf. A058527.
%Y A253316 Number of possible rows equals A177790.
%K A253316 nonn,more,changed
%O A253316 0,2
%A A253316 _Brian Kell_, Dec 30 2014
%E A253316 a(5) from _Frans J. Faase_, Dec 14 2015
%E A253316 a(6) from _Kirill Khoruzhii_, Sep 27 2026
```

a(6) entered OEIS on the same day our run finished. Our value is identical:
**4972955852202492**.

### Why a confirmation is worth submitting

- **Different algorithm.** Khoruzhii's program enumerates the first 9 rows (up to
  0↔1 swap and left↔right reflection) and counts the last 3 rows from a bitmap
  index. Its header quotes **about 650 single-thread hours** for a(6). Ours splits
  the grid between rows 6 and 7 and counts compatible top/bottom halves for each
  symmetry class of middle-row pairs. It took **31 minutes on 14 threads (about 7
  CPU-hours)**, roughly 100× less compute.
- **Checked internally twice.** Two full runs used two independent pair-counting
  methods (bitset and inclusion–exclusion). All 5548 work units gave identical
  counts (`compare_ckpt.sh`). The self-test reproduces n = 1..5, and symmetry
  orbit checks pass on n = 6.

### Proposed edit (draft)

Replace `Your Name` / `_Your Name_` with your OEIS registered name. OEIS writes
the plain name in `%H` lines and `_Name_` in `%E` lines. Upload `a253316.c` as an attachment in the edit form; OEIS
assigns the final file name (probably `a253316_3.c.txt`).

```
%H A253316 Your Name, <a href="/A253316/a253316_3.c.txt">C program</a> (meet-in-the-middle at the middle rows, multithreaded, with checkpoints; a(6) in about 7 CPU-hours)
%E A253316 a(6) confirmed by _Your Name_, Sep 28 2026
```

"a(n) confirmed by" lines are a common OEIS convention. Editors may still choose
to keep only one of the two lines.

---

## 2. New sequence: no distinctness rule

**Searched OEIS for:** `2,90,11222,12413918`, `1,2,90,11222`, `90,11222,12413918`,
`12413918` alone. **No results** (the last search also rules out the terms
appearing inside a table).

### Draft entry

```
%N Number of 2n X 2n 0-1 matrices with n ones in each row and each column and no three consecutive equal entries in any row or column.
%C These are the Takuzu grids of A253316 without the requirement that all rows and all columns be distinct.
%C a(n) <= A058527(n), with equality for n <= 2 (a run of three is impossible in a balanced line of length 4).
%S 1,2,90,11222,12413918
%Y Cf. A058527, A177790 (number of possible rows), A253316.
%K nonn,more
%O 0,2
%A _Your Name_, <date>
```

**Provenance.** n ≤ 4 by brute force (`tools/variants_bruteforce.c`). M(4) was also
reproduced by an independent row-by-row dynamic program during development.

---

## 3 and 4. New sequences: Takuzu grids up to symmetry

**Searched OEIS for:** `1,13,589,516478`, `13,589,516478`, `516478`,
`1,10,321,259400`, `10,321,259400`, `259400 589`. **No results.** The keyword search
`takuzu` returns only A253316, A398298 and A177790.

### Draft entries

```
%N Number of 2n X 2n Takuzu grids up to rotations, reflections and complementation (interchanging 0 and 1).
%C See A253316 for the definition of a Takuzu grid.
%S 1,1,10,321,259400
%Y Cf. A253316.
%K nonn,more
%O 0,3
```

```
%N Number of 2n X 2n Takuzu grids up to rotations and reflections.
%S 1,1,13,589,516478
%Y Cf. A253316.
%K nonn,more
%O 0,3
```

**Provenance.** n ≤ 4 by canonical forms over all grids, cross-checked with
Burnside's lemma (`tools/symmetry_classes.c`). Both methods agree for every n.

### Burnside fixed-point counts (for comments, and for extending to n = 5, 6)

| n | identity | rot 90° | rot 180° | rot 270° | flip top↔bottom | flip left↔right | transpose | anti-transpose |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | 2 | 0 | 0 | 0 | 2 | 2 |
| 2 | 72 | 0 | 8 | 0 | 0 | 0 | 12 | 12 |
| 3 | 4140 | 0 | 156 | 0 | 0 | 0 | 208 | 208 |
| 4 | 4111116 | 8 | 2996 | 8 | 0 | 0 | 8848 | 8848 |

The same symmetries composed with complementation:

| n | complement | rot 90° | rot 180° | rot 270° | flip top↔bottom | flip left↔right | transpose | anti-transpose |
|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 2 | 0 | 2 | 2 | 2 | 0 | 0 |
| 2 | 0 | 4 | 0 | 4 | 24 | 24 | 0 | 0 |
| 3 | 0 | 8 | 0 | 8 | 204 | 204 | 0 | 0 |
| 4 | 0 | 80 | 32 | 80 | 9192 | 9192 | 0 | 0 |

Five columns are always 0. A flip would repeat rows (top↔bottom) or columns
(left↔right). Complementation alone cannot fix a grid. Complement∘transpose and
complement∘anti-transpose would need an entry on the diagonal (or anti-diagonal)
equal to its own complement.

---

## 5. New sequence: distinct rows only

**Searched OEIS for:** `2,72,5868,6417612`, `1,2,72,5868`, `72,5868,6417612`, `6417612`.
**No results.**

```
%N Number of 2n X 2n 0-1 matrices with n ones in each row and each column, no three consecutive equal entries in any row or column, and all rows distinct.
%C By transposition, also the number of such matrices with all columns distinct.
%S 1,2,72,5868,6417612
%Y Cf. A253316, A058527.
%K nonn,more
%O 0,2
```

This is lower priority, since editors may consider it too close to #2 and A253316.
For reference, the brute-force split for n = 4 (`tools/variants_bruteforce.c`):

| n = 4 | count |
|---|---|
| all balanced, no-three-in-a-row grids (#2) | 12413918 |
| rows and columns distinct (A253316) | 4111116 |
| rows distinct, some columns repeated | 2306496 |
| columns distinct, some rows repeated | 2306496 |
| both a repeated row and a repeated column | 3689810 |

---

## 6. Minor edit: A177790 cross-references

A177790 (number of valid rows: 1, 2, 6, 14, 34, 84, 208, 518, 1296, …) is cited by
A253316 and A398298. Its own cross-reference line is only
`Equals twice A003440 …`. Optional edit:

```
%Y A177790 Cf. A253316, A398298.
```

Related entries, checked and needing nothing:
- **A398298** (valid row patterns, Stef Bordo, Jul 2026) already cites A177790 and A253316.
- **A058527** (balanced 0-1 matrices) already cites A253316.
- **A003440** equals A177790 / 2. A253316 now cites A177790 rather than A003440, which is correct.

If #2 is accepted, A058527 and A253316 should get a "Cf." to the new number.
Editors usually add these back-references.

---

## 7 and 8. Not submittable yet

- **a(7)** (14 × 14). Estimated around 10²¹–10²². It needs about 1.1×10¹⁵
  half-grids of work, roughly 9–12 CPU-years, and up to about 5 TB per work unit.
  That is plausible only on the 4 TB machine, over weeks, after porting the
  program to n = 7.
- **Rectangular 2m × 2n Takuzu arrays.** No such table exists in OEIS. The seam
  method generalizes to rectangles, but `a253316.c` is square-only today.

---

## Extending the new sequences before submitting

OEIS favours more than 4–5 terms. Rough effort with what we have:

| Sequence | n = 5 | n = 6 | n = 7 |
|---|---|---|---|
| #2 (no distinctness) | minutes | under an hour (no pair checks needed) | about 10¹⁵ half-grids to enumerate but only small histograms to keep: roughly a week on this Mac, about a day on a large server |
| #3 / #4 (up to symmetry) | new code for the fixed-point counts; minutes of compute | same code; the 180°-rotation and complement∘flip counts need one pass over all ~10¹⁰ top halves (minutes), the transpose counts a separate small search | no |
| #5 (distinct rows only) | minutes | similar to the a(6) run or less | no |

#2 and #5 need a small "variant" mode in `a253316.c`: allow repeated rows and
equal middle rows, and drop one or both distinctness tests.

## How to submit

1. **Account.** Log in or request an account at <https://oeis.org/login>. OEIS
   requires a real name.
2. **Edit an existing entry** (#1, #6). Open the entry, click **edit**, add the
   lines, attach the program file, and submit for review.
3. **New sequence** (#2–#5). Use <https://oeis.org/Submit.html>. Check the format
   against the [Style Sheet](https://oeis.org/wiki/Style_Sheet): name, terms,
   offset, comments, cross-references, keywords and author.
4. **Programs.** `a253316.c` is about 1600 lines, so attach it as a file rather
   than an inline `%o` program. The small tools in `tools/` are short enough for
   inline `%o (C)` lines or attachments.
5. **Review.** Editors review each change, usually within days to a few weeks;
   expect questions or edits. Re-check the live entries just before submitting.

## Files backing these numbers

| File | Supports |
|---|---|
| `a253316.c`, `a253316_n6.ckpt`, `a253316_n6_ie.ckpt`, `compare_ckpt.sh` | #1: a(6), two independent full runs |
| `tools/variants_bruteforce.c` | #2, #5: n ≤ 4 |
| `tools/symmetry_classes.c` | #3, #4: n ≤ 4, with the Burnside tables |
