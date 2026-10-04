# OEIS submissions from the A034822 project

Prepared 2026-10-03. All data below comes from Patrick De Geest's table of square palindromes
(worldofnumbers.com, "Plain Text Squares.txt", update of Nov 16 2025, complete through 69 digits; the 68- and
69-digit squares were found with Robert Xiao's CUDA program). Nothing here depends on our own L = 70 search,
which has not been run. Credit the data to De Geest / Xiao; we are only transferring it into b-files.

**Current OEIS state, checked 2026-10-03:** no new revisions since 2026-09-29, and no pending drafts on any of the
six sequences: A034822 #35 (Jun 04 2026), A263618 #29 (Apr 08 2025), A016113 #34 (Apr 02 2025),
A002778 #100 (Nov 05 2025), A002779 #117 (Sep 16 2026), A027829 #51 (Jul 05 2025). An OEIS search for the
68-digit root 7256171055736382499839982033184475 returns no results.

Files are in [`oeis/`](oeis/): proposed b-files `b*.txt`, the current OEIS b-files `orig_b*.txt`, De Geest's
table `degeest_squares_2025-11-16.txt`, and `verify_bfiles.py`, which checks everything (`python3 verify_bfiles.py`
prints ALL OK). Every proposed b-file extends the current one unchanged.

## Submission order

Submit the b-files in this order, so the cross-references hold.

### 1. A263618 — Number of palindromic squares with exactly n digits

* **New terms:** a(68) = 1, a(69) = 1443. Data line becomes `..., 578, 0, 1199, 2, 701, 1, 1443`.
* **b-file:** `oeis/b263618.txt` (n = 1..69). There is currently no %H b-file line; the OEIS file is synthesized
  from the data.
* %H: `Jeff Sponaugle, <a href="/A263618/b263618.txt">Table of n, a(n) for n = 1..69</a> (from Patrick De Geest's table)`
* %E: `a(68)-a(69) from Patrick De Geest's table (computed with Robert Xiao's CUDA program) added by _Jeff Sponaugle_, Oct 2026`
* Optional %C: `a(68) = 1: 7256171055736382499839982033184475^2 (Patrick De Geest, Apr 14 2025).`

### 2. A002778 — Numbers whose square is a palindrome

* **b-file:** `oeis/b002778.txt`, n = 1..10173 (current 1..8729). Adds n = 8730 = 7256171055736382499839982033184475
  (the 68-digit square) and n = 8731..10173 (the 1443 roots of the 69-digit squares).
* %H: replace the current b-file line with
  `Patrick De Geest, <a href="/A002778/b002778.txt">Table of n, a(n) for n = 1..10173</a> (all terms with squares of at most 69 digits; terms 8730..10173 added by _Jeff Sponaugle_)`
  (or keep the "communicated by Max Alekseyev" form and only change the range; follow the editor's preference).

### 3. A002779 — Palindromic squares

* **b-file:** `oeis/b002779.txt`, n = 1..10173 (current 1..1940, i.e. squares up to 45 digits). These are the
  squares of the A002778 terms, all palindromic squares with at most 69 digits. Size about 630 KB.
* %H: `Jeff Sponaugle, <a href="/A002779/b002779.txt">Table of n, a(n) for n = 1..10173</a> (from Patrick De Geest's table)`
  (keep the existing Havermann/Noe credit for the first 1940 terms if the editor prefers).

### 4. A016113 — Numbers whose square is a palindrome with an even number of digits

* **New term:** a(23) = 7256171055736382499839982033184475.
* **b-file:** `oeis/b016113.txt` (n = 1..23, current 1..22).
* %H: change `n = 1..22` to `n = 1..23`.
* %E: `a(23) from Patrick De Geest's table (found Apr 14 2025 with Robert Xiao's CUDA program) added by _Jeff Sponaugle_, Oct 2026`
* Optional %C: `The search of De Geest and Xiao is exhaustive through 69 digits, so a(24) has at least 70 digits.`

### 5. A027829 — Palindromic squares with an even number of digits

* **b-file:** `oeis/b027829.txt`, n = 1..23 (current 1..14, i.e. up to 44 digits). These are the squares of the
  A016113 terms. a(23) = 52652018390106447787037098707897055079870789073078774460109381025625 (68 digits).
* %H: `Jeff Sponaugle, <a href="/A027829/b027829.txt">Table of n, a(n) for n = 1..23</a> (from Patrick De Geest's table)`

### 6. A034822 — Lengths with no palindromic square (comment only, no new term)

* Proposed %C: `All lengths up to 69 have been searched exhaustively (Patrick De Geest with Robert Xiao's CUDA
  program); there is exactly one palindromic square of length 68 (7256171055736382499839982033184475^2), so
  a(17) >= 70. - _Jeff Sponaugle_, Oct 2026`
* Optional %H: `Patrick De Geest, <a href="https://www.worldofnumbers.com/Plain%20Text%20Squares.txt">Table of square palindromes</a>`
  (already linked from A002778).

## Not yet submittable

* **A034822 a(17).** It is either 70 or larger and needs the exhaustive L = 70 search (`cuda/a034822_cuda`,
  ~10 days on two DGX Sparks; not run, see README). If the run finds no 70-digit palindromic square, submit
  a(17) = 70 with A263618 a(70) = 0. If it finds one, submit it to A016113/A027829/A002778/A002779/A263618
  (a(70) >= 1) and note a(17) >= 72. Partial runs give nothing submittable unless they find a square.
