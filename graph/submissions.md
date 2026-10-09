# OEIS submissions from the planar-graph enumeration

Status as of 2026-10-03, checked against the live OEIS entries and b-files on that date.

The source data is the pruned-geng enumeration of unlabeled connected planar graphs:
a(14) = 117,243,822,184 and a(15) = 2,385,519,362,107, together with the full edge
distributions (rows 14 and 15 of the connected triangle). Everything below is either
already in the OEIS from those runs, derivable from them by a formula, or a correction
found while cross-checking.

## 1. Summary

| Sequence | Counts | In the OEIS today | Still to do |
|---|---|---|---|
| A003094 | connected planar graphs, by nodes | a(14), a(15) approved (Sep 11 / Sep 21, 2026) | b-file stops at n = 14; upload `b003094.txt` (0..15), fix the credit line |
| A005470 | all planar graphs, by nodes | a(13), a(14), a(15) approved | nothing required; `b005470.txt` replaces the synthesized b-file if wanted |
| A049334 | connected planar graphs, by nodes and edges | b-file rows 1..15 (289 terms) approved | optional formulas (section 3) |
| A039735 | all planar graphs, by nodes and edges | b-file rows 1..15 (289 terms) approved | nothing |
| A145270 | connected nonplanar graphs | ends at a(13) | **submit a(14), a(15)** |
| A145269 | nonplanar graphs | ends at a(12) | **submit a(13), a(14), a(15)** |
| A173422 | partial sums of A003094 | ends at a(12) | **submit a(13)..a(15)** |
| A173310 | partial sums of A005470 | ends at a(12) | **submit a(13)..a(15)** |
| A049369 | planar graphs with no isolated node | ends at a(12) | **submit a(13)..a(15)** plus a formula |
| A174573 | partial sums of A049369 | ends at a(11) | **submit a(12)..a(15)** plus a formula |
| A343873 | connected planar graphs, by edges and nodes | b-file claims rows 0..22 but is malformed | **submit corrected b-file** `b343873_corrected.txt` |
| A046091 | connected planar graphs, by edges | ends at 22 edges | nothing new; our data agrees for 0..14 edges |
| A021103, A049336, A126201, subclasses | see section 5 | not derivable from counts alone | future runs |

Files in this directory ready to upload: `b003094.txt`, `b005470.txt`, `b343873_corrected.txt`.
The September drafts `b049334txt` and `b039735.txt` (249 terms each) are superseded by the
289-term b-files now live on the OEIS.

## 2. Cross-checks of what is live

All of the following were recomputed on 2026-10-03 from the live OEIS files. Every check passed.

| # | Check | Result |
|---|---|---|
| 1 | Row sums of the live A049334 b-file (rows 1..15) equal A003094(1..15) | pass |
| 2 | Row sums of the live A039735 b-file (rows 1..15) equal A005470(1..15) | pass |
| 3 | Euler transform of A003094(0..15) equals A005470(0..15), so a(13), a(14), a(15) of A005470 are consistent | pass |
| 4 | Two-variable Euler transform of A049334 rows 1..15 reproduces the live A039735 rows 1..15 term for term | pass |
| 5 | Sparse end of rows 14 and 15: T(n,k) = A054924(n,k) for k = n-1..n+2 (cyclomatic number at most 3 forces planarity), and the first deficit appears exactly at k = n+3 | pass |
| 6 | Dense end: T(14,36) = A000109(14) = 339,722 and T(15,39) = A000109(15) = 2,406,841 (triangulations) | pass |
| 7 | Row 14, k = 33..36 against the 2-connected triangle A049336 row 14: excesses 10,378,594; 633,015; 0; 0, each explained exactly by cut-vertex graphs (plantri + nauty, September) | pass |
| 8 | Row 15, k = 38: must equal the number of distinct graphs obtained by deleting one edge from a 15-node triangulation (no cut vertex is possible at 38 edges); plantri + deledgeg + shortg give 46,427,839 | pass |
| 9 | Rows 14 and 15 against Grasegger's transposed triangle A343873: 10 cells of row 14 (k = 13..22) and 9 cells of row 15 (k = 14..22) agree | pass |
| 10 | Column sums of A049334 for 0..14 edges equal A046091(0..14) | pass |
| 11 | Growth: a(n)/a(n-1) for n = 11..15 is 16.57, 17.96, 18.96, 19.73, 20.35, with smoothly shrinking increments | consistent |

Check 8 passes: deleting each of the 39 edges from each of the 2,406,841 triangulations on
15 nodes gives 93,866,799 graphs, of which 46,427,839 are pairwise non-isomorphic, exactly
the live T(15,38). This is a new independent reference (A049336 has no row 15), so the
dense end of row 15 is now pinned at k = 38 and 39.

What is **not** independently verified: the middle of the two new rows, k = 23..32 of row 14
and k = 23..37 of row 15. Those cells hold almost all of a(14) and a(15) and rest on the
single pruned-geng computation. Both ends of each row, the totals' consistency across four
sequences, and the two transposed-triangle overlaps are all pinned by independent data.

A by-product of check 8: the vertex-rooted 14-node triangulations number R(14) = 4,717,745,
so if anyone computes the 2-connected row 15 (A049336), its k = 37 entry must be
426,720,176 - 4,717,745 = 422,002,431. That closes the last dense-end cell of row 15.

## 3. Entry-by-entry detail and proposed text

Author credit below uses the OEIS user name shown in the edit histories, _Jeff Sponaugle_.
Dates in EXTENSIONS lines should be the actual submission date.

### A003094, connected planar graphs

Live: DATA through a(15) = 2385519362107. EXTENSIONS lines for a(14) (Sep 11 2026) and
a(15) (Sep 21 2026). COMMENT: "a(14) calculated using nauty with enhanced PREPRUNE functions
to decrease search space by ~100x." LINK: verification file `a003094.txt` for a(14)
(the contents of `Compute14.txt`). PROG line is still the old `geng -c $n | planarg -q | countg -q`.

Problems:

- The b-file link reads "Table of n, a(n) for n = 0..14" and the file indeed stops at
  n = 14 (it also carries 14 trailing blank lines). The data field has a(15).
- The credit string is garbled: "David Wasserman, Brendan McKay and Georg Grasegger,Jeff Sponaugle,".

Submit:

1. Upload `b003094.txt` (n = 0..15). Link text:
   `David Wasserman, Brendan McKay, Georg Grasegger, and Jeff Sponaugle, <a href="/A003094/b003094.txt">Table of n, a(n) for n = 0..15</a>`
2. Replace the a(14) comment with one that covers both terms and names the method, for example:
   `a(14) and a(15) were computed with nauty's geng using a PREPRUNE hook that rejects any intermediate graph that is nonplanar (planarity is closed under vertex deletion), with a left-right planarity test; candidates with more than 3n-6 edges are rejected without testing and a new vertex of degree <= 1 (or degree 2 with adjacent neighbors) is accepted without testing. Compared with generating all connected graphs and filtering, this reduces the work by a factor of roughly 100 at n = 14. - _Jeff Sponaugle_, Oct 03 2026`
3. Optionally a verification file for a(15) in the same format as `a003094.txt` (the n = 15
   aggregate log is not in this directory).

### A005470, all planar graphs

Live: DATA through a(15) = 2516143483576; EXTENSIONS "a(13)-a(14) from _Jeff Sponaugle_, Sep 11 2026"
and "a(15) from _Jeff Sponaugle_, Sep 21 2026". The b-file is synthesized from the data field (0..15).

Checks 2 and 3 confirm all three new terms. Nothing is required. `b005470.txt` (0..15) can
replace the synthesized file; editors generally accept either.

### A049334, connected planar graphs by nodes and edges

Live: b-file n = 1..289 (rows 1..15), credited to Jeff Sponaugle with rows 1..11 from
Andrew Howroyd and rows 12..13 from Georg Grasegger. PROG line is the old pipeline.
Formulas present: T(n,n-1) = A000055(n), row sums = A003094, and the log relation to A039735.

Checks 1, 4, 5, 6, 7, 8, 9, 10 all bear on this entry. Optional FORMULA additions that the
data now supports and that are not in the entry:

- `T(n,k) = A054924(n,k) for n >= 4 and k <= n+2, since a connected graph with cyclomatic number at most 3 is planar (subdivisions of K_{3,3} have cyclomatic number 4).`
- `T(n,n) = A001429(n) for n >= 3 (unicyclic graphs).`
- `T(n,3n-6) = A000109(n) for n >= 3 (triangulations).`
- `T(n,3n-7) = A049336(n,3n-7) for n >= 4, because a connected graph with a cut vertex has at most 3n-8 edges.`

### A039735, all planar graphs by nodes and edges

Live: b-file n = 1..289 (rows 1..15), credited to Jeff Sponaugle; it replaced a synthesized
76-term file. Check 4 reproduces every term from A049334, and check 2 confirms the row sums.
Nothing to do.

### A145270, connected nonplanar graphs (a(n) = A001349(n) - A003094(n))

Live: ends at a(13) = 50329965610911 (Max Alekseyev, Feb 2025). Offset 0, keyword `more`.

| n | A001349(n) | A003094(n) | a(n) |
|---|---|---|---|
| 14 | 29003487462848061 | 117243822184 | **29003370219025877** |
| 15 | 31397381142761241960 | 2385519362107 | **31397378757241879853** |

EXTENSIONS: `a(14)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A145269, nonplanar graphs (a(n) = A000088(n) - A005470(n))

Live: ends at a(12) = 164757860141 (Max Alekseyev, Aug 2015). Offset 0, keyword `more`.

| n | A000088(n) | A005470(n) | a(n) |
|---|---|---|---|
| 13 | 50502031367952 | 6295835195 | **50495735532757** |
| 14 | 29054155657235488 | 123897775473 | **29054031759460015** |
| 15 | 31426485969804308768 | 2516143483576 | **31426483453660825192** |

EXTENSIONS: `a(13)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A173422, partial sums of A003094

Live: ends at a(12) = 331953037. Offset 0, keywords `hard,nonn`. The name has a typo
("Partials sums") that could be fixed in the same edit.

| n | a(n) |
|---|---|
| 13 | **6274211345** |
| 14 | **123518033529** |
| 15 | **2509037395636** |

EXTENSIONS: `a(13)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A173310, partial sums of A005470

Live: ends at a(12) = 353222210. Offset 0, keywords `nonn,hard,more`.

| n | a(n) |
|---|---|
| 13 | **6649057405** |
| 14 | **130546832878** |
| 15 | **2646690316454** |

EXTENSIONS: `a(13)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A049369, planar graphs with minimum degree at least 1

Live: ends at a(12) = 314631443 (Sean A. Irvine, 2021). Offset 1, keywords `nonn,nice,hard`.
The entry has no formula. A planar graph with an isolated node is a planar graph on n-1 nodes
plus that node, so a(n) = A005470(n) - A005470(n-1); the published terms satisfy this.

| n | a(n) |
|---|---|
| 13 | **5962522744** |
| 14 | **117601940278** |
| 15 | **2392245708103** |

FORMULA: `a(n) = A005470(n) - A005470(n-1).`
EXTENSIONS: `a(13)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A174573, partial sums of A049369

Live: ends at a(11) = 18681007, even though a(12) has been derivable since 2014. Offset 1.
The sum telescopes: a(n) = A005470(n) - A005470(0) = A005470(n) - 1.

| n | a(n) |
|---|---|
| 12 | **333312450** |
| 13 | **6295835194** |
| 14 | **123897775472** |
| 15 | **2516143483575** |

FORMULA: `a(n) = A005470(n) - 1.`
EXTENSIONS: `a(12)-a(15) from _Jeff Sponaugle_, Oct 03 2026`

### A343873, connected planar graphs by edges and nodes (transpose of A049334)

Live: b-file by Georg Grasegger, "Table of n, a(n) for n = 0..285 (rows 0..22)", with rows
0..19 from Andrew Howroyd. Rows 0..22 of this triangle hold 276 terms, not 286, so the
link text and the file disagree. Aligning the file against A049334, A000055, A001429,
A001435 and A054924 shows exactly what happened:

- Rows 0..20 agree with every available reference: all cells with at most 15 nodes match
  A049334, and each row ends with A001435(n-1), A001429(n), A000055(n+1).
- Row 21 is missing its last entry T(21,22) = A000055(22) = 5,623,756 (trees on 22 nodes).
- Row 22 is missing its last two entries T(22,22) = A001429(22) = 164,235,501 and
  T(22,23) = A000055(23) = 14,828,074. Its other entries check out: T(22,21) = 911,325,798
  = A001435(21), T(22,20) = 3,288,208,176 = A054924(20,22), and T(22,10..15) match A049334.
- The last 13 terms (indices 273..285) are a partial row 23 for 1..13 nodes
  (0 x 9, 2306, 962057, 52417175, 894803870), which matches column 23 of A049334.

Because of the missing entries, every term from index 252 on is at the wrong index.

`b343873_corrected.txt` contains rows 0..22 (n = 0..275): the original terms with the three
missing values inserted and the partial row removed. After the fix every row n >= 4 ends
with A001435(n-1), A001429(n), A000055(n+1), as it must. Row 23 cannot be completed from
existing data: it needs connected planar graphs with 23 edges on 16..20 nodes. (Our rows
supply 14 and 15 nodes, 6,737,246,119 and 26,778,167,163, and A001435, A001429, A000055
supply 21..24 nodes, but 16..20 are unknown.)

Submit the corrected b-file with link text `Table of n, a(n) for n = 0..275 (rows 0..22)`,
keeping Grasegger's and Howroyd's credits, and a short note in the edit discussion
describing the three missing terms. This is a correction to someone else's upload, so it
is worth being explicit about the evidence above.

### A046091, connected planar graphs by edges

Live: ends at 22 edges. A graph with k edges has at most k+1 nodes, so rows 1..15 of
A049334 determine this sequence only through 14 edges, all of which agree (check 10).
Nothing to submit.

## 4. How the checks were run

Tools: nauty 2.8.9 (geng, deledgeg, shortg, countg, planarg, and a 20-line orbit counter
linked against nauty.a) and plantri 5.8, built from source. Reference values came from the
OEIS b-files of A000055, A000088, A000109, A001349, A001429, A001435, A046091, A049336,
A054924 and A343873. The Euler transform used the divisor-sum recurrence
b(n) = (1/n) * Sum_{k=1..n} c(k) b(n-k) with c(k) = Sum_{d|k} d a(d); the two-variable
version was done directly as a product of (1 - x^n y^k)^(-T(n,k)) truncated at 15 nodes
and 39 edges.

The dense-end identities use three facts: a planar graph with n >= 3 nodes has at most
3n-6 edges, with equality only for triangulations; a 2-connected planar graph with 3n-7
edges is a triangulation minus one edge (its single non-triangular face is a 4-cycle, and
at most one of its two diagonals can already be an edge); and a connected graph with a cut
vertex has at most 3n-8 edges (blocks of sizes n1 + n2 = n + 1 give 3(n+1) - 12 edges,
or 3(n-1) - 6 + 1 for a pendant vertex). The excess of T(n,3n-8) over the 2-connected count
is therefore the number of vertex-rooted triangulations on n-1 nodes.

## 5. Not derivable from the present counts, but within reach of the same generator

- **A021103 and A049336** (2-connected planar graphs, by nodes and by nodes and edges):
  both end at n = 14, from Gagarin, Labelle, Leroux and Walsh (2009), who never generated
  the graphs. Running the pruned generator with geng's `-C` switch gives row 15 directly,
  independently confirms rows 12..14, and settles the predicted T(15,37) value 422,002,431.
  The 2-connected graphs are a subset of the connected ones, so the cost is below the
  connected run.
- **A126201** (rooted connected planar graphs, ends at n = 10): a(n) is the total number of
  vertex orbits over the graphs counted by A003094(n). Calling nauty on each generated graph
  in an output procedure gives n = 11..13 cheaply and n = 14 at a few times the cost of the
  original run.
- **Subclasses countable by a filter at output time**, all currently well short of n = 15:
  A049339 and A049340 (all degrees even, end at n = 14), A049371 (minimum degree 3, n = 12),
  A243337 (K4-free connected planar, n = 12), A243787 (chordal planar, n = 14). Each adds a
  cheap test per output graph to a run that already exists.
- **A003094(16)**: the n = 13 and n = 14 logs show the work growing about 21-fold per step
  in planarity tests (4.73e9 to 9.83e10) and about 25-fold in CPU (1.81 to 44.4 core-hours).
  Extrapolating, n = 16 is on the order of 25,000 to 30,000 core-hours with an expected
  value near 5e13. Feasible on a cluster, not on one machine.
- **A046091 and A343872** (edge-indexed): extending them needs sparse graphs on up to 24
  nodes, a different regime where the planarity prune helps far less.

## 6. Files

| File | Contents |
|---|---|
| `b003094.txt` | A003094, n = 0..15, ready to upload |
| `b005470.txt` | A005470, n = 0..15, optional replacement for the synthesized b-file |
| `b343873_corrected.txt` | A343873 rows 0..22 (n = 0..275) with the three missing terms inserted |
| `b049334txt`, `b039735.txt` | September drafts through row 14; superseded by the live 289-term files |
| `Compute13.txt`, `Compute14.txt` | run logs; `Compute14.txt` is what is posted as `a003094.txt` on the OEIS |
| `planarprune.c` | prototype geng plugin using nauty's own planarity test; not the production code |
