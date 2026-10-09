# OEIS submissions — space-group counts (A004029 family)

Results of the dimension 1–6 recount in this directory (see README.md),
checked against the live OEIS entries as of **2026-10-04**.

## Status summary

| Sequence | Definition | OEIS last term | Our result | Status |
|---|---|---|---|---|
| [A395859](https://oeis.org/A395859) | enantiomorphic pairs among n-dim space groups | a(5) = 79; comment "7045 <= a(6) <= 7052" | **a(6) = 7052** | 📝 to submit |
| [A006227](https://oeis.org/A006227) | n-dim space groups incl. enantiomorphs | a(5) = 222097 | **a(6) = 28934967** | 📝 to submit |
| [A307288](https://oeis.org/A307288) | arithmetic classes incl. enantiomorphs | a(5) = 6079 | **a(6) = 85338** | 📝 to submit |
| [A004029](https://oeis.org/A004029) | n-dim space groups | a(6) = 28927915 | all terms re-verified; **a(7) ≥ 12540299110** (and a(n) ≥ A000273(n)) | 📝 comment |
| [A000273](https://oeis.org/A000273) | unlabeled simple digraphs | — | crystallographic interpretation | 📝 comment |
| A006227 / A395859 / A004027 / A004028 | 7-dim terms | unknown | a(7) ≥ 12540769790 / 470680 / 1431317 / 38074 | 📝 comments (optional) |

All dims 0–5 values of A004027, A004029, A006227, A395859 and A307288 are
reproduced exactly, and all 7103 six-dimensional geometric classes match
CARAT's stored Z-class and affine-class counts. Without a public link, the
"computed by" credits rest on the description alone. Consider pushing this
directory to GitHub and adding it as a `%H` link to each edit.

## 1. A395859 — a(6) = 7052

**Data:** `0, 0, 0, 11, 111, 79, 7052`

**Remove** the comment `7045 <= a(6) <= 7052.` and **add**:

> %C The CARAT tables used by Souvignier (2003) to obtain a(6) = 7052 listed
> the geometric class of the hyperoctahedral group (order 46080 = 2^6*6!)
> twice. That duplicate accounts for the 7 extra affine classes in
> A114268(6) - A004029(6). The point group contains reflections, so none of
> those 7 classes is enantiomorphic and a(6) = 7052 stands. - _Jeff Sponaugle_, Oct 04 2026

> %E a(6) = 7052 confirmed by _Jeff Sponaugle_, Oct 04 2026: recomputed for
> all 85308 arithmetic classes of the corrected CARAT tables, with the
> cohomology H^1(G, Q^6/Z^6) and the normalizer action computed independently
> of CARAT.

One-line justification for the editors: *"Diffing CARAT 2.0 (2003) against
CARAT 2.1 shows the only removed geometric class is group.7477 (family 6-1,
order 46080, 3 Z-classes, 7 affine classes), a copy of max.22 = W(B6). Its
generators include a coordinate swap (det -1), so the 7 duplicate types are
achiral."*

## 2. A006227 — a(6) = 28934967

**Data:** `1, 2, 17, 230, 4894, 222097, 28934967`

> %F a(n) = A004029(n) + A395859(n).

> %E a(6) = 28934967 from _Jeff Sponaugle_, Oct 04 2026: counted directly over
> all 85308 arithmetic classes of the corrected CARAT tables, as the orbits of
> N cap SL_6(Z) on H^1(G, Q^6/Z^6), where N is the normalizer of G in GL_6(Z);
> this equals A004029(6) + A395859(6) = 28927915 + 7052. (Souvignier's 2003
> total 28934974 = 28927922 + 7052 was based on the old CARAT tables, which
> double-counted 7 achiral affine classes.)

## 3. A307288 — a(6) = 85338

**Data:** `1, 2, 13, 73, 780, 6079, 85338`

> %C a(n) = A004027(n) + (number of arithmetic classes whose normalizer in
> GL_n(Z) lies in SL_n(Z)). For odd n, -I has determinant -1, so
> a(n) = A004027(n). - _Jeff Sponaugle_, Oct 04 2026

> %E a(6) = 85338 from _Jeff Sponaugle_, Oct 04 2026: 30 of the 85308
> six-dimensional arithmetic classes split into enantiomorphic pairs. This
> agrees with Souvignier (2003), Table 2, which lists 85311 (+30) arithmetic
> classes using the old CARAT tables; the 3 extra classes come from the
> duplicated hyperoctahedral geometric class (see A395859), which contains
> reflections and so is achiral.

> %H B. Souvignier, <a href="https://doi.org/10.1107/S0108767303004161">Enantiomorphism of crystallographic groups in higher dimensions with results in dimensions up to 6</a>, Acta Cryst., A59 (2003), 210-220.

Supporting detail (not for the entry): Souvignier's remarks name the three
geometric classes in which every arithmetic class is enantiomorphic,
C2 x S5, PSL(2,7) and C2 x PSL(2,7). In our run these are group.7457
(order 240), group.7621 (168) and group.7624 (336), with 3/3, 2/2 and 2/2
chiral arithmetic classes.

## 4. A004029 — comments (lower bounds for a(7))

The computational bound:

> %C a(7) >= 12540299110. This counts the space groups whose point group's
> rational representation has a one-dimensional constituent. There are at
> least 38074 such geometric classes in dimension 7, built from the 7103
> six-dimensional ones; they split into 1431317 arithmetic classes, computed
> with CARAT and counted by an independent program. In dimensions 5 and 6
> such point groups give 98.3% and 98.5% of all space groups, which suggests
> a(7) is about 1.27*10^10. - _Jeff Sponaugle_, Oct 04 2026

and the structural one:

> %C a(n) >= A000273(n): the space groups whose point group is the group of
> all diagonal matrices with entries +-1 acting on Z^n form one arithmetic
> class, and their types correspond to unlabeled digraphs on n nodes. Here
> H^1(G, R^n/Z^n) is the space of binary n X n matrices with zero diagonal,
> and the normalizer (signed permutation matrices) acts on it by
> simultaneous permutation of rows and columns. This class gives 3 of the 17
> plane groups (p2mm, p2mg, p2gg) and 16 of the 219 space groups (Pmmm, ...,
> Pnma). In particular a(7) >= 882033440. - _Jeff Sponaugle_, Oct 04 2026

> %Y add A000273.

## 4b. Optional 7-dim comments for the sibling entries

Same computation (README, "Dimension 7"); put a `%H` link to the code next
to each:

> A006227: %C a(7) >= 12540769790. - _Jeff Sponaugle_, Oct 04 2026
>
> A395859: %C a(7) >= 470680. - _Jeff Sponaugle_, Oct 04 2026
>
> A004027: %C a(7) >= 1431317 (arithmetic classes of the 7-dimensional point
> groups with a one-dimensional rational constituent). - _Jeff Sponaugle_, Oct 04 2026
>
> A004028: %C a(7) >= 38074: there are at least that many geometric classes
> of 7-dimensional point groups with a one-dimensional rational constituent,
> pairwise distinguished by conjugacy-class invariants. - _Jeff Sponaugle_, Oct 04 2026

The A004028 bound rests only on our own code. The others also rely on
CARAT's Q→Z splitting and normalizers in dimension 7, cross-checked as
described in the README (`Z_equiv` certificates for all 580,256 look-alike
pairs; normalizers recomputed for 95% of the mass with 0 differences).

## 5. A000273 — comment

> %C a(n) is the number of types of n-dimensional space groups whose point
> group is the group of all diagonal sign matrices acting on the lattice Z^n
> (the n-dimensional analogs of Pmmm); see A004029. - _Jeff Sponaugle_, Oct 04 2026

## 6. (Optional) CARAT bug report — draft for https://github.com/lbfm-rwth/carat/issues

> **Wrong H^1 when the presentation is computed on the fly (10 six-dimensional Z-classes)**
>
> For 10 of the 85308 Z-classes of dimension 6, `Extensions file -H` (no
> presentation file, so `pres()` is used) returns a wrong
> H^1(G,Q^n/Z^n). With the stored `pres.*` file from the Q-catalog, the
> result is correct. Example: QtoZ splits `dim6/dir.4-1;2-2/.../group.5428`
> into 10 Z-classes; for the 3rd one, `Extensions g.3.1 -H` gives
> `C_2 X C_6`, while `Extensions pres.group.5428 g.3.1 -H` gives
> `C_2 X C_2` (correct, confirmed by two independent implementations). In the
> wrong case, normalop() maps one normalizer generator to `[1 0; 0 3]` on
> Z/2 x Z/6, which is not an automorphism, and the orbit lengths sum to 8,
> not 12. Complete list: group.1447 Z24, group.3846 Z14/Z15/Z23/Z24,
> group.3873 Z7, min.442 Z20/Z32/Z33, group.5428 Z3 (in QtoZ order). The
> catalog counts are unaffected. Also: `QtoZ`/`q2z` segfaults on 1x1 groups,
> and in dimension 7 `q2z` fails with "matrix is singular. Can't invert" on
> 6-2'/group.7619 (+) 1, and `Z_equiv` errors ("matrix is singular",
> "Normal: Error: divide by zero") on 2-2;2-2;1,1/group.2704 and 2705
> extended by a sign character (Z-classes 9 vs 31).
