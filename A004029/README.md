# A004029 — number of n-dimensional space groups

**A004029:** 1, 2, 17, 219, 4783, 222018, 28927915 (n = 0…6).

This project rebuilds the full count of crystallographic (space) groups in
dimensions 1–6 with two independent C/C++ implementations. Along the way it
fills in three open terms of neighbouring OEIS sequences:

| Sequence | Meaning | OEIS before | Result here |
|---|---|---|---|
| [A395859](https://oeis.org/A395859)(6) | enantiomorphic pairs of 6-dim space groups | "7045 ≤ a(6) ≤ 7052" | **7052** |
| [A006227](https://oeis.org/A006227)(6) | 6-dim space groups including enantiomorphs | unknown | **28934967** |
| [A307288](https://oeis.org/A307288)(6) | 6-dim arithmetic classes including enantiomorphs | unknown | **85338** |

For dimension 7, where every term is unknown, it proves lower bounds:

| Sequence | 7-dim lower bound | How |
|---|---|---|
| [A004029](https://oeis.org/A004029)(7) space groups | **≥ 12,540,299,110** (probably ≈ 1.27×10¹⁰) | computation over 1.43 M arithmetic classes (below) |
|  | ≥ 882,033,440 | proof: one arithmetic class ↔ digraphs (A000273) |
| [A006227](https://oeis.org/A006227)(7) incl. enantiomorphs | ≥ 12,540,769,790 | same computation |
| [A395859](https://oeis.org/A395859)(7) enantiomorphic pairs | ≥ 470,680 | same computation |
| [A004027](https://oeis.org/A004027)(7) arithmetic classes | ≥ 1,431,317 | certified pairwise non-conjugate |
| [A004028](https://oeis.org/A004028)(7) geometric classes | ≥ 38,074 | distinct conjugacy invariants |

It also re-verifies every known term of A004029, A004027, A006227, A395859
and A307288. Ready-to-paste OEIS text is in [submissions.md](submissions.md).

## What the sequence counts

A *space group* (crystallographic group) is a discrete group R of isometries
of Rⁿ with a compact fundamental domain, i.e. the symmetry group of an
n-dimensional crystal pattern. Space groups are counted up to affine
conjugacy (equivalently, by Bieberbach's second theorem, up to abstract
isomorphism). Under this equivalence an enantiomorphic pair, a group and its
mirror image (P4₁ / P4₃), counts **once**. A006227 counts the pair twice.

- n = 1: 2 (translations only; translations plus reflections)
- n = 2: 17 wallpaper groups
- n = 3: 219 (the crystallographers' 230 counts 11 enantiomorphic pairs twice)
- n = 4: 4783 (Brown, Bülow, Neubüser, Wondratschek & Zassenhaus 1978, corrected 2002)
- n = 5, 6: 222018 and 28927915 (Plesken & Schulz 2000, with the CARAT package; a(6) corrected in 2009)

**Structure (Bieberbach).** The translations in R form a lattice L ≅ Zⁿ, and
G = R/L is a finite group acting on L, i.e. G ≤ GL_n(Z). The classification
has three layers:

1. **Geometric crystal classes (Q-classes, A004028):** G up to conjugacy in GL_n(Q).
2. **Arithmetic classes (Z-classes, A004027):** G up to conjugacy in GL_n(Z).
3. **Space-group types (A004029):** within one Z-class, R is determined by
   a class in H = H¹(G, Rⁿ/Zⁿ), a finite abelian group. Two classes give the
   same type iff they lie in the same orbit of the normalizer
   N = N_GL_n(Z)(G) (Zassenhaus' algorithm). So
   **a(n) = Σ over Z-classes of #(N-orbits on H)**.

**Enantiomorphism.** Orientation-preserving types correspond to orbits of
N⁺ = N ∩ SL_n(Z). Both counts come from one orbit computation: let N act on
H × {+1, −1} by n·(x, s) = (n·x, det(n)·s). The number of orbits is
#N⁺-orbits on H (or 2·#N-orbits if N ≤ SL_n(Z)), so

    proper  = #orbits on H × {±1}        (A006227 summand)
    affine  = #orbits on H               (A004029 summand)
    pairs   = proper − affine            (A395859 summand)

A Z-class is chiral (it splits into two SL_n(Z)-classes, which is what
A307288 counts) iff proper = 2·affine. The split extension is fixed by all
of N, so it splits exactly when N ≤ SL_n(Z).

## Results

`scripts/run_all.sh` rebuilds everything from scratch in about 2 minutes on
a 14-core M4 Max:

```
 n      Q       Z  chiralZ  Z+chiral      affine      proper  pairs  largest single Z-class
 0      1       1        0         1           1           1      0  1 groups (H^1 = 1)
 1      2       2        0         2           2           2      0  1 groups (H^1 = 1)
 2     10      13        0        13          17          17      0  3 groups (H^1 = 2x2)
 3     32      73        0        73         219         230     11  16 groups (H^1 = 2x2x2x2x2x2)
 4    227     710       70       780        4783        4894    111  272 groups (H^1 = 2^9)
 5    955    6079        0      6079      222018      222097     79  11456 groups (H^1 = 2^16)
 6   7103   85308       30     85338    28927915    28934967   7052  1540944 groups (H^1 = 2^30)

A004027 (Z-classes)                    ok, ok, ok, ok, ok, ok, ok
A307288 (Z-classes incl. enant.)       ok, ok, ok, ok, ok, ok, NEW a(6)=85338
A004029 (space groups)                 ok, ok, ok, ok, ok, ok, ok
A006227 (space groups incl. enant.)    ok, ok, ok, ok, ok, ok, NEW a(6)=28934967
A395859 (enantiomorphic pairs)         ok, ok, ok, ok, ok, ok, NEW a(6)=7052
```

Odd dimensions have no chiral Z-classes, because −I is in every normalizer
and has det −1.

In dimension 6, the 7052 pairs come from 969 Z-classes in 197 Q-classes. The
richest Q-class (`2-1;2-1;1;1 / group.1505`, 18 Z-classes) alone has 883
pairs. The 30 chiral Z-classes lie in 17 Q-classes: 12 of family `3;2-2;1`
(orders 12–96), plus point groups of order 32 (`4-1;1;1`), 64 (`4-1;2-1`),
168 and 336 (`6-3`) and 240 (`6-1`).

### Why the OEIS bound was 7045 ≤ a(6) ≤ 7052, and why the answer is 7052

Souvignier (2003) found 7052 pairs, but with the CARAT tables of the time
(28927922 affine classes, the erroneous A114268). The 2009 correction removed
7 affine classes, and nobody had checked whether any of them were
enantiomorphic. Diffing CARAT 2.0 (git `fb4b267302`, 2003) against the
current tables shows exactly one removed entry:

    group.7477 6 6-1 46080 (1)^16(2)^20 3 7 0 65 1      (CARAT 2.0 only)
    max.22     6 6-1 46080 (1)^16(2)^20 3 7 0 65 1      (both versions)

It is a second copy of the hyperoctahedral group W(B₆) = C₂ ≀ S₆ (order
46080 = 2⁶·6!). Its first generator swaps two coordinates (det −1), so all 7
duplicated types are achiral and Souvignier's 7052 stands. The direct
recomputation here also gives 7052.

### A rigorous lower bound for a(7), and digraphs

Take the point group G = {all diagonal ±1 matrices} acting on Zⁿ (the
n-dimensional "Pmmm"). Then H¹(G, Rⁿ/Zⁿ) ≅ (Z/2)^{n(n−1)}: a cocycle is a
binary n×n matrix with zero diagonal, whose entry (i, j) is the half-translation
in coordinate i attached to the reflection in coordinate j. The normalizer is
the signed permutation group W(B_n). Its sign part is G, which acts
trivially, and S_n permutes rows and columns simultaneously. So the types in
this one arithmetic class correspond exactly to **unlabeled simple digraphs
on n nodes** ([A000273](https://oeis.org/A000273)): 1, 3, 16 (Pmmm…Pnma),
218, 9608, 1540944, 882033440, … Hence

    a(n) ≥ A000273(n),   in particular   a(7) ≥ 882,033,440   (≈ 30 × a(6)).

`zverify` confirms this directly for n = 3…7 (`work/dim7/P_mmm*`; n = 7 has
|H| = 2⁴² and takes 0.6 s). For n = 6 it matches the catalog-validated count
of CARAT's Z-class `group.1080`.

## Dimension 7: lower bounds by computation

A full classification of 7-dim point groups is out of reach (see the next
section). But almost all space groups come from point groups whose rational
representation has a **one-dimensional constituent**: 84% of a(3) and a(4),
98.3% of a(5), 98.5% of a(6). Those geometric classes can all be built from
the complete 6-dim catalog:

1. **Q-classes** (`extend1`). Every finite G ≤ GL₇(Q) with a 1-dim rational
   constituent is, up to Q-conjugacy, one of the following, for some 6-dim
   Q-class G₆ (all 7103 of them):
   - G₆ ⊕ 1 (trivial on the new coordinate);
   - G₆ × C₂ (an extra −1 on it);
   - the graph of a sign character χ of G₆, i.e. gᵢ ⊕ χ(gᵢ).

   That gives 67,706 constructions. One per *fingerprint* is kept: |G| plus
   the multisets of per-conjugacy-class data (size, order, traces of
   powers), power maps, and the multiplicity and kernel of every sign
   character. Equal Q-classes always have equal fingerprints, so nothing is
   counted twice. This yields **38,074** pairwise distinct 7-dim Q-classes.
2. **Z-classes and normalizers.** CARAT's `q2z` (it runs fine in dimension 7)
   splits each Q-class, with normalizer generators for every Z-class.
3. **Counting** with the independent `zverify`.
4. **Distinctness certificate.** Z-classes of one Q-class with different
   *Z-fingerprints* are certainly not GL₇(Z)-conjugate. The Z-fingerprint
   adds the elementary divisors of g − 1, g + 1 and of the norm of ⟨g⟩ for
   each class, plus the coinvariants of L and of L*. The 580,256 look-alike
   pairs that remain were all tested with CARAT's `Z_equiv`: none is
   conjugate. Two pairs where `Z_equiv` errored are dropped from the count.
5. **Normalizer check.** A normalizer that is too small would *over*-count.
   For the 40,000 largest contributors (94.9% of the total), the normalizer
   was recomputed from scratch with CARAT's standalone `Normalizer` and the
   count redone with the union of both generating sets: 0 differences.

The run took about 20 minutes of wall-clock time on this M4 Max (14 cores) plus an M5 Pro (18 cores): 4.6 CPU-hours to split, 2.2 to count.
Of the 38,074 Q-classes, 7 were excluded. Three extensions of the
irreducible 6-dim group `6-2'/group.7619` fail inside CARAT's `q2z`. Four
tiny ones (⟨−I₆⊕1⟩, a single reflection, their product, and one hexagonal
case) were stopped after 15–27 min, because their huge form spaces make the
Voronoi normalizer step crawl; they hold only a handful of space groups.
Excluding them only lowers the bound. Result (`results/dim7_bound.txt`):

```
Z-classes counted:             1431317   (dropped: 2 unclear Z_equiv, 0 errors, 0 orphans)
space-group types (affine):    12540299110
incl. enantiomorphs (proper):  12540769790
enantiomorphic pairs:          470680
normalizer cross-check: 40000 Z-classes recomputed, 0 differences, covering 94.94% of the affine total
```

**Validation.** The same construction one dimension lower (5 → 6)
reproduces *exactly* the 3620 six-dimensional Q-classes with a 1-dim
constituent: none missing, none spurious. Built from these constructed
groups, the Q→Z splitting and counting reproduce the catalog's Z-class and
space-group numbers for all 3620 (67,362 Z-classes, 28,504,167 space
groups). The fingerprint separates all 955 five- and all 7103 six-dimensional
catalog Q-classes. In dimensions 3–6, the diagonal Q-class's counts from
this pipeline match the catalog: (4, 28), (8, 419), (16, 19325),
(36, 3164168). In 7D that class has 80 Z-classes and 1,846,316,952 space
groups.

**Rigor.** The Q-class distinctness (A004028(7) ≥ 38074) rests only on this
project's code. The other bounds also rely on CARAT's algorithms: the Q→Z
splitting (cross-checked by the Z-fingerprints and `Z_equiv`) and
normalizer completeness (cross-checked on 95% of the mass). These are the
same algorithms the published 5- and 6-dim values rest on.

**Estimate.** If, as in dimensions 5 and 6, about 98.5% of all space groups
have a point group with a 1-dim constituent, then a(7) ≈ 1.27×10¹⁰.

## Can we compute a(7)?

**Not exactly, yet.** The counting step is now cheap: `zverify`
handles any Z-class in milliseconds to seconds, including H = (Z/2)⁴². The
missing input is the classification itself, the complete list of finite
subgroups of GL₇(Z) up to conjugacy. It has never been computed: even
A004028(7) (Q-classes) and A004027(7) (Z-classes) are unknown. The section
above covers every Q-class with a 1-dim constituent. The rest, roughly
1.5% of the space groups, have point groups whose rational representation
splits into constituents of dimension ≥ 2: 7, 5+2, 4+3, 3+2+2. Computing
a(7) exactly would require:

1. **Q-classes** of that remaining type, i.e. the corresponding finite
   subgroups of the maximal finite subgroups of GL₇(Q), up to
   GL₇(Q)-conjugacy. The irreducible maximal ones are W(E₇)
   (order 2,903,040), C₂ × S₈ (80,640) and C₂ ≀ S₇ (645,120); the reducible
   ones are products over each partition of 7. This needs subgroup lattices
   and rational-character equivalence (GAP/Magma scale) plus subdirect
   products for the decomposable classes.
2. **Z-classes and normalizers**, the step that works well in 7D (above).
   CARAT's `MAXDIM 6` only affects its catalog lookups, but it still fails
   on a few groups and does not trap integer overflow.
3. **Counting.** Done (`zverify`).

Step 1 for those constituent types (irreducible groups in W(E₇), C₂×S₈,
C₂≀S₇ and their subdirect products) still needs a group-theory system such
as GAP or Magma. It is the remaining research-level piece.

## How it works

| File | Role |
|---|---|
| `carat/` | CARAT (Plesken, Schulz et al., GPL), built with `carat/Makefile.local` (no autotools needed). Supplies the Q-class catalog, the Q→Z splitting (`q2z`) and normalizer generators. |
| `src/enantio.c` | C tool on top of CARAT's library: splits each Q-class into Z-classes and dumps them (`-k`, with `ENANTIO_DUMP_DIR`). Can also count with CARAT's own cohomology code (`-w`, `-e` debug dump, `-b` print CARAT's H¹). |
| `src/zverify.cpp` | **Independent** C++ checker: shares no code with CARAT. Enumerates G, builds H¹ from the full Cayley graph (no presentation needed), computes each normalizer element's action on H¹ exactly (GMP), and counts orbits on H and H × {±1}. Small H: explicit BFS. Large H: Burnside over the image of N in Aut(H) × {±1}. |
| `scripts/verify_zclass.py` | Third, pure-Python exact implementation, used to adjudicate disagreements. |
| `scripts/run_all.sh` | Full reproduction, dims 1–6, then `scripts/summarize.py`. |
| `scripts/renorm.sh` | Recompute normalizers from scratch with CARAT's standalone `Normalizer`, then recount. |
| `src/extend1.cpp` | Builds the (n+1)-dim groups with a 1-dim constituent from n-dim Q-classes (`-o`); Q-class fingerprints (`-f`) and Z-class fingerprints (`-z`). |
| `scripts/run_dim7.sh` | The whole dimension-7 computation on one machine (`NREPS=` to test on a subset). |
| `scripts/dedup_q7.py`, `run7.py`, `zdedup7.py`, `renorm7.py`, `finalize7.py` | Its steps: one group per Q-fingerprint; Q→Z split + count (restartable, parallel, time limits); Z_equiv distinctness certificate; normalizer cross-check; final bounds. |
| `results/` | Per-Q-class and per-Z-class tables (`zverify_dim*.tsv`: file, \|G\|, H¹, affine, proper, pairs, improper, method). |
| `work/` | Z-class dumps (~500 MB; regenerated by `run_all.sh`) and the case studies. |

Self-checks inside `zverify`:
- every normalizer generator conjugates G into itself;
- every computed action is an automorphism of H (well-defined and bijective);
- G itself acts trivially on H;
- Burnside sums are divisible by the group order.

## Verification

- **Known terms:** dimensions 1–5 reproduce every OEIS value (A004027,
  A004029, A006227, A395859, A307288). For dimension 6, all 7103 Q-classes
  match CARAT's stored Z-class and affine-class counts exactly.
- **Two implementations:** `enantio` (CARAT cohomology) and `zverify`
  (independent) agree on all Z-classes except the 10 CARAT-bug cases below,
  and in each of those the Python reference sides with `zverify`.
- **Normalizer completeness:** a missing det −1 normalizer element would
  leave the affine counts intact but inflate the pair count. So for all 969
  Z-classes that contribute pairs in dimension 6, the normalizer was
  recomputed from scratch with CARAT's standalone `Normalizer` program.
  Recount: identical, 7052.
- **Literature:** the dimension-3/4/5 pair counts 11, 111, 79 match
  Neubüser–Souvignier–Wondratschek (2002) and Souvignier (2003). The
  dimension-6 value matches Souvignier's 7052 once the duplicate W(B₆) class
  is accounted for (see above).

## A bug found in CARAT

When no presentation file is given, CARAT computes a presentation of the
point group on the fly (`pres()`) and passes it to `cohomology()`. On that
path, H¹(G, Q⁶/Z⁶) comes out wrong for **10 of the 85308** six-dimensional
Z-classes. Sometimes it is too large and sometimes too small, so the fault is
not simply a missing relator:

| Z-class (family / Q-class / Z#) | CARAT H¹ | true H¹ |
|---|---|---|
| 2-1;2-1;1;1 / group.1447 / 24 | 2⁴ | 2⁴×4 |
| 3;1;1;1 / group.3846 / 14 | 2³×6 | 2⁴ |
| 3;1;1;1 / group.3846 / 15 | 6 | 2 |
| 3;1;1;1 / group.3846 / 23, 24 | 2²×6 | 2³ |
| 3;1;1;1 / group.3873 / 7 | 2⁵ | 2⁴ |
| 3;1;1;1 / min.442 / 20 | 6×6 | 6 |
| 3;1;1;1 / min.442 / 32, 33 | 6×6 | 2×6 |
| 4-1;2-2 / group.5428 / 3 | 2×6 | 2×2 |

The catalog numbers are unaffected, because they were computed with the
stored presentations (`pres.*` files). Given those, `Extensions` returns the
true H¹ and the right count in all 10 cases. Reproduce with `work/g5428/`:

    carat/bin/Extensions g.3.1 -H -n                      # C_2 X C_6  (wrong)
    carat/bin/Extensions <qcatalog>/.../pres.group.5428 g.3.1 -H -n   # C_2 X C_2

In the bad case, one normalizer generator's induced map on the wrong H¹ is
`[1 0; 0 3]` on Z/2×Z/6, which is not an automorphism. That is how the bug
was first noticed. None of the 10 classes contributes enantiomorphic pairs.
(CARAT also segfaults on 1×1 input; dimension 1 is handled by hand.)

## Reproduce

```bash
scripts/run_all.sh
```

and for the dimension-7 bounds (about 7 CPU-hours):

```bash
scripts/run_dim7.sh
```

Requirements: clang/clang++, GMP (`/opt/homebrew`), python3. CARAT's tables
are found through `CARAT_DIR` (set by the script).

## References

- W. Plesken, T. Schulz, *Counting crystallographic groups in low dimensions*, Experimental Math. 9 (2000) 407–411.
- B. Souvignier, *Enantiomorphism of crystallographic groups in higher dimensions with results in dimensions up to 6*, Acta Cryst. A59 (2003) 210–220.
- H. Brown, R. Bülow, J. Neubüser, H. Wondratschek, H. Zassenhaus, *Crystallographic Groups of Four-Dimensional Space*, Wiley 1978; corrections Acta Cryst. A58 (2002) 301.
- CARAT: https://github.com/lbfm-rwth/carat (tables at git `65619d6`; old tables at `fb4b267302`).
- A. Hoshi, M. Kang, A. Yamasaki, *Multiplicative invariant fields of dimension ≤ 6* (notes the 7104 → 7103 Q-class correction).
