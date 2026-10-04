# A112389 — LEGO buildings with n 2×4 bricks

[OEIS A112389](https://oeis.org/A112389): number of contiguous buildings made from
n 2×4 LEGO bricks, counted up to translation and rotation about the vertical axis.

    1, 24, 1560, 119580, 10166403, 915103765, 85747377755, 8274075616387,
    816630819554486, 82052796578652749          (a(1)..a(10); a(10) Simon 2018)

Status (2026-10-03): **paused.** a(11) is still unknown. Lasse Deleuran
([github.com/LasseD/A112389](https://github.com/LasseD/A112389)) has every size-11
refinement except 9. We computed two of those 9 (<32321> and <353>) and the
symmetric counts of all 9. The related 2×4 entries A123829, A123831 and A112390 can
be extended or corrected now; see [submission.md](submission.md).

## Results

| what | value | how it was checked |
|---|---|---|
| A123829(10) (180°-symmetric buildings) | 145036229 | also equals Deleuran's (unsubmitted) value |
| **A123829(11)** | **2556043967** | Deleuran's partial sum plus our counts for his 9 missing refinements gives the same total |
| **A123831(9)** (symmetric, on a fixed base brick) | **34006742** | also follows from Deleuran's table alone |
| **A123831(10)** | **471305645** | Deleuran's table plus our values for the missing refinements |
| A112390 T(5,2) | **298688**, not 248688 as in the OEIS | row 5 must sum to A112389(5) = 10166403; recomputed |
| A112390 rows 7 and 8 | see `b112390.txt` | recomputed refinement by refinement; each row sums to A112389(n) |
| <32321> (size-11 refinement) | 16226885812568655 (861639) | two runs (with and without D4 reduction) |
| <353> (size-11 refinement) | 12004326444628855 (2459880) | same-seed samples identical across code paths |

All 192 symmetric counts in Deleuran's table (every bottleneck-free refinement with
n ≤ 11) agree with ours (`sym_vs_lasse.py`). Every refinement of size ≤ 8 agrees
with his table too (`results/rows_a112390.txt`).

Symmetric counts of the 9 refinements still missing from a(11) (classes fixed by
a 180° rotation): <353> 2459880, <4421> 29863601, <3431> 3571587,
<2531> 3560118, <2441> 64268008, <2432> 121440336, <1541> 993142,
<32321> 861639, <13421> 1559051.

## Tool

`lego.c` is a single file with pthreads and no external libraries.

    make
    ./lego -s                 # self test: a(1..7) and a180(1..7) vs OEIS
    ./lego -S 11              # A123829(1..11), A123831(1..10): symmetric counts only (~45 min)
    ./lego -SV 11             # the same, plus S180/S90 for every layer profile
    ./lego 7                  # a(1..7) by plain enumeration plus refinement table
    ./lego -r 353 -a 02       # one refinement by layer aggregation
    ./lego -r 3431 -a 02 -p 0.005   # sampled estimate of count and run time
    ./lego -b 1421            # one refinement by plain enumeration (cross-check)

Helper scripts:

* `rows_a112390.py SYM NMAX [DELEURAN]` recomputes A112390 rows from `./lego`:
  aggregation for bottleneck-free profiles and composition for the others.
  Optionally it compares every profile with Deleuran's table.
* `sym_vs_lasse.py SYM DELEURAN` compares the symmetric counts with Deleuran's.
* `verify_lasse.py DELEURAN SIZE [secs] [profiles]` recomputes chosen refinements
  and compares them with Deleuran's.

`results/deleuran-sum-for-size.py` is a copy of Deleuran's table, unchanged since
2026-09-18; his code is public domain.

Per-refinement output follows Deleuran's convention: `<profile> count (symmetric)`.
Here `count` is the number of classes up to rotation and `symmetric` is the number
of classes fixed by a 180° rotation. The `F=`, `S180=` and `S90=` fields are the
Burnside ingredients: buildings up to translation, and those fixed by a 180° or
90° rotation, so count = (F + S180 + 2·S90)/4 and symmetric = (S180 + S90)/2.

### Methods

* **Plain enumeration** (`./lego N`, `-b`): Redelmeier over the brick adjacency
  graph, rooted at the least brick, with leaf counting. Profiles with a bottleneck
  (an interior layer holding one brick) are composed from smaller ones:
  F(<R1 1 R2>) = F(<R1 1>)·F(<1 R2>)/2 and S180 likewise. Runs at about 10⁸
  buildings/s/core, so it is only practical up to n ≈ 8.
* **Layer aggregation** (`-r`): choose pairwise non-adjacent layers A and
  enumerate only the remaining layers R. A brick in an A-layer touches only R-bricks,
  so for fixed R-bricks with components K₁..K_r the A-layers are counted by Möbius
  inversion over partitions σ of the components:

      #buildings = Σ_σ μ(σ,1̂) · Π_{l∈A} I_{z_l}(X^l_σ)

  X^l_σ is the set of placements in layer l that touch R-bricks lying inside a
  single block of σ, and I_k(X) counts the k-subsets of X with no two bricks
  overlapping. R is enumerated by Redelmeier on an augmented graph: real overlaps,
  plus "virtual" edges between two R-bricks that one A-brick could overlap. Each
  R-set is weighted by its D4 orbit, so only one image per orbit is evaluated.
  I₁ to I₃ come from per-class tallies (O(q³) per partition, q = number of touch-mask
  classes). I₄ uses inclusion–exclusion over edge sets on four vertices, with
  cliques counted over grid cells, which works because overlapping 2×4 bricks are
  axis-parallel boxes and so have the Helly property.
* **Symmetric counts** (`-S`): Redelmeier over brick orbits under a 180° or 90°
  rotation about a fixed centre class (4 classes for 180°, 2 for 90°), rooted at the
  least orbit, keeping only unions that are actually connected. No F computation is
  needed, which is why A123829(11) and A123831(10) are within reach while
  A112389(11) is not.

### Validation

* `./lego -s`: a(1..7) and a180(1..7) (A123829) match OEIS.
* `./lego -S 11` reproduces A123829(1..9) and A123831(1..8) exactly, including the
  122 90°-symmetric buildings at n = 8 (all for profile <44>).
* The aggregated counter reproduces Deleuran's refinement values, including size-10
  values such as <22222> = 697608586669144 (10421527) and <22321> = 1209535848675777
  (2034360). Different choices of aggregated layers give the same count.
* The I₄ formula and the class tables were checked against brute-force bitset/DFS
  counts on sampled runs (`-x 1`, `-c 1000000`).
* Bug found and fixed on 2026-10-03: the memo key for I_k packed the layer tag and
  k into the top 6 bits of the key. Layer 1 with k=1 then collided with layer 0
  with k=4, and the only plan affected was <1421> with `-a 13`. It was caught
  because row 8 of A112390 did not sum to a(8). The fixed build reproduces
  Deleuran's value 60442092848, and same-seed sampled sums for <32321>, <353> and
  <3431> are unchanged, so those results were not affected.

## Missing size-11 refinements: costs

From sampled runs (`-p`), 10 threads on an M1 Pro:

| refinement | plan | estimate |
|---|---|---|
| <3431> | A={0,2} (R = 4+1 bricks) | ≈4.4·10¹⁶ classes, ~1.3 h |
| <4421>, <2441> | A={0,2} with I₄ | several hours each (to be measured) |
| <2432>, <13421> | I₄ on a layer that touches two R-layers | ~40 h each (rough) |
| <2531>, <1541> | would need I₅ or ~10¹⁰ R-sets | not feasible with the current method |

To resume: `./lego -r 3431 -a 02`, then `results/followup.sh`, which checks two
known size-10 values and then compares 25 of Deleuran's size-11 values.

## Files

* `results/missing11.txt` — full runs of the missing refinements.
* `results/sym11.txt` — `./lego -SV 11`: per-n totals and S180/S90 for every
  profile with n ≤ 11 (n = 11 took 42 min).
* `results/rows_a112390.txt` — A112390 rows 2..8 recomputed, with a per-profile
  comparison against Deleuran.
* `b112390.txt` — b-file for A112390 rows 2..8 (index from 2, as in %O 2,1).
  `results/b112390_rows2-10.txt` also has rows 9..10, which come from
  Deleuran's table.
