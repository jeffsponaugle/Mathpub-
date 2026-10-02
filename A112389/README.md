# A112389 — LEGO buildings with n 2×4 bricks

[OEIS A112389](https://oeis.org/A112389): number of contiguous buildings made from
n 2×4 LEGO bricks, counted up to translation and rotation about the vertical axis.

    1, 24, 1560, 119580, 10166403, 915103765, 85747377755, 8274075616387,
    816630819554486, 82052796578652749          (a(1)..a(10); a(10) Simon 2018)

a(11) is not known exactly. Gunning–Guttmann–Nilsson (arXiv:2605.07380, 2026)
extrapolate a(11) ≈ 8.3649·10¹⁸ ± 2.4·10¹⁴ (1σ). Lasse Deleuran
([github.com/LasseD/A112389](https://github.com/LasseD/A112389)) is computing a(11)
refinement by refinement. As of 2026-09-18 his table has everything except 9
bottleneck-free refinements, and his partial sum is 7718602988330167504.

## Tool

`lego.c` is a single file with pthreads and no external libraries.

    make
    ./lego -s                 # self test: a(1..7) and a180(1..7) vs OEIS
    ./lego 7                  # a(1..7) by plain enumeration plus refinement table
    ./lego -r 353 -a 02       # one refinement by layer aggregation
    ./lego -r 3431 -a 02 -p 0.005   # sampled estimate of count and run time

Per-refinement output follows Deleuran's convention: `<profile> count (symmetric)`.
Here `count` is the number of classes up to rotation and `symmetric` is the number
of classes fixed by a 180° rotation. The `F=`, `S180=` and `S90=` fields are the
Burnside ingredients: buildings up to translation, and those fixed by a 180° or
90° rotation.

### Methods

* **Plain enumeration** (`./lego N`, `-b`): Redelmeier over the brick adjacency
  graph, rooted at the least brick, with leaf counting. Profiles with a bottleneck
  (an interior layer holding one brick) are composed from smaller ones:
  F(<R1 1 R2>) = F(<R1 1>)·F(<1 R2>)/2. Runs at about 10⁸ buildings/s/core, so it is
  only practical up to n ≈ 8.
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
  I₁ to I₃ use closed formulas. I₄ uses inclusion–exclusion over edge sets on four
  vertices, with cliques counted over grid cells, which works because overlapping
  2×4 bricks are axis-parallel boxes and so have the Helly property.
* **Symmetric counts**: Redelmeier over brick orbits under a 180° or 90° rotation
  about a fixed centre class, keeping only unions that are actually connected.

### Validation

* `./lego -s`: a(1..7) and a180(1..7) (A123829) match OEIS.
* The aggregated counter reproduces Deleuran's refinement values, including
  size-10 values such as <22222> = 697608586669144 (10421527) and
  <22321> = 1209535848675777 (2034360). Every choice of aggregated layers gives
  the same count.
* The I₄ formula was checked against brute-force DFS on every evaluation, using
  sampled runs of <44>, <45> and <141>.

## Results: missing size-11 refinements

See `results/missing11.txt`.

| refinement | count (symmetric) | method | time (10 threads, M1 Pro) | checks |
|---|---|---|---|---|
| <32321> | 16226885812568655 (861639) | A={0,2,4} | 17 s | identical F without D4 reduction (`-n`) |
| <353> | 12004326444628855 (2459880) | A={0,2} | 566 s | sampled sums identical with class tables on and off |

In both rows, F + S180 + 2·S90 is divisible by 4, and (S180 + S90) is even.
Raw values: <32321> F=64907543248551342, S180=1723278;
<353> F=48017305773595660, S180=4919760.

Cost estimates for the rest, taken from sampled runs (`-p`):

| refinement | plan | estimate |
|---|---|---|
| <3431> | A={0,2} (R = 4+1 bricks) | ≈4.4·10¹⁶ classes, ~1.3 h |
| <4421>, <2441> | A={0,2} with I₄ | several hours each (to be measured) |
| <2432>, <13421> | I₄ on a layer that touches two R-layers | ~40 h each (rough) |
| <2531>, <1541> | would need I₅ or ~10¹⁰ R-sets | not feasible with the current method |
