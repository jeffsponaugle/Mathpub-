# OEIS submission: A285873–A285878

Proposed edits from this project, as of **Oct 03 2026**. Each sequence is a separate edit on
oeis.org: sign in, choose **Edit** on the sequence, change the **Data** field, add the
**%E** line, and paste the explanation into the edit's discussion box. Replace
`_<your OEIS name>_` with your OEIS user name in the usual `_First Last_` form.

**Ready** means two independent computations agree: different Zobrist hash keys, different
task decompositions, and usually different machines. **Pending** means one computation has
finished and its confirming run is queued or in progress.

## Summary

| Sequence | Change | Terms | Status |
|---|---|---|---|
| [A285874](https://oeis.org/A285874) no rooks | extend | a(10)–a(13) | **ready** |
| [A285877](https://oeis.org/A285877) no pawns | **correct** | a(6), a(7) | **ready**: three independent programs agree, including Stockfish 19 |
| [A285877](https://oeis.org/A285877) no pawns | extend | a(8)–a(10) | **ready** |
| [A285878](https://oeis.org/A285878) pawns and king | extend | a(11)–a(13) | **ready** |
| [A285875](https://oeis.org/A285875) no knights | extend | a(10)–a(12) | **ready** |
| [A285873](https://oeis.org/A285873) no queens | extend | a(10)–a(12) | **ready** |
| [A285876](https://oeis.org/A285876) no bishops | extend | a(10)–a(12) | **ready** |

All existing OEIS terms of these six sequences are reproduced exactly, except A285877
a(6) and a(7).

---

## A285874 — chess games after n plies, starting without rooks

Ready now:

```
DATA  1, 20, 400, 8702, 188473, 4505624, 106770421, 2770746488, 71151220765,
      1969755500063, 53938077183794, 1568917594686774, 45144112430506979,
      1363086189786054605

%C A285874 a(0)-a(9) were confirmed and a(10)-a(13) computed with a perft program that uses
   a transposition table. Each new term was computed twice independently, with different hash
   keys and different subdivisions of the game tree, and the results agreed.
   - _<your OEIS name>_, Oct 06 2026

%E A285874 a(10)-a(13) from _<your OEIS name>_, Oct 04 2026
```

The %C line is optional but recommended: the entry's existing comment says the initial terms
came from the Chesspresso engine, so this one says where the new terms came from. In the OEIS
editor a comment is a single line; the wrapping above is only for readability.

| n | a(n) | run 1 | run 2 |
|---|---|---|---|
| 10 | 53938077183794 | about 1 min, M4 Max | 53 s, M4 Max |
| 11 | 1568917594686774 | 796 s, M4 Max | 747 s, M4 Max |
| 12 | 45144112430506979 | 3.9 h, M4 Max | 4.0 h, M4 Max |
| 13 | 1363086189786054605 | 9 h, 3 x86 servers | 17 h, 5 machines |

**a(13)** comes from two independent computations that agree exactly:

- **Run 1:** seed 1, split ply 5 (734,593 tasks), divided among three x86 servers (g++). The
  merge checked that every task was covered exactly once.
- **Run 2:** seed 2, split ply 6 (8,095,245 tasks). The Mac Studio (ARM, Apple clang) worked
  from the front of the task list. The three servers and an M5 Pro (ARM) worked from the end.
  Where the two fronts met, 193,017 tasks were computed twice, mostly on different
  architectures, and every one of them agreed.

## A285877 — chess games after n plies, starting without pawns (**correction**)

```
DATA  1, 50, 2125, 96062, 4200525, 191462298, 8509434052, 390020558283,
      17528222547068, 806635564568344, 36538968387272899

%E A285877 a(6)-a(7) corrected and a(8)-a(10) added by _<your OEIS name>_, Oct 03 2026
```

| n | OEIS now | corrected | difference |
|---|---|---|---|
| 6 | 8509434855 | **8509434052** | OEIS is 803 too high |
| 7 | 390020597683 | **390020558283** | OEIS is 39400 too high |
| 8 | — | 17528222547068 | new |
| 9 | — | 806635564568344 | new (three runs, ARM and x86) |
| 10 | — | 36538968387272899 | new (ARM and x86) |

### Discussion text to paste with this edit

> a(6) and a(7) were too large. They were computed with the Chesspresso library, which has
> a castling bug. Its castling-rights update is an if/else-if chain, so after 1. Rxh8
> (rook from h1 to h8) it clears only White's kingside right, and Black keeps its own even
> though its h8 rook was captured. It also allows castling without checking that the rook
> is present, e.g. 1. Rxh8 Nf6 2. Rh1 Ba3 3. Ra2 O-O. Chesspresso itself reproduces the old
> a(6), and a program imitating only this bug reproduces the old a(5)-a(7). The corrected
> values were confirmed with Stockfish 19 and with two other independent programs, with
> identical counts after every first move. Stockfish's "go perft 6" gives 8509434052, and
> its perft 6 after each of White's 50 first moves totals 390020558283.

### Evidence

Three independent programs agree on the corrected values move by move: Stockfish 19, and
the two programs written for this project.

1. **Two independent programs agree on a(6) = 8509434052**, including all 50
   per-first-move counts (appendix):
   - `perft.cpp`: bitboard legal-move generator, bulk counting, hash tables;
   - `naive.c`: 0x88 mailbox board, pseudo-legal moves with make-and-test legality, no
     hashing, no bulk counting, no shared code.
2. **The same engine reproduces every other existing term exactly:** all of A285873,
   A285874, A285875, A285876 and A285878, A285877 a(0)–a(5), standard chess (A048987)
   through a(9), and the Chess Programming Wiki perft test positions. That is 232 of 234
   checks, the two exceptions being A285877 a(6)–a(7). The same result was obtained on six
   machines covering ARM and x86 with Apple clang and g++.
3. **The two programs also agree on a(7) = 390020558283**, again with identical per-first-move
   counts for all 50 moves (`runs/nopawns_d7_divide_*.txt`), and `perft.cpp` gave the same
   value in two runs with independent hash keys. The a(7) gap, 39400 ≈ 803 × 49.07, is about
   one ply of no-pawn branching. That fits the 803 spurious games at ply 6 being carried forward,
   not a separate error.
4. **The old values come from a bug in Chesspresso, the library they were computed with.**
   - **It reproduces the old numbers.** Running Chesspresso from its author's repository
     with the same calls as the original computation gives 8509434855, as OEIS does.
   - **All 803 extra games start 1.Rxh8.** Chesspresso updates castling rights with an
     `else if` chain, so only one corner square is checked per move. For h1→h8 it removes
     White's kingside right and never reaches the test that should remove Black's, whose
     rook was just captured.
   - **The missing rook isn't checked.** Chesspresso doesn't verify that the rook is still
     present before castling, so it lets Black castle with no rook at ply 6, as in
     `1. Rxh8 Nf6 2. Rh1 Ba3 3. Ra2 O-O`. Its castling code then flips the f8 and h8 bits,
     so the castle creates two black rooks out of nothing.
   - **The bug only shows up when a whole file is empty.** A rook can only travel straight
     from one corner to the opposite one along an empty file, which is why the bug is
     invisible in the variants that keep their pawns.
   - **Imitating this one bug reproduces OEIS exactly.** Our program with just this
     behaviour added gives OEIS's a(5) = 191462298, a(6) = 8509434855 and
     a(7) = 390020597683. Files are in `chesspresso/` (`perft_cpbug.cpp`, per-move counts).
5. **Stockfish 19 confirms both values.** `position fen rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1`
   then `go perft 6` prints `Nodes searched: 8509434052` (35 s). `perft 6` from each of the 50
   positions after White's first move totals 390020558283 = a(7). Stockfish's per-first-move
   counts match both of our programs for all 50 moves at both depths
   (`runs/stockfish_nopawns_d6.txt`, `runs/stockfish_nopawns_d7_by_move.txt`).

## A285878 — pawns and king

Ready now:

```
DATA  1, 18, 324, 5658, 98766, 1683597, 28677387, 479763588, 8014917042,
      132060434889, 2170519816231, 35270612954472, 570918858679054, 9146467706177090

%E A285878 a(11)-a(13) from _<your OEIS name>_, Oct 03 2026
```

a(13) was confirmed on two machines: run 1 on the Mac Studio, run 2 on the M5 Pro.

## A285875 — no knights

Ready now:

```
DATA  1, 18, 324, 6572, 132640, 3030492, 68633066, 1733220521, 43321058602,
      1182486223832, 31936882857782, 924460957058722, 26485630391196432

%E A285875 a(10)-a(12) from _<your OEIS name>_, Oct 03 2026
```

a(12) was confirmed on two machines: run 1 on the Mac Studio (81 min), run 2 on the M5 Pro.

## A285873 — no queens

Ready now:

```
DATA  1, 21, 441, 9872, 220447, 5247292, 124278971, 3113440755, 77520962327,
      2024021927610, 52504446988449, 1417924743230683, 38045182847789814

%E A285873 a(10)-a(12) from _<your OEIS name>_, Oct 03 2026
```

All three new terms fit the sequence's comment that A285873(n) ≤ A048987(n) for n ≥ 7.

## A285876 — no bishops

Ready now:

```
DATA  1, 22, 484, 11248, 260904, 6434922, 158069690, 4126252938, 107097735673,
      2940365284820, 80165663096785, 2297455042438951, 65331730056646121

%E A285876 a(10)-a(12) from _<your OEIS name>_, Oct 03 2026
```

All three new terms still exceed A048987, as the sequence's comment anticipates for small n.
a(12) was confirmed on two architectures: run 1 on mathd (x86, g++) and run 2 on the M5 Pro
(ARM, Apple clang).

---

## Method note to paste into each edit's discussion

> Computed with a perft program written for this purpose. It uses a bitboard legal-move
> generator with bulk counting at the leaves. The unique positions at a split ply are
> enumerated exactly with their path multiplicities, and the count from each is computed
> using a transposition table keyed by two independent 64-bit Zobrist keys. Every new term
> was computed twice with independent hash keys and different split plies, and both results
> agree. The program reproduces A048987 through a(9), the Chess Programming Wiki perft
> positions, and every existing term of A285873-A285878 except A285877 a(6)-a(7), which are
> corrected here. The corrected a(6) and a(7) were confirmed independently with Stockfish 19
> (go perft) and with a second, unrelated program, with identical counts for every first move.

If you publish the code, for example on GitHub, add a link line:
`%H A285874 <your name>, <a href="...">C++ perft program and logs</a>` (likewise for the
other sequences).

## Appendix: no-pawn perft(6) by first move

Identical from `perft.cpp` and `naive.c` (`runs/nopawns_d6_divide_*.txt`); the total is
8509434052.

| move | count | move | count | move | count | move | count | move | count |
|---|---|---|---|---|---|---|---|---|---|
| a1a2 | 227718968 | a1a7 | 155647583 | c1a3 | 132435882 | c1g5 | 195724388 | d1d2 | 215729499 |
| a1a3 | 204412866 | a1a8 | 128748164 | c1b2 | 202996145 | c1h6 | 139652209 | d1d3 | 227070044 |
| a1a4 | 208494166 | b1a3 | 149187624 | c1d2 | 173107076 | d1a4 | 24093846 | d1d4 | 257564302 |
| a1a5 | 197549616 | b1c3 | 240344662 | c1e3 | 231217206 | d1b3 | 249719440 | d1d5 | 219033851 |
| a1a6 | 174995210 | b1d2 | 151978323 | c1f4 | 226523114 | d1c2 | 251375545 | d1d6 | 184293349 |
| d1d7 | 9211919 | d1h5 | 8071501 | f1c4 | 206083864 | g1e2 | 167058302 | h1h4 | 208964544 |
| d1d8 | 4290788 | e1e2 | 147442491 | f1d3 | 190071855 | g1f3 | 210745651 | h1h5 | 196698750 |
| d1e2 | 19202634 | e1f2 | 200216325 | f1e2 | 174227785 | g1h3 | 149062914 | h1h6 | 175291099 |
| d1f3 | 235632049 | f1a6 | 140350286 | f1g2 | 205381457 | h1h2 | 227803183 | h1h7 | 154546973 |
| d1g4 | 221237079 | f1b5 | 20617550 | f1h3 | 138168990 | h1h3 | 204491837 | h1h8 | 124951148 |
