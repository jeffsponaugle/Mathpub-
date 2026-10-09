# Extending OEIS A285874 and its siblings

This project computes the number of chess games that last exactly n plies (*perft*) from
six modified starting positions. Proposed OEIS edits, including a correction to A285877,
are in **[submission.md](submission.md)**.

| Sequence | Start position | FEN |
|---|---|---|
| [A285874](https://oeis.org/A285874) | no rooks | `1nbqkbn1/pppppppp/8/8/8/8/PPPPPPPP/1NBQKBN1 w - - 0 1` |
| [A285873](https://oeis.org/A285873) | no queens | `rnb1kbnr/pppppppp/8/8/8/8/PPPPPPPP/RNB1KBNR w KQkq - 0 1` |
| [A285875](https://oeis.org/A285875) | no knights | `r1bqkb1r/pppppppp/8/8/8/8/PPPPPPPP/R1BQKB1R w KQkq - 0 1` |
| [A285876](https://oeis.org/A285876) | no bishops | `rn1qk1nr/pppppppp/8/8/8/8/PPPPPPPP/RN1QK1NR w KQkq - 0 1` |
| [A285877](https://oeis.org/A285877) | no pawns | `rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1` |
| [A285878](https://oeis.org/A285878) | pawns and king | `4k3/pppppppp/8/8/8/8/PPPPPPPP/4K3 w - - 0 1` |

## Results (Oct 04 2026)

**Verified** means two independent runs agree; see `runs/results.txt`. **One run** means
the confirming run is queued.

| Sequence | n | a(n) | Status |
|---|---|---|---|
| A285874 | 10 | 53938077183794 | verified |
| A285874 | 11 | 1568917594686774 | verified |
| A285874 | 12 | 45144112430506979 | verified |
| A285874 | 13 | 1363086189786054605 | verified: run 1 on 3 x86 servers; run 2 on the Mac Studio, 3 servers and M5 Pro; 193,017 overlapping tasks agree |
| A285877 | 6 | 8509434052 | verified **correction** (OEIS: 8509434855) |
| A285877 | 7 | 390020558283 | verified **correction** (OEIS: 390020597683) |
| A285877 | 8 | 17528222547068 | verified |
| A285877 | 9 | 806635564568344 | verified (3 runs, ARM and x86) |
| A285878 | 11 | 35270612954472 | verified |
| A285878 | 12 | 570918858679054 | verified |
| A285878 | 13 | 9146467706177090 | verified |
| A285875 | 10 | 31936882857782 | verified |
| A285875 | 11 | 924460957058722 | verified |
| A285875 | 12 | 26485630391196432 | verified |
| A285873 | 10 | 52504446988449 | verified |
| A285873 | 11 | 1417924743230683 | verified |
| A285876 | 10 | 80165663096785 | verified |
| A285876 | 11 | 2297455042438951 | verified |
| A285873 | 12 | 38045182847789814 | verified |
| A285876 | 12 | 65331730056646121 | verified (x86 and ARM) |
| A285877 | 10 | 36538968387272899 | verified (x86 and ARM) |

Every existing OEIS term of the six sequences is reproduced exactly, except A285877 a(6)
and a(7). The evidence for that correction is in [submission.md](submission.md#a285877--chess-games-after-n-plies-starting-without-pawns-correction).
All computations are finished. Every term in the table is verified.

**Why the old A285877 values are wrong.** They were computed with the Chesspresso library.
Running Chesspresso itself (`chesspresso/`) reproduces them and pins down the bug:
- **Rights not cleared:** Chesspresso updates castling rights with an `else if` chain. A rook
  capturing from one corner to the opposite corner (1.Rxh8) therefore clears only the
  mover's right, and Black keeps kingside castling after losing its h8 rook.
- **Rook not checked:** it also castles without checking that the rook is still there.
- **Confirmed:** our engine with exactly that behaviour added reproduces OEIS's
  a(5)–a(7). The bug needs an empty a- or h-file, which is why only the no-pawn variant is
  affected.

## How the counts are computed (`perft.cpp`)

- **Move generation.** A bitboard legal-move generator using magic bitboards and pin and
  check masks. Nodes one ply from the end are bulk-counted: their legal moves are counted,
  not made.
- **Split ply.** The unique positions at a split ply k are enumerated exactly
  (full-position equality) with their path multiplicities. Each becomes a task computing
  perft(n−k), and the answer is Σ multiplicity × count, summed in 128-bit arithmetic.
- **Hash table.**
  - A lockless transposition table shared by all threads. An entry is located by one
    64-bit Zobrist key and verified by a second, independent 64-bit key plus the depth.
    A single 64-bit key is not safe at this scale: about 10¹³ probes would give hundreds of
    false hits.
  - An XOR lock rejects reads that catch an entry mid-write.
  - Buckets match the cache line: 128 bytes on Apple Silicon, 64 bytes on x86.
- **Huge pages on Linux.** The table is touched with writes right after allocation so it
  lands on 2 MB huge pages. If the first touch is a read, Linux maps its shared zero page
  and later splits it into 4 KB pages. With huge pages, the same job ran 1.63× faster on a
  200 GB table.
- **Checkpoints.** Every run writes a checkpoint (`--ckpt`) and resumes from it. A final
  line cut short by a crash is ignored and trimmed, never parsed as a shorter number.
- **Distributed runs.**
  - `--tasks FILE` or `--range A:B` makes a process compute only part of a run.
    `--reverse` walks its list from the end.
  - `--merge` combines the checkpoints of all parts. It re-verifies every entry against
    the task list, requires every task to be present, and requires tasks computed by more
    than one machine to agree.
  - Task cost is far from uniform. Tasks are ordered by White's first move, and the
    pawn moves come last. Those tasks (after moves like 1.e4 and 1.d4) cost about 3×
    more than earlier ones. Machines that start early in the list therefore look
    fast at first. Balance on measured rates, not on task counts.
- **Scope.** As in standard perft, the halfmove clock, the fullmove number and repetition
  are ignored. An en passant square is hashed only when a capture is pseudo-legal.

## Verification

1. **Self-test** (`./perft --test --max 3000000000000`, 234 checks): perft of the Chess
   Programming Wiki positions 1–6 including mirrored position 4, A048987 through a(9), and
   every known term of the six sequences. Each depth is computed three ways where
   affordable:
   - a reference perft that makes every move and cross-checks the generator against the
     counter, incremental hash keys against recomputed ones, and move legality;
   - bulk counting without the hash table;
   - bulk counting with it.

   Every build was checked on every machine that used it for results. The full suite
   (to 3·10¹² nodes) passed on all six machines (ARM and x86, Apple clang and g++ 13),
   each giving 232 of 234, the two exceptions being A285877 a(6)–a(7). Build
   `f938ceff9324` also passed the full suite on mathb, and `7de4dcb9ad85` on the M5 Pro. Both later x86 builds
   (`f938ceff9324`, `7de4dcb9ad85`) passed a reduced suite (to 10¹⁰) on each server, with
   identical results everywhere.
2. **Two independent runs per term**: different Zobrist keys (`--seed`), different split
   plies (so different task lists), and different table sizes, often on different machines.
3. **`naive.c`**: a second implementation sharing no code with `perft.cpp`. It uses a 0x88
   mailbox board, pseudo-legal moves with make-and-test legality, and no hashing. It gave
   identical per-first-move counts for A285877 a(6) and a(7) (`runs/nopawns_d*_divide_*.txt`).
4. **Stockfish 19** (third party, `brew install stockfish` on the Mac Studio) gives the same
   A285877 a(6) and a(7), with the same count for every first move
   (`runs/stockfish_nopawns_*`; `sf_move.sh` runs one first move).
5. **Distributed-merge tests**: the merge reproduced A285874 a(9) from three overlapping
   parts, one of them computed in reverse. It refused an incomplete set, flagged a
   deliberately altered count (`CONFLICT`), and a part interrupted mid-run resumed
   correctly.

## Machines

| Name | Hardware | RAM | Table | Role |
|---|---|---|---|---|
| M4 Max (laptop) | MacBook Pro, 14 cores | 36 GB | 8 GB | A285874 a(10)–a(12), early sibling terms; now kept idle |
| M5 Pro | MacBook Pro, 18 cores | 64 GB | 32 GB | sibling run-2 queue (`run_queue.sh`, seeds 2), then A285873 a(12) run 1; now helps the A285874 a(13) verify run from the end of its task list |
| Mac Studio | M2 Ultra, 24 cores | 64 GB | 24–32 GB | sibling confirmations, naive a(7) check, then the A285874 a(13) verify run (`studio_chain.sh`) |
| mathd | 2× Xeon Platinum 8168, 48 cores / 96 threads | 754 GB | 600 GB | A285874 a(13) run 1: its share, then phase 2 (`phase2.sh`); runs the merge (`merge_a13.sh`); then A285876 a(12) run 1 |
| mathb | 2× Xeon Platinum 8268, 48 cores / 96 threads | 754 GB | 600 GB | A285874 a(13) run 1: its share, then phase 2; then A285877 a(10) run 1 |
| mathg | 2× Xeon Gold 5315Y, 16 cores / 32 threads | 1.5 TB | 1.2 TB | A285874 a(13) run 1: a reduced share (it runs at about 0.3× mathd's task rate), then phase 2 |

Approximate run times for A285874: a(10) about 1 min, a(11) about 13 min and a(12) about
4 h on the M4 Max with an 8 GB table.

## Program versions

Every log starts with `build: <source hash> <compiler> <arch>`.

| Source hash | Change | Used for |
|---|---|---|
| `71656dc13d5c` (`runs/perft_v1.cpp`) | first version; power-of-two tables | M4 Max runs; M5 Pro runs before 09:45 Oct 03 |
| `7c75e38c1115` | portable (x86 and ARM), any table size, build stamp | Mac Studio; first x86 validation |
| `452923c2df94` | huge-page fix on Linux | validation on mathb |
| `f938ceff9324` | task lists, `--range`, `--reverse`, `--merge` | A285874 a(13) shares (`perft3`) |
| `7de4dcb9ad85` | hardened checkpoint reading | current `perft.cpp`; the merge and phase-2 runs on x86; M5 Pro runs from 09:45 Oct 03 |

The move generator and counting logic are unchanged since v1. Later versions change only
memory handling, the driver and checkpoint I/O.

## Files

- `perft.cpp`: the engine. `naive.c`: the independent cross-check.
- `build.sh`: builds both, stamping the source hash, compiler and architecture.
- `extend.sh`: runs a job file (`variant depth split1 split2 [split3]`) for the seeds in
  `$SEEDS`.
- `run_queue.sh`: one machine's share of the sibling queue.
- `collect.sh`: pulls logs from all machines into `runs/<machine>/` and appends agreeing
  results to `runs/results.txt`.
- `assign.py`: splits a run's remaining tasks into interleaved blocks weighted by machine
  speed. `rebalance.py`: moves the unfinished tail of one machine's list to others' phase-2
  lists, which were later re-split three ways from measured rates (originals kept in `runs/`).
  `phase2.sh`: runs a machine's phase-2 list after its share finishes. `after_a13.sh`: starts a
  server's next job (`jobs_after_<machine>.txt`) once its a(13) work is done.
- `merge_a13.sh`: mathd's script that retries the merge every 10 minutes until the three
  checkpoint files cover every task.
- `studio_chain.sh`, `verify_a13.sh`: the Mac Studio's queue and the a(13) verify run.
  `coord_a13v.sh` (on the laptop): coordinates the verify run. The Studio works from the
  front of the task list. A back group (mathd, mathb, mathg and the M5 Pro, lists in
  `tasks_v13_*.txt`, rebalanced once from measured rates; the first lists are in
  `runs/v13/lists_round1/`) works from the end toward it. The script merges all five
  checkpoints on mathd when the fronts meet or the back group finishes, then stops everything.
- `jobs_*.txt`: job lists. `run_a12.sh`: the script that computed A285874 a(12).
- `runs/`: logs, checkpoints, self-test logs, `results.txt`, and A285877 per-move
  breakdowns (`nopawns_d6_divide_*.txt`).

## Usage

```
./build.sh
./perft --test --max 3000000000000 --hash 8192                     # validation suite
./perft --variant norooks --depth 11 --hash 8192 --split 5 --seed 1 --ckpt runs/x.ckpt
./perft --variant nopawns --depth 6 --divide                        # per-move breakdown
./perft --fen "<FEN>" --depth 7                                     # any position
# a distributed run: each machine computes its share, then the checkpoints are merged
python3 assign.py 734593 done.ckpt 2000 mathd=1.0 mathb=1.05 mathg=0.5
./perft --variant norooks --depth 13 --split 5 --seed 1 --tasks tasks_mathd.txt --ckpt part.ckpt
./perft --variant norooks --depth 13 --split 5 --seed 1 --merge a.ckpt,b.ckpt,c.ckpt
```

Variants: `std noqueens norooks noknights nobishops nopawns pawnsking`. Choose a `--hash`
size (MB) that fits in free RAM. If the table spills into compressed or swapped memory, the
run slows down by about 10×. On Linux, run large tables under `numactl --interleave=all`.
