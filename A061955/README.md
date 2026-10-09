# Extending A061955 and the "n divides concatenation of 1..n" families

Searching for new terms of OEIS [A061955](https://oeis.org/A061955) and its
cross-referenced families, using an O(b·log²n)-per-n evaluator instead of the O(n)
direct computation used for the published data (most of which stops at 10^7, 2011).

New terms and draft OEIS edits: **[submissions.md](submissions.md)**; b-files in `bfiles/`.

**Result (Oct 4 2026):** all 96 reversed-digit sequences searched completely to 10^11, giving
1,215 new terms. Every new term was confirmed by direct evaluation, and the search
reproduces every published term (including the four existing OEIS b-files) exactly.
56 sequences now need a b-file because their full list no longer fits in the data field.

## The sequences

Each family has 24 members, base b = 2..25 (A-number = first + b − 2).
"Left" concatenation writes n … 2 1, "right" writes 1 2 … n.

| tag | family | concatenation | digits of each k | published search limit |
|---|---|---|---|---|
| Rn | A029447–A029470 | right | normal | (no `more`; full lists) |
| Ln | A029471–A029494 | left | normal | 3·10^10 (2024); A029471 10^11 (Aug 2026) |
| Rk | A029495–A029518 | right | reversed, zeros kept | 10^7 (2011); A029495 10^8 (GPU, Apr 2026) |
| Ld | A029519–A029542 | left | reversed, least significant zeros dropped | 10^7 (2011) |
| Rd | A061931–A061954 | right | reversed, zeros dropped | 10^7 (2011) |
| Lk | A061955–A061978 | left | reversed, zeros kept | 7·10^6–10^7 (2011) |

The four reversed-digit families (96 sequences) are the targets; the normal-digit
families are only used to cross-check the machinery against well-studied data.

## Technique

Write V(n) for the concatenation as an integer. We need V(n) mod n for every n, and
the modulus changes with n, so nothing carries over from one n to the next. Direct
evaluation is O(n) per n, which caps a search near 10^7.

* **String monoid.** A digit string is kept mod q as (V, X) = (value, b^length).
  Concatenation is (V1,X1)·(V2,X2) = (V1·X2 + V2, X1·X2).
* **Blocks and aligned ranges.** Words are grouped by digit count D. In an aligned range
  k = P + m, 0 ≤ m < b^j, the values are val(P) + s·v_j(m). Here v_j is the j-digit
  reversal (s = b^(D−j)), or m itself for normal digits. The lengths depend on m only,
  except for m = 0 in the zeros-dropped families, where the count of trailing zeros
  comes from P.
* **Digit recurrence.** The words m = 1..b^j−1 (the "interior" I_j) are summarized by
  A = Σ b^position, B = Σ v_j(m)·b^position, X = b^length and G = A + X. Splitting m on
  its top digit gives a closed recurrence using the powers of μ = X·b^(word length):
  `X' = μ^(b−1)X, G' = P·G, B' = P(bB+X) + QG − bX'` (left; right is `B' = bPB + SG`),
  with P = Σμ^t, Q = Σtμ^t, S = Σ_{t≥1} t·μ^(b−1−t). So a complete block of D-digit
  words costs O(b·D) multiplications, and V(n) mod n costs about (b+3)·D²/2, with
  D = log_b n.
* **Zeros kept.** All words in a block have the same length, so the full sets
  m = 0..b^j−1 obey a cheaper recurrence in y_j = b^(D·b^j). In base 2 that is
  3 multiplications per step (G' = G(1+y), R' = 2R(1+y) + yG, y' = y²). Zeros dropped
  costs 4 per step in base 2, with μ_j = 2·(2^(D−1))^(2^j) taken from a squaring chain.
* **Last block.** The final partial block [b^(D−1), n] is cut along the digits of n+1
  into aligned ranges. At each digit position the ranges t ≥ 1 are summed in closed
  form from the same powers; the range starting at Q_j itself (t = 0) is added
  separately.
* **Arithmetic.** Odd moduli use 64-bit Montgomery multiplication. Even n = 2^s·q are
  tested mod 2^64 first, which is cheap and rejects most of them, then mod q. Four
  moduli are evaluated in lockstep for instruction-level parallelism. Base-2 kernels
  and general-base kernels (with a shallow product tree for the powers) are both
  hand-specialized.
* **AVX-512 floating-point engine (default on x86).** For n < 2^43 residues are held as
  exact integers in doubles: h = a·b, l = fma(a, b, −h) is the exact low part,
  q = round(h·(1/n)), r = fma(−q, n, h) + l. For |a|, |b| ≤ 4n the quotient is never off
  by more than 1, so |r| < 0.52n and every intermediate is an exact integer below 2^53;
  sums are re-reduced before they can exceed 4n. Eight moduli go in one zmm register
  (NEON: 2), with per-lane masks where the digits of n differ in the last block. It runs
  2.6–6× faster than the 64-bit Montgomery integer engine on the Xeons (which do only one
  64×64→128 multiply per cycle); `-x 0` selects the integer engine.
* **Prefilters** from the last word. Left: V ≡ 1 (mod b), so gcd(n, b) = 1.
  Right, reversed: V ≡ rev(n) (mod b^len(n)), so gcd(n, b^len(n)) must divide rev(n).

Cost per n at 10^10 is about 0.5–2 µs on one Apple performance core, depending on
family and base, versus about 10^10 operations for the direct method.

## Files

| file | purpose |
|---|---|
| `concat.c` | the search tool: all 6 families, bases 2..36, multithreaded; integer and AVX-512/NEON engines |
| `brute2.c` | independent reference from the definition; `-v N` verifies one large N on all cores |
| `run_jobs.sh` | runs a queue file of `A-number LO HI` jobs; resumable, skips finished jobs (`done_LO_HI` markers) |
| `run_pool.sh` | runs jobs from a pool shared by all boxes (atomic claims on mathb) |
| `work.sh`, `restart.sh` | a box's whole workload (own queue, then the pool), and how to (re)start it |
| `jobs_base2_1e12.txt` | the shared pool: base-2 sequences over [1e11, 1e12) |
| `validate.sh` | checks a build: brute2 vs concat (`res`, `resl`, `resv`) for all 144 sequences, n < N |
| `collect.py` | pulls results from all machines, checks them against OEIS, verifies new terms, writes submissions.md |
| `jobs_<machine>_all.txt` | the queue each x86 box is running |
| `data/oeis/` | the 144 OEIS entries (text format) the results are compared against |
| `data/oeis_bfiles/` | the 4 existing OEIS b-files (A029512, A029515, A029518, A061951; Lars Blomberg, terms below 10^7) |
| `bfiles/bNNNNNN.txt` | new b-files for all 96 reversed-digit sequences: every term below 10^11, OEIS format |
| `results/final_collect.txt` | the last complete collection (all machines, Oct 4 19:51); `collect.py --snapshot` rebuilds from it |
| `results/verified.json` | brute-force verification record for every new term |
| `a061955.c`, `brute.c`, `check_small.py` | first A061955-only tools (special-purpose searcher, O(n) C and big-integer Python references) |
| `gpu_concat.metal` | Metal kernel prototype for the base-2 families (not yet used) |

```
cc -O3 -mcpu=native -ffp-contract=off -o concat concat.c -lpthread     # Apple silicon / ARM
gcc -O3 -march=native -ffp-contract=off -o concat concat.c -pthread    # x86-64 Linux (AVX-512 if available)
cc -O3 -o brute2 brute2.c -lpthread

./concat search A029519 1 10000000000 -t 14 -o runs/A029519   # SEQ LO HI
./concat search Ld:2 10000000000 100000000000 -o runs/A029519   # TAG:BASE form, next decade
./concat res A061955 2323794751 2323794752                   # V(n) mod n
./brute2 Lk 2 -v 2323794751                                  # independent check
```
`search` appends hits to `<dir>/hits_<A-number>.txt` and keeps
`hits_<A-number>.txt.progress` = the contiguous range finished so far (restart from there).
On a new machine: build, run `./validate.sh` (must print 144/144), then
`setsid nohup sh run_jobs.sh JOBFILE THREADS &`. Collect with `./collect.py --verify`.

## Validation

* `brute2` vs `concat` (scalar and 4-lane paths): every residue V(n) mod n equal for all
  n ≤ 1500, all 6 families × bases 2..25 (144 sequences), re-run after every kernel change.
* 7,800+ structured spot checks up to 3·10^6 (b^k ± 1, c·b^k, many trailing zeros, random).
* Lane path vs scalar path: residues equal for 400 consecutive n near 10^10, all 96 sequences.
* Search mode vs OEIS: hits on [1, 2·10^6] equal the published terms for all 144 sequences.
* Every new term is re-checked by direct O(n) evaluation (`brute2 -v`, `brute`).
* The same validation (`validate.sh`, now including the AVX-512/NEON path) passed on every
  machine before it ran searches.
* Every new term is confirmed by `brute2 -v` (results/verified.json); a sample was also
  re-checked through brute2's plain per-word path. `collect.py` reports any contradiction
  with OEIS (a published term missed inside the searched range, or an extra hit below the
  largest published term); there have been none.

## Machines and runs

All searching runs on the three Linux x86 boxes (user `jbs`, `~/src/math/A061955`, gcc 13,
`-O3 -march=native`, AVX-512 engine, validated 144/144 before use). Each box runs
`work.sh MACHINE THREADS` at `nice 10`, which does two things in order:

1. `run_jobs.sh jobs_<machine>_all.txt` -- the box's own queue (static split): the rest of
   the 92 general-base sequences to 10^10, then to 10^11 (mathb also runs the base-2
   sequences to 10^11 first).
2. `run_pool.sh jobs_base2_1e12.txt` -- the shared base-2 pool (disabled, see the plan
   section): A061955, A029519,
   A029495, A061931 over [10^11, 10^12) in 36 pieces of 10^11. A box claims the next piece by
   creating `~/src/math/A061955/pool/claims/<A>_<LO>_<HI>` on mathb (atomic `mkdir`;
   mathd and mathg do it over ssh), so whichever box frees up first takes the next piece.

| machine | hardware | own queue |
|---|---|---|
| mathb `10.1.30.21` | 2x Xeon Platinum 8268, 48c/96t | base-2 to 1e11; 18 general [1, 1e10); 40 general [1e10, 1e11) |
| mathd `10.1.30.23` | 2x Xeon Platinum 8168, 48c/96t | 53 general [1, 1e10); 37 general [1e10, 1e11) |
| mathg `10.1.30.26` | 2x Xeon Gold 5315Y, 16c/32t | 20 general [1, 1e10); 15 general [1e10, 1e11) |

Finished jobs leave `runs/<A>/done_<LO>_<HI>`; `runs/jobs.log` records starts, ends and
pool claims. To stop or restart a box: `sh restart.sh MACHINE THREADS` (restarts `work.sh`;
the interrupted job resumes from its progress file, finished jobs are skipped).
To install a new `concat`: build it elsewhere, run `validate.sh`, then `mv` it over
`concat` (running jobs keep the old binary until their next job).

Both Macs are free. Mac 1 (M4 Max) ran A061955 and A029519 to 10^10; Mac 2 (M5 Pro,
`jbs@10.1.1.43`) ran A029495, A061931 and A029520 to 10^10. Their output stays in their
own `runs/` directories and `collect.py` still reads it.

## Cost model and plan

Work grows like N·log²N per sequence. Engine history on the Xeons (µs per n, one
thread, while the boxes were busy, so only the ratios matter):

| sequence | 64-bit Montgomery | AVX-512 FP | + even n to lanes | + no divisions |
|---|---|---|---|---|
| Lk:2 (A061955) | 3.09 | 1.21 | 1.21 | 1.08 |
| Ld:3 (A029520) | 5.45 | 3.13 | 1.61 | 1.31 |
| Lk:10 (A061963) | 2.13 | 0.51 | 0.51 | 0.38 |
| Rk:25 (A029518) | 3.84 | 2.37 | 1.23 | 0.94 |

Overall 3–5x over the integer engine. The current queues ([1, 10^10) for the rest of
the 92 general-base sequences, base 2 to 10^11, then all 92 general-base sequences to
10^11) were estimated at ~61 h with the integer engine; at current speeds they should
finish in roughly 8–10 h. A [1, 10^10) job now takes 1.5–3 min on mathd; a
[10^10, 10^11) job ~10–15 min.

The base-2 stage to 10^12 (shared pool, see above) was not started: on Oct 4 we chose
to stop once every sequence is complete to 10^11. The pool file is parked as
`jobs_base2_1e12.txt.disabled` on each box. To run it later, rename it back to
`jobs_base2_1e12.txt` and start `sh run_pool.sh jobs_base2_1e12.txt THREADS` on each box
(or `sh restart.sh MACHINE THREADS`, which skips the finished queue). Estimate: ~18
mathb-hours, about 8 h on all three Xeons.

Check progress any time with `./collect.py` (add `--verify` to confirm new terms); it
rewrites submissions.md and `bfiles/`. It caches each machine's report in `results/raw/` and
won't rewrite anything from partial data. `./collect.py --snapshot results/final_collect.txt`
rebuilds from the final complete run without contacting any machine.

## GPU (Metal)

Not used. A kernel for the base-2 families (`gpu_concat.metal`, runtime-compiled, no host
program yet) was started before the x86 machines became available. The AVX-512 engine
on the Xeons made it a lower priority: Apple GPUs have no FP64 and only 32-bit integer
multipliers, so a 64-bit modular multiply costs ~16–20 multiply instructions, and the
estimated gain was ~2–4x over each Mac's own CPU, against the Xeons' combined ~4x.
