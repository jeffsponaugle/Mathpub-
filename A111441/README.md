# Integer sequence searches — status & results

*Last updated: 2026-10-01. a(18)–a(21) are **published in OEIS**
(A007357 revision 64, edited 2026-08-31). a(22) and a(23) are proven
here and not yet submitted — see [submission.md](submission.md).*

## A007357 — infinitary perfect numbers (`ipn_search.c`)

[OEIS A007357](https://oeis.org/A007357): N is infinitary perfect iff
σ∞(N) = 2N. Writing N as a product of *distinct* Fermi–Dirac primes
q = p^(2^k), the condition is ∏(q+1) = 2·∏q.

Until Aug 2026 OEIS listed **17 terms** (exhaustive to ~3.3×10^16,
circa 2004); it now lists **21**, the last four from this project. ~187
larger examples are known non-exhaustively (Cohen 1990, Moews & Moews
1998, Ball 2000–2007; David Moews' aliquot database at
https://djm.cc/aliquot-database/aliquot-database-.1.txt).

### RESULTS — PUBLISHED IN OEIS ✅

All four new terms were submitted and accepted; A007357 now lists 21
terms (revision 64, 2026-08-31), with:
- DATA: a(18)–a(21) appended.
- COMMENT: "No other infinitary perfect number below 10^28.
  - _Jeff Sponaugle_, Aug 30 2026"
- LINK: "Jeff Sponaugle, Factorization for a(18)-a(21)"
  (https://oeis.org/A007357/a007357.txt)
- EXTENSIONS: "a(18)-a(21) from _Jeff Sponaugle_, Aug 30 2026"
  (alongside "More terms from _Eric W. Weisstein_, Jan 27 2004")
- The "more" keyword was removed by the editors (keyword now "nonn").

| term | value | factorization | run | status |
|---|---|---|---|---|
| a(18) | 1342989770695372800 | 2^12·3^5·5^2·7^2·11·13·17·41·43·257 | 3 s, 8.2e6 nodes | published |
| a(19) | 20144846560430592000 | 2^12·3^6·5^3·7^2·11·13·17·41·43·257 | 7 s, 2.3e7 nodes | published |
| a(20) | 1441560489708423064780800 | 2^16·3^5·5^2·7^3·11·13·41·83·331·65537 | 451 s, 1.34e9 nodes | published |
| a(21) | 21623407345626345971712000 | 2^16·3^6·5^3·7^3·11·13·41·83·331·65537 | 1173 s, 3.6e9 nodes | published |

a(18) = Moews & Moews 1998 example; a(19) = Cohen 1990; a(20) also
appears in the Moews database (FD components
3·7·11·13·25·41·49·81·83·331·65536·65537).

The searches prove the *gaps are empty*: no infinitary perfect number
exists between consecutive terms from a(17) up through a(21), nor in
(a(21), 10^28] — intervals that had never been exhaustively swept.

Validation: the program reproduces the 17 known OEIS terms exactly
(X = a(17): 1.1 s, 2.0e6 nodes) before any new claims.

### Certified bounds

- **Published**: no infinitary perfect number exists in (a(21), 10^28]
  (sweep 2026-08-30: 3.1 h, 3.28e10 nodes, warnings 0; `sweep_1e28.log`).
  So a(22) > 10^28, as recorded in the OEIS comment.

### ✅ a(23) = 9224733293334342541129236480000 — PROVEN (2026-10-01)

Exhaustive sweep of (a(22), 9.22×10^30] completed on the 96-core +
32-core pair. X = 9224733293334342541129236480000, 311775-task depth-5
frontier, sharded 2 ways (`-M 2 -K 0/1`).

- a(23) = 2^17·3^6·5^4·7^2·11·41·79·83·157·313·331·65537
  (FD components 2·9·11·41·49·79·81·83·157·313·331·625·65536·65537) —
  a Moews & Moews (1998) example, now proven to be the 23rd term. It is
  1.5×a(22) and shares a(22)'s large components (79, 157, 313, 331,
  625, 65537); the two differ only in the 2- and 3-parts.
- Certificate: shard K=1 complete on .26 (34386 s, 1.90e11 nodes,
  warnings 0). Shard K=0's tail was split between .23's own run and
  five prefix-restricted runs on .26 (`-P` task lists); the union of
  all checkpoints leaves **0 of 311775 tasks uncovered**. Every run
  reports `wn 0` throughout; the one deliberately interrupted run's
  completed tasks are counted only from its checkpoint bitmap.
- Solutions found across all checkpoints: exactly 23, equal to
  a(1)..a(23). Nothing new below 9.22×10^30.
- Total ≈ 4.6e11 nodes across all runs (includes the ~1845 tasks
  done twice before the overlap was caught). a(24) > 9.22×10^30; the next known example is
  9.22×10^31 (within the code's audited 10^32 arithmetic ceiling).

### ✅ a(22) = 6149822195556228360752824320000 — PROVEN (2026-09-30)

Exhaustive sweep of (10^28, 6.15×10^30] completed across two machines
(96-core + 32-core, checkpoint-partitioned). Certificates: warnings 0
on both final runs; all 308792 subtree tasks complete in both
checkpoints; 22 distinct solutions = a(1)..a(22) exactly.
a(22) = 2^18·3^5·5^4·7^2·11·41·79·83·157·313·331·65537 (a Moews &
Moews 1998 example, now proven to be the 22nd term; earlier notes here
wrongly credited it to Ball). a(23) > 6.15×10^30; run toward the
next known example at 9.22×10^30 in progress (`a23_*.log`).

Engineering log for the final phase (v12→v15): work-stealing queue for
subtree-level parallelism; skipped redundant per-candidate primality
proofs; **batched interval factoring** (sieve small primes across each
scan block instead of per-candidate Pollard rho) — turned subtrees that
had absorbed 500+ core-hours into 80-second work; children reuse the
precomputed factorization of v+1. Every revision gated on exact
node-count reproduction of the a(20) run (1340384551 nodes) plus the
17 published terms.

### Superseded frontier notes — a(22)

The smallest known infinitary perfect number above 10^28 is (Moews &
Moews 1998)

```
6149822195556228360752824320000  ≈ 6.15×10^30
```

so **a(22) ≤ 6149822195556228360752824320000**, and confirming it
requires exhaustively emptying (10^28, 6.15×10^30]. That run
(~3.4×10^11 nodes projected from the observed nodes ≈ X^0.365 scaling)
was started 2026-08-30, reached ~17% (5.9×10^10 nodes, warnings 0) on
the M4 Max, and was **checkpointed and migrated to another server**
on 2026-08-30 (`ipn_a22.ckpt`, 35062/308792 tasks done at handoff):

```sh
./ipn_search_v2 -x 6149822195556228360752824320000 -t <cores> -c ipn_a22.ckpt
```

**2026-09-29**: run resumed on `mathd` (10.1.30.23, 96 cores, ~4.0M
nodes/s — 2.5× the M4 Max) and then split across a second box
(10.1.30.26, 32 cores) by checkpoint partition: remaining tasks divided
3:1 interleaved, each machine's checkpoint pre-marking the other's
share as done; exhaustiveness = union of both runs' completions, both
ending with warnings 0. Fixed en route: ckpt_apply read the done-bitmap
line through a 64 KB buffer, truncating resumes past 262k tasks
(caused only redone work, never missed work; buffer now sized to the
task count).

Beyond a(22), the next known examples (Moews database) are close
together and then jump:
- 9224733293334342541129236480000 ≈ 9.22×10^30 — only 1.5× above the
  a(22) candidate, so the incremental a(23) run is cheap once a(22)
  is done;
- 92247332933343425411292364800000 ≈ 9.22×10^31 (a(24) candidate) —
  still within the code's audited 10^32 arithmetic ceiling;
- then ≈ 3.76×10^33 — **beyond the current 10^32 ceiling**; extending
  it needs a fresh overflow audit (mainly compute_ab guard, den/q1
  bounds, and record verification widths).

### History (2026-08-29 → 08-30 session)

Single overnight session: built the branch-and-bound searcher,
validated against the 17 published terms, proved a(18)–a(21), swept to
10^28, upgraded arithmetic to a 10^32 ceiling (`ipn_search_v2`:
256-bit mul/divmod in the two-component solver, 128-bit Miller–Rabin
FD tests for candidates > 2^64, bit-identical regression at a(20)
scale), and started the a(22) sweep.

### Algorithm (soundness sketch)

State = chosen components (strictly increasing), exact needed ratio
r = a/b kept as a signed prime-exponent list (a | 2∏q ⇒ exact in
integers), budget B = X/∏q. Rules, each proven sound:
- a==b ⇒ solution; a<b ⇒ dead.
- Closure: (a−b) | b ⇒ q* = b/(a−b) is the unique 1-component finish.
- Every prime P | b forces a future component P^(2^k) ∈ (f, B]; product
  of these minima over distinct P must fit in B.
- Every prime P | a forces a future q ≡ −1 (mod P) in (f, B].
- qmin = max(f+1, ⌈b/(a−b)⌉); qmin³ > B ⇒ ≤2 components remain ⇒ solved
  exactly (closure + divisor identity ((a−b)q₁−b)((a−b)q₂−b) = ab).
- Else scan the smallest next component up to Q ≤ ⌊log_qmin B⌋·a/(a−b).
- Any candidate component that would exceed 2^64 is counted in the
  `wn` (warn) statistic rather than silently skipped: **wn = 0 in the
  final line is the exhaustiveness certificate.**

Parallelism: deterministic BFS frontier of subtree tasks, worker pool,
checkpoint = done-bitmap + solutions every 30 s (`-c file`, resume by
rerunning the same command).

```sh
cc -O3 -I/opt/homebrew/include -L/opt/homebrew/lib -o ipn_search ipn_search.c -lprimesieve -lpthread
./ipn_search -x 21623407345626345971712000 -t 14   # a(21) run
```

## A111441 — sum of squares of first k primes divisible by k

Earlier session work, kept for reference. a(17) > 8×10^15 (Dyson);
linear prime-sieve problem, ~200 days to the frontier on this machine.
- `a111441.c` — single-thread reference (verifies a(1)-a(11) in ~5 s).
- `a111441_mt.c` — multithreaded two-pass searcher, ~700M primes/s on
  12 cores, checkpointed (`test.ckpt` holds progress to 5.5×10^12).
Feasibility notes in git history of this README; parked in favor of
A007357, which has proven far more productive.
