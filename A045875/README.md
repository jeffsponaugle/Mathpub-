# A045875 — hunting digit runs in powers of two

Search tools for [OEIS A045875](https://oeis.org/A045875):

> **a(n) is the smallest m for which the decimal representation of 2^m
> contains n consecutive identical digits.**

For example a(2) = 16 because 2^16 = 65**55**36 is the first power of 2 with
two equal adjacent digits, and a(8) = 8554 because 2^8554 is the first with
eight (a run of `00000000` deep inside a 2576-digit number).

The sequence is tagged `hard, more` on OEIS: each term requires scanning the
decimal digits of *every* power of 2 up to the answer, and the answers grow
fast. When this project started, the sequence was known through
a(17) = 356,677,212 — a number with 107 million digits. The terms found here
required searching powers of two with **over a billion digits each**.

## New results from this project (September 2026)

| term | value | the run | digits of 2^m |
|---|---|---|---|
| **a(18)** | **1,973,572,199** | 18 consecutive `4`s | 594,104,431 |
| **a(19)** | **3,238,684,956** | 19 consecutive `9`s | 974,941,319 |

Both are minimality-proven (every smaller exponent searched, in logged,
checkpointed, self-verifying runs) and independently confirmed by a direct
GMP recomputation of 2^m through a separate code path. Because each run has
*exactly* the target length, each result also bounds the next term for free.
Full provenance and coverage tables are in [RESULTS.md](RESULTS.md).

**Both terms are now published in OEIS** (A045875 revision #54, 2026-09-26),
along with the run digits ('4', '9') in the companion sequence A215732.
The one remaining submittable item is the a(20) lower-bound comment — see
[SUBMISSION.md](SUBMISSION.md) for the verified OEIS state and the exact
draft line.

**Current state (as of 2026-10-03): the a(20) search is PAUSED**, all state
checkpointed (the Sparks are temporarily on another job; resume commands
are in [RESUME-A20.md](RESUME-A20.md)). No 20-run found yet. The proven,
contiguous bound is

> **a(20) > 4,033,865,312**

### a(20) coverage ledger (what has been searched)

| range (m) | status | source |
|---|---|---|
| ≤ 3,238,684,956 | clear of 20-runs | a(19) proof (the a(19) run is exactly 19; a 20-run contains a 19-run) |
| 3,238,684,957 – 3,416,000,000 | searched, empty | atom1 leg A1 (a20-a1.log) |
| 3,416,000,000 – 3,585,000,000 | searched, empty | atom2 leg B1 (a20-b1.log) |
| 3,585,000,000 – 3,725,000,000 | searched, empty | atom1 leg A2 (a20-a2.log) |
| 3,725,000,000 – 3,904,306,680 | clear of 20-runs | cleared of 19-runs during the a(19) tail (a19-r3725.log) |
| 3,904,306,680 – 4,033,865,312 | searched, empty | atom2 leg B2, paused (a20-b2.ckpt) |
| 4,033,865,312 – 4,046,000,000 | **NOT searched** | the gap to close first on resume |
| 4,046,000,000 – 4,065,471,256 | searched, empty | atom1 leg A3, paused (a20-a3.ckpt) — splices in once the gap closes |
| > 4,065,471,256 | not searched | — |

### Where the next search starts

1. **atom2 resumes `a20-b2.ckpt` (m = 4,033,865,312) with `-e 4046000000`**
   — closing the 12M gap makes coverage contiguous through 4,065,471,256.
2. **atom1 resumes `a20-a3.ckpt` (m = 4,065,471,256) with `-e 4190000000`.**
3. Beyond that, continue interleaved legs from **m = 4,190,000,000**
   (atom2's chained leg B3 in RESUME-A20.md starts there, open-ended).

Exact relaunch commands: [RESUME-A20.md](RESUME-A20.md) (also on both
atoms). The heuristic median for the actual find is around m ≈ 4.9×10⁹.

## How the search works

### The state is the decimal expansion itself

You cannot afford to recompute 2^m from scratch for each m — at half a
billion digits even one binary-to-decimal conversion takes minutes. Instead
the tools keep the full decimal expansion of 2^m as an array of base-10⁹
limbs (`uint32`, nine digits each, least significant first) and **double it
in place** once per step, scanning digits as they go.

Two structural facts make this fast:

1. **Doubling has purely local carries.** In base 10⁹, the carry into limb
   *i* when doubling is exactly `(old a[i-1] >= 5*10^8)` — independent of
   anything below limb *i−1*, because 2·a[i] is even and can never become
   10⁹−1 by itself. So there is no carry chain: every limb's new value
   depends only on itself and its lower neighbor, and the whole pass
   parallelizes perfectly.

2. **After j doublings, limb i depends only on original limbs i−j..i.**
   (Induction from fact 1.) This "dependency triangle" is what allows
   batching many doublings into one memory pass, recomputing small windows
   at intermediate steps, and splitting work across threads and machines
   with only tiny overlaps.

### Run detection: the aligned-repdigit filter

Scanning all N digits of every intermediate number would dominate the cost.
The key observation: **any run of ≥ 17 identical digits must fully contain
one aligned 9-digit limb**, i.e. a limb whose value is d·111111111 for some
digit d (worst case, a 17-run starts one digit into a limb: 8 digits + a
full limb + 0 more). Checking `v % 111111111 == 0` is one multiply and one
compare — the compiler reduces it to exactly that — fused into the doubling
loop at nearly zero cost.

Candidates are rare (about one limb in 10⁸), and each is confirmed by an
exact run-length walk into its neighboring limbs, respecting the top limb's
true digit count so leading zeros are never counted. The filter is complete
only for n ≥ 17; for smaller n the CPU tool falls back to a generic
digit-stream scan (used to verify the tool against all the small known
terms).

### Getting to the starting line, and staying honest

- **Initialization**: GMP computes 2^start once (a binary shift plus
  divide-and-conquer radix conversion — ~3 minutes for a billion digits),
  parsed into limbs. Every init and checkpoint resume is self-checked
  against an independently computed `2^m mod 10^9` by modular
  exponentiation, so a corrupt state can never silently start a search.
- **Checkpointing**: state is saved every 10 minutes (and on SIGINT/SIGTERM)
  via atomic rename, in a format shared by the CPU and GPU tools — a search
  can migrate between machines and between tools mid-flight.
- **On a hit** the tool prints a FINAL REPORT: the exponent, the digit, the
  exact run length and its position measured from both ends of the number,
  plus search statistics.

## The CPU version (`a045875.c`)

Portable C + pthreads (`-t` sets thread count). Each thread owns a
contiguous block of limbs; per step the threads double their blocks in
place (descending, so the neighbor read happens before the overwrite) and
scan for candidates in the same fused pass. macOS lacks
`pthread_barrier_t`, so synchronization is a custom sense-reversing spin
barrier — three barrier crossings per doubling.

The CPU tool is **memory-bandwidth-bound**: each doubling streams the whole
number through memory (read + write). Measured ~3,000 steps/s at 107M
digits on an Apple M-series (10 threads), ~770 steps/s on a 16-core x86,
~1,265 steps/s on the DGX Spark's 20 Grace cores. It supports any n ≥ 1 and
doubles as the reference implementation.

## The GPU version (`a045875_gpu.cu`)

CUDA, single file, no dependencies beyond the toolkit (GMP linked for
init). Three ideas take it ~50× past the CPU tool:

1. **k-batching (the big one).** Because of the dependency triangle, a
   thread that owns a span of limbs plus a K-limb *halo* below it can run
   **K = 28 doublings entirely in registers** — carries kept as 28-bit
   masks passed limb to limb — checking every intermediate value for
   repdigit limbs, and touching memory only once per 28 steps. Memory
   traffic per doubling drops 28×, converting a bandwidth-bound problem
   into a compute-bound one. The halo (28 extra limbs per 64-limb span) is
   recomputed redundantly; that 37% compute overhead buys the absence of
   any cross-thread communication inside a batch.

2. **Tile-transposed memory layout.** A thread walking its own span makes
   every warp load touch 32 different cache lines. The limb array is
   therefore stored permanently in a tile-transposed order (logical limb
   `tile·TILE + t·SPAN + i` lives at physical index `tile·TILE + i·BT + t`)
   so warp accesses are fully coalesced. Checkpoints stay canonical —
   they're permuted on load/save, preserving CPU-tool compatibility.

3. **Host-side exact confirmation.** The kernel only *flags* candidates
   `(limb, step)` into a small atomic buffer. After each batch the host
   recomputes a ±10-limb window at that exact intermediate step from the
   pre-batch state (ping-pong buffers keep it intact) using the dependency
   triangle — a few hundred operations — and measures the true run length,
   including per-step digit counts at the top of the number.

The GPU tool requires n ≥ 17 (it relies on the aligned-limb filter).
Measured on one GB10: **~38,300 steps/s at 107M digits**, ~7,300 at 590M,
~3,900 at 1.1B — throughput scales inversely with digit count, pinned near
the memory-bandwidth ceiling of the chip.

| | CPU tool | GPU tool |
|---|---|---|
| doublings per memory pass | 1 | 28 |
| parallelism | pthreads, block per thread | 1 thread per 64-limb span |
| run lengths supported | any n | n ≥ 17 |
| candidate confirmation | inline, same pass | host, from pre-batch state |
| checkpoint format | shared | shared |
| 107M digits | ~770–3,000 steps/s | ~38,300 steps/s |
| 1.1B digits | ~120 steps/s (est.) | ~3,900 steps/s |

### Verification discipline

- The GPU tool is validated **bit-exactly** against the CPU tool:
  md5-identical checkpoints after identical doubling ranges (including the
  single-step tail path, which caught a real halo-stride bug during
  development).
- Both tools reproduce a(16) = 106,892,452 and a(17) = 356,677,212 —
  digit, run length, and position — when started just below them; the CPU
  tool also reproduces a(1)..a(13) from scratch.
- Every found term is re-verified by `verify_a18.c`: a direct GMP
  recomputation of 2^m (no shared code with the search) checking digit
  count, every digit of the run, and that the run's length is exact.
- Mid-search sanity: during the a(18) hunt, an n=17 probe at 530M digits
  found a 17-run within statistical expectation, confirming the pipeline
  healthy at scales far beyond the unit tests.

### Splitting work across machines

Minimality demands that *every* exponent below the answer be searched, so
disjoint range splits waste nothing. The a(20) search runs interleaved
~12-hour legs (atom1 takes the lowest uncovered range, atom2 the next), so
the lower machine's live frontier is always a valid contiguous bound —
useful for publishing `a(20) > x` before the term itself is found. Legs
chain automatically on clean "nothing found" completions (shell scripts
keyed on exit codes: 0 = found, 1 = range exhausted, 130 = interrupted).
A bonus: ranges cleared of 19-runs during the a(19) search are free for
a(20), since a 20-run contains a 19-run — the schedule skips
[3,725M, 3,904.3M] outright.

## Building and running

macOS (Apple Silicon, Homebrew GMP):

```bash
make                    # builds the CPU tool
```

DGX Spark / Linux aarch64 (GB10 = sm_121, CUDA ≥ 12.9; works without
libgmp-dev by declaring the stable GMP ABI and linking libgmp.so.10):

```bash
make -f Makefile.spark  # builds a045875 (CPU) and a045875gpu (GPU)
```

Usage (flags are shared; `-t` is CPU-only):

```bash
# reproduce a(17) as a smoke test (~15 s init + instant find)
./a045875gpu -n 17 -s 356677200

# production search with checkpointing (resumes if the file exists)
./a045875gpu -n 20 -s 3238684957 -e 3416000000 -c a20.ckpt

# CPU tool, 10 threads, any n
./a045875 -n 12 -s 0 -t 10
```

`-s` start exponent (GMP init), `-e` optional end bound, `-c` checkpoint
file (written every `-i` seconds, default 600, and on stop/interrupt).
Progress streams to stderr on one line; results print a FINAL REPORT to
stdout.

## Files

| file | what |
|---|---|
| `a045875.c` | CPU search tool (any n) |
| `a045875_gpu.cu` | CUDA search tool (n ≥ 17) |
| `verify_a18.c` | independent GMP verifier for a found term (on the Sparks) |
| `Makefile` / `Makefile.spark` | macOS / DGX Spark builds |
|  proven results, coverage tables, run details | proven results, coverage tables, OEIS submission lines |

Search logs and checkpoints for the live campaigns live in
`/home/jbs/A045875` on atom1 (10.1.30.36) and atom2 (10.1.30.37).
| `SUBMISSION.md` | verified OEIS state and the remaining submittable item |
| `RESUME-A20.md` | paused a(20) search: checkpoints and relaunch commands |
