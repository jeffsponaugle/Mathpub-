# prime_runs_ps

Search for **runs of consecutive primes that share the same digit sum**, at
scale, using the [primesieve](https://github.com/kimwalisch/primesieve)
library as the prime source.

---

## What it searches for

Take the primes in increasing order and compute the digit sum of each one
(base 10 by default):

```
  ...  22193  22229  22247  22279  ...
        17     17     17     28
```

A **run of length L** is a block of `L` *consecutive* primes — no primes
skipped — whose digit sums are all equal. The block above contains a run of
length 3 with digit sum 17, starting at 22193, which is the smallest such run
in base 10.

The program scans a half-open interval `[start, end)` and reports either the
first run of the requested length or, with `--all`, every one it finds.

A few details that matter when interpreting the output:

* **Every window counts.** A *maximal* run of length `M` contains `M - L + 1`
  distinct windows of length `L`, and `--all` reports each of them. So a
  maximal run of 5 primes reported at `-n 3` yields three lines.
* **Ownership by first member.** A run is attributed to the interval
  containing its *first* prime. This makes the output independent of the
  thread count and of where a resumed job happened to stop.
* **Other bases.** `--base B` computes digit sums in base `B` (2..36). The
  cheap pre-filter is congruence mod `B-1`, since digit sum ≡ n (mod B-1) in
  any base; the expensive digit sum is only computed when that test passes.
* **Emirps.** With `--emirps`, every member of the run must also be an emirp:
  a prime whose base-`B` digit reversal is a *different* prime. Palindromic
  primes are excluded by the usual convention. The reversal test is only
  applied once a window already shares a digit sum, since it is by far the
  most expensive check in the pipeline.

### How the digit sums are computed

Primes arrive in increasing order, so digit sums are never computed from
scratch. The program keeps the low `k` digits of the current prime in a
window (`k` chosen so that `base^k` fits in 1 MB — 10⁶ for base 10) plus the
digit sum of everything above the window. Sliding the window forward by the
prime gap costs one byte load from an L2-resident table and one add; the
upper half only has to be recomputed when the window carries, which at base
10 happens about once every 25,000 primes.

This gives the exact digit sum of every prime for roughly the cost of the
mod-9 pre-filter it replaces, so no pre-filter is needed at all. Measured on
a pre-materialized prime stream near 10¹², it is about 1.7× faster than
computing digit sums lazily from a 400 KB lookup table.

### Statistics

Alongside the search, the program builds a histogram of **maximal** run
lengths — how many runs of length 1, 2, 3, … it saw. This is essentially free
(the run tracker already knows every run boundary) and is the interesting
by-product: it shows how fast runs of each length thin out. Typical base-10
output:

```
  length        count      share    prev/this   smallest run starts at
    1          417442  95.3530%          --   2
    2           19849   4.5340%      21.031   523
    3             482   0.1101%      41.180   22193
    4              13   0.0030%      37.077   1442173
```

* `count` — number of maximal runs of that exact length
* `share` — percentage of all maximal runs
* `prev/this` — ratio to the previous length, i.e. the cost of one more term
* last column — the smallest prime at which a run of that length begins

Lengths of 64 or more are pooled into a single `>=64` bucket. The reported
`primes counted` is exactly `sum(length × count)` over the histogram, which
doubles as a consistency check.

Two caveats the program also prints for itself:

* Without `--all`, each thread stops at its first hit, so the histogram only
  covers the range scanned up to that point.
* With `--emirps`, the histogram is still over *digit-sum* runs. Testing every
  prime for emirp-ness would cost far more than the search itself.

---

## Requirements

| | |
|---|---|
| Compiler | any C99 compiler (gcc, clang) |
| RAM | ~1 MB of shared tables plus primesieve's per-thread segments; negligible |
| Library | **primesieve** ≥ 7.0 (headers + shared library) |
| Threads | POSIX threads |
| OS | Linux, macOS, or any POSIX system (uses `getline`, `truncate`, `fsync`, `sigaction`) |

primesieve does the heavy lifting: it generates primes in order, in the tens
to hundreds of millions per second, with a small memory footprint.

---

## Installing primesieve

**Debian / Ubuntu**

```sh
sudo apt install libprimesieve-dev
```

**Fedora / RHEL**

```sh
sudo dnf install primesieve-devel
```

**Arch**

```sh
sudo pacman -S primesieve
```

**macOS (Homebrew)**

```sh
brew install primesieve
```

**From source** (any platform, if the packaged version is old)

```sh
git clone https://github.com/kimwalisch/primesieve
cd primesieve
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
sudo cmake --install build
```

---

## Compiling

### Linux

```sh
cc -O3 -march=native -pthread -o prime_runs_ps prime_runs_ps.c -lprimesieve -lm
```

If you installed primesieve from source into `/usr/local` and the loader can't
find it at run time:

```sh
sudo ldconfig
```

### macOS

Homebrew does not install into the default search paths, so point the compiler
at its prefix (`/opt/homebrew` on Apple Silicon, `/usr/local` on Intel):

```sh
cc -O3 -pthread \
   -I"$(brew --prefix)/include" -L"$(brew --prefix)/lib" \
   -o prime_runs_ps prime_runs_ps.c -lprimesieve -lm
```

Notes for macOS:

* Apple clang on Apple Silicon may reject `-march=native`; use `-mcpu=native`
  or simply drop the flag — `-O3` alone gets nearly all of the benefit here.
* `-lm` is harmless but not strictly required.
* Linking against a Homebrew-only primesieve produces a binary that needs
  Homebrew present at run time. Use `-Wl,-rpath,"$(brew --prefix)/lib"` if you
  intend to move the binary around.

### Quick check

```sh
./prime_runs_ps -n 3 -E 6e4 --all
```

should report the runs starting at 22193 and 25373.

---

## Command-line options

Everything is optional except `--length`.

### Search definition

| Option | Default | Meaning |
|---|---|---|
| `-n`, `--length N` | *required* | Run length: find `N` consecutive primes sharing a digit sum. `N = 1` matches every prime, which is mostly useful for exercising the statistics. |
| `-S`, `--start V` | `0` | Lower bound of the search range, inclusive. |
| `-E`, `--end V` | `1e9` | Upper bound, **exclusive**. `-L` / `--limit` is an alias kept for compatibility. |
| `-b`, `--base B` | `10` | Base in which digit sums (and, with `--emirps`, digit reversals) are computed. Valid range 2..36. Base 10 uses a five-digit lookup table and a hard-coded mod 9; other bases take a slightly slower generic path. |

**Number format for `-S` and `-E`.** Values may be written as:

| Form | Example | Value |
|---|---|---|
| plain integer | `1000000000000` | exact, parsed as a 64-bit integer |
| scientific | `3e5`, `1e10`, `1.5e12` | 300000, 10000000000, 1500000000000 |
| decimal | `2.5e9` | 2500000000 |
| suffixed | `12K`, `0.5M`, `2.5B`, `1.5T` | ×10³, ×10⁶, ×10⁹, ×10¹² |
| separated | `1_000_000`, `1,000,000` | 1000000 |

Suffixes are case-insensitive and may be combined with decimals. Plain
integers never pass through floating point, so full 64-bit precision is
preserved; scientific and decimal forms carry double precision (53 bits),
which is exact for every power of ten up to 10¹⁹. Values above 2⁶⁴ clamp.

### Search behaviour

| Option | Default | Meaning |
|---|---|---|
| `--all` | off | Report every run found rather than stopping at the first. Without it, each worker thread stops as soon as it finds a hit in its own sub-range, and only the earliest is printed. |
| `--emirps` | off | Require every prime in the run to be an emirp — its base-`B` digit reversal must be a different prime. Palindromic primes do not qualify. Reversal primality is tested with a deterministic Miller–Rabin over the 12 standard bases, valid for all n < 2⁶⁴, and results are cached per prime so overlapping windows never retest. |
| `-t`, `--threads T` | `1` | Number of worker threads. The range is split into `T` contiguous chunks; each thread scans its chunk plus a 1 MiB (2²⁰ integers) overlap on each side so that runs straddling a boundary are seen whole. Ownership by first member means results and statistics are **identical for any `T`**. If the range is too small to split (fewer than 8 overlaps per chunk), the thread count is reduced automatically. |

### Checkpointing and interruption

| Option | Default | Meaning |
|---|---|---|
| `-c`, `--checkpoint F` | off | Enable crash-safe checkpointing to file `F`. Hits are streamed to `F.runs` as they are found. Rerunning the identical command line resumes from where it stopped. |
| `--ckpt-secs S` | `30` | Seconds between checkpoint writes. Minimum 1. |
| `--restart` | off | Delete an existing checkpoint and its results file, and start over. |

**How resume works.** Each worker records the *head of its in-progress run*,
not its current position, so on resume the rescan rebuilds both the run state
and its statistics contribution exactly. A per-thread high-water mark of
already-emitted hits prevents duplicates in the rescanned region. The
checkpoint records the byte offset of `F.runs` at snapshot time and resume
truncates the file to it, so the two files can never disagree no matter when
the process died. The checkpoint itself is written to a temp file and renamed,
so it is never observed half-written.

The thread partition is stored in the checkpoint and restored from it, so `-t`
may differ between sessions — it will simply be ignored on resume. Every other
parameter (`-n`, `-b`, `-S`, `-E`, `--all`, `--emirps`) must match; if it does
not, the program explains the mismatch and exits rather than producing
nonsense. Use a different checkpoint name or `--restart`.

**Interruption.** `SIGINT` (Ctrl-C) and `SIGTERM` stop the workers cleanly at
the next check, flush a final checkpoint, close the results file, print
whatever has been found so far along with the statistics, and exit with
status **130**. The output is explicitly marked as partial. Without a
checkpoint this still gives you a clean summary of the work done — it just
can't be resumed.

### Output

| Option | Default | Meaning |
|---|---|---|
| `--no-stats` | off | Suppress the run-length histogram. |
| `--sieve-size K` | auto | Override primesieve's segment size, in KB. primesieve normally picks this from detected cache sizes; if that detection is wrong (it can be on some ARM parts) a manual value near half the per-core L2 is often worth 20–50%. Sweep 64/128/256/512/1024 on a fixed sub-range and keep the winner. |
| `--stream-only` | off | Generate primes but skip the search entirely. Timing baseline: the difference against a normal run is the true cost of the search layer, which tells you whether to optimize the sieve or the tracker. |
| `--quiet` | off | Suppress the progress line on stderr. Results still go to stdout, and checkpointing is unaffected. |
| `-h`, `--help` | | Print a usage summary. |

The progress line (stderr) shows percentage of the value range covered,
current throughput, elapsed time, ETA, and hits so far. Results go to stdout,
so `--quiet` plus a redirect gives a clean capture:

```sh
./prime_runs_ps -n 8 -E 1e12 -t 16 --all --quiet > hits.txt
```

### Exit status

| Code | Meaning |
|---|---|
| 0 | completed normally |
| 2 | bad arguments |
| 3 | checkpoint parameter mismatch |
| 4 | could not open the results file |
| 130 | interrupted by SIGINT/SIGTERM |

---

## Performance notes

Two costs are in play: generating the primes, and tracking runs over them.
`--stream-only` separates them. On a typical machine near 10¹² the split is
roughly 60/40 in favour of the sieve, so:

* **If the sieve dominates**, tune `--sieve-size` first (see the option
  table), then check that you aren't splitting the range so finely that each
  worker pays the fixed startup cost of building sieving primes up to
  `sqrt(end)`. That cost is a fraction of a second per thread near 10¹⁷ but
  grows with the range, and it is paid again on every resume — favour long
  runs between restarts and a thread count that matches physical cores.
* **If the tracker dominates**, the digit-sum window above is already the
  cheap path; the remaining per-prime work is a compare and a branch.

Larger `--sieve-size` is not monotonically better: throughput peaks when the
segment fits comfortably in L2 and falls off sharply beyond it.

Runs of length `L` require every gap in the run to be divisible by 18 — equal
digit sums force equality mod 9, and both primes are odd — so at large `p`,
where the mean gap approaches and exceeds 18, candidate regions are unusually
dense prime clusters. That structure is not currently exploited but is the
obvious lever for a future version.

## Examples

Smallest run of 3 consecutive primes with equal digit sum:

```sh
./prime_runs_ps -n 3 -E 1e5
```

Every run of 4 below 10¹⁰, on 16 threads, capturing to a file:

```sh
./prime_runs_ps -n 4 -E 1e10 -t 16 --all --quiet > runs4.txt
```

A long hunt with checkpointing every two minutes — safe to Ctrl-C and rerun:

```sh
./prime_runs_ps -n 9 -E 1e14 -t 16 --checkpoint run9.ck --ckpt-secs 120
# ...later, same command line resumes:
./prime_runs_ps -n 9 -E 1e14 -t 16 --checkpoint run9.ck --ckpt-secs 120
```

Emirp runs in base 10:

```sh
./prime_runs_ps -n 3 -E 1e9 --emirps --all
```

Digit sums in base 7, statistics only:

```sh
./prime_runs_ps -n 64 -b 7 -E 1e9 -t 8
```

Resume a specific segment of a large search (each segment can use its own
checkpoint file and run on a different machine):

```sh
./prime_runs_ps -n 8 -S 4e12 -E 5e12 -t 32 --all --checkpoint seg4.ck
```

---

## Files written

| File | Contents |
|---|---|
| `F` | checkpoint state: parameters, per-thread positions, histogram, results-file offset |
| `F.runs` | one line per hit: digit sum followed by the `L` primes, appended as found |
| `F.tmp` | transient; renamed over `F` on each checkpoint |

`F.runs` is plain text and safe to read while the job runs, though the tail
may be a partial line. On completion, the program re-reads it, sorts by first
member, and prints the full ordered list.
