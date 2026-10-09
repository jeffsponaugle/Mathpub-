# dschain and dspattern — digit-sum prime chains (OEIS A090009)

Two search programs for **digit-sum prime chains**: sequences of primes where
each term is the previous term plus the sum of its own decimal digits,

```
p  ->  p + digitsum(p)  ->  p + digitsum(p) + digitsum(p + digitsum(p))  ->  ...
```

with every term along the way prime. For example, starting at 516493:

```
516493 (+28) -> 516521 (+20) -> 516541 (+22) -> 516563 (+26) -> 516589 (+34) -> 516623
```

is a chain of length 6 — six consecutive primes under the map `p -> p + digitsum(p)`
(the step after 516623 has an odd digit sum, so it lands on an even number and the
chain ends).

- **`dschain`** scans every prime in a range, follows the chain from each, and
  keeps a histogram of chain lengths. It is the reference: simple, and it sees
  everything.
- **`dspattern`** tests only numbers whose *digits* allow a chain of the
  requested length at all — for 10-prime chains, one number in 10¹¹ below 2⁶⁴ —
  and works in 128 bits up to 10³⁰. It is exhaustive, and for long chains it is
  faster than `dschain` by a factor of about 10⁹. On a Mac, `--gpu` adds a
  Metal GPU search, about ten times the speed of the CPU alone.

## New results

- **a(10) = 63,782,998,969,989,799,429.** The first chain of ten primes:

  ```
  63782998969989799429 (+142) -> 63782998969989799571 (+140) -> 63782998969989799711 (+136)
  -> 63782998969989799847 (+146) -> 63782998969989799993 (+148) -> 63782998969989800141 (+116)
  -> 63782998969989800257 (+124) -> 63782998969989800381 (+122) -> 63782998969989800503 (+118)
  -> 63782998969989800621 (+119)
  ```

  Found by `dspattern` in 5.3 seconds on a laptop. Every term is proven prime
  (13-base Miller–Rabin and an independent Pratt certificate), the chain is a
  true start, and every smaller candidate above a(9) was tested — twice, by two
  independent programs. See [*Finding a(10)*](#finding-a10) for the evidence.
- **a(11) = 8,968,999,974,999,983,898,865,998,461.** The first chain of eleven
  primes:

  ```
  8968999974999983898865998461 (+206) -> 8968999974999983898865998667 (+214)
  -> 8968999974999983898865998881 (+212) -> 8968999974999983898865999093 (+208)
  -> 8968999974999983898865999301 (+200) -> 8968999974999983898865999501 (+202)
  -> 8968999974999983898865999703 (+206) -> 8968999974999983898865999909 (+214)
  -> 8968999974999983898866000123 (+176) -> 8968999974999983898866000299 (+190)
  -> 8968999974999983898866000489 (+191)
  ```

  Found in 53 minutes on a Mac Studio (M2 Ultra), 23 of them on its GPU. All
  eleven terms are proven prime by Pratt certificates, it is a true start, and
  an independent re-search with a different GPU kernel tested all 18.4 trillion
  candidates below it and found no other chain. See
  [*Finding a(11)*](#finding-a11).
- **No chain of 12 or more primes starts below 10³⁰.** Below 7.99989×10²⁴ no
  number's digits even allow one; below 10³⁰ there is a single digit pattern,
  and all 379 million numbers that fit it fail, each by a proof of
  compositeness: all but 711 have a term divisible by a prime below 500, and
  those 711 each have a term that fails the base-2 strong test.
- **A long chain is a digit pattern.** Every chain of 10 primes below 10²¹ ends
  in one of five four-digit patterns, …9743, …9429, …9329, …9419 or …9319,
  each with a prescribed digit sum and number of trailing 9s above them. See
  [*Why long chains are rare*](#why-long-chains-are-rare-digit-patterns).

## OEIS submissions

[`submission.md`](submission.md) lists every submission this project's data
supports, each checked against the live OEIS (last on October 3, 2026), with
draft text; the data files and their evidence are in [`oeis/`](oeis/README.md).

| sequence | submission | status (October 3, 2026) |
|---|---|---|
| A090009 | a(10) = 63782998969989799429, a(11) = 8968999974999983898865998461, and a comment | approved |
| A320879 | b-file of 4,379 terms < 10¹⁶, restoring 16 terms missing from the old one | approved |
| A320880 | b-file of 886 terms < 10¹⁸ (was 15) | approved |
| A320878 | b-file of 10,000 terms (was 7,626) | approved |
| new | chains of ≥ 10 primes: 665 terms < 10²² | ready |
| A048527, A048523, A048519, A320881 | b-files of 10,000 terms (were 3,000, 1,000, 1,000 and 63) | ready |
| A048524, A048525, A048526 | b-files of 20,000 terms (were 10,000) | optional |
| A090009, A320880 | cross-references to the new entry | after it is approved |
| Puzzle 163 (primepuzzles.net) | the earliest chains of 10 and 11 primes | ready, by email |

## The sequence: OEIS A090009

[A090009](https://oeis.org/A090009) records, for each n, the smallest prime that
starts a chain of n such primes (offset 1):

| n | a(n) — smallest prime starting a length-n chain |
|---|---|
| 1 | 2 |
| 2 | 11 |
| 3 | 11 |
| 4 | 277 |
| 5 | 37,783 |
| 6 | 516,493 |
| 7 | 286,330,897 |
| 8 | 286,330,897 |
| 9 | 56,676,324,799 |
| 10 | **63,782,998,969,989,799,429** (new here; submitted to the OEIS, September 2026) |
| 11 | **8,968,999,974,999,983,898,865,998,461** (new here; submitted to the OEIS, September 2026) |
| 12 | above 10³⁰, if it exists |

(11 begins the chain 11 → 13 → 17, so it is the earliest start for both length 2
and length 3; likewise 286330897 actually starts a length-8 chain, covering n=7
and n=8 at once.)

The sequence was contributed by Joseph L. Pe (2004); a(7)–a(8) are due to Donovan
Johnson and a(9) to Giovanni Resta (2013). As of September 2026 the OEIS entry
lists a(1)–a(9) with the keyword "more". Related OEIS entries: A062028 (the map
n + digitsum(n)), A048519/A048523–A048527, A320878–A320880.

## Why long chains are rare: digit patterns

A prime p > 5 is odd and prime to 3 and 5, and so must every term of a chain be.
Those three conditions already do almost all the work:

- **Parity.** p is odd, so p + digitsum(p) is odd only if digitsum(p) is even.
  Every step of a chain needs an even digit sum.
- **3.** digitsum(p) ≡ p (mod 3), so each step doubles p mod 3; if p is prime
  to 3, so is every later term. Only the first term needs checking.
- **5.** No term may end in 0 or 5.

The digit sums along a chain are not independent coin flips, though, because a
chain barely moves. Write p = H·10⁴ + x, with x the last four digits. Each step
adds at most a couple of hundred, so a chain of up to 40 terms stays within 10⁴
of its start: it carries into H + 1 at most once. Every term's digit sum is then
digitsum(H) plus the digit sum of its own last four digits — or digitsum(H + 1)
plus that, after the carry — and digitsum(H + 1) = digitsum(H) + 1 − 9t, where t
is the number of trailing 9s of H.

So whether the first L terms are all odd and prime to 3 and 5 depends only on
**(S, t, x)**: the digit sum of the high part, its trailing 9s, and the last four
digits. That is a small, finite set, and checking all of it is instant. The
triples that pass are the **patterns** for length L. Below 2⁶⁴:

| L | patterns | note |
|---|---|---|
| 5 | 9,265 | |
| 6 | 2,801 | |
| 7 | 937 | |
| 8 | 410 | |
| 9 | 167 | 6 of them always have a term divisible by a prime ≤ 61 |
| 10 | **5** | every one needs the carry through trailing 9s |
| 11 | **0** | the first appears at p ≈ 8.0×10²³ (digit sum 196) |

The five patterns for ten primes (these are all of them below 10²¹, too):

| last four digits | digit sum of H | trailing 9s of H | step sizes | smallest p |
|---|---|---|---|---|
| 9743 | 53 | 4 | +76 +80 +88 +86 +28 +20 +22 +26 +34 | 9,899,999,743 |
| 9429 | 118 | 1 | +142 +140 +136 +146 +148 +116 +124 +122 +118 | 299,999,999,999,899,429 |
| 9329 | 119 | 1 | the same | 399,999,999,999,899,329 |
| 9419 | 119 | 11 | +142 +140 +136 +146 +148 +26 +34 +32 +28 | 398,999,999,999,999,419 |
| 9319 | 120 | 11 | the same | 498,999,999,999,999,319 |

The large jump in each step list (+88 → +28, +148 → +116, +148 → +26) is the carry
running through the trailing 9s of H, which lowers the digit sum by a multiple
of 9 at exactly the right moment. Nothing else keeps ten digit sums even. This is
why `dschain`'s long runs to 2×10¹⁴ never saw a chain of 10: in that range only a
few thousand numbers have the right digits, and none of them works.

`dspattern --patterns` prints the patterns and candidate counts for any length
and range.

## The two programs

### dschain

`dschain` scans a range `[--start, --end]` of chain **starting points** and reports
every chain of at least `--length` terms.

- **Candidate generation**: primes in the range are enumerated with
  [primesieve](https://github.com/kimwalisch/primesieve) if available, otherwise
  a built-in segmented Eratosthenes sieve (about 2–3× slower).
- **Chain following**: every chain *continuation* is checked with a deterministic
  12-base Miller–Rabin test, valid for all 64-bit integers. Chains are therefore
  never truncated at `--end` — a chain starting just inside the range is followed
  as far as it truly goes.
- **Maximal chains only (default)**: a hit is reported only if no prime r exists
  with `r + digitsum(r) == p`, i.e. p is the true start of the chain, matching
  A090009's definition. Use `--all` to report mid-chain suffixes too.
- **Parallelism**: the range is split into chunks handed out in order to OpenMP
  worker threads.
- **Checkpointing**: with `-k FILE`, progress is saved periodically, on Ctrl-C,
  and at the end; re-running the identical command resumes where it stopped.
- **Statistics**: every run ends with a histogram of chain lengths over all
  primes examined, plus an extrapolation of how much more work each extra link
  costs (empirically ~20–25× per link).

Everything fits in 64-bit arithmetic, so the search is sound up to
`--end` ≈ 1.8×10¹⁹ (the program guards against overflow of chain terms as well).

### dspattern

`dspattern` classifies the patterns for `--length` at startup, then, for each
pattern, generates every number in range that fits it — every H with the right
digit sum (and trailing 9s), in increasing order, without scanning — and tests
it:

1. **Sieve.** All L terms at once against the primes 7..499: the term offsets
   are fixed by the pattern, so each prime is one bitmask lookup. The loop never
   forms p itself. p is linear in the digits, so p mod 7·11·13·17·19·23 is a
   constant plus the digits times per-pattern weights, and it is updated as the
   digits change (they usually change in the last two places). That is tested
   for every candidate; 29..61 and 67..499 are worked out from the digits only
   for the few percent that pass, in 32-bit arithmetic. About 1 candidate in
   300,000 survives for chains of 11.
2. **Screen.** A base-2 strong probable-prime test on each term, in 64-bit or
   two-word Montgomery arithmetic.
3. **Proof.** A start whose L terms all pass is proven with deterministic
   Miller–Rabin (12 bases below 2⁶⁴, 13 above, proven below 3.3×10²⁴) and
   followed to its true length exactly as `dschain` does, including the maximal
   check.

Every rejection is a proof of compositeness, so the search is exhaustive:
`--first` returns the confirmed smallest start in the range, and "no chain" means
there is none. Positions are 128-bit, up to 10³⁰ (less a margin for chain
extensions). Miller–Rabin with 13 bases is only a *proof* of primality below
3.3×10²⁴ (Sorenson and Webster, 2015); above it the search is still exhaustive,
since rejections are still proofs, but the terms of a hit are strong probable
primes, and the output says so. `verify_chain.py` then certifies them with Pratt
certificates, which work at any size. Starts below 10⁴ are checked by brute
force. The Montgomery arithmetic, primality tests and number parser come from
`gapsieve.c` in the sibling `gapchain` project, and the GPU search follows
`gapsieve`'s design; see [*On the GPU*](#on-the-gpu-macos).

## Building

The easiest route on both platforms is the included Makefile, which autodetects
OpenMP and primesieve for `dschain` and prints what it found:

```bash
make
```

```bash
make config    # show detection results without building
```

```bash
make test      # quick correctness check of both programs
```

`dspattern` has no dependencies — any 64-bit gcc or clang with pthreads and
`__int128` — and builds on its own with:

```bash
make dspattern
```

On a Mac that includes the Metal GPU search (`--gpu`): `dspattern_gpu.m`
includes `dspattern.c`, and the Makefile embeds `dspattern.metal` with `xxd`, so
the binary stands alone and compiles its kernel at startup (no Xcode Metal
toolchain needed). `make GPU=0` leaves it out. Anywhere, the CPU-only program
is just:

```bash
cc -O3 -pthread -o dspattern dspattern.c
```

For `dschain` neither dependency is strictly required, but you want both for
real searches: **OpenMP** for multithreading and **primesieve** for fast prime
generation.

### macOS

Apple's stock clang ships without OpenMP support (`-fopenmp` is rejected), so
install Homebrew's `libomp`; `primesieve` is also in Homebrew:

```bash
brew install libomp primesieve
```

```bash
make
```

Manual build without the Makefile:

```bash
clang -O3 -march=native -Xpreprocessor -fopenmp \
  -I"$(brew --prefix libomp)/include" -I"$(brew --prefix primesieve)/include" \
  -L"$(brew --prefix libomp)/lib" -L"$(brew --prefix primesieve)/lib" \
  -lomp -lprimesieve -o dschain dschain.c
```

### Linux

```bash
sudo apt install build-essential libprimesieve-dev   # Debian/Ubuntu
```

```bash
make
```

Manual build:

```bash
cc -O3 -march=native -fopenmp -o dschain dschain.c -lprimesieve
```

### Without primesieve

If primesieve is not installed, the Makefile falls back automatically, or build
by hand with the built-in segmented sieve:

```bash
cc -O3 -march=native -fopenmp -DNO_PRIMESIEVE -o dschain dschain.c -lm
```

The built-in sieve first tabulates all base primes up to √end (using about √end
bytes of memory while doing so — ~1 GB at `--end` 1e18) and refuses to go past
√end = 4×10⁹, i.e. `--end` ≈ 1.6e19; use primesieve for the extreme top of the
64-bit range.

If the binary was built without OpenMP, `--threads` is ignored and the program
prints a prominent warning with per-platform rebuild instructions.

## dschain usage

```bash
./dschain [options]
```

| Flag | Long form | Default | Meaning |
|---|---|---|---|
| `-l N` | `--length N` | 6 | Minimum chain length to report |
| `-s N` | `--start N` | 2 | Lowest starting prime to consider |
| `-e N` | `--end N` | 1e9 | Highest starting prime to consider |
| `-t N` | `--threads N` | all cores | Number of OpenMP worker threads |
| `-c N` | `--chunk N` | 1e7 | Range size per work unit (auto-adjusted to keep all threads busy and the chunk map small) |
| `-a` | `--all` | off | Also report chains that are suffixes of longer chains (default: maximal chains only) |
| `-F` | `--first` | off | Find the **earliest** chain of length ≥ `--length`: hits print as they are found, work above the best hit is dropped, and the sweep below it runs to completion so the smallest starting prime is confirmed |
| `-k F` | `--checkpoint F` | none | Save progress to file F; resume from it if it exists |
| `-i S` | `--interval S` | 60 | Seconds between periodic checkpoint writes |
| | `--restart` | | Ignore an existing checkpoint and start the range over |
| `-q` | `--quiet` | off | Suppress the live progress line; print results only |
| `-h` | `--help` | | Usage message |

### Number formats

`--start`, `--end`, `--chunk`, `--threads`, and `--interval` all accept:

- plain integers, with optional `_`, `,`, `'` or space separators: `1_000_000`, `1,000,000`
- scientific notation: `1e9`, `2.5e8`, `3.2e14`
- magnitude suffixes (case-insensitive): `900K`, `12M`, `4B` (= `4G`), `2.5T`, `3P`, `1E`
  — thousand, million, billion, trillion, quadrillion, quintillion. Fractions
  work (`1.5B` = 1,500,000,000). A trailing `E` means exa (`3E` = 3×10¹⁸), but
  `3E9` still parses as scientific notation. Maximum is the 64-bit limit, ~18.4E.

Anything unparseable is a hard error, not silently ignored.

### Checkpointing

```bash
./dschain -l 8 -s 1e10 -e 3.2e14 -k len8.ckpt
```

With `-k`, progress is written every `--interval` seconds, when the run finishes,
and on the first Ctrl-C / SIGTERM (a second Ctrl-C aborts immediately without
saving). Re-run the **same command** with the same `-k` file to continue:
`--start`/`--end`/`--length`/`--all` must match what's in the file, and the
chunk size stored in the checkpoint wins over `-c`. `--restart` throws the
checkpoint away and starts over.

The checkpoint is a small, human-readable text file (chunk watermark, a list of
out-of-order completed chunks, partial-chunk positions, and the accumulated
statistics). Chunks that were mid-flight when a checkpoint landed are redone on
resume, so a resumed run may re-print a hit it already reported — dedupe output
if you're collecting hits across restarts.

Chunks are handed to the threads in increasing order (`schedule(monotonic:
dynamic)`). Earlier builds used a plain dynamic schedule, which LLVM's libomp
runs as work stealing: each thread started with its own contiguous 1/T of the
range, so the watermark crept along with thread 0 while the checkpoint listed
millions of finished chunks above it — 55 MB at 46% of a 10¹⁴ range. The
results were never affected, only the file size and the progress display.
Those older checkpoints still resume correctly.

### Progress and interruption

When stdout is a terminal, a single self-updating status line shows percent
done, current position, primes examined, throughput, elapsed/ETA, the longest
chain seen, and the hit count. When piped or redirected, the output is plain
line-oriented text instead. Exit status: 0 for a completed run (or a `--first`
hit), 130 when stopped by a signal, 2 on usage errors.

## dspattern usage

```bash
./dspattern [options]
```

| Flag | Long form | Default | Meaning |
|---|---|---|---|
| `-l N` | `--length N` | 10 | Minimum chain length to find (4–40; `dschain` is the tool for shorter chains) |
| `-s N` | `--start N` | 2 | Lowest starting prime to consider |
| `-e N` | `--end N` | max | Highest starting prime; `max` is 10³⁰ less a margin (a larger `--end` is lowered to it, with a note). Above 3.3×10²⁴ hits are reported as probable primes, to certify with `verify_chain.py` |
| `-t N` | `--threads N` | all cores | CPU worker threads. With `--gpu` the default is all cores but two (the GPU is fed by one more thread), and `-t 0` runs the GPU alone |
| | `--gpu` | off | Also search on the GPU (macOS builds). The GPU takes work from the same queue as the CPU threads, so checkpoints, `--first` and the statistics are unchanged |
| `-a` | `--all` | off | Also report chains that are suffixes of longer chains (default: maximal chains only) |
| `-F` | `--first` | off | Report only the earliest chain. Work items run in order of their smallest start, work above a hit is dropped, and everything below it finishes, so the answer is the confirmed smallest start in the range |
| `-k F` | `--checkpoint F` | none | Save progress to file F; resume from it if it exists |
| `-i S` | `--interval S` | 60 | Seconds between checkpoint writes |
| | `--restart` | | Ignore an existing checkpoint and start over |
| `-p` | `--patterns` | | List the digit patterns for `--length`, with the candidates each has in range and the smallest, then exit |
| `-q` | `--quiet` | off | No progress output; results and summary only |
| `-h` | `--help` | | Usage message |

Numbers are converted **exactly**, in 128 bits: `1e20`, `2.5e21`, `300T`,
`18_446_744_073_709_551_616`, plain digits, and `max`. The suffixes are
dschain's (K, M, B/G, T, P, E).

Hits print in `dschain`'s format as they are found, then a summary: candidates
(numbers that fit a pattern), how many survived the sieve, how many base-2 tests
ran, and the chains found. On a terminal a status line shows percent done,
candidates per second, the position below which the search is complete,
elapsed time and ETA. With output redirected, the status line appears once a
minute instead.

Checkpoints work like `dschain`'s — the same command with the same `-k` file
resumes — and additionally record the `--first` hit, so a resumed `--first` run
never searches past its answer. `--length`, `--start`, `--end`, `--all` and
`--first` must match, as must the way the work was cut (checked with a hash, so
a checkpoint from a different build is refused rather than misread). Items in
flight at an interrupt are simply redone. Exit status as for `dschain`.

## Examples and expected results

Times are for an Apple M4 Max laptop (14 CPU cores).

### Quick sanity check (`make test`)

There are exactly **12 maximal chains of length ≥ 6 below 10⁸**:

```bash
./dschain -l 6 -s 2 -e 100M -q
```

```
[len 6] 516493(+28) -> 516521(+20) -> 516541(+22) -> 516563(+26) -> 516589(+34) -> 516623(+23)
[len 6] 1056493(+28) -> ...
[len 6] 1427383(+28) -> ...
[len 6] 1885943(+38) -> ...
[len 6] 3166183(+28) -> ...
[len 6] 3805183(+28) -> ...
[len 6] 4241593(+28) -> ...
[len 6] 6621283(+28) -> ...
[len 6] 7646953(+40) -> ...
[len 6] 12912283(+28) -> ...
[len 6] 17987839(+52) -> ...
[len 6] 32106493(+28) -> ...

done
  primes examined  : 5761455
  chains >= 6      : 12
  longest seen     : 6 (starting at 516493)
```

`./dspattern -l 6 -s 2 -e 100M -q` prints the same 12. (Hits arrive in whatever
order the worker threads finish, so the printed order varies run to run; the set
does not. Each term is printed with the digit sum it adds in parentheses. The
trailing `(+23)` on the last term shows why the chain stops: an odd digit sum
means the next number is even.)

The earliest of these is 516493 = A090009(6). Note the striking repetition:
nine of the twelve chains end in `...283/...493/...383/...593` and step through
the identical digit-sum pattern `+28, +20, +22, +26, +34` — these are digit
patterns at work, as described above.

### Reproducing the known sequence terms

```bash
./dschain -l 7 -s 2 -e 1B --first -q
```

```
[len 8] 286330897(+46) -> 286330943(+38) -> 286330981(+40) -> 286331021(+26) -> 286331047(+34) -> 286331081(+32) -> 286331113(+28) -> 286331141(+29)

earliest chain of length >= 7 starts at 286330897
(confirmed smallest: every prime below it in range was examined)
```

confirming a(7) = a(8) = 286330897 (the first length-7 chain happens to extend
to length 8). Runs in well under a second with primesieve. Similarly:

```bash
./dschain -l 9 -s 2 -e 60B --first -k len9.ckpt
```

finds a(9) = 56676324799 in a few minutes:

```
[len 9] 56676324799(+64) -> 56676324863(+56) -> 56676324919(+58) -> 56676324977(+62) -> 56676325039(+52) -> 56676325091(+50) -> 56676325141(+46) -> 56676325187(+56) -> 56676325243(+49)
```

`dspattern` finds it in a quarter of a second, testing 7 million candidates
instead of 2.4 billion primes:

```bash
./dspattern -l 9 -s 2 -e 1e11 --first
```

`--first` is race-safe in both programs despite the parallel work: when a hit
is printed, work *above* it is abandoned, but every in-flight region *below* it
runs to completion. If a thread sweeping a lower region finds an earlier chain,
that hit is printed too and becomes the new bound (hits therefore appear in
decreasing order of starting prime). The final confirmation line is only
printed once everything below the winner has actually been examined; a run
interrupted by Ctrl-C reports its best hit as *not* confirmed instead.

### Finding a(10)

`dschain` searched 10¹² to 2×10¹⁴ for a chain of ten primes (at about
2×10¹⁰ integers per second, ~4 hours) and found none; in that range only a few
thousand numbers have digits that allow one. The whole 64-bit range holds just
176,288,301 such numbers, and `dspattern` tests them all in half a second:

```bash
./dspattern -l 10 -e 18446744073709551615
```

```
no chain of length >= 10 starts in this range (exhaustive)
  candidates       : 176,288,301
  after sieve      : 173,444
  base-2 tests     : 210,385
  chains >= 10     : 0
```

So a(10) > 2⁶⁴. Above 2⁶⁴, the command that found it:

```bash
./dspattern -l 10 -s 56676324800 -e 1e21 --first -k a10.ckpt
```

```
dspattern: p -> p + digitsum(p), chains of length >= 10
  range     : 56,676,324,800 .. 1,000,000,000,000,000,000,000
  patterns  : 5 (5 of them carry past the low four digits)
  candidates: 2,009,904,995,948 numbers in range fit a pattern (1 in 4.98e+08)
  threads   : 14
  mode      : maximal chains only, earliest only
  checkpoint: a10.ckpt every 60s

[len 10] 63782998969989799429(+142) -> 63782998969989799571(+140) -> 63782998969989799711(+136) -> 63782998969989799847(+146) -> 63782998969989799993(+148) -> 63782998969989800141(+116) -> 63782998969989800257(+124) -> 63782998969989800381(+122) -> 63782998969989800503(+118) -> 63782998969989800621(+119)

earliest chain of length >= 10 starts at 63782998969989799429
(confirmed smallest: every candidate below it in range was tested)
  elapsed          : 5.32s
  candidates       : 3,581,358,838
  after sieve      : 3,909,550
  base-2 tests     : 4,703,523
  chains >= 10     : 1
  longest found    : 10 (starting at 63782998969989799429)
```

It starts just above a(9), since a chain of ten primes is also a chain of nine.
The log of the run is in `a10.log`.

**The evidence.**

- **Primality.** `verify_chain.py` (pure Python, sharing no code with either
  program) proves all ten terms prime twice: by 13-base Miller–Rabin, which is a
  proof below 3.3×10²⁴, and by a Pratt certificate for each term (n − 1 fully
  factored, a Lucas witness, and every factor proven the same way). It also
  confirms the digit sums, that the tenth digit sum (119) is odd so the chain
  ends there, and that no prime r has r + digitsum(r) = a(10).
- **Minimality.** A separate program, sharing no code with `dspattern`,
  re-derived the five patterns, enumerated candidates by recursive digit search
  instead of `dspattern`'s generator, and tested each with GMP's BPSW test,
  following chains by the real digit sums of the full numbers. Over the whole
  64-bit range it found the same 176,288,301 candidates, and none starts a chain
  longer than 8. Over (a(9), a(10)] it found the same 3,577,235,998 candidates
  as `dspattern` and exactly one chain of ten, a(10) itself.
- **The patterns.** For every digit sum and trailing-9 count possible below
  2⁶⁴, and every one of the 10⁴ possible last four digits, 92 million chain
  prefixes computed from the decomposition matched the same prefixes walked with
  full-number digit sums, with no disagreements. The classification gives the
  same answer whether the number is split at 10⁴, 10⁵ or 10⁶.
- **Above 2⁶⁴.** An independent Python brute force (every odd number, 13-base
  Miller–Rabin, real digit sums) matched `dspattern --all` exactly in windows at
  2⁶⁴, at 10²⁰, across 10²¹, across 3.3×10²⁴ (where primality switches from
  proven to probable) and just below the old 3.3×10²⁴ ceiling. Checkpoints were
  interrupted and resumed, including a `--first` checkpoint with a planted fake
  hit, which the real a(10) overturned, and resumed runs give the same counts
  as uninterrupted ones.
- **The heuristic.** A Hardy–Littlewood estimate built from each pattern's
  residues mod every prime matched the observed counts of partial chains among
  the 64-bit candidates (11,540 predicted vs 11,481 observed starts with four
  prime terms; 799 vs 757 with five; 75 vs 66 with six). It predicts 0.37 chains
  of ten below 10²⁰ and 20 below 10²¹, so a(10) at 6.4×10¹⁹ is right where it
  should be.

To check it yourself:

```bash
python3 verify_chain.py 63782998969989799429
```

### Listing patterns

```bash
./dspattern -l 10 --patterns
```

```
dspattern: p -> p + digitsum(p), chains of length >= 10
  range     : 2 .. 999,999,999,999,999,999,999,999,961,119
  patterns  : ... (... of them carry past the low four digits)
  candidates: ... numbers in range fit a pattern

  low digits 9743  ds(H)=53  t=4    ds(p)=76  steps +76 +80 +88 +86 +28 +20 +22 +26 +34
      ... candidates in range, the smallest 9899999743
  ...
```

`./dspattern -l 11 -e 3.3e24 --patterns` shows the single 11-prime pattern
below 3.3×10²⁴ (low digits 8911, fourteen trailing 9s above them), with its 203
candidates, none of which survives the sieve.

### Finding a(11)

No number below 8×10²³ has digits that allow 11 terms, and the heuristic put the
first chain far above the 3.3×10²⁴ limit of proven Miller–Rabin. It expected it
between about 5×10²⁷ and 3×10²⁸. The search ran on the Mac Studio, in
`~/A090009` there:

```bash
nohup caffeinate -i ./dspattern -l 11 -s 2 -e 5e28 --first -k a11.ckpt -i 300 --gpu > a11.log 2>&1 &
```

It started on the CPU alone (24 threads, 1.06×10⁹ candidates a second) and,
after half an hour at 5×10²⁷, was stopped and resumed from its checkpoint with
the first GPU kernel (about 1.2×10¹⁰ a second). 23 minutes later it reported
the chain above, at 8.97×10²⁷, and finished the ranges below it:

```
earliest chain of length >= 11 starts at 8968999974999983898865998461
(confirmed smallest: every candidate below it in range was tested)
(its terms reach 3.3e24, beyond proven Miller-Rabin: certify them with verify_chain.py)
  elapsed          : 3166.76s (including earlier runs)
  candidates       : 18,415,991,319,744
```

(18.4 trillion candidates, out of the 7×10¹⁴ in the whole range to 5×10²⁸;
`caffeinate -i` keeps the machine from idle-sleeping during such a run.) The
chain is the length-11 pattern with low digits …8461 below a high part whose
digit sum is 187 and which ends in two 9s: its eighth step, +214, carries
through …99909 into …00123, and the digit sum drops to 176.

**The evidence.**

- **Primality.** `python3 verify_chain.py 8968999974999983898865998461`
  proves all eleven terms prime with Pratt certificates (each n − 1 fully
  factored, a Lucas witness, every factor proven the same way) in 1.5 seconds.
  It also confirms the digit sums, that the eleventh (191) is odd, and that no
  prime r has r + digitsum(r) = a(11).
- **Minimality.** The laptop re-searched everything from 2 to a(11) with the
  second GPU kernel (`list_kernel`), which enumerates candidates differently
  from both the CPU code and the first kernel that found it, in `--all` mode:

  ```bash
  ./dspattern --gpu -t 0 -l 11 -s 2 -e 8968999974999983898865998461 -a
  ```

  It tested all 18,411,673,898,715 candidates in 402 seconds and found exactly
  one chain of 11 or more, a(11). The log is in `a11-recheck.log`, and the
  Studio's is `a11.log`.
- **The candidates, counted independently.** `oeis/indep_count.py`, written
  from scratch in Python with its own pattern classification and digit-sum
  counting, finds exactly the same 18,411,673,898,715 candidates below a(11)
  as `dspattern` (and 176,288,301 for chains of 10 below 2⁶⁴).
- **The patterns at this height.** For every high-part digit sum up to 216 and
  up to 24 trailing 9s — everything a 28-digit number can have — and every
  one of the 10⁴ low parts, 140.9 million chain prefixes from the classification
  matched the same prefixes walked on real 128-bit numbers, with no
  disagreements.
- **The sieve at this height.** On the slice 7×10²⁷..7.2×10²⁷ the GPU and the
  CPU produced the same set of sieve survivors, number for number.

To see the chain in a few seconds:

```bash
./dspattern --gpu -l 11 -s 8.968e27 -e 8.97e27
```

## Performance

Chains of 11, the decade 10²⁶..10²⁷ (68.8 billion candidates, none of which
starts a chain), with the first GPU kernel (`enum_kernel`, one GPU thread per
run of candidates):

| machine | CPU alone | GPU alone | GPU + CPU threads |
|---|---|---|---|
| M4 Max laptop (14 cores, 32-core GPU) | 80.0 s (8.6×10⁸/s) | 8.2 s (8.4×10⁹/s) | 7.5 s with 12 threads (9.2×10⁹/s) |
| M2 Ultra Mac Studio (24 cores, 60-core GPU) | 64.7 s (1.06×10⁹/s) | 7.5 s (9.1×10⁹/s) | 6.7 s with 22 threads (1.02×10¹⁰/s) |

The second kernel (`list_kernel`, the default now: a SIMD group shares the
high digits and its lanes split a list of low ones) on the laptop's GPU:

| range | `enum_kernel` | `list_kernel` |
|---|---|---|
| 10²⁶..10²⁷ (6.9×10¹⁰ candidates) | 8.2 s | 2.8 s (2.5×10¹⁰/s) |
| 7×10²⁷..7.8×10²⁷, in the a(11) range (1.7×10¹² candidates) | 184.4 s | 41.4 s (4.2×10¹⁰/s) |

For the a(10) search, `dschain` managed 2×10¹⁰ *integers* a second, which is
about 100 years to reach 6.4×10¹⁹. `dspattern` on the laptop's CPU took 5.3
seconds, and its GPU 0.5.

## On the GPU (macOS)

A GPU would have given `dschain` perhaps a 10–50× speedup, sieving as `gapsieve`
does with Metal and CUDA in the sibling `gapchain` project. The digit patterns
gave about 10⁹×, so a(10) needed no GPU at all. a(11) is a search of hours to
days on the CPU, though, and there `--gpu` gives about ten times more.

**How it works.** A work item is an aligned block of Q — a fixed prefix and r
free digits with a given digit sum — so the number of candidates in it is known
exactly from the counting table. The GPU takes a unit of items (about 2³¹
candidates) at a time. There are two kernels, which find exactly the same
survivors.

`enum_kernel`, the first, gives each GPU thread K = 1024 consecutive
candidates of one item (`DSPATTERN_GPU_KERNEL=1` selects it):

1. **Unranking.** The thread finds the digits of its first candidate directly
   from the counting table, choosing each digit from the top by skipping the
   blocks of candidates with smaller digits.
2. **Enumeration and sieve.** It then runs `run_item`'s loop: the next digit
   string with the same digit sum, p mod 7·11·13·17·19·23 updated from the
   pattern's digit weights, and the stage-1 test; the later stages are worked
   out from the digits for the few that pass. Everything is 32-bit, which is
   what Apple GPUs do natively. The digits live in four 32-bit registers, as
   nibbles, rather than in an array, which would spill to memory. The step to
   the next digit string needs no loop: with z the lowest nonzero digit and i
   the lowest digit above z that is not 9, every digit in between is 9, so
   the digit sum below i is d_z + 9(i − 1 − z) and the weighted sum comes from
   a prefix table. Both positions come from bit tricks on the nibbles.
3. **Survivors.** The rare starts that pass all three stages (about 1 in
   300,000 for chains of 11) go back as (item, rank) pairs. The CPU rebuilds
   each by unranking, and screens and proves it exactly as `run_item` does.

Two units are in flight, so the GPU sieves one while the CPU finishes the
other. Every command buffer's status is checked: a failed unit is redone, and
after three failures its items are searched on the CPU; a survivor list that
overflows is enlarged and the unit redone.

Getting the digits out of a per-thread array and into registers took the
laptop from 12.4 s to 8.8 s on the benchmark decade, and the loop-free
successor to 8.2 s. Timing stripped-down variants of the kernel showed where
the rest went: about 70% was the general successor step. The simple step
(units digit down one, tens digit up one) covers most candidates, but a
32-lane SIMD group can take it only when all 32 lanes can, and in the a(11)
range, where numbers are full of 9s, that is rare.

`list_kernel`, the default, keeps the lanes in lockstep instead. Write the
free digits as y = U·10⁷ + V. For a given U, the candidates are exactly the
7-digit V whose digit sum is need − ds(U) (and, when the chain carries, which
do not end in 9): one of the precomputed sorted lists, the same list for every
U with that digit sum. So each SIMD group takes K = 16,384 consecutive
candidates, shares one U at a time, and its 32 lanes stride through U's list
together:

1. **Starting point.** The group unranks its first candidate as above, then
   finds V's place in its list from the same counting tables.
2. **Stage 1 without division.** For each trailing-9 count t (which fixes
   the digit weights), a table holds each listed V's share of p mod 7, 11, 13,
   17, 19 and 23, packed into 26 bits. Once per U the six prime masks are
   rotated by U's share, so a candidate costs one coalesced load and six shifts.
3. **The next U.** When the list runs out, the group moves to the next U whose
   digit sum leaves a nonempty list: usually U + 1, otherwise the lowest digit
   that can rise is raised and the digits below it rebuilt, with U's residues
   moved by the difference. Every lane does the same, so the only divergence
   left is the rare candidate that passes stage 1.

The length of V matters because the numbers are full of 9s. With 5-digit V
the lists in the a(11) range averaged 17 to 45 entries, which left lanes idle
and made the group change U every iteration, so the kernel was no faster
than `enum_kernel`. 7 digits gives lists of 85 to 277 entries on average, and
tables of 40 MB per trailing-9 count. (`DSPATTERN_GPU_VS` sets it from 4 to 7;
8 would need a 32-bit shift and 400 MB tables.) The later stages now run to
499 rather than 199, which spares the CPU four in five of the starts it would
otherwise have to test; before that the feeder thread was nearly saturated.

**What is left.** On the a(11)-range slice, stripped-down variants put about
40% of the time in the stage-2 and stage-3 code: when any lane of 32 passes
stage 1 the whole group waits for it. Stage-2 residues tabulated like stage
1's, and stage-3 residues for U kept in threadgroup memory, would cut most of
that.

**Tuning and test hooks** (environment variables; none changes the results):
`DSPATTERN_GPU_KERNEL` (2, `list_kernel`, or 1, `enum_kernel`),
`DSPATTERN_GPU_VS` (digits of V, 7), `DSPATTERN_GPU_K` (candidates per SIMD
group, 16,384, or per thread for `enum_kernel`, 1024), `DSPATTERN_GPU_UNIT`
(candidates per unit, in millions, 2048), `DSPATTERN_GPU_DEPTH` (units in
flight, 2), `DSPATTERN_GPU_FAIL=N` (treat every Nth command buffer as failed),
`DSPATTERN_GPU_CAP=N` (start with N-entry survivor lists), and
`DSPATTERN_SURVIVORS=file` (write every start that passes the sieve to file,
from the CPU or the GPU, for comparing the two).

**Validation.** On the laptop's GPU (M4 Max), for each kernel:

- **Survivor sets.** The GPU alone produced exactly the CPU's set of sieve
  survivors, every start that passes the sieve, compared as sorted lists.
  `enum_kernel` (sieving to 199 at the time): chains of 6 at 10⁹, 8 at 10¹², 9
  at 10¹⁴ (781,740 survivors), 10 over all of 64 bits, and 11 over
  10²⁶..10²⁷. `list_kernel` (sieving to 499): the same five, plus chains of 10
  at 10²⁰..3×10²⁰ (2,045,464 survivors, numbers with few nonzero digits
  rather than many 9s) and 11 at 7×10²⁷..7.2×10²⁷ in the a(11) range.
- **Hits.** The same chains as `dschain` for lengths 4 to 9, including the top
  of the 64-bit range, and a(9) and a(10) under `--first`.
- **Failure paths.** Identical results with every third command buffer failed
  (redone), with every one failed (the CPU took over), with 100-entry survivor
  lists (enlarged and redone), with 7, 33, 100,000 and 1,000,000 candidates
  per GPU thread or SIMD group, 1M-candidate units, one unit in flight, CPU
  threads beside the GPU, and `list_kernel` with V of 4, 5, 6 and 7 digits.
- **Checkpoints.** A GPU run interrupted and resumed gave the same totals as an
  uninterrupted run. (The a(11) search moved from the CPU to the GPU the same
  way.)

On the Studio's GPU (M2 Ultra): with `enum_kernel`, before the a(11) search
moved to it, the same survivor set as its CPU for chains of 9 at 10¹⁴, the same
counts as the CPU on the benchmark decade, and a(10); with `list_kernel`, the
same survivor sets as its CPU for chains of 9 at 10¹⁴ (201,998 survivors) and
for 11 at 7×10²⁷..7.1×10²⁷ (38,221), and a(10).
