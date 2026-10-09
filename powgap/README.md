# powgap

Search for chains of primes whose gaps are successive powers of two:

```
p,  p+2,  p+6,  p+14,  p+30,  p+62,  p+126,  p+254, ...
```

The gap sequence is 2, 4, 8, 16, 32, ..., so the k-th member sits at offset
`2^(k+1) - 2` from the starting prime `p`. Example of a length-10 chain:

```
len 10  p = 100,217,229,301,697
        100217229301697 100217229301699 100217229301703 100217229301711
        100217229301727 100217229301759 100217229301823 100217229301951
        100217229302207 100217229302719
```

## Building

```bash
gcc -O3 -march=native -o powgap powgap.c -lpthread -lm
```

```bash
clang -O3 -mcpu=native -o powgap powgap.c -lpthread -lm   # Apple silicon
```

No dependencies beyond pthreads and libm. Requires a compiler with
`__uint128_t` (gcc or clang on any 64-bit target).

## How it works

A naive search would test every odd number with a primality test. powgap
instead layers three filters, each much cheaper than the next:

**1. Wheel of admissible residues.** For a small prime `q`, the chain
survives mod `q` only if no member is divisible by `q`, i.e. `p` avoids the
residues `-offset mod q` for every offset in the pattern. The offsets
`2^(k+1) - 2 = 2(2^k - 1)` take at most `ord_q(2)` distinct values mod `q`
(the multiplicative order of 2), so even long patterns forbid few residues
and the pattern is admissible for every length — but the allowed fraction is
tiny: 1 of 2 residues mod 2, 1 of 3 mod 3, 1 of 5 mod 5 (once the length is
at least 4), and so on. powgap multiplies small primes together into a wheel
modulus `W` (by default at most 8192 surviving residues; see `--wheel`) and only ever
looks at numbers congruent to a surviving residue mod `W`. For length 6 the
wheel alone discards over 99.9% of all integers.

**2. Constellation sieve.** The survivors are indexed as `c = base + r + W·j`
for each wheel residue `r`. For every sieving prime `q` up to `--bound`
(beyond the wheel primes), the `j` for which *some* chain member is divisible
by `q` form a few arithmetic progressions mod `q`; those bits are marked in a
per-residue-class bitmap. Because the offsets collapse to at most `ord_q(2)`
distinct classes, each prime needs only a handful of passes.

**3. Verification.** Whatever survives the sieve is checked member by member
with a deterministic 7-base Miller–Rabin test, valid for the entire 64-bit
range. The chain is extended as far as it goes (not just to `--len`), so the
report shows the true full length.

The range is cut into blocks of `W × JLEN` integers handed out to worker
threads via an atomic counter. Output is reassembled in ascending order
regardless of which thread finishes first. Numbers below the sieve bound are
handled by brute force, since a chain member there could equal a sieving
prime.

The startup header prints a first-Hardy–Littlewood-style estimate (singular
series × ∫ dx/ln^L x) of how many chains to expect in the range; it is
usually within a factor of ~2.

## Usage

```
powgap [options] <start> <end>
```

`start` and `end` accept `K M B/G T P E` suffixes, decimals, scientific
notation, and digit separators: `250B`, `1.5T`, `3e15`, `1_000_000`.

```bash
# length >= 6 chains up to 100 billion
./powgap -l 6 0 100B
```

```bash
# long-running search with checkpointing; ctrl-c stops cleanly, rerunning resumes
./powgap -l 9 -C search.ckpt 1T 100T
```

```bash
# chains of >= 5 CONSECUTIVE primes (no other prime between the members)
./powgap -l 5 -c 0 1e12
```

### Relation to OEIS A090807

With `-c`, powgap searches exactly the objects catalogued in OEIS
[A090807](https://oeis.org/A090807), *"First prime in the earliest chain of n
consecutive primes with gaps 2^1, 2^2, ..., 2^(n-1)"* (the same sequence as
[A079014](https://oeis.org/A079014) under different indexing). Every
published term has been re-confirmed by exhaustive search with this program —
for each n, a sweep from 1 with `powgap -l n -c` reports a(n) as the first
chain found:

| n | a(n) | verification command |
|---|---|---|
| 2 | 3 | `./powgap -l 2 -c 1 100` |
| 3 | 5 | `./powgap -l 3 -c 1 100` |
| 4 | 1,997 | `./powgap -l 4 -c 1 3000` |
| 5 | 2,237 | `./powgap -l 5 -c 1 10000` |
| 6 | 6,824,897 | `./powgap -l 6 -c 1 7e6` |
| 7 | 1,356,705,137 | `./powgap -l 7 -c 1 1.4e9` |
| 8 | 3,637,803,390,827 | `./powgap -l 8 -c 1 3.65e12` (~25 s on an M4 Max) |
| 9 | 14,014,732,040,120,297 | `./powgap -l 9 -c -b 4096 3.6e12 1.5e16` (~32 min) |

**The a(9) term was found with this program and is now published in the
OEIS** ("a(9) from Jeff Sponaugle, Aug 19 2026"; the same value also extends
[A079014](https://oeis.org/A079014)): an exhaustive sweep of
[3.6 × 10¹², 1.40148 × 10¹⁶] contains exactly one chain of 9 consecutive
primes, starting at p = 14,014,732,040,120,297. Starting the sweep at a(8)
is rigorous because a 9-chain's first 8 members form an 8-consecutive chain,
and the from-1 sweep above proves none exists below a(8). The chain:

```
14014732040120297 +2 14014732040120299 +4 14014732040120303 +8 14014732040120311
+16 14014732040120327 +32 14014732040120359 +64 14014732040120423
+128 14014732040120551 +256 14014732040120807
```

The would-be 10th member, 14,014,732,040,121,319 (gap 512), is divisible
by 13. The sequence remains marked `hard, more`, with a(10) as the open
term; sweeps with this program have so far established
a(10) > 4.1 × 10¹⁸ (October 2026). Extending the estimate methodology below to length 10 (the new
512-gap adds ~510/ln x to the exponential cost) puts the median for a(10)
near 10²², far beyond this program's 2⁶⁴ ceiling; the chance it lies within
reach (< 1.8 × 10¹⁹) is only about 1 in 2,500. Hunting a(10) needs 128-bit
arithmetic.

### Where to expect a(9) — prediction and outcome

This estimate was made *before* a(9) was found, calibrated from measured
runs of this program (Apple M4 Max) rather than pure theory. It predicted a
median first occurrence of 2 × 10¹⁶ with an 80% window of
[4 × 10¹⁵, 7 × 10¹⁶]; a(9) then turned up at 1.4 × 10¹⁶ — near the middle
of the window. The methodology, kept for the a(10)-and-beyond record:

**Density of raw length-9 chains** (consecutiveness ignored) follows the
Hardy–Littlewood form `A₉ / ln⁹ x`. Calibrating `A₉` from two windows four
orders of magnitude apart gives the same constant within 7%:

| window | length ≥ 9 chains | implied A₉ |
|---|---|---|
| `powgap -l 9 1e12 2e12` | 185 | 1943 |
| `powgap -l 9 1e16 1.005e16` | 726 | 1817 |

**The consecutiveness penalty** — the probability that no *other* prime
interrupts the chain — is what actually drives this sequence. Modelling each
gap of size g as costing a factor `exp(−α(g−2)/ln x)` and measuring the
fraction of chains that are fully consecutive:

| chain length | top gap | consecutive fraction | implied α per gap |
|---|---|---|---|
| 6 | 32 | 22.2% (1804 / 8129 at 10¹²) | — |
| 7 | 64 | 2.38% (416 / 17464 at 10¹²) | 1.014 |
| 8 | 128 | 0.049% (16 / 32344 at ~10¹³) | 0.993 |

α is flat at ≈ 1.0 for the two largest measurable gaps, so extrapolating
α ≈ 1.0 for the new 256-gap that length 9 introduces is well-grounded. As a
cross-check, the model retrodicts a(8): it predicts ~10–16 consecutive-8
chains in [4 × 10¹², 4 × 10¹³]; the measured count is 16.

Integrating the resulting density `A₉/ln⁹x · exp(−C/ln x)` (C ≈ 480–510)
until the expected count crosses the Poisson thresholds:

- 10% chance a(9) lies below ~4 × 10¹⁵
- **median first occurrence ≈ 2 × 10¹⁶**
- 90% chance found by ~7 × 10¹⁶ (a conservative +10% on the 256-gap α
  stretches this to ~1.2 × 10¹⁷)

That is a jump of ~6,000× over a(8), continuing the sequence's accelerating
growth (a(7) → a(8) was ~2,700×).

### Tuning the search for a(9)

The one optimization that matters is **lowering the sieve bound**. The
default auto-bound for `--len 9` (2¹⁷) is tuned for reporting chains, not for
this regime: at length 9 the wheel alone leaves only 1 in ~3,400 integers, so
deep sieving saves fewer Miller–Rabin tests than it costs in sieve passes.
Measured on a 2 × 10¹³ slice at 2 × 10¹⁶ (all bounds produce byte-identical
results — a useful correctness check):

| `-b` | candidates tested | time | throughput |
|---|---|---|---|
| 131072 (default) | 46,152 | 13 s | 1.5 × 10¹²/s |
| 16384 | 261,768 | 4 s | 5 × 10¹²/s |
| **2048–8192** | ~0.5–2.3 M | **3 s** | **~6.7 × 10¹²/s** |

Also measured: keep the default `-m 4M` (16 MB bitmaps ran *slower* — cache
effects), and start at 3.6 × 10¹² since a(9) ≥ a(8) by the prefix property.
With the tuned bound the hunt is tractable on a single desktop — about
**1 hour to the median estimate, 3 hours to 90% confidence, 8 hours to sweep
clear out to 2 × 10¹⁷**:

```bash
./powgap -l 9 -c -b 4096 -C a9.ckpt 3.6e12 1e17
```

The first `consec  9` line printed is a(9): output is emitted in ascending
order, so the first hit is the earliest, and checkpointing makes the
multi-hour run interruption-safe.

Two further options exist but are not needed at this scale: raising the
wheel's residue cap (8192 in `build_wheel`) would let 23 and 29 join the
wheel — for this pattern 23 admits 14 of 23 residues and 29 admits 20 of 29,
roughly another 2.4× — and the range splits trivially across machines.

### Options

| option | meaning |
|---|---|
| `-l, --len N` | Minimum chain length to report (default 6, max 24). Chains longer than N are still reported at their full length. |
| `-t, --threads N` | Worker threads (default: all online CPUs). |
| `-c, --consecutive` | Require the members to be *consecutive* primes: a chain only counts if no other prime falls between its members. Each report shows both the full length and the consecutive-prefix length. |
| `-b, --bound N` | Small-prime sieve bound (default: auto, from 2¹⁷ to 2²¹ depending on `--len`; smaller lengths leave more survivors so they are sieved harder). Accepts suffixes. |
| `-m, --mem N` | Sieve bitmap bytes per thread (default 4M). Larger blocks amortize per-block overhead over long ranges. |
| `-w, --wheel N` | Maximum wheel residues (default 8192). Larger values fold more small primes into the wheel, shrinking the bitmap every sieve prime must process, at the cost of more memory per thread (see *Memory vs speed* below). |
| `-C, --checkpoint F` | Save progress to file `F` every 10 s; if `F` exists and matches the parameters, resume from where it left off. |
| `--brute` | Naive mode: primality-test everything, no wheel or sieve. For testing/verification only. |
| `-q, --quiet` | Suppress the header and progress lines (final statistics still print). |
| `-h, --help` | Usage summary. |

Search results go to **stdout**; the header, progress line, and statistics go
to **stderr**, so `./powgap ... > chains.txt` keeps them separate.

### Progress, interruption, and statistics

While running, a status line on stderr updates once per second:

```
# 28/56 blocks  50.0%  14,748 chains  2s elapsed, ~2s left
```

Pressing **ctrl-c** (or sending SIGTERM) stops cleanly: workers finish their
current block, a checkpoint is saved if `-C` was given, and statistics print.
A second ctrl-c aborts immediately. On completion or interruption you get a
histogram of every chain length found:

```
# finished -- 42,760 chains found, 318,450 candidates confirmed, 12s elapsed
# chains by length:
#   len  6 : 38,753
#   len  7 : 3,691
#   len  8 : 295
#   len  9 : 21
```

With `-c` a second histogram of consecutive-prefix lengths is printed. Exit
codes: 0 = range completed, 130 = interrupted (resumable), 1 = error.

### Checkpointing details

The checkpoint is a small human-readable text file recording the search
parameters, the first number not yet searched, cumulative statistics, and
elapsed time. It is written atomically (temp file + rename) every 10 seconds
and deleted when the range completes. Statistics are counted only for blocks
whose output has actually been emitted, so a resumed run's combined output is
byte-identical to an uninterrupted run — no duplicates, no gaps. You may
resume with a different `-t`, `-b`, or `-m`; changing the range, `--len`, or
`-c` is rejected (delete the file to start over). A checkpoint written by an
interrupted `--brute` run resumes the same way.

## powcheck: verifying a result

`powcheck` is a small standalone companion tool that verifies a single chain
by hand: given a starting number, it walks the pattern member by member,
printing each value, its primality, the gap it should satisfy, and any
intervening primes that break consecutiveness, until the chain ends.

```bash
gcc -O3 -march=native -o powcheck powcheck.c        # or clang -mcpu=native
```

```bash
./powcheck 1997
```

```
p = 1,997
   1:                  1,997                   prime
   2:                  1,999   gap   2 = 2^1   prime
   3:                  2,003   gap   4 = 2^2   prime
   4:                  2,011   gap   8 = 2^3   prime
   5:                  2,027   gap  16 = 2^4   prime   [intervening prime: 2,017]
   6:                  2,059   gap  32 = 2^5   composite (29 x 71) -- chain ends
  -> chain length 5 (consecutive-prime prefix 4)  [powgap: len 5 (consec 4)]
```

The trailing `[powgap: ...]` line is in the same format powgap reports, so
the two are directly comparable. Composite members are shown with a small
factor when one exists (trial division to 10⁶). Numbers accept commas,
underscores, and powgap-style suffixes, so results can be pasted verbatim
from powgap output:

```bash
./powcheck 3,637,803,390,827        # a(8) -- prints all 8 members, then the
                                    # 256-gap failure that stops a(9) here
```

Multiple numbers can be checked in one invocation; the exit code is nonzero
only if an argument fails to parse. Like powgap, it uses the deterministic
7-base Miller–Rabin test, so verification is exact for all n < 2⁶⁴.

## Performance

Measured on an Apple M4 Max (14 cores, default settings). Throughput is the
width of the searched interval divided by wall time; it is nearly independent
of where the interval sits, because the sieve cost per block is constant and
only the (rare) surviving candidates pay for primality tests.

| search | interval width | time | throughput |
|---|---|---|---|
| `-l 6` near 10¹² | 100 × 10⁹ | 4 s | ~25 × 10⁹/s |
| `-l 6` near 10¹⁵ | 10¹² | 12 s | ~83 × 10⁹/s |
| `-l 6` near 10¹⁸ | 10¹² | 12 s | ~83 × 10⁹/s |
| `-l 8` near 10¹² | 10¹² | 7 s | ~140 × 10⁹/s |
| `-l 10` near 10¹⁴ | 10¹² | 1 s | ~1 × 10¹²/s |

Longer patterns are *faster* per integer: each extra member forbids more
residues, so the wheel and sieve leave fewer survivors. Wider intervals are
also faster per integer, since block sizes auto-scale and startup costs
amortize. As rough guidance on this class of machine: a `-l 6` sweep of a
10¹⁵-wide interval takes on the order of 3–4 hours; `-l 8` about half that;
`-l 10` well under an hour.

### Memory vs speed

The sieve's cost is bit-marking, already near the hardware limit (~2 cycles
per mark), and the survivors' Miller–Rabin tests are only a few percent of
the run. The one lever that trades memory for speed is the wheel: each
prime folded in shrinks the bitmap (and thus every sieve prime's work) by
the fraction of residues it admits, while multiplying the residue count —
for the length-10 pattern, 23 admits 13/23 (13× residues), 29 admits 19/29
(19×), 31 admits 26/31 (26×). `--wheel` raises the residue cap; pair it with
`--mem` so rows do not collapse to the 1024-bit minimum:

| setting | wheel | residues | bitmap/thread | measured gain (`-l 10`) |
|---|---|---|---|---|
| default | primes ≤ 19 | 972 | 4 MB | — |
| `-w 16K -m 16M` | ≤ 23 | 12,636 | 16 MB | ~15–20% |
| `-w 300K -m 64M` | ≤ 29 | 240,084 | 64 MB + 288 MB shared | none measured |

Gains fall well short of the 1.8×/2.7× that the density reduction alone
predicts, because per-(residue, prime) overhead grows with the residue
count. All settings produce identical output, and checkpoints are
compatible across them, so it is safe to benchmark a short slice on the
target machine and switch a running sweep to the winner.

Chains themselves thin out with height like `1/ln^L x` — the header's
expected-count estimate is the thing to consult before committing to a big
range. Memory use is modest: roughly `-m` bytes per thread (4 MB default)
plus the sieve-prime tables.

## Notes and limits

- Maximum pattern length is 24 (`MAXLEN`); offsets then reach 2²⁵ − 2.
- `end` must stay below 2⁶⁴ minus the largest offset (~1.8 × 10¹⁹).
- All arithmetic is exact 64-bit; primality is deterministic Miller–Rabin
  with the 7-base set valid for all n < 2⁶⁴ — no probabilistic results.
- `--consecutive` decides consecutiveness by testing every odd number in the
  gaps, which is expensive for long chains (gaps up to 2²⁴); it only runs on
  chains that already passed the length filter, so the cost is usually
  negligible.
- `--brute` exists to validate the sieve: `powgap -l 5 1 3e7 -q` and
  `powgap -l 5 1 3e7 -q --brute` produce identical output.
