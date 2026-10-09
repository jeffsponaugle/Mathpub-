# gapsearch — prime chains following an arbitrary gap pattern

Finds primes `p` such that `p, p+g1, p+g1+g2, …` are all prime, for a chosen
gap pattern. `powgap.c` (doubling gaps only) is the earlier, fixed-pattern
version; `gapsearch.c` supersedes it and reproduces it exactly with the
default `--gaps 2^n`.

## Gap patterns

```
-g 2^n         2,4,8,16,32,...      doubling
-g 3^n         2,6,18,54,162,...    2*r^k for any r (5^n, 10^n, ...)
-g 2n          2,4,6,8,10,...
-g 2p          4,6,10,14,22,...     twice each prime
-g 2p@3        6,10,14,22,26,...    twice each prime >= 3
-g p           2,3,5,7,11,...
-g 30x5        30,30,30,30,30       constant gap (AP / CPAP)
-g kt12        densest admissible 12-tuple (kt2 .. kt18)
-g 4,6,10,14   explicit list
```

## Pattern catalogue

Admissibility ceiling and expected chain count for `p < 2^64`, at the ceiling
length. "consec cost" is roughly what the `-c` requirement multiplies the
yield by.

| pattern | max length | expected chains < 10^18 | consec cost | notes |
|---|---|---|---|---|
| `2n`   | unbounded | ~10^4 at len 13 | severe past len 8 | OEIS A016045 |
| `2^n`  | unbounded | 1568 at len 12  | severe past len 8 | doubling |
| `3^n`  | unbounded | 9.9e4 at len 11 | hopeless past len 6 | wide open |
| `2p`   | **3** | — | — | blocked mod 3; only `3,7,13,23,37,59` |
| `2p@3` | 9 | plentiful | mild | |
| `2,6,30,210,2310,...` | 8 | 6.3e8 | hopeless | primorial gaps |
| `30x5` (CPAP-6) | 6 | plentiful | moderate | smallest is 121174811 |
| `210x9` | 10 | 2.2e6 | brutal | CPAP-10 record is 93 digits |
| `kt12` | 12 | ~550 | **~2x only** | prime 12-tuplet |
| `kt14` | 14 | 1.8 | ~2x | genuine 64-bit frontier |
| `kt16` | 16 | 0.017 | ~2x | needs bignum |

Patterns whose gaps grow fast (Fibonacci, squares, triangular, factorial-based)
almost all die at length 3-6 on a mod-3 or mod-5 obstruction. Run them with
`-l 2` first and read the admissibility line before committing to a search.

The k-tuplet patterns are the sweet spot for `-c`: the whole chain spans 42-60
integers, so demanding consecutive primes only costs a factor of about two,
versus factors of 10^-20 for wide patterns.

## Admissibility — read this first

Before searching, the program checks whether the pattern *can* produce chains
of the requested length. If the first `L` offsets cover every residue class
mod some prime `q`, then one member of every chain is divisible by `q`, and no
chains of length `L` exist at all — apart from the degenerate ones where a
member *is* `q`. The program prints the largest admissible length and, when
you ask for more than that, checks the degenerate cases directly instead of
burning CPU on an impossible search.

Worked example — `--gaps 2p` (gaps 4, 6, 10, 14, …, offsets 0, 4, 10, 20, …):

```
offsets mod 3:  0, 1, 1, 2, ...
```

By the fourth member the offsets hit all of 0, 1, 2 mod 3, so one of
`p, p+4, p+10, p+20` is always a multiple of 3. Maximum admissible length: 3.
The single exception is `p = 3` itself:

```
3, 7, 13, 23, 37, 59        (the next term, 85 = 5·17, is composite)
```

and that is the *only* chain of length ≥ 4 in this pattern, anywhere. Nudging
the pattern to `2p@3` (gaps 6, 10, 14, 22, …) dodges the obstruction and is
admissible to length 9.

## Build

```
gcc   -O3 -march=native -o gapsearch gapsearch.c -lpthread -lm   # Linux / x86
clang -O3 -mcpu=native  -o gapsearch gapsearch.c -lpthread -lm   # macOS
```

On macOS `-march=native` is not accepted by Apple clang on arm64 — use
`-mcpu=native`, or just `-O3`.

No external libraries. (primesieve isn't needed — the search never enumerates
primes, it sieves the constellation directly, which is far faster here.)

## Usage

```
gapsearch [options] <start> <end>

  -g, --gaps SPEC    gap pattern (default 2^n) -- see above
  -l, --len N        minimum chain length to report (default 6)
  -t, --threads N    worker threads (default: online CPUs)
  -c, --consecutive  require members to be CONSECUTIVE primes
  -b, --bound N      small-prime sieve bound (default: auto, by length)
  -m, --mem N        bitmap bytes per thread (default 4M)
      --brute        naive verification mode, no sieve (for testing)
  -q, --quiet        no header/progress on stderr
```

`start` and `end` accept `K M B/G T P E` suffixes, decimals and scientific
notation, and `,`/`_` separators: `250B`, `1.5T`, `3e15`, `1_000_000`.

```
./gapsearch -l 10 -t 16 1T 10T             # doubling gaps, length 10+
./gapsearch -g 2p@3 -l 8 -c 1 100B         # 2*prime gaps, consecutive primes
./gapsearch -g 2n -l 12 -t 16 1P 2P        # the 2,4,6,8,... chain search
./gapsearch -g 2p -l 4 1 1T                # reports the obstruction, finds p=3
./gapsearch -l 5 --brute 1 3M              # cross-check the sieve
```

Results stream to stdout in ascending order of `p` regardless of thread count;
progress and stats go to stderr, so `> chains.txt` gives a clean file.

## How it works

For the default doubling pattern the offsets are `2^(k+1) − 2`, so their
residues mod a prime `q` are governed by the multiplicative order of 2 — which
makes the wheel unusually effective:

| q  | offsets mod q          | allowed residues for p |
|----|------------------------|------------------------|
| 2  | {0}                    | 1 of 2                 |
| 3  | {0,2}                  | 1 of 3                 |
| 5  | {0,2,1,4}              | 1 of 5 (len ≥ 4)       |
| 7  | {0,2,6}                | 4 of 7                 |
| 11 | 10 values              | 1 of 11 (len ≥ 10)     |

Since `ord_q(2) ≤ q−1`, at least one residue always survives — the doubling
pattern is admissible at *every* length, so arbitrarily long chains should
exist. Other patterns are not so lucky; see the admissibility section above.

1. **Wheel.** Residues mod `W = 2·3·5·7·11·13·17·19 = 9,699,690` that survive
   the pattern are enumerated by CRT (972 of them for length 10 — a density of
   1 in ~10,000). Everything downstream only ever touches those.
2. **Segmented constellation sieve.** Each block is a bitmap of `nres × JLEN`
   bits, one row per wheel residue class. For each small prime `q` and each
   offset `o_k`, the killed positions inside a row form an arithmetic
   progression of stride `q`. The start index is
   `j ≡ −(base + r_i)·W⁻¹ − o_k·W⁻¹ (mod q)`, so the `−o_k·W⁻¹` terms are
   precomputed once per prime (deduplicated, since offsets collide mod small
   `q`) and the per-row work is one modular multiply plus an add per offset.
3. **Confirmation.** Survivors get deterministic 64-bit Miller–Rabin
   (7-base set, valid to 2⁶⁴), members tested in order so failures bail early.
   Confirmed chains are then extended past the requested length for free.

Throughput on one core: ~1×10¹¹ integers/second at length 10 near 10¹². It
scales linearly with threads (blocks are handed out by an atomic counter).

## `--consecutive`

By default, a chain only requires the listed members to be prime; other primes
may sit between them. `-c` additionally requires that no prime lies inside any
gap, so the gaps are genuine prime gaps (the stricter definition).

Be aware this hits a wall fast. Length `L` requires roughly `2^L` interior odd
numbers to all be composite, and the probability of that near `N` is about
`exp(−2^L / ln N)`. Length 7 is easy, length 8 is findable, length 9 is
already brutal, and length 10 needs a maximal-gap-scale event of 512. Length
7 chains show up every few billion; a length-8 consecutive chain is a
worthwhile hunt.

Without `-c` the frontier is much further out — expected counts near `N`:

| length | first expected around |
|--------|-----------------------|
| 10     | 10¹¹                  |
| 11     | 10¹³                  |
| 12     | 10¹⁵                  |
| 13     | 10¹⁷–10¹⁸             |

The header prints the singular series and an expected-count estimate for the
range you asked for, so you can size a run before starting it.

## Tuning

* `--bound` is the main knob. Higher = fewer survivors to Miller–Rabin, more
  sieve setup per block. The default scales with `--len` (2²¹ at length ≤ 4
  down to 2¹⁷ at length ≥ 9); raise it if the "candidates confirmed" count in
  the summary is large relative to runtime.
* `--mem` sets the per-thread bitmap (default 4 MB). Bigger blocks amortize
  per-prime setup better; smaller blocks keep the bitmap in L2/L3 and give
  finer load balancing. The block size is automatically reduced when the range
  is too small to keep every thread fed.
* Chains found are extended beyond `--len` automatically, so searching with
  `-l 10` will still report a length-12 chain as length 12 — no need to run
  multiple lengths.

## Verification

`--brute` runs the same predicate with no sieve at all. Outputs match exactly
over tested ranges, including across the low-range/sieve boundary:

```
./gapsearch -g 2p -l 3 -q --brute 1 500K > b.txt
./gapsearch -g 2p -l 3 -q -t 8      1 500K > s.txt
diff b.txt s.txt
```

Numbers below the sieve bound are always handled by the brute path, since a
chain member could otherwise coincide with one of the sieving primes.
