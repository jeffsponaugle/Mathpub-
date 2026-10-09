# lemverify

Verify y-cruncher digit files of the **Lemniscate constant** — as far as
mathematics allows.

y-cruncher's "Lemniscate" is the arc length of the unit lemniscate,

```
L = 2*varpi = 2*pi / AGM(1, sqrt 2) = Gamma(1/4)^2 / sqrt(2*pi)
  = 5.24411510858423962...  (hex 5.3E7E53E7D42C1B9E6...)
```

Files holding `varpi = 2.62205755429...` itself are auto-detected too.

## The fundamental limit (read this first)

**There is no BBP-type digit-extraction formula for the lemniscate
constant.** All known digit-extraction formulas live in the polylogarithm
class (pi, log 2, Catalan, ...); this constant is essentially
`Gamma(1/4)^2`, which is outside it. So unlike pi — where `piverify` can
independently compute hex digits at any position and spot-check a file
anywhere — **digits of L at an arbitrary position cannot be checked without
computing the whole prefix.**

Consequently, the only complete verification of a 3-trillion-digit
computation is the method every published record uses: **run the
computation again with a mathematically independent formula and compare.**
y-cruncher ships several independent formulas for this constant (Zuniga's
hypergeometric series in v0.8.x, plus the classic custom formulas
Gauss `8 ArcSinlemn(1/2) + 4 ArcSinlemn(7/23)`, Sebah, AGM-Pi, Series-Pi).
If the first run used one Zuniga series, verify with a *different* one
(or Gauss/Sebah) and compare outputs — `lemverify compare` does the
comparison and reports the first mismatch.

Everything short of that is partial evidence. `lemverify` implements all
the partial verification that *is* possible:

| mode | what it proves | cost on a 3 TB file |
|---|---|---|
| `verify` | first K digits match an **independent computation** (MPFR, `2*pi/AGM(1,sqrt 2)` — different library *and* formula than y-cruncher/Zuniga); plus 10-digit windows at up to **64 published checkpoint positions reaching 600 billion** (decimal files; values from independent record computations) | seconds–minutes (seek reads; K=10^7 ≈ 22 s, K=10^8 ≈ 4 min of MPFR time) |
| `cross` | the DEC and HEX files **encode the same real number** over their first N digits, verified in both directions by modular base-conversion anchors — one wrong digit anywhere in either covered prefix scrambles the anchors (miss probability ~16^-12 per anchor) | N=2.5×10^8 (default) ≈ 4 min; N=10^9 ≈ 25 min, ~5 GB RAM |
| `compare` | two files are digit-for-digit identical — **the** full verification when the two files come from independent formulas | streaming, I/O-bound (~200 Mdigit/s + disk) |
| `scan` | files are well-formed; digit counts, frequency chi-square, SHA-256 fingerprints | one streaming pass, I/O-bound |

What `verify` + `cross` together do **not** prove: that y-cruncher's
computed *binary* value is correct beyond the independently recomputed
prefix and the checkpoint windows. Both the dec and hex files derive from
the same internal value, so their mutual consistency checks
conversion/storage, not the computation. (y-cruncher also performs internal
modular error-checking during the computation and radix conversion — keep
its `Validation File`.)

## Build

Needs GMP and MPFR (`brew install gmp mpfr` / `apt install libgmp-dev
libmpfr-dev`). On Windows, build under WSL or MSYS2.

```
make
```

## Usage on a 3-trillion-digit y-cruncher output

```bash
# 1. Independent check of the head + all 64 published checkpoints (to 600B):
./lemverify verify "Lemniscate - Dec - Zuniga.txt" --prefix 100000000 -o cert_dec.txt

# 2. Hex file head (no published hex table exists; cross ties it to the dec file):
./lemverify verify "Lemniscate - Hex - Zuniga.txt" --prefix 50000000 -o cert_hex.txt

# 3. Mutual dec<->hex consistency over the first 10^9 digits (~25 min, ~5 GB RAM):
./lemverify cross "Lemniscate - Dec - Zuniga.txt" "Lemniscate - Hex - Zuniga.txt" \
    --span 1000000000 -o cert_cross.txt

# 4. Whole-file integrity + fingerprints (streams all 3 TB once):
./lemverify scan "Lemniscate - Dec - Zuniga.txt"

# 5. THE real verification: recompute with a different y-cruncher formula, then
./lemverify compare "Lemniscate - Dec - Zuniga.txt" "Lemniscate - Dec - Gauss.txt" -o cert_final.txt
```

`--no-hash` skips SHA-256 in certificates (saves a full read of a 3 TB
file). Multiple chunk files are accepted and read in series:
`lemverify verify part0.txt part1.txt ...`.

Other commands:

```bash
./lemverify gen 1000000 --out ref.txt --base hex      # reference digits (MPFR)
./lemverify digits 1000000 --count 10                 # digits at a position
                                                      # (computes the whole
                                                      #  prefix — no shortcut
                                                      #  exists; fine to ~10^8)
```

## How the cross-check works

The decimal file defines the rational `x = 5 + D/10^N` (D = first N
fractional digits as one integer). Then

```
frac(16^m * x) = ((2^(4m) mod 10^N) * D mod 10^N) / 10^N
```

is computable by modular arithmetic, and its leading hex digits must equal
the digits sitting at position m+1 of the **hex file** (an O(1) file seek).
Because multiplication mod 10^N mixes *every* decimal digit into the
result, a single wrong digit anywhere in the covered range changes the
anchor value essentially at random. Anchors are placed near the maximum
usable position `m ≈ N/1.20412` (beyond that, the file-truncation error
would reach the compared digits; 20 guard digits are kept).

The mirror direction does the same from the hex file (`frac(10^r * x_hex)`
mod `2^(4N)`, where the reductions are just bit masks) against decimal
digits read from the dec file. A failure localizes which file is bad: a
corrupt dec file fails only dec→hex anchors, a corrupt hex file only
hex→dec.

This is the same "anchor" machinery as `piverify`'s decimal mode — but
where piverify compared anchors against BBP (independent truth), here the
two files vouch only for *each other*.

## Trust model / self-checks

- Reference digits are computed with MPFR (correct rounding) at 96 guard
  bits + 12 guard digits, truncated like y-cruncher files, with an explicit
  truncation-boundary ambiguity check.
- Every reference computation is automatically cross-checked against the
  built-in first-100-digit constants **and** y-cruncher's distributed
  reference table (64 checkpoint values contributed by independent record
  computations: Watkins, Trueb, Kim, Cutress, et al., up to position
  600,000,000,000). The tool refuses to proceed if its own computation
  disagrees with the table.
- Negative tests: one flipped digit at an arbitrary interior position
  fails `verify` (when inside the prefix), fails `cross` in the matching
  direction, and fails `compare` at the exact position.
- `--seed` makes anchor placement reproducible; `-o cert.txt` writes a
  certificate with SHA-256 fingerprints, every check performed, and an
  interpretation of exactly what was and was not proven.

## Performance (M-series Mac, single-threaded GMP/MPFR)

| operation | size | time |
|---|---|---|
| `verify --prefix` (MPFR recompute) | 10^7 digits | ~22 s |
| | 10^8 digits | ~4 min |
| checkpoint windows on a 3 TB file | 64 seeks | seconds |
| `cross --span` | 5×10^7 | ~40 s |
| | 2.5×10^8 (default) | ~4 min |
| | 10^9 | ~25 min, ~5 GB RAM |
| | 10^10 | hours, ~50 GB RAM |
| `compare` / `scan` | any | disk speed (~0.2–2 GB/s) |
