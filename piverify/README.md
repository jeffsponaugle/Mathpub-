# piverify

Verify a file of pi digits using the **BBP digit-extraction formula** — the
Bailey–Borwein–Plouffe formula (1996), which computes digits of pi at an
arbitrary position without computing any of the earlier digits:

```
pi = sum_{k>=0} 16^-k * ( 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) )
```

Multiplying by `16^n` and reducing each term with modular exponentiation
yields `frac(16^n * pi)` — the hex digits starting at position `n+1` — in
O(n) small operations and O(1) memory.

The tool is written in C. `make` builds two binaries:

- **`piverify`** — the full tool (`verify` / `digit` / `gen`). Needs GMP
  for the decimal-anchor big-integer step and for `gen`
  (`brew install gmp` / `apt install libgmp-dev`).
- **`bbp`** — a minimal standalone BBP digit extractor (no dependencies
  beyond pthreads), useful on machines without GMP.

Layout: [bbp.c](bbp.c)/[bbp.h](bbp.h) is the BBP engine — multithreaded
(pthreads) with Montgomery modular multiplication on native 64/128-bit
integers and exact 128-bit fixed-point accumulation (no float round-off).
[piverify.c](piverify.c) is the tool; [bbp_main.c](bbp_main.c) the thin
CLI. [piverify.py](piverify.py) is the original Python reference
implementation — no longer required, kept because it independently
cross-validates the C code (both produce identical results).

BBP is embarrassingly parallel — every term of the sum is independent —
so threads pull chunks of the term range from an atomic counter and
scaling is essentially linear in physical cores (7.8x on 14 threads of a
14-core M-series Mac, which has a mix of performance/efficiency cores).

## The catch, and how this tool handles it

**BBP extracts *hexadecimal* (base-16) digits, not decimal.** No practical
decimal digit-extraction formula is known (Plouffe published a base-10
formula in 2022, but it costs more than recomputing pi, so it is not usable
for spot checks). This forces two different verification strategies:

### Hex files → true random spot checks

Digits at random positions are compared directly against BBP. Each check
verifies exactly the digits sampled. The tool always includes a check at
the **end of the file**, because that is the decisive one for files produced
by a pi computation: any arithmetic error propagates to every later digit,
so a correct tail plus random interior samples is precisely how world-record
computations are verified in practice.

### Decimal files → whole-file "anchor" verification

You cannot spot-check decimal digits with BBP, but you can do something
stronger. The file defines a rational number `x = 3 + D/10^N` (D = all N
digits as one integer). The tool computes

```
frac(16^m * x)  =  ( (2^(4m) mod 10^N) * D mod 10^N ) / 10^N
```

and compares its leading hex digits against BBP's `frac(16^m * pi)`, at
anchor positions `m` chosen randomly near the maximum usable value
`~0.8305*N` (the base-conversion sensitivity horizon: `log10/log16`).

Because base conversion mixes **every** decimal digit into the result, a
single wrong digit at *any* covered position changes the anchor value by an
essentially random amount mod 1 — a mismatch with probability
`~1 - 16^-t` per anchor (t = 12 compared hex digits by default, 3
independent anchors by default). This reads the whole file once, but unlike
sampling it is a genuine whole-file integrity check: it catches even one
flipped bit from storage corruption. The only blind spot is the last
~few-hundred digits, which sit beyond the anchors' resolution (the tool
reports the exact covered range).

## Usage

```
make
./piverify verify pi_digits.txt               # auto-detects dec/hex
./piverify verify pi.hex --checks 16 --seed 1
./piverify verify pi.txt --anchors 5 --digits 16
./piverify digit 1000000 --count 8            # BBP at any position (1-based)
./piverify digit 1000000000 -t 8              # limit to 8 threads
./piverify gen 1000000 --base dec --out pi_1m.txt  # test data (Chudnovsky)
./bbp -d 12 -t 14 -v 1000000000               # standalone digit extraction
```

`verify` and `digit` take `-t N` / `--threads N` to set the BBP worker
thread count (default: all cores).

### Verification certificates

`verify --out cert.txt` writes a plain-text certificate documenting the
specific verification performed:

- tool version, UTC timestamp, host, and the exact command line
- the RNG seed used (recorded even when auto-generated, so the identical
  run can be reproduced with `--seed`)
- every input file with its size and **SHA-256 fingerprint**, binding
  the certificate to the exact bytes that were verified (skip the
  hashing pass over huge files with `--no-hash`)
- each file's digit range in the series and the total digit count
- every check performed with both computed values (file vs BBP) and its
  PASS/FAIL status, including anchor positions and file-boundary checks
- an interpretation section stating precisely what the result does and
  does not guarantee (e.g. the anchor coverage range and the uncovered
  tail digits), and the overall RESULT

The SHA-256 implementation is built in (validated against `shasum -a
256`), so no new dependencies are introduced.

### Multi-file digit sets

`verify` accepts multiple files and reads them in series as one
continuous digit stream, in the order given on the command line — the
usual layout for very large digit sets stored as chunks:

```
./piverify verify pi_part_00.txt pi_part_01.txt pi_part_02.txt
```

The first file carries the `3.` prefix; later files are raw digit
chunks (any formatting; files may split mid-line, since whitespace is
ignored). Hex mode automatically adds a spot check at every file
boundary — the place where a missing, truncated, or out-of-order chunk
shows up — and labels those `[file boundary]` in the output. Decimal
mode needs nothing extra: the anchors cover the concatenated stream, so
any chunk swap or gap scrambles them.

Accepted file formats: raw digits with or without a leading `3.` or `3`;
whitespace/newlines/commas ignored. Contiguous files are read with O(1)
seeks; formatted files fall back to streaming. Exit code 0 = pass,
1 = fail, 2 = error.

The first ~100 digits are always checked against a built-in constant first,
which catches format and off-by-one problems immediately.

`gen` computes reference digits with the Chudnovsky series — a completely
independent algorithm — so `gen` + `verify` cross-validate each other.

## Performance (M-series Mac, 14 cores)

BBP cost is linear in position; the whole-file anchor also needs one
modular exponentiation over the full file.

| position | Python engine (14 procs) | C engine (14 threads) |
|---|---|---|
| 10^6 | 0.3 s | < 0.05 s |
| 10^8 | 29 s | 5.6 s (46 s on 1 thread) |
| 10^9 | ~5 min | 70 s |
| 10^12 (one trillion) | ~4–5 days | **~20 hours** (extrapolated; linear) |

So a trillion-digit hex spot check, or the BBP side of a trillion-digit
decimal anchor, is an overnight run on this machine with the C engine.
(Dedicated GPU implementations do the same in under an hour.)

The GMP big-integer side is far from the bottleneck: the decimal anchor's
convert step on a 1,000,000-digit file takes 0.1 s (vs 13 s with Python
ints), and `gen` produces 1,000,000 reference digits in 0.2 s (vs 7 s).
For a trillion-digit decimal file the anchor still needs the whole value
in memory (~420 GB) — at that scale verify the hex representation
instead, which streams with near-zero memory.

### Very large decimal files (`--prefix`)

Exhaustive decimal anchor verification must hold the entire verified
span in memory as one big integer — roughly **3.2 bytes of working
memory per digit** (so a 10-trillion-digit file would need ~33 TB, and
GMP additionally caps a single integer at ~38 billion digits). When a
file is too large, `verify` now says so up front and suggests options
instead of failing in `malloc`.

For such files use `--prefix N` (suffixes `k`/`m`/`g`/`t` accepted):

```
./piverify verify /bigoutput/pi-10t.txt --prefix 9g
```

This anchor-verifies the first N digits — still *exhaustive within the
prefix* (any single wrong digit there is caught) — reads nothing beyond
it, and reports the unverified remainder honestly on stdout and in the
certificate. Pick N ≈ RAM/4. Expect the anchors to take a while at this
scale: at 9 billion digits each anchor is a ~7.5-billion-position BBP
computation (minutes to tens of minutes) plus a GMP powm of similar
wall time.

There is no low-memory way to exhaustively verify the *tail* of a huge
decimal file — base conversion fundamentally mixes all leading digits
into any checkable quantity. For full coverage at multi-terabyte scale,
verify a **hexadecimal** representation of the same computation
instead: BBP spot checks (including the decisive end-of-stream check)
stream at any size with near-zero memory. If the file came out of
y-cruncher, its hex output (or the computation's own BBP check) is the
natural companion.

Tip: with `--out`, the certificate's SHA-256 pass reads every byte of
every input file — on a multi-TB file use `--no-hash` unless you want
that.

## Trust model

- BBP digits are computed with exact integer fixed-point arithmetic
  (96 guard bits), not floats, so results are reliable at any position.
- Verified four ways: against published constants, against the original
  BBP paper's value at position 10^6 (`26C65E52...`), against
  independently generated Chudnovsky digits at random anchors, and
  C-vs-Python cross-checks (identical digits, anchors, and gen output).
- Negative tests: a single flipped digit mid-file (decimal) or at a sampled
  position (hex) produces FAIL.
