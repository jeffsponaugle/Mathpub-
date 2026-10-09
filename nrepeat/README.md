# nrepeat — extending OEIS A331881 / A331882

- [A331881](https://oeis.org/A331881)(n): the first n-digit substring to occur
  n times in the decimal expansion of the fractional part of Pi.
- [A331882](https://oeis.org/A331882)(n): the number of digits of the
  fractional part of Pi needed to contain those n occurrences (the 1-based end
  position of the n-th occurrence).

Known terms (n = 1..8):

```
A331881: 1, 26, 446, 2796, 86538, 872117, 1591292, 66416662
A331882: 1, 22, 219, 1805, 25499, 168882, 3566679, 29325629
```

All window start positions count (overlapping occurrences included), matching
the OEIS examples ("26" at positions 6 and 21, so a(2)=22).

## New results

n = 9 (computed 2026-08-31 on 600M MPFR-generated digits; independently
confirmed by `tools/verify_exhaustive.py` — numpy count of all 481,426,663
windows in the prefix shows exactly one 9-digit substring reaching 9
occurrences):

```
A331881(9) = 858275944
A331882(9) = 481426671
occurrence starts: 8364516 108721072 118923087 209578974 347121287
                   365354466 397701128 426812266 481426663
```

n = 10 (computed 2026-09-01 on the 10T-digit file, house1, -t 30, 50s;
single candidate):

```
A331881(10) = 5762246467
A331882(10) = 4519908297
occurrence starts: 62083470 334368659 734076310 901367000 1611279195
                   2368347516 3645033154 4282821718 4472382985 4519908288
```

n = 11 (computed 2026-09-01 on the 10T-digit file, house1, -t 30 -c 4,
10m43s; single candidate):

```
A331881(11) = 47982672136
A331882(11) = 53506468331
occurrence starts: 10501601989 12403732092 13269952849 15713942497
                   16424571944 22902234759 25685796054 31751231100
                   42257680350 51109092955 53506468321
```

n = 12 (computed 2026-08-31 on the 10T-digit file, house1, -t 30 -c 4 -P 2,
2h32m; single candidate, 12 positions verified by the tool's exact
verification pass):

```
A331881(12) = 100085029093
A331882(12) = 579330690478
occurrence starts: 51738833731 68671221941 75604597502 160932526320
                   231367395165 314009479133 324889612204 445320004727
                   447112780643 468314848852 526877014657 579330690467
```

All four terms n = 9..12 are computed; n = 9 was reproduced identically on
both independent pi computations (MPFR locally and the 10T y-cruncher-style
file on house1). See `oeis-submission.md` for the prepared OEIS edit text.

## Build

```
make            # c++ -O3 -march=native -std=c++17 -pthread
```

Single source file, no dependencies. Linux and macOS.

## Usage

```
nrepeat -n N [options] <pi-digit-file>

  -n N            compute term N (required unless --selftest)
  -o FILE         append the summary report to FILE
  -t T            worker threads (default: all hardware threads)
  -c 4|8          counter bits (default 8; 4 halves table memory, needs n<=15)
  -P K            split the value space into K partitions scanned in K passes;
                  divides table memory by K at the cost of K reads of the
                  digit prefix (default 1)
  -b MDIGITS      max block size in Mdigits (default 256)
  -m DIGITS       cap the number of fractional digits scanned
  --selftest      recompute n=1..8 and compare with the known terms
  --progress SEC  progress interval on stderr (default 1.0)
  -q              quiet (no progress output)
```

The pi file is plain text. Every non-digit byte (newlines, spaces, the decimal
point) is ignored, and a leading `3.` (or a bare leading `3` before `141...`)
is auto-detected and skipped so the stream is the fractional part
`1415926535...`. The first 20 digits are sanity-checked against Pi and a
warning is printed on mismatch. The file must be seekable (the verification
phase rescans it), so no pipes/stdin.

Progress goes to stderr (percent done, digits/s, elapsed, ETA, candidate
count); the summary report goes to stdout and, with `-o`, is appended to the
given file. Exit codes: 0 found, 2 not found within the file, 130 interrupted
(Ctrl-C stops cleanly at the next block boundary).

Typical runs:

```
./nrepeat --selftest pi.txt                       # verify n=1..8 first
./nrepeat -n 9  -o results.txt pi.txt
./nrepeat -n 10 -o results.txt pi.txt
./nrepeat -n 11 -c 4 -o results.txt pi.txt        # 50 GB table
./nrepeat -n 12 -c 4 -P 8 -o results.txt pi.txt   # 62.5 GB table, 8 passes
./nrepeat -n 13 -c 4 -P 32 -o results.txt pi.txt  # 156 GB table, 32 passes
```

## How it works

**Detection pass.** The digit stream is processed in blocks (adaptive: 1 Md
doubling up to 256 Md), each block split across threads. Every n-digit window
value indexes directly into an exact table of 10^n saturating atomic counters
— no hashing, so no collisions and no false positives. The moment any counter
reaches n, that substring becomes a candidate; the scan stops at the end of
that block. Any substring whose n-th occurrence lies in the scanned prefix
necessarily has count ≥ n there, so the true answer is always among the
recorded candidates.

**Partitioning (`-P K`).** When the full table doesn't fit in RAM, the value
space [0, 10^n) is split into K contiguous ranges and the file is scanned K
times, each pass counting only the windows whose value falls in its range
with a full-resolution table of 10^n/K slots. Each pass is exact, so no
information is lost; the answer is the minimum completion position across
partitions. Once one partition has produced a completion at position E,
every later partition only scans to E-1 (a later partition only matters if
it beats the best so far), so passes after the first are usually cheaper.
Extra cost is K sequential reads of the prefix — I/O, not compute.

(Why not a low-bit prefilter — e.g. 2-bit counters, dump every substring seen
3+ times, then rescan for those? At n=12 the scan runs to ~5.3e11 digits, so
each 12-digit string is seen λ≈0.53 times on average; P(Poisson(0.53) >= 3)
≈ 1.7%, i.e. ~1.7e10 survivors — a ~130 GB candidate list that the rescan
would have to test every window against. A threshold that low filters almost
nothing, and 2-bit counters can't hold a higher one. Partitioning reaches the
same memory target with exact counting and no intermediate files.)

**Verification pass.** The table is freed and the prefix is rescanned,
tracking only the (few) candidates and keeping each one's n smallest
occurrence end positions. The winner is the candidate whose n-th occurrence
ends first; that end position is A331882(n) and the substring is A331881(n).
The report also lists the n occurrence start positions for the OEIS example
line.

The inner loop is branch-light: the window value rolls with
`v = (v - d_out * 10^(n-1)) * 10 + d_in` (no division), and counter lines are
software-prefetched a few thousand windows ahead. Throughput is dominated by
random access into the counter table: roughly 300–400 Mdigits/s on an M-series
laptop while the table fits in RAM.

## Resource planning

Expected completion position from the Poisson/birthday estimate
N(n) ≈ (n!)^(1/n) · 10^(n-1) (matches the known terms within a factor ~1.4):

| n  | est. digits needed | table (-c 4, P=1) | with -P                  |
|----|--------------------|-------------------|--------------------------|
| 9  | ~4.2e8             | 0.5 GB            | —                        |
| 10 | ~4.5e9             | 5 GB              | —                        |
| 11 | ~4.9e10            | 50 GB             | -P 4: 12.5 GB, 4 passes  |
| 12 | ~5.3e11            | 500 GB            | -P 8: 62.5 GB, 8 passes  |
| 13 | ~5.7e12            | 5 TB              | -P 32: 156 GB, 32 passes |

Pick K so the per-partition table fits comfortably in RAM. Extra cost is K
sequential reads of the ~N-digit prefix (and usually less: later partitions
are capped at the best completion found so far). For n=13 with -P 32 that is
up to ~32 x 6 TB ≈ 190 TB of sequential reads — a day or two at a few GB/s,
with a 10T-digit file covering the expected completion point.

The table is allocated with `MAP_NORESERVE` and touched lazily, but a run to
the expected completion point touches essentially all of it — plan real RAM
for the per-partition table size, and prefer `-c 4` for n ≥ 11.

## Verification

`--selftest` recomputes n=1..8 and compares both sequences against the known
OEIS values (needs ≥ 30M digits of input). Verified passing against 40M MPFR-
computed digits, with identical results in 8-bit and 4-bit counter modes.
