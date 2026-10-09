# pi_missing

A multithreaded scanner for the decimal digits of Pi that computes terms of
five OEIS sequences in a single pass over a digit file:

| Sequence | Definition | Flag |
|---|---|---|
| [A228988](https://oeis.org/A228988) | a(n) = smallest missing number in the first 10^n digits after the decimal point | always on |
| [A260627](https://oeis.org/A260627) | a(n) = largest n-digit number missing in the first 10^n digits after the decimal point | `-L` |
| [A153221](https://oeis.org/A153221) | numbers k whose digit string occurs at (1-based) position k−3 | `-S` |
| [A153223](https://oeis.org/A153223) | same, at position k−4 | `-S` |
| [A153224](https://oeis.org/A153224) | same, at position k−5 | `-S` |

"Missing" means the number's decimal string (no leading zeros) does not occur
as a substring of the digit stream. Example: a(1) of A228988 is 0 because the
first 10 digits after the decimal point, `1415926535`, contain every digit
string `1`–`9` but no `0`.

## How it works

The digit file is memory-mapped and a window slides over every position. At
each position that does not start with `0`, the values of the 1-, 2-, 3-, …
digit windows starting there are computed incrementally (`v = v*10 + digit`)
and each value sets one bit in a shared bit array. At every 10^k boundary the
lowest clear bit is A228988(k). Correctness details:

- **Leading zeros don't count**, but any occurrence of the canonical digit
  string does: an occurrence of `014523` contains `14523` starting one digit
  later, so skipping zero-led windows loses nothing.
- **Windows never cross a boundary they are being counted for.** Windows are
  clamped at the current segment end; the last few positions before a
  boundary are re-scanned in the next segment so their longer windows count
  toward later boundaries only.
- **Bounded tracking, guarded.** Only numbers below a limit (`-m`) get a bit.
  The limit defaults to max(digits/25, 2e6) — about 4x the answer, since
  empirically a(n) ≈ 1.0000x · 10^(n−2). If every tracked number turns out to
  be present, the program refuses to answer and asks for a larger `-m`; it
  never reports a wrong term.
- **A260627** (`-L`) rides the same loop: for each window length l it keeps a
  small bit window covering only the top `-W` values below 10^l (default
  2^20; the answers are empirically within ~10 of 10^n − 1). At boundary 10^k
  the highest clear bit of the length-k window is the term, with the same
  refuse-don't-guess behavior if the window saturates.
- **Self-locating searches** (`-S`) need no memory at all: each thread
  carries the decimal strings of i+4, i+5, i+6 (the candidates for offsets
  k−3, k−4, k−5 at 0-based position i) as ripple-incremented counters and
  compares them against the stream at every position. Matches print
  immediately with their position. Print order can race slightly between
  threads; sort by position before publishing.

### Performance design

- **Threads** partition each ~1-gigadigit slice of the stream; bits are set
  with an atomic OR guarded by a plain read (almost all bits are already set,
  so the atomic rarely executes).
- **Smallest-missing lower bound:** once a boundary reports smallest-missing
  = s, every value below s is provably present, so later segments skip the
  bit-array access for window values < s. In the final segment of a big run
  this removes ~90% of the random memory traffic — the dominant cost — and is
  worth ~25x wall-clock on a 10^10 run.
- The workload is memory-latency-bound: throughput scales nearly linearly
  with threads until DRAM parallelism saturates. Use all cores (the default).

### Checkpointing

`-c FILE` snapshots the full state (bit arrays + position + loop state) every
`-i` seconds (default 1800), written atomically via `FILE.tmp` + rename. With
`-c` active, **Ctrl-C / SIGTERM write a final checkpoint and exit cleanly**
(status 3). Resume with the same command plus `-r` — on the same machine or
another one (copy the checkpoint and use the same digit file and identical
`-n/-m/-W/-L/-S` options; a header fingerprint refuses mismatches). At most
one checkpoint interval of work is repeated; re-marking is idempotent, and an
interrupted+resumed run has been verified to produce output identical to an
uninterrupted one.

## Input file format

Plain text digits of Pi. A leading `3.` is detected and skipped
automatically (the sequences count digits *after* the decimal point);
trailing whitespace at end of file is ignored. Everything else must be pure
digit bytes — the program aborts with the file offset if it hits anything
else. If your file has line breaks, strip them first:

```bash
tr -cd '0-9' < pi_raw.txt > pi_clean.txt
```

(then pass `-o 1` if that leaves a leading `3`).

## Building

```bash
make            # or: cc -O3 -march=native -pthread -o pi_missing pi_missing.c
```

Plain C11 + POSIX (pthreads, mmap); builds on macOS and Linux.

## Running

```bash
./pi_missing [-L] [-S] [-n digits] [-m limit] [-W window] [-t threads]
             [-o skip] [-c ckpt [-i secs] [-r]] pifile
```

| Option | Meaning |
|---|---|
| `-L` | also compute A260627 at each power-of-10 boundary |
| `-S` | also search for A153221 / A153223 / A153224 terms |
| `-n digits` | use only the first so-many digits (default: whole file). Sizes accept `k/M/G/T` suffixes and `1e12` notation everywhere |
| `-m limit` | track numbers 0..limit−1 for A228988; memory = limit/8 bytes (default max(digits/25, 2e6)) |
| `-W window` | A260627 per-length window below each 10^l (default 2^20) |
| `-t threads` | worker threads (default: all CPUs) |
| `-o skip` | skip leading bytes (a leading `3.` is auto-skipped) |
| `-c file` | write checkpoints to `file`; SIGINT/SIGTERM save and exit |
| `-i secs` | checkpoint interval (default 1800) |
| `-r` | resume from the `-c` checkpoint |

Results go to stdout (flushed line by line — safe to redirect). Progress goes
to stderr: on a terminal a self-updating line; redirected to a log, one line
per minute:

```
 29.9%  3082/10300G  25.3 Md/s  elapsed 1d 10:12  10^13 in 2d 03:14  ETA 3d 08:45 (Sep 05 14:32)
```

(position, smoothed current rate, elapsed, time to the next boundary report,
and overall ETA with projected finish time.)

### Examples

Quick verification against the known OEIS terms (any pi file with ≥ 10^7
digits):

```bash
./pi_missing -L -S -n 10M pi.txt
```

Full production run on a 10-trillion-digit file, checkpointed:

```bash
./pi_missing -L -S -m 2e11 -c ck.bin pi-10t.txt > results.txt
```

`-m 2e11` uses a 25 GB bit array (checkpoints are the same size — keep ~50 GB
free for the atomic swap). A228988's answer at 10^13 digits is expected near
1.00000x·10^11, so 2e11 is a 2x margin; the default (digits/25 = 4e11, 50 GB)
is fine too if RAM allows, and marks somewhat more than it needs to.

### Sizing and runtime intuition

- a(n) of A228988 ≈ 1.0000x · 10^(n−2), so `-m` around 2·10^(n−2) is a
  comfortable margin at any scale.
- Memory: `-m`/8 bytes. Never let the bit array exceed physical RAM.
- Runtime is dominated by the last decade of digits and by random DRAM
  access. Reference points: 10^10 digits ≈ 26 s on an Apple-silicon laptop
  (14 threads, 50 MB array); 10^13 digits with a 25 GB array lands in the
  0.5–5 day range depending on core count and memory system — thread count
  is the lever.

## Verification

Every known OEIS term of all five sequences reproduces exactly: A228988
a(1)–a(11), A260627 a(1)–a(7) and a(9), A153221 a(1)–a(8), A153223
a(1)–a(10), A153224 a(1)–a(6). Digit data was cross-checked against an
independent Chudnovsky computation, and interrupted+resumed runs were
verified bit-identical to clean runs.

Results from runs with this tool (as of Sep 2026, first 10^12 digits, runs to
10^13 in progress):

- **A260627(8) = 99999991** — corrects the published 99999992, which in fact
  occurs at position 3389380; 99999991 first occurs at position 209136149,
  and 99999993–99999999 all occur within the first 10^8 digits.
- New terms: A228988(12) = 10000031203; A260627(10)–(12) = 9999999999,
  99999999998, 999999999998; A153221 a(9)–a(11) = 44159961449, 73514522243,
  249973441618; A153223 a(11)–a(12) = 4927138295, 29478491757; A153224
  a(7)–a(9) = 2544423437, 7730967174, 870163762036.
