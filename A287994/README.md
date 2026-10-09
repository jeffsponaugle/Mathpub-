# A287994 / A290977 tool: a287994

One multithreaded C program that finds terms of
[OEIS A287994](https://oeis.org/A287994) (position of the first time an n-digit
number appears twice in a row after the decimal point of pi) and its companion
[A290977](https://oeis.org/A290977) (the number itself), for a whole range of
lengths in a single pass.

    make            # builds a287994 (cc -O3 -march=native -pthread)
    make test       # checks a(1)..a(6) against ../A197123/pi-1b.txt

## Usage

    ./a287994 -p pi-10b.txt -l 9-12 [-n 10B] [-t 14] [-o results.txt]

| option | meaning |
|---|---|
| `-p FILE` | pi text file (`3.14159...`; non-digits are ignored) |
| `-l SPEC` | block length n, or an inclusive range: `-l 9` or `-l 1-12` (1..64) |
| `-n COUNT` | digits to scan: `500M`, `10B`, `all` (default all) |
| `-t N` | worker threads (default: all cores) |
| `-o FILE` | append confirmed results and the final summary to FILE |
| `-P SECS` | progress interval (default 0.5; min 30 when stderr is not a tty) |
| `-C COUNT` | work chunk size in digits (default 8M) |

Positions are 1-based counting from the first digit after the decimal point
(the usual OEIS convention: a(1) = 24 because "33" starts 24 digits in).
Status goes to stderr on one self-overwriting line (%done, digits searched,
rate, ETA, which lengths are still open); confirmed results print to stdout
the moment they are proven minimal, so a long run shows the small-n answers
within seconds.

## Algorithm

A repeat at start p means `d[i] == d[i+n]` for all i in `[p, p+n)`.  Each
candidate window is checked **back to front**: a mismatch at i rules out every
start in `[i-n+1, i]`, so the scan jumps straight to `i+1`.  On random digits
that is ~1.24 byte-compares per **n digits** advanced, so searching a whole
range of lengths costs only a few compares per digit and runs at memory speed
(~24 GB/s across 14 threads on cached digits; a cold run is disk-bound).

**Loading (built for multi-TB files).** The pi file is mmap'd read-only and
never read up front: the size comes from `fstat`, only the first 64MB are
sniffed to classify the format, and trailing whitespace is trimmed by reading
a few KB backwards.  A pure-digit file (y-cruncher / MIT format) is then
searched in place, each worker validating its chunk's bytes as it scans, so
the file streams from disk exactly once - during the timed search, with live
progress/ETA from the first second.  Pages behind the confirmed frontier are
released (`MADV_DONTNEED`), so a 10T-digit file passes through bounded RAM
instead of filling it.  A non-digit byte past the sniff is a hard error with
the exact digit position.  A file with junk near the start (line wraps etc.)
is compacted into RAM instead (needs ~1 byte per searched digit, honors `-n`),
which requires it to fit.

**Threading / exactness.** Threads pull fixed-size chunks of start positions
in increasing order from a shared counter.  A per-length candidate is only
declared *the* first occurrence once every chunk before it has completed, so
results are exact regardless of thread scheduling; a length drops out of the
scan as soon as it is resolved.  Ctrl-C stops cleanly and reports how far the
completed frontier got (results confirmed before the interrupt remain exact).

**Leading zeros.** A block like `012012` is arguably not a repeated "3-digit
number".  The programs behind the published terms do not exclude leading
zeros (for every known term the two readings agree), so the tool tracks both
the first repeat of *any* n-digit block and the first whose block does not
start with `0`, and reports both if they ever differ.

## Results

Verified against the previously published terms, a(1)..a(8), on
pi-1b.txt / pi-10b.txt.

Terms a(9)-a(12) of both sequences were found with this tool and are now
published in OEIS ("a(9)-a(12) from Jeff Sponaugle, Sep 2026"): a(9) on the
10^10-digit file (2026-08-30), a(10)-a(12) on the 10^13-digit file.

|  n | A287994 (position) | A290977 (value) |
|---:|---:|---:|
|  9 | 1,267,602,413 | 669,926,997 |
| 10 | 11,732,857,329 | 6,024,801,436 |
| 11 | 82,401,465,264 | 73,533,771,380 |
| 12 | 1,701,747,216,801 | 405,287,536,357 |

a(9) was additionally confirmed by a single-threaded rerun and by reading the
raw bytes at file offset 2 + 1267602413 - 1 ("...0016120594 669926997
669926997 5965270717...").

The expected first position for a given n is ~10^n, so a(13) (expected
~10^13, a coin flip in the 10T file if it wasn't already ruled out) and
beyond need a bigger digit file - the scan itself is not the bottleneck,
reading the file is.
