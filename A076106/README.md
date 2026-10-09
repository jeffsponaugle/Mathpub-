# a076106 — Pi Prime Search (OEIS A076106 / A076130)

Computes two related OEIS sequences from a file of Pi digits:

- **[A076106](https://oeis.org/A076106)**: of all n-digit primes, the one whose
  first occurrence in the decimal digits of Pi (ignoring the initial 3) is
  latest: `7, 73, 373, 9337, 35569, 805289, 9271903, ...`
- **[A076130](https://oeis.org/A076130)**: the 1-indexed position where that
  prime first appears: `13, 299, 5229, 75961, 715492, 11137824, 135224164, ...`

Example: of the 2-digit primes (11..97), the last to show up in Pi is 73, at
position 299.

## Build

No dependencies:

```
make
```

## Usage

```
a076106 [-n maxN] [-d maxdigits] [-t threads] [-c statefile] [-i ckptdigits] pifile

  -n maxN       compute a(1)..a(maxN), 1..12 (default 8)
  -d maxdigits  stop after scanning this many digits (default: whole file)
  -t threads    threads for prime sieving (default: all CPUs)
  -c statefile  checkpoint file (resume if it exists; save progress)
  -i ckptdigits digits between periodic checkpoints (default 100e9, 0 = off)
  -o outfile    append discoveries (timestamped, as they happen) and the
                final results table to this file
  pifile        text file of Pi digits
```

## Checkpointing

With `-c state.ckpt`, progress is saved to that file periodically (every
`-i` digits), at end of input, and on Ctrl-C/SIGTERM (saves are atomic:
tmp file + rename). Rerunning with the same `-c` file resumes: the tool
fast-skips the digits file to the saved position and continues -- so a scan
that exhausted one digits file can be continued later against a longer file
without redoing the work. Requirements on resume: the same `-n`, and a
digits file whose digit stream matches the original from the beginning
(formatting/whitespace may differ). The state file holds each n's
seen-bitmap plus counters: ~1.4 GB for maxN=10, ~14 GB for maxN=11,
~140 GB for maxN=12. Consistency (prime counts, bitmap popcounts) is
verified on load.

Typical long-run workflow:

```
./a076106 -n 11 -c state11.ckpt pi-1t.txt     # runs out of digits, saves
./a076106 -n 11 -c state11.ckpt pi-4t.txt     # resumes at the 1T mark
```

The digits file is **streamed**, never loaded into memory, so files of any
size work (tested design point: multi-terabyte y-cruncher output). Non-digit
bytes (whitespace, newlines, the `.` in `3.`) are ignored, and a leading `3`
is dropped so positions count from the first digit after the decimal point.
A warning is printed if the stream doesn't begin `3.14159265`.

All requested n are searched in **one pass**, and the tool stops reading the
file as soon as every n-digit prime has been found for every n — e.g. with
`-n 7` it reads only ~135.3M digits no matter how big the file is. Results
are printed to stderr as each n completes (with a progress line every 10^9
digits scanned), and a summary table goes to stdout at the end. If the file
ends before every prime has appeared, the affected n report how many primes
remain unseen instead of guessing.

Digits needed per n (set by A076130(n)):

| n | digits needed  |
|---|----------------|
| 5 | ~0.72M         |
| 6 | ~11.2M         |
| 7 | ~135.3M        |
| 8 | > 10^9 (a(8) was still unknown to OEIS as of 2021) |

Example:

```
./a076106 -n 8 pi-10t.txt
```

Memory use is dominated by bitmaps totaling ~2.8*10^maxN bits (prime sieve
plus seen-flags for every n): ~28 MB for maxN=8, ~2.8 GB for maxN=10,
~28 GB for maxN=11, ~280 GB for maxN=12 (a large-memory server). The tool
prints the estimate at startup. The sieve is segmented and multithreaded:
sieving to 10^12 takes on the order of minutes on a many-core server.

Expected digit requirements (coupon-collector estimates beyond n=8):
A076130(11) is likely around 2-4 trillion digits and A076130(12) around
25-40 trillion. maxN=13 is out of reach on both axes (a ~1.25 TB sieve
plus a ~1.25 TB seen bitmap, and an estimated ~3*10^14 digits), hence the
cap at 12.

## How it works

1. A segmented, multithreaded bitmap sieve of Eratosthenes marks all primes
   up to 10^maxN; per-length prime counts come from a popcount pass.
2. The file is read in 4 MB chunks. For every digit position, each still-active
   n maintains a rolling n-digit window value, updated with a small ring
   buffer of recent digits (no division in the inner loop); windows with a
   leading zero are skipped.
3. Cheap composite filters (last digit must be 1/3/7/9, value coprime to 3
   and 7) reject ~73% of windows before the big bitmaps are touched. The
   survivors go through a 1024-entry candidate queue whose sieve/seen bytes
   are software-prefetched, hiding DRAM latency when the bitmaps are far
   larger than cache (they are 125 GB each at maxN=12).
4. The first time a window equals an unseen n-digit prime, that prime is
   marked along with its position. The last prime to be marked is A076106(n)
   and its position is A076130(n); each n drops out of the scan once all of
   its primes have appeared, and the scan ends when none remain.
