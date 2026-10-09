# A197123 tools

C tools for finding a(n) of [OEIS A197123](https://oeis.org/A197123) (the first
n-digit substring to repeat in the decimal expansion of pi) and for exploring pi
generally:

* **pibloom** - pass 1: bloom-filter candidate finder (single-threaded)
* **pibloomt** - pass 1 with 10 worker threads sharing one bloom filter
* **pisearch** - pass 2: exact verification of the candidates
* **pifind** - find given numbers in pi, with position + context
* **piself** - find self-locating numbers (14 OEIS sequences in one pass)
* **pimatchscan / pimatchsort / pimatchverify** - the minimizer + partition + sort
  repeat finder: one read of pi, scratch on local disks, exact results for every
  length at once (see *Minimizer/sort repeat finder* below)
* **ycd library** - `ycd.h` + `ycd.c`, a standalone C99/POSIX reader and writer
  for y-cruncher `.ycd` files and directories of them, with no dependency on
  the pi tools, for use in other projects; `ycdtool` (info / dump / totext /
  check / encode) and `ycd_test` sit on top of it.  The format and the API are
  described in [ycd.md](ycd.md).

They replace `seqbuildbloom.c` / `sequsebloom.c` (kept for reference).

## Digit sources (all tools)

`-p` (or piself's file argument) takes a comma-separated list of files and/or
directories, read as one continuous digit stream.  Two formats are auto-detected
per file:

* **Text** - `3.14159...`; the leading `3.` is skipped and newlines/spaces are
  ignored.  Continuation files must be raw digits (a second file starting with
  `3.` is rejected).  Files are read in the order given.
* **y-cruncher `.ycd`** (`#Compressed Digit File`, v1.0.0/1.1.0, base 10) - read
  natively, no conversion needed.  Each file is block `BlockID` of `Blocksize`
  digits; a directory contributes every `*.ycd` inside it.  Files are ordered by
  `BlockID` (command-line order does not matter), must be contiguous and share a
  `Blocksize`, and a set may start at any block - positions are then absolute
  (block *b* starts at position *b*·Blocksize+1), so a single 100B block from the
  middle of a set scans correctly on its own.  Digits are bounded by the header,
  never by file size, so the padded tail of a block never leaks in.  When block 0
  is present the first word is checked against the known digits of pi (and a
  stream that carries the integer `3` at offset 0 is detected); otherwise the
  offset-0-after-the-point convention is assumed and the banner says so.

  ycd is **0.42 bytes/digit** (19 digits per little-endian 64-bit word) versus 1
  byte/digit as text - 2.4x less I/O per pass, which matters more than anything
  else on a network-mounted 20T-digit set.  Words are decoded with an
  ILP-friendly split (two constant divisions, then two-digit table lookups)
  and, on machines with 8+ CPUs, in parallel by a small pool (3 helper threads +
  the reader; override with `PI_DECODE_THREADS=N`), so the reader stage delivers
  ~900M digits/s and is never the pass's bottleneck.

Positions are 1-based from the first digit after the decimal point (OEIS
convention) in every tool and both formats.  Size arguments take suffixes
(`512M`, `64G`; digit counts `10B`, `1T`, `all`); status lines update in place on a
terminal and log one line per interval when redirected; Ctrl-C stops cleanly with
statistics and flushed output.  Every tool prints a `source format` line
describing what it opened (format, files/blocks, position range).

## Pass 1: pibloom (bloom-filter candidate finder)

    ./pibloom -p pi-10b.txt -l 19 -b 8G -f 1 [-B] [-n 10B] [-o candidates.txt]

| option | meaning |
|---|---|
| `-p SRC` | text file(s) or `.ycd` file(s)/directories, comma-separated (see *Digit sources*) |
| `-l N` | sequence length, 1..38 |
| `-b SIZE` | bloom size in bytes, `1M`..`16T` (1024-based: `512M`, `64G`, `3.5T`) |
| `-n COUNT` | digits to scan: `500M`, `10B`, `1T`, `all` (default all) |
| `-f SPEC` | only test windows whose 1-3 digit prefix is SPEC: a value (`3`, `24`), an inclusive range (`0-4`, `10-15`, `05-09`, `334-666`), or **`K/N`** = the K-th (1-based) of N equal partitions: `-f 2/20` = `[05-09]`, `-f 1/10` = `[0]`, `-f 2/3` = `[333-665]`. K/N resolves to the coarsest prefix length that divides N exactly, else to 3 digits (equal to within 0.1%). The banner shows the resolved range and its percentage |
| `-B` | blocked bloom: 4 bits in one cache line. ~1.7x faster, ~2x more spurious candidates |
| `-o FILE` | candidate file (default `candidates_L<len>_N<count>[_F<d>].txt`) |
| `-P SECS` | progress interval (default 30) |
| `-c SIZE` | read chunk size (default 64M) |
| `-N` | don't pre-fault the bloom memory |
| `-T` (pibloomt) | dry run: no bloom, nothing inserted; measures source + pipeline throughput |

Before scanning it prints the expected fill and false-positive rate for the chosen size so a
too-small table is obvious up front (warns above 50% fill, refuses > physical RAM).
Rule of thumb: bits per insertion ~ 20-70. With `-f`, insertions = digits/10, e.g. 1T digits with
`-f` = 100G insertions → `-b 512G` gives 41 bits/insert, fill ~9%, FP ~0.007% (~7M candidates).

Output lines are `<digits>:<position>` (position = 1-based index after the decimal point,
the same convention as OEIS and the old tools); `#` lines are metadata.

On completion it reports windows scanned/tested, candidates, rate, peak RSS, bits set
(tracked exactly during the run and re-verified with a full popcount), fill vs. theoretical fill
(ratio ≈ 1.0 means the hashes are behaving), and the false-positive rate.

### Hashing
The window is encoded injectively as a 128-bit key (`lo` = last 19 digits as an integer,
`hi` = the preceding digits), kept up to date with O(1) rolling arithmetic, then mixed with
splitmix64 into two 64-bit hashes and expanded to 4 bit positions by double hashing with a
division-free range reduction (works for any table size, not just powers of two). The same
code handles every length 1..38 with no per-length hash functions; measured fill matches
theory to 4 decimals. Bit lookups are software-prefetched 32 windows ahead so the
random DRAM accesses overlap instead of serialising.

## Pass 1, multithreaded: pibloomt

    ./pibloomt -p pi-10b.txt -l 20 -b 16G [-f 8] [-B] [-t 30] ...

Same options and output format as `pibloom` (plus `-t N`, worker threads, default
10, up to 256), but 1 reader thread + N worker threads sharing one bloom filter with
atomic (relaxed `fetch_or`) bit test-and-set.

Work is partitioned by digit *value*, never by file region: the three digits right
after the `-f` prefix form a number v in 0..999 and worker w owns the windows with
v in [w·1000/N, (w+1)·1000/N).  With `-t 10` that is exactly "split by the next
digit", with `-t 100` the next two; any N balances to within 0.1%.  Each worker
scans every chunk with a 16-lane SIMD prefilter (per-byte range tests on the prefix
digits and the leading digits of v - a superset test) and runs the exact check only
on lanes that pass, so scanning costs tens of ms per 64 MB chunk per worker and the
file is read and decoded once.  `-l` must be at least the prefix length + 3.

A campaign of N equal partitions is therefore `-f 1/N` ... `-f N/N` (each run gets
its own candidate file, `pisearch` auto-detects each range); the union covers every
window exactly once, with no range arithmetic to get wrong.

Why this is safe: both occurrences of any given window share the same prefix and
the same v, so they always land on the *same* worker, which sees them in file order
- a true repeat can never be missed to thread interleaving.  Only the few spurious
(false-positive) candidates can differ from a single-threaded run; `pisearch`
verifies exactly either way.

When to raise `-t`: a huge bloom is memory-latency-bound and throughput scales
with the number of independent miss chains in flight.  Measured on a 32-CPU box:
a DRAM-resident table went 45 -> 68 -> 85M tested/s at `-t 10/16/24` (near-linear);
a 3.42 TiB table in Optane Memory Mode went 18 -> ~26M/s from 10 to 24 workers and
then flattens, because the PMem DIMMs' internal concurrency is the ceiling (~26M
random line accesses/s).  On a machine with fewer cores than workers it only adds
contention.  Before the SIMD scanner (v1.0.0) the workers' `memchr` restarts on a
10%-dense digit capped a filtered pass at ~215M digits/s regardless of `-t`; that
is gone.

Validated: `-f 1` L=19 on 10B digits gives the identical a(19); unfiltered L=12
runs on 1B digits with `-t 1`, `10`, `16`, `30` and `100` each produced a verified
repeat set byte-identical to the single-threaded tool (499,299 repeats,
a(12)=756130190263 matches OEIS), with per-worker balance within 0.1%.  On this laptop it is
~2x pibloom (memory-bandwidth-bound); on big-RAM machines where the table dwarfs the
caches the latency hiding across 10 threads should scale further.

## Pass 2: pisearch (exact check)

    ./pisearch -p pi-10b.txt -l 19 -c candidates.txt [-n 10B] [-o results.txt] [-1]

| option | meaning |
|---|---|
| `-p SRC` | text or `.ycd` source(s), comma-separated (see *Digit sources*) |
| `-l N` | sequence length (must match the candidate file) |
| `-c FILES` | candidate file(s) from pibloom/pibloomt (`digits:pos` lines, `#` comments); repeat `-c` or give a comma list to verify several partitions in **one** pass over the digits (the table is the union, the auto filter the union range) |
| `-n COUNT` | digits to scan (default all) |
| `-f X` | test only windows whose 1-3 digit prefix is a value, range or `K/N` (`3`, `24`, `0-4`, `10-15`, `2/20`); `auto` (default) uses the tightest range covering every candidate (up to 3 digits); `none` disables |
| `-o FILE` | results file (default `results_L<len>.txt`) |
| `-1` | stop at the first repeat found |
| `-P SECS` | progress interval (default 30) |
| `-C SIZE` | read chunk size (default 64M) |

Startup shows progress for each phase (counting the candidate file, pre-faulting
the table, loading); the scan status line overwrites itself on a terminal.

Loads the candidates into an open-addressing hash table (dedups, validates every line),
rescans pi, records the first occurrence of each candidate and reports every second
occurrence. The filter prefix range is detected automatically from the candidate file. A small
bitmap prefilter rejects >99.8% of windows before touching the table. `-1` stops at the
first repeat found.

Results file: `<digits>:<first>:<second>`. The summary names the repeat with the
**earliest second occurrence** (that is a(n)) and, if several, lists them all and the
numerically smallest. Leading zeros are preserved and flagged.

Both programs: Ctrl-C prints the statistics gathered so far and flushes output (a second
Ctrl-C hard-exits); all file/allocation/parse errors are fatal with a clear message; a
pisearch run that saw fewer digits than pibloom warns about never-seen candidates.

## Minimizer/sort repeat finder: pimatchscan, pimatchsort, pimatchverify

An alternative to the bloom passes that reads pi **once** and moves the work to
scratch disks.  For every window of w consecutive k-mers, `pimatchscan` keeps the
k-mer with the smallest random rank (the *minimizer*) and writes a 14-byte record
(k-mer value, position) into one of P partition files chosen by a hash of the
value.  Two identical strings of length >= **L = k + w − 1** contain the same
window of k-mers, hence the same minimizer, hence equal keys in the same
partition - so every repeat of at least L digits is found; shorter ones only by
luck.  With the default k = 19 the key is the k-mer itself (exact, no hash
collisions), and about 2/(w+1) of positions produce a record.

    ./pimatchscan   -p /data/pi_ycd -d /scratch/run26 -L 26 -P 128 -t 24
    ./pimatchsort   -d /scratch/run26 -m 900G -t 24
    ./pimatchverify -d /scratch/run26

**pimatchscan** (pass 1): `-p` source (text or ycd, as everywhere), `-d` scratch
dir or comma-separated list of dirs (see below), `-L` guaranteed length (`-k`
1..19 [19], `-w` derived), `-n` digits, `-P` partitions (size them to <= 30% of
the sorting machine's RAM budget, i.e. a third of `-m`; small values make small
files for tests), `-H K/R` hash-range pass (see below), `-t` workers, `-B` write block per partition [4M] - writes
are buffered and never synced except at checkpoints - `-C` checkpoint seconds
[600], `-R` resume, `-f` overwrite, `-T` dry run (count records, write nothing),
`-I` progress, `-c` chunk.  The banner states the guarantee, the sampling
fraction, expected records and scratch bytes per partition, the RAM pimatchsort
will need, expected false pairs, the scratch write rate per 100M digits/s of
source, and free space (refuses to start if the scratch will not fit).
Checkpoints record the resume position and every partition's length; `-R`
truncates and continues (duplicate records around chunk/resume boundaries are
dropped by the sort).  Pipeline: reader thread -> chunk ring -> workers (rolling
k-mer value, random rank, monotonic-deque sliding minimum, O(1)/digit).

**Multi-directory scratch and hash-range passes.**  `-d /a,/b,/c` spreads the
partition files over several directories (different pools, mounts or NFS
exports) in proportion to their free space at scan start, with a 5% margin; the
first directory is the *primary* and holds the manifest, checkpoint, candidate
files and results, and the manifest records which partitions live where, so
`pimatchsort -d /a` and `pimatchverify -d /a` find everything (the full list is
accepted there too and ignored beyond the first entry).  The dry run `-T` prints
the layout.  When even the whole set cannot hold the records, `-H K/R` runs the
scan as R independent hash-range passes: the partition hash is spread over R x P
global partitions and pass K writes only its P of them, so each pass produces
1/R of the records (and of the false pairs) and every repeat lands wholly inside
one pass (both occurrences share the key, hence the hash).  Sort and verify each
pass as usual (pass K's files are identical in content to partitions (K-1)P ..
KP-1 of a single R x P scan); the answer is the earliest second occurrence over
the R matches.txt files.  `pimatchsort -D` deletes each partition file once its
candidates are written, so the next pass can reuse the same space; give every
pass its own primary directory (e.g. `/scratch/a27/p1`, `/scratch/a27/p2`).

**pimatchsort** (pass 2): loads each partition into RAM, parallel LSD radix
sorts the records by key (`-t` threads, `-m` RAM budget [half of physical] for
three buffers: sort source, sort destination and the next partition being
loaded; a partition larger than a third of the budget is sorted in sub-passes
that re-read the file keeping one hash slice each), and writes every key seen
at >= 2 distinct positions as a candidate group `key pos1 pos2 ...` to
`cand_NNNN.txt` with a `.done` marker (`-R` skips finished partitions, also
those whose file `-D` already deleted); `candidates.txt` is the concatenation.  A loader thread reads the next unit
(partition or sub-pass) straight from disk into the free buffer while the main
thread sorts and emits the current one, so disk and sort time overlap; the
per-partition status line shows both (`load Ns sort Ns`) and the MB/s figure
is the combined throughput of this run.

**pimatchverify** (pass 3): reads the groups, collects the digits around every
candidate position (`-E` context either side [100]) in one streaming pass over pi
- or by seeking with `-S` - and extends every pair to its full common length.
Seek mode reads only each snippet's bytes with `pread()` from `-j` parallel
threads [32] (no reopen, readahead or prefetcher per seek), so a disk array
delivers its aggregate random-read IOPS: use it whenever the candidate set is
much smaller than the source (millions of positions against terabytes); raise
`-j` to 64-128 on a wide array or NFS, and fall back to streaming (one
sequential pass at the source's read rate) only when the positions are dense.  Repeats of at least `-L` [the scan's L] are written
as `len pos1 pos2 digits` sorted by second occurrence, and the summary lists the
earliest second occurrence for every length the run covers: those are the a(n)
candidates, e.g. one L=14 run over 1B digits reports a(14)..a(17).  Pairs whose
match is shorter than L are the "false pairs" the banner predicted.

### Scratch pool and calibration (house1)

The scratch pool is 12 x 18 TB Exos as **2 x raidz1** on ZFS (`ashift=12`,
`recordsize=1M`, `compression=off`, `atime=off`, `logbias=throughput`, ARC capped
at 64 GB), ~160 TB usable, ~2.5 GB/s streaming.  A SLOG, special vdev or L2ARC
buys nothing for this workload (large sequential appends, no sync writes, every
byte read once); the one place SSDs would matter is holding the ycd *source*.
Scratch is regenerable, so raidz1 is the balance between a stripe (one dead
disk = restart) and raidz2.

Calibration on 200B digits of the local text source (`-P 512 -t 16 -B 16M`):

| step | result |
|---|---|
| pimatchscan | 282M digits/s (source read 283 MB/s = the limiter), 878 MB/s of records written (1.09 GB/s on the pool with parity), records exactly 22.22% = 2/9 of positions, 579 GiB in 512 partitions of 1.13 GiB each, RSS 27 GiB |
| pimatchsort | 17:53 for 579 GiB (579 MB/s with v1.3, where load and sort did not overlap; v1.4 overlaps them), 301 candidate pairs, 2,323 boundary duplicates dropped |
| pimatchverify -S | 14 s; all 301 pairs shorter than 26 (expected: 200B digits has a 0.02% chance of a 26-digit repeat) |

Projection for the full 21.1T digits (88% confidence for a(26)): pass 1 ~21 h,
66 TB of scratch (512 x ~128 GB; `-m 400G` on the 500 GB machine = 3 x 133 GiB
buffers, enough for one pass per partition), pass 2 ~32 h at the v1.3 sort rate
(v1.4's overlapped loading brings it toward the slower of the pool read rate and
the sort rate), pass 3 minutes.  The a(26) run:

    zfs create scratch/a26 && chown jbs:jbs /scratch/a26
    nohup ./pimatchscan -p /bigoutput/pi-20t.txt -d /scratch/a26 -L 26 -P 512 -t 16 -B 16M -C 900 > scan26.log 2>&1 &
    nohup ./pimatchsort -d /scratch/a26 -m 400G -t 16 > sort26.log 2>&1 &
    ./pimatchverify -d /scratch/a26 -S -j 64

Both long passes resume with `-R` after an interruption.  The verify summary's
`n=26` line (and `n=27`... for longer repeats found along the way) is the answer.

The a(27) plan (68T digits of the Google 100T ycd set = 680 blocks, ~29 TB, 90%
confidence): 68T x 2/9 x 14 B = ~212 TB of records, more than the 160 TB pool, so
two hash-range passes of ~106 TB each (P=1024 -> ~104 GiB partitions, sorted in
one pass under `-m 400G`), each in its own primary directory on the pool, the
first sorted with `-D` before the second scan starts:

    zfs create scratch/a27 && mkdir /scratch/a27/p1 /scratch/a27/p2 && chown -R jbs:jbs /scratch/a27
    nohup ./pimatchscan -p /scratch/pi100t -n 68T -d /scratch/a27/p1 -L 27 -P 1024 -t 16 -B 16M -C 900 -H 1/2 > scan27_1.log 2>&1 &
    nohup ./pimatchsort -d /scratch/a27/p1 -m 400G -t 16 -D > sort27_1.log 2>&1 &
    ./pimatchverify -d /scratch/a27/p1 -S -j 64
    nohup ./pimatchscan -p /scratch/pi100t -n 68T -d /scratch/a27/p2 -L 27 -P 1024 -t 16 -B 16M -C 900 -H 2/2 > scan27_2.log 2>&1 &
    nohup ./pimatchsort -d /scratch/a27/p2 -m 400G -t 16 -D > sort27_2.log 2>&1 &
    ./pimatchverify -d /scratch/a27/p2 -S -j 64

a(27) is the earlier of the two `n=27` lines.  With a second pool or export
mounted, `-d /scratch/a27/p1,/mnt/other` spreads one pass over both; if the set
holds ~225 TB, a single pass without `-H` does the whole job.

Validation: L=14 (k=12) over 1B digits reproduces a(14)-a(17) exactly, and every
one of the 5,282 distinct 14-digit repeats found by pibloom+pisearch is covered by
a maximal match (the window counts agree exactly); ycd and text sources give
identical output; an interrupted + resumed scan, a sub-pass sort, a resumed sort
and seek-mode verify all give identical results.  Throughput on a laptop: scan
323M digits/s writing 2.3 GB/s of records; the whole 1B-digit run takes ~10 s.
Scaling (90% confidence, sampled 18-22%): a(26) ≈ 55-66 TB of scratch, a(27) ≈
170-200 TB; scratch write bandwidth, not the source, is the pass-1 limiter.

## Utility: pifind

    ./pifind -p pi-10b.txt 1350168131352524443 84756845106452435773
    ./pifind -p /data/pi_ycd_dir -s 100B 5528194976825473781      # ycd set, seek to 100B
    Match 1 at position 4,750,253,204:  ...741<u>1350168131352524443</u>181...
    Match 2 at position 8,858,170,606:  ...593<u>1350168131352524443</u>778...

Prints every occurrence of one or more digit strings (1..30 digits each, up to 64
numbers, leading zeros allowed) - all searched in a single pass over the file - with
3 digits of context each side; the match itself is underlined (ANSI) on a terminal,
`[bracketed]` when piped.  `^` / `$` mark the start / end of the file when there is
less than 3 digits of context.  `-n` limits the scan, `-s` skips to an absolute position before
searching (`-s 1T` seeks instantly on `.ycd` and pure-digit text, across file/block
boundaries, and falls back to a counting skip if a text file contains newlines;
reported positions stay absolute), `-m` caps how many matches are
printed per number (default 1000; counting always continues), `-c` sets the read
chunk size.  With multiple numbers
each match line is labeled `Match N of <number>`, and matches are grouped per number
within each 64MB read chunk rather than globally position-sorted.  Runs at ~700M digits/s (memmem).
A live status line (%% done, digits searched, rate, ETA) ticks on stderr when run in a
terminal.  Exit status 0 = every number found at least once, 1 otherwise (missing
numbers are flagged NOT FOUND in the summary).

## Utility: piself (self-locating numbers)

    ./piself pi-10b.txt [-n 1B] [-o hits.txt] [flags]

The source must start at position 1 (a text file from `3.`, or a ycd set that
includes block 0): every self-locating check is relative to the start of pi, so a
partial block set is refused rather than producing meaningless hits.

Finds every number in pi that describes its own location, in one streaming pass.
The literature counts positions several ways; each convention is a separate OEIS
sequence, and every hit is printed (and written to the `-o` file) with its label
and OEIS number:

    HIT (k-3, A153221): 51 begins at position 48
    HIT (pithy-3, A109514): 9696 is the 9,696th 4-tuple

**Sliding-window family** - the digit window starting at 1-based after-decimal
position p spells p+d:

| flag | label | rule | OEIS |
|---|---|---|---|
| `-0` (= `-d -1`) | `0-based` | k at position k, first decimal = position 0 | [A057680](https://oeis.org/A057680) |
| `-1` (= `-d 0`) | `1-based` | k at position k, first decimal = position 1 | [A064810](https://oeis.org/A064810) |
| `-3` (= `-d 1`) | `3-first` | k at digit k, the leading 3 = digit 1 | [A057679](https://oeis.org/A057679) |
| `-d 2` | `k-2` | k at position k-2 (the 3 and the point count) | [A153220](https://oeis.org/A153220) |
| `-d 3` | `k-3` | k at after-decimal position k-3 | [A153221](https://oeis.org/A153221) |
| `-d 4` | `k-4` | k at after-decimal position k-4 | [A153223](https://oeis.org/A153223) |
| `-d 5` | `k-5` | k at after-decimal position k-5 | [A153224](https://oeis.org/A153224) |

`-d` takes any offset (negative allowed, comma lists OK), so `-d 6,7` searches
conventions nobody has cataloged.

**Multiplier family** (`-x N`, comma lists OK) - k found at position N*k:
`-x 2` = [A153227](https://oeis.org/A153227), `-x 3` = [A153228](https://oeis.org/A153228).

**Reversed family** (`-r N`, negative and comma lists OK) - the digit-REVERSAL of
k found at position k-N ("0161" starting at position 1610): `-r 0` =
[A366831](https://oeis.org/A366831) (after decimal), `-r 1` =
[A366830](https://oeis.org/A366830) (the 3 is digit 1). Internally the rolling
counter is simply kept least-significant-digit first, so the reversed compare
costs the same as a forward one.

**Pithy numbers** (`--pithy`) - chop pi into consecutive non-overlapping m-tuples;
an m-digit k is pithy when the k-th m-tuple is exactly k. Both conventions are
checked: tuples starting after the decimal point (`pithy`,
[A109513](https://oeis.org/A109513)) and starting at the 3 (`pithy-3`,
[A109514](https://oeis.org/A109514)).

**Residue family** (`-m N`, comma lists OK) - only k mod N has to appear at
position k: `-m 100` = [A153225](https://oeis.org/A153225), `-m 1000` =
[A153226](https://oeis.org/A153226). These are DENSE (~1.9% of all positions hit
for `-m 100`, ~0.28% for `-m 1000` - the sequences are effectively infinite), so
they are opt-in and excluded from the default; use with `-n` or a redirect.

With no selection flags everything except `-m` runs: offsets -1..5, reversed
offsets 0 and 1, multipliers 2 and 3, and both pithy variants - 14 sequences in
one pass, ~4 min for the 10B file. All checks share the single file read: each offset keeps a rolling
decimal counter (one compare per digit), multiplier/pithy checks fire only at
their sparse positions via countdowns / a next-event table, and `-m` uses a
precomputed residue-string table.

Each sparse convention expects ~0.9 hits per order of magnitude (an m-digit
match has probability 10^-m and there are 9*10^(m-1) candidates). The 10B scan
reproduces every published term of all 14 sequences and continues past them;
terms found here that OEIS does not yet list (verified against the raw digits):

| sequence | published up to | new terms found (position) |
|---|---|---|
| A153223 (`k-4`) | 783,609,697 | 4,927,138,295 (4,927,138,291) |
| A153224 (`k-5`) | 469,113,068 | 2,544,423,437 (2,544,423,432); 7,730,967,174 (7,730,967,169) |
| A153227 (`2k`) | 9,552,919 | 781,117,981 (1,562,235,962); 3,798,368,934 (7,596,737,868) |
| A153228 (`3k`) | 447,263 | 998,701,890 (2,996,105,670) |

Output file lines are `<label> <number> <OEIS>` (`-` when the convention has no
cataloged sequence); the file is flushed per hit, so an interrupted run keeps
everything found so far.

## Performance notes (measured on the campaign machines)

* **I/O is the budget.** Every pass reads the whole digit set; the tools process
  700-900M digits/s from cache, so a pass is normally bound by the source.  The
  readers issue kernel readahead (`fadvise WILLNEED`) and additionally run a
  prefetch thread that keeps 8 chunks (512 MB) of the file in the page cache ahead
  of the consumer; the consumer waits for it rather than reading the disk itself, so
  the source sees exactly one sequential stream (two streams at different offsets
  make a spinning-disk array seek itself to a crawl) while the large reads keep
  several NFS requests in flight.  Set `PI_NO_PREFETCH=1` to disable it.  For NFS beyond that:
  `read_ahead_kb` on the mount's bdi (live, no remount), and `rsize=1048576,
  nconnect=8` (remount) - a single NFS TCP stream tends to plateau near 200 MB/s.
  ycd sources need 2.4x fewer bytes than text for the same digits.
* **Pass 1 on a 3.42 TiB bloom in Optane Memory Mode** (1 TiB DRAM cache, 4 TiB
  PMem): ~26M tested/s with `pibloomt -B -t 24` (18M at `-t 10`) versus 45-85M/s
  for a table that fits the DRAM cache - a real penalty, not a cliff, and it is
  the PMem's random-access ceiling, so more workers do not help past ~24.  With
  1/10 partitions of 20T digits the pass is Optane-bound (~22h at ~250M digits/s);
  with 1/20 partitions it is network-bound instead (~10h at ~570M digits/s) and
  the bloom sits at 12% fill with ~200M candidates per partition, so pass 2 is
  4x lighter - the two plans cost about the same pass-1 time (20 x 10h vs
  10 x 22h) but 1/20 is the better shape, especially with multi-file pass 2.
  Pre-faulting the 3.42 TiB table takes ~6 min at ~9.5 GiB/s with the parallel
  first-touch (was 28 min single-threaded); the banner shows the thread count and
  rate.  Random bloom access is Optane's worst case, so benchmark before
  committing a campaign to it.
* **Sizing rule:** spurious candidates scale as fill^4 integrated over the run.
  700 GiB / 400B inserts (23% fill) gave ~366M candidates; 1 TiB gives ~90M,
  1.5 TiB ~20M, 3.5 TiB ~1.5M.  Prefer a longer `-f` prefix (or a range) over a
  fuller table: 1T inserts into 3.5 TiB (`-f 20-24` on 20T digits) is 12% fill
  and ~70M candidates.
* **Pass 2 memory:** 40 bytes per hash-table slot at <= 50% load (so 80-160
  bytes per candidate, power-of-two sized) plus a prefilter bitmap of 64 bits per
  candidate capped at 8 GiB.  65M candidates ~ 11 GiB; 366M ~ 44 GiB; 1.5G ~ 170
  GiB.  The scan is a three-stage prefetch pipeline (hash -> prefilter bit ->
  table probe); it doubled throughput on a laptop and should gain more on servers
  with higher memory latency, typically making pass 2 disk-bound again.
* **Pass 2 across partitions:** `pisearch -c a,b,c` verifies several partitions in
  one read of the digits.  Reads are shared but tests are not (three 1/20
  partitions = 15% of windows), and pisearch is single-threaded at ~50-70M
  tests/s, so 3-4 partitions per run stays close to I/O-bound; once pass 1 has
  freed the RAM, several such runs can share the page cache concurrently.
* **Shortcut for later partitions:** once a repeat is known, the remaining
  partitions only need scanning to its second-occurrence position (`-n`), since
  any repeat that beats it has both occurrences before that point.

## Versioning

All tools share one version, `PI_TOOLS_VERSION` in `pi_common.h`, printed in each
banner and `-h`, by `-V`, and stamped into candidate/result file headers (so an
output file records which build produced it).  It is bumped on every change:
patch for fixes/docs, minor for new options or behaviour, major for incompatible
output or file formats.

### Changelog

* **ycd library 1.0.0** (2026-09-19; the pi tools are unchanged at 1.6.0) -
  new standalone `ycd.h`/`ycd.c` with `ycdtool` and `ycd_test`, and `ycd.md`
  documenting the format.  Findings while building it: the bytes after a
  block's last word in a real file are valid-looking non-zero words, not
  padding, so digits must be bounded by the header; and a short final block
  with `TotalDigits: 0` is unrecoverable (the writer refuses to create one).
  `make test` gains the library self-test and a ycdtool <-> pifind cross-check.

* **1.6.0** - scratch sets and hash-range passes.  `pimatchscan -d a,b,c`
  spreads partitions over several directories in proportion to free space
  (layout recorded in the manifest; the first directory is the primary that
  holds manifest, checkpoint, candidates and results; `-T` shows the layout;
  old single-directory manifests still read).  `-H K/R` restricts a scan to
  hash-range pass K of R (each pass 1/R of the records, every repeat wholly
  inside one pass; pass K equals partitions (K-1)P..KP-1 of a single R x P
  scan).  `pimatchsort -D` deletes partition files once their candidates are
  written; `pimatchsort`/`pimatchverify -d` accept the list (first entry used).
  Verified: multi-directory and hash-pass candidate sets identical to a
  single-directory scan; interrupted multi-directory scan resumed to identical
  candidates; `make test` gains a two-directory hash-pass check.

* **1.5.0** - pimatchverify `-S` rewritten: `-j` parallel seek threads [32]
  issue small `pread()`s of just the snippet bytes (the old path reopened the
  file, issued two 4 MB readahead hints and woke the sequential prefetcher per
  position, which held the a(26) verify on house1 to ~73 snippets/s = the
  pool's sequential bandwidth); the status line shows snippets/s and an ETA.
  Fix: ycd block-0 sets that carry the leading "3" (lead3) reported every
  streamed position one too low (first_pos was 0); seek mode, pifind and
  piself on such a set now agree with the text source exactly.  Verified: seek
  and streaming results identical on text and on a synthetic two-block ycd.

* **1.4.0** - pimatchsort overlaps loading with sorting: a loader thread reads
  the next partition (or sub-pass) straight from disk into a third buffer while
  the main thread radix sorts and emits the current one (`-m` now covers three
  buffers, so a partition must fit in a third of the budget to sort in one
  pass); unfiltered loads no longer go through a staging copy; the status line
  shows per-partition load and sort seconds; MB/s and ETA after `-R` count only
  this run's bytes.  Output verified identical to 1.3.3 (single-pass, sub-pass,
  and interrupted+resumed runs).

* **1.3.3** - pimatchscan's status "records (% of digits)" now divides by windows
  the workers have finished rather than chunks the reader has published, so it
  reads the true density from the first tick instead of converging to it.

* **1.3.2** - pimatchscan's false-pair estimate scaled by the measured 3x
  (minimizer sampling favours low-rank k-mers, which collide with each other more
  than a uniform sample); README gains the house1 scratch-pool layout, calibration
  numbers and the a(26) run commands.

* **1.3.1** - gcc warning fixes in the pimatch tools (temp-file name buffers,
  redundant unsigned test); no behaviour change.

* **1.3.0** - new `pimatchscan` / `pimatchsort` / `pimatchverify`: the minimizer +
  partition + radix-sort repeat finder (one read of pi, scratch on local disks,
  exact results for all lengths >= L in one run, checkpoint/resume, dry run,
  sub-pass sorting, seek-mode verification).  Text-source seeks cache the
  pure-digit check (seek-mode verify was re-sampling the file per position).

* **1.2.1** - dry-run banner no longer prints bloom size / fill estimates for the
  placeholder table.

* **1.2.0** - `pibloomt -T` dry run (read, decode, dispatch and hash everything, but
  allocate no bloom and insert nothing) to measure a source + pipeline rate without
  the pre-fault or the bloom; status lines of both bloom tools now also show
  M digits/s and MB/s, since tested/s scales with the filter fraction and is not
  comparable between partition sizes.

* **1.1.2** - the reader's prefetch thread is now the only disk reader (the consumer
  waits for it and reads from page cache), so a source sees one sequential stream;
  two readers at different offsets made spinning-disk NFS exports crawl (30-40 MB/s
  where `dd` got 230).

* **1.1.1** - `pibloom` and `pibloomt` status lines overwrite in place on a terminal
  (one line per interval when redirected), like `pisearch` and `pifind`.

* **1.1.0** - `pisearch -c` accepts several candidate files (repeat `-c` or a comma
  list) and verifies all of them in one pass over the digits; the auto filter
  becomes the union range and the results file lists every input.

* **1.0.0** - first versioned release: text + y-cruncher `.ycd` sources (multi-file,
  directories, absolute positions), prefix filters as values / ranges / `K/N`
  fractions up to 3 digits, `pibloomt -t` worker count with SIMD window scanning,
  three-stage prefetch pipeline in `pisearch`, parallel pre-faulting, readahead +
  prefetch thread + parallel ycd decode in the reader, `pifind` multi-number / `-s`
  seek, `piself`, in-place status lines.

## Validation (this machine, Apple Silicon, single thread)

| test | result | time |
|---|---|---|
| a(7) on 10K digits | 8530614 at 4,167 / 4,601 ✓ | <1 s |
| a(10) on 1B digits | 4392366484 at 182,105 / 240,479 ✓ | 17 s + 2:18 (46.8M repeats written) |
| a(19) on 10B digits, `-f 1 -b 8G` | 1350168131352524443 at 4,750,253,204 / 8,858,170,606 ✓ | 26 s + 16 s |
| synthetic 30-digit repeat, L=25 (two-word key), classic and `-B` | all expected positions ✓ | |
| ycd: synthetic block set built from pi-1b.txt (1 block, and a 2-block directory set) | pibloom candidates byte-identical to text; pibloomt+pisearch verified repeats identical; seam window across blocks, cross-block seek, EOF marker ✓ | |
| **a(25) on 10T digits** (pibloom/pibloomt + pisearch, ten prefix partitions on mathh/house1) | 2420215682379102240706463 at 399,369,008,413 / 2,564,331,708,881 ✓ (pifind-confirmed) | days |
| **a(26) on 21.1T digits** (pimatchscan/sort/verify on house1, 2026-09-11) | 21381728990849390313328992 at 1,803,801,456,180 / 11,975,559,503,988; runner-up 95798890137711775569397580 at 3,905,852,113,873 / 15,807,462,740,009; 2 repeats vs 2.2 expected, 3,339,972 false pairs vs ~3.3M expected | 20:42 + 1d 05:27 + 1:20 |
| ycd: real 42GB y-cruncher block 1 (positions 100B+1..200B) | first word at 100,000,000,001; padded tail bounded exactly; 3 in-block 19-digit repeats found by pibloomt+pisearch and confirmed by pifind seeks ✓ | 725M digits/s decode |

Throughput: ~60-70M inserted windows/s (classic), ~115M/s (`-B`) when the table fits in RAM;
pisearch ~600M windows/s with the filter. For the 10T-digit file with `-f`, expect roughly
1T insertions → 4-5 h for pass 1 per filter digit if the table fits in RAM.
