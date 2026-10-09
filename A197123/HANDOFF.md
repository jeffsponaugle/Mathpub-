# HANDOFF - A197123 pi repeat search: state of the work

Read this first in a new session. It summarizes the problem, the machines, the
tools, the design decisions, what was measured, what went wrong, and what is next.
The README has the per-tool reference; this file has the *why* and the numbers.

## 1. The problem

OEIS A197123: a(n) is the first n-digit substring to repeat in the decimal
expansion of pi. "First" means the repeat whose **second occurrence** comes
earliest. Positions are **1-based from the first digit after the decimal point**
("1 appears at positions 1 and 3"); the leading 3 is never part of a window;
leading zeros are preserved (a(4) = 0582, a(18) = 013724950651727463).

Known terms: a(1)..a(24) are on OEIS. **a(25) = 2420215682379102240706463 at
positions 399,369,008,413 and 2,564,331,708,881** was found here (bloom pipeline
over 4T digits, all ten first-digit partitions scanned to the threshold position,
both positions spot-checked with pifind on the raw file). Runner-up:
2209805007189562893228030 at 1,607,917,458,042 / 2,903,162,347,272. OEIS
submission text (data, %C comment with positions, %E extension) was drafted in
conversation; Jeff (jeff@sponaugle.com) submits.

**a(26) = 21381728990849390313328992 at positions 1,803,801,456,180 and
11,975,559,503,988** was found 2026-09-11 with the pimatch pipeline over the
21.1T-digit text file on house1 (exact: every window >= 26 digits covered).
Runner-up: 95798890137711775569397580 at 3,905,852,113,873 / 15,807,462,740,009.
Two 26-digit repeats vs 2.2 expected; no 27-digit repeat (0.2 expected), so
a(27) needs the 68T ycd run. pifind confirmation of both positions pending at
the time of writing; OEIS entry drafted (same shape as a(25)'s).

Birthday model: P(repeat of length n by N digits) = 1 - exp(-N^2 / (2*10^n)); the
median first repeat is ~1.18*10^(n/2) digits. a(26): 86% within 20T, 88% within
21.1T, median 11.8T. a(27): 18% in 20T, ~68T for 90%. a(28)+ are beyond any
published digit set (215T for 90%).

## 2. Machines and data

| machine | role | notes |
|---|---|---|
| laptop (this repo, macOS/clang, Apple Silicon, 36 GB) | development | `pi-1b.txt`, `pi-10b.txt` (text), `pi.ycd` (y-cruncher block 1 = positions 100B+1..200B, 42 GB) |
| mathstore (10.20.30.33) | file server, old 8-core Xeon | `/databig/pi/google/` = ycd set, 312 blocks x 100B digits = 31.2T digits (`dd` 592 MB/s); NFS export `databig` |
| house1 (10.20.30.29) | file server, ~500 GB RAM | `/bigoutput/pi-20t.txt` = 21.1T digits text (2 x raidz2 of 1.8 TB 10K SAS, local read ~283 MB/s); **`/scratch`** = new 2 x raidz1 of 12 x 18 TB Exos, ~160 TB, for pimatch scratch; `bigwork` was destroyed to build it |
| mathh (10.20.30.27) | big-memory compute, 32 CPUs | 1 TiB DRAM + 4 TiB Optane in **Memory Mode** (3.9 TiB visible, DRAM is a direct-mapped cache); NFS client with `nconnect=8,rsize=1M` |
| mathc | older compute | ran pass-2 verifications |

NFS from mathh: `/mnt/databig` (ycd), `/mnt/bigoutput` and `/mnt/bigwork` (text).

## 3. The tools (all C, `make`, `make test`, one shared version `PI_TOOLS_VERSION`)

Bloom pipeline (RAM-bound, many passes):
* **pibloom** - pass 1 single-threaded: every window into a k=4 bloom; windows whose 4 bits are already set are candidates (`digits:position`). `-B` blocked mode (1 cache line/window). Any table size 1M..16T.
* **pibloomt** - pass 1 with a reader thread + `-t N` workers sharing one bloom (atomic bit ops). Workers own windows by the 3 digits after the `-f` prefix (value v in 0..999, worker w owns [w*1000/N,(w+1)*1000/N)), scanning with a 16-lane SIMD range prefilter. **Invariant: both occurrences of a repeat share prefix and v, so they hit the same worker in file order - no false negatives from threading.** `-T` dry run.
* **pisearch** - pass 2: candidates -> open-addressing table + prefilter bitmap, three-stage FIFO prefetch pipeline (hash -> bit -> probe); several `-c` files verified in one pass; reports every repeat and the earliest second occurrence.
* `-f` everywhere: 1-3 digit prefix value/range or `K/N` fraction (`2/20` = `[05-09]`); pisearch auto-detects the tightest range from its candidates.

Minimizer/sort pipeline (I/O-bound, ONE read of pi + scratch; exact; all lengths >= L at once):
* **ycd library** (`ycd.h`, `ycd.c`, `ycdtool.c`, `ycd_test.c`, `ycd.md`; v1.0.0, 2026-09-19) - standalone reader/writer for y-cruncher files for Jeff's other projects; independent of pi_common.h (the pi tools keep their own reader with prefetch + decode pool). ycd.md is the format spec with a verified/unverified table. Open item: confirm the partial-last-word convention at a real block seam on house1 (`ycdtool dump DIR 199999999987 14` vs the text file).
* **pimatchscan** - minimizers over windows of w k-mers (L = k+w-1 guaranteed), 14-byte records (8-byte exact k-mer value for k<=19, 6-byte position) into P partition files by `fastrange(mix64(value+C), P)`; block-buffered writers; checkpoint/resume; `-T` dry run; `-d a,b,c` scratch set split by free space (primary = first dir, layout in the manifest); `-H K/R` hash-range pass (hash over R x P global partitions, pass K keeps its P).
* **pimatchsort** - per partition: load, parallel LSD radix sort, emit `key pos1 pos2...` groups; loader thread + third buffer overlap the next load with the sort (`-m` covers 3 buffers); sub-passes by low bits of the same hash if a partition > a third of the RAM budget; `.done` markers, `-R`; `-D` deletes each partition file after its candidates.
* **pimatchverify** - snippets around every candidate position (streaming, or `-S` direct preads from `-j` threads [32] - ~1,350/s on house1's pool), extend each pair to its maximal length, dedupe, report by second occurrence + the earliest per length (= a(n) candidates); with `-H` passes, take the earliest across the passes' matches.txt.

Utilities: **pifind** (multi-number search with context, `-s` seek), **piself** (self-locating numbers, 14 OEIS sequences; needs a source starting at position 1).

Shared code: `pi_common.h` (digit sources, hashing, window keys, filters, timing, formatting, huge allocations, signals, version) and `pimatch.h` (records, manifest, checkpoint).

## 4. Design facts worth knowing before touching the code

* **pistream** (in pi_common.h) reads text or y-cruncher `.ycd`, auto-detected per file; `-p` takes comma lists and directories (all `*.ycd`, ordered by BlockID, contiguity checked; a set may start at any block -> absolute positions). ycd: text header ending in `EndHeader` + NUL, then little-endian u64 words of 19 digits (MSD first), `Blocksize`/`BlockID` bound the digits (files are padded); offset 0 = first digit after the point (verified against text at position 100,000,000,001). Decoder: 8/8/3-digit split + pair table; a decode pool of 3 helper threads.
* **Reader I/O**: `fadvise(WILLNEED)` readahead plus a prefetch thread that is the ONLY disk reader (the consumer waits for it) so the source sees one sequential stream - two readers at different offsets made a spinning-disk NFS export crawl (30-40 MB/s vs 230). `PI_NO_PREFETCH=1` disables it. The consumer's EOF wait is capped at the file end (that was a deadlock once).
* **pireader** hands out chunks with an L-1 digit overlap so windows tile exactly; `first_pos`/`digits_est` feed banners and % done; never roll the window key past the last window (reads d[L]).
* **Window key**: (hi, lo) = leading L-19 digits and last 19 digits as integers - injective, rolled in O(1); hashed with splitmix64 into two 64-bit hashes; bloom bits by double hashing + fastrange (any table size). Measured fill matches theory to 4 decimals.
* **pisearch memory**: 40-byte entries at <=50% load, power-of-two slots (80-160 B/candidate) + 64 bits/candidate prefilter capped at 8 GiB. 65M candidates ~ 11 GiB; 366M ~ 44 GiB; 1.5G ~ 170 GiB.
* **pimatch**: guarantee L = k+w-1; density ~2/(w+1) (measured exactly 2/9 at w=8); k=19 keys are exact so the only false pairs are equal 19-mers that do not extend (measured ~3x the uniform estimate N^2/(2*10^k), because minimizers favour low-rank k-mers); duplicate records only at chunk/resume boundaries, dropped by the sort. Checkpoint = resume position (start of the first not-fully-processed chunk) + every partition's byte length; the reader keeps in-flight chunks within nslots of `done_contig` so per-slot arrays are unambiguous.
* **Status lines** overwrite in place on a terminal (`isatty`), one line per interval when redirected. All tools print `M digits/s (MB/s)`; **tested/s scales with the filter fraction and is not comparable across partition sizes** (this caused confusion once).
* Every code change: bump `PI_TOOLS_VERSION`, add a README changelog line (also in memory).

## 5. Measured performance (the numbers decisions were made on)

| what | result |
|---|---|
| bloom pass 1, DRAM-resident table | 45 / 68 / 85M tested/s at `-t 10/16/24` (near-linear) |
| bloom pass 1, 3.42 TiB table in Optane Memory Mode | 18M (10 workers) -> ~26M tested/s (24); PMem random-access ceiling, more workers do not help |
| bloom sizing | 700 GiB / 400B inserts -> 23% fill, 366M candidates (measured); 1 TiB -> ~90M; 3.5 TiB / 1T inserts -> 12%, ~70M |
| pisearch | 8.1M tested/s unpipelined (DRAM latency) -> pipelined ~2x on laptop, disk-bound on servers |
| pre-fault | 32 threads: 3.42 TiB in ~6 min at ~9.5 GiB/s (single-thread was 28 min) |
| sources from mathh | databig ycd ~570M digits/s; bigoutput text 409 MB/s (our reader beats `dd`'s 230 thanks to 8 NFS connections); bigwork 121; databig text 170 |
| pimatchscan on house1 (local) | 282M digits/s = source-bound (283 MB/s local read); 878 MB/s of records; RSS 27 GiB with `-t 16 -P 512` |
| pimatchsort on house1 | 579-619 MB/s with v1.3 (load and sort serial) over the full a(26) set; v1.4 overlaps them with a loader thread + third buffer - not yet measured at scale (`load Ns sort Ns` in the status line shows which side limits) |
| a(26) via pimatch, actual | pass 1: 20:41:42 (283M digits/s, 881 MB/s), 4,688,890,321,625 records = 22.22%, 59.70 TiB in 512 x 119.4 GiB; pass 2 (v1.3): 579 MB/s, ~6,500 groups per partition = the 3x false-pair model (~3.3M total) |

## 6. Bugs found and lessons (do not repeat them)

* Several threads claimed the same chunk (counter read before the condvar wait) -> 3x records. Claim after the wait.
* Prefetch gating deadlocked at EOF (last read wants more than the file has). Cap the wait at total bytes.
* Two `fmt_u64()` calls sharing one buffer in a single printf print the same value twice - always use distinct buffers.
* gcc (not clang) flags: `%0*d` with unbounded args (bound them explicitly), `-Wmisleading-indentation` for `if (a) x; if (b) y;` on one line, `-Wmaybe-uninitialized` when a callee's early return skips outputs (check the return), temp-file buffers the same size as the path. Jeff builds with gcc on Linux; clang here is quieter.
* macOS has no `timeout`; use `perl -e 'alarm shift; exec @ARGV' N cmd`. Killing a hung background batch can let its remaining commands run later and race the next test.
* Do not compare 500M-record sets in Python (50 GB); compare via the tools' outputs (candidates/matches) instead.
* A tool's rate limit can hide behind another: the workers' `memchr` restarts capped filtered passes at ~215M digits/s regardless of reader speed until the SIMD scan replaced them; measure with `-T` dry runs and `pifind` (single-thread read rate) to isolate.
* Old binaries on servers: check `-V`.
* Random access through a streaming abstraction is a trap: pistream_seek reopened the file, issued 4 MB `WILLNEED` hints and woke the prefetcher per position, so 200-digit snippets cost ~4 MB of disk traffic each (73/s on a 290 MB/s pool). Seek workloads need small `pread()`s, no readahead, and enough threads in flight to use the array's IOPS.
* Raw partition files are not comparable between runs (record order depends on worker interleaving); compare candidates.txt bodies (or sorted unions across hash passes) instead.
* ycd block 0 with the leading "3": first_pos must be 1 (the 3 is dropped), not `start + 1 - lead3` = 0; the synthetic ycd test recipe must include the lead 3 to catch this.

## 7. Test recipes (fast, exact)

Known terms for quick checks (positions after the point): a(7) 8530614 @ 4,167/4,601; a(10) 4392366484 @ 182,105/240,479; a(12) 756130190263 @ 447,673/857,982; a(13) 3186120489507 @ 3,143,598/5,563,712; a(14) 18220874234996 @ 4,821,309/9,289,694; a(15) 276854551127715 @ 19,552,188/28,048,917; a(16) 8230687217052243 @ 67,889,102/129,440,088; a(17) 93415455347042966 @ 109,303,021/262,527,966; a(19) 1350168131352524443 @ 4,750,253,204/8,858,170,606 (10B file, `-f 1` or `-f 13`).

* `make test` covers a(7) (both bloom tools), pifind, piself, and pimatch a(10).
* Exactness of any pass-1 variant: run pibloom+pisearch unfiltered on pi-1b.txt L=12 (499,299 verified repeats), then the variant, and diff the sorted verified-repeat sets (or the exact prefix subset for a filtered run).
* pimatch vs bloom: the L=14 pisearch set (5,282 repeats) must each be covered by a maximal match; window counts agree exactly.
* ycd: build a synthetic block set from pi-1b.txt with a 40-line Python encoder (header with `\r\n`, `TotalDigits: 0`, 4 KiB padding) and diff against the text run; real `pi.ycd` checks: first word at 100,000,000,001, EOF marker at 199,999,999,994.
* CR status lines: run under `pty.spawn` in Python and count `\r`.

## 8. Where things stand / next steps

* **a(26) run on house1**: pass 1 finished 2026-09-10 (20:41:42); pass 2 finished 2026-09-11 with v1.3.3 (1d 05:27, 619 MB/s, 3,339,973 groups = the 3x false-pair model, largest group 3, RSS 239 GiB). Pass 3 with the v1.3/1.4 seek path ran at ~73 snippets/s (25 h for 6,679,947 positions) - replaced in 1.5.0 by parallel direct preads: `./pimatchverify -d /scratch/a26 -S -j 64` then ran at ~1,350 snippets/s (80 min). **Result (2026-09-11): a(26) = 21381728990849390313328992 at positions 1,803,801,456,180 and 11,975,559,503,988** (exactly 26 digits, not truncated; second occurrence at 12.0T is the 51st percentile of the birthday model). pifind confirmation of both positions was in progress; OEIS submission pending.
* Bloom partition runs for a(26) on mathh (`-f 1/10`, 3.5T bloom) were in progress before the pivot; partitions done there still count.
* **a(27) next**: 1.6.0 (2026-09-11) added everything the run needs - multi-directory scratch (`-d a,b`, split by free space, layout in the manifest), hash-range passes (`-H K/R`), `pimatchsort -D`, parallel seek verify (`-j`). Commands are in the README ("The a(27) plan"): two passes of ~106 TB on the 160 TB pool, P=1024, `-m 400G`, source = first 680 blocks of the Google 100T ycd set (to be downloaded to a `scratch/pi100t` dataset, ~29 TB). Estimate ~10 days end to end.
* Pending improvements, in value order: (1) DONE in 1.4.0 - pimatchsort overlaps loading with sorting (third buffer + loader thread); (2) several parallel `pread`s in the reader's prefetcher for local ZFS single-stream sources (~283 -> 400-500 MB/s); (3) multi-threaded pisearch (partition by prefix classes, as pibloomt); (4) a(27) needs ~68T digits and ~190 TB scratch - more spindles.
* A friend is building a similar sort-based approach; the analysis of its I/O (sampling makes scratch traffic ~ the bloom plan's source traffic; the win is moving bytes to local disks and exactness) is in the README's performance notes.
