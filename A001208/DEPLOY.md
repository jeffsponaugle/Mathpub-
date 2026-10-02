# Deploying the postage-stamp searches on other machines

**Status 2026-10-01 08:35:** deployed on `mathd` (jbs@10.1.30.23, Ubuntu 24.04, 2 x Xeon Platinum 8168 = 48 cores /
96 threads, 754 GB) in `~/A0001208/psph`. gcc 13 `make linux` (AVX2 path) built cleanly; both `psph` and `psph16`
pass the 422-check selftest there. Calibration: the (7,7) selftest case takes 7.2 s wall on 96 threads vs ~15 s on
the M1 Pro, so the whole box is only ~2x one M1 Pro for this code (per thread ~2.2x slower than an M1 core).
`jobs.sh` (nohup) runs job 1 (n(9,8) items 1200..11003, chunks of 1000, target 5522) and then job 2
(A001209(303), psph16, depth 2) back to back; logs in `~/A0001208/psph/logs/<tag>/` and `jobs.log`.
Revised ETAs on this box: job 1 ~1-1.5 days, job 2 ~4-7 days.

**08:45 — second box `mathg`** (jbs@10.1.30.26, 2 x Xeon Gold 5315Y = 16 cores / 32 threads at 3.2 GHz, 1.5 TB;
~0.4 x mathd) deployed in `~/A001208/psph`, both binaries pass the 422-check selftest. Job 1 was re-split:
mathd runs items [1200, 8000) (chunks of 1000) then job 2 (A001209(303)); mathg runs items [8000, 11003)
(chunks of 500) then job 3 (A001211(27), depth 4, target 176381, chunks of 100000). Each box has its own
`jobs.sh`/`jobs.log`. Rebalance later by changing HI/LO before a box reaches the other's items.

**09:25 — third box `orin1`** (jbs@10.1.30.40, Jetson Orin Nano Super: 6 x Cortex-A78AE 1.7 GHz, 7.5 GB, CUDA 13.2):
NEON build, smoke tests pass, selftest running; measured 5.5e7 candidates/s on n(9,8) = 6 % of mathd, 10 % of mathg.
Assigned the cheap tail of the item list, items [10500, 11003) (mathg's chunks 5 and 6 were marked .done so mathg
skips them), chunks of 100, starting automatically once its selftest passes (`jobs.sh` waits for it). Its main
value is as a CUDA development/test box for a GPU port of the psph leaf, which would be the way to use the two
DGX Sparks' GPUs on these searches.

Measured same-target throughput (candidates/s, n(9,8), TGT 5600): mathd 8.7e8, mathg 5.6e8, M1 Pro ~5e8, orin1 5.5e7.

**11:00 — fourth box `mathb`** (jbs@10.1.30.21, 2 x Xeon Platinum 8268 = 48 cores / 96 threads at 2.9 GHz, 754 GB; ~mathd
speed) in `~/A001208/psph`. Rebalance of the n(9,8) search (measured chunk rates at 10:53: mathd 5 h/1000 items,
mathg 5.5 h/500, orin1 6 h/100): mathb takes items [5200, 8000) from mathd (chunks 4-6 marked done there),
[10000, 10500) from mathg (its chunk 4) and [10700, 11003) from orin1 (its chunks 2-5), at target 5650 (best known 5649).
Expected finishes: mathd ~Oct 2 05:00 PDT, mathb ~Oct 2 06:00, mathg ~Oct 2 08:00, orin1 ~Oct 1 23:00.
mathb then chains job 4: A001210(91) = n(91,5), target 8897043, depth 3, 133730 items, chunks of 5000.

Everything needed is in `psph/` (single C file + Makefile + driver script). Copy that directory (or the
tarball `psph_deploy.tgz` in the project root) to each machine.

## Build and validate (once per machine)

```bash
cd psph
make linux            # Linux x86-64: gcc -march=native, AVX2 path; Linux aarch64 (Sparks): NEON path
# make linux LCC=clang  (to use clang);   make linux NOSIMD=1  (scalar fallback if the SIMD build misbehaves)
# macOS Apple Silicon: just  make
./psph -selftest quick     # ~4 min on 10 threads, 422 checks against the published tables; MUST print 0 failures
./psph16 -h 9 -k 5 -t 797 -j 4 -p 0   # smoke test of the 16-bit build: must print SOLUTION ... 1 9 23 108 181
```

The AVX2 path was written on the Mac and cannot be tested there, so run the selftest on the x86 box before
any production run. If it fails, rebuild with `NOSIMD=1` and run the selftest again.

## Jobs, in priority order

All runs use `run_range.sh TAG H K TGT D LO HI CHUNK THREADS [BINARY]`; the item range `[LO,HI)` is what
partitions a search across machines (item order is deterministic; the split depth D must be identical
everywhere). The final answer is the maximum `n_h` over all machines' `logs/TAG/solutions.log` once every
range is complete; every chunk is resumable (re-run the same command).

### 1. n(9,8) = A053348(9) = A005344(8)  — remainder of the overnight run

* Best known: **6082** (basis 1 8 27 88 197 521 1226 1461, mathg, items 9000..9250, 23:00); earlier 5945, 5844, 5822, 5809, 5781, 5728, 5708, 5674, 5649, 5599, 5521. mathd/mathg/orin1 run target 5600, mathb 5709 (targets only affect speed; the answer is the max over all solutions).
* Measured chunk times (13:40 PDT Oct 1): mathd 3.8 h per 1000 items (chunk 0 done, 0 solutions in items 1200-2200),
  mathb ~6.2 h per 1000, mathg 4.8 h per 500 (20 bases >= 5600 in items 8000-8500!), orin1 ~10 h per 100.
  Projected finishes: mathg ~Oct 2 04:30, orin1 ~05:30, mathd ~08:00, mathb ~11:00 PDT. When mathg/mathd finish,
  hand them mathb's remaining chunks (mark the chunk .done on mathb before it starts; run it elsewhere with a new tag).
* 16:36 PDT Oct 1: the Mac rejoined: ranges [10700,11003) and [10000,10500) released by mathb (logs/h9k8c, h9k8b marked done
  there) and [4200,5200) released by mathd (its chunk 3 marked done), target 5823, 10 threads (psph/mac_jobs.sh).
* 18:00 PDT Oct 1 rebalance #2 (items 8500-10000 cost ~2x the earlier mathg chunk; the Mac's third range was too big):
  mathg's chunk 3 [9500,10000) -> mathd (runs after its own range, before A001209(303)); the Mac's third range cut to
  [4200,4700) and [4700,5200) -> mathd as well (mathd wrapper `jobs2.sh`, Mac wrapper `mac_jobs2.sh`, target 5845).
  Projected: mathg ~Oct 2 07:00, mathd ~08:00-11:00, mathb ~08:30, Mac ~09:00, orin1 ~03:30 -> n(9,8) settled ~Oct 2 late morning.
* 22:40 PDT Oct 1 rebalance #3: mathg's chunk [9000,9500) projected ~20 h, so mathg now runs only [9000,9250) (then job 3)
  and mathd adds [9250,9500) (mathd wrapper `jobs3.sh`: own range -> [9500,10000) -> [9250,9500) -> [4700,5200) -> A001209(303)).
  Projected: mathb ~Oct 2 08:20 (then A001210(91)), orin1 ~08:20, mathg ~09:00, Mac ~11:00, mathd ~12:00 -> n(9,8) settled ~Oct 2 noon.
  Solutions >= target found so far: mathg 266, mathd 60, mathb 49, orin1 13, Mac 1; best 5945 (mathb).
* 02:00 PDT Oct 2 rebalance #4 (mathg finished [9000,9250) in 3 h, not 20): final queue —
  mathd: [9500,10000), [9250,9500) then A001209(303); mathb: chunk [7200,8000), then [10250,10500) (from the Mac), then
  A001210(91); mathg: [4700,5200) (from mathd), then [4450,4700) (from the Mac), then A001211(27); Mac: [10000,10250) then
  [4200,4450); orin1: [10600,10700). Projected n(9,8) completion ~Oct 2 09:30 PDT. Best 6082 (mathg).
* 04:05 PDT Oct 2 rebalance #5: mathb's chunk [7200,8000) is slower than planned (ends ~09:45), so [10250,10500) moved from mathb to mathd (after [9250,9500)); mathb goes straight to A001210(91). Projected n(9,8) completion ~09:45 PDT.
* 07:47 PDT Oct 2: mathd finished all its n(9,8) ranges and started A001209(303) (psph16, 301 items).
* 09:10 PDT Oct 2: Mac paused on request; its range [4200,4450) moved to mathd (jobs7.sh: [4200,4450) then A001209(303) restarted from scratch, ~45 min lost).
* 08:05 PDT Oct 2: mathg finished all its n(9,8) ranges and started A001211(27).
* 03:50 PDT Oct 2: orin1 finished its range [10500,10700) (29+ bases >= 5600, best 5782); orin1 now idle.
* CUDA port (orin1): VALIDATED — all completed differential cases identical (incl. (8,8) items 3000-3002, the -s2 0 cases with every
  counter matching, and int64+paranoid cases); only the two 15-minute-timeout cases are inconclusive and are being redone
  with small ranges (tests/supp.sh). GPU 2.3-3.3x the Orin's 6 cores on these cases.
* Items at depth 4: 11003. Items 0..1200 are DONE (Mac, logs/overnight). Remaining: [1200, 11003).
* Observed cost: 400 items took 50 min, 66 min and 5.4 h on 10 M1 Pro threads (the expensive region had
  a_2 = 3); expect 60-150 M1-core-hours for the remainder, i.e. roughly one to two hours on 96 cores.
* Suggested split (adjust to core counts): 96-core box items [1200, 8000), Mac Studio [8000, 11003).
  Keep CHUNK well above the thread count (a chunk ends when its slowest item ends):

```bash
./run_range.sh h9k8_box    9 8 5522 4 1200 8000  1000 96
./run_range.sh h9k8_studio 9 8 5522 4 8000 11003  600 24
```

* If a solution with n_h > 5521 appears anywhere, later chunks can be restarted with TGT = n_h + 1
  (delete nothing; just change TGT — the maximum over all solutions is still the answer).

### 2. A001209(303) = n(303,4)  — new term, needs psph16

* Lower bound **71148327** from the formula basis 1 228 17657 912312 (range verified). Target 71148328.
* Expected cost 2000-3300 core-hours (full checks dominate); ~1 day on 96 cores.
* Work items: 301 at depth 2 (one per a_2). Depth 3 would give ~3e6 items but generating them costs
  ~10 min of single-thread DP on every machine, so use depth 2; the ~100 heaviest items (a_2 between
  ~100 and ~230) each carry ~1 % of the work, which balances acceptably over 96 threads.

```bash
./run_range.sh h303k4_box 303 4 71148328 2 0 301 301 96 ./psph16     # one chunk = whole job (~1 day)
# or in two resumable halves:  ... 2 0 150 150 96 ./psph16   and   ... 2 150 301 151 96 ./psph16
```

* RAM: psph16 tables are 3 x (TGT+1) x 2 bytes = 430 MB per thread -> 41 GB at 96 threads. Fine on the
  2 TB box; use -j 48 elsewhere.

* After h=303: the same basis gives the formula values for h=304..306 (72060639, 72972951, 73885263); each
  is a separate run of similar cost.

### 3. A001211(27) = n(27,6)

* Lower bound **176380** (basis 1 19 177 1016 6649 22876, heuristic); target 176381.
* Items: 4005 at depth 3, 2,964,437 at depth 4 (generated in 2 s). Expected 300-3000 core-hours (growth
  exponent still rising at h=14), so start it only when a machine is otherwise idle; the ETA printed after
  the first chunks tells the truth. Depth 4 balances best:

```bash
./run_range.sh h27k6_box 27 6 176381 4 0 2964437 100000 96
```
* The h=26 extremal basis 1 19 177 816 6708 18060 has 27-range 174804 (a weaker but certain bound).

### 4. A001210(91) = n(91,5)

* Lower bound: the h=90 extremal basis 1 70 1412 32085 371775 has 91-range **8897042** (certain); the true
  value is probably ~9.0e6 (ratio ~1.06 per step). Target 8897043 (safe) — a better heuristic bound would
  speed the run a little. Expected 1600-8000 core-hours. Items: 133,730 at depth 3.

```bash
./run_range.sh h91k5_box 91 5 8897043 3 0 133730 5000 96
```

### Not for now
A001212(25) (4e17-8e17 nodes; needs the CUDA port and/or Challis's pruning), k=9 h=8 (~1e5 core-h),
the h=3/4/5 row proofs.

## Notes

* NITEMS for a given (H,K,TGT,D): `./psph -h H -k K -t TGT -d D -i 0:0` prints "split depth D: N work items".
  Item generation runs single-threaded on every machine before its range starts; it is instantaneous
  except for k=4 at depth 3 (see job 2).
* n(15,7) = A053346(15): the h=14 basis 1 12 52 225 546 3033 5464 has 15-range 29939 (certain lower bound);
  expected 1300-38000 core-hours — only after the jobs above.
* `-p 300` prints progress with an ETA every 5 minutes to the chunk's .err file.
* Memory per thread is about (k-1) x TGT bytes (x2 for psph16): negligible except for k=4 h=303
  (3 x 71 MB x 2 = 430 MB per thread, 41 GB for 96 threads — use fewer threads if RAM is short, e.g. -j 48).
* Do not mix item ranges from different D or different TGT in one logs/TAG directory.
