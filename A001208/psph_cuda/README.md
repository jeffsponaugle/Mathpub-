# psph_cuda — CUDA leaf for the postage-stamp exhaustive search

`gpu_leaf.cu` + `gpu_leaf.h` move the **leaf** of `../psph/psph.c` (the test of every candidate `a_k` of a
surviving `(k-1)`-prefix: Challis X test, deep-hole / reach / top-block targets, exact full check) onto an
NVIDIA GPU. The prefix DFS, the lazily extended min-stamps tables, the work-item partitioning, the
`SOLUTION` verification (`report_solution`, an independent exact DP) and all statistics stay in `psph.c`;
a `-gpu` mode (default in the GPU build) swaps the CPU candidate loop for `gpu_submit_prefix()`.

Tested on a Jetson Orin Nano (orin1: sm_87, 8 SMs / 1024 cores, 7.5 GB shared, CUDA 13.2); built for the
DGX Spark (GB10) with `ARCH=-arch=sm_121`.

## Build and run

```
cd psph_cuda
make                       # psph_gpu (8-bit cells, h <= 254) and psph_gpu16 (-DWIDE, 16-bit cells, h <= 65534)
make psph_cpu psph_cpu16   # reference CPU binaries from the same psph.c (no GPU code), for differential tests
make ARCH="-arch=sm_121"   # DGX Spark GB10 (Blackwell); "-arch=native" also works with CUDA >= 12.8
```

`psph.c` is compiled with `gcc -O3 -march=native -DGPU_LEAF` (NEON on aarch64, AVX2 on x86 — identical to
the CPU build), `gpu_leaf.cu` with `nvcc -O3 $(ARCH)`, and the two are linked by nvcc. The Makefile uses
`../psph/psph.c` unless a `psph.c` copy exists in the directory (deployment layout, as on orin1).

```
./psph_gpu -h 9 -k 8 -t 5600 -j 6 -d 4 -i 10900:11003      # same options as psph; leaf on the GPU
./psph_gpu -cpu ...                                        # same binary, CPU leaf (PSPH_GPU=0 does the same)
./psph_gpu -selftest quick -j 6                            # 422 checks against the published tables, GPU leaf
./psph_gpu16 -h 303 -k 4 -t 71148328 -j 6 -d 2 ...         # 16-bit cells for h >= 255
```

`-j` is the number of CPU threads; each owns a batch and a CUDA stream, so 4–8 threads keep the GPU fed
while the others build tables. All `psph` options (`-d`, `-i`, `-b`, `-s2`, `-s2cap`, `-p`, `-v`) work
unchanged; `-s2 0` disables the deep-hole targets on both sides (useful for exact differential comparisons).
GPU knobs are environment variables (defaults in parentheses):

| variable | meaning |
|---|---|
| `PSPH_GPU_BATCH_MB` (64) | table bytes per batch and thread (pinned host + device copy) |
| `PSPH_GPU_BATCH_PREFIXES` (8192) | prefixes per batch |
| `PSPH_GPU_MAX_CAND` (2^26) | candidates per batch (a larger prefix range is split, exactly) |
| `PSPH_GPU_SCRATCH_MB` (64) | full-check scratch per thread; grows automatically if one survivor needs more |
| `PSPH_GPU_FULL` (auto) | `thread` / `warp`: full-check kernel flavour (auto: warp when < 2048 survivors in a chunk or `a_k` > 65536) |
| `PSPH_GPU_STRATA` (4) | deep-hole strata per prefix (see below) |
| `PSPH_GPU_INT64=1` | force the 64-bit arithmetic filter kernel (default: 32-bit when provably overflow-free) |
| `PSPH_GPU_PARANOID=1` | re-run every batch on the host with the same (shared) kernel code and abort on any mismatch |
| `PSPH_GPU_VERBOSE=1` | one log line per batch |

## Design

```
psph.c worker thread t                     gpu_leaf.cu, stream t
----------------------                     ------------------------------------------
prefix DFS (unchanged)                     batch_t: pinned host buffers + device mirrors
leaf(): extend T_{k-1} to E,                 tables  T_p[0..E_p] (+32 B zero pad, 16 B aligned)
        reach r[0..h], Z, deep holes  --->   meta    {tab_off, E, lo, hi, a_{k-1}, Z, deep_off, ndeep, reach_off}
        gpu_submit_prefix(...)               cst[p]  exclusive prefix sum of candidate counts
                                             deep[], reach[], a[1..k-1] copies (host only)
   batch full (64 MB / 8192 prefixes / 2^26 candidates) or gpu_flush():
                                           H2D copy  ->  kernel 1 (filter)  ->  D2H survivors  ->  sort
                                           -> chunks that fit the scratch -> kernel 2 (full check) -> D2H gaps
report_solution_basis()  <---  callback    for every survivor with gap == 0
```

**Kernel 1 — filter, one thread per candidate.** The flat thread index `t` in `[0, total)` is mapped to
`(p, j)` by a binary search over `cst[]` (`cst[p] <= t < cst[p+1]`, 13 steps for 8192 prefixes, served from
L1/L2 because the array is tiny and hot); the candidate is `a = lo_p + j`. This was chosen over "one block
per prefix with a grid-stride loop" because the candidate counts per prefix vary from 1 to ~6·10^5
(k=5, h=91): block-per-prefix leaves most lanes idle on small prefixes and serialises the large ones,
whereas the flat mapping is perfectly balanced, and consecutive threads still handle consecutive `a` of the
same prefix, so a warp is almost always prefix-uniform (same table, similar `C_k`, little divergence).
Each thread runs `filter_one()`: `C_k = TGT/a`, `C_{k-1} = a/a_{k-1}`,
`X = (C_k-1)a + (C_{k-1}-1)a_{k-1} + Z`, the `c` loop from `C_k-1` down to
`cmin = ceil((X - h a_{k-1})/(a - a_{k-1}))` (clipped at 0) exactly as the stage-1 code of `psph.c`; then
`leaf_rep`-identical tests of `(C_k-1)a + y` and (if `C_k >= 3`) `(C_k-2)a + y` for the prefix's deep holes
with `y < a`; the reach test (`ma` = min `m` with `r[m] >= a-1` by scanning `r[]`, `c* = h-ma+1`, if
`c* <= C_k` then `xs = c* a + r[ma-1] + 1` must be representable when `xs <= TGT`, otherwise `direct++`);
the top-block deep holes `C_k a + y <= TGT`; and `TGT` itself. Survivors are appended with `atomicAdd`
(if the list overflows, it is enlarged and the kernel is simply re-run — a survivor is never dropped).
Counters (X-pass, deep-hole-pass, stage-2-pass, direct, probes) are warp-reduced with shuffles and added
once per warp. The kernel is a template on the integer type: 32-bit arithmetic is used when every
intermediate provably fits (`TGT + max hi < 2^31` and `h·max a_{k-1} + TGT < 2^31`, checked per batch —
true for every case in the tables including k=4 h=303 and k=5 h=91), otherwise 64-bit.

**Kernel 2 — full check.** `fullcheck_core<NL>()` is a line-by-line port of `full_check()`:
the streaming block DP `u[x] = min(T[x], u[x-a]+1)` over `[0, E]` runs **in place** in a buffer of `a`
cells (block `b` overwrites block `b-1` at the same offsets, so a ring of exactly `a` cells needs no wrap
handling), 4 bytes per step with the CUDA SIMD-in-a-word intrinsics (`__vaddus4/__vminu4/__vcmpgtu4`, or the
`…2` variants for 16-bit cells); the unaligned read `T[b·a + i]` is assembled from two aligned words with
`__funnelshift_r`. The last word of a partial block merges the previous block's cells back in (the window
formula needs them), and gap detection masks cells beyond the block. Early exit at the first gap (the
minimum over the block, i.e. the same position the CPU returns). After the DP, the first gap above `E` is
`min_{z in (E-a, E]} z + (h+1-u[z]) a`, taken at the smallest `z` carrying the maximum `u` — the window is
exactly the buffer read in rotated order (`i >= len_last` first, from the block before the last, or from `T`
when only one block exists; then `i < len_last`), and when `a > E` it is `T[1..E]` itself.
Two flavours share this code: `NL=1` (one thread per survivor, lane 0 only) and `NL=32` (one **warp** per
survivor; the 32 lanes stride over the words of a block, which is dependency-free within a block, with
`__shfl_xor_sync` reductions for the first gap and the window maximum). Survivors are sorted (so neighbouring
threads share tables and the output order is deterministic per batch) and processed in chunks whose scratch
(`roundup(a·CS, 4) + 16` bytes each, 16-byte aligned offsets) fits the per-thread scratch buffer; a single
survivor larger than the buffer makes the buffer grow.

**Host pipeline.** Per CPU thread and batch: `cudaMemcpyAsync` of tables/meta/deep/reach (pinned),
kernel 1, survivor count back, (re-run on overflow), survivor list back, `std::sort`, chunked kernel 2 with
results back, callbacks for `gap == 0`, statistics. Everything is on the thread's own non-blocking stream
and `cudaDeviceScheduleBlockingSync` makes waiting threads sleep, so the CPU threads overlap their table
building with other threads' GPU work and the process coexists with other CPU jobs.

### Memory layout (per CPU thread, defaults)

| buffer | size | contents |
|---|---|---|
| `h_tab` / `d_tab` | 64 MB (grows if one table is larger) | concatenated tables `T_p[0..E_p]`, each padded with 32 zero bytes to a 16-byte boundary |
| `h_meta` / `d_meta` | 8192 × 64 B | `gpu_prefix` records |
| `h_cst` / `d_cst` | 8193 × 8 B | candidate prefix sums |
| `h_deep` / `d_deep` | variable (uint32 per hole) | deep-hole remainders, ascending per prefix |
| `h_reach` / `d_reach` | 8192 × (h+1) × 8 B | reach sequences |
| `h_a` | 8192 × 32 × 8 B (host only) | prefix copies for the callback |
| `d_surv` | grows | survivor flat indices (uint32) |
| `d_scratch` | 64 MB (grows) | full-check ring buffers |
| chunk arrays | survivors × 24 B | `(p, j, offset)` per survivor + results |

For k=8 h=9 (E ≤ 5600) a 64 MB batch holds 8192 prefixes ≈ 1–4·10^6 candidates; for k=5 h=91 (E = TGT ≈ 8.9·10^6)
it holds 7 prefixes (≈ 10^5–6·10^5 candidates each); for k=4 h=303 with 16-bit cells one table is 142 MB and
the buffer grows to hold one prefix per batch (≈ 10^6 candidates). The batch is always worth a launch.

## Why the results are exact

Every GPU rejection is justified by the **exact unrepresentability of one number t ≤ TGT** or by the exact
full check, exactly as in the CPU program:

* `T[0..E]` is the exact min-stamps table of the prefix (built by `psph.c`, copied verbatim; the extension to
  `E = min(TGT, h·a_{k-1})` is the CPU's own `extend()`). A probe `T[t - c·a] + c <= h` with `t - c·a <= E`
  proves representability; indices above `E` are treated as unrepresentable, which is correct because
  `T[y] >= y/a_{k-1} > h` for `y > h·a_{k-1}` and `t - c·a <= t <= TGT = E` otherwise. The `cmin` bound only
  skips probes that cannot succeed (`t - c·a > (h-c)·a_{k-1}` implies `T[t-c·a] + c > h`). A target is
  declared unrepresentable only when **all** admissible `c` fail — a false "unrepresentable" would need a
  wrong table entry or a wrong index, never a wrong choice of targets.
* The X test, the reach test and the TGT test use the same formulas and loop bounds as `psph.c`, so the
  **X-pass counts match the CPU exactly** (verified below), and with `-s2 0` every counter (stage-2 pass,
  full-check fails, direct, probes) matches exactly.
* Deep holes: the set of tested targets may differ from the CPU (chosen per prefix instead of per candidate),
  which changes only *which* numbers are tested, not the exactness of each test. The GPU list is the union
  over strata `[1, lo)` and `GPU_STRATA` equal parts of `[lo, hi)` of the `NDH` positions with the largest
  `T[y]` (one per residue class mod `a_{k-1}`, same `deep_insert` rule and `NDH = clamp(E/1024, 8, s2cap)`
  as the CPU), sorted ascending; a candidate uses the entries with `y < a`. This is at least as strong as
  the CPU's list for most candidates (e.g. k=5 h=9: 55 full checks instead of 1705, at +14 % probes).
* Full check: exact DP over `[1, E]` with saturating 8/16-bit adds (values > h are never trusted as finite),
  exact formula above `E` (proof in `psph.c`); `PSPH_GPU_PARANOID=1` re-runs the same `__host__ __device__`
  code sequentially on the host for every candidate and every survivor and compares survivor lists, gap
  positions and all counters, which guards against indexing, batching and race errors in the GPU
  execution (logic errors are caught by the differential tests against the real CPU program).
* Every reported basis is still verified by `report_solution_basis()` (`hrange_exact`, independent DP); a
  basis with `n_h < TGT` aborts the program, as before.
* Overflow: positions are 64-bit on the host and in kernel 2; kernel 1 uses 32-bit arithmetic only under the
  per-batch bound above (the 64-bit instantiation is selected otherwise and can be forced for testing).
* Candidate ranges `[lo, hi]`, `E`, `Z`, `r[]` and all bounds are computed by the unchanged `psph.c` code.

## Validation

### Validation (orin1, 2026-10-01)

* `./psph_gpu -selftest quick -j 6` with the GPU leaf active: **422 tests, 0 failures** (81 search cases run at TGT = n and
  TGT = n+1, 341 formula/table checks). Log `logs/selftest_gpu.log`.
* Differential GPU vs CPU (same `psph.c`, CPU build without GPU code). "X-pass identical" means the stage-1 Challis test made
  the identical decision on every candidate; solution sets are compared as sorted files. Stage-2 counts differ by design
  (per-prefix stratified deep-hole targets on the GPU reject 20-100x more candidates before the full check).

| case | items | solutions | X-pass decisions compared | GPU vs CPU wall (shared CPU) |
|---|---|---|---|---|
| (12,6) TGT 5118 | full search | identical (1) | 1,683,140,899 identical | 52 s vs 122 s |
| (12,6) TGT 5118, -s2 0 (all counters must match) | full | identical (1), all counters equal | 1,683,140,899 | 253 s vs 306 s |
| (9,5) TGT 797, 16-bit build | full | identical (1) | 217,557 | ~0 vs 1 s |
| (9,5) TGT 797, -s2 0, 16-bit; int64+PARANOID variants | full | identical, paranoid host re-check silent | 217,557 | |
| (12,5) TGT 2047, int64+PARANOID | full | identical | | |
| (9,7) TGT 3191 / 3192 | -d 3 items 144-152 (extremal prefix {1,7,30}) | identical (1 / 0) | 2.81e9 / 2.80e9 | 57 vs 184 s / 56 vs 179 s |
| (14,6) TGT 9748 | -d 3 items 440-600 (extremal prefix {1,11,49}) | identical (1) | 2.92e9 | 67 vs 220 s |
| (8,8) TGT 3485 | -d 4 items 3000-3002 | identical (0) | 3.95e9 | 153 vs 456 s |
| (9,8) TGT 5600 | -d 4 items 10995-11003 | identical (0) | 1.47e10 | |
| (9,8) TGT 5600 | -d 4 items 10980-10995 | identical (0) | 4.21e10 | |
| (9,8) TGT 5521 | -d 4 item 1200 | identical (0) | 6.20e9 | |
| (9,8) TGT 5781 | -d 4 items 5203-5205 | identical (0) | 2.58e9 | 177 vs 254 s |

Two further (9,8) cases in the original script hit its 15-minute per-run timeout on both binaries (items 1195-1205 and
10900-11003); they are not comparisons and were replaced by the small-range cases above. No completed case differs.


## Throughput

### Throughput (orin1: 1024-core Ampere GPU at ~1 GHz, 6 x Cortex-A78AE host cores shared with a production job)

| case | CPU binary wall (6 thr) | GPU binary wall | in-kernel filter rate | note |
|---|---|---|---|---|
| (12,6) TGT 5118 | 122 s | 52 s | 2.8e8 cand/s | GPU time 33 s = H2D 3 + filter 22 + full check 8; host pipeline is the bottleneck |
| (9,7) TGT 3191 items 144-152 | 184 s | 57 s | | |
| (14,6) TGT 9748 items 440-600 | 220 s | 67 s | | |
| (8,8) TGT 3485 items 3000-3002 | 456 s | 153 s | | |
| (9,8) TGT 5600 items 10900-11003 | > 900 s (timeout, 0 bases) | > 900 s (timeout, 6 bases found) | | GPU covered several times more items in the same time |
| (27,6) TGT 176381 -d 3 items 0-20 | 94 s | 93 s | 7.3e8 cand/s | host-bound: the 6 Orin cores cannot build prefix tables fast enough |
| (25,5) TGT 31108 | 79 s | 23 s (warp full check) / 26 s (thread) | 2.1e8 / 1.6e8 cand/s | warp-per-survivor full check 13 s vs 32 s |
| (60,4) TGT 143814 | 6 s | 9 s (warp) / 41 s (thread) | 6.0e7 / 3.4e7 | full-check-bound at small E: the Orin GPU's bandwidth loses to 6 CPU cores with NEON |

Reading: the filter kernel runs at 2-7e8 candidates/s on this small GPU (5-13x the Orin's own CPU), so on the Orin the
end-to-end speed is set by the host side (prefix DFS, table extension to E, deep-hole/reach preparation), which the six slow
cores deliver at roughly 1-3e8 candidates/s worth of prefixes. A DGX Spark (20 fast ARM cores, 6144-core GPU, ~10x the
bandwidth) should remove both limits; expect several x mathd per Spark on the probe-dominated searches (k=6..9, h near k).


## What limits scaling

### Scaling and limits

* Batch = 64 MB of tables by default: 8192 prefixes at k=8 h=9 (E ~ 5.6 KB each), 7 prefixes at k=5 h=91 (E = TGT ~ 8.9 MB
  each), 1 prefix at k=4 h=303 (142 MB 16-bit; the buffer grows automatically). Candidate counts per prefix (1..6e5) are
  balanced by the flat index -> (prefix, j) mapping.
* Full-check scratch (64 MB default) holds a ring of a_k cells per survivor: ~95-640 survivors at a_k ~ 1e5-7e5, so for
  k=5 h=91 the warp-per-survivor flavour is required (implemented, auto-selected) and PSPH_GPU_SCRATCH_MB should be 256-512
  on a Spark. Full checks at E ~ 8.9e6 cells move ~27 MB each and are bandwidth-bound there, as on the CPU.
* GPU mode requires TGT < 2^31 (uint32 deep holes / survivor indices; batches of <= 2^26 candidates, larger ranges split).
* Host side: keep one CPU thread per stream; on a weak host (Orin) the CPU prefix work limits throughput, on a strong host the
  GPU does. `-j` = number of host threads = number of streams.


## Files

* `gpu_leaf.cu`, `gpu_leaf.h` — the CUDA leaf (single translation unit + C API).
* `Makefile` — `psph_gpu`, `psph_gpu16`, `psph_cpu`, `psph_cpu16`.
* `tests/validate.sh`, `tests/throughput.sh` — the scripts that produced the tables above (run on orin1).
* `../psph/psph.c` — the search program; the `-gpu` additions are all inside `#ifdef GPU_LEAF`.
