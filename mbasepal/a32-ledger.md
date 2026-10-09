# a(32) search ledger — n=32, OEIS A171740

## MATCH FOUND 2026-09-28 (Spark 2, CUDA):
## a(32)-designate = 48849507047691135113700496464447321657547821120
## bases 29 (mi0koes8osf75lpddpl57fso8seok0im) and 31 (2qm7gfg5h42ochl44lhco24h5gfg7mq2)
## ~4.885e46, 156 bits. Independently verified (bases 2-400 scan).
## SEALED 2026-09-28: ALL certificates complete. s29 CLEAR (upper +
## absorbed lower), s30 CLEAR (Ultra, 45.8T, cert "below 1706917413..."),
## s31 CLEAR below match (Spark1 lower cert "below 3819427559..." +
## Spark2 monotone sweep), s32 sliver CLEAR ([32^31, match), Spark2).
## a(32) = 48849507047691135113700496464447321657547821120 is FINAL.


All runs: MITM binary (w192 + two-level filter), `-M -n 32 -T <cap>`.
Stage windows are disjoint; a match in stage S is a(32) once all stages
below S are certified empty. Save each run's final `done:`/banner lines.

| stage | partners | steps (yh) | machine | status |
|---|---|---|---|---|
| ≤24 | — | small | M1 Max (first local run) | **CLEAR** (run reached stage 25) |
| 25 | {23} | 2.80T | house1 | **CLEAR** (empty, 2026-09-24) |
| 26 | {24} | 5.47T | mathf (10.1.7.25) | **CLEAR** (empty, 2026-09-24). Floor: a(32) > 27^31 ≈ 2.36e44 |
| 27 | {25} | 9.92T | house1 (lower) + mathf (upper) | **CLEAR** (both slices empty, 2026-09-25; certificates: house1 "below 4101734...", mathf "below 7276105..." — valid jointly). Floor: a(32) > 28^31 ≈ 7.28e44 |
| 28 | {26} | 17.07T | Studio CPU 0→14.8% + Ultra GPU 14.8→100% | **CLEAR** (empty, 2026-09-27; certificate "no ... below 2159424..." in a32-s28-gpu.log) |
| 29 | {27} | 28.24T | laptop GPU (upper, **CLEAR** 2026-09-28, cert in a32-s29hi-gpu.log) + 10.1.1.74 (lower 0→35.1%, killed) + laptop GPU absorbing [902.1G, 2569.7G) (a32-s29-absorb.log) | lower absorption in flight (~12 min) — stage certified when it completes |
| 30 | {27,28} | ~46T (t=6) | **Studio Ultra GPU** (fresh full stage, launched after s28 clear; log a32-s30-gpu.log) | in flight. 10.1.10.20's CPU run confirmed stopped (was 1.7% of p28); box retired from campaign. |
| 31 | {28,29} | 75.7T (t=6) | **Spark1 + Spark2 (CUDA, GB10)** | launched 2026-09-27 @ ~1.27G/s each, ~8h. Spark1 (10.1.30.36): p28 + p29-lower, `-B 31 -E 32 -L 38194275598097520417729325221884039804056785392 -T 900000000` (a32-s31lo.log). Spark2 (10.1.30.37): p29-upper, `-B 31 -H 32722128008401 -P 29 -E 32 -T 900000000` (a32-s31hi.log). CUDA worker validated 24/24. |
| 29+ | ... | ~2× per stage | — | pending |

## Commands

- New stage N on a box: `./a171740 -M -n 32 -B N -T 700000000 2>&1 | tee a32-sN.log`
  (`-T 300000000` on 16GB machines.)
- Mid-stage resume (new binary required): `-B N -H <done-counter>`, plus
  `-P <b1>` on multi-partner stages to name the partner phase the counter
  came from (see the last `MITM b1=X:` banner in the log).
- To cap a run so it doesn't roll into a stage another box owns: add
  `-E <next stage b2>` (new binary; e.g. house1's stage-27 run should be
  `-B 27 -E 28`). The old `-L <lo of next stage>` works on old binaries.

## Certified floor

Everything below stage 25's window. Updates as stages complete; each
empty stage's final banner line is the certificate — keep the logs.
