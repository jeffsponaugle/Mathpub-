# a(24) search ledger — n=24, OEIS A171740

## RESOLVED 2026-09-23: a(24) = 1205824352934861911999927836207125

Found by the new meet-in-the-middle scanner (`-M`) in 50 seconds on one
machine: 24-digit palindrome in bases 24 (`lgi53dd79enccne97dd35igl`) and
26 (`3bd6fe8i3cbbbbc3i8ef6db3`), ≈1.21×10^33 (110 bits). It sits in stage
26 — stages 23, 24 and 25 are all EMPTY, so the fleet's stage-23 chunks
were never going to find it. The MITM run is a self-contained minimality
proof (it scanned stages 18–27 itself) and independently reproduces the
fleet's stage 18–22 clears exactly. Independently re-verified over bases
2–300. Fleet runs can be stopped; letting the stage-23 chunks finish adds
an independent-implementation cross-check of "stage 23 empty" but is not
required. Historical record of the fleet campaign below.


Minimality of any match requires every stage below it exhaustively cleared.
Stage windows are disjoint through stage 24, so cross-machine results
combine by simple minimum. Update this file as runs complete.

| stage | partners | candidates | machine | status |
|---|---|---|---|---|
| 18 | {16} | 4.2T | buildtest | **CLEAR** (empty) |
| 19 | {17} | 36.9T | buildtest | **CLEAR** (empty) |
| 20 | {18} | 121.9T | house1 (migrated from buildtest) | **CLEAR** (empty, 2026-09-04) |
| 21 | {19} | 315.7T | house1 | in flight — 89.8% @ 3.7G/s, ETA ~2.4h (2026-09-23) |
| 22 | {20} | 720.8T | M1 Max 0→78.5%, Studio 78.5%→100% | **CLEAR** (empty; composite coverage, Studio finished 2026-09-03) |
| 23 | {21} | 1516.2T | Studio (low), house1 (high, queued) | **low chunk CLEAR** (empty, 2026-09-23: no match below 35533897096292198492445107690441); mid/high split below |
| 24 | {22} | 3003.2T | — | pending (only if 23 empty; ~42h split 3 ways) |

## Stage 23 split (rebalanced 2026-09-03 for Studio + house1; M1 Max out)

```
low  (Studio, RUNNING): ./a171740 -n 24 -B 23 -L 35533897096292198492445107690441
mid  (Studio, next):    ./a171740 -n 24 -B 23 -H 668659832245794 -L 48374975418702820372293333762879
high (house1, after 20/21): ./a171740 -n 24 -B 23 -H 1254619147325272 -L 54108198377272584130510593262881
```

Mid = 586T (~21h Studio), high = 262T (~21h house1). A match in a higher
chunk is only a(24) if everything below is clear — take the minimum
across all matches.

## Notes

- Measured rates: M1 Max ~8.8G/s, Studio ~5.9–7.5G/s, buildtest ~3.7G/s,
  house1 TBD.
- On any interruption: restart with `-B <stage>` plus `-H <done-counter>`
  pasted from the last status line (margin is auto-subtracted).
- Migration caution: if a run was stopped mid-stage and restarted on
  another machine *without* `-H`, it redoes the stage from the start —
  correct, just slower. Gaps only occur if a `-H` value exceeds what was
  actually scanned.
- Cumulative find probability (Poisson heuristic): ~51% by end of 22,
  ~67% by 23, ~74% by 24. If 23 clears empty, consider building the
  meet-in-the-middle algorithm before committing to stage 24+.
