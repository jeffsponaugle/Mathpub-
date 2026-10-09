#!/bin/bash
# Overnight campaign, M4 Max.  Starts when results/queue.sh has finished.  Resumable.
# Split lengths of bases 2,9: L=51 parts 0/2 (Studio 1/2); L=53 parts 0,1 of 5 (Studio 2,3,4);
# L=55 parts 0..7 of 20 (Studio 8..19).  All parts of a length must use the same plan (fingerprinted).
cd "$(dirname "$0")/.."
while pgrep -f "results/queue.sh" >/dev/null; do sleep 30; done
C="python3 tools/campaign.py run"
$C --bases 2,9  --below 9^51 --Lmin 50 --part 0/2   # A393014 / A259385
$C --bases 8,10 --below 10^52                       # A046477 / A029804 (re-run with the base-8 fix)
for i in 0 1; do $C --bases 2,9 --below 9^53 --Lmin 52 --part $i/5; done
$C --bases 9,10 --below 10^46                       # A046478 / A029965
$C --bases 3,10 --below 10^46                       # A046473 / A007633
$C --bases 2,10 --below 10^46                       # A046472 / A007632 cross-check
for i in 0 1 2 3 4 5 6 7; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done
echo OVERNIGHT DONE
