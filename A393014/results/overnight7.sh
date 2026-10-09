#!/bin/bash
# M4 Max (rev. 7, 2026-10-02): takes over from overnight6.sh at its next step boundary.  A060792
# (bases 2,3) has the best remaining odds per GPU-hour, so its next lengths come first.  L=113 and
# L=115 are split with the Studio: this machine takes parts in ascending order, the Studio in
# descending order, and parts already logged by the other machine are skipped.  A393014 L=55 last.
cd "$(dirname "$0")/.."
while pgrep -f "campaign.py run" >/dev/null; do sleep 20; done
C="python3 tools/campaign.py run"
$C --bases 2,3 --below 3^112                                                     # A060792 L=111 (+112 skip)
for i in 0 1 2 3 4 5 6 7; do $C --bases 2,3 --below 3^113 --Lmin 113 --part $i/8; done
$C --bases 3,10 --below 10^49                                                    # A046473 / A007633
for i in $(seq 0 15); do $C --bases 2,3 --below 3^115 --Lmin 114 --part $i/16; done
for i in 2 3 4 5 6 7; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done   # A393014
echo OVERNIGHT DONE
