#!/bin/bash
# M4 Max (rev. 8, 2026-10-03): takes over from overnight7.sh at its next step boundary.  Finishes the
# A060792 L=115 parts with the Studio, then the best remaining lengths per GPU-hour (user-approved):
# A046479 L=46-47 and the non-coprime L=55 steps (this machine is the faster one at those), then
# A060792 L=117 split with the Studio (ascending here, descending there; claims prevent duplicates).
# A393014 (bases 2,9) L=55 last.
cd "$(dirname "$0")/.."
while pgrep -f "campaign.py run" >/dev/null; do sleep 20; done
C="python3 tools/campaign.py run"
for i in $(seq 1 15); do $C --bases 2,3 --below 3^115 --Lmin 115 --part $i/16; done      # A060792 L=115
$C --bases 10,11 --below 11^47                                                      # A046479 / A029966
$C --bases 2,10  --below 10^55                                                      # A046472 / A007632
$C --bases 6,10  --below 10^55                                                      # A046475 / A029963
$C --bases 8,10  --below 10^55                                                      # A046477 / A029804
$C --bases 4,10  --below 10^55                                                      # A046474 / A029961
for i in $(seq 0 31); do $C --bases 2,3 --below 3^117 --Lmin 116 --part $i/32; done     # A060792 L=117
for i in 2 3 4 5 6 7; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done      # A393014
echo OVERNIGHT DONE
