#!/bin/bash
# M4 Max (rev. 10, 2026-10-03): takes over from overnight9.sh at its next step boundary.  The M4 Max
# runs bases 2,3 (v2) ~1.6x faster than a DGX Spark while the Sparks match it on the table method,
# so it takes A060792 L=117 parts (ascending from 0) and atom1 takes the non-coprime L=55 steps.
# A393014 L=55 parts afterwards.  Claims keep machines from running the same part.
cd "$(dirname "$0")/.."
while pgrep -f "campaign.py run" >/dev/null; do sleep 20; done
C="python3 tools/campaign.py run"
for i in $(seq 0 31); do $C --bases 2,3 --below 3^117 --Lmin 116 --part $i/32; done     # A060792 L=117
for i in $(seq 2 19); do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done      # A393014 L=55
echo OVERNIGHT DONE
