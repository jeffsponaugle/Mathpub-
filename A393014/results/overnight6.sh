#!/bin/bash
# M4 Max (rev. 6, 2026-10-01): the prime-sequence lengths with the best expected new terms per GPU-hour
# (heuristic expected count from tools/stats.py / dbpal plan estimate), roughly best-first.  These pairs
# are listed in the Studio's results/remote_pairs.txt so it skips them.  Bases 2,9 L=55 moved to the end.
cd "$(dirname "$0")/.."
C="python3 tools/campaign.py run"
$C --bases 2,10  --below 10^52     # A046472 / A007632
$C --bases 10,11 --below 10^46     # A046479 / A029966   (base-11 L=41..45)
$C --bases 7,10  --below 10^46     # A046476 / A029964
$C --bases 10,13 --below 13^41     # A046481 / A029968
$C --bases 9,10  --below 10^47     # A046478 / A029965
$C --bases 2,3   --below 3^109     # A060792
$C --bases 3,10  --below 10^47     # A046473 / A007633
$C --bases 7,10  --below 10^47
$C --bases 2,10  --below 10^54
$C --bases 9,10  --below 10^49
$C --bases 2,3   --below 3^111
$C --bases 3,10  --below 10^49
$C --bases 2,3   --below 3^113
for i in 2 3 4 5 6 7; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done
echo OVERNIGHT DONE
