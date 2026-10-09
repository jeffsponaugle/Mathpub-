#!/bin/bash
# M4 Max (rev. 5): term-productive family pairs first (OEIS data ends near 10^18), then bases 2,9 L=55.
cd "$(dirname "$0")/.."
C="python3 tools/campaign.py run"
for b in 10,20 10,24 10,27 10,32 10,36 7,29; do $C --bases $b --below 10^40; done
for i in 1 2 3 4 5 6 7; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done
for i in 19 18 17 16 15 14 13 12 11 10 9 8; do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done
echo OVERNIGHT DONE
