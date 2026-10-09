#!/bin/bash
# M4 Max (rev. 9, 2026-10-03): four machines (M4 Max, Studio, two DGX Sparks).  Its own claimed
# A060792 L=115 part, then the non-coprime L=55 steps (the table method is no faster on a Spark),
# then A060792 L=117 parts ascending from 0, then A393014 L=55 parts.  Split lengths are shared by
# all machines in different orders; claims keep them from running the same part.
cd "$(dirname "$0")/.."
C="python3 tools/campaign.py run"
for i in $(seq 1 15); do $C --bases 2,3 --below 3^115 --Lmin 115 --part $i/16; done      # A060792 L=115
$C --bases 2,10  --below 10^55                                                      # A046472 / A007632
$C --bases 6,10  --below 10^55                                                      # A046475 / A029963
$C --bases 8,10  --below 10^55                                                      # A046477 / A029804
$C --bases 4,10  --below 10^55                                                      # A046474 / A029961
for i in $(seq 0 31); do $C --bases 2,3 --below 3^117 --Lmin 116 --part $i/32; done     # A060792 L=117
for i in $(seq 2 19); do $C --bases 2,9 --below 9^55 --Lmin 54 --part $i/20; done      # A393014 L=55
echo OVERNIGHT DONE
