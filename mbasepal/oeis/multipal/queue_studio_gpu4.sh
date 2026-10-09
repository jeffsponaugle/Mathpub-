#!/bin/sh
# Studio GPU queue, stage 4: A171703 toward a(23) (upper part, to 1965600^3)
cd "$(dirname "$0")" || exit 1
while pgrep -f "multipalg2 -k 4 -l 2164304937550479360" > /dev/null; do sleep 60; done
./multipalg2 -k 4 -l 5200000000000000000 -u 7594259452416000001 -r 21 -H 8 > runs/k4_gpu_5.2e18_to_7.594e18.txt 2> runs/k4_gpu_5.2e18_to_7.594e18.err < /dev/null
