#!/bin/sh
# Studio GPU queue, stage 3: A171703 toward a(22) (upper part, to 1441440^3)
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_studio_gpu2.sh" > /dev/null || pgrep -f "multipalg2 -k 4 -l 374370997248000001" > /dev/null; do sleep 60; done
./multipalg2 -k 4 -l 2100000000000000000 -u 2994950912937984001 -r 19 -H 2 > runs/k4_gpu_2.1e18_to_2.995e18.txt 2> runs/k4_gpu_2.1e18_to_2.995e18.err < /dev/null
