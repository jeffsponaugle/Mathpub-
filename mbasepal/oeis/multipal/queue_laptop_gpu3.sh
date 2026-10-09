#!/bin/sh
# laptop GPU queue, stage 3: A171703 extension toward a(21) (upper part)
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_laptop_gpu2.sh" > /dev/null || pgrep -f "multipalg2 -k 4 -l 150000000000000000" > /dev/null; do sleep 60; done
./multipalg2 -k 4 -l 870000000000000000 -u 1168773800173056001 -r 15 -H 4 > runs/k4_gpu_8.7e17_to_1.169e18.txt 2> runs/k4_gpu_8.7e17_to_1.169e18.err < /dev/null
