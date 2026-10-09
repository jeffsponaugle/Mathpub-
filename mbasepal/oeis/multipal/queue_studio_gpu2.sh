#!/bin/sh
# Studio GPU queue, stage 2: A171703 extension toward a(21) (lower part)
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_studio_gpu.sh" > /dev/null || pgrep -f "multipalg2 -k 4 -l 8400090327552001" > /dev/null; do sleep 60; done
./multipalg2 -k 4 -l 374370997248000001 -u 870000000000000000 -r 15 -H 2 > runs/k4_gpu_3.744e17_to_8.7e17.txt 2> runs/k4_gpu_3.744e17_to_8.7e17.err < /dev/null
