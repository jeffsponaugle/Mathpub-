#!/bin/sh
# laptop GPU queue, stage 2: A171703 extension (upper half) after the A171704(5) run
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_laptop_gpu.sh" > /dev/null || pgrep -f "multipalg -k 5" > /dev/null; do sleep 30; done
./multipalg2 -k 4 -l 150000000000000000 -u 374370997248000001 -r 15 -H 4 > runs/k4_gpu_1.5e17_to_3.744e17.txt 2> runs/k4_gpu_1.5e17_to_3.744e17.err < /dev/null
