#!/bin/sh
# Studio GPU queue: A171703 extension (lower half) after the A171704(5) run
cd "$(dirname "$0")" || exit 1
while pgrep -f "multipalg -k 5" > /dev/null; do sleep 30; done
./multipalg2 -k 4 -l 8400090327552001 -u 150000000000000000 -r 15 -H 2 > runs/k4_gpu_8.4e15_to_1.5e17.txt 2> runs/k4_gpu_8.4e15_to_1.5e17.err < /dev/null
