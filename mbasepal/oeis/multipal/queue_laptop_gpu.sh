#!/bin/sh
# laptop GPU queue: A171704(5) lower range after the A171703 proof run
cd "$(dirname "$0")" || exit 1
while pgrep -f "multipalg -k 4" > /dev/null; do sleep 20; done
./multipalg -k 5 -l 1000000000000 -u 35000000000000000 -r 5 > runs/k5_gpu_1e12_to_3.5e16.txt 2> runs/k5_gpu_1e12_to_3.5e16.err < /dev/null
