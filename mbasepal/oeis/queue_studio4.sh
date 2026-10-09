#!/bin/sh
# Studio CPU: independent CPU confirmation of A171703 a(16)-a(18) (GPU-found), after queue_studio3
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_studio3.sh" > /dev/null; do sleep 60; done
cd multipal || exit 1
./multipal -k 4 -l 8400090327552001 -u 205045005215232001 -t 24 -r 16 > runs/k4_cpu_8.4e15_to_a18.txt 2> runs/k4_cpu_8.4e15_to_a18.err
echo "queue_studio4 finished" >> ../queue_studio.log
