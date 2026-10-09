#!/bin/sh
# Studio CPU queue 2: independent CPU confirmation of A171704(5) and (6),
# after the A216 L=16/17 batch in queue_studio.sh
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_studio.sh" > /dev/null; do sleep 60; done
cd multipal || exit 1
./multipal -k 5 -l 1000000000000 -u 4922057407205377 -t 24 -r 5 > runs/k5_cpu_1e12_to_a5.txt 2> runs/k5_cpu_1e12_to_a5.err
./multipal -k 5 -l 4922057407205377 -u 69554307029409793 -t 24 -r 5 > runs/k5_cpu_a5_to_a6.txt 2> runs/k5_cpu_a5_to_a6.err
echo "queue_studio2 finished" >> ../queue_studio.log
