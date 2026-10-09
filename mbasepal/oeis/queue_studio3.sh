#!/bin/sh
# Studio CPU: A171704(6)-range CPU confirmation, after the A216 batch and the a(5) confirmation
cd "$(dirname "$0")" || exit 1
while pgrep -f "queue_studio.sh" > /dev/null || pgrep -f "multipal -k 5 -l 1000000000000 " > /dev/null; do sleep 60; done
cd multipal || exit 1
./multipal -k 5 -l 4922057407205377 -u 34777153514704897 -t 24 -r 5 > runs/k5_cpu_a5_to_a6.txt 2> runs/k5_cpu_a5_to_a6.err
echo "queue_studio3 finished" >> ../queue_studio.log
