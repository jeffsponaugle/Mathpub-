#!/bin/sh
# Studio CPU: A171703 a(19)-a(21) CPU confirmation, after the Studio GPU a(22) range finishes
cd "$(dirname "$0")" || exit 1
while pgrep -f "multipalg2 -k 4 -l 2164304937550479360" > /dev/null; do sleep 60; done
./multipal -k 4 -l 205045005215232001 -u 1168773800173056001 -t 14 -r 19 > runs/k4_cpu_a18_to_a21.txt 2> runs/k4_cpu_a18_to_a21.err
echo "queue_studio5 finished" >> ../queue_studio.log
