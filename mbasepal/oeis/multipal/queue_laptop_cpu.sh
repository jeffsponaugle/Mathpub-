#!/bin/sh
# laptop CPU queue: A171702 (3-digit palindromes in >= n bases) beyond Resta's a(100) = 13967553601
cd "$(dirname "$0")" || exit 1
while pgrep -f "batch.sh 10" > /dev/null; do sleep 60; done
./multipal -k 3 -l 13967553601 -u 1000000000000 -t 10 -r 101 > runs/k3_1.4e10_to_1e12.txt 2> runs/k3_1.4e10_to_1e12.err < /dev/null
