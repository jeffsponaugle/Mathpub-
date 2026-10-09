#!/bin/bash
# Studio queue: seed-1 confirmations, the independent naive-engine check of A285877 a(7),
# then the full A285874 a(13) verification run.
cd "$(dirname "$0")" || exit 1
SEEDS=1 HASH=24576 ./extend.sh jobs_confirm.txt
while pgrep -x perft > /dev/null; do sleep 30; done
grep -q "^total:" runs/nopawns_d7_divide_fast.txt 2>/dev/null ||
  nice -n 10 ./perft --variant nopawns --depth 7 --divide --hash 8192 > runs/nopawns_d7_divide_fast.txt
grep -q "^total:" runs/nopawns_d7_divide_naive.txt 2>/dev/null ||
  nice -n 10 ./naive "rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1" 7 24 > runs/nopawns_d7_divide_naive.txt
./verify_a13.sh
