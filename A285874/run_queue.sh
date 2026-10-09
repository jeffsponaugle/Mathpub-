#!/bin/bash
# One machine's share of the work: SEEDS=1 on the M4 Max, SEEDS=2 on the M5 Pro.
# Start:  SEEDS=1 HASH=8192 nohup caffeinate -i ./run_queue.sh > runs/queue.out 2>&1 &
# Stop:   pkill -f run_queue.sh; pkill -f extend.sh; pkill -x perft; pkill -x naive
cd "$(dirname "$0")" || exit 1
./extend.sh jobs_quick.txt || exit 1
if [ "${NAIVE:-0}" = 1 ] && ! grep -q "^total:" runs/nopawns_d7_divide_naive.txt 2>/dev/null; then
  nice -n 10 ./perft --variant nopawns --depth 7 --divide --hash "${HASH:-8192}" > runs/nopawns_d7_divide_fast.txt
  nice -n 10 ./naive "rnbqkbnr/8/8/8/8/8/8/RNBQKBNR w KQkq - 0 1" 7 "$(getconf _NPROCESSORS_ONLN)" > runs/nopawns_d7_divide_naive.txt
fi
./extend.sh jobs_long.txt
