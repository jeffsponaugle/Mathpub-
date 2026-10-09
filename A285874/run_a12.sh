#!/bin/bash
# Computes A285874 a(12) twice, with independent Zobrist keys and different task splits,
# then prints both results. Resumable: rerun after an interruption and finished tasks
# are read back from the checkpoints instead of recomputed.
#
# Start:  nohup caffeinate -i ./run_a12.sh > runs/run_a12.out 2>&1 &
# Stop:   pkill -f run_a12.sh; pkill -x perft
cd "$(dirname "$0")" || exit 1
run() { nice -n 10 ./perft --variant norooks --depth 12 --hash 8192 --report 300 "$@"; }
run --split 5 --seed 1 --ckpt runs/norooks_d12_seed1.ckpt >> runs/norooks_d12_seed1.log 2>&1 || exit 1
run --split 4 --seed 2 --ckpt runs/norooks_d12_seed2.ckpt >> runs/norooks_d12_seed2.log 2>&1 || exit 1
grep -h "perft(12)" runs/norooks_d12_seed1.log runs/norooks_d12_seed2.log
