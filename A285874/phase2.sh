#!/bin/bash
# On a server in the distributed A285874 a(13) run: after this machine's current part
# finishes, computes its phase-2 tasks, appending to the same checkpoint file.
# Usage: HASH=<MB> setsid nohup ./phase2.sh NAME CKPT > runs/phase2.out 2>&1 < /dev/null &
cd ~/src/math/A285874 || exit 1
while pgrep -f "[p]erft3* --variant norooks --depth 13" > /dev/null; do sleep 60; done
numactl --interleave=all nice -n 10 ./perft --variant norooks --depth 13 --split 5 --seed 1 \
  --hash "${HASH:-614400}" --report 300 --tasks "tasks_$1.txt" --ckpt "$2" < /dev/null \
  >> "runs/norooks_d13_seed1.${1%2}.log" 2>&1
