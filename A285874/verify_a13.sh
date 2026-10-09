#!/bin/bash
# On the Mac Studio: the independent verification of A285874 a(13) -- a full single run with
# seed 2 and split 6 (the x86 machines compute seed 1, split 5, distributed).
# Start:  nohup caffeinate -i ./verify_a13.sh > runs/verify_a13.out 2>&1 &
# Stop:   pkill -f verify_a13.sh; pkill -x perft      (resumes from its checkpoint)
cd "$(dirname "$0")" || exit 1
while pgrep -x perft > /dev/null; do sleep 30; done
nice -n 10 ./perft --variant norooks --depth 13 --split 6 --seed 2 --hash "${HASH:-32768}" --report 600 \
  --ckpt runs/norooks_d13_seed2.ckpt < /dev/null >> runs/norooks_d13_seed2.log 2>&1
