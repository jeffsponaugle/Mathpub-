#!/bin/bash
# On mathd: merges the distributed A285874 a(13) run (seed 1, split 5) once the checkpoints of
# its parts together cover every task: mathd (share + phase 2), mathb (share + phase 2) and
# mathg. Retries every 10 minutes; a partial final line in a file still being written is
# ignored by the reader. The merge re-verifies every entry and stops on any conflict.
# Start:  setsid nohup ./merge_a13.sh > runs/merge_a13.out 2>&1 < /dev/null &
cd ~/src/math/A285874 || exit 1
while true; do
  scp -q -o BatchMode=yes jbs@10.1.30.21:src/math/A285874/runs/norooks_d13_seed1.mathb.ckpt runs/
  scp -q -o BatchMode=yes jbs@10.1.30.26:src/math/A285874/runs/norooks_d13_seed1.mathg.ckpt runs/
  ./perft --variant norooks --depth 13 --split 5 --seed 1 \
    --merge runs/norooks_d13_seed1.ckpt,runs/norooks_d13_seed1.mathb.ckpt,runs/norooks_d13_seed1.mathg.ckpt \
    > runs/merge_a13.last 2>&1
  if grep -q -E "^perft\(13\) = |CONFLICT" runs/merge_a13.last; then
    cat runs/merge_a13.last >> runs/norooks_d13_seed1.log
    break
  fi
  sleep 600
done
