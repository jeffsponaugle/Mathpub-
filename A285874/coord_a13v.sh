#!/bin/bash
# Coordinates the A285874 a(13) verification (run 2: seed 2, split ply 6, N tasks). The Mac
# Studio works through the task list from the front (verify_a13.sh). A back group (mathd,
# mathb, mathg, M5 Pro) works from the end toward it, each on its own tasks_v13_<m>.txt. The
# lists were rebalanced once from measured rates; the work finished before that is in each
# machine's checkpoint. Runs on the laptop, the only machine with SSH access to all of them.
#
# Every 5 minutes it checks two triggers:
#   - everyone's progress counts add up to N (the two fronts have met), or
#   - no back-group process is still running (their lists cover everything right of the Studio).
# Either way it gathers all checkpoints on mathd and merges them. The merge re-verifies every
# entry, requires every task, and requires tasks computed twice (where the fronts meet) to
# agree. After a successful merge it records the result in mathd's runs/norooks_d13_seed2.log
# and stops all five machines. A conflict is reported and nothing is stopped.
STUDIO=jeffsponauglejbs@10.1.10.6
M5=jbs@10.1.1.43
MD=jbs@10.1.30.23
MB=jbs@10.1.30.21
MG=jbs@10.1.30.26
D=src/math/A285874
N=8095245
DONE0=985120  # tasks the back group had finished before the rebalance (excluded from new lists)
RUN='[p]erft --variant norooks --depth 13 --split 6'
last() { ssh -o BatchMode=yes -o ConnectTimeout=8 "$1" "grep 'tasks (' ~/$D/runs/$2 | tail -1" < /dev/null 2>/dev/null | awk '{split($1,a,"/"); print a[1]+0}'; }
group_running() {  # "no" only if every back-group machine answered and none is running
  local h r
  for h in "$MD" "$MB" "$MG" "$M5"; do
    r=$(ssh -o BatchMode=yes -o ConnectTimeout=8 "$h" "pgrep -f '$RUN' > /dev/null && echo yes || echo no" < /dev/null 2>/dev/null)
    [ "$r" = no ] || { echo yes; return; }
  done
  echo no
}
while true; do
  v=$(last "$STUDIO" norooks_d13_seed2.log)
  c=$(( ${v:-0} + DONE0 ))
  for hm in "$MD mathd" "$MB mathb" "$MG mathg" "$M5 m5pro"; do
    set -- $hm
    v=$(last "$1" "norooks_d13_seed2.$2.r2.log")
    c=$(( c + ${v:-0} ))
  done
  if [ "$c" -ge "$N" ] || [ "$(group_running)" = no ]; then
    ok=1
    scp -3 -q -o BatchMode=yes "$STUDIO:$D/runs/norooks_d13_seed2.ckpt" "$MD:$D/runs/v13_studio.ckpt" || ok=0
    scp -3 -q -o BatchMode=yes "$M5:$D/runs/norooks_d13_seed2.m5pro.ckpt" "$MD:$D/runs/v13_m5pro.ckpt" || ok=0
    scp -3 -q -o BatchMode=yes "$MB:$D/runs/norooks_d13_seed2.mathb.ckpt" "$MD:$D/runs/v13_mathb.ckpt" || ok=0
    scp -3 -q -o BatchMode=yes "$MG:$D/runs/norooks_d13_seed2.mathg.ckpt" "$MD:$D/runs/v13_mathg.ckpt" || ok=0
    if [ $ok = 1 ]; then
      out=$(ssh -o BatchMode=yes "$MD" "cd ~/$D && ./perft --variant norooks --depth 13 --split 6 --seed 2 --merge \
        runs/v13_studio.ckpt,runs/v13_m5pro.ckpt,runs/norooks_d13_seed2.mathd.ckpt,runs/v13_mathb.ckpt,runs/v13_mathg.ckpt 2>&1 | tee runs/merge_a13v.last" < /dev/null)
      if echo "$out" | grep -q CONFLICT; then
        echo "CONFLICT in verification merge:"; echo "$out"; exit 1
      fi
      if echo "$out" | grep -q "^perft(13) = "; then
        ssh -o BatchMode=yes "$MD" "cd ~/$D && { echo '# merged: Mac Studio (forward) + mathd, mathb, mathg, M5 Pro (from the end); seed 2, split 6'; cat runs/merge_a13v.last; } >> runs/norooks_d13_seed2.log" < /dev/null
        ssh -o BatchMode=yes "$STUDIO" 'pkill -f "[v]erify_a13"; pkill -x perft' < /dev/null
        ssh -o BatchMode=yes "$M5" 'pkill -x perft' < /dev/null
        for h in "$MD" "$MB" "$MG"; do ssh -o BatchMode=yes "$h" "pkill -f '$RUN'" < /dev/null; done
        echo "$out" | tail -2
        exit 0
      fi
    fi
  fi
  sleep 300
done
