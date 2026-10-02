#!/bin/bash
# stop_a20.sh X PID - on atom1: SIGTERM the a(20) run (process PID) once its chunk log contains
# chunk X+1.  atom2 counts chunks [X, 92569) (a20b.*), so chunks X and X+1 end up computed by both
# GPUs, which gives a direct cross-check; the merge uses atom1's lines below X only.
#
#   nohup bash stop_a20.sh 53416 1003403 > /dev/null 2>&1 < /dev/null &
X=$1
PID=$2
cd ~/A350826 || exit 1
echo "$(date) stop_a20: will stop PID $PID after chunk $((X + 1)) (split at chunk $X)" >> runs/queue.log
while kill -0 "$PID" 2> /dev/null; do
    last=$(tail -1 runs/a20.chunks | awk '{print $2}')
    if [ -n "$last" ] && [ "$last" -ge $((X + 1)) ]; then
        kill -TERM "$PID"
        echo "$(date) stop_a20: sent SIGTERM, last logged chunk $last" >> runs/queue.log
        exit 0
    fi
    sleep 20
done
echo "$(date) stop_a20: PID $PID exited before reaching chunk $((X + 1))" >> runs/queue.log
