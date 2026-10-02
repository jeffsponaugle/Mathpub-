#!/bin/bash
# stop_at_chunk.sh LOG X PID - SIGTERM process PID (an a350826_cuda count that writes chunk log LOG,
# relative to ~/A350826) once LOG contains chunk X+1.  Another run covers chunks X.. ; chunks X and
# X+1 then end up in both logs, as a cross-check.  The stopped run saves its checkpoint.
#
#   nohup bash stop_at_chunk.sh runs/a20b.chunks 91842 4149926 > /dev/null 2>&1 < /dev/null &
LOG=$1
X=$2
PID=$3
cd ~/A350826 || exit 1
echo "$(date) stop_at_chunk: will stop PID $PID after chunk $((X + 1)) of $LOG" >> runs/queue.log
while kill -0 "$PID" 2> /dev/null; do
    last=$(tail -1 "$LOG" | awk '{print $2}')
    if [ -n "$last" ] && [ "$last" -ge $((X + 1)) ]; then
        kill -TERM "$PID"
        echo "$(date) stop_at_chunk: sent SIGTERM to $PID, last logged chunk $last" >> runs/queue.log
        exit 0
    fi
    sleep 20
done
echo "$(date) stop_at_chunk: PID $PID exited before chunk $((X + 1))" >> runs/queue.log
