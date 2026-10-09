#!/bin/sh
# restart.sh MACHINE THREADS -- stop this machine's workload and start it again with
# work.sh; finished jobs are skipped and the interrupted one resumes.
cd "$(dirname "$0")" || exit 1
pkill -f '^sh work.sh'; pkill -f '^sh -c while kill -0'; pkill -f '^sh run_'; sleep 1
pkill -f 'concat search'; sleep 2
(setsid nohup sh work.sh "$1" "$2" > /dev/null 2>&1 < /dev/null &)
sleep 3
tail -2 runs/jobs.log
