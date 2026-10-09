#!/bin/sh
# work.sh MACHINE THREADS -- this machine's whole workload: its own queue file, then the
# shared base-2 pool.  Both stages skip finished jobs and resume interrupted ones, so a
# (re)start is always safe:
#   pkill -f '^sh work.sh'; pkill -f '^sh run_'; sleep 1; pkill -f 'concat search'
#   setsid nohup sh work.sh MACHINE THREADS > /dev/null 2>&1 < /dev/null &
cd "$(dirname "$0")" || exit 1
sh run_jobs.sh "jobs_$1_all.txt" "$2"
sh run_pool.sh jobs_base2_1e12.txt "$2"
