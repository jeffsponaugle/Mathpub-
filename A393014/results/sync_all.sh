#!/bin/bash
# Hub sync for a multi-machine campaign (runs on the M4 Max).  Every $2 seconds (default 60), for up to $1 hours
# (default 72): merge each worker's results and log into results/, then push the merged log to every
# worker, so all machines see each other's finished lengths, parts and claims.  results/sync.log
# gets one line per host per round.
cd "$(dirname "$0")/.."
HOSTS="jeffsponauglejbs@10.1.10.6 jbs@10.1.30.36 jbs@10.1.30.37"
end=$(( $(date +%s) + ${1:-72}*3600 ))
while [ "$(date +%s)" -lt "$end" ]; do
  for h in $HOSTS; do
    printf "%s merge %s: %s\n" "$(date '+%H:%M:%S')" "$h" "$(python3 tools/merge_remote.py $h 2>&1 | tail -1)" >> results/sync.log
  done
  for h in $HOSTS; do
    printf "%s push  %s: %s\n" "$(date '+%H:%M:%S')" "$h" "$(bash tools/push_log.sh $h m4max 2>&1 | tail -1)" >> results/sync.log
  done
  sleep ${2:-60}
done
