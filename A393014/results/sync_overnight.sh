#!/bin/bash
# Every 10 minutes (for up to $1 hours, default 14): push this machine's log to the Studio and merge the Studio's
# results and log into results/.  Writes a timestamped line per round to results/sync.log.
cd "$(dirname "$0")/.."
end=$(( $(date +%s) + ${1:-14}*3600 ))
while [ "$(date +%s)" -lt "$end" ]; do
  {
    printf "%s " "$(date '+%H:%M:%S')"
    bash tools/push_log.sh jeffsponauglejbs@10.1.10.6 m4max 2>&1 | tr '\n' ' '
    python3 tools/merge_remote.py jeffsponauglejbs@10.1.10.6 2>&1 | tail -1
  } >> results/sync.log
  sleep 600
done
