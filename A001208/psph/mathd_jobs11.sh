#!/bin/sh
# mathd phase 2, v2 (Oct 3 07:50 PDT): A001209(303) = n(303,4), target 71148328 (formula value 71148327 + 1), split depth 3
# (4158140 work items; generation costs ~8.6 min per psph16 invocation, so chunks are 200000 items = 21 chunks of ~1-2 h).
# Depth 2 (301 items in chunks of 10) was abandoned: at most 10 of 96 threads could be busy.
cd ~/A0001208/psph
echo "=== JOB2b n(303,4) depth 3, 4158140 items in chunks of 200000 start $(date)" >> jobs.log
./run_range.sh h303k4b 303 4 71148328 3 0 4158140 200000 96 ./psph16 >> jobs.log 2>&1
echo "=== JOB2b done $(date)" >> jobs.log
echo "=== ALL JOBS DONE $(date)" >> jobs.log
