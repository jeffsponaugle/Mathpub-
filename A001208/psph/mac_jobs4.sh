#!/bin/sh
# Mac: after [10000,10250), run [4200,4450) only ([4450,4700) -> mathg, [10250,10500) -> mathb).
cd /Users/Jeff.Sponaugle/src/math/A001208/psph
while pgrep -f "run_range.sh mac_h9k8b" > /dev/null; do sleep 60; done
echo "=== MAC JOB1b done $(date)" >> mac_jobs.log
echo "=== MAC JOB1d items [4200,4450) start $(date)" >> mac_jobs.log
./run_range.sh mac_h9k8d 9 8 6083 4 4200 4450 250 10 ./psph >> mac_jobs.log 2>&1
echo "=== MAC JOB1 (all Mac ranges) done $(date)" >> mac_jobs.log
