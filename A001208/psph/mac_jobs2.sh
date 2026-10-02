#!/bin/sh
# Mac: after the current range [10700,11003), run [10000,10500) and the reduced [4200,4700); target 5845.
cd /Users/Jeff.Sponaugle/src/math/A001208/psph
while pgrep -f "run_range.sh mac_h9k8c" > /dev/null; do sleep 60; done
echo "=== MAC JOB1c done $(date)" >> mac_jobs.log
echo "=== MAC JOB1b items [10000,10500) start $(date)" >> mac_jobs.log
./run_range.sh mac_h9k8b 9 8 5946 4 10000 10500 250 10 ./psph >> mac_jobs.log 2>&1
echo "=== MAC JOB1d items [4200,4700) start $(date)" >> mac_jobs.log
./run_range.sh mac_h9k8d 9 8 5946 4 4200 4700 250 10 ./psph >> mac_jobs.log 2>&1
echo "=== MAC JOB1 (all Mac ranges) done $(date)" >> mac_jobs.log
