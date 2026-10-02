#!/bin/sh
# Mac: tail of the A001211(27) item list (mathg's chunks 20-29), target 176381.
cd /Users/Jeff.Sponaugle/src/math/A001208/psph
echo "=== MAC JOB3 n(27,6) items [2000000,2964437) start $(date)" >> mac_jobs.log
./run_range.sh mac_h27k6 27 6 176381 4 2000000 2964437 100000 10 ./psph >> mac_jobs.log 2>&1
echo "=== MAC JOB3 done $(date)" >> mac_jobs.log
