#!/bin/sh
# mathd overnight v2 (Oct 3 01:10 PDT): let the h27k6d driver finish its current chunk, then the remaining A001211(27)
# chunks at TGT 186943 (13/14 remainder, 27-29, then 25 taken from mathb), then A001209(303) in 10-item chunks.
cd ~/A0001208/psph
while pgrep -f "run_range.sh h27k6d" > /dev/null; do sleep 60; done
echo "=== JOB3d done (old driver exited) $(date)" >> jobs.log
echo "=== JOB3i n(27,6) items [1300000,1500000) remainder, TGT 186943 start $(date)" >> jobs.log
./run_range.sh h27k6i 27 6 186943 4 1300000 1500000 100000 96 ./psph >> jobs.log 2>&1
echo "=== JOB3i done $(date)" >> jobs.log
echo "=== JOB3e n(27,6) items [2700000,2964437) (chunks 27-29), TGT 186943 start $(date)" >> jobs.log
./run_range.sh h27k6e 27 6 186943 4 2700000 2964437 100000 96 ./psph >> jobs.log 2>&1
echo "=== JOB3e done $(date)" >> jobs.log
echo "=== JOB3g n(27,6) items [2500000,2600000) (chunk 25, from mathb), TGT 186943 start $(date)" >> jobs.log
./run_range.sh h27k6g 27 6 186943 4 2500000 2600000 100000 96 ./psph >> jobs.log 2>&1
echo "=== JOB3g done $(date)" >> jobs.log
echo "=== JOB2 n(303,4) restart in chunks of 10 items $(date)" >> jobs.log
./run_range.sh h303k4 303 4 71148328 2 0 301 10 96 ./psph16 >> jobs.log 2>&1
echo "=== JOB2 done $(date)" >> jobs.log
echo "=== ALL JOBS DONE $(date)" >> jobs.log
