#!/bin/bash
# queue_a20.sh - on atom1: wait for the a(19) count to finish, then count a(20).
#
# Bins: [1e19, 2^64) reproduces pi_6(2^64) - pi_6(10^19) = 48629687343 - 28722086297
# = 19907601046 (Desfontaines); the k*10^19 and 2^65, 2^66 boundaries cost nothing and
# give pi_6 at those points for later comparisons.  Resume after an interruption by
# rerunning the a(20) command below (same arguments, same state file).
#
#   nohup bash queue_a20.sh > /dev/null 2>&1 < /dev/null &
cd ~/A350826 || exit 1
log=runs/queue.log
BINS=2^64,2e19,3e19,2^65,4e19,5e19,6e19,7e19,2^66,8e19,9e19
echo "$(date) queue: waiting for the a(19) run" >> $log
while pgrep -x a350826_cuda > /dev/null; do sleep 60; done
if ! grep -q "^RESULT" runs/a19.out; then
    echo "$(date) queue: a(19) did not end with a RESULT line; not starting a(20)" >> $log
    exit 1
fi
echo "$(date) queue: a(19) done: $(grep '^RESULT \[' runs/a19.out)" >> $log
echo "$(date) queue: starting a(20)" >> $log
./a350826_cuda count 1e19 1e20 -b $BINS -S runs/a20.state -L runs/a20.chunks -i 300 \
    >> runs/a20.out 2>> runs/a20.err < /dev/null
echo "$(date) queue: a(20) exited with status $?" >> $log
