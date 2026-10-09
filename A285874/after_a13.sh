#!/bin/bash
# On an x86 server: once this machine's part of the distributed A285874 a(13) run has
# finished (its share and phase 2), runs the run-1 (seed 1) job in the given job file.
# Usage: HASH=<MB> setsid nohup ./after_a13.sh JOBFILE > runs/after_a13.out 2>&1 < /dev/null &
cd ~/src/math/A285874 || exit 1
while pgrep -f "[p]hase2.sh|[p]erft3* --variant norooks --depth 13" > /dev/null; do sleep 60; done
SEEDS=1 numactl --interleave=all ./extend.sh "$1"
