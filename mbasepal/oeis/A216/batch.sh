#!/bin/sh
# usage: ./batch.sh JOBS L:DMAX [L:DMAX ...]   (runs each length in turn, resumable)
cd "$(dirname "$0")" || exit 1
J=$1; shift
for spec in "$@"; do
    L=${spec%%:*}; D=${spec##*:}
    ./run_gap.py -L "$L" --dmax "$D" -j "$J" >> "run_L$L.log" 2>&1
done
echo "batch finished: $*" >> batch.log
