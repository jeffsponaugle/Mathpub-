#!/bin/sh
# launch_n14.sh -- detached full run for a(14): 1000 pieces, P processes,
# low priority, survives the terminal.  Output goes to ../runs/n14, which
# also holds frozen copies of the binary and driver used.
#   usage: ./launch_n14.sh [P]         (default P = 8)
# Progress:  ls ../runs/n14/part_*_of_1000.txt | wc -l
# Stop:      pkill -f 'geng_coplanar_lr -u -v 14'      (pieces are restartable:
#            rerun this script and completed pieces are skipped)
# Finish:    python3 aggregate.py 14 1000 ../runs/n14
#            python3 euler.py ../runs/n14/rows_with_n14.json 14
#            python3 bfile.py ../runs/n14/rows_with_n14.json 14
set -eu
P=${1:-8}
DIR=$(cd "$(dirname "$0")" && pwd)
OUT="$DIR/../runs/n14"
mkdir -p "$OUT"
cp "$DIR/geng_coplanar_lr" "$DIR/run_split.sh" "$OUT/"
shasum -a 256 "$OUT/geng_coplanar_lr" > "$OUT/geng_coplanar_lr.sha256"
cd "$OUT"
nohup env BIN=geng_coplanar_lr nice -n 15 ./run_split.sh 14 1000 "$P" "$OUT" > "$OUT/run.log" 2>&1 &
echo "started pid $! at $(date); log: $OUT/run.log"
# watcher: aggregates, verifies and derives the other sequences when all pieces are done
nohup "$DIR/finish_n14.sh" > /dev/null 2>&1 &
echo "result will appear in $OUT/RESULT.txt"
