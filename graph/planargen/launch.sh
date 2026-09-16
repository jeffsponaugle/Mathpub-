#!/bin/sh
# launch.sh -- detached run for a(n): MOD pieces, P processes at a time, low
# priority, survives closing the terminal.  Output in ../runs/n<n>, with frozen
# copies of the binary and driver.  When the full piece range is run on this
# host, finish.sh is started too and writes ../runs/n<n>/RESULT.txt (aggregate,
# verification against OEIS, Euler transform, b-file lines) at the end.
#
#   usage: ./launch.sh n MOD P [FIRST LAST]
#   e.g.   ./launch.sh 14 1000 8            # one machine, 8 processes
#          ./launch.sh 15 4000 20 0 1999    # host A   (then host B: 2000 3999)
#
# Progress: ls ../runs/n<n>/part_*_of_<MOD>.txt | wc -l
# Stop:     pkill -f "geng_coplanar_lr -u -v <n> "   (pieces are restartable:
#           rerunning the same launch skips completed pieces)
set -eu
n=$1; MOD=$2; P=$3; FIRST=${4:-0}; LAST=${5:-$((MOD-1))}
DIR=$(cd "$(dirname "$0")" && pwd)
OUT="$DIR/../runs/n$n"
mkdir -p "$OUT"
cp "$DIR/geng_coplanar_lr" "$DIR/run_split.sh" "$OUT/"
( shasum -a 256 "$OUT/geng_coplanar_lr" 2>/dev/null || sha256sum "$OUT/geng_coplanar_lr" ) > "$OUT/geng_coplanar_lr.sha256"
cd "$OUT"
nohup env BIN=geng_coplanar_lr nice -n 15 ./run_split.sh "$n" "$MOD" "$P" "$OUT" "$FIRST" "$LAST" > "$OUT/run.log" 2>&1 &
echo "started pid $! at $(date): n=$n, pieces $FIRST..$LAST of $MOD, $P at a time; log $OUT/run.log"
if [ "$FIRST" -eq 0 ] && [ "$LAST" -eq $((MOD-1)) ]; then
    nohup "$DIR/finish.sh" "$n" "$MOD" "$OUT" > /dev/null 2>&1 &
    echo "result will appear in $OUT/RESULT.txt"
else
    echo "partial range: when every host is done, copy all part files into one directory and run"
    echo "  python3 $DIR/aggregate.py $n $MOD <dir> && python3 $DIR/euler.py <dir>/rows_with_n$n.json $n && python3 $DIR/bfile.py <dir>/rows_with_n$n.json $n"
fi
