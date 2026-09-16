#!/bin/sh
# check_backend.sh -- recompute selected production pieces of the n=14 run with
# the nauty-planarity backend (geng_coplanar) and diff against the main run's
# output for the same pieces.  Writes ../runs/n14_check_backend/RESULT_backend.txt
DIR=$(cd "$(dirname "$0")" && pwd)
MAIN="$DIR/../runs/n14"
OUT="$DIR/../runs/n14_check_backend"
PIECES="100 250 333 421 500 587 666 713 800 876 950 999"
mkdir -p "$OUT"; cd "$OUT"
for r in $PIECES; do echo $r; done | xargs -P 2 -I{} sh -c '
  [ -f nauty_$0.txt ] && grep -q ">Z" nauty_$0.txt && exit 0
  nice -n 19 "'"$DIR"'/geng_coplanar" -u -v 14 55:78 $0/1000 > nauty_$0.tmp 2>&1 && mv nauty_$0.tmp nauty_$0.txt' {}
{
  echo "nauty-backend recomputation finished $(date)"
  fail=0
  for r in $PIECES; do
    while [ ! -f "$MAIN/part_${r}_of_1000.txt" ]; do sleep 120; done
    a=$(grep "^>C [0-9]* graphs with\|^>Z" nauty_$r.txt | sed 's/ in .*//')
    b=$(grep "^>C [0-9]* graphs with\|^>Z" "$MAIN/part_${r}_of_1000.txt" | sed 's/ in .*//')
    if [ "$a" = "$b" ]; then echo "piece $r/1000: identical ($(grep '>Z' nauty_$r.txt | awk '{print $2}') graphs)"
    else echo "piece $r/1000: DIFFERENT"; fail=1; fi
  done
  echo "RESULT: $([ $fail = 0 ] && echo PASS || echo FAIL)"
} > "$OUT/RESULT_backend.txt" 2>&1
