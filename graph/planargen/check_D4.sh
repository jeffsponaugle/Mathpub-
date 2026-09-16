#!/bin/sh
# check_D4.sh -- independent n=14 cross-check on the subclass "maximum degree <= 4":
#   naive:  geng -c -D4 14 r/4 | planarg -q | countg --e      (nauty tools only)
#   ours:   geng_coplanar_lr -u -v -d9 14 55:78 r/4           (complement min degree >= 9)
# Runs 3 jobs at a time at lowest priority; writes RESULT_D4.txt when done.
DIR=$(cd "$(dirname "$0")" && pwd)
NAUTY="$DIR/../nauty2_9_3"
OUT="$DIR/../runs/n14_check_D4"
mkdir -p "$OUT"
cd "$OUT"
{
  for r in 0 1 2 3; do
    echo "naive $r"
    echo "ours $r"
  done
} | xargs -P 3 -L 1 sh -c '
  kind=$0; r=$1
  if [ "$kind" = naive ]; then
    [ -f naive_$r.txt ] && grep -q altogether naive_$r.txt && exit 0
    nice -n 19 "'"$NAUTY"'/geng" -c -D4 -q 14 $r/4 | nice -n 19 "'"$NAUTY"'/planarg" -q | nice -n 19 "'"$NAUTY"'/countg" --e > naive_$r.tmp 2>&1 && mv naive_$r.tmp naive_$r.txt
  else
    [ -f ours_$r.txt ] && grep -q ">Z" ours_$r.txt && exit 0
    nice -n 19 "'"$DIR"'/geng_coplanar_lr" -u -v -d9 14 55:78 $r/4 > ours_$r.tmp 2>&1 && mv ours_$r.tmp ours_$r.txt
  fi
'
{ echo "finished $(date)"; python3 "$DIR/compare_naive.py" 14 4 "$OUT"; } > "$OUT/RESULT_D4.txt" 2>&1
