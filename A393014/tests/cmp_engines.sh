#!/bin/bash
# Compare CPU and GPU engines on several base pairs (outputs must be identical).
S=${TMPDIR:-/tmp}/dbpal_cmp.$$
mkdir -p $S
status=0
while read -r b lo hi; do
  [ -z "$b" ] && continue
  ./dbpal-cpu search -b $b --Lmin $lo --Lmax $hi > $S/c.out 2>$S/c.err
  ./dbpal search -b $b --Lmin $lo --Lmax $hi > $S/g.out 2>$S/g.err
  nc=$(wc -l < $S/c.out); ng=$(wc -l < $S/g.out)
  tc=$(grep -o "search [0-9.]*s" $S/c.err | awk '{s+=$2} END {printf "%.1f", s}')
  tg=$(grep -o "search [0-9.]*s" $S/g.err | awk '{s+=$2} END {printf "%.1f", s}')
  if diff -q <(cut -f1 $S/c.out) <(cut -f1 $S/g.out) >/dev/null; then r=SAME; else r=DIFF; status=1; fi
  ov=$(grep -c OVERFLOW $S/g.err $S/c.err | awk -F: '{s+=$2} END {print s}')
  echo "bases $b L=$lo..$hi: cpu $nc hits (${tc}s)  gpu $ng hits (${tg}s)  $r  overflow:$ov"
done <<'LIST'
2,10 1 28
3,10 1 26
8,10 1 26
9,10 1 24
2,3 1 45
7,10 1 24
10,15 1 22
6,10 1 24
10,11 1 24
2,5 1 45
4,10 1 26
10,16 1 24
10,12 1 24
2,9 1 33
LIST
rm -rf $S
exit $status
