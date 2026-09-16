#!/bin/sh
# run_split.sh -- run geng_coplanar for n vertices split into MOD pieces,
# P at a time, writing one output file per piece.  Pieces are independent, so
# the same script can be run on several hosts with disjoint piece ranges
# (FIRST..LAST) and the output directories merged before aggregation.
#
#   usage: [BIN=geng_coplanar_lr] run_split.sh n MOD P outdir [FIRST [LAST]]
#
# BIN selects the generator binary (default geng_coplanar_lr, the fast
# left-right-planarity build; geng_coplanar uses nauty's planarity tester).
#
# Each piece writes outdir/part_<res>_of_<MOD>.txt; a piece is complete when
# its file contains a ">Z" line.  Rerunning skips completed pieces, so an
# interrupted run can simply be restarted.
set -eu
n=$1; MOD=$2; P=$3; OUT=$4; FIRST=${5:-0}; LAST=${6:-$((MOD-1))}
E=$((n*(n-1)/2)); LO=$((E-(3*n-6))); HI=$((E-(n-1)))
DIR=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-geng_coplanar_lr}
mkdir -p "$OUT"
echo "$BIN -u -v $n $LO:$HI res/$MOD  for res=$FIRST..$LAST, $P at a time -> $OUT" >&2
seq "$FIRST" "$LAST" | xargs -P "$P" -I{} sh -c '
  f="$1/part_$2_of_$3.txt"
  if [ -f "$f" ] && grep -q "^>Z" "$f"; then exit 0; fi
  "$4/$8" -u -v "$5" "$6:$7" "$2/$3" > "$f.tmp" 2>&1 && mv "$f.tmp" "$f"
' _ "$OUT" {} "$MOD" "$DIR" "$n" "$LO" "$HI" "$BIN"
echo "done: $(grep -l "^>Z" "$OUT"/part_*_of_"$MOD".txt 2>/dev/null | wc -l | tr -d " ") of $((LAST-FIRST+1)) pieces complete" >&2
