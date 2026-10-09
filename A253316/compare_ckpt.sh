#!/bin/sh
# Compare the per-unit counts of two A253316 checkpoint files and recompute
# the total (sum of weight * count) of each with exact integer arithmetic (bc).
#   usage: ./compare_ckpt.sh run1.ckpt run2.ckpt
set -e
[ $# -eq 2 ] || { echo "usage: $0 A.ckpt B.ckpt" >&2; exit 2; }
for f in "$1" "$2"; do
    printf '%s: %s units, total = ' "$f" "$(grep -c '^U' "$f")"
    { printf '0'; awk '/^U/ { printf "+%s*%s", $4, $5 }' "$f"; echo; } | BC_LINE_LENGTH=0 bc
done
awk 'NR == FNR { if (/^U/) a[$2 " " $3] = $5; next }
     /^U/ { k = $2 " " $3; if (k in a) { common++; if (a[k] != $5) { bad++; print "MISMATCH seam " k ": " a[k] " vs " $5 } } }
     END { printf "%d units in both files, %d mismatches\n", common, bad + 0; exit(bad > 0) }' "$1" "$2"
