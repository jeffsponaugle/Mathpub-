#!/bin/sh
# finish.sh n MOD outdir -- wait until run.log reports completion, then
# aggregate, verify against OEIS, derive A005470/A039735 and b-file lines.
DIR=$(cd "$(dirname "$0")" && pwd)
n=$1; MOD=$2; OUT=$3
while ! grep -q "^done:" "$OUT/run.log" 2>/dev/null; do sleep 30; done
{
  echo "n=$n run finished at $(date)"
  tail -1 "$OUT/run.log"
  echo "=== aggregate + verification ==="
  python3 "$DIR/aggregate.py" "$n" "$MOD" "$OUT"
  echo "=== Euler transform (A005470; A039735 rows) ==="
  python3 "$DIR/euler.py" "$OUT/rows_with_n$n.json" "$n" | grep -E "^A003094|^A005470|^A039735 row ($((n-1))|$n) "
  echo "=== b-file lines ==="
  python3 "$DIR/bfile.py" "$OUT/rows_with_n$n.json" "$n"
} > "$OUT/RESULT.txt" 2>&1
