#!/bin/sh
# finish_n14.sh -- wait for the n=14 run to complete, then aggregate, verify
# (against A049334's partial row 14 and A000109(14)), and derive A005470 /
# A039735 by Euler transform.  Results: ../runs/n14/RESULT.txt
DIR=$(cd "$(dirname "$0")" && pwd)
OUT="$DIR/../runs/n14"
while ! grep -q "^done:" "$OUT/run.log" 2>/dev/null; do sleep 60; done
{
  echo "n=14 run finished at $(date)"
  tail -1 "$OUT/run.log"
  echo "=== aggregate + verification ==="
  python3 "$DIR/aggregate.py" 14 1000 "$OUT"
  echo "=== Euler transform (A005470, A039735 rows) ==="
  python3 "$DIR/euler.py" "$OUT/rows_with_n14.json" 14 | grep -v "^A039735 row [0-9] \|^A039735 row 1[0-2] "
  echo "=== b-file lines ==="
  python3 "$DIR/bfile.py" "$OUT/rows_with_n14.json" 14
} > "$OUT/RESULT.txt" 2>&1
