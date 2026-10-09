#!/bin/bash
# Regression suite: brute-force self-tests (v1 and v2, CPU and GPU), CPU/GPU engine comparison,
# reproduction of the A259385 b-file and of the known A393014 terms.
set -u
cd "$(dirname "$0")/.."
fail=0
run() { echo "== $*"; "$@" || { echo "FAILED: $*"; fail=1; }; }
summ() { awk '/configurations|OK$|MISMATCH/ {n++} /MISMATCH/ {bad++} END {printf "   %d lines, %d mismatches\n", n, bad+0; exit bad>0}'; }

for b in 2,9 2,10 3,10 8,10 9,10 2,3 10,15 6,10; do
  echo "== selftest (v1, cpu) bases $b"; ./dbpal-cpu selftest -b $b --Lmin 5 --Lmax 13 --engine cpu 2>&1 | summ || fail=1
done
for b in 2,9 3,10 9,10 7,10 2,3 10,11; do
  echo "== selftest2 (v2, cpu) bases $b"; ./dbpal-cpu selftest2 -b $b --Lmin 7 --Lmax 15 --engine cpu 2>&1 | summ || fail=1
  echo "== selftest2 (v2, gpu) bases $b"; ./dbpal selftest2 -b $b --Lmin 7 --Lmax 15 2>&1 | summ || fail=1
done
run bash tests/cmp_engines.sh

echo "== multi-limb arithmetic (small limbs, dbpal and dbpal4)"
python3 tests/limbs.py || fail=1

echo "== OEIS families (independent check)"
python3 tests/oeis_families.py || fail=1

echo "== A259385 b-file (bases 2,9 below 9^33) with v1 and v2"
T=$(mktemp -d)
python3 - "$T" <<'PY'
import json, sys, os
sys.path.insert(0, "tools")
import oeis
open(os.path.join(sys.argv[1], "bfile.txt"), "w").write("\n".join(map(str, oeis.all_terms(oeis.entry("A259385")))) + "\n")
PY
for algo in v1 v2; do
  ./dbpal search -b 2,9 --Lmin 1 --Lmax 33 --algo $algo 2>/dev/null | cut -f1 | sort > $T/$algo.txt
  if diff -q <(sort $T/bfile.txt) $T/$algo.txt >/dev/null; then echo "   $algo: identical (69 terms)"; else echo "   $algo: DIFFERENT"; fail=1; fi
done
rm -rf "$T"
[ $fail = 0 ] && echo "ALL TESTS PASSED" || echo "SOME TESTS FAILED"
exit $fail
