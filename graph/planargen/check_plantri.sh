#!/bin/sh
# check_plantri.sh -- independent check of the dense end of A049334 rows.
#
# Every planar graph is a spanning subgraph of a simple triangulation, and a
# triangulation (3-connected) stays connected after deleting up to two edges.
# So the connected planar graphs with 3n-7 (resp. 3n-8) edges are exactly the
# distinct graphs obtained by deleting one (resp. two) edges from the
# triangulations on n vertices.  plantri (Brinkmann & McKay) generates the
# triangulations by a completely different method (face additions), so this
# is independent of the geng-based generator.
#
#   usage: check_plantri.sh n [2]        (2 = also do two-edge deletions)
# Needs plantri in $PLANTRI and the nauty tools deledgeg / shortg in $NAUTY.
set -eu
n=$1; two=${2:-}
DIR=$(cd "$(dirname "$0")" && pwd)
NAUTY=${NAUTY:-$DIR/../nauty2_9_3}
PLANTRI=${PLANTRI:-$DIR/../plantri/plantri}
TMP=${TMPDIR:-/tmp}
echo "n=$n: triangulations: $($PLANTRI -u $n 2>&1 | grep -o '[0-9]* triangulations')  (A000109)"
echo "T($n,$((3*n-7))) = $($PLANTRI -g $n 2>/dev/null | $NAUTY/deledgeg -q | $NAUTY/shortg -u -T$TMP 2>&1 | grep produced | awk '{print $2}')   (distinct triangulation-minus-one-edge graphs)"
if [ -n "$two" ]; then
  echo "T($n,$((3*n-8))) = $($PLANTRI -g $n 2>/dev/null | $NAUTY/deledgeg -q | $NAUTY/shortg -q -T$TMP | $NAUTY/deledgeg -q | $NAUTY/shortg -u -T$TMP 2>&1 | grep produced | awk '{print $2}')   (distinct triangulation-minus-two-edges graphs)"
fi
