#!/bin/sh
# recompute the normalizer of a Z-class from scratch with CARAT's Normalizer
# program, then recount with zverify.  usage: renorm.sh srcdir dstdir file...
src=$1; dst=$2; shift 2
for f in "$@"; do
  "$CARAT_DIR/bin/Normalizer" "$src/$f" > "$dst/$f" 2>/dev/null || { echo "NORMFAIL $f" >&2; continue; }
  (cd "$dst" && "$CARAT_DIR/../zverify" "$f")
done
