#!/bin/sh
# Reproduce everything from scratch:
#   1. build CARAT (C) and the two tools (enantio: C on top of CARAT;
#      zverify: independent C++),
#   2. split every Q-class of CARAT's catalog (dims 2..6) into Z-classes with
#      normalizers (CARAT q2z, via `enantio -k`),
#   3. count space-group types, proper types and enantiomorphic pairs for every
#      Z-class with the independent checker zverify,
#   4. print the per-dimension table and compare with the OEIS.
# Takes about 3 minutes on a 14-core M4 Max.  J=<jobs> overrides parallelism.
set -e
cd "$(dirname "$0")/.."
ROOT=$PWD
export CARAT_DIR=$ROOT/carat
J=${J:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}

make -C carat -f Makefile.local -j"$J" all > /dev/null
[ -d carat/tables/qcatalog ] || (cd carat/tables && tar xzf qcatalog.tar.gz)
make enantio zverify > /dev/null
mkdir -p results work

# dimension 1 by hand (CARAT's q2z crashes on 1x1 groups):
# G = {1} and G = {+-1}, normalizer GL_1(Z) = {+-1} in both cases
rm -rf work/z1 && mkdir -p work/z1
printf '#g1 n1\n1\n1\n1\n-1\n1 = 1\n' > work/z1/trivial__Z1
printf '#g1\n1\n-1\n2^1 = 2\n' > work/z1/minus1__Z1
(cd work/z1 && "$ROOT/zverify" trivial__Z1 minus1__Z1 > "$ROOT/results/zverify_dim1.tsv")

for d in 2 3 4 5 6; do
    echo "dimension $d ..." >&2
    find carat/tables/qcatalog/dim$d -type f \( -name 'group.*' -o -name 'min.*' -o -name 'max.*' \) \
        | sort > results/qfiles_dim$d.txt
    rm -rf work/z$d && mkdir -p work/z$d
    tr '\n' '\0' < results/qfiles_dim$d.txt \
        | ENANTIO_DUMP_DIR=$ROOT/work/z$d xargs -0 -P "$J" -n 1 ./enantio -k > /dev/null 2> results/split_dim$d.err
    (cd work/z$d && ls | tr '\n' '\0' | xargs -0 -P "$J" -n 10 "$ROOT/zverify" > "$ROOT/results/zverify_dim$d.tsv")
done

python3 scripts/summarize.py
