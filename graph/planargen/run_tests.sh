#!/bin/sh
# run_tests.sh -- test suite for lrplanar_sg.c (sparsegraph left-right
# planarity test, no size limit), with the original lrplanar.c and nauty's
# planarity.c as references.
#
#   usage: ./run_tests.sh [quick|full]      (default: full)
#
# Tiers:
#   1  exhaustive: every graph (connected or not) on 1..10 vertices, three
#      testers compared (nauty, sparse LR, bitmask LR)
#   2  random simple graphs near the planarity threshold, n = 11..64
#   3  random multigraphs with loops and parallel edges (nauty vs sparse LR)
#   4  random sparse graphs with thousands of vertices
#   5  known-answer families at 10^4..10^5 vertices, three-way where possible
#   6  huge known-answer graphs (10^6..10^7 vertices), sparse LR alone for the
#      high-degree families (nauty's tester is quadratic in vertex degree),
#      with nauty as reference otherwise; reports time and peak memory
#   7  regression: the enumeration of connected planar graphs with the sparse
#      tester as geng backend reproduces OEIS A003094 for n = 9, 10, 11
#   8  sanitizer build (AddressSanitizer + UBSan) on a sample of tiers 1-5
#   9  digraph input (digraph = TRUE): every orientation of every graph on
#      <= 6 vertices (directg), random digraphs (genrang -z), an exact count
#      comparison with McKay's reference pipeline "underlyingg | planarg",
#      and the huge families oriented inside the harness (-D; digraph6 is a
#      dense format, so huge digraphs cannot come from files)
#
# "quick" skips tiers 6-8, shortens the others and runs only the small
# part of tier 9 (about one minute).
# Exit status 0 iff every test passed.
set -u
MODE=${1:-full}
DIR=$(cd "$(dirname "$0")" && pwd)
NAUTY=${NAUTY:-$DIR/../nauty2_9_3}
TMP=${TMPDIR:-/tmp}/lrtests.$$
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT
cd "$DIR"
PASS=0; FAIL=0

check() {   # check NAME EXPECTED_RC  <- runs "$@" after those two
    name=$1; shift
    if "$@" > "$TMP/out" 2>&1; then
        PASS=$((PASS+1)); printf '  ok   %s\n' "$name"
        grep -h "^lrtest_sg:" "$TMP/out" | sed 's/^/         /'
    else
        FAIL=$((FAIL+1)); printf '  FAIL %s\n' "$name"; sed 's/^/         /' "$TMP/out" | head -12
    fi
}

# pipe helper: check NAME "producer command" "consumer command"
pcheck() { name=$1; prod=$2; cons=$3; check "$name" sh -c "$prod | $cons"; }

# family helper: FAM p|n args...  (three-way for small n, nauty+sg otherwise)
fam() { f=$1; e=$2; shift 2; pcheck "$f $* (expected $e)" "./biggraphs $f $* 2>/dev/null" "./lrtest_sg -e $e -q"; }
famN() { f=$1; e=$2; shift 2; pcheck "$f $* (expected $e, sg only)" "./biggraphs $f $* 2>/dev/null" "./lrtest_sg -N -B -e $e -q"; }

echo "== building =="
make -s lrtest_sg biggraphs geng_coplanar_lrsg lrtest 2>&1 | grep -i "error" && exit 1
[ "$MODE" = full ] && { make -s lrtest_sg_asan 2>&1 | grep -i "error" && exit 1; }

echo "== tier 1: exhaustive, all graphs on 1..10 vertices (three testers) =="
NMAX=9; [ "$MODE" = full ] && NMAX=10
for n in 1 2 3 4 5 6 7 8 9 10; do
    [ $n -gt $NMAX ] && break
    pcheck "all graphs n=$n" "$NAUTY/geng -q $n" "./lrtest_sg -q"
done

echo "== tier 2: random simple graphs near the planarity threshold =="
NUM=100000; [ "$MODE" = quick ] && NUM=20000
for n in 11 12 14 16 20 24 32; do
    pcheck "random n=$n e=$((2*n+3)) x$NUM" "$NAUTY/genrang -q -e$((2*n+3)) $n $NUM" "./lrtest_sg -q"
done
pcheck "random n=64 e=130 x$((NUM/5))" "$NAUTY/genrang -q -e130 64 $((NUM/5))" "./lrtest_sg -q"
pcheck "random sparse n=16 e=17 x$NUM (many disconnected)" "$NAUTY/genrang -q -e17 16 $NUM" "./lrtest_sg -q"

echo "== tier 3: random multigraphs with loops and parallel edges =="
for n in 6 8 10 14 20 40; do
    pcheck "3-regular multigraphs n=$n loops<=2 mult<=3 x$((NUM/2))" "$NAUTY/genrang -q -r3 -l2 -m3 $n $((NUM/2))" "./lrtest_sg -q"
done
pcheck "4-regular multigraphs n=12 loops<=1 mult<=2 x$((NUM/2))" "$NAUTY/genrang -q -r4 -l1 -m2 12 $((NUM/2))" "./lrtest_sg -q"

echo "== tier 4: random sparse graphs with thousands of vertices =="
for n in 100 1000 10000; do
    pcheck "random n=$n e=$((n+n/2)) x$((200000/n))" "$NAUTY/genrang -q -e$((n+n/2)) $n $((200000/n))" "./lrtest_sg -q"
done

echo "== tier 5: known-answer families, 10^4..10^5 vertices =="
fam grid p 300 400;        fam cylinder p 200 300;    fam ladder p 50000
fam path p 100000;         fam cycle p 100000;        fam star p 100000
fam wheel p 100000;        fam k2n p 100000;          fam tree p 100000 7
fam forest p 100000 3;     fam apollonian p 100000 1; fam apollodel p 100000 50000 2
fam multi p 50000 5
fam torus n 300 400;       fam mobius n 50000;        fam k3n n 100000
fam k5sub n 100005;        fam apollok5 n 100000 4;   fam apollok33 n 100000 9
fam k5far n 100000;        fam multink33 n 100000 11

if [ "$MODE" = full ]; then
    echo "== tier 6: huge graphs (10^6..10^7 vertices) =="
    fam grid p 1000 1000;          fam apollonian p 1000000 1;    fam apollodel p 2000000 1000000 3
    fam multi p 1000000 4;         fam path p 10000000;           fam forest p 2000000 5
    fam torus n 1000 1000;         fam mobius n 1000000;          fam apollok5 n 1000000 6
    fam apollok33 n 1000000 7;     fam k5far n 3000000;           fam k5sub n 2000005
    fam multink33 n 1000000 8
    famN star p 1000000;           famN k2n p 1000000;            famN k3n n 1000000
    famN wheel p 1000000;          famN grid p 2000 2000;         famN apollonian p 10000000 3
    famN apollok5 n 10000000 6
    echo "   time and peak memory of the sparse tester alone (includes reading the input):"
    for spec in "apollonian 10000000 3" "path 10000000" "grid 2000 2000" "star 1000000"; do
        ./biggraphs $spec > "$TMP/g.s6" 2>/dev/null
        /usr/bin/time -l ./lrtest_sg -N -B -q < "$TMP/g.s6" > "$TMP/t" 2>&1 || /usr/bin/time -v ./lrtest_sg -N -B -q < "$TMP/g.s6" > "$TMP/t" 2>&1
        printf '     %-24s wall %ss  peak RSS %s MB\n' "$spec" "$(grep real "$TMP/t" | awk '{print $1}')" \
            "$(grep -i "maximum resident" "$TMP/t" | awk '{v=$1; if (v > 1e7) v/=1048576; else v/=1024; printf "%.0f", v}')"
    done

    echo "== tier 7: enumeration regression (OEIS A003094 with the sparse tester as backend) =="
    for n in 9 10 11; do
        e=$((n*(n-1)/2))
        check "a($n) via geng_coplanar_lrsg" sh -c "./geng_coplanar_lrsg -u -v $n $((e-(3*n-6))):$((e-(n-1))) > '$TMP/g$n' 2>&1 && python3 verify.py $n coplanar '$TMP/g$n' | tail -1 | grep -q PASS"
    done

    echo "== tier 8: sanitizer build on a sample =="
    export ASAN_OPTIONS=detect_leaks=0
    pcheck "asan: all graphs n=8" "$NAUTY/geng -q 8" "./lrtest_sg_asan -q"
    pcheck "asan: 3-regular multigraphs n=14" "$NAUTY/genrang -q -r3 -l2 -m3 14 20000" "./lrtest_sg_asan -q"
    pcheck "asan: random n=16 e=30" "$NAUTY/genrang -q -e30 16 50000" "./lrtest_sg_asan -q"
    for spec in "grid p 100 150" "apollonian p 20000 1" "multi p 20000 2" "torus n 50 60" "apollok33 n 20000 3" \
                "multink33 n 20000 4" "star p 50000" "k5far n 20000" "forest p 30000 5" "k5sub n 20005"; do
        set -- $spec; f=$1; e=$2; shift 2
        pcheck "asan: $f $* (expected $e)" "./biggraphs $f $* 2>/dev/null" "./lrtest_sg_asan -e $e -q"
    done
fi

echo "== tier 9: digraph input =="
for n in 3 4 5 6; do
    pcheck "all orientations of all graphs n=$n (directg)" "$NAUTY/geng -q $n | $NAUTY/directg -q" "./lrtest_sg -q"
done
pcheck "random digraphs n=8 arcs=12 x$NUM" "$NAUTY/genrang -z -q -e12 8 $NUM" "./lrtest_sg -q"
pcheck "random digraphs n=14 arcs=28 x$NUM" "$NAUTY/genrang -z -q -e28 14 $NUM" "./lrtest_sg -q"
pcheck "random digraphs n=30 arcs=70 x$((NUM/5))" "$NAUTY/genrang -z -q -e70 30 $((NUM/5))" "./lrtest_sg -q"
pcheck "-D: all graphs n=8 oriented in the harness" "$NAUTY/geng -q 8" "./lrtest_sg -D -q"
pcheck "-D: 3-regular multigraphs n=14 oriented" "$NAUTY/genrang -q -r3 -l2 -m3 14 $((NUM/2))" "./lrtest_sg -D -q"
$NAUTY/genrang -z -q -e20 12 $NUM > "$TMP/dig.d6"
c1=$(./lrtest_sg -N -B < "$TMP/dig.d6" 2>&1 | sed -n 's/.*graphs (digraphs).*, \([0-9]*\) planar.*/\1/p')
c2=$($NAUTY/underlyingg -q "$TMP/dig.d6" | $NAUTY/planarg -u 2>&1 | sed -n 's/.* \([0-9]*\) graphs planar.*/\1/p')
if [ -n "$c1" ] && [ "$c1" = "$c2" ]; then PASS=$((PASS+1)); printf '  ok   planar count of %s random digraphs equals underlyingg | planarg (%s)\n' "$NUM" "$c1"
else FAIL=$((FAIL+1)); printf '  FAIL planar count vs underlyingg | planarg: %s vs %s\n' "$c1" "$c2"; fi
if [ "$MODE" = full ]; then
    famD() { f=$1; e=$2; shift 2; pcheck "-D $f $* (expected $e)" "./biggraphs $f $* 2>/dev/null" "./lrtest_sg -D -e $e -q"; }
    famDN() { f=$1; e=$2; shift 2; pcheck "-D $f $* (expected $e, sg only)" "./biggraphs $f $* 2>/dev/null" "./lrtest_sg -D -N -B -e $e -q"; }
    famD grid p 1000 1000;      famD apollonian p 1000000 1;   famD multi p 500000 3
    famD path p 3000000;        famD torus n 700 700;          famD apollok33 n 1000000 2
    famD k5far n 1000000;       famD multink33 n 1000000 8;    famDN star p 1000000
    famDN k3n n 1000000
fi

echo "== summary: $PASS passed, $FAIL failed =="
[ $FAIL -eq 0 ]
