#!/bin/sh
# A001211(27) chunk-boundary windows. The depth-4 work-item enumeration depends on the target through the
# element-wise lower bounds (2964437 items at TGT 176381, 2964430 at TGT 177915: the 7 extra prefixes are
# proven unable to reach 177915). Chunks run at different targets therefore use index spaces that differ by
# a shift of at most 7, so up to 7 items could fall through at a boundary where the lower chunk used TGT
# 176381 and the upper one 177915. Re-running items [B-8,B+8) at TGT 176381 (the larger enumeration) for
# every chunk boundary B = 100000*i covers every such item regardless of which side used which target.
cd ~/A001208/psph
L=logs/h27k6_bnd; mkdir -p $L
for i in $(seq 1 29); do
  B=$((i*100000)); lo=$((B-8)); hi=$((B+8))
  [ -f $L/w_$i.done ] && continue
  echo "=== window $i items [$lo,$hi) h=27 k=6 TGT=176381 start $(date)" >> $L/progress.log
  nice -n 19 ./psph -h 27 -k 6 -t 176381 -j 4 -p 300 -d 4 -i $lo:$hi > $L/w_$i.out 2> $L/w_$i.err
  grep SOLUTION $L/w_$i.out >> $L/solutions.log
  grep -q "worker CPU split" $L/w_$i.err && touch $L/w_$i.done
  echo "=== window $i done $(date): $(grep -c SOLUTION $L/w_$i.out) solutions; $(grep -m1 'processing work items' $L/w_$i.err)" >> $L/progress.log
done
echo "=== BOUNDARY WINDOWS COMPLETE $(date)" >> $L/progress.log
