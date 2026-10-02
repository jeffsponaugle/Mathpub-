#!/bin/sh
# Resumable, chunked driver for psph over a range of work items (for splitting a search across machines).
#
#   run_range.sh TAG H K TGT D LO HI CHUNK THREADS [BINARY]
#
#   TAG      name for the log directory (e.g. h9k8_box1)
#   H K TGT  the search (TGT = best known lower bound + 1)
#   D        split depth (must be the same on every machine so that item indices agree)
#   LO HI    this machine's item range [LO,HI)  (item counts: run  psph -h H -k K -t TGT -d D -i 0:0)
#   CHUNK    items per psph invocation (a finished chunk is marked by a .done file -> re-run to resume)
#   THREADS  worker threads (number of cores)
#   BINARY   psph (h <= 254) or psph16 (h >= 255); default ./psph
#
# Output: logs/TAG/progress.log, logs/TAG/solutions.log (SOLUTION lines), per-chunk .out/.err/.done.
# Result of the whole search = max n_h over ALL machines' solutions.log once every chunk everywhere is done.
TAG=$1; H=$2; K=$3; TGT=$4; D=$5; LO=$6; HI=$7; C=$8; J=$9; P=${10:-./psph}
[ -z "$J" ] && { echo "usage: run_range.sh TAG H K TGT D LO HI CHUNK THREADS [BINARY]"; exit 1; }
L=logs/$TAG; mkdir -p $L
i=0
while :; do
  lo=$((LO + i*C)); [ $lo -ge $HI ] && break
  hi=$((lo + C)); [ $hi -gt $HI ] && hi=$HI
  if [ -f $L/chunk_$i.done ]; then i=$((i+1)); continue; fi
  echo "=== chunk $i items [$lo,$hi) h=$H k=$K TGT=$TGT start $(date)" >> $L/progress.log
  $P -h $H -k $K -t $TGT -j $J -p 300 -d $D -i $lo:$hi > $L/chunk_$i.out 2> $L/chunk_$i.err
  rc=$?
  grep SOLUTION $L/chunk_$i.out >> $L/solutions.log
  grep -q "worker CPU split" $L/chunk_$i.err && touch $L/chunk_$i.done
  echo "=== chunk $i done rc=$rc $(date): $(grep -c SOLUTION $L/chunk_$i.out) solutions; $(grep -m1 'leaf candidates' $L/chunk_$i.err | cut -c1-110)" >> $L/progress.log
  i=$((i+1))
done
echo "=== RANGE [$LO,$HI) COMPLETE $(date)" >> $L/progress.log
sort -u $L/solutions.log 2>/dev/null
