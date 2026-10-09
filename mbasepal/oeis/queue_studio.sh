#!/bin/sh
# Studio work queue: waits for the running A216 batch, then continues.
cd "$(dirname "$0")" || exit 1
while pgrep -f "batch.sh 22" > /dev/null; do sleep 30; done
cd multipal && mkdir -p runs && \
  ./multipal -k 4 -u 8400090327552001 -t 24 -r 10 > runs/k4_to_8.4e15.txt 2> runs/k4_to_8.4e15.err
cd ../A216 && ./batch.sh 22 16:60 17:60
echo "queue_studio finished" >> ../queue_studio.log
