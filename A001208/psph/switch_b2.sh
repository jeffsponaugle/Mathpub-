#!/bin/sh
cd ~/A001208/psph
for p in $(pgrep -f "sh [j]obs7.sh"); do kill $p; done
sleep 1
nohup sh jobs8.sh > /dev/null 2>&1 < /dev/null &
sleep 2
echo "wrappers: $(pgrep -af '[j]obs[0-9]*.sh' | cut -c1-60 | tr '\n' ' ')"; echo "drivers: $(pgrep -af '[r]un_range' | cut -c1-90)"
echo "=== REBALANCE7 $(date): after chunk 24, chunks 25 and 29 (from mathd) at TGT 186943, then h91k5b" >> jobs.log
tail -1 logs/h27k6h/chunk_2.err | cut -c1-110
