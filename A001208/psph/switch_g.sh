#!/bin/sh
cd ~/A001208/psph
for p in $(pgrep -f "sh [j]obs5.sh"); do kill $p; done
sleep 1
nohup sh jobs7.sh > /dev/null 2>&1 < /dev/null &
echo "bnd v1 status: $(ls logs/h27k6_bnd/*.done | wc -l)/29 done, $(tail -1 logs/h27k6_bnd/progress.log | cut -c1-60)"
if pgrep -f "sh [b]nd_h27k6.sh" > /dev/null; then echo "bnd v1 still running; bnd2 will start after it"; (while pgrep -f "sh [b]nd_h27k6.sh" > /dev/null; do sleep 30; done; nohup sh bnd2_h27k6.sh > bnd2.log 2>&1 < /dev/null) > /dev/null 2>&1 < /dev/null &
else nohup sh bnd2_h27k6.sh > bnd2.log 2>&1 < /dev/null & fi
sleep 2
echo "wrappers: $(pgrep -af '[j]obs[0-9]*.sh' | cut -c1-60 | tr '\n' ' ')"; echo "drivers: $(pgrep -af '[r]un_range' | cut -c1-90)"; echo "bnd: $(pgrep -af '[b]nd' | cut -c1-40 | tr '\n' ' ')"
echo "=== SWITCH2 $(date): after chunk 19, chunk 26 at TGT 186943 (h27k6f), then h91k5c; boundary windows v2 (half-width 24)" >> jobs.log
tail -1 logs/h27k6/chunk_19.err | cut -c1-120
