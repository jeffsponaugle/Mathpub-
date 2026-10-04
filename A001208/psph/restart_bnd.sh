#!/bin/sh
# restart the boundary-window run at normal priority (resumable via .done markers)
cd ~/A001208/psph
for p in $(pgrep -f "sh [b]nd_h27k6.sh"); do kill $p; done
sleep 1
for p in $(pgrep -f "psph -h 27 -k 6 -t 176381 -j [4] -p"); do kill $p; done
sleep 1
sed -i 's/^  nice -n 19 \.\/psph/  .\/psph/' bnd_h27k6.sh
grep -n "psph -h 27" bnd_h27k6.sh
nohup sh bnd_h27k6.sh > bnd.log 2>&1 < /dev/null &
sleep 3
echo "running: $(pgrep -af '[b]nd_h27k6' | wc -l) script, $(pgrep -af 'psph -h 27 -k 6 -t 176381 -j [4] -p' | wc -l) psph"
echo "main chunk run still up: $(pgrep -af 'psph -h 27 -k 6 -t 176381 -j 3[2] ' | wc -l)"
