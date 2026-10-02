#!/bin/sh
# Stop every postage-stamp job on this box and report what is complete (resumable state = .done markers).
cd "$(dirname "$0")"
for p in $(pgrep -f "jobs[0-9]*\.sh"); do [ "$p" != "$$" ] && kill $p 2>/dev/null; done
pkill -f "run_range.sh" 2>/dev/null; sleep 1
pkill -x psph 2>/dev/null; pkill -x psph16 2>/dev/null; pkill -x psph_gpu 2>/dev/null; sleep 2
echo "=== STOPPED ALL by stop_all.sh $(date)" >> jobs.log
echo "remaining processes: $(pgrep -fa 'psph|run_range' | grep -v stop_all | wc -l)"
for t in logs/*/; do n=$(ls $t/chunk_*.done 2>/dev/null | wc -l); echo "$t: $n chunks done; solutions: $(sort -u $t/solutions.log 2>/dev/null | wc -l)"; done
