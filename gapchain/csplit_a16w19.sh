#!/bin/sh
# Confirming A016045's candidate at 307223174679659356577: stop this
# machine's a16w19 run once every start below 258536912857616220160 is searched
# (20164987 segments); a DGX Spark searches the rest, 258536912857616220160..265600000000000000000
# (a16w67), so that 2e20..3e20 is finished at about 16:52.
cd "/Users/jbs/src/math/gapchain" || exit 1
need=20164987
while pgrep -f "gapsieve.*-c a16w19.ck" > /dev/null; do
    d=$(awk '$1 == "done" { print $2 }' a16w19.ck)
    if [ -n "$d" ] && [ "$d" -ge "$need" ]; then
        pkill -INT -f "^./gapsieve.*-c a16w19.ck"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  handoff: $d segments done, every start below 258536912857616220160 searched; stopped, a Spark has the rest (a16w67)" >> a16w19.handoff
        exit 0
    fi
    sleep 30
done
echo "$(date '+%Y-%m-%d %H:%M:%S')  handoff: the a16w19 run ended on its own" >> a16w19.handoff
