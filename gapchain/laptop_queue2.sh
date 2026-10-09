#!/bin/sh
# After the laptop's d16w14 run finishes, search the top of the M1 Max's
# d16w13 (d16w20) and of the Intel Mac's d16w19 (d16w21), whose watchers stop
# those runs at the split points, then the laptop's share of 2e20..3e20
# (d16w22), one after another, on the GPU.  A run that does not finish cleanly
# or reports a candidate chain starts nothing further.
cd /Users/jbs/src/math/gapchain || exit 1
while pgrep -f "gapsieve.*-c d16w14.ck" > /dev/null; do
    sleep 30
done
ok() {    # name: finished cleanly, with no candidate chain
    grep -qE 'candidate(: | [0-9])' "$1.log" && return 1
    case "$(tail -1 "$1.log")" in *"exit 0") return 0 ;; esac
    return 1
}
run() {    # name, start, end
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: starting $1, $2..$3" >> laptop_queue.log
    caffeinate -i nice -n 19 ./gapsieve -l 17 -s "$2" -e "$3" --first --gaps dec \
        -c "$1.ck" -r --log "$1.log" -q --gpu > "$1.out" 2>&1 < /dev/null
    ok "$1" && return 0
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: $1 did not finish cleanly or reported a candidate; stopping" >> laptop_queue.log
    return 1
}
if ! ok d16w14; then
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: d16w14 did not finish cleanly or reported a candidate; nothing started" >> laptop_queue.log
    exit 1
fi
run d16w20 98191564376758026240 107000000000000000000 &&
run d16w21 199943325641973104640 200000000000000000000 &&
run d16w22 2e20 2.645e20 &&
echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: d16w20, d16w21 and d16w22 finished" >> laptop_queue.log
