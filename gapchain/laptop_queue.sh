#!/bin/sh
# After the laptop's d16w4 run finishes (at 5.75e19), search the tops of the
# other A263049 pieces, whose watchers (split_d16w*.sh) stop those runs at the
# split points, so that every piece finishes together: d16w8 (mathc's top),
# d16w9 (the M1 Pro's), d16w10 (math4's), d16w11 (the WSL laptop's) and
# d16w12 (the Intel Mac's), one after another, on the GPU.  A run that does
# not finish cleanly starts nothing further.  Usage: laptop_queue.sh PID
cd /Users/jbs/src/math/gapchain || exit 1
pid=$1
while kill -0 "$pid" 2>/dev/null; do
    sleep 30
done
last=$(tail -1 d16w4.log)
case "$last" in
*"exit 0") ;;
*)
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: d16w4 did not finish cleanly ($last); nothing started" >> laptop_queue.log
    exit 1
    ;;
esac
run() {    # name, start, end
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: starting $1, $2..$3" >> laptop_queue.log
    caffeinate -i nice -n 19 ./gapsieve -l 17 -s "$2" -e "$3" --first --gaps dec \
        -c "$1.ck" -r --log "$1.log" -q --gpu > "$1.out" 2>&1 < /dev/null
    last=$(tail -1 "$1.log")
    case "$last" in
    *"exit 0") return 0 ;;
    esac
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: $1 did not finish cleanly ($last); stopping" >> laptop_queue.log
    return 1
}
run d16w8 28637618584061050880 31000000000000000000 &&
run d16w9 35957542650153533440 38000000000000000000 &&
run d16w10 59937403886031994880 60500000000000000000 &&
run d16w11 63084449125441863680 64000000000000000000 &&
run d16w12 65531703868629975040 66000000000000000000 &&
echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: d16w8 to d16w12 finished" >> laptop_queue.log
