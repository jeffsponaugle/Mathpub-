#!/bin/sh
# laptop: after a16w19, its piece of 2e20..3e20, search its piece of 3e20..5e20,
# 4.443e20..4.733e20 (a16w43), then its piece of 5e20..7e20, 6.252e20..6.643e20 (a16w56), one
# after the other.  A run that does not finish cleanly or reports a candidate
# chain starts nothing further.
cd "/Users/jbs/src/math/gapchain" || exit 1
ok() {    # name: finished cleanly, with no candidate chain
    grep -qE 'candidate(: | [0-9])' "$1.log" && return 1
    case "$(tail -1 "$1.log")" in *"exit 0") return 0 ;; esac
    return 1
}
run() {    # name, start, end
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: starting $1, $2..$3" >> laptop_queue.log
    caffeinate -i nice -n 19 ./gapsieve -l 17 -s "$2" -e "$3" --first -c "$1.ck" -r --log "$1.log" -q --gpu \
        > "$1.out" 2>&1 < /dev/null
    ok "$1" && return 0
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: $1 did not finish cleanly or reported a candidate; stopping" >> laptop_queue.log
    return 1
}
until pgrep -f "gapsieve.*-c a16w19.ck" > /dev/null; do
    [ -e a16w19.log ] && tail -1 a16w19.log | grep -q "  exit " && break
    sleep 30
done
while pgrep -f "gapsieve.*-c a16w19.ck" > /dev/null; do
    sleep 30
done
if ! ok a16w19; then
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: a16w19 did not finish cleanly or reported a candidate; nothing started" >> laptop_queue.log
    exit 1
fi
run a16w43 4.443e20 4.733e20 &&
run a16w56 6.252e20 6.643e20 &&
echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: a16w43 and a16w56 finished" >> laptop_queue.log
