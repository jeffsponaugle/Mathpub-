#!/bin/sh
# After d16w20 (the last piece below A263049's find), the laptop moves to
# A016045: the top of the Studio's a16w9, 166246697682959073280..1.73e20 (a16w17), then
# its piece of 2e20..3e20, 2.329e20..2.656e20 (a16w19), on the GPU.  A run that
# does not finish cleanly or reports a candidate chain starts nothing further.
cd /Users/jbs/src/math/gapchain || exit 1
while pgrep -f "gapsieve.*-c d16w20.ck" > /dev/null; do
    sleep 30
done
echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: d16w20 ended: $(tail -1 d16w20.log)" >> laptop_queue.log
run() {    # name, start, end
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: starting $1, $2..$3" >> laptop_queue.log
    caffeinate -i nice -n 19 ./gapsieve -l 17 -s "$2" -e "$3" --first \
        -c "$1.ck" -r --log "$1.log" -q --gpu > "$1.out" 2>&1 < /dev/null
    if grep -qE 'candidate(: | [0-9])' "$1.log"; then
        echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: $1 reported a candidate chain; stopping" >> laptop_queue.log
        return 1
    fi
    case "$(tail -1 "$1.log")" in *"exit 0") return 0 ;; esac
    echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: $1 did not finish cleanly; stopping" >> laptop_queue.log
    return 1
}
run a16w17 166246697682959073280 173000000000000000000 &&
run a16w19 2.329e20 2.656e20 &&
echo "$(date '+%Y-%m-%d %H:%M:%S')  queue: a16w17 and a16w19 finished" >> laptop_queue.log
