#!/bin/sh
# When the laptop's a16w19 run (2.329e20..2.656e20) finishes, start its piece of
# 3e20..5e20: 3.69e20..4.373e20 (a16w30).  A run that did not finish cleanly, or
# that reported a candidate chain, starts nothing.
cd "/Users/jbs/src/math/gapchain" || exit 1
until pgrep -f "gapsieve.*-c a16w19.ck" > /dev/null; do
    [ -e a16w19.log ] && tail -1 a16w19.log | grep -q "  exit " && break
    sleep 30
done
while pgrep -f "gapsieve.*-c a16w19.ck" > /dev/null; do
    sleep 30
done
if grep -qE 'candidate(: | [0-9])' a16w19.log; then
    echo "$(date '+%Y-%m-%d %H:%M:%S')  next: a16w19 reported a candidate chain; a16w30 not started" >> laptop_queue.log
    exit 1
fi
case "$(tail -1 a16w19.log)" in
*"exit 0") ;;
*)
    echo "$(date '+%Y-%m-%d %H:%M:%S')  next: a16w19 did not finish cleanly; a16w30 not started" >> laptop_queue.log
    exit 1
    ;;
esac
echo "$(date '+%Y-%m-%d %H:%M:%S')  next: a16w19 finished; starting a16w30, 3.69e20..4.373e20" >> laptop_queue.log
exec caffeinate -i nice -n 19 ./gapsieve -l 17 -s 3.69e20 -e 4.373e20 --first -c a16w30.ck -r --log a16w30.log -q --gpu > a16w30.out 2>&1 < /dev/null
