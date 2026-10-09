#!/bin/sh
# When the laptop's queue (laptop_queue.sh) has finished d16w8 to d16w12,
# start its next A263049 piece on the GPU, 1.07e20..1.85e20 (d16w14).  A queue
# that stopped early starts nothing.  Usage: next_laptop.sh QUEUE_PID
cd /Users/jbs/src/math/gapchain || exit 1
pid=$1
while kill -0 "$pid" 2>/dev/null; do
    sleep 30
done
case "$(tail -1 laptop_queue.log)" in
*"d16w8 to d16w12 finished") ;;
*)
    echo "$(date '+%Y-%m-%d %H:%M:%S')  next: the queue did not finish; d16w14 not started" >> laptop_queue.log
    exit 1
    ;;
esac
echo "$(date '+%Y-%m-%d %H:%M:%S')  next: queue finished; starting d16w14, 1.07e20..1.85e20" >> laptop_queue.log
exec caffeinate -i nice -n 19 ./gapsieve -l 17 -s 1.07e20 -e 1.85e20 --first --gaps dec \
    -c d16w14.ck -r --log d16w14.log -q --gpu > d16w14.out 2>&1 < /dev/null
