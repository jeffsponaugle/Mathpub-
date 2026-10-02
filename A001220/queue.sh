#!/bin/bash
# queue.sh QUEUEFILE
#
# Run the scans listed in QUEUEFILE one after another on this Spark.  Lines are
#     START END TAG BASES          (decimal START/END; '#' starts a comment)
# and run in the directory of QUEUEFILE with state/log files TAG.*.  An item is done when
# TAG.txt has the SUMMARY line of its full range [START, END); the first item that is not
# done is run next, resuming from TAG.state if it exists.  The file is re-read after every
# scan, so items can be added, removed or moved to another Spark while the queue runs.
# If a scan ends without completing its range (stopped by hand, crash), the queue stops.
# Every scan gets a light CPU double-checker (dcheck.sh: 4 threads pinned to efficiency
# cores 10-13, every 96th chunk).
set -u
Q=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
cd "$(dirname "$Q")"
done_item() { grep -q "^SUMMARY \[$1, $2)" "$3.txt" 2> /dev/null; }
while true; do
    item=""
    while read -r s e tag bases rest; do
        case "${s:-#}" in \#*) continue ;; esac
        done_item "$s" "$e" "$tag" && continue
        item="$s $e $tag $bases"
        break
    done < "$Q"
    if [ -z "$item" ]; then echo "$(date +%FT%T) queue finished"; exit 0; fi
    set -- $item
    s=$1; e=$2; tag=$3; bases=$4
    echo "$(date +%FT%T) running $tag: scan $s $e -b $bases"
    /home/jbs/A001220/cuda/wieferich_cuda scan "$s" "$e" -b "$bases" -S "$tag.state" -L "$tag.chunks" -i 60 \
        >> "$tag.txt" 2>> "$tag.log" < /dev/null &
    gp=$!
    taskset -c 10-13 /home/jbs/A001220/dcheck.sh "$tag.chunks" "$bases" 96 4 "$tag.dcheck" "$gp" \
        >> "$tag.dcheck.log" 2>&1 < /dev/null &
    wait "$gp"
    if ! done_item "$s" "$e" "$tag"; then
        echo "$(date +%FT%T) $tag ended before $e: queue stops"
        exit 0
    fi
    echo "$(date +%FT%T) $tag complete"
done
