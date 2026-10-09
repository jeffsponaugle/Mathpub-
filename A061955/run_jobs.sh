#!/bin/sh
# run_jobs.sh JOBFILE [THREADS] -- run `concat search` jobs one after another.
#
# JOBFILE lines:  A-number LO HI [# comment]     (blank lines and # comments ignored)
# Results go to runs/<A-number>/ (hits_<A>.txt, its .progress file, search.out/.log);
# starts/ends are logged to runs/jobs.log.  A finished job leaves runs/<A>/done_<LO>_<HI>
# and is skipped from then on; an interrupted job resumes from the contiguous point in
# its .progress file.  Safe to stop (kill the runner, then concat) and restart at any time.
cd "$(dirname "$0")" || exit 1
JOBS=$1
TH=${2:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
[ -f "$JOBS" ] || { echo "usage: $0 JOBFILE [THREADS]" >&2; exit 1; }
mkdir -p runs
grep -v '^[[:space:]]*\(#\|$\)' "$JOBS" | while read -r seq lo hi _rest; do
    dir=runs/$seq
    prog=$dir/hits_$seq.txt.progress
    mark=$dir/done_${lo}_${hi}
    mkdir -p "$dir"
    [ -f "$mark" ] && continue
    start=$lo
    if [ -f "$prog" ]; then             # "<tag> base <b> complete [a, x)"
        a=$(sed -n 's/.*complete \[\([0-9]*\), \([0-9]*\)).*/\1/p' "$prog")
        x=$(sed -n 's/.*complete \[\([0-9]*\), \([0-9]*\)).*/\2/p' "$prog")
        if [ -n "$x" ] && [ "$a" -ge "$lo" ] && [ "$a" -lt "$hi" ] && [ "$x" -gt "$start" ] && [ "$x" -le "$hi" ]; then
            start=$x
        fi
    fi
    if [ "$start" -lt "$hi" ]; then
        echo "$(date '+%F %T') $(hostname -s) start $seq [$start, $hi)" >> runs/jobs.log
        nice -n 10 ./concat search "$seq" "$start" "$hi" -t "$TH" -o "$dir" >> "$dir/search.out" 2>> "$dir/search.log" < /dev/null
        rc=$?
        echo "$(date '+%F %T') $(hostname -s) end   $seq [$start, $hi) exit=$rc" >> runs/jobs.log
        x=$(sed -n 's/.*complete \[\([0-9]*\), \([0-9]*\)).*/\2/p' "$prog" 2>/dev/null)
        [ "$rc" -eq 0 ] && [ "$x" = "$hi" ] && touch "$mark"
    else
        touch "$mark"
    fi
done
