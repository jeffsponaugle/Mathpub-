#!/bin/sh
# run_pool.sh POOLFILE [THREADS] -- run jobs from a pool shared by several machines.
#
# POOLFILE lines: A-number LO HI [# comment]   (the same file on every machine).
# A job is claimed by atomically creating claims/<A>_<LO>_<HI> on the coordinator
# (mkdir succeeds for exactly one caller); the claiming machine runs it through
# run_jobs.sh, so it gets the usual resume and runs/<A>/done_<LO>_<HI> marker.  Jobs this
# machine claimed earlier but did not finish (pool/mine.txt) are resumed first.  Jobs
# are claimed in file order, so earlier lines start first.
#   POOL_COORD   coordinator ssh target   (default jbs@10.1.30.21 = mathb)
#   POOL_CLAIMS  claims directory under the coordinator's home
cd "$(dirname "$0")" || exit 1
POOL=$1
TH=${2:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
COORD=${POOL_COORD:-jbs@10.1.30.21}
CDIR=${POOL_CLAIMS:-src/math/A061955/pool/claims}
HOST=$(hostname -s)
[ -f "$POOL" ] || { echo "usage: $0 POOLFILE [THREADS]" >&2; exit 1; }
mkdir -p pool runs
touch pool/mine.txt

is_coord() { hostname -I 2>/dev/null | tr ' ' '\n' | grep -qx "${COORD#*@}"; }

claim() {                       # prints CLAIMED if this machine now owns job $1
    for try in 1 2 3 4 5 6; do
        if is_coord; then
            mkdir -p "$HOME/$CDIR" || return 1
            if mkdir "$HOME/$CDIR/$1" 2>/dev/null; then echo "$HOST" > "$HOME/$CDIR/$1/host"; echo CLAIMED; fi
            return 0
        fi
        out=$(ssh -n -o BatchMode=yes -o ConnectTimeout=10 "$COORD" \
            "mkdir -p ~/$CDIR && if mkdir ~/$CDIR/$1 2>/dev/null; then echo $HOST > ~/$CDIR/$1/host; echo CLAIMED; else echo TAKEN; fi") \
            && { [ "$out" = CLAIMED ] && echo CLAIMED; return 0; }
        sleep 30
    done
    return 1
}

run_one() {
    printf '%s %s %s\n' "$1" "$2" "$3" > "pool/job_$HOST.txt"
    sh run_jobs.sh "pool/job_$HOST.txt" "$TH" < /dev/null
}

while read -r a lo hi; do       # resume jobs claimed earlier
    [ -f "runs/$a/done_${lo}_${hi}" ] || run_one "$a" "$lo" "$hi"
done < pool/mine.txt

grep -v '^[[:space:]]*\(#\|$\)' "$POOL" | while read -r a lo hi _rest; do
    [ -f "runs/$a/done_${lo}_${hi}" ] && continue
    r=$(claim "${a}_${lo}_${hi}") || { echo "$(date '+%F %T') $HOST pool: cannot reach $COORD, stopping" >> runs/jobs.log; exit 1; }
    [ "$r" = CLAIMED ] || continue
    echo "$a $lo $hi" >> pool/mine.txt
    echo "$(date '+%F %T') $HOST pool: claimed $a [$lo, $hi)" >> runs/jobs.log
    run_one "$a" "$lo" "$hi"
done
