#!/bin/bash
# dcheck.sh CHUNKLOG BASES [STEP] [THREADS] [OUT] [PID]
#
# Double-check a GPU run: recompute every STEP-th chunk of CHUNKLOG (a -L log of
# wieferich_cuda) with the CPU tool ../wieferich (independent code: 64-bit limbs,
# double-precision quotient estimate, primesieve on the CPU) and compare prime count
# and per-base checksums, which must be identical.  The CPU tool splits the chunk
# into 50 sub-chunks to use THREADS cores; checksums are sums mod 2^64, so its
# SUMMARY totals are compared with the GPU line.  Appends "MATCH lo hi" or
# "MISMATCH ..." to OUT (default CHUNKLOG.dcheck).  Resumable: chunks already in OUT
# are skipped.  Waits for new chunks while process PID (the GPU run) is alive.
set -u
LOG=$1; BASES=$2; STEP=${3:-16}; T=${4:-12}; OUT=${5:-$1.dcheck}; PID=${6:-0}
CPU="$(dirname "$0")/wieferich"
touch "$OUT"
k=0
while true; do
    progressed=0
    while IFS= read -r line; do
        k=$((k + 1))
        [ $((k % STEP)) -eq 0 ] || continue
        set -- $line
        lo=$2; hi=$3
        grep -q " $lo $hi " "$OUT" && continue
        len=$((hi - lo)); sub=$len
        [ $((len % 50)) -eq 0 ] && sub=$((len / 50))
        cpu=$("$CPU" scan "$lo" "$hi" -c "$sub" -b "$BASES" -q -t "$T" | awk -v lo="$lo" -v hi="$hi" '
            /^SUMMARY/ { primes = $5; flt = $9; rng = $2 " " $3 }
            /^base/    { for (i = 1; i < NF; i++) if ($i == "cks") { s = s " " $2 ":" $(i + 1) } }
            END { if (rng != "[" lo ", " hi ")" || flt != 0) print "BAD range " rng " flt " flt;
                  else print "C " lo " " hi " " primes s }')
        if [ "$cpu" = "$line" ]; then
            echo "MATCH $lo $hi $(date +%FT%T)" >> "$OUT"
        else
            echo "MISMATCH $lo $hi gpu: $line cpu: $cpu" >> "$OUT"
        fi
        progressed=1
    done < <(sed -n "$((k + 1)),\$p" "$LOG")
    if [ "$PID" -gt 0 ] && kill -0 "$PID" 2> /dev/null; then
        [ $progressed -eq 1 ] || sleep 60
    else
        [ $progressed -eq 1 ] || break
    fi
done
