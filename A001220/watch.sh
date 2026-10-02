#!/bin/bash
# watch.sh [SECS [HOST...]] -- poll the Sparks (default: atom1 and atom2) every SECS seconds
# (default 600) and exit with a report when a NEW solution (FOUND), double-check MISMATCH or
# FLT-ERROR line appears (lines present at the first poll are the baseline and are ignored),
# or when a host has had no wieferich_cuda scan running on two consecutive polls (finished,
# crashed or stopped; hand-overs between queued scans take seconds).
# Works with the bash 3.2 that ships with macOS.
SECS=${1:-600}
shift
HOSTS=("$@")
[ ${#HOSTS[@]} -gt 0 ] || HOSTS=(10.1.30.36 10.1.30.37)
IDLE=(); BASE=(); HAVE=()
for i in "${!HOSTS[@]}"; do IDLE[$i]=0; BASE[$i]=""; HAVE[$i]=0; done
while true; do
    for i in "${!HOSTS[@]}"; do
        H=${HOSTS[$i]}
        out=$(ssh -o ConnectTimeout=20 -o BatchMode=yes jbs@$H 'cd /home/jbs/A001220/runs &&
            { grep -H "^FOUND\|^FLT-ERROR" *.txt 2>/dev/null; grep -H "^MISMATCH" *.dcheck 2>/dev/null;
              pgrep -f "^/home/jbs/A001220/cuda/[w]ieferich_cuda scan" > /dev/null || echo NOSCAN; }') ||
            { echo "$(date +%FT%T) ssh to $H failed, retrying"; continue; }
        events=$(echo "$out" | grep -v '^NOSCAN$')
        if [ "${HAVE[$i]}" -eq 0 ]; then
            BASE[$i]=$events; HAVE[$i]=1
            [ -n "$events" ] && { echo "$(date +%FT%T) $H baseline (ignored):"; echo "$events"; }
        else
            new=$(echo "$events" | grep -vxF -f <(echo "${BASE[$i]}"))
            if [ -n "$new" ]; then echo "$(date +%FT%T) $H:"; echo "$new"; exit 0; fi
        fi
        if echo "$out" | grep -q '^NOSCAN$'; then
            IDLE[$i]=$(( ${IDLE[$i]} + 1 ))
            if [ "${IDLE[$i]}" -ge 2 ]; then echo "$(date +%FT%T) $H: no scan running"; exit 0; fi
        else
            IDLE[$i]=0
        fi
    done
    sleep "$SECS"
done
