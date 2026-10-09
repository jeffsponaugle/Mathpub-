#!/bin/bash
# Pulls run logs from the other machines into runs/<machine>/ and appends a line to
# runs/results.txt for every job finished under at least two different seeds, listing
# which machine ran each seed and whether all results agree. Logs directly in runs/
# come from this machine (m4max).
cd "$(dirname "$0")" || exit 1
HOSTS=${HOSTS:-"m5pro=jbs@10.1.1.43 mathd=jbs@10.1.30.23 mathb=jbs@10.1.30.21 mathg=jbs@10.1.30.26 m2ultra=jeffsponauglejbs@10.1.10.6"}
for hm in $HOSTS; do
  name=${hm%%=*}
  mkdir -p "runs/$name"
  scp -q -o BatchMode=yes -o ConnectTimeout=8 "${hm#*=}:src/math/A285874/runs/*.log" "runs/$name/" 2>/dev/null
done
cat jobs_*.txt | while read -r variant depth _; do
  [[ -z "$variant" || "$variant" == \#* ]] && continue
  vals=() where=() seeds=()
  for f in runs/"${variant}_d${depth}"_seed*.log runs/*/"${variant}_d${depth}"_seed*.log; do
    [ -f "$f" ] || continue
    v=$(grep -h "^perft($depth) = " "$f" | tail -1 | awk '{print $3}')
    [ -n "$v" ] || continue
    s=${f##*_seed}
    m=$(dirname "$f"); m=${m#runs}; m=${m#/}
    vals+=("$v") where+=("seed ${s%.log}@${m:-m4max}") seeds+=("${s%.log}")
  done
  [ "$(printf '%s\n' "${seeds[@]}" | sort -u | grep -c .)" -ge 2 ] || continue
  status=AGREE
  for v in "${vals[@]}"; do [ "$v" = "${vals[0]}" ] || status=MISMATCH; done
  line="$variant a($depth) = ${vals[0]} $status [$(printf '%s, ' "${where[@]}" | sed 's/, $//')]"
  [ "$status" = MISMATCH ] && line="$line values: ${vals[*]}"
  grep -qF "$line" runs/results.txt 2>/dev/null || echo "$(date '+%F %T') $line" >> runs/results.txt
done
