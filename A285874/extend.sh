#!/bin/bash
# Runs the perft jobs in a job file. For each job and each seed in $SEEDS, computes perft
# with that seed and that seed's split (split1 for seed 1, split2 for seed 2), writing
# runs/<variant>_d<depth>_seed<seed>.{log,ckpt}. The two seeds use independent Zobrist keys
# and different task decompositions; collect.sh compares them.
#
# Job file lines: variant depth split1 split2 [split3]   ('#' starts a comment)
# Seed 3 (a third confirmation) uses split3, or split2 if absent.
# Env: SEEDS (default "1 2"), HASH table size in MB (default 8192).
# Resumable: finished runs are skipped and interrupted runs continue from their checkpoints.
cd "$(dirname "$0")" || exit 1
HASH=${HASH:-8192}
SEEDS=${SEEDS:-1 2}
while pgrep -x perft > /dev/null; do sleep 30; done  # never run two computations at once
while read -r variant depth s1 s2 s3 _; do
  [[ -z "$variant" || "$variant" == \#* ]] && continue
  for seed in $SEEDS; do
    case "$seed" in 1) split=$s1 ;; 2) split=$s2 ;; *) split=${s3:-$s2} ;; esac
    log="runs/${variant}_d${depth}_seed${seed}.log"
    grep -q "^perft($depth) = " "$log" 2>/dev/null && continue
    nice -n 10 ./perft --variant "$variant" --depth "$depth" --hash "$HASH" --split "$split" --seed "$seed" \
      --report 300 --ckpt "runs/${variant}_d${depth}_seed${seed}.ckpt" < /dev/null >> "$log" 2>&1 || exit 1
  done
done < "$1"
