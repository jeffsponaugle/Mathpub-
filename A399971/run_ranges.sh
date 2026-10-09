#!/bin/bash
# Resumable chunked run of a399971 for one n.
#
#   run_ranges.sh N CHUNKS [FIRST_CHUNK [LAST_CHUNK]] [-- extra a399971 options]
#
# Splits the work units into CHUNKS equal ranges and computes chunks
# FIRST_CHUNK..LAST_CHUNK (default: all), writing results/nN/chunk_K.txt.
# Chunks that already have a RANGE line are skipped, so the script can be
# stopped and restarted, and different machines can take different chunk
# ranges (copy the result files together afterwards).  Then:
#
#   ./merge_ranges.py results/nN/chunk_*.txt
set -euo pipefail
cd "$(dirname "$0")"
pos=()
while [ $# -gt 0 ] && [ "$1" != "--" ]; do pos+=("$1"); shift; done
[ "${1:-}" = "--" ] && shift
[ ${#pos[@]} -ge 2 ] || { sed -n '2,13p' "$0" >&2; exit 1; }
n=${pos[0]} chunks=${pos[1]}
first=${pos[2]:-0}
last=${pos[3]:-$((chunks - 1))}
units=$(./a399971 -q -u 0:0 "$@" "$n" | sed -n 's/.* units=\([0-9]*\) .*/\1/p')
[ -n "$units" ] || { echo "could not determine the number of work units" >&2; exit 1; }
out=results/n$n
mkdir -p "$out"
echo "n=$n: $units work units in $chunks chunks; doing chunks $first..$last" >&2
for ((k = first; k <= last; k++)); do
    f=$out/chunk_$k.txt
    if [ -s "$f" ] && grep -q '^RANGE ' "$f"; then continue; fi
    a=$((units * k / chunks)) b=$((units * (k + 1) / chunks))
    [ "$a" -lt "$b" ] || continue
    echo "$(date '+%F %T') chunk $k: units $a..$b" >&2
    ./a399971 -q -u "$a:$b" "$@" "$n" > "$f.tmp" && mv "$f.tmp" "$f"
done
