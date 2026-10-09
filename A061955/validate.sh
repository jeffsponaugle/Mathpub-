#!/bin/sh
# validate.sh [N] [JOBS] -- check this machine's build: for all 6 families x bases 2..25,
# every residue V(n) mod n for n < N (default 1500) must agree between brute2 (direct
# definition) and concat (scalar `res`, 4-lane `resl`, and floating-point vector `resv`
# paths; resv only for the reversed-digit families).  Prints FAIL lines and
# a summary; exit status 0 only if everything matches.
cd "$(dirname "$0")" || exit 1
N=${1:-1500}
J=${2:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
T=$(mktemp -d)
for fam in Ln Rn Lk Rk Ld Rd; do for b in $(seq 2 25); do echo "$fam $b"; done; done |
    xargs -P "$J" -L 1 sh -c '
        ./brute2 $0 $1 1 '"$N"' > '"$T"'/$0_$1.b
        ./concat res $0:$1 1 '"$N"' > '"$T"'/$0_$1.f
        ./concat resl $0:$1 1 '"$N"' > '"$T"'/$0_$1.l
        case $0 in ?n) cp '"$T"'/$0_$1.b '"$T"'/$0_$1.v ;; *) ./concat resv $0:$1 1 '"$N"' > '"$T"'/$0_$1.v ;; esac
        if cmp -s '"$T"'/$0_$1.b '"$T"'/$0_$1.f && cmp -s '"$T"'/$0_$1.b '"$T"'/$0_$1.l && cmp -s '"$T"'/$0_$1.b '"$T"'/$0_$1.v
        then echo ok > '"$T"'/$0_$1.ok; else echo "FAIL $0 base $1"; fi'
ok=$(ls "$T" | grep -c '\.ok$')
rm -rf "$T"
echo "validate: $ok/144 family-base combinations match for n < $N on $(hostname)"
[ "$ok" -eq 144 ]
