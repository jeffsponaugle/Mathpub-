#!/bin/sh
# Lower bounds for dimension 7 (see README, "Dimension 7"):
#   1. extend every 6-dim Q-class of CARAT's catalog by a 7th coordinate
#      (trivial, sign character, or extra -1) and keep one group per Q-class
#      fingerprint -> all 7-dim Q-classes with a 1-dim rational constituent;
#   2. split each into Z-classes (CARAT q2z) and count space groups (zverify);
#   3. certify the Z-classes pairwise non-conjugate (fingerprints + Z_equiv);
#   4. recompute the normalizers of the TOP largest contributors from scratch;
#   5. print the bounds.
# About 7 CPU-hours for step 2 (heavy tail: a few groups with huge form
# spaces take far longer; SPLIT_LIMIT caps them, which only lowers the bound).
# Optional: NREPS=<n> to run on the first n Q-classes only (testing);
# W=<dir> for the work files (default work/), R=<dir> for results.
set -e
cd "$(dirname "$0")/.."
ROOT=$PWD
J=${J:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
R=${R:-results/dim7}
TOP=${TOP:-20000}
W=${W:-work}
mkdir -p "$R" "$W/ext7"
make enantio zverify extend1 > /dev/null
[ -f results/qfiles_dim6.txt ] || find carat/tables/qcatalog/dim6 -type f \( -name 'group.*' -o -name 'min.*' -o -name 'max.*' \) | sort > results/qfiles_dim6.txt

tr '\n' '\0' < results/qfiles_dim6.txt | xargs -0 -P "$J" -n 10 ./extend1 -o "$W/ext7" > "$R/ext7_from6.tsv"
python3 scripts/dedup_q7.py "$R/ext7_from6.tsv" "$R/ext7_reps.tsv"
if [ -n "$NREPS" ]; then head -n "$NREPS" "$R/ext7_reps.tsv" > "$R/reps_run.tsv"; else cp "$R/ext7_reps.tsv" "$R/reps_run.tsv"; fi

python3 scripts/run7.py "$R/reps_run.tsv" "$W/z7" "$R" run -j "$J" --split-limit "${SPLIT_LIMIT:-3600}"

find "$W/z7" -type f -name 'x__*__Z*' ! -name '*.N2' -print0 | xargs -0 -P "$J" -n 200 ./extend1 -z > "$R/zfp7_run.tsv"
python3 scripts/zdedup7.py "$R/zfp7_run.tsv" "$R/zdups_run.tsv" -j "$J"

sort -t "$(printf '\t')" -k4,4nr "$R/z_run.tsv" | head -n "$TOP" \
    | awk -F'\t' '{z=$1; q=substr(z,4); sub(/__Z[0-9]+$/,"",q); print q "/" z}' > "$R/renorm_list.txt"
python3 scripts/renorm7.py "$R/renorm_list.txt" "$W/z7" "$R/renorm_run.tsv" -j "$J" --limit 900

python3 scripts/finalize7.py "$R/reps_run.tsv" "$R" | tee "$R/bound.txt"
