#!/usr/bin/env python3
"""Aggregate geng -v outputs of split pieces and verify against OEIS.

usage: aggregate.py n MOD outdir [more outdirs...]

Sums the '>C <count> graphs with <k> edges' lines and '>Z' totals over all
complete pieces (res = 0..MOD-1), reports missing pieces, writes
<first outdir>/combined_n<n>.txt in geng -v format, and runs verify.py on
it (complement edge counts are translated there).
"""
import glob, json, os, re, subprocess, sys

def main():
    n, mod = int(sys.argv[1]), int(sys.argv[2])
    dirs = sys.argv[3:]
    counts, total, cpu, done = {}, 0, 0.0, set()
    stats = {"calls": 0, "tests": 0, "nonplanar": 0, "trivial": 0}
    for d in dirs:
        for f in glob.glob(os.path.join(d, f"part_*_of_{mod}.txt")):
            text = open(f).read()
            m = re.search(r">Z (\d+) graphs generated in ([\d.]+) sec", text)
            if not m:
                continue
            res = int(re.search(r"part_(\d+)_of_", os.path.basename(f)).group(1))
            if res in done:
                print(f"duplicate piece {res} in {f}, skipped"); continue
            done.add(res)
            total += int(m.group(1)); cpu += float(m.group(2))
            for c, k in re.findall(r">C (\d+) graphs with (\d+) edges", text):
                counts[int(k)] = counts.get(int(k), 0) + int(c)
            m2 = re.search(r"preprune calls=(\d+).*planarity_tests=(\d+) nonplanar=(\d+)", text)
            if m2:
                stats["calls"] += int(m2.group(1)); stats["tests"] += int(m2.group(2))
                stats["nonplanar"] += int(m2.group(3))
            m3 = re.search(r"skipped test\)=(\d+)", text)
            if m3: stats["trivial"] += int(m3.group(1))
    missing = sorted(set(range(mod)) - done)
    print(f"pieces complete: {len(done)}/{mod}" + (f"  MISSING: {missing[:20]}{'...' if len(missing) > 20 else ''}" if missing else ""))
    print(f"total graphs: {total}   total cpu: {cpu:.1f} s = {cpu/3600:.2f} core-hours")
    print(f"preprune calls {stats['calls']}, planarity tests {stats['tests']}, "
          f"nonplanar {stats['nonplanar']}, trivially planar {stats['trivial']}")
    out = os.path.join(dirs[0], f"combined_n{n}.txt")
    with open(out, "w") as fh:
        for k in sorted(counts):
            fh.write(f">C {counts[k]} graphs with {k} edges\n")
        fh.write(f">Z {total} graphs generated in {cpu:.2f} sec\n")
    print("wrote", out)
    # rows JSON (planar edge counts) for euler.py: known rows + this one
    e = n*(n-1)//2
    here = os.path.dirname(os.path.abspath(__file__))
    rows_path = os.path.join(here, "a049334_rows.json")
    rows = json.load(open(rows_path)) if os.path.exists(rows_path) else {}
    row = [0]*(3*n-5)
    for k, v in counts.items():
        kk = e - k
        if 0 <= kk < len(row): row[kk] = v
    rows[str(n)] = row
    rows_out = os.path.join(dirs[0], f"rows_with_n{n}.json")
    json.dump(rows, open(rows_out, "w"))
    print("wrote", rows_out, "(input for euler.py)")
    if missing:
        print("RESULT: INCOMPLETE (not verifying)"); sys.exit(2)
    here = os.path.dirname(os.path.abspath(__file__))
    sys.exit(subprocess.call([sys.executable, os.path.join(here, "verify.py"), str(n), "coplanar", out]))

if __name__ == "__main__":
    main()
