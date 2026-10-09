#!/usr/bin/env python3
"""Write OEIS b-files (lines "n a(n)") for the sequences matching a dbpal results file.

  bfile.py results/b2_9.tsv            -> results/bfiles/b259385.txt, b393014.txt, ...

Only terms below the certified search bound (from results/search_log.jsonl) are written, so the
b-file is complete as far as it goes.  The first index is the sequence's OEIS offset, and the
existing OEIS terms are checked to be an exact prefix of the new list.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import oeis  # noqa: E402
from report import bound_for  # noqa: E402


def offset(anum):
    txt = oeis._get(f"https://oeis.org/search?q=id:{anum}&fmt=text")
    m = re.search(r"%O " + anum + r" (-?\d+)", txt)
    return int(m.group(1)) if m else 1


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    path = sys.argv[1]
    b1, b2 = map(int, os.path.basename(path)[1:-4].split("_"))
    P, L, lim, _ = bound_for(b1, b2)
    terms = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2 and p[0].isdigit() and int(p[0]) < lim:
            terms[int(p[0])] = p[1] == "1"
    fam = oeis.family()
    outdir = os.path.join(os.path.dirname(os.path.abspath(path)), "bfiles")
    os.makedirs(outdir, exist_ok=True)
    for primes in (False, True):
        e = fam.get(((min(b1, b2), max(b1, b2)), primes))
        if not e:
            continue
        seq = sorted(n for n, isp in terms.items() if isp or not primes)
        known = [n for n in oeis.all_terms(e) if n < lim]
        if seq[: len(known)] != known:
            print(f"{e['anum']}: OEIS terms are not a prefix of the computed list -- not writing")
            continue
        off = offset(e["anum"])
        fn = os.path.join(outdir, "b" + e["anum"][1:] + ".txt")
        with open(fn, "w") as f:
            f.write(f"# {e['anum']}: {e['name']}\n")
            f.write(f"# complete for all terms < {P}^{L} (dbpal exhaustive search)\n")
            for i, n in enumerate(seq):
                f.write(f"{i + off} {n}\n")
        print(f"{fn}: {len(seq)} terms (OEIS has {len(known)}), complete below {P}^{L}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
