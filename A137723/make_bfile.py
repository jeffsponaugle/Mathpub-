#!/usr/bin/env python3
"""
make_bfile.py [scan_table] [structural_terms]

Merge the table written by `a137723 scan`/`next` (default next_1e11.txt) with
structural_terms.txt (terms found by `a137723 search`) and write

  b137723.txt  the OEIS b-file (n a(n) per line, header comments start with #)
  DATA.txt     the consecutive terms as an OEIS DATA line

Only the consecutive run a(1..k) is written; the first missing n is reported.
Every structural term must agree with the scan table where both know it, and
the forms given in structural_terms.txt (like 2^54-32) are re-evaluated.
"""
import re
import sys
import datetime

scan_file = sys.argv[1] if len(sys.argv) > 1 else "next_1e11.txt"
struct_file = sys.argv[2] if len(sys.argv) > 2 else "structural_terms.txt"

rows = {}
scan_limit = None
for line in open(scan_file):
    m = re.match(r"\s*(\d+)\s+(\d+)\s", line)
    if m:
        rows[int(m.group(1))] = int(m.group(2))
    m = re.search(r"scan to (\d+)", line)
    if m:
        scan_limit = int(m.group(1))

for line in open(struct_file):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    n, a, form = line.split()[:3]
    n, a = int(n), int(a)
    if form != "-":
        val = eval(form.replace("^", "**"), {"__builtins__": {}})
        if val != a:
            sys.exit(f"form {form} does not evaluate to a({n}) = {a}")
    if n in rows and rows[n] != a:
        sys.exit(f"conflict for n = {n}: scan {rows[n]} vs structural {a}")
    rows[n] = a

k = 0
while k + 1 in rows:
    k += 1
missing = [n for n in sorted(set(range(1, max(rows) + 1)) - set(rows))]

today = datetime.date.today().strftime("%b %d %Y")
with open("b137723.txt", "w") as f:
    f.write("# A137723: First occurrence of a set of n consecutive numbers having at least one prime gap\n")
    f.write("# in their factorization: a(n) = smallest number of this set.\n")
    f.write(f"# Table of n, a(n) for n = 1..{k}.\n")
    f.write(f"# a(1)-a(31) as in the OEIS entry; a(32)-a({k}) computed with a137723.c ({today}):\n")
    f.write(f"# exhaustive scan of 1..{scan_limit} plus the structural search for even n (see README.md).\n")
    for n in range(1, k + 1):
        f.write(f"{n} {rows[n]}\n")
with open("DATA.txt", "w") as f:
    f.write(", ".join(str(rows[n]) for n in range(1, k + 1)) + "\n")

print(f"b137723.txt: n = 1..{k}; DATA.txt: {k} terms")
print(f"first missing n: {missing[:8]}")
print(f"terms above 10^11: {sum(1 for n in range(1, k + 1) if rows[n] > 10**11)}, "
      f"largest: a({max(range(1, k + 1), key=lambda n: rows[n])}) with "
      f"{len(str(max(rows[n] for n in range(1, k + 1))))} digits")
