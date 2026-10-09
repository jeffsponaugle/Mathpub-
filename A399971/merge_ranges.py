#!/usr/bin/env python3
"""Sum the RANGE lines printed by `a399971 -u first:last n` into a(n).

Usage: merge_ranges.py FILE... (or RANGE lines on stdin)

Checks that the ranges of each run (same n, split, order, unit count) cover
every work unit exactly once before printing the totals.
"""
import re
import sys
from collections import defaultdict

FIELDS = ("compact", "labeled", "puzzles", "multi", "covers", "errors", "seconds")


def main():
    texts = [open(f).read() for f in sys.argv[1:]] or [sys.stdin.read()]
    runs = defaultdict(list)
    for text in texts:
        for line in text.splitlines():
            if not line.startswith("RANGE "):
                continue
            kv = dict(re.findall(r"(\w+)=([\d.]+)", line))
            key = (int(kv["n"]), int(kv["split"]), int(kv["order"]), int(kv["units"]))
            runs[key].append(kv)
    if not runs:
        sys.exit("no RANGE lines found")
    status = 0
    for (n, split, order, units), parts in sorted(runs.items()):
        parts.sort(key=lambda kv: int(kv["first"]))
        covered, problems = 0, []
        for kv in parts:
            first, last = int(kv["first"]), int(kv["last"])
            if first > covered:
                problems.append(f"gap {covered}..{first}")
            elif first < covered:
                problems.append(f"overlap at {first}..{min(last, covered)}")
            covered = max(covered, last)
        if covered < units:
            problems.append(f"gap {covered}..{units}")
        tot = {f: sum(float(kv[f]) if f == "seconds" else int(kv[f]) for kv in parts) for f in FIELDS}
        print(f"n={n} (split {split}, order {order}): {len(parts)} ranges, {units} work units")
        if problems:
            status = 1
            print("  INCOMPLETE: " + ", ".join(problems))
            print(f"  partial sum so far: a(n) >= {tot['compact']}")
            continue
        print(f"a({n}) = {tot['compact']}")
        print(f"  compact solutions, not reduced for symmetry: {tot['labeled']}")
        print(f"  compact starting positions, reduced: {tot['puzzles']} "
              f"({tot['multi']} solutions share a puzzle with another)")
        print(f"  canonical covers checked: {tot['covers']}, compute time {tot['seconds'] / 3600:.1f} h")
        if tot["errors"]:
            status = 1
            print(f"  INTERNAL ERROR count: {tot['errors']}")
    sys.exit(status)


if __name__ == "__main__":
    main()
