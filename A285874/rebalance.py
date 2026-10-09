#!/usr/bin/env python3
"""Moves part of one machine's unfinished task list to others' phase-2 lists.

Usage: rebalance.py tasks_SRC.txt SRC.ckpt KEEP block name=weight ...

Reads SRC's task list (ranges, in processing order) and its checkpoint, keeps the first
KEEP unfinished tasks for SRC (rewriting tasks_SRC.txt), and deals the rest out in blocks,
weighted, to tasks_<name>2.txt files.
"""
import sys

src, ckpt, keep, block = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
weights = {k: float(v) for k, v in (a.split("=") for a in sys.argv[5:])}
done = {int(l.split()[0]) for l in open(ckpt) if not l.startswith("#")}
order = [i for line in open(src) for i in range(*map(int, line.split()))]
rest = [i for i in order if i not in done]


def ranges(ids):
    out, start = [], None
    for k, i in enumerate(ids):
        if start is None:
            start = i
        if k + 1 == len(ids) or ids[k + 1] != i + 1:
            out.append((start, i + 1))
            start = None
    return out


def write(name, ids):
    with open(name, "w") as fh:
        fh.writelines(f"{a} {b}\n" for a, b in ranges(ids))
    print(f"{name}: {len(ids)} tasks")


write(src, rest[:keep])
load = dict.fromkeys(weights, 0)
share = {m: [] for m in weights}
for s in range(keep, len(rest), block):
    chunk = rest[s:s + block]
    m = min(weights, key=lambda m: (load[m] + len(chunk)) / weights[m])
    load[m] += len(chunk)
    share[m] += chunk
for m in weights:
    write(f"tasks_{m}2.txt", sorted(share[m]))
print(f"done before: {len(done & set(order))}, remaining: {len(rest)}")
