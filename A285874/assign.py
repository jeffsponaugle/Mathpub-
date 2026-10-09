#!/usr/bin/env python3
"""Splits the remaining tasks of one distributed run among machines.

Usage: assign.py N done.ckpt[,more.ckpt] block name=weight [name=weight ...]

Tasks already present in the given checkpoints are skipped. The rest are cut into blocks
of consecutive tasks (neighbours share many positions, so a block keeps one machine's hash
table useful) and dealt out greedily in proportion to the weights, which interleaves the
blocks so each machine gets a similar mix of cheap and expensive positions. Writes
tasks_<name>.txt files ("first last+1" per line) for perft --tasks.
"""
import sys

n, ckpts, block = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
weights = {k: float(v) for k, v in (a.split("=") for a in sys.argv[4:])}
done = set()
for f in filter(None, ckpts.split(",")):
    with open(f) as fh:
        done.update(int(line.split()[0]) for line in fh if not line.startswith("#"))
load = dict.fromkeys(weights, 0)
ranges = {m: [] for m in weights}
for start in range(0, n, block):
    todo = [i for i in range(start, min(n, start + block)) if i not in done]
    if not todo:
        continue
    m = min(weights, key=lambda m: (load[m] + len(todo)) / weights[m])
    load[m] += len(todo)
    run_start = prev = todo[0]
    for i in todo[1:] + [None]:
        if i != prev + 1 if i is not None else True:
            ranges[m].append((run_start, prev + 1))
            run_start = i
        prev = i if i is not None else prev
for m, rs in ranges.items():
    with open(f"tasks_{m}.txt", "w") as fh:
        fh.writelines(f"{a} {b}\n" for a, b in rs)
    print(f"{m}: {load[m]} tasks in {len(rs)} ranges")
print(f"already done: {len(done)}, total: {n}")
