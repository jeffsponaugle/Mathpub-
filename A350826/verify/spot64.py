#!/usr/bin/env python3
"""
spot64.py - recount the [1e19, 2^64) part of randomly chosen a(20) chunks with the CPU tool
(../a350826, independent of the CUDA code) and compare with the GPU logs' first bin.

Usage: spot64.py K ORIG_SEGMENTS...      (segments as for combine_a20.py, paths relative to ..)
Appends results to spot64.txt.
"""
import random
import subprocess
import sys
import time

K = int(sys.argv[1])
orig = {}
for seg in sys.argv[2:]:
    path, lo, hi = seg.rsplit(':', 2)
    for raw in open(path):
        f = raw.split()
        if len(f) >= 6 and f[0] == 'C' and int(lo) <= int(f[1]) < int(hi):
            orig[int(f[1])] = f[5]
C, CH = 1516640125, 16384
out = open('spot64.txt', 'a')
rng = random.Random(time.time_ns())
picks = rng.sample(sorted(orig), K)
bad = 0
t0 = time.time()
for n, k in enumerate(picks, 1):
    a, b = k * CH, min((k + 1) * CH, C)
    log = f'/tmp/spot64_{k}.chunks'
    subprocess.run(['rm', '-f', log])
    subprocess.run(['nice', '-n', '10', '../a350826', 'count', '1e19', '2^64', '-w', '37', '-c', str(CH),
                    '-C', f'{a}:{b}', '-q', '-L', log], capture_output=True, check=True)
    line = open(log).read().split()
    mine = line[5]
    ok = mine == orig[k]
    bad += not ok
    msg = f"{'OK      ' if ok else 'MISMATCH'} chunk {k}: GPU {orig[k]}  CPU {mine}"
    print(msg, flush=True)
    out.write(msg + '\n')
    out.flush()
summary = f"spot64: {K} random chunks of [1e19, 2^64) recounted on the CPU, {bad} mismatches ({time.time() - t0:.0f} s)"
print(summary)
out.write(summary + '\n')
