#!/usr/bin/env python3
"""
combine_a20.py - merge the a(20) chunk logs of several runs (two Sparks).

All runs are `a350826_cuda count 1e19 1e20 -b BINS` with the same wheel (37), bound and chunk
size (16384), over different class ranges (-C); each writes one line per chunk with absolute
chunk indices.  Every argument LOG:LO:HI takes the lines of chunks LO <= k < HI from LOG.  The
merge checks that every chunk 0..92568 is covered exactly once, that chunks present in several
logs have identical lines (cross-check between GPUs), and adds up survivors, candidates,
pseudoprimes and the per-bin counts and checksums (mod 2^64).

For a(20): atom1 ran all classes and was stopped after chunk 53417 (runs/a20.chunks), atom2 ran
chunks 53416.. (runs/a20b.chunks) and was stopped after Y+1 when atom1 took the tail Y..
(runs/a20c.chunks):

    combine_a20.py runs/a20.chunks:0:53416 runs/a20b.chunks:53416:Y runs/a20c.chunks:Y:92569
"""
import sys

NCHUNKS = 92569
PI6_1E19 = 28722086297                 # pi_6(10^19), Desfontaines; confirmed by our a(19) run
PI6_2P64 = 48629687343                 # pi_6(2^64), Desfontaines
BOUNDS = ['1e19', '2^64', '2e19', '3e19', '2^65', '4e19', '5e19', '6e19', '7e19', '2^66', '8e19', '9e19', '1e20']
MASK = (1 << 64) - 1


def read(path):
    lines = {}
    dup = 0
    for raw in open(path):
        if not raw.startswith('C '):
            continue
        raw = raw.rstrip('\n')
        f = raw.split()
        if len(f) != 5 + len(BOUNDS) - 1:
            continue                    # partial line at the end of a log being written
        k = int(f[1])
        if k in lines and lines[k] != raw:
            dup += 1
        lines[k] = raw
    return lines, dup


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    logs, segs = {}, []
    for arg in sys.argv[1:]:
        path, lo, hi = arg.rsplit(':', 2)
        if path not in logs:
            logs[path], dups = read(path)
            if dups:
                print(f"CONFLICTING DUPLICATE LINES in {path}: {dups}")
            ks = logs[path]
            print(f"{path}: {len(ks)} chunks ({min(ks) if ks else '-'}..{max(ks) if ks else '-'})")
        segs.append((path, int(lo), int(hi)))
    # cross-check: chunks present in more than one log
    seen, both, same = {}, 0, 0
    for path, ks in logs.items():
        for k, line in ks.items():
            if k in seen:
                both += 1
                same += seen[k] == line
                if seen[k] != line:
                    print(f"MISMATCH chunk {k}:\n  {seen[k]}\n  {line}")
            else:
                seen[k] = line
    print(f"chunks computed in more than one run: {both}, identical: {same}" + ("" if same == both else "  MISMATCH!"))
    merged = {}
    for path, lo, hi in segs:
        for k, line in logs[path].items():
            if lo <= k < hi:
                if k in merged:
                    print(f"segments overlap at chunk {k}")
                merged[k] = line
    missing = [k for k in range(NCHUNKS) if k not in merged]
    if missing:
        print(f"INCOMPLETE: {len(missing)} chunks missing (first {missing[:5]})")
    nb = len(BOUNDS) - 1
    cnt, cks = [0] * nb, [0] * nb
    surv = cand = spsp = 0
    for k in sorted(merged):
        f = merged[k].split()
        surv += int(f[2])
        cand += int(f[3])
        spsp += int(f[4])
        for i, field in enumerate(f[5:]):
            c, h = field.split(':')
            cnt[i] += int(c)
            cks[i] = (cks[i] + int(h, 16)) & MASK
    total = sum(cnt)
    tag = 'PARTIAL' if missing else 'RESULT'
    pi = PI6_1E19
    for i in range(nb):
        pi += cnt[i]
        print(f"{tag} [{BOUNDS[i]}, {BOUNDS[i + 1]}) count {cnt[i]} cks {cks[i]:016x}   pi_6({BOUNDS[i + 1]}) = {pi}")
    print(f"{tag} [1e19, 1e20) count {total} cks {sum(cks) & MASK:016x}  A350826(20)")
    print(f"{tag} survivors {surv} candidates {cand} spsp {spsp} (candidates - count = {cand - total})")
    print(f"{tag} A063501(20) = pi_6(10^20) = {PI6_1E19 + total}")
    if not missing:
        ok = PI6_1E19 + cnt[0] == PI6_2P64
        print(f"check pi_6(2^64) = {PI6_1E19 + cnt[0]} vs Desfontaines {PI6_2P64}: {'OK' if ok else 'MISMATCH'}")
        print(f"check candidates - count = spsp: {'OK' if cand - total == spsp else 'MISMATCH'}")


if __name__ == '__main__':
    main()
