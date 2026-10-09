#!/usr/bin/env python3
"""Keep one construction per Q-class fingerprint.

usage: dedup_q7.py EXT_TSV OUT_TSV      (EXT_TSV: extend1 -o output)
OUT_TSV: representative \t fingerprint \t |G| \t #constructions with this fingerprint
"""
import sys
from collections import defaultdict

cons, order = defaultdict(list), {}
for l in open(sys.argv[1]):
    f = l.rstrip('\n').split('\t')
    if len(f) != 5:
        sys.exit(f'malformed line: {l!r}')
    cons[f[2]].append(f[0])
    order[f[2]] = f[1]
with open(sys.argv[2], 'w') as out:
    for h in sorted(cons):
        out.write(f'{sorted(cons[h])[0]}\t{h}\t{order[h]}\t{len(cons[h])}\n')
print(f'{sum(len(v) for v in cons.values())} constructions, {len(cons)} distinct Q-class fingerprints')
