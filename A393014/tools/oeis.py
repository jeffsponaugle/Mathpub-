#!/usr/bin/env python3
"""OEIS helper for the double-base palindrome family.

  oeis.py catalog                      list the related sequences (bases, primes-only?, #terms, max term)
  oeis.py compare RESULTS.tsv          compare dbpal results with every matching OEIS sequence
  oeis.py show A393014                 print the cached entry

Sequence entries and b-files are cached in tools/cache/ (delete it to refresh).
"""
import json
import os
import re
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache")

# The family around A393014: (bases, primes only) -> A-number.
# Built from the cross references of A393014, A259385 and A046473 (names checked by `catalog`).
SEEDS = ["A393014", "A259385", "A046473"]


def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "dbpal-oeis-check/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")


def entry(anum, refresh=False):
    """Return dict(name, terms, bfile, xrefs) for an A-number (cached)."""
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, anum + ".json")
    if os.path.exists(path) and not refresh:
        with open(path) as f:
            return json.load(f)
    txt = _get(f"https://oeis.org/search?q=id:{anum}&fmt=text")
    name, data, xrefs = "", "", []
    for line in txt.splitlines():
        m = re.match(r"%(\w) (A\d+) ?(.*)", line)
        if not m or m.group(2) != anum:
            continue
        tag, rest = m.group(1), m.group(3)
        if tag == "N":
            name = rest
        elif tag in "STU":
            data += rest
        elif tag == "Y":
            xrefs += re.findall(r"A\d{6}", rest)
    terms = [int(x) for x in data.split(",") if x.strip()]
    bfile = []
    try:
        num = anum[1:]
        bt = _get(f"https://oeis.org/{anum}/b{num}.txt")
        for line in bt.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                bfile.append(int(parts[1]))
    except Exception:
        pass
    e = {"anum": anum, "name": name, "terms": terms, "bfile": bfile, "xrefs": sorted(set(xrefs))}
    with open(path, "w") as f:
        json.dump(e, f)
    time.sleep(0.3)
    return e


def all_terms(e):
    return e["bfile"] if len(e["bfile"]) >= len(e["terms"]) else e["terms"]


NAME_PATTERNS = [
    (re.compile(r"^Primes that are palindromic in (?:both )?bases? (\d+) and (?:base )?(\d+)", re.I), True),
    (re.compile(r"^Palindromic primes in bases (\d+) and (\d+)", re.I), True),
    (re.compile(r"^Palindromic numbers in bases (\d+) and (\d+)", re.I), False),
    (re.compile(r"^Numbers that are palindromic in bases (\d+) and (\d+)", re.I), False),
    (re.compile(r"^Palindromic in bases (\d+) and (\d+)", re.I), False),
    (re.compile(r"^Base[- ]10 numbers that are palindromic in bases (\d+) and (\d+)", re.I), False),
    (re.compile(r"^Numbers that are palindromes in bases (\d+) and (\d+)", re.I), False),
]


def classify(name):
    for pat, primes in NAME_PATTERNS:
        m = pat.search(name)
        if m:
            b = tuple(sorted((int(m.group(1)), int(m.group(2)))))
            return b, primes
    return None, None


def family(refresh=False):
    """Collect (bases, primes) -> entry for the palindrome family reachable from the seeds."""
    seen, fam, todo = set(), {}, list(SEEDS)
    while todo:
        a = todo.pop(0)
        if a in seen:
            continue
        seen.add(a)
        try:
            e = entry(a, refresh)
        except Exception as ex:
            print(f"  ({a}: {ex})", file=sys.stderr)
            continue
        b, pr = classify(e["name"])
        if b is None:
            continue
        fam[(b, pr)] = e
        # follow cross references only from family members
        for x in e["xrefs"]:
            if x not in seen:
                todo.append(x)
    return fam


def load_results(path):
    """dbpal TSV: n, isprime, b1, b2, repr1, repr2"""
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or not p[0].isdigit():
                continue
            b = tuple(sorted((int(p[2]), int(p[3]))))
            out.setdefault(b, {})[int(p[0])] = p[1] == "1"
    return out


def searched_bound(logpath, bases):
    """Largest N such that every base-P length below it is recorded as done/skipped in the log."""
    if not os.path.exists(logpath):
        return None
    sys.path.insert(0, HERE)
    from campaign import complete_lengths  # only finished lengths count: not parts in progress or claims
    recs = {}
    with open(logpath) as f:
        for line in f:
            try:
                r = json.loads(line)
            except Exception:
                continue
            if tuple(sorted(r["bases"])) != bases:
                continue
            recs.setdefault((r["P"], r["primesOnly"]), []).append(r)
    done = {key: complete_lengths(rs) for key, rs in recs.items()}
    res = {}
    for (P, primes), Ls in done.items():
        L = 1
        while L in Ls:
            L += 1
        res[primes] = (P, L - 1, P ** (L - 1))
    return res


def cmd_catalog(args):
    fam = family("--refresh" in args)
    print(f"{'A-number':9} {'bases':8} {'primes':6} {'#terms':>6}  {'largest known':>45}  name")
    for (b, pr), e in sorted(fam.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        t = all_terms(e)
        mx = str(max(t)) if t else "-"
        print(f"{e['anum']:9} {b[0]:>2},{b[1]:<5} {'yes' if pr else 'no':6} {len(t):6}  {mx:>45}  {e['name'][:70]}")


def cmd_compare(args):
    path = args[0]
    logpath = os.path.join(os.path.dirname(os.path.abspath(path)), "search_log.jsonl")
    res = load_results(path)
    fam = family()
    rc = 0
    for b, found in sorted(res.items()):
        bounds = searched_bound(logpath, b) or {}
        for primes in (False, True):
            key = (b, primes)
            mine = sorted(n for n, isp in found.items() if isp or not primes)
            bnd = bounds.get(False) or (bounds.get(True) if primes else None)
            if key not in fam:
                print(f"bases {b} primes={primes}: no OEIS sequence found; {len(mine)} terms in results")
                continue
            e = fam[key]
            known = sorted(set(all_terms(e)))
            if bnd:
                P, L, lim = bnd
                mine_b = [n for n in mine if n < lim]
                known_b = [n for n in known if n < lim]
                missing = sorted(set(known_b) - set(mine_b))
                extra = sorted(set(mine_b) - set(known_b))
                # terms the OEIS lacks: beyond its largest listed term
                newer = [n for n in extra if not known or n > max(known)]
                gaps = [n for n in extra if known and n <= max(known)]
                status = "OK" if not missing else "MISSING TERMS"
                if missing:
                    rc = 1
                print(f"{e['anum']} bases {b[0]},{b[1]}{' primes' if primes else ''}: searched all N < {P}^{L} "
                      f"(~{float(lim):.3g}); OEIS {len(known)} terms, {len(known_b)} below bound; "
                      f"mine {len(mine_b)}  [{status}]")
                for n in missing:
                    print(f"    missing (in OEIS, not found): {n}")
                for n in gaps:
                    print(f"    not in OEIS data (gap):       {n}")
                for n in newer:
                    print(f"    NEW (beyond OEIS data):       {n}")
            else:
                print(f"{e['anum']} bases {b}: no completed search log")
    return rc


def cmd_show(args):
    e = entry(args[0], "--refresh" in args)
    print(e["anum"], e["name"])
    t = all_terms(e)
    print(f"{len(t)} terms; last: {t[-3:] if t else []}")
    print("xrefs:", " ".join(e["xrefs"]))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    cmd, args = sys.argv[1], sys.argv[2:]
    if cmd == "catalog":
        return cmd_catalog(args)
    if cmd == "compare":
        return cmd_compare(args)
    if cmd == "show":
        return cmd_show(args)
    print(__doc__)
    return 1


if __name__ == "__main__":
    sys.exit(main() or 0)
