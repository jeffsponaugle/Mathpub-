#!/usr/bin/env python3
"""Resumable search campaigns.

  campaign.py run --bases 3,10 --below 10^40 [--primes] [--engine gpu|cpu] [--max-seconds S]
      search every length (in dbpal's enumerated base) below the bound that is not yet recorded
      as done in results/search_log.jsonl; results go to results/b<B1>_<B2>.tsv
  campaign.py status
      searched bound and term counts per base pair

Lengths are run one at a time, so an interrupted campaign resumes where it stopped.
"""
import argparse
import functools
import json
import os
import subprocess
import sys
import time

print = functools.partial(print, flush=True)  # noqa: A001 (unbuffered progress for logs)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
LOG = os.path.join(RES, "search_log.jsonl")
DBPAL = os.path.join(ROOT, "dbpal")
DBPAL4 = os.path.join(ROOT, "dbpal4")  # NLIMB=4 build (`make dbpal4`), used for lengths beyond the 3-limb range


def parse_num(s):
    if "^" in s:
        b, e = s.split("^")
        return int(b) ** int(e)
    if "e" in s.lower():
        m, e = s.lower().split("e")
        return int(m) * 10 ** int(e)
    return int(s)


def roles(b1, b2):
    """(P, Q) exactly as dbpal chooses them."""
    pow2 = lambda b: b & (b - 1) == 0
    if pow2(b1) != pow2(b2):
        return (b2, b1) if pow2(b1) else (b1, b2)
    return max(b1, b2), min(b1, b2)


def ndigits(n, b):
    d = 0
    while n:
        n //= b
        d += 1
    return max(d, 1)


def log_entries():
    """Records of results/search_log.jsonl plus logs synced from other machines
    (results/remote_*/search_log.jsonl), de-duplicated."""
    paths = [LOG]
    if os.path.isdir(RES):
        paths += [os.path.join(RES, d, "search_log.jsonl") for d in sorted(os.listdir(RES)) if d.startswith("remote_")]
    out, seen = [], set()
    for path in paths:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except Exception:
                    continue
                key = json.dumps({k: v for k, v in r.items() if k != "host"}, sort_keys=True)
                if key not in seen:
                    seen.add(key)
                    out.append(r)
    return out


def complete_lengths(records):
    """Lengths searched completely: a 'done'/'skipped' record, or 'part' records i/n for every i."""
    s, parts = set(), {}
    for r in records:
        if r["status"] in ("done", "skipped"):
            s.add(r["L"])
        elif r["status"] == "part" and "plan" in r:  # parts only combine when they used the same plan
            i, n = map(int, r["part"].split("/"))
            parts.setdefault((r["L"], n, r["plan"]), set()).add(i)
    for (L, n, plan), got in parts.items():
        if len(got) == n:
            s.add(L)
    return s


CLAIM_HOURS = 12  # a part claimed by another machine this recently is left to it
ME = __import__("socket").gethostname().split(".")[0]


def claimed_elsewhere(b1, b2, L, primes, part):
    """Another machine's recent claim on this part (machines take split lengths from opposite ends,
    and a claim keeps them from both running the part where they meet).  Returns (host, time) or None."""
    now = time.time()
    for r in log_entries():
        if (r.get("status") == "claim" and sorted(r["bases"]) == sorted([b1, b2]) and r["L"] == L
                and r.get("part") == part and r["primesOnly"] == primes and r.get("by") != ME):
            t = time.mktime(time.strptime(r["time"], "%Y-%m-%dT%H:%M:%S"))
            if now - t < CLAIM_HOURS * 3600:
                return r.get("by"), r["time"]
    return None


def claim(b1, b2, P, Q, L, primes, part):
    with open(LOG, "a") as f:
        f.write(json.dumps({"bases": [b1, b2], "P": P, "Q": Q, "L": L, "primesOnly": primes, "status": "claim",
                            "part": part, "by": ME, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}) + "\n")


def done_parts(b1, b2, primes, part):
    """Lengths whose part 'i/n' is logged."""
    P, Q = roles(b1, b2)
    return {r["L"] for r in log_entries() if sorted(r["bases"]) == sorted([b1, b2]) and r["P"] == P
            and r["status"] == "part" and r.get("part") == part and (primes or not r["primesOnly"])}


def done_lengths(b1, b2, primes):
    """Lengths fully searched; an all-palindromes search also covers the prime search."""
    P, Q = roles(b1, b2)
    rs = [r for r in log_entries() if sorted(r["bases"]) == sorted([b1, b2]) and r["P"] == P
          and (primes or not r["primesOnly"])]
    return P, Q, complete_lengths(rs)


def assigned_elsewhere(b1, b2):
    """results/remote_pairs.txt: lines 'B1,B2 host' -- pairs searched on another machine."""
    path = os.path.join(RES, "remote_pairs.txt")
    if not os.path.exists(path):
        return None
    for line in open(path):
        p = line.split()
        if p and not p[0].startswith("#") and sorted(map(int, p[0].split(","))) == sorted([b1, b2]):
            return p[1] if len(p) > 1 else "elsewhere"
    return None


def predict_seconds(took, L):
    """Cost of length L from this run's measured lengths.  Odd and even lengths can differ a lot, so
    scale the last length of the same parity by the growth seen over two lengths; 0 if unknown."""
    same = sorted(l for l in took if l % 2 == L % 2 and l < L)
    if len(same) >= 2:
        return took[same[-1]] * max(1.0, took[same[-1]] / took[same[-2]]) ** ((L - same[-1]) / (same[-1] - same[-2]))
    ls = sorted(took)
    if len(ls) >= 2:  # one length of each parity: assume the growth between them applies again
        return took[ls[-1]] * max(1.0, took[ls[-1]] / took[ls[-2]]) ** (L - ls[-1])
    return 0.0


def cmd_run(a):
    b1, b2 = map(int, a.bases.split(","))
    host = assigned_elsewhere(b1, b2)
    if host:
        print(f"bases ({b1},{b2}): assigned to {host}, skipping here")
        return 0
    bound = parse_num(a.below)
    P, Q, done = done_lengths(b1, b2, a.primes)
    if a.part:
        done |= done_parts(b1, b2, a.primes, a.part)
    Lmax = ndigits(bound - 1, P)
    todo = [L for L in range(max(1, a.Lmin), Lmax + 1) if L not in done]
    out = os.path.join(RES, f"b{min(b1, b2)}_{max(b1, b2)}.tsv")
    print(f"bases ({b1},{b2}): P={P} Q={Q}; lengths 1..{Lmax}; {len(todo)} to do: {todo[:12]}{'...' if len(todo) > 12 else ''}")
    t0 = time.time()
    took = {}  # L -> seconds for lengths searched in this run (> 5 s), to predict the next length's cost
    for L in todo:
        spent = time.time() - t0
        if a.max_seconds and spent > a.max_seconds:
            print("time budget exhausted")
            break
        guess = predict_seconds(took, L)
        if a.max_seconds and guess and spent + guess > a.max_seconds:
            print(f"time budget: L={L} would take ~{guess:.0f}s, stopping")
            break
        if a.part:
            c = claimed_elsewhere(b1, b2, L, a.primes, a.part)
            if c:
                print(f"  L={L}: part {a.part} claimed by {c[0]} at {c[1]}, skipping")
                continue
            claim(b1, b2, P, Q, L, a.primes, a.part)
        cmd = [DBPAL, "search", "-b", f"{b1},{b2}", "--Lmin", str(L), "--Lmax", str(L), "--out", out, "--log", LOG]
        if a.primes:
            cmd.append("--primes")
        if a.engine:
            cmd += ["--engine", a.engine]
        if a.part:
            cmd += ["--part", a.part]
        cmd += a.extra
        t1 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0 and "too large for 3 limbs" in r.stderr and os.path.exists(DBPAL4):
            r = subprocess.run([DBPAL4] + cmd[1:], capture_output=True, text=True)
        hits = [l for l in r.stdout.splitlines() if l and l[0].isdigit()]
        dt = time.time() - t1
        if dt > 5:
            took[L] = dt
        print(f"  L={L}: {len(hits)} hits in {dt:.1f}s" + ("" if r.returncode == 0 else f"  (exit {r.returncode})"))
        for h in hits:
            print("     ", h.split("\t")[0], h.split("\t")[1])
        if r.returncode != 0:
            print(r.stderr[-2000:])
            return 1
    return 0


def cmd_status(a):
    pairs = {}
    for r in log_entries():
        key = tuple(sorted(r["bases"]))
        pairs.setdefault(key, []).append(r)
    for key, rs in sorted(pairs.items()):
        P = rs[0]["P"]
        for primes in (False, True):
            Ls = complete_lengths([r for r in rs if r["primesOnly"] == primes or not r["primesOnly"]])
            L = 0
            while L + 1 in Ls:
                L += 1
            if L == 0:
                continue
            path = os.path.join(RES, f"b{key[0]}_{key[1]}.tsv")
            nt = np_ = 0
            if os.path.exists(path):
                for line in open(path):
                    p = line.split("\t")
                    if len(p) > 1 and int(p[0]) < P ** L:
                        nt += 1
                        np_ += p[1] == "1"
            secs = sum(r.get("seconds", 0) + r.get("tableSeconds", 0) for r in rs)
            tag = "primes only" if primes else "all"
            print(f"bases {key[0]},{key[1]} ({tag}): complete below {P}^{L} ~ {float(P**L):.3g}; "
                  f"{nt} double palindromes ({np_} prime); {secs/3600:.2f} h of search logged")
            if not primes:
                break
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("run")
    r.add_argument("--bases", required=True)
    r.add_argument("--below", required=True)
    r.add_argument("--primes", action="store_true")
    r.add_argument("--engine", default=None)
    r.add_argument("--max-seconds", type=float, default=0)
    r.add_argument("--part", default=None, help="i/n: this machine searches part i of n of every length")
    r.add_argument("--Lmin", type=int, default=1, help="only lengths >= Lmin (in the enumerated base)")
    r.add_argument("extra", nargs="*", help="extra dbpal options after --")
    sub.add_parser("status")
    a = ap.parse_args()
    if a.cmd == "run":
        return cmd_run(a)
    if a.cmd == "status":
        return cmd_status(a)
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
