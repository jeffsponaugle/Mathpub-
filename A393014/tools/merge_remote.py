#!/usr/bin/env python3
"""Pull results from another machine running dbpal campaigns and merge them into results/.

  merge_remote.py user@host[:path]      (default path ~/dbpal)
  merge_remote.py DIR                   (a results directory copied here beforehand)

Copies the remote results/ into results/remote_<host>/ (rsync), then appends to the local
search_log.jsonl every remote log line that is not already present, and to each local
b<B1>_<B2>.tsv every remote result line whose number is not already listed.  Idempotent.
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    target = sys.argv[1]
    if os.path.isdir(target):  # a results directory already copied here (e.g. from a machine that can't be reached)
        stage = host = target.rstrip("/")
    else:
        host, _, path = target.partition(":")
        path = path or "~/dbpal"
        stage = os.path.join(RES, "remote_" + host.split("@")[-1].replace(".", "_"))
        os.makedirs(stage, exist_ok=True)
        subprocess.run(["rsync", "-az", f"{host}:{path}/results/", stage + "/"], check=True)
    tag = os.path.basename(stage)[len("remote_"):] if os.path.isdir(target) else host.split("@")[-1]
    # log lines
    local_log = os.path.join(RES, "search_log.jsonl")
    def key(r):  # identity of a log record, ignoring the host tag added here
        return json.dumps({k: v for k, v in r.items() if k != "host"}, sort_keys=True)

    have = set()
    if os.path.exists(local_log):
        for line in open(local_log):
            if line.strip():
                have.add(key(json.loads(line)))
    added = 0
    rlog = os.path.join(stage, "search_log.jsonl")
    if os.path.exists(rlog):
        with open(local_log, "a") as out:
            for line in open(rlog):
                if not line.strip():
                    continue
                r = json.loads(line)
                if key(r) in have:
                    continue
                have.add(key(r))
                r.setdefault("host", tag)  # keep the original machine of records merged before
                out.write(json.dumps(r) + "\n")
                added += 1
    # result files
    nres = 0
    for fn in sorted(os.listdir(stage)):
        if not (fn.startswith("b") and fn.endswith(".tsv")):
            continue
        lpath = os.path.join(RES, fn)
        known = set()
        if os.path.exists(lpath):
            for line in open(lpath):
                known.add(line.split("\t")[0])
        with open(lpath, "a") as out:
            for line in open(os.path.join(stage, fn)):
                n = line.split("\t")[0]
                if n and n[0].isdigit() and n not in known:
                    out.write(line if line.endswith("\n") else line + "\n")
                    known.add(n)
                    nres += 1
    print(f"merged from {host}: {added} log lines, {nres} result lines")
    return 0


if __name__ == "__main__":
    sys.exit(main())
