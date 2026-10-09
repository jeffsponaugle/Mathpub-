#!/usr/bin/env python3
"""collect.py -- gather search results from every machine, verify new terms, and
regenerate submissions.md.

    ./collect.py             fetch, print a report, rewrite submissions.md
    ./collect.py --verify    also confirm unverified new terms with brute2 -v (direct
                             O(n) evaluation from the definition), spread over VERIFY_HOSTS
    ./collect.py --snapshot results/final_collect.txt
                             rebuild from a saved complete console report (no ssh)

Each machine's raw report is cached in results/raw/; an unreachable machine's cached report
is used instead.  If a machine has neither, outputs are not rewritten (unless --partial).

Coverage of a sequence = union over machines of completed ranges, read from .progress
files ("complete [a, x)"), done_<lo>_<hi> markers, and runs/jobs.log (a job restarted at
s had [first start, s) complete; "end ... exit=0" means the job's range completed).
Only the contiguous part starting at 1 is claimed in submissions.md.
"""
import datetime, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
MACHINES = [                      # name, ssh target (None = this machine), project dir
    ("mac1", None, HERE),
    ("mac2", "jbs@10.1.1.43", "~/src/math/A061955"),
    ("mathb", "jbs@10.1.30.21", "~/src/math/A061955"),
    ("mathd", "jbs@10.1.30.23", "~/src/math/A061955"),
    ("mathg", "jbs@10.1.30.26", "~/src/math/A061955"),
]
VERIFY_HOSTS = ["jbs@10.1.30.21", "jbs@10.1.30.23", "jbs@10.1.30.26"]   # brute2 -v runs here, in parallel
VERIFY_DIR = "~/src/math/A061955"
FAMS = [(29447, "Rn"), (29471, "Ln"), (29495, "Rk"), (29519, "Ld"), (61931, "Rd"), (61955, "Lk")]
DESCR = {"Lk": "left concatenation, digits reversed, zeros kept",
         "Ld": "left concatenation, digits reversed, least significant zeros dropped",
         "Rk": "right concatenation, digits reversed, zeros kept",
         "Rd": "right concatenation, digits reversed, least significant zeros dropped",
         "Ln": "left concatenation, normal digits", "Rn": "right concatenation, normal digits"}
NAME = "_[name]_"                 # replace with the submitter's OEIS name

REPORT = r'''cd {path} 2>/dev/null || exit 0
for d in runs/A[0-9]*; do
  [ -d "$d" ] || continue
  a=${{d#runs/}}
  case $a in *.*) continue;; esac
  echo "@DIR $a"
  for f in "$d"/hits_*.txt "$d"/hits.txt; do [ -f "$f" ] && grep -v '^#' "$f" | sed 's/^/@HIT /'; done
  for f in "$d"/*.progress; do [ -f "$f" ] && sed 's/^/@PROG /' "$f"; done
  [ -f "$d/progress.txt" ] && sed 's/^/@OLDPROG /' "$d/progress.txt"
  ls "$d" | grep '^done_' | sed 's/^/@DONE /'
done
[ -f runs/jobs.log ] && sed 's/^/@LOG /' runs/jobs.log
exit 0
'''


def family(anum):
    for first, tag in FAMS:
        if first <= anum < first + 24:
            return tag, anum - first + 2
    raise ValueError(anum)


def run_report(target, path):
    script = REPORT.format(path=path)
    cmd = ["sh", "-s"] if target is None else ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", target, "sh", "-s"]
    try:
        r = subprocess.run(cmd, input=script, capture_output=True, text=True, timeout=120)
        return r.stdout if r.returncode == 0 else None
    except subprocess.TimeoutExpired:
        return None


def oeis_terms(anum):
    """Published terms (the existing b-file if OEIS has one, else the data section),
    the name, the offset, and the existing b-file's %H line (or None)."""
    path = os.path.join(HERE, "data", "oeis", f"A{anum:06d}.txt")
    t = open(path).read()
    data = "".join(re.findall(rf"^%[STU] A{anum:06d} (.*)$", t, re.M)).replace(" ", "")
    name = re.search(rf"^%N A{anum:06d} (.*)$", t, re.M).group(1)
    offset = int(re.search(rf"^%O A{anum:06d} (-?\d+)", t, re.M).group(1))
    terms = [int(x) for x in data.split(",") if x]
    bpath = os.path.join(HERE, "data", "oeis_bfiles", f"b{anum:06d}.txt")
    hline = None
    if os.path.exists(bpath):
        terms = [int(l.split()[1]) for l in open(bpath) if l.strip() and not l.startswith("#")]
        m = re.search(rf"^%H A{anum:06d} (.*b{anum:06d}\.txt.*)$", t, re.M)
        hline = re.sub(r"<[^>]*>", "", m.group(1)) if m else "existing b-file"
    return terms, name, offset, hline


def data_terms(anum):
    t = open(os.path.join(HERE, "data", "oeis", f"A{anum:06d}.txt")).read()
    data = "".join(re.findall(rf"^%[STU] A{anum:06d} (.*)$", t, re.M)).replace(" ", "")
    return [int(x) for x in data.split(",") if x]


def unfmt(x):
    m = re.fullmatch(r"(?:(\d+)\*)?10\^(\d+)", x)
    if m: return int(m.group(1) or 1) * 10 ** int(m.group(2))
    m = re.fullmatch(r"2\^(\d+)", x)
    return 2 ** int(m.group(1)) if m else int(x)


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def fmt_bound(x):
    for k in range(1, 20):
        if x == 10 ** k: return f"10^{k}"
        if x == 2 ** (k + 20): return f"2^{k + 20}"
    for m in range(2, 10):
        for k in range(1, 19):
            if x == m * 10 ** k: return f"{m}*10^{k}"
    return str(x)


def nice_floor(x):
    """Largest bound <= x of the form m*10^k (m = 1..9) or 2^k, for "No other terms < X"."""
    best = 0
    for k in range(0, 19):
        for m in range(1, 10):
            if m * 10 ** k <= x: best = max(best, m * 10 ** k)
    for k in range(20, 63):
        if 2 ** k <= x and 2 ** k > best and 2 ** k * 0.999 > best: best = max(best, 2 ** k)
    return best


def wrap_data(anum, terms):
    """%S/%T/%U lines (at most 3 of <= 69 chars), the number of terms they hold, and
    whether the list had to be cut (then the full list goes in the b-file)."""
    lines, cur, cnt = [], "", []
    for t in map(str, terms):
        piece = (cur + "," if cur else "") + t
        if len(piece) > 69:
            lines.append(cur + ",")
            cnt.append(len(cur.split(",")))
            cur = t
        else:
            cur = piece
    lines.append(cur)
    cnt.append(len(cur.split(",")))
    tags = ["%S", "%T", "%U"]
    out = [f"{tags[i]} A{anum:06d} {l}" for i, l in enumerate(lines[:3])]
    if len(lines) > 3: out[-1] = out[-1].rstrip(",")
    return out, sum(cnt[:3]), len(lines) > 3


def main():
    verify = "--verify" in sys.argv
    hits, cover, starts, seen = {}, {}, {}, set()
    rawdir = os.path.join(HERE, "results", "raw")
    os.makedirs(rawdir, exist_ok=True)
    missing = []
    if "--snapshot" in sys.argv:              # rebuild from a saved complete console report
        snap = sys.argv[sys.argv.index("--snapshot") + 1]
        for line in open(snap):
            m = re.match(r"A(\d+) \S+\s+contiguous<(\S+)\s+new=\d+(.*)$", line)
            if not m: continue
            a = int(m.group(1)); seen.add(a)
            hits.setdefault(a, set()).update(data_terms(a))
            hits[a].update(int(x) for x in re.findall(r"(\d+)\(", m.group(3)))
            if unfmt(m.group(2)) > 1: cover.setdefault(a, []).append((1, unfmt(m.group(2))))
        sources = []
    else:
        sources = MACHINES
    for name, target, path in sources:
        out = run_report(target, path)
        cache = os.path.join(rawdir, f"{name}.txt")
        if out is not None:
            open(cache, "w").write(out)
        elif os.path.exists(cache):
            print(f"warning: {name} unreachable, using its cached report from "
                  f"{datetime.datetime.fromtimestamp(os.path.getmtime(cache)):%Y-%m-%d %H:%M}", file=sys.stderr)
            out = open(cache).read()
        else:
            print(f"warning: {name} unreachable and no cached report", file=sys.stderr)
            missing.append(name)
            continue
        cur = None
        for line in out.splitlines():
            tag, _, rest = line.partition(" ")
            if tag == "@DIR":
                cur = int(rest[1:]); seen.add(cur)
            elif tag == "@HIT":
                hits.setdefault(cur, set()).add(int(rest))
            elif tag == "@PROG":
                m = re.search(r"complete \[(\d+), (\d+)\)", rest)
                if m: cover.setdefault(cur, []).append((int(m.group(1)), int(m.group(2))))
            elif tag == "@OLDPROG":
                cover.setdefault(cur, []).append((1, int(rest)))
            elif tag == "@DONE":
                lo, hi = map(int, rest.split("_")[1:3])
                cover.setdefault(cur, []).append((lo, hi))
            elif tag == "@LOG":
                m = re.search(r" (start|end  ) A(\d+) \[(\d+), (\d+)\)(?: exit=(\d+))?", rest)
                if not m: continue
                kind, a, s, hi = m.group(1).strip(), int(m.group(2)), int(m.group(3)), int(m.group(4))
                if kind == "start":
                    starts.setdefault((name, a, hi), []).append(s)
                elif m.group(5) == "0":
                    cover.setdefault(a, []).append((s, hi))
    for (name, a, hi), ss in starts.items():
        s0 = min(ss)
        for s in ss:
            if s > s0: cover.setdefault(a, []).append((s0, s))

    vpath = os.path.join(HERE, "results", "verified.json")
    os.makedirs(os.path.dirname(vpath), exist_ok=True)
    verified = json.load(open(vpath)) if os.path.exists(vpath) else {}

    report, problems = [], []
    for a in sorted(seen | set(hits)):
        tag, b = family(a)
        known, name, offset, hline = oeis_terms(a)
        iv = merge(cover.get(a, []))
        contig = iv[0][1] if iv and iv[0][0] <= 1 else 0
        h = sorted(hits.get(a, set()))
        kmax = max(known) if known else 0
        new = [x for x in h if x not in known]
        for x in new:
            if x <= kmax: problems.append(f"A{a:06d}: hit {x} below the largest published term but not published")
        for x in known:
            if x < contig and x not in h: problems.append(f"A{a:06d}: published term {x} not found inside covered range")
        report.append(dict(a=a, tag=tag, b=b, name=name, known=known, new=new, contig=contig, iv=iv,
                           offset=offset, hline=hline, hits=h))

    if verify:                    # one worker per host, largest n first
        import threading
        todo = sorted(((x, f"A{r['a']:06d}", r["tag"], r["b"]) for r in report for x in r["new"]
                       if str(x) not in verified.get(f"A{r['a']:06d}", {})), reverse=True)
        lock = threading.Lock()
        def worker(host):
            while True:
                with lock:
                    if not todo: return
                    x, key, tag, b = todo.pop(0)
                cmd = f"cd {VERIFY_DIR} && ./brute2 {tag} {b} -v {x}"
                out = subprocess.run(["ssh", "-o", "BatchMode=yes", host, cmd], capture_output=True, text=True).stdout.split()
                if not (len(out) == 2 and out[0] == str(x)):
                    print(f"verify {key} {x}: no answer from {host}, will retry later")
                    continue
                with lock:
                    verified.setdefault(key, {})[str(x)] = int(out[1])
                    print(f"verify {key} {x}: residue {out[1]} ({host.split('@')[1]})", flush=True)
                    json.dump(verified, open(vpath, "w"), indent=1, sort_keys=True)
        th = [threading.Thread(target=worker, args=(h,)) for h in VERIFY_HOSTS]
        for t in th: t.start()
        for t in th: t.join()

    # ---- console report
    for r in report:
        key = f"A{r['a']:06d}"
        st = [("ok" if verified.get(key, {}).get(str(x)) == 0 else "BAD" if str(x) in verified.get(key, {}) else "?") for x in r["new"]]
        cov = " ".join(f"[{fmt_bound(lo)},{fmt_bound(hi)})" for lo, hi in r["iv"])
        print(f"{key} {r['tag']}:{r['b']:<2} contiguous<{fmt_bound(r['contig']):<8} new={len(r['new'])} "
              + " ".join(f"{x}({s})" for x, s in zip(r["new"], st)) + (f"   ranges {cov}" if len(r["iv"]) > 1 else ""))
    for p in problems: print("PROBLEM", p)

    if missing and "--partial" not in sys.argv:
        print(f"not rewriting submissions.md or bfiles/: no data from {', '.join(missing)} "
              f"(bring them up, or use --snapshot / --partial)", file=sys.stderr)
        return

    # ---- b-files: all terms below the contiguous bound, from the search itself
    bdir = os.path.join(HERE, "bfiles")
    os.makedirs(bdir, exist_ok=True)
    for r in report:
        bound = nice_floor(r["contig"]) if r["contig"] > 1 else 0
        terms = sorted(set(r["hits"]) | set(r["known"]))
        r["bterms"] = [x for x in terms if x < bound]
        if not bound or problems: continue
        key = f"A{r['a']:06d}"
        with open(os.path.join(bdir, f"b{key[1:]}.txt"), "w") as fh:
            fh.write(f"# {key}: all terms below {fmt_bound(bound)}.\n")
            fh.write(f"# Every n < {fmt_bound(bound)} was tested; terms beyond the previously published ones were\n")
            fh.write(f"# each confirmed by direct evaluation of the concatenation mod n. {datetime.date.today():%b %d %Y}.\n")
            for i, x in enumerate(r["bterms"]):
                fh.write(f"{r['offset'] + i} {x}\n")

    # ---- submissions.md
    today = datetime.date.today().strftime("%b %d %Y")
    L = ["# Draft OEIS submissions", "",
         "Generated by `collect.py` from the search results on all machines; do not edit by hand.",
         f"Last update: {datetime.datetime.now():%Y-%m-%d %H:%M}.", "",
         "New terms come from `concat` (O(b log^2 n) evaluator, see README.md) and are confirmed by",
         "direct O(n) evaluation from the definition (`brute2 FAM B -v N`, an independent word-by-word",
         "implementation): **verified** = residue 0 confirmed, **unverified** = not yet run.",
         "\"No other terms < X\" uses only the contiguous range [1, X) that every n was tested on.",
         f"Author placeholder: `{NAME}` -- replace with your OEIS name before submitting.", ""]
    if problems:
        L += ["## Problems (must be resolved before submitting)", ""] + [f"- {p}" for p in problems] + [""]
    with_new = [r for r in report if r["new"]]
    extended = [r for r in report if not r["new"] and r["contig"] >= 10 ** 8]
    L += ["## Summary", "", "| sequence | family | base | searched (contiguous) | new terms |", "|---|---|---|---|---|"]
    for r in sorted(with_new, key=lambda r: (r["b"] != 2, r["a"])) + extended:
        srch = f"< {fmt_bound(nice_floor(r['contig']))}" if r["contig"] > 1 else "in progress"
        L.append(f"| [A{r['a']:06d}](https://oeis.org/A{r['a']:06d}) | {r['tag']} | {r['b']} | {srch} | {len(r['new'])} |")
    L.append("")
    for r in sorted(with_new, key=lambda r: (r["b"] != 2, r["a"])):
        key = f"A{r['a']:06d}"
        k0 = len(r["known"])
        L += ["---", "", f"## {key} — base {r['b']}, {DESCR[r['tag']]}", "", r["name"], "",
              "| term | value | status |", "|---|---|---|"]
        for i, x in enumerate(r["new"]):
            v = verified.get(key, {}).get(str(x), "missing")
            st = "verified" if v == 0 else "unverified" if v == "missing" else f"FAILED (residue {v})"
            L.append(f"| a({k0 + i + 1}) | {x} | {st} |")
        lines, shown, overflow = wrap_data(r["a"], r["known"] + r["new"])
        bound = nice_floor(r["contig"]) if r["contig"] > 1 else 0
        L += ["", f"Searched: every n < {fmt_bound(bound)}." if bound else "Searched: first range still in progress (terms so far).",
              "", "Proposed edit:", "```"] + lines
        if bound:
            L.append(f"%C {key} No other terms < {fmt_bound(bound)}. - {NAME}, {today}")
        last = k0 + len(r["new"])
        if bound and (overflow or r["hline"]):
            n_b = r["offset"] + len(r["bterms"]) - 1
            credit = f"; terms 1..{len(r['known'])} from Lars Blomberg" if r["hline"] else ""
            L.append(f'%H {key} {NAME}, <a href="/{key}/b{key[1:]}.txt">Table of n, a(n) for n = {r["offset"]}..{n_b}</a> '
                     f'(all terms below {fmt_bound(bound)}{credit})')
        L.append(f"%E {key} a({k0 + 1}){'-a(%d)' % last if last > k0 + 1 else ''} from {NAME}, {today}")
        L.append("```")
        if overflow:
            L.append(f"Data shows the first {shown} of {len(r['known']) + len(r['new'])} terms; the full list is "
                     f"`bfiles/b{key[1:]}.txt` (upload it with the edit).")
        elif r["hline"]:
            L.append(f"OEIS already has a b-file ({r['hline']}); replace it with `bfiles/b{key[1:]}.txt`.")
        else:
            L.append(f"A b-file is optional here; `bfiles/b{key[1:]}.txt` is available.")
        L.append("")
    if extended:
        L += ["---", "", "## Search limits extended, no new terms", "",
              "Proposed comment for each: `%C Axxxxxx No other terms < X.`", ""]
        for r in extended:
            L.append(f"- A{r['a']:06d} (base {r['b']}, {r['tag']}): no other terms < {fmt_bound(nice_floor(r['contig']))}")
        L.append("")
    open(os.path.join(HERE, "submissions.md"), "w").write("\n".join(L))
    print(f"submissions.md: {len(with_new)} sequences with new terms, {len(extended)} extended without new terms")


if __name__ == "__main__":
    main()
