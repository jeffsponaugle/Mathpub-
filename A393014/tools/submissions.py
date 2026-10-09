#!/usr/bin/env python3
"""Write submissions.md: concrete OEIS edit proposals for every sequence the searches touched.

  submissions.py [--who "Jeff Sponaugle"] [--date "Oct 01 2026"] [--refresh]

Uses results/ (terms, certified bounds from the search log, b-files written by tools/bfile.py),
the OEIS entries and their revision histories (cached in tools/cache/full/, --refresh to re-fetch),
and tools/certify.py for primality certificates of new prime terms.  Edits already made on the OEIS
(published, or pending review) are listed in their own section and not proposed again.  Re-run
after more searching or after OEIS edits; the file is regenerated.
"""
import argparse
import datetime
import html
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import certify  # noqa: E402
import oeis  # noqa: E402
from campaign import roles  # noqa: E402
from report import bound_for  # noqa: E402
from verify import to_str  # noqa: E402

RES = os.path.join(ROOT, "results")
FULL = os.path.join(HERE, "cache", "full")
DATA_LIMIT = 260  # practical length of the OEIS DATA field (%S/%T/%U)


def entry_text(anum, refresh):
    os.makedirs(FULL, exist_ok=True)
    path = os.path.join(FULL, anum + ".txt")
    if refresh or not os.path.exists(path):
        open(path, "w").write(oeis._get(f"https://oeis.org/search?q=id:{anum}&fmt=text"))
        time.sleep(0.4)
    return open(path).read()


def parse(anum, txt):
    lines = [l for l in txt.splitlines() if re.match(r"%\w " + anum + r"\b", l)]
    get = lambda tag: [l[11:] for l in lines if l.startswith("%" + tag)]
    data = "".join(get("S") + get("T") + get("U"))
    rev = re.search(r"#(\d+) (\w{3} \d{2} \d{4})", (get("I") or [""])[0])
    return {
        "name": (get("N") or [""])[0],
        "data": data,
        "terms": [int(x) for x in data.split(",") if x.strip()],
        "keywords": (get("K") or [""])[0],
        "offset": int((get("O") or ["1"])[0].split(",")[0]),
        "bfile": next((h for h in get("H") if "Table of n" in h), None),
        "bound_comments": [c for c in get("C") if re.search(r">|exist", c)],
        "rev": int(rev.group(1)) if rev else 0,  # published revision
        "revdate": rev.group(2) if rev else "",
    }


def history_html(anum, refresh):
    path = os.path.join(FULL, anum + ".history.html")
    if refresh or not os.path.exists(path):
        open(path, "w").write(oeis._get(f"https://oeis.org/history?seq={anum}"))
        time.sleep(0.4)
    return open(path).read()


def _text(h):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", h))).strip()


def revisions(page):
    """Revisions on an OEIS history page, newest first: rev, who, date, and per changed section
    (DATA, COMMENTS, ..., STATUS) the inserted and deleted texts; plus discussion notes."""
    out = []
    for b in page.split('<div class="revbar">')[1:]:
        m = re.search(r'v=(\d+)">#\d+</a>\s*by\s*<a[^>]*>([^<]+)</a>\s*at\s*\w{3} (\w{3} \d{2}) [\d:]+ \w+ (\d{4})', b)
        if not m:
            continue
        secs = {}
        for name, body in re.findall(r"<div class=sectname>(\w+)</div>(.*?)(?=<div class=sectname>|<div class=discussbar>|$)", b, re.S):
            secs[name] = ([_text(x) for x in re.findall(r"<ins>(.*?)</ins>", body, re.S)],
                          [_text(x) for x in re.findall(r"<del>(.*?)</del>", body, re.S)])
        notes = [(u.strip(), _text(t).lstrip(": ")) for u, t in re.findall(r"<span class=user>([^<]+)</span>(.*?)</pre>", b, re.S)]
        out.append(dict(rev=int(m.group(1)), who=m.group(2).strip(), date=f"{m.group(3)} {m.group(4)}", secs=secs, notes=notes))
    return out


def oeis_state(info, page, who):
    """The entry's edits beyond the published text (pending review), and the edits by `who`."""
    revs = revisions(page)
    pending = [r for r in revs if r["rev"] > info["rev"]]
    ins = lambda rs, sec: [t for r in rs for t in r["secs"].get(sec, ([], []))[0] if t.strip(", ")]
    status = next((r["secs"]["STATUS"][0][-1] for r in revs if r["secs"].get("STATUS", ([], []))[0]), "approved")
    mine = [r for r in revs if r["who"] == who and set(r["secs"]) - {"STATUS"}]
    return {
        "pending": bool(pending),
        "status": status if pending else "approved",
        "pending_terms": [int(x) for t in ins(pending, "DATA") for x in re.findall(r"\d+", t)],
        "pending_comments": ins(pending, "COMMENTS"),
        "pending_links": ins(pending, "LINKS"),
        "mine": [dict(r, published=r["rev"] <= info["rev"]) for r in mine],
        "notes": [(u, t, r["date"]) for r in revs if r["rev"] > min([x["rev"] for x in mine] or [10**9]) for u, t in r["notes"] if u != who],
    }


def edit_summary(r):
    """Short description of one revision's content: the extension (new terms), comments and b-file."""
    parts = []
    for t in r["secs"].get("EXTENSIONS", ([], []))[0]:
        parts.append(t.split(" from ")[0])
    for t in r["secs"].get("COMMENTS", ([], []))[0]:
        parts.append("comment `" + t.split(" - ")[0] + "`")
    for t in r["secs"].get("LINKS", ([], []))[0]:
        m = re.search(r"for n = (\d+\.\.\d+)", t)
        parts.append(f"b-file {m.group(1)}" if m else "link")
    return "; ".join(parts) or ", ".join(sorted(set(r["secs"]) - {"STATUS"})).lower()


def sci(n):
    m, e = f"{float(n):.2e}".split("e")
    return f"{m}*10^{int(e)}"


def nice_bound(P, L):
    """P^L as text, plus a round power of 10 not above it."""
    lim = P ** L
    k = len(str(lim)) - 1
    return f"{P}^{L}" if P != 10 else f"10^{L}", f"10^{k}", lim


def credit_suffix(old_bfile_line, n_old):
    """Previous contributors of a b-file, for the '(terms 1..m from ...)' part of the new %H line.
    The old line's own credits are kept, flattened into a ';'-separated list (no nested parentheses)."""
    if not old_bfile_line:
        return ""
    who = old_bfile_line.split(",")[0].strip()
    extra = re.search(r"</a>\s*,?\s*\(?(.*?)\)?\s*$", old_bfile_line)
    extra = extra.group(1).strip().rstrip(".") if extra else ""
    m = re.fullmatch(r"terms < (\S+)", extra)
    if m:  # e.g. "(terms < 10^18)": fold into the range
        return f"terms 1..{n_old} (those < {m.group(1)}) from {who}"
    if extra:
        return f"terms 1..{n_old} from {who}; {extra}"
    return f"terms 1..{n_old} from {who}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--who", default="Jeff Sponaugle")
    ap.add_argument("--date", default=datetime.date.today().strftime("%b %d %Y"))
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()
    who, date = a.who, a.date
    sig = f"_{who}_, {date}"

    fam = oeis.family()
    items = []  # one per OEIS sequence touched
    newseq = []
    for fn in sorted(os.listdir(RES)):
        if not (fn.startswith("b") and fn.endswith(".tsv")):
            continue
        b1, b2 = map(int, fn[1:-4].split("_"))
        bd = bound_for(b1, b2)
        if not bd or bd[1] == 0:
            continue
        P, L, lim, hours = bd
        terms = {}
        for line in open(os.path.join(RES, fn)):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].isdigit() and int(p[0]) < lim:
                terms[int(p[0])] = p[1] == "1"
        for primes in (False, True):
            e = fam.get(((b1, b2), primes))
            mine = sorted(n for n, isp in terms.items() if isp or not primes)
            if not e:
                if primes and len([x for x in mine if x > max(b1, b2)]) >= 2:
                    newseq.append((b1, b2, P, L, mine))
                continue
            if a.refresh:  # terms and b-file as published now
                e = oeis.entry(e["anum"], refresh=True)
            info = parse(e["anum"], entry_text(e["anum"], a.refresh))
            state = oeis_state(info, history_html(e["anum"], a.refresh), who)
            # terms in a pending edit count as known: they are submitted already
            published = set(oeis.all_terms(e))
            known = sorted(published | set(state["pending_terms"]))
            kset = set(known)
            new = [n for n in mine if n not in kset]
            pending_new = [n for n in mine if n in kset and n not in published]
            missing = [n for n in known if n < lim and n not in set(mine)]
            items.append(dict(anum=e["anum"], b1=b1, b2=b2, P=P, L=L, lim=lim, hours=hours, primes=primes,
                              info=info, state=state, known=known, mine=mine, new=new, missing=missing,
                              pending_new=pending_new))

    out = []
    w = out.append
    w("# OEIS submission candidates")
    w("")
    w(f"Generated {date} by `tools/submissions.py` from the dbpal results (`results/`), the search log")
    w("(`results/search_log.jsonl`, which certifies each bound length by length) and the current OEIS entries.")
    w(f"Attribution below uses _{who}_; change it if your OEIS user name differs. Lines are in OEIS internal")
    w("format (`%S` data, `%C` comment, `%H` link, `%E` extension) so they can be pasted into the edit form.")
    w("")
    w("General notes for every submission:")
    w("")
    w("* **Completeness.** Every bound below means *all* numbers under it were searched: each length")
    w("  (in the enumerated base) is logged as done only after the whole length finished. Lengths split")
    w("  across two machines count only when all parts used an identical, fingerprinted search plan.")
    w("* **Verification.** Every reported number was re-checked independently (GMP and a separate")
    w("  pure-Python checker: palindromic in both bases, Baillie-PSW for primality). For each sequence,")
    w("  every term already in the OEIS below the bound was re-found (none missing). The engines are")
    w("  cross-checked against brute force and against 14 full OEIS b-files (`make test`).")
    w("* **Method** (for a comment or a program link, if wanted): exhaustive meet-in-the-middle search.")
    w("  Fix the outer digits of a palindrome in one base; that fixes the leading digits in the other")
    w("  base and hence (palindrome) its trailing digits; matching middles are found by table lookup,")
    w("  ~N^(1/4) work for coprime bases. GPU (Apple Metal and NVIDIA CUDA) + CPU implementation, 'dbpal'.")
    w("* **b-files** are in `results/bfiles/` (format `n a(n)`, with `#` header lines that you may delete).")
    w("")
    w("How to submit (one edit per sequence; editors review each):")
    w("")
    w("1. Sign in at oeis.org, open the sequence, choose **Edit**.")
    w("2. For a b-file: upload it first (the b-file link on the edit page), then change the `%H` line.")
    w("3. Paste the `%S` / `%C` / `%E` lines below; put a short note in the edit summary, e.g.")
    w("   \"exhaustive search to X with a GPU meet-in-the-middle program; new term(s) verified, primality")
    w("   proven (Pocklington)\".")
    w("4. Start with section 1 (new terms), then the bound comments, then the b-files. Several")
    w("   sequences of one family can be submitted the same day, but keep each edit self-contained.")
    w("5. If you publish the program (e.g. on GitHub), a `%H` link to it can be added to A393014 and the")
    w("   entries whose terms it found, e.g. `%H A393014 " + who + ", <a href=\"https://...\">dbpal</a>, GPU")
    w("   search for numbers palindromic in two bases.`")
    w("")
    w("Bounds are those of the search log when this file was generated. Re-run")
    w("`python3 tools/submissions.py --refresh` after more searching or after OEIS edits, so that this")
    w("file proposes only what is still missing from the OEIS.")
    w("")

    # ---------------------------------------------------------------- bound comments vs the OEIS
    for it in items:
        inf, st = it["info"], it["state"]
        bexact, bround, lim = nice_bound(it["P"], it["L"])
        nxt = inf["offset"] + len(it["known"])
        btxt = bexact if it["P"] == 10 else f"{bexact} (about {sci(lim)})"
        it["core"] = f"a({nxt}) > {btxt}, if it exists."
        existing = inf["bound_comments"] + st["pending_comments"]
        it["bound_done"] = any(c.startswith(it["core"]) for c in existing)
        it["own_bound"] = next((c for c in existing if who in c and not c.startswith(it["core"])), None)

    # ---------------------------------------------------------------- already submitted
    tier1 = [it for it in items if it["new"] and (it["primes"] or len(it["known"]) <= 20)]
    tier2 = [it for it in items if it["primes"] and not it["new"] and not it["bound_done"]]
    tier3 = [it for it in items if not it["primes"] and it["new"] and it not in tier1]
    done = [it for it in items if it["state"]["mine"]]
    w("## 0. Edits already made on the OEIS")
    w("")
    if not done:
        w(f"None found: no revision by {who} in the OEIS histories of these sequences.")
    else:
        fetched = min(os.path.getmtime(os.path.join(FULL, it["anum"] + ".history.html")) for it in done)
        w(f"From the OEIS revision histories (fetched {time.strftime('%b %d %Y %H:%M', time.localtime(fetched))}).")
        w("What is published or pending review is not proposed again below.")
        w("")
        w("| Sequence | Your edit | OEIS status | Still to do |")
        w("|---|---|---|---|")
        for it in sorted(done, key=lambda x: x["anum"]):
            st, inf = it["state"], it["info"]
            edits = "<br>".join(f"{edit_summary(r)} ({r['date']})" for r in reversed(st["mine"]))
            if st["pending"]:
                status = {"proposed": "proposed, awaiting review", "reviewed": "reviewed, awaiting approval",
                          "editing": "draft (status editing): not yet submitted"}.get(st["status"], st["status"])
            else:
                status = f"approved, published {inf['revdate']} (#{inf['rev']})"
            todo = []
            if it["new"]:
                todo.append(f"{len(it['new'])} new term(s) (section {'1' if it not in tier3 else '3'})")
            if it in tier2:
                todo.append("new bound comment `" + it["core"] + "`" + (" replacing yours" if it["own_bound"] else "") + " (section 2)")
            w(f"| {it['anum']} | {edits} | {status} | {'; '.join(todo) or 'nothing'} |")
        notes = [(it["anum"], u, t, d) for it in done for u, t, d in it["state"]["notes"]]
        if notes:
            w("")
            w("Editor notes:")
            w("")
            for anum, u, t, d in notes:
                w(f"* {anum}, {u} ({d}): {t}")
    w("")

    # ---------------------------------------------------------------- new terms
    w("## 1. New terms")
    w("")
    for it in tier1:
        inf, anum = it["info"], it["anum"]
        P, L = it["P"], it["L"]
        bexact, bround, lim = nice_bound(P, L)
        n_old = len(it["known"])
        w(f"### {anum} — {inf['name']}")
        w("")
        w(f"Currently {n_old} terms (largest {it['known'][-1]}); keywords `{inf['keywords']}`.")
        for c in inf["bound_comments"]:
            w(f"Existing bound comment: \"{c}\"")
        w("")
        idx0 = inf["offset"] + n_old
        newtxt = ",".join(str(n) for n in it["new"])
        fits = len(inf["data"]) + 1 + len(newtxt) <= DATA_LIMIT
        w("Proposed edits:")
        w("")
        w("```")
        if fits:
            w(f"%S {anum} {inf['data']},{newtxt}")
            w(f"   (append to the data; the field stays within ~{DATA_LIMIT} characters)")
        else:
            w(f"   DATA is already {len(inf['data'])} characters: the new term(s) go in the b-file only")
        rng = f"a({idx0})" if len(it["new"]) == 1 else f"a({idx0})-a({idx0 + len(it['new']) - 1})"
        w(f"%E {anum} {rng} from {sig}")
        nxt = idx0 + len(it["new"])
        w(f"%C {anum} a({nxt}) > {bexact}" + (f" (about {sci(lim)})" if P != 10 else "") + f", if it exists. - {sig}")
        if inf["bfile"] or not fits:
            n_new = inf["offset"] + len(it["mine"]) - 1
            cred = credit_suffix(inf["bfile"], n_old) if inf["bfile"] else ""
            cred = "; ".join(p for p in (cred, f"all terms < {bexact}") if p)
            w(f"%H {anum} {who}, <a href=\"/{anum}/b{anum[1:]}.txt\">Table of n, a(n) for n = {inf['offset']}..{n_new}</a>"
              + (f" ({cred})" if cred else ""))
            w(f"   (upload results/bfiles/b{anum[1:]}.txt)")
        w("```")
        w("")
        for n in it["new"]:
            w(f"* {n}")
            w(f"  * base {it['b1']}: `{to_str(n, it['b1'])}`")
            w(f"  * base {it['b2']}: `{to_str(n, it['b2'])}`")
            if it["primes"]:
                c = certify.certify(n)
                w(f"  * prime: " + ("proven (Pocklington certificate, Appendix B)" if c else "probable prime (BPSW)"))
        if inf["bound_comments"]:
            w("")
            w("The existing bound comment(s) are superseded by the new term; editors may delete or keep them.")
        if not it["primes"] and anum == "A060792":
            w("")
            w("Notable: first new term of A060792 since 2014 (a(18) > 3^93 was the last published progress, 2016).")
        if it["primes"] and "more" not in inf["keywords"]:
            w("")
            w("Keyword: consider adding `more` (further terms would be welcome; the sequence is not known to be finite).")
        w("")

    # ---------------------------------------------------------------- bounds
    w("## 2. New search bounds for prime sequences (no new term)")
    w("")
    w("Each line is a proposed `%C` comment; where an older bound comment exists it can be replaced.")
    w("\"(yours)\" marks your own earlier comment: edit it in place rather than adding a second one.")
    w("")
    w("| Sequence | Bases | Terms | Existing bound | Proposed comment |")
    w("|---|---|---|---|---|")
    for it in sorted(tier2, key=lambda x: x["anum"]):
        inf = it["info"]
        old = "; ".join(c.split(" - ")[0] + (" (yours)" if who in c else "")
                        for c in inf["bound_comments"] + it["state"]["pending_comments"]) or "—"
        w(f"| {it['anum']} | {it['b1']}, {it['b2']} | {len(it['known'])} | {old} | `{it['core']} - {sig}` |")
    w("")
    uptodate = [it for it in items if it["primes"] and not it["new"] and it["bound_done"]]
    if uptodate:
        w("Already in the OEIS (published or pending), nothing to add: " +
          ", ".join(f"{it['anum']} (`{it['core']}`)" for it in sorted(uptodate, key=lambda x: x["anum"])) + ".")
        w("")
    w("Notes:")
    w("")
    w("* A393014: the largest known term is a(9) = 2073722573406568398837217962245683 (about 2.07*10^33).")
    w("  The search covered all numbers below 9^54 (every base-9 length up to 53; even base-9 lengths cannot")
    w("  occur, since an even-length base-9 palindrome is a multiple of 10). Worth adding to A393014 also:")
    w(f"  `%C A393014 a(10) > 9^54 (about 3.38*10^51), if it exists. - {sig}`.")
    w("* A046472: the bound above is from this search. In addition, none of the 183 known terms of")
    w("  A007632 (its b-file, complete to 53 decimal digits according to E. Schacham, plus eight further")
    w("  terms with 54-55 digits) is prime, which would give a(7) > 10^53 if that b-file is complete.")
    w("")

    # ---------------------------------------------------------------- b-files
    w("## 3. Extended b-files (numbers palindromic in two bases)")
    w("")
    w("Each b-file is complete up to the stated bound and starts with the existing terms (checked equal).")
    w("Upload the file from `results/bfiles/` and replace the `%H` b-file line with the one given.")
    w("")
    w("| Sequence | Bases | Existing b-file | New b-file | New terms |")
    w("|---|---|---|---|---|")
    for it in sorted(tier3, key=lambda x: x["anum"]):
        inf = it["info"]
        bexact, bround, lim = nice_bound(it["P"], it["L"])
        n_new = inf["offset"] + len(it["mine"]) - 1
        old = re.search(r"for n = (\d+\.\.\d+)", inf["bfile"] or "")
        w(f"| {it['anum']} | {it['b1']}, {it['b2']} | {old.group(1) if old else '—'} | {inf['offset']}..{n_new} (all terms < {bexact}) | {len(it['new'])} |")
    w("")
    w("Proposed `%H` lines:")
    w("")
    w("```")
    for it in sorted(tier3, key=lambda x: x["anum"]):
        inf, anum = it["info"], it["anum"]
        bexact, bround, lim = nice_bound(it["P"], it["L"])
        n_old = len(it["known"])
        n_new = inf["offset"] + len(it["mine"]) - 1
        cred = credit_suffix(inf["bfile"], n_old)
        parts = [p for p in (cred, f"all terms < {bexact}") if p]
        w(f"%H {anum} {who}, <a href=\"/{anum}/b{anum[1:]}.txt\">Table of n, a(n) for n = {inf['offset']}..{n_new}</a> ({'; '.join(parts)})")
    w("```")
    w("")
    w("Where the data section is short of the limit, a few more terms could also be appended to `%S`;")
    w("usually the b-file is enough.")
    w("")

    # ---------------------------------------------------------------- unchanged
    same = [it for it in items if not it["primes"] and not it["new"]]
    if same:
        w("Searched with no new terms (b-file not extended): " +
          ", ".join(f"{it['anum']} ({it['b1']},{it['b2']}; all terms < {nice_bound(it['P'], it['L'])[0]} confirmed)"
                    for it in same) + ".")
        w("")

    # ---------------------------------------------------------------- new sequences
    w("## 4. Possible new sequences (optional)")
    w("")
    w("Primes that are palindromic in two bases where the OEIS has no such sequence. An OEIS search for")
    w("the first terms of each list returned no match (Oct 2026). The (10, b) ones would continue the")
    w("family A046472-A046484 (which stops at b = 16); the others are the primes in an existing")
    w("\"palindromic in bases b1 and b2\" sequence (cross-reference it). Editors may consider some of these")
    w("too sparse or too arbitrary; submit only the ones you find worthwhile.")
    w("")
    for b1, b2, P, L, pr in newseq:
        allseq = fam.get(((b1, b2), False))
        bexact, bround, lim = nice_bound(P, L)
        w(f"* **Primes that are palindromic in bases {b1} and {b2}** "
          + (f"(primes in {allseq['anum']}); " if allseq else "; ")
          + f"complete below {bexact}" + (f" (about {sci(lim)})" if P != 10 else "") + f": {', '.join(map(str, pr))}")
    w("")

    # ---------------------------------------------------------------- appendix A
    w("## Appendix A. Certified bounds per base pair")
    w("")
    w("| Bases | Searched below | GPU hours logged |")
    w("|---|---|---|")
    seen = set()
    for it in sorted(items, key=lambda x: (x["b1"], x["b2"])):
        k = (it["b1"], it["b2"])
        if k in seen:
            continue
        seen.add(k)
        bexact, bround, lim = nice_bound(it["P"], it["L"])
        w(f"| {it['b1']}, {it['b2']} | {bexact} (about {sci(lim)}) | {it['hours']:.2f} |")
    w("")

    # ---------------------------------------------------------------- appendix B
    w("## Appendix B. Primality certificates of the new prime terms")
    w("")
    w("Pocklington's theorem: if N-1 = F*R with F > sqrt(N), and for every prime q dividing F some base a")
    w("satisfies a^(N-1) = 1 (mod N) and gcd(a^((N-1)/q) - 1, N) = 1, then N is prime. Here N-1 is factored")
    w("completely; prime factors above 10^12 are certified the same way (recursively), smaller ones by")
    w("trial division. Regenerate with `python3 tools/certify.py N`.")
    w("")
    for it in sorted(items, key=lambda x: (x not in tier1, x["anum"])):
        if not it["primes"]:
            continue
        for n in it["new"] + it["pending_new"]:
            c = certify.certify(n)
            w(f"### {it['anum']}: {n}" + (" (submitted, pending review)" if n in it["pending_new"] else ""))
            w("")
            w("```")
            lines = []

            def emit(cc, ind=0):
                pad = "  " * ind
                if cc["method"] == "trial division":
                    lines.append(f"{pad}{cc['N']} prime (trial division)")
                    return
                lines.append(f"{pad}N = {cc['N']}: N-1 fully factored, witnesses:")
                for f in cc["factors"]:
                    lines.append(f"{pad}  q = {f['q']}" + (f"^{f['e']}" if f["e"] > 1 else "") + f"  a = {f['a']}")
                    if f["cert"]["method"] != "trial division":
                        emit(f["cert"], ind + 2)

            if c:
                emit(c)
            else:
                lines.append("no certificate produced (probable prime by BPSW)")
            for l in lines:
                w(l)
            w("```")
            w("")
    path = os.path.join(ROOT, "submissions.md")
    open(path, "w").write("\n".join(out) + "\n")
    print(path)


if __name__ == "__main__":
    main()
