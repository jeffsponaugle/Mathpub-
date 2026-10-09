#!/usr/bin/env python3
"""Build a self-contained HTML results page from results/ (writes results/site/index.html)."""
import html
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import oeis  # noqa: E402
from campaign import complete_lengths, log_entries, roles  # noqa: E402
from stats import expected, ndig  # noqa: E402
from verify import to_str  # noqa: E402

RES = os.path.join(ROOT, "results")


def bound(b1, b2):
    rs = [r for r in log_entries() if sorted(r["bases"]) == sorted([b1, b2]) and not r["primesOnly"]]
    if not rs:
        return None
    P = rs[0]["P"]
    Ls = complete_lengths(rs)
    L = 0
    while L + 1 in Ls:
        L += 1
    secs = sum(r.get("seconds", 0) + r.get("tableSeconds", 0) for r in rs if r["L"] <= L)
    hosts = sorted({r.get("host", "M4 Max") for r in rs})
    return P, L, secs, hosts


def collect():
    fam = oeis.family()
    pairs = []
    for fn in sorted(os.listdir(RES)):
        if not (fn.startswith("b") and fn.endswith(".tsv")):
            continue
        b1, b2 = map(int, fn[1:-4].split("_"))
        bd = bound(b1, b2)
        if not bd or bd[1] == 0:
            continue
        P, L, secs, hosts = bd
        Q = b1 if P == b2 else b2
        lim = P ** L
        terms = {}
        for line in open(os.path.join(RES, fn)):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].isdigit() and int(p[0]) < lim:
                terms[int(p[0])] = p[1] == "1"
        seqs = []
        for primes in (False, True):
            e = fam.get(((min(b1, b2), max(b1, b2)), primes))
            mine = sorted(n for n, isp in terms.items() if isp or not primes)
            if not e:
                continue
            known = sorted(set(oeis.all_terms(e)))
            kset = set(known)
            new = [n for n in mine if n not in kset]
            missing = [n for n in known if n < lim and n not in set(mine)]
            seqs.append({"anum": e["anum"], "name": e["name"], "primes": primes, "known": len(known),
                         "maxKnown": str(max(known)) if known else None, "found": len(mine),
                         "new": [str(n) for n in new], "missing": [str(n) for n in missing]})
        # observed vs expected (cumulative, positive terms)
        rows = []
        co = ce = 0.0
        obs = {}
        for n in terms:
            if n > 0:
                obs[ndig(n, P)] = obs.get(ndig(n, P), 0) + 1
        for l in range(1, L + 1):
            e, _ = expected(P, Q, l)
            co += obs.get(l, 0)
            ce += e
            rows.append([l, co, round(ce, 3)])
        newterms = []
        for s in seqs:
            for n in s["new"]:
                newterms.append(n)
        detail = {}
        for n in sorted(set(int(x) for x in newterms)):
            detail[str(n)] = {"prime": terms.get(n, False), "r1": to_str(n, min(b1, b2)), "r2": to_str(n, max(b1, b2))}
        pairs.append({"b1": min(b1, b2), "b2": max(b1, b2), "P": P, "Q": Q, "L": L, "bound": str(lim),
                      "boundLog10": math.log10(lim), "hours": secs / 3600, "hosts": hosts, "count": len(terms),
                      "primes": sum(terms.values()), "seqs": seqs, "rows": rows, "detail": detail})
    return pairs


CSS = r"""
:root{
  --paper:#f2f4f1; --surface:#fbfcfa; --ink:#16201c; --muted:#56655f; --rule:#d2dad5;
  --accent:#0b6b66; --accent2:#9a5a14; --hl:#e4eee9; --good:#1f6f3c; --warn:#9a3b12;
  --display:"STIX Two Text", "STIX Two Math", Georgia, "Times New Roman", serif;
  --body:"IBM Plex Sans", "Helvetica Neue", Arial, sans-serif;
  --mono:"IBM Plex Mono", ui-monospace, "SF Mono", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){ color-scheme:dark;
    --paper:#0e1412; --surface:#141c19; --ink:#e1e9e5; --muted:#92a39c; --rule:#2a3531;
    --accent:#4dc0b6; --accent2:#e0a456; --hl:#182621; --good:#6fcf8e; --warn:#f08a5d; }
}
:root[data-theme="dark"]{ color-scheme:dark;
  --paper:#0e1412; --surface:#141c19; --ink:#e1e9e5; --muted:#92a39c; --rule:#2a3531;
  --accent:#4dc0b6; --accent2:#e0a456; --hl:#182621; --good:#6fcf8e; --warn:#f08a5d; }
body{ background:var(--paper); color:var(--ink); font:15px/1.55 var(--body); }
.wrap{ max-width:1080px; margin:0 auto; padding-inline:clamp(16px,4vw,40px); padding-block:28px 64px; }
h1,h2,h3{ font-family:var(--display); font-weight:600; text-wrap:balance; margin:0; }
h1{ font-size:clamp(30px,5vw,46px); line-height:1.08; letter-spacing:-.01em; }
h2{ font-size:26px; margin-top:56px; padding-bottom:8px; border-bottom:1px solid var(--rule); }
h3{ font-size:19px; }
p{ max-width:68ch; margin:0; }
.eyebrow{ font:600 11px/1 var(--body); letter-spacing:.14em; text-transform:uppercase; color:var(--muted); }
.lede{ display:grid; gap:14px; margin-top:10px; }
.lede p{ color:var(--muted); font-size:16px; }
.mono{ font-family:var(--mono); font-variant-numeric:tabular-nums; }
.finding{ margin-top:28px; background:var(--surface); border:1px solid var(--rule); border-radius:6px; padding:22px clamp(16px,3vw,28px); display:grid; gap:14px; }
.finding .label{ display:flex; flex-wrap:wrap; gap:8px 14px; align-items:baseline; }
.finding .label b{ font-family:var(--display); font-size:21px; font-weight:600; }
.digits{ font-family:var(--mono); font-size:clamp(14px,2.2vw,19px); line-height:1.5; word-break:break-all; }
.digits .base{ display:inline-block; min-width:5.5em; color:var(--muted); font-size:12px; letter-spacing:.06em; text-transform:uppercase; font-family:var(--body); }
.m1{ color:var(--accent); } .m2{ color:var(--accent2); } .mc{ font-weight:700; text-decoration:underline; text-underline-offset:3px; }
.facts{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:12px; margin-top:16px; }
.fact{ border-top:2px solid var(--accent); padding-top:8px; display:grid; gap:2px; }
.fact .v{ font-family:var(--display); font-size:22px; }
.fact .k{ color:var(--muted); font-size:13px; }
.tbl{ overflow-x:auto; margin-top:16px; border:1px solid var(--rule); border-radius:6px; background:var(--surface); }
table{ border-collapse:collapse; width:100%; min-width:760px; font-size:14px; }
th,td{ text-align:left; padding:9px 12px; border-bottom:1px solid var(--rule); vertical-align:top; }
th{ font:600 11px/1.2 var(--body); letter-spacing:.08em; text-transform:uppercase; color:var(--muted); background:var(--hl); }
td.num{ font-family:var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap; }
tr:last-child td{ border-bottom:0; }
.pill{ display:inline-block; padding:1px 8px; border-radius:999px; font-size:12px; font-weight:600; white-space:nowrap; }
.pill.new{ background:color-mix(in srgb,var(--good) 16%,transparent); color:var(--good); }
.pill.bound{ background:color-mix(in srgb,var(--accent) 14%,transparent); color:var(--accent); }
.pill.miss{ background:color-mix(in srgb,var(--warn) 16%,transparent); color:var(--warn); }
a{ color:var(--accent); text-underline-offset:2px; }
a:focus-visible, summary:focus-visible, button:focus-visible{ outline:2px solid var(--accent); outline-offset:2px; }
.pairs{ display:grid; gap:14px; margin-top:16px; }
details.pair{ background:var(--surface); border:1px solid var(--rule); border-radius:6px; }
details.pair > summary{ cursor:pointer; padding:14px 18px; display:flex; flex-wrap:wrap; gap:6px 18px; align-items:baseline; list-style:none; }
details.pair > summary::-webkit-details-marker{ display:none; }
details.pair > summary::before{ content:"\25B8"; color:var(--muted); width:1em; transition:transform .15s; }
details.pair[open] > summary::before{ transform:rotate(90deg); }
@media (prefers-reduced-motion: reduce){ details.pair > summary::before{ transition:none; } }
.pair .sumk{ color:var(--muted); font-size:13px; }
.pairbody{ padding:0 18px 18px; display:grid; gap:18px; }
.cols{ display:grid; grid-template-columns:minmax(0,1fr) minmax(0,1fr); gap:18px; }
@media (max-width:760px){ .cols{ grid-template-columns:1fr; } }
.terms{ display:grid; gap:10px; max-height:420px; overflow:auto; padding-right:4px; }
.term{ border-left:2px solid var(--rule); padding-left:10px; display:grid; gap:2px; }
.term.prime{ border-left-color:var(--good); }
.term .n{ font-family:var(--mono); font-size:13.5px; word-break:break-all; }
.term .r{ font-family:var(--mono); font-size:11.5px; color:var(--muted); word-break:break-all; }
.chart svg{ width:100%; height:auto; display:block; }
.chart .cap{ font-size:12.5px; color:var(--muted); margin-top:6px; }
.method{ display:grid; gap:12px; margin-top:16px; }
.method p{ color:var(--ink); }
pre{ background:var(--surface); border:1px solid var(--rule); border-radius:6px; padding:12px 14px; overflow-x:auto; font:13px/1.5 var(--mono); }
.foot{ margin-top:48px; color:var(--muted); font-size:13px; }
"""

JS = r"""
const D = JSON.parse(document.getElementById('data').textContent);
const $ = (t, a = {}, ...c) => { const e = document.createElement(t); for (const [k, v] of Object.entries(a)) { if (k === 'class') e.className = v; else if (k === 'html') e.innerHTML = v; else e.setAttribute(k, v); } for (const x of c) e.append(x); return e; };
function sci(s) { const n = BigInt(s); const len = s.length; if (len <= 6) return s; const m = (Number(s.slice(0, 4)) / 1000).toFixed(2); return m + '·10' + sup(len - 1); }
function sup(k) { const m = '⁰¹²³⁴⁵⁶⁷⁸⁹'; return String(k).split('').map(d => m[+d]).join(''); }
function mirror(str) {  // colour the two halves of a palindromic digit string
  const n = str.length, h = Math.floor(n / 2);
  const a = str.slice(0, h), c = n % 2 ? str[h] : '', b = str.slice(n - h);
  return `<span class="m1">${a}</span>${c ? `<span class="mc">${c}</span>` : ''}<span class="m2">${b}</span>`;
}
function chart(rows, P) {
  const W = 520, H = 230, L = 46, R = 14, T = 14, B = 34;
  const xmax = rows.length, ymax = Math.max(1, ...rows.map(r => Math.max(r[1], r[2])));
  const step = ymax > 60 ? 20 : ymax > 24 ? 10 : 5;
  const ytop = Math.ceil(ymax / step) * step;
  const x = l => L + (l - 1) / Math.max(1, xmax - 1) * (W - L - R), y = v => H - B - v / ytop * (H - T - B);
  let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="cumulative count of double palindromes by length, observed and expected">`;
  for (let v = 0; v <= ytop; v += step) s += `<line x1="${L}" x2="${W - R}" y1="${y(v)}" y2="${y(v)}" stroke="var(--rule)" stroke-width="1"/><text x="${L - 6}" y="${y(v) + 4}" text-anchor="end" font-size="11" fill="var(--muted)" font-family="var(--mono)">${v}</text>`;
  const xt = xmax > 60 ? 20 : xmax > 30 ? 10 : 5;
  for (let l = xt; l <= xmax; l += xt) s += `<text x="${x(l)}" y="${H - B + 16}" text-anchor="middle" font-size="11" fill="var(--muted)" font-family="var(--mono)">${l}</text>`;
  s += `<text x="${(L + W - R) / 2}" y="${H - 4}" text-anchor="middle" font-size="11" fill="var(--muted)">length in base ${P}</text>`;
  const pe = rows.map((r, i) => `${i ? 'L' : 'M'}${x(r[0]).toFixed(1)},${y(r[2]).toFixed(1)}`).join('');
  s += `<path d="${pe}" fill="none" stroke="var(--accent2)" stroke-width="1.6" stroke-dasharray="4 3"/>`;
  let po = ''; rows.forEach((r, i) => { const xx = x(r[0]).toFixed(1), yy = y(r[1]).toFixed(1); po += i ? `H${xx}V${yy}` : `M${xx},${yy}`; });
  s += `<path d="${po}" fill="none" stroke="var(--accent)" stroke-width="2"/>`;
  const last = rows[rows.length - 1];
  s += `<circle cx="${x(last[0])}" cy="${y(last[1])}" r="3.5" fill="var(--accent)"/>`;
  s += `<text x="${x(last[0]) - 6}" y="${y(last[1]) - 8}" text-anchor="end" font-size="11.5" fill="var(--accent)" font-family="var(--mono)">${last[1]} found</text>`;
  s += `<text x="${x(last[0]) - 6}" y="${y(last[2]) + 16}" text-anchor="end" font-size="11.5" fill="var(--accent2)" font-family="var(--mono)">${last[2].toFixed(1)} expected</text>`;
  return s + '</svg>';
}
// summary table
const tb = document.getElementById('summary');
for (const p of D.pairs) for (const q of p.seqs) {
  const status = q.missing.length ? `<span class="pill miss">${q.missing.length} missing</span>` : q.new.length ? `<span class="pill new">${q.new.length} new</span>` : `<span class="pill bound">bound raised</span>`;
  tb.append($('tr', { html: `<td><a href="https://oeis.org/${q.anum}">${q.anum}</a></td><td>${q.primes ? 'primes' : 'all'} in bases ${p.b1} &amp; ${p.b2}</td><td class="num">${q.known}</td><td class="num">${q.maxKnown ? sci(q.maxKnown) : '–'}</td><td class="num">${p.P}${sup(p.L)} ≈ 10${sup(Math.floor(p.boundLog10))}</td><td class="num">${q.found}</td><td>${status}</td>` }));
}
// per pair details
const host = document.getElementById('pairs');
for (const p of D.pairs) {
  const d = $('details', { class: 'pair', id: `b${p.b1}-${p.b2}` });
  const newCount = new Set(p.seqs.flatMap(q => q.new)).size;
  d.append($('summary', { html: `<h3>Bases ${p.b1} and ${p.b2}</h3><span class="sumk">complete below ${p.P}${sup(p.L)} ≈ 10${sup(Math.floor(p.boundLog10))} · ${p.count} double palindromes, ${p.primes} prime · ${newCount} not in OEIS · ${p.hours < 0.1 ? Math.round(p.hours * 60) + ' min' : p.hours.toFixed(1) + ' h'} of search</span>` }));
  const body = $('div', { class: 'pairbody' });
  const cols = $('div', { class: 'cols' });
  const terms = $('div', { class: 'terms' });
  const nk = Object.keys(p.detail);
  if (!nk.length) terms.append($('p', { class: 'sumk' }, 'No terms beyond the OEIS data below this bound.'));
  for (const n of nk) {
    const t = p.detail[n];
    terms.append($('div', { class: 'term' + (t.prime ? ' prime' : ''), html: `<span class="n">${n}${t.prime ? ' <span class="pill new">prime</span>' : ''}</span><span class="r">base ${p.b1}: ${mirror(t.r1)}</span><span class="r">base ${p.b2}: ${mirror(t.r2)}</span>` }));
  }
  const left = $('div', {}, $('div', { class: 'eyebrow' }, 'Terms not in the OEIS data'), terms);
  left.style.display = 'grid'; left.style.gap = '10px';
  const ch = $('div', { class: 'chart', html: `<div class="eyebrow" style="margin-bottom:8px">Found vs. expected (cumulative)</div>${chart(p.rows, p.P)}<p class="cap">Solid: double palindromes found up to each base-${p.P} length. Dashed: naive expectation, treating the two palindrome conditions as independent.</p>` });
  cols.append(left, ch); body.append(cols); d.append(body); host.append(d);
}
if (location.hash) { const el = document.getElementById(location.hash.slice(1)); if (el && el.tagName === 'DETAILS') el.open = true; }
"""


def page(pairs, highlights):
    data = json.dumps({"pairs": pairs})
    hl = ""
    if highlights:
        cards = []
        for n, b1, b2, anum, idx, note, isprime in highlights:
            r1, r2 = to_str(n, b1), to_str(n, b2)
            pill = '<span class="pill new">prime</span>' if isprime else '<span class="pill bound">first new term since 2014</span>' if anum == "A060792" else ''
            cards.append(f"""
<section class="finding" aria-label="new term {anum} a({idx})">
  <div class="label"><span class="eyebrow">New term</span><b><a href="https://oeis.org/{anum}">{anum}</a>, a({idx})</b>{pill}</div>
  <div class="digits"><span class="base">base {b2}</span><span data-mirror>{r2}</span></div>
  <div class="digits"><span class="base">base {b1}</span><span data-mirror>{r1}</span></div>
  <p class="sumk" style="color:var(--muted)">{html.escape(note)}</p>
</section>""")
        hl = "\n".join(cards)
    return f"""<title>Two-Base Palindrome Census</title>
<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;600&family=STIX+Two+Text:wght@500;600&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="lede">
  <span class="eyebrow">OEIS A393014 and relatives &middot; exhaustive search with dbpal</span>
  <h1>Numbers that read the same in two bases</h1>
  <p>Every number below each bound was searched for palindromes in both bases at once, on Apple-silicon GPUs. Bounds are certified length by length in the search log, so a bound means no further term exists below it.</p>
</header>
{hl}
<div class="facts" id="facts"></div>
<h2>Sequences</h2>
<div class="tbl"><table><thead><tr><th>OEIS</th><th>sequence</th><th>terms in OEIS</th><th>largest in OEIS</th><th>searched below</th><th>found below bound</th><th>status</th></tr></thead><tbody id="summary"></tbody></table></div>
<h2>By base pair</h2>
<div class="pairs" id="pairs"></div>
<h2>How the search works</h2>
<div class="method">
<p>A palindrome in base P of length L is cut into k outer digit pairs and a middle palindrome M: N = C + M&middot;P<sup>k</sup>. Fixing the outer digits pins N to an interval, which fixes its top s digits in the other base Q. If N is also a base-Q palindrome, its bottom s digits are those top digits reversed, so M&middot;P<sup>k</sup> &equiv; rev(T) &minus; C (mod Q<sup>s</sup>). All middles are bucketed by this residue in advance, and each node then costs one lookup instead of P<sup>m/2</sup> palindrome tests. Balanced, that is about N<sup>1/4</sup> work, against N<sup>1/2</sup> for testing the palindromes of one base.</p>
<p>For very large N the middles are generated class by class (residue mod Q<sup>r</sup>) and matched against the outer digits of the same class, which keeps every lookup table small enough for the GPU cache. When the bases share a factor (2 &amp; 10, 8 &amp; 10, 10 &amp; 15), the two residues must also agree modulo the common factor, which prunes the search tree at every level.</p>
<p>Every reported number is re-checked with GMP and by an independent Python verifier (Baillie&ndash;PSW for primality). The engines are cross-checked against brute force on thousands of parameter choices and against the existing OEIS b-files.</p>
</div>
<h2>Reproduce</h2>
<pre>make &amp;&amp; make test
./dbpal search -b 2,9 --below 9^47 --out results/b2_9.tsv --log results/search_log.jsonl
python3 tools/campaign.py run --bases 9,10 --below 10^44
python3 tools/report.py          # RESULTS.md
python3 tools/bfile.py results/b9_10.tsv</pre>
<p class="foot">Generated by tools/site.py from the dbpal results directory.</p>
</div>
<script type="application/json" id="data">{data}</script>
<script>{JS}
for (const el of document.querySelectorAll('[data-mirror]')) el.innerHTML = mirror(el.textContent);
const F = document.getElementById('facts');
const facts = {json.dumps(FACTS)};
for (const f of facts) F.append($('div', {{class:'fact', html:`<span class="v">${{f[0]}}</span><span class="k">${{f[1]}}</span>`}}));
</script>
"""


FACTS = []


def main():
    pairs = collect()
    # headline: every prime term that is not in the OEIS data
    highlights = []
    for p in pairs:
        for q in p["seqs"]:
            # prime sequences, and sparse sequences (few terms known) where a new term is notable
            sparse = (not q["primes"]) and q["known"] <= 20 and 0 < len(q["new"]) <= 3
            if not (q["primes"] or sparse):
                continue
            prev = q["maxKnown"]
            for i, ns in enumerate(q["new"]):
                n = int(ns)
                kind = "Prime, and a palindrome" if q["primes"] else "A palindrome"
                highlights.append((n, p["b1"], p["b2"], q["anum"], q["known"] + 1 + i,
                                   f"{kind} in base {p['b1']} and in base {p['b2']}. It follows "
                                   f"{prev} with no other term in between; the search is complete below "
                                   f"{p['P']}^{p['L']}.", q["primes"]))
                prev = ns
    totnew = len({(p["b1"], p["b2"], n) for p in pairs for s in p["seqs"] for n in s["new"]})
    a393 = next((p for p in pairs if (p["b1"], p["b2"]) == (2, 9)), None)
    if a393:
        FACTS.append([f"{a393['P']}<sup>{a393['L']}</sup> ≈ 10<sup>{int(a393['boundLog10'])}</sup>",
                      "A393014 searched to here: no prime beyond a(9) = 2.07·10<sup>33</sup>"])
    FACTS.append([str(len(highlights)), "new terms of sparse sequences, each the next term after the OEIS data"])
    FACTS.append([str(totnew), "terms not in the OEIS data, across all sequences searched"])
    FACTS.append([str(len(pairs)), "base pairs searched exhaustively"])
    out = os.path.join(RES, "site")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "index.html"), "w") as f:
        f.write(page(pairs, highlights))
    print(os.path.join(out, "index.html"))


if __name__ == "__main__":
    main()
