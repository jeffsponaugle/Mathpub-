# Draft OEIS submission: A001211(27) = n(27,6)

STATUS: COMPLETE 2026-10-03 08:13 PDT. **n(27,6) = 186942.** Final audit: logs/n27k6_final_audit.txt (30/30 chunks of 100000
work items, no gaps/overlaps; all 4 mixed-target chunk boundaries covered by re-run windows; 1147 reported bases).

## A001211 (postage stamp problem with 6 denominations and n stamps): new term a(27)

Current OEIS data ends at a(25) = 129783. The July 2013 Challis-Robinson addendum value a(26) = 156744 is not on OEIS yet
(see bfiles/bA001211.txt, which also corrects a(20) = 45754, misprinted as 45745 on OEIS).

Data (to append): ..., 106888, 129783, 156744, **186942**

Comment / example:
  a(27) = 186942: the extremal 27-basis with 6 denominations is {1, 19, 194, 1095, 7370, 27669} (every integer 1..186942 is a
  sum of at most 27 of these denominations, 186943 is not), found by an exhaustive search over all admissible bases
  (Challis's H-program method: admissibility a_{j+1} <= n_27(A_j)+1, element-wise bounds from the proven values
  n_27(1..5) = 27, 223, 1524, 8913, 42744, the Challis-Robinson "difficult target" test, deep-hole targets and an exact
  full check; 2964437 prefix work items at depth 4, 1.11e+14 leaf candidates, 7.35e+10 full checks, 3326 core-hours
  on three Xeon servers, Oct 2-3 2026).

Extension line: a(26) from the Challis-Robinson addendum (2013), a(27) from Jeff Sponaugle, Oct 03 2026.

Also affects A084192 / A084193 / A196416 (arrays): the (h=27, k=6) entry (beyond the current b-files' antidiagonals).

## Verification notes (for the editors / our records)

* Exactness: every rejected candidate failed an exact representability test of a number <= target, or the exact full
  check; the admissibility bounds used were the proven extremal values n_27(k), k = 1..5 (A014616(27), A001208(27),
  A001209(27), A001210(27)).
* Targets: chunks 0-5, 15-19 ran at target 176381 (heuristic bound 176380 + 1, mathg), chunks 6-13, 20, 21 at 177915
  (mathd, mathb), chunks 14, 22-29 at 186943 (after the maximum appeared in chunk 21). All targets are <= 186942 except
  186943 = maximum + 1, so every basis with range >= its chunk's target was reported, and the maximum over all reported
  bases is n(27,6). Coverage of chunks 0..29 exactly once is checked by `tools/collect_h27k6.py`.
* Enumeration subtlety: the depth-4 work-item enumeration has 2964437 items at target 176381, 2964430 at 177915 and
  2964418 at 186943 (prefixes pruned at a higher target are proven unable to reach it). Chunk index spaces at different
  targets differ by a shift <= 19, so every chunk boundary B was re-run as a window [B-24, B+24) at target 176381
  (psph/bnd2_h27k6.sh; an earlier pass with [B-8, B+8) also completed); the audit script checks that each boundary whose
  lower chunk used a lower target than its upper chunk is covered.
* Independent check of the winning basis: `tools/hrange 27 1 19 194 1095 7370 27669` = 186942 (16-bit min-stamp DP).
* Uniqueness: {1, 19, 194, 1095, 7370, 27669} is the only basis with 27-range 186942 among the chunks searched with
  targets <= 186942 (chunks 0-13, 15-21); chunks 14 and 22-29 were searched with target 186943 and would not have
  reported a second basis with range exactly 186942, so uniqueness is not fully established (re-run those 10 chunks at
  target 186942 to settle it, ~15 h on a 96-thread box). The VALUE 186942 is fully established.
* Machines: mathd (2x Xeon 8168, 96 thr; chunks 6-14, 27, 28), mathb (2x Xeon 8268, 96 thr; chunks 20-25, 29), mathg
  (2x Xeon Gold 5315Y, 32 thr; chunks 0-5, 15-19, 26, boundary windows); an M1 Pro briefly ran chunk 20 (abandoned,
  redone on mathb).
* Lower-bound history during the run: 176380 (heuristic) -> 176863 -> 177568 -> 177914 -> 180223 -> 182481 -> 182504
  -> 182526 -> 183362 -> 183614 -> 184129 -> 186942 (mathb chunk 21, Oct 3 00:45 PDT).
