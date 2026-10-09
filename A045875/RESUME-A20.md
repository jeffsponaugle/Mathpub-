# a(20) search — PAUSED 2026-09-27

Paused so the Sparks can run another sequence. All state checkpointed.

## Proven bound at pause

**a(20) > 4,033,865,312** — contiguous coverage:
[3,238,684,957 → 3,725,000,000] (legs A1+B1+A2, logged empty) +
[3,725,000,000 → 3,904,306,680] (free: cleared of 19-runs during a(19)) +
[3,904,306,680 → 4,033,865,312] (leg B2 checkpoint, atom2).

atom1's leg A3 additionally cleared [4,046,000,000 → 4,065,471,256]
(a20-a3.ckpt) — it splices into the contiguous bound once B2 reaches 4,046M.

## To resume (exact commands)

atom2 (finish leg B2, then chained B3 starts automatically):
```bash
ssh jbs@10.1.30.37 'cd A045875 && nohup ./run_a20_atom2_resume.sh > /dev/null 2>&1 < /dev/null &'
```
where run_a20_atom2_resume.sh is:
```bash
#!/bin/bash
cd /home/jbs/A045875
./a045875gpu -n 20 -e 4046000000 -c a20-b2.ckpt >> a20-b2.log 2>&1; [ $? -ne 1 ] && exit 0
exec ./a045875gpu -n 20 -s 4190000000 -c a20-b3.ckpt >> a20-b3.log 2>&1
```

atom1 (finish leg A3; next leg beyond 4,190M to be planned at resume time):
```bash
ssh jbs@10.1.30.36 'cd A045875 && nohup ./a045875gpu -n 20 -e 4190000000 -c a20-a3.ckpt >> a20-a3.log 2>&1 < /dev/null &'
```

Or just ask Claude to "resume the a(20) search" — this file plus the
checkpoints carry everything needed (it will also rebalance legs for
whatever frontier exists then).

## Checkpoint inventory (do not delete)

- atom1: a20-a3.ckpt (m=4,065,471,256), logs a20-a1/a2/a3.log
- atom2: a20-b2.ckpt (m=4,033,865,312), logs a20-b1/b2.log
- a(18)/a(19) provenance checkpoints and logs per RESULTS.md
