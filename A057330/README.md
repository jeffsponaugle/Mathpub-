# cchain — smallest Cunningham chain search (OEIS A057330)

Exhaustive search for the smallest prime starting a Cunningham chain of
length k.

- **kind 2** (`p -> 2p-1`, members `2^i(p-1)+1`): [A057330](https://oeis.org/A057330).
  16 terms known; **a(17) is open**, with a(16) = 3203000719597029781 (≈3.2e18)
  and the a(16) chain having length exactly 16, so a(17) > a(16).
- **kind 1** (`p -> 2p+1`, members `2^i(p+1)-1`): "length ≥ n" analogue derived
  from [A005602](https://oeis.org/A005602).

Both sequences use "chain of length **at least** k" semantics (no maximality
condition), which is exactly what the search tests.

## Method

1. **Wheel**: for each wheel prime q in {2,3,5,7,11,13,17,19,23}, p must avoid
   the residues r = ±(2^-i − 1) mod q (i < min(k, ord_2(q))) that would make a
   chain member divisible by q. CRT gives the admissible classes mod
   M = 223092870 (for k = 17, kind 2: 864 classes — a 258000× reduction).
2. **Segmented sieve**: candidates p = base + r + M·t; for each sieve prime
   q ≤ B, forbidden t values are marked in a bitmap per (class, window).
3. **PRP tests**: survivors get GMP Miller–Rabin chain-prefix tests; full-length
   hits are re-verified with 40 rounds and appended to the results file.

Perhaps surprisingly, small `-B` wins (default 1024): marking with mid-size
primes is cache-hostile, while GMP PRP tests on the extra survivors are cheap.
On an M-series Mac (14 threads) throughput is ≈3.5e14 p/sec.

**Exhaustiveness caveat**: a chain member equal to a sieve prime q would be
wrongly sieved out; that requires p ≤ B, so results are exhaustive for
p > B. All real searches here start far above B. (The selftest ranges start
at 3, relying on the known values it is checking being > B.)

## Build & test

    make            # needs GMP (brew install gmp) and pthreads
    make test       # --verify-known + --selftest

`--verify-known` checks all 32 known terms (both kinds) by direct primality.
`--selftest` re-discovers five known terms from scratch through the full
sieve+PRP pipeline (kind 2: k = 8, 10, 12; kind 1: k = 8, 10).

## The a(17) campaign

    ./cchain -k 17 -K 2 -s 3203000719597029781 -e 1e21 \
             -c cchain.ckpt -o found.txt >> campaign.log 2>&1 &

- Status lines print every `-u` seconds (Ctrl-T / SIGINFO for an immediate one).
- Checkpoint written every `-i` seconds (default 60) and on Ctrl-C/SIGTERM;
  resume with the same arguments plus `-r`. The frontier only advances over
  fully-completed superblocks, so a crash costs at most a few minutes of work.
- When a chain is found the search continues until the frontier passes it,
  proving minimality; it then prints `RESULT` and exits (unless `-A`).
- Growth of the sequence suggests a(17) ≈ 1e20 (could be 1e19–1e21).
  At ≈3.5e14 p/s that is about 3 days to 1e20, about a month to 1e21.

Before announcing a result: certify all 17 members with an independent tool
(e.g. PARI/GP `isprime`, which is a proving algorithm at this size), and
re-verify the empty range on a second machine if practical.

## Files

- `cchain.c` — the whole tool
- `found.txt` — every full-length chain found (append-only)
- `cchain.ckpt` — checkpoint (config echo + frontier + stats + finds)
- `campaign.log` — status output of the running campaign
