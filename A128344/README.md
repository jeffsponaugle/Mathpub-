# a128344

Multithreaded search tool for [OEIS A128344](https://oeis.org/A128344):
numbers k such that (7^k − 5^k)/2 is prime.

```
3, 5, 7, 113, 397, 577, 7573, 14561, 58543, 100019, 123407, 136559,
208283, 210761, 457871, 608347, 636043
```

## Build

Requires [primesieve](https://github.com/kimwalisch/primesieve) and
[GMP](https://gmplib.org/) (`brew install primesieve gmp`), then:

```
make
```

On Intel Macs or Linux set `PREFIX` to your install prefix, e.g.
`make PREFIX=/usr/local`.

## Usage

```
a128344 [options] nterms          # print the first nterms terms
a128344 -N [options]              # hunt the next term, print "HIT: k"
a128344 -F file [options]         # trial-factor only, export survivors
```

- `nterms` — how many terms to print (one per line, on stdout).
- `-N` — search for the next term only: stop at the first hit and print
  it as `HIT: k` plus search statistics. Unless `-s` or a checkpoint says
  otherwise, the search starts at 636044, just past a(17) = 636043.
- `-F file` — skip the bignum PRP tests: append trial-factoring
  survivors to `file` in PFGW ABC format (header `ABC (7^$a-5^$a)/2`).
  Run the PRP tests elsewhere with `pfgw -f0 file` — PFGW/gwnum is
  typically 3–6× faster than GMP at these sizes. The file is written in
  k order and appends are deduplicated on resume.
- `-t threads` — worker threads (default: online CPU count).
- `-s start` — first candidate k to consider (default 2).
- `-e end` — stop after every candidate k ≤ end is resolved. Use
  disjoint `-s`/`-e` slices to split a search across machines.
- `-c file` — checkpoint: progress is saved every few seconds, and a
  killed or crashed run restarted with the same `-c` resumes where the
  frontier left off (a checkpoint overrides `-s`). Counters (terms
  found, candidates checked) accumulate across sessions.
- `-d mmax` — fixed trial-factor depth, q = 2km+1 for m ≤ mmax
  (default: auto-scales as ~4096·(k/1000)², capped at 2^33).
- `-v` — log per-candidate results (factors found, PRP timings) to
  stderr.

A status line goes to stderr showing where the search is in k, the
rates, and what the workers are doing. On a terminal it redraws in place
on a single line (CR, no LF) once a second; when stderr is redirected to
a file it writes one normal line every 10 seconds instead:

```
status: k=636247 (dispatched to 636263), 6 checked, 0.20 cand/s, 6.8 k/s, TF: 10, PRP: 0, lead k=636193 TF 31% (~66s left), elapsed 30s
```

`TF:`/`PRP:` count workers in each stage. The lead entry shows the
furthest-along check: during trial factoring the fraction of the m-range
swept, during a PRP test the true position inside the modular
exponentiation, each with an ETA. A lead shown as `confirming (BPSW)`
has passed the Fermat scan and is being verified.

In `-N` mode a hit prints statistics (cumulative across resumed
sessions):

```
HIT: 14561
a128344: HIT at k=14561 after 25.4s: (7^k-5^k)/2 is a 12306-digit (probable) prime
  searched k=13000..14561: 223 candidates checked (135 eliminated by trial factoring, 88 PRP-tested)
  session rates: 8.79 cand/s, 61.5 k/s
```

Examples:

```
./a128344 8                             # first 8 terms (~1 min on an M-series)
./a128344 -N -c a18.ckpt                # hunt a(18), resumable
./a128344 -N -s 1000000 -e 1100000 -c m1.ckpt        # one machine's slice
./a128344 -s 1000000 -e 1100000 -F m1.abc -c m1.ckpt # TF here, PRP in PFGW
```

## How it works

- Only prime k can be terms (for composite k = ab, (7^a − 5^a)/2 properly
  divides (7^k − 5^k)/2), so candidates are the primes, generated with
  primesieve.
- Any prime q dividing (7^k − 5^k)/2 satisfies q ≡ 1 (mod 2k): the
  multiplicative order of 7·5⁻¹ mod q divides the prime k and can't be 1.
  Each candidate is trial-factored over q = 2km+1: a bitmap sieve first
  removes the ~90% of m where q has a prime factor below 2^16 (any prime
  divisor that matters is itself ≡ 1 (mod 2k) and appears at its own
  smaller m), and the survivors get a 64/128-bit check of
  7^k ≡ 5^k (mod q). The depth auto-scales ~quadratically with k, because
  that's how the cost of the PRP test it might save grows.
- Survivors get a base-2 Fermat test, computed in 16-bit windows of the
  exponent so real progress can be reported while it runs. Composites
  essentially never pass a Fermat test, so this one pass settles them;
  the rare number that does pass is confirmed with GMP's BPSW test
  (`mpz_probab_prime_p`). Terms are strong probable primes — the same
  status as the published OEIS terms.
- Worker threads pull candidates from a shared queue; the main thread
  prints confirmed terms, writes export lines, and updates the
  checkpoint in sequence order as soon as every earlier candidate has
  been resolved. The checkpoint stores that frontier, so after a crash
  at most the work in flight above it is redone.

## pfgw_split.sh: running the survivors on a PFGW box

`pfgw_split.sh` turns a survivor file into a ready-to-run PFGW workload
for a multi-socket machine:

```
./pfgw_split.sh -n 24 -T 1 -S 2 -p pfgw64 survivors.abc
```

splits the survivors round-robin into 24 chunks (balanced mix of k sizes)
and writes `pfgw-run/` containing one directory per chunk plus:

- `run_all.sh` — launches every unfinished chunk, one PFGW instance per
  chunk, pinned to alternating sockets with
  `numactl --physcpubind --membind`, each instance on its own physical
  core (the default: 24 single-threaded instances on a dual 12-core
  Xeon, SMT siblings idle — measured to run each test at full solo
  speed, beating every multithreaded config by 2× or more in
  throughput). Run it under screen/tmux or nohup.
- `status.sh` — per-chunk done/left table and any PRP hits found.
- `resume.sh` — after a crash or reboot: removes already-tested k from
  each chunk (parsed from the captured PFGW output), then `./run_all.sh` continues where
  it left off.
- `manifest.txt` — chunk → socket/threads/count/k-range map.

The defaults (`-n 24 -T 1 -S 2`) fit a dual 12-core Xeon; set instances
to the physical core count and leave `-T 1`. Measured on a dual Xeon
8158 (wall-clock, k≈10^6): gwnum's threaded FFT scales poorly and its
extra coordinator thread oversubscribes an exact-T cpuset, so 24
single-threaded instances outran every multithreaded layout by 2× or
more, and scaled linearly with no memory-bandwidth penalty. Running a
second instance on each core's SMT sibling (48×) was a net 16% loss.
Note the pinning reads the real topology from `lscpu`: on many dual
sockets logical CPUs alternate between sockets (even = socket 0, odd =
socket 1), so hand-written contiguous ranges would straddle sockets and
cripple the FFT. (The trial-factoring stage is the opposite: its integer
loop gains ~37% from SMT, so run `a128344 -t` with all logical CPUs.)

## Searching beyond 10^6

Only 17 terms are known; Grantham & Granville report no others below
10^6 (2023), so the hunt for a(18) starts there. At k ≈ 10^6 each value
has ~845,000 digits and one Fermat test costs hours of a single core, so
for a serious search: split the range across machines with `-s`/`-e`,
give each slice its own `-c` checkpoint, and preferably use `-F` to do
only the (cheap, deep) trial factoring here and feed the survivor files
to PFGW on your biggest hardware. RAM is not a factor — each test's
working set is a few MB; core count and SIMD FFT throughput are what
matter. Note that many concurrent big-FFT tests can saturate memory
bandwidth before they run out of cores.
