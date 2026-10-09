# gapsieve — the constellation sieve

`gapsieve` searches for runs of consecutive primes whose gaps increase by two
(2, 4, 6, ... — see [README.md](README.md) and
[OEIS A016045](https://oeis.org/A016045)) without ever enumerating primes.
On a 14-core Apple Silicon machine where the enumeration tool `gapchain_ps`
manages ~26 billion integers/second near 2×10¹⁷, `gapsieve` sustains
**~145 trillion integers/second** in about 17 MB of memory. That turns
the hunt for `a(15)` from a ~15-year enumeration into roughly a day and a
half for the *entire* remaining 64-bit range. Positions are 128-bit, so the
search carries on past 2⁶⁴, up to 3.3×10²⁴, at about 70% of that speed. On a
Mac, `--gpu` adds a Metal GPU search that runs about eleven times as fast as the
CPU version on an M4 Max laptop, and nine times on an M2 Ultra Mac Studio.

This document explains the algorithm, the optimizations and the measurements
behind them, how to build and run it, and what to expect from a long search.

## Why enumeration loses at large heights

`gapchain_ps` must look at every prime, because it detects chains by watching
gaps in the prime stream. Near 2×10¹⁷ that means primesieve does the full work
of sieving every integer for primality, just so that ~2.5% of integers (the
primes) can each fail a gap test almost immediately. Nearly all of that work
answers a question nobody asked.

The chain condition is far more restrictive than primality, and the
restriction is *visible to small primes*. A chain head `p` needs

```
p + m(m+1)   prime, for m = 0 .. L-1     (offsets 0, 2, 6, 12, 20, 30, ...)
```

Reduce that system mod 3: the offsets hit the residues {0, 2}, so unless
`p ≡ 2 (mod 3)`, one of the first three members is divisible by 3. Two thirds
of all integers are eliminated by the number 3 alone. Every small prime adds
its own constraint, and they multiply. The constellation sieve takes this
observation seriously at industrial scale.

## The pipeline

A candidate `p` passes four stages, each more expensive and each seeing far
fewer candidates than the one before. Every candidate the sieve *eliminates*
is eliminated by an exact divisibility proof, and every chain the tool
*reports* is certified by a deterministic primality test. The probabilistic
tests in between only decide where the *work* goes, never what is true.

### Stage 1: the wheel — 99.87% gone before any work

For chain length 16, reduce the 16 offsets modulo each small prime and count
which residues of `p` survive:

| mod q | surviving residues | fraction |
| ----- | ------------------ | -------- |
| 2     | 1 of 2  (p odd)    | 0.500 |
| 3     | 1 of 3  (p ≡ 2)    | 0.333 |
| 5     | 2 of 5             | 0.400 |
| 7     | 3 of 7             | 0.429 |
| 11    | 5 of 11            | 0.455 |
| 13    | 6 of 13            | 0.462 |
| 17    | 8 of 17            | 0.471 |
| 19    | 9 of 19            | 0.474 |

By the Chinese Remainder Theorem these combine into a wheel of modulus
`W = 2·3·5·7·11·13·17·19 = 9,699,690` with only `12,960` admissible residue
classes — **0.134% of all integers**, or one candidate per 748. Candidates are
addressed *by class*, so the other 99.87% cost nothing, not even a skipped
loop iteration.

The wheel adapts to `--length`: it extends prime by prime while the class
count stays under 24,000 and the modulus under 3×10⁷.

### Stage 2: the pattern pre-sieve — every member, whole words at a time

The range is cut into segments of `2¹⁷ × W ≈ 1.27×10¹²` integers. Within a
segment each admissible class is an arithmetic progression

```
p = base + s0 + k·W,    k = 0 .. 2¹⁷-1
```

held as a 2¹⁷-bit (16 KB) bitmap. The job of stage 2 is to clear bit `k`
whenever some member `p + o` has a factor `q` between 23 and 199.

The key fact is that **the set of `k` a prime kills is the same in every class,
up to a rotation**. Member `p + o` is divisible by `q` exactly when

```
k ≡ W⁻¹·(−(base + s0) − o)   (mod q)
```

Substitute `x = k + r` with `r = W⁻¹·(base + s0) mod q` and the condition
becomes `x ≡ −W⁻¹·o (mod q)`, which no longer mentions the class at all. So
each prime gets one fixed bit pattern — bit `x` clear iff `x ≡ −W⁻¹·o` for any
of the L offsets `o` — and a class just reads that pattern starting at its own
rotation `r`. Several primes combine into one pattern with period
`P = q₁·q₂·q₃·…`, read at `R = CRT(r₁, r₂, r₃, …)`. With the default cap of
2²⁴ bits per pattern, primes 23..199 fit in **12 patterns** (~4 MB total,
shared read-only by all threads).

Building a class bitmap is then 12 pattern reads per 64-bit word: shift,
AND, done. All 12 patterns are ANDed into an 8-word accumulator held in
registers, so the bitmap is written once. There are no per-bit operations at
all.

Because a pattern costs the same no matter how many bits it clears, the
pre-sieve checks **all 16 members** rather than the handful a bit-striking
sieve can afford. That is where the power comes from: prime 23 now eliminates
12/23 = 52% of candidates instead of 6/23 = 26%, prime 31 eliminates 16/31
instead of 6/31, and so on down the line. After the 38 primes from 23 to 199,
about **13 candidates survive per 131,072-candidate class**, one per ~7.3
million integers. The old sieve struck 1,900 primes one bit at a time and
still left 164.

### Stage 3: batched screening — keep the multipliers busy

A survivor is known to have all 16 members free of factors up to 199. Each
member is prime with probability ~1/4.5, so the screen tests one member at a
time across all of a class's survivors, dropping each candidate at its first
composite member:

- **No trial division.** Every member is already 199-rough, so dividing by
  small primes would find nothing.
- **A specialised base-2 strong test.** With base 2, "multiply by the base"
  is a modular doubling, and 1 in Montgomery form is just `2⁶⁴ mod n`. So there
  is no `R²` constant and no 128-bit division, and the per-bit branch becomes a
  select, since the exponent's bits are unpredictable.
- **Four tests interleaved.** One strong test is a chain of ~60 *dependent*
  Montgomery squarings, which is latency-bound and leaves the multipliers idle
  most of the time. Running four independent tests in lockstep overlaps their
  chains. Four lanes measured best: two leave throughput on the table, and
  eight was slower again.

A "composite" verdict from the strong test is a proof. Only a candidate whose
16 members all pass moves on, and for L = 16 that means essentially only
genuine prime 16-tuples.

**Above 2⁶⁴** the same test runs in two-word Montgomery arithmetic (CIOS,
R = 2¹²⁸). A segment whose numbers all stay below 2⁶⁴ uses the one-word code
exactly as before, so nothing changes there. Candidates are stored as 64-bit
offsets from their segment's start, and the sieve only ever needs the start's
remainders, so it is unchanged at any height. Two details matter for speed.
1 in Montgomery form is now `2¹²⁸ mod n`, and a 128-bit remainder is a slow
library loop, so it comes from a floating-point estimate of the quotient plus
an exact correction. And four interleaved tests are still best.

### Stage 4: exact verification — nothing probabilistic survives it

A candidate that passes the screen gets `chain_verify`:

- **Interloper check.** Every odd number strictly between consecutive members
  must be composite, or the real consecutive-prime gaps are not 2, 4, 6, ....
  Each is screened with trial division and then base-2, both exact when they
  say "composite". Only a number surviving both gets the full deterministic
  test, and if it really is prime the chain ends at that gap.
- **Extension.** The walk continues past length L until a member fails or an
  interloper appears, up to L+80 steps, so the reported length is the true
  length.
- **Certification.** Every member is finally re-proved prime with the
  deterministic 12-base Miller-Rabin test (bases 2..37, correct for all
  n < 3.18×10²³ ≫ 2⁶⁴), truncating the chain if a base-2 pseudoprime ever
  slipped into the provisional walk. Above 2⁶⁴ the test uses the 13 smallest
  primes as bases (2..41). Sorenson and Webster (2015) proved that correct for
  every n below ψ₁₃ = 3,317,044,064,679,887,385,961,981, so that is the
  search's ceiling. Every number a chain check touches must stay below it.

## Correctness

**The sieve never kills a real chain.** A bit is cleared only when
`q | p + o` for a sieve prime `q` and a real member offset `o`. The sieved
region starts above `depth + 512`, so every member `p + o > q`. That makes
`p + o` a genuine composite and `p` genuinely not a chain head. A member
*equal* to a sieve prime can only occur below the sieved region, and that
strip is searched by `direct_scan`, a brute-force local prime sieve that walks
the actual consecutive-prime list.

**Tuning never changes results.** `--presieve` and `--depth` only move work
between stages whose rejections are all exact.

**Everything printed is deterministically certified.** A "prime" verdict from
a base-2 test is always re-established by the 12-base test (13 bases above
2⁶⁴) before anything is reported.

### Validation performed

- **Bit-exact pre-sieve.** For chain lengths 2, 4, 8, 15 and 16, 600 random
  classes each (393 million bits in all) at 10⁶, a(14), 8×10¹⁸ and 1.8×10¹⁹,
  every bitmap bit matched an independent brute-force predicate: plain
  `p mod q` arithmetic, with no CRT, patterns or inverses. There were 0
  mismatches. Survivor density matched the model exactly (13.4 vs 13.3 per
  class for L = 16).
- **Batched primality.** 12 million random odd numbers, one third of them
  above 2⁶³ to exercise the 129-bit carry and doubling paths, plus the
  largest 64-bit primes and base-2 strong pseudoprimes (2047, 3215031751,
  3825123056546413051). The 4-lane test agreed with the scalar test and never
  rejected a true prime.
- **Byte-identical output to `gapchain_ps`** on `-l 2 3..2e6` (14,871
  chains), `-l 4 3..1e6` (straddling the direct/sieve boundary),
  `-l 5 3..3e8`, `-l 7 3..2e11` (7,471 chains), `-l 6` at 10¹², and `-l 8` at
  2.2×10¹⁷. Also under deliberately odd tuning (`-p 0`, `-p 37 -d 40000`,
  `-p 2000`, ...).
- **A016045 from scratch.** `--first` reproduces every term through a(13).
  Most tellingly, `-l 15 -s 484511389338941 -e 2.3e17 --first` searched the
  entire range from a(13) upward and returned **a(14) = 221860944705726407**,
  the last known term, originally found with distributed computing. It took
  28.6 minutes at ~128 T/s and produced exactly one candidate, the correct
  one. That run exercises every stage at the real altitude and confirms that
  nothing below the true answer was missed or falsely reported.
- **Checkpoints.** Mid-run interrupt and resume gave output identical to an
  uninterrupted run (31,079 chains at 2.2×10¹⁷). A checkpoint in the older
  format, written for the exact a(15) command line, resumes correctly.
  `--first` restarts reach the right answer in every stop scenario (see
  *Stopping and restarting*).
- **Up to 2⁶⁴.** In the last 10⁹ below 2⁶⁴, every chain matched a brute
  force over primesieve's exact list of primes: 311 chains of length 5 of both
  shapes, 2,038 increasing chains of length 4, and a decreasing-only run right
  up to 2⁶⁴ − 21. Built with clang's unsigned-overflow sanitizer, searches
  there wrap only in the Montgomery arithmetic, where wrapping mod 2⁶⁴ is
  intended. The candidate positions and the chain checks never wrap.
- **Above 2⁶⁴: the arithmetic.** About 25,000 cases were checked against
  Python's exact integers:
  - Montgomery products for moduli up to 2¹²⁷.
  - `2¹²⁸ mod n`, concentrated around 2⁶⁴, where the quotient estimate is
    largest.
  - Exact parsing of numbers of any size.
  - The base-2 test and the full primality test on 4,533 numbers. These
    included 29 base-2 strong pseudoprimes and Carmichael numbers above 2⁶⁴,
    and an independent Baillie–PSW test agreed on every one.

  The cases passed under clang and GCC, and under the address and
  undefined-behaviour sanitizers.
- **Above 2⁶⁴: whole searches.** Chains matched a brute force (a small-prime
  sieve plus 13-base Miller–Rabin) in windows straddling 2⁶⁴, at 10²⁰ and
  10²², and right up to the ceiling. They also matched across segment
  boundaries, including a one-word segment followed by a two-word one. The
  checks covered `--first` for both shapes and an interrupted, resumed run.
  Below 2⁶⁴ the 128-bit build reports exactly the same chains as the 64-bit
  build did, and old checkpoints resume.
- ThreadSanitizer is clean on the counting and `--first` paths, above 2⁶⁴ too.

## How it got fast: the optimization history

Each step was measured single-threaded on one full segment at 2.2×10¹⁷,
`-l 16`:

| Version | s / segment | vs. start |
| ------- | ----------- | --------- |
| Bit-striking sieve, 6 offsets, primes to 16,384 | 5.40 | 1× |
| + pattern pre-sieve, all 16 offsets, primes to 199 | 0.159 | 34× |
| + 4-lane batched base-2 screen, larger patterns (2²⁴-bit cap) | 0.09 | 60× |

The full machine went from 1.97 T/s to **142–149 T/s** (~73×). Multi-core scales
slightly better than single-core here, because the old sieve's scattered
memory writes contended in shared caches and the patterns do not.

**What the profile said.** The original sieve spent 90% of its time (4.37 of
4.86 s) doing ~11 billion single-bit strikes per segment. Screening was only
10%. So the goal was never to make each strike cheaper. It was to stop
striking.

**Pre-sieve limit.** Once patterns carry all 16 members, a prime's marginal
value is set by the survivors it removes, measured against one more pattern
read per word. Per-segment time at a 2²⁴ cap:

| `--presieve` | s / segment |
| ------------ | ----------- |
| 151 | 0.103 |
| **199** | **0.091** |
| 263 | 0.103 |
| 331 | 0.119 |
| 509 | 0.211 |

**Bit-striking beyond the pre-sieve is now a loss.** Adding a strike stage
for primes 211..1024 on top of the default patterns raises per-segment time
from 0.091 s to 0.949 s. Striking costs far more than screening the few survivors it
removes, so the default `--depth` equals `--presieve` and the strike stage is
empty. It remains in the code for tuning experiments.

**Pattern size.** Larger patterns mean fewer pattern reads per word, at the
cost of cache footprint. On the full machine: 2²⁰ cap → 125 T/s, 2²² → 127
T/s, **2²⁴ → 142 T/s** (~4 MB of patterns, which fits comfortably in the
shared L2 cache).

**Fusing the pattern loop** (one accumulator in registers instead of an
in-place AND per pattern) measured within noise. The loop is bound by loads
and shifts, not by stores. It stays because it is no slower and scales to any
number of patterns.

**Earlier rounds** (still in effect): segments are sized so setup amortizes;
the exact verification screens interlopers with cheap proofs before the
12-base test (this turned a timeout into 17.6 s for `-l 8` at 2.2×10¹⁷); and
the Montgomery arithmetic handles n > 2⁶³ correctly. Testing caught two
silent-failure bugs there: a signed `2⁶⁴ mod n` that disabled the primality
test entirely, and a REDC sum that needed 129 bits.

**Above 2⁶⁴**, measured at `-l 17` on the laptop, where the one-word path
runs at ~160 T/s. The first two-word version ran at 103 T/s. A profile showed
13% of the time in the library's 128-bit division, computing `2¹²⁸ mod n`
once per test. Replacing it with a floating-point quotient estimate and an
exact correction brought the rate to ~115 T/s. Other changes measured within
noise:

- a dedicated squaring routine, since the compiler already shares the cross
  product;
- two interleaved tests instead of four (eight was 17% slower);
- a deeper pre-sieve (229 was within noise, and 263 and above were slower).

## Building

No dependencies. It needs plain C with pthreads and `__uint128_t`, so any
modern gcc or clang on macOS or Linux works:

```bash
cc -O3 -pthread gapsieve.c -o gapsieve
```

On a Mac, `make` also builds in the GPU search. See
[*On the GPU (macOS)*](#on-the-gpu-macos) below.

## Running

The CLI is the same as `gapchain_ps`, plus two tuning knobs:

```
gapsieve -l LEN -e END [options]

  -l, --length N      chain length (primes per run), 2..1000   [required]
  -e, --end N         highest chain head to consider, or max   [required]
  -s, --start N       lowest chain head (default 3)
  -t, --threads N     worker threads (default: online CPUs)
  -g, --gaps WHICH    inc (2, 4, 6, ...; A016045; default), dec (..., 6, 4, 2;
                      A263049) or both
  -f, --first         report only the smallest qualifying chain (one of each
                      shape with --gaps both), then stop
  -c, --checkpoint F  checkpoint to F periodically and on SIGINT/SIGTERM
  -i, --interval S    checkpoint seconds, 0 = exit only (default 60)
  -r, --resume        resume from F (required if F exists)
  -p, --presieve N    primes up to N sieved against all L members with
                      bit patterns (default 199)
  -d, --depth N       largest sieve prime; primes above --presieve are
                      struck bit by bit (default 199, i.e. none)
      --log FILE      also append every message, timestamped, to FILE,
                      with a status line once a minute
      --gpu           also search on the GPU (macOS builds); -t then
                      sets the CPU threads beside it (default half
                      the cores, 0 for the GPU alone)
  -q, --quiet         suppress progress output
```

`-p` and `-d` are tuning only; output is identical at any values. Numbers
accept `1e14`, `2.5e15`, `300T`, `1_000_000`, or plain digits, all converted
exactly, so `1e20` is exactly 10²⁰. Chains print to stdout, progress and
statistics to stderr.

`-e max` searches as high as the proven primality test allows. Every number
checked for a start p must stay below ψ₁₃ ≈ 3.317×10²⁴, so gapsieve lowers
the end to the last start where that holds, and says so. That is ψ₁₃ − 9,507
at `-l 17`, because the increasing search follows a chain up to 80 members
past `-l` to report its full length; `--gaps dec` alone reaches ψ₁₃ − 273.
Any end above that limit is lowered the same way.

### Status line and log

On a terminal, a status line repaints itself four times a second:

```
0.06%  eta 1d10:14:58  no hit yet  144.10T/s  at 2.33223e+17  elapsed 00:01:09
```

The fields are percent done, time remaining, the hit status (under `--first`,
`no hit yet` or `candidate <p>`; otherwise the chain count), speed, position
reached, and elapsed time. The line never exceeds the terminal width, because
a wrapped line would scroll instead of overwriting itself. On a narrow window
the least important fields are dropped whole rather than cut off, so a prime is
never shown truncated.

`--log FILE` keeps that live line on the terminal and appends a complete,
timestamped record to `FILE`: the command line and process ID, resume
details, the settings, a status line once a minute, every candidate the moment
it is found, the winning chain, the final summary, and the exit code:

```
2026-09-23 10:49:24  started (pid 84817): ./gapsieve -l 16 ... --log a15.log
2026-09-23 10:49:24  resuming from a15.ck at 1005/13983585 segments (0.01%)
2026-09-23 10:50:24  0.06%  eta 1d10:14:35  no hit yet  144.13T/s  at 2.31811e+17  elapsed 00:01:00
2026-09-23 10:50:34  interrupt: finishing segments in flight...
2026-09-23 10:50:34  exit 130
```

The file is opened in append mode, so each restart continues the same log.
With stderr redirected instead (`2> file`), there is no live line: stderr is
then not a terminal, so it gets one plain status line a minute.

### The a(15) and a(16) searches

A016045's a(15) was found with this command, starting at a(14) because the
sequence is non-decreasing:

```bash
./gapsieve -l 16 -s 221860944705726407 -e 1.8e19 --first -c a15.ck --log a15.log
```

The answer is 1397398433200922807, found after 6.6% of the range. A263049's
a(15), 253253149671986953, came out of validating the decreasing search.
README.md's [*New results*](README.md#new-results) has the details and the
evidence for both.

`--first` stops at the smallest hit, with the guarantee that in-flight earlier
segments finish before a winner is declared. A segment is about 0.1 s of one
thread's work, so checkpoints are fine-grained.

Both sequences' next terms need 17-prime chains, each starting from its own
a(15). A263049's can dip by one gap, so its search starts 32 lower:

```bash
./gapsieve -l 17 -s 1397398433200922807 -e 1.8e19 --first -c a16.ck --log a16.log
```

```bash
./gapsieve -l 17 -s 253253149671986921 -e 1.8e19 --first --gaps dec -c d16.ck --log d16.log
```

Neither search found a chain. Follow-up runs took both to 2⁶⁴ (the last with
the 64-bit build's `-e max`, plus a separate check of the few thousand starts
at the very top), so both a(16) terms are above 2⁶⁴. README.md's
[*The a(16) searches*](README.md#the-a16-searches) has the details.

With 128-bit positions the searches continue above 2⁶⁴. Starting at
1.8446744×10¹⁹ overlaps the finished range by 7×10¹⁰, which costs under a
second:

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 1e20 --first -c a16w.ck -r --log a16w.log
```

```bash
./gapsieve -l 17 -s 1.8446744e19 -e 1e20 --first --gaps dec -c d16w.ck -r --log d16w.log
```

`--gaps both` could run these as one search now that the lengths match. But it
would have to start at the lower of the two starts, 2.5×10¹⁷. That makes the
increasing half re-scan 10¹⁸ where the a(15) search has already ruled out any
chain. Two separate runs, one per machine, are better.

### Stopping and restarting

The checkpoint is rewritten every 60 seconds (`-i` to change) and on Ctrl-C
or SIGTERM. Each write goes to `a15.ck.tmp`, is fsynced, and is then renamed
into place, so a crash or power cut never leaves a damaged file. The most that
is lost is the last interval of work. To restart, run the same command with
`-r` added:

```bash
./gapsieve -l 16 -s 221860944705726407 -e 1.8e19 --first -c a15.ck -r --log a15.log
```

Under `--first` the checkpoint also records the best candidate found so far
(a `best` line). Resuming then does the right thing wherever the run was
stopped:

- **Before any hit:** the search continues where it left off.
- **After a hit, with earlier segments still unfinished:** only those earlier
  segments are rescanned. If one of them holds a smaller chain, that chain
  replaces the recorded one.
- **After the answer was settled:** the checkpoint is reported as already
  complete and the answer is printed again, with no further searching.

Resuming a finished `--first` checkpoint used to be unsafe in `gapsieve`. The
segment holding the hit was marked done but the hit itself was not saved, so
a restart searched on past the answer and could report a larger prime as the
first. Tests covered each case above: resuming finished runs; checkpoints
stopped with a hit recorded, including one planted with a fake larger
candidate, which the real earlier chain overturned; and runs hard-killed
(SIGKILL) the moment the hit appeared. `gapchain_ps` never had the problem,
because it abandons the segment holding the hit, so its checkpoint always
stops before the answer.

**If a run is already in progress with an older build**, stop it with
Ctrl-C, rebuild, and resume with `-r` as above. The segment grid for `-l 16`
is unchanged, older checkpoint files are read correctly, and the work already
done is kept.

When a candidate appears, the log shows it immediately:

```
candidate: <p> (length 16, segment <s>); finishing earlier segments
```

Cross-check any hit with the independent enumeration tool:

```bash
./gapchain_ps -l 16 -s <p-1000> -e <p+1000>
```

## Expected performance

Measured on 14-core Apple Silicon, `-l 16`, defaults:

| Metric | Value |
| ------ | ----- |
| Single thread | ~0.09 s per 1.27×10¹² segment (~14 T/s) |
| Full machine (14 threads) | **142–149 T/s** sustained |
| Peak memory | 16.6 MB (patterns ~4 MB, shared) |
| vs. previous `gapsieve` | ~73× |
| vs. `gapchain_ps` at 2.2×10¹⁷ | ~5,000× |

Throughput is nearly flat with height (measured 149 T/s at 2.2×10¹⁷, 145 at
8×10¹⁸, 137 at 1.7×10¹⁹, where numbers above 2⁶³ take a slightly longer
arithmetic path), so ETAs extrapolate honestly. (`-l 15` runs ~128 T/s, because
15-member patterns leave more survivors.) Above 2⁶⁴ the two-word arithmetic
costs about 30%: at `-l 17` the laptop measured 163 T/s at 1.8×10¹⁹,
115 T/s at 2×10¹⁹ and 110 T/s at 10²⁴. At ~145 T/s:

| If a(15) is near | time to find it |
| ---------------- | --------------- |
| 1×10¹⁸ | ~1.5 hours |
| 8×10¹⁸ (growth-trend estimate) | ~15 hours |
| 1.8×10¹⁹ (end of the range) | ~35 hours |

In the event, a(15) was at 1.4×10¹⁸. The search found it after 6.6% of the
range; house1, an x86 machine running at about 112 T/s, finished it off in
under three hours.

Both a(16) searches came up empty below 2⁶⁴, which is why gapsieve now
works in 128 bits. With all three machines at about 70% of their 64-bit
speed, roughly 350 T/s together, reaching 10²⁰ takes about three days per
sequence.

## On the GPU (macOS)

On a Mac, `--gpu` adds a search on the GPU, through Apple's Metal. At `-l 17`
above 2⁶⁴, near 10²⁰ where the a(16) searches run:

| Configuration | M4 Max laptop (32-core GPU) | M2 Ultra Mac Studio (60-core GPU) |
| ------------- | --------------------------- | --------------------------------- |
| CPU only, all cores | 108–111 T/s | 133–138 T/s |
| GPU alone (`--gpu -t 0`) | — | 1,118 T/s |
| GPU + all cores but two (the default) | ~1,280 T/s | ~1,280 T/s |

That's about eleven times the CPU version on the laptop and nine on the
Studio. Both shapes run at the same speed. (The laptop went straight from its
tests into its search, so its GPU alone hasn't been timed.)

**How it works.** The GPU takes units of 16 consecutive segments from the
same queue as the CPU threads, so checkpoints, `--first`, the statistics and
the status line work unchanged. A segment is exactly 2¹⁷ turns of the wheel,
so consecutive segments share their classes, and a unit is just a longer
stretch of each class. Each unit and shape runs in three steps:

1. **Offsets.** One thread per class, pattern group and segment works out
   where the class starts in that group's pattern, as `class_bitmap` does.
2. **Sieve.** A threadgroup takes 1024 32-bit words of one class. First
   each thread ANDs the first three rotated patterns into its four words.
   Those patterns hold the smallest primes and strike the most: afterwards
   about 90% of words are zero. The words still standing are gathered in
   threadgroup memory (one atomic per 32 threads, through SIMD prefix sums),
   and each gets a thread of its own for the other nine patterns, so no lane
   works on a zero word. Survivors go on the candidate list. Apple GPUs are
   32-bit machines, so 32-bit words sieve much faster than 64-bit ones did:
   1.69 ms per segment instead of 2.87. The two phases took the Studio's
   sieve from 0.97 ms per segment to 0.84. `GAPSIEVE_GPU_DENSE` sets how many
   patterns come before the gathering: 3 measured best on the Studio (1,118
   T/s, against 1,056 with 2 and 1,098 with 5).
3. **Screen.** Three rounds each test one member of every remaining candidate
   with the base-2 strong test, and keep the candidates that pass. Then one
   kernel tests all the remaining members of the few candidates left (under
   1%), stopping at the first failure. The GPU runs threads in groups
   of 32, and about a fifth of candidates pass each test. Testing each
   candidate's members in one thread from the start left most of each group
   idle, waiting on the few that kept passing: member by member took the
   screen from 1.99 ms per segment to 1.46. But each round costs a dispatch
   and a barrier whatever its size. Once few candidates remain, the idle
   lanes cost less than more rounds would: three rounds then the tail
   measured best (Studio: 773 T/s, against 687 with one round and 763 with
   all seventeen). The GPU sizes each round itself (indirect dispatch), so a
   unit runs as one command buffer with no round trips to the CPU.

The CPU then proves the rare survivors with `chain_verify`, as it does its
own. The GPU's base-2 test is `gapsieve.c`'s Montgomery arithmetic, redone
for the GPU. The first version ported the CPU's two 64-bit words. But Apple
GPUs multiply 32 bits at a time, so each 64-bit product costs several
multiplies. Every number searched is below 2⁸², so the test now works in
three 32-bit limbs, whose 64-bit sums can't overflow, with about half the
instructions. And with R = 2⁹⁶ far above 4n, Montgomery products of values
below 2n stay below 2n, so they're reduced only for the final comparisons.
That took the Studio's screening from 0.66 ms per segment to 0.26. Apple GPUs
have no 64-bit floating point, so `R mod n` comes from repeated doubling
instead of the CPU's floating-point estimate. In a first benchmark, the
two-word version already ran the base-2 test about five times as fast as all
14 CPU cores of the laptop: 178 million tests a second against 38.

**Units, one at a time.** The first GPU version sent single segments, two in
flight. It reached 592 T/s on the laptop, but only 274–411 T/s on the Studio,
varying from run to run, and adding CPU threads more than doubled the
Studio's GPU. A profile (`GAPSIEVE_GPU_PROF=1`) showed why. A segment's
screening took almost as long as its sieve (1.20 ms against 1.30 on the
laptop), though it does a fraction of the arithmetic, because its 17 rounds
of two dispatches each spend most of their time waiting on barriers. The
two-die M2 Ultra waited longer. Units of 16 segments give every dispatch 16
times the work, and the tail kernel cuts 34 dispatches to 8.

A second finding: units running side by side slow each other badly (on the
Studio, two at once ran at 263 T/s and one at a time at 735), and the first
version's two segments in flight had been overlapping too. So units now run
strictly one after another: every unit writes one small shared buffer, and
Metal's hazard tracking starts each only once the one before has finished.
The next unit is always queued, so the GPU never waits on the CPU. GPU
alone, at `-l 17` near 10²⁰:

| Version | Laptop | Studio |
| ------- | ------ | ------ |
| single segments, two in flight | 592 T/s | 371 T/s |
| 8-segment units, one in flight | 593 T/s | 735 T/s |
| 8-segment units, two side by side | 417 T/s | 263 T/s |
| 8-segment units, next one queued | 604 T/s | 754 T/s |
| 16-segment units, next one queued | 604 T/s | 773 T/s |
| three-limb arithmetic, two-phase sieve (now) | — | 1,118 T/s |

The laptop gains little: its GPU was less starved, and its time now goes to
the sieve and the arithmetic. The Studio doubled. With its GPU fed, CPU
threads add to the total instead of competing: 805 T/s with 6 beside it,
855 with 12, 891 with 18 and 899 with 21. Hence the new default of all cores
but two, leaving room for the thread that feeds the GPU. (On the laptop,
whose CPU and GPU share one power budget, 6 to 15 threads all gave about
640 T/s.) The three-limb arithmetic and the two-phase sieve came next, and
with them both machines search at about 1,280 T/s with their CPU threads:
the Studio up from 909, the laptop up from 640.

Short chains leave many more candidates: `-l 7` leaves about 68 million per
segment, against 103,000 at `-l 17`. So the unit shrinks to keep its lists
near 8 million entries, sized from the patterns in advance. The estimate is
exact on average; at `-l 17` it predicted 102,968 candidates per segment,
and 102,958 were measured.

**Failures are handled, not hoped away.** Every command buffer's status is
checked. A failed unit is redone, and after three failures the CPU searches
its segments instead. A candidate list that overflows is enlarged and the
unit redone. Two environment variables exercise those paths for testing:
`GAPSIEVE_GPU_FAIL=N` treats every Nth command buffer as failed, and
`GAPSIEVE_GPU_CAP=N` starts with N-entry lists. Four more override the
tuning: `GAPSIEVE_GPU_UNIT` (segments per unit, 16), `GAPSIEVE_GPU_SPLIT`
(rounds before the tail kernel, 3), `GAPSIEVE_GPU_DEPTH` (units in flight,
2), `GAPSIEVE_GPU_DENSE` (patterns before the sieve's gathering, 3) and
`GAPSIEVE_GPU_ORDER` (0 lets units overlap, 2 orders only their sieves).
None of them changes the results.

**Validation.** Before its first real run, the GPU build was checked against
the CPU build and against brute force:
- **The base-2 test:** identical answers on 8.4 million numbers near 5×10¹⁹
  and 10²¹.
- **Segment by segment:** identical sieve candidate lists and screen
  survivors, at `-l 17` for both shapes, and at `-l 8` and `-l 9`, where
  survivors come through.
- **Whole runs:** the same chains as the CPU build, with the GPU alone and
  with CPU threads. That held with every third command buffer failed on
  purpose, with every one failed so that the CPU took over, and with tiny
  lists. The runs included all 441,634 chains of `-l 7` from 10¹⁵ to
  1.05×10¹⁵ (below 2⁶⁴, which the GPU searches in two words), all 23,234
  chains of `-l 9` at 1.2×10¹⁹, and the known terms a(11) and a(10) under
  `--first`.
- **As for the 128-bit build:** the brute-force windows below and above 2⁶⁴,
  across segment boundaries and up to the ceiling; `--first` above 2⁶⁴; and
  an interrupted, resumed run.
- **Units:** all of the above again, plus the regression runs with units of
  1, 2, 3, 5 and 64 segments, with no screening rounds and with all
  seventeen, with one and five units in flight, with units free to overlap
  and with only their sieves ordered, with forced failures and tiny lists in
  odd-sized units, and an interrupted run in 3-segment units beside two CPU
  threads. The Studio reproduced the laptop's results for every test.
- **Three limbs and two phases:** the three-limb test gave the same answer as
  the CPU's two-word test on 18.9 million numbers. They were spread over nine
  ranges from 701 up to 2⁹³, across the 32- and 64-bit limb boundaries, plus
  24 special cases: primes, and base-2 strong pseudoprimes such as
  3,825,123,056,546,413,051 and 3,317,044,064,679,887,385,961,981, which the
  test must pass. Then the whole suite ran again on the laptop, adding 1, 2
  and 12 patterns before the gathering (84 runs identical, every brute-force
  window matching), and on the Studio (42 of 42 identical).

**Other Macs.** Any Apple Silicon Mac should work. Before using one, check it
with the fingerprint command in README.md, run with `--gpu -t 0`. Besides the
laptop's M4 Max and the Studio's M2 Ultra, an M1 Max MacBook Pro (32-core GPU)
reproduced the laptop's results in 19 of 19 checks. It searched at 614 T/s
with 8 CPU threads beside its GPU, against 59 on its CPU alone, but only for
its first 20 minutes: once warm it settled at about 438 T/s. The GPU slows
itself; with 2 CPU threads instead of 8, to cut the heat, it gave 428. The
AMD GPUs in Intel Macs do 64-bit integer arithmetic slowly and are untested.

## On NVIDIA GPUs (CUDA)

`make CUDA=1` builds the same GPU search for Linux with an NVIDIA GPU. It was
written for two DGX Sparks, whose GB10 pairs a 20-core ARM CPU (10
Cortex-X925, 10 Cortex-A725) with a Blackwell GPU (compute capability 12.1)
and 121 GB of memory shared between them.

**How it is built.** `gapsieve_cuda.cu` holds the kernels, ported from
`gapsieve.metal`: the same units of 16 segments, the two-phase sieve (the
words still standing after the first three patterns are gathered in shared
memory, one atomic per 32-thread warp through shuffle-based prefix sums), and
the three-limb base-2 test. `gapsieve_cuda_host.c` is `gapsieve_gpu.m`'s
host side in C. It includes `gapsieve.c` unchanged and feeds the GPU from the
same queue of segments as the CPU threads. `gapsieve_cuda.h` is the small C
interface between the two, so the kernels compile with `nvcc` and the host
with the C compiler.

Two things differ from the Metal version. The screening kernels read their
candidate count from GPU memory and loop over it, so there is no indirect
dispatch and no round trip to the host. And every unit goes on one CUDA
stream, which runs units strictly in order with the next one queued; on
Metal that took a shared buffer and hazard tracking.

**Speed.** At `-l 17` near 3×10²⁰:

| Configuration | atom1 | atom2 |
| ------------- | ----- | ----- |
| GPU alone (`--gpu -t 0`) | 1,428 T/s | 1,338 T/s |
| GPU + 2 CPU threads | 1,426 T/s | 1,340 T/s |
| GPU + 6 CPU threads | 1,319 T/s | 1,282 T/s |
| GPU + 18 CPU threads (the default) | 1,288 T/s | 1,254 T/s |
| GPU alone, 32-segment units | 1,378 T/s | 1,324 T/s |

On the Spark the CPU and GPU share one power budget and one memory, so CPU
threads slow the GPU by more than they add. Run it with `--gpu -t 0`. In
production the two ran at 1,470 and 1,360 T/s: each about as fast as the M4
Max laptop, and faster than the Mac Studio. atom2 was consistently about 6%
slower than atom1.

**Validation.** Before its first real run, on both Sparks:
- **The fingerprint** (all 234 chains of length 7 in a window above 2⁶⁴):
  `71fe408bdec8f5c6`, as on every other machine.
- **The base-2 test:** the same 18.9 million numbers the Metal version was
  checked on, across nine ranges from 701 to 2⁹³, the 32- and 64-bit limb
  boundaries, strong pseudoprimes, and the members of A263049's a(16). All
  agreed with the CPU's `sprp2w`.
- **Whole runs:** all 14 regression ranges the Macs were checked on, with the
  GPU alone and with CPU threads, forced failures, tiny lists, and different
  unit sizes, rounds, depths and dense patterns. 42 of 42 results on each
  Spark were identical to the CPU build's.

## Searching for both shapes

`--gaps dec` searches for A263049's decreasing chains, and `--gaps both`
searches for both kinds in one pass.

**How it works.** A decreasing chain of L primes from `p` has members
`p + m(2L-1-m)`, the increasing pattern reflected end to end. The
reflection has a useful consequence. At every wheel prime the two shapes rule
out exactly as many residues, so they need the same wheel, have the same
12,960 classes, and leave the same density of survivors. The decreasing
search reuses everything: the wheel, the segment grid, the pattern groups and
their CRT constants. Only its bit patterns, its residue list and its final
check are its own. That check verifies the exact L-prime pattern, then
extends the chain **backwards**. A decreasing chain always ends in a gap of 2,
so it can only grow at the front. `--first` ranks by the head of the exact
pattern, since that is what A263049 records.

**What it costs.** Measured on the full machine at 2.2×10¹⁷, `-l 16`:

| `--gaps` | Throughput |
| -------- | ---------- |
| `inc` | 148 T/s |
| `dec` | 142 T/s |
| `both` | 70 T/s |

The increasing search is unchanged by the addition. `both` covers the range
at half speed because it does twice the work: every class is sieved once per
shape. The two searches share almost nothing that costs anything. The
per-class rotation arithmetic is common, but it is negligible; the sieve
passes need different survivor sets and cannot be shared. So `both` takes as
long as the two searches run back to back. For 16-prime chains from 2.2×10¹⁷,
that is about 70 hours to cover the whole range, against about 35 hours for
either shape alone.

What `both` changes is *when* each answer arrives: each shape's answer comes
at half its solo speed. Run `both` for one command and one log. Run the shapes
separately if one answer matters more.

**Validation.**
- **A263049 from scratch.** `--first --gaps dec` reproduces every known term,
  a(1)–a(13). That includes the a(9)/a(10) pair, where a(10) is *smaller*: the
  10-prime pattern is the tail of the 11-prime chain. The run reports that
  the pattern ends a longer chain.
- **A263049's a(14), re-derived.** Searching `-l 15` from a(13) − 30 upward
  took 32.5 minutes at ~128 T/s. It returned exactly one candidate,
  **253253149671986983**, the last known term. It also reported that this
  pattern ends a longer chain of 16 primes from 253253149671986953, which
  makes that prime **A263049's a(15)**. The claim was confirmed with
  primesieve's exact sieve, a separate Python primality check and
  `gapchain_ps`.
- **Independent cross-check.** `gapchain_ps` gained `--gaps` too. It detects
  decreasing chains by walking the real prime sequence, keeping a ring of
  earlier primes to extend backwards, including past the start of a segment.
  Its output matched `gapsieve` exactly: about 100,000 decreasing chains, from
  `-l 2` up to 3×10⁹ through `-l 8` at 2.2×10¹⁷, in `dec` and `both` modes.
- **Bit-exact sieve.** 524 million bitmap bits matched a brute-force predicate,
  covering both shapes, five lengths, and heights up to 1.8×10¹⁹.
- **`--first` with both shapes.** Terms for both sequences came out of single
  runs, including one where a larger decreasing candidate appeared first and
  was then overturned.
- **Checkpoints.** Both shapes' candidates are stored and restored, and
  interrupted runs resume to the same answers. An increasing-only checkpoint
  refuses to resume as a `both` search. ThreadSanitizer is clean.
