# pdigits — digit counts in the first 10^N primes

Computes the next term of **ten OEIS sequences at once**: A119291 (digit 0)
through A119300 (digit 9), "total number of digit-d occurrences in the first
10^n primes". All ten currently stop at n = 12; one run at `-n 13` produces
a(13) for every one of them. Exact integer arithmetic throughout — nothing
probabilistic.

## OEIS sequences

**Sequences this tool can extend (new terms):**

| Sequence | Definition |
|----------|------------|
| [A119291](https://oeis.org/A119291) | number of 0 digits in the first 10^n primes |
| [A119292](https://oeis.org/A119292) | number of 1 digits in the first 10^n primes |
| [A119293](https://oeis.org/A119293) | number of 2 digits in the first 10^n primes |
| [A119294](https://oeis.org/A119294) | number of 3 digits in the first 10^n primes |
| [A119295](https://oeis.org/A119295) | number of 4 digits in the first 10^n primes |
| [A119296](https://oeis.org/A119296) | number of 5 digits in the first 10^n primes |
| [A119297](https://oeis.org/A119297) | number of 6 digits in the first 10^n primes |
| [A119298](https://oeis.org/A119298) | number of 7 digits in the first 10^n primes |
| [A119299](https://oeis.org/A119299) | number of 8 digits in the first 10^n primes |
| [A119300](https://oeis.org/A119300) | number of 9 digits in the first 10^n primes |

All ten have 12 published terms (n = 1..12); `-n 13` yields a(13) for each.

**Sequences used for verification (embedded reference tables):**

| Sequence | Definition | Role |
|----------|------------|------|
| [A006988](https://oeis.org/A006988) | 10^n-th prime | every index boundary must land on it |
| [A119290](https://oeis.org/A119290) | total digits in first 10^n primes | known to n=21, so it independently checks the digit-sum of the NEW terms |
| [A006880](https://oeis.org/A006880) | pi(10^n) | prime counts at every power of ten passed |
| [A091634](https://oeis.org/A091634)–[A091643](https://oeis.org/A091643) | primes < 10^n with no digit d (d = 0..9) | verified at every 10^n the sweep covers (n <= 14 at N=13) |

The A091634 family is known to n = 19, so at our scopes those are pure
cross-verification, not new terms (extending them would need enumeration to
10^20 — far out of reach). The complement family A091644-style values follow
from a(n) + complement = pi(10^n) and are not printed separately.

## Build

Requires [primesieve](https://github.com/kimwalisch/primesieve) and pthreads.

    brew install primesieve      # once
    make                         # builds ./pdigits
    make test                    # selftest: full pipeline at N=8, ~1 s

## Run

    ./pdigits -n 13 -c pdigits.ckpt -o report13.txt        # the campaign (~hours)
    ./pdigits -n 13 -c pdigits.ckpt -o report13.txt -r     # resume after a stop

Ctrl-C is always safe: completed chunks are already in the checkpoint log,
and `-r` redoes only what is missing. Ctrl-T (SIGINFO) prints an immediate
status line. Exit code: 0 = complete with all checks passed, 2 = complete
with mismatches, 130 = interrupted.

## Command line options

| Option | Default | Meaning |
|--------|---------|---------|
| `-n N` | 13 | target: count the first 10^N primes (N in 1..16) |
| `-t THREADS` | online cores | number of worker threads |
| `-s E` | 9 | chunk span 10^E (E in 9..12; use `-s 10` for N=14+ to keep the chunk table small) |
| `-c FILE` | `pdigits.ckpt` | checkpoint log; `none` disables checkpointing |
| `-o FILE` | (none) | also write the final report/summary to FILE |
| `-r` | off | resume from the checkpoint (must match `-n` and `-s`) |
| `-u SECS` | 10 | status line interval |
| `-q` | off | quiet (suppress status lines) |
| `--selftest` | — | run the full pipeline at N=8 and verify against all published values |
| `--benchmark[=S]` | 45 s | measure this machine's throughput and print a comparable SCORE; honors `-t` (use `=` form for a custom duration, e.g. `--benchmark=60`) |
| `--chunk-min I` | 0 | only process chunks >= I (multi-machine split) |
| `--chunk-max I` | all | only process chunks < I (multi-machine split) |
| `-h` | — | help |

## Benchmarking machines

    ./pdigits --benchmark          # 45 s, all cores
    ./pdigits --benchmark=60 -t 8  # custom duration / thread count

Every machine tallies the identical fixed region (stripes at 10^17,
N=16-representative height) for the given wall time; the printed
**SCORE (Gprimes/s)** is directly comparable across machines and `-t`
settings, and drives the rough full-run projections. Benchmark on an
otherwise-idle machine — a competing CPU load deflates the score
proportionally (the same laptop read 0.385 while a ~10-core job ran).
Reference scores, 14 threads:

    M1 Max laptop (M1MaxJS): 0.656 Gprimes/s   (N=16 solo: ~180 days)

Use the score ratio between machines to split `--chunk-min/--chunk-max`
ranges proportionally.

## Multi-machine runs

Chunk tallies are independent, so a big run can be split across machines:
give each machine the same `-n`/`-s` plus a disjoint `--chunk-min`/`--chunk-max`
range and its own checkpoint file. When all ranges are complete, merge the
logs on one machine and rerun with `-r` for the final report:

    # machine A                                # machine B
    ./pdigits -n 16 -s 12 -c a.ckpt \          ./pdigits -n 16 -s 12 -c b.ckpt \
        --chunk-max 200000                         --chunk-min 200000

    # then on either machine:
    grep '^C ' b.ckpt >> a.ckpt
    ./pdigits -n 16 -s 12 -c a.ckpt -o report16.txt -r

Balance ranges by remembering later chunks are somewhat slower (higher sieve
region); give the faster machine the tail. Interrupt/resume works per machine
as usual.

## Method

The range up to ~p(10^N) is cut into span-aligned chunks. Worker threads pull
chunks and tally prime count + per-digit counts + digit-presence masks using a
primesieve iterator. Digit counting is SWAR: two lookup tables cover the low
7 digits (4+3 split, byte-lane packed counts flushed every 30 primes), and the
shared high-digit prefix of each 10^7 window is counted once and multiplied by
the window's prime count. A 1024-bin presence-mask histogram per chunk yields
the "contains no digit d" counts for the A091634 family. Chunk 0 and all
boundary rescans use a plain exact per-digit loop. Completed chunks append to
the checkpoint log; exact 10^n index boundaries are located afterward by
rescanning the single chunk containing each one.

## Verification at completion

The final report (stdout, plus `-o FILE`) prints, for every n up to N:
the boundary prime vs A006988, the total digit count vs A119290, all ten
per-digit counts vs A119291..A119300 (marked NEW past n=12), the
A091634..A091643 row for every 10^n value boundary covered, and pi(10^d)
vs A006880 — ending with an overall ALL CHECKS PASSED / FAILURES verdict.
A new term is only trustworthy if every check says OK.

## Timing (M-series Mac, 14 threads)

- N=10: ~12 s (used as a full-verification shakedown)
- N=13: a few hours — the campaign target
- N=14: ~1.5 days (completed 2026-08-21); `-s 10`
- N=15: ~2 weeks-scale (completed 2026-09 on bigger hardware); `-s 11`
- N=16: ~13-14x the N=15 work — months; `-s 12`, and consider a
  multi-machine split (above). This is the last N with an independent
  A119290 digit-sum check (A119290 is published through n=21).
