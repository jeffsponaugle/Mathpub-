# digitprimes

Extends the OEIS family of sequences counting the first 10^n primes by decimal
digit content:

| digit d | no digit d | at least one digit d |
|---------|------------|----------------------|
| 0 | [A231412](https://oeis.org/A231412) | [A231726](https://oeis.org/A231726) |
| 1 | [A228413](https://oeis.org/A228413) | [A231787](https://oeis.org/A231787) |
| 2 | [A228414](https://oeis.org/A228414) | [A231788](https://oeis.org/A231788) |
| 3 | [A228415](https://oeis.org/A228415) | [A231789](https://oeis.org/A231789) |
| 4 | [A228416](https://oeis.org/A228416) | [A231790](https://oeis.org/A231790) |
| 5 | [A228417](https://oeis.org/A228417) | [A231792](https://oeis.org/A231792) |
| 6 | [A228418](https://oeis.org/A228418) | [A231793](https://oeis.org/A231793) |
| 7 | [A228419](https://oeis.org/A228419) | [A231794](https://oeis.org/A231794) |
| 8 | [A228420](https://oeis.org/A228420) | [A231795](https://oeis.org/A231795) |
| 9 | [A228421](https://oeis.org/A228421) | [A231796](https://oeis.org/A231796) |

The "at least one d" sequences are exact complements: a(n) = 10^n − (no-d
count), so a single pass produces all 20 sequences at once. All were known
through a(12); the next term a(13) requires digit statistics over the first
10^13 primes (all primes up to ~3.24 × 10^14).

## How it works

- [primesieve](https://github.com/kimwalisch/primesieve) iterates every prime
  in order.
- Each prime gets a 10-bit digit-presence mask (table lookup, 4 digits at a
  time); a 1024-entry histogram per work chunk is reduced to per-digit
  "missing" counts, so the hot loop is one lookup chain + one increment per
  prime.
- Work is split into fixed 10^9-wide chunks processed by a thread pool, but
  results are **merged strictly in order**, so the running prime count is
  exact. When a merge would cross a power of 10, that one chunk is re-scanned
  serially to pin down the counts at exactly the 10^n-th prime (and the
  re-scan totals are cross-checked against the parallel result).
- Every crossing prints the counts with an `OK` / `MISMATCH` verdict against
  the embedded known OEIS values (n ≤ 12), so the run self-validates on the
  way to the new term.
- Progress checkpoints to `digitprimes.ckpt` every minute; SIGINT/SIGTERM
  shut down cleanly and a rerun resumes where it left off.

## Usage

```
brew install primesieve
make
./digitprimes --target 13          # the full a(13) run (resumable)
```

Options: `--target N` (default 13), `--threads T`, `--chunk C`,
`--checkpoint FILE`, `--out FILE`, `--fresh` (ignore existing checkpoint).

Results are written to `digitprimes_results.txt` (CSV of all crossings plus
OEIS-formatted term lists) and updated at every crossing, so partial results
survive interruption.

On an Apple M4 Max (14 threads) the tool verifies all known terms through
a(10) in ~6 s and runs at roughly 1.7–2 billion primes/s; the full a(13) pass
takes on the order of 2 hours.
