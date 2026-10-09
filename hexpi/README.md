# hexpi — hex-digit frequencies of Pi (OEIS A099333–A099348)

`hexcount` streams a file of hexadecimal digits and, at each 10^n digit
boundary, prints the cumulative count of each of the 16 hex digits. Row
`10^n`, column `d` is the term a(n) of OEIS A099333 (d=0) through
A099348 (d=f): "Frequency of the hexadecimal d in the first 10^n
hexadecimal digits of Pi."

Conventions match the OEIS sequences: the integer part is excluded (a
leading `3.` is auto-detected and skipped), counting starts at the first
fractional digit (`243f6a88...`). Whitespace (newlines, spaces, tabs) is
ignored anywhere; any other byte aborts with its file offset. Upper- and
lowercase hex both work.

## Build

```
cc -O2 -o hexcount hexcount.c
```

## Usage

```
./hexcount [-b MB] [--csv FILE] [--check-pi] [--no-skip-prefix] FILE
```

- `FILE` — the hex-digit file (`-` reads stdin, so you can pipe from
  `zstdcat`/`gzcat`).
- `--csv FILE` — also write rows as CSV (`digits,0,1,...,f`), flushed as
  each boundary is reached.
- `--check-pi` — compare each row against the published OEIS terms
  (n = 1..12) and tag it `OEIS:OK` / `OEIS:MISMATCH`; exits 2 on any
  mismatch.
- `--no-skip-prefix` — count a leading `3.` literally instead of
  skipping it.
- `-b MB` — read buffer size (default 16).

A `final` row with the total is printed when the file length is not an
exact power of ten. Progress (digits, GB, MB/s) goes to stderr every 5 s
when stderr is a terminal. Throughput is ~1.8 GB/s, so a 4-trillion-digit
file is disk-bound (~35 min of CPU).

## Test data

`gen_pihex.py N outfile [--plain]` generates the first N fractional hex
digits of Pi (Machin's formula, pure Python, ~7 s for N=100000). Default
format has the `3.` prefix and 64-char lines; `--plain` is the bare
digit stream.

Verified: `python3 gen_pihex.py 100000 pihex_1e5.txt` then
`./hexcount --check-pi pihex_1e5.txt` reproduces all 16 sequences
exactly for n = 1..5.
