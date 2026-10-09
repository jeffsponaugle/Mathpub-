# The y-cruncher `.ycd` compressed digit format, and the `ycd` C library

`.ycd` is the file format [y-cruncher](http://www.numberworld.org/y-cruncher/)
uses to store the digits of a constant compactly: 19 decimal digits in every
8 bytes (2.37 digits per byte, against 1 for text). Large results are split
into numbered block files. The 100-trillion-digit pi computation published by
Google is distributed this way, as 1,000 files of 100 billion digits each.

This document describes the format as established by reading real y-cruncher
output, and the small C library in this directory (`ycd.h`, `ycd.c`) that reads
and writes it. Each statement about the format is marked in the last section
as verified against real files or not.

## 1. File layout

```
+---------------------------+
| text header               |   ASCII, CRLF line ends, ~200 bytes
| NUL byte (0x00)           |   ends the header
+---------------------------+
| word 0                    |   8 bytes each, little-endian unsigned 64-bit
| word 1                    |
| ...                       |
| last word of the block    |   may hold fewer valid digits (see 1.3)
+---------------------------+
| trailing bytes            |   NOT digits of this block; ignore (see 1.4)
+---------------------------+
```

### 1.1 Header

The header of a real file, byte for byte (`\t` is a tab, every line ends in
CR LF, and blank lines are part of the layout):

```
#Compressed Digit File

FileVersion:	1.1.0

Base:	10

FirstDigits:	3.14159265358979323846264338327950288419716939937510

TotalDigits:	0

Blocksize:	100000000000
BlockID:	1

EndHeader

<NUL>
```

| field | meaning |
|---|---|
| `#Compressed Digit File` | magic; the first bytes of every file |
| `FileVersion` | `1.0.0` or `1.1.0`. The data layout is the same in both. |
| `Base` | `10` or `16` |
| `FirstDigits` | the constant's first 50 digits as text, integer part and radix point included. Informational, but useful as a check (section 2.2). |
| `TotalDigits` | digits after the radix point in the whole set, identical in every file of the set. `0` means not recorded. |
| `Blocksize` | digits per block file |
| `BlockID` | index of this block, from 0 |
| `EndHeader` | last header line, followed by CR LF CR LF and one NUL byte |

Parse it as `Name:<tab>Value` lines. Skip blank lines and unknown names, stop
at `EndHeader`, and take the data to begin at the byte after the first NUL in
the file. Do not hard-code the header length. It was 200 bytes in the file
above, and it changes with the number of digits in every value.

### 1.2 Words

Each word is an unsigned 64-bit integer stored little-endian. In base 10 it
holds 19 decimal digits, the most significant digit first in stream order:

```
bytes (LE)  f5 fa c7 b2 0a 18 b8 4c   =  5528194976825473781
digits      5 5 2 8 1 9 4 9 7 6 8 2 5 4 7 3 7 8 1
```

So digit `i` of a word (0 = first) is `(w / 10^(18-i)) % 10`, and every valid
word is below 10^19. A word at or above 10^19 means corruption or a wrong file
offset, which makes a cheap integrity check. In base 16 a word holds 16 hex
digits, most significant nibble first.

Words start at the beginning of each block. Block files are self-contained:
word 0 of block `b` begins with digit number `b * Blocksize` of the stream,
whatever `Blocksize mod 19` is.

### 1.3 The last word of a block

A block holds `Blocksize` digits, or fewer in the final block of a set where
`TotalDigits` cuts it short. That takes `ceil(n / 19)` words. When `n` is not
a multiple of 19 the last word is partial. Its valid digits are the leading
ones, and the remaining low-order digits of that word are not part of the
block. For a 100,000,000,000-digit block that is 5,263,157,895 words, the last
carrying 14 valid digits.

### 1.4 Trailing bytes and file size

Files are longer than their data. y-cruncher sizes them to a multiple of 4,096
bytes, and the real 100-billion-digit block examined here carries 352 more
words after its last block word. Those words are not zero and all decode as
valid 19-digit values, so a reader that trusts the file size will deliver
digits that do not belong to the block.

**The number of digits in a file comes from the header, never from its size:**

```
start   = BlockID * Blocksize
end     = start + Blocksize            (capped at TotalDigits when that is non-zero)
ndigits = end - start
```

The file size is only good for detecting truncation, when fewer than
`ceil(ndigits / 19)` words are present.

One consequence matters for anyone writing files. If `TotalDigits` is 0 and the
final block is short, its length cannot be recovered at all, because the
padding is indistinguishable from digits. A set with `TotalDigits: 0` must end
on a block boundary, as the 1,000-block Google set does. Otherwise record
`TotalDigits`. The library's writer refuses to produce the ambiguous case.

## 2. Positions

### 2.1 Offset 0 is the first digit after the point

The digit stream does not include the integer part. Word 0 of block 0 of pi is
`1415926535897932384`, not `3141592653589793238`. With positions counted from 1
at the first digit after the radix point, which is the usual convention for
digit positions in pi and the one the OEIS uses:

```
position p  ->  stream offset o = p - 1
            ->  block  b = o / Blocksize
            ->  word   k = (o % Blocksize) / 19,  digit (o % Blocksize) % 19 within it
byte offset of that word in block b's file = data_offset + 8 * k
```

Block `b` therefore covers positions `b * Blocksize + 1` through
`(b + 1) * Blocksize`. Random access costs one `pread` of a few words.

### 2.2 Checking the convention

When block 0 is present the library compares its first word with the
`FirstDigits` header. If the word equals the first 19 digits after the point,
the layout above is confirmed. If it equals the integer part followed by the
fractional digits, the stream carries the integer part, as files from some
other encoder might, and the library shifts every position so that position 1
is still the first fractional digit. When block 0 is absent nothing can be
checked, and the standard layout is assumed.

## 3. Sets of files

- A set is one or more files with the same `Base` and `Blocksize` and
  contiguous `BlockID`s. Order them by `BlockID`. File names carry no meaning
  to a reader. y-cruncher names them `<Constant> - <Dec|Hex> - <Algorithm> - <BlockID>.ycd`,
  which sorts wrongly as text (`10` before `2`).
- A set need not start at block 0. Positions stay absolute.
- Every block except the last must be full.

## 4. The library

Two files, C99 and POSIX, no other dependencies: [ycd.h](ycd.h) and
[ycd.c](ycd.c). Add them to a project and compile.

```c
#include "ycd.h"

char err[256];
ycd_set *s = ycd_open("/data/pi_ycd", err, sizeof err);   /* file, directory, or a,b,c list */
if (!s) { fprintf(stderr, "%s\n", err); return 1; }
printf("positions %llu..%llu\n", (unsigned long long)ycd_first_pos(s), (unsigned long long)ycd_last_pos(s));

ycd_reader *r = ycd_reader_new(s, YCD_ASCII);             /* 0 = digit values 0..9 */
uint8_t buf[101];
int64_t n = ycd_reader_read_at(r, 1000000, buf, 100);     /* 100 digits from position 1,000,000 */
if (n < 0) fprintf(stderr, "%s\n", ycd_reader_error(r));
buf[n] = 0; puts((char *)buf);

while ((n = ycd_reader_read(r, buf, 100)) > 0) { /* ...continues sequentially... */ }

ycd_reader_free(r);
ycd_close(s);
```

| call | purpose |
|---|---|
| `ycd_open(spec, err, n)` / `ycd_open_paths` | open a set from files and directories (a directory contributes every `*.ycd` in it); validates base, block size, contiguity, truncation |
| `ycd_first_pos`, `ycd_last_pos`, `ycd_ndigits`, `ycd_base`, `ycd_nfiles`, `ycd_get_file` | what the set covers, and each file's parsed header (`ycd_file`) |
| `ycd_integer_digits`, `ycd_convention_checked`, `ycd_length_uncertain` | results of the checks in 2.2 and 1.4 |
| `ycd_reader_new(set, flags)` | a cursor; `YCD_ASCII` delivers text instead of values |
| `ycd_reader_seek`, `ycd_reader_tell`, `ycd_reader_read`, `ycd_reader_read_at` | sequential and random access across block boundaries; a short count only at the end of the set |
| `ycd_reader_set_buffer` | raw buffer size: about 4 to 64 KiB for random access, several MiB for streaming [1 MiB] |
| `ycd_writer_create`, `ycd_writer_append`, `ycd_writer_close` | write one block file, with the header laid out exactly as y-cruncher's |
| `ycd_decode_dec`, `ycd_encode_dec`, `ycd_decode_hex`, `ycd_encode_hex`, `ycd_load_le64`, `ycd_store_le64` | the inline word codec, usable on its own buffers |
| `ycd_is_ycd`, `ycd_read_header` | classify and parse a single file |

Behaviour worth knowing:

- **Errors.** Functions return -1 (or NULL) with a message, and never exit or
  print. Open errors go to the caller's buffer, read errors to
  `ycd_reader_error`.
- **Threads.** A `ycd_set` is read-only after opening and can be shared. Give
  each thread its own `ycd_reader`. Readers use `pread`, so they never disturb
  each other. For millions of random lookups on a disk array, run many readers
  in parallel with small buffers: in this project that change took a
  6.7-million-position lookup against a 21 TB source from 25 hours to 80
  minutes.
- **Integrity.** Every decimal word is checked to be below 10^19 as it is
  decoded, and a read fails with the file name and word index otherwise.
- **Speed.** The decoder splits a word into 8 + 8 + 3 digits with two constant
  divisions and a two-digit lookup table rather than 19 divisions in series.
  One thread streams about 1.4 billion digits/s from cache on an Apple M-series
  core, and 5 billion digits of the real block in 2.0 s including text output.
  Across threads, split the position range and give each its own reader.
- **Portability.** Byte order is handled explicitly, so big-endian hosts work.
  POSIX calls (`open`, `pread`, `opendir`) are used, so Windows needs a shim.

### ycdtool

[ycdtool.c](ycdtool.c) is a command-line front end and a worked example of the
API.

```
ycdtool info   SRC                       headers, block table, position range, convention check
ycdtool dump   SRC POS [N]               N digits [50] from position POS
ycdtool totext SRC [-s POS] [-n COUNT]   digits as text on stdout
ycdtool check  SRC                       read everything, validate every word, digit frequencies, speed
ycdtool encode IN.txt OUTDIR -b BLOCKSIZE [-n COUNT] [-N NAME] [-F FIRSTDIGITS]
```

`SRC` is a file, a directory, or a comma-separated list. Counts accept K, M, B
and T suffixes. `encode` takes text digits, splits a leading `3.` off into
`FirstDigits`, counts the digits first so that `TotalDigits` is recorded, and
writes `OUTDIR/NAME - <id>.ycd`.

### Building and testing

```
make ycdtool ycd_test
./ycd_test
```

Or directly: `cc -O2 -o ycdtool ycdtool.c ycd.c`. The self-test covers the
codec against a reference for 2 million words, write and read round trips in
both bases with block sizes that are and are not multiples of the word length
(including blocks shorter than one word), 4,000 random reads per set aimed at
block seams, files named out of block order, sub-sets without block 0,
integer-part detection, and the error paths.

## 5. What is verified, and how

| statement | status |
|---|---|
| Header layout and field names (1.1) | verified on a real y-cruncher file (`FileVersion 1.1.0`, block 1 of a 100-billion-digit-per-block pi set); the writer reproduces it byte for byte |
| Little-endian 64-bit words of 19 digits, most significant first (1.2) | verified: block 1's first word decodes to `5528194976825473781`, the digits of pi at positions 100,000,000,001 to 100,000,000,019 |
| Blocks start word-aligned at position `BlockID * Blocksize + 1`, and offset 0 is the first digit after the point (1.2, 2.1) | verified by the same check |
| Digits are bounded by the header, and the trailing bytes are not padding zeros (1.4) | verified: 352 non-zero, valid-looking words follow the block's last word in the real file |
| Files sized to a multiple of 4,096 bytes (1.4) | observed on the one real file; the writer does the same |
| Partial last word carries its valid digits first (1.3) | consistent with the most-significant-first layout and is what the library reads and writes, but not yet confirmed at a real block boundary. To confirm, compare the last 14 digits of a real block with a text copy of the same positions: `ycdtool dump DIR 199999999987 14` |
| `FileVersion 1.0.0` | accepted on the assumption that the layout is unchanged; no 1.0.0 file has been examined |
| `Base 16`, 16 hex digits per word, most significant first | implemented from y-cruncher's documented behaviour and covered by round-trip tests, but not checked against a real hex file |
| A set with `TotalDigits: 0` ends on a block boundary (1.4) | a consequence of the format rather than a documented rule; holds for the real set examined |

Cross-checks between independent implementations: files written by this
library are read correctly by the pi tools in this directory (`pifind` finds
known digit strings at their known positions, including across block seams),
and the whole first billion digits of pi round-trip through `ycdtool encode`
and `ycdtool totext` with identical checksums.
