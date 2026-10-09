/*
 * ycd.h - read and write y-cruncher compressed digit files (.ycd)
 *
 * A small, dependency-free C99/POSIX library.  Add ycd.h and ycd.c to a
 * project; nothing else is needed.  See ycd.md for the file format.
 *
 *   ycd_set *s = ycd_open("/data/pi_ycd", err, sizeof err);   // file(s) and/or directories
 *   ycd_reader *r = ycd_reader_new(s, 0);
 *   ycd_reader_seek(r, 1000000);                  // position 1 = first digit after the point
 *   uint8_t d[50];
 *   int64_t n = ycd_reader_read(r, d, 50);        // digit values 0..9 (or text with YCD_ASCII)
 *
 * Conventions
 *   - Positions are 1-based and count digits after the radix point: for pi,
 *     position 1 is the "1" of 3.14159...  (ycd_first_pos/ycd_last_pos give
 *     the range a set covers; a set need not start at block 0.)
 *   - A "set" is one or more .ycd files of the same constant: they are
 *     ordered by their BlockID header field (file names are irrelevant), and
 *     must have equal Blocksize and contiguous BlockIDs.
 *   - A ycd_set is immutable after ycd_open and may be shared by threads;
 *     each thread uses its own ycd_reader (readers hold the file descriptor
 *     and buffers).
 *   - Functions return 0 (or a count) on success and -1 on error; the message
 *     is in the err buffer passed to ycd_open / ycd_read_header, or in
 *     ycd_reader_error() / ycd_writer_error().
 *
 * This code is released into the public domain (or under the MIT license,
 * at your option).
 */
#ifndef YCD_H
#define YCD_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define YCD_VERSION "1.0.0"
#define YCD_MAGIC "#Compressed Digit File"
#define YCD_DEC_DIGITS_PER_WORD 19      /* Base 10: 19 decimal digits per 64-bit word */
#define YCD_HEX_DIGITS_PER_WORD 16      /* Base 16: 16 hex digits per 64-bit word     */
#define YCD_MAX_DIGITS_PER_WORD 19
#define YCD_FILE_ALIGN 4096             /* y-cruncher sizes files to a multiple of this */

/* reader/writer flags */
#define YCD_ASCII 1                     /* digits as text ('0'..'9', 'a'..'f') instead of values 0..15 */

/* ---------------------------------------------------------------- */
/* One file                                                          */
/* ---------------------------------------------------------------- */
typedef struct {
    char *path;
    char version[16];           /* FileVersion, "1.0.0" or "1.1.0"                   */
    int base;                   /* 10 or 16                                          */
    int digits_per_word;        /* 19 or 16                                          */
    char first_digits[128];     /* FirstDigits header, e.g. "3.14159265358979..."    */
    uint64_t total_digits;      /* TotalDigits header (0 = not recorded)             */
    uint64_t blocksize;         /* digits per block                                  */
    uint64_t blockid;           /* index of this block, 0-based                      */
    uint64_t data_offset;       /* byte offset of the first data word                */
    uint64_t file_size;         /* bytes                                             */
    uint64_t start;             /* 0-based offset of this file's first digit in the  */
                                /* whole digit stream = blockid * blocksize          */
    uint64_t ndigits;           /* digits this file holds (from the header fields,   */
                                /* never from the file size: files are over-long)    */
    int length_uncertain;       /* 1: the file is shorter than Blocksize and         */
                                /* TotalDigits is 0, so ndigits is only an upper     */
                                /* bound (whole words present, padding included)     */
} ycd_file;

/* ycd_is_ycd - 1 if the file starts with the ycd magic line, 0 if not, -1 if unreadable */
int ycd_is_ycd(const char *path);

/* ycd_read_header - parse one file's header into *f (f->path is strdup'd; free with ycd_file_free) */
int ycd_read_header(const char *path, ycd_file *f, char *err, size_t errlen);
void ycd_file_free(ycd_file *f);

/* ---------------------------------------------------------------- */
/* A set of files (one constant, contiguous blocks)                  */
/* ---------------------------------------------------------------- */
typedef struct ycd_set ycd_set;

/*
 * ycd_open - open a set from a comma-separated list of .ycd files and/or
 * directories (a directory contributes every *.ycd directly inside it).
 * ycd_open_paths takes the list as an array instead (paths may then contain
 * commas).  Returns NULL and fills err on failure.
 */
ycd_set *ycd_open(const char *spec, char *err, size_t errlen);
ycd_set *ycd_open_paths(const char *const *paths, int npaths, char *err, size_t errlen);
void ycd_close(ycd_set *s);

int ycd_nfiles(const ycd_set *s);
const ycd_file *ycd_get_file(const ycd_set *s, int i);     /* ordered by BlockID */
int ycd_base(const ycd_set *s);                            /* 10 or 16 */
uint64_t ycd_first_pos(const ycd_set *s);                  /* first position available (1 when block 0 is present) */
uint64_t ycd_last_pos(const ycd_set *s);                   /* last position available */
uint64_t ycd_ndigits(const ycd_set *s);                    /* ycd_last_pos - ycd_first_pos + 1 */

/*
 * ycd_integer_digits - how many integer-part digits the stream carries before
 * the first fractional digit.  y-cruncher's files carry none (0): word 0 of
 * block 0 starts with the first digit after the point.  The library checks
 * block 0's first word against the FirstDigits header and also recognises a
 * stream that includes the integer part (e.g. files from another encoder);
 * positions are adjusted so that position 1 is always the first fractional
 * digit.  ycd_convention_checked is 1 when block 0 was present and its first
 * word matched one of the two layouts, 0 when it could not be checked (no
 * block 0, or no FirstDigits): then 0 integer digits are assumed.
 */
int ycd_integer_digits(const ycd_set *s);
int ycd_convention_checked(const ycd_set *s);

/* ycd_length_uncertain - 1 if the last block's digit count is only an upper bound (see ycd_file) */
int ycd_length_uncertain(const ycd_set *s);

/* ---------------------------------------------------------------- */
/* Reading                                                           */
/* ---------------------------------------------------------------- */
typedef struct ycd_reader ycd_reader;

/*
 * ycd_reader_new - a sequential/random-access cursor over a set, positioned
 * at ycd_first_pos.  flags: 0 or YCD_ASCII.  The buffer size can be changed
 * before the first read (default 1 MiB of words; use 64 KiB or less for
 * random access, several MiB for streaming).
 */
ycd_reader *ycd_reader_new(const ycd_set *s, int flags);
void ycd_reader_free(ycd_reader *r);
int ycd_reader_set_buffer(ycd_reader *r, size_t bytes);

/* ycd_reader_seek - move to a position (first_pos .. last_pos + 1); -1 if outside */
int ycd_reader_seek(ycd_reader *r, uint64_t pos);
uint64_t ycd_reader_tell(const ycd_reader *r);             /* position of the next digit */

/*
 * ycd_reader_read - read up to max digits at the cursor, crossing file
 * boundaries; returns the number delivered (0 at the end of the set, fewer
 * than max only there) or -1 on error.
 */
int64_t ycd_reader_read(ycd_reader *r, uint8_t *dst, size_t max);

/* ycd_reader_read_at - seek + read in one call */
int64_t ycd_reader_read_at(ycd_reader *r, uint64_t pos, uint8_t *dst, size_t n);

const char *ycd_reader_error(const ycd_reader *r);

/* ---------------------------------------------------------------- */
/* Writing                                                           */
/* ---------------------------------------------------------------- */
typedef struct {
    int base;                   /* 10 (default when 0) or 16                         */
    uint64_t blocksize;         /* digits per block (required)                       */
    uint64_t blockid;           /* this file's block                                 */
    uint64_t total_digits;      /* TotalDigits header: digits in the WHOLE set, the  */
                                /* same in every file.  0 = not recorded; then the   */
                                /* set must end on a block boundary, because a short */
                                /* final block's length cannot be recovered          */
    const char *first_digits;   /* FirstDigits header, e.g. "3.1415926535..." (may be NULL) */
    const char *version;        /* FileVersion [1.1.0]                               */
} ycd_write_opts;

typedef struct ycd_writer ycd_writer;

/*
 * ycd_writer_create - start one block file; append its digits in any number
 * of calls (values 0..15, or text when flags has YCD_ASCII; at most
 * blocksize digits in total), then close.  close writes the final partial
 * word left-aligned and pads the file with zero bytes to a multiple of 4096.
 */
ycd_writer *ycd_writer_create(const char *path, const ycd_write_opts *opts, int flags, char *err, size_t errlen);
int ycd_writer_append(ycd_writer *w, const uint8_t *digits, size_t n);
/*
 * ycd_writer_close - finish the file and free w.  Returns -1 on a write
 * error, and also when fewer than blocksize digits were appended to a file
 * whose opts.total_digits was 0: such a file is ambiguous (its digits cannot
 * be told from the padding), so always set total_digits for a set that does
 * not end on a block boundary.
 */
int ycd_writer_close(ycd_writer *w);
const char *ycd_writer_error(const ycd_writer *w);

/* ---------------------------------------------------------------- */
/* Word codec (inline; usable on its own)                            */
/* ---------------------------------------------------------------- */

/* ycd_load_le64 / ycd_store_le64 - unaligned little-endian 64-bit access */
static inline uint64_t ycd_load_le64(const uint8_t *p)
{
    uint64_t w = 0;
    for (int i = 7; i >= 0; i--) w = (w << 8) | p[i];
    return w;
}
static inline void ycd_store_le64(uint8_t *p, uint64_t w)
{
    for (int i = 0; i < 8; i++) { p[i] = (uint8_t)w; w >>= 8; }
}

#define YCD_DEC_WORD_LIMIT 10000000000000000000ULL          /* 10^19: a decimal word is below this */

extern const uint8_t ycd_pairs[200];                        /* digits of 00..99 as values */

static inline void ycd_put8_(uint32_t v, uint8_t *out)      /* 8 digits of v < 10^8 */
{
    uint32_t a = v / 10000u, b = v % 10000u;
    uint32_t a1 = a / 100u, a2 = a % 100u, b1 = b / 100u, b2 = b % 100u;
    out[0] = ycd_pairs[2 * a1]; out[1] = ycd_pairs[2 * a1 + 1];
    out[2] = ycd_pairs[2 * a2]; out[3] = ycd_pairs[2 * a2 + 1];
    out[4] = ycd_pairs[2 * b1]; out[5] = ycd_pairs[2 * b1 + 1];
    out[6] = ycd_pairs[2 * b2]; out[7] = ycd_pairs[2 * b2 + 1];
}

/*
 * ycd_decode_dec - the 19 decimal digits of a word (< 10^19), most
 * significant first, as values 0..9.  Two constant divisions split the word
 * into 8 + 8 + 3 digits, finished with 32-bit arithmetic and a pair table:
 * several times faster than a chain of 19 divisions.
 */
static inline void ycd_decode_dec(uint64_t w, uint8_t *out)
{
    uint64_t top = w / 100000000000ULL;
    uint64_t rest = w - top * 100000000000ULL;
    uint32_t mid = (uint32_t)(rest / 1000u);
    uint32_t last = (uint32_t)(rest - (uint64_t)mid * 1000u);
    ycd_put8_((uint32_t)top, out);
    ycd_put8_(mid, out + 8);
    uint32_t l2 = last % 100u;
    out[16] = (uint8_t)(last / 100u);
    out[17] = ycd_pairs[2 * l2]; out[18] = ycd_pairs[2 * l2 + 1];
}

/* ycd_encode_dec - pack 19 digit values (most significant first) into a word */
static inline uint64_t ycd_encode_dec(const uint8_t *digits)
{
    uint64_t w = 0;
    for (int i = 0; i < YCD_DEC_DIGITS_PER_WORD; i++) w = w * 10 + digits[i];
    return w;
}

/* ycd_decode_hex / ycd_encode_hex - 16 hex digits per word, most significant first */
static inline void ycd_decode_hex(uint64_t w, uint8_t *out)
{
    for (int i = 0; i < YCD_HEX_DIGITS_PER_WORD; i++) out[i] = (uint8_t)((w >> (60 - 4 * i)) & 15u);
}
static inline uint64_t ycd_encode_hex(const uint8_t *digits)
{
    uint64_t w = 0;
    for (int i = 0; i < YCD_HEX_DIGITS_PER_WORD; i++) w = (w << 4) | (digits[i] & 15u);
    return w;
}

#ifdef __cplusplus
}
#endif
#endif /* YCD_H */
