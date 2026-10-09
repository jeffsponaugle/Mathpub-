/*
 * pi_common.h - shared helpers for pibloom / pisearch (A197123 tools)
 *
 * Conventions (match OEIS A197123 and the original tools):
 *   - Digit positions are 1-based, counting from the first digit AFTER the
 *     decimal point.  "1 appears at positions 1 and 3", "26 at 6 and 21".
 *   - The leading "3." in the source file is skipped and never part of a window.
 *   - A "window" of length L at position p is the L digits starting at p.
 */
#ifndef PI_COMMON_H
#define PI_COMMON_H

#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>              /* printf/fprintf/snprintf, FILE for candidate & result files */
#include <stdlib.h>             /* malloc/free/realloc, strtod/strtoull, qsort, exit */
#include <stdint.h>             /* uint64_t/uint8_t fixed-width integers */
#include <string.h>             /* memcpy/memmove/memset/memchr/strlen/strsep */
#include <errno.h>              /* errno + strerror for error messages */
#include <fcntl.h>              /* open() flags, fcntl(F_RDAHEAD) on macOS */
#include <unistd.h>             /* read/lseek/close/pread, sysconf, isatty, getopt */
#include <time.h>               /* clock_gettime for timing/progress */
#include <signal.h>             /* sigaction for clean Ctrl-C handling */
#include <math.h>               /* exp/pow for bloom fill and false-positive estimates */
#include <inttypes.h>           /* PRIu64 printf formats for 64-bit integers */
#include <sys/mman.h>           /* mmap/munmap/madvise for the huge bloom and table allocations */
#include <sys/stat.h>           /* stat/fstat to size input files */
#include <sys/resource.h>       /* getrusage for peak RSS reporting */
#include <stdarg.h>             /* va_list for the die()/warn() printf-style helpers */
#include <pthread.h>            /* parallel pre-faulting of huge allocations */
#include <stdatomic.h>          /* progress counter shared by the pre-fault threads */

/* Suite version: bump on EVERY change to any tool or this header (see the
 * CHANGELOG in README.md).  Printed by each tool's banner, -h and -V, and
 * stamped into candidate/result file headers. */
#define PI_TOOLS_VERSION "1.6.0"

#define MAX_SEQ_LEN 38          /* two uint64 decimal halves: 19 + 19 digits */
#define LO_DIGITS   19          /* digits held in the 'lo' word              */

/* ------------------------------------------------------------------ */
/* Error / logging                                                     */
/* ------------------------------------------------------------------ */
/*
 * die - print a formatted FATAL message to stderr and exit(2).
 * Used for every unrecoverable condition (bad arguments, I/O errors,
 * allocation failures, corrupt input) so that a run never continues on
 * bad state.  Flushes stdout first so partial progress output is not lost.
 */
static void die(const char *fmt, ...) __attribute__((noreturn, format(printf, 1, 2)));
static void die(const char *fmt, ...)
{
    va_list ap;
    fflush(stdout);
    fprintf(stderr, "\nFATAL: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    exit(2);
}
/*
 * warn - print a formatted WARNING to stderr and continue.
 * For conditions the run can survive but the user should know about
 * (e.g. a too-small bloom filter, an unusual file header).
 */
static void warn(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void warn(const char *fmt, ...)
{
    va_list ap;
    fflush(stdout);
    fprintf(stderr, "WARNING: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
}

/*
 * print_version - "<tool> v<version> (built <date>)" on stdout; used by -V.
 */
static void print_version(const char *tool)
{
    printf("%s v%s (built %s %s)\n", tool, PI_TOOLS_VERSION, __DATE__, __TIME__);
}

/* ------------------------------------------------------------------ */
/* Time                                                                */
/* ------------------------------------------------------------------ */
/*
 * now_sec - monotonic wall-clock time in seconds (CLOCK_MONOTONIC).
 * Used for all rate, ETA and progress-interval calculations; immune to
 * system clock changes during multi-day runs.
 */
static inline double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/*
 * fmt_duration - format seconds as "m:ss", "h:mm:ss" or "Nd hh:mm:ss"
 * into the caller-supplied buffer; returns buf for use inside printf.
 */
static const char *fmt_duration(double s, char *buf, size_t n)
{
    if (s < 0) s = 0;
    long t = (long)s;
    long d = t / 86400, h = (t / 3600) % 24, m = (t / 60) % 60, sec = t % 60;
    if (d > 0)      snprintf(buf, n, "%ldd %02ld:%02ld:%02ld", d, h, m, sec);
    else if (h > 0) snprintf(buf, n, "%ld:%02ld:%02ld", h, m, sec);
    else            snprintf(buf, n, "%ld:%02ld", m, sec);
    return buf;
}

/*
 * fmt_u64 - format an integer with thousands separators ("1,234,567")
 * into the caller-supplied buffer; returns buf.  Callers must use a
 * distinct buffer per value within one printf call.
 */
static const char *fmt_u64(uint64_t v, char *buf, size_t n)
{
    char tmp[32];
    int len = snprintf(tmp, sizeof tmp, "%" PRIu64, v);
    int commas = (len - 1) / 3;
    int out = len + commas;
    if ((size_t)out + 1 > n) { snprintf(buf, n, "%" PRIu64, v); return buf; }
    buf[out] = 0;
    int j = out - 1;
    for (int i = len - 1, c = 0; i >= 0; i--) {
        buf[j--] = tmp[i];
        if (++c == 3 && i > 0) { buf[j--] = ','; c = 0; }
    }
    return buf;
}

/*
 * fmt_bytes - format a byte count as a human-readable binary size
 * ("700.00 GiB"); returns buf.
 */
static const char *fmt_bytes(double b, char *buf, size_t n)
{
    const char *u[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
    int i = 0;
    while (b >= 1024.0 && i < 5) { b /= 1024.0; i++; }
    snprintf(buf, n, "%.2f %s", b, u[i]);
    return buf;
}

/* ------------------------------------------------------------------ */
/* Size parsing                                                        */
/* ------------------------------------------------------------------ */
/*
 * parse_mem_size - parse a memory size with optional K/M/G/T suffix
 * (1024-based; "MB"/"GiB" spellings accepted; fractions allowed, e.g.
 * "1.5T").  Returns 0 and stores the byte count, or -1 on a bad string.
 */
static int parse_mem_size(const char *s, uint64_t *out)
{
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (end == s || errno != 0 || v < 0) return -1;
    double mult = 1;
    if (*end) {
        switch (*end) {
        case 'k': case 'K': mult = 1024.0; break;
        case 'm': case 'M': mult = 1024.0 * 1024; break;
        case 'g': case 'G': mult = 1024.0 * 1024 * 1024; break;
        case 't': case 'T': mult = 1024.0 * 1024 * 1024 * 1024; break;
        default: return -1;
        }
        end++;
        if (*end == 'B' || *end == 'b') end++;     /* allow "MB", "GiB" */
        if (*end == 'i') end++;
        if (*end == 'B' || *end == 'b') end++;
        if (*end) return -1;
    }
    double r = v * mult;
    if (r > 1.8e19) return -1;
    *out = (uint64_t)r;
    return 0;
}

/*
 * parse_count - parse a digit count with optional K/M/B(G)/T suffix
 * (1000-based), or the word "all" (= UINT64_MAX).  Returns 0 on success,
 * -1 on a bad string.
 */
static int parse_count(const char *s, uint64_t *out)
{
    if (!strcmp(s, "all") || !strcmp(s, "ALL")) { *out = UINT64_MAX; return 0; }
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (end == s || errno != 0 || v < 0) return -1;
    double mult = 1;
    if (*end) {
        switch (*end) {
        case 'k': case 'K': mult = 1e3; break;
        case 'm': case 'M': mult = 1e6; break;
        case 'b': case 'B': case 'g': case 'G': mult = 1e9; break;
        case 't': case 'T': mult = 1e12; break;
        default: return -1;
        }
        end++;
        if (*end) return -1;
    }
    double r = v * mult;
    if (r > 1.8e19) return -1;
    *out = (uint64_t)r;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Hashing                                                             */
/* ------------------------------------------------------------------ */
/*
 * mix64 - splitmix64 finalizer: a fast 64-bit avalanche mixer where every
 * output bit depends on every input bit.  This is what makes a fixed
 * leading-digit prefix (the -f filter) still spread uniformly over the
 * whole bloom table.
 */
static inline uint64_t mix64(uint64_t x)
{
    /* splitmix64 finalizer */
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27; x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
}

/*
 * hash_window - derive two independent 64-bit hashes (h1, h2) from a
 * window key (hi, lo).  The k bloom positions are then h1 + i*h2 (Kirsch-
 * Mitzenmacher double hashing); h2 is forced odd so the stride never
 * degenerates.  pisearch uses h1 for its hash table and prefilter.
 */
static inline void hash_window(uint64_t lo, uint64_t hi, uint64_t *h1, uint64_t *h2)
{
    uint64_t a = mix64(lo + 0x9E3779B97F4A7C15ULL);
    uint64_t b = mix64(hi ^ a ^ 0xD1B54A32D192ED03ULL);
    *h1 = a ^ (b << 1);
    *h2 = mix64(b + a) | 1;          /* odd: good stride for double hashing */
}

/*
 * fastrange64 - map a uniform 64-bit hash into [0, range) with a 128-bit
 * multiply instead of a modulo (Lemire).  Lets the bloom table be ANY size
 * in bytes, not just a power of two, with no division in the hot loop.
 */
static inline uint64_t fastrange64(uint64_t h, uint64_t range)
{
    return (uint64_t)(((__uint128_t)h * (__uint128_t)range) >> 64);
}

/* ------------------------------------------------------------------ */
/* Prefix-range filter: windows whose first `len` digits, read as a      */
/* number, lie in [lo, hi].  Specs: a value or inclusive range of 1..3   */
/* digits ("3", "24", "0-4", "10-15", "05-09", "334-666"), or K/N = the  */
/* K-th (1-based) of N equal partitions, resolved to the coarsest prefix */
/* length that divides N exactly (1/10 -> [1], 2/20 -> [05-09]) or to 3  */
/* digits otherwise (2/3 -> [334-666], equal to within 0.1%).            */
/* ------------------------------------------------------------------ */
#define PFILTER_MAX_LEN 3

typedef struct {
    int len;            /* 0 = no filter, else 1..3 digits */
    int lo, hi;         /* inclusive prefix range          */
    int frac_k, frac_n; /* set when given as K/N, for display */
    char str[16];       /* canonical text, e.g. "10-15"    */
} pfilter;

/*
 * pfilter_set - fill a filter from an explicit (len, lo, hi) and build its
 * canonical display string.
 */
static void pfilter_set(pfilter *f, int len, int lo, int hi)
{
    if (len < 1 || len > PFILTER_MAX_LEN || lo < 0 || hi < lo || hi >= 1000)
        die("internal: bad prefix filter (len %d, %d-%d)", len, lo, hi);
    /* explicitly bounded copies so the formatted length is provably < 16 */
    int w = len < 1 ? 1 : len > PFILTER_MAX_LEN ? PFILTER_MAX_LEN : len;
    unsigned ulo = (unsigned)lo % 1000u, uhi = (unsigned)hi % 1000u;
    f->len = len; f->lo = lo; f->hi = hi;
    if (lo == hi) snprintf(f->str, sizeof f->str, "%0*u", w, ulo);
    else          snprintf(f->str, sizeof f->str, "%0*u-%0*u", w, ulo, w, uhi);
}

/*
 * pfilter_parse - parse a -f prefix spec: a 1-3 digit value ("3", "24"),
 * an inclusive range with equal-length ends ("0-4", "10-15", "05-09"), or
 * a fraction "K/N" selecting the K-th of N equal partitions (1 <= K <= N
 * <= 1000).  Fills len/lo/hi and a canonical display string; returns -1
 * on a malformed spec.
 */
static int pfilter_parse(const char *spec, pfilter *f)
{
    memset(f, 0, sizeof *f);
    const char *slash = strchr(spec, '/');
    if (slash) {
        char *e1, *e2;
        long k = strtol(spec, &e1, 10);
        long n = strtol(slash + 1, &e2, 10);
        if (e1 != slash || *e2 || k < 1 || n < 1 || k > n || n > 1000) return -1;
        int len = 3, scale = 1000;
        if (10 % n == 0)        { len = 1; scale = 10; }
        else if (100 % n == 0)  { len = 2; scale = 100; }
        int lo = (int)((k - 1) * scale / n);
        int hi = (int)(k * scale / n) - 1;
        pfilter_set(f, len, lo, hi);
        f->frac_k = (int)k; f->frac_n = (int)n;
        return 0;
    }
    const char *dash = strchr(spec, '-');
    const char *a = spec, *b = dash ? dash + 1 : spec;
    size_t la = dash ? (size_t)(dash - spec) : strlen(spec), lb = strlen(b);
    if (la < 1 || la > PFILTER_MAX_LEN || lb != la) return -1;
    int lo = 0, hi = 0;
    for (size_t i = 0; i < la; i++) {
        if (a[i] < '0' || a[i] > '9' || b[i] < '0' || b[i] > '9') return -1;
        lo = lo * 10 + (a[i] - '0');
        hi = hi * 10 + (b[i] - '0');
    }
    if (lo > hi) return -1;
    pfilter_set(f, (int)la, lo, hi);
    return 0;
}

/*
 * pfilter_match - true if the window starting at w (digit values 0..9) has
 * its first f->len digits, read as a number, within [lo, hi].  Callers
 * guarantee f->len < window length so every digit read is inside the window.
 */
static inline int pfilter_match(const pfilter *f, const uint8_t *w)
{
    int v = w[0];
    for (int i = 1; i < f->len; i++) v = v * 10 + w[i];
    return v >= f->lo && v <= f->hi;
}

/*
 * pfilter_fraction - the fraction of all windows the filter admits
 * ((hi-lo+1)/10^len), used to scale insertion, fill and repeat estimates.
 */
static inline double pfilter_fraction(const pfilter *f)
{
    if (f->len == 0) return 1.0;
    return (double)(f->hi - f->lo + 1) / pow(10.0, f->len);
}

/*
 * pfilter_describe - "[05-09]" or "[05-09] (2/20)" for banners.
 */
static const char *pfilter_describe(const pfilter *f, char *buf, size_t n)
{
    if (f->frac_n) snprintf(buf, n, "[%s] (%d/%d)", f->str, f->frac_k, f->frac_n);
    else           snprintf(buf, n, "[%s]", f->str);
    return buf;
}

/* ------------------------------------------------------------------ */
/* Window (hi,lo) rolling state                                        */
/* ------------------------------------------------------------------ */
static const uint64_t POW10[20] = {
    1ULL, 10ULL, 100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 10000000ULL,
    100000000ULL, 1000000000ULL, 10000000000ULL, 100000000000ULL, 1000000000000ULL,
    10000000000000ULL, 100000000000000ULL, 1000000000000000ULL, 10000000000000000ULL,
    100000000000000000ULL, 1000000000000000000ULL, 10000000000000000000ULL };

typedef struct {
    int L;          /* sequence length                */
    int Lh;         /* digits in hi (0 if L <= 19)     */
    int Ll;         /* digits in lo (min(L,19))        */
    uint64_t hi_top;/* POW10[Lh-1] (0 if Lh == 0)      */
    uint64_t lo_top;/* POW10[Ll-1]                     */
} winspec;

/*
 * winspec_init - set up the (hi, lo) split for sequence length L: lo holds
 * the last min(L,19) digits, hi the preceding L-19 (0 for L <= 19), plus
 * the powers of ten needed by win_roll.  Dies if L is outside 1..38.
 */
static void winspec_init(winspec *w, int L)
{
    if (L < 1 || L > MAX_SEQ_LEN) die("sequence length must be 1..%d", MAX_SEQ_LEN);
    w->L = L;
    w->Ll = L < LO_DIGITS ? L : LO_DIGITS;
    w->Lh = L - w->Ll;
    w->lo_top = POW10[w->Ll - 1];
    w->hi_top = w->Lh ? POW10[w->Lh - 1] : 0;
}

/*
 * win_compute - compute the (hi, lo) key of the window d[0..L-1] from
 * scratch (L multiply-adds).  Used at the start of each chunk and by
 * pibloomt, which only touches ~1/10 of the windows so rolling would not
 * pay off.
 */
static inline void win_compute(const winspec *w, const uint8_t *d, uint64_t *hi, uint64_t *lo)
{
    uint64_t h = 0, l = 0;
    int i = 0;
    for (; i < w->Lh; i++) h = h * 10 + d[i];
    for (; i < w->L; i++)  l = l * 10 + d[i];
    *hi = h; *lo = l;
}

/*
 * win_roll - advance the (hi, lo) key from the window at d to the window
 * at d+1 in O(1): drop the leading digit of each half, shift, append the
 * digit that enters.  Reads d[L], so callers must not roll past the last
 * window of a buffer.
 */
static inline void win_roll(const winspec *w, const uint8_t *d, uint64_t *hi, uint64_t *lo)
{
    if (w->Lh) {
        *hi = (*hi - (uint64_t)d[0] * w->hi_top) * 10 + d[w->Lh];
    }
    *lo = (*lo - (uint64_t)d[w->Lh] * w->lo_top) * 10 + d[w->L];
}

/*
 * win_parse_ascii - build the (hi, lo) key from L ASCII digits (used when
 * loading candidate files).  Returns -1 without touching the outputs if a
 * non-digit is found, so callers must check the result.
 */
static int win_parse_ascii(const winspec *w, const char *s, uint64_t *hi, uint64_t *lo)
{
    uint8_t d[MAX_SEQ_LEN];
    for (int i = 0; i < w->L; i++) {
        if (s[i] < '0' || s[i] > '9') return -1;
        d[i] = (uint8_t)(s[i] - '0');
    }
    win_compute(w, d, hi, lo);
    return 0;
}

/*
 * win_to_ascii - render a (hi, lo) key back into its L-digit string
 * (leading zeros preserved), NUL-terminated.  Inverse of win_compute.
 */
static void win_to_ascii(const winspec *w, uint64_t hi, uint64_t lo, char *out)
{
    for (int i = w->L - 1; i >= w->Lh; i--) { out[i] = (char)('0' + lo % 10); lo /= 10; }
    for (int i = w->Lh - 1; i >= 0; i--)   { out[i] = (char)('0' + hi % 10); hi /= 10; }
    out[w->L] = 0;
}

/* ------------------------------------------------------------------ */
/* Pi digit sources                                                    */
/* ------------------------------------------------------------------ */
/*
 * Two formats are read natively, auto-detected per file:
 *
 *   TEXT  "3.14159..." (the leading "3." is skipped; newlines/spaces are
 *         ignored).  Continuation files hold raw digits only.
 *   YCD   y-cruncher compressed digit files (#Compressed Digit File,
 *         version 1.0.0/1.1.0, Base 10): a text header terminated by
 *         "EndHeader" and a NUL byte, then little-endian 64-bit words each
 *         holding 19 decimal digits (most-significant digit first).  A file
 *         is block BlockID of Blocksize digits; the last word may be partial
 *         and files are padded, so digits are bounded by the header, never
 *         by file size.  Digit offset 0 is the first digit AFTER the decimal
 *         point (verified against the known digits whenever block 0 is
 *         present; a stream that carries the integer "3" is detected and
 *         handled).
 *
 * A source spec is a comma-separated list of files and/or directories (a
 * directory contributes every *.ycd inside it).  ycd files are ordered by
 * BlockID and must be contiguous; the set may start at any block, in which
 * case positions are absolute (block b starts at position b*Blocksize+1).
 * Text files are read in the order given.
 *
 * Positions everywhere are 1-based from the first digit after the decimal
 * point (OEIS convention).
 */
#include <dirent.h>             /* opendir/readdir to expand a directory of .ycd files */

#define YCD_DIGITS_PER_WORD 19
#define YCD_MAGIC "#Compressed Digit File"

typedef enum { PIF_TEXT = 0, PIF_YCD = 1 } pifile_kind;

typedef struct {
    char *path;
    pifile_kind kind;
    uint64_t data_off;      /* byte offset of the first digit-bearing byte      */
    uint64_t start;         /* 0-based digit offset of the file's first digit   */
    uint64_t ndigits;       /* digits in the file (ycd exact, text = byte count) */
    uint64_t blocksize, blockid, totaldigits;   /* ycd header fields            */
    uint64_t data_bytes;    /* bytes from data_off to end of file (for prefetching) */
} pifile;

typedef struct {
    pifile *files;
    int nfiles, cur;
    pifile_kind kind;
    int fd;
    size_t chunk;           /* max raw bytes per read                           */
    uint8_t *raw;           /* raw read buffer (keeps a partial ycd word)        */
    size_t rawlen;
    uint64_t file_off;      /* byte offset of the next unread byte, current file */
    uint64_t file_left;     /* digits left in current file (text: UINT64_MAX)    */
    uint64_t drop;          /* digits to discard before delivering (seek/lead3)  */
    uint64_t next_pos;      /* 1-based position of the next digit to deliver     */
    uint64_t first_pos;     /* position of the set's first digit                 */
    uint64_t digits_total;  /* digits in the whole set (text: estimate)          */
    uint64_t total_raw, total_skipped;
    int lead3;              /* ycd stream carries the integer digit 3 at offset 0 */
    int text_pure;          /* -1 unknown, else result of pistream_text_pure()   */
    int eof;
    char *pathbuf;
    /* background prefetcher: keeps PREFETCH_DEPTH chunks of the file in the
     * page cache ahead of the consumer so a high-latency source (NFS, slow
     * array) always has several requests in flight */
    uint64_t *cum_bytes;    /* cum_bytes[i] = data bytes of files 0..i-1; [nfiles] = total */
    _Atomic uint64_t consumer_global;   /* consumer's position in that byte space */
    _Atomic uint64_t pf_pos;            /* bytes the prefetcher has pulled into the page cache */
    _Atomic int pf_stop;
    int pf_running;
    pthread_t pf_thread;
} pistream;

#define PREFETCH_DEPTH 8        /* chunks kept ahead of the consumer */

/*
 * ycd_load_le64 - read a little-endian 64-bit word from unaligned bytes
 * (portable; compiles to a plain load on little-endian hosts).
 */
static inline uint64_t ycd_load_le64(const uint8_t *p)
{
    uint64_t w = 0;
    for (int i = 7; i >= 0; i--) w = (w << 8) | p[i];
    return w;
}

/* two-digit lookup: PAIRS[2*v], PAIRS[2*v+1] are the digits of v (0..99) as values 0..9 */
static const uint8_t YCD_PAIRS[200] = {
    0,0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9, 1,0,1,1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,1,9,
    2,0,2,1,2,2,2,3,2,4,2,5,2,6,2,7,2,8,2,9, 3,0,3,1,3,2,3,3,3,4,3,5,3,6,3,7,3,8,3,9,
    4,0,4,1,4,2,4,3,4,4,4,5,4,6,4,7,4,8,4,9, 5,0,5,1,5,2,5,3,5,4,5,5,5,6,5,7,5,8,5,9,
    6,0,6,1,6,2,6,3,6,4,6,5,6,6,6,7,6,8,6,9, 7,0,7,1,7,2,7,3,7,4,7,5,7,6,7,7,7,8,7,9,
    8,0,8,1,8,2,8,3,8,4,8,5,8,6,8,7,8,8,8,9, 9,0,9,1,9,2,9,3,9,4,9,5,9,6,9,7,9,8,9,9 };

/* write the 8 digits of v (< 10^8) to out[0..7] */
static inline void ycd_put8(uint32_t v, uint8_t *out)
{
    uint32_t a = v / 10000u, b = v % 10000u;          /* two 4-digit halves */
    uint32_t a1 = a / 100u, a2 = a % 100u, b1 = b / 100u, b2 = b % 100u;
    out[0] = YCD_PAIRS[2 * a1]; out[1] = YCD_PAIRS[2 * a1 + 1];
    out[2] = YCD_PAIRS[2 * a2]; out[3] = YCD_PAIRS[2 * a2 + 1];
    out[4] = YCD_PAIRS[2 * b1]; out[5] = YCD_PAIRS[2 * b1 + 1];
    out[6] = YCD_PAIRS[2 * b2]; out[7] = YCD_PAIRS[2 * b2 + 1];
}

/*
 * ycd_decode_word - expand one ycd word (< 10^19) into its 19 decimal
 * digits, most significant first, as values 0..9.  Splits the word with two
 * constant divisions (compiled to multiply-shift) into 8 + 8 + 3 digit
 * pieces and finishes with 32-bit arithmetic and two-digit table lookups,
 * so the work is ~20 mostly independent operations instead of a serial
 * chain of 19 divisions - several times faster on a modest server core,
 * where the single reader thread's decode rate bounds the whole pass.
 */
static inline void ycd_decode_word(uint64_t w, uint8_t *out)
{
    uint64_t top = w / 100000000000ULL;                /* digits 0..7   (8 digits) */
    uint64_t rest = w - top * 100000000000ULL;         /* 11 digits              */
    uint32_t mid = (uint32_t)(rest / 1000u);           /* digits 8..15  (8 digits) */
    uint32_t last = (uint32_t)(rest - (uint64_t)mid * 1000u);   /* digits 16..18 */
    ycd_put8((uint32_t)top, out);
    ycd_put8(mid, out + 8);
    uint32_t l1 = last / 100u, l2 = last % 100u;
    out[16] = (uint8_t)l1;
    out[17] = YCD_PAIRS[2 * l2]; out[18] = YCD_PAIRS[2 * l2 + 1];
}

/* reference decoder (simple division chain), kept for self-tests */
static inline void ycd_decode_word_ref(uint64_t w, uint8_t *out)
{
    for (int i = YCD_DIGITS_PER_WORD - 1; i >= 0; i--) { out[i] = (uint8_t)(w % 10); w /= 10; }
}

/* ------------------------------------------------------------------ */
/* Parallel ycd decode pool                                            */
/* ------------------------------------------------------------------ */
/* A 64 MB chunk is 8.4M words / 160M digits; a single reader thread that
 * decodes it all can bound a whole pass on a modest server core.  The pool
 * splits each chunk's full words across a few threads; the main thread
 * dispatches, takes the last share itself, and waits for the rest. */
#define DECODE_MAX_THREADS 8
#define DECODE_MIN_WORDS   (1u << 16)      /* below this, decode inline */

typedef struct {
    pthread_mutex_t m;
    pthread_cond_t start, done;
    int nthreads, gen, pending, stop, started;
    const uint8_t *src;
    uint8_t *dst;
    size_t nwords;
    pthread_t th[DECODE_MAX_THREADS];
    int ids[DECODE_MAX_THREADS];
} decode_pool;

static decode_pool g_dpool = { PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER, PTHREAD_COND_INITIALIZER,
                               0, 0, 0, 0, 0, NULL, NULL, 0, {0}, {0} };

/* decode words [w0, w1) of the current job */
static inline void decode_range(const decode_pool *dp, size_t w0, size_t w1)
{
    const uint8_t *p = dp->src + w0 * 8;
    uint8_t *o = dp->dst + w0 * YCD_DIGITS_PER_WORD;
    for (size_t i = w0; i < w1; i++, p += 8, o += YCD_DIGITS_PER_WORD)
        ycd_decode_word(ycd_load_le64(p), o);
}

/*
 * decode_worker_main - pool thread: wait for a job generation, decode its
 * share of the words, report completion.
 */
static void *decode_worker_main(void *arg)
{
    int id = *(int *)arg;
    decode_pool *dp = &g_dpool;
    int seen = 0;
    for (;;) {
        pthread_mutex_lock(&dp->m);
        while (dp->gen == seen && !dp->stop) pthread_cond_wait(&dp->start, &dp->m);
        if (dp->stop) { pthread_mutex_unlock(&dp->m); return NULL; }
        seen = dp->gen;
        size_t n = dp->nwords, k = (size_t)dp->nthreads + 1;      /* +1: main thread's share */
        pthread_mutex_unlock(&dp->m);
        decode_range(dp, n * (size_t)id / k, n * (size_t)(id + 1) / k);
        pthread_mutex_lock(&dp->m);
        if (--dp->pending == 0) pthread_cond_signal(&dp->done);
        pthread_mutex_unlock(&dp->m);
    }
}

/*
 * decode_threads_wanted - helper threads for decoding: PI_DECODE_THREADS if
 * set, else 3 on machines with 8+ CPUs (4 decoders with the caller), else 0.
 */
static int decode_threads_wanted(void)
{
    const char *e = getenv("PI_DECODE_THREADS");
    if (e) { int v = atoi(e); return v < 1 ? 0 : v > DECODE_MAX_THREADS ? DECODE_MAX_THREADS : v - 1; }
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n >= 8 ? 3 : 0;
}

/*
 * decode_pool_start - create the helper threads once per process.
 */
static void decode_pool_start(void)
{
    decode_pool *dp = &g_dpool;
    if (dp->started) return;
    dp->started = 1;
    dp->nthreads = decode_threads_wanted();
    for (int i = 0; i < dp->nthreads; i++) {
        dp->ids[i] = i;
        if (pthread_create(&dp->th[i], NULL, decode_worker_main, &dp->ids[i])) { dp->nthreads = i; break; }
    }
}

/*
 * decode_words_parallel - decode nwords consecutive full words from src into
 * dst using the pool (or inline when the job is small / no helpers).
 */
static void decode_words_parallel(const uint8_t *src, uint8_t *dst, size_t nwords)
{
    decode_pool *dp = &g_dpool;
    if (dp->nthreads == 0 || nwords < DECODE_MIN_WORDS) {
        dp->src = src; dp->dst = dst;
        decode_range(dp, 0, nwords);
        return;
    }
    pthread_mutex_lock(&dp->m);
    dp->src = src; dp->dst = dst; dp->nwords = nwords;
    dp->pending = dp->nthreads;
    dp->gen++;
    pthread_cond_broadcast(&dp->start);
    pthread_mutex_unlock(&dp->m);
    size_t k = (size_t)dp->nthreads + 1;
    decode_range(dp, nwords * (size_t)dp->nthreads / k, nwords);       /* main thread's share */
    pthread_mutex_lock(&dp->m);
    while (dp->pending) pthread_cond_wait(&dp->done, &dp->m);
    pthread_mutex_unlock(&dp->m);
}

/*
 * pistream_readahead - ask the kernel to start fetching the next few
 * chunks (POSIX_FADV_WILLNEED) so a high-latency source such as an NFS
 * mount keeps streaming while the current chunk is processed.  Took the
 * network read rate from 67 MB/s to the mount ceiling of ~200 MB/s.
 */
static void pistream_readahead(pistream *ps)
{
#ifdef POSIX_FADV_WILLNEED
    posix_fadvise(ps->fd, (off_t)ps->file_off, (off_t)(ps->chunk * 4), POSIX_FADV_WILLNEED);
#else
    (void)ps;
#endif
}

/*
 * pistream_note_position - publish the consumer's current byte position
 * (in the concatenated data-byte space) for the prefetch thread.
 */
static inline void pistream_note_position(pistream *ps)
{
    if (!ps->cum_bytes) return;
    const pifile *f = &ps->files[ps->cur];
    uint64_t in_file = ps->file_off > f->data_off ? ps->file_off - f->data_off : 0;
    atomic_store_explicit(&ps->consumer_global, ps->cum_bytes[ps->cur] + in_file, memory_order_relaxed);
}

/*
 * pistream_prefetch_main - background thread: reads the file data (into a
 * throwaway buffer, so it lands in the page cache) up to PREFETCH_DEPTH
 * chunks ahead of the consumer, following it across files and seeks, and
 * publishes how far it has got (pf_pos).  The consumer waits for pf_pos
 * before each of its own reads, so the source ever sees ONE sequential
 * stream (two readers at different offsets make a spinning-disk array seek
 * itself to a crawl) while the large reads keep several NFS requests in
 * flight.  If the prefetcher cannot read a region it skips it and the
 * consumer reads it itself.
 */
static void *pistream_prefetch_main(void *arg)
{
    pistream *ps = arg;
    uint8_t *buf = malloc(ps->chunk);
    if (!buf) return NULL;
    const uint64_t total = ps->cum_bytes[ps->nfiles];
    const uint64_t depth = (uint64_t)ps->chunk * PREFETCH_DEPTH;
    uint64_t g = atomic_load_explicit(&ps->consumer_global, memory_order_relaxed);
    int fd = -1, fidx = -1;
    while (!atomic_load_explicit(&ps->pf_stop, memory_order_relaxed)) {
        uint64_t c = atomic_load_explicit(&ps->consumer_global, memory_order_relaxed);
        if (g < c || g > c + depth + (uint64_t)ps->chunk) {           /* consumer moved (seek): follow it */
            g = c;
            atomic_store_explicit(&ps->pf_pos, g, memory_order_release);
        }
        if (g >= total) { struct timespec ts = { 0, 20000000 }; nanosleep(&ts, NULL); continue; }
        if (g >= c + depth) { struct timespec ts = { 0, 2000000 }; nanosleep(&ts, NULL); continue; }
        /* locate the file holding byte g */
        int i = fidx >= 0 && g >= ps->cum_bytes[fidx] && g < ps->cum_bytes[fidx + 1] ? fidx : -1;
        if (i < 0) for (i = 0; i < ps->nfiles; i++) if (g < ps->cum_bytes[i + 1]) break;
        if (i >= ps->nfiles) break;
        if (i != fidx) {
            if (fd >= 0) close(fd);
            fd = open(ps->files[i].path, O_RDONLY);
            fidx = i;
            if (fd < 0) { g = ps->cum_bytes[i + 1]; continue; }
        }
        uint64_t in_file = g - ps->cum_bytes[i];
        uint64_t left = ps->cum_bytes[i + 1] - g;
        size_t want = left < ps->chunk ? (size_t)left : ps->chunk;
        ssize_t n = pread(fd, buf, want, (off_t)(ps->files[i].data_off + in_file));
        if (n <= 0) g = ps->cum_bytes[i + 1];                 /* skip on error/EOF: consumer reads it itself */
        else g += (uint64_t)n;
        atomic_store_explicit(&ps->pf_pos, g, memory_order_release);
    }
    atomic_store_explicit(&ps->pf_pos, UINT64_MAX, memory_order_release);   /* gone: consumer never waits */
    if (fd >= 0) close(fd);
    free(buf);
    return NULL;
}

/*
 * pistream_start_prefetch - build the cumulative byte table and start the
 * prefetch thread (disabled when the environment variable PI_NO_PREFETCH is
 * set, for benchmarking the raw reader).
 */
static void pistream_start_prefetch(pistream *ps)
{
    if (getenv("PI_NO_PREFETCH")) return;
    ps->cum_bytes = calloc((size_t)ps->nfiles + 1, sizeof *ps->cum_bytes);
    if (!ps->cum_bytes) die("out of memory");
    for (int i = 0; i < ps->nfiles; i++) ps->cum_bytes[i + 1] = ps->cum_bytes[i] + ps->files[i].data_bytes;
    pistream_note_position(ps);
    atomic_store(&ps->pf_stop, 0);
    atomic_store(&ps->pf_pos, atomic_load(&ps->consumer_global));
    if (pthread_create(&ps->pf_thread, NULL, pistream_prefetch_main, ps) == 0) ps->pf_running = 1;
    else { free(ps->cum_bytes); ps->cum_bytes = NULL; }
}

/*
 * pistream_open_fd - make file idx the current file: open it, position at
 * its first digit-bearing byte, reset the per-file digit budget and the
 * partial-word carry, and prime readahead.
 */
static void pistream_open_fd(pistream *ps, int idx)
{
    if (ps->fd >= 0) close(ps->fd);
    ps->cur = idx;
    pifile *f = &ps->files[idx];
    ps->fd = open(f->path, O_RDONLY);
    if (ps->fd < 0) die("cannot open pi source '%s': %s", f->path, strerror(errno));
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(ps->fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
    fcntl(ps->fd, F_RDAHEAD, 1);
#endif
    ps->file_off = f->data_off;
    ps->file_left = f->kind == PIF_YCD ? f->ndigits : UINT64_MAX;
    ps->rawlen = 0;
    if (lseek(ps->fd, (off_t)f->data_off, SEEK_SET) < 0) die("lseek failed on '%s': %s", f->path, strerror(errno));
    pistream_readahead(ps);
    pistream_note_position(ps);
}

/*
 * ycd_parse_header - parse a y-cruncher .ycd text header (FileVersion,
 * Base, TotalDigits, Blocksize, BlockID, terminated by "EndHeader" and a
 * NUL byte) and derive the file's absolute digit range.  Refuses versions
 * other than 1.0.0/1.1.0 and bases other than 10, and checks the file is
 * not truncated (files are padded, so only a shortfall is detectable).
 */
static void ycd_parse_header(pifile *f)
{
    int fd = open(f->path, O_RDONLY);
    if (fd < 0) die("cannot open '%s': %s", f->path, strerror(errno));
    char hdr[8192];
    ssize_t n = read(fd, hdr, sizeof hdr - 1);
    if (n <= 0) die("cannot read header of '%s'", f->path);
    hdr[n] = 0;
    struct stat st;
    if (fstat(fd, &st) != 0) die("cannot stat '%s'", f->path);
    close(fd);

    /* find the NUL that ends the header */
    ssize_t nul = -1;
    for (ssize_t i = 0; i < n; i++) if (hdr[i] == 0) { nul = i; break; }
    if (nul < 0) die("'%s': ycd header not terminated within %zu bytes", f->path, sizeof hdr);
    f->data_off = (uint64_t)nul + 1;

    char version[32] = "";
    int base = 0, seen_end = 0;
    f->blocksize = 0; f->blockid = 0; f->totaldigits = 0;
    char *save = hdr, *line;
    while ((line = strsep(&save, "\n")) != NULL) {
        size_t len = strlen(line);
        while (len && (line[len - 1] == '\r' || line[len - 1] == ' ')) line[--len] = 0;
        if (!len) continue;
        if (!strcmp(line, "EndHeader")) { seen_end = 1; break; }
        char *colon = strchr(line, ':');
        if (!colon) continue;
        *colon = 0;
        char *val = colon + 1;
        while (*val == '\t' || *val == ' ') val++;
        if (!strcmp(line, "FileVersion")) snprintf(version, sizeof version, "%s", val);
        else if (!strcmp(line, "Base")) base = atoi(val);
        else if (!strcmp(line, "TotalDigits")) f->totaldigits = strtoull(val, NULL, 10);
        else if (!strcasecmp(line, "Blocksize")) f->blocksize = strtoull(val, NULL, 10);
        else if (!strcmp(line, "BlockID")) f->blockid = strtoull(val, NULL, 10);
    }
    if (!seen_end) die("'%s': ycd header has no EndHeader", f->path);
    if (strcmp(version, "1.0.0") && strcmp(version, "1.1.0"))
        die("'%s': unsupported ycd FileVersion '%s'", f->path, version);
    if (base != 10) die("'%s': ycd Base %d is not supported (need 10)", f->path, base);
    if (f->blocksize < 100) die("'%s': invalid ycd Blocksize", f->path);

    uint64_t block_start = f->blockid * f->blocksize;
    uint64_t block_end = block_start + f->blocksize;
    if (f->totaldigits) {
        if (f->totaldigits <= block_start) die("'%s': BlockID %" PRIu64 " lies beyond TotalDigits", f->path, f->blockid);
        if (block_end > f->totaldigits) block_end = f->totaldigits;
    }
    f->start = block_start;
    f->ndigits = block_end - block_start;

    /* sanity against the file size (files are padded, so only a shortfall is checkable) */
    uint64_t avail_words = (uint64_t)st.st_size > f->data_off ? ((uint64_t)st.st_size - f->data_off) / 8 : 0;
    uint64_t need_words = (f->ndigits + YCD_DIGITS_PER_WORD - 1) / YCD_DIGITS_PER_WORD;
    if (avail_words < need_words) {
        if (f->totaldigits)
            die("'%s' is truncated: holds %" PRIu64 " words, block needs %" PRIu64, f->path, avail_words, need_words);
        f->ndigits = avail_words * YCD_DIGITS_PER_WORD;
        warn("'%s' is shorter than Blocksize and TotalDigits is unknown - treating it as the final block with %" PRIu64
             " digits (the last few thousand may be padding)", f->path, f->ndigits);
    }
    f->kind = PIF_YCD;
}

/*
 * pifile_cmp_block - qsort comparator ordering ycd files by BlockID, so
 * command-line/directory order does not matter.
 */
static int pifile_cmp_block(const void *a, const void *b)
{
    const pifile *x = a, *y = b;
    return x->blockid < y->blockid ? -1 : x->blockid > y->blockid;
}

/*
 * pistream_open - open a digit source spec: a comma-separated list of files
 * and/or directories (a directory contributes every *.ycd inside).  Each
 * file is auto-detected as text or ycd (mixing is refused).  ycd sets are
 * ordered by BlockID, checked for contiguity and equal Blocksize, and may
 * start at any block (positions become absolute).  The offset-0 convention
 * is verified against the known digits when block 0 is present.  Text sets
 * skip a leading "3." on the first file only.  Sets first_pos and
 * digits_total and opens the first file.
 */
static void pistream_open(pistream *ps, const char *spec, size_t chunk)
{
    memset(ps, 0, sizeof *ps);
    ps->fd = -1;
    ps->text_pure = -1;
    ps->chunk = chunk;
    ps->raw = malloc(chunk + 16);
    if (!ps->raw) die("cannot allocate %zu byte read buffer", chunk);

    /* expand the spec into a file list */
    int cap = 64;
    ps->files = calloc((size_t)cap, sizeof *ps->files);
    ps->pathbuf = strdup(spec);
    if (!ps->files || !ps->pathbuf) die("out of memory");
    char *save = ps->pathbuf, *tok;
    while ((tok = strsep(&save, ",")) != NULL) {
        if (!*tok) continue;
        struct stat st;
        if (stat(tok, &st) != 0) die("cannot stat pi source '%s': %s", tok, strerror(errno));
        if (S_ISDIR(st.st_mode)) {
            DIR *d = opendir(tok);
            if (!d) die("cannot open directory '%s': %s", tok, strerror(errno));
            struct dirent *de;
            int found = 0;
            while ((de = readdir(d)) != NULL) {
                size_t l = strlen(de->d_name);
                if (l < 5 || strcasecmp(de->d_name + l - 4, ".ycd")) continue;
                if (ps->nfiles == cap) { cap *= 2; ps->files = realloc(ps->files, (size_t)cap * sizeof *ps->files); if (!ps->files) die("out of memory"); }
                size_t pl = strlen(tok) + 1 + l + 1;
                char *p = malloc(pl);
                if (!p) die("out of memory");
                snprintf(p, pl, "%s/%s", tok, de->d_name);
                memset(&ps->files[ps->nfiles], 0, sizeof ps->files[0]);
                ps->files[ps->nfiles++].path = p;
                found++;
            }
            closedir(d);
            if (!found) die("directory '%s' contains no .ycd files", tok);
        } else {
            if (ps->nfiles == cap) { cap *= 2; ps->files = realloc(ps->files, (size_t)cap * sizeof *ps->files); if (!ps->files) die("out of memory"); }
            memset(&ps->files[ps->nfiles], 0, sizeof ps->files[0]);
            ps->files[ps->nfiles++].path = tok;
        }
    }
    if (ps->nfiles == 0) die("empty pi source spec");

    /* classify each file */
    int nycd = 0;
    for (int i = 0; i < ps->nfiles; i++) {
        pifile *f = &ps->files[i];
        int fd = open(f->path, O_RDONLY);
        if (fd < 0) die("cannot open pi source '%s': %s", f->path, strerror(errno));
        char head[32];
        ssize_t n = read(fd, head, sizeof head);
        struct stat st;
        fstat(fd, &st);
        close(fd);
        if (n <= 0) die("pi source '%s' is empty", f->path);
        if ((size_t)n >= strlen(YCD_MAGIC) && !memcmp(head, YCD_MAGIC, strlen(YCD_MAGIC))) {
            ycd_parse_header(f);
            nycd++;
        } else {
            f->kind = PIF_TEXT;
            f->data_off = 0;
            f->ndigits = (uint64_t)st.st_size;
        }
        f->data_bytes = (uint64_t)st.st_size > f->data_off ? (uint64_t)st.st_size - f->data_off : 0;
    }
    if (nycd && nycd != ps->nfiles) die("pi sources mix text and .ycd files - use one format per run");
    ps->kind = nycd ? PIF_YCD : PIF_TEXT;

    if (ps->kind == PIF_YCD) {
        qsort(ps->files, (size_t)ps->nfiles, sizeof *ps->files, pifile_cmp_block);
        for (int i = 1; i < ps->nfiles; i++) {
            pifile *a = &ps->files[i - 1], *b = &ps->files[i];
            if (b->blocksize != a->blocksize)
                die("ycd files '%s' and '%s' have different Blocksize", a->path, b->path);
            if (b->blockid == a->blockid) die("duplicate ycd BlockID %" PRIu64 " ('%s' and '%s')", b->blockid, a->path, b->path);
            if (b->blockid != a->blockid + 1)
                die("ycd set is not contiguous: BlockID %" PRIu64 " ('%s') is followed by %" PRIu64 " ('%s')",
                    a->blockid, a->path, b->blockid, b->path);
            if (a->ndigits != a->blocksize)
                die("ycd file '%s' is short but is not the last block", a->path);
        }
        /* offset-0 convention check when block 0 is present */
        ps->lead3 = 0;
        if (ps->files[0].blockid == 0) {
            int fd = open(ps->files[0].path, O_RDONLY);
            uint8_t w8[8];
            if (fd >= 0 && pread(fd, w8, 8, (off_t)ps->files[0].data_off) == 8) {
                uint64_t w = ycd_load_le64(w8);
                if (w == 1415926535897932384ULL) ps->lead3 = 0;
                else if (w == 3141592653589793238ULL) ps->lead3 = 1;
                else warn("ycd block 0 does not start with the known digits of pi (word=%" PRIu64 ") - is this pi?", w);
            }
            if (fd >= 0) close(fd);
        }
        for (int i = 0; i < ps->nfiles; i++) ps->digits_total += ps->files[i].ndigits;
        if (ps->lead3) ps->digits_total--;                     /* the "3" is not a position */
        ps->first_pos = ps->files[0].start + 1 - (uint64_t)ps->lead3;
        if (ps->lead3 && ps->files[0].blockid == 0) { ps->drop = 1; ps->first_pos = 1; }   /* the "3" is dropped: first digit delivered is position 1 */
    } else {
        /* text: first file may start with "3." / "3"; continuation files are raw digits */
        pifile *f0 = &ps->files[0];
        int fd = open(f0->path, O_RDONLY);
        unsigned char head[4];
        ssize_t n = fd >= 0 ? read(fd, head, sizeof head) : -1;
        if (fd >= 0) close(fd);
        uint64_t skip = 0;
        if (n >= 2 && head[0] == '3' && head[1] == '.') skip = 2;
        else if (n >= 4 && head[0] == '3' && head[1] == '1' && head[2] == '4' && head[3] == '1') skip = 1;
        else if (!(n >= 3 && head[0] == '1' && head[1] == '4' && head[2] == '1'))
            warn("pi source '%s' does not start with '3.1415' or '1415' - assuming it begins at the first digit after the decimal point", f0->path);
        f0->data_off = skip;
        f0->ndigits = f0->ndigits > skip ? f0->ndigits - skip : 0;
        f0->data_bytes = f0->data_bytes > skip ? f0->data_bytes - skip : 0;
        uint64_t acc = 0;
        for (int i = 0; i < ps->nfiles; i++) {
            pifile *f = &ps->files[i];
            if (i > 0) {
                int fd2 = open(f->path, O_RDONLY);
                unsigned char h2[2];
                ssize_t n2 = fd2 >= 0 ? read(fd2, h2, 2) : -1;
                if (fd2 >= 0) close(fd2);
                if (n2 == 2 && h2[0] == '3' && h2[1] == '.')
                    die("continuation file '%s' starts with \"3.\" - it looks like a from-the-start pi file, not a continuation", f->path);
            }
            f->start = acc;
            acc += f->ndigits;
        }
        ps->digits_total = acc;      /* estimate: newlines would overcount slightly */
        ps->first_pos = 1;
    }
    ps->next_pos = ps->first_pos;
    pistream_open_fd(ps, 0);
    ps->total_raw = ps->files[0].data_off;
    if (ps->kind == PIF_YCD) decode_pool_start();
    pistream_start_prefetch(ps);
}

/*
 * pistream_wait_prefetch - block until the prefetcher has pulled the bytes
 * the consumer is about to read into the page cache (or has given up),
 * so the consumer never opens a second disk stream.
 */
static inline void pistream_wait_prefetch(pistream *ps, uint64_t bytes)
{
    if (!ps->pf_running) return;
    uint64_t need = atomic_load_explicit(&ps->consumer_global, memory_order_relaxed) + bytes;
    uint64_t total = ps->cum_bytes[ps->nfiles];
    if (need > total) need = total;         /* the last read wants more than the file has */
    while (atomic_load_explicit(&ps->pf_pos, memory_order_acquire) < need &&
           !atomic_load_explicit(&ps->pf_stop, memory_order_relaxed)) {
        struct timespec ts = { 0, 1000000 };
        nanosleep(&ts, NULL);
    }
}

/*
 * pistream_read - deliver up to max digits (values 0..9) into dst,
 * chaining across files transparently.  Text: reads bytes and drops
 * anything that is not a digit.  ycd: reads whole words, decodes them,
 * keeps a partial trailing word for the next call, and stops exactly at
 * the block's digit count so padding never leaks.  Honours the "drop"
 * counter used by seeks and the leading-3 case.  Returns 0 at end of data.
 * ycd callers must pass max >= 19.
 */
static size_t pistream_read(pistream *ps, uint8_t *dst, size_t max)
{
    size_t produced = 0;
    if (ps->kind == PIF_YCD && max < YCD_DIGITS_PER_WORD) die("internal: pistream_read needs room for a full ycd word");
    while (produced == 0 && !ps->eof) {
        pifile *f = &ps->files[ps->cur];
        if (f->kind == PIF_YCD && ps->file_left == 0) {
            if (ps->cur + 1 < ps->nfiles) { pistream_open_fd(ps, ps->cur + 1); continue; }
            ps->eof = 1; break;
        }
        size_t want;
        if (f->kind == PIF_YCD) {
            size_t words = max / YCD_DIGITS_PER_WORD;
            if (words == 0) words = 1;
            if (words * 8 > ps->chunk) words = ps->chunk / 8;
            want = words * 8;
            if (want > ps->rawlen) want -= ps->rawlen; else want = 8;
        } else {
            want = max < ps->chunk ? max : ps->chunk;
        }
        pistream_wait_prefetch(ps, want);
        ssize_t n = read(ps->fd, ps->raw + ps->rawlen, want);
        if (n < 0) {
            if (errno == EINTR) continue;
            die("read error on pi source '%s': %s", f->path, strerror(errno));
        }
        if (n == 0) {
            if (f->kind == PIF_YCD) die("'%s' ended before its block was complete (%" PRIu64 " digits missing)", f->path, ps->file_left);
            if (ps->cur + 1 < ps->nfiles) { pistream_open_fd(ps, ps->cur + 1); continue; }
            ps->eof = 1; break;
        }
        ps->total_raw += (uint64_t)n;
        ps->file_off += (uint64_t)n;
        pistream_readahead(ps);
        pistream_note_position(ps);

        if (f->kind == PIF_YCD) {
            size_t total = ps->rawlen + (size_t)n;
            size_t nwords = total / 8;
            const uint8_t *p = ps->raw;
            size_t i = 0;
            /* fast path: all the words that decode whole (no drop, no partial
             * word, room in dst) go through the parallel decoder */
            if (ps->drop == 0) {
                size_t full = nwords;
                if (full > ps->file_left / YCD_DIGITS_PER_WORD) full = (size_t)(ps->file_left / YCD_DIGITS_PER_WORD);
                if (full > (max - produced) / YCD_DIGITS_PER_WORD) full = (max - produced) / YCD_DIGITS_PER_WORD;
                if (full) {
                    decode_words_parallel(p, dst + produced, full);
                    produced += full * YCD_DIGITS_PER_WORD;
                    ps->file_left -= full * YCD_DIGITS_PER_WORD;
                    p += full * 8;
                    i = full;
                }
            }
            for (; i < nwords && ps->file_left > 0; i++, p += 8) {
                uint64_t w = ycd_load_le64(p);
                if (ps->file_left >= YCD_DIGITS_PER_WORD && ps->drop == 0 && produced + YCD_DIGITS_PER_WORD <= max) {
                    ycd_decode_word(w, dst + produced);
                    produced += YCD_DIGITS_PER_WORD;
                    ps->file_left -= YCD_DIGITS_PER_WORD;
                } else {
                    uint8_t tmp[YCD_DIGITS_PER_WORD];
                    ycd_decode_word(w, tmp);
                    size_t take = ps->file_left < YCD_DIGITS_PER_WORD ? (size_t)ps->file_left : YCD_DIGITS_PER_WORD;
                    ps->file_left -= take;
                    size_t k = 0;
                    if (ps->drop) { size_t d = ps->drop < take ? (size_t)ps->drop : take; ps->drop -= d; k = d; }
                    for (; k < take && produced < max; k++) dst[produced++] = tmp[k];
                    if (k < take) {
                        /* dst full mid-word: cannot happen, max >= 19 is guaranteed by callers */
                        die("internal: ycd decode overflow");
                    }
                }
            }
            /* keep any partial trailing word */
            size_t used = (size_t)(p - ps->raw);
            size_t left = total - used;
            if (ps->file_left == 0) left = 0;                     /* rest is padding */
            if (left) memmove(ps->raw, ps->raw + used, left);
            ps->rawlen = left;
        } else {
            const uint8_t *src = ps->raw;
            for (ssize_t i = 0; i < n; i++) {
                unsigned c = (unsigned)src[i] - '0';
                if (c < 10u) {
                    if (ps->drop) { ps->drop--; continue; }
                    dst[produced++] = (uint8_t)c;
                } else ps->total_skipped++;
            }
            if (produced == 0 && ps->drop == 0 && ps->total_skipped > (uint64_t)ps->chunk * 4)
                die("pi source '%s' contains no digits after %" PRIu64 " bytes - wrong file?", f->path, ps->total_raw);
        }
    }
    ps->next_pos += produced;
    return produced;
}

/*
 * pistream_text_pure - sample the first 4 MiB of every text file; true if
 * all bytes are digits, meaning byte offset == digit offset and the source
 * can be seeked instead of scanned.  Always true for ycd.
 */
static int pistream_text_pure(pistream *ps)
{
    if (ps->kind != PIF_TEXT) return 1;
    if (ps->text_pure >= 0) return ps->text_pure;       /* cached: seeks are called per position */
    char *buf = malloc(1 << 22);
    if (!buf) die("out of memory");
    int pure = 1;
    for (int i = 0; i < ps->nfiles && pure; i++) {
        int fd = open(ps->files[i].path, O_RDONLY);
        if (fd < 0) die("cannot open '%s': %s", ps->files[i].path, strerror(errno));
        ssize_t n = pread(fd, buf, 1 << 22, (off_t)ps->files[i].data_off);
        close(fd);
        for (ssize_t k = 0; k < n; k++) if (buf[k] < '0' || buf[k] > '9') { pure = 0; break; }
    }
    free(buf);
    ps->text_pure = pure;
    return pure;
}

/*
 * pistream_seek - position the stream so the next digit delivered is at
 * 1-based position pos: finds the containing file, lseeks to the byte
 * (ycd: to the word, discarding the intra-word remainder).  Returns -1 if
 * the source is text with non-digit bytes (caller must read-and-discard).
 * Seeking past the end leaves the stream at EOF.
 */
static int pistream_seek(pistream *ps, uint64_t pos)
{
    if (pos < ps->first_pos) pos = ps->first_pos;
    if (ps->kind == PIF_TEXT && !pistream_text_pure(ps)) return -1;
    uint64_t o = ps->kind == PIF_YCD ? (pos - 1) + (uint64_t)ps->lead3 : pos - 1;   /* 0-based stream offset */
    int idx = -1;
    for (int i = 0; i < ps->nfiles; i++) {
        pifile *f = &ps->files[i];
        if (o >= f->start && o < f->start + f->ndigits) { idx = i; break; }
    }
    if (idx < 0) { ps->eof = 1; ps->next_pos = pos; return 0; }
    pifile *f = &ps->files[idx];
    uint64_t local = o - f->start;
    ps->eof = 0;
    ps->drop = 0;
    pistream_open_fd(ps, idx);
    if (f->kind == PIF_YCD) {
        uint64_t word = local / YCD_DIGITS_PER_WORD;
        ps->file_off = f->data_off + word * 8;
        ps->file_left = f->ndigits - word * YCD_DIGITS_PER_WORD;
        ps->drop = local % YCD_DIGITS_PER_WORD;
    } else {
        ps->file_off = f->data_off + local;
    }
    if (lseek(ps->fd, (off_t)ps->file_off, SEEK_SET) < 0) die("lseek failed on '%s': %s", f->path, strerror(errno));
    pistream_readahead(ps);
    pistream_note_position(ps);
    ps->next_pos = pos;
    return 0;
}

/*
 * pistream_describe - one-line human description of what was opened
 * (format, file/block count, BlockID range, position range) for banners.
 */
static void pistream_describe(const pistream *ps, char *out, size_t n)
{
    char b1[32], b2[32], b3[32];
    if (ps->kind == PIF_YCD)
        snprintf(out, n, "ycd, %d block%s (BlockID %" PRIu64 "-%" PRIu64 ", Blocksize %s), positions %s..%s%s",
                 ps->nfiles, ps->nfiles == 1 ? "" : "s", ps->files[0].blockid, ps->files[ps->nfiles - 1].blockid,
                 fmt_u64(ps->files[0].blocksize, b1, sizeof b1), fmt_u64(ps->first_pos, b2, sizeof b2),
                 fmt_u64(ps->first_pos + ps->digits_total - 1, b3, sizeof b3),
                 ps->files[0].blockid ? "  (block 0 absent: offset-0 convention assumed, cross-check with pifind)" : "");
    else
        snprintf(out, n, "text, %d file%s, ~%s digits", ps->nfiles, ps->nfiles == 1 ? "" : "s",
                 fmt_u64(ps->digits_total, b1, sizeof b1));
}

/*
 * pistream_close - close the current file and free all buffers and path
 * strings (directory-expanded paths were malloc'd individually).
 */
static void pistream_close(pistream *ps)
{
    if (ps->pf_running) {
        atomic_store(&ps->pf_stop, 1);
        pthread_join(ps->pf_thread, NULL);
        ps->pf_running = 0;
    }
    free(ps->cum_bytes); ps->cum_bytes = NULL;
    if (ps->fd >= 0) close(ps->fd);
    ps->fd = -1;
    free(ps->raw); ps->raw = NULL;
    /* paths inside pathbuf are freed with it; directory-expanded paths were malloc'd */
    for (int i = 0; i < ps->nfiles; i++) {
        const char *p = ps->files[i].path;
        if (p < ps->pathbuf || p >= ps->pathbuf + strlen(ps->pathbuf) + 1) free((void *)p);
    }
    free(ps->files); ps->files = NULL;
    free(ps->pathbuf); ps->pathbuf = NULL;
}

/* ------------------------------------------------------------------ */
/* Window chunker on top of a digit stream                             */
/* ------------------------------------------------------------------ */
/*
 * After pireader_fill():
 *    r->buf[0 .. r->ndigits-1]   digits (values 0..9), starting with the last
 *                                L-1 digits of the previous chunk
 *    r->base                     1-based position of r->buf[0]
 *    windows this chunk = r->ndigits - (L-1); window i is at position base+i
 */
typedef struct {
    pistream ps;
    int L;
    size_t chunk;           /* digits per fill                            */
    uint8_t *buf;
    size_t ndigits;
    uint64_t base;
    uint64_t next_pos;
    uint64_t first_pos;     /* position of the source's first digit       */
    uint64_t digits_est;    /* digits in the source (exact for ycd)       */
    uint64_t total_skipped; /* non-digit bytes dropped (text)             */
    uint64_t limit;         /* max digits to deliver (UINT64_MAX = all)   */
    uint64_t delivered;
    int eof, started;
} pireader;

/*
 * pireader_open - open a source for windowed scanning: chunk digits per
 * fill, windows of length L, at most limit_digits delivered in total.
 * Exposes first_pos / digits_est for banners and percent-done.
 */
static void pireader_open(pireader *r, const char *spec, int L, size_t chunk, uint64_t limit_digits)
{
    memset(r, 0, sizeof *r);
    r->L = L;
    r->chunk = chunk;
    r->limit = limit_digits;
    pistream_open(&r->ps, spec, chunk);
    r->buf = malloc(chunk + (size_t)L + YCD_DIGITS_PER_WORD + 64);
    if (!r->buf) die("cannot allocate %zu byte digit buffer", chunk);
    r->first_pos = r->ps.first_pos;
    r->digits_est = r->ps.digits_total;
    r->base = r->next_pos = r->first_pos;
}

/*
 * pireader_fill - load the next chunk.  Carries the last L-1 digits of the
 * previous chunk to the front of the buffer so windows straddling a chunk
 * boundary are complete, then fills up to chunk new digits from the
 * stream (respecting the digit limit).  Returns the number of windows
 * available (window i starts at buf[i], position base+i), 0 at the end.
 */
static size_t pireader_fill(pireader *r)
{
    const int L = r->L;
    size_t keep = 0;
    if (r->started && r->ndigits >= (size_t)(L - 1)) {
        keep = (size_t)(L - 1);
        memmove(r->buf, r->buf + r->ndigits - keep, keep);
        r->base = r->next_pos - keep;
    } else if (r->started) {
        keep = r->ndigits;
        r->base = r->next_pos - keep;
    }
    r->started = 1;
    r->ndigits = keep;

    while (r->ndigits - keep < r->chunk && !r->eof) {
        if (r->delivered >= r->limit) { r->eof = 1; break; }
        size_t room = r->chunk - (r->ndigits - keep);
        uint64_t lim = r->limit - r->delivered;
        if ((uint64_t)room > lim) room = (size_t)lim;
        if (room < YCD_DIGITS_PER_WORD) room = YCD_DIGITS_PER_WORD;   /* decoder needs a full word of space */
        size_t got = pistream_read(&r->ps, r->buf + r->ndigits, room);
        if (got == 0) { r->eof = 1; break; }
        if ((uint64_t)got > lim) got = (size_t)lim;                    /* trim overshoot from the word slack */
        r->ndigits += got;
        r->delivered += got;
        r->next_pos += got;
    }
    r->total_skipped = r->ps.total_skipped;
    if (r->ndigits < (size_t)L) return 0;
    return r->ndigits - (size_t)(L - 1);
}

/*
 * pireader_close - release the stream and the digit buffer.
 */
static void pireader_close(pireader *r)
{
    pistream_close(&r->ps);
    free(r->buf); r->buf = NULL;
}

/* ------------------------------------------------------------------ */
/* Large memory allocation                                             */
/* ------------------------------------------------------------------ */
/* pre-fault worker: first-touch (zero) one slice of the region */
typedef struct {
    char *base;
    uint64_t from, to;
    _Atomic uint64_t *done;
} prefault_arg;

/*
 * prefault_main - thread body for alloc_huge: memsets its slice in 256 MiB
 * steps, publishing progress through the shared atomic counter.
 */
static void *prefault_main(void *a)
{
    prefault_arg *pa = a;
    const uint64_t step = 256ULL << 20;
    for (uint64_t off = pa->from; off < pa->to; ) {
        uint64_t n = pa->to - off < step ? pa->to - off : step;
        memset(pa->base + off, 0, n);
        off += n;
        atomic_fetch_add_explicit(pa->done, n, memory_order_relaxed);
    }
    return NULL;
}

/*
 * prefault_threads - how many threads to use for first-touching: the online
 * CPU count capped at 32 (page-fault handling scales across cores until
 * memory or Optane write bandwidth saturates).
 */
static int prefault_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) n = 1;
    if (n > 32) n = 32;
    return (int)n;
}

/*
 * alloc_huge - allocate a large zeroed region with mmap (requesting
 * transparent huge pages).  With prefault set, every page is touched up
 * front - in parallel across prefault_threads() threads, since the cost is
 * per-page kernel fault handling rather than memset bandwidth - with live
 * percent/rate/ETA output, so page faults never land in the hot loop and an
 * over-committed allocation fails at startup rather than days in.  Dies on
 * failure.
 */
static void *alloc_huge(uint64_t bytes, int prefault, const char *what)
{
    if (bytes == 0) return NULL;
    void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS
#ifdef MAP_NORESERVE
                   | (prefault ? 0 : MAP_NORESERVE)
#endif
                   , -1, 0);
    if (p == MAP_FAILED) {
        char b[32];
        die("cannot allocate %s for %s: %s", fmt_bytes((double)bytes, b, sizeof b), what, strerror(errno));
    }
#ifdef MADV_HUGEPAGE
    madvise(p, bytes, MADV_HUGEPAGE);
#endif
    if (prefault) {
        char b[32], b2[32], b3[32], d1[32];
        int nth = prefault_threads();
        if (bytes < (1ULL << 30)) nth = 1;                   /* not worth threads for small tables */
        printf("Pre-faulting %s for %s (%d thread%s) ...", fmt_bytes((double)bytes, b, sizeof b), what, nth, nth == 1 ? "" : "s");
        fflush(stdout);
        int tty = isatty(1);
        double t0 = now_sec(), tlast = t0;
        int started_progress = 0, last_decile = 0;

        _Atomic uint64_t done = 0;
        pthread_t th[32];
        prefault_arg args[32];
        for (int i = 0; i < nth; i++) {
            args[i].base = (char *)p;
            args[i].from = bytes * (uint64_t)i / (uint64_t)nth;
            args[i].to   = bytes * (uint64_t)(i + 1) / (uint64_t)nth;
            args[i].done = &done;
            if (pthread_create(&th[i], NULL, prefault_main, &args[i])) die("pthread_create(prefault) failed");
        }
        /* progress display while the threads work */
        for (;;) {
            struct timespec ts = { 0, 100 * 1000 * 1000 };
            nanosleep(&ts, NULL);
            uint64_t off = atomic_load_explicit(&done, memory_order_relaxed);
            double t = now_sec();
            if (off >= bytes) break;
            if (t - t0 < 1.0) continue;                       /* only bother once it is taking a while */
            double pct = 100.0 * (double)off / (double)bytes;
            double rate = (double)off / (t - t0);
            double eta = rate > 0 ? (double)(bytes - off) / rate : 0;
            if (tty) {
                if (t - tlast < 0.25) continue;
                if (!started_progress) { printf("\n"); started_progress = 1; }
                printf("\r\033[K  %5.1f%%  %s of %s   %s/s   ETA %s ", pct,
                       fmt_bytes((double)off, b, sizeof b), fmt_bytes((double)bytes, b2, sizeof b2),
                       fmt_bytes(rate, b3, sizeof b3), fmt_duration(eta, d1, sizeof d1));
                fflush(stdout);
                tlast = t;
            } else {
                int decile = (int)(pct / 10.0);
                if (decile > last_decile) {
                    if (!started_progress) { printf("\n"); started_progress = 1; }
                    printf("  %3.0f%%  %s   %s/s   ETA %s\n", 10.0 * decile,
                           fmt_bytes((double)off, b, sizeof b), fmt_bytes(rate, b3, sizeof b3),
                           fmt_duration(eta, d1, sizeof d1));
                    last_decile = decile;
                }
            }
        }
        for (int i = 0; i < nth; i++) pthread_join(th[i], NULL);
        double dt = now_sec() - t0;
        if (started_progress && tty) printf("\r\033[K");
        if (started_progress)
            printf("  pre-fault done: %s in %.1fs (%s/s)\n", fmt_bytes((double)bytes, b, sizeof b), dt,
                   fmt_bytes(dt > 0 ? (double)bytes / dt : 0, b2, sizeof b2));
        else
            printf(" done in %.1fs\n", dt);
    }
    return p;
}

/*
 * phys_mem_bytes - physical RAM size from sysconf (0 if unknown).
 */
static uint64_t phys_mem_bytes(void)
{
    long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE);
    if (pages <= 0 || psz <= 0) return 0;
    return (uint64_t)pages * (uint64_t)psz;
}

/*
 * check_fits_in_ram - refuse an allocation larger than physical RAM and
 * warn above 85%: random access into swapped memory is ~1000x slower
 * than RAM and such a run would effectively never finish.
 */
static void check_fits_in_ram(uint64_t bytes, const char *what)
{
    uint64_t phys = phys_mem_bytes();
    if (!phys) return;
    char b1[32], b2[32];
    if (bytes > phys)
        die("%s (%s) is larger than physical RAM (%s)", what, fmt_bytes((double)bytes, b1, sizeof b1), fmt_bytes((double)phys, b2, sizeof b2));
    if (bytes > phys * 0.85)
        warn("%s (%s) is more than 85%% of physical RAM (%s); if anything else is using memory the run will swap and crawl",
             what, fmt_bytes((double)bytes, b1, sizeof b1), fmt_bytes((double)phys, b2, sizeof b2));
}

/*
 * peak_rss_mb - peak resident set size of this process in MiB (getrusage;
 * units differ between Linux and macOS).
 */
static double peak_rss_mb(void)
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return (double)ru.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double)ru.ru_maxrss / 1024.0;
#endif
}

/* ------------------------------------------------------------------ */
/* Signals                                                             */
/* ------------------------------------------------------------------ */
static volatile sig_atomic_t g_stop = 0;
/*
 * on_signal - SIGINT/SIGTERM handler: the first signal sets g_stop so the
 * scan loops finish the current chunk, flush output and print statistics;
 * a second signal exits immediately.
 */
static void on_signal(int sig)
{
    (void)sig;
    if (g_stop) _exit(130);   /* second Ctrl-C: hard exit */
    g_stop = 1;
}
/*
 * install_signals - route SIGINT and SIGTERM to on_signal.
 */
static void install_signals(void)
{
    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
}

#endif /* PI_COMMON_H */
