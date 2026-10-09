/*
 * pipal.c — search the decimal expansion of Pi for the terms of OEIS A279885:
 *           "First n-digit palindrome in the decimal expansion of Pi."
 *
 * Digits are counted with the leading '3' as digit 0 (the OEIS convention:
 * a(1) = 3).  Terms may begin with '0' (e.g. a(10) = 0136776310), so results
 * are reported as digit strings, not numbers.
 *
 * Designed for very large digit files (multiple terabytes): the file is
 * streamed once, single-threaded, with constant memory.  Any non-digit bytes
 * (the "3." decimal point, newlines, spaces) are ignored, so raw y-cruncher
 * output, "3.14159..." text files, and wrapped-line files all work.
 *
 * Algorithm: for each palindrome center (both odd- and even-length centers),
 * expand outward while digits match, capped at MAXN/2.  Random digits
 * mismatch with probability 0.9, so the expected cost is ~2.2 comparisons
 * per digit regardless of file size.  Centers are processed strictly left to
 * right, so the first center whose maximal palindrome covers a still-missing
 * length n yields the first n-digit palindrome in the expansion.
 *
 * Usage: pipal [-o outfile] [-n maxlen] pifile
 *
 * Build:  cc -O2 -o pipal pipal.c
 */

#define _FILE_OFFSET_BITS 64

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MAXN_CAP 128            /* hard cap on -n */
#ifndef DBUF_CAP
#define DBUF_CAP (4u << 20)     /* digit window: 4 MiB */
#endif
#ifndef RAW_CAP
#define RAW_CAP  (8u << 20)     /* raw read buffer: 8 MiB */
#endif
#define REPORT_EVERY 1000000000ull

static int      maxn = 64;      /* longest palindrome length tracked (-n) */
static long     maxr;           /* maxn / 2 */
static FILE    *out;
static int      stderr_tty;

static char     found[MAXN_CAP + 1];
static int      remaining;      /* lengths 1..maxn not yet found */
static int      min_odd, min_even; /* smallest missing length per parity; 0 = none missing */

static void update_mins(void)
{
    min_odd = min_even = 0;
    for (int n = 1; n <= maxn; n += 2)
        if (!found[n]) { min_odd = n; break; }
    for (int n = 2; n <= maxn; n += 2)
        if (!found[n]) { min_even = n; break; }
}

static void report(int n, uint64_t start, const uint8_t *digits)
{
    found[n] = 1;
    remaining--;
    fprintf(out, "%d\t%" PRIu64 "\t%.*s\n", n, start, n, digits);
    fflush(out);
    if (out != stdout || !isatty(STDOUT_FILENO))
        fprintf(stderr, "%sa(%d) = %.*s  (digit index %" PRIu64 ")\n",
                stderr_tty ? "\r\033[K" : "", n, n, digits, start);
}

/* Maximal odd palindrome centered at buffer index i (absolute index c). */
static void center_odd(const uint8_t *d, size_t i, size_t dlen, uint64_t c)
{
    long rmax = maxr;
    if ((long)i < rmax)              rmax = (long)i;
    if ((long)(dlen - 1 - i) < rmax) rmax = (long)(dlen - 1 - i);
    long r = 0;
    while (r < rmax && d[i - r - 1] == d[i + r + 1])
        r++;
    long len = 2 * r + 1;
    if (min_odd == 0 || len < min_odd)
        return;
    for (int n = min_odd; n <= len; n += 2)
        if (!found[n])
            report(n, c - (uint64_t)(n - 1) / 2, d + i - (n - 1) / 2);
    update_mins();
}

/* Maximal even palindrome centered between buffer indices i and i+1. */
static void center_even(const uint8_t *d, size_t i, size_t dlen, uint64_t c)
{
    long rmax = maxr;
    if ((long)i + 1 < rmax)          rmax = (long)i + 1;
    if ((long)(dlen - 1 - i) < rmax) rmax = (long)(dlen - 1 - i);
    long r = 0;
    while (r < rmax && d[i - r] == d[i + 1 + r])
        r++;
    long len = 2 * r;
    if (min_even == 0 || len < min_even)
        return;
    for (int n = min_even; n <= len; n += 2)
        if (!found[n])
            report(n, c - (uint64_t)n / 2 + 1, d + i - n / 2 + 1);
    update_mins();
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "usage: %s [-o outfile] [-n maxlen] pifile\n"
        "  pifile      file of decimal digits of Pi (non-digit bytes ignored)\n"
        "  -o outfile  write results there instead of stdout\n"
        "  -n maxlen   longest palindrome length to track (default %d, max %d)\n",
        prog, maxn, MAXN_CAP);
    exit(2);
}

int main(int argc, char **argv)
{
    const char *outpath = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "o:n:h")) != -1) {
        switch (opt) {
        case 'o': outpath = optarg; break;
        case 'n':
            maxn = atoi(optarg);
            if (maxn < 1 || maxn > MAXN_CAP) {
                fprintf(stderr, "pipal: -n must be 1..%d\n", MAXN_CAP);
                return 2;
            }
            break;
        default: usage(argv[0]);
        }
    }
    if (optind != argc - 1)
        usage(argv[0]);
    const char *pipath = argv[optind];

    maxr = maxn / 2;
    remaining = maxn;
    update_mins();
    stderr_tty = isatty(STDERR_FILENO);

    int fd = open(pipath, O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "pipal: cannot open %s: %s\n", pipath, strerror(errno));
        return 1;
    }
#ifdef POSIX_FADV_SEQUENTIAL
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
#ifdef F_RDAHEAD
    fcntl(fd, F_RDAHEAD, 1);
#endif

    out = stdout;
    if (outpath) {
        out = fopen(outpath, "w");
        if (!out) {
            fprintf(stderr, "pipal: cannot open %s: %s\n", outpath, strerror(errno));
            return 1;
        }
    }
    fprintf(out, "# A279885: first n-digit palindrome in the decimal expansion of Pi\n");
    fprintf(out, "# start_index is 0-based; digit 0 is the leading 3\n");
    fprintf(out, "# n\tstart_index\tpalindrome\n");
    fflush(out);

    uint8_t *dbuf = malloc(DBUF_CAP);
    uint8_t *raw  = malloc(RAW_CAP);
    if (!dbuf || !raw) {
        fprintf(stderr, "pipal: out of memory\n");
        return 1;
    }

    uint64_t base = 0;          /* absolute digit index of dbuf[0] */
    size_t   dlen = 0;          /* digits currently buffered */
    uint64_t next_c = 0;        /* next center (absolute index) to process */
    uint64_t next_report = REPORT_EVERY;
    int      eof = 0;
    struct timespec t0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    for (;;) {
        /* Refill the digit window, filtering out non-digit bytes. */
        while (!eof && dlen < DBUF_CAP) {
            size_t want = DBUF_CAP - dlen;
            if (want > RAW_CAP)
                want = RAW_CAP;
            ssize_t got = read(fd, raw, want);
            if (got < 0) {
                if (errno == EINTR)
                    continue;
                fprintf(stderr, "pipal: read error: %s\n", strerror(errno));
                return 1;
            }
            if (got == 0) { eof = 1; break; }
            for (ssize_t k = 0; k < got; k++) {
                uint8_t ch = raw[k];
                if ((uint8_t)(ch - '0') <= 9)
                    dbuf[dlen++] = ch;
            }
        }

        uint64_t total = base + dlen;
        if (stderr_tty && total >= next_report) {
            struct timespec t1;
            clock_gettime(CLOCK_MONOTONIC, &t1);
            double secs = (double)(t1.tv_sec - t0.tv_sec) +
                          (double)(t1.tv_nsec - t0.tv_nsec) * 1e-9;
            fprintf(stderr, "\r\033[K%.1fG digits, %.0f Mdigit/s, %d length(s) pending",
                    (double)total / 1e9, (double)total / secs / 1e6, remaining);
            next_report = (total / REPORT_EVERY + 1) * REPORT_EVERY;
        }

        /* Centers are complete once maxr digits of right context exist
         * (even centers additionally read one digit past the center,
         * hence maxr + 1); at EOF everything left is processable. */
        uint64_t limit = total;
        if (!eof)
            limit = (total > (uint64_t)maxr + 1) ? total - (uint64_t)maxr - 1 : 0;

        while (next_c < limit) {
            size_t i = (size_t)(next_c - base);
            if (min_odd)
                center_odd(dbuf, i, dlen, next_c);
            if (min_even && i + 1 < dlen)
                center_even(dbuf, i, dlen, next_c);
            next_c++;
        }

        if (remaining == 0 || (eof && next_c >= total))
            break;

        /* Slide the window: unprocessed centers need maxr digits of
         * left context. */
        size_t keep_from = 0;
        if (next_c - base > (uint64_t)maxr)
            keep_from = (size_t)(next_c - base) - (size_t)maxr;
        if (keep_from > 0) {
            memmove(dbuf, dbuf + keep_from, dlen - keep_from);
            dlen -= keep_from;
            base += keep_from;
        }
    }

    uint64_t total = base + dlen;
    if (stderr_tty)
        fprintf(stderr, "\r\033[K");
    fprintf(stderr, "done: scanned %" PRIu64 " digits, %d of %d lengths found\n",
            total, maxn - remaining, maxn);

    fprintf(out, "# scanned %" PRIu64 " digits of Pi\n", total);
    if (remaining) {
        fprintf(out, "# no n-digit palindrome found for n =");
        for (int n = 1; n <= maxn; n++)
            if (!found[n])
                fprintf(out, " %d", n);
        fprintf(out, "\n");
    }
    if (out != stdout)
        fclose(out);
    return 0;
}
