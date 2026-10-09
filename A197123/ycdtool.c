/*
 * ycdtool - inspect, extract, check and create y-cruncher .ycd digit files.
 * A small front end to the ycd library (ycd.h / ycd.c); see ycd.md.
 *
 *   ycdtool info    SRC                      headers, block ranges, position range
 *   ycdtool dump    SRC POS [N]              N digits [50] starting at position POS
 *   ycdtool totext  SRC [-s POS] [-n COUNT]  digits as text on stdout
 *   ycdtool check   SRC                      read every word, validate, report speed
 *   ycdtool encode  IN.txt OUTDIR -b BLOCKSIZE [-n COUNT] [-N NAME] [-F FIRSTDIGITS]
 *                                            text digits -> OUTDIR/NAME - <id>.ycd
 *
 * SRC is a .ycd file, a directory of them, or a comma-separated list.
 * Positions are 1-based from the first digit after the radix point.
 *
 * Build: cc -O2 -o ycdtool ycdtool.c ycd.c
 */
#define _POSIX_C_SOURCE 200809L
#include "ycd.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

static void die(const char *msg, const char *detail)
{
    fprintf(stderr, "ycdtool: %s%s%s\n", msg, detail ? ": " : "", detail ? detail : "");
    exit(1);
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void usage(void)
{
    fprintf(stderr,
"ycdtool (ycd library v" YCD_VERSION ") - y-cruncher .ycd digit files\n"
"  ycdtool info    SRC\n"
"  ycdtool dump    SRC POS [N]\n"
"  ycdtool totext  SRC [-s POS] [-n COUNT]\n"
"  ycdtool check   SRC\n"
"  ycdtool encode  IN.txt OUTDIR -b BLOCKSIZE [-n COUNT] [-N NAME] [-F FIRSTDIGITS]\n"
"SRC: a .ycd file, a directory of them, or a comma-separated list.\n"
"Positions are 1-based from the first digit after the radix point.\n");
    exit(1);
}

static ycd_set *open_set(const char *spec)
{
    char err[512];
    ycd_set *s = ycd_open(spec, err, sizeof err);
    if (!s) die(err, NULL);
    return s;
}

static int cmd_info(const char *spec)
{
    ycd_set *s = open_set(spec);
    int n = ycd_nfiles(s);
    const ycd_file *f0 = ycd_get_file(s, 0);
    printf("files        : %d (BlockID %" PRIu64 "..%" PRIu64 ")\n", n, f0->blockid, ycd_get_file(s, n - 1)->blockid);
    printf("format       : FileVersion %s, Base %d, %d digits per 64-bit word\n", f0->version, f0->base, f0->digits_per_word);
    printf("blocksize    : %" PRIu64 " digits\n", f0->blocksize);
    printf("total digits : %" PRIu64 "%s\n", f0->total_digits, f0->total_digits ? "" : " (not recorded in the header)");
    if (f0->first_digits[0]) printf("first digits : %s\n", f0->first_digits);
    printf("positions    : %" PRIu64 " .. %" PRIu64 " (%" PRIu64 " digits)\n", ycd_first_pos(s), ycd_last_pos(s), ycd_ndigits(s));
    if (ycd_convention_checked(s))
        printf("layout       : verified against FirstDigits - %s\n",
               ycd_integer_digits(s) ? "stream carries the integer part (adjusted)" : "word 0 starts with the first digit after the point");
    else
        printf("layout       : not verifiable (block 0 or FirstDigits absent) - assuming word 0 of block 0 is position 1\n");
    if (ycd_length_uncertain(s))
        printf("WARNING      : the last block is shorter than Blocksize and TotalDigits is 0 - its true length is\n"
               "               unknown; the last position above is an upper bound and may include padding\n");
    for (int i = 0; i < n; i++) {
        const ycd_file *f = ycd_get_file(s, i);
        printf("  block %-6" PRIu64 " %" PRIu64 " digits, data at byte %" PRIu64 ", %" PRIu64 " bytes  %s\n",
               f->blockid, f->ndigits, f->data_offset, f->file_size, f->path);
        if (i == 7 && n > 12) { printf("  ... %d more ...\n", n - 10); i = n - 3; }
    }
    ycd_close(s);
    return 0;
}

static int cmd_dump(const char *spec, uint64_t pos, size_t count)
{
    ycd_set *s = open_set(spec);
    ycd_reader *r = ycd_reader_new(s, YCD_ASCII);
    if (!r) die("out of memory", NULL);
    ycd_reader_set_buffer(r, 4096);
    uint8_t *buf = malloc(count + 1);
    if (!buf) die("out of memory", NULL);
    int64_t n = ycd_reader_read_at(r, pos, buf, count);
    if (n < 0) die(ycd_reader_error(r), NULL);
    buf[n] = 0;
    printf("%" PRIu64 ": %s\n", pos, (char *)buf);
    free(buf); ycd_reader_free(r); ycd_close(s);
    return 0;
}

static int cmd_totext(const char *spec, uint64_t pos, uint64_t count)
{
    ycd_set *s = open_set(spec);
    ycd_reader *r = ycd_reader_new(s, YCD_ASCII);
    if (!r) die("out of memory", NULL);
    ycd_reader_set_buffer(r, 8u << 20);
    if (pos && ycd_reader_seek(r, pos)) die(ycd_reader_error(r), NULL);
    size_t cap = 16u << 20;
    uint8_t *buf = malloc(cap);
    if (!buf) die("out of memory", NULL);
    while (count) {
        size_t want = count < cap ? (size_t)count : cap;
        int64_t n = ycd_reader_read(r, buf, want);
        if (n < 0) die(ycd_reader_error(r), NULL);
        if (n == 0) break;
        if (fwrite(buf, 1, (size_t)n, stdout) != (size_t)n) die("write error", NULL);
        count -= (uint64_t)n;
    }
    free(buf); ycd_reader_free(r); ycd_close(s);
    return 0;
}

static int cmd_check(const char *spec)
{
    ycd_set *s = open_set(spec);
    ycd_reader *r = ycd_reader_new(s, 0);
    if (!r) die("out of memory", NULL);
    ycd_reader_set_buffer(r, 8u << 20);
    size_t cap = 32u << 20;
    uint8_t *buf = malloc(cap);
    if (!buf) die("out of memory", NULL);
    uint64_t total = 0, hist[16] = { 0 };
    double t0 = now_sec(), tl = t0;
    int tty = isatty(2);
    for (;;) {
        int64_t n = ycd_reader_read(r, buf, cap);
        if (n < 0) die(ycd_reader_error(r), NULL);
        if (n == 0) break;
        for (int64_t i = 0; i < n; i++) hist[buf[i] & 15]++;
        total += (uint64_t)n;
        double t = now_sec();
        if (tty && t - tl >= 1.0) {
            fprintf(stderr, "\r  %.1f%%  %.0f M digits/s", 100.0 * (double)total / (double)ycd_ndigits(s), (double)total / (t - t0) / 1e6);
            tl = t;
        }
    }
    if (tty) fprintf(stderr, "\r\033[K");
    double el = now_sec() - t0;
    printf("read %" PRIu64 " digits (positions %" PRIu64 "..%" PRIu64 ") in %.1f s = %.0f M digits/s; every word valid\n",
           total, ycd_first_pos(s), ycd_last_pos(s), el, el > 0 ? (double)total / el / 1e6 : 0);
    if (total != ycd_ndigits(s)) die("digit count does not match the headers", NULL);
    printf("digit frequencies:");
    for (int d = 0; d < ycd_base(s); d++) printf(" %c=%.4f%%", d < 10 ? '0' + d : 'a' + d - 10, 100.0 * (double)hist[d] / (double)total);
    printf("\n");
    free(buf); ycd_reader_free(r); ycd_close(s);
    return 0;
}

/*
 * cmd_encode - read decimal digits as text ("3." prefix on the first line is
 * split off into FirstDigits; other non-digit bytes are skipped) and write
 * them as a set of block files.
 */
static int cmd_encode(const char *in, const char *outdir, uint64_t blocksize, uint64_t count, const char *name, const char *firstdigits)
{
    FILE *f = fopen(in, "rb");
    if (!f) die("cannot open input", in);
    size_t cap = 8u << 20;
    uint8_t *buf = malloc(cap), *dig = malloc(cap);
    if (!buf || !dig) die("out of memory", NULL);

    /* integer part: everything before the first '.', if one appears in the first 32 bytes */
    char fd_auto[64] = "";
    size_t n = fread(buf, 1, cap, f);
    size_t skip = 0;
    for (size_t i = 0; i < n && i < 32; i++) if (buf[i] == '.') { skip = i + 1; break; }
    if (skip && !firstdigits) {
        size_t k = 0;
        for (size_t i = 0; i < n && k < 52 + skip; i++) if ((buf[i] >= '0' && buf[i] <= '9') || buf[i] == '.') fd_auto[k++] = (char)buf[i];
        fd_auto[k] = 0;
        firstdigits = fd_auto;
    }

    /* count the digits first: TotalDigits must be recorded, or a short final
     * block's length could not be recovered by any reader */
    uint64_t avail = 0;
    {
        size_t m = n, o2 = skip;
        while (m > 0) {
            for (size_t i = o2; i < m; i++) if (buf[i] >= '0' && buf[i] <= '9') avail++;
            o2 = 0;
            if (avail >= count) break;
            m = fread(buf, 1, cap, f);
        }
        if (count > avail) count = avail;
        rewind(f);
        n = fread(buf, 1, cap, f);
    }

    ycd_write_opts o;
    memset(&o, 0, sizeof o);
    o.base = 10; o.blocksize = blocksize; o.first_digits = firstdigits; o.total_digits = count;
    ycd_writer *w = NULL;
    uint64_t in_block = 0, total = 0;
    char err[512], path[4096];
    size_t off = skip;
    while (n > 0 && total < count) {
        size_t nd = 0;
        for (size_t i = off; i < n; i++) if (buf[i] >= '0' && buf[i] <= '9') dig[nd++] = (uint8_t)(buf[i] - '0');
        off = 0;
        size_t p = 0;
        while (p < nd && total < count) {
            if (!w) {
                snprintf(path, sizeof path, "%s/%s - %" PRIu64 ".ycd", outdir, name, o.blockid);
                w = ycd_writer_create(path, &o, 0, err, sizeof err);
                if (!w) die(err, NULL);
                in_block = 0;
            }
            uint64_t room = blocksize - in_block;
            size_t take = nd - p;
            if (take > room) take = (size_t)room;
            if (take > count - total) take = (size_t)(count - total);
            if (ycd_writer_append(w, dig + p, take)) die(ycd_writer_error(w), NULL);
            p += take; in_block += take; total += take;
            if (in_block == blocksize) {
                if (ycd_writer_close(w)) die("write error", path);
                w = NULL; o.blockid++;
            }
        }
        n = fread(buf, 1, cap, f);
    }
    if (w && ycd_writer_close(w)) die("write error", path);
    fclose(f);
    printf("encoded %" PRIu64 " digits into %" PRIu64 " block file(s) of %" PRIu64 " in %s\n",
           total, o.blockid + (in_block && in_block != blocksize ? 1 : 0), blocksize, outdir);
    free(buf); free(dig);
    return 0;
}

/* parse_count - "123", "10K", "5M", "2B"/"2G", "1T" (powers of 1000) */
static uint64_t parse_count(const char *s)
{
    char *e;
    double v = strtod(s, &e);
    uint64_t mul = 1;
    switch (*e) {
    case 'k': case 'K': mul = 1000ULL; break;
    case 'm': case 'M': mul = 1000000ULL; break;
    case 'b': case 'B': case 'g': case 'G': mul = 1000000000ULL; break;
    case 't': case 'T': mul = 1000000000000ULL; break;
    case 0: break;
    default: die("bad number", s);
    }
    if (v < 0) die("bad number", s);
    return (uint64_t)(v * (double)mul + 0.5);
}

int main(int argc, char **argv)
{
    if (argc < 3) usage();
    const char *cmd = argv[1];
    if (!strcmp(cmd, "info")) return cmd_info(argv[2]);
    if (!strcmp(cmd, "check")) return cmd_check(argv[2]);
    if (!strcmp(cmd, "dump")) {
        if (argc < 4) usage();
        return cmd_dump(argv[2], parse_count(argv[3]), argc > 4 ? (size_t)parse_count(argv[4]) : 50);
    }
    if (!strcmp(cmd, "totext")) {
        uint64_t pos = 0, count = UINT64_MAX;
        for (int i = 3; i + 1 < argc; i += 2) {
            if (!strcmp(argv[i], "-s")) pos = parse_count(argv[i + 1]);
            else if (!strcmp(argv[i], "-n")) count = parse_count(argv[i + 1]);
            else usage();
        }
        return cmd_totext(argv[2], pos, count);
    }
    if (!strcmp(cmd, "encode")) {
        if (argc < 4) usage();
        uint64_t bs = 0, count = UINT64_MAX;
        const char *name = "Digits - Dec", *fdig = NULL;
        for (int i = 4; i + 1 < argc; i += 2) {
            if (!strcmp(argv[i], "-b")) bs = parse_count(argv[i + 1]);
            else if (!strcmp(argv[i], "-n")) count = parse_count(argv[i + 1]);
            else if (!strcmp(argv[i], "-N")) name = argv[i + 1];
            else if (!strcmp(argv[i], "-F")) fdig = argv[i + 1];
            else usage();
        }
        if (!bs) die("encode needs -b BLOCKSIZE", NULL);
        return cmd_encode(argv[2], argv[3], bs, count, name, fdig);
    }
    usage();
    return 1;
}
