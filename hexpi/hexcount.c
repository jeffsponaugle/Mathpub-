/*
 * hexcount - frequency of each hexadecimal digit in the first 10^n digits
 *            of a hex-digit file (e.g. the hexadecimal expansion of Pi).
 *
 * Emits one row per power of ten (n = 1, 2, 3, ...) containing the 16
 * cumulative counts, matching OEIS A099333..A099348 ("Frequency of the
 * hexadecimal <d> in the first 10^n hexadecimal digits of Pi").
 *
 * Per those sequences the integer part is excluded: a leading "3." (or
 * lone "3" followed by a non-hex separator) is detected and skipped, so
 * the count starts at the first fractional digit. Whitespace anywhere in
 * the file is ignored. Any other byte is a fatal error.
 *
 * Usage: hexcount [-b MB] [--csv FILE] [--check-pi] [--no-skip-prefix] FILE
 *        FILE may be "-" for stdin.
 *
 * Build:  cc -O2 -o hexcount hexcount.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

#define WS      16   /* whitespace class */
#define BAD     17   /* invalid byte class */
#define NCLASS  18

static uint8_t lut[256];

static void build_lut(void)
{
    memset(lut, BAD, sizeof lut);
    for (int c = '0'; c <= '9'; c++) lut[c] = (uint8_t)(c - '0');
    for (int c = 'a'; c <= 'f'; c++) lut[c] = (uint8_t)(c - 'a' + 10);
    for (int c = 'A'; c <= 'F'; c++) lut[c] = (uint8_t)(c - 'A' + 10);
    lut[' '] = lut['\t'] = lut['\n'] = lut['\r'] = lut['\v'] = lut['\f'] = WS;
}

/* OEIS A099333..A099348, terms n = 1..12 (rows: digit 0..f, cols: n). */
#define NKNOWN 12
static const uint64_t oeis[16][NKNOWN] = {
    {0,9,59,634,6296,62522,624597,6250690,62501979,625011206,6249979329,62499881108},
    {0,5,68,627,6325,62385,624342,6246592,62498560,624994420,6249991124,62500212206},
    {1,8,62,642,6355,62644,625896,6255492,62508519,624984447,6249917131,62499924780},
    {1,11,69,630,6283,62432,623853,6250592,62514233,624983935,6250077541,62500188844},
    {1,8,60,582,6171,62235,623635,6246996,62482731,625014007,6249925117,62499807368},
    {1,3,64,636,6291,62649,625491,6253329,62501692,625011152,6250052283,62500007205},
    {1,6,71,646,6287,62545,625489,6251102,62501242,624989166,6250025693,62499925426},
    {0,4,53,585,6274,62515,625369,6251231,62504587,625010556,6249958503,62499878794},
    {3,11,81,614,6244,62459,624959,6249740,62495786,625047187,6250044644,62500216752},
    {0,7,65,618,6234,62699,624584,6246654,62494524,624982085,6249968714,62500120671},
    {1,6,65,632,6228,62070,625563,6248822,62497832,624987694,6250055390,62500266095},
    {0,1,60,639,6278,62345,625077,6249493,62489514,624962532,6250019064,62499955595},
    {0,7,56,609,6193,62621,624618,6250491,62500734,625010592,6250006866,62500188610},
    {0,3,50,596,6241,62829,626014,6253342,62504195,625010529,6249933819,62499613666},
    {0,7,59,674,6155,62415,625462,6248121,62501388,624987488,6249989678,62499875079},
    {1,4,58,636,6145,62635,625051,6247313,62502484,625013004,6250055104,62499937801},
};

static uint64_t counts[16];
static uint64_t digits_done;     /* hex digits tallied so far            */
static uint64_t bytes_done;      /* raw bytes consumed (for error msgs)  */
static uint64_t next_boundary = 10;
static int      boundary_n = 1;  /* n of next_boundary = 10^n            */
static FILE    *csv;
static int      check_pi;
static int      check_failures;

static void emit_row(int n, uint64_t ndigits, int final)
{
    if (final)
        printf("%-6s %15llu", "final", (unsigned long long)ndigits);
    else
        printf("10^%-3d %15llu", n, (unsigned long long)ndigits);
    for (int d = 0; d < 16; d++)
        printf(" %12llu", (unsigned long long)counts[d]);
    if (check_pi && !final && n >= 1 && n <= NKNOWN) {
        int ok = 1;
        for (int d = 0; d < 16; d++)
            if (counts[d] != oeis[d][n - 1]) ok = 0;
        printf("  %s", ok ? "OEIS:OK" : "OEIS:MISMATCH");
        if (!ok) check_failures++;
    }
    putchar('\n');
    fflush(stdout);

    if (csv) {
        fprintf(csv, "%llu", (unsigned long long)ndigits);
        for (int d = 0; d < 16; d++)
            fprintf(csv, ",%llu", (unsigned long long)counts[d]);
        fputc('\n', csv);
        fflush(csv);
    }
}

static void header(void)
{
    printf("%-6s %15s", "n", "digits");
    for (int d = 0; d < 16; d++)
        printf(" %12x", d);
    putchar('\n');
}

static void die_bad_byte(unsigned c, uint64_t off)
{
    fprintf(stderr,
        "hexcount: invalid byte 0x%02x at file offset %llu "
        "(only hex digits and whitespace allowed)\n",
        c, (unsigned long long)off);
    exit(1);
}

int main(int argc, char **argv)
{
    size_t bufmb = 16;
    const char *csvpath = NULL, *path = NULL;
    int skip_prefix = 1;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-b") && i + 1 < argc) {
            bufmb = (size_t)strtoul(argv[++i], NULL, 10);
            if (bufmb < 1) bufmb = 1;
        } else if (!strcmp(argv[i], "--csv") && i + 1 < argc) {
            csvpath = argv[++i];
        } else if (!strcmp(argv[i], "--check-pi")) {
            check_pi = 1;
        } else if (!strcmp(argv[i], "--no-skip-prefix")) {
            skip_prefix = 0;
        } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            printf("usage: hexcount [-b MB] [--csv FILE] [--check-pi] "
                   "[--no-skip-prefix] FILE\n");
            return 0;
        } else if (!path) {
            path = argv[i];
        } else {
            fprintf(stderr, "hexcount: unexpected argument '%s'\n", argv[i]);
            return 1;
        }
    }
    if (!path) {
        fprintf(stderr, "usage: hexcount [-b MB] [--csv FILE] [--check-pi] "
                        "[--no-skip-prefix] FILE\n");
        return 1;
    }

    FILE *in = strcmp(path, "-") ? fopen(path, "rb") : stdin;
    if (!in) { perror(path); return 1; }
    if (csvpath) {
        csv = fopen(csvpath, "w");
        if (!csv) { perror(csvpath); return 1; }
        fprintf(csv, "digits,0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f\n");
    }

    build_lut();

    size_t bufsz = bufmb << 20;
    uint8_t *buf = malloc(bufsz);
    if (!buf) { fprintf(stderr, "hexcount: out of memory\n"); return 1; }

    /* Detect and skip an integer-part prefix "3." so counting starts at
     * the first fractional digit, as in OEIS A099333..A099348. */
    if (skip_prefix) {
        int c0 = fgetc(in);
        if (c0 == '3') {
            int c1 = fgetc(in);
            if (c1 == '.') {
                fprintf(stderr, "hexcount: skipped leading \"3.\" "
                                "(integer part of Pi)\n");
                bytes_done = 2;
            } else {
                if (c1 != EOF) ungetc(c1, in);
                ungetc('3', in);
            }
        } else if (c0 != EOF) {
            ungetc(c0, in);
        }
    }

    header();

    time_t t0 = time(NULL), tlast = t0;
    int progress = isatty(fileno(stderr));

    for (;;) {
        size_t got = fread(buf, 1, bufsz, in);
        if (got == 0) break;

        size_t pos = 0;
        while (pos < got) {
            uint64_t room = next_boundary - digits_done;  /* digits to boundary */
            size_t   span = got - pos;
            if (room > span) {
                /* Fast path: boundary cannot fall inside this span. */
                uint64_t tally[NCLASS] = {0};
                const uint8_t *p = buf + pos, *end = buf + got;
                while (p < end)
                    tally[lut[*p++]]++;
                if (tally[BAD]) {
                    for (size_t j = pos; j < got; j++)
                        if (lut[buf[j]] == BAD)
                            die_bad_byte(buf[j], bytes_done + j);
                }
                uint64_t nd = 0;
                for (int d = 0; d < 16; d++) {
                    counts[d] += tally[d];
                    nd += tally[d];
                }
                digits_done += nd;
                pos = got;
            } else {
                /* Careful path: walk bytes until the boundary digit. */
                while (pos < got && digits_done < next_boundary) {
                    uint8_t v = lut[buf[pos]];
                    if (v < 16) {
                        counts[v]++;
                        digits_done++;
                    } else if (v == BAD) {
                        die_bad_byte(buf[pos], bytes_done + pos);
                    }
                    pos++;
                }
                if (digits_done == next_boundary) {
                    emit_row(boundary_n, next_boundary, 0);
                    next_boundary *= 10;
                    boundary_n++;
                }
            }
        }
        bytes_done += got;

        if (progress) {
            time_t now = time(NULL);
            if (now - tlast >= 5) {
                double el = (double)(now - t0);
                fprintf(stderr,
                    "\r%llu digits (%.1f GB, %.0f MB/s)   ",
                    (unsigned long long)digits_done,
                    (double)bytes_done / 1e9,
                    el > 0 ? (double)bytes_done / 1e6 / el : 0.0);
                tlast = now;
            }
        }
    }
    if (progress) fputc('\n', stderr);

    if (ferror(in)) { perror(path); return 1; }

    /* Partial tail (file length not an exact power of ten). */
    if (digits_done != next_boundary / 10 || digits_done == 0)
        emit_row(0, digits_done, 1);

    fprintf(stderr, "hexcount: %llu hex digits processed\n",
            (unsigned long long)digits_done);

    if (in != stdin) fclose(in);
    if (csv) fclose(csv);
    if (check_pi && check_failures) {
        fprintf(stderr, "hexcount: %d row(s) mismatched OEIS reference data\n",
                check_failures);
        return 2;
    }
    return 0;
}
