/*
 * ycd_test - self-test of the ycd library: word codec against a reference,
 * write/read round trips in base 10 and 16 with awkward block sizes, seeks
 * at every kind of boundary, a set given in shuffled order, integer-part
 * detection, and the error paths.  Exit status 0 = all passed.
 *
 * Build: cc -O2 -o ycd_test ycd_test.c ycd.c      Run: ./ycd_test [tmpdir]
 */
#define _POSIX_C_SOURCE 200809L
#include "ycd.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failures;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__); printf("\n"); } } while (0)

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { uint64_t x = rng_state; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return rng_state = x; }

static void test_codec(void)
{
    uint8_t a[19], b[19];
    uint64_t edge[] = { 0, 1, 9, 10, 99999999ULL, 100000000ULL, 99999999999ULL, 100000000000ULL, 999ULL, 1000ULL,
                        9999999999999999999ULL, 1415926535897932384ULL, 5000000000000000000ULL };
    for (size_t i = 0; i < sizeof edge / sizeof *edge + 2000000; i++) {
        uint64_t w = i < sizeof edge / sizeof *edge ? edge[i] : rnd() % YCD_DEC_WORD_LIMIT;
        ycd_decode_dec(w, a);
        uint64_t v = w;
        for (int k = 18; k >= 0; k--) { b[k] = (uint8_t)(v % 10); v /= 10; }
        if (memcmp(a, b, 19)) { CHECK(0, "decimal decode of %" PRIu64, w); return; }
        if (ycd_encode_dec(a) != w) { CHECK(0, "decimal encode of %" PRIu64, w); return; }
    }
    for (int i = 0; i < 100000; i++) {
        uint64_t w = rnd();
        ycd_decode_hex(w, a);
        CHECK(ycd_encode_hex(a) == w, "hex round trip");
    }
    uint8_t le[8]; ycd_store_le64(le, 0x0102030405060708ULL);
    CHECK(le[0] == 8 && le[7] == 1 && ycd_load_le64(le) == 0x0102030405060708ULL, "little-endian helpers");
}

/* write `total` digits as blocks of `bs`, returning the digits */
static uint8_t *make_set(const char *dir, int base, uint64_t bs, uint64_t total, int with_total, const char *first, const uint8_t *prefix, size_t nprefix)
{
    uint8_t *d = malloc(total);
    for (uint64_t i = 0; i < total; i++) d[i] = i < nprefix ? prefix[i] : (uint8_t)(rnd() % (uint64_t)base);
    char err[256], path[1024];
    uint64_t nblocks = (total + bs - 1) / bs;
    for (uint64_t b = 0; b < nblocks; b++) {
        ycd_write_opts o; memset(&o, 0, sizeof o);
        o.base = base; o.blocksize = bs; o.blockid = b; o.total_digits = with_total ? total : 0; o.first_digits = first;
        snprintf(path, sizeof path, "%s/blk%03d.ycd", dir, (int)((b * 7 + 3) % nblocks));   /* names unrelated to block order */
        ycd_writer *w = ycd_writer_create(path, &o, 0, err, sizeof err);
        CHECK(w != NULL, "create %s: %s", path, err);
        if (!w) continue;
        uint64_t n = total - b * bs < bs ? total - b * bs : bs;
        /* append in uneven pieces */
        uint64_t off = 0;
        while (off < n) { uint64_t k = 1 + rnd() % 977; if (k > n - off) k = n - off; CHECK(ycd_writer_append(w, d + b * bs + off, (size_t)k) == 0, "append"); off += k; }
        if (b + 1 == nblocks && n == bs) CHECK(ycd_writer_append(w, d, 1) == -1, "overflow must be refused");
        int rc = ycd_writer_close(w);
        CHECK(rc == ((n < bs && !with_total) ? -1 : 0), "close status for block %d", (int)b);
    }
    return d;
}

static void rmdir_all(const char *dir)
{
    char cmd[1200]; snprintf(cmd, sizeof cmd, "rm -rf '%s'", dir);
    if (system(cmd)) { /* best effort */ }
}

static void test_roundtrip(const char *tmp, int base, uint64_t bs, uint64_t total)
{
    char dir[1024], err[256];
    snprintf(dir, sizeof dir, "%s/ycdtest_b%d_%" PRIu64, tmp, base, bs);
    rmdir_all(dir);
    char cmd[1200]; snprintf(cmd, sizeof cmd, "mkdir -p '%s'", dir);
    if (system(cmd)) { CHECK(0, "mkdir"); return; }
    uint8_t *d = make_set(dir, base, bs, total, 1, NULL, NULL, 0);
    ycd_set *s = ycd_open(dir, err, sizeof err);
    CHECK(s != NULL, "open %s: %s", dir, err);
    if (!s) { free(d); return; }
    CHECK(ycd_first_pos(s) == 1 && ycd_last_pos(s) == total && ycd_ndigits(s) == total, "range 1..%" PRIu64 " got %" PRIu64 "..%" PRIu64, total, ycd_first_pos(s), ycd_last_pos(s));
    CHECK(ycd_base(s) == base && !ycd_length_uncertain(s), "base/length flags");
    for (int i = 1; i < ycd_nfiles(s); i++) CHECK(ycd_get_file(s, i)->blockid == ycd_get_file(s, i - 1)->blockid + 1, "block order");

    /* sequential, with odd request sizes and a small buffer */
    ycd_reader *r = ycd_reader_new(s, 0);
    ycd_reader_set_buffer(r, 256);
    uint8_t *got = malloc(total + 64);
    uint64_t n = 0;
    for (;;) { int64_t k = ycd_reader_read(r, got + n, (size_t)(1 + rnd() % 300)); CHECK(k >= 0, "read: %s", ycd_reader_error(r)); if (k <= 0) break; n += (uint64_t)k; }
    CHECK(n == total && !memcmp(got, d, total), "sequential read of base %d set (%" PRIu64 " of %" PRIu64 ")", base, n, total);
    CHECK(ycd_reader_tell(r) == total + 1, "tell at end");

    /* random access, including block seams and the end */
    int dpw = base == 10 ? 19 : 16;
    for (int i = 0; i < 4000; i++) {
        uint64_t pos = 1 + rnd() % total;
        if (i % 4 == 0) { uint64_t b = 1 + rnd() % ((total + bs - 1) / bs); pos = b * bs - (rnd() % (uint64_t)(2 * dpw)) + (uint64_t)dpw; if (pos < 1 || pos > total) pos = total; }
        size_t want = (size_t)(1 + rnd() % 80);
        int64_t k = ycd_reader_read_at(r, pos, got, want);
        uint64_t exp = total - pos + 1 < want ? total - pos + 1 : want;
        if (k != (int64_t)exp || memcmp(got, d + pos - 1, (size_t)exp)) { CHECK(0, "read_at %" PRIu64 " x %zu (got %" PRId64 ")", pos, want, k); break; }
    }
    CHECK(ycd_reader_read_at(r, total + 1, got, 5) == 0, "read at end returns 0");
    CHECK(ycd_reader_seek(r, total + 2) == -1 && ycd_reader_seek(r, 0) == -1, "seek outside refused");

    /* text output */
    ycd_reader *ra = ycd_reader_new(s, YCD_ASCII);
    ycd_reader_read_at(ra, 1, got, 16);
    for (int i = 0; i < 16; i++) { char c = (char)(d[i] < 10 ? '0' + d[i] : 'a' + d[i] - 10); CHECK(got[i] == (uint8_t)c, "ascii digit %d", i); }
    ycd_reader_free(ra);
    ycd_reader_free(r);

    /* a sub-set without block 0 has absolute positions */
    if (ycd_nfiles(s) >= 3) {
        const char *two[2] = { ycd_get_file(s, 2)->path, ycd_get_file(s, 1)->path };
        ycd_set *sub = ycd_open_paths(two, 2, err, sizeof err);
        CHECK(sub != NULL, "open subset: %s", err);
        if (sub) {
            CHECK(ycd_first_pos(sub) == bs + 1, "subset first_pos");
            ycd_reader *rs = ycd_reader_new(sub, 0);
            size_t want = 2 * bs < 40 ? (size_t)(2 * bs) : 40;
            CHECK(ycd_reader_read_at(rs, bs + 1, got, want) == (int64_t)want && !memcmp(got, d + bs, want), "subset read");
            CHECK(ycd_reader_seek(rs, bs) == -1, "subset refuses earlier positions");
            ycd_reader_free(rs); ycd_close(sub);
        }
        const char *gap[2] = { ycd_get_file(s, 0)->path, ycd_get_file(s, 2)->path };
        CHECK(ycd_open_paths(gap, 2, err, sizeof err) == NULL, "gap in blocks must be refused");
    }
    ycd_close(s); free(got); free(d);
    rmdir_all(dir);
}

static void test_integer_part(const char *tmp)
{
    /* a foreign encoder that put the "3" of pi into the stream: positions must still start at the "1" */
    char dir[1024], err[256], cmd[1200];
    snprintf(dir, sizeof dir, "%s/ycdtest_int", tmp); rmdir_all(dir);
    snprintf(cmd, sizeof cmd, "mkdir -p '%s'", dir); if (system(cmd)) { CHECK(0, "mkdir"); return; }
    const uint8_t pi[] = { 3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3,2,3,8,4,6,2,6,4,3,3,8,3,2,7,9,5 };
    uint8_t *d = make_set(dir, 10, 500, 1200, 1, "3.14159265358979323846264338327950288419716939937510", pi, sizeof pi);
    ycd_set *s = ycd_open(dir, err, sizeof err);
    CHECK(s != NULL, "open: %s", err);
    if (s) {
        CHECK(ycd_convention_checked(s) && ycd_integer_digits(s) == 1, "integer digit detected");
        CHECK(ycd_first_pos(s) == 1 && ycd_last_pos(s) == 1199, "positions exclude the integer digit");
        uint8_t got[64]; ycd_reader *r = ycd_reader_new(s, 0);
        CHECK(ycd_reader_read_at(r, 1, got, 5) == 5 && got[0] == 1 && got[1] == 4 && got[4] == 9, "position 1 is the first fractional digit");
        CHECK(ycd_reader_read_at(r, 500, got, 10) == 10 && !memcmp(got, d + 500, 10), "seam with integer offset");
        ycd_reader_free(r); ycd_close(s);
    }
    free(d); rmdir_all(dir);

    /* the normal layout is recognised too */
    snprintf(cmd, sizeof cmd, "mkdir -p '%s'", dir); if (system(cmd)) return;
    d = make_set(dir, 10, 500, 1000, 0, "3.14159265358979323846264338327950288419716939937510", pi + 1, sizeof pi - 1);
    s = ycd_open(dir, err, sizeof err);
    CHECK(s && ycd_convention_checked(s) && ycd_integer_digits(s) == 0 && ycd_last_pos(s) == 1000, "fractional layout verified");
    ycd_close(s); free(d); rmdir_all(dir);
}

static void test_errors(const char *tmp)
{
    char dir[1024], path[1100], err[256], cmd[1200];
    snprintf(dir, sizeof dir, "%s/ycdtest_err", tmp); rmdir_all(dir);
    snprintf(cmd, sizeof cmd, "mkdir -p '%s'", dir); if (system(cmd)) { CHECK(0, "mkdir"); return; }
    CHECK(ycd_open(dir, err, sizeof err) == NULL && strstr(err, "no .ycd"), "empty directory: %s", err);
    snprintf(path, sizeof path, "%s/nope.ycd", dir);
    CHECK(ycd_open(path, err, sizeof err) == NULL, "missing file");
    FILE *f = fopen(path, "w"); fputs("3.14159", f); fclose(f);
    CHECK(ycd_is_ycd(path) == 0 && ycd_open(path, err, sizeof err) == NULL, "text file refused");

    /* a short final block without TotalDigits: the writer refuses it (checked in
     * make_set); a reader that meets one anyway flags the length as uncertain
     * when the file is visibly shorter than a block */
    unlink(path);
    uint8_t *d = make_set(dir, 10, 1000000, 1000050, 0, NULL, NULL, 0);
    ycd_set *s = ycd_open(dir, err, sizeof err);
    CHECK(s && ycd_length_uncertain(s) && ycd_last_pos(s) >= 1000050 && ycd_last_pos(s) < 2000000, "uncertain length flagged");
    ycd_close(s);

    /* truncation is detected when TotalDigits is known */
    rmdir_all(dir); if (system(cmd)) return;
    free(d); d = make_set(dir, 10, 100000, 100000, 1, NULL, NULL, 0);
    snprintf(path, sizeof path, "%s/blk000.ycd", dir);
    CHECK(truncate(path, 20000) == 0, "truncate");
    CHECK(ycd_open(dir, err, sizeof err) == NULL && strstr(err, "truncated"), "truncated file: %s", err);
    free(d); rmdir_all(dir);
}

int main(int argc, char **argv)
{
    const char *tmp = argc > 1 ? argv[1] : "/tmp";
    test_codec();
    test_roundtrip(tmp, 10, 1000, 4321);          /* blocks not a multiple of 19, short last block */
    test_roundtrip(tmp, 10, 19 * 50, 19 * 200);   /* everything word-aligned */
    test_roundtrip(tmp, 10, 7, 100);              /* blocks shorter than a word */
    test_roundtrip(tmp, 10, 100003, 1000000);
    test_roundtrip(tmp, 16, 1000, 4321);
    test_roundtrip(tmp, 16, 64, 640);
    test_integer_part(tmp);
    test_errors(tmp);
    printf(failures ? "ycd_test: %d FAILURE(S)\n" : "ycd_test: all passed\n", failures);
    return failures ? 1 : 0;
}
