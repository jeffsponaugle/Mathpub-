/*
 * piverify.c — verify a file of pi digits using the BBP formula.
 *
 * Full C port of piverify.py; see README.md for the math.
 *
 *   hex files  -> true random spot checks against BBP (bbp.c engine).
 *   dec files  -> whole-file "anchor" verification: base 10 has no
 *                 practical digit-extraction formula, so the file's full
 *                 value x = 3 + D/10^N is base-converted at anchor points:
 *                     frac(16^m * x) = ((2^(4m) mod 10^N) * D mod 10^N)/10^N
 *                 and compared with BBP's frac(16^m * pi). Base conversion
 *                 mixes every decimal digit into the result, so one wrong
 *                 digit almost anywhere scrambles the anchors.
 *
 * The big-integer arithmetic for the decimal anchor and for `gen`
 * (Chudnovsky reference digits) uses GMP.
 *
 * Usage:
 *   piverify verify FILE [FILE...] [--base auto|dec|hex] [--checks K]
 *                        [--anchors A] [--digits t<=16] [--seed S]
 *                        [-t|--threads T] [-v]
 *     Multiple FILEs are read in series as one continuous digit stream in
 *     the order given; the first carries the "3." prefix, the rest are raw
 *     digit chunks. Hex mode adds a spot check at every file boundary.
 *   piverify digit POS   [--count t<=16] [-t|--threads T] [-v]
 *   piverify gen N       --out FILE [--base dec|hex] [--line-width W]
 *
 * Exit codes: 0 pass, 1 fail, 2 error.
 */

#include <ctype.h>
#include <getopt.h>
#include <gmp.h>
#include <inttypes.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "bbp.h"

/* First fractional digits of pi, for format/offset sanity checks. */
static const char KNOWN_DEC[] =
    "1415926535897932384626433832795028841971693993751"
    "058209749445923078164062862089986280348253421170679";
static const char KNOWN_HEX[] =
    "243f6a8885a308d313198a2e03707344a4093822299f31d0082efa98ec4e6c89";

/* log16(10) = 1.2041199826559248 as rational upper/lower bounds */
#define L16_HI 12041199827ULL
#define L16_LO 12041199826ULL
#define L16_DEN 10000000000ULL

#define CHUNK_IO (4u << 20)

enum { BASE_AUTO, BASE_DEC, BASE_HEX };

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "error: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(2);
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* thousands-separated u64 for readable output */
static const char *fmt_u64(u64 v) {
    static char bufs[8][32];
    static int slot;
    char raw[24], *out = bufs[slot = (slot + 1) & 7];
    int len = snprintf(raw, sizeof raw, "%" PRIu64, v);
    int oi = 0;
    for (int i = 0; i < len; i++) {
        if (i && (len - i) % 3 == 0)
            out[oi++] = ',';
        out[oi++] = raw[i];
    }
    out[oi] = '\0';
    return out;
}

/* GMP stores one integer in at most 2^31-1 64-bit limbs (~38e9 decimal
 * digits) — the hard ceiling for the decimal anchor's verified span. */
#define GMP_DIGIT_CAP 38000000000ULL

static u64 phys_ram_bytes(void) {
    long pages = sysconf(_SC_PHYS_PAGES), psize = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && psize > 0)
        return (u64)pages * (u64)psize;
    return 0;
}

/* "8000000000", "8g", "500m", "2t" — k/m/g/t = 1e3/1e6/1e9/1e12 */
static u64 parse_size(const char *s) {
    char *end;
    unsigned long long v = strtoull(s, &end, 10);
    u64 mult = 1;
    if (*end) {
        switch (tolower((unsigned char)*end)) {
        case 'k': mult = 1000ULL; break;
        case 'm': mult = 1000000ULL; break;
        case 'g': case 'b': mult = 1000000000ULL; break;
        case 't': mult = 1000000000000ULL; break;
        default: die("bad size \"%s\" (examples: 8000000000, 8g, 2t)", s);
        }
        if (end[1])
            die("bad size \"%s\" (examples: 8000000000, 8g, 2t)", s);
    }
    return (u64)v * mult;
}

/* ------------------------------------------------------------------ */
/* SHA-256 (certificate file fingerprints; no external dependency)     */
/* ------------------------------------------------------------------ */

#define PIVERIFY_VERSION "1.1"

static const uint32_t SHA_K[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1,
    0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786,
    0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147,
    0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
    0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a,
    0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

typedef struct {
    uint32_t h[8];
    uint64_t nbytes;
    unsigned char buf[64];
    size_t fill;
} sha256_t;

static uint32_t ror32(uint32_t x, int r) {
    return (x >> r) | (x << (32 - r));
}

static void sha256_init(sha256_t *s) {
    static const uint32_t h0[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372,
                                   0xa54ff53a, 0x510e527f, 0x9b05688c,
                                   0x1f83d9ab, 0x5be0cd19};
    memcpy(s->h, h0, sizeof h0);
    s->nbytes = 0;
    s->fill = 0;
}

static void sha256_block(sha256_t *s, const unsigned char *p) {
    uint32_t w[64];
    for (int i = 0; i < 16; i++)
        w[i] = (uint32_t)p[4 * i] << 24 | (uint32_t)p[4 * i + 1] << 16 |
               (uint32_t)p[4 * i + 2] << 8 | p[4 * i + 3];
    for (int i = 16; i < 64; i++) {
        uint32_t s0 = ror32(w[i - 15], 7) ^ ror32(w[i - 15], 18) ^
                      (w[i - 15] >> 3);
        uint32_t s1 = ror32(w[i - 2], 17) ^ ror32(w[i - 2], 19) ^
                      (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = s->h[0], b = s->h[1], c = s->h[2], d = s->h[3];
    uint32_t e = s->h[4], f = s->h[5], g = s->h[6], h = s->h[7];
    for (int i = 0; i < 64; i++) {
        uint32_t S1 = ror32(e, 6) ^ ror32(e, 11) ^ ror32(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = h + S1 + ch + SHA_K[i] + w[i];
        uint32_t S0 = ror32(a, 2) ^ ror32(a, 13) ^ ror32(a, 22);
        uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + mj;
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    s->h[0] += a; s->h[1] += b; s->h[2] += c; s->h[3] += d;
    s->h[4] += e; s->h[5] += f; s->h[6] += g; s->h[7] += h;
}

static void sha256_update(sha256_t *s, const unsigned char *p, size_t n) {
    s->nbytes += n;
    if (s->fill) {
        while (n && s->fill < 64) {
            s->buf[s->fill++] = *p++;
            n--;
        }
        if (s->fill == 64) {
            sha256_block(s, s->buf);
            s->fill = 0;
        }
    }
    while (n >= 64) {
        sha256_block(s, p);
        p += 64;
        n -= 64;
    }
    while (n--)
        s->buf[s->fill++] = *p++;
}

static void sha256_hex(sha256_t *s, char out[65]) {
    u64 bits = s->nbytes * 8;
    unsigned char pad[72] = {0x80};
    size_t padlen = (s->fill < 56 ? 56 : 120) - s->fill;
    unsigned char lenb[8];
    for (int i = 0; i < 8; i++)
        lenb[i] = (unsigned char)(bits >> (56 - 8 * i));
    sha256_update(s, pad, padlen);
    sha256_update(s, lenb, 8);
    for (int i = 0; i < 8; i++)
        snprintf(out + 8 * i, 9, "%08x", s->h[i]);
}

static void sha256_file(const char *path, char out[65]) {
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open %s", path);
    unsigned char *buf = malloc(CHUNK_IO);
    sha256_t s;
    sha256_init(&s);
    size_t got;
    while ((got = fread(buf, 1, CHUNK_IO, f)) > 0)
        sha256_update(&s, buf, got);
    free(buf);
    fclose(f);
    sha256_hex(&s, out);
}

/* ------------------------------------------------------------------ */
/* verification certificate                                            */
/* ------------------------------------------------------------------ */

static FILE *g_cert;

/* append one line to the certificate (no-op when none requested) */
static void cert(const char *fmt, ...) {
    if (!g_cert)
        return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(g_cert, fmt, ap);
    va_end(ap);
    fputc('\n', g_cert);
}

/* ------------------------------------------------------------------ */
/* deterministic RNG (splitmix64)                                      */
/* ------------------------------------------------------------------ */

static u64 rng_state;

static u64 rng_next(void) {
    u64 z = (rng_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static u64 rng_below(u64 n) {
    u64 lim = UINT64_MAX - UINT64_MAX % n, r;
    do
        r = rng_next();
    while (r >= lim);
    return r % n;
}

/* ------------------------------------------------------------------ */
/* digit file access                                                   */
/* ------------------------------------------------------------------ */

/* per-byte classes: 1 digit, 0 skip (whitespace . , _), -1 invalid */
static signed char CLS[2][256];

static void cls_init(void) {
    for (int b = 0; b < 2; b++) {
        for (int i = 0; i < 256; i++)
            CLS[b][i] = -1;
        const char *skip = " \t\r\n\v\f.,_";
        for (const char *p = skip; *p; p++)
            CLS[b][(unsigned char)*p] = 0;
        for (const char *p = "0123456789"; *p; p++)
            CLS[b][(unsigned char)*p] = 1;
        if (b == 1)
            for (const char *p = "abcdefABCDEF"; *p; p++)
                CLS[b][(unsigned char)*p] = 1;
    }
}

typedef struct {
    const char *path;
    int base;      /* BASE_DEC or BASE_HEX */
    int skip;      /* leading '3' present in the digit stream */
    int pure;      /* contiguous digits: O(1) seeks possible */
    off_t data0;   /* pure mode: byte offset of first fractional digit */
    u64 n;         /* fractional digit count */
} digitfile;

/* strip skip-chars from head into out (lowercased), return count */
static size_t strip_head(const unsigned char *head, size_t len, int bi,
                         char *out, size_t cap) {
    size_t o = 0;
    for (size_t i = 0; i < len && o < cap; i++) {
        signed char c = CLS[bi][head[i]];
        if (c == 1)
            out[o++] = (char)tolower(head[i]);
        else if (c < 0)
            break; /* purity handled elsewhere; stop at junk */
    }
    return o;
}

static int detect_base(const char *path) {
    unsigned char head[8192];
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open %s", path);
    size_t got = fread(head, 1, sizeof head, f);
    fclose(f);
    if (!got)
        die("file is empty");
    char s[100];
    size_t sl = strip_head(head, got, 1, s, sizeof s - 1); /* hex superset */
    s[sl] = '\0';
    if (!strncmp(s, "3243f6a88", 9) || !strncmp(s, "243f6a88", 8))
        return BASE_HEX;
    if (!strncmp(s, "31415926", 8) || !strncmp(s, "1415926", 7))
        return BASE_DEC;
    die("file does not begin with pi in base 10 or 16 (got \"%.24s...\")", s);
    return 0;
}

static u64 df_count(digitfile *df) {
    fprintf(stderr, "  (formatted file: counting digits with a full pass...)\n");
    FILE *f = fopen(df->path, "rb");
    if (!f)
        die("cannot open %s", df->path);
    unsigned char *buf = malloc(CHUNK_IO);
    int bi = df->base == BASE_HEX;
    u64 n = 0;
    size_t got;
    while ((got = fread(buf, 1, CHUNK_IO, f)) > 0)
        for (size_t i = 0; i < got; i++) {
            signed char c = CLS[bi][buf[i]];
            if (c < 0)
                die("unexpected character 0x%02x in digit file", buf[i]);
            n += (u64)(c == 1);
        }
    free(buf);
    fclose(f);
    return n - (u64)df->skip;
}

/* first=1: the file must open with pi's leading digits (optionally "3.");
 * first=0: a continuation chunk — raw digits, no prefix expected. */
static void df_open(digitfile *df, const char *path, int base, int first) {
    memset(df, 0, sizeof *df);
    df->path = path;
    df->base = base;
    int bi = base == BASE_HEX;
    const char *known = bi ? KNOWN_HEX : KNOWN_DEC;

    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open %s", path);
    unsigned char head[8192];
    size_t got = fread(head, 1, sizeof head, f);
    if (!got)
        die("%s is empty", path);

    if (first) {
        char s[100];
        size_t sl = strip_head(head, got, bi, s, sizeof s - 1);
        s[sl] = '\0';
        if (s[0] == '3' && !strncmp(s + 1, known, 8))
            df->skip = 1;
        else if (!strncmp(s, known, 8))
            df->skip = 0;
        else
            die("%s does not start with the digits of pi (base %s)", path,
                bi ? "16" : "10");
    } else {
        df->skip = 0;
    }

    /* byte-exact prefix, then purity check on the sample */
    size_t i = 0;
    while (i < got && CLS[bi][head[i]] == 0 && head[i] != '.')
        i++;
    off_t data0;
    if (first && i + 1 < got && head[i] == '3' && head[i + 1] == '.')
        data0 = (off_t)i + 2;
    else if (first && df->skip && head[i] == '3')
        data0 = (off_t)i + 1;
    else
        data0 = (off_t)i;
    df->pure = 1;
    for (size_t j = (size_t)data0; j < got; j++)
        if (CLS[bi][head[j]] != 1) {
            df->pure = 0;
            break;
        }

    if (df->pure) {
        if (fseeko(f, 0, SEEK_END))
            die("seek failed on %s", path);
        off_t size = ftello(f);
        /* trailing whitespace */
        off_t tail0 = size > 4096 ? size - 4096 : 0;
        unsigned char tail[4096];
        fseeko(f, tail0, SEEK_SET);
        size_t tlen = fread(tail, 1, sizeof tail, f);
        while (tlen > 0 && CLS[bi][tail[tlen - 1]] == 0 &&
               tail[tlen - 1] != '.')
            tlen--;
        if (tlen > 0 && CLS[bi][tail[tlen - 1]] != 1)
            die("file tail contains non-digit data");
        off_t trail = (off_t)(size - tail0) - (off_t)tlen;
        df->data0 = data0;
        df->n = (u64)(size - data0 - trail);
        fclose(f);
    } else {
        fclose(f);
        df->data0 = -1;
        df->n = df_count(df);
    }
}

typedef struct {
    u64 start, len;
    char *out;  /* len+1 bytes, filled lowercased + NUL */
    char *dest; /* multifile stitching only (see mf_fetch) */
} drange;

static int drange_cmp(const void *a, const void *b) {
    const drange *x = a, *y = b;
    return x->start < y->start ? -1 : x->start > y->start;
}

/* fetch all ranges (sorted by caller via df_fetch) */
static void df_fetch(digitfile *df, drange *rs, size_t nr) {
    qsort(rs, nr, sizeof *rs, drange_cmp);
    int bi = df->base == BASE_HEX;
    if (df->pure) {
        FILE *f = fopen(df->path, "rb");
        if (!f)
            die("cannot open %s", df->path);
        for (size_t i = 0; i < nr; i++) {
            fseeko(f, df->data0 + (off_t)rs[i].start, SEEK_SET);
            if (fread(rs[i].out, 1, rs[i].len, f) != rs[i].len)
                die("short read at digit %s", fmt_u64(rs[i].start + 1));
            for (u64 j = 0; j < rs[i].len; j++) {
                if (CLS[bi][(unsigned char)rs[i].out[j]] != 1)
                    die("non-digit byte at position %s",
                        fmt_u64(rs[i].start + j + 1));
                rs[i].out[j] = (char)tolower((unsigned char)rs[i].out[j]);
            }
            rs[i].out[rs[i].len] = '\0';
        }
        fclose(f);
        return;
    }
    /* streaming scan */
    FILE *f = fopen(df->path, "rb");
    if (!f)
        die("cannot open %s", df->path);
    unsigned char *buf = malloc(CHUNK_IO);
    u64 idx = 0;
    int skip_left = df->skip;
    size_t ri = 0;
    size_t got;
    while (ri < nr && (got = fread(buf, 1, CHUNK_IO, f)) > 0) {
        for (size_t i = 0; i < got && ri < nr; i++) {
            signed char c = CLS[bi][buf[i]];
            if (c < 0)
                die("unexpected character 0x%02x in digit file", buf[i]);
            if (c != 1)
                continue;
            if (skip_left) {
                skip_left--;
                continue;
            }
            for (size_t j = ri; j < nr && rs[j].start <= idx; j++)
                if (idx < rs[j].start + rs[j].len)
                    rs[j].out[idx - rs[j].start] =
                        (char)tolower(buf[i]);
            while (ri < nr && idx + 1 >= rs[ri].start + rs[ri].len)
                if (idx + 1 == rs[ri].start + rs[ri].len)
                    ri++;
                else
                    break;
            idx++;
        }
    }
    free(buf);
    fclose(f);
    if (ri < nr)
        die("EOF before requested digit position");
    for (size_t i = 0; i < nr; i++)
        rs[i].out[rs[i].len] = '\0';
}

/* ------------------------------------------------------------------ */
/* multifile: several chunk files read in series as one digit stream   */
/* ------------------------------------------------------------------ */

typedef struct {
    digitfile *segs;
    u64 *cum; /* nseg+1 entries; cum[s] = global index where seg s starts */
    int nseg;
    u64 n;
} multifile;

static void mf_open(multifile *mf, char **paths, int nf, int base) {
    mf->segs = malloc((size_t)nf * sizeof *mf->segs);
    mf->cum = malloc(((size_t)nf + 1) * sizeof *mf->cum);
    mf->nseg = nf;
    u64 tot = 0;
    for (int i = 0; i < nf; i++) {
        df_open(&mf->segs[i], paths[i], base, i == 0);
        if (mf->segs[i].n == 0)
            die("%s contains no digits", paths[i]);
        mf->cum[i] = tot;
        tot += mf->segs[i].n;
    }
    mf->cum[nf] = tot;
    mf->n = tot;
}

/* fetch global ranges, splitting any range that spans a file boundary */
static void mf_fetch(multifile *mf, drange *rs, size_t nr) {
    drange *sub = malloc(nr * sizeof *sub);
    for (int s = 0; s < mf->nseg; s++) {
        u64 s0 = mf->cum[s], s1 = mf->cum[s + 1];
        size_t ns = 0;
        for (size_t i = 0; i < nr; i++) {
            u64 lo = rs[i].start > s0 ? rs[i].start : s0;
            u64 end = rs[i].start + rs[i].len;
            u64 hi = end < s1 ? end : s1;
            if (lo >= hi)
                continue;
            sub[ns].start = lo - s0;
            sub[ns].len = hi - lo;
            sub[ns].out = malloc((size_t)(hi - lo) + 1);
            sub[ns].dest = rs[i].out + (lo - rs[i].start);
            ns++;
        }
        if (!ns)
            continue;
        df_fetch(&mf->segs[s], sub, ns);
        for (size_t i = 0; i < ns; i++) {
            memcpy(sub[i].dest, sub[i].out, sub[i].len);
            free(sub[i].out);
        }
    }
    free(sub);
    for (size_t i = 0; i < nr; i++)
        rs[i].out[rs[i].len] = '\0';
}

/* load fractional decimal digits from the file series; caller frees.
 * limit = 0 loads everything; otherwise loading stops after `limit`
 * fractional digits and *rest_out receives a byte-based estimate of the
 * unread digits. fcounts[i] receives file i's loaded digit count. */
static char *load_dec_digits(char **paths, int nf, u64 limit, u64 *n_out,
                             u64 *rest_out, u64 *fcounts) {
    off_t *sizes = malloc((size_t)nf * sizeof *sizes);
    off_t total = 0;
    for (int i = 0; i < nf; i++) {
        fcounts[i] = 0;
        FILE *f = fopen(paths[i], "rb");
        if (!f)
            die("cannot open %s", paths[i]);
        fseeko(f, 0, SEEK_END);
        sizes[i] = ftello(f);
        total += sizes[i];
        fclose(f);
    }
    u64 store_max = limit ? limit + 1 : (u64)total; /* +1: leading '3' */
    char *digits = malloc((size_t)store_max + 1);
    unsigned char *buf = malloc(CHUNK_IO);
    if (!digits || !buf)
        die("out of memory (%llu bytes needed for the digit buffer)\n"
            "use --prefix to verify only as many leading digits as fit"
            " in RAM",
            (unsigned long long)store_max + 1);
    u64 n = 0, rest = 0;
    int last_file = 0, stopped = 0;
    for (int i = 0; i < nf && !stopped; i++) {
        FILE *f = fopen(paths[i], "rb");
        if (!f)
            die("cannot open %s", paths[i]);
        u64 n_before = n;
        off_t fpos = 0;
        size_t got;
        while (!stopped && (got = fread(buf, 1, CHUNK_IO, f)) > 0) {
            for (size_t j = 0; j < got; j++) {
                signed char c = CLS[0][buf[j]];
                if (c < 0)
                    die("unexpected character 0x%02x in %s", buf[j],
                        paths[i]);
                if (c != 1)
                    continue;
                if (n == store_max) { /* prefix reached; estimate the rest */
                    rest = (u64)(sizes[i] - fpos - (off_t)j);
                    for (int k = i + 1; k < nf; k++)
                        rest += (u64)sizes[k];
                    stopped = 1;
                    break;
                }
                digits[n++] = (char)buf[j];
            }
            fpos += (off_t)got;
        }
        fclose(f);
        if (n > n_before)
            last_file = i;
        else if (!stopped)
            die("%s contains no digits", paths[i]);
        fcounts[i] = n - n_before;
    }
    free(buf);
    free(sizes);
    digits[n] = '\0';
    if (n > 9 && digits[0] == '3' && !strncmp(digits + 1, KNOWN_DEC, 8)) {
        memmove(digits, digits + 1, n); /* includes NUL */
        n--;
        fcounts[0]--;
    }
    if (n < 9 || strncmp(digits, KNOWN_DEC, 8))
        die("file does not start with 3.14159265...");
    if (limit && n > limit) {
        fcounts[last_file] -= n - limit;
        rest += n - limit;
        n = limit;
        digits[n] = '\0';
    }
    *n_out = n;
    *rest_out = rest;
    return digits;
}

/* ------------------------------------------------------------------ */
/* options                                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    int base;
    int checks;
    int anchors;
    int digits;   /* 0 = default per mode */
    int threads;
    int verbose;
    int has_seed;
    u64 seed;
    int line_width;
    const char *out;
    int no_hash;
    u64 prefix;
} opts_t;

/* ------------------------------------------------------------------ */
/* verify: hex files (random spot checks)                              */
/* ------------------------------------------------------------------ */

static int verify_hex(char **files, int nf, const opts_t *o) {
    int t = o->digits ? o->digits : 8;
    printf("mode: HEX file -> random spot checks with the BBP formula\n\n");

    printf("[1/3] scanning %d file(s)...\n", nf);
    multifile mf;
    mf_open(&mf, files, nf, BASE_HEX);
    cert("MODE: hexadecimal random spot checks (%d hex digits compared"
         " per check)", t);
    cert("");
    for (int i = 0; i < nf; i++)
        cert("  %s: hex digits %s..%s", mf.segs[i].path,
             fmt_u64(mf.cum[i] + 1), fmt_u64(mf.cum[i + 1]));
    cert("  total fractional hex digits: %s", fmt_u64(mf.n));
    cert("");
    for (int i = 0; i < nf && nf > 1; i++)
        printf("      %s: hex digits %s..%s (%s access)\n",
               mf.segs[i].path, fmt_u64(mf.cum[i] + 1),
               fmt_u64(mf.cum[i + 1]),
               mf.segs[i].pure ? "seek" : "streaming");
    printf("      %s fractional hex digits total%s\n", fmt_u64(mf.n),
           nf == 1 ? (mf.segs[0].pure ? " (seek access)"
                                      : " (streaming access)")
                   : "");
    if (mf.n < (u64)(2 * t))
        die("file too small");

    u64 lead_n = mf.n < 64 ? mf.n : 64;
    printf("[2/3] checking first %" PRIu64 " digits against known value...\n",
           lead_n);
    char leadbuf[65];
    drange lead = {0, lead_n, leadbuf, NULL};
    mf_fetch(&mf, &lead, 1);
    if (strncmp(leadbuf, KNOWN_HEX, lead_n)) {
        printf("      FAIL: leading digits are wrong\n");
        cert("CHECKS\n  leading %" PRIu64 " digits vs reference: FAIL",
             lead_n);
        return 1;
    }
    printf("      ok\n");
    cert("CHECKS");
    cert("  leading %" PRIu64 " hex digits vs built-in reference: PASS",
         lead_n);

    /* checks: end of stream (computation errors propagate), every file
     * boundary (stitching/order errors show up there), plus K random */
    size_t np = 0, cap = (size_t)o->checks + (size_t)nf + 1;
    u64 *pos = malloc(cap * sizeof *pos);
    pos[np++] = mf.n - (u64)t;
    int nbound = 0;
    for (int s = 1; s < nf; s++) {
        u64 p = mf.cum[s] <= mf.n - (u64)t ? mf.cum[s] : mf.n - (u64)t;
        int dup = 0;
        for (size_t j = 0; j < np; j++)
            dup |= pos[j] == p;
        if (!dup) {
            pos[np++] = p;
            nbound++;
        }
    }
    for (int i = 0; i < o->checks; i++) {
        u64 p = rng_below(mf.n - (u64)t + 1);
        int dup = 0;
        for (size_t j = 0; j < np; j++)
            dup |= pos[j] == p;
        if (!dup)
            pos[np++] = p;
    }
    drange *rs = calloc(np, sizeof *rs);
    for (size_t i = 0; i < np; i++) {
        rs[i].start = pos[i];
        rs[i].len = (u64)t;
        rs[i].out = malloc((size_t)t + 1);
    }
    printf("[3/3] running %zu spot checks (random%s%s)...\n", np,
           o->has_seed ? ", seeded" : "",
           nbound ? " + file boundaries" : "");
    mf_fetch(&mf, rs, np);

    int ok = 1;
    for (size_t i = 0; i < np; i++) {
        double t0 = now();
        char expect[17];
        bbp_hex_digits(rs[i].start, t, o->threads, o->verbose, expect);
        int match = !strcmp(rs[i].out, expect);
        ok &= match;
        const char *tag = "";
        for (int s = 1; s < nf; s++)
            if (rs[i].start == mf.cum[s])
                tag = "  [file boundary]";
        printf("      hex digits %s..%s: file->%s  bbp->%s  %s  (%.1fs)%s\n",
               fmt_u64(rs[i].start + 1), fmt_u64(rs[i].start + (u64)t),
               rs[i].out, expect, match ? "PASS" : "FAIL", now() - t0, tag);
        cert("  hex digits %s..%s: file value %s, BBP value %s: %s%s",
             fmt_u64(rs[i].start + 1), fmt_u64(rs[i].start + (u64)t),
             rs[i].out, expect, match ? "PASS" : "FAIL", tag);
        fflush(stdout);
    }

    printf("\n");
    if (ok) {
        printf("RESULT: PASS — all sampled positions match BBP.\n");
        printf("  Note: spot checks verify the sampled digits. The end-of-file\n"
               "  check is the strong one: any error in a pi *computation*\n"
               "  propagates to all later digits, so a matching tail plus\n"
               "  random interior samples is the standard verification method.\n");
    } else {
        printf("RESULT: FAIL — at least one sampled position does not match pi.\n");
    }
    cert("");
    cert("INTERPRETATION");
    cert("  Each spot check compares hex digits read from the file series"
         " with the");
    cert("  same digits computed independently by the BBP digit-extraction"
         " formula.");
    cert("  Checks performed: the end of the stream (an error in a pi"
         " computation");
    cert("  propagates to all later digits, so a correct tail is strong"
         " evidence");
    cert("  the whole computation was correct)%s, and %d random"
         " position(s).",
         nbound ? ", every file boundary (where a\n  missing, truncated,"
                  " or out-of-order chunk would show)"
                : "", o->checks);
    cert("  Spot checks prove the sampled digits only; unsampled digits"
         " could still");
    cert("  differ (e.g. from storage corruption). For an exhaustive check,"
         " verify");
    cert("  a decimal representation with anchor mode.");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* verify: decimal files (whole-file anchor verification)              */
/* ------------------------------------------------------------------ */

static int verify_dec(char **files, int nf, const opts_t *o) {
    int t = o->digits ? o->digits : 12;
    printf("mode: DECIMAL file -> whole-file anchor verification against BBP\n"
           "  (base 10 has no practical digit-extraction formula; instead the\n"
           "   file's full value is base-converted at random hex anchor points\n"
           "   and compared with BBP, which catches an error at ANY position)\n\n");

    /* memory feasibility: the anchor needs the whole verified span in
     * memory as one big integer — roughly 3.2 bytes of working memory
     * per verified digit (digit buffer + GMP operands and scratch). */
    off_t total_bytes = 0;
    for (int i = 0; i < nf; i++) {
        FILE *f = fopen(files[i], "rb");
        if (!f)
            die("cannot open %s", files[i]);
        fseeko(f, 0, SEEK_END);
        total_bytes += ftello(f);
        fclose(f);
    }
    u64 want = o->prefix ? o->prefix : (u64)total_bytes;
    u64 ram = phys_ram_bytes();
    u64 need = want * 16 / 5;
    u64 sugg = ram ? ram / 4 : 8000000000ULL;
    sugg -= sugg % 1000000000ULL;
    if (sugg < 1000000000ULL)
        sugg = 1000000000ULL;
    if (sugg > GMP_DIGIT_CAP)
        sugg = GMP_DIGIT_CAP;
    if (want > GMP_DIGIT_CAP && o->prefix)
        die("--prefix %s exceeds GMP's single-integer ceiling (~%s digits)",
            fmt_u64(o->prefix), fmt_u64(GMP_DIGIT_CAP));
    if (!o->prefix && (want > GMP_DIGIT_CAP ||
                       (ram && need > ram - ram / 5)))
        die("exhaustive decimal anchor verification needs the ENTIRE value in\n"
            "memory as one big integer: ~%s digits here means ~%s GB of\n"
            "working memory (this machine has %s GB%s).\n"
            "Options:\n"
            "  --prefix %sg   anchor-verify the first %s billion digits —\n"
            "                 exhaustive within that prefix, fits this"
            " machine\n"
            "  or verify a HEXADECIMAL representation of the same"
            " computation:\n"
            "  BBP spot checks stream at any size with near-zero memory"
            " (see README)",
            fmt_u64(want), fmt_u64(need / 1000000000ULL),
            fmt_u64(ram / 1000000000ULL),
            want > GMP_DIGIT_CAP ? "; GMP also caps one integer at ~38e9"
                                   " digits" : "",
            fmt_u64(sugg / 1000000000ULL), fmt_u64(sugg / 1000000000ULL));
    if (o->prefix && ram && need > ram - ram / 5)
        printf("warning: --prefix %s wants ~%s GB working memory; this"
               " machine has %s GB\n\n",
               fmt_u64(o->prefix), fmt_u64(need / 1000000000ULL),
               fmt_u64(ram / 1000000000ULL));

    double t0 = now();
    printf("[1/4] reading %d file(s)...\n", nf);
    u64 n, rest;
    u64 *fcounts = malloc((size_t)nf * sizeof *fcounts);
    char *digits = load_dec_digits(files, nf, o->prefix, &n, &rest, fcounts);
    cert("MODE: decimal anchor verification (%d anchors, %d hex digits"
         " compared each)", o->anchors, t);
    if (rest)
        cert("      PREFIX MODE: only the first %s digits were read and"
             " verified", fmt_u64(n));
    cert("");
    u64 at = 0;
    for (int i = 0; i < nf; i++) {
        if (!fcounts[i])
            continue; /* beyond the prefix; never read */
        if (nf > 1)
            printf("      %s: %s digits\n", files[i], fmt_u64(fcounts[i]));
        cert("  %s: decimal digits %s..%s", files[i], fmt_u64(at + 1),
             fmt_u64(at + fcounts[i]));
        at += fcounts[i];
    }
    free(fcounts);
    printf("      %s fractional decimal digits loaded (%.1fs)\n",
           fmt_u64(n), now() - t0);
    cert("  fractional decimal digits loaded and verified: %s", fmt_u64(n));
    if (rest) {
        printf("      prefix mode: ~%s further digits were NOT read\n",
               fmt_u64(rest));
        cert("  further digits in the series (NOT read): ~%s", fmt_u64(rest));
    }
    cert("");
    if (n < 2000)
        die("need at least 2000 digits for anchor verification");

    u64 lead_n = n < 100 ? n : 100;
    printf("[2/4] checking first %" PRIu64 " digits against known value...\n",
           lead_n);
    if (strncmp(digits, KNOWN_DEC, lead_n)) {
        printf("      FAIL: leading digits are wrong\n");
        cert("CHECKS\n  leading %" PRIu64 " digits vs reference: FAIL",
             lead_n);
        return 1;
    }
    printf("      ok\n");
    cert("CHECKS");
    cert("  leading %" PRIu64 " decimal digits vs built-in reference: PASS",
         lead_n);

    printf("[3/4] converting file to a %s-digit integer (GMP)...\n",
           fmt_u64(n));
    t0 = now();
    mpz_t D, mod, y, q;
    mpz_inits(D, mod, y, q, NULL);
    if (mpz_set_str(D, digits, 10))
        die("GMP could not parse the digit string");
    free(digits);
    mpz_ui_pow_ui(mod, 10, n);
    printf("      done (%.1fs)\n", now() - t0);

    /* highest anchor with 16^(m+t+20) <= 10^n (file truncation error
     * stays far below the compared digits) */
    u64 m_cap = (u64)((u128)n * L16_DEN / L16_HI) - (u64)t - 20;
    int na = o->anchors > 0 ? o->anchors : 1;
    u64 *anchors = malloc((size_t)na * sizeof *anchors);
    anchors[0] = m_cap;
    for (int i = 1; i < na; i++) {
        u64 m;
        int dup;
        do {
            m = m_cap - 1 - rng_below(511);
            dup = 0;
            for (int j = 0; j < i; j++)
                dup |= anchors[j] == m;
        } while (dup);
        anchors[i] = m;
    }

    printf("[4/4] checking %d anchor(s) near hex position %s...\n", na,
           fmt_u64(m_cap));
    int ok = 1;
    u64 m_min = m_cap;
    for (int i = 0; i < na; i++) {
        u64 m = anchors[i];
        if (m < m_min)
            m_min = m;
        t0 = now();
        mpz_set_ui(y, 2);
        mpz_t e;
        mpz_init_set_ui(e, 4 * m);
        mpz_powm(y, y, e, mod);      /* 2^(4m) mod 10^n */
        mpz_clear(e);
        mpz_mul(y, y, D);
        mpz_mod(y, y, mod);          /* * D mod 10^n */
        mpz_mul_2exp(q, y, 4 * (unsigned)t);
        mpz_tdiv_q(q, q, mod);       /* top t hex digits of the fraction */
        char raw[24], file_hex[20];
        mpz_get_str(raw, 16, q);       /* q < 16^t, so at most t digits */
        size_t rl = strlen(raw);
        size_t pad = rl < (size_t)t ? (size_t)t - rl : 0;
        memset(file_hex, '0', pad);
        memcpy(file_hex + pad, raw, rl + 1);
        double t1 = now();
        char bbp[17];
        bbp_hex_digits(m, t, o->threads, o->verbose, bbp);
        int match = !strcmp(file_hex, bbp);
        ok &= match;
        printf("      anchor hex pos %s: file->%s  bbp->%s  %s"
               "  (convert %.1fs, bbp %.1fs)\n",
               fmt_u64(m), file_hex, bbp, match ? "PASS" : "FAIL", t1 - t0,
               now() - t1);
        cert("  anchor at hex position %s: file value %s, BBP value %s: %s",
             fmt_u64(m), file_hex, bbp, match ? "PASS" : "FAIL");
        fflush(stdout);
    }
    mpz_clears(D, mod, y, q, NULL);

    u64 covered = (u64)((u128)(m_min + (u64)t) * L16_LO / L16_DEN);
    if (covered > n)
        covered = n;
    printf("\n");
    if (ok) {
        printf("RESULT: PASS\n");
        printf("  A single wrong digit at any decimal position 1..%s would have\n"
               "  scrambled every anchor (miss probability ~16^-%d per anchor).\n",
               fmt_u64(covered), t);
        if (covered < n)
            printf("  Note: the final %s digits of the verified span are"
                   " beyond the\n  anchors' resolution and are not covered"
                   " by this check.\n",
                   fmt_u64(n - covered));
        if (rest)
            printf("  PREFIX MODE: only the first %s digits were verified;"
                   " the\n  remaining ~%s digits were not read (an"
                   " exhaustive decimal check\n  needs the whole value in"
                   " memory — see README for options).\n",
                   fmt_u64(n), fmt_u64(rest));
    } else {
        printf("RESULT: FAIL — the file's digits do not match pi.\n"
               "  (An anchor mismatch means at least one digit somewhere in the\n"
               "   covered range is wrong, or digits are missing/inserted.)\n");
    }
    cert("");
    cert("INTERPRETATION");
    cert("  The file series encodes x = 3.d(1) d(2) ... d(N) with"
         " N = %s.", fmt_u64(n));
    cert("  Each anchor compares the leading hexadecimal digits of"
         " frac(16^m * x),");
    cert("  computed from ALL N file digits by modular arithmetic, against");
    cert("  frac(16^m * pi) computed independently with the BBP"
         " digit-extraction");
    cert("  formula. Base conversion mixes every decimal digit into the"
         " result, so");
    if (ok) {
        cert("  a single wrong digit at any decimal position 1..%s would"
             " have", fmt_u64(covered));
        cert("  scrambled every anchor (miss probability ~16^-%d per"
             " anchor).", t);
        if (covered < n)
            cert("  The final %s digits of the verified span lie beyond the"
                 " anchors'\n  resolution and are NOT covered by this"
                 " verification.",
                 fmt_u64(n - covered));
        if (rest)
            cert("  PREFIX MODE: this verification covers only the first %s"
                 " digits\n  of the series; the remaining ~%s digits were"
                 " not read.",
                 fmt_u64(n), fmt_u64(rest));
    } else {
        cert("  at least one digit in decimal positions 1..%s is wrong, or",
             fmt_u64(covered));
        cert("  digits are missing, inserted, or out of order.");
    }
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* gen: Chudnovsky reference digits (independent of BBP)               */
/* ------------------------------------------------------------------ */

#define CHUD_C3_24 10939058860032000ULL /* 640320^3 / 24 */

static void chud_bs(u64 a, u64 b, mpz_t P, mpz_t Q, mpz_t T) {
    if (b - a == 1) {
        if (a == 0) {
            mpz_set_ui(P, 1);
            mpz_set_ui(Q, 1);
        } else {
            mpz_set_ui(P, 6 * a - 5);
            mpz_mul_ui(P, P, 2 * a - 1);
            mpz_mul_ui(P, P, 6 * a - 1);
            mpz_set_ui(Q, a);
            mpz_mul_ui(Q, Q, a);
            mpz_mul_ui(Q, Q, a);
            mpz_mul_ui(Q, Q, CHUD_C3_24);
        }
        mpz_mul_ui(T, P, 13591409UL + 545140134UL * a);
        if (a & 1)
            mpz_neg(T, T);
        return;
    }
    u64 m = (a + b) / 2;
    mpz_t P2, Q2, T2, tmp;
    mpz_inits(P2, Q2, T2, tmp, NULL);
    chud_bs(a, m, P, Q, T);
    chud_bs(m, b, P2, Q2, T2);
    mpz_mul(T, T, Q2);      /* T = T1*Q2 + P1*T2 */
    mpz_mul(tmp, P, T2);
    mpz_add(T, T, tmp);
    mpz_mul(P, P, P2);
    mpz_mul(Q, Q, Q2);
    mpz_clears(P2, Q2, T2, tmp, NULL);
}

/* result = floor(pi * F) */
static void pi_scaled(mpz_t result, const mpz_t F) {
    double prec_digits = (double)mpz_sizeinbase(F, 2) * 0.30103 + 10;
    u64 terms = (u64)(prec_digits / 14.181647462725477) + 2;
    mpz_t P, Q, T, sq;
    mpz_inits(P, Q, T, sq, NULL);
    chud_bs(0, terms, P, Q, T);
    mpz_abs(T, T);
    mpz_mul_2exp(sq, F, 64);        /* sq = sqrt(10005 * (F*2^64)^2) */
    mpz_mul(sq, sq, sq);
    mpz_mul_ui(sq, sq, 10005);
    mpz_sqrt(sq, sq);
    mpz_mul(result, sq, Q);
    mpz_mul_ui(result, result, 426880);
    mpz_mul_2exp(T, T, 64);
    mpz_tdiv_q(result, result, T);
    if (mpz_sgn(result) <= 0)
        die("internal error: Chudnovsky sign");
    mpz_clears(P, Q, T, sq, NULL);
}

static int cmd_gen(u64 n, const opts_t *o) {
    if (n < 100)
        die("gen needs at least 100 digits");
    if (!o->out)
        die("gen requires --out FILE");
    double t0 = now();
    mpz_t F, v;
    mpz_inits(F, v, NULL);
    char *s;
    const char *check;
    size_t check_n;
    if (o->base == BASE_HEX) {
        mpz_set_ui(F, 1);
        mpz_mul_2exp(F, F, 4 * n + 32);
        pi_scaled(v, F);
        mpz_tdiv_q_2exp(v, v, 32);
        s = mpz_get_str(NULL, 16, v);
        check = KNOWN_HEX;
        check_n = n < 64 ? n : 64;
    } else {
        mpz_ui_pow_ui(F, 10, n + 8);
        pi_scaled(v, F);
        mpz_ui_pow_ui(F, 10, 8);
        mpz_tdiv_q(v, v, F);
        s = mpz_get_str(NULL, 10, v);
        check = KNOWN_DEC;
        check_n = n < 100 ? n : 100;
    }
    if (strlen(s) != n + 1 || s[0] != '3' || strncmp(s + 1, check, check_n))
        die("internal self-check failed generating digits");
    FILE *f = fopen(o->out, "w");
    if (!f)
        die("cannot write %s", o->out);
    fputs("3.", f);
    if (o->line_width <= 0)
        fputs(s + 1, f);
    else {
        fputc('\n', f);
        for (u64 i = 1; i <= n; i += (u64)o->line_width) {
            u64 len = n + 1 - i < (u64)o->line_width ? n + 1 - i
                                                     : (u64)o->line_width;
            fwrite(s + i, 1, len, f);
            fputc('\n', f);
        }
    }
    fclose(f);
    printf("wrote %s %s digits of pi to %s (%.1fs)\n", fmt_u64(n),
           o->base == BASE_HEX ? "hex" : "dec", o->out, now() - t0);
    free(s);
    mpz_clears(F, v, NULL);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(void) {
    fprintf(stderr,
        "usage: piverify verify FILE [FILE...] [--base auto|dec|hex]\n"
        "                        [--checks K] [--anchors A] [--digits t<=16]\n"
        "                        [--seed S] [-t|--threads T] [-v]\n"
        "                        [-o|--out CERT.txt] [--no-hash]\n"
        "                        [--prefix N[kmgt]]\n"
        "         (multiple FILEs are read in series as one digit stream,\n"
        "          in the order given; later files are raw digit chunks.\n"
        "          --out writes a verification certificate with sha-256\n"
        "          fingerprints of the inputs; --no-hash skips the hashes.\n"
        "          --prefix limits decimal anchor verification to the first\n"
        "          N digits — for files too large to hold in memory)\n"
        "       piverify digit POS [--count t<=16] [-t|--threads T] [-v]\n"
        "       piverify gen N --out FILE [--base dec|hex] [--line-width W]\n");
    exit(2);
}

int main(int argc, char **argv) {
    cls_init();
    if (argc < 3)
        usage();
    const char *cmd = argv[1];

    opts_t o = {BASE_AUTO, 8, 3, 0, (int)sysconf(_SC_NPROCESSORS_ONLN),
                0, 0, 0, 0, NULL, 0, 0};
    static struct option lo[] = {
        {"base", required_argument, 0, 'b'},
        {"checks", required_argument, 0, 'c'},
        {"anchors", required_argument, 0, 'a'},
        {"digits", required_argument, 0, 'd'},
        {"count", required_argument, 0, 'd'},
        {"seed", required_argument, 0, 's'},
        {"threads", required_argument, 0, 't'},
        {"out", required_argument, 0, 'o'},
        {"line-width", required_argument, 0, 'w'},
        {"no-hash", no_argument, 0, 'H'},
        {"prefix", required_argument, 0, 'P'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}};
    optind = 2;
    int c;
    while ((c = getopt_long(argc, argv, "b:c:a:d:s:t:T:o:w:v", lo, NULL)) != -1)
        switch (c) {
        case 'b':
            if (!strcmp(optarg, "dec")) o.base = BASE_DEC;
            else if (!strcmp(optarg, "hex")) o.base = BASE_HEX;
            else if (!strcmp(optarg, "auto")) o.base = BASE_AUTO;
            else die("bad --base %s", optarg);
            break;
        case 'c': o.checks = atoi(optarg); break;
        case 'a': o.anchors = atoi(optarg); break;
        case 'd': o.digits = atoi(optarg); break;
        case 's': o.has_seed = 1; o.seed = strtoull(optarg, NULL, 10); break;
        case 't':
        case 'T': o.threads = atoi(optarg); break;
        case 'o': o.out = optarg; break;
        case 'w': o.line_width = atoi(optarg); break;
        case 'H': o.no_hash = 1; break;
        case 'P': o.prefix = parse_size(optarg); break;
        case 'v': o.verbose = 1; break;
        default: usage();
        }
    if (optind >= argc)
        usage();
    const char *arg = argv[optind];
    if (o.digits < 0 || o.digits > 16)
        die("--digits/--count must be 1..16");
    rng_state = o.has_seed
                    ? o.seed
                    : (u64)time(NULL) ^ ((u64)getpid() << 32) ^ (u64)clock();

    if (!strcmp(cmd, "verify")) {
        char **files = &argv[optind];
        int nf = argc - optind;
        int base = o.base == BASE_AUTO ? detect_base(files[0]) : o.base;
        u64 seed_used = rng_state;
        printf("piverify: %s%s\n", files[0], nf > 1 ? " (+ more)" : "");
        printf("engine: C | threads: %d | GMP: %s\n\n", o.threads,
               gmp_version);
        if (o.out) {
            g_cert = fopen(o.out, "w");
            if (!g_cert)
                die("cannot write certificate to %s", o.out);
            char ts[64], host[256] = "?";
            time_t tt = time(NULL);
            struct tm tmv;
            gmtime_r(&tt, &tmv);
            strftime(ts, sizeof ts, "%Y-%m-%dT%H:%M:%SZ", &tmv);
            gethostname(host, sizeof host - 1);
            cert("================================================================");
            cert("              PIVERIFY VERIFICATION CERTIFICATE");
            cert("================================================================");
            cert("tool:        piverify %s (BBP digit-extraction"
                 " verification of pi)", PIVERIFY_VERSION);
            cert("date (UTC):  %s", ts);
            cert("host:        %s", host);
            fputs("command:     ", g_cert);
            for (int i = 0; i < argc; i++)
                fprintf(g_cert, "%s%s", i ? " " : "", argv[i]);
            fputc('\n', g_cert);
            cert("threads:     %d | GMP: %s", o.threads, gmp_version);
            cert("seed:        %" PRIu64 "%s", seed_used,
                 o.has_seed ? "" : "  (auto; rerun with --seed to reproduce)");
            cert("");
            cert("INPUT FILES (%d, read in series)", nf);
            if (!o.no_hash)
                printf("computing sha-256 fingerprints of %d file(s)...\n\n",
                       nf);
            for (int i = 0; i < nf; i++) {
                FILE *sf = fopen(files[i], "rb");
                if (!sf)
                    die("cannot open %s", files[i]);
                fseeko(sf, 0, SEEK_END);
                off_t sz = ftello(sf);
                fclose(sf);
                cert("  %d. %s (%s bytes)", i + 1, files[i],
                     fmt_u64((u64)sz));
                if (!o.no_hash) {
                    char hex[65];
                    sha256_file(files[i], hex);
                    cert("     sha-256: %s", hex);
                }
            }
            cert("");
        }
        int rc = base == BASE_DEC ? verify_dec(files, nf, &o)
                                  : verify_hex(files, nf, &o);
        if (g_cert) {
            cert("");
            cert("----------------------------------------------------------------");
            cert("RESULT: %s", rc == 0 ? "PASS" : "FAIL");
            cert("================================================================");
            fclose(g_cert);
            g_cert = NULL;
            printf("\ncertificate written to %s\n", o.out);
        }
        return rc;
    }
    if (argc - optind > 1)
        die("unexpected extra argument \"%s\"", argv[optind + 1]);
    if (!strcmp(cmd, "digit")) {
        char *end;
        unsigned long long p = strtoull(arg, &end, 10);
        if (*end || p < 1 || p > (1ULL << 59))
            die("bad position %s", arg);
        int t = o.digits ? o.digits : 8;
        double t0 = now();
        char buf[17];
        bbp_hex_digits((u64)p - 1, t, o.threads, o.verbose, buf);
        printf("hex digits %s..%s of pi (after the point): %s   [%.1fs]\n",
               fmt_u64(p), fmt_u64(p + (u64)t - 1), buf, now() - t0);
        return 0;
    }
    if (!strcmp(cmd, "gen")) {
        char *end;
        unsigned long long ngen = strtoull(arg, &end, 10);
        if (*end || ngen < 1)
            die("bad digit count %s", arg);
        if (o.base == BASE_AUTO)
            o.base = BASE_DEC;
        return cmd_gen((u64)ngen, &o);
    }
    usage();
    return 2;
}
