/*
 * lemverify.c — verify files of Lemniscate-constant digits (y-cruncher
 * output) as far as mathematics allows.
 *
 * THE FUNDAMENTAL LIMIT: unlike pi, the lemniscate constant has NO known
 * BBP-type digit-extraction formula (it is essentially Gamma(1/4)^2, which
 * lies outside the polylogarithm class where all known BBP formulas live).
 * Digits at an arbitrary position can NOT be spot-checked independently.
 * Complete verification of a large computation requires a second full
 * computation with a mathematically independent formula (y-cruncher offers
 * several: Zuniga, Gauss, Sebah, AGM-Pi, Series-Pi) — that is how the
 * published records are verified, and `compare` mode is built for it.
 *
 * What CAN be verified without a second full computation:
 *
 *   verify  — (a) the first K digits against an independent computation
 *             done here with MPFR via  2*pi / AGM(1, sqrt 2)  (independent
 *             library and formula); (b) for decimal files, 10-digit windows
 *             at up to 64 known checkpoint positions reaching 600 billion,
 *             from y-cruncher's distributed reference table (values from
 *             independent record computations by different people/machines).
 *   cross   — mutual base-conversion verification between the DEC and HEX
 *             files over their first N digits (RAM-limited), in BOTH
 *             directions:
 *                 dec->hex:  frac(16^m * x_dec) computed from ALL N decimal
 *                            digits by modular arithmetic must equal the hex
 *                            digits read at file position m+1;
 *                 hex->dec:  frac(10^r * x_hex) likewise vs the dec file.
 *             Base conversion mixes every source digit into the anchor, so
 *             one wrong digit anywhere in the covered prefix of either file
 *             scrambles it (miss probability ~ base^-t per anchor).
 *   compare — stream-compare two digit files (e.g. the outputs of two
 *             y-cruncher runs with different formulas): the actual
 *             full-verification method. Reports the first mismatch.
 *   scan    — whole-file streaming integrity: digit counts, frequency
 *             chi-square, invalid bytes, sha-256.
 *   gen     — reference digits via MPFR (test data / small cross-checks).
 *   digits  — print digits at a position (computes the whole prefix; cheap
 *             only for positions up to ~10^8).
 *
 * The constant: y-cruncher's "Lemniscate" is the arc length of the unit
 * lemniscate, L = 2*varpi = 5.24411510858423962... ; files holding
 * varpi = 2.62205755429211981... are auto-detected and handled too.
 *
 * Exit codes: 0 pass, 1 fail, 2 error.
 */

#include <ctype.h>
#include <getopt.h>
#include <gmp.h>
#include <mpfr.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

#define LEMVERIFY_VERSION "1.0"

/* ------------------------------------------------------------------ */
/* the constant: known leading fractional digits (both variants/bases) */
/* ------------------------------------------------------------------ */

enum { VAR_2PV, VAR_PV };            /* 5.244... (2*varpi) | 2.622... (varpi) */
enum { BASE_AUTO, BASE_DEC, BASE_HEX };

/* first 100 fractional digits (verified against y-cruncher's reference
 * table; the gen/verify self-checks recompute these with MPFR) */
static const char KNOWN_DEC_2PV[] =
    "2441151085842396209296791797822388273655099028632463256336434076"
    "015811741408285004605910659228581868";
static const char KNOWN_HEX_2PV[] =
    "3e7e53e7d42c1b9e61ddaeb0236785edae2e5c91c86e6b5ed68d8f5f225a21ce"
    "002f689ed621841cceabd03b24e29a059c51";
static const char KNOWN_DEC_PV[] =
    "6220575542921198104648395898911194136827549514316231628168217038"
    "007905870704142502302955329614290934";
static const char KNOWN_HEX_PV[] =
    "9f3f29f3ea160dcf30eed75811b3c2f6d7172e48e43735af6b46c7af912d10e7"
    "0017b44f6b10c20e6755e81d92714d02ce28";

static const char INT_DIGIT[2] = {'5', '2'};

static const char *known_frac(int base, int variant) {
    if (base == BASE_HEX)
        return variant == VAR_2PV ? KNOWN_HEX_2PV : KNOWN_HEX_PV;
    return variant == VAR_2PV ? KNOWN_DEC_2PV : KNOWN_DEC_PV;
}

/* y-cruncher reference table for L = 2*varpi, base 10: the 10 decimal
 * digits ENDING at each (1-based) fractional position. Deep entries come
 * from independent record computations (Watkins, Trueb, Kim, ...). */
typedef struct {
    u64 pos;
    const char d[11];
} checkpoint;

static const checkpoint CHECKPOINTS[] = {
    {50ULL, "0990286324"}, {100ULL, "9228581868"},
    {10000ULL, "1514152075"}, {12000ULL, "9665417773"},
    {15000ULL, "5067751464"}, {20000ULL, "8990743126"},
    {25000ULL, "8674057327"}, {40000ULL, "7469851918"},
    {50000ULL, "8774039057"}, {75000ULL, "3462648045"},
    {100000ULL, "0323336121"}, {120000ULL, "8953301096"},
    {150000ULL, "3757098246"}, {200000ULL, "6656564573"},
    {250000ULL, "4027883484"}, {400000ULL, "1346890248"},
    {500000ULL, "9038045486"}, {750000ULL, "5570453931"},
    {1000000ULL, "3648089791"}, {1200000ULL, "5122997377"},
    {1500000ULL, "3475609603"}, {2000000ULL, "5217087931"},
    {2500000ULL, "6665268558"}, {4000000ULL, "9934557873"},
    {5000000ULL, "6554023086"}, {7500000ULL, "2313396129"},
    {10000000ULL, "8237413280"}, {12000000ULL, "3699127022"},
    {15000000ULL, "1058024895"}, {20000000ULL, "7529885129"},
    {25000000ULL, "0627119064"}, {40000000ULL, "3773193407"},
    {50000000ULL, "6781634645"}, {75000000ULL, "6903021122"},
    {100000000ULL, "6241334147"}, {120000000ULL, "2777622176"},
    {150000000ULL, "5334181771"}, {200000000ULL, "6450042477"},
    {250000000ULL, "9054554181"}, {400000000ULL, "8212955239"},
    {500000000ULL, "4406345766"}, {750000000ULL, "2841568662"},
    {1000000000ULL, "6405298443"}, {1200000000ULL, "5413489027"},
    {1500000000ULL, "1077041491"}, {2000000000ULL, "6535367085"},
    {2500000000ULL, "4068558793"}, {4000000000ULL, "2299235074"},
    {5000000000ULL, "2248748574"}, {7500000000ULL, "5547124594"},
    {10000000000ULL, "6284950888"}, {12000000000ULL, "4655560628"},
    {15000000000ULL, "4594177681"}, {20000000000ULL, "4909466364"},
    {25000000000ULL, "3297267317"}, {40000000000ULL, "9745814812"},
    {50000000000ULL, "8394614199"}, {55000000000ULL, "0113703931"},
    {75000000000ULL, "2065064605"}, {80000000000ULL, "2153808528"},
    {100000000000ULL, "2821098536"}, {120000000000ULL, "8791164665"},
    {125000000000ULL, "9474352718"}, {130000000000ULL, "7591436736"},
    {190000000000ULL, "9568426919"}, {200000000000ULL, "5223459045"},
    {250000000000ULL, "5116943536"}, {300000000000ULL, "6955064161"},
    {310000000000ULL, "0895680415"}, {600000000000ULL, "4549395116"},
};
#define NCHECKPOINTS (sizeof CHECKPOINTS / sizeof CHECKPOINTS[0])

/* log10(16) = 1.2041199826559248 as rational upper/lower bounds
 * (decimal digits per hex digit) */
#define L16_HI 12041199827ULL
#define L16_LO 12041199826ULL
#define L16_DEN 10000000000ULL

#define CHUNK_IO (4u << 20)

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

/* portable u64 -> mpz (unsigned long is 32-bit on some platforms) */
static void mpz_set_u64(mpz_t z, u64 v) {
    char buf[24];
    snprintf(buf, sizeof buf, "%" PRIu64, v);
    mpz_set_str(z, buf, 10);
}

/* ------------------------------------------------------------------ */
/* SHA-256 (certificate file fingerprints; no external dependency)     */
/* ------------------------------------------------------------------ */

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
/* digit file access (y-cruncher .txt digit files, or any digit dump)  */
/* ------------------------------------------------------------------ */

/* per-byte classes: 1 digit, 0 skip (whitespace . , _ :), -1 invalid */
static signed char CLS[2][256];

static void cls_init(void) {
    for (int b = 0; b < 2; b++) {
        for (int i = 0; i < 256; i++)
            CLS[b][i] = -1;
        const char *skip = " \t\r\n\v\f.,_:";
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
    int skip;      /* leading integer digit present in the digit stream */
    int pure;      /* contiguous digits: O(1) seeks possible */
    off_t data0;   /* pure mode: byte offset of first fractional digit */
    u64 n;         /* fractional digit count */
} digitfile;

static size_t strip_head(const unsigned char *head, size_t len, int bi,
                         char *out, size_t cap) {
    size_t o = 0;
    for (size_t i = 0; i < len && o < cap; i++) {
        signed char c = CLS[bi][head[i]];
        if (c == 1)
            out[o++] = (char)tolower(head[i]);
        else if (c < 0)
            break;
    }
    return o;
}

/* Identify base and variant from the file head. Returns 1 on match. */
static int head_match(const char *s, int base, int variant) {
    const char *k = known_frac(base, variant);
    char c0 = INT_DIGIT[variant];
    if (s[0] == c0 && !strncmp(s + 1, k, 12))
        return 1;
    if (!strncmp(s, k, 12))
        return 1;
    return 0;
}

static void detect_file(const char *path, int *base, int *variant) {
    unsigned char head[8192];
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open %s", path);
    size_t got = fread(head, 1, sizeof head, f);
    fclose(f);
    if (!got)
        die("%s is empty", path);
    char s[128];
    size_t sl = strip_head(head, got, 1, s, sizeof s - 1); /* hex superset */
    s[sl] = '\0';
    for (int b = BASE_DEC; b <= BASE_HEX; b++)
        for (int v = VAR_2PV; v <= VAR_PV; v++)
            if ((*base == BASE_AUTO || *base == b) && head_match(s, b, v)) {
                *base = b;
                *variant = v;
                return;
            }
    die("%s does not begin with the lemniscate constant in base 10 or 16\n"
        "       (expected 5.2441151085... / 5.3e7e53e7d... for 2*varpi, or\n"
        "        2.6220575542... / 2.9f3f29f3ea... for varpi; got \"%.24s...\")",
        path, s);
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

/* first=1: the file must open with the constant's leading digits
 * (optionally "5." / "2."); first=0: a raw continuation chunk. */
static void df_open(digitfile *df, const char *path, int base, int variant,
                    int first) {
    memset(df, 0, sizeof *df);
    df->path = path;
    df->base = base;
    int bi = base == BASE_HEX;
    const char *known = known_frac(base, variant);
    char c0 = INT_DIGIT[variant];

    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open %s", path);
    unsigned char head[8192];
    size_t got = fread(head, 1, sizeof head, f);
    if (!got)
        die("%s is empty", path);

    if (first) {
        char s[128];
        size_t sl = strip_head(head, got, bi, s, sizeof s - 1);
        s[sl] = '\0';
        if (s[0] == c0 && !strncmp(s + 1, known, 12))
            df->skip = 1;
        else if (!strncmp(s, known, 12))
            df->skip = 0;
        else
            die("%s does not start with the lemniscate constant (base %s)",
                path, bi ? "16" : "10");
    } else {
        df->skip = 0;
    }

    /* byte-exact prefix, then purity check on the sample */
    size_t i = 0;
    while (i < got && CLS[bi][head[i]] == 0 && head[i] != '.')
        i++;
    off_t data0;
    if (first && i + 1 < got && head[i] == c0 && head[i + 1] == '.')
        data0 = (off_t)i + 2;
    else if (first && df->skip && head[i] == c0)
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

static void df_fetch(digitfile *df, drange *rs, size_t nr) {
    qsort(rs, nr, sizeof *rs, drange_cmp);
    int bi = df->base == BASE_HEX;
    if (df->pure) {
        FILE *f = fopen(df->path, "rb");
        if (!f)
            die("cannot open %s", df->path);
        for (size_t i = 0; i < nr; i++) {
            fseeko(f, df->data0 + (off_t)rs[i].start, SEEK_SET);
            size_t want = (size_t)rs[i].len, done = 0;
            while (done < want) {
                size_t chunk = want - done < CHUNK_IO ? want - done : CHUNK_IO;
                if (fread(rs[i].out + done, 1, chunk, f) != chunk)
                    die("short read at digit %s", fmt_u64(rs[i].start + 1));
                done += chunk;
            }
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
                    rs[j].out[idx - rs[j].start] = (char)tolower(buf[i]);
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
    u64 *cum;
    int nseg;
    u64 n;
} multifile;

static void mf_open(multifile *mf, char **paths, int nf, int base,
                    int variant) {
    mf->segs = malloc((size_t)nf * sizeof *mf->segs);
    mf->cum = malloc(((size_t)nf + 1) * sizeof *mf->cum);
    mf->nseg = nf;
    u64 tot = 0;
    for (int i = 0; i < nf; i++) {
        df_open(&mf->segs[i], paths[i], base, variant, i == 0);
        if (mf->segs[i].n == 0)
            die("%s contains no digits", paths[i]);
        mf->cum[i] = tot;
        tot += mf->segs[i].n;
    }
    mf->cum[nf] = tot;
    mf->n = tot;
}

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

/* ------------------------------------------------------------------ */
/* independent reference digits: L = 2*pi / AGM(1, sqrt 2)  via MPFR   */
/* ------------------------------------------------------------------ */

/* Compute the first nd fractional digits of the constant in `base`.
 * Returns a malloc'd string of exactly nd lowercase digit chars.
 * Truncated (not rounded), matching y-cruncher digit files. */
static char *lem_digits(u64 nd, int base, int variant, int verbose) {
    if (nd < 1 || nd > 2000000000ULL)
        die("reference computation limited to 2e9 digits");
    int extra = 12; /* guard digits for the truncation-boundary check */
    double bits_per = base == BASE_HEX ? 4.0 : 3.3219280948873626;
    mpfr_prec_t prec = (mpfr_prec_t)((double)(nd + (u64)extra + 2) * bits_per)
                       + 96;
    double t0 = now();
    mpfr_t pi, s2, ag, one, L;
    mpfr_inits2(prec, pi, s2, ag, one, L, (mpfr_ptr)0);
    if (verbose)
        fprintf(stderr, "  [mpfr] pi to %s bits...\n", fmt_u64((u64)prec));
    mpfr_const_pi(pi, MPFR_RNDN);
    if (verbose)
        fprintf(stderr, "  [mpfr] pi done (%.1fs); agm(1, sqrt 2)...\n",
                now() - t0);
    mpfr_sqrt_ui(s2, 2, MPFR_RNDN);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_agm(ag, one, s2, MPFR_RNDN);
    mpfr_mul_2ui(L, pi, 1, MPFR_RNDN);
    mpfr_div(L, L, ag, MPFR_RNDN);       /* 2*pi / agm(1, sqrt 2) */
    if (variant == VAR_PV)
        mpfr_div_2ui(L, L, 1, MPFR_RNDN);
    if (verbose)
        fprintf(stderr, "  [mpfr] value ready (%.1fs); extracting digits...\n",
                now() - t0);

    size_t want = (size_t)nd + (size_t)extra + 1; /* + integer digit */
    mpfr_exp_t e;
    char *s = mpfr_get_str(NULL, &e, base == BASE_HEX ? 16 : 10, want, L,
                           MPFR_RNDZ);
    if (!s || e != 1 || s[0] != INT_DIGIT[variant])
        die("internal error: reference value has unexpected leading digits");
    /* if the guard digits are a run of 0s or top-digits, the truncation at
     * nd could sit on a digit boundary within our error bound — never seen
     * in practice, but refuse to guess */
    char topd = base == BASE_HEX ? 'f' : '9';
    int all0 = 1, allt = 1;
    for (size_t i = (size_t)nd + 1; i < want; i++) {
        all0 &= s[i] == '0';
        allt &= s[i] == topd;
    }
    if (all0 || allt)
        die("truncation-boundary ambiguity at digit %s — rerun with a "
            "different --prefix", fmt_u64(nd));
    char *out = malloc((size_t)nd + 1);
    memcpy(out, s + 1, (size_t)nd);
    out[(size_t)nd] = '\0';
    mpfr_free_str(s);
    mpfr_clears(pi, s2, ag, one, L, (mpfr_ptr)0);
    mpfr_free_cache();

    /* self-checks: built-in leading digits + y-cruncher reference table */
    u64 kn = nd < 100 ? nd : 100;
    if (strncmp(out, known_frac(base, variant), kn))
        die("internal error: computed reference digits disagree with the "
            "built-in constant");
    if (base == BASE_DEC && variant == VAR_2PV) {
        int used = 0;
        for (size_t i = 0; i < NCHECKPOINTS && CHECKPOINTS[i].pos <= nd; i++) {
            if (strncmp(out + CHECKPOINTS[i].pos - 10, CHECKPOINTS[i].d, 10))
                die("internal error: computed reference digits disagree with "
                    "the y-cruncher reference table at position %s",
                    fmt_u64(CHECKPOINTS[i].pos));
            used++;
        }
        if (verbose && used)
            fprintf(stderr,
                    "  [mpfr] reference digits cross-validated against %d "
                    "y-cruncher table value(s)\n", used);
    }
    if (verbose)
        fprintf(stderr, "  [mpfr] %s digits in %.1fs\n", fmt_u64(nd),
                now() - t0);
    return out;
}

/* ------------------------------------------------------------------ */
/* options                                                             */
/* ------------------------------------------------------------------ */

typedef struct {
    int base;
    int variant;      /* -1 = auto */
    u64 prefix;       /* verify: independently recomputed digits */
    u64 span;         /* cross: digits loaded from each file */
    int anchors;
    int digits;       /* compared digits per anchor/window */
    int verbose;
    int has_seed;
    u64 seed;
    int line_width;
    const char *out;
    int no_hash;
} opts_t;

static const char *variant_name(int v) {
    return v == VAR_2PV ? "L = 2*varpi = 5.2441151085..."
                        : "varpi = 2.6220575542...";
}

/* ------------------------------------------------------------------ */
/* verify: prefix vs MPFR + reference-table checkpoints                */
/* ------------------------------------------------------------------ */

static int verify_files(char **files, int nf, const opts_t *o) {
    int base = o->base, variant = o->variant;
    if (base == BASE_AUTO || variant < 0)
        detect_file(files[0], &base, &variant);
    int bi = base == BASE_HEX;

    printf("mode: VERIFY — independent prefix recomputation%s\n",
           !bi && variant == VAR_2PV
               ? " + reference-table checkpoints" : "");
    printf("constant: %s (base %s)\n\n", variant_name(variant),
           bi ? "16" : "10");

    printf("[1/4] scanning %d file(s)...\n", nf);
    multifile mf;
    mf_open(&mf, files, nf, base, variant);
    for (int i = 0; i < nf && nf > 1; i++)
        printf("      %s: digits %s..%s (%s access)\n", mf.segs[i].path,
               fmt_u64(mf.cum[i] + 1), fmt_u64(mf.cum[i + 1]),
               mf.segs[i].pure ? "seek" : "streaming");
    printf("      %s fractional digits total%s\n", fmt_u64(mf.n),
           nf == 1 ? (mf.segs[0].pure ? " (seek access)"
                                      : " (streaming access)") : "");
    cert("MODE: verify (independent prefix + checkpoint windows)");
    cert("CONSTANT: %s, base %s", variant_name(variant), bi ? "16" : "10");
    cert("");
    for (int i = 0; i < nf; i++)
        cert("  %s: digits %s..%s", mf.segs[i].path, fmt_u64(mf.cum[i] + 1),
             fmt_u64(mf.cum[i + 1]));
    cert("  total fractional digits: %s", fmt_u64(mf.n));
    cert("");
    cert("CHECKS");

    u64 K = o->prefix;
    if (K > mf.n)
        K = mf.n;
    if (K < 100)
        die("file too small (or --prefix too small): need at least 100 digits");

    int ok = 1;

    printf("[2/4] recomputing the first %s digits independently\n"
           "      (MPFR: 2*pi / AGM(1, sqrt 2))...\n", fmt_u64(K));
    double t0 = now();
    char *ref = lem_digits(K, base, variant, o->verbose);
    printf("      done (%.1fs)\n", now() - t0);

    printf("[3/4] comparing file prefix against the independent value...\n");
    char *filepfx = malloc((size_t)K + 1);
    drange r = {0, K, filepfx, NULL};
    mf_fetch(&mf, &r, 1);
    u64 mism = K;
    for (u64 i = 0; i < K; i++)
        if (filepfx[i] != ref[i]) {
            mism = i;
            break;
        }
    if (mism == K) {
        printf("      PASS: all %s digits match\n", fmt_u64(K));
        cert("  digits 1..%s vs independent MPFR computation "
             "(2*pi/AGM(1,sqrt 2)): PASS", fmt_u64(K));
    } else {
        ok = 0;
        printf("      FAIL: first mismatch at digit %s (file '%c', "
               "computed '%c')\n", fmt_u64(mism + 1), filepfx[mism],
               ref[mism]);
        cert("  digits 1..%s vs independent MPFR computation: FAIL at "
             "digit %s (file '%c', computed '%c')", fmt_u64(K),
             fmt_u64(mism + 1), filepfx[mism], ref[mism]);
    }
    free(filepfx);
    free(ref);

    /* reference-table windows (decimal 2*varpi only) */
    int used = 0, deep = 0;
    if (!bi && variant == VAR_2PV) {
        size_t nc = 0;
        for (size_t i = 0; i < NCHECKPOINTS; i++)
            if (CHECKPOINTS[i].pos <= mf.n)
                nc++;
        printf("[4/4] checking %zu reference-table windows (10 digits each, "
               "up to position %s)...\n", nc,
               nc ? fmt_u64(CHECKPOINTS[nc - 1].pos) : "0");
        drange *rs = calloc(nc, sizeof *rs);
        for (size_t i = 0; i < nc; i++) {
            rs[i].start = CHECKPOINTS[i].pos - 10;
            rs[i].len = 10;
            rs[i].out = malloc(11);
        }
        mf_fetch(&mf, rs, nc);
        /* mf_fetch sorts by position; CHECKPOINTS is sorted too */
        for (size_t i = 0; i < nc; i++) {
            int match = !strncmp(rs[i].out, CHECKPOINTS[i].d, 10);
            ok &= match;
            used++;
            if (CHECKPOINTS[i].pos > K)
                deep++;
            if (!match || o->verbose)
                printf("      position %s: file->%s  table->%s  %s\n",
                       fmt_u64(CHECKPOINTS[i].pos), rs[i].out,
                       CHECKPOINTS[i].d, match ? "PASS" : "FAIL");
            cert("  digits %s..%s vs y-cruncher reference table: file %s, "
                 "table %s: %s", fmt_u64(CHECKPOINTS[i].pos - 9),
                 fmt_u64(CHECKPOINTS[i].pos), rs[i].out, CHECKPOINTS[i].d,
                 match ? "PASS" : "FAIL");
            free(rs[i].out);
        }
        free(rs);
        if (ok)
            printf("      all %d windows match (%d beyond the recomputed "
                   "prefix)\n", used, deep);
    } else {
        printf("[4/4] reference-table windows: n/a for this base/variant\n"
               "      (the published table covers decimal 2*varpi only; use\n"
               "       `cross` to tie this file to a verified decimal file)\n");
    }

    /* tail display: informational, cannot be verified without recompute */
    u64 tn = mf.n < 20 ? mf.n : 20;
    char tailbuf[21];
    drange tr = {mf.n - tn, tn, tailbuf, NULL};
    mf_fetch(&mf, &tr, 1);
    printf("\n      final %s digits (positions %s..%s): %s\n", fmt_u64(tn),
           fmt_u64(mf.n - tn + 1), fmt_u64(mf.n), tailbuf);
    cert("");
    cert("  final digits %s..%s (informational, NOT verified): %s",
         fmt_u64(mf.n - tn + 1), fmt_u64(mf.n), tailbuf);

    printf("\nRESULT: %s\n", ok ? "PASS" : "FAIL");
    if (ok) {
        printf("  Verified: digits 1..%s (independent recomputation)", fmt_u64(K));
        if (used)
            printf(" and %d\n  10-digit windows at published positions up "
                   "to %s", used, fmt_u64(CHECKPOINTS[0].pos <= mf.n
                       ? CHECKPOINTS[(size_t)used - 1].pos : 0));
        printf(".\n"
           "  NOT verified: everything else. The lemniscate constant has no\n"
           "  digit-extraction (BBP-type) formula, so unlike pi there is no\n"
           "  independent spot check at arbitrary positions. Digits beyond the\n"
           "  covered points can only be verified by a second full computation\n"
           "  with a different formula (then: lemverify compare FILE1 FILE2),\n"
           "  and dec/hex mutual consistency by `lemverify cross`.\n");
    }
    cert("");
    cert("INTERPRETATION");
    cert("  The first %s digits were recomputed with MPFR (arbitrary-"
         "precision,", fmt_u64(K));
    cert("  correctly rounded) via 2*pi/AGM(1, sqrt 2) — an implementation "
         "and");
    cert("  formula independent of y-cruncher's — and matched digit for "
         "digit.");
    if (used) {
        cert("  %d ten-digit windows at published reference positions "
             "(sourced from", used);
        cert("  independent record computations, up to position "
             "600,000,000,000) were");
        cert("  compared by direct file reads.");
    }
    cert("  No digit-extraction formula exists for this constant, so "
         "positions");
    cert("  outside the recomputed prefix and the sampled windows are NOT");
    cert("  covered. Full verification requires an independent recomputation");
    cert("  (different formula) compared with `lemverify compare`.");
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* cross: mutual dec<->hex base-conversion anchors                     */
/* ------------------------------------------------------------------ */

/* load the first `count` fractional digits of a file series into an mpz
 * (as an integer, most significant digit first) */
static void load_prefix_mpz(multifile *mf, u64 count, int base, mpz_t out,
                            const char *what) {
    char *buf = malloc((size_t)count + 1);
    if (!buf)
        die("out of memory loading %s prefix (%s digits)", what,
            fmt_u64(count));
    drange r = {0, count, buf, NULL};
    mf_fetch(mf, &r, 1);
    double t0 = now();
    if (mpz_set_str(out, buf, base == BASE_HEX ? 16 : 10))
        die("GMP could not parse the %s digit prefix", what);
    free(buf);
    (void)t0;
}

static int cross_mode(char **files, int nf, const opts_t *o) {
    if (nf != 2)
        die("cross needs exactly two files: DEC_FILE HEX_FILE (either order)");
    int t = o->digits ? o->digits : 12;
    if (t > 16)
        die("--digits must be 1..16");

    /* identify which file is which; variants must agree */
    int b0 = BASE_AUTO, v0 = -1, b1 = BASE_AUTO, v1 = -1;
    detect_file(files[0], &b0, &v0);
    detect_file(files[1], &b1, &v1);
    if (b0 == b1)
        die("cross needs one decimal and one hexadecimal file "
            "(both are base %s)", b0 == BASE_HEX ? "16" : "10");
    if (v0 != v1)
        die("the files hold different constants (%s vs %s)",
            variant_name(v0), variant_name(v1));
    const char *decpath = b0 == BASE_DEC ? files[0] : files[1];
    const char *hexpath = b0 == BASE_DEC ? files[1] : files[0];
    int variant = v0;

    printf("mode: CROSS — mutual dec<->hex base-conversion verification\n");
    printf("constant: %s\n\n", variant_name(variant));

    printf("[1/5] scanning files...\n");
    multifile mfd, mfh;
    char *dp = (char *)decpath, *hp = (char *)hexpath;
    mf_open(&mfd, &dp, 1, BASE_DEC, variant);
    mf_open(&mfh, &hp, 1, BASE_HEX, variant);
    printf("      dec: %s — %s digits\n", decpath, fmt_u64(mfd.n));
    printf("      hex: %s — %s digits\n", hexpath, fmt_u64(mfh.n));

    /* choose matched spans: Nd dec digits ~ Nh hex digits of information */
    u64 Nd = o->span;
    if (Nd > mfd.n)
        Nd = mfd.n;
    u64 Nh = (u64)((u128)Nd * L16_DEN / L16_HI); /* dec -> hex equivalent */
    if (Nh > mfh.n) {
        Nh = mfh.n;
        u64 nd2 = (u64)((u128)Nh * L16_LO / L16_DEN);
        if (nd2 < Nd)
            Nd = nd2;
    }
    if (Nd < 2000 || Nh < 2000)
        die("need at least 2000 digits in both files for cross verification");
    printf("      spans: first %s dec digits <-> first %s hex digits\n",
           fmt_u64(Nd), fmt_u64(Nh));

    cert("MODE: cross (mutual dec<->hex base-conversion anchors, %d digits "
         "compared per anchor)", t);
    cert("CONSTANT: %s", variant_name(variant));
    cert("");
    cert("  dec: %s (%s fractional digits, %s used)", decpath,
         fmt_u64(mfd.n), fmt_u64(Nd));
    cert("  hex: %s (%s fractional digits, %s used)", hexpath,
         fmt_u64(mfh.n), fmt_u64(Nh));
    cert("");
    cert("CHECKS");

    int na = o->anchors > 0 ? o->anchors : 3;
    int ok = 1;
    double t0;

    printf("[2/5] loading digit prefixes into big integers (GMP)...\n");
    t0 = now();
    mpz_t D, H;
    mpz_inits(D, H, NULL);
    load_prefix_mpz(&mfd, Nd, BASE_DEC, D, "decimal");
    load_prefix_mpz(&mfh, Nh, BASE_HEX, H, "hexadecimal");
    printf("      done (%.1fs)\n", now() - t0);

    /* ---- direction 1: dec -> hex ----------------------------------- */
    /* frac(16^m * (a + D/10^Nd)) = (2^(4m) mod 10^Nd) * D mod 10^Nd / 10^Nd
     * The top t hex digits must equal hex-file digits m+1..m+t.
     * Sensitivity: covers dec digits 1..~1.204*(m+t); anchor cap keeps the
     * file-truncation error 20 hex digits below the compared window. */
    u64 m_cap = (u64)((u128)Nd * L16_DEN / L16_HI) - (u64)t - 20;
    if (m_cap > Nh - (u64)t)
        m_cap = Nh - (u64)t;
    u64 *am = malloc((size_t)na * sizeof *am);
    am[0] = m_cap;
    for (int i = 1; i < na; i++) {
        u64 m;
        int dup;
        do {
            m = m_cap - 1 - rng_below(511);
            dup = 0;
            for (int j = 0; j < i; j++)
                dup |= am[j] == m;
        } while (dup);
        am[i] = m;
    }
    /* one powmod at the smallest anchor, then cheap shifts upward */
    u64 m_min = am[0];
    for (int i = 1; i < na; i++)
        if (am[i] < m_min)
            m_min = am[i];

    printf("[3/5] dec->hex: %d anchor(s) near hex position %s\n"
           "      (one 2^(4m) mod 10^N powmod, then shifts)...\n", na,
           fmt_u64(m_cap));
    t0 = now();
    mpz_t mod10, P, y, q, e;
    mpz_inits(mod10, P, y, q, e, NULL);
    mpz_ui_pow_ui(mod10, 10, (unsigned long)Nd);
    mpz_set_ui(P, 2);
    mpz_set_u64(e, 4 * m_min);
    mpz_powm(P, P, e, mod10);                /* 2^(4*m_min) mod 10^Nd */
    printf("      powmod done (%.1fs)\n", now() - t0);

    /* process anchors in increasing m so shifting only goes up */
    for (int pass = 0; pass < na; pass++) {
        /* next smallest unprocessed anchor */
        int bi_ = -1;
        for (int i = 0; i < na; i++)
            if (am[i] != UINT64_MAX && (bi_ < 0 || am[i] < am[bi_]))
                bi_ = i;
        u64 m = am[bi_];
        am[bi_] = UINT64_MAX;
        if (m > m_min) {
            mpz_mul_2exp(P, P, (mp_bitcnt_t)(4 * (m - m_min)));
            mpz_mod(P, P, mod10);
            m_min = m;
        }
        t0 = now();
        mpz_mul(y, P, D);
        mpz_mod(y, y, mod10);                /* frac numerator */
        mpz_mul_2exp(q, y, (mp_bitcnt_t)(4 * (unsigned)t));
        mpz_tdiv_q(q, q, mod10);             /* top t hex digits */
        char raw[20], conv[20];
        mpz_get_str(raw, 16, q);
        size_t rl = strlen(raw);
        memset(conv, '0', (size_t)t);
        memcpy(conv + ((size_t)t - rl), raw, rl + 1);
        conv[t] = '\0';
        char filebuf[20];
        drange fr = {m, (u64)t, filebuf, NULL};
        mf_fetch(&mfh, &fr, 1);
        int match = !strcmp(conv, filebuf);
        ok &= match;
        printf("      hex position %s..%s: hexfile->%s  from-dec->%s  %s "
               "(%.1fs)\n", fmt_u64(m + 1), fmt_u64(m + (u64)t), filebuf,
               conv, match ? "PASS" : "FAIL", now() - t0);
        cert("  dec->hex anchor, hex digits %s..%s: hex file %s, computed "
             "from dec file %s: %s", fmt_u64(m + 1), fmt_u64(m + (u64)t),
             filebuf, conv, match ? "PASS" : "FAIL");
        fflush(stdout);
    }
    u64 dec_covered = (u64)((u128)(m_cap + (u64)t) * L16_LO / L16_DEN);
    if (dec_covered > Nd)
        dec_covered = Nd;
    mpz_clears(P, y, q, NULL);

    /* ---- direction 2: hex -> dec ----------------------------------- */
    /* frac(10^r * (a + H/2^(4Nh))) -> top t dec digits must equal dec-file
     * digits r+1..r+t. Modulus 2^(4Nh): reductions are bit masks. */
    u64 r_cap = (u64)((u128)Nh * L16_LO / L16_DEN) - (u64)t - 24;
    if (r_cap > Nd - (u64)t)
        r_cap = Nd - (u64)t;
    u64 *ar = malloc((size_t)na * sizeof *ar);
    ar[0] = r_cap;
    for (int i = 1; i < na; i++) {
        u64 r_;
        int dup;
        do {
            r_ = r_cap - 1 - rng_below(511);
            dup = 0;
            for (int j = 0; j < i; j++)
                dup |= ar[j] == r_;
        } while (dup);
        ar[i] = r_;
    }
    u64 r_min = ar[0];
    for (int i = 1; i < na; i++)
        if (ar[i] < r_min)
            r_min = ar[i];

    printf("[4/5] hex->dec: %d anchor(s) near dec position %s\n"
           "      (one 10^r mod 2^(4N) powmod, then small multiplies)...\n",
           na, fmt_u64(r_cap));
    t0 = now();
    mp_bitcnt_t hbits = (mp_bitcnt_t)(4 * Nh);
    mpz_t P2, y2, q2, mod16;
    mpz_inits(P2, y2, q2, mod16, NULL);
    mpz_set_ui(mod16, 1);
    mpz_mul_2exp(mod16, mod16, hbits);       /* 2^(4Nh) */
    mpz_set_ui(P2, 10);
    mpz_set_u64(e, r_min);
    mpz_powm(P2, P2, e, mod16);              /* 10^r_min mod 2^(4Nh) */
    printf("      powmod done (%.1fs)\n", now() - t0);

    for (int pass = 0; pass < na; pass++) {
        int bi_ = -1;
        for (int i = 0; i < na; i++)
            if (ar[i] != UINT64_MAX && (bi_ < 0 || ar[i] < ar[bi_]))
                bi_ = i;
        u64 r_ = ar[bi_];
        ar[bi_] = UINT64_MAX;
        if (r_ > r_min) {
            mpz_ui_pow_ui(y2, 10, (unsigned long)(r_ - r_min));
            mpz_mul(P2, P2, y2);
            mpz_fdiv_r_2exp(P2, P2, hbits);
            r_min = r_;
        }
        t0 = now();
        mpz_mul(y2, P2, H);
        mpz_fdiv_r_2exp(y2, y2, hbits);      /* frac numerator */
        mpz_ui_pow_ui(q2, 10, (unsigned)t);
        mpz_mul(q2, q2, y2);
        mpz_fdiv_q_2exp(q2, q2, hbits);      /* top t dec digits */
        char raw[20], conv[20];
        mpz_get_str(raw, 10, q2);
        size_t rl = strlen(raw);
        memset(conv, '0', (size_t)t);
        memcpy(conv + ((size_t)t - rl), raw, rl + 1);
        conv[t] = '\0';
        char filebuf[20];
        drange fr = {r_, (u64)t, filebuf, NULL};
        mf_fetch(&mfd, &fr, 1);
        int match = !strcmp(conv, filebuf);
        ok &= match;
        printf("      dec position %s..%s: decfile->%s  from-hex->%s  %s "
               "(%.1fs)\n", fmt_u64(r_ + 1), fmt_u64(r_ + (u64)t), filebuf,
               conv, match ? "PASS" : "FAIL", now() - t0);
        cert("  hex->dec anchor, dec digits %s..%s: dec file %s, computed "
             "from hex file %s: %s", fmt_u64(r_ + 1), fmt_u64(r_ + (u64)t),
             filebuf, conv, match ? "PASS" : "FAIL");
        fflush(stdout);
    }
    u64 hex_covered = (u64)((u128)(r_cap + (u64)t) * L16_DEN / L16_HI);
    if (hex_covered > Nh)
        hex_covered = Nh;
    mpz_clears(D, H, mod10, mod16, P2, y2, q2, e, NULL);
    free(am);
    free(ar);

    printf("[5/5] summary\n");
    printf("\nRESULT: %s\n", ok ? "PASS" : "FAIL");
    if (ok) {
        printf(
        "  The two files agree as the same real number over the checked "
        "spans:\n"
        "    dec digits 1..%s verified against the hex file (one wrong "
        "digit\n"
        "      anywhere in that range scrambles every dec->hex anchor),\n"
        "    hex digits 1..%s verified against the dec file likewise.\n"
        "  Miss probability per anchor ~ %d^-%d; %d anchors per direction.\n"
        "  This checks conversion/storage integrity of both files over the\n"
        "  span; it does NOT verify y-cruncher's computation itself (both\n"
        "  files derive from the same binary value). Raise --span (RAM\n"
        "  permitting) to cover more digits.\n",
        fmt_u64(dec_covered), fmt_u64(hex_covered), 16, t, na);
    } else {
        printf("  The dec and hex files DISAGREE within the checked spans —\n"
               "  at least one of them is corrupted (or they are from "
               "different runs).\n");
    }
    cert("");
    cert("INTERPRETATION");
    cert("  dec->hex: frac(16^m * x_dec) is computed from ALL %s loaded "
         "decimal", fmt_u64(Nd));
    cert("  digits by modular arithmetic and compared with digits read "
         "directly");
    cert("  from the hex file; base conversion mixes every source digit "
         "into the");
    cert("  anchor value. Covered: dec digits 1..%s.", fmt_u64(dec_covered));
    cert("  hex->dec: the mirror check. Covered: hex digits 1..%s.",
         fmt_u64(hex_covered));
    cert("  A PASS means both files encode the same real number over these");
    cert("  prefixes (miss probability ~16^-%d per anchor, %d anchors per");
    cert("  direction). It does not verify the underlying computation.", t, na);
    return ok ? 0 : 1;
}

/* ------------------------------------------------------------------ */
/* compare: stream two digit files against each other                  */
/* ------------------------------------------------------------------ */

static int compare_mode(char **files, int nf, const opts_t *o) {
    if (nf != 2)
        die("compare needs exactly two files");
    int base = o->base, variant = o->variant;
    if (base == BASE_AUTO || variant < 0)
        detect_file(files[0], &base, &variant);
    {
        int b2 = BASE_AUTO, v2 = -1;
        detect_file(files[1], &b2, &v2);
        if (b2 != base || v2 != variant)
            die("files differ in base or constant (%s base %s vs %s base %s)",
                variant_name(variant), base == BASE_HEX ? "16" : "10",
                variant_name(v2), b2 == BASE_HEX ? "16" : "10");
    }
    int bi = base == BASE_HEX;
    printf("mode: COMPARE — streaming digit-for-digit comparison\n");
    printf("constant: %s (base %s)\n\n", variant_name(variant),
           bi ? "16" : "10");
    printf("  a: %s\n  b: %s\n\n", files[0], files[1]);

    FILE *fa = fopen(files[0], "rb"), *fb = fopen(files[1], "rb");
    if (!fa || !fb)
        die("cannot open input files");
    unsigned char *ba = malloc(CHUNK_IO), *bb = malloc(CHUNK_IO);
    char *da = malloc(CHUNK_IO), *db = malloc(CHUNK_IO);
    size_t la = 0, lb = 0;   /* digits buffered */
    int skipa = -1, skipb = -1; /* leading integer digit not yet dropped */
    u64 pos = 0;             /* fractional digits compared */
    u64 mismatch_at = 0;
    char mma = 0, mmb = 0;
    int done_a = 0, done_b = 0, mismatched = 0;
    double t0 = now(), tlast = t0;

    while (!mismatched) {
        /* refill a */
        while (la == 0 && !done_a) {
            size_t got = fread(ba, 1, CHUNK_IO, fa);
            if (!got) {
                done_a = 1;
                break;
            }
            for (size_t i = 0; i < got; i++) {
                signed char c = CLS[bi][ba[i]];
                if (c < 0)
                    die("unexpected character 0x%02x in %s", ba[i], files[0]);
                if (c == 1) {
                    if (skipa) {
                        skipa = 0; /* first digit = integer part */
                        continue;
                    }
                    da[la++] = (char)tolower(ba[i]);
                }
            }
        }
        while (lb == 0 && !done_b) {
            size_t got = fread(bb, 1, CHUNK_IO, fb);
            if (!got) {
                done_b = 1;
                break;
            }
            for (size_t i = 0; i < got; i++) {
                signed char c = CLS[bi][bb[i]];
                if (c < 0)
                    die("unexpected character 0x%02x in %s", bb[i], files[1]);
                if (c == 1) {
                    if (skipb) {
                        skipb = 0;
                        continue;
                    }
                    db[lb++] = (char)tolower(bb[i]);
                }
            }
        }
        if (la == 0 || lb == 0)
            break;
        size_t n = la < lb ? la : lb;
        if (memcmp(da, db, n)) {
            for (size_t i = 0; i < n; i++)
                if (da[i] != db[i]) {
                    mismatched = 1;
                    mismatch_at = pos + i + 1;
                    mma = da[i];
                    mmb = db[i];
                    break;
                }
        }
        pos += n;
        memmove(da, da + n, la - n);
        memmove(db, db + n, lb - n);
        la -= n;
        lb -= n;
        if (now() - tlast > 5.0) {
            fprintf(stderr, "  ... %s digits compared (%.0f Mdigit/s)\r",
                    fmt_u64(pos), (double)pos / (now() - t0) / 1e6);
            tlast = now();
        }
    }
    fprintf(stderr, "%-60s\r", "");
    free(ba); free(bb); free(da); free(db);
    fclose(fa); fclose(fb);

    cert("MODE: compare (streaming digit-for-digit)");
    cert("  a: %s", files[0]);
    cert("  b: %s", files[1]);
    if (mismatched) {
        printf("RESULT: FAIL — first mismatch at fractional digit %s "
               "('%c' vs '%c')\n", fmt_u64(mismatch_at), mma, mmb);
        printf("  %s digits agreed before the mismatch.\n",
               fmt_u64(mismatch_at - 1));
        cert("  FAIL: first mismatch at digit %s ('%c' vs '%c')",
             fmt_u64(mismatch_at), mma, mmb);
        return 1;
    }
    printf("RESULT: PASS — %s fractional digits identical (%.1fs, "
           "%.0f Mdigit/s)\n", fmt_u64(pos), now() - t0,
           (double)pos / (now() - t0) / 1e6);
    if (la || lb || !done_a || !done_b)
        printf("  note: %s is longer; extra digits beyond %s were not "
               "compared.\n", (la || !done_a) ? files[0] : files[1],
               fmt_u64(pos));
    printf("  If the two files come from independent computations with\n"
           "  different formulas, this IS the standard full verification of\n"
           "  every compared digit.\n");
    cert("  PASS: %s digits identical", fmt_u64(pos));
    return 0;
}

/* ------------------------------------------------------------------ */
/* scan: whole-file streaming integrity + statistics                   */
/* ------------------------------------------------------------------ */

static int scan_mode(char **files, int nf, const opts_t *o) {
    int base = o->base, variant = o->variant;
    if (base == BASE_AUTO || variant < 0)
        detect_file(files[0], &base, &variant);
    int bi = base == BASE_HEX;
    int nb = bi ? 16 : 10;
    printf("mode: SCAN — whole-file streaming integrity/statistics\n");
    printf("constant: %s (base %s)\n\n", variant_name(variant),
           bi ? "16" : "10");
    u64 counts[16] = {0};
    u64 total = 0;
    double t0 = now(), tlast = t0;
    for (int fi = 0; fi < nf; fi++) {
        FILE *f = fopen(files[fi], "rb");
        if (!f)
            die("cannot open %s", files[fi]);
        unsigned char *buf = malloc(CHUNK_IO);
        sha256_t sh;
        sha256_init(&sh);
        u64 fdigits = 0, fbytes = 0;
        size_t got;
        while ((got = fread(buf, 1, CHUNK_IO, f)) > 0) {
            if (!o->no_hash)
                sha256_update(&sh, buf, got);
            for (size_t i = 0; i < got; i++) {
                signed char c = CLS[bi][buf[i]];
                if (c < 0)
                    die("unexpected character 0x%02x at byte %s of %s",
                        buf[i], fmt_u64(fbytes + i), files[fi]);
                if (c == 1) {
                    unsigned char ch = buf[i];
                    int d = ch <= '9' ? ch - '0'
                                      : (int)(tolower(ch) - 'a') + 10;
                    counts[d]++;
                    fdigits++;
                }
            }
            fbytes += got;
            if (now() - tlast > 5.0) {
                fprintf(stderr, "  ... %s digits scanned\r",
                        fmt_u64(total + fdigits));
                tlast = now();
            }
        }
        free(buf);
        fclose(f);
        fprintf(stderr, "%-50s\r", "");
        printf("  %s: %s digits", files[fi], fmt_u64(fdigits));
        if (!o->no_hash) {
            char hex[65];
            sha256_hex(&sh, hex);
            printf("\n      sha-256: %s", hex);
        }
        printf("\n");
        total += fdigits;
        if (fi == 0 && fdigits)
            counts[INT_DIGIT[variant] - '0']--; /* drop the integer digit */
    }
    total--; /* integer digit */
    printf("\n  total fractional digits: %s (%.1fs)\n", fmt_u64(total),
           now() - t0);
    double exp_ = (double)total / nb, chi2 = 0;
    printf("  digit frequencies:\n");
    for (int d = 0; d < nb; d++) {
        double dev = (double)counts[d] - exp_;
        chi2 += dev * dev / exp_;
        printf("    %x: %s (%+.6f%%)\n", d, fmt_u64(counts[d]),
               100.0 * dev / exp_);
    }
    printf("  chi-square (%d dof): %.2f  (expected ~%d +- %.1f for random "
           "digits)\n", nb - 1, chi2, nb - 1, sqrt(2.0 * (nb - 1)));
    printf("\n  A clean scan proves the files are well-formed and gives "
           "their\n  fingerprints; it does not verify digit values.\n");
    return 0;
}

/* ------------------------------------------------------------------ */
/* gen / digits                                                        */
/* ------------------------------------------------------------------ */

static int cmd_gen(u64 n, const opts_t *o) {
    if (n < 100)
        die("gen needs at least 100 digits");
    if (!o->out)
        die("gen requires --out FILE");
    int base = o->base == BASE_AUTO ? BASE_DEC : o->base;
    int variant = o->variant < 0 ? VAR_2PV : o->variant;
    double t0 = now();
    char *s = lem_digits(n, base, variant, o->verbose);
    FILE *f = fopen(o->out, "w");
    if (!f)
        die("cannot write %s", o->out);
    fprintf(f, "%c.", INT_DIGIT[variant]);
    if (o->line_width <= 0)
        fputs(s, f);
    else {
        fputc('\n', f);
        for (u64 i = 0; i < n; i += (u64)o->line_width) {
            u64 len = n - i < (u64)o->line_width ? n - i : (u64)o->line_width;
            fwrite(s + i, 1, len, f);
            fputc('\n', f);
        }
    }
    fclose(f);
    printf("wrote %s %s digits of %s to %s (%.1fs)\n", fmt_u64(n),
           base == BASE_HEX ? "hex" : "dec", variant_name(variant), o->out,
           now() - t0);
    free(s);
    return 0;
}

static int cmd_digits(u64 pos, const opts_t *o) {
    int t = o->digits ? o->digits : 10;
    int base = o->base == BASE_AUTO ? BASE_DEC : o->base;
    int variant = o->variant < 0 ? VAR_2PV : o->variant;
    u64 need = pos + (u64)t - 1;
    if (need > 200000000ULL)
        fprintf(stderr,
                "note: no digit-extraction formula exists for this constant; "
                "computing\nthe whole prefix (%s digits) — this may take a "
                "long time...\n", fmt_u64(need));
    double t0 = now();
    char *s = lem_digits(need, base, variant, o->verbose);
    printf("%s digits %s..%s of %s (after the point): %.*s   [%.1fs]\n",
           base == BASE_HEX ? "hex" : "dec", fmt_u64(pos),
           fmt_u64(pos + (u64)t - 1), variant_name(variant), t,
           s + pos - 1, now() - t0);
    free(s);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(void) {
    fprintf(stderr,
        "usage: lemverify verify FILE [FILE...] [--prefix K] [--base auto|dec|hex]\n"
        "                        [--constant auto|2pv|pv] [-o CERT.txt] [--no-hash] [-v]\n"
        "         independently recompute the first K digits (default 1,000,000; MPFR\n"
        "         2*pi/AGM(1,sqrt 2)) + check y-cruncher reference-table windows at up\n"
        "         to 64 positions reaching 600 billion (decimal files).\n"
        "         Multiple FILEs are read in series as one digit stream.\n"
        "       lemverify cross DEC_FILE HEX_FILE [--span N] [--anchors A]\n"
        "                        [--digits t<=16] [--seed S] [-o CERT.txt] [-v]\n"
        "         mutual dec<->hex base-conversion verification over the first N\n"
        "         dec digits (default 250,000,000; RAM- and CPU-bound).\n"
        "       lemverify compare FILE_A FILE_B [-o CERT.txt]\n"
        "         stream-compare two runs (use different formulas for a real proof).\n"
        "       lemverify scan FILE [FILE...] [--no-hash]\n"
        "         whole-file integrity: counts, frequencies, sha-256.\n"
        "       lemverify gen N --out FILE [--base dec|hex] [--constant 2pv|pv]\n"
        "                        [--line-width W]\n"
        "       lemverify digits POS [--count t<=16] [--base dec|hex]\n"
        "\n"
        "constant: 2pv = y-cruncher \"Lemniscate\" = 5.2441151085... (default),\n"
        "          pv  = varpi = 2.6220575542...; auto-detected from files.\n");
    exit(2);
}

int main(int argc, char **argv) {
    cls_init();
    if (argc < 3)
        usage();
    const char *cmd = argv[1];

    opts_t o = {BASE_AUTO, -1, 1000000ULL, 250000000ULL, 3, 0, 0, 0, 0, 0,
                NULL, 0};
    static struct option lo[] = {
        {"base", required_argument, 0, 'b'},
        {"constant", required_argument, 0, 'C'},
        {"prefix", required_argument, 0, 'p'},
        {"span", required_argument, 0, 'n'},
        {"anchors", required_argument, 0, 'a'},
        {"digits", required_argument, 0, 'd'},
        {"count", required_argument, 0, 'd'},
        {"seed", required_argument, 0, 's'},
        {"out", required_argument, 0, 'o'},
        {"line-width", required_argument, 0, 'w'},
        {"no-hash", no_argument, 0, 'H'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}};
    optind = 2;
    int c;
    while ((c = getopt_long(argc, argv, "b:C:p:n:a:d:s:o:w:Hv", lo, NULL))
           != -1)
        switch (c) {
        case 'b':
            if (!strcmp(optarg, "dec")) o.base = BASE_DEC;
            else if (!strcmp(optarg, "hex")) o.base = BASE_HEX;
            else if (!strcmp(optarg, "auto")) o.base = BASE_AUTO;
            else die("bad --base %s", optarg);
            break;
        case 'C':
            if (!strcmp(optarg, "2pv")) o.variant = VAR_2PV;
            else if (!strcmp(optarg, "pv")) o.variant = VAR_PV;
            else if (!strcmp(optarg, "auto")) o.variant = -1;
            else die("bad --constant %s (use 2pv, pv, or auto)", optarg);
            break;
        case 'p': o.prefix = strtoull(optarg, NULL, 10); break;
        case 'n': o.span = strtoull(optarg, NULL, 10); break;
        case 'a': o.anchors = atoi(optarg); break;
        case 'd': o.digits = atoi(optarg); break;
        case 's': o.has_seed = 1; o.seed = strtoull(optarg, NULL, 10); break;
        case 'o': o.out = optarg; break;
        case 'w': o.line_width = atoi(optarg); break;
        case 'H': o.no_hash = 1; break;
        case 'v': o.verbose = 1; break;
        default: usage();
        }
    if (optind >= argc)
        usage();
    if (o.digits < 0 || o.digits > 16)
        die("--digits/--count must be 1..16");
    if (o.span > 20000000000ULL)
        die("--span above 20e9 digits exceeds GMP/RAM practicality");
    if (sizeof(unsigned long) < 8 &&
        (o.span > 500000000ULL || o.prefix > 500000000ULL))
        die("this platform has 32-bit unsigned long (native Windows?): "
            "keep --span/--prefix below 5e8, or build under WSL");
    rng_state = o.has_seed
                    ? o.seed
                    : (u64)time(NULL) ^ ((u64)getpid() << 32) ^ (u64)clock();

    int is_verify = !strcmp(cmd, "verify");
    int is_cross = !strcmp(cmd, "cross");
    int is_compare = !strcmp(cmd, "compare");
    int is_scan = !strcmp(cmd, "scan");
    if (is_verify || is_cross || is_compare || is_scan) {
        char **files = &argv[optind];
        int nf = argc - optind;
        u64 seed_used = rng_state;
        printf("lemverify %s: %s%s\n", cmd, files[0],
               nf > 1 ? " (+ more)" : "");
        printf("engine: C | GMP %s | MPFR %s\n\n", gmp_version,
               mpfr_get_version());
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
            cert("             LEMVERIFY VERIFICATION CERTIFICATE");
            cert("================================================================");
            cert("tool:        lemverify %s (lemniscate-constant digit"
                 " verification)", LEMVERIFY_VERSION);
            cert("date (UTC):  %s", ts);
            cert("host:        %s", host);
            fputs("command:     ", g_cert);
            for (int i = 0; i < argc; i++)
                fprintf(g_cert, "%s%s", i ? " " : "", argv[i]);
            fputc('\n', g_cert);
            cert("libs:        GMP %s | MPFR %s", gmp_version,
                 mpfr_get_version());
            cert("seed:        %" PRIu64 "%s", seed_used,
                 o.has_seed ? "" : "  (auto; rerun with --seed to reproduce)");
            cert("");
            cert("INPUT FILES (%d)", nf);
            if (!o.no_hash && !is_scan)
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
                if (!o.no_hash && !is_scan) {
                    char hex[65];
                    sha256_file(files[i], hex);
                    cert("     sha-256: %s", hex);
                }
            }
            cert("");
        }
        int rc;
        if (is_verify)
            rc = verify_files(files, nf, &o);
        else if (is_cross)
            rc = cross_mode(files, nf, &o);
        else if (is_compare)
            rc = compare_mode(files, nf, &o);
        else
            rc = scan_mode(files, nf, &o);
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
    const char *arg = argv[optind];
    if (!strcmp(cmd, "gen")) {
        char *end;
        unsigned long long n = strtoull(arg, &end, 10);
        if (*end || n < 1)
            die("bad digit count %s", arg);
        return cmd_gen((u64)n, &o);
    }
    if (!strcmp(cmd, "digits")) {
        char *end;
        unsigned long long p = strtoull(arg, &end, 10);
        if (*end || p < 1)
            die("bad position %s", arg);
        return cmd_digits((u64)p, &o);
    }
    usage();
    return 2;
}
