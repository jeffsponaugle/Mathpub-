/*
 * piself.c - find self-locating numbers in the digits of pi
 *
 * A number k is "self-locating" when the decimal digits of k appear in pi
 * starting at the position that k itself names.  The literature counts
 * positions several ways; all reduce to one rule: the digit window starting
 * at 1-based after-decimal position p equals p + d for a fixed offset d.
 *
 *   d = -1  0-based: position 0 is the first digit after the decimal point
 *           (OEIS A057680: 6, 27, 13598, 43611, 24643510, ...)
 *   d =  0  1-based: position 1 is the first digit after the decimal point
 *           (OEIS A064810: 1, 16470, 44899, 79873884, ...)
 *   d =  1  3-first: the leading 3 itself is digit 1, so "5" at digit 5
 *           (OEIS A057679: 5, 242424, 271070, 9292071, ...)
 *   d =  2  "k-2": the 3 and the decimal point both count as characters
 *           (OEIS A153220: 4, 315, 360, 384, 47696, ...)
 *   d =  3  "k-3": k found at after-decimal position k-3
 *           (OEIS A153221: 51, 875, 62843, 242424, ...)
 *   d =  4  "k-4" (OEIS A153223: 9, 233, 1614, 9218, ...)
 *   d =  5  "k-5" (OEIS A153224: 26, 32, 41, 86, 2799, ...)
 *
 * Any other offset can be requested with -d N (negative N allowed).
 *
 * The multiplicative variants place k at a multiple of itself: the window
 * at 1-based after-decimal position m*k equals k.
 *   m = 2  "2k" (OEIS A153227: 5, 95, 171, 529, 1913, 2753, ...)
 *   m = 3  "3k" (OEIS A153228: 1, 2, 3, 85, 200, 447263, ...)
 * Other multipliers can be requested with -x N.
 *
 * The reversed variants (-r) place the digit-reversal of k at the position
 * k names ("0161" starting at position 1610):
 *   d = 0  "rev-1-based" (OEIS A366831: 1, 50, 576, 242424, 746074, ...)
 *   d = 1  "rev-3-first" (OEIS A366830: 5, 1610, 5833, 82699856, ...)
 * Any other offset can be requested with -r N.
 *
 * The residue variants (-m) only require k mod M at position k:
 *   M = 100   "mod-100"  (OEIS A153225: 1, 102, 104, 189, 193, ...)
 *   M = 1000  "mod-1000" (OEIS A153226: 1, 1005, 1053, 1255, ...)
 * These are dense (~1.9% of positions hit for M=100, ~0.28% for M=1000),
 * so they are opt-in and not part of the default scan.
 *
 * "Pithy numbers" (--pithy) are the non-overlapping variant: chop pi into
 * consecutive m-tuples of digits; an m-digit k is pithy when the k-th
 * m-tuple is exactly k.  Again two conventions:
 *   pithy:   tuples start after the decimal point
 *            (OEIS A109513: 1, 19, 94, 3542, 7295, 318320, ...)
 *   pithy-3: tuples start at the leading 3
 *            (OEIS A109514: 5, 9696, 19781, 199898, ...)
 *
 * Usage:  piself [-p] <pi-file> [-n digits] [-1] [-0] [-3] [-d N] [--pithy]
 */
#include "pi_common.h"          /* digit sources (text/ycd), POW10, timing and formatting helpers */

/* Overlap window: the longest number we ever match is ndigits(position),
 * 20 digits covers any file that fits on a disk. */
#define SELF_L 20
#define MAX_OFF 32

/*
 * usage - print the command-line help to stderr and exit(1).
 */
static void usage(const char *argv0)
{
    fprintf(stderr,
        "piself v" PI_TOOLS_VERSION " - self-locating numbers in pi\n"
        "Usage: %s [-p] <pi-file> [options]\n"
        "\n"
        "Scan pi for self-locating numbers: numbers k whose decimal digits\n"
        "appear in pi starting at the position k names.\n"
        "\n"
        "Options:\n"
        "  -p FILE   pi digit source (may also be given as a bare argument)\n"
        "  -n N      scan only the first N digits (500M, 1B, all; default all)\n"
        "  -o FILE   also write hits to FILE (one per line: convention number)\n"
        "  -1        1-based hits, A064810 (position 1 = first decimal)\n"
        "  -0        0-based hits, A057680 (position 0 = first decimal)\n"
        "  -3        3-first hits, A057679 (the leading 3 is digit 1)\n"
        "  -d N      k at after-decimal position k-N; N may be negative and\n"
        "            may be a comma list.  -d 2,3,4,5 = A153220/21/23/24;\n"
        "            -d 0/-d -1/-d 1 are the same as -1/-0/-3\n"
        "  -x N      k at after-decimal position N*k; comma list allowed.\n"
        "            -x 2 = A153227, -x 3 = A153228\n"
        "  -r N      k REVERSED at after-decimal position k-N; N may be negative\n"
        "            and may be a comma list.  -r 0 = A366831, -r 1 = A366830\n"
        "  -m N      k mod N at after-decimal position k; comma list allowed.\n"
        "            -m 100 = A153225, -m 1000 = A153226.  Dense: expect about\n"
        "            2%% of positions to hit for -m 100 - use with -n\n"
        "  --pithy   pithy hits: the k-th non-overlapping m-tuple equals the\n"
        "            m-digit k; A109513 (after decimal) and A109514 (3 first)\n"
        "\n"
        "Selection flags may be combined; with none given, everything except\n"
        "the dense -m checks is run (offsets -1..5, reversed offsets 0 and 1,\n"
        "multipliers 2 and 3, both pithy variants).\n"
        "  -h        this help\n", argv0);
    exit(1);
}

/* Decimal counter held as digit values 0..9, most significant first. */
typedef struct {
    uint8_t d[SELF_L + 4];
    int len;
} deccnt;

/*
 * deccnt_inc - increment a decimal counter kept as a digit array (most
 * significant digit first), growing it when it carries out of the top.
 * Keeping candidates as digit strings lets a window compare be a plain
 * byte comparison instead of a conversion per digit.
 */
static inline void deccnt_inc(deccnt *c)
{
    int i = c->len - 1;
    while (i >= 0 && ++c->d[i] == 10) c->d[i--] = 0;
    if (i < 0) {
        memmove(c->d + 1, c->d, (size_t)c->len);
        c->d[0] = 1;
        c->len++;
        if (c->len > SELF_L) die("position exceeds %d digits", SELF_L);
    }
}

/*
 * deccnt_set - initialise a decimal counter from an integer.
 */
static void deccnt_set(deccnt *c, uint64_t v)
{
    char tmp[24];
    int n = snprintf(tmp, sizeof tmp, "%" PRIu64, v);
    c->len = n;
    for (int i = 0; i < n; i++) c->d[i] = (uint8_t)(tmp[i] - '0');
}

/*
 * deccnt_inc_rev - increment a counter stored least-significant digit
 * first, i.e. the array IS the reversed decimal string, so the reversed
 * conventions (-r) compare as cheaply as the forward ones.
 */
static inline void deccnt_inc_rev(deccnt *c)
{
    int i = 0;
    while (i < c->len && ++c->d[i] == 10) c->d[i++] = 0;
    if (i == c->len) {
        if (c->len >= SELF_L) die("position exceeds %d digits", SELF_L);
        c->d[c->len++] = 1;
    }
}

/*
 * deccnt_set_rev - initialise a least-significant-first counter.
 */
static void deccnt_set_rev(deccnt *c, uint64_t v)
{
    int n = 0;
    do { c->d[n++] = (uint8_t)(v % 10); v /= 10; } while (v);
    c->len = n;
}

/*
 * digits_match - true if the n digits at w equal the n digits at d.
 */
static inline int digits_match(const uint8_t *w, const uint8_t *d, int n)
{
    for (int i = 0; i < n; i++)
        if (w[i] != d[i]) return 0;
    return 1;
}

static uint64_t nhits;
static FILE *hits_fp;

/*
 * report - print and log a hit: convention label, OEIS id, the number and
 * the position it was found at; flushes the hits file per hit so an
 * interrupted run keeps everything found so far.
 */
static void report(const char *conv, const char *oeis, uint64_t pos, const deccnt *c)
{
    char b[40], s[SELF_L + 1];
    for (int i = 0; i < c->len; i++) s[i] = (char)('0' + c->d[i]);
    s[c->len] = 0;
    if (oeis && *oeis)
        printf("HIT (%s, %s): %s begins at position %s\n", conv, oeis, s, fmt_u64(pos, b, sizeof b));
    else
        printf("HIT (%s): %s begins at position %s\n", conv, s, fmt_u64(pos, b, sizeof b));
    fflush(stdout);
    if (hits_fp) {
        fprintf(hits_fp, "%s %s %s\n", conv, s, (oeis && *oeis) ? oeis : "-");
        fflush(hits_fp);           /* hits are rare; keep the file current */
    }
    nhits++;
}

/* ------------------------------------------------------------------ */
/* Sliding-window offset checks: window at position p equals p + d     */
/* ------------------------------------------------------------------ */
typedef struct {
    int64_t d;
    int rev;            /* compare the REVERSED digits of p + d        */
    char label[24];
    char oeis[12];
    deccnt c;           /* holds p + d once that is >= 0; digits are   */
                        /* LSB-first when rev is set                   */
    uint64_t skip;      /* positions to wait before counting (d < 0)   */
} offchk;

static offchk off[MAX_OFF];
static int noff;

/*
 * add_offset - register an offset convention (window at position p must
 * equal p+d), forward or reversed, with its label and OEIS id, and seed
 * its rolling counter so it reads p+d when the scan reaches p=1.
 */
static void add_offset(int64_t d, int rev)
{
    if (d < -1000000000000LL || d > 1000000000000LL) die("offset %lld out of range", (long long)d);
    for (int j = 0; j < noff; j++) if (off[j].d == d && off[j].rev == rev) return;
    if (noff >= MAX_OFF) die("too many -d/-r offsets (max %d)", MAX_OFF);
    offchk *oc = &off[noff++];
    memset(oc, 0, sizeof *oc);
    oc->d = d;
    oc->rev = rev;
    const char *pre = rev ? "rev-" : "";
    if (d == -1)     snprintf(oc->label, sizeof oc->label, "%s0-based", pre);
    else if (d == 0) snprintf(oc->label, sizeof oc->label, "%s1-based", pre);
    else if (d == 1) snprintf(oc->label, sizeof oc->label, "%s3-first", pre);
    else if (d > 0)  snprintf(oc->label, sizeof oc->label, "%sk-%lld", pre, (long long)d);
    else             snprintf(oc->label, sizeof oc->label, "%sk+%lld", pre, (long long)-d);
    if (rev)
        strcpy(oc->oeis, d == 0 ? "A366831" : d == 1 ? "A366830" : "");
    else
        strcpy(oc->oeis, d == -1 ? "A057680" : d == 0 ? "A064810" : d == 1 ? "A057679" :
                         d == 2  ? "A153220" : d == 3 ? "A153221" : d == 4 ? "A153223" :
                         d == 5  ? "A153224" : "");
    if (d >= 0) {                                  /* incremented to p+d at p=1 */
        if (rev) deccnt_set_rev(&oc->c, (uint64_t)d);
        else     deccnt_set(&oc->c, (uint64_t)d);
    } else {
        deccnt_set(&oc->c, 0);                     /* "0" is the same both ways */
        oc->skip = (uint64_t)-d;
    }
}

/*
 * parse_offsets - parse a comma list of -d/-r offsets into add_offset.
 */
static void parse_offsets(const char *arg, int rev)
{
    const char *s = arg;
    while (*s) {
        char *end;
        errno = 0;
        long long v = strtoll(s, &end, 10);
        if (end == s || errno != 0) die("bad offset list '%s'", arg);
        add_offset((int64_t)v, rev);
        s = end;
        if (*s == ',') s++;
        else if (*s) die("bad offset list '%s'", arg);
    }
}

/* ------------------------------------------------------------------ */
/* Multiplier checks: window at position m*k equals k (A153227/28)     */
/* ------------------------------------------------------------------ */
typedef struct {
    uint64_t mult;
    char label[24];
    char oeis[12];
    deccnt c;           /* holds k = p / mult when the check fires     */
    uint64_t countdown; /* positions until the next multiple of mult   */
} multchk;

static multchk mchk[MAX_OFF];
static int nmult;

/*
 * add_mult - register a multiplier convention (k found at position m*k):
 * a countdown fires the check only at multiples of m.
 */
static void add_mult(uint64_t m)
{
    if (m < 2 || m > 1000000000000ULL) die("multiplier %" PRIu64 " out of range (2..1e12)", m);
    for (int j = 0; j < nmult; j++) if (mchk[j].mult == m) return;
    if (nmult >= MAX_OFF) die("too many -x multipliers (max %d)", MAX_OFF);
    multchk *mc = &mchk[nmult++];
    memset(mc, 0, sizeof *mc);
    mc->mult = m;
    snprintf(mc->label, sizeof mc->label, "%" PRIu64 "k", m);
    strcpy(mc->oeis, m == 2 ? "A153227" : m == 3 ? "A153228" : "");
    deccnt_set(&mc->c, 0);          /* incremented to k=1 at p=m */
    mc->countdown = m;
}

/*
 * parse_mults - parse a comma list of -x multipliers into add_mult.
 */
static void parse_mults(const char *arg)
{
    const char *s = arg;
    while (*s) {
        char *end;
        errno = 0;
        unsigned long long v = strtoull(s, &end, 10);
        if (end == s || errno != 0) die("bad multiplier list '%s'", arg);
        add_mult((uint64_t)v);
        s = end;
        if (*s == ',') s++;
        else if (*s) die("bad multiplier list '%s'", arg);
    }
}

/* ------------------------------------------------------------------ */
/* Residue checks: window at position k equals k mod M (A153225/26)    */
/* ------------------------------------------------------------------ */
typedef struct {
    uint64_t mod;
    char label[24];
    char oeis[12];
    uint8_t *tab;       /* per residue: [len, digit, digit, ...]       */
    int esz;            /* entry size = 1 + ndigits(mod-1)             */
    uint64_t r;         /* p mod mod                                   */
} modchk;

static modchk rchk[MAX_OFF];
static int nmod;

/*
 * add_mod - register a residue convention (only k mod m must appear at
 * position k) with a precomputed table of residue digit strings.
 */
static void add_mod(uint64_t m)
{
    if (m < 2 || m > 10000000) die("modulus %" PRIu64 " out of range (2..1e7)", m);
    for (int j = 0; j < nmod; j++) if (rchk[j].mod == m) return;
    if (nmod >= MAX_OFF) die("too many -m moduli (max %d)", MAX_OFF);
    modchk *rc = &rchk[nmod++];
    memset(rc, 0, sizeof *rc);
    rc->mod = m;
    snprintf(rc->label, sizeof rc->label, "mod-%" PRIu64, m);
    strcpy(rc->oeis, m == 100 ? "A153225" : m == 1000 ? "A153226" : "");
    int dl = 1;
    for (uint64_t v = m - 1; v >= 10; v /= 10) dl++;
    rc->esz = 1 + dl;
    rc->tab = malloc((size_t)m * (size_t)rc->esz);
    if (!rc->tab) die("cannot allocate residue table for mod %" PRIu64, m);
    for (uint64_t v = 0; v < m; v++) {
        uint8_t *e = rc->tab + v * (uint64_t)rc->esz;
        char tmp[24];
        int n = snprintf(tmp, sizeof tmp, "%" PRIu64, v);
        e[0] = (uint8_t)n;
        for (int i = 0; i < n; i++) e[1 + i] = (uint8_t)(tmp[i] - '0');
    }
    rc->r = 0;
}

/*
 * parse_mods - parse a comma list of -m moduli into add_mod.
 */
static void parse_mods(const char *arg)
{
    const char *s = arg;
    while (*s) {
        char *end;
        errno = 0;
        unsigned long long v = strtoull(s, &end, 10);
        if (end == s || errno != 0) die("bad modulus list '%s'", arg);
        add_mod((uint64_t)v);
        s = end;
        if (*s == ',') s++;
        else if (*s) die("bad modulus list '%s'", arg);
    }
}

/* ------------------------------------------------------------------ */
/* Pithy numbers (A109513 / A109514): non-overlapping m-tuples         */
/* ------------------------------------------------------------------ */
/* For an m-digit k, the k-th m-tuple starts at exactly one position, so
 * instead of testing every window we keep, per convention and per tuple
 * length m, the next position where a check is due.  The main loop only
 * pays one comparison (p == pithy_next) per digit. */
#define PITHY_MAXM 18   /* (10^18-1)*18 still fits in a uint64 position */

typedef struct {
    deccnt k;           /* current tuple index = candidate number, m digits */
    uint64_t pos;       /* 1-based after-decimal position of that tuple    */
    int m;
    int active;
} pithy_state;

static pithy_state pithy[2][PITHY_MAXM + 1];  /* [0]=A109513, [1]=A109514 */
static uint64_t pithy_next = UINT64_MAX;
static const char *pithy_name[2] = { "pithy", "pithy-3" };
static const char *pithy_oeis[2] = { "A109513", "A109514" };

/*
 * pithy_recompute_next - recompute the earliest position at which any
 * pithy tuple check is due, so the main loop pays one compare per digit.
 */
static void pithy_recompute_next(void)
{
    uint64_t ne = UINT64_MAX;
    for (int v = 0; v < 2; v++)
        for (int m = 1; m <= PITHY_MAXM; m++)
            if (pithy[v][m].active && pithy[v][m].pos < ne)
                ne = pithy[v][m].pos;
    pithy_next = ne;
}

/*
 * pithy_init - set up the pithy-number states for tuple lengths 1..18 in
 * both conventions (tuples starting after the point, or at the 3): the
 * first m-digit candidate k = 10^(m-1) and the position of the k-th tuple.
 */
static void pithy_init(void)
{
    for (int v = 0; v < 2; v++)
        for (int m = 1; m <= PITHY_MAXM; m++) {
            pithy_state *t = &pithy[v][m];
            memset(t, 0, sizeof *t);
            t->m = m;
            t->k.len = m;
            t->k.d[0] = 1;                     /* k = 10^(m-1), first m-digit */
            uint64_t k0 = POW10[m - 1];
            /* A109513: tuple k covers after-decimal positions (k-1)m+1..km.
             * A109514: tuple k covers overall digits (k-1)m+1..km with the 3
             * as digit 1, i.e. after-decimal positions (k-1)m..km-1. */
            t->pos = (k0 - 1) * (uint64_t)m + (v == 0 ? 1 : 0);
            t->active = 1;
            if (t->pos == 0) {                 /* v=1, m=1: tuple 1 is the '3' */
                deccnt_inc(&t->k);
                t->pos = 1;
            }
        }
    pithy_recompute_next();
}

/*
 * pithy_events - handle every pithy check due at position p: compare the
 * window with the current tuple index, report hits, advance to the next
 * tuple (retiring a length once its index outgrows m digits) and
 * recompute the next due position.
 */
static void pithy_events(uint64_t p, const uint8_t *w, int room)
{
    for (int v = 0; v < 2; v++)
        for (int m = 1; m <= PITHY_MAXM; m++) {
            pithy_state *t = &pithy[v][m];
            if (!t->active || t->pos != p) continue;
            if (m <= room && digits_match(w, t->k.d, m)) {
                char b[40], s[SELF_L + 1];
                uint64_t kv = 0;
                for (int j = 0; j < m; j++) {
                    s[j] = (char)('0' + t->k.d[j]);
                    kv = kv * 10 + t->k.d[j];
                }
                s[m] = 0;
                unsigned last2 = (unsigned)(kv % 100), last1 = (unsigned)(kv % 10);
                const char *suf = (last2 >= 11 && last2 <= 13) ? "th" :
                                  last1 == 1 ? "st" : last1 == 2 ? "nd" :
                                  last1 == 3 ? "rd" : "th";
                printf("HIT (%s, %s): %s is the %s%s %d-tuple\n",
                       pithy_name[v], pithy_oeis[v], s, fmt_u64(kv, b, sizeof b), suf, m);
                fflush(stdout);
                if (hits_fp) { fprintf(hits_fp, "%s %s %s\n", pithy_name[v], s, pithy_oeis[v]); fflush(hits_fp); }
                nhits++;
            }
            deccnt_inc(&t->k);
            if (t->k.len > m) { t->active = 0; continue; }
            t->pos += (uint64_t)m;
        }
    pithy_recompute_next();
}

/*
 * main - parse options (defaulting to all sparse conventions), open the
 * digit stream (must start at position 1, since every check is relative
 * to the start of pi), and stream the digits once: per position, advance
 * each offset counter and compare it with the window, fire multiplier /
 * pithy / residue checks at their due positions, and report hits.
 */
int main(int argc, char **argv)
{
    const char *path = NULL, *out_path = NULL;
    uint64_t limit = UINT64_MAX;
    int wantp = 0, only = 0;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-p") && i + 1 < argc)      path = argv[++i];
        else if (!strcmp(argv[i], "-n") && i + 1 < argc) {
            if (parse_count(argv[++i], &limit) != 0) die("bad digit count '%s'", argv[i]);
        }
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) out_path = argv[++i];
        else if (!strcmp(argv[i], "-1")) { only = 1; add_offset(0, 0); }
        else if (!strcmp(argv[i], "-0")) { only = 1; add_offset(-1, 0); }
        else if (!strcmp(argv[i], "-3")) { only = 1; add_offset(1, 0); }
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) { only = 1; parse_offsets(argv[++i], 0); }
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) { only = 1; parse_offsets(argv[++i], 1); }
        else if (!strcmp(argv[i], "-x") && i + 1 < argc) { only = 1; parse_mults(argv[++i]); }
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) { only = 1; parse_mods(argv[++i]); }
        else if (!strcmp(argv[i], "--pithy")) { only = 1; wantp = 1; }
        else if (!strcmp(argv[i], "-h")) usage(argv[0]);
        else if (!strcmp(argv[i], "-V")) { print_version("piself"); return 0; }
        else if (argv[i][0] != '-' && !path) path = argv[i];
        else usage(argv[0]);
    }
    if (!path) usage(argv[0]);
    printf("piself v%s\n", PI_TOOLS_VERSION);
    if (!only) {                       /* default: check everything known */
        for (int64_t d = -1; d <= 5; d++) add_offset(d, 0);
        add_offset(0, 1);
        add_offset(1, 1);
        add_mult(2);
        add_mult(3);
        wantp = 1;
    }

    if (out_path) {
        hits_fp = fopen(out_path, "w");
        if (!hits_fp) die("cannot open output file '%s': %s", out_path, strerror(errno));
    }

    install_signals();

    pireader rd;
    pireader_open(&rd, path, SELF_L, 64 << 20, limit);
    if (rd.first_pos != 1)
        die("piself needs a source that starts at position 1 (text from \"3.\", or a ycd set including block 0); "
            "this source starts at position %" PRIu64, rd.first_pos);

    char b1[40], b2[40], d1[32], convs[512] = "";
    for (int j = 0; j < noff; j++) {
        if (*convs) strcat(convs, ", ");
        strcat(convs, off[j].label);
    }
    for (int j = 0; j < nmult; j++) {
        if (*convs) strcat(convs, ", ");
        strcat(convs, mchk[j].label);
    }
    for (int j = 0; j < nmod; j++) {
        if (*convs) strcat(convs, ", ");
        strcat(convs, rchk[j].label);
    }
    if (wantp) { if (*convs) strcat(convs, ", "); strcat(convs, "pithy"); pithy_init(); }
    printf("Scanning %s for self-locating numbers (%s; up to %s digits)\n",
           path, convs, limit == UINT64_MAX ? "all" : fmt_u64(limit, b1, sizeof b1));

    uint64_t scanned = 0, next_tick = 100000000;
    double t0 = now_sec();

    for (;;) {
        size_t nwin = pireader_fill(&rd);
        size_t avail = rd.ndigits;                 /* tail handling at EOF */
        size_t nchk = nwin ? nwin : avail;
        if (nchk == 0) break;
        const uint8_t *buf = rd.buf;
        uint64_t base = rd.base;

        for (size_t i = 0; i < nchk; i++) {
            uint64_t p = base + i;
            if (p <= scanned)   /* overlap region already checked; counters in step */
                continue;
            int room = nwin ? SELF_L : (int)(avail - i);

            /* Sliding-window checks: does the window here spell p + d
             * (or its digit reversal, for the -r family)? */
            for (int j = 0; j < noff; j++) {
                offchk *oc = &off[j];
                if (oc->skip) {            /* p + d still negative */
                    if (--oc->skip) continue;
                    /* p + d == 0 this position: c is already 0 */
                } else if (oc->rev) {
                    deccnt_inc_rev(&oc->c);
                } else {
                    deccnt_inc(&oc->c);    /* c = p + d */
                }
                if (oc->c.len <= room && digits_match(buf + i, oc->c.d, oc->c.len)) {
                    /* the classic conventions number positions so that
                     * position == the number itself; for the k-N family
                     * report where in pi (after decimal) it starts */
                    uint64_t pos = (oc->d >= -1 && oc->d <= 1)
                                   ? (uint64_t)((int64_t)p + oc->d) : p;
                    if (oc->rev) {         /* c is LSB-first: flip for display */
                        deccnt fwd;
                        fwd.len = oc->c.len;
                        for (int q = 0; q < fwd.len; q++)
                            fwd.d[q] = oc->c.d[fwd.len - 1 - q];
                        report(oc->label, oc->oeis, pos, &fwd);
                    } else {
                        report(oc->label, oc->oeis, pos, &oc->c);
                    }
                }
            }

            /* Multiplier checks: at p = m*k, does the window spell k? */
            for (int j = 0; j < nmult; j++) {
                multchk *mc = &mchk[j];
                if (--mc->countdown) continue;
                mc->countdown = mc->mult;
                deccnt_inc(&mc->c);        /* c = k = p / mult */
                if (mc->c.len <= room && digits_match(buf + i, mc->c.d, mc->c.len))
                    report(mc->label, mc->oeis, p, &mc->c);
            }

            /* Residue checks: does the window here spell p mod M? */
            for (int j = 0; j < nmod; j++) {
                modchk *rc = &rchk[j];
                if (++rc->r == rc->mod) rc->r = 0;
                const uint8_t *e = rc->tab + rc->r * (uint64_t)rc->esz;
                int rl = e[0];
                if (rl <= room && digits_match(buf + i, e + 1, rl)) {
                    char b[40], b2m[40], rs[24];
                    for (int q = 0; q < rl; q++) rs[q] = (char)('0' + e[1 + q]);
                    rs[rl] = 0;
                    if (rc->oeis[0])
                        printf("HIT (%s, %s): %s (%s begins at position %s)\n",
                               rc->label, rc->oeis, fmt_u64(p, b, sizeof b), rs,
                               fmt_u64(p, b2m, sizeof b2m));
                    else
                        printf("HIT (%s): %s (%s begins at position %s)\n",
                               rc->label, fmt_u64(p, b, sizeof b), rs,
                               fmt_u64(p, b2m, sizeof b2m));
                    fflush(stdout);
                    if (hits_fp) {
                        fprintf(hits_fp, "%s %" PRIu64 " %s\n", rc->label, p,
                                rc->oeis[0] ? rc->oeis : "-");
                        fflush(hits_fp);
                    }
                    nhits++;
                }
            }

            /* pithy: any tuple checks due at this position? */
            if (p == pithy_next)
                pithy_events(p, buf + i, room);

            scanned = p;
            if (p >= next_tick) {
                double dt = now_sec() - t0;
                fprintf(stderr, "  ... %s digits, %.1fs (%s digits/s)\n",
                        fmt_u64(p, b1, sizeof b1), dt,
                        fmt_u64((uint64_t)((double)p / dt), b2, sizeof b2));
                next_tick += 100000000;
            }
            if (g_stop) break;
        }
        if (g_stop || !nwin) break;
    }

    double dt = now_sec() - t0;
    char b3[40];
    printf("%s: scanned %s digits in %s (%s digits/s), %s hit%s\n",
           g_stop ? "Interrupted" : "Done",
           fmt_u64(scanned, b1, sizeof b1), fmt_duration(dt, d1, sizeof d1),
           fmt_u64(dt > 0 ? (uint64_t)((double)scanned / dt) : 0, b2, sizeof b2),
           fmt_u64(nhits, b3, sizeof b3), nhits == 1 ? "" : "s");

    if (hits_fp) fclose(hits_fp);
    pireader_close(&rd);
    return 0;
}
