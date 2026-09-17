/*
 * a053686.c
 *
 * Extend OEIS A053686:
 *
 *   "Record gaps between consecutive primes that repeat at least once before
 *    a new record occurs."
 *
 * Known terms (12):  2, 4, 6, 14, 34, 36, 52, 86, 132, 154, 250, 336.
 * They first occur after the primes 3, 7, 23, 113, 1327, 9551, 19609, 155921,
 * 1357201, 4652353, 387096133, 3842610773 (A133788).
 *
 * Definition
 * ----------
 * Let g_1 = 1 < g_2 = 2 < g_3 = 4 < ... be the record ("maximal") prime gaps
 * (A005250) and P_n the prime after which g_n first occurs (A002386): P_n + g_n
 * is prime and no gap between primes below P_(n+1) exceeds g_n.  The record g_n
 * is a term of A053686 iff it occurs a second time before the next record, i.e.
 * iff some prime p with P_n < p < P_(n+1) is followed by the gap g_n.
 * Equivalently the terms are the values that occur at least twice in A085237
 * (nondecreasing gaps).  Example: 336 first occurs after 3842610773 and again
 * after 4275912661 (4275912661 + 336 = 4275912997 is the next prime), before
 * the record 354 after 4302407359, so 336 = a(12).
 *
 * What is known
 * -------------
 * All maximal gaps below 2^64 are known (Nicely, Oliveira e Silva, Jacobsen
 * and others; REC[] below holds the 80 records up to 1550 after 1.836e19).
 * The OEIS entry decides the records through g_34 = 336.  A085237, computed
 * to its term 778 by Donovan Johnson (2008) and Charles Greathouse (2011),
 * lists 354, 382, ..., 766 once each, i.e. records 35..57 do not repeat, so
 * a(13) >= 778; that inference was never written into A053686, and nothing
 * beyond record 58 (778 after 4.28e13) has been examined.  Deciding records
 * 58..63 (778, 804, 806, 906, 916, 924) needs every gap up to
 * P_64 = 1693182318746371; record 64 (1132) would need 4.4e16.
 *
 * Method
 * ------
 * Deciding record n means looking at every prime gap in [P_n, P_(n+1)), so
 * the cost is one pass of a prime sieve over the interval.  Nothing cleverer
 * applies: the target is one exact gap size in the extreme tail, and a
 * partial sieve (small primes only) leaves ~5% survivors against ~3% primes,
 * so it would not reduce the primality work enough to pay for itself.
 * primesieve sieves at ~1.8e9 numbers/s per core near 10^15 on an M1 Pro,
 * and the per-prime bookkeeping here is one subtraction and one compare.
 *
 * The range [START, END) is cut into chunks handed to pthreads through an
 * atomic counter.  A worker sieves its chunk [lo, hi) with a primesieve
 * iterator, reads the primes straight out of the iterator's buffer and looks
 * at every gap (p, q) with lo <= p < hi (the gap that crosses hi belongs to
 * the chunk that owns p, so nothing is missed or counted twice).  Gaps of at
 * least g_n - 198, where n is the record interval containing p, are rare and
 * get a closer look: the gap g_(n+1) after P_(n+1) is checked against the
 * table (every record in the range must turn up exactly where A002386 says,
 * a strong check on the sieve), a gap equal to g_n is a repeat (a term), a
 * gap above g_n would be an unlisted maximal gap or a sieve error (reported
 * as ANOMALY), and all of them are counted in a per-record histogram of gaps
 * g_n - 2k (k < 100), which shows how close the near misses come.
 *
 * Completed chunks are folded into the totals in order (the frontier), which
 * is when their repeats are printed and when the record interval that ends
 * in the chunk is announced as DECIDED.  With -S FILE the position scanned
 * to, the prime count, the histogram and the repeat positions are
 * checkpointed; rerunning with the same START resumes from that position
 * (END and the chunk size may be changed between runs).
 *
 * Usage
 * -----
 *   a053686 scan [START] END [-t T] [-c CHUNK] [-s KIB] [-n NEAR] [-S FILE] [-i SECS] [-H] [-q]
 *       Examine every gap (p, q) with START <= p < END (START defaults to 0;
 *       START may be written rN for the prime P_N that starts record N).  Prints REPEAT / RECORD / DECIDED
 *       lines as they become final and a table at the end.  -n NEAR also
 *       prints gaps within NEAR of the record ("near" lines, NEAR <= 198).
 *       -H prints the histogram as machine-readable "hist N GAP COUNT" lines.
 *       -S FILE keeps a checkpoint (every SECS seconds, default 60, and on
 *       Ctrl-C); rerunning with the same START resumes from it, END may be
 *       larger than before.  -t T threads
 *       (default: all cores), -c CHUNK chunk size (default 10^10), -s KIB
 *       primesieve sieve size (default: library choice).
 *   a053686 verify P G
 *       Check with deterministic Miller-Rabin (and GMP when compiled in) that
 *       P and P+G are consecutive primes, and say which record interval the
 *       gap falls in and whether it is a repeat of that record.
 *   a053686 bench [N] [-d SPAN] [-t T] [-c CHUNK] [-s KIB]
 *       Time the scan of [N, N+SPAN) (defaults 10^15 and 10^11) and project
 *       when each of the next record decisions would be reached.
 *   a053686 selftest [LIMIT] [-t T]
 *       Check the table, the chunking against a brute-force pass, and the
 *       decisions for all records below LIMIT (default 10^11) against OEIS.
 *
 * Numbers may be written as decimal, 2^k, 10^k, 1e12, 1.7e15, or X+Y / X-Y
 * of those.
 *
 * Build:  cc -O2 -std=gnu11 -pthread -I/opt/homebrew/include a053686.c \
 *            -L/opt/homebrew/lib -lprimesieve -lm [-DHAVE_GMP -lgmp] -o a053686
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <limits.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>
#include <inttypes.h>
#include <primesieve.h>
#ifdef HAVE_GMP
#include <gmp.h>
#endif

typedef uint64_t u64;
typedef unsigned __int128 u128;

#define NREC      80        /* maximal gaps below 2^64 */
#define HK        100       /* histogram of gaps g_n - 2k, k = 0..HK-1, per record interval */
#define REPMAX    16        /* repeat positions kept per record (the count is in hist[n][0]) */
#define RINGC     4096      /* chunk completion ring: bounds how far workers run ahead of the frontier */
#define EVMAX     4096      /* events kept per chunk */
#define LOOKAHEAD 65536     /* sieve hint past a chunk end, for the gap that crosses it */
#define NKNOWN    12

/* A005250 / A002386: the record gap REC[i].g first occurs after the prime REC[i].p.
 * Index i is 0-based here; the OEIS index is n = i + 1.  Beyond 2^64 the known
 * records continue 1552 (after 18470057946260698231), 1572, 1676, 1724, 1854. */
static const struct { unsigned g; u64 p; } REC[NREC] = {
    {    1,                     2ULL},   /*  1 */
    {    2,                     3ULL},   /*  2 */
    {    4,                     7ULL},   /*  3 */
    {    6,                    23ULL},   /*  4 */
    {    8,                    89ULL},   /*  5 */
    {   14,                   113ULL},   /*  6 */
    {   18,                   523ULL},   /*  7 */
    {   20,                   887ULL},   /*  8 */
    {   22,                  1129ULL},   /*  9 */
    {   34,                  1327ULL},   /* 10 */
    {   36,                  9551ULL},   /* 11 */
    {   44,                 15683ULL},   /* 12 */
    {   52,                 19609ULL},   /* 13 */
    {   72,                 31397ULL},   /* 14 */
    {   86,                155921ULL},   /* 15 */
    {   96,                360653ULL},   /* 16 */
    {  112,                370261ULL},   /* 17 */
    {  114,                492113ULL},   /* 18 */
    {  118,               1349533ULL},   /* 19 */
    {  132,               1357201ULL},   /* 20 */
    {  148,               2010733ULL},   /* 21 */
    {  154,               4652353ULL},   /* 22 */
    {  180,              17051707ULL},   /* 23 */
    {  210,              20831323ULL},   /* 24 */
    {  220,              47326693ULL},   /* 25 */
    {  222,             122164747ULL},   /* 26 */
    {  234,             189695659ULL},   /* 27 */
    {  248,             191912783ULL},   /* 28 */
    {  250,             387096133ULL},   /* 29 */
    {  282,             436273009ULL},   /* 30 */
    {  288,            1294268491ULL},   /* 31 */
    {  292,            1453168141ULL},   /* 32 */
    {  320,            2300942549ULL},   /* 33 */
    {  336,            3842610773ULL},   /* 34 */
    {  354,            4302407359ULL},   /* 35 */
    {  382,           10726904659ULL},   /* 36 */
    {  384,           20678048297ULL},   /* 37 */
    {  394,           22367084959ULL},   /* 38 */
    {  456,           25056082087ULL},   /* 39 */
    {  464,           42652618343ULL},   /* 40 */
    {  468,          127976334671ULL},   /* 41 */
    {  474,          182226896239ULL},   /* 42 */
    {  486,          241160624143ULL},   /* 43 */
    {  490,          297501075799ULL},   /* 44 */
    {  500,          303371455241ULL},   /* 45 */
    {  514,          304599508537ULL},   /* 46 */
    {  516,          416608695821ULL},   /* 47 */
    {  532,          461690510011ULL},   /* 48 */
    {  534,          614487453523ULL},   /* 49 */
    {  540,          738832927927ULL},   /* 50 */
    {  582,         1346294310749ULL},   /* 51 */
    {  588,         1408695493609ULL},   /* 52 */
    {  602,         1968188556461ULL},   /* 53 */
    {  652,         2614941710599ULL},   /* 54 */
    {  674,         7177162611713ULL},   /* 55 */
    {  716,        13829048559701ULL},   /* 56 */
    {  766,        19581334192423ULL},   /* 57 */
    {  778,        42842283925351ULL},   /* 58 */
    {  804,        90874329411493ULL},   /* 59 */
    {  806,       171231342420521ULL},   /* 60 */
    {  906,       218209405436543ULL},   /* 61 */
    {  916,      1189459969825483ULL},   /* 62 */
    {  924,      1686994940955803ULL},   /* 63 */
    { 1132,      1693182318746371ULL},   /* 64 */
    { 1184,     43841547845541059ULL},   /* 65 */
    { 1198,     55350776431903243ULL},   /* 66 */
    { 1220,     80873624627234849ULL},   /* 67 */
    { 1224,    203986478517455989ULL},   /* 68 */
    { 1248,    218034721194214273ULL},   /* 69 */
    { 1272,    305405826521087869ULL},   /* 70 */
    { 1328,    352521223451364323ULL},   /* 71 */
    { 1356,    401429925999153707ULL},   /* 72 */
    { 1370,    418032645936712127ULL},   /* 73 */
    { 1442,    804212830686677669ULL},   /* 74 */
    { 1476,   1425172824437699411ULL},   /* 75 */
    { 1488,   5733241593241196731ULL},   /* 76 */
    { 1510,   6787988999657777797ULL},   /* 77 */
    { 1526,  15570628755536096243ULL},   /* 78 */
    { 1530,  17678654157568189057ULL},   /* 79 */
    { 1550,  18361375334787046697ULL},   /* 80 */
};

/* the records (1-based OEIS index n) that repeat, i.e. A053686 = { REC[n-1].g } */
static const int KNOWN_N[NKNOWN] = { 2, 3, 4, 6, 10, 11, 13, 15, 20, 22, 29, 34 };
#define A085237_LAST 57     /* A085237 (to its term 778) shows records 35..57 do not repeat */

static bool stderr_tty;
static volatile sig_atomic_t g_stop = 0;

/* ------------------------------------------------------------------ */
/* Generic helpers                                                     */
/* ------------------------------------------------------------------ */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void die(const char *fmt, ...)
{
    va_list ap;
    if (stderr_tty) fputs("\r\033[K", stderr);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

static const char *fmt_dur(double s, char *buf, size_t n)
{
    if (s < 0 || s != s || s > 1e9) { snprintf(buf, n, "?"); return buf; }
    long t = (long)s;
    if (t >= 86400) snprintf(buf, n, "%ldd%02ldh%02ldm", t / 86400, (t / 3600) % 24, (t / 60) % 60);
    else if (t >= 3600) snprintf(buf, n, "%ldh%02ldm%02lds", t / 3600, (t / 60) % 60, t % 60);
    else if (t >= 60) snprintf(buf, n, "%ldm%02lds", t / 60, t % 60);
    else snprintf(buf, n, "%.1fs", s);
    return buf;
}

static bool mul_ok(u64 a, u64 b, u64 *r)
{
    u128 x = (u128)a * b;
    if (x >> 64) return false;
    *r = (u64)x;
    return true;
}

static bool all_digits(const char *s, size_t n)
{
    if (n == 0) return false;
    for (size_t i = 0; i < n; i++)
        if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

/* decimal, B^E, MeE (M may have a decimal point) */
static bool parse_term(const char *s, size_t n, u64 *out)
{
    char buf[64], *c;
    if (n == 0 || n >= sizeof buf) return false;
    memcpy(buf, s, n);
    buf[n] = 0;
    if ((c = strchr(buf, '^'))) {
        *c = 0;
        if (!all_digits(buf, strlen(buf)) || !all_digits(c + 1, strlen(c + 1))) return false;
        u64 b = strtoull(buf, 0, 10), e = strtoull(c + 1, 0, 10), r = 1;
        while (e--) if (!mul_ok(r, b, &r)) return false;
        *out = r;
        return true;
    }
    if ((c = strpbrk(buf, "eE"))) {
        *c = 0;
        if (!all_digits(c + 1, strlen(c + 1))) return false;
        long e = strtol(c + 1, 0, 10);
        char mant[64];
        int md = 0, frac = -1;
        for (const char *q = buf; *q; q++) {
            if (*q == '.') { if (frac >= 0) return false; frac = 0; continue; }
            if (!isdigit((unsigned char)*q)) return false;
            mant[md++] = *q;
            if (frac >= 0) frac++;
        }
        if (md == 0 || md > 19) return false;
        mant[md] = 0;
        if (frac > 0) e -= frac;
        if (e < 0) return false;
        u64 m = strtoull(mant, 0, 10);
        while (e-- > 0) if (!mul_ok(m, 10, &m)) return false;
        *out = m;
        return true;
    }
    if (!all_digits(buf, n) || n > 20) return false;
    errno = 0;
    u64 v = strtoull(buf, 0, 10);
    if (errno) return false;
    *out = v;
    return true;
}

static u64 parse_num(const char *s)
{
    size_t n = strlen(s);
    if ((s[0] == 'r' || s[0] == 'R') && all_digits(s + 1, n - 1)) {   /* rN = start prime of record N */
        u64 k = strtoull(s + 1, 0, 10);
        if (k < 1 || k > NREC) die("record number in '%s' must be 1..%d", s, NREC);
        return REC[k - 1].p;
    }
    for (size_t i = n; i-- > 1; ) {           /* split at the last top-level + or - */
        if (s[i] != '+' && s[i] != '-') continue;
        u64 a, b;
        if (!parse_term(s, i, &a) || !parse_term(s + i + 1, n - i - 1, &b)) die("bad number '%s'", s);
        if (s[i] == '+') {
            if (a + b < a) die("overflow in '%s'", s);
            return a + b;
        }
        if (b > a) die("negative result in '%s'", s);
        return a - b;
    }
    u64 v;
    if (!parse_term(s, n, &v)) die("bad number '%s'", s);
    return v;
}

static int parse_int(const char *s, int lo, int hi, const char *what)
{
    u64 v = parse_num(s);
    if (v < (u64)lo || v > (u64)hi) die("%s must be between %d and %d", what, lo, hi);
    return (int)v;
}

/* ------------------------------------------------------------------ */
/* Independent primality: deterministic Miller-Rabin below 2^64         */
/* ------------------------------------------------------------------ */

static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }

static u64 powmod(u64 a, u64 e, u64 m)
{
    u64 r = 1;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod(r, a, m);
        a = mulmod(a, a, m);
        e >>= 1;
    }
    return r;
}

static bool is_prime64(u64 n)
{
    static const u64 small[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
    static const u64 bases[] = { 2, 325, 9375, 28178, 450775, 9780504, 1795265022 };
    if (n < 2) return false;
    for (size_t i = 0; i < sizeof small / sizeof *small; i++)
        if (n % small[i] == 0) return n == small[i];
    u64 d = n - 1;
    int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    for (size_t i = 0; i < sizeof bases / sizeof *bases; i++) {
        u64 a = bases[i] % n;
        if (a == 0) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        bool composite = true;
        for (int j = 1; j < r && composite; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) composite = false;
        }
        if (composite) return false;
    }
    return true;
}

static u64 next_prime64(u64 n)
{
    if (n < 2) return 2;
    u64 m = n + 1;
    if (!(m & 1)) m++;
    while (!is_prime64(m)) m += 2;
    return m;
}

/* ------------------------------------------------------------------ */
/* Record table helpers                                                */
/* ------------------------------------------------------------------ */

/* index i with REC[i].p <= x < REC[i+1].p (0 for x < 3) */
static int rec_index(u64 x)
{
    int lo = 0, hi = NREC - 1;
    while (lo < hi) {
        int mid = (lo + hi + 1) / 2;
        if (REC[mid].p <= x) lo = mid; else hi = mid - 1;
    }
    return lo;
}

static u64 rec_end(int i) { return i + 1 < NREC ? REC[i + 1].p : UINT64_MAX; }

static int known_term_index(int n)     /* a(k) = g_n?  returns k or 0 */
{
    for (int k = 0; k < NKNOWN; k++)
        if (KNOWN_N[k] == n) return k + 1;
    return 0;
}

/* how a decision for record n (1-based) compares with what OEIS has */
static const char *oeis_tag(int n, bool repeats)
{
    static char buf[96];
    if (n <= KNOWN_N[NKNOWN - 1]) {
        int k = known_term_index(n);
        if (repeats && k) snprintf(buf, sizeof buf, "= OEIS a(%d)", k);
        else if (!repeats && !k) snprintf(buf, sizeof buf, "as in OEIS");
        else snprintf(buf, sizeof buf, "DIFFERS FROM OEIS (%s)", k ? "listed as a term" : "not listed");
    } else if (n <= A085237_LAST) {
        snprintf(buf, sizeof buf, repeats ? "NEW, contradicts A085237" : "agrees with A085237");
    } else {
        snprintf(buf, sizeof buf, repeats ? "NEW" : "new information, undecided in OEIS");
    }
    return buf;
}

/* ------------------------------------------------------------------ */
/* Scan engine                                                         */
/* ------------------------------------------------------------------ */

enum { EV_REPEAT, EV_RECORD, EV_ANOMALY, EV_NEAR };

typedef struct { unsigned char type, i; u64 p, g; } event_t;
typedef struct { unsigned char i, k; u64 c; } hent_t;

typedef struct {                /* per-chunk result, held until the frontier folds it */
    bool done;
    u64 primes;
    unsigned nh, nev, dropped;
    hent_t *h;
    event_t *ev;
} slot_t;

static struct {
    u64 start, end, chunk, nchunks;
    u64 grid0;                      /* chunk i of this run is [grid0 + i*chunk, ...): start, or the resumed position */
    _Atomic u64 next;               /* next chunk index to hand out */
    /* protected by g_mu: */
    u64 frontier;                   /* chunks [0, frontier) of this run are complete and folded */
    slot_t slot[RINGC];
    u64 primes;                     /* primes in [start, frontier position) */
    u64 hist[NREC][HK];             /* hist[i][k] = gaps of size REC[i].g - 2k after primes in interval i (not at REC[i].p) */
    u64 rep[NREC][REPMAX];          /* the first repeat positions, ascending */
    unsigned nrep[NREC];
    bool recseen[NREC];             /* the gap after REC[i].p was checked */
    u64 anomalies;
    int last_rec;                   /* highest 1-based record confirmed so far, 0 = none */
} S;
static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_sieve_kib = 0;         /* primesieve sieve size in KiB, 0 = library default */
static unsigned g_near = 0;         /* print gaps within this of the record */
static bool g_hist_lines = false;
static bool g_quiet_events = false;  /* bookkeeping only, no REPEAT/RECORD/DECIDED lines */

typedef struct {
    int ci;                         /* record interval containing the current p */
    u64 (*lhist)[HK];               /* NREC x HK */
    event_t *ev;
    unsigned nev, dropped;
} cstate_t;

typedef struct {
    int id;
    primesieve_iterator it;
    cstate_t c;
} worker_t;

static void add_event(cstate_t *c, int type, int i, u64 p, u64 g)
{
    if (c->nev >= EVMAX) {
        if (type == EV_NEAR) { c->dropped++; return; }
        /* never lose a repeat, record or anomaly: overwrite the last near miss */
        for (unsigned j = c->nev; j-- > 0; )
            if (c->ev[j].type == EV_NEAR) { c->dropped++; c->nev = j; break; }
        if (c->nev >= EVMAX) { c->dropped++; return; }
    }
    event_t *e = &c->ev[c->nev++];
    e->type = (unsigned char)type;
    e->i = (unsigned char)i;
    e->p = p;
    e->g = g;
}

/* A gap g after the prime p that is large for its record interval. */
static inline void big_gap(cstate_t *c, u64 p, u64 g)
{
    while (c->ci + 1 < NREC && p >= REC[c->ci + 1].p) c->ci++;
    const int ci = c->ci;
    const unsigned G = REC[ci].g;
    if (p == REC[ci].p) {                       /* the record itself */
        add_event(c, g == G ? EV_RECORD : EV_ANOMALY, ci, p, g);
        return;
    }
    if (g > G) {                                /* larger than the record of its interval */
        add_event(c, EV_ANOMALY, ci, p, g);
        return;
    }
    const u64 k = (G - g) >> 1;                 /* both even: ci >= 1 here (interval 0 is just p = 2) */
    if (k < HK) c->lhist[ci][k]++;
    if (g == G) add_event(c, EV_REPEAT, ci, p, g);
    else if (G - g <= g_near) add_event(c, EV_NEAR, ci, p, g);
}

static void fold_locked(void);

static void scan_chunk(worker_t *w, u64 idx)
{
    const u64 lo = S.grid0 + idx * S.chunk;
    u64 hi = lo + S.chunk;
    if (hi > S.end || hi < lo) hi = S.end;
    cstate_t *c = &w->c;

    memset(c->lhist, 0, sizeof(u64) * NREC * HK);
    c->nev = c->dropped = 0;
    c->ci = rec_index(lo);
    const u64 thr = REC[c->ci].g > 2 * (HK - 1) ? REC[c->ci].g - 2 * (HK - 1) : 1;
    u64 count = 0;

    primesieve_jump_to(&w->it, lo, hi + LOOKAHEAD);
    /* The primes are read straight out of the iterator's buffer, which
     * primesieve_generate_next_primes() fills with primes[0..size). */
    primesieve_generate_next_primes(&w->it);
    if (w->it.is_error || w->it.size == 0) die("primesieve error in chunk %" PRIu64, idx);
    const u64 *P = w->it.primes;
    size_t np = w->it.size, i = 1;
    u64 prev = P[0];
    if (prev >= hi) goto finish;                    /* no prime in this chunk: the gap over it belongs to the previous one */
    count = 1;
    for (;;) {
        for (; i < np; i++) {
            const u64 q = P[i];
            const u64 g = q - prev;
            if (g >= thr) big_gap(c, prev, g);
            if (q >= hi) goto finish;               /* the gap crossing hi is done; q belongs to the next chunk */
            count++;
            prev = q;
        }
        primesieve_generate_next_primes(&w->it);
        if (w->it.is_error || w->it.size == 0) die("primesieve error in chunk %" PRIu64, idx);
        P = w->it.primes;
        np = w->it.size;
        i = 0;
    }

finish:;
    /* compact the histogram */
    unsigned nh = 0;
    for (int a = 0; a < NREC; a++)
        for (int k = 0; k < HK; k++)
            if (c->lhist[a][k]) nh++;
    hent_t *h = nh ? malloc(nh * sizeof *h) : NULL;
    event_t *ev = c->nev ? malloc(c->nev * sizeof *ev) : NULL;
    if ((nh && !h) || (c->nev && !ev)) die("out of memory");
    for (int a = 0, j = 0; a < NREC; a++)
        for (int k = 0; k < HK; k++)
            if (c->lhist[a][k]) { h[j].i = (unsigned char)a; h[j].k = (unsigned char)k; h[j].c = c->lhist[a][k]; j++; }
    if (c->nev) memcpy(ev, c->ev, c->nev * sizeof *ev);

    pthread_mutex_lock(&g_mu);
    slot_t *s = &S.slot[idx % RINGC];
    s->done = true;
    s->primes = count;
    s->nh = nh;
    s->h = h;
    s->nev = c->nev;
    s->ev = ev;
    s->dropped = c->dropped;
    fold_locked();
    pthread_mutex_unlock(&g_mu);
}

static void clear_line(void) { if (stderr_tty) fputs("\r\033[K", stderr); }

static void announce_decision_locked(int i, u64 upto)   /* interval i (0-based) has been scanned up to `upto` */
{
    const int n = i + 1;
    const u64 reps = S.hist[i][0];
    const bool full = S.start <= REC[i].p && upto >= rec_end(i);
    if (g_quiet_events) return;
    clear_line();
    if (full) {
        if (reps)
            printf("DECIDED  #%-2d gap %-5u repeats %" PRIu64 " time%s before the next record, first after %" PRIu64
                   "  -> TERM  [%s]\n", n, REC[i].g, reps, reps == 1 ? "" : "s", S.rep[i][0], oeis_tag(n, true));
        else
            printf("DECIDED  #%-2d gap %-5u no second occurrence in (%" PRIu64 ", %" PRIu64 ")  -> not a term  [%s]\n",
                   n, REC[i].g, REC[i].p, rec_end(i), oeis_tag(n, false));
    } else {
        u64 from = S.start > REC[i].p ? S.start : REC[i].p, to = upto < rec_end(i) ? upto : rec_end(i);
        if (reps)
            printf("PARTIAL  #%-2d gap %-5u repeats %" PRIu64 " time%s in the scanned part [%" PRIu64 ", %" PRIu64
                   "), first after %" PRIu64 "  -> TERM  [%s]\n", n, REC[i].g, reps, reps == 1 ? "" : "s", from, to,
                   S.rep[i][0], oeis_tag(n, true));
        else
            printf("PARTIAL  #%-2d gap %-5u no second occurrence in the scanned part [%" PRIu64 ", %" PRIu64 ") of (%" PRIu64
                   ", %" PRIu64 ")\n", n, REC[i].g, from, to, REC[i].p, rec_end(i));
    }
    fflush(stdout);
}

static void handle_event_locked(const event_t *e)
{
    const int i = e->i, n = i + 1;
    switch (e->type) {
    case EV_REPEAT:
        if (S.nrep[i] < REPMAX) S.rep[i][S.nrep[i]++] = e->p;
        if (g_quiet_events) return;
        clear_line();
        printf("REPEAT   #%-2d gap %-5u after %" PRIu64 " (next prime %" PRIu64 "); first occurrence after %" PRIu64
               "  [%s]\n", n, REC[i].g, e->p, e->p + e->g, REC[i].p, oeis_tag(n, true));
        break;
    case EV_RECORD:
        S.recseen[i] = true;
        if (n > S.last_rec) S.last_rec = n;
        if (g_quiet_events) return;
        clear_line();
        printf("RECORD   #%-2d gap %-5u after %" PRIu64 " confirmed where A002386 lists it\n", n, REC[i].g, e->p);
        if (i >= 1 && S.start < REC[i].p) announce_decision_locked(i - 1, REC[i].p);
        break;
    case EV_ANOMALY:
        S.anomalies++;
        clear_line();
        if (e->p == REC[i].p)
            printf("ANOMALY  gap after %" PRIu64 " is %" PRIu64 " but A005250 says record #%d is %u\n", e->p, e->g, n, REC[i].g);
        else
            printf("ANOMALY  gap %" PRIu64 " after %" PRIu64 " exceeds the record %u (#%d) of its interval: unlisted maximal gap or sieve error\n",
                   e->g, e->p, REC[i].g, n);
        break;
    case EV_NEAR:
        if (g_quiet_events) return;
        clear_line();
        printf("near     #%-2d gap %-5" PRIu64 " after %" PRIu64 " (record %u, short by %" PRIu64 ")\n", n, e->g, e->p, REC[i].g, REC[i].g - e->g);
        break;
    }
    fflush(stdout);
}

/* Fold complete chunks in frontier order (caller holds g_mu). */
static void fold_locked(void)
{
    while (S.frontier < S.nchunks && S.slot[S.frontier % RINGC].done) {
        slot_t *s = &S.slot[S.frontier % RINGC];
        s->done = false;
        S.primes += s->primes;
        for (unsigned j = 0; j < s->nh; j++) S.hist[s->h[j].i][s->h[j].k] += s->h[j].c;
        for (unsigned j = 0; j < s->nev; j++) handle_event_locked(&s->ev[j]);
        if (s->dropped) {
            clear_line();
            printf("note     %u near-miss lines dropped in chunk %" PRIu64 " (lower -n)\n", s->dropped, S.frontier);
        }
        free(s->h);
        free(s->ev);
        s->h = NULL;
        s->ev = NULL;
        s->nh = s->nev = s->dropped = 0;
        S.frontier++;
    }
}

static void *worker(void *arg)
{
    worker_t *w = arg;
    primesieve_init(&w->it);
    w->c.lhist = malloc(sizeof(u64) * NREC * HK);
    w->c.ev = malloc(sizeof(event_t) * EVMAX);
    if (!w->c.lhist || !w->c.ev) die("out of memory");
    while (!g_stop) {
        u64 idx = atomic_fetch_add(&S.next, 1);
        if (idx >= S.nchunks) break;
        for (;;) {                                   /* stay within RINGC of the frontier */
            pthread_mutex_lock(&g_mu);
            u64 f = S.frontier;
            pthread_mutex_unlock(&g_mu);
            if (idx < f + RINGC - 1 || g_stop) break;
            usleep(2000);
        }
        if (g_stop) break;
        scan_chunk(w, idx);
    }
    free(w->c.lhist);
    free(w->c.ev);
    primesieve_free_iterator(&w->it);
    return NULL;
}

/* ---- checkpoint file ---- */

static void save_state(const char *path)
{
    char tmp[PATH_MAX];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "\nwarning: cannot write %s: %s\n", tmp, strerror(errno)); return; }
    pthread_mutex_lock(&g_mu);
    u64 pos = S.grid0 + S.frontier * S.chunk;
    if (pos > S.end) pos = S.end;
    fprintf(f, "A053686 checkpoint 2\nstart %" PRIu64 "\npos %" PRIu64 "\nchunk %" PRIu64 "\nend %" PRIu64
               "\nprimes %" PRIu64 "\nanomalies %" PRIu64 "\nlastrec %d\n",
            S.start, pos, S.chunk, S.end, S.primes, S.anomalies, S.last_rec);
    for (int i = 0; i < NREC; i++)
        if (S.recseen[i]) fprintf(f, "rec %d\n", i);
    for (int i = 0; i < NREC; i++)
        for (unsigned j = 0; j < S.nrep[i]; j++) fprintf(f, "rep %d %" PRIu64 "\n", i, S.rep[i][j]);
    for (int i = 0; i < NREC; i++)
        for (int k = 0; k < HK; k++)
            if (S.hist[i][k]) fprintf(f, "hist %d %d %" PRIu64 "\n", i, k, S.hist[i][k]);
    pthread_mutex_unlock(&g_mu);
    fputs("end\n", f);
    if (fclose(f) != 0 || rename(tmp, path) != 0)
        fprintf(stderr, "\nwarning: cannot update %s: %s\n", path, strerror(errno));
}

/* Returns false if there is no checkpoint; dies on a mismatch. */
static bool load_state(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        if (errno == ENOENT) return false;
        die("cannot read %s: %s", path, strerror(errno));
    }
    char line[256];
    if (!fgets(line, sizeof line, f) || strncmp(line, "A053686 checkpoint 2", 20) != 0)
        die("%s is not an A053686 checkpoint (format 2)", path);
    u64 start = 0, pos = 0, primes = 0, anomalies = 0;
    int lastrec = 0;
    bool complete = false, have_pos = false;
    while (fgets(line, sizeof line, f)) {
        int i, k;
        u64 v;
        if (sscanf(line, "start %" SCNu64, &start) == 1) continue;
        if (sscanf(line, "pos %" SCNu64, &pos) == 1) { have_pos = true; continue; }
        if (sscanf(line, "chunk %" SCNu64, &v) == 1) continue;      /* informational */
        if (sscanf(line, "end %" SCNu64, &v) == 1) continue;        /* informational */
        if (sscanf(line, "primes %" SCNu64, &primes) == 1) continue;
        if (sscanf(line, "anomalies %" SCNu64, &anomalies) == 1) continue;
        if (sscanf(line, "lastrec %d", &lastrec) == 1) continue;
        if (sscanf(line, "rec %d", &i) == 1) {
            if (i < 0 || i >= NREC) die("%s: bad rec line", path);
            S.recseen[i] = true;
            continue;
        }
        if (sscanf(line, "rep %d %" SCNu64, &i, &v) == 2) {
            if (i < 0 || i >= NREC) die("%s: bad rep line", path);
            if (S.nrep[i] < REPMAX) S.rep[i][S.nrep[i]++] = v;
            continue;
        }
        if (sscanf(line, "hist %d %d %" SCNu64, &i, &k, &v) == 3) {
            if (i < 0 || i >= NREC || k < 0 || k >= HK) die("%s: bad hist line", path);
            S.hist[i][k] = v;
            continue;
        }
        if (strncmp(line, "end", 3) == 0) { complete = true; break; }
    }
    fclose(f);
    if (!complete || !have_pos) die("%s is truncated", path);
    if (start != S.start)
        die("%s was written for START=%" PRIu64 ", not %" PRIu64, path, start, S.start);
    if (pos < S.start) die("%s: bad position", path);
    if (pos > S.end)
        die("%s already covers [%" PRIu64 ", %" PRIu64 "), beyond END=%" PRIu64, path, start, pos, S.end);
    /* lay out this run's chunks from the resumed position (END and CHUNK may differ from the earlier run) */
    S.grid0 = pos;
    S.nchunks = (S.end - pos + S.chunk - 1) / S.chunk;
    S.frontier = 0;
    S.primes = primes;
    S.anomalies = anomalies;
    S.last_rec = lastrec;
    return true;
}

static void on_sigint(int sig)
{
    (void)sig;
    g_stop = 1;
}

typedef struct {
    double seconds;
    u64 primes;
    u64 numbers;
} scan_stats_t;

static u64 frontier_pos(void)
{
    u64 pos = S.grid0 + S.frontier * S.chunk;
    return pos > S.end ? S.end : pos;
}

static void reset_state(u64 start, u64 end, u64 chunk)
{
    if (end <= start) die("END must be larger than START");
    if (chunk == 0) die("CHUNK must be positive");
    S.start = start;
    S.grid0 = start;
    S.chunk = chunk;
    S.nchunks = (end - start + chunk - 1) / chunk;      /* the last chunk may be short */
    if (S.nchunks == 0 || (S.nchunks - 1) > (UINT64_MAX - start) / chunk) die("range too large");
    S.end = end;
    if (S.end > primesieve_get_max_stop() - 2 * LOOKAHEAD) die("END beyond primesieve's limit");
    S.frontier = 0;
    S.primes = 0;
    S.anomalies = 0;
    S.last_rec = 0;
    for (u64 i = 0; i < RINGC; i++) { free(S.slot[i].h); free(S.slot[i].ev); }
    memset(S.slot, 0, sizeof S.slot);
    memset(S.hist, 0, sizeof S.hist);
    memset(S.rep, 0, sizeof S.rep);
    memset(S.nrep, 0, sizeof S.nrep);
    memset(S.recseen, 0, sizeof S.recseen);
}

/* Scan [start, end) with the given number of threads. */
static scan_stats_t run_scan(u64 start, u64 end, u64 chunk, int threads, const char *state, int interval, bool quiet)
{
    reset_state(start, end, chunk);
    if (start > 0 && start != REC[rec_index(start)].p && !quiet)
        fprintf(stderr, "note: START is inside record interval #%d (gap %u after %" PRIu64 "), which will only be partially scanned\n",
                rec_index(start) + 1, REC[rec_index(start)].g, REC[rec_index(start)].p);

    if (state && load_state(state))
        fprintf(stderr, "resuming from %s: [%" PRIu64 ", %" PRIu64 ") done, %" PRIu64 " primes examined; %" PRIu64
                        " chunks to go\n", state, S.start, S.grid0, S.primes, S.nchunks);
    atomic_store(&S.next, S.frontier);

    if (!quiet)
        fprintf(stderr, "scanning gaps after the primes in [%" PRIu64 ", %" PRIu64 ") in %" PRIu64 " chunks of %" PRIu64
                        " with %d threads\n", S.start, S.end, S.nchunks, S.chunk, threads);

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_sigint;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    worker_t *ws = calloc((size_t)threads, sizeof *ws);
    pthread_t *th = calloc((size_t)threads, sizeof *th);
    if (!ws || !th) die("out of memory");
    for (int i = 0; i < threads; i++) {
        ws[i].id = i;
        if (pthread_create(&th[i], NULL, worker, &ws[i]) != 0) die("pthread_create failed");
    }

    /* rate window: one sample per second, ~2 minutes deep */
    enum { NS = 128 };
    double ts[NS];
    u64 ps[NS];
    int ns = 0;
    const double t0 = now();
    const u64 pos0 = frontier_pos();
    double last_status = 0, last_ckpt = t0, last_sample = 0;
    ts[0] = t0; ps[0] = pos0; ns = 1;
    char b1[32], b2[32];

    for (;;) {
        usleep(200000);
        double t = now();
        pthread_mutex_lock(&g_mu);
        u64 f = S.frontier, primes = S.primes, pos = frontier_pos(), reps = 0;
        int lastrec = S.last_rec;
        for (int i = 0; i < NREC; i++) reps += S.hist[i][0];
        pthread_mutex_unlock(&g_mu);
        if (f >= S.nchunks || g_stop) break;

        if (t - last_sample >= 1.0) {
            last_sample = t;
            if (ns < NS) { ts[ns] = t; ps[ns] = pos; ns++; }
            else { memmove(ts, ts + 1, (NS - 1) * sizeof *ts); memmove(ps, ps + 1, (NS - 1) * sizeof *ps); ts[NS - 1] = t; ps[NS - 1] = pos; }
        }
        double rate = (t - ts[0] > 0.5) ? (double)(pos - ps[0]) / (t - ts[0]) : 0;
        double eta = rate > 0 ? (double)(S.end - pos) / rate : -1;
        bool show = stderr_tty ? (t - last_status >= 1.0) : (t - last_status >= interval);
        if (show && !quiet) {
            last_status = t;
            fprintf(stderr, "%s%s  pos %.6g (%.2f%%)  %.3g/s  %.3g primes  ETA %s  records to #%d  repeats %" PRIu64,
                    stderr_tty ? "\r\033[K" : "progress: ",
                    fmt_dur(t - t0, b1, sizeof b1), (double)pos,
                    100.0 * (double)(pos - S.start) / (double)(S.end - S.start),
                    rate, (double)primes, fmt_dur(eta, b2, sizeof b2), lastrec, reps);
            if (!stderr_tty) fputc('\n', stderr);
            fflush(stderr);
        }
        if (state && t - last_ckpt >= interval) {
            last_ckpt = t;
            save_state(state);
        }
    }
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    double t1 = now();
    if (stderr_tty && !quiet) fputs("\r\033[K", stderr);

    if (state) save_state(state);
    if (g_stop && S.frontier < S.nchunks && !quiet)
        fprintf(stderr, "stopped at position %" PRIu64 " (%" PRIu64 " of %" PRIu64 " chunks complete)%s\n",
                frontier_pos(), S.frontier, S.nchunks, state ? ", checkpoint saved" : "");

    scan_stats_t st;
    st.seconds = t1 - t0;
    st.primes = S.primes;
    st.numbers = frontier_pos() - pos0;
    free(ws);
    free(th);
    return st;
}

/* pi(10^k) for k = 6..18 (A006880) */
static const u64 PI10[19] = {
    0, 0, 0, 0, 0, 0, 78498ULL, 664579ULL, 5761455ULL, 50847534ULL, 455052511ULL,
    4118054813ULL, 37607912018ULL, 346065536839ULL, 3204941750802ULL,
    29844570422669ULL, 279238341033925ULL, 2623557157654233ULL, 24739954287740860ULL
};

/* Summary of the scan [S.start, frontier position): one line per record interval touched. */
static int print_summary(void)
{
    const u64 pos = frontier_pos();
    int fails = 0;
    printf("\nRecord intervals touched by the scan [%" PRIu64 ", %" PRIu64 "):\n", S.start, pos);
    printf("  n   gap   first after            next record after      repeats   runner-up gap (count)  gaps within 20 / 100 of record   status\n");
    for (int i = 0; i < NREC; i++) {
        const u64 lo = REC[i].p, hi = rec_end(i);
        if (hi <= S.start || lo >= pos) continue;
        const bool full = S.start <= lo && pos >= hi;
        const u64 reps = S.hist[i][0];
        u64 w20 = 0, w100 = 0;
        int krun = 0;
        for (int k = 1; k < HK; k++) {
            if (S.hist[i][k] && !krun) krun = k;
            if (k <= 10) w20 += S.hist[i][k];
            if (k <= 50) w100 += S.hist[i][k];
        }
        char run[40], nxt[24], status[160];
        if (krun) snprintf(run, sizeof run, "%u (%" PRIu64 ")", REC[i].g - 2 * krun, S.hist[i][krun]);
        else if ((int)REC[i].g - 2 * (HK - 1) <= 2) snprintf(run, sizeof run, "none");
        else snprintf(run, sizeof run, "< %d", (int)REC[i].g - 2 * (HK - 1));
        if (hi == UINT64_MAX) snprintf(nxt, sizeof nxt, "(beyond 2^64)");
        else snprintf(nxt, sizeof nxt, "%" PRIu64, hi);
        if (full) {
            snprintf(status, sizeof status, "%s, %s", reps ? "TERM" : "not a term", oeis_tag(i + 1, reps > 0));
            if (i + 1 <= KNOWN_N[NKNOWN - 1] && (reps > 0) != (known_term_index(i + 1) > 0)) fails++;
        } else {
            u64 from = S.start > lo ? S.start : lo, to = pos < hi ? pos : hi;
            snprintf(status, sizeof status, "%s scanned [%" PRIu64 ", %" PRIu64 ") only", reps ? "TERM (repeat found);" : "undecided:", from, to);
        }
        printf(" %2d  %5u  %-21" PRIu64 "  %-21s  %-8" PRIu64 "  %-21s  %8" PRIu64 " / %-8" PRIu64 "  %s\n",
               i + 1, REC[i].g, lo, nxt, reps, run, w20, w100, status);
        if (reps) {
            printf("             repeats after:");
            for (unsigned j = 0; j < S.nrep[i]; j++) printf(" %" PRIu64, S.rep[i][j]);
            if (reps > S.nrep[i]) printf(" ... (%" PRIu64 " more)", reps - S.nrep[i]);
            putchar('\n');
        }
    }

    /* integrity: every record start in the range must have been seen */
    int nrec = 0, missing = 0;
    for (int i = 0; i < NREC; i++) {
        if (REC[i].p < S.start || REC[i].p >= pos) continue;
        nrec++;
        if (!S.recseen[i]) {
            missing++;
            printf("MISSING  record #%d (gap %u after %" PRIu64 ") was not seen where A002386 lists it\n", i + 1, REC[i].g, REC[i].p);
        }
    }
    printf("\n%d record gap%s in the range%s; %" PRIu64 " anomal%s\n", nrec, nrec == 1 ? "" : "s",
           missing ? "" : (nrec ? ", all confirmed where A002386 lists them" : ""), S.anomalies, S.anomalies == 1 ? "y" : "ies");
    fails += missing + (S.anomalies > 0);

    printf("%" PRIu64 " primes in [%" PRIu64 ", %" PRIu64 ")", S.primes, S.start, pos);
    for (int k = 6; k <= 18; k++) {
        if (pos != (u64)pow(10, k)) continue;
        u64 expect = PI10[k];
        if (S.start > 0) {
            if (S.start > 100000000000000ULL) break;          /* pi(START) would cost too much */
            expect -= primesieve_count_primes(0, S.start - 1);
        }
        printf(S.primes == expect ? "  = pi(10^%d)%s" : "  BUT pi(10^%d)%s = %" PRIu64 " -- MISMATCH", k,
               S.start ? " - pi(START-1)" : "", expect);
        if (S.primes != expect) fails++;
        break;
    }
    putchar('\n');

    if (g_hist_lines) {
        printf("\n# hist N GAP COUNT: gaps of size GAP after primes strictly inside record interval N (gap g_N first after A002386(N))\n");
        for (int i = 0; i < NREC; i++)
            for (int k = 0; k < HK; k++)
                if (S.hist[i][k]) printf("hist %d %u %" PRIu64 "\n", i + 1, REC[i].g - 2 * k, S.hist[i][k]);
    }
    fflush(stdout);
    return fails;
}

/* Brute-force reference for [0, limit): the same statistics from a plain loop
 * over an array of primes. */
typedef struct {
    u64 primes;
    u64 hist[NREC][HK];
    u64 rep[NREC][REPMAX];
    unsigned nrep[NREC];
    bool recseen[NREC];
    u64 anomalies;
} ref_t;

static void reference_scan(u64 limit, ref_t *r)
{
    size_t n;
    u64 *pr = primesieve_generate_primes(0, limit + LOOKAHEAD, &n, UINT64_PRIMES);
    if (!pr || n < 2) die("primesieve_generate_primes failed");
    memset(r, 0, sizeof *r);
    for (size_t j = 0; j + 1 < n && pr[j] < limit; j++) {
        const u64 p = pr[j], g = pr[j + 1] - p;
        r->primes++;
        int i = 0;
        while (i + 1 < NREC && REC[i + 1].p <= p) i++;
        if (p == REC[i].p) {
            if (g == REC[i].g) r->recseen[i] = true; else r->anomalies++;
            continue;
        }
        if (g > REC[i].g) { r->anomalies++; continue; }
        u64 k = (REC[i].g - g) / 2;
        if (k < HK) r->hist[i][k]++;
        if (g == REC[i].g && r->nrep[i] < REPMAX) r->rep[i][r->nrep[i]++] = p;
    }
    primesieve_free(pr);
}

static bool same_as_reference(const ref_t *r)
{
    if (r->primes != S.primes || r->anomalies != S.anomalies) return false;
    if (memcmp(r->hist, S.hist, sizeof S.hist) != 0) return false;
    if (memcmp(r->recseen, S.recseen, sizeof S.recseen) != 0) return false;
    if (memcmp(r->nrep, S.nrep, sizeof S.nrep) != 0) return false;
    for (int i = 0; i < NREC; i++)
        if (memcmp(r->rep[i], S.rep[i], S.nrep[i] * sizeof(u64)) != 0) return false;
    return true;
}

/* ------------------------------------------------------------------ */
/* Commands                                                            */
/* ------------------------------------------------------------------ */

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return n < 1 ? 1 : (n > 1024 ? 1024 : (int)n);
}

static void usage(void)
{
    fprintf(stderr,
        "usage: a053686 scan [START] END [-t T] [-c CHUNK] [-s KIB] [-n NEAR] [-S FILE] [-i SECS] [-H] [-q]\n"
        "       a053686 verify P G\n"
        "       a053686 bench [N] [-d SPAN] [-t T] [-c CHUNK] [-s KIB]\n"
        "       a053686 selftest [LIMIT] [-t T]\n"
        "START may be rN, the prime after which record gap N first occurs (e.g. r35 = 4302407359).\n");
    exit(2);
}

static int cmd_scan(int argc, char **argv)
{
    u64 pos[2];
    int npos = 0, threads = default_threads(), interval = 60;
    u64 chunk = 10000000000ULL;
    const char *state = NULL;
    bool quiet = false;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(argv[i], "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) g_sieve_kib = parse_int(argv[++i], 16, 65536, "sieve KiB");
        else if (!strcmp(argv[i], "-n") && i + 1 < argc) g_near = (unsigned)parse_int(argv[++i], 0, 2 * (HK - 1), "NEAR");
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = parse_int(argv[++i], 1, 86400, "SECS");
        else if (!strcmp(argv[i], "-H")) g_hist_lines = true;
        else if (!strcmp(argv[i], "-q")) quiet = true;
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos < 2) pos[npos++] = parse_num(argv[i]);
        else usage();
    }
    if (npos == 0) usage();
    u64 start = npos == 2 ? pos[0] : 0, end = pos[npos - 1];
    if (chunk < 1000000) die("CHUNK must be at least 10^6");
    if (g_sieve_kib) primesieve_set_sieve_size(g_sieve_kib);
    scan_stats_t st = run_scan(start, end, chunk, threads, state, interval, quiet);
    int fails = print_summary();
    char b[32];
    fprintf(stderr, "%" PRIu64 " numbers, %" PRIu64 " primes in %s (%.3g numbers/s)\n",
            st.numbers, S.primes, fmt_dur(st.seconds, b, sizeof b),
            st.seconds > 0 ? (double)st.numbers / st.seconds : 0);
    return fails ? 1 : 0;
}

static int cmd_verify(int argc, char **argv)
{
    if (argc < 2) usage();
    u64 p = parse_num(argv[0]), g = parse_num(argv[1]);
    if (!is_prime64(p)) die("%" PRIu64 " is not prime", p);
    u64 q = next_prime64(p);
    printf("%" PRIu64 " is prime; the next prime is %" PRIu64 ", gap %" PRIu64 "%s\n", p, q, q - p,
           q - p == g ? "" : "  -- NOT the claimed gap");
    bool ok = q - p == g;
#ifdef HAVE_GMP
    {
        mpz_t z;
        mpz_init_set_ui(z, p);
        mpz_nextprime(z, z);
        bool agree = mpz_cmp_ui(z, q) == 0;
        mpz_clear(z);
        printf("GMP mpz_nextprime %s\n", agree ? "agrees" : "DISAGREES");
        ok = ok && agree;
    }
#endif
    int i = rec_index(p);
    printf("record interval #%d: gap %u first after %" PRIu64 ", next record after ", i + 1, REC[i].g, REC[i].p);
    if (rec_end(i) == UINT64_MAX) printf("2^64 or beyond\n"); else printf("%" PRIu64 "\n", rec_end(i));
    if (p == REC[i].p)
        printf("this is the first occurrence of record #%d itself%s\n", i + 1, q - p == REC[i].g ? "" : " -- but the gap does not match the table!");
    else if (q - p == REC[i].g)
        printf("REPEAT: gap %u occurs again before the next record, so %u is a term of A053686  [%s]\n",
               REC[i].g, REC[i].g, oeis_tag(i + 1, true));
    else if (q - p > REC[i].g)
        printf("the gap exceeds the record of its interval: this would be an unlisted maximal gap (check the table)\n");
    else
        printf("not a repeat: the record of this interval is %u, short by %" PRIu64 "\n", REC[i].g, REC[i].g - (q - p));
    return ok ? 0 : 1;
}

static int cmd_bench(int argc, char **argv)
{
    u64 N = 1000000000000000ULL, span = 100000000000ULL, chunk = 10000000000ULL;
    int threads = default_threads(), npos = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) span = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-c") && i + 1 < argc) chunk = parse_num(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) g_sieve_kib = parse_int(argv[++i], 16, 65536, "sieve KiB");
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos++ == 0) N = parse_num(argv[i]);
        else usage();
    }
    if (chunk < 1000000) die("CHUNK must be at least 10^6");
    if (g_sieve_kib) primesieve_set_sieve_size(g_sieve_kib);
    g_quiet_events = true;
    scan_stats_t st = run_scan(N, N + span, chunk, threads, NULL, 60, true);
    g_quiet_events = false;
    double rate = (double)st.numbers / st.seconds;
    printf("bench: [%" PRIu64 ", %" PRIu64 ") with %d threads, chunk %" PRIu64 ": %.2f s\n",
           N, N + span, threads, chunk, st.seconds);
    printf("       %.4g numbers/s, %.4g primes/s (%" PRIu64 " primes)\n", rate, (double)st.primes / st.seconds, st.primes);
    const u64 s0 = REC[34].p;    /* 4302407359, where the OEIS entry stops */
    printf("       time to decide each record with a scan from %" PRIu64 " at this rate\n"
           "       (the sieve is faster below %.3g and slower above):\n", s0, (double)N);
    char b[32];
    for (int i = 34; i + 1 < NREC; i++) {
        double t = (double)(rec_end(i) - s0) / rate;
        if (t > 400 * 86400) break;
        printf("         #%-2d gap %-5u needs the scan to %-21" PRIu64 " %s\n", i + 1, REC[i].g, rec_end(i), fmt_dur(t, b, sizeof b));
    }
    return 0;
}

static int cmd_selftest(int argc, char **argv)
{
    u64 limit = 100000000000ULL;
    int threads = default_threads(), npos = 0, fails = 0;
    for (int i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) threads = parse_int(argv[++i], 1, 1024, "threads");
        else if (argv[i][0] == '-' && !isdigit((unsigned char)argv[i][1])) usage();
        else if (npos++ == 0) limit = parse_num(argv[i]);
        else usage();
    }
#define CHECK(cond, ...) do { if (cond) printf("ok    " __VA_ARGS__); else { printf("FAIL  " __VA_ARGS__); fails++; } putchar('\n'); } while (0)

    CHECK(parse_num("2^54-32") == 18014398509481952ULL && parse_num("1e12") == 1000000000000ULL &&
          parse_num("1.7e15") == 1700000000000000ULL && parse_num("10^15+7") == 1000000000000007ULL &&
          parse_num("r35") == 4302407359ULL,
          "number parsing");
    CHECK(is_prime64(2) && is_prime64(3) && !is_prime64(1) && !is_prime64(9) && is_prime64(1000000007ULL) &&
          !is_prime64(3215031751ULL) && is_prime64(18446744073709551557ULL) && !is_prime64(18446744073709551555ULL) &&
          is_prime64(4302407359ULL) && is_prime64(4302407359ULL + 354) && !is_prime64(4302407359ULL + 2),
          "Miller-Rabin spot checks");

    /* the table: increasing, each entry a genuine prime gap of the listed size */
    {
        bool mono = true, gaps = true;
        for (int i = 0; i < NREC; i++) {
            if (i && (REC[i].g <= REC[i - 1].g || REC[i].p <= REC[i - 1].p)) mono = false;
            if (!is_prime64(REC[i].p) || next_prime64(REC[i].p) != REC[i].p + REC[i].g) {
                gaps = false;
                printf("      table entry #%d (%u after %" PRIu64 ") is not a prime gap of that size\n", i + 1, REC[i].g, REC[i].p);
            }
        }
        CHECK(mono, "A005250/A002386 table is increasing in both columns");
        CHECK(gaps, "all %d table entries are prime gaps of the listed size (Miller-Rabin next-prime search)", NREC);
        CHECK(rec_index(2) == 0 && rec_index(3) == 1 && rec_index(6) == 1 && rec_index(7) == 2 && rec_index(22) == 2 &&
              rec_index(23) == 3 && rec_index(4302407358ULL) == 33 && rec_index(4302407359ULL) == 34 &&
              rec_index(UINT64_MAX) == NREC - 1, "record interval lookup");
    }

    /* chunk boundaries: odd chunk sizes and thread counts against the brute-force reference */
    {
        const u64 lim = 200000000ULL;
        ref_t *ref = malloc(sizeof *ref);
        if (!ref) die("out of memory");
        static const struct { u64 chunk; int threads; } cfg[] = {
            { 1234567ULL, 0 }, { 7000001ULL, 3 }, { 50000000ULL, 1 }, { 200000000ULL, 1 }, { 1000000ULL, 0 }
        };
        unsigned save_near = g_near;
        g_near = 0;
        for (size_t ci = 0; ci < sizeof cfg / sizeof *cfg; ci++) {
            int t = cfg[ci].threads ? cfg[ci].threads : threads;
            g_quiet_events = true;                  /* no REPEAT/RECORD lines for these small scans */
            run_scan(0, lim, cfg[ci].chunk, t, NULL, 60, true);
            g_quiet_events = false;
            reference_scan(S.end, ref);
            u64 tot = 0;
            for (int i = 0; i < NREC; i++) for (int k = 0; k < HK; k++) tot += S.hist[i][k];
            CHECK(same_as_reference(ref), "chunked scan (chunk %" PRIu64 ", %d threads) matches brute force on [0, %" PRIu64
                  "): %" PRIu64 " primes, %" PRIu64 " histogram entries, %d records", cfg[ci].chunk, t, S.end, S.primes, tot, S.last_rec);
            if (S.end == lim) CHECK(S.primes == 11078937ULL, "prime count is pi(2*10^8) = 11078937");
        }
        g_near = save_near;
        free(ref);
    }

    printf("sieve scan of [0, %" PRIu64 ") with %d threads ...\n", limit, threads);
    fflush(stdout);
    u64 chunk = limit >= 100000000000ULL ? 10000000000ULL : (limit >= 10000000000ULL ? 1000000000ULL : 100000000ULL);
    scan_stats_t st = run_scan(0, limit, chunk, threads, NULL, 60, true);
    char b[32];
    printf("      %" PRIu64 " primes in %s (%.3g numbers/s)\n", st.primes, fmt_dur(st.seconds, b, sizeof b),
           st.seconds > 0 ? (double)st.numbers / st.seconds : 0);
    int ndecided = 0;
    for (int i = 0; i + 1 < NREC && REC[i + 1].p <= S.end; i++) {
        const int n = i + 1;
        const bool reps = S.hist[i][0] > 0;
        ndecided++;
        if (n <= KNOWN_N[NKNOWN - 1]) {
            int k = known_term_index(n);
            CHECK(reps == (k > 0), "record #%d = %u %s before the next record (OEIS: %s)", n, REC[i].g,
                  reps ? "repeats" : "does not repeat", k ? "a term" : "not a term");
        } else if (n <= A085237_LAST) {
            CHECK(!reps, "record #%d = %u does not repeat (A085237 lists it once)", n, REC[i].g);
        } else {
            printf("      record #%d = %u %s before the next record (not decided in OEIS)\n", n, REC[i].g, reps ? "REPEATS" : "does not repeat");
        }
        if (reps)
            CHECK(next_prime64(S.rep[i][0]) == S.rep[i][0] + REC[i].g, "first repeat of %u after %" PRIu64 " confirmed by next-prime search",
                  REC[i].g, S.rep[i][0]);
    }
    CHECK(ndecided >= 34 || limit < REC[34].p, "%d records decided below %" PRIu64, ndecided, S.end);
    {
        int missing = 0;
        for (int i = 0; i < NREC && REC[i].p < S.end; i++) if (!S.recseen[i]) missing++;
        CHECK(missing == 0 && S.anomalies == 0, "every record below the limit found where A002386 lists it, no anomalies");
    }
    for (int k = 6; k <= 18; k++)
        if (S.end == (u64)pow(10, k)) CHECK(S.primes == PI10[k], "prime count %" PRIu64 " = pi(10^%d)", S.primes, k);
#undef CHECK
    printf(fails ? "SELFTEST FAILED (%d)\n" : "selftest passed\n", fails);
    return fails ? 1 : 0;
}

int main(int argc, char **argv)
{
    stderr_tty = isatty(2);
    setvbuf(stdout, NULL, _IOLBF, 0);
    if (argc < 2) usage();
    if (!strcmp(argv[1], "scan")) return cmd_scan(argc - 2, argv + 2);
    if (!strcmp(argv[1], "verify")) return cmd_verify(argc - 2, argv + 2);
    if (!strcmp(argv[1], "bench")) return cmd_bench(argc - 2, argv + 2);
    if (!strcmp(argv[1], "selftest")) return cmd_selftest(argc - 2, argv + 2);
    usage();
    return 2;
}
