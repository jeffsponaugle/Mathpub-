/*
 * pdigits.c — count decimal digit occurrences in the first 10^N primes.
 *
 * Produces the next term(s) of OEIS A119291..A119300 (digits 0..9) in a
 * single run, with cross-verification against:
 *   - A119291..A119300: all published per-digit counts (n = 1..12)
 *   - A006988: the 10^n-th prime at every boundary
 *   - A119290: total digits in first 10^n primes (known through n=21,
 *              so it also checks the NEW n=13/14 digit sums)
 *   - A006880: pi(10^d) prime counts passed along the way
 *
 * Method: the range [0, ~p_{10^N}] is cut into span-aligned chunks
 * (default 10^9). Worker threads pull chunk indices and tally prime count
 * plus per-digit counts per chunk with a primesieve iterator. Digit
 * tallying uses SWAR lookup tables over the low 7 digits (3+4 split) with
 * the shared high-digit prefix counted once per 10^7 window. Completed
 * chunks are appended to the checkpoint log; resume re-reads it and redoes
 * only missing chunks. When the contiguous chunk prefix contains 10^N
 * primes, exact boundary crossings are located by rescanning the (single)
 * chunk containing each 10^n-th prime.
 *
 * Everything is exact 64-bit integer arithmetic; no probabilistic steps.
 *
 * Build: make   (needs primesieve; brew install primesieve)
 */

#include <inttypes.h>
#include <pthread.h>
#include <stdarg.h>
#include <signal.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

#include <primesieve.h>

typedef uint32_t u32;
typedef uint64_t u64;

/* ------------------------------------------------------------------ */
/* reference data from OEIS                                            */
/* ------------------------------------------------------------------ */

/* A006988: 10^n-th prime, n = 0..18 */
static const u64 P10N[] = { 2, 29, 541, 7919, 104729, 1299709, 15485863,
    179424673, 2038074743ULL, 22801763489ULL, 252097800623ULL,
    2760727302517ULL, 29996224275833ULL, 323780508946331ULL,
    3475385758524527ULL, 37124508045065437ULL, 394906913903735329ULL,
    4185296581467695669ULL, 0 /* 44211790234832169331 > 2^64? no: fits, but unused */ };
#define NP10N 18

/* A006880: pi(10^d), d = 1..15 */
static const u64 PI10D[] = { 4, 25, 168, 1229, 9592, 78498, 664579,
    5761455, 50847534, 455052511ULL, 4118054813ULL, 37607912018ULL,
    346065536839ULL, 3204941750802ULL, 29844570422669ULL,
    279238341033925ULL, 2623557157654233ULL, 24739954287740860ULL };
#define NPI10D 18

/* A119290: total digits in first 10^n primes, n = 0..15 (u64-safe part) */
static const u64 TOTD[] = { 1, 16, 271, 3803, 48982, 610484, 7245905,
    83484450, 942636916ULL, 10487584405ULL, 115369529592ULL,
    1257761617574ULL, 13611696080735ULL, 146406754329933ULL,
    1566562183907264ULL, 16687323842873339ULL, 177063766685219106ULL };
#define NTOTD 17

/* A119291..A119300: count of digit d in first 10^n primes, n = 1..12 */
static const u64 DIGN[10][12] = {
 {0,9,191,3303,46188,557005,6481183,76292782,881025347ULL,9763247930ULL,106864564286ULL,1162019145892ULL},
 {5,55,574,7043,102370,1222003,13896979,152844768,1611113572ULL,16626301609ULL,171878734261ULL,1781706224877ULL},
 {3,26,339,4070,55213,632418,7133747,82051293,1041785731ULL,12182327373ULL,137771632675ULL,1525041000685ULL},
 {3,54,569,6500,72717,880415,9616078,105052677,1139295891ULL,12676601300ULL,135996846047ULL,1450287630358ULL},
 {0,27,311,3919,47647,628061,7098709,79587350,888504625ULL,10169426048ULL,110936566287ULL,1199775766812ULL},
 {1,15,327,3904,47525,590450,7087898,79504457,887852429ULL,9862623240ULL,110885914721ULL,1199346258292ULL},
 {0,10,315,3824,47269,560677,7079879,79433407,887292243ULL,9778050835ULL,110518283071ULL,1198946387631ULL},
 {2,34,551,6338,72319,809360,9543704,104376285,1136782466ULL,12273965395ULL,134080968533ULL,1448607569210ULL},
 {0,8,195,3763,47174,558842,6541723,79273100,882188472ULL,9770202402ULL,106927984586ULL,1198303968170ULL},
 {2,33,431,6318,72062,806674,9004550,104220797,1131743629ULL,12266783460ULL,131900123107ULL,1447662128808ULL},
};
#define NDIGN 12

static const char *SEQID[10] = { "A119291","A119292","A119293","A119294",
    "A119295","A119296","A119297","A119298","A119299","A119300" };

/* A091634..A091643: primes < 10^n NOT containing digit d, n = 1..19 */
static const u64 NOD[10][19] = {
 {4,25,153,1010,7122,52313,397866,3103348,24649318ULL,198536215ULL,1616808581ULL,13287264748ULL,110033428309ULL,917072930187ULL,7685458755706ULL,64714772654036ULL,547198653032787ULL,4643886890826448ULL,39539902821647035ULL},
 {4,17,101,670,4675,34425,262549,2051466,16312743ULL,131464721ULL,1071368863ULL,8809580516ULL,72986908554ULL,608542410004ULL,5101540409471ULL,42969395328824ULL,363420769121973ULL,3084910521042831ULL,26271330977102835ULL},
 {3,22,139,877,6235,46105,352155,2747284,21831323ULL,175881412ULL,1432781905ULL,11778245565ULL,97558533214ULL,813253056497ULL,6816503587850ULL,57405945577835ULL,485459545302540ULL,4120383980056923ULL,35086052977411238ULL},
 {3,16,102,668,4715,34813,265015,2067152,16413535ULL,132200223ULL,1076692515ULL,8849480283ULL,73288053795ULL,610860050965ULL,5119536757406ULL,43110730041621ULL,364540518319108ULL,3093850266657393ULL,26343191757392667ULL},
 {4,22,136,903,6361,46545,354123,2761106,21925170ULL,176544507ULL,1437663500ULL,11814853749ULL,97837428598ULL,815398741896ULL,6833236940218ULL,57537517496756ULL,486502895023133ULL,4128723230358420ULL,35153146157359046ULL},
 {3,22,136,905,6310,46549,354910,2765749,21955845ULL,176781643ULL,1439380189ULL,11827571824ULL,97933795005ULL,816144146010ULL,6839035710792ULL,57583214837091ULL,486865837479563ULL,4131626542889840ULL,35176529911304840ULL},
 {4,23,136,897,6367,46706,355148,2770239,21984207ULL,176966593ULL,1440765209ULL,11838096715ULL,98014747908ULL,816769206831ULL,6843920386961ULL,57621676585577ULL,487171557272477ULL,4134074824797519ULL,35196259497882848ULL},
 {3,16,100,680,4773,34992,266823,2079512,16503238ULL,132852644ULL,1081509855ULL,8885472675ULL,73563855306ULL,612982476612ULL,5136111809120ULL,43241276000455ULL,365577774983773ULL,3102150168920636ULL,26410045901566231ULL},
 {4,23,141,915,6375,46799,355805,2774348,22023132ULL,177273427ULL,1443074791ULL,11855541525ULL,98146301284ULL,817786989282ULL,6851887572613ULL,57684575826536ULL,487671951877316ULL,4138082681079015ULL,35228572737080252ULL},
 {4,19,108,687,4766,35139,267486,2083814,16531372ULL,133059504ULL,1082995490ULL,8896945667ULL,73651718719ULL,613664827254ULL,5141432448369ULL,43283348254700ULL,365912414930465ULL,3104832036379344ULL,26431670649900716ULL},
};
#define NNOD 19

/* ------------------------------------------------------------------ */
/* SWAR digit tables: byte lanes 0..7 = digits 0..7 in u64,            */
/* bits 0-7 / 16-23 = digits 8,9 in u32                                */
/* ------------------------------------------------------------------ */

static u64 t3lo[1000], t4lo[10000];
static u32 t3hi[1000], t4hi[10000];
static uint16_t t3mask[1000], t4mask[10000];   /* digit-presence bitmasks */

static void build_tables(void)
{
    for (int v = 0; v < 10000; v++) {
        u64 lo = 0; u32 hi = 0;
        int a = v / 1000, b = v / 100 % 10, c = v / 10 % 10, d = v % 10;
        int dig[4] = { a, b, c, d };
        for (int i = 0; i < 4; i++) {
            if (dig[i] < 8) lo += 1ULL << (8 * dig[i]);
            else            hi += 1u << (16 * (dig[i] - 8));
        }
        t4lo[v] = lo; t4hi[v] = hi;
        t4mask[v] = (uint16_t)((1u << a) | (1u << b) | (1u << c) | (1u << d));
        if (v < 1000) {
            lo = 0; hi = 0;
            for (int i = 1; i < 4; i++) {   /* skip thousands digit */
                if (dig[i] < 8) lo += 1ULL << (8 * dig[i]);
                else            hi += 1u << (16 * (dig[i] - 8));
            }
            t3lo[v] = lo; t3hi[v] = hi;
            t3mask[v] = (uint16_t)((1u << b) | (1u << c) | (1u << d));
        }
    }
}

static void slow_digits(u64 v, u64 d[10])   /* no leading zeros */
{
    if (v == 0) return;
    while (v) { d[v % 10]++; v /= 10; }
}

/* ------------------------------------------------------------------ */
/* configuration & state                                               */
/* ------------------------------------------------------------------ */

typedef struct { u64 np; u64 d[10]; u64 nod[10]; } tally_t;

/* fold a presence-mask histogram into "lacks digit d" counts */
static void hist_to_nod(const u64 *hist, u64 nod[10])
{
    for (int m = 0; m < 1024; m++) {
        u64 c = hist[m];
        if (!c) continue;
        for (int i = 0; i < 10; i++)
            if (!(m & (1 << i))) nod[i] += c;
    }
}

static struct {
    int N;                       /* count first 10^N primes */
    int span_e;                  /* chunk span = 10^span_e  */
    u64 span;
    int threads;
    double status_iv;
    const char *ckpt_file;
    const char *report_file;
    FILE *rf;
    int resume, quiet, selftest;

    u64 target;                  /* 10^N */
    u64 pmax;                    /* planned sieve limit */
    u64 nplan;                   /* planned chunk count */
    u64 chunk_min, chunk_max;    /* this machine's assignment [min,max) */
    tally_t *res;
    unsigned char *done;
    _Atomic u64 next_chunk;
    _Atomic u64 chunks_done;
    _Atomic u64 done_in_range;
    _Atomic int stop;
    pthread_mutex_t lock;        /* log file + frontier */
    FILE *log;
    u64 frontier;                /* chunks [0,frontier) all done */
    u64 cum_np_frontier;
    u64 pi_small[20];            /* pi(10^d) for 10^d inside chunk 0 */
    u64 nod_small[12][10];       /* no-digit-d counts at 10^d in chunk 0 */
    int pi_small_max;
    double t_start;
} S;

static volatile sig_atomic_t g_sigint = 0, g_siginfo = 0;
static void on_sigint(int sig) { (void)sig; g_sigint = 1; }
#ifdef SIGINFO
static void on_siginfo(int sig) { (void)sig; g_siginfo = 1; }
#endif

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static u64 pow10u(int e) { u64 v = 1; while (e-- > 0) v *= 10; return v; }

/* thousands separators for readability */
static char *sep_u64(u64 v, char *buf)
{
    char tmp[32]; int n = 0;
    if (!v) tmp[n++] = '0';
    while (v) { tmp[n++] = (char)('0' + v % 10); v /= 10; }
    int j = 0;
    for (int i = n - 1; i >= 0; i--) {
        buf[j++] = tmp[i];
        if (i && i % 3 == 0) buf[j++] = ',';
    }
    buf[j] = 0;
    return buf;
}

/* ------------------------------------------------------------------ */
/* chunk tallying                                                      */
/* ------------------------------------------------------------------ */

/* exact per-prime path; also records pi(10^d) and no-digit-d counts at
 * the 10^d crossings (chunk 0) */
static void tally_slow(primesieve_iterator *it, u64 start, u64 stop,
                       tally_t *t, int record_pi)
{
    int pd = 1;
    u64 next_pow = 10;
    u64 *hist = calloc(1024, sizeof(u64));
    primesieve_jump_to(it, start < 2 ? 2 : start, stop);
    for (u64 p = primesieve_next_prime(it); p < stop;
         p = primesieve_next_prime(it)) {
        if (record_pi) {
            while (p > next_pow && pd <= 19) {
                if (pd <= (int)(sizeof S.pi_small / sizeof S.pi_small[0])) {
                    S.pi_small[pd - 1] = t->np;
                    if (pd > S.pi_small_max) S.pi_small_max = pd;
                }
                if (pd <= 12) {
                    memset(S.nod_small[pd - 1], 0, sizeof S.nod_small[0]);
                    hist_to_nod(hist, S.nod_small[pd - 1]);
                }
                pd++;
                next_pow = (pd >= 20) ? UINT64_MAX : next_pow * 10;
            }
        }
        t->np++;
        u64 v = p; unsigned m = 0;
        while (v) { t->d[v % 10]++; m |= 1u << (v % 10); v /= 10; }
        hist[m]++;
    }
    hist_to_nod(hist, t->nod);
    free(hist);
}

/* fast SWAR path; requires start >= 10^7 so every prefix is nonzero.
 * deadline (monotonic seconds, 0 = none) stops after the current 10^7
 * window — used only by --benchmark. */
static void tally_fast(primesieve_iterator *it, u64 start, u64 stop,
                       tally_t *t, u64 *hist /* zeroed u64[1024] */,
                       double deadline)
{
    const u64 W = 10000000ULL;
    primesieve_jump_to(it, start, stop);
    u64 p = primesieve_next_prime(it);
    while (p < stop) {
        if (deadline != 0 && now_s() >= deadline) break;
        u64 pref = p / W;
        u64 base = pref * W;
        u64 wend = base + W;
        if (wend > stop) wend = stop;
        u64 pd[10] = {0};
        slow_digits(pref, pd);
        unsigned pmask = 0;
        { u64 v = pref; while (v) { pmask |= 1u << (v % 10); v /= 10; } }
        u64 npw = 0;
        u64 alo = 0; u32 ahi = 0; int pending = 0;
        while (p < wend) {
            u32 v = (u32)(p - base);           /* < 10^7 */
            u32 a = v / 1000, b = v - a * 1000;
            alo += t4lo[a] + t3lo[b];
            ahi += t4hi[a] + t3hi[b];
            hist[pmask | t4mask[a] | t3mask[b]]++;
            npw++;
            if (++pending == 30) {
                for (int i = 0; i < 8; i++) t->d[i] += (alo >> (8 * i)) & 0xFF;
                t->d[8] += ahi & 0xFF; t->d[9] += (ahi >> 16) & 0xFF;
                alo = 0; ahi = 0; pending = 0;
            }
            p = primesieve_next_prime(it);
        }
        if (pending) {
            for (int i = 0; i < 8; i++) t->d[i] += (alo >> (8 * i)) & 0xFF;
            t->d[8] += ahi & 0xFF; t->d[9] += (ahi >> 16) & 0xFF;
        }
        for (int i = 0; i < 10; i++) t->d[i] += pd[i] * npw;
        t->np += npw;
    }
    hist_to_nod(hist, t->nod);
}

static void *worker(void *arg)
{
    (void)arg;
    primesieve_iterator it;
    primesieve_init(&it);
    u64 *hist = malloc(1024 * sizeof(u64));
    for (;;) {
        if (atomic_load_explicit(&S.stop, memory_order_relaxed)) break;
        u64 idx = atomic_fetch_add(&S.next_chunk, 1);
        if (idx >= S.chunk_max) break;
        if (S.done[idx]) continue;

        tally_t t; memset(&t, 0, sizeof t);
        u64 start = idx * S.span, stop = (idx + 1) * S.span;
        if (idx == 0) tally_slow(&it, 0, stop, &t, 1);
        else {
            memset(hist, 0, 1024 * sizeof(u64));
            tally_fast(&it, start, stop, &t, hist, 0);
        }

        pthread_mutex_lock(&S.lock);
        S.res[idx] = t;
        S.done[idx] = 1;
        if (S.log) {
            fprintf(S.log, "C %" PRIu64 " %" PRIu64, idx, t.np);
            for (int i = 0; i < 10; i++)
                fprintf(S.log, " %" PRIu64, t.d[i]);
            for (int i = 0; i < 10; i++)
                fprintf(S.log, " %" PRIu64, t.nod[i]);
            fprintf(S.log, "\n");
            fflush(S.log);
        }
        while (S.frontier < S.nplan && S.done[S.frontier]) {
            S.cum_np_frontier += S.res[S.frontier].np;
            S.frontier++;
        }
        if (S.cum_np_frontier >= S.target)
            atomic_store(&S.stop, 1);
        pthread_mutex_unlock(&S.lock);
        atomic_fetch_add(&S.chunks_done, 1);
        atomic_fetch_add(&S.done_in_range, 1);
    }
    free(hist);
    primesieve_free_iterator(&it);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* checkpoint                                                          */
/* ------------------------------------------------------------------ */

static int load_checkpoint(void)
{
    FILE *f = fopen(S.ckpt_file, "r");
    if (!f) { fprintf(stderr, "error: cannot open %s\n", S.ckpt_file); return 0; }
    char line[512];
    int N = 0, span_e = 0;
    if (!fgets(line, sizeof line, f) ||
        sscanf(line, "pdigits-v2 N=%d span_e=%d", &N, &span_e) != 2 ||
        N != S.N || span_e != S.span_e) {
        fprintf(stderr, "error: checkpoint mismatch (want N=%d span_e=%d)\n",
                S.N, S.span_e);
        fclose(f);
        return 0;
    }
    u64 loaded = 0;
    while (fgets(line, sizeof line, f)) {
        u64 idx, np, d[10], nod[10];
        int k = sscanf(line, "C %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64
                       " %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64
                       " %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64
                       " %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64
                       " %" SCNu64 " %" SCNu64 " %" SCNu64 " %" SCNu64
                       " %" SCNu64 " %" SCNu64,
                       &idx, &np, &d[0], &d[1], &d[2], &d[3], &d[4],
                       &d[5], &d[6], &d[7], &d[8], &d[9],
                       &nod[0], &nod[1], &nod[2], &nod[3], &nod[4],
                       &nod[5], &nod[6], &nod[7], &nod[8], &nod[9]);
        if (k != 22 || idx >= S.nplan) continue;
        S.res[idx].np = np;
        memcpy(S.res[idx].d, d, sizeof d);
        memcpy(S.res[idx].nod, nod, sizeof nod);
        S.done[idx] = 1;
        loaded++;
    }
    fclose(f);
    while (S.frontier < S.nplan && S.done[S.frontier]) {
        S.cum_np_frontier += S.res[S.frontier].np;
        S.frontier++;
    }
    atomic_store(&S.chunks_done, loaded);
    printf("resumed: %" PRIu64 " chunks loaded, frontier at chunk %" PRIu64
           " (%.3e primes)\n", loaded, S.frontier, (double)S.cum_np_frontier);
    /* chunk 0 pi_small was not persisted; recompute quickly if needed */
    if (S.done[0]) {
        primesieve_iterator it; primesieve_init(&it);
        tally_t t; memset(&t, 0, sizeof t);
        tally_slow(&it, 0, S.span, &t, 1);
        primesieve_free_iterator(&it);
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* final report                                                        */
/* ------------------------------------------------------------------ */

/* print to stdout and, if configured, to the report file */
static void rp(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    if (S.rf) {
        va_start(ap, fmt);
        vfprintf(S.rf, fmt, ap);
        va_end(ap);
    }
}

static int report(void)
{
    char b1[32], b2[32];
    int all_ok = 1;

    if (S.report_file) {
        S.rf = fopen(S.report_file, "w");
        if (!S.rf)
            fprintf(stderr, "warn: cannot open report file %s\n", S.report_file);
        else {
            time_t tt = time(NULL);
            char ts[64];
            strftime(ts, sizeof ts, "%Y-%m-%d %H:%M:%S", localtime(&tt));
            fprintf(S.rf, "pdigits report  N=%d  %s\n", S.N, ts);
        }
    }

    /* prefix sums over chunks */
    u64 nch = S.frontier;
    u64 *cumnp = malloc((nch + 1) * sizeof(u64));
    tally_t cum; memset(&cum, 0, sizeof cum);
    cumnp[0] = 0;
    for (u64 i = 0; i < nch; i++) cumnp[i + 1] = cumnp[i] + S.res[i].np;

    if (cumnp[nch] < S.target) {
        rp("incomplete: %.3e of %.3e primes counted; no report.\n",
           (double)cumnp[nch], (double)S.target);
        free(cumnp);
        if (S.rf) fclose(S.rf);
        return 0;
    }

    /* pi(10^d) verification */
    rp("\npi(10^d) checks vs A006880:\n");
    for (int d = 1; d <= NPI10D; d++) {
        u64 pd = pow10u(d), got;
        if (d <= S.pi_small_max) got = S.pi_small[d - 1];
        else if (pd % S.span == 0 && pd / S.span <= nch)
            got = cumnp[pd / S.span];
        else continue;
        int ok = got == PI10D[d - 1];
        if (!ok) all_ok = 0;
        rp("  pi(10^%-2d) = %-18s %s\n", d, sep_u64(got, b1),
           ok ? "OK" : "** MISMATCH **");
    }

    /* value-boundary family: primes < 10^n lacking digit d (A091634..43) */
    rp("\nprimes < 10^n not containing digit d, vs A091634..A091643:\n");
    for (int n = 1; n <= 19; n++) {
        u64 pn = pow10u(n);
        u64 nod[10]; int have = 0;
        if (n <= S.pi_small_max && pn < S.span) {
            memcpy(nod, S.nod_small[n - 1], sizeof nod);
            have = 1;
        } else if (pn % S.span == 0 && pn / S.span <= nch) {
            memset(nod, 0, sizeof nod);
            for (u64 j = 0; j < pn / S.span; j++)
                for (int i = 0; i < 10; i++) nod[i] += S.res[j].nod[i];
            have = 1;
        }
        if (!have) break;
        int row_ok = 1, row_new = 0;
        char row[512]; int len = 0;
        for (int i = 0; i < 10; i++) {
            int ok = n <= NNOD ? nod[i] == NOD[i][n - 1] : -1;
            if (ok == 0) { row_ok = 0; all_ok = 0; }
            if (ok < 0) row_new = 1;
            len += snprintf(row + len, sizeof row - len, "%s%" PRIu64 "%s",
                            i ? " " : "", nod[i], ok == 0 ? "!" : "");
        }
        rp("n=%-2d %s  [%s]\n", n, row,
           !row_ok ? "** MISMATCH (! marks) **"
                   : row_new ? "NEW TERMS" : "all OK vs A091634..A091643");
    }

    /* boundary rescans, in increasing n: walk chunks with running tally */
    rp("\nper-boundary results:\n");
    primesieve_iterator it; primesieve_init(&it);
    for (int n = 1; n <= S.N; n++) {
        u64 T = pow10u(n);
        /* find chunk c with cumnp[c] < T <= cumnp[c+1] */
        u64 lo = 0, hi = nch;
        while (lo + 1 < hi) {
            u64 mid = (lo + hi) / 2;
            if (cumnp[mid] < T) lo = mid; else hi = mid;
        }
        u64 c = lo;
        tally_t part; memset(&part, 0, sizeof part);
        u64 need = T - cumnp[c], np = 0, bp = 0;
        primesieve_jump_to(&it, c ? c * S.span : 2, (c + 1) * S.span);
        for (u64 p = primesieve_next_prime(&it); np < need;
             p = primesieve_next_prime(&it)) {
            np++;
            slow_digits(p, part.d);
            bp = p;
        }
        /* totals at T = full chunks below c + partial */
        u64 dig[10], tot = 0;
        for (int i = 0; i < 10; i++) {
            dig[i] = part.d[i];
            for (u64 j = 0; j < c; j++) dig[i] += S.res[j].d[i];
            tot += dig[i];
        }
        /* (re-walk is O(nch) per digit per n; fine at these sizes) */

        int okp = (n < NP10N && P10N[n]) ? bp == P10N[n] : -1;
        int okt = n < NTOTD ? tot == TOTD[n] : -1;
        if (okp == 0 || okt == 0) all_ok = 0;
        rp("n=%-2d p(10^n) = %-22s %s\n", n, sep_u64(bp, b1),
           okp < 0 ? "(beyond A006988 table)" : okp ? "OK A006988" : "** MISMATCH A006988 **");
        rp("     digits  = %-22s %s\n", sep_u64(tot, b2),
           okt < 0 ? "(beyond A119290 table)" : okt ? "OK A119290" : "** MISMATCH A119290 **");
        for (int i = 0; i < 10; i++) {
            int ok = n <= NDIGN ? dig[i] == DIGN[i][n - 1] : -1;
            if (!ok) all_ok = 0;
            rp("     d%d = %-22s %s%s\n", i, sep_u64(dig[i], b1),
               SEQID[i],
               ok < 0 ? "  ** NEW TERM **" : ok ? " OK" : " ** MISMATCH **");
        }
    }
    primesieve_free_iterator(&it);

    rp("\noverall verification: %s\n",
       all_ok ? "ALL CHECKS PASSED" : "** FAILURES PRESENT **");
    free(cumnp);
    if (S.rf) {
        fclose(S.rf);
        S.rf = NULL;
        printf("report written to %s\n", S.report_file);
    }
    return all_ok;
}

/* ------------------------------------------------------------------ */
/* benchmark                                                           */
/* ------------------------------------------------------------------ */

/* Fixed workload for cross-machine comparison: every thread tallies its
 * own contiguous stripe starting at 10^17 (representative of an N=16
 * run) until the shared deadline. Score = primes tallied per second. */

#define BENCH_BASE   100000000000000000ULL   /* 10^17 */
#define BENCH_STRIDE 1000000000000ULL        /* 10^12 per thread */

typedef struct { int id; double deadline; u64 primes; } bworker_t;

static void *bench_worker(void *arg)
{
    bworker_t *b = arg;
    primesieve_iterator it;
    primesieve_init(&it);
    u64 *hist = calloc(1024, sizeof(u64));
    tally_t t;
    memset(&t, 0, sizeof t);
    u64 start = BENCH_BASE + (u64)b->id * BENCH_STRIDE;
    tally_fast(&it, start, start + BENCH_STRIDE, &t, hist, b->deadline);
    b->primes = t.np;
    free(hist);
    primesieve_free_iterator(&it);
    return NULL;
}

static int benchmark(double secs)
{
    printf("pdigits benchmark: %d threads, %.0f s, region 10^17 "
           "(representative of an N=16 run)\n", S.threads, secs);
    fflush(stdout);
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)S.threads);
    bworker_t *bw = calloc((size_t)S.threads, sizeof(bworker_t));
    double t0 = now_s();
    for (int i = 0; i < S.threads; i++) {
        bw[i].id = i;
        bw[i].deadline = t0 + secs;
        pthread_create(&th[i], NULL, bench_worker, &bw[i]);
    }
    u64 tot = 0;
    for (int i = 0; i < S.threads; i++) {
        pthread_join(th[i], NULL);
        tot += bw[i].primes;
    }
    double el = now_s() - t0;
    double rate = (double)tot / el;
    printf("  %.4e primes tallied in %.1f s (%.3e primes/s/thread)\n",
           (double)tot, el, rate / S.threads);
    printf("\nSCORE: %.3f Gprimes/s  (%d threads)\n", rate / 1e9, S.threads);
    printf("\nrough full-run projections at this rate:\n");
    printf("  N=13 (10^13 primes): %6.1f hours\n", 1e13 / rate / 3600.0);
    printf("  N=15 (10^15 primes): %6.1f days\n",  1e15 / rate / 86400.0);
    printf("  N=16 (10^16 primes): %6.1f days\n",  1e16 / rate / 86400.0);
    printf("(projections ignore height variation, ~+-20%%; the SCORE line is\n"
           " the number to compare between machines / -t settings)\n");
    free(th); free(bw);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void usage(void)
{
    printf(
"pdigits — digit counts in the first 10^N primes (OEIS A119291..A119300)\n\n"
"usage: pdigits [options]\n"
"  -n N            target: first 10^N primes (default 13, max 16)\n"
"  -t THREADS      worker threads (default: online cores)\n"
"  -s E            chunk span 10^E (default 9; 10 recommended for N=14+)\n"
"  -c FILE         checkpoint log (default pdigits.ckpt; 'none' disables)\n"
"  -o FILE         also write the final report/summary to FILE\n"
"  -r              resume from checkpoint\n"
"  -u SECS         status interval (default 10)\n"
"  -q              quiet\n"
"  --selftest      full pipeline at N=8 vs published values\n"
"  --benchmark[=S] measure this machine's speed for S seconds (default 45)\n"
"                  and print a SCORE in Gprimes/s; honors -t\n"
"  --chunk-min I   only process chunks >= I (multi-machine split)\n"
"  --chunk-max I   only process chunks < I; merge logs afterward and\n"
"                  rerun with -r on the combined file for the report\n"
"  -h              help\n\n"
"Ctrl-C exits cleanly (checkpoint log is already on disk); resume with -r.\n"
"Ctrl-T (SIGINFO) prints an immediate status line.\n");
}

int main(int argc, char **argv)
{
    memset(&S, 0, sizeof S);
    S.N = 13; S.span_e = 9;
    S.threads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    S.status_iv = 10;
    S.ckpt_file = "pdigits.ckpt";

    static struct option lo[] = {
        { "selftest", 0, 0, 1 }, { "chunk-min", 1, 0, 2 },
        { "chunk-max", 1, 0, 3 }, { "benchmark", 2, 0, 4 },
        { "help", 0, 0, 'h' }, { 0, 0, 0, 0 } };
    double bench_secs = 0;
    int c;
    while ((c = getopt_long(argc, argv, "n:t:s:c:o:ru:qh", lo, NULL)) != -1) {
        switch (c) {
        case 'n': S.N = atoi(optarg); break;
        case 't': S.threads = atoi(optarg); break;
        case 's': S.span_e = atoi(optarg); break;
        case 'c': S.ckpt_file = strcmp(optarg, "none") ? optarg : NULL; break;
        case 'o': S.report_file = optarg; break;
        case 'r': S.resume = 1; break;
        case 'u': S.status_iv = atof(optarg); break;
        case 'q': S.quiet = 1; break;
        case 1: S.selftest = 1; break;
        case 2: S.chunk_min = strtoull(optarg, NULL, 10); break;
        case 3: S.chunk_max = strtoull(optarg, NULL, 10); break;
        case 4: bench_secs = optarg ? atof(optarg) : 45; break;
        case 'h': usage(); return 0;
        default: usage(); return 1;
        }
    }
    if (bench_secs != 0) {
        if (bench_secs < 10 || bench_secs > 600) {
            fprintf(stderr, "error: benchmark duration in [10,600] seconds\n");
            return 1;
        }
        if (S.threads < 1) S.threads = 1;
        build_tables();
        return benchmark(bench_secs);
    }
    if (S.selftest) { S.N = 8; S.ckpt_file = NULL; S.resume = 0; }
    if (S.N < 1 || S.N > 16) { fprintf(stderr, "error: N in [1,16]\n"); return 1; }
    if (S.span_e < 9 || S.span_e > 12) { fprintf(stderr, "error: span 10^E, E in [9,12]\n"); return 1; }
    if (S.threads < 1) S.threads = 1;

    build_tables();
    S.span = pow10u(S.span_e);
    S.target = pow10u(S.N);
    S.pmax = (u64)((double)P10N[S.N] * 1.0005) + S.span;
    S.nplan = (S.pmax + S.span - 1) / S.span;
    if (S.chunk_max == 0 || S.chunk_max > S.nplan) S.chunk_max = S.nplan;
    if (S.chunk_min >= S.chunk_max) {
        fprintf(stderr, "error: chunk range [%" PRIu64 ",%" PRIu64
                ") is empty (%" PRIu64 " chunks planned)\n",
                S.chunk_min, S.chunk_max, S.nplan);
        return 1;
    }
    S.res = calloc(S.nplan, sizeof(tally_t));
    S.done = calloc(S.nplan, 1);
    pthread_mutex_init(&S.lock, NULL);

    if (!S.quiet) {
        printf("pdigits: first 10^%d primes (to ~%.4e), %" PRIu64
               " chunks of 10^%d, %d threads\n",
               S.N, (double)P10N[S.N], S.nplan, S.span_e, S.threads);
        if (S.chunk_min > 0 || S.chunk_max < S.nplan)
            printf("  this machine: chunks [%" PRIu64 ", %" PRIu64 ")\n",
                   S.chunk_min, S.chunk_max);
    }

    if (S.resume) {
        if (!load_checkpoint()) return 1;
    }
    {
        u64 dir = 0;
        for (u64 i = S.chunk_min; i < S.chunk_max; i++)
            if (S.done[i]) dir++;
        atomic_store(&S.done_in_range, dir);
    }
    atomic_store(&S.next_chunk, S.chunk_min);
    if (S.ckpt_file) {
        S.log = fopen(S.ckpt_file, S.resume ? "a" : "w");
        if (!S.log) { fprintf(stderr, "error: cannot open %s\n", S.ckpt_file); return 1; }
        if (!S.resume) {
            fprintf(S.log, "pdigits-v2 N=%d span_e=%d\n", S.N, S.span_e);
            fflush(S.log);
        }
    }

    signal(SIGINT, on_sigint);
    signal(SIGTERM, on_sigint);
#ifdef SIGINFO
    signal(SIGINFO, on_siginfo);
#endif

    S.t_start = now_s();
    u64 done0 = atomic_load(&S.chunks_done);
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)S.threads);
    for (int i = 0; i < S.threads; i++)
        pthread_create(&th[i], NULL, worker, NULL);

    double last_status = now_s();
    int interrupted = 0;
    for (;;) {
        usleep(200000);
        if (g_sigint) { atomic_store(&S.stop, 1); interrupted = 1; break; }
        if (atomic_load(&S.stop)) break;
        if (atomic_load(&S.done_in_range) >= S.chunk_max - S.chunk_min) break;
        double t = now_s();
        if (!S.quiet && (t - last_status >= S.status_iv || g_siginfo)) {
            g_siginfo = 0; last_status = t;
            pthread_mutex_lock(&S.lock);
            u64 fr = S.frontier, np = S.cum_np_frontier;
            pthread_mutex_unlock(&S.lock);
            u64 cd = atomic_load(&S.chunks_done);
            double el = t - S.t_start;
            double rate = (double)(cd - done0) / (el + 1e-9);
            u64 remain = S.chunk_max - S.chunk_min - atomic_load(&S.done_in_range);
            double eta = rate > 0 ? (double)remain / rate : 0;
            printf("[%7.0fs] chunks %" PRIu64 "/%" PRIu64
                   " (frontier %" PRIu64 ")  %.3e primes  %.1f chunks/s  ETA %.1fh\n",
                   el, cd, S.nplan, fr, (double)np, rate, eta / 3600.0);
            fflush(stdout);
        }
    }
    for (int i = 0; i < S.threads; i++) pthread_join(th[i], NULL);
    if (S.log) fclose(S.log);

    if (interrupted) {
        pthread_mutex_lock(&S.lock);
        printf("\ninterrupted: %" PRIu64 " chunks done, frontier %" PRIu64
               " (%.3e primes); resume with -r\n",
               atomic_load(&S.chunks_done), S.frontier,
               (double)S.cum_np_frontier);
        pthread_mutex_unlock(&S.lock);
        return 130;
    }

    if (S.chunk_min > 0 || S.chunk_max < S.nplan) {
        pthread_mutex_lock(&S.lock);
        u64 np = S.cum_np_frontier;
        pthread_mutex_unlock(&S.lock);
        if (np < S.target) {
            printf("\nchunk range [%" PRIu64 ", %" PRIu64 ") complete.\n"
                   "To finish: collect every machine's log into one file, e.g.\n"
                   "    grep '^C ' other.ckpt >> %s\n"
                   "then run the same command on it with -r for the final report.\n",
                   S.chunk_min, S.chunk_max,
                   S.ckpt_file ? S.ckpt_file : "merged.ckpt");
            return 0;
        }
    }

    int ok = report();
    if (S.selftest)
        printf("selftest: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 2;
}
