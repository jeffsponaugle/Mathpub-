/*
 * gapchain_ps.c -- same search as gapchain.c, but built on primesieve.
 *
 * No bitmap and no candidate rescanning: primesieve hands us the primes in
 * order, and we slide a window of W = k+EXTRA consecutive primes across the
 * stream.  A prime heads a chain iff the gaps forward from it are 2, 4, 6,
 * ...; a prime that sits *inside* a chain fails at the first test (its next
 * gap is 4 or more), so each chain is reported exactly once and no
 * backward context is needed -- which is what makes the segments
 * independent and the whole thing trivially threadable.
 *
 * That same independence is what makes checkpointing cheap.  A checkpoint is
 * the number of *leading* segments that are completely finished plus the
 * chain-length histogram for exactly those segments, so it always describes a
 * consistent prefix of the range.  Segments that finished out of order past
 * the frontier are simply redone on resume; their chains are printed twice,
 * so treat a resumed output file as a multiset if you need exact counts.
 *
 * --first wants the *smallest* qualifying prime, which needs more care than
 * "stop on the first hit": threads run out of order, so a hit in segment 7
 * proves nothing while segments 3 and 5 are still being scanned.  Segments are
 * handed out in increasing order, so once any segment b produces a hit, every
 * unclaimed segment is above b and can be skipped, every in-flight segment
 * above b is irrelevant and is abandoned, and every in-flight segment below b
 * must be finished before b's hit can be declared the winner.
 *
 * Build (see README.md):
 *   cc -O3 -pthread gapchain_ps.c $(pkg-config --cflags --libs primesieve) \
 *      -o gapchain_ps
 * or, against a static primesieve:
 *   cc -O3 -pthread -I<primesieve>/include gapchain_ps.c \
 *      -L<primesieve> -lprimesieve -lstdc++ -o gapchain_ps
 */

#define _GNU_SOURCE
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <primesieve.h>

#define SEG_MIN (1ULL << 24)
#define SEG_MAX (1ULL << 32)
#define EXTRA 8
#define NBUCKET (EXTRA + 1)          /* reportable lengths are k .. k+EXTRA */
#define CKPT_MAGIC "gapchain_ps-checkpoint 1"
#define NO_SEG UINT64_MAX

static uint64_t g_lo, g_hi, g_scanoff, g_segsize, g_nseg, g_next_seg;
static int g_k, g_W;
static int g_progress, g_quiet, g_first;

/*
 * --gaps: INC is gaps 2, 4, 6, ... (A016045); DEC is the mirror image,
 * ..., 6, 4, 2 (A263049).  A DEC chain is found by its exact L-prime pattern
 * and extended backwards, since the pattern ends every chain it is part of.
 */
enum { INC, DEC };
static int g_on[2] = { 1, 0 };

/*
 * Completion bookkeeping.  g_frontier counts the leading segments that are
 * done; g_hist/g_found cover those segments only.  A segment that finishes
 * early parks its histogram in g_pend until the frontier reaches it.
 */
static pthread_mutex_t g_lock = PTHREAD_MUTEX_INITIALIZER;
static uint64_t g_frontier, g_found[2], g_hist[2][NBUCKET];
struct pend { uint64_t seg, hist[2][NBUCKET]; };
static struct pend *g_pend;
static size_t g_npend, g_cpend;

/* --first: best candidate so far.  g_best_seg is read lock-free by workers. */
static uint64_t g_best_seg = NO_SEG;
static uint64_t *g_best_chain;
static int g_best_len;

static pthread_mutex_t g_out = PTHREAD_MUTEX_INITIALIZER;
static volatile sig_atomic_t g_stop;
static int g_active;

static const char *g_ckfile;
static char *g_cktmp;
static unsigned g_ckint = 60;
static double g_t0;
static uint64_t g_base;              /* frontier we started this run from */

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* stderr message that does not leave half a progress line behind */
static void note(const char *fmt, ...)
{
    va_list ap;

    pthread_mutex_lock(&g_out);
    if (g_progress)
        fprintf(stderr, "\r%78s\r", "");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    pthread_mutex_unlock(&g_out);
}

static void fmt_dur(char *buf, size_t n, double s)
{
    if (!(s >= 0.0) || s > 3.0e8) {
        snprintf(buf, n, "--:--:--");
        return;
    }
    uint64_t t = (uint64_t)s, d = t / 86400;
    t %= 86400;
    if (d)
        snprintf(buf, n, "%" PRIu64 "d%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
                 d, t / 3600, (t / 60) % 60, t % 60);
    else
        snprintf(buf, n, "%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64,
                 t / 3600, (t / 60) % 60, t % 60);
}

static void fmt_rate(char *buf, size_t n, double v)
{
    static const char *unit[] = { "", "K", "M", "G", "T", "P" };
    int i = 0;
    while (v >= 1000.0 && i < 5) {
        v /= 1000.0;
        i++;
    }
    snprintf(buf, n, "%.2f%s/s", v, unit[i]);
}

/*
 * Accepts plain digits, e notation (1e14, 2.5e15), a K/M/B/G/T/P suffix, and
 * ,/_/space as digit separators.  Digit-only strings go through integer
 * arithmetic so nothing is lost above 2^53, where long double is only a
 * double on some platforms.
 */
static int parse_u64(const char *arg, uint64_t *out)
{
    char buf[128];
    size_t n = 0;
    for (const char *s = arg; *s; s++) {
        if (*s == ',' || *s == '_' || *s == ' ')
            continue;
        if (n + 1 >= sizeof buf)
            return -1;
        buf[n++] = *s;
    }
    buf[n] = '\0';
    if (!n)
        return -1;

    size_t digits = strspn(buf, "0123456789");
    if (digits && !buf[digits]) {
        uint64_t v = 0;
        for (size_t i = 0; i < digits; i++) {
            unsigned d = (unsigned)(buf[i] - '0');
            if (v > (UINT64_MAX - d) / 10)
                return -1;
            v = v * 10 + d;
        }
        if (v > 18440000000000000000ULL)    /* same cap as below */
            return -1;
        *out = v;
        return 0;
    }

    if (!digits && buf[0] != '.')    /* keep strtold off inf/nan/0x... */
        return -1;

    char *end;
    errno = 0;
    long double v = strtold(buf, &end);
    if (end == buf || !isfinite((double)v) || v < 0.0L)
        return -1;

    long double mul = 1.0L;
    if (*end) {
        switch (toupper((unsigned char)*end)) {
        case 'K': mul = 1e3L;  break;
        case 'M': mul = 1e6L;  break;
        case 'B':
        case 'G': mul = 1e9L;  break;
        case 'T': mul = 1e12L; break;
        case 'P': mul = 1e15L; break;
        default:  return -1;
        }
        if (*++end)
            return -1;
    }
    long double r = v * mul;
    if (r > 1.844e19L)
        return -1;
    *out = (uint64_t)(r + 0.5L);
    return 0;
}

static void seg_complete(uint64_t seg, uint64_t (*hist)[NBUCKET])
{
    pthread_mutex_lock(&g_lock);

    if (g_npend == g_cpend) {
        size_t c = g_cpend ? g_cpend * 2 : 32;
        struct pend *p = realloc(g_pend, c * sizeof *p);
        if (!p) {
            fprintf(stderr, "gapchain_ps: out of memory tracking segments\n");
            _exit(1);
        }
        g_pend = p;
        g_cpend = c;
    }
    g_pend[g_npend].seg = seg;
    memcpy(g_pend[g_npend].hist, hist, sizeof g_pend[g_npend].hist);
    g_npend++;

    for (int moved = 1; moved; ) {
        moved = 0;
        for (size_t i = 0; i < g_npend; i++) {
            if (g_pend[i].seg != g_frontier)
                continue;
            for (int d = 0; d < 2; d++)
                for (int b = 0; b < NBUCKET; b++) {
                    g_hist[d][b] += g_pend[i].hist[d][b];
                    g_found[d] += g_pend[i].hist[d][b];
                }
            g_pend[i] = g_pend[--g_npend];
            g_frontier++;
            moved = 1;
            break;
        }
    }

    pthread_mutex_unlock(&g_lock);
}

/*
 * --first: remember this chain if it beats the current best.  A worker scans
 * its segment in increasing order, so its first hit is that segment's minimum
 * and lower segment index always wins.
 */
static void record_candidate(uint64_t seg, const uint64_t *win, int head,
                             int len, int W)
{
    uint64_t frontier;
    int better;

    pthread_mutex_lock(&g_lock);
    better = seg < g_best_seg;
    frontier = g_frontier;
    if (better) {
        for (int j = 0, q = head; j < len; j++, q = q + 1 == W ? 0 : q + 1)
            g_best_chain[j] = win[q];
        g_best_len = len;
        __atomic_store_n(&g_best_seg, seg, __ATOMIC_RELEASE);
    }
    pthread_mutex_unlock(&g_lock);

    if (better)
        note("candidate: %" PRIu64 " (length %d, segment %" PRIu64 ")%s\n",
             win[head], len, seg,
             seg > frontier ? "; finishing earlier segments" : "");
}

static void *worker(void *arg)
{
    (void)arg;
    const int W = g_W;
    uint64_t *win = malloc((size_t)W * sizeof *win);
    uint64_t back[EXTRA];            /* DEC: the primes just before the head */

    while (win && !g_stop) {
        /*
         * Segments are handed out in ascending order, so once a candidate
         * exists every segment we could still claim lies above it.
         */
        if (g_first && __atomic_load_n(&g_best_seg, __ATOMIC_ACQUIRE) != NO_SEG)
            break;

        uint64_t s = __atomic_fetch_add(&g_next_seg, 1, __ATOMIC_RELAXED);
        if (s >= g_nseg)
            break;

        uint64_t seg_lo = g_lo + s * g_segsize;
        uint64_t seg_hi = seg_lo + g_segsize - 1;
        if (seg_hi > g_hi)
            seg_hi = g_hi;

        uint64_t hist[2][NBUCKET] = { { 0 } };
        unsigned tick = 0;
        int aborted = 0, bpos = 0;

        primesieve_iterator it;
        primesieve_init(&it);
        primesieve_jump_to(&it, seg_lo, seg_hi + g_scanoff);

        for (int i = 0; i < W; i++)
            win[i] = primesieve_next_prime(&it);

        /*
         * back[(bpos - j) mod EXTRA] is the j-th prime before the head, 0 if
         * there is none; seeded here with the primes below the segment.
         */
        if (g_on[DEC]) {
            primesieve_iterator bi;
            primesieve_init(&bi);
            primesieve_jump_to(&bi, seg_lo - 1, 0);
            for (int j = 0; j < EXTRA; j++)
                back[(EXTRA - j) % EXTRA] = primesieve_prev_prime(&bi);
            primesieve_free_iterator(&bi);
        }

        int head = 0;
        while (win[head] <= seg_hi) {
            uint64_t p = win[head];

            /* someone found a hit below us: nothing here can win any more */
            if (g_first && !(++tick & 0x3FF) &&
                __atomic_load_n(&g_best_seg, __ATOMIC_ACQUIRE) < s) {
                aborted = 1;
                break;
            }

            int len = 1, i = head;
            for (int j = 1; g_on[INC] && j < W; j++) {
                int nx = i + 1 == W ? 0 : i + 1;
                if (win[nx] - win[i] != (uint64_t)(2 * j))
                    break;
                len++;
                i = nx;
            }

            if (g_on[INC] && len >= g_k) {
                if (g_first) {
                    record_candidate(s, win, head, len, W);
                    aborted = 1;      /* rest of this segment is above p */
                    break;
                }
                hist[INC][len - g_k]++;
                pthread_mutex_lock(&g_out);
                printf("%s%2d: %" PRIu64, g_on[DEC] ? "inc " : "", len, p);
                int q = head;
                for (int j = 1; j < len; j++) {
                    q = q + 1 == W ? 0 : q + 1;
                    printf(" %" PRIu64, win[q]);
                }
                if (len == W)
                    fputs(" +", stdout);
                putchar('\n');
                pthread_mutex_unlock(&g_out);
            }

            if (g_on[DEC]) {
                /* the exact pattern from p: gaps 2(L-1), ..., 4, 2 */
                int ok = 1;
                i = head;
                for (int j = 1; ok && j < g_k; j++) {
                    int nx = i + 1 == W ? 0 : i + 1;
                    ok = win[nx] - win[i] == 2 * (uint64_t)(g_k - j);
                    i = nx;
                }
                if (ok) {
                    uint64_t h = p;
                    int dl = g_k;
                    for (int j = 0; j < EXTRA; j++) {
                        uint64_t pb = back[(bpos - j + EXTRA) % EXTRA];
                        if (!pb || h - pb != 2 * (uint64_t)dl)
                            break;
                        h = pb;
                        dl++;
                    }
                    hist[DEC][dl - g_k]++;
                    pthread_mutex_lock(&g_out);
                    printf("dec %2d: %" PRIu64, dl, h);
                    for (int m = 1; m < dl; m++)
                        printf(" %" PRIu64,
                               h + (uint64_t)m * (uint64_t)(2 * dl - 1 - m));
                    if (dl == W)
                        fputs(" +", stdout);
                    putchar('\n');
                    pthread_mutex_unlock(&g_out);
                }
                bpos = (bpos + 1) % EXTRA;
                back[bpos] = p;
            }

            win[head] = primesieve_next_prime(&it);
            head = head + 1 == W ? 0 : head + 1;
        }
        primesieve_free_iterator(&it);

        /* only a fully scanned segment may advance the checkpoint frontier */
        if (!aborted)
            seg_complete(s, hist);
    }

    free(win);
    __atomic_sub_fetch(&g_active, 1, __ATOMIC_RELEASE);
    return NULL;
}

static void show_progress(int newline)
{
    uint64_t f, found;
    pthread_mutex_lock(&g_lock);
    f = g_frontier;
    found = g_found[INC] + g_found[DEC];
    pthread_mutex_unlock(&g_lock);

    double el = now_s() - g_t0;
    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    double rate = el > 0.0 ? (double)(f - g_base) * (double)g_segsize / el
                           : 0.0;
    double eta = rate > 0.0 ? (double)(g_nseg - f) * (double)g_segsize / rate
                            : -1.0;

    char sr[32], se[32], sl[32];
    fmt_rate(sr, sizeof sr, rate);
    fmt_dur(se, sizeof se, eta);
    fmt_dur(sl, sizeof sl, el);

    pthread_mutex_lock(&g_out);
    fprintf(stderr,
            "%s%6.2f%%  at %.6g  %s  eta %s  elapsed %s  %" PRIu64 " chains%s",
            newline ? "" : "\r", pct,
            (double)(g_lo + (f < g_nseg ? f * g_segsize : g_hi - g_lo)),
            sr, se, sl, found, newline ? "\n" : "        ");
    fflush(stderr);
    pthread_mutex_unlock(&g_out);
}

static const char *mode_name(void)
{
    return g_on[INC] && g_on[DEC] ? "both" : g_on[DEC] ? "dec" : "inc";
}

static void write_checkpoint(void)
{
    uint64_t frontier, found[2], hist[2][NBUCKET];

    pthread_mutex_lock(&g_lock);
    frontier = g_frontier;
    memcpy(found, g_found, sizeof found);
    memcpy(hist, g_hist, sizeof hist);
    pthread_mutex_unlock(&g_lock);

    /*
     * Sample the frontier first, then flush: everything the checkpoint claims
     * is finished is guaranteed on disk, and the flush can only carry extra
     * output past the frontier, which resume merely repeats.
     */
    pthread_mutex_lock(&g_out);
    fflush(stdout);
    pthread_mutex_unlock(&g_out);

    FILE *f = fopen(g_cktmp, "w");
    if (!f) {
        note("gapchain_ps: %s: %s\n", g_cktmp, strerror(errno));
        return;
    }
    fprintf(f, "%s\n", CKPT_MAGIC);
    fprintf(f, "k %d\n", g_k);
    fprintf(f, "lo %" PRIu64 "\n", g_lo);
    fprintf(f, "hi %" PRIu64 "\n", g_hi);
    fprintf(f, "segsize %" PRIu64 "\n", g_segsize);
    fprintf(f, "nseg %" PRIu64 "\n", g_nseg);
    fprintf(f, "first %d\n", g_first ? 1 : 0);
    fprintf(f, "gaps %s\n", mode_name());
    fprintf(f, "done %" PRIu64 "\n", frontier);
    /* INC keeps the original line names; DEC's carry a "dec" prefix */
    for (int d = 0; d < 2; d++) {
        const char *pre = d == DEC ? "dec" : "";
        if (!g_on[d])
            continue;
        fprintf(f, "%sfound %" PRIu64 "\n", pre, found[d]);
        for (int b = 0; b < NBUCKET; b++)
            fprintf(f, "%slen %d %" PRIu64 "\n", pre, g_k + b, hist[d][b]);
    }

    if (fflush(f) || fsync(fileno(f)) || fclose(f)) {
        note("gapchain_ps: %s: %s\n", g_cktmp, strerror(errno));
        return;
    }
    if (rename(g_cktmp, g_ckfile))
        note("gapchain_ps: %s: %s\n", g_ckfile, strerror(errno));
}

/* 0 = state loaded, -1 = no such file.  A mismatched checkpoint is fatal. */
static int read_checkpoint(void)
{
    FILE *f = fopen(g_ckfile, "r");
    if (!f) {
        if (errno == ENOENT)
            return -1;
        fprintf(stderr, "gapchain_ps: %s: %s\n", g_ckfile, strerror(errno));
        exit(1);
    }

    char line[256];
    uint64_t lo = 0, hi = 0, segsize = 0, nseg = 0, done = NO_SEG;
    uint64_t found[2] = { 0 }, hist[2][NBUCKET] = { { 0 } };
    int k = 0, magic = 0, first = 0;   /* pre-0.2 files have no "first" line */
    char gaps[16] = "inc";             /* ... nor "gaps" */

    while (fgets(line, sizeof line, f)) {
        uint64_t v;
        int len;
        if (!strncmp(line, CKPT_MAGIC, strlen(CKPT_MAGIC)))
            magic = 1;
        else if (sscanf(line, "k %d", &k) == 1)
            ;
        else if (sscanf(line, "lo %" SCNu64, &lo) == 1)
            ;
        else if (sscanf(line, "hi %" SCNu64, &hi) == 1)
            ;
        else if (sscanf(line, "segsize %" SCNu64, &segsize) == 1)
            ;
        else if (sscanf(line, "nseg %" SCNu64, &nseg) == 1)
            ;
        else if (sscanf(line, "first %d", &first) == 1)
            ;
        else if (sscanf(line, "done %" SCNu64, &done) == 1)
            ;
        else if (sscanf(line, "gaps %15s", gaps) == 1)
            ;
        else if (sscanf(line, "found %" SCNu64, &found[INC]) == 1)
            ;
        else if (sscanf(line, "decfound %" SCNu64, &found[DEC]) == 1)
            ;
        else if (sscanf(line, "len %d %" SCNu64, &len, &v) == 2) {
            if (len >= k && len - k < NBUCKET)
                hist[INC][len - k] = v;
        }
        else if (sscanf(line, "declen %d %" SCNu64, &len, &v) == 2) {
            if (len >= k && len - k < NBUCKET)
                hist[DEC][len - k] = v;
        }
    }
    fclose(f);

    if (!magic || done == NO_SEG) {
        fprintf(stderr, "gapchain_ps: %s: not a checkpoint file\n", g_ckfile);
        exit(1);
    }
    if (k != g_k || lo != g_lo || hi != g_hi ||
        first != (g_first ? 1 : 0) || strcmp(gaps, mode_name())) {
        fprintf(stderr,
                "gapchain_ps: %s describes a different search\n"
                "  checkpoint: -l %d -s %" PRIu64 " -e %" PRIu64 " --gaps %s%s\n"
                "  requested:  -l %d -s %" PRIu64 " -e %" PRIu64 " --gaps %s%s\n",
                g_ckfile, k, lo, hi, gaps, first ? " --first" : "",
                g_k, g_lo, g_hi, mode_name(), g_first ? " --first" : "");
        exit(1);
    }

    /*
     * The file's segment size defines what "done" counts, so adopt it even if
     * this binary would pick a different one for a fresh run.
     */
    if (!segsize || nseg != (g_hi - g_lo) / segsize + 1 || done > nseg) {
        fprintf(stderr, "gapchain_ps: %s: inconsistent segment data\n",
                g_ckfile);
        exit(1);
    }
    g_segsize = segsize;
    g_nseg = nseg;

    g_frontier = g_next_seg = done;
    memcpy(g_found, found, sizeof g_found);
    memcpy(g_hist, hist, sizeof g_hist);
    return 0;
}

static void print_first_result(void)
{
    uint64_t f = g_frontier;
    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    char el[32];
    fmt_dur(el, sizeof el, now_s() - g_t0);

    if (g_best_seg != NO_SEG) {
        printf("%2d: %" PRIu64, g_best_len, g_best_chain[0]);
        for (int j = 1; j < g_best_len; j++)
            printf(" %" PRIu64, g_best_chain[j]);
        if (g_best_len == g_W)
            fputs(" +", stdout);
        putchar('\n');
        fflush(stdout);

        fprintf(stderr, "\nfirst chain of length >= %d is at %" PRIu64
                ", after searching %.2f%% of the range in %s\n",
                g_k, g_best_chain[0], pct, el);
        return;
    }

    if (g_stop)
        fprintf(stderr, "\nstopped at %.2f%% of %" PRIu64 "..%" PRIu64
                " with no chain of length >= %d yet, %s elapsed\n",
                pct, g_lo, g_hi, g_k, el);
    else
        fprintf(stderr, "\nno chain of length >= %d in %" PRIu64 "..%" PRIu64
                ", %s elapsed\n", g_k, g_lo, g_hi, el);

    if (g_stop && g_ckfile)
        fprintf(stderr, "resume with: -c %s -r\n", g_ckfile);
}

static void print_stats(void)
{
    uint64_t f, found[2], hist[2][NBUCKET];

    if (g_first) {
        print_first_result();
        return;
    }

    pthread_mutex_lock(&g_lock);
    f = g_frontier;
    memcpy(found, g_found, sizeof found);
    memcpy(hist, g_hist, sizeof hist);
    pthread_mutex_unlock(&g_lock);

    double pct = g_nseg ? 100.0 * (double)f / (double)g_nseg : 100.0;
    char el[32];
    fmt_dur(el, sizeof el, now_s() - g_t0);

    fprintf(stderr, "\n%s %.2f%% of %" PRIu64 "..%" PRIu64
            " (%" PRIu64 "/%" PRIu64 " segments), %s elapsed\n",
            f >= g_nseg ? "completed" : "stopped at", pct, g_lo, g_hi,
            f, g_nseg, el);

    for (int d = 0; d < 2; d++) {
        if (!g_on[d])
            continue;
        if (g_on[DEC])
            fprintf(stderr, "%s gaps (%s):\n",
                    d == INC ? "increasing" : "decreasing",
                    d == INC ? "A016045" : "A263049");
        if (!found[d]) {
            fprintf(stderr, "no chains found\n");
            continue;
        }
        fprintf(stderr, "chains by length:\n");
        for (int b = 0; b < NBUCKET; b++) {
            if (!hist[d][b])
                continue;
            char tag[16];
            snprintf(tag, sizeof tag, "%d%s", g_k + b,
                     b == NBUCKET - 1 ? "+" : "");
            fprintf(stderr, "  %6s  %14" PRIu64 "  %6.2f%%%s\n", tag,
                    hist[d][b], 100.0 * (double)hist[d][b] / (double)found[d],
                    b == NBUCKET - 1 ? "   (filled the window)" : "");
        }
        fprintf(stderr, "  %6s  %14" PRIu64 "\n", "total", found[d]);
    }

    if (f < g_nseg && g_ckfile)
        fprintf(stderr, "resume with: -c %s -r\n", g_ckfile);
    else if (f < g_nseg)
        fprintf(stderr, "(no checkpoint file; rerun with -c FILE to make "
                        "interrupted runs resumable)\n");
}

static void on_signal(int sig)
{
    (void)sig;
    if (g_stop)                       /* second one: give up immediately */
        _exit(130);
    g_stop = 1;
}

static void monitor(void)
{
    double last_ck = now_s(), last_log = now_s();

    for (;;) {
        struct timespec ts = { 0, 250000000L };
        nanosleep(&ts, NULL);         /* a signal cuts this short */

        int live = __atomic_load_n(&g_active, __ATOMIC_ACQUIRE);
        double t = now_s();

        if (g_progress)
            show_progress(0);
        else if (!g_quiet && t - last_log >= 60.0) {
            show_progress(1);
            last_log = t;
        }
        if (g_ckfile && g_ckint && t - last_ck >= (double)g_ckint) {
            write_checkpoint();
            last_ck = t;
        }
        if (!live || g_stop)
            break;
    }
}

/* never returns, so option cases that call it need no break */
static __attribute__((noreturn)) void usage(FILE *f, const char *prog,
                                            int status)
{
    fprintf(f,
        "usage: %s -l LEN -e END [options]\n"
        "\n"
        "Finds runs of consecutive primes whose gaps grow by two:\n"
        "2, 4, 6, ... (OEIS A016045) or, mirrored, ..., 6, 4, 2 (A263049).\n"
        "See README.md.\n"
        "\n"
        "required:\n"
        "  -l, --length N      chain length: least number of primes in a run\n"
        "  -e, --end N         highest prime to consider as a chain head\n"
        "\n"
        "options:\n"
        "  -s, --start N       lowest prime to consider (default 3)\n"
        "  -t, --threads N     worker threads (default: online CPUs)\n"
        "  -g, --gaps WHICH    inc (2, 4, 6, ...; default), dec (..., 6, 4, 2)\n"
        "                      or both\n"
        "  -f, --first         report only the smallest qualifying chain,\n"
        "                      then stop (inc only; gapsieve does all three)\n"
        "  -c, --checkpoint F  checkpoint to F periodically and on SIGINT\n"
        "  -i, --interval S    checkpoint seconds, 0 = on exit (default 60)\n"
        "  -r, --resume        resume from F (required if F exists)\n"
        "  -q, --quiet         suppress progress output\n"
        "  -h, --help          this message\n"
        "\n"
        "N accepts 1e14, 2.5e15, 300T, 1_000_000 or plain digits.\n",
        prog);
    exit(status);
}

int main(int argc, char **argv)
{
    static const struct option longopts[] = {
        { "length",     required_argument, NULL, 'l' },
        { "start",      required_argument, NULL, 's' },
        { "end",        required_argument, NULL, 'e' },
        { "threads",    required_argument, NULL, 't' },
        { "first",      no_argument,       NULL, 'f' },
        { "gaps",       required_argument, NULL, 'g' },
        { "checkpoint", required_argument, NULL, 'c' },
        { "interval",   required_argument, NULL, 'i' },
        { "resume",     no_argument,       NULL, 'r' },
        { "quiet",      no_argument,       NULL, 'q' },
        { "help",       no_argument,       NULL, 'h' },
        { NULL,         0,                 NULL,  0  },
    };
    const char *prog = argv[0];
    const char *arg_end = NULL;
    long k = 0, nthreads = 0;
    int resume = 0, opt;

    g_lo = 3;
    while ((opt = getopt_long(argc, argv, "l:s:e:t:fg:c:i:rqh",
                              longopts, NULL)) != -1) {
        switch (opt) {
        case 'l':
            k = strtol(optarg, NULL, 10);
            if (k < 2 || k > 1000) {
                fprintf(stderr, "gapchain_ps: --length must be 2..1000\n");
                return 2;
            }
            break;
        case 's':
            if (parse_u64(optarg, &g_lo)) {
                fprintf(stderr, "gapchain_ps: bad --start value '%s'\n", optarg);
                return 2;
            }
            break;
        case 'e':
            arg_end = optarg;
            if (parse_u64(optarg, &g_hi)) {
                fprintf(stderr, "gapchain_ps: bad --end value '%s'\n", optarg);
                return 2;
            }
            break;
        case 't':
            nthreads = strtol(optarg, NULL, 10);
            if (nthreads < 1) {
                fprintf(stderr, "gapchain_ps: --threads must be >= 1\n");
                return 2;
            }
            break;
        case 'f': g_first = 1; break;
        case 'g':
            g_on[INC] = !strcmp(optarg, "inc") || !strcmp(optarg, "both");
            g_on[DEC] = !strcmp(optarg, "dec") || !strcmp(optarg, "both");
            if (!g_on[INC] && !g_on[DEC]) {
                fprintf(stderr, "gapchain_ps: --gaps must be inc, dec or both\n");
                return 2;
            }
            break;
        case 'c': g_ckfile = optarg; break;
        case 'i': g_ckint = (unsigned)strtoul(optarg, NULL, 10); break;
        case 'r': resume = 1; break;
        case 'q': g_quiet = 1; break;
        case 'h': usage(stdout, prog, 0);
        default:  usage(stderr, prog, 2);
        }
    }

    if (optind < argc) {
        fprintf(stderr,
                "gapchain_ps: unexpected argument '%s'\n"
                "  arguments are no longer positional; use -l LEN -s START "
                "-e END\n", argv[optind]);
        return 2;
    }
    if (!k || !arg_end) {
        fprintf(stderr, "gapchain_ps: %s is required\n",
                !k ? "--length" : "--end");
        usage(stderr, prog, 2);
    }
    if (resume && !g_ckfile) {
        fprintf(stderr, "gapchain_ps: --resume needs --checkpoint FILE\n");
        return 2;
    }
    if (g_first && g_on[DEC]) {
        fprintf(stderr, "gapchain_ps: --first supports --gaps inc only; "
                        "gapsieve --first handles dec and both\n");
        return 2;
    }

    g_k = (int)k;
    g_W = g_k + EXTRA;
    g_scanoff = (uint64_t)g_W * (uint64_t)(g_W - 1);

    if (g_lo < 3)
        g_lo = 3;
    if (g_hi < g_lo) {
        fprintf(stderr, "gapchain_ps: --end is below --start\n");
        return 2;
    }
    if (!nthreads)
        nthreads = sysconf(_SC_NPROCESSORS_ONLN);
    if (nthreads < 1)
        nthreads = 1;

    /*
     * Each primesieve_jump_to costs O(sqrt(hi)) to rebuild sieving primes, so
     * high searches need segments large enough to amortize it: at 2.2e17,
     * 2^24 segments spend 80% of their time initializing.  Scale the segment
     * with sqrt(hi), but keep at least ~4 segments per thread for load
     * balance and checkpoint granularity.
     */
    g_segsize = SEG_MIN;
    while (g_segsize < SEG_MAX &&
           (long double)g_segsize < 8.0L * sqrtl((long double)g_hi))
        g_segsize <<= 1;
    while (g_segsize > SEG_MIN &&
           (g_hi - g_lo) / g_segsize + 1 < (uint64_t)(4 * nthreads))
        g_segsize >>= 1;

    g_nseg = (g_hi - g_lo) / g_segsize + 1;
    g_progress = isatty(STDERR_FILENO) && !g_quiet;

    g_best_chain = malloc((size_t)g_W * sizeof *g_best_chain);
    if (!g_best_chain) {
        fprintf(stderr, "gapchain_ps: out of memory\n");
        return 1;
    }

    if (g_ckfile) {
        size_t n = strlen(g_ckfile) + 5;
        g_cktmp = malloc(n);
        if (!g_cktmp) {
            fprintf(stderr, "gapchain_ps: out of memory\n");
            return 1;
        }
        snprintf(g_cktmp, n, "%s.tmp", g_ckfile);

        if (!resume && !access(g_ckfile, F_OK)) {
            fprintf(stderr, "gapchain_ps: %s exists; pass --resume to resume "
                            "from it or remove it first\n", g_ckfile);
            return 2;
        }
        if (resume && read_checkpoint() == 0) {
            if (g_frontier >= g_nseg) {
                fprintf(stderr, "gapchain_ps: %s is already complete\n",
                        g_ckfile);
                g_t0 = now_s();
                print_stats();
                return 0;
            }
            fprintf(stderr, "resuming from %s at %" PRIu64 "/%" PRIu64
                    " segments (%.2f%%)", g_ckfile, g_frontier, g_nseg,
                    100.0 * (double)g_frontier / (double)g_nseg);
            if (g_first)
                fputc('\n', stderr);
            else
                fprintf(stderr, ", %" PRIu64 " chains so far\n",
                        g_found[INC] + g_found[DEC]);
        }
    }

    fprintf(stderr,
        "consecutive primes%s, chain length >= %d, range %" PRIu64 "..%" PRIu64
        ", %ld thread%s%s (primesieve)\n",
        !g_on[DEC] ? "" : g_on[INC] ? ", increasing and decreasing gaps"
                                    : ", decreasing gaps",
        g_k, g_lo, g_hi, nthreads, nthreads == 1 ? "" : "s",
        g_first ? ", first match only" : "");

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    /* Keep the handler off the workers so their stdio is never interrupted. */
    sigset_t block, prev;
    sigemptyset(&block);
    sigaddset(&block, SIGINT);
    sigaddset(&block, SIGTERM);
    pthread_sigmask(SIG_BLOCK, &block, &prev);

    g_t0 = now_s();
    g_base = g_frontier;

    pthread_t *t = calloc((size_t)nthreads, sizeof *t);
    if (!t) {
        fprintf(stderr, "gapchain_ps: out of memory\n");
        return 1;
    }
    long started = 0;
    __atomic_store_n(&g_active, (int)nthreads, __ATOMIC_RELEASE);
    for (long i = 0; i < nthreads; i++) {
        if (pthread_create(&t[i], NULL, worker, NULL)) {
            __atomic_sub_fetch(&g_active, (int)(nthreads - i), __ATOMIC_RELEASE);
            break;
        }
        started++;
    }
    if (!started) {
        fprintf(stderr, "gapchain_ps: could not start any worker threads\n");
        return 1;
    }

    pthread_sigmask(SIG_SETMASK, &prev, NULL);
    monitor();

    if (g_stop && g_progress)
        note("interrupt: finishing segments in flight...\n");
    for (long i = 0; i < started; i++)
        pthread_join(t[i], NULL);
    free(t);

    if (g_ckfile)
        write_checkpoint();

    if (g_progress)
        fprintf(stderr, "\r%78s\r", "");
    fflush(stdout);
    print_stats();
    free(g_pend);
    free(g_cktmp);
    free(g_best_chain);
    return g_stop ? 130 : 0;
}
