/*
 * prime_runs_ps.c -- Runs of CONSECUTIVE primes sharing a digit sum,
 *                    powered by the primesieve library.
 *
 * The prime stream comes from primesieve's iterator (state-of-the-art,
 * memory-light, generates primes in order very fast).  On top of it sits the
 * same machinery as before:
 *   - mod-(base-1) lazy reject: digit sum computed only when a run could extend
 *   - table-based digit sum for base 10
 *   - pthreads with overlap-and-ownership so output is identical for any -t
 *
 * Build:
 *   cc -O3 -march=native -pthread -o prime_runs_ps prime_runs_ps.c -lprimesieve -lm
 *
 * Usage:
 *   ./prime_runs_ps -n 8 -L 3000000000 -t 8
 *   ./prime_runs_ps -n 3 -L 60000 --all
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <primesieve.h>

typedef uint64_t u64;
typedef unsigned long long ull;

static inline double mono(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
}

/* Parse a count with optional suffix: T=1e12, B=1e9, M=1e6, K=1e3.
 * With a suffix the numeric part may be fractional ("1.5B" -> 1500000000).
 * Without a suffix it is read as an exact integer. */
static u64 parse_num(const char *s) {
    size_t n = strlen(s);
    if (!n) return 0;
    double mult = 0;
    switch (s[n - 1]) {
        case 'T': case 't': mult = 1e12; break;
        case 'B': case 'b': mult = 1e9;  break;
        case 'M': case 'm': mult = 1e6;  break;
        case 'K': case 'k': mult = 1e3;  break;
    }
    if (mult > 0) return (u64)(strtod(s, NULL) * mult + 0.5);
    return strtoull(s, NULL, 10);
}

static int  g_L, g_base, g_find_all, g_quiet, g_emirps;
static u64  g_M, g_start, g_end;

#define OVERLAP (1u << 20)        /* boundary margin (integers) */

/* ---- digit sum: base-10 table, generic fallback ---------------------- */
static int ds5[100000];
static inline int dsum(u64 n) {
    if (g_base == 10) { int s = 0; while (n) { s += ds5[n % 100000]; n /= 100000; } return s; }
    int s = 0; while (n) { s += (int)(n % (u64)g_base); n /= (u64)g_base; } return s;
}

/* ---- deterministic Miller-Rabin for all 64-bit n ---------------------- */
static inline u64 mulmod(u64 a, u64 b, u64 m) {
    return (u64)((__uint128_t)a * b % m);
}
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1; a %= m;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_u64(u64 n) {
    if (n < 2) return 0;
    static const int sm[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n == (u64)sm[i]) return 1;
        if (n % (u64)sm[i] == 0) return 0;
    }
    u64 d = n - 1; int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }
    for (int i = 0; i < 12; i++) {                 /* these bases are proven  */
        u64 x = powmod((u64)sm[i], d, n);          /* sufficient for n < 2^64 */
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int j = 1; j < r; j++) {
            x = mulmod(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

/* Reverse the digits of n in base g_base; returns 0 on overflow. */
static u64 revdigits(u64 n) {
    u64 b = (u64)g_base, r = 0;
    while (n) {
        u64 dg = n % b; n /= b;
        if (r > (UINT64_MAX - dg) / b) return 0;   /* would overflow */
        r = r * b + dg;
    }
    return r;
}

/* Leading digit of n in base 10, O(1) via a power table. */
static const u64 PW10[20] = {
    1ULL,10ULL,100ULL,1000ULL,10000ULL,100000ULL,1000000ULL,10000000ULL,
    100000000ULL,1000000000ULL,10000000000ULL,100000000000ULL,
    1000000000000ULL,10000000000000ULL,100000000000000ULL,
    1000000000000000ULL,10000000000000000ULL,100000000000000000ULL,
    1000000000000000000ULL,10000000000000000000ULL };
static inline int lead_digit10(u64 n) {
    int d = (int)(((63 - __builtin_clzll(n)) * 1233) >> 12) + 1;   /* ~log10 */
    if (d < 20 && n >= PW10[d]) d++;                               /* fix under-estimate */
    return (int)(n / PW10[d - 1]);
}
/* A prime whose LEADING digit is 2,4,5,6 or 8 can never be an emirp: the
 * reversal then ENDS in that digit, so it is even or a multiple of 5.
 * (Single-digit 2 and 5 reverse to themselves and are palindromes anyway.)
 * Index = leading digit; 1 means "cannot be an emirp". */
static const char BAD_LEAD10[10] = {1,0,1,0,1,1,1,0,1,0};

/* An emirp: prime whose digit-reversal is a DIFFERENT prime.
 * (Palindromic primes are excluded, per the usual convention.)
 * n is already known to be prime when this is called. */
static int is_emirp(u64 n) {
    if (g_base == 10 && n >= 10 && BAD_LEAD10[lead_digit10(n)])
        return 0;                                  /* cheap reject, no reversal */
    u64 r = revdigits(n);
    if (r == 0 || r == n) return 0;                /* overflow or palindrome */
    return is_prime_u64(r);
}

/* ---- per-thread tracker + results ------------------------------------ */
typedef struct {
    u64 lo, hi, A, B; int tid;
    u64 *first; int *dsv; u64 *members; size_t count, cap;
    int len; u64 res; int ds; u64 *ring; size_t rpos; u64 head; int done;
    signed char *eflag;  /* per ring slot: -1 unknown, 0 not emirp, 1 emirp */
    u64 covered;         /* good-measure scanned so far (for progress) */
} Targ;

static void record(Targ *t) {
    if (t->count == t->cap) {
        t->cap = t->cap ? t->cap * 2 : 16;
        t->first   = realloc(t->first,   t->cap * sizeof(u64));
        t->dsv     = realloc(t->dsv,     t->cap * sizeof(int));
        t->members = realloc(t->members, t->cap * (size_t)g_L * sizeof(u64));
    }
    /* the current window is the last g_L primes pushed into the ring */
    size_t base = t->rpos - (size_t)g_L;
    u64 *dst = &t->members[t->count * (size_t)g_L];
    for (int j = 0; j < g_L; j++)
        dst[j] = t->ring[(base + (size_t)j) % (size_t)g_L];
    t->first[t->count] = dst[0];
    t->dsv[t->count]   = t->ds;
    t->count++;
}

static inline void feed(Targ *t, u64 q) {
    if (t->done) return;
    /* For base 10 the divisor is the literal 9, which the compiler turns into
     * a multiply-shift; the runtime g_M path (other bases) keeps a real mod. */
    u64 rq = (g_base == 10) ? (q % 9u) : (q % g_M);
    int extended = 0;
    if (t->len >= 1 && rq == t->res) {                 /* cheap reject first */
        if (t->ds < 0) t->ds = dsum(t->head);
        if (dsum(q) == t->ds) {
            t->eflag[t->rpos % (size_t)g_L] = -1;
            t->ring[t->rpos % (size_t)g_L] = q; t->rpos++;
            t->len++; extended = 1;
        }
    }
    if (!extended) {                                   /* start a new run */
        t->len = 1; t->res = rq; t->ds = -1; t->head = q;
        t->rpos = 0;
        t->ring[0] = q; t->eflag[0] = -1; t->rpos = 1;
    }

    /* Every window of g_L consecutive primes inside a maximal run counts, so
     * emit once for each prime that brings the run to length >= g_L.  A
     * maximal run of length M therefore yields M - g_L + 1 windows. */
    if (t->len >= g_L) {
        if (g_emirps) {
            /* Only now -- once L primes actually share a digit sum -- is the
             * expensive reversal-primality test worth doing.  Flags are cached
             * per ring slot, so overlapping windows never retest a prime, and
             * we scan newest-first since that slot is the only new one. */
            size_t wbase = t->rpos - (size_t)g_L;
            for (int j = g_L - 1; j >= 0; j--) {
                size_t sl = (wbase + (size_t)j) % (size_t)g_L;
                if (t->eflag[sl] < 0)
                    t->eflag[sl] = is_emirp(t->ring[sl]) ? 1 : 0;
                if (!t->eflag[sl]) return;             /* not all emirps */
            }
        }
        u64 f = t->ring[(t->rpos - (size_t)g_L) % (size_t)g_L];
        if (f >= t->lo && f < t->hi) {
            if (t->ds < 0) t->ds = dsum(t->head);
            record(t);
            if (!g_find_all) t->done = 1;
        }
    }
}

/* ---- progress -------------------------------------------------------- */
static pthread_mutex_t g_pmx = PTHREAD_MUTEX_INITIALIZER;
static u64 g_primes_done;
static double g_t_start, g_t_last, g_rate;   /* timing + windowed rate */
static u64 g_count_last;
static int g_nthreads;                       /* for summing progress */
static u64 *g_ivlo, *g_ivhi; static int g_niv;   /* scan intervals (merged) */
static u64 *g_pos;                           /* per-thread value-progress (p - A) */
static u64 g_total_span;                      /* sum of (B - A) over threads */

static void fmt_hms(double secs, char *out, size_t n) {
    if (secs < 0) secs = 0;
    if (secs > 3.15e9) { snprintf(out, n, ">100y"); return; }
    long s = (long)(secs + 0.5), h = s/3600; s %= 3600; long m = s/60; s %= 60;
    if (h > 0) snprintf(out, n, "%ld:%02ld:%02ld", h, m, s);
    else       snprintf(out, n, "%ld:%02ld", m, s);
}

/* Reset run state.  Required at every interval boundary: two primes on
 * opposite sides of a skipped region are NOT consecutive in the surviving
 * stream, so carrying state across a gap would fabricate a false window. */
static void reset_tracker(Targ *t) {
    t->len = 0; t->res = 0; t->ds = -1; t->rpos = 0; t->head = 0;
    memset(t->eflag, -1, (size_t)g_L);
}

/* ---- worker: stream primes over this thread's pieces of the intervals -- */
static void *worker(void *arg) {
    Targ *t = (Targ *)arg;
    t->ring = (u64 *)malloc((size_t)g_L * sizeof(u64));
    t->eflag = (signed char *)malloc((size_t)g_L);
    t->covered = 0;
    reset_tracker(t);

    u64 local = 0;
    for (int j = 0; j < g_niv && !t->done; j++) {
        u64 a = g_ivlo[j] > t->A ? g_ivlo[j] : t->A;
        u64 b = g_ivhi[j] < t->B ? g_ivhi[j] : t->B;
        if (a >= b) continue;

        reset_tracker(t);                 /* gap => primes are not consecutive */

        primesieve_iterator it;
        primesieve_init(&it);
        primesieve_jump_to(&it, a, b);
        u64 p;
        while ((p = primesieve_next_prime(&it)) < b) {
            feed(t, p);
            if (t->done) break;
            if (!g_quiet && ((++local & ((1u << 22) - 1)) == 0)) {
                pthread_mutex_lock(&g_pmx);
                g_primes_done += (1u << 22);
                g_pos[t->tid] = t->covered + (p - a);
                double now = mono(), since = now - g_t_last;
                if (since >= 0.2) {
                    g_rate = (double)(g_primes_done - g_count_last) / since;
                    g_t_last = now; g_count_last = g_primes_done;
                }
                double elapsed = now - g_t_start;
                double avg = elapsed > 0 ? (double)g_primes_done / elapsed : 0;
                double cur = g_rate > 0 ? g_rate : avg;
                u64 covered = 0; for (int i = 0; i < g_nthreads; i++) covered += g_pos[i];
                double frac = g_total_span ? (double)covered / (double)g_total_span : 0;
                if (frac > 1) frac = 1;
                char el[32], eta[32];
                fmt_hms(elapsed, el, sizeof el);
                if (frac > 1e-6) fmt_hms(elapsed * (1.0 - frac) / frac, eta, sizeof eta);
                else             snprintf(eta, sizeof eta, "--:--");
                fprintf(stderr,
                    "\r%5.1f%% | %.1f M/s | elapsed %s | ETA %s        ",
                    frac * 100.0, cur / 1e6, el, eta);
                fflush(stderr);
                pthread_mutex_unlock(&g_pmx);
            }
        }
        primesieve_free_iterator(&it);
        t->covered += b - a;
        g_pos[t->tid] = t->covered;
    }
    return NULL;
}

static u64 align30(u64 x) { return x - (x % 30); }

typedef struct { u64 first; int ds; const u64 *mem; } Rec;
static int reccmp(const void *a, const void *b) {
    u64 x = ((const Rec *)a)->first, y = ((const Rec *)b)->first;
    return (x > y) - (x < y);
}

/* Build the list of value-intervals that must actually be scanned.
 *
 * Without --emirps that is simply [start,end).  With --emirps we may skip any
 * prime whose LEADING digit is 2,4,5,6,8 (its reversal would end in that digit
 * and so be even or a multiple of 5).  Those primes can never be emirps, and a
 * valid window needs ALL of its primes to be emirps, so no valid window can
 * contain one -- every valid window lies wholly inside a "good" band.
 * Good bands are [d*10^k,(d+1)*10^k) for d in {1,3,7,9}, which is 4/9 of each
 * decade.  Bands [9*10^k,10^(k+1)) and [10^(k+1),2*10^(k+1)) are contiguous, so
 * bands are MERGED -- otherwise a legitimate window straddling that boundary
 * would be split.  Single-digit primes are skipped: they reverse to themselves
 * (palindromes) and are excluded by the emirp definition. */
static void build_intervals(void) {
    if (!g_emirps) {
        g_niv = 1;
        g_ivlo = (u64 *)malloc(sizeof(u64)); g_ivhi = (u64 *)malloc(sizeof(u64));
        g_ivlo[0] = g_start; g_ivhi[0] = g_end;
        return;
    }
    int cap = 4 * 20 + 4;
    u64 *lo = (u64 *)malloc((size_t)cap * sizeof(u64));
    u64 *hi = (u64 *)malloc((size_t)cap * sizeof(u64));
    int n = 0;
    static const int GOOD[4] = {1, 3, 7, 9};
    for (int k = 1; k <= 19; k++) {
        for (int gi = 0; gi < 4; gi++) {
            int d = GOOD[gi];
            __uint128_t blo = (__uint128_t)d * PW10[k];
            __uint128_t bhi = (__uint128_t)(d + 1) * PW10[k];
            if (blo > (__uint128_t)UINT64_MAX) continue;
            if (bhi > (__uint128_t)UINT64_MAX) bhi = (__uint128_t)UINT64_MAX;
            u64 a = (u64)blo, b = (u64)bhi;
            if (a < g_start) a = g_start;
            if (b > g_end)   b = g_end;
            if (a >= b) continue;
            lo[n] = a; hi[n] = b; n++;
        }
    }
    /* merge contiguous / touching bands (generated in ascending order) */
    int m = 0;
    for (int i = 0; i < n; i++) {
        if (m > 0 && lo[i] <= hi[m-1]) { if (hi[i] > hi[m-1]) hi[m-1] = hi[i]; }
        else { lo[m] = lo[i]; hi[m] = hi[i]; m++; }
    }
    g_ivlo = lo; g_ivhi = hi; g_niv = m;
}

int main(int argc, char **argv) {
    g_L = 0; g_base = 10; g_find_all = 0; g_quiet = 0; g_emirps = 0;
    g_start = 0; g_end = 1000000000ULL;
    int nthreads = 1;
    for (int i = 1; i < argc; i++) {
        if      ((!strcmp(argv[i],"--length") ||!strcmp(argv[i],"-n")) && i+1<argc) g_L=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--start")  ||!strcmp(argv[i],"-S")) && i+1<argc) g_start=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--end")    ||!strcmp(argv[i],"-E")) && i+1<argc) g_end=parse_num(argv[++i]);
        else if ((!strcmp(argv[i],"--limit")  ||!strcmp(argv[i],"-L")) && i+1<argc) g_end=parse_num(argv[++i]); /* alias for -E */
        else if ((!strcmp(argv[i],"--base")   ||!strcmp(argv[i],"-b")) && i+1<argc) g_base=atoi(argv[++i]);
        else if ((!strcmp(argv[i],"--threads")||!strcmp(argv[i],"-t")) && i+1<argc) nthreads=atoi(argv[++i]);
        else if ( !strcmp(argv[i],"--all"))    g_find_all=1;
        else if ( !strcmp(argv[i],"--emirps")) g_emirps=1;
        else if ( !strcmp(argv[i],"--quiet")) g_quiet=1;
        else { fprintf(stderr,"bad arg: %s\n", argv[i]); return 2; }
    }
    if (g_L < 1) { fprintf(stderr,"need --length >= 1\n"); return 2; }
    if (g_base < 2 || g_base > 36) { fprintf(stderr,"base 2..36\n"); return 2; }
    if (g_end <= g_start) { fprintf(stderr,"need end (-E) > start (-S)\n"); return 2; }
    if (nthreads < 1) nthreads = 1;
    g_M = (u64)g_base - 1; if (g_M == 0) g_M = 1;
    if (g_base == 10) for (int i = 1; i < 100000; i++) ds5[i] = ds5[i/10] + i%10;

    u64 span = g_end - g_start;
    while (nthreads > 1 && span / (u64)nthreads < (u64)OVERLAP * 8) nthreads--;

    build_intervals();
    u64 measure = 0;
    for (int j = 0; j < g_niv; j++) measure += g_ivhi[j] - g_ivlo[j];
    if (measure == 0) {
        printf("Nothing to scan in [%llu, %llu)%s.\n", (ull)g_start, (ull)g_end,
               g_emirps ? " after skipping disqualified leading digits" : "");
        return 0;
    }
    while (nthreads > 1 && measure / (u64)nthreads < (u64)OVERLAP * 8) nthreads--;

    /* Cut the range so each thread gets an equal share of the SCANNED measure
     * (not of the raw value span), otherwise threads landing in skipped bands
     * would idle while others do all the work. */
    u64 *cut = (u64 *)malloc((size_t)(nthreads + 1) * sizeof(u64));
    cut[0] = g_start; cut[nthreads] = g_end;
    {
        double per = (double)measure / (double)nthreads;
        int i = 1; u64 acc = 0;
        for (int j = 0; j < g_niv && i < nthreads; j++) {
            u64 w = g_ivhi[j] - g_ivlo[j];
            while (i < nthreads && (double)(acc + w) >= per * (double)i) {
                u64 need = (u64)(per * (double)i - (double)acc);
                u64 c = g_ivlo[j] + (need < w ? need : w);
                if (c < cut[i-1]) c = cut[i-1];
                cut[i] = c; i++;
            }
            acc += w;
        }
        while (i < nthreads) { cut[i] = cut[i-1]; i++; }
    }

    Targ *ta = (Targ *)calloc(nthreads, sizeof(Targ));
    g_total_span = 0;
    for (int i = 0; i < nthreads; i++) {
        u64 lo = cut[i], hi = cut[i+1];
        ta[i].lo = lo; ta[i].hi = hi; ta[i].tid = i;
        ta[i].A = (lo > OVERLAP) ? lo - OVERLAP : 0;
        ta[i].B = (hi + OVERLAP < g_end) ? hi + OVERLAP : g_end;
        for (int j = 0; j < g_niv; j++) {      /* measure this thread scans */
            u64 a = g_ivlo[j] > ta[i].A ? g_ivlo[j] : ta[i].A;
            u64 b = g_ivhi[j] < ta[i].B ? g_ivhi[j] : ta[i].B;
            if (a < b) g_total_span += b - a;
        }
    }
    g_nthreads = nthreads;
    g_pos = (u64 *)calloc(nthreads, sizeof(u64));

    g_t_start = g_t_last = mono();
    g_count_last = 0; g_rate = 0;
    pthread_t *th = (pthread_t *)malloc(nthreads * sizeof(pthread_t));
    for (int i = 0; i < nthreads; i++) pthread_create(&th[i], NULL, worker, &ta[i]);
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    if (!g_quiet) {
        char el[32]; fmt_hms(mono() - g_t_start, el, sizeof el);
        fprintf(stderr, "\rcompleted in %s%40s\n", el, "");
    }

    size_t total = 0;
    for (int i = 0; i < nthreads; i++) total += ta[i].count;
    Rec *recs = (Rec *)malloc((total ? total : 1) * sizeof(Rec));
    size_t k = 0;
    for (int i = 0; i < nthreads; i++)
        for (size_t j = 0; j < ta[i].count; j++)
            recs[k++] = (Rec){ ta[i].first[j], ta[i].dsv[j], &ta[i].members[j*(size_t)g_L] };
    qsort(recs, total, sizeof(Rec), reccmp);

    if (total == 0) {
        printf("No run of %d consecutive %sprimes with equal digit sum "
               "(base %d) in [%llu, %llu).\n",
               g_L, g_emirps ? "emirp " : "", g_base, (ull)g_start, (ull)g_end);
    } else {
        size_t show = g_find_all ? total : 1;
        printf("%s of %d consecutive %sprimes sharing a digit sum (base %d), "
               "in [%llu, %llu):\n",
               g_find_all ? "All runs" : "First run", g_L,
               g_emirps ? "emirp " : "", g_base,
               (ull)g_start, (ull)g_end);
        int idxw = 1; for (size_t tt = total; tt >= 10; tt /= 10) idxw++;
        for (size_t r = 0; r < show; r++) {
            if (g_find_all) printf("  #%*zu  digit sum %d: [", idxw, r + 1, recs[r].ds);
            else            printf("  digit sum %d: [", recs[r].ds);
            for (int j = 0; j < g_L; j++)
                printf("%llu%s", (ull)recs[r].mem[j], j+1 < g_L ? ", " : "");
            printf("]\n");
        }
        if (g_find_all)
            printf("Total sequences discovered: %zu\n", total);
    }
    return 0;
}
