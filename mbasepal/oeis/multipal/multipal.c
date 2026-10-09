/*
 * multipal.c — exhaustive search for numbers that are k-digit palindromes
 * in many bases at once (OEIS A171701-A171706, A171741, A171742).
 *
 * For a digit length k, every k-digit palindrome in every base b >= 2 is
 * enumerated over consecutive value chunks [S, E). Within a chunk the
 * number of bases producing each value is counted, either with a u16
 * array (dense mode, used when palindromes are at least ~1/8 as numerous
 * as values, i.e. k = 2, 3) or by radix-sorting the generated values
 * (sparse mode). The search is exhaustive over [lo, hi): for every
 * multiplicity c it reports the least v with exactly c representations,
 * and from those the least v with at least n (the A1717xx definition).
 *
 * A k-digit palindrome in base b with digits d0 d1 .. d1 d0 is
 *   v = sum_i d_i * P_i(b),   P_i = b^(k-1-i) + b^i
 * (P = b^i for the middle digit of odd k), d0 in [1, b-1], others in
 * [0, b-1]. v increases lexicographically in (d0, d1, ...), because the
 * inner digits contribute less than one step of the outer one, so the
 * digits landing in a chunk are found level by level with one division
 * per bound, and the innermost digit sweeps an arithmetic progression.
 *
 *   ./multipal -k 4 -u 27653197824001 [-l 1] [-t threads] [-N keys/chunk]
 *              [-r report_min] [-q]
 *   -l/-u  search lo <= v < hi (hi must be < 2^64)
 *   -r c   print "HIT v c" for every v with >= c representations
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint16_t u16;

#define MAXC 4096
#define MAXLV 5

static int K, LV;                      /* digit length, ceil(K/2) levels */
static u64 LO = 1, HI = 0;
static u64 NCAP = 1ULL << 26;          /* sparse: max keys per chunk */
static u64 DCAP = 1ULL << 27;          /* dense: values per chunk */
static int REPORT = 0, QUIET = 0, NT = 0;

static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static u64 g_min[MAXC + 1];            /* least v with exactly c reps */
static u64 g_hist[MAXC + 1];           /* number of v with exactly c reps */
static _Atomic u64 g_keys = 0;         /* palindromes generated */
static _Atomic u64 g_covered = 0;      /* width of completed chunks */
static _Atomic int g_done = 0;

/* ---------------- chunk queue ---------------- */

typedef struct { u64 S, E; } Chunk;
static Chunk *g_q = NULL;
static size_t g_qn = 0, g_qcap = 0, g_qnext = 0;

static void q_push(u64 S, u64 E){
    if (g_qn == g_qcap){
        g_qcap = g_qcap ? 2 * g_qcap : 1024;
        g_q = realloc(g_q, g_qcap * sizeof *g_q);
        if (!g_q){ fprintf(stderr, "out of memory (queue)\n"); exit(1); }
    }
    g_q[g_qn].S = S; g_q[g_qn].E = E; g_qn++;
}

/* LIFO over the tail so split halves are processed next (they are the
 * lowest unfinished values); the initial list is pushed in reverse. */
static int q_pop(Chunk *c){
    pthread_mutex_lock(&g_mu);
    int ok = 0;
    if (g_qn > g_qnext){ *c = g_q[--g_qn]; ok = 1; }
    pthread_mutex_unlock(&g_mu);
    return ok;
}

/* ---------------- arithmetic ---------------- */

static u128 powsat(u64 b, int e){
    const u128 CAP = (u128)1 << 100;
    u128 r = 1;
    while (e-- > 0){ if (r > CAP / b) return CAP; r *= b; }
    return r;
}

static u64 iroot(u64 x, int k){       /* largest r with r^k <= x */
    if (k <= 1) return x;
    u64 r = (u64)pow((double)x, 1.0 / k);
    while (r > 0 && powsat(r, k) > x) r--;
    while (powsat(r + 1, k) <= x) r++;
    return r;
}

/* palindromes per unit value near v: sum over bases b with
 * b^(K-1) <= v < b^K of b^-(K-LV), via the integral approximation */
static double density(u64 v){
    double x = (double)iroot(v, K) + 1.0, y = (double)iroot(v, K - 1);
    if (y < 2) y = 2;
    if (x < 2) x = 2;
    if (y < x) return 0.0;
    int m = K - LV;
    if (m <= 1) return log((y + 0.5) / (x - 0.5));
    return (pow(x - 0.5, 1 - m) - pow(y + 0.5, 1 - m)) / (m - 1);
}

/* ---------------- generation ---------------- */

typedef struct {
    u64 S, E;
    int dense;
    u64 *buf; u64 n, cap; int overflow;     /* sparse */
    u16 *cnt;                               /* dense */
    u64 P[MAXLV], R[MAXLV], b;
} Ctx;

static void emit(Ctx *C, u64 v0, u64 step, u64 m){
    if (C->dense){
        u16 *c = C->cnt + (v0 - C->S);
        for (u64 j = 0; j < m; j++, c += step) (*c)++;
    } else {
        if (C->n + m > C->cap){ C->overflow = 1; return; }
        u64 *o = C->buf + C->n, k = v0 - C->S;
        for (u64 j = 0; j < m; j++, k += step) o[j] = k;
        C->n += m;
    }
}

static void gen(Ctx *C, int i, u64 acc){
    if (C->overflow) return;
    const u64 P = C->P[i], R = C->R[i];
    u64 dmin = (i == 0) ? 1 : 0, dmax = C->b - 1;
    u128 reach = (u128)acc + R;               /* largest v below this node */
    if ((u128)C->S > reach){
        u64 t = (u64)(((u128)C->S - reach + P - 1) / P);
        if (t > dmin) dmin = t;
    }
    if (acc > C->E - 1) return;
    u64 t2 = (C->E - 1 - acc) / P;
    if (t2 < dmax) dmax = t2;
    if (dmin > dmax) return;
    if (i == LV - 1){
        emit(C, acc + dmin * P, P, dmax - dmin + 1);
        return;
    }
    for (u64 d = dmin; d <= dmax && !C->overflow; d++) gen(C, i + 1, acc + d * P);
}

static void gen_chunk(Ctx *C){
    u64 blo = iroot(C->S, K) + 1;             /* smallest b with b^K > S */
    u64 bhi = iroot(C->E - 1, K - 1);         /* largest b with b^(K-1) <= E-1 */
    if (blo < 2) blo = 2;
    for (u64 b = blo; b <= bhi && !C->overflow; b++){
        C->b = b;
        for (int i = 0; i < LV; i++){
            int j = K - 1 - i;
            u128 p = (i == j) ? powsat(b, i) : powsat(b, j) + powsat(b, i);
            C->P[i] = (u64)p;
        }
        u128 r = 0;
        for (int i = LV - 1; i >= 0; i--){ C->R[i] = (u64)r; r += (u128)(b - 1) * C->P[i]; }
        gen(C, 0, 0);
    }
}

/* LSD radix sort, 11-bit digits over the low `bits`; returns sorted array */
static u64 *radix(u64 *a, u64 *t, u64 n, int bits){
    static _Thread_local u64 cnt[2048];
    for (int sh = 0; sh < bits; sh += 11){
        memset(cnt, 0, sizeof cnt);
        for (u64 i = 0; i < n; i++) cnt[(a[i] >> sh) & 2047]++;
        u64 s = 0;
        for (int j = 0; j < 2048; j++){ u64 c = cnt[j]; cnt[j] = s; s += c; }
        for (u64 i = 0; i < n; i++) t[cnt[(a[i] >> sh) & 2047]++] = a[i];
        u64 *x = a; a = t; t = x;
    }
    return a;
}

static void record(u64 *tmin, u64 *thist, u64 v, unsigned c){
    if (c > MAXC) c = MAXC;
    thist[c]++;
    if (v < tmin[c]) tmin[c] = v;
    if (REPORT && (int)c >= REPORT){
        pthread_mutex_lock(&g_mu);
        printf("HIT %llu %u\n", (unsigned long long)v, c);
        fflush(stdout);
        pthread_mutex_unlock(&g_mu);
    }
}

static void *worker(void *arg){
    (void)arg;
    u64 *buf = NULL, *tmp = NULL; u16 *cnt = NULL;
    u64 *tmin = malloc((MAXC + 1) * sizeof *tmin);
    u64 *thist = malloc((MAXC + 1) * sizeof *thist);
    Chunk ch;
    while (q_pop(&ch)){
        Ctx C; memset(&C, 0, sizeof C);
        C.S = ch.S; C.E = ch.E;
        C.dense = density(ch.S) >= 0.125;
        if (C.dense && ch.E - ch.S > DCAP){
            pthread_mutex_lock(&g_mu);
            q_push(ch.S + DCAP, ch.E); q_push(ch.S, ch.S + DCAP);
            pthread_mutex_unlock(&g_mu);
            continue;
        }
        if (C.dense){
            if (!cnt && !(cnt = malloc(DCAP * sizeof *cnt))){ fprintf(stderr, "oom\n"); exit(1); }
            C.cnt = cnt;
            memset(cnt, 0, (ch.E - ch.S) * sizeof *cnt);
        } else {
            if (!buf){
                buf = malloc(NCAP * sizeof *buf); tmp = malloc(NCAP * sizeof *tmp);
                if (!buf || !tmp){ fprintf(stderr, "oom\n"); exit(1); }
            }
            C.buf = buf; C.cap = NCAP;
        }
        gen_chunk(&C);
        if (C.overflow){
            /* density estimate was low here: split and retry */
            u64 mid = ch.S + (ch.E - ch.S) / 2;
            pthread_mutex_lock(&g_mu);
            q_push(mid, ch.E); q_push(ch.S, mid);
            pthread_mutex_unlock(&g_mu);
            continue;
        }
        for (int c = 0; c <= MAXC; c++){ tmin[c] = UINT64_MAX; thist[c] = 0; }
        u64 keys;
        if (C.dense){
            u64 w = ch.E - ch.S; keys = 0;
            for (u64 i = 0; i < w; i++){
                unsigned c = cnt[i];
                if (c){ keys += c; record(tmin, thist, ch.S + i, c); }
            }
        } else {
            keys = C.n;
            int bits = 64 - __builtin_clzll((ch.E - ch.S - 1) | 1);
            u64 *s = radix(buf, tmp, C.n, bits);
            for (u64 i = 0; i < C.n; ){
                u64 j = i + 1;
                while (j < C.n && s[j] == s[i]) j++;
                record(tmin, thist, ch.S + s[i], (unsigned)(j - i));
                i = j;
            }
        }
        pthread_mutex_lock(&g_mu);
        for (int c = 1; c <= MAXC; c++){
            g_hist[c] += thist[c];
            if (tmin[c] < g_min[c]) g_min[c] = tmin[c];
        }
        pthread_mutex_unlock(&g_mu);
        atomic_fetch_add(&g_keys, keys);
        atomic_fetch_add(&g_covered, ch.E - ch.S);
    }
    free(buf); free(tmp); free(cnt); free(tmin); free(thist);
    return NULL;
}

static double now_s(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void *status(void *arg){
    double t0 = *(double *)arg;
    while (!atomic_load(&g_done)){
        for (int i = 0; i < 100 && !atomic_load(&g_done); i++) usleep(100000);
        if (atomic_load(&g_done)) break;
        double el = now_s() - t0;
        double cov = (double)atomic_load(&g_covered) / (double)(HI - LO);
        u64 keys = atomic_load(&g_keys);
        int top = 0; u64 topv = 0;
        pthread_mutex_lock(&g_mu);
        for (int c = MAXC; c >= 1; c--) if (g_hist[c]){ top = c; topv = g_min[c]; break; }
        pthread_mutex_unlock(&g_mu);
        fprintf(stderr, "[k=%d] %.3f%% of values | %.3g palindromes, %.1fM/s | max reps %d (least %llu) | %.0fs%s\n",
                K, 100.0 * cov, (double)keys, keys / el / 1e6, top,
                (unsigned long long)topv, el, "");
    }
    return NULL;
}

int main(int argc, char **argv){
    int opt;
    while ((opt = getopt(argc, argv, "k:l:u:t:N:D:r:q")) != -1){
        switch (opt){
        case 'k': K = atoi(optarg); break;
        case 'l': LO = strtoull(optarg, NULL, 10); break;
        case 'u': HI = strtoull(optarg, NULL, 10); break;
        case 't': NT = atoi(optarg); break;
        case 'N': NCAP = strtoull(optarg, NULL, 10); break;
        case 'D': DCAP = strtoull(optarg, NULL, 10); break;
        case 'r': REPORT = atoi(optarg); break;
        case 'q': QUIET = 1; break;
        default:
            fprintf(stderr, "usage: %s -k digits -u hi [-l lo] [-t threads] [-N keys] [-D values] [-r c] [-q]\n", argv[0]);
            return 2;
        }
    }
    if (K < 2 || K > 2 * MAXLV || HI <= LO || LO < 1){
        fprintf(stderr, "need 2 <= k <= %d and 1 <= lo < hi\n", 2 * MAXLV);
        return 2;
    }
    LV = (K + 1) / 2;
    if (NT <= 0) NT = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (NT > 256) NT = 256;
    for (int c = 0; c <= MAXC; c++){ g_min[c] = UINT64_MAX; g_hist[c] = 0; }

    /* chunk list, pushed in reverse so the lowest values pop first */
    {
        Chunk *tmpq = NULL; size_t n = 0, cap = 0;
        for (u64 S = LO; S < HI; ){
            double d = density(S);
            u64 W;
            if (d >= 0.125) W = DCAP;
            else {
                double w = 0.8 * (double)NCAP / (d > 0 ? d : 1e-30);
                W = w > 1e18 ? (u64)1e18 : (u64)w;
                if (W < 1024) W = 1024;
            }
            u64 E = (HI - S > W) ? S + W : HI;
            if (n == cap){ cap = cap ? 2 * cap : 1024; tmpq = realloc(tmpq, cap * sizeof *tmpq); }
            tmpq[n].S = S; tmpq[n].E = E; n++;
            S = E;
        }
        for (size_t i = n; i-- > 0; ) q_push(tmpq[i].S, tmpq[i].E);
        free(tmpq);
        fprintf(stderr, "multipal: k=%d, v in [%llu, %llu), %zu chunks, %d threads\n",
                K, (unsigned long long)LO, (unsigned long long)HI, n, NT);
    }

    double t0 = now_s();
    pthread_t st, tid[256];
    if (!QUIET) pthread_create(&st, NULL, status, &t0);
    for (int i = 0; i < NT; i++) pthread_create(&tid[i], NULL, worker, NULL);
    for (int i = 0; i < NT; i++) pthread_join(tid[i], NULL);
    atomic_store(&g_done, 1);
    if (!QUIET) pthread_join(st, NULL);

    double el = now_s() - t0;
    int top = 0;
    for (int c = MAXC; c >= 1; c--) if (g_hist[c]){ top = c; break; }
    printf("# multipal k=%d: all v in [%llu, %llu) searched, %llu palindromes, %.1fs\n",
           K, (unsigned long long)LO, (unsigned long long)HI,
           (unsigned long long)atomic_load(&g_keys), el);
    printf("# c  least_v_with_exactly_c  number_of_such_v\n");
    for (int c = 1; c <= top; c++)
        if (g_hist[c])
            printf("EXACT %d %llu %llu\n", c, (unsigned long long)g_min[c], (unsigned long long)g_hist[c]);
    printf("# n  least_v_with_at_least_n  (a(n) if lo = 1)\n");
    u64 run = UINT64_MAX;
    u64 *atl = malloc((top + 2) * sizeof *atl);
    for (int c = top; c >= 1; c--){ if (g_min[c] < run) run = g_min[c]; atl[c] = run; }
    for (int c = 1; c <= top; c++) printf("ATLEAST %d %llu\n", c, (unsigned long long)atl[c]);
    free(atl);
    return 0;
}
