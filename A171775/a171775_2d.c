/*
 * a171775_2d.c -- faster search for OEIS A171775 ("swapped 2D" method).
 *
 *   a(n) = smallest M such that for every k = 2..n there is a base b_k in
 *   which M is a k-digit palindrome.  Known: a(7..10) = 2^30, 2^42, 2^56,
 *   2^72 (a(10) by a171775.c, Sep 2026); conjecture a(n) = 2^((n-1)(n-2)).
 *
 * a171775.c enumerates the n-digit palindromes and solves one digit.  For odd
 * n the (n-1)-digit palindromes are the smaller set (below 2^90 there are
 * 7.0e15 ten-digit ones and 6.1e16 eleven-digit ones), and solving two digits
 * at once removes a factor c^2 instead of B.  This tool therefore
 *
 *   * enumerates the Le = (n-1)-digit palindromes M in base c:
 *     free digits e_0..e_{De-1} with weights w_i = c^(Le-1-i) + c^i
 *     (w_mid = c^(Le/2) for odd Le); x = e_{De-3}, y = e_{De-2}, z = e_{De-1},
 *     block = e_0..e_{De-4}, M = M1 + y w_y + z w_z with M1 = Mblock + x w_x;
 *   * for every base B < c in which M can have Lc = n digits
 *     (B^(n-1) <= M < B^n) imposes "the k lowest base-B digits mirror the k
 *     highest" as a linear congruence in (y, z) and finds all its solutions
 *     in the box [0, c)^2 with a lookup table.
 *
 * Per base pair (c, B) choose the largest k <= 3 with max(span, w_x) < P =
 * B^(n-k), span = (c-1)(w_y + w_z).  Then over the (y, z) box q = floor(M/P)
 * (the k leading base-B digits of M) takes at most the two values q1, q1+1 of
 * M1, and x -> x+1 raises q1 by at most one.  M must satisfy M = T(q) mod N,
 * N = B^k, T = the digits of q reversed.  With g = gcd(w_z, N), N' = N/g,
 * u = inverse of a unit s'' = w_z/g (mod N') and R = w_y u mod N:
 *
 *   y w_y + z w_z = T - M1 (mod N)   <=>   V(y) = (W - y R) mod N,  g | V(y),
 *   z = V(y)/g + i N',   where W = (T - M1) u mod N.
 *
 * So the y that give some z < c are those with y R mod N in the cyclic window
 * (W - min(g c, N), W]; a table of y R mod N sorted into ~c buckets finds them
 * in O(1).  W, M1 mod P and the digits of q1 are updated with additions only
 * as x runs (a carry stopping at digit m moves T by B^(m+1) + B^m mod B^k, as
 * in a171775.c); each x costs one lookup, or two when the box straddles q1+1.
 *
 * A match is filtered by the next digit pair (digit k from M mod B^(k+1),
 * digit n-1-k from the remainder with a floating-point quotient and an exact
 * fallback), then checked in full.  Survivors (n-digit palindrome in B and
 * (n-1)-digit palindrome in c) go through the same canonical-base rule and
 * length n-2..3 tests as a171775.c, so the statistics (survivor count,
 * checksum, failures by length) are directly comparable with it.
 *
 * Usage
 *   a171775_2d selftest [full]        naive comparisons; n = 8, 9 [, 10] full searches vs a171775.c
 *   a171775_2d search n LO HI [-t T] [-S state] [-i sec] [-v]
 *   a171775_2d paranoid [reps]        (build with -DPARANOID) state checked at every step
 * Numbers accept 2^90, 1e27, 2^90-1, ...
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <math.h>
#include <pthread.h>
#include <signal.h>
#ifndef A171775_2D_NO_DRIVER
#include <stdatomic.h>
#endif
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef int64_t i64;
typedef uint32_t u32;
typedef uint16_t u16;

#define MAXN 24
#define MAXK 4
#define KMAX 3
#define U128_MAX (~(u128)0)

/* ------------------------------------------------------------------ */
/* integer helpers (as in a171775.c)                                   */
/* ------------------------------------------------------------------ */

static u128 pow_sat(u64 b, int e)
{
    u128 r = 1;
    for (int i = 0; i < e; i++) {
        if (b && r > U128_MAX / b) return U128_MAX;
        r *= b;
    }
    return r;
}

static u64 iroot(u128 x, int k)
{
    if (x == 0) return 0;
    double a = pow((double)x, 1.0 / k);
    u64 r = a >= 1.8e19 ? (u64)1.8e19 : (u64)a;
    if (r == 0) r = 1;
    while (r > 1 && pow_sat(r, k) > x) r--;
    while (pow_sat(r + 1, k) <= x) r++;
    return r;
}

static char *u128s(u128 x, char *buf)
{
    char t[48];
    int i = 0;
    do { t[i++] = (char)('0' + (int)(x % 10)); x /= 10; } while (x);
    for (int j = 0; j < i; j++) buf[j] = t[i - 1 - j];
    buf[i] = 0;
    return buf;
}

static int parse_term(const char **ps, u128 *out)
{
    const char *s = *ps;
    u128 v = 0;
    int nd = 0;
    while (*s >= '0' && *s <= '9') { v = v * 10 + (u128)(*s - '0'); s++; nd++; }
    if (!nd) return 0;
    if (*s == '^' || *s == 'e' || *s == 'E') {
        char op = *s++;
        char *e;
        unsigned long ex = strtoul(s, &e, 10);
        if (e == s) return 0;
        s = e;
        if (op == '^') {
            if (v > UINT64_MAX) return 0;
            v = pow_sat((u64)v, (int)ex);
        } else v *= pow_sat(10, (int)ex);
    }
    *ps = s;
    *out = v;
    return 1;
}

static u128 parse_or_die(const char *s0)
{
    const char *s = s0;
    u128 acc, t;
    if (!parse_term(&s, &acc)) goto bad;
    while (*s == '+' || *s == '-') {
        char op = *s++;
        if (!parse_term(&s, &t)) goto bad;
        acc = op == '+' ? acc + t : acc - t;
    }
    if (*s) goto bad;
    return acc;
bad:
    fprintf(stderr, "bad number: %s\n", s0);
    exit(1);
}

static u64 gcd64(u64 a, u64 b)
{
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

static u64 modinv64(u64 a, u64 m)
{
    if (m == 1) return 0;
    __int128 t = 0, nt = 1, r = m, nr = a % m;
    while (nr) {
        __int128 q = r / nr, tmp;
        tmp = t - q * nt; t = nt; nt = tmp;
        tmp = r - q * nr; r = nr; nr = tmp;
    }
    if (r != 1) return 0;
    if (t < 0) t += m;
    return (u64)t;
}

static inline u64 mulmod64(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }

static inline u64 mod128_64(u128 M, u64 b)
{
    u64 hi = (u64)(M >> 64), lo = (u64)M;
    if (!hi) return lo % b;
    if (b <= 0xffffffffULL) {
        u64 r = hi % b;
        r = ((r << 32) | (lo >> 32)) % b;
        r = ((r << 32) | (lo & 0xffffffffULL)) % b;
        return r;
    }
    return (u64)(M % b);
}

static inline double u128d(u128 x) { return (double)(u64)(x >> 64) * 18446744073709551616.0 + (double)(u64)x; }

static int is_pal_base(u128 M, u64 b, int L)
{
    u64 d[130];
    int len = 0;
    if (b < 2) return 0;
    if ((M >> 64) == 0) {
        u64 x = (u64)M;
        while (x) { if (len >= L) return 0; d[len++] = x % b; x /= b; }
    } else {
        u128 x = M;
        while (x) { if (len >= L) return 0; d[len++] = (u64)(x % b); x /= b; }
    }
    if (len != L) return 0;
    for (int i = 0; i < L / 2; i++)
        if (d[i] != d[L - 1 - i]) return 0;
    return 1;
}

static u64 find_pal_base(u128 M, int L)
{
    if (L < 2 || M == 0) return 0;
    if (L == 2) {
        if (M < 3) return 0;
        return (M - 1 > UINT64_MAX) ? UINT64_MAX : (u64)(M - 1);
    }
    u64 blo = iroot(M, L) + 1, bhi = iroot(M, L - 1);
    if (blo < 2) blo = 2;
    double Md = (double)M;
    for (u64 b = blo; b <= bhi; b++) {
        u64 r = mod128_64(M, b);
        if (r == 0) continue;
        double la = Md / pow((double)b, (double)(L - 1));
        if (la < (double)r - 1.5 || la > (double)r + 1.5) continue;
        u128 P = pow_sat(b, L - 1);
        if (M / P != r) continue;
        if (is_pal_base(M, b, L)) return b;
    }
    return 0;
}

static inline u64 surv_hash(u128 M)
{
    u64 a = (u64)M * 0x9E3779B97F4A7C15ULL, b = (u64)(M >> 64) * 0xC2B2AE3D27D4EB4FULL;
    a ^= a >> 29;
    return a * 0xBF58476D1CE4E5B9ULL + b;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ------------------------------------------------------------------ */
/* problem and per-base constants                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, Le, Lc, he, De, bd, ud; /* check length Lc = n, enumerated length Le = n-1 */
    u128 lo, hi;
    u64 cmin, cmax;                /* enumeration bases */
    u64 *uoff;                     /* uoff[c-cmin] = first unit of base c */
    u64 nunits;
} prob_t;

typedef struct {                   /* per (c, B) */
    u64 B, k, N, g, Np, u, lim, R, Wx, Wc[MAXK];
    u128 P, Rc;                    /* P = B^(n-k); crossing iff rem >= Rc = P - span */
    u64 BK1, oxK, oyK, ozK;        /* B^(k+1) and w_x, w_y, w_z mod B^(k+1) */
    u64 Bk;                        /* B^k */
    u128 Ptop;                     /* B^(n-1-k): digit n-1-k = floor(r / Ptop) */
    double invPtop;
    u32 *val;                      /* table: key(y R mod N) sorted, with y */
    u16 *ys;
    u32 *bst;                      /* bucket starts, nb + 1 entries */
    int sh;
    u32 nb;
    u32 *ay_lo, *az_lo;            /* t w_y, t w_z mod B^(k+1) for t < c, split as (mod B^k, div B^k) */
    u16 *ay_hi, *az_hi;
    u64 ox_lo, ox_hi;              /* w_x mod B^(k+1), split the same way */
} bconst_t;

typedef struct {                   /* per enumeration base c */
    u64 c;
    u128 w[MAXN];
    u128 wx, wy, wz, span, span_block, span_unit;
    u64 Bbase, Bcount;
    bconst_t *tab;
    size_t tabcap;
    char *mem;                     /* table storage */
    size_t memcap;
} cctx_t;

typedef struct { u128 M; u64 B; } pair_t;

typedef struct {
    const prob_t *pr;
    u64 steps, lookups, hits, f1, st2, canon, csum, sol;
    u64 fail[MAXN + 1];
    int record;                    /* test hook: 1 = full n-digit palindromes, 2 = all k-digit matches */
    pair_t *pairs;
    size_t np, pcap;
    int verbose;
#ifdef PARANOID
    u64 pchk, pz;                  /* brute-force steps checked, (y,z) hits checked */
#endif
} wstat_t;

static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
static FILE *sol_out;

static void prob_init(prob_t *pr, int n, u128 lo, u128 hi)
{
    memset(pr, 0, sizeof *pr);
    pr->n = n;
    pr->Lc = n;
    pr->Le = n - 1;
    pr->he = pr->Le / 2;
    pr->De = (pr->Le + 1) / 2;
    pr->bd = pr->De - 3;
    pr->ud = pr->bd - 1 >= 1 ? pr->bd - 1 : 1;
    if (lo < 1) lo = 1;
    pr->lo = lo;
    pr->hi = hi;
    /* enumeration bases c with c^(Le-1) <= hi and c^Le - 1 >= lo */
    u64 cmax = iroot(hi, pr->Le - 1), cmin = 3;
    while (cmin <= cmax && pow_sat(cmin, pr->Le) - 1 < lo) cmin++;
    pr->cmin = cmin;
    pr->cmax = cmax;
    u64 nc = cmax >= cmin ? cmax - cmin + 1 : 0;
    pr->uoff = (u64 *)calloc(nc + 1, sizeof(u64));
    u64 tot = 0;
    for (u64 i = 0; i < nc; i++) {
        u64 c = cmin + i;
        pr->uoff[i] = tot;
        tot += (u64)((u128)(c - 1) * pow_sat(c, pr->ud - 1));
    }
    pr->uoff[nc] = tot;
    pr->nunits = tot;
}

static int u32cmp_pairs(const void *a, const void *b)
{
    const u64 *x = (const u64 *)a, *y = (const u64 *)b;   /* packed (val << 16 | y) */
    return *x < *y ? -1 : *x > *y;
}

/* constants and lookup table for the base pair (c, B) */
static void bconst_set(bconst_t *K, const prob_t *pr, const cctx_t *cc, u64 B, char **memp)
{
    const u64 c = cc->c;
    const int Lc = pr->Lc;
    u128 need = cc->span > cc->wx ? cc->span : cc->wx;
    int k = 0;
    while (k + 1 <= KMAX && need < pow_sat(B, Lc - (k + 1))) k++;
    if (k < 1) {
        fprintf(stderr, "internal: no usable k for c=%llu B=%llu\n", (unsigned long long)c, (unsigned long long)B);
        exit(1);
    }
    memset(K, 0, sizeof *K);
    K->B = B;
    K->k = (u64)k;
    u128 N128 = pow_sat(B, k);
    if (N128 >= ((u128)1 << 32)) {
        fprintf(stderr, "B^k too large for 32-bit tables (B=%llu k=%d)\n", (unsigned long long)B, k);
        exit(1);
    }
    u64 N = (u64)N128;
    K->N = N;
    K->Bk = N;
    K->P = pow_sat(B, Lc - k);
    K->Rc = K->P - cc->span;
    u64 sN = (u64)(cc->wz % N), g = gcd64(sN, N), Np = N / g;
    u64 sp = Np > 1 ? (sN / g) % Np : 0, spp = sp;
    while (gcd64(spp % N, N) != 1) spp += Np;
    K->g = g;
    K->Np = Np;
    K->u = modinv64(spp % N, N);
    K->lim = g * c < N ? g * c : N;
    K->R = mulmod64((u64)(cc->wy % N), K->u, N);
    K->Wx = mulmod64((u64)(cc->wx % N), K->u, N);
    u64 bj = 1;
    for (int j = 0; j < k; j++) { K->Wc[j] = mulmod64(bj % N, K->u, N); bj *= B; }
    K->BK1 = (u64)pow_sat(B, k + 1);
    K->oxK = (u64)(cc->wx % K->BK1);
    K->oyK = (u64)(cc->wy % K->BK1);
    K->ozK = (u64)(cc->wz % K->BK1);
    K->Ptop = pow_sat(B, Lc - 1 - k);
    K->invPtop = 1.0 / u128d(K->Ptop);
    /* table */
    K->val = (u32 *)*memp; *memp += c * sizeof(u32);
    K->ys = (u16 *)*memp; *memp += ((c * sizeof(u16) + 7) & ~(size_t)7);
    int sh = 0;
    while ((N >> sh) > c) sh++;
    K->sh = sh;
    K->nb = (u32)(((N - 1) >> sh) + 1);
    K->bst = (u32 *)*memp; *memp += ((K->nb + 1) * sizeof(u32) + 7) & ~(size_t)7;
    K->ay_lo = (u32 *)*memp; *memp += c * sizeof(u32);
    K->az_lo = (u32 *)*memp; *memp += c * sizeof(u32);
    K->ay_hi = (u16 *)*memp; *memp += ((c * sizeof(u16) + 7) & ~(size_t)7);
    K->az_hi = (u16 *)*memp; *memp += ((c * sizeof(u16) + 7) & ~(size_t)7);
    for (u64 t = 0; t < c; t++) {
        u64 vy = (u64)((u128)t * K->oyK % K->BK1), vz = (u64)((u128)t * K->ozK % K->BK1);
        K->ay_lo[t] = (u32)(vy % K->Bk); K->ay_hi[t] = (u16)(vy / K->Bk);
        K->az_lo[t] = (u32)(vz % K->Bk); K->az_hi[t] = (u16)(vz / K->Bk);
    }
    K->ox_lo = K->oxK % K->Bk;
    K->ox_hi = K->oxK / K->Bk;
    /* key(v) = (v mod g) N' + v div g: the y whose v = y R mod N gives a solution z < c for target W
       are exactly those with v = W (mod g) and (W div g - v div g) mod N' < c, i.e. a window of width
       min(c, N') inside one residue class, whatever g is */
    u64 *tmp = (u64 *)malloc(c * sizeof(u64));
    u64 v = 0;
    for (u64 y = 0; y < c; y++) {
        u64 key = (v % g) * Np + v / g;
        tmp[y] = (key << 16) | y;
        v += K->R;
        if (v >= N) v -= N;
    }
    qsort(tmp, c, sizeof(u64), u32cmp_pairs);
    for (u64 i = 0; i < c; i++) { K->val[i] = (u32)(tmp[i] >> 16); K->ys[i] = (u16)(tmp[i] & 0xffff); }
    free(tmp);
    u32 b = 0;
    for (u64 i = 0; i < c; i++) {
        u32 bi = K->val[i] >> sh;
        while (b <= bi) K->bst[b++] = (u32)i;
    }
    while (b <= K->nb) K->bst[b++] = (u32)c;
}

static void cctx_set(cctx_t *cc, const prob_t *pr, u64 c)
{
    int Le = pr->Le, he = pr->he, De = pr->De, Lc = pr->Lc;
    cc->c = c;
    for (int i = 0; i < he; i++) cc->w[i] = pow_sat(c, Le - 1 - i) + pow_sat(c, i);
    if (Le & 1) cc->w[he] = pow_sat(c, he);
    cc->wx = cc->w[De - 3];
    cc->wy = cc->w[De - 2];
    cc->wz = cc->w[De - 1];
    cc->span = (u128)(c - 1) * (cc->wy + cc->wz);
    cc->span_block = (u128)(c - 1) * (cc->wx + cc->wy + cc->wz);
    u128 su = 0;
    for (int i = pr->ud; i < De; i++) su += cc->w[i];
    cc->span_unit = (u128)(c - 1) * su;
    if (c > 65535) { fprintf(stderr, "c too large\n"); exit(1); }
    /* B range for this c */
    u128 Mlo = pow_sat(c, Le - 1), Mhi = pow_sat(c, Le) - 1;
    if (Mlo < pr->lo) Mlo = pr->lo;
    if (Mhi > pr->hi) Mhi = pr->hi;
    u64 Bmin = iroot(Mlo, Lc) + 1, Bmax = iroot(Mhi, Lc - 1);
    if (Bmin < 2) Bmin = 2;
    if (Bmax >= c) Bmax = c - 1;
    if (Bmax < Bmin) { cc->Bbase = Bmin; cc->Bcount = 0; return; }
    cc->Bbase = Bmin;
    cc->Bcount = Bmax - Bmin + 1;
    if (cc->Bcount > cc->tabcap) {
        free(cc->tab);
        cc->tabcap = cc->Bcount;
        cc->tab = (bconst_t *)malloc(cc->tabcap * sizeof(bconst_t));
    }
    size_t need = cc->Bcount * (c * 4 + c * 2 + 8 + (c + 2) * 4 + 16 + c * 12 + 32) + 64;
    if (need > cc->memcap) {
        free(cc->mem);
        cc->memcap = need;
        cc->mem = (char *)malloc(need);
    }
    char *mp = cc->mem;
    for (u64 i = 0; i < cc->Bcount; i++) bconst_set(&cc->tab[i], pr, cc, Bmin + i, &mp);
}

/* ------------------------------------------------------------------ */
/* survivors                                                           */
/* ------------------------------------------------------------------ */

static void report_solution(wstat_t *st, u128 M, u64 B, u64 c)
{
    const prob_t *pr = st->pr;
    char buf[48];
    u64 bases[MAXN + 1];
    for (int L = 3; L <= pr->n; L++) bases[L] = find_pal_base(M, L);
    pthread_mutex_lock(&out_lock);
    fprintf(sol_out, "SOLUTION n=%d M=%s", pr->n, u128s(M, buf));
    fprintf(sol_out, " bases(L=2..%d)=%s", pr->n, u128s(M - 1, buf));
    for (int L = 3; L <= pr->n; L++) fprintf(sol_out, ",%llu", (unsigned long long)bases[L]);
    fprintf(sol_out, "  [B=%llu c=%llu]\n", (unsigned long long)B, (unsigned long long)c);
    fflush(sol_out);
    pthread_mutex_unlock(&out_lock);
}

static void record_pair(wstat_t *st, u128 M, u64 B)
{
    if (st->np == st->pcap) {
        st->pcap = st->pcap ? 2 * st->pcap : 1024;
        st->pairs = (pair_t *)realloc(st->pairs, st->pcap * sizeof(pair_t));
    }
    st->pairs[st->np].M = M;
    st->pairs[st->np].B = B;
    st->np++;
}

static void survivor(wstat_t *st, u128 M, u64 B, u64 c);

/* M is an n-digit palindrome in base B and an (n-1)-digit palindrome in base c */
static void survivor(wstat_t *st, u128 M, u64 B, u64 c)
{
    const prob_t *pr = st->pr;
    const int n = pr->n;
    st->st2++;
    if (st->record == 1) { record_pair(st, M, B); return; }
    if (find_pal_base(M, n - 1) != c) return;          /* canonical: smallest bases of both lengths */
    if (find_pal_base(M, n) != B) return;
    st->canon++;
    st->csum += surv_hash(M);
    for (int L = n - 2; L >= 3; L--) {
        if (!find_pal_base(M, L)) {
            st->fail[L]++;
            if (st->verbose && L <= n - 3) {
                char buf[48];
                pthread_mutex_lock(&out_lock);
                fprintf(sol_out, "DEEP n=%d M=%s fails at L=%d [B=%llu c=%llu]\n", n, u128s(M, buf), L,
                        (unsigned long long)B, (unsigned long long)c);
                fflush(sol_out);
                pthread_mutex_unlock(&out_lock);
            }
            return;
        }
    }
    st->sol++;
    report_solution(st, M, B, c);
}

/* ------------------------------------------------------------------ */
/* one block, one check base B                                         */
/* ------------------------------------------------------------------ */

#ifdef PARANOID
static int u64cmp(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

static u64 T_of(u64 q, u64 B, int k)
{
    u64 E[MAXK], T = 0;
    for (int j = k - 1; j >= 0; j--) { E[j] = q % B; q /= B; }
    for (int j = k - 1; j >= 0; j--) T = T * B + E[j];
    return T;
}
#endif

static void block_B(wstat_t *st, const cctx_t *cc, const bconst_t *K, u128 Mblock, int inside)
{
    const prob_t *pr = st->pr;
    const int n = pr->n;
    const u64 c = cc->c, B = K->B, N = K->N, g = K->g, Np = K->Np, Wx = K->Wx;
    const int k = (int)K->k;
    const u128 P = K->P, Rc = K->Rc, wx = cc->wx, wy = cc->wy, wz = cc->wz;
    const u64 BK1 = K->BK1, Bk = K->Bk, ox_lo = K->ox_lo, ox_hi = K->ox_hi;
    const u128 Ptop = K->Ptop;
    const double invPtop = K->invPtop;
    const u32 *ay_lo = K->ay_lo, *az_lo = K->az_lo;
    const u16 *ay_hi = K->ay_hi, *az_hi = K->az_hi;
    const u32 *val = K->val, *bst = K->bst;
    const u16 *ys = K->ys;
    const int sh = K->sh;

    u128 q128 = Mblock / P;
    if (q128 >= N) return;                 /* M >= B^n for the whole block */
    u128 rem = Mblock % P;
    u64 q = (u64)q128, E[MAXK];
    for (int j = k - 1; j >= 0; j--) { E[j] = q % B; q /= B; }
    u64 T = 0;
    for (int j = k - 1; j >= 0; j--) T = T * B + E[j];
    u64 MN = (u64)(Mblock % N);
    u64 W = mulmod64(T >= MN ? T - MN : T + N - MN, K->u, N);
    u64 mk_lo, mk_hi;                      /* M1 mod B^(k+1) = mk_hi B^k + mk_lo */
    {
        u64 m = (u64)(Mblock % BK1);
        mk_lo = m % Bk;
        mk_hi = m / Bk;
    }
    u128 M1 = Mblock;
    u64 wnext;
    int last;
    /* q -> q+1 with the carry stopping at digit m changes T by B^(m+1) + B^m (mod B^k), or by B^(k-1) */
#define RECOMPUTE_NEXT()                                                        \
    do {                                                                        \
        int jj = k - 1;                                                         \
        while (jj >= 0 && E[jj] == B - 1) jj--;                                 \
        last = jj < 0;                                                          \
        if (jj == k - 1) wnext = K->Wc[k - 1];                                  \
        else if (jj >= 0) {                                                     \
            wnext = K->Wc[jj + 1] + K->Wc[jj];                                  \
            if (wnext >= N) wnext -= N;                                         \
        } else wnext = 0;                                                       \
    } while (0)
    RECOMPUTE_NEXT();

    u64 nlook = 0;
#ifdef PARANOID
    u64 *pfound = malloc((c * c + 1) * sizeof(u64)), *want = malloc((c * c + 1) * sizeof(u64));
    const int pcap = (int)(c * c + 1);
#endif
    for (u64 x = 0; x < c; x++) {
#ifdef PARANOID
        {
            u128 M1p = Mblock + (u128)x * wx;
            u128 qp = M1p / P;
            int bad = 0;
            if (qp >= N) bad |= 1;
            else {
                if (M1p % P != rem) bad |= 2;
                u64 qq = (u64)qp, Et[MAXK];
                for (int j = k - 1; j >= 0; j--) { Et[j] = qq % B; qq /= B; }
                for (int j = 0; j < k; j++) if (Et[j] != E[j]) bad |= 4;
                u64 Tt = T_of((u64)qp, B, k), MNt = (u64)(M1p % N);
                if (mulmod64(Tt >= MNt ? Tt - MNt : Tt + N - MNt, K->u, N) != W) bad |= 8;
                if ((u64)qp + 1 < N) {
                    u64 T1 = T_of((u64)qp + 1, B, k);
                    if (mulmod64(T1 >= Tt ? T1 - Tt : T1 + N - Tt, K->u, N) != wnext) bad |= 16;
                    if (last) bad |= 32;
                } else if (!last) bad |= 32;
                if ((u64)(M1p % BK1) != mk_lo + mk_hi * Bk || M1p != M1) bad |= 64;
            }
            if (bad) {
                char buf[48];
                fprintf(stderr, "PARANOID state mismatch %#x: n=%d c=%llu B=%llu x=%llu M1=%s\n", bad, st->pr->n,
                        (unsigned long long)c, (unsigned long long)B, (unsigned long long)x, u128s(M1p, buf));
                abort();
            }
        }
        int npf = 0;
#endif
        /* the two lookups: target T(q1) (region r < P) and, if the box reaches it, T(q1+1) */
        for (int reg = 0; reg < 2; reg++) {
            u64 Wt;
            if (reg == 0) Wt = W;
            else {
                if (last || rem < Rc) break;
                Wt = W + wnext;
                if (Wt >= N) Wt -= N;
            }
            nlook++;
            u64 r = 0, mW = Wt;
            if (g > 1) { r = Wt % g; mW = Wt / g; }
            const u64 base = r * Np;
            u64 m_lo, m_hi = mW;
            int wrap;
            if (c >= Np) { m_lo = 0; m_hi = Np - 1; wrap = 0; }
            else {
                m_lo = mW + 1 >= c ? mW + 1 - c : mW + 1 + Np - c;
                wrap = m_lo > m_hi;
            }
            for (int part = 0; part <= wrap; part++) {
                u64 a = base + (part == 0 ? m_lo : 0), b = base + ((wrap && part == 0) ? Np - 1 : m_hi);
                u32 i = bst[a >> sh], e = bst[(b >> sh) + 1];
                for (; i < e; i++) {
                    u64 v = val[i];
                    if (v < a) continue;
                    if (v > b) break;
                    u64 m = v - base;
                    u64 y = ys[i];
                    for (u64 z = mW >= m ? mW - m : mW + Np - m; z < c; z += Np) {
                        u128 d = (u128)y * wy + (u128)z * wz;
                        u128 r2 = rem + d;
                        if (reg == 0) {
                            if (r2 >= P) break;
                        } else {
                            if (r2 < P) continue;
                            r2 -= P;
                        }
#ifdef PARANOID
                        if (npf < pcap) pfound[npf++] = (y << 16) | z;
#endif
                        /* same steps as candidate(): count, range, next digit pair, full test */
                        st->hits++;
                        u128 M = M1 + d;
                        if (!inside && (M < pr->lo || M > pr->hi)) continue;
                        if (st->record == 2) {
                            if (pow_sat(B, n - 1) <= M && M < pow_sat(B, n)) record_pair(st, M, B);
                            continue;
                        }
                        u64 dlo = mk_lo + ay_lo[y] + az_lo[z], dhi2 = mk_hi + ay_hi[y] + az_hi[z];
                        if (dlo >= Bk) { dlo -= Bk; dhi2++; }
                        if (dlo >= Bk) { dlo -= Bk; dhi2++; }
                        while (dhi2 >= B) dhi2 -= B;           /* digit k of M */
                        double dq = u128d(r2) * invPtop;       /* digit n-1-k of M: top digit of r2 */
                        u64 dtop = (u64)dq;
                        double fr = dq - (double)dtop;
                        if (fr < 1e-9 || fr > 1 - 1e-9) dtop = (u64)(r2 / Ptop);
                        if (dhi2 != dtop) continue;
                        st->f1++;
                        if (!is_pal_base(M, B, n)) continue;
                        survivor(st, M, B, c);
                    }
                }
            }
        }
#ifdef PARANOID
        if (c <= 48 || x % 29 == 0) {       /* brute force over the (y, z) box */
            u128 M1 = Mblock + (u128)x * wx;
            int nw = 0;
            for (u64 y = 0; y < c; y++)
                for (u64 z = 0; z < c; z++) {
                    u128 M = M1 + (u128)y * wy + (u128)z * wz;
                    u128 qz = M / P;
                    if (qz >= N) continue;
                    if ((u64)(M % N) == T_of((u64)qz, B, k) && nw < pcap) want[nw++] = (y << 16) | z;
                }
            /* both lists are small; compare as sorted sets */
            qsort(pfound, npf, sizeof(u64), u64cmp);
            int same = nw == npf;
            for (int i = 0; same && i < nw; i++) same = want[i] == pfound[i];
            if (!same) {
                char buf[48];
                fprintf(stderr, "PARANOID (y,z) mismatch: n=%d c=%llu B=%llu x=%llu M1=%s fast %d brute %d\n",
                        st->pr->n, (unsigned long long)c, (unsigned long long)B, (unsigned long long)x,
                        u128s(M1, buf), npf, nw);
                abort();
            }
            st->pchk++;
            st->pz += nw;
        }
#endif
        /* advance x */
        rem += wx;
        if (rem >= P) {
            rem -= P;
            if (last) break;                   /* q1 reaches B^k: M >= B^n from here on */
            W += wnext;
            if (W >= N) W -= N;
            int j = k - 1;
            while (j >= 0) {
                E[j]++;
                if (E[j] < B) break;
                E[j] = 0;
                j--;
            }
            RECOMPUTE_NEXT();
        }
        W = W >= Wx ? W - Wx : W + N - Wx;
        mk_lo += ox_lo;
        mk_hi += ox_hi;
        if (mk_lo >= Bk) { mk_lo -= Bk; mk_hi++; }
        if (mk_hi >= B) mk_hi -= B;
        M1 += wx;
    }
#undef RECOMPUTE_NEXT
#ifdef PARANOID
    free(pfound);
    free(want);
#endif
    st->lookups += nlook;
}

static void do_block(wstat_t *st, const cctx_t *cc, u128 Mblock)
{
    const prob_t *pr = st->pr;
    if (!cc->Bcount) return;
    u128 Mmin = Mblock, Mmax = Mblock + cc->span_block;
    if (Mmax < pr->lo || Mmin > pr->hi) return;
    int inside = Mmin >= pr->lo && Mmax <= pr->hi;
    if (Mmin < pr->lo) Mmin = pr->lo;
    if (Mmax > pr->hi) Mmax = pr->hi;
    u64 Blo = iroot(Mmin, pr->Lc) + 1, Bhi = iroot(Mmax, pr->Lc - 1);
    if (Blo < cc->Bbase) Blo = cc->Bbase;
    if (Bhi > cc->Bbase + cc->Bcount - 1) Bhi = cc->Bbase + cc->Bcount - 1;
    for (u64 B = Blo; B <= Bhi; B++) block_B(st, cc, &cc->tab[B - cc->Bbase], Mblock, inside);
    if (Bhi >= Blo) st->steps += (Bhi - Blo + 1) * cc->c;
}

static u64 unit_decode(const prob_t *pr, u64 u, u64 *pref)
{
    u64 lo = 0, hi = pr->cmax - pr->cmin;
    while (lo < hi) {
        u64 mid = (lo + hi + 1) / 2;
        if (pr->uoff[mid] <= u) lo = mid; else hi = mid - 1;
    }
    u64 c = pr->cmin + lo, local = u - pr->uoff[lo];
    for (int i = pr->ud - 1; i >= 1; i--) { pref[i] = local % c; local /= c; }
    pref[0] = 1 + local;
    return c;
}

static void do_unit(wstat_t *st, cctx_t *cc, u64 u)
{
    const prob_t *pr = st->pr;
    u64 pref[MAXN];
    u64 c = unit_decode(pr, u, pref);
    if (cc->c != c) cctx_set(cc, pr, c);
    u128 Mu = 0;
    for (int i = 0; i < pr->ud; i++) Mu += (u128)pref[i] * cc->w[i];
    if (Mu + cc->span_unit < pr->lo || Mu > pr->hi) return;
    int nb = pr->bd - pr->ud;
    if (nb == 0) do_block(st, cc, Mu);
    else {
        u128 wb = cc->w[pr->ud];
        for (u64 d = 0; d < c; d++) {
            u128 Mb = Mu + (u128)d * wb;
            if (Mb > pr->hi) break;
            do_block(st, cc, Mb);
        }
    }
}

/* ------------------------------------------------------------------ */
/* threaded search                                                     */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t g_stop;
static u64 g_end_unit = UINT64_MAX;       /* -E: process units below this index only */
static void on_sig(int s) { (void)s; g_stop = 1; }

static void cctx_free(cctx_t *cc) { free(cc->tab); free(cc->mem); }

static void add_stats(wstat_t *t, const wstat_t *s)
{
    t->steps += s->steps; t->lookups += s->lookups; t->hits += s->hits; t->f1 += s->f1;
    t->st2 += s->st2; t->canon += s->canon; t->csum += s->csum; t->sol += s->sol;
    for (int L = 0; L <= MAXN; L++) t->fail[L] += s->fail[L];
}

static int load_state(const char *path, const prob_t *pr, u64 *frontier, wstat_t *acc)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char lo[64], hi[64];
    int n;
    unsigned long long fr, a[9];
    int ok = fscanf(f, "n=%d lo=%63s hi=%63s frontier=%llu steps=%llu lookups=%llu hits=%llu f1=%llu st2=%llu "
                       "canon=%llu csum=%llu sol=%llu",
                    &n, lo, hi, &fr, &a[0], &a[1], &a[2], &a[3], &a[4], &a[5], &a[6], &a[7]) == 12;
    unsigned long long fl[MAXN + 1] = {0};
    if (ok) {
        for (int L = 3; L <= MAXN; L++) {
            int LL;
            unsigned long long v;
            if (fscanf(f, " fail%d=%llu", &LL, &v) != 2) break;
            if (LL >= 0 && LL <= MAXN) fl[LL] = v;
        }
    }
    fclose(f);
    if (!ok) { fprintf(stderr, "cannot parse state file %s\n", path); exit(1); }
    if (n != pr->n || parse_or_die(lo) != pr->lo || parse_or_die(hi) != pr->hi) {
        fprintf(stderr, "state file %s is for a different search\n", path);
        exit(1);
    }
    *frontier = fr;
    acc->steps = a[0]; acc->lookups = a[1]; acc->hits = a[2]; acc->f1 = a[3]; acc->st2 = a[4];
    acc->canon = a[5]; acc->csum = a[6]; acc->sol = a[7];
    for (int L = 0; L <= MAXN; L++) acc->fail[L] = fl[L];
    return 1;
}

static void save_state(const char *path, const prob_t *pr, u64 frontier, const wstat_t *t)
{
    char tmp[4096], a[48], b[48];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror(tmp); return; }
    fprintf(f, "n=%d lo=%s hi=%s frontier=%llu steps=%llu lookups=%llu hits=%llu f1=%llu st2=%llu canon=%llu "
               "csum=%llu sol=%llu",
            pr->n, u128s(pr->lo, a), u128s(pr->hi, b), (unsigned long long)frontier, (unsigned long long)t->steps,
            (unsigned long long)t->lookups, (unsigned long long)t->hits, (unsigned long long)t->f1,
            (unsigned long long)t->st2, (unsigned long long)t->canon, (unsigned long long)t->csum,
            (unsigned long long)t->sol);
    for (int L = pr->n - 2; L >= 3; L--) fprintf(f, " fail%d=%llu", L, (unsigned long long)t->fail[L]);
    fprintf(f, "\n");
    fclose(f);
    rename(tmp, path);
}

static void print_summary(FILE *o, const prob_t *pr, const wstat_t *t, int complete, double el)
{
    char a[48], b[48];
    fprintf(o, "# n=%d [%s, %s] %s: %.1fs, steps %llu, lookups %llu, hits %llu, f1 %llu, (n-1)-palindromes %llu, "
               "canonical %llu, checksum %016llx, solutions %llu\n",
            pr->n, u128s(pr->lo, a), u128s(pr->hi, b), complete ? "complete" : "INTERRUPTED", el,
            (unsigned long long)t->steps, (unsigned long long)t->lookups, (unsigned long long)t->hits,
            (unsigned long long)t->f1, (unsigned long long)t->st2, (unsigned long long)t->canon,
            (unsigned long long)t->csum, (unsigned long long)t->sol);
    fprintf(o, "# failed at length:");
    for (int L = pr->n - 2; L >= 3; L--) fprintf(o, " L=%d:%llu", L, (unsigned long long)t->fail[L]);
    fprintf(o, "\n");
    fflush(o);
}

#ifndef A171775_2D_NO_DRIVER
typedef struct {
    prob_t *pr;
    _Atomic u64 next;
    unsigned char *done;
    u64 frontier;
    pthread_mutex_t lock;
    wstat_t *st;
    int nthreads;
} search_t;

typedef struct { search_t *S; int id; } warg_t;

static void *worker(void *arg)
{
    warg_t *a = arg;
    search_t *S = a->S;
    wstat_t *st = &S->st[a->id];
    cctx_t cc;
    memset(&cc, 0, sizeof cc);
    while (!g_stop) {
        u64 u = atomic_fetch_add(&S->next, 1);
        if (u >= S->pr->nunits) break;
        do_unit(st, &cc, u);
        pthread_mutex_lock(&S->lock);
        S->done[u] = 1;
        while (S->frontier < S->pr->nunits && S->done[S->frontier]) S->frontier++;
        pthread_mutex_unlock(&S->lock);
    }
    cctx_free(&cc);
    return NULL;
}

static void sum_stats(search_t *S, wstat_t *tot)
{
    memset(tot, 0, sizeof *tot);
    for (int i = 0; i < S->nthreads; i++) add_stats(tot, &S->st[i]);
}

/* rough work of a unit, for progress only */
static double unit_weight(const prob_t *pr, u64 c, u64 d0)
{
    int Le = pr->Le, Lc = pr->Lc;
    u128 Mmin = (u128)d0 * pow_sat(c, Le - 1), Mmax = (u128)(d0 + 1) * pow_sat(c, Le - 1);
    if (Mmax < pr->lo || Mmin > pr->hi) return 0;
    if (Mmin < pr->lo) Mmin = pr->lo;
    if (Mmax > pr->hi) Mmax = pr->hi;
    double bl = pow((double)Mmin, 1.0 / Lc), bh = pow((double)Mmax, 1.0 / (Lc - 1));
    double nB = bh - bl + 1;
    if (nB < 0) nB = 0;
    double frac = (double)(Mmax - Mmin) / (double)pow_sat(c, Le - 1);
    return nB * pow((double)c, pr->De - 3) * frac;
}

/* returns 0 if complete; result (optional) receives the totals */
static int run_search(int n, u128 lo, u128 hi, int nthreads, const char *state, double interval, int verbose,
                      int quiet, wstat_t *result)
{
    char a[48], b[48];
    prob_t pr;
    prob_init(&pr, n, lo, hi);
    if (g_end_unit < pr.nunits) pr.nunits = g_end_unit;
    search_t S;
    memset(&S, 0, sizeof S);
    S.pr = &pr;
    S.nthreads = nthreads;
    S.done = calloc(pr.nunits + 1, 1);
    S.st = calloc(nthreads, sizeof(wstat_t));
    pthread_mutex_init(&S.lock, NULL);
    wstat_t base;
    memset(&base, 0, sizeof base);
    u64 start = 0;
    if (state && load_state(state, &pr, &start, &base) && !quiet)
        fprintf(stderr, "resuming at unit %llu of %llu\n", (unsigned long long)start, (unsigned long long)pr.nunits);
    for (u64 u = 0; u < start; u++) S.done[u] = 1;
    S.frontier = start;
    atomic_store(&S.next, start);
    for (int i = 0; i < nthreads; i++) { S.st[i].pr = &pr; S.st[i].verbose = verbose; }

    double *cw = malloc((pr.nunits + 1) * sizeof(double));
    cw[0] = 0;
    for (u64 u = 0; u < pr.nunits; u++) {
        u64 pref[MAXN];
        u64 c = unit_decode(&pr, u, pref);
        double wgt = unit_weight(&pr, c, pref[0]);
        if (pr.ud > 1) wgt /= (double)c;
        cw[u + 1] = cw[u] + wgt;
    }
    double wtot = cw[pr.nunits] > 0 ? cw[pr.nunits] : 1;
    if (!quiet)
        fprintf(stderr, "search2 n=%d [%s, %s]: enumerate %d-digit palindromes in bases c in [%llu, %llu], %llu units, "
                        "%d threads\n",
                n, u128s(lo, a), u128s(hi, b), pr.Le, (unsigned long long)pr.cmin, (unsigned long long)pr.cmax,
                (unsigned long long)pr.nunits, nthreads);
    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);
    pthread_t *th = malloc(nthreads * sizeof(pthread_t));
    warg_t *wa = malloc(nthreads * sizeof(warg_t));
    double t0 = now_sec(), tlast = t0, tsave = t0, w0 = cw[start];
    for (int i = 0; i < nthreads; i++) {
        wa[i].S = &S;
        wa[i].id = i;
        pthread_create(&th[i], NULL, worker, &wa[i]);
    }
    for (;;) {
        usleep(100000);
        pthread_mutex_lock(&S.lock);
        u64 fr = S.frontier;
        pthread_mutex_unlock(&S.lock);
        u64 nx = atomic_load(&S.next);
        int finished = fr >= pr.nunits;
        double t = now_sec();
        if (!quiet && (t - tlast >= interval || finished || g_stop)) {
            wstat_t tot;
            sum_stats(&S, &tot);
            double el = t - t0, wd = cw[fr < pr.nunits ? fr : pr.nunits];
            double eta = (wd - w0) > 0 ? (wtot - wd) * el / (wd - w0) : 0;
            u64 pref[MAXN];
            u64 cc_ = unit_decode(&pr, fr < pr.nunits ? fr : pr.nunits - 1, pref);
            fprintf(stderr,
                    "[%7.0fs] frontier %llu/%llu (c=%llu) %.3f%%  %.3e steps/s %.3e lookups/s  st2 %llu canon %llu sol %llu  "
                    "ETA %.0fs\n",
                    el, (unsigned long long)fr, (unsigned long long)pr.nunits, (unsigned long long)cc_, 100 * wd / wtot,
                    tot.steps / (el > 0 ? el : 1), tot.lookups / (el > 0 ? el : 1),
                    (unsigned long long)(base.st2 + tot.st2), (unsigned long long)(base.canon + tot.canon),
                    (unsigned long long)(base.sol + tot.sol), eta);
            tlast = t;
        }
        if (state && (t - tsave >= 60 || finished || g_stop)) {
            wstat_t tot;
            sum_stats(&S, &tot);
            add_stats(&tot, &base);
            save_state(state, &pr, fr, &tot);
            tsave = t;
        }
        if (finished || (g_stop && nx >= pr.nunits) || g_stop) break;
    }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    wstat_t tot;
    sum_stats(&S, &tot);
    add_stats(&tot, &base);
    if (state) save_state(state, &pr, S.frontier, &tot);
    int complete = S.frontier >= pr.nunits;
    if (!quiet) print_summary(stdout, &pr, &tot, complete, now_sec() - t0);
    if (result) *result = tot;
    free(th); free(wa); free(cw); free(S.done); free(S.st); free(pr.uoff);
    return complete ? 0 : 3;
}

/* ------------------------------------------------------------------ */
/* selftest                                                            */
/* ------------------------------------------------------------------ */

static u64 rng_state = 0x2545F4914F6CDD1DULL;
static u64 rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

static int pair_cmp(const void *a, const void *b)
{
    const pair_t *x = a, *y = b;
    if (x->M != y->M) return x->M < y->M ? -1 : 1;
    if (x->B != y->B) return x->B < y->B ? -1 : 1;
    return 0;
}

static int ends_mirror(u128 M, u64 B, int L, int k)
{
    u64 d[130];
    int len = 0;
    u128 x = M;
    while (x) { d[len++] = (u64)(x % B); x /= B; }
    if (len != L) return 0;
    for (int i = 0; i < k; i++)
        if (d[i] != d[L - 1 - i]) return 0;
    return 1;
}

/* all (M, B) for M in the block: mode 1 = n-digit palindrome in B, mode 2 = the k(c,B) ends mirror */
static void naive_block(const prob_t *pr, const cctx_t *cc, u128 Mblock, int mode, pair_t **out, size_t *np)
{
    size_t cap = 1024, cnt = 0;
    pair_t *p = malloc(cap * sizeof(pair_t));
    u64 c = cc->c;
    for (u64 x = 0; x < c; x++)
        for (u64 y = 0; y < c; y++)
            for (u64 z = 0; z < c; z++) {
                u128 M = Mblock + (u128)x * cc->wx + (u128)y * cc->wy + (u128)z * cc->wz;
                if (M < pr->lo || M > pr->hi) continue;
                u64 Blo = iroot(M, pr->Lc) + 1, Bhi = iroot(M, pr->Lc - 1);
                if (Blo < 2) Blo = 2;
                for (u64 B = Blo; B <= Bhi; B++) {
                    if (B < cc->Bbase || B >= cc->Bbase + cc->Bcount) {
                        fprintf(stderr, "naive: B=%llu outside the table range\n", (unsigned long long)B);
                        exit(1);
                    }
                    int ok = mode == 1 ? is_pal_base(M, B, pr->Lc)
                                       : ends_mirror(M, B, pr->Lc, (int)cc->tab[B - cc->Bbase].k);
                    if (!ok) continue;
                    if (cnt == cap) { cap *= 2; p = realloc(p, cap * sizeof(pair_t)); }
                    p[cnt].M = M;
                    p[cnt].B = B;
                    cnt++;
                }
            }
    *out = p;
    *np = cnt;
}

static int compare_block(int n, u64 c, u64 e0, u64 e1, int mode, size_t *npairs)
{
    prob_t pr;
    prob_init(&pr, n, 1, U128_MAX >> 2);
    cctx_t cc;
    memset(&cc, 0, sizeof cc);
    cctx_set(&cc, &pr, c);
    u128 Mb = (u128)e0 * cc.w[0] + (pr.bd >= 2 ? (u128)e1 * cc.w[1] : 0);
    wstat_t st;
    memset(&st, 0, sizeof st);
    st.pr = &pr;
    st.record = mode;
    do_block(&st, &cc, Mb);
    pair_t *np;
    size_t nn;
    naive_block(&pr, &cc, Mb, mode, &np, &nn);
    qsort(st.pairs, st.np, sizeof(pair_t), pair_cmp);
    qsort(np, nn, sizeof(pair_t), pair_cmp);
    int ok = st.np == nn;
    for (size_t i = 0; ok && i < nn; i++) ok = st.pairs[i].M == np[i].M && st.pairs[i].B == np[i].B;
    if (!ok) {
        char buf[48];
        printf("  MISMATCH n=%d c=%llu e0=%llu e1=%llu mode %d (Mblock=%s): fast %zu, naive %zu\n", n,
               (unsigned long long)c, (unsigned long long)e0, (unsigned long long)e1, mode, u128s(Mb, buf), st.np, nn);
    }
    *npairs += nn;
    free(st.pairs);
    free(np);
    cctx_free(&cc);
    free(pr.uoff);
    return ok;
}

static int cmd_selftest(int full, int nthreads)
{
    int bad = 0;
    printf("1. fast lookup vs naive enumeration on random blocks\n");
    struct { int n; u64 cmin, cmax; int reps; } cfg[] = {
        {8, 4, 40, 60}, {9, 4, 36, 60}, {10, 4, 26, 50}, {11, 4, 22, 40}, {12, 4, 12, 20},
    };
    for (size_t ci = 0; ci < sizeof cfg / sizeof cfg[0]; ci++) {
        int n = cfg[ci].n, okc = 0, tot = 0;
        size_t p1 = 0, p2 = 0;
        for (int r = 0; r < cfg[ci].reps; r++) {
            u64 c = cfg[ci].cmin + rnd() % (cfg[ci].cmax - cfg[ci].cmin + 1);
            u64 e0 = 1 + rnd() % (c - 1), e1 = rnd() % c;
            if (r == 0) { e0 = 1; e1 = 0; }
            if (r == 1) { e0 = c - 1; e1 = c - 1; }
            int ok = compare_block(n, c, e0, e1, 1, &p1) & compare_block(n, c, e0, e1, 2, &p2);
            okc += ok;
            tot++;
        }
        printf("  n=%d c in [%llu,%llu]: %d/%d blocks agree (%zu n-digit palindrome pairs, %zu k-digit matches)\n", n,
               (unsigned long long)cfg[ci].cmin, (unsigned long long)cfg[ci].cmax, okc, tot, p1, p2);
        fflush(stdout);
        bad += okc != tot;
    }
    printf("2. full searches vs a171775.c (1D) statistics\n");
    struct { int n; const char *hi; u64 st2, canon, csum; u64 fails[3]; } ref[] = {
        {8, "2^42", 4725, 4724, 0x13747afea2768e7aULL, {4723, 0, 0}},
        {9, "2^56", 19140, 19140, 0x380d9946bb59f02eULL, {19138, 1, 0}},
        {10, "2^72", 135382, 135382, 0xe9eed6127c51a50fULL, {135381, 0, 0}},
    };
    int nref = full ? 3 : 2;
    for (int i = 0; i < nref; i++) {
        FILE *mem = tmpfile();
        FILE *save = sol_out;
        sol_out = mem;
        wstat_t t;
        double t0 = now_sec();
        run_search(ref[i].n, 1, parse_or_die(ref[i].hi), nthreads, NULL, 1e9, 0, 1, &t);
        sol_out = save;
        rewind(mem);
        char line[1024];
        int nsol = 0;
        u128 best = 0;
        while (fgets(line, sizeof line, mem)) {
            char *p = strstr(line, "M=");
            if (!p) continue;
            char num[64];
            sscanf(p + 2, "%63[0-9]", num);
            u128 M = parse_or_die(num);
            if (!best || M < best) best = M;
            nsol++;
        }
        fclose(mem);
        int n = ref[i].n;
        int ok = nsol == 1 && best == parse_or_die(ref[i].hi) && t.st2 == ref[i].st2 && t.canon == ref[i].canon &&
                 t.csum == ref[i].csum && t.fail[n - 2] == ref[i].fails[0] && t.fail[n - 3] == ref[i].fails[1] &&
                 t.fail[n - 4] == ref[i].fails[2];
        char buf[48];
        printf("  n=%d [1, %s]: %.2fs, %.3e steps, %.3e lookups, hits %llu, f1 %llu, st2 %llu, canonical %llu, "
               "checksum %016llx, fails %llu/%llu/%llu, %d solution(s), smallest %s  %s\n",
               n, ref[i].hi, now_sec() - t0, (double)t.steps, (double)t.lookups, (unsigned long long)t.hits,
               (unsigned long long)t.f1, (unsigned long long)t.st2, (unsigned long long)t.canon,
               (unsigned long long)t.csum, (unsigned long long)t.fail[n - 2], (unsigned long long)t.fail[n - 3],
               (unsigned long long)t.fail[n - 4], nsol, best ? u128s(best, buf) : "none", ok ? "ok" : "WRONG");
        fflush(stdout);
        bad += !ok;
    }
    printf(bad ? "SELFTEST FAILED (%d)\n" : "selftest passed\n", bad);
    return bad ? 1 : 0;
}

#ifdef PARANOID
static int cmd_paranoid(int reps)
{
    struct { int n; u64 cmin, cmax; u128 hi; } cfg[] = {
        {8, 4, 64, (u128)1 << 60}, {9, 4, 64, (u128)1 << 62}, {10, 4, 64, (u128)1 << 72},
        {11, 4, 64, (u128)1 << 90}, {12, 4, 40, (u128)1 << 100},
        {10, 200, 512, (u128)1 << 72}, {11, 400, 1024, (u128)1 << 90},
    };
    FILE *devnull = fopen("/dev/null", "w");
    FILE *save = sol_out;
    sol_out = devnull;
    for (size_t ci = 0; ci < sizeof cfg / sizeof cfg[0]; ci++) {
        int n = cfg[ci].n;
        u64 blocks = 0, chk = 0, zh = 0, steps = 0;
        double t0 = now_sec();
        for (int r = 0; r < reps && now_sec() - t0 < 60; r++) {
            u64 c = cfg[ci].cmin + rnd() % (cfg[ci].cmax - cfg[ci].cmin + 1);
            prob_t pr;
            prob_init(&pr, n, 1, cfg[ci].hi);
            cctx_t cc;
            memset(&cc, 0, sizeof cc);
            cctx_set(&cc, &pr, c);
            u128 dm = cfg[ci].hi / cc.w[0];
            u64 dmax = dm < c - 1 ? (u64)dm : c - 1;
            if (dmax >= 1 && cc.Bcount) {
                u64 e0 = 1 + rnd() % dmax, e1 = rnd() % c;
                u128 Mb = (u128)e0 * cc.w[0] + (pr.bd >= 2 ? (u128)e1 * cc.w[1] : 0);
                wstat_t st;
                memset(&st, 0, sizeof st);
                st.pr = &pr;
                do_block(&st, &cc, Mb);
                blocks++;
                chk += st.pchk;
                zh += st.pz;
                steps += st.steps;
            }
            cctx_free(&cc);
            free(pr.uoff);
        }
        printf("  paranoid n=%d c in [%llu,%llu]: %llu blocks, %llu x-steps state-checked, %llu brute-force (y,z) boxes, "
               "%llu matches: ok\n",
               n, (unsigned long long)cfg[ci].cmin, (unsigned long long)cfg[ci].cmax, (unsigned long long)blocks,
               (unsigned long long)steps, (unsigned long long)chk, (unsigned long long)zh);
        fflush(stdout);
    }
    sol_out = save;
    return 0;
}
#endif

/* ------------------------------------------------------------------ */

static void usage(void)
{
    fprintf(stderr,
            "usage: a171775_2d selftest [full] [-t T]\n"
            "       a171775_2d search n LO HI [-t T] [-S state] [-i sec] [-v]\n"
#ifdef PARANOID
            "       a171775_2d paranoid [reps]\n"
#endif
    );
    exit(1);
}

#ifndef A171775_2D_LIB
int main(int argc, char **argv)
{
    sol_out = stdout;
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    const char *state = NULL;
    double interval = 10;
    int verbose = 0;
    char *pos[8];
    int npos = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-t") && i + 1 < argc) nthreads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-S") && i + 1 < argc) state = argv[++i];
        else if (!strcmp(argv[i], "-i") && i + 1 < argc) interval = atof(argv[++i]);
        else if (!strcmp(argv[i], "-v")) verbose = 1;
        else if (!strcmp(argv[i], "-E") && i + 1 < argc) g_end_unit = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-u") && i + 1 < argc) {   /* print the first unit of base c */
            int n = atoi(argv[++i]);
            u64 c = strtoull(argv[++i], NULL, 10);
            prob_t pr;
            prob_init(&pr, n, 1, parse_or_die(argv[++i]));
            printf("%llu\n", (unsigned long long)pr.uoff[c - pr.cmin]);
            return 0;
        }
        else if (npos < 8) pos[npos++] = argv[i];
        else usage();
    }
    if (nthreads < 1) nthreads = 1;
    if (npos < 1) usage();
    if (!strcmp(pos[0], "selftest")) return cmd_selftest(npos > 1 && !strcmp(pos[1], "full"), nthreads);
#ifdef PARANOID
    if (!strcmp(pos[0], "paranoid")) return cmd_paranoid(npos > 1 ? atoi(pos[1]) : 300);
#endif
    if (!strcmp(pos[0], "search") && npos == 4) {
        int n = atoi(pos[1]);
        if (n < 8 || n > 20) { fprintf(stderr, "n must be 8..20 (use a171775 for smaller n)\n"); return 1; }
        return run_search(n, parse_or_die(pos[2]), parse_or_die(pos[3]), nthreads, state, interval, verbose, 0, NULL);
    }
    usage();
    return 1;
}
#endif /* A171775_2D_LIB */
#endif /* A171775_2D_NO_DRIVER */
