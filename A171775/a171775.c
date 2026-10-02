/*
 * a171775.c
 *
 * Compute and extend OEIS A171775:
 *
 *   "a(n) = smallest number M such that there exist bases b_2, b_3, ..., b_n
 *    with the property that M written in base b_k is a k-digit palindrome
 *    for all k = 2..n."
 *
 * Known terms (n = 1..9):
 *   1, 3, 5, 52, 130, 1885, 1073741824 = 2^30, 4398046511104 = 2^42,
 *   72057594037927936 = 2^56 (a(9): Max Alekseyev, Jun 2026).
 * Bound: a(n) <= 2^((n-1)(n-2)) (James G. Merickel); conjecture (Alekseyev):
 * equality for n >= 7.  So a(10) <= 2^72 = 4722366482869645213696.
 *
 * Why 2^((n-1)(n-2)) works: in base b = 2^e - 1, 2^r (b+1)^(L-1) has the
 * digits 2^r C(L-1, j), a palindrome of L digits whenever 2^r C(L-1, j) < b.
 * With N = (n-1)(n-2) = e(L-1) + r every length L = 2..n finds such (e, r).
 *
 * Method
 * ------
 * L = 2 is free (M = "11" in base M-1 for M >= 3).  Every solution is an
 * n-digit palindrome in some base B and an (n-1)-digit palindrome in some
 * base c.  The n-digit palindromes are the smallest of the sets (about
 * 3.7e12 below 2^72 for n = 10), so they are enumerated, and the (n-1)-digit
 * condition is imposed by solving a congruence for the innermost digit
 * instead of trying every value of it.
 *
 * n-digit palindrome in base B: h = n/2 mirrored pairs, D = ceil(n/2) free
 * digits d_0..d_{D-1} (d_0 >= 1) with weights w_i = B^(n-1-i) + B^i
 * (i < h) and, for odd n, w_h = B^h for the middle digit.  Write the three
 * innermost free digits as x = d_{D-3}, y = d_{D-2}, z = d_{D-1} and the
 * rest as the "block" (d_0..d_{D-4}), so M = Mblock + x w_x + y w_y + z s.
 * For every block and every base c that can give M exactly n-1 digits
 * (c^(n-2) <= M < c^(n-1); always c > B), the tool loops over (x, y) and
 * solves for z.
 *
 * Let k = n-2-h and P = c^(h+1), N = c^k.  For fixed (x, y), z moves M by
 * less than (B-1) s < B^(h+1) < P, so q = floor(M / P) takes at most two
 * values, q and q+1.  q holds the k leading base-c digits (E_0, ..., E_{k-1})
 * of M and, M being a palindrome, the k trailing digits must be the same in
 * reverse: M = T(q) (mod N) with T = E_0 + E_1 c + ... + E_{k-1} c^(k-1).
 * That is  z s = T - M0 (mod N),  M0 = Mblock + x w_x + y w_y.
 *
 * With g = gcd(s, N), N' = N/g and s'' = s/g + j N' coprime to N, put
 * u = s''^(-1) mod N and w = (T - M0) u mod N.  Then a solution exists iff
 * g | w, and the smallest one is z0 = w / g (the others are z0 + i N').
 * So a single comparison, w < g B, is the whole test per (x, y); w, the
 * remainder M0 mod P and the digits of q are all updated incrementally with
 * additions only (y -> y+1 adds w_y to M0; the wrap y = B-1 -> 0, x -> x+1
 * adds w_x - (B-1) w_y; q grows by at most one per step, which shifts T by
 * c^(k-1), or by c^(m+1) + c^m (mod c^k) when the carry stops at digit m).  The rare hits (probability
 * about B/c^k) are checked in full: all n-1 digits of M in base c.
 *
 * Survivors (M an n-digit and an (n-1)-digit palindrome) are tested for
 * L = n-2, ..., 3 by scanning the bases b with b^(L-1) <= M < b^L: last
 * digit M mod b, a floating-point filter on the leading digit, exact check.
 * Duplicates (M a palindrome of the same length in two bases) are removed by
 * keeping only the pair whose B and c are the smallest possible bases.
 *
 * Work is cut into units (B, d_0) (plus d_1 for n >= 11) handed to pthreads
 * through an atomic counter; the lowest unfinished unit is the checkpoint
 * frontier (-S state file).  Everything is exact integer arithmetic
 * (unsigned __int128 for M); n <= 6 uses a plain enumeration.
 *
 * Usage
 *   a171775 selftest [full]        known terms, fast filter vs naive, a(7), a(8) [, a(9)]
 *   a171775 search n LO HI [-t T] [-S state] [-i sec] [-v]
 *   a171775 check M [n]            bases for each length 2..n
 *   a171775 naive n LO HI          reference enumeration (small ranges)
 * Numbers accept 2^72, 1e21, 2^72-1, ...
 */
#define _GNU_SOURCE
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>
#include <stdatomic.h>
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

#define MAXN 24
#define MAXK 8
#define U128_MAX (~(u128)0)

/* ------------------------------------------------------------------ */
/* integer helpers                                                     */
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

/* floor(x^(1/k)), k >= 2 */
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
    if (*s == '^') {
        s++;
        char *e;
        unsigned long ex = strtoul(s, &e, 10);
        if (e == s) return 0;
        s = e;
        if (v > UINT64_MAX) return 0;
        v = pow_sat((u64)v, (int)ex);
    } else if (*s == 'e' || *s == 'E') {
        s++;
        char *e;
        unsigned long ex = strtoul(s, &e, 10);
        if (e == s) return 0;
        s = e;
        v *= pow_sat(10, (int)ex);
    }
    *ps = s;
    *out = v;
    return 1;
}

static int parse_u128(const char *s, u128 *out)
{
    u128 acc, t;
    if (!parse_term(&s, &acc)) return 0;
    while (*s == '+' || *s == '-') {
        char op = *s++;
        if (!parse_term(&s, &t)) return 0;
        acc = op == '+' ? acc + t : acc - t;
    }
    if (*s) return 0;
    *out = acc;
    return 1;
}

static u128 parse_or_die(const char *s)
{
    u128 v;
    if (!parse_u128(s, &v)) { fprintf(stderr, "bad number: %s\n", s); exit(1); }
    return v;
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

/* ------------------------------------------------------------------ */
/* palindrome tests                                                    */
/* ------------------------------------------------------------------ */

/* M written in base b has exactly L digits and is a palindrome */
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

/* smallest base b >= 2 in which M is an L-digit palindrome, 0 if none */
static u64 find_pal_base(u128 M, int L)
{
    if (L < 2 || M == 0) return 0;
    if (L == 2) {                       /* "11" in base M-1; only existence matters */
        if (M < 3) return 0;
        return (M - 1 > UINT64_MAX) ? UINT64_MAX : (u64)(M - 1);
    }
    u64 blo = iroot(M, L) + 1, bhi = iroot(M, L - 1);
    if (blo < 2) blo = 2;
    double Md = (double)M;
    for (u64 b = blo; b <= bhi; b++) {
        u64 r = mod128_64(M, b);
        if (r == 0) continue;
        double Pd = pow((double)b, (double)(L - 1));
        double la = Md / Pd;
        if (la < (double)r - 1.5 || la > (double)r + 1.5) continue;
        u128 P = pow_sat(b, L - 1);
        if (M / P != r) continue;
        if (is_pal_base(M, b, L)) return b;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* problem setup                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    int n, h, D, k, bd, ud;     /* length, pairs, free digits, base-c digits used, block digits, unit digits */
    u128 lo, hi;
    u64 Bmin, Bmax;
    u64 *uoff;                  /* uoff[B-Bmin] = first unit of base B */
    u64 nunits;
} prob_t;

typedef struct {
    u64 c, N, P, g, Np, inv, limit, Wy, NWy, Wx, NWx, Rcross;
    u64 Wc[MAXK];
} cconst_t;

typedef struct {
    u64 B;
    u128 w[MAXN];
    u64 s, wy, dx;              /* innermost weight, y weight, x-wrap delta */
    u128 wx;
    u128 span_block, span_unit;
    u64 cbase, ccount;          /* table of c constants */
    cconst_t *tab;
    size_t tabcap;
} bctx_t;

typedef struct { u128 M; u64 c; } pair_t;

typedef struct {
    const prob_t *pr;
    u64 tuples, hits, st2, canon;
    u64 csum;                   /* order-independent checksum of the canonical survivors */
    u64 fail[MAXN + 1];         /* fail[L]: canonical survivors that failed at length L */
    u64 sol;
    int record;                 /* test hook: record (M, c) pairs with M an (n-1)-palindrome in base c */
    pair_t *pairs;
    size_t np, pcap;
    int verbose;
} wstat_t;

static pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
static FILE *sol_out;           /* solutions go here (stdout) */

static void prob_init(prob_t *pr, int n, u128 lo, u128 hi)
{
    memset(pr, 0, sizeof *pr);
    pr->n = n;
    pr->h = n / 2;
    pr->D = (n + 1) / 2;
    pr->k = n - 2 - pr->h;
    pr->bd = pr->D - 3;
    pr->ud = pr->bd - 1 >= 1 ? pr->bd - 1 : 1;
    if (lo < 1) lo = 1;
    pr->lo = lo;
    pr->hi = hi;
    /* B with B^(n-1) <= hi and B^n - 1 >= lo */
    u64 Bmax = iroot(hi, n - 1);
    u64 Bmin = 2;
    while (Bmin <= Bmax && pow_sat(Bmin, n) - 1 < lo) Bmin++;
    if (Bmin < 2) Bmin = 2;
    pr->Bmin = Bmin;
    pr->Bmax = Bmax;
    u64 nb = Bmax >= Bmin ? Bmax - Bmin + 1 : 0;
    pr->uoff = calloc(nb + 1, sizeof(u64));
    u64 tot = 0;
    for (u64 i = 0; i < nb; i++) {
        u64 B = Bmin + i;
        pr->uoff[i] = tot;
        u128 cnt = (u128)(B - 1) * pow_sat(B, pr->ud - 1);
        tot += (u64)cnt;
    }
    pr->uoff[nb] = tot;
    pr->nunits = tot;
}

static void bctx_set(bctx_t *bc, const prob_t *pr, u64 B)
{
    int n = pr->n, h = pr->h, D = pr->D;
    bc->B = B;
    for (int i = 0; i < h; i++) bc->w[i] = pow_sat(B, n - 1 - i) + pow_sat(B, i);
    if (n & 1) bc->w[h] = pow_sat(B, h);
    bc->s = (u64)bc->w[D - 1];
    bc->wy = (u64)bc->w[D - 2];
    bc->wx = bc->w[D - 3];
    bc->dx = (u64)(bc->w[D - 3] - (u128)(B - 1) * bc->w[D - 2]);
    bc->span_block = (u128)(B - 1) * (bc->w[D - 3] + bc->w[D - 2] + bc->w[D - 1]);
    u128 su = 0;
    for (int i = pr->ud; i < D; i++) su += bc->w[i];
    bc->span_unit = (u128)(B - 1) * su;

    /* c range for this B */
    u128 lo = pr->lo, hi = pr->hi;
    u128 Mlo = pow_sat(B, n - 1), Mhi = pow_sat(B, n) - 1;
    if (Mlo < lo) Mlo = lo;
    if (Mhi > hi) Mhi = hi;
    u64 cmin = iroot(Mlo, n - 1) + 1, cmax = iroot(Mhi, n - 2);
    if (cmax < cmin) cmax = cmin;
    bc->cbase = cmin;
    bc->ccount = cmax - cmin + 1;
    if (bc->ccount > bc->tabcap) {
        free(bc->tab);
        bc->tabcap = bc->ccount;
        bc->tab = malloc(bc->tabcap * sizeof(cconst_t));
    }
    int k = pr->k;
    for (u64 i = 0; i < bc->ccount; i++) {
        u64 c = cmin + i;
        cconst_t *K = &bc->tab[i];
        u128 P128 = pow_sat(c, h + 1), N128 = pow_sat(c, k);
        if (P128 >= ((u128)1 << 62) || N128 >= ((u128)1 << 42)) {
            fprintf(stderr, "c = %llu too large for 64-bit state (n = %d)\n", (unsigned long long)c, n);
            exit(1);
        }
        u64 P = (u64)P128, N = (u64)N128;
        if ((u128)(B - 1) * bc->s >= P || bc->wy >= P || bc->dx >= P) {
            fprintf(stderr, "internal: weight bound violated (B=%llu c=%llu)\n",
                    (unsigned long long)B, (unsigned long long)c);
            exit(1);
        }
        K->c = c;
        K->N = N;
        K->P = P;
        u64 sN = bc->s % N;
        u64 g = gcd64(sN, N);          /* gcd(0, N) = N */
        u64 Np = N / g;
        u64 sp = Np > 1 ? (sN / g) % Np : 0;
        u64 spp = sp;
        while (gcd64(spp % N, N) != 1) spp += Np;
        K->g = g;
        K->Np = Np;
        K->inv = modinv64(spp % N, N);
        K->limit = g * B;
        K->Wy = mulmod64(bc->wy % N, K->inv, N);
        K->NWy = N - K->Wy;
        K->Wx = mulmod64(bc->dx % N, K->inv, N);
        K->NWx = N - K->Wx;
        u64 cj = 1;
        for (int j = 0; j < k; j++) { K->Wc[j] = mulmod64(cj % N, K->inv, N); cj *= c; }
        K->Rcross = P - (B - 1) * bc->s;
    }
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

static inline u64 surv_hash(u128 M)
{
    u64 a = (u64)M * 0x9E3779B97F4A7C15ULL, b = (u64)(M >> 64) * 0xC2B2AE3D27D4EB4FULL;
    a ^= a >> 29;
    return a * 0xBF58476D1CE4E5B9ULL + b;
}

static void candidate(wstat_t *st, u128 M, u64 B, u64 c)
{
    const prob_t *pr = st->pr;
    st->hits++;
    if (M < pr->lo || M > pr->hi) return;
    int n = pr->n;
    if (st->record == 2) {              /* test hook: every k-digit match with n-1 digits */
        if (pow_sat(c, n - 2) <= M && M < pow_sat(c, n - 1)) {
            if (st->np == st->pcap) {
                st->pcap = st->pcap ? 2 * st->pcap : 1024;
                st->pairs = realloc(st->pairs, st->pcap * sizeof(pair_t));
            }
            st->pairs[st->np].M = M;
            st->pairs[st->np].c = c;
            st->np++;
        }
        return;
    }
    if (!is_pal_base(M, c, n - 1)) return;
    st->st2++;
    if (st->record) {
        if (st->np == st->pcap) {
            st->pcap = st->pcap ? 2 * st->pcap : 1024;
            st->pairs = realloc(st->pairs, st->pcap * sizeof(pair_t));
        }
        st->pairs[st->np].M = M;
        st->pairs[st->np].c = c;
        st->np++;
        return;
    }
    /* canonical: c and B are the smallest bases of their lengths */
    if (find_pal_base(M, n - 1) != c) return;
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
/* the fast filter: one block, one base c                              */
/* ------------------------------------------------------------------ */

#ifdef PARANOID
/* debug build: recompute the incremental state from scratch at every step */
static u64 T_of(u64 q, u64 c, int k)
{
    u64 E[MAXK], T = 0;
    for (int j = k - 1; j >= 0; j--) { E[j] = q % c; q /= c; }
    for (int j = k - 1; j >= 0; j--) T = T * c + E[j];
    return T;
}

static void paranoid_state(const prob_t *pr, const bctx_t *bc, const cconst_t *K, u128 M0, u64 w, u64 rem,
                           const u64 *E, u64 elast, u64 wnext, int last)
{
    const int k = pr->k;
    const u64 c = K->c, N = K->N, P = K->P;
    char buf[48];
    u128 q = M0 / P;
    int bad = 0;
    if ((u64)(M0 % P) != rem) bad |= 1;
    if (q >= N) bad |= 2;
    else {
        u64 qq = (u64)q, Et[MAXK];
        for (int j = k - 1; j >= 0; j--) { Et[j] = qq % c; qq /= c; }
        if (Et[k - 1] != elast) bad |= 4;
        for (int j = 0; j < k - 1; j++)
            if (Et[j] != E[j]) bad |= 8;
        u64 T = T_of((u64)q, c, k), MN = (u64)(M0 % N);
        u64 wt = mulmod64(T >= MN ? T - MN : T + N - MN, K->inv, N);
        if (wt != w) bad |= 16;
        if ((u64)q + 1 < N) {
            u64 T1 = T_of((u64)q + 1, c, k);
            u64 wn = mulmod64(T1 >= T ? T1 - T : T1 + N - T, K->inv, N);
            if (wn != wnext) bad |= 32;
            if (last) bad |= 64;
        } else if (!last) bad |= 64;
    }
    if (bad) {
        fprintf(stderr, "PARANOID state mismatch %#x: n=%d B=%llu c=%llu M0=%s\n", bad, pr->n,
                (unsigned long long)bc->B, (unsigned long long)c, u128s(M0, buf));
        abort();
    }
}

static int u64cmp(const void *a, const void *b)
{
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}
#endif

static void block_c(wstat_t *st, const bctx_t *bc, const cconst_t *K, u128 Mblock)
{
    const prob_t *pr = st->pr;
    const int k = pr->k;
    const u64 B = bc->B, c = K->c, N = K->N, P = K->P, limit = K->limit;
    const u64 Wy = K->Wy, NWy = K->NWy, Wx = K->Wx, NWx = K->NWx;
    const u64 wy = bc->wy, dx = bc->dx, s = bc->s, Rcross = K->Rcross;
    const u64 g = K->g, Np = K->Np;

    u64 rem = (u64)(Mblock % P);
    u128 q128 = Mblock / P;
    if (q128 >= N) return;             /* M >= c^(n-1) already */
    u64 q = (u64)q128;
    u64 E[MAXK];
    for (int j = k - 1; j >= 0; j--) { E[j] = q % c; q /= c; }
    u64 T = 0;
    for (int j = k - 1; j >= 0; j--) T = T * c + E[j];
    u64 MN = (u64)(Mblock % N);
    u64 t = T >= MN ? T - MN : T + N - MN;
    u64 w = mulmod64(t, K->inv, N);

    u64 wnext;
    int last;
    u64 elast, trig;
    /* q -> q+1 with the carry stopping at digit m (E[m+1..k-1] = c-1 wrap to 0,
       E[m] += 1) changes T by -(c-1)(c^(m+1) + ... + c^(k-1)) + c^m, which is
       c^(m+1) + c^m (mod c^k) for m < k-1 and c^(k-1) for m = k-1. */
#define RECOMPUTE_NEXT()                                                        \
    do {                                                                        \
        int jj = k - 1;                                                         \
        while (jj >= 0 && E[jj] == c - 1) jj--;                                 \
        last = jj < 0;                                                          \
        if (jj == k - 1) wnext = K->Wc[k - 1];                                  \
        else if (jj >= 0) {                                                     \
            wnext = K->Wc[jj + 1] + K->Wc[jj];                                  \
            if (wnext >= N) wnext -= N;                                         \
        } else wnext = 0;                                                       \
    } while (0)
    RECOMPUTE_NEXT();
    elast = E[k - 1];
    trig = (elast == c - 1) ? c : c - 1;

#ifdef PARANOID
#define CAND(zv)                                                                \
    do {                                                                        \
        if (nfound < 256) found[nfound++] = (zv);                               \
        candidate(st, M0 + (u128)(zv) * s, B, c);                               \
    } while (0)
#else
#define CAND(zv) candidate(st, M0 + (u128)(zv) * s, B, c)
#endif
    for (u64 x = 0; x < B; x++) {
        for (u64 y = 0; y < B; y++) {
#ifdef PARANOID
            u64 found[256];
            int nfound = 0;
            {
                u128 M0p = Mblock + (u128)x * bc->wx + (u128)y * wy;
                paranoid_state(pr, bc, K, M0p, w, rem, E, elast, wnext, last);
            }
#endif
            u64 w1 = w + wnext;
            w1 = (w1 >= N) ? w1 - N : w1;
            int hit = (w < limit) | ((rem >= Rcross) & (w1 < limit));
            if (__builtin_expect(hit, 0)) {
                u128 M0 = Mblock + (u128)x * bc->wx + (u128)y * wy;
                if (w < limit && (g == 1 || w % g == 0)) {
                    for (u64 z = w / g; z < B; z += Np) {
                        if (rem + z * s >= P) break;
                        CAND(z);
                    }
                }
                if (!last && rem >= Rcross && w1 < limit && (g == 1 || w1 % g == 0)) {
                    u64 z = w1 / g;
                    u64 zs = (P - rem + s - 1) / s;    /* first z in the q+1 range */
                    if (z < zs) z += ((zs - z + Np - 1) / Np) * Np;
                    for (; z < B; z += Np) CAND(z);
                }
            }
#ifdef PARANOID
            if (B <= 64 || (x * B + y) % 97 == 0) {   /* brute force over z */
                u128 M0 = Mblock + (u128)x * bc->wx + (u128)y * wy;
                u64 want[256];
                int nw = 0;
                for (u64 z = 0; z < B; z++) {
                    u128 M = M0 + (u128)z * s;
                    u128 qz = M / P;
                    if (qz >= N) continue;
                    if ((u64)(M % N) == T_of((u64)qz, c, k) && nw < 256) want[nw++] = z;
                }
                qsort(found, nfound, sizeof(u64), u64cmp);
                int same = nw == nfound;
                for (int i = 0; same && i < nw; i++) same = want[i] == found[i];
                if (!same) {
                    char buf[48];
                    fprintf(stderr, "PARANOID z mismatch: n=%d B=%llu c=%llu M0=%s fast %d brute %d\n", pr->n,
                            (unsigned long long)B, (unsigned long long)c, u128s(M0, buf), nfound, nw);
                    abort();
                }
                st->pcap++;                        /* count checked steps */
                st->np += nw;                      /* count checked z hits */
            }
#endif
            /* advance to the next (x, y) */
            u64 dW, nW, dR;
            if (y + 1 < B) { dW = Wy; nW = NWy; dR = wy; }
            else { dW = Wx; nW = NWx; dR = dx; }
            w = (w >= dW) ? w - dW : w + nW;
            rem += dR;
            u64 inc = rem >= P;
            rem -= inc ? P : 0;
            u64 wa = w + (inc ? wnext : 0);
            w = (wa >= N) ? wa - N : wa;
            elast += inc;
            if (__builtin_expect(elast >= trig, 0)) {
                if (elast == c - 1) {          /* next increment carries */
                    E[k - 1] = c - 1;
                    RECOMPUTE_NEXT();
                    trig = c;
                } else {                       /* elast == c: carry happened */
                    if (last) goto done;       /* q reached c^k: M >= c^(n-1) from here on */
                    E[k - 1] = 0;
                    int j = k - 2;
                    while (j >= 0) {
                        E[j]++;
                        if (E[j] < c) break;
                        E[j] = 0;
                        j--;
                    }
                    if (j < 0) goto done;
                    elast = 0;
                    RECOMPUTE_NEXT();
                    trig = c - 1;
                }
            }
        }
    }
done:
#undef RECOMPUTE_NEXT
#undef CAND
    return;
}

/* process one block for all its bases c */
static void do_block(wstat_t *st, const bctx_t *bc, u128 Mblock)
{
    const prob_t *pr = st->pr;
    u128 Mmin = Mblock, Mmax = Mblock + bc->span_block;
    if (Mmax < pr->lo || Mmin > pr->hi) return;
    if (Mmin < pr->lo) Mmin = pr->lo;
    if (Mmax > pr->hi) Mmax = pr->hi;
    u64 clo = iroot(Mmin, pr->n - 1) + 1, chi = iroot(Mmax, pr->n - 2);
    if (clo < bc->cbase) clo = bc->cbase;
    if (chi > bc->cbase + bc->ccount - 1) chi = bc->cbase + bc->ccount - 1;
    for (u64 c = clo; c <= chi; c++) block_c(st, bc, &bc->tab[c - bc->cbase], Mblock);
    if (chi >= clo) st->tuples += (chi - clo + 1) * bc->B * bc->B;
}

/* unit u -> base and prefix digits */
static u64 unit_decode(const prob_t *pr, u64 u, u64 *pref)
{
    u64 lo = 0, hi = pr->Bmax - pr->Bmin;   /* find i with uoff[i] <= u < uoff[i+1] */
    while (lo < hi) {
        u64 mid = (lo + hi + 1) / 2;
        if (pr->uoff[mid] <= u) lo = mid; else hi = mid - 1;
    }
    u64 B = pr->Bmin + lo;
    u64 local = u - pr->uoff[lo];
    for (int i = pr->ud - 1; i >= 1; i--) { pref[i] = local % B; local /= B; }
    pref[0] = 1 + local;
    return B;
}

static void do_unit(wstat_t *st, bctx_t *bc, u64 u)
{
    const prob_t *pr = st->pr;
    u64 pref[MAXN];
    u64 B = unit_decode(pr, u, pref);
    if (bc->B != B) bctx_set(bc, pr, B);
    u128 Mu = 0;
    for (int i = 0; i < pr->ud; i++) Mu += (u128)pref[i] * bc->w[i];
    if (Mu + bc->span_unit < pr->lo || Mu > pr->hi) return;
    int nb = pr->bd - pr->ud;          /* extra block digits: 0 or 1 */
    if (nb == 0) {
        do_block(st, bc, Mu);
    } else {
        u128 wb = bc->w[pr->ud];
        for (u64 d = 0; d < B; d++) do_block(st, bc, Mu + (u128)d * wb);
    }
}

/* ------------------------------------------------------------------ */
/* naive references                                                    */
/* ------------------------------------------------------------------ */

/* all (M, c) with M in the block an (n-1)-digit palindrome in base c */
static int ends_mirror(u128 M, u64 c, int L, int k)
{
    u64 d[130];
    int len = 0;
    u128 x = M;
    while (x) { d[len++] = (u64)(x % c); x /= c; }
    if (len != L) return 0;
    for (int i = 0; i < k; i++)
        if (d[i] != d[L - 1 - i]) return 0;
    return 1;
}

static void naive_block_pairs(const prob_t *pr, const bctx_t *bc, u128 Mblock, pair_t **out, size_t *np, int mode)
{
    size_t cap = 1024, cnt = 0;
    pair_t *p = malloc(cap * sizeof(pair_t));
    u64 B = bc->B;
    u128 xspan = (u128)(B - 1) * (bc->wy + bc->s);
    for (u64 x = 0; x < B; x++) {
        u128 Mx = Mblock + (u128)x * bc->wx;
        if (Mx + xspan < pr->lo || Mx > pr->hi) continue;
        for (u64 y = 0; y < B; y++)
            for (u64 z = 0; z < B; z++) {
                u128 M = Mblock + (u128)x * bc->wx + (u128)y * bc->wy + (u128)z * bc->s;
                if (M < pr->lo || M > pr->hi) continue;
                u64 clo = iroot(M, pr->n - 1) + 1, chi = iroot(M, pr->n - 2);
                for (u64 c = clo; c <= chi; c++) {
                    u64 r = mod128_64(M, c);
                    if (!r) continue;
                    u128 Pc = pow_sat(c, pr->n - 2);
                    if (M / Pc != r) continue;
                    if (mode == 2 ? !ends_mirror(M, c, pr->n - 1, pr->k) : !is_pal_base(M, c, pr->n - 1)) continue;
                    if (cnt == cap) { cap *= 2; p = realloc(p, cap * sizeof(pair_t)); }
                    p[cnt].M = M;
                    p[cnt].c = c;
                    cnt++;
                }
            }
    }
    *out = p;
    *np = cnt;
}

static int pair_cmp(const void *a, const void *b)
{
    const pair_t *x = a, *y = b;
    if (x->M != y->M) return x->M < y->M ? -1 : 1;
    if (x->c != y->c) return x->c < y->c ? -1 : 1;
    return 0;
}

/* plain enumeration of n-digit palindromes in [lo, hi], all lengths checked */
static u128 naive_search(int n, u128 lo, u128 hi, int print)
{
    u128 best = 0;
    char buf[48];
    if (n <= 1) return lo <= 1 ? 1 : lo;
    if (n == 2) return lo <= 3 ? 3 : lo;
    u64 Bmax = iroot(hi, n - 1);
    int D = (n + 1) / 2, h = n / 2;
    for (u64 B = 2; B <= Bmax; B++) {
        u128 w[MAXN];
        for (int i = 0; i < h; i++) w[i] = pow_sat(B, n - 1 - i) + pow_sat(B, i);
        if (n & 1) w[h] = pow_sat(B, h);
        u64 d[MAXN] = {0};
        d[0] = 1;
        for (;;) {
            u128 M = 0;
            for (int i = 0; i < D; i++) M += (u128)d[i] * w[i];
            if (M >= lo && M <= hi && (!best || M < best)) {
                int ok = 1;
                for (int L = n - 1; L >= 3 && ok; L--) ok = find_pal_base(M, L) != 0;
                if (ok) {
                    best = M;
                    if (print) printf("  naive n=%d: %s (B=%llu)\n", n, u128s(M, buf), (unsigned long long)B);
                }
            }
            int i = D - 1;
            while (i >= 0) {
                d[i]++;
                if (d[i] < B) break;
                d[i] = (i == 0) ? 1 : 0;
                i--;
            }
            if (i < 0) break;
        }
    }
    return best;
}

/* ------------------------------------------------------------------ */
/* threaded search                                                     */
/* ------------------------------------------------------------------ */

static volatile sig_atomic_t g_stop;
static void on_sig(int s) { (void)s; g_stop = 1; }

typedef struct {
    prob_t *pr;
    _Atomic u64 next;
    unsigned char *done;
    u64 frontier;
    pthread_mutex_t lock;
    wstat_t *st;                /* per thread */
    int nthreads;
    int verbose;
} search_t;

typedef struct { search_t *S; int id; } warg_t;

static void *worker(void *arg)
{
    warg_t *a = arg;
    search_t *S = a->S;
    wstat_t *st = &S->st[a->id];
    bctx_t bc;
    memset(&bc, 0, sizeof bc);
    while (!g_stop) {
        u64 u = atomic_fetch_add(&S->next, 1);
        if (u >= S->pr->nunits) break;
        do_unit(st, &bc, u);
        pthread_mutex_lock(&S->lock);
        S->done[u] = 1;
        while (S->frontier < S->pr->nunits && S->done[S->frontier]) S->frontier++;
        pthread_mutex_unlock(&S->lock);
    }
    free(bc.tab);
    return NULL;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void sum_stats(search_t *S, wstat_t *tot)
{
    memset(tot, 0, sizeof *tot);
    for (int i = 0; i < S->nthreads; i++) {
        wstat_t *s = &S->st[i];
        tot->tuples += s->tuples;
        tot->hits += s->hits;
        tot->st2 += s->st2;
        tot->canon += s->canon;
        tot->csum += s->csum;
        tot->sol += s->sol;
        for (int L = 0; L <= MAXN; L++) tot->fail[L] += s->fail[L];
    }
}

/* estimated tuples of a unit (for progress only) */
static double unit_weight(const prob_t *pr, u64 B, u64 d0)
{
    int n = pr->n;
    u128 Mmin = (u128)d0 * pow_sat(B, n - 1), Mmax = (u128)(d0 + 1) * pow_sat(B, n - 1);
    if (Mmax < pr->lo || Mmin > pr->hi) return 0;
    if (Mmin < pr->lo) Mmin = pr->lo;
    if (Mmax > pr->hi) Mmax = pr->hi;
    double cl = pow((double)Mmin, 1.0 / (n - 1)), ch = pow((double)Mmax, 1.0 / (n - 2));
    double nc = ch - cl + 1;
    if (nc < 0) nc = 0;
    double frac = (double)(Mmax - Mmin) / (double)pow_sat(B, n - 1);
    return nc * pow((double)B, pr->D - 2) * frac;   /* blocks (d0 fixed) x B^2 */
}

static int load_state(const char *path, const prob_t *pr, u64 *frontier, wstat_t *acc)
{
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    char lo[64], hi[64];
    int n;
    unsigned long long fr, tu, hi_, s2, ca, so, cs = 0;
    int ok = fscanf(f, "n=%d lo=%63s hi=%63s frontier=%llu tuples=%llu hits=%llu st2=%llu canon=%llu sol=%llu",
                    &n, lo, hi, &fr, &tu, &hi_, &s2, &ca, &so) == 9;
    if (ok && fscanf(f, " csum=%llu", &cs) != 1) cs = 0;
    fclose(f);
    if (!ok) { fprintf(stderr, "cannot parse state file %s\n", path); exit(1); }
    if (n != pr->n || parse_or_die(lo) != pr->lo || parse_or_die(hi) != pr->hi) {
        fprintf(stderr, "state file %s is for a different search\n", path);
        exit(1);
    }
    *frontier = fr;
    acc->tuples = tu; acc->hits = hi_; acc->st2 = s2; acc->canon = ca; acc->sol = so; acc->csum = cs;
    return 1;
}

static void save_state(const char *path, const prob_t *pr, u64 frontier, const wstat_t *t)
{
    char tmp[4096], a[48], b[48];
    snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror(tmp); return; }
    fprintf(f, "n=%d lo=%s hi=%s frontier=%llu tuples=%llu hits=%llu st2=%llu canon=%llu sol=%llu csum=%llu\n",
            pr->n, u128s(pr->lo, a), u128s(pr->hi, b), (unsigned long long)frontier, (unsigned long long)t->tuples,
            (unsigned long long)t->hits, (unsigned long long)t->st2, (unsigned long long)t->canon,
            (unsigned long long)t->sol, (unsigned long long)t->csum);
    fclose(f);
    rename(tmp, path);
}

static int cmd_search(int n, u128 lo, u128 hi, int nthreads, const char *state, double interval, int verbose)
{
    char a[48], b[48];
    if (n <= 6) {
        u128 m = naive_search(n, lo, hi, 1);
        printf("naive: smallest solution %s\n", m ? u128s(m, a) : "none");
        return 0;
    }
    prob_t pr;
    prob_init(&pr, n, lo, hi);
    search_t S;
    memset(&S, 0, sizeof S);
    S.pr = &pr;
    S.nthreads = nthreads;
    S.verbose = verbose;
    S.done = calloc(pr.nunits + 1, 1);
    S.st = calloc(nthreads, sizeof(wstat_t));
    pthread_mutex_init(&S.lock, NULL);
    wstat_t base;
    memset(&base, 0, sizeof base);
    u64 start = 0;
    if (state && load_state(state, &pr, &start, &base))
        fprintf(stderr, "resuming at unit %llu of %llu\n", (unsigned long long)start, (unsigned long long)pr.nunits);
    for (u64 u = 0; u < start; u++) S.done[u] = 1;
    S.frontier = start;
    atomic_store(&S.next, start);
    for (int i = 0; i < nthreads; i++) { S.st[i].pr = &pr; S.st[i].verbose = verbose; }

    /* progress weights */
    double *cw = malloc((pr.nunits + 1) * sizeof(double));
    cw[0] = 0;
    for (u64 u = 0; u < pr.nunits; u++) {
        u64 pref[MAXN];
        u64 B = unit_decode(&pr, u, pref);
        double wgt = unit_weight(&pr, B, pref[0]);
        if (pr.ud > 1) wgt /= (double)B;
        cw[u + 1] = cw[u] + wgt;
    }
    double wtot = cw[pr.nunits] > 0 ? cw[pr.nunits] : 1;

    fprintf(stderr, "search n=%d [%s, %s]: B in [%llu, %llu], %llu units, k=%d, %d threads\n", n, u128s(lo, a),
            u128s(hi, b), (unsigned long long)pr.Bmin, (unsigned long long)pr.Bmax,
            (unsigned long long)pr.nunits, pr.k, nthreads);
    signal(SIGINT, on_sig);
    signal(SIGTERM, on_sig);
    pthread_t *th = malloc(nthreads * sizeof(pthread_t));
    warg_t *wa = malloc(nthreads * sizeof(warg_t));
    double t0 = now_sec(), tlast = t0, tsave = t0;
    for (int i = 0; i < nthreads; i++) {
        wa[i].S = &S;
        wa[i].id = i;
        pthread_create(&th[i], NULL, worker, &wa[i]);
    }
    double w0 = cw[start];
    for (;;) {
        usleep(200000);
        pthread_mutex_lock(&S.lock);
        u64 fr = S.frontier;
        pthread_mutex_unlock(&S.lock);
        u64 nx = atomic_load(&S.next);
        int finished = fr >= pr.nunits;
        double t = now_sec();
        if (t - tlast >= interval || finished || g_stop) {
            wstat_t tot;
            sum_stats(&S, &tot);
            double el = t - t0;
            double wd = cw[fr < pr.nunits ? fr : pr.nunits];
            double prog = wd / wtot;
            double rate = tot.tuples / (el > 0 ? el : 1);
            double eta = (wd - w0) > 0 ? (wtot - wd) * el / (wd - w0) : 0;
            u64 pref[MAXN];
            u64 Bc = unit_decode(&pr, fr < pr.nunits ? fr : pr.nunits - 1, pref);
            fprintf(stderr,
                    "[%7.0fs] frontier %llu/%llu (B=%llu) %.3f%%  %.3e tuples/s  st2 %llu canon %llu sol %llu  ETA %.0fs\n",
                    el, (unsigned long long)fr, (unsigned long long)pr.nunits, (unsigned long long)Bc, 100 * prog,
                    rate, (unsigned long long)(base.st2 + tot.st2), (unsigned long long)(base.canon + tot.canon),
                    (unsigned long long)(base.sol + tot.sol), eta);
            tlast = t;
            if (state && (t - tsave >= 60 || finished || g_stop)) {
                wstat_t all = tot;
                all.tuples += base.tuples; all.hits += base.hits; all.st2 += base.st2;
                all.canon += base.canon; all.sol += base.sol; all.csum += base.csum;
                save_state(state, &pr, fr, &all);
                tsave = t;
            }
        }
        if (finished || (g_stop && nx >= pr.nunits)) break;
        if (g_stop) break;
    }
    for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
    wstat_t tot;
    sum_stats(&S, &tot);
    if (state) {
        wstat_t all = tot;
        all.tuples += base.tuples; all.hits += base.hits; all.st2 += base.st2;
        all.canon += base.canon; all.sol += base.sol; all.csum += base.csum;
        save_state(state, &pr, S.frontier, &all);
    }
    double el = now_sec() - t0;
    printf("# n=%d [%s, %s] %s: %.1fs, tuples %llu, hits %llu, (n-1)-palindromes %llu, canonical %llu, "
           "checksum %016llx, solutions %llu\n",
           n, u128s(lo, a), u128s(hi, b), S.frontier >= pr.nunits ? "complete" : "INTERRUPTED", el,
           (unsigned long long)(base.tuples + tot.tuples), (unsigned long long)(base.hits + tot.hits),
           (unsigned long long)(base.st2 + tot.st2), (unsigned long long)(base.canon + tot.canon),
           (unsigned long long)(base.csum + tot.csum), (unsigned long long)(base.sol + tot.sol));
    printf("# failed at length:");
    for (int L = n - 2; L >= 3; L--) printf(" L=%d:%llu", L, (unsigned long long)tot.fail[L]);
    printf("\n");
    return 0;
}

/* ------------------------------------------------------------------ */
/* selftest                                                            */
/* ------------------------------------------------------------------ */

static u64 rng_state = 88172645463325252ULL;
static u64 rnd(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}

/* compare fast (M, c) pairs with naive ones on one block */
/* win > 0: restrict to M in [Mblock + x0 w_x, Mblock + (x0 + win) w_x] for a random x0 */
static int compare_block(int n, u64 B, u64 blockidx, u128 lo, u128 hi, int quiet, int win, int mode, size_t *npairs)
{
    prob_t pr;
    prob_init(&pr, n, lo, hi);
    bctx_t bc;
    memset(&bc, 0, sizeof bc);
    bctx_set(&bc, &pr, B);
    /* block digits d_0..d_{bd-1} from blockidx: d_0 = 1 + idx / B^(bd-1) */
    u64 d[MAXN];
    u64 t = blockidx;
    for (int i = pr.bd - 1; i >= 1; i--) { d[i] = t % B; t /= B; }
    d[0] = 1 + t;
    if (d[0] >= B) { free(bc.tab); free(pr.uoff); return 1; }
    u128 Mb = 0;
    for (int i = 0; i < pr.bd; i++) Mb += (u128)d[i] * bc.w[i];
    if (win) {
        u64 x0 = rnd() % B;
        u128 wlo = Mb + (u128)x0 * bc.wx, whi = wlo + (u128)win * bc.wx;
        if (wlo < lo) wlo = lo;
        if (whi > hi) whi = hi;
        if (wlo > whi) { wlo = Mb; whi = Mb + (u128)win * bc.wx; }
        free(pr.uoff);
        prob_init(&pr, n, wlo, whi);
        bctx_set(&bc, &pr, B);
    }
    wstat_t st;
    memset(&st, 0, sizeof st);
    st.pr = &pr;
    st.record = mode;
    do_block(&st, &bc, Mb);
    pair_t *np;
    size_t nn;
    naive_block_pairs(&pr, &bc, Mb, &np, &nn, mode);
    if (npairs) *npairs += nn;
    qsort(st.pairs, st.np, sizeof(pair_t), pair_cmp);
    qsort(np, nn, sizeof(pair_t), pair_cmp);
    int ok = st.np == nn;
    for (size_t i = 0; ok && i < nn; i++) ok = st.pairs[i].M == np[i].M && st.pairs[i].c == np[i].c;
    char buf[48];
    if (!ok || !quiet)
        printf("  n=%d B=%llu block %llu (Mblock=%s) mode %d: fast %zu pairs, naive %zu pairs, hits %llu  %s\n", n,
               (unsigned long long)B, (unsigned long long)blockidx, u128s(Mb, buf), mode, st.np, nn,
               (unsigned long long)st.hits, ok ? "ok" : "MISMATCH");
    free(st.pairs);
    free(np);
    free(bc.tab);
    free(pr.uoff);
    return ok;
}

static int check_value(u128 M, int n, int print)
{
    char buf[48];
    int ok = 1;
    if (print) printf("M = %s\n", u128s(M, buf));
    for (int L = 2; L <= n; L++) {
        u64 b = find_pal_base(M, L);
        if (print) {
            if (b) {
                printf("  L=%2d base %llu digits:", L, (unsigned long long)b);
                u128 x = M;
                u64 dg[130];
                int len = 0;
                while (x) { dg[len++] = (u64)(x % b); x /= b; }
                if (L <= 12)
                    for (int i = len - 1; i >= 0; i--) printf(" %llu", (unsigned long long)dg[i]);
                printf("\n");
            } else printf("  L=%2d none\n", L);
        }
        if (!b) ok = 0;
    }
    return ok;
}

static int cmd_selftest(int full, int nthreads)
{
    int bad = 0;
    char buf[48];
    const char *known[] = {"1", "3", "5", "52", "130", "1885", "2^30", "2^42", "2^56"};
    printf("1. naive a(1..6)\n");
    for (int n = 1; n <= 6; n++) {
        u128 m = naive_search(n, 1, 100000, 0);
        u128 want = parse_or_die(known[n - 1]);
        printf("  a(%d) = %s %s\n", n, u128s(m, buf), m == want ? "ok" : "WRONG");
        bad += m != want;
    }
    printf("2. 2^((n-1)(n-2)) satisfies lengths 2..n (n = 7..12), and fails n+1\n");
    for (int n = 7; n <= 12; n++) {
        u128 M = (u128)1 << ((n - 1) * (n - 2));
        int ok = check_value(M, n, 0);
        printf("  n=%d 2^%d: %s\n", n, (n - 1) * (n - 2), ok ? "ok" : "FAIL");
        bad += !ok;
    }
    printf("3. fast filter vs naive on random blocks\n");
    struct { int n; u64 Bmin, Bmax; int reps, win; } cfg[] = {
        {7, 3, 40, 60, 0}, {8, 3, 30, 60, 0}, {9, 3, 24, 60, 0}, {10, 3, 20, 60, 0}, {11, 3, 14, 40, 0},
        {12, 3, 10, 20, 0}, {9, 60, 90, 4, 0}, {10, 40, 60, 3, 0}, {11, 30, 40, 2, 0},
        {10, 100, 255, 8, 2}, {11, 200, 511, 4, 1}, {9, 100, 127, 8, 2}, {8, 50, 63, 8, 0},
    };
    for (size_t ci = 0; ci < sizeof cfg / sizeof cfg[0]; ci++) {
        int n = cfg[ci].n, okc = 0, tot = 0;
        size_t pairs1 = 0, pairs2 = 0;
        for (int r = 0; r < cfg[ci].reps; r++) {
            u64 B = cfg[ci].Bmin + rnd() % (cfg[ci].Bmax - cfg[ci].Bmin + 1);
            u128 hi = U128_MAX >> 1;
            if (n == 10) hi = (u128)1 << 72;     /* production limits: c^k and c^(h+1) fit */
            if (n == 11) hi = (u128)1 << 90;
            if (n == 12) hi = (u128)1 << 80;
            u128 dm = hi / pow_sat(B, n - 1);    /* leading digits with Mblock <= hi */
            u64 dmax = dm < B - 1 ? (u64)dm : B - 1;
            if (dmax < 1) { r--; continue; }
            int bdg = (n + 1) / 2 - 3;
            u64 nblocks = dmax * (u64)pow_sat(B, bdg - 1);
            u64 bi = rnd() % nblocks;
            if (r == 0) bi = 0;
            if (r == 1) bi = nblocks - 1;
            u64 save = rng_state;
            int ok = compare_block(n, B, bi, 1, hi, 1, cfg[ci].win, 1, &pairs1);
            rng_state = save;                    /* same window for mode 2 */
            ok &= compare_block(n, B, bi, 1, hi, 1, cfg[ci].win, 2, &pairs2);
            okc += ok;
            tot++;
        }
        printf("  n=%d B in [%llu,%llu]%s: %d/%d blocks agree (%zu palindrome pairs, %zu k-digit matches)\n", n,
               (unsigned long long)cfg[ci].Bmin, (unsigned long long)cfg[ci].Bmax, cfg[ci].win ? " windowed" : "",
               okc, tot, pairs1, pairs2);
        bad += okc != tot;
    }
    printf("4. full searches reproduce a(7), a(8)%s\n", full ? ", a(9)" : "");
    int nmax = full ? 9 : 8;
    for (int n = 7; n <= nmax; n++) {
        u128 want = parse_or_die(known[n - 1]);
        prob_t pr;
        prob_init(&pr, n, 1, want);
        search_t S;
        memset(&S, 0, sizeof S);
        S.pr = &pr;
        S.nthreads = nthreads;
        S.done = calloc(pr.nunits + 1, 1);
        S.st = calloc(nthreads, sizeof(wstat_t));
        pthread_mutex_init(&S.lock, NULL);
        for (int i = 0; i < nthreads; i++) S.st[i].pr = &pr;
        /* capture solutions */
        FILE *mem = tmpfile();
        FILE *save = sol_out;
        sol_out = mem;
        double t0 = now_sec();
        pthread_t *th = malloc(nthreads * sizeof(pthread_t));
        warg_t *wa = malloc(nthreads * sizeof(warg_t));
        for (int i = 0; i < nthreads; i++) {
            wa[i].S = &S;
            wa[i].id = i;
            pthread_create(&th[i], NULL, worker, &wa[i]);
        }
        for (int i = 0; i < nthreads; i++) pthread_join(th[i], NULL);
        sol_out = save;
        wstat_t tot;
        sum_stats(&S, &tot);
        rewind(mem);
        char line[1024];
        u128 best = 0;
        int nsol = 0;
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
        int ok = best == want && nsol == 1;
        printf("  n=%d [1, %s]: %.2fs, %.3e tuples, %llu (n-1)-pals, %llu canonical, %d solution(s), smallest %s  %s\n", n,
               known[n - 1], now_sec() - t0, (double)tot.tuples, (unsigned long long)tot.st2,
               (unsigned long long)tot.canon, nsol, best ? u128s(best, buf) : "none", ok ? "ok" : "WRONG");
        bad += !ok;
        free(th); free(wa); free(S.done); free(S.st); free(pr.uoff);
    }
    printf(bad ? "SELFTEST FAILED (%d)\n" : "selftest passed\n", bad);
    return bad ? 1 : 0;
}

/* ------------------------------------------------------------------ */

#ifdef PARANOID
/* run the real block code with full state checks on many random blocks */
static int cmd_paranoid(int reps)
{
    struct { int n; u64 Bmin, Bmax; u128 hi; } cfg[] = {
        {7, 3, 60, (u128)1 << 62}, {8, 3, 60, (u128)1 << 62}, {9, 3, 60, (u128)1 << 62},
        {10, 3, 60, (u128)1 << 72}, {11, 3, 60, (u128)1 << 90}, {12, 3, 40, (u128)1 << 80},
        {10, 100, 255, (u128)1 << 72}, {11, 200, 511, (u128)1 << 90}, {9, 64, 127, (u128)1 << 62},
    };
    FILE *devnull = fopen("/dev/null", "w");
    FILE *save = sol_out;
    sol_out = devnull;
    for (size_t ci = 0; ci < sizeof cfg / sizeof cfg[0]; ci++) {
        int n = cfg[ci].n;
        u64 steps = 0, zh = 0, blocks = 0;
        double t0 = now_sec();
        int r;
        for (r = 0; r < reps && now_sec() - t0 < 60; r++) {
            u64 B = cfg[ci].Bmin + rnd() % (cfg[ci].Bmax - cfg[ci].Bmin + 1);
            u128 hi = cfg[ci].hi;
            u128 dm = hi / pow_sat(B, n - 1);
            u64 dmax = dm < B - 1 ? (u64)dm : B - 1;
            if (dmax < 1) continue;
            prob_t pr;
            prob_init(&pr, n, 1, hi);
            bctx_t bc;
            memset(&bc, 0, sizeof bc);
            bctx_set(&bc, &pr, B);
            u64 d[MAXN];
            d[0] = 1 + rnd() % dmax;
            for (int i = 1; i < pr.bd; i++) d[i] = rnd() % B;
            u128 Mb = 0;
            for (int i = 0; i < pr.bd; i++) Mb += (u128)d[i] * bc.w[i];
            wstat_t st;
            memset(&st, 0, sizeof st);
            st.pr = &pr;
            do_block(&st, &bc, Mb);
            steps += st.pcap;
            zh += st.np;
            blocks++;
            free(bc.tab);
            free(pr.uoff);
        }
        printf("  paranoid n=%d B in [%llu,%llu]: %llu blocks, every step state-checked, %llu brute-force steps, %llu z hits: ok\n",
               n, (unsigned long long)cfg[ci].Bmin, (unsigned long long)cfg[ci].Bmax, (unsigned long long)blocks,
               (unsigned long long)steps, (unsigned long long)zh);
        fflush(stdout);
    }
    sol_out = save;
    return 0;
}
#endif

static void usage(void)
{
    fprintf(stderr,
            "usage: a171775 selftest [full] [-t T]\n"
            "       a171775 search n LO HI [-t T] [-S state] [-i sec] [-v]\n"
            "       a171775 check M [n]\n"
            "       a171775 naive n LO HI\n");
    exit(1);
}

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
        else if (npos < 8) pos[npos++] = argv[i];
        else usage();
    }
    if (nthreads < 1) nthreads = 1;
    if (npos < 1) usage();
    if (!strcmp(pos[0], "selftest")) return cmd_selftest(npos > 1 && !strcmp(pos[1], "full"), nthreads);
#ifdef PARANOID
    if (!strcmp(pos[0], "paranoid")) return cmd_paranoid(npos > 1 ? atoi(pos[1]) : 200);
#endif
    if (!strcmp(pos[0], "search") && npos == 4) {
        int n = atoi(pos[1]);
        if (n < 3 || n > 20) usage();
        return cmd_search(n, parse_or_die(pos[2]), parse_or_die(pos[3]), nthreads, state, interval, verbose);
    }
    if (!strcmp(pos[0], "check") && npos >= 2) {
        u128 M = parse_or_die(pos[1]);
        int n = npos > 2 ? atoi(pos[2]) : 12;
        return check_value(M, n, 1) ? 0 : 2;
    }
    if (!strcmp(pos[0], "naive") && npos == 4) {
        char buf[48];
        u128 m = naive_search(atoi(pos[1]), parse_or_die(pos[2]), parse_or_die(pos[3]), 1);
        printf("smallest: %s\n", m ? u128s(m, buf) : "none");
        return 0;
    }
    usage();
    return 1;
}
