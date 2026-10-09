/*
 * a171740.c — search terms of OEIS A171740:
 *   a(n) = the first number that is palindromic with exactly n digits
 *          in more than one base.
 *
 * Method: a number N has n digits in base b iff b^(n-1) <= N < b^n.
 * For two bases b1 < b2 both to give n digits we need b2^(n-1) < b1^n.
 * For each "stage" b2 = 3, 4, 5, ... (increasing, so stage lower bounds
 * b2^(n-1) increase monotonically) we enumerate the n-digit palindromes
 * of base b2 in increasing order via their leading half, and test each
 * value for n-digit palindromicity in every valid smaller base b1.
 * The first match, minimized across all stages whose lower bound is
 * below it, is a(n).
 *
 * Enumerating in the *larger* base of each pair minimizes candidates:
 * there are only ~(range)/b2^floor(n/2) palindromes of base b2 in a range.
 *
 * Build:  cc -O3 -o a171740 a171740.c -lpthread
 * Usage:  ./a171740 -n 18 [-t threads] [-L decimal_limit] [-q]
 *   -n N   digit count (required)
 *   -t T   worker threads (default: all cores)
 *   -L X   only search values < X (exclusive); useful to confirm a
 *          known candidate is minimal. Default: unbounded.
 *   -q     suppress periodic status lines
 * Status/progress to stderr; result to stdout.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include <time.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

#define U128MAX (~(u128)0)
#define MAXB1 24

static int g_n, g_k, g_m, g_odd;    /* n, floor(n/2), ceil(n/2), n&1 */
static int g_nthreads = 0;
static u128 g_limit = 0;            /* 0 = unbounded */
static int g_quiet = 0;
static int g_mitm = 0;              /* -M: meet-in-the-middle scanner */
static u64 g_tablecap = 300000000;  /* max MITM table entries (yl side) */
static int g_limitw_set = 0;        /* -L given (wide) */

/* ---------------- u128 utilities ---------------- */

/* b^e saturating at U128MAX */
static u128 pow_sat(u64 b, int e){
    u128 r = 1;
    while (e-- > 0){
        if (r > U128MAX / b) return U128MAX;
        r *= b;
    }
    return r;
}

static u64 pw64(u64 b, int e){ u64 r = 1; while (e-- > 0) r *= b; return r; }

static unsigned gcd_u(unsigned a, unsigned b){
    while (b){ unsigned t = a % b; a = b; b = t; }
    return a;
}

/* inverse of a mod m (gcd(a,m)=1); returns 0 for m==1 */
static unsigned modinv_u(unsigned a, unsigned m){
    if (m == 1) return 0;
    long long t = 0, nt = 1, r = m, nr = a % m;
    while (nr){
        long long qq = r / nr, tmp;
        tmp = t - qq*nt; t = nt; nt = tmp;
        tmp = r - qq*nr; r = nr; nr = tmp;
    }
    return (unsigned)((t % m + m) % m);
}

/* divide *v by small d, return remainder; fast path for 64-bit values,
 * limb long-division otherwise (compiles to multiplies when d is a
 * compile-time constant via the specialized scan instantiations). */
static inline unsigned udm128(u128 *v, unsigned d){
    u64 hi = (u64)(*v >> 64), lo = (u64)*v;
    if (!hi){ *v = lo / d; return (unsigned)(lo % d); }
    u64 qh = hi / d, r = hi % d;
    u64 a  = (r << 32) | (lo >> 32);
    u64 q1 = a / d;  r = a % d;
    u64 b  = (r << 32) | (lo & 0xffffffffu);
    u64 q0 = b / d;  r = b % d;
    *v = ((u128)qh << 64) | (q1 << 32) | q0;
    return (unsigned)r;
}

/* ---------------- w192: 3-limb wide values ----------------
 * Values beyond u128 (needed from n ~ 28) live in a little-endian
 * 3x64-bit struct. Only the operations the search needs are provided;
 * range tops out near 6.3e57 (b^n for n ~ 40). The legacy scanners stay
 * u128 internally and only run on stages whose window fits u128. */

typedef struct { u64 d[3]; } w192;

static inline w192 w_from_u128(u128 x){
    w192 r = {{ (u64)x, (u64)(x >> 64), 0 }}; return r;
}
static inline int w_fits_u128(w192 a){ return a.d[2] == 0; }
static inline u128 w_to_u128(w192 a){ return ((u128)a.d[1] << 64) | a.d[0]; }
static inline w192 w_max(void){ w192 r = {{ ~0ULL, ~0ULL, ~0ULL }}; return r; }
static inline int w_is_max(w192 a){
    return a.d[0] == ~0ULL && a.d[1] == ~0ULL && a.d[2] == ~0ULL;
}
static inline int w_cmp(w192 a, w192 b){
    for (int i = 2; i >= 0; i--){
        if (a.d[i] < b.d[i]) return -1;
        if (a.d[i] > b.d[i]) return 1;
    }
    return 0;
}
static inline int w_is_zero(w192 a){ return !(a.d[0] | a.d[1] | a.d[2]); }
static inline w192 w_add(w192 a, w192 b){
    w192 r; unsigned c = 0;
    for (int i = 0; i < 3; i++){
        u64 s = a.d[i] + b.d[i];
        u64 s2 = s + c;
        c = (s < a.d[i]) || (s2 < s);
        r.d[i] = s2;
    }
    return r;
}
static inline w192 w_addu(w192 a, u64 b){
    w192 r = a;
    r.d[0] += b;
    if (r.d[0] < b){ if (++r.d[1] == 0) ++r.d[2]; }
    return r;
}
static inline w192 w_sub(w192 a, w192 b){   /* requires a >= b */
    w192 r; unsigned br = 0;
    for (int i = 0; i < 3; i++){
        u64 s = a.d[i] - b.d[i];
        u64 s2 = s - br;
        br = (a.d[i] < b.d[i]) || (s < br);
        r.d[i] = s2;
    }
    return r;
}
static inline w192 w_subu(w192 a, u64 b){
    w192 r = a;
    u64 s = r.d[0] - b;
    if (r.d[0] < b){ if (r.d[1]-- == 0) r.d[2]--; }
    r.d[0] = s;
    return r;
}
/* a*m; *ov set if the true product exceeds 192 bits */
static inline w192 w_mulu_ov(w192 a, u64 m, int *ov){
    u128 t0 = (u128)a.d[0] * m;
    u128 t1 = (u128)a.d[1] * m + (u64)(t0 >> 64);
    u128 t2 = (u128)a.d[2] * m + (u64)(t1 >> 64);
    if (t2 >> 64) *ov = 1;
    w192 r = {{ (u64)t0, (u64)t1, (u64)t2 }};
    return r;
}
static inline w192 w_mulu(w192 a, u64 m){ int ov = 0; return w_mulu_ov(a, m, &ov); }
static w192 w_pow_sat(u64 b, int e){
    w192 r = w_from_u128(1);
    while (e-- > 0){
        int ov = 0;
        r = w_mulu_ov(r, b, &ov);
        if (ov) return w_max();
    }
    return r;
}
/* divide *v by small d (< 2^32), return remainder */
static inline unsigned w_divmod_small(w192 *v, unsigned d){
    u64 r = 0;
    for (int i = 2; i >= 0; i--){
        u64 limb = v->d[i];
        u64 a = (r << 32) | (limb >> 32);
        u64 q1 = a / d; r = a % d;
        u64 b = (r << 32) | (limb & 0xffffffffu);
        u64 q0 = b / d; r = b % d;
        v->d[i] = (q1 << 32) | q0;
    }
    return (unsigned)r;
}
/* num / den with quotient known to fit u64; remainder via *rem */
static u64 w_div_big(w192 num, w192 den, w192 *rem){
    w192 r = {{0,0,0}}; u64 q = 0;
    for (int i = 191; i >= 0; i--){
        r.d[2] = (r.d[2] << 1) | (r.d[1] >> 63);
        r.d[1] = (r.d[1] << 1) | (r.d[0] >> 63);
        r.d[0] = (r.d[0] << 1) | ((num.d[i >> 6] >> (i & 63)) & 1);
        if (w_cmp(r, den) >= 0){
            r = w_sub(r, den);
            if (i < 64) q |= 1ULL << i;
        }
    }
    if (rem) *rem = r;
    return q;
}
/* v mod d for u64 divisors up to 2^63 (bitwise; setup-time use only) */
static u64 w_mod_u64(w192 v, u64 d){
    u64 r = 0;
    for (int i = 191; i >= 0; i--){
        int hib = (int)(r >> 63);
        r = (r << 1) | ((v.d[i >> 6] >> (i & 63)) & 1);
        if (hib || r >= d) r -= d;
    }
    return r;
}

static char *w_str(w192 v, char *buf /* >= 64 bytes */){
    char *p = buf + 63; *--p = 0;
    if (w_is_zero(v)) { *--p = '0'; return p; }
    while (!w_is_zero(v)){
        unsigned r = w_divmod_small(&v, 1000000000u);
        for (int i = 0; i < 9 && !(w_is_zero(v) && r == 0); i++){
            *--p = (char)('0' + r % 10); r /= 10;
        }
    }
    return p;
}
static w192 w_parse(const char *s){
    w192 v = {{0,0,0}};
    if (!*s){ fprintf(stderr, "empty number\n"); exit(1); }
    for (; *s; s++){
        if (*s < '0' || *s > '9'){ fprintf(stderr, "bad number '%s'\n", s); exit(1); }
        v = w_addu(w_mulu(v, 10), (u64)(*s - '0'));
    }
    return v;
}

static w192 g_limitW;               /* -L, wide (g_limitw_set says whether given) */

/* ---------------- best-so-far ---------------- */

static pthread_mutex_t g_mu = PTHREAD_MUTEX_INITIALIZER;
static w192 g_best; static int g_have = 0; static unsigned g_bb1, g_bb2;

static w192 best_getw(void){
    pthread_mutex_lock(&g_mu);
    w192 r = g_have ? g_best : w_max();
    pthread_mutex_unlock(&g_mu);
    return r;
}

/* legacy u128 scanners: best clamped into u128 range (monotone-safe) */
static u128 best_get(void){
    w192 b = best_getw();
    return w_fits_u128(b) ? w_to_u128(b) : U128MAX;
}

static void report_matchw(w192 v, unsigned b1, unsigned b2){
    pthread_mutex_lock(&g_mu);
    if (!g_have || w_cmp(v, g_best) < 0){
        g_best = v; g_have = 1; g_bb1 = b1; g_bb2 = b2;
        char buf[64];
        fprintf(stderr, "\n*** MATCH: %s is an %d-digit palindrome in bases %u and %u\n",
                w_str(v, buf), g_n, b1, b2);
    }
    pthread_mutex_unlock(&g_mu);
}

static void report_match(u128 v, unsigned b1, unsigned b2){
    report_matchw(w_from_u128(v), b1, b2);
}

/* ---------------- stage state ---------------- */

typedef struct {
    unsigned b2;
    int nb1;
    unsigned b1[MAXB1];
    u128 p_lo[MAXB1], p_hi[MAXB1];   /* b1^(n-1), b1^n */
    u128 lo, hi;                     /* value window [b2^(n-1), min(max b1^n, limit)) */
    u128 mulk;                       /* b2^k */
    u64 h_start, h_end;              /* half-value range in base b2 (exclusive end) */
    int mode;                        /* 0 generic, 1 fast adjacent (odd n), 2 even-n AP, 3 MITM */
    u128 step;                       /* v increment per trailing half digit d0 */
    u64 T0, T1;                      /* t = h/b2 range for structured modes */
    int wide;                        /* stage window exceeds u128 (MITM only) */
    w192 loW, hiW;                   /* exact wide window */
    w192 p_hiW[MAXB1];               /* exact per-partner b1^n */
} Stage;

static char g_b1str[256];            /* partner list for banners/status */

static Stage g_st;
static atomic_ullong g_next;         /* absolute t (fast) or h (generic) cursor */
static atomic_ullong g_done;         /* halves processed in current stage */
static atomic_int g_cut;
static atomic_int g_stage_active = 0;
static atomic_int g_phasegen = 0;    /* bumped at every scan-phase start */
static u64 g_prev_done = 0;          /* halves from completed stages */
static volatile int g_finished = 0;

/* full palindrome test of v in base B1 (exactly-n-digits already ensured
 * by the caller's range bookkeeping): reverse low k digits, compare with
 * the leading k digits. */
static inline int fullpal(u128 v, unsigned B1){
    u128 w = v; u64 rr = 0;
    for (int i = 0; i < g_k; i++) rr = rr*B1 + udm128(&w, B1);
    if (g_odd) udm128(&w, B1);
    return w == (u128)rr;
}

/* wide variant for the MITM path (rr can exceed u64 from n ~ 28) */
static inline int fullpal_w(w192 v, unsigned B1){
    w192 w = v; u128 rr = 0;
    for (int i = 0; i < g_k; i++) rr = rr*B1 + w_divmod_small(&w, B1);
    if (g_odd) w_divmod_small(&w, B1);
    return w.d[2] == 0 && w_to_u128(w) == rr;
}

/* ---------------- fast structured scan (single b1 = b2-1) ----------------
 * Half h = t*B2 + d0.  With r = reverse of t's (m-1) digits:
 *   v(h) = t*B2*mulk + r + d0*step,  step = mulk (n odd) or mulk + mulk/B2.
 * Per candidate we maintain incrementally:
 *   bm   = v mod B1          (last digit in base B1)
 *   lead = v div B1^(n-1)    (first digit in base B1, via t2 = v - lead*p_lo)
 * and only run the full check when bm == lead. */

#define TBLOCK 4096ULL

static inline void scan_fast(const unsigned B1, const unsigned B2){
    const Stage *S = &g_st;
    const int mm1 = g_m - 1;
    const u128 mulk = S->mulk;
    const u128 step = g_odd ? mulk : mulk + mulk / B2;
    const u128 p_lo = S->p_lo[0];
    unsigned sm;
    { u128 tmp = step; sm = udm128(&tmp, B1); }
    u64 localdone = 0;

    for (;;){
        if (atomic_load_explicit(&g_cut, memory_order_relaxed)) break;
        u64 t0 = atomic_fetch_add(&g_next, TBLOCK);
        if (t0 >= S->T1) break;
        u64 t1 = t0 + TBLOCK; if (t1 > S->T1) t1 = S->T1;
        u128 stop = best_get(); if (S->hi < stop) stop = S->hi;

        for (u64 t = t0; t < t1; t++){
            u64 x = t, r = 0;
            for (int i = 0; i < mm1; i++){ r = r*B2 + x % B2; x /= B2; }
            u128 v = (u128)(t * B2) * mulk + r;
            if (v >= stop){ atomic_store(&g_cut, 1); goto out; }
            u128 tmp = v; unsigned bm = udm128(&tmp, B1);
            u128 t2 = v; unsigned lead = 0;
            while (t2 >= p_lo){ t2 -= p_lo; lead++; }
            u64 d0max = B2;
            if (t*B2 + B2 > S->h_end) d0max = S->h_end - t*B2;

            for (u64 d0 = 0; d0 < d0max; d0++){
                if (v >= stop){
                    atomic_store(&g_cut, 1); localdone += d0; goto out;
                }
                if (bm == lead && fullpal(v, B1)){
                    report_match(v, B1, B2);
                    stop = best_get(); if (S->hi < stop) stop = S->hi;
                }
                v += step;
                bm += sm; if (bm >= B1) bm -= B1;
                t2 += step;
                while (t2 >= p_lo){ t2 -= p_lo; lead++; }
            }
            localdone += d0max;
        }
        atomic_fetch_add(&g_done, localdone); localdone = 0;
    }
out:
    atomic_fetch_add(&g_done, localdone);
}

/* ---------------- generic scan (small n / multiple b1) ---------------- */

#define HBLOCK 65536ULL

static void scan_gen(void){
    const Stage *S = &g_st;
    const unsigned B2 = S->b2;
    u64 localdone = 0;
    for (;;){
        if (atomic_load(&g_cut)) break;
        u64 h0 = atomic_fetch_add(&g_next, HBLOCK);
        if (h0 >= S->h_end) break;
        u64 h1 = h0 + HBLOCK; if (h1 > S->h_end) h1 = S->h_end;
        u128 stopall = best_get(); if (S->hi < stopall) stopall = S->hi;
        for (u64 h = h0; h < h1; h++){
            u64 x = g_odd ? h / B2 : h;
            u64 r = 0;
            for (int i = 0; i < g_k; i++){ r = r*B2 + x % B2; x /= B2; }
            u128 v = (u128)h * S->mulk + r;
            if (v >= stopall){
                atomic_store(&g_cut, 1); localdone += (h - h0); goto out;
            }
            for (int i = 0; i < S->nb1; i++){
                if (v >= S->p_hi[i]) continue;
                unsigned B1 = S->b1[i];
                u128 tmp = v; unsigned last = udm128(&tmp, B1);
                u128 t2 = v; unsigned lead = 0;
                while (t2 >= S->p_lo[i]){ t2 -= S->p_lo[i]; lead++; }
                if (last == lead && fullpal(v, B1)) report_match(v, B1, B2);
            }
        }
        localdone += (h1 - h0);
        atomic_fetch_add(&g_done, localdone); localdone = 0;
    }
out:
    atomic_fetch_add(&g_done, localdone);
}

/* ---------------- even-n scanner with divisibility skipping ----------------
 * Even-length palindromes in base b are divisible by b+1 (pair mirrored
 * digits under b = -1).  So a candidate v = base(t) + d0*step can be a
 * base-b1 palindrome only if v = 0 mod q, q = b1+1.  For each row t and
 * partner b1, the valid d0 form an arithmetic progression with common
 * difference q/gcd(step,q) — visit only those (a ~q/g-fold reduction).
 * Partners with b1 = b2-1 never appear (excluded in stage_setup). */

static void scan_even(void){
    const Stage *S = &g_st;
    const unsigned B2 = S->b2;
    const int mm1 = g_m - 1;
    const u128 mulk = S->mulk;
    const u128 step = S->step;
    unsigned q[MAXB1], gg[MAXB1], stride[MAXB1], inv[MAXB1];
    u128 vstep[MAXB1];
    for (int i = 0; i < S->nb1; i++){
        q[i] = S->b1[i] + 1;
        u128 tmp = step; unsigned s = udm128(&tmp, q[i]);
        gg[i] = s ? gcd_u(s, q[i]) : q[i];
        stride[i] = q[i] / gg[i];
        inv[i] = modinv_u((s / gg[i]) % stride[i], stride[i]);
        vstep[i] = (u128)stride[i] * step;
    }
    u64 localdone = 0;
    for (;;){
        if (atomic_load_explicit(&g_cut, memory_order_relaxed)) break;
        u64 t0 = atomic_fetch_add(&g_next, TBLOCK);
        if (t0 >= S->T1) break;
        u64 t1 = t0 + TBLOCK; if (t1 > S->T1) t1 = S->T1;
        u128 stop = best_get(); if (S->hi < stop) stop = S->hi;
        for (u64 t = t0; t < t1; t++){
            u64 x = t, r = 0;
            for (int i = 0; i < mm1; i++){ r = r*B2 + x % B2; x /= B2; }
            u128 base = (u128)(t * B2) * mulk + r;
            if (base >= stop){ atomic_store(&g_cut, 1); goto out; }
            u64 d0max = B2;
            if (t*B2 + B2 > S->h_end) d0max = S->h_end - t*B2;
            for (int i = 0; i < S->nb1; i++){
                u128 tmp = base; unsigned bq = udm128(&tmp, q[i]);
                unsigned need = bq ? q[i] - bq : 0;   /* -base mod q */
                if (need % gg[i]) continue;           /* no solutions */
                u64 d0 = ((u64)(need / gg[i]) * inv[i]) % stride[i];
                const unsigned B1 = S->b1[i];
                u128 lim = S->p_hi[i]; if (stop < lim) lim = stop;
                u128 v = base + (u128)d0 * step;
                for (; d0 < d0max; d0 += stride[i], v += vstep[i]){
                    if (v >= lim) break;
                    u128 tv = v; unsigned last = udm128(&tv, B1);
                    u128 t2 = v; unsigned lead = 0;
                    while (t2 >= S->p_lo[i]){ t2 -= S->p_lo[i]; lead++; }
                    if (last == lead && fullpal(v, B1)){
                        report_match(v, B1, B2);
                        stop = best_get(); if (S->hi < stop) stop = S->hi;
                        if (stop < lim) lim = stop;
                    }
                }
            }
            localdone += d0max;
        }
        atomic_fetch_add(&g_done, localdone); localdone = 0;
    }
out:
    atomic_fetch_add(&g_done, localdone);
}

static void *even_entry(void *a){ (void)a; scan_even(); return NULL; }

/* ---------------- meet-in-the-middle scanner (even n, mode 3) ----------------
 * For even n with m = n/2, a base-b2 palindrome is v = y*b2^m + rev(y).
 * Split the half y = yh*b2^t + yl (s high digits, t = m-s low digits):
 *
 *   v = yh*b2^(m+t) + yl*b2^m + rev_t(yl)*b2^s + rev_s(yh)
 *
 * b1-palindromicity requires v mod b1^j = rev_j(p), where p = the top j
 * base-b1 digits of v — and p is determined by yh alone (up to carry,
 * a range of at most a few values), since the yl-terms are < b2^(m+t).
 * Working mod Q = b1^j, the condition splits:
 *
 *   B(yl) := (yl*b2^m + rev_t(yl)*b2^s) mod Q  ==  rev_j(p) - A(yh)  (mod Q)
 *   A(yh) := (yh*b2^(m+t) + rev_s(yh)) mod Q
 *
 * So: bucket all yl by B once (b2^t entries ~ W^(1/4)), then sweep yh in
 * increasing order (increasing v, preserving first-match semantics),
 * looking up matching yl in O(1) and full-checking the few survivors.
 * Total cost ~ W^(1/4) instead of the W^(1/2) of half-enumeration.
 * Every doubly palindromic v is found: its (yh, yl) split is unique and
 * its true prefix p lies in the scanned carry range.
 *
 * Odd n works the same way with m = k+1 half digits (middle digit not
 * mirrored): v = y*b2^k + rev_k(y div b2), and with y = yh*b2^t + yl,
 *   v = yh*b2^(k+t) + yl*b2^k + rev_(t-1)(yl div b2)*b2^s + rev_s(yh)
 * — the same separation with shift exponent k instead of m and the low
 * reverse taken over t-1 digits of yl div b2. */

typedef struct {
    unsigned b1, b2;
    int s, t, j;
    int je;                /* j+e: extended prefix depth for in-bucket filtering */
    u64 Q, QF;             /* b1^j (bucket count), b1^je (full filter modulus) */
    u64 c1f, c2f, c3f;     /* b2^K, b2^s, b2^(K+t) mod QF */
    w192 PM, PS, PMT;      /* b2^K, b2^s, b2^(K+t) (K = m even / k odd) */
    w192 PB1NJE;           /* b1^(n-je): granularity of the extended prefix */
    w192 hi_p;             /* min(b1^n, stage hi) */
    u64 yh_start, yh_end;
    u64 *off;              /* bucket ENDS: bucket k = [k?off[k-1]:0, off[k]) */
    u64 *ent;              /* per entry: yl<<30 | (B_full/Q), sorted by B_full%Q */
} Mitm;
static Mitm g_mm;

static inline u64 revdig(u64 x, unsigned b, int d){
    u64 r = 0;
    for (int i = 0; i < d; i++){ r = r*b + x % b; x /= b; }
    return r;
}

static inline u64 mitm_bfull(const Mitm *M, u64 yl){
    u64 r = g_odd ? revdig(yl / M->b2, M->b2, M->t - 1)
                  : revdig(yl, M->b2, M->t);
    return (u64)(((u128)yl * M->c1f + (u128)r * M->c2f) % M->QF);
}

/* memory-lean parallel build: counts and fill cursors live in off[]
 * itself; after the fill pass off[k] holds bucket k's END, so bucket k
 * is read as [k ? off[k-1] : 0, off[k]). Entries pack yl<<30 | bhi
 * (t-selection keeps yl < 2^34, je-selection keeps bhi < 2^30). */
typedef struct { const Mitm *M; u64 lo, hi; int pass; } MBArg;
static void *mitm_build_worker(void *ap){
    MBArg *A = (MBArg*)ap; const Mitm *M = A->M;
    for (u64 yl = A->lo; yl < A->hi; yl++){
        u64 Bf = mitm_bfull(M, yl);
        u64 key = Bf % M->Q;
        if (A->pass == 1){
            __atomic_fetch_add(&M->off[key], 1, __ATOMIC_RELAXED);
        } else {
            u64 slot = __atomic_fetch_add(&M->off[key], 1, __ATOMIC_RELAXED);
            M->ent[slot] = (yl << 30) | (Bf / M->Q);
        }
    }
    return NULL;
}

static int entcmp(const void *a, const void *b){
    u64 x = *(const u64*)a & 0x3fffffffu, y = *(const u64*)b & 0x3fffffffu;
    return x < y ? -1 : x > y ? 1 : 0;
}
#define SORT_CHUNK 65536ULL
static void *mitm_sort_worker(void *ap){
    MBArg *A = (MBArg*)ap; const Mitm *M = A->M;
    _Atomic(u64) *cur = (_Atomic(u64)*)(uintptr_t)A->lo;
    for (;;){
        u64 q0 = atomic_fetch_add(cur, SORT_CHUNK);
        if (q0 >= M->Q) break;
        u64 q1 = q0 + SORT_CHUNK; if (q1 > M->Q) q1 = M->Q;
        for (u64 q = q0; q < q1; q++){
            u64 s = q ? M->off[q-1] : 0, e = M->off[q];
            if (e > s + 1) qsort(M->ent + s, e - s, 8, entcmp);
        }
    }
    return NULL;
}

static int mitm_build(Mitm *M){
    u64 Nyl = pw64(M->b2, M->t);
    M->off = calloc(M->Q + 1, sizeof *M->off);
    M->ent = malloc(Nyl * sizeof *M->ent);
    if (!M->off || !M->ent){
        fprintf(stderr, "MITM: table allocation failed (%llu entries)\n",
                (unsigned long long)Nyl);
        free(M->off); free(M->ent);
        return -1;
    }
    int nth = g_nthreads; if (nth > 64) nth = 64;
    pthread_t tids[64]; MBArg args[64];
    for (int pass = 1; pass <= 2; pass++){
        if (pass == 2){
            u64 run = 0;
            for (u64 q = 0; q < M->Q; q++){ u64 c = M->off[q]; M->off[q] = run; run += c; }
            M->off[M->Q] = run;
        }
        u64 chunk = (Nyl + nth - 1) / nth;
        for (int i = 0; i < nth; i++){
            args[i].M = M; args[i].lo = (u64)i * chunk;
            args[i].hi = (u64)(i + 1) * chunk > Nyl ? Nyl : (u64)(i + 1) * chunk;
            args[i].pass = pass;
            pthread_create(&tids[i], NULL, mitm_build_worker, &args[i]);
        }
        for (int i = 0; i < nth; i++) pthread_join(tids[i], NULL);
    }
    /* pass 3: sort each bucket by its low-30-bit residue so the sweep can
     * binary-search instead of scanning (defuses degenerate mega-buckets) */
    {
        _Atomic(u64) qcursor = 0;
        MBArg sargs[64];
        void *(*sfn)(void*) = mitm_sort_worker;
        for (int i = 0; i < nth; i++){
            sargs[i].M = M; sargs[i].lo = (u64)(uintptr_t)&qcursor; sargs[i].hi = 0; sargs[i].pass = 3;
            pthread_create(&tids[i], NULL, sfn, &sargs[i]);
        }
        for (int i = 0; i < nth; i++) pthread_join(tids[i], NULL);
    }
    return 0;
}

static void *mitm_worker(void *a){
    (void)a;
    const Mitm *M = &g_mm;
    const unsigned B1 = M->b1, B2 = M->b2;
    u64 localdone = 0;
    for (;;){
        if (atomic_load_explicit(&g_cut, memory_order_relaxed)) break;
        u64 y0 = atomic_fetch_add(&g_next, TBLOCK);
        if (y0 >= M->yh_end) break;
        u64 y1 = y0 + TBLOCK; if (y1 > M->yh_end) y1 = M->yh_end;
        w192 stop = best_getw();
        if (w_cmp(M->hi_p, stop) < 0) stop = M->hi_p;
        /* extended prefix p (je base-b1 digits) and its residual,
         * maintained incrementally across yh */
        w192 vbase = w_mulu(M->PMT, y0);
        w192 prem;
        u64 p = w_div_big(vbase, M->PB1NJE, &prem);
        /* yh*c3f mod QF: seeded once per block, then incremental */
        u64 term3 = (u64)(((u128)(y0 % M->QF) * M->c3f) % M->QF);
        for (u64 yh = y0; yh < y1; yh++){
            if (w_cmp(vbase, stop) >= 0){ atomic_store(&g_cut, 1); goto out; }
            u64 rs = revdig(yh, B2, M->s);
            u64 Af = term3 + rs % M->QF;
            if (Af >= M->QF) Af -= M->QF;
            w192 vtop = w_subu(w_add(vbase, M->PMT), 1);
            if (w_cmp(vtop, stop) >= 0) vtop = w_subu(stop, 1);
            for (u64 pp = p;; pp++){
                w192 plo = w_mulu(M->PB1NJE, pp);
                if (w_cmp(plo, vtop) > 0) break;
                w192 phi = w_add(plo, M->PB1NJE);
                u64 rjf = revdig(pp, B1, M->je);
                u64 R = rjf >= Af ? rjf - Af : rjf + M->QF - Af;
                u64 key = R % M->Q;
                u64 needhi = R / M->Q;
                u64 ib0 = key ? M->off[key-1] : 0, ib1 = M->off[key];
                { u64 blo = ib0, bhi2 = ib1;   /* buckets sorted: find needhi run */
                  while (blo < bhi2){ u64 mid = (blo + bhi2) >> 1;
                      if ((M->ent[mid] & 0x3fffffffu) < needhi) blo = mid + 1; else bhi2 = mid; }
                  ib0 = blo; }
                for (u64 idx = ib0; idx < ib1; idx++){
                    u64 en = M->ent[idx];
                    if ((en & 0x3fffffffu) != needhi) break;   /* sorted: past the run */
                    u64 yl = en >> 30;
                    u64 rt = g_odd ? revdig(yl / B2, B2, M->t - 1)
                                   : revdig(yl, B2, M->t);
                    w192 v = w_addu(w_add(w_add(vbase, w_mulu(M->PM, yl)),
                                          w_mulu(M->PS, rt)), rs);
                    if (w_cmp(v, plo) < 0 || w_cmp(v, phi) >= 0) continue;  /* prefix mismatch */
                    if (w_cmp(v, stop) >= 0) continue;
                    if (fullpal_w(v, B1)){
                        report_matchw(v, B1, B2);
                        stop = best_getw();
                        if (w_cmp(M->hi_p, stop) < 0) stop = M->hi_p;
                    }
                }
            }
            localdone++;
            vbase = w_add(vbase, M->PMT);
            prem = w_add(prem, M->PMT);
            while (w_cmp(prem, M->PB1NJE) >= 0){ prem = w_sub(prem, M->PB1NJE); p++; }
            term3 += M->c3f; if (term3 >= M->QF) term3 -= M->QF;
        }
        atomic_fetch_add(&g_done, localdone); localdone = 0;
    }
out:
    atomic_fetch_add(&g_done, localdone);
    return NULL;
}

static u64 g_mitm_skip = 0;         /* -H in MITM mode: yh done-counter to resume at */
static unsigned g_mitm_skip_b1 = 0; /* -P: partner phase the counter belongs to */

/* returns 0 ok / -1 alloc failure; scans every partner of the stage */
static int run_mitm_stage(Stage *S){
    const int m = g_m;
    if (g_mitm_skip && g_mitm_skip_b1){
        int found = 0;
        for (int i = 0; i < S->nb1; i++) if (S->b1[i] == g_mitm_skip_b1) found = 1;
        if (!found){
            fprintf(stderr, "warning: -P %u is not a partner of stage b2=%u; ignoring -H\n",
                    g_mitm_skip_b1, S->b2);
            g_mitm_skip = 0;
        }
    }
    for (int i = 0; i < S->nb1; i++){
        Mitm *M = &g_mm;
        memset(M, 0, sizeof *M);
        M->b1 = S->b1[i]; M->b2 = S->b2;
        if (g_mitm_skip){
            unsigned target = g_mitm_skip_b1 ? g_mitm_skip_b1 : S->b1[0];
            if (M->b1 < target){
                fprintf(stderr, "  MITM b1=%u: skipped (resume asserts this partner completed)\n", M->b1);
                continue;
            }
        }
        /* pick t: largest with b2^t <= tablecap, capped at m/2 */
        int t = 0; u64 e = 1;
        while (t + 1 <= m/2 && e * M->b2 <= g_tablecap && e * M->b2 < ((u64)1 << 34)){ e *= M->b2; t++; }
        M->t = t; M->s = m - t;
        /* pick j: largest with b1^j <= b2^t (keeps buckets ~1 and the
         * prefix carry range small) */
        int j = 0; u64 q = 1;
        while (q * M->b1 <= e){ q *= M->b1; j++; }
        if (j < 1){ j = 1; q = M->b1; }
        M->j = j; M->Q = q;
        int K = g_odd ? g_k : m;          /* shift exponent: k (odd) / m (even) */
        M->PM = w_pow_sat(M->b2, K);
        M->PS = w_pow_sat(M->b2, M->s);
        M->PMT = w_pow_sat(M->b2, K + t);
        if (w_is_max(M->PMT)){
            fprintf(stderr, "  MITM b1=%u: values exceed w192 range, skipping\n", M->b1);
            continue;
        }
        /* extend the prefix depth j -> je while: the mirror constraint
         * still applies (je <= k), the per-yh prefix carry range stays
         * small (PMT / b1^(n-je) <= ~8), the packed high residue fits
         * u32, and QF fits comfortably in u64 */
        {
            int je = j; u64 qf = M->Q;
            while (je + 1 <= g_k && je + 1 <= g_n - 1 &&
                   (u128)qf * M->b1 <= ((u64)1 << 40) &&
                   (qf / M->Q) * (u64)M->b1 <= ((u64)1 << 30)){
                w192 pnje = w_pow_sat(M->b1, g_n - (je + 1));
                if (w_is_max(pnje)) break;
                if (w_div_big(M->PMT, pnje, NULL) > 8) break;
                je++; qf *= M->b1;
            }
            M->je = je; M->QF = qf;
        }
        M->PB1NJE = w_pow_sat(M->b1, g_n - M->je);
        if (w_is_max(M->PB1NJE)){
            fprintf(stderr, "  MITM b1=%u: values exceed w192 range, skipping\n", M->b1);
            continue;
        }
        M->hi_p = S->p_hiW[i];
        if (w_cmp(S->hiW, M->hi_p) < 0) M->hi_p = S->hiW;
        M->c1f = w_mod_u64(M->PM,  M->QF);
        M->c2f = w_mod_u64(M->PS,  M->QF);
        M->c3f = w_mod_u64(M->PMT, M->QF);
        M->yh_start = pw64(M->b2, M->s - 1);
        M->yh_end = w_div_big(M->hi_p, M->PMT, NULL) + 1;
        u64 cap = pw64(M->b2, M->s);
        if (M->yh_end > cap) M->yh_end = cap;
        if (M->yh_end <= M->yh_start) continue;
        w192 bg = best_getw();
        if (w_cmp(w_mulu(M->PMT, M->yh_start), bg) >= 0) continue;

        fprintf(stderr, "  MITM b1=%u: s=%d t=%d j=%d je=%d Q=%llu table=%llu yh=%llu\n",
                M->b1, M->s, M->t, M->j, M->je, (unsigned long long)M->Q,
                (unsigned long long)pw64(M->b2, t),
                (unsigned long long)(M->yh_end - M->yh_start));
        if (mitm_build(M)) return -1;

        u64 skip0 = 0;
        if (g_mitm_skip){
            /* same in-flight-block margin logic as the legacy -H */
            u64 margin = (u64)g_nthreads * TBLOCK;
            skip0 = g_mitm_skip > margin ? g_mitm_skip - margin : 0;
            u64 span = M->yh_end - M->yh_start;
            if (skip0 > span) skip0 = span;
            fprintf(stderr, "  resuming partner b1=%u at +%llu of %llu yh (margin rescan %llu)\n",
                    M->b1, (unsigned long long)skip0,
                    (unsigned long long)(M->yh_end - M->yh_start),
                    (unsigned long long)(g_mitm_skip > skip0 ? g_mitm_skip - skip0 : 0));
            g_mitm_skip = 0;
        }
        S->h_start = M->yh_start; S->h_end = M->yh_end;  /* for status */
        atomic_store(&g_next, M->yh_start + skip0);
        atomic_store(&g_done, skip0);
        atomic_store(&g_cut, 0);
        atomic_fetch_add(&g_phasegen, 1);
        pthread_t tids[256];
        for (int w = 0; w < g_nthreads; w++)
            pthread_create(&tids[w], NULL, mitm_worker, NULL);
        for (int w = 0; w < g_nthreads; w++) pthread_join(tids[w], NULL);
        g_prev_done += atomic_load(&g_done);
        free(M->off); free(M->ent);
        M->off = NULL; M->ent = NULL;
    }
    return 0;
}

/* specialized entry points so base divisions become multiply-by-magic */
#define DEFFAST(B) static void *fast_##B(void *a){ (void)a; scan_fast(B, (B)+1); return NULL; }
DEFFAST(3)  DEFFAST(4)  DEFFAST(5)  DEFFAST(6)  DEFFAST(7)  DEFFAST(8)
DEFFAST(9)  DEFFAST(10) DEFFAST(11) DEFFAST(12) DEFFAST(13) DEFFAST(14)
DEFFAST(15) DEFFAST(16) DEFFAST(17) DEFFAST(18) DEFFAST(19) DEFFAST(20)
DEFFAST(21) DEFFAST(22) DEFFAST(23) DEFFAST(24) DEFFAST(25) DEFFAST(26)
DEFFAST(27) DEFFAST(28) DEFFAST(29) DEFFAST(30)

static void *gen_entry(void *a){ (void)a; scan_gen(); return NULL; }

static void *(*fast_tab(unsigned b1))(void *){
    switch (b1){
    case 3: return fast_3;   case 4: return fast_4;   case 5: return fast_5;
    case 6: return fast_6;   case 7: return fast_7;   case 8: return fast_8;
    case 9: return fast_9;   case 10: return fast_10; case 11: return fast_11;
    case 12: return fast_12; case 13: return fast_13; case 14: return fast_14;
    case 15: return fast_15; case 16: return fast_16; case 17: return fast_17;
    case 18: return fast_18; case 19: return fast_19; case 20: return fast_20;
    case 21: return fast_21; case 22: return fast_22; case 23: return fast_23;
    case 24: return fast_24; case 25: return fast_25; case 26: return fast_26;
    case 27: return fast_27; case 28: return fast_28; case 29: return fast_29;
    case 30: return fast_30;
    }
    return NULL;
}

/* ---------------- stage setup ---------------- */

/* returns 0 = run it, 1 = skip (no valid b1), -1 = past bound: stop */
static int stage_setup(unsigned b2, Stage *S, w192 bound){
    w192 loW = w_pow_sat(b2, g_n - 1);
    if (w_is_max(loW) || w_cmp(loW, bound) >= 0) return -1;
    memset(S, 0, sizeof *S);
    S->b2 = b2;
    w192 himaxW = {{0,0,0}};
    for (unsigned b1 = 2; b1 < b2; b1++){
        /* even n: an even-length palindrome in base b1 is divisible by
         * b1+1; if b1+1 == b2, then b2 | v, contradicting v's nonzero
         * trailing digit in base b2 (= its leading digit). Impossible. */
        if (!(g_n & 1) && b1 == b2 - 1) continue;
        w192 h = w_pow_sat(b1, g_n);
        if (w_cmp(h, loW) > 0){
            if (S->nb1 >= MAXB1){
                fprintf(stderr, "warning: >%d partner bases at b2=%u; truncating\n", MAXB1, b2);
                break;
            }
            S->b1[S->nb1]   = b1;
            S->p_hiW[S->nb1] = h;
            S->p_lo[S->nb1] = pow_sat(b1, g_n - 1);
            S->p_hi[S->nb1] = w_fits_u128(h) ? w_to_u128(h) : U128MAX;
            S->nb1++;
            if (w_cmp(h, himaxW) > 0) himaxW = h;
        }
    }
    if (!S->nb1) return 1;
    w192 hiW = himaxW;
    if (g_limitw_set && w_cmp(g_limitW, hiW) < 0) hiW = g_limitW;
    S->loW = loW; S->hiW = hiW;
    S->wide = !w_fits_u128(hiW);
    S->lo = w_fits_u128(loW) ? w_to_u128(loW) : U128MAX;
    S->hi = w_fits_u128(hiW) ? w_to_u128(hiW) : U128MAX;
    S->mulk = pow_sat(b2, g_k);
    S->step = g_odd ? S->mulk : S->mulk + S->mulk / b2;
    if (g_mitm && g_n >= 8)
        S->mode = 3;
    else if (S->wide){
        fprintf(stderr, "stage b2=%u: window exceeds u128 — rerun with -M (MITM); skipping\n", b2);
        return 1;
    }
    else if (!(g_n & 1) && g_n >= 6)
        S->mode = 2;
    else if (S->nb1 == 1 && S->b1[0] == b2 - 1 && g_n >= 6 &&
             fast_tab(S->b1[0]) != NULL)
        S->mode = 1;
    else
        S->mode = 0;
    if (!S->wide){
        u128 hs128 = pow_sat(b2, g_m - 1);
        u128 hecap128 = pow_sat(b2, g_m);
        u128 he128 = S->hi / S->mulk + 1;
        if (he128 > hecap128) he128 = hecap128;
        if (hs128 >= ((u128)1 << 63) || he128 >= ((u128)1 << 63)){
            if (S->mode != 3){
                fprintf(stderr, "stage b2=%u: half range exceeds 2^63, cannot scan\n", b2);
                return 1;
            }
        } else {
            u64 hs = (u64)hs128, he = (u64)he128;
            if (he <= hs) return 1;
            S->h_start = hs; S->h_end = he;
        }
    }
    if (S->mode == 1 || S->mode == 2){
        if (!S->h_end) return 1;
        S->T0 = pw64(b2, g_m - 2);
        S->T1 = (S->h_end + b2 - 1) / b2;
    }
    {
        size_t off = 0; g_b1str[0] = 0;
        for (int i = 0; i < S->nb1; i++)
            off += (size_t)snprintf(g_b1str + off, sizeof g_b1str - off,
                                    "%s%u", i ? "," : "", S->b1[i]);
    }
    return 0;
}

/* half-count of stage b2 clamped to bound, for overall ETA */
static u64 stage_total_bounded(unsigned b2, u128 bound){
    u128 lo = pow_sat(b2, g_n - 1);
    if (lo >= bound) return 0;
    u128 himax = 0;
    for (unsigned b1 = 2; b1 < b2; b1++){
        if (!(g_n & 1) && b1 == b2 - 1) continue;
        u128 h = pow_sat(b1, g_n);
        if (h > lo && h > himax) himax = h;
    }
    if (!himax) return 0;
    if (bound < himax) himax = bound;
    u128 mulk = pow_sat(b2, g_k);
    u64 hs = pw64(b2, g_m - 1);
    u64 he = (u64)(himax / mulk) + 1;
    u64 cap = pw64(b2, g_m);
    if (he > cap) he = cap;
    return he > hs ? he - hs : 0;
}

/* ---------------- status thread ---------------- */

static double now_s(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9*ts.tv_nsec;
}

static void fmt_eta(double s, char *out){
    if (s < 0 || s > 86400.0*365) { strcpy(out, "?"); return; }
    long t = (long)s;
    if (t >= 3600) sprintf(out, "%ld:%02ld:%02ld", t/3600, (t/60)%60, t%60);
    else sprintf(out, "%ld:%02ld", t/60, t%60);
}

static void *status_thread(void *a){
    (void)a;
    double t_prev = now_s(); u64 d_prev = 0; int prev_stage = -1;
    double rate = 0.0;
    while (!g_finished){
        struct timespec ts = {2, 0}; nanosleep(&ts, NULL);
        if (g_finished) break;
        if (!atomic_load(&g_stage_active)) { prev_stage = -1; continue; }
        u64 done = atomic_load(&g_done);
        double t = now_s();
        int st = atomic_load(&g_phasegen);
        if (st != prev_stage || done < d_prev){
            /* first tick of a scan phase (new stage, new partner, or a
             * resume that jumped the counter): seed the sample at the
             * current done count and skip one interval so the first
             * printed rate comes from a clean delta */
            d_prev = done; t_prev = t; prev_stage = st; rate = 0.0;
            continue;
        }
        double inst = (done > d_prev && t > t_prev) ? (double)(done - d_prev)/(t - t_prev) : 0.0;
        if (inst > 0) rate = rate > 0 ? 0.5*rate + 0.5*inst : inst;
        d_prev = done; t_prev = t;

        u64 stot = g_st.h_end - g_st.h_start;
        double pct = stot ? 100.0*(double)done/(double)stot : 100.0;
        char eta1[32] = "?", eta2[32] = "?", bbuf[64];
        if (rate > 0) fmt_eta((double)(stot - (done < stot ? done : stot))/rate, eta1);

        /* overall vs. bound (limit and/or best found) */
        u128 bound = best_get();
        if (g_limit && g_limit < bound) bound = g_limit;
        char opct[32] = "?";
        if (bound != U128MAX){
            u64 total = 0, before = 0;
            for (unsigned b2 = 3; b2 <= g_st.b2 + 64; b2++){
                u64 c = stage_total_bounded(b2, bound);
                if (!c && pow_sat(b2, g_n-1) >= bound) break;
                total += c;
                if (b2 < g_st.b2) before += c;
            }
            u64 odone = before + (done < stot ? done : stot);
            if (total){
                snprintf(opct, sizeof opct, "%.1f%%", 100.0*(double)odone/(double)total);
                if (rate > 0) fmt_eta((double)(total - (odone < total ? odone : total))/rate, eta2);
            }
        }
        fprintf(stderr, "[b2=%u b1={%s}] stage %5.1f%% (%llu/%llu) | %.1fM/s | stage ETA %s | overall %s ETA %s | best %s\n",
                g_st.b2, g_b1str, pct,
                (unsigned long long)done, (unsigned long long)stot,
                rate/1e6, eta1, opct, eta2,
                g_have ? w_str(g_best, bbuf) : "-");
    }
    return NULL;
}

/* ---------------- result verification / printing ---------------- */

static const char *DIGCH = "0123456789abcdefghijklmnopqrstuvwxyz";

static int digits_of(w192 v, unsigned b, unsigned *d, int maxd){
    int i = 0;
    while (!w_is_zero(v) && i < maxd){ d[i++] = w_divmod_small(&v, b); }
    return w_is_zero(v) ? i : -1;   /* little-endian */
}

static void print_all_bases(w192 v){
    unsigned d[200];
    for (unsigned b = 2; b <= 128; b++){
        w192 plo = w_pow_sat(b, g_n - 1), phi = w_pow_sat(b, g_n);
        if (w_cmp(v, plo) < 0 || w_cmp(v, phi) >= 0) continue;
        int nd = digits_of(v, b, d, 200);
        if (nd != g_n) continue;
        int pal = 1;
        for (int i = 0, j = nd-1; i < j; i++, j--) if (d[i] != d[j]) { pal = 0; break; }
        if (!pal) continue;
        printf("  base %-3u: ", b);
        for (int i = nd-1; i >= 0; i--){
            if (b <= 36) putchar(DIGCH[d[i]]);
            else printf("%u%s", d[i], i ? "." : "");
        }
        printf("\n");
    }
}

/* ---------------- main ---------------- */

int main(int argc, char **argv){
    int opt;
    unsigned b2start = 3, b2end = 0;
    u64 hskip = 0;
    while ((opt = getopt(argc, argv, "n:t:L:B:E:H:P:T:Mqh")) != -1){
        switch (opt){
        case 'n': g_n = atoi(optarg); break;
        case 't': g_nthreads = atoi(optarg); break;
        case 'L':
            g_limitW = w_parse(optarg); g_limitw_set = 1;
            g_limit = w_fits_u128(g_limitW) ? w_to_u128(g_limitW) : U128MAX;
            break;
        case 'B': b2start = (unsigned)atoi(optarg); break;
        case 'E': b2end = (unsigned)atoi(optarg); break;
        case 'H': hskip = strtoull(optarg, NULL, 10); break;
        case 'M': g_mitm = 1; break;
        case 'P': g_mitm_skip_b1 = (unsigned)atoi(optarg); break;
        case 'T': g_tablecap = strtoull(optarg, NULL, 10); break;
        case 'q': g_quiet = 1; break;
        default:
            fprintf(stderr, "usage: %s -n digits [-t threads] [-L limit] [-B b2start] [-H done] [-M] [-T entries] [-q]\n"
                            "  -M: meet-in-the-middle scanner for even n (~W^(1/4) work)\n"
                            "  -T: MITM table entry cap (default 3e8, ~2.5GB peak)\n"
                            "  -B: resume at stage b2 (asserts earlier stages are empty)\n"
                            "  -E: stop before stage b2 (fleet mode: keep off stages other machines own)\n"
                            "  -H: with -B, resume that stage from a status-line done-counter\n"
                            "      (paste the first number of the last status line; a safety\n"
                            "      margin for in-flight work is subtracted automatically)\n"
                            "  -P: (MITM, multi-partner stages) the b1 whose partner phase the\n"
                            "      -H counter belongs to; earlier partners are assumed complete\n", argv[0]);
            return 2;
        }
    }
    if (hskip && b2start == 3){
        fprintf(stderr, "-H requires -B <b2> naming the stage it applies to\n");
        return 2;
    }
    if (g_n < 1 || g_n > 36){ fprintf(stderr, "need -n in 1..36\n"); return 2; }
    if (g_n == 1){ printf("a(1) = 1  (single digit '1' in every base >= 2)\n"); return 0; }
    if (g_nthreads <= 0) g_nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN);
    if (g_nthreads < 1) g_nthreads = 1;
    if (g_nthreads > 256) g_nthreads = 256;

    g_k = g_n / 2; g_m = (g_n + 1) / 2; g_odd = g_n & 1;

    char lb[64];
    fprintf(stderr, "A171740 search: n=%d, threads=%d, limit=%s\n",
            g_n, g_nthreads,
            g_limitw_set ? w_str(g_limitW, lb) : "none");

    double t_start = now_s();
    pthread_t stat_tid;
    if (!g_quiet) pthread_create(&stat_tid, NULL, status_thread, NULL);

    for (unsigned b2 = b2start;; b2++){
        if (b2end && b2 >= b2end){
            fprintf(stderr, "reached -E %u: stopping before stage b2=%u\n", b2end, b2);
            if (!g_limitw_set){   /* report the true covered bound */
                g_limitW = w_pow_sat(b2, g_n - 1);
                g_limitw_set = 1;
            }
            break;
        }
        w192 bound = best_getw();
        if (g_limitw_set && w_cmp(g_limitW, bound) < 0) bound = g_limitW;
        int rc = stage_setup(b2, &g_st, bound);
        if (rc < 0) break;
        if (rc > 0) continue;

        {
            char a[64], b[64];
            static const char *modename[] = { " (generic path)", "", " (even-AP path)", " (MITM path)" };
            fprintf(stderr, "stage b2=%u vs b1={%s}: values [%s, %s), %llu candidates%s\n",
                    b2, g_b1str, w_str(g_st.loW, a), w_str(g_st.hiW, b),
                    (unsigned long long)(g_st.h_end - g_st.h_start),
                    modename[g_st.mode]);
        }

        if (g_st.mode == 3){
            if (hskip){
                if (b2 == b2start){ g_mitm_skip = hskip; }
                else fprintf(stderr, "warning: -H ignored (stage b2=%u from -B was not runnable)\n", b2start);
                hskip = 0;
            }
            atomic_store(&g_stage_active, 1);
            if (run_mitm_stage(&g_st)){ g_finished = 1; return 1; }
            atomic_store(&g_stage_active, 0);
            continue;
        }

        u64 skipped = 0;
        if (hskip){
            if (b2 == b2start){
                /* the done counter is a sum over completed blocks, not a
                 * contiguous frontier: at kill time up to nthreads blocks
                 * below it may have been in flight. Rescan that margin. */
                u64 blk = TBLOCK * (u64)g_st.b2;
                if (blk < HBLOCK) blk = HBLOCK;
                u64 margin = (u64)g_nthreads * blk;
                u64 eff = hskip > margin ? hskip - margin : 0;
                u64 total = g_st.h_end - g_st.h_start;
                if (eff > total) eff = total;
                if (g_st.mode){
                    u64 tskip = eff / g_st.b2;
                    atomic_store(&g_next, g_st.T0 + tskip);
                    skipped = tskip * g_st.b2;
                } else {
                    atomic_store(&g_next, g_st.h_start + eff);
                    skipped = eff;
                }
                fprintf(stderr, "resuming stage b2=%u at +%llu of %llu halves (margin rescan %llu)\n",
                        b2, (unsigned long long)skipped, (unsigned long long)total,
                        (unsigned long long)(hskip > skipped ? hskip - skipped : 0));
            } else {
                fprintf(stderr, "warning: -H ignored (stage b2=%u from -B was not runnable)\n", b2start);
            }
            hskip = 0;
        }
        if (!skipped)
            atomic_store(&g_next, g_st.mode ? g_st.T0 : g_st.h_start);
        atomic_store(&g_done, skipped);
        atomic_store(&g_cut, 0);
        atomic_fetch_add(&g_phasegen, 1);
        atomic_store(&g_stage_active, 1);

        void *(*fn)(void *) =
            g_st.mode == 1 ? fast_tab(g_st.b1[0]) :
            g_st.mode == 2 ? even_entry : gen_entry;
        pthread_t tids[256];
        for (int i = 0; i < g_nthreads; i++) pthread_create(&tids[i], NULL, fn, NULL);
        for (int i = 0; i < g_nthreads; i++) pthread_join(tids[i], NULL);

        atomic_store(&g_stage_active, 0);
        g_prev_done += atomic_load(&g_done);
    }

    g_finished = 1;
    if (!g_quiet) pthread_join(stat_tid, NULL);

    double el = now_s() - t_start;
    fprintf(stderr, "done: %.2fs, %llu candidates tested\n",
            el, (unsigned long long)(g_prev_done));

    if (g_have){
        char buf[64];
        printf("a(%d) = %s   (palindromic with %d digits in bases %u and %u)\n",
               g_n, w_str(g_best, buf), g_n, g_bb1, g_bb2);
        print_all_bases(g_best);
        return 0;
    }
    char b[64];
    printf("no doubly %d-digit palindrome found below %s\n",
           g_n, g_limitw_set ? w_str(g_limitW, b) : "(unbounded?!)");
    return 1;
}
