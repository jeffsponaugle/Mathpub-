/*
 * enantio.c -- count space-group types per geometric crystal class (Q-class),
 *              both up to affine conjugacy (A004029) and up to orientation-
 *              preserving affine conjugacy (A006227); the difference is the
 *              number of enantiomorphic pairs (A395859).
 *
 * Input: one or more CARAT bravais_TYP files, each a representative of a
 * Q-class (e.g. the files group.N / min.N / max.N of CARAT's Q-catalog).
 *
 * For each Q-class:
 *   1. CARAT's q2z() splits it into Z-classes (arithmetic classes) and
 *      returns, for each representative G <= GL_n(Z), generators of the
 *      normalizer N = N_{GL_n(Z)}(G).
 *   2. A presentation of G and the cohomology group
 *      H = H^1(G, Q^n/Z^n) = Z/d_1 x ... x Z/d_m are computed with CARAT
 *      (Zassenhaus algorithm); normalop() gives the action of each
 *      normalizer generator on H as an m x m matrix.
 *   3. Space-group types in this Z-class  <->  N-orbits on H.
 *      Orientation-preserving types       <->  N+-orbits on H,
 *      N+ = N cap SL_n(Z).  We count N-orbits on H x {+1,-1}, where
 *      n.(x,s) = (n.x, det(n) s); these are in bijection with N+-orbits on
 *      H when N != N+, and number 2*(N-orbits) when N = N+ (then the whole
 *      arithmetic class is chiral).  Either way
 *          proper types  = #orbits on H x {+-1}
 *          pairs         = proper - affine.
 *   4. Orbits are counted either explicitly (BFS over H x {+-1}) when |H| is
 *      small, or with Burnside's lemma over the image of N in
 *      Aut(H) x {+-1}:
 *          #orbits(H)        = (1/|I|) sum_{(A,s) in I} |Fix(A)|
 *          #orbits(H x +-1)  = (1/|I|) sum_{(A,+1) in I} 2 |Fix(A)|
 *      where |Fix(A)| = |ker(A - 1)| is the product of the elementary
 *      divisors of [A - 1 | diag(d)].
 *
 * Output (stdout), one tab-separated line per input file:
 *   file  Zclasses  affine  proper  pairs  improper_pointgroup  badZ  seconds
 * badZ = number of Z-classes where some normalizer generator's computed
 * action on H^1 is not an automorphism (counts there are unreliable).
 *
 * Options (letters chosen to avoid those read inside CARAT's library):
 *   -x=k   explicit enumeration when |H| <= 2^k        (default 22)
 *   -y=n   cap on |image of N| for Burnside            (default 3000000)
 *   -z     cross-check: run both methods when feasible and compare
 *   -m     minimal: skip Q-classes whose point group contains det -1
 *          elements (they contribute no enantiomorphs); prints affine = -1
 *   -w     write one line per Z-class to stderr
 *   -j     inputs are Z-class representatives (with normalizer), not Q-classes
 *   -e     dump H^1 and the normalizer action matrices to stderr
 *   -k     only split into Z-classes and write them to $ENANTIO_DUMP_DIR
 *   -b     only print CARAT's H^1(G,Q^n/Z^n) per Z-class (line "H <file> <Z> <d1xd2..>")
 */

#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <gmp.h>
#include <zass.h>
#include <longtools.h>
#include <presentation.h>
#include <base.h>
#include <bravais.h>
#include <ZZ.h>
#include <datei.h>
#include <tools.h>

#include <stdint.h>
#include <sys/time.h>

typedef unsigned __int128 u128;

static int EXPLICIT_LOG2 = 22;
static long MAX_IMAGE = 3000000;
static int CROSSCHECK = 0;
static int MINIMAL = 0;
static int ZLINES = 0;
static int DUMP = 0;
static int ZFILES = 0;
static int DUMPONLY = 0;
static int COHOONLY = 0;

static double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

static void die(const char *msg, const char *file)
{
    fprintf(stderr, "enantio: FATAL: %s (%s)\n", msg, file ? file : "");
    exit(2);
}

static char *u128str(u128 x, char *buf)
{
    char tmp[64];
    int n = 0;
    if (x == 0) tmp[n++] = '0';
    while (x) { tmp[n++] = '0' + (int)(x % 10); x /= 10; }
    for (int i = 0; i < n; i++) buf[i] = tmp[n - 1 - i];
    buf[n] = 0;
    return buf;
}

/* determinant of a small integer matrix (Bareiss, exact) */
static long long det_ll(matrix_TYP *A)
{
    int n = A->rows;
    __int128 M[16][16], prev = 1;
    int sign = 1;
    if (n > 16 || A->cols != n) die("det_ll: bad matrix", NULL);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) M[i][j] = A->array.SZ[i][j];
    for (int k = 0; k < n - 1; k++) {
        if (M[k][k] == 0) {
            int p = -1;
            for (int i = k + 1; i < n; i++) if (M[i][k] != 0) { p = i; break; }
            if (p < 0) return 0;
            for (int j = 0; j < n; j++) { __int128 t = M[k][j]; M[k][j] = M[p][j]; M[p][j] = t; }
            sign = -sign;
        }
        for (int i = k + 1; i < n; i++)
            for (int j = k + 1; j < n; j++)
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev;
        prev = M[k][k];
    }
    return sign * (long long)M[n - 1][n - 1];
}

static int det_sign(matrix_TYP *A)
{
    long long d = det_ll(A);
    if (d != 1 && d != -1) die("matrix is not unimodular", NULL);
    return (int)d;
}

/* ------------------------------------------------------------------------ */
/* the cohomology group H = Z/d[0] x ... x Z/d[m-1], d[i] | d[i+1]           */

#define MAXM 96

typedef struct {
    int m;
    int d[MAXM];
    int L;              /* exponent = d[m-1] */
    double log2h;       /* log2 |H| */
} hgrp;

/* generators of the image of N in Aut(H) x {+-1} */
typedef struct {
    int n;
    int cap;
    int **A;            /* m*m row-major, entries in [0,d[row]) */
    int *s;             /* determinant sign */
} genset;

static int gs_has(const genset *g, const int *A, int s, int mm)
{
    for (int i = 0; i < g->n; i++)
        if (g->s[i] == s && memcmp(g->A[i], A, sizeof(int) * mm) == 0) return 1;
    return 0;
}

static void gs_add(genset *g, int *A, int s)
{
    if (g->n == g->cap) {
        g->cap = g->cap ? 2 * g->cap : 16;
        g->A = realloc(g->A, g->cap * sizeof(int *));
        g->s = realloc(g->s, g->cap * sizeof(int));
    }
    g->A[g->n] = A;
    g->s[g->n] = s;
    g->n++;
}

static void gs_free(genset *g)
{
    for (int i = 0; i < g->n; i++) free(g->A[i]);
    free(g->A);
    free(g->s);
}

/* |ker(A - 1)| on H, as the lattice index [Z^m : (A-1)Z^m + diag(d)Z^m]   */
static u128 fixcount(const hgrp *H, const int *A)
{
    int m = H->m;
    long long L = H->L;
    /* columns: m of (A-1), m of diag(d), plus up to m extra */
    int maxc = 3 * m + 1;
    long long *C = malloc(sizeof(long long) * maxc * m);
    int *active = malloc(sizeof(int) * maxc);
    int nc = 0;
#define COL(c, r) C[(c) * m + (r)]
    for (int c = 0; c < m; c++, nc++) {
        for (int r = 0; r < m; r++) {
            long long v = A[r * m + c] - (r == c ? 1 : 0);
            v %= L; if (v < 0) v += L;
            COL(nc, r) = v;
        }
        active[nc] = 1;
    }
    for (int c = 0; c < m; c++, nc++) {
        for (int r = 0; r < m; r++) COL(nc, r) = (r == c) ? (H->d[c] % L) : 0;
        active[nc] = 1;
    }
    u128 index = 1;
    for (int r = 0; r < m; r++) {
        int piv = -1;
        for (int c = 0; c < nc; c++) {
            if (!active[c] || COL(c, r) == 0) continue;
            if (piv < 0) { piv = c; continue; }
            /* combine columns piv and c so that c gets 0 in row r */
            long long a = COL(piv, r), b = COL(c, r);
            long long g, u, v;
            {   /* extended gcd */
                long long r0 = a, r1 = b, s0 = 1, s1 = 0, t0 = 0, t1 = 1;
                while (r1) {
                    long long q = r0 / r1, tmp;
                    tmp = r0 - q * r1; r0 = r1; r1 = tmp;
                    tmp = s0 - q * s1; s0 = s1; s1 = tmp;
                    tmp = t0 - q * t1; t0 = t1; t1 = tmp;
                }
                g = r0; u = s0; v = t0;
            }
            long long ag = a / g, bg = b / g;
            for (int k = r; k < m; k++) {
                __int128 p = COL(piv, k), q = COL(c, k);
                __int128 np = (u * p + v * q) % L;
                __int128 nq = (-bg * p + ag * q) % L;
                if (np < 0) np += L;
                if (nq < 0) nq += L;
                COL(piv, k) = (long long)np;
                COL(c, k) = (long long)nq;
            }
        }
        if (piv < 0) {
            /* row r is 0 mod L in every active column: L e_r is in the lattice */
            index *= (u128)L;
            continue;
        }
        long long a = COL(piv, r);
        /* include L e_r: (a, L) -> (g, 0) */
        long long g, u;
        {
            long long r0 = a, r1 = L, s0 = 1, s1 = 0;
            while (r1) {
                long long q = r0 / r1, tmp;
                tmp = r0 - q * r1; r0 = r1; r1 = tmp;
                tmp = s0 - q * s1; s0 = s1; s1 = tmp;
            }
            g = r0; u = s0;
        }
        (void)u;
        index *= (u128)g;
        /* Q = (-L/g) * P + (a/g) * L e_r has zero in row r; lower rows
           are (-L/g) * P_lower.  Add it if nonzero. */
        long long Lg = L / g;
        int nz = 0;
        if (nc >= maxc) die("fixcount: column overflow", NULL);
        for (int k = 0; k < m; k++) COL(nc, k) = 0;
        for (int k = r + 1; k < m; k++) {
            __int128 q = ((__int128)(-Lg) * COL(piv, k)) % L;
            if (q < 0) q += L;
            COL(nc, k) = (long long)q;
            if (q) nz = 1;
        }
        active[piv] = 0;
        if (nz) { active[nc] = 1; nc++; }
    }
#undef COL
    free(C);
    free(active);
    return index;
}

/* brute force |ker(A-1)|, for cross-checking fixcount on small H */
static u128 fixcount_brute(const hgrp *H, const int *A)
{
    int m = H->m;
    u128 h = 1;
    for (int i = 0; i < m; i++) h *= H->d[i];
    int x[MAXM];
    u128 cnt = 0;
    for (u128 v = 0; v < h; v++) {
        u128 t = v;
        for (int i = 0; i < m; i++) { x[i] = (int)(t % H->d[i]); t /= H->d[i]; }
        int ok = 1;
        for (int r = 0; r < m && ok; r++) {
            long long s = 0;
            for (int c = 0; c < m; c++) s += (long long)A[r * m + c] * x[c];
            s -= x[r];
            if (s % H->d[r] != 0) ok = 0;
        }
        cnt += ok;
    }
    return cnt;
}

/* Is x -> A x an automorphism of H?  Well-defined: d_r | A[r][c] d_c.
   Bijective: |ker A| = [Z^m : A Z^m + diag(d) Z^m] = 1. */
static int check_aut(const hgrp *H, const int *A, char *why)
{
    int m = H->m;
    for (int r = 0; r < m; r++)
        for (int c = 0; c < m; c++)
            if (((long long)A[r * m + c] * H->d[c]) % H->d[r] != 0) {
                sprintf(why, "not well-defined at (%d,%d)", r, c);
                return 0;
            }
    int *B = malloc(sizeof(int) * m * m);
    for (int r = 0; r < m; r++)
        for (int c = 0; c < m; c++) B[r * m + c] = A[r * m + c] + (r == c);
    u128 k = fixcount(H, B);     /* = |ker((B - 1))| = |ker A| */
    free(B);
    if (k != 1) {
        sprintf(why, "kernel of order %llu", (unsigned long long)k);
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------------ */
/* Burnside over the image I of N in Aut(H) x {+-1}                          */

typedef struct {
    int m, mm;
    long n, cap;
    uint16_t *E;        /* n * mm */
    int8_t *S;
    int32_t *tab;       /* open addressing, -1 empty */
    long tabsize;
} imgset;

static uint64_t hash_elt(const uint16_t *e, int mm, int s)
{
    uint64_t h = 1469598103934665603ULL ^ (uint64_t)(s + 2);
    for (int i = 0; i < mm; i++) { h ^= e[i]; h *= 1099511628211ULL; }
    h ^= h >> 29;
    return h;
}

static long img_find_or_add(imgset *I, const uint16_t *e, int s, int *added)
{
    uint64_t h = hash_elt(e, I->mm, s);
    long pos = (long)(h & (I->tabsize - 1));
    for (;;) {
        int32_t k = I->tab[pos];
        if (k < 0) break;
        if (I->S[k] == s && memcmp(I->E + (size_t)k * I->mm, e, I->mm * sizeof(uint16_t)) == 0) {
            *added = 0;
            return k;
        }
        pos = (pos + 1) & (I->tabsize - 1);
    }
    if (I->n >= I->cap) return -1;
    memcpy(I->E + (size_t)I->n * I->mm, e, I->mm * sizeof(uint16_t));
    I->S[I->n] = (int8_t)s;
    I->tab[pos] = (int32_t)I->n;
    *added = 1;
    return I->n++;
}

/* returns 1 on success, 0 if the image is larger than MAX_IMAGE */
static int burnside(const hgrp *H, const genset *G, u128 *orbH, u128 *orbHS, long *imgsize)
{
    int m = H->m, mm = m * m;
    imgset I;
    I.m = m; I.mm = mm; I.n = 0; I.cap = MAX_IMAGE;
    {   /* keep the element store under ~768 MB per process */
        long memcap = (768L << 20) / (long)(mm * sizeof(uint16_t) + 9);
        if (I.cap > memcap) I.cap = memcap;
    }
    I.tabsize = 1;
    while (I.tabsize < 2 * I.cap) I.tabsize <<= 1;
    I.E = malloc((size_t)I.cap * mm * sizeof(uint16_t));
    I.S = malloc((size_t)I.cap);
    I.tab = malloc((size_t)I.tabsize * sizeof(int32_t));
    if (!I.E || !I.S || !I.tab) die("burnside: out of memory", NULL);
    memset(I.tab, 0xff, (size_t)I.tabsize * sizeof(int32_t));

    uint16_t *e = calloc(mm, sizeof(uint16_t));
    for (int i = 0; i < m; i++) e[i * m + i] = (uint16_t)(1 % H->d[i]);
    int added;
    img_find_or_add(&I, e, 1, &added);

    int ok = 1;
    for (long k = 0; k < I.n && ok; k++) {
        for (int j = 0; j < G->n && ok; j++) {
            const int *A = G->A[j];
            const uint16_t *B = I.E + (size_t)k * mm;
            for (int r = 0; r < m; r++) {
                for (int c = 0; c < m; c++) {
                    long long s = 0;
                    for (int t = 0; t < m; t++) s += (long long)A[r * m + t] * B[t * m + c];
                    e[r * m + c] = (uint16_t)(s % H->d[r]);
                }
            }
            if (img_find_or_add(&I, e, G->s[j] * I.S[k], &added) < 0) ok = 0;
        }
    }
    free(e);
    if (!ok) {
        free(I.E); free(I.S); free(I.tab);
        *imgsize = -1;
        return 0;
    }
    *imgsize = I.n;

    u128 sumAll = 0, sumPlus = 0;
    int *A = malloc(sizeof(int) * mm);
    for (long k = 0; k < I.n; k++) {
        for (int i = 0; i < mm; i++) A[i] = I.E[(size_t)k * mm + i];
        u128 f = fixcount(H, A);
        if (CROSSCHECK && H->log2h <= 14.0) {
            u128 fb = fixcount_brute(H, A);
            if (fb != f) die("fixcount mismatch", NULL);
        }
        sumAll += f;
        if (I.S[k] > 0) sumPlus += 2 * f;
    }
    free(A);
    if (sumAll % (u128)I.n || sumPlus % (u128)I.n) die("Burnside sums not divisible", NULL);
    *orbH = sumAll / (u128)I.n;
    *orbHS = sumPlus / (u128)I.n;
    free(I.E); free(I.S); free(I.tab);
    return 1;
}

/* ------------------------------------------------------------------------ */
/* explicit orbits on H x {+-1}                                              */

static void explicit_orbits(const hgrp *H, const genset *G, u128 *orbH, u128 *orbHS)
{
    int m = H->m;
    uint64_t h = 1;
    for (int i = 0; i < m; i++) h *= (uint64_t)H->d[i];
    uint64_t nbits = 2 * h;
    uint64_t *vis = calloc((nbits + 63) / 64, sizeof(uint64_t));
    if (!vis) die("explicit: out of memory", NULL);
#define VIS(st) ((vis[(st) >> 6] >> ((st) & 63)) & 1)
#define SETVIS(st) (vis[(st) >> 6] |= (1ULL << ((st) & 63)))
    uint64_t radix[MAXM];
    radix[0] = 1;
    for (int i = 1; i < m; i++) radix[i] = radix[i - 1] * H->d[i - 1];

    size_t qcap = 1 << 16;
    uint64_t *Q = malloc(qcap * sizeof(uint64_t));
    int x[MAXM];
    u128 oH = 0, oHS = 0;

    for (uint64_t v0 = 0; v0 < h; v0++) {
        for (int pass = 0; pass < 2; pass++) {
            uint64_t st0 = 2 * v0 + (uint64_t)pass;   /* pass 1: sign -1 */
            if (VIS(st0)) continue;
            if (pass == 0) oH++;
            oHS++;
            size_t qh = 0, qt = 0;
            SETVIS(st0);
            Q[qt++] = st0;
            while (qh < qt) {
                uint64_t st = Q[qh++];
                uint64_t v = st >> 1;
                int neg = (int)(st & 1);
                for (int i = 0; i < m; i++) { x[i] = (int)(v % H->d[i]); v /= H->d[i]; }
                for (int j = 0; j < G->n; j++) {
                    const int *A = G->A[j];
                    uint64_t w = 0;
                    for (int r = 0; r < m; r++) {
                        long long s = 0;
                        const int *row = A + r * m;
                        for (int c = 0; c < m; c++) s += (long long)row[c] * x[c];
                        w += (uint64_t)(s % H->d[r]) * radix[r];
                    }
                    int nneg = neg ^ (G->s[j] < 0);
                    uint64_t nst = 2 * w + (uint64_t)nneg;
                    if (!VIS(nst)) {
                        SETVIS(nst);
                        if (qt == qcap) { qcap *= 2; Q = realloc(Q, qcap * sizeof(uint64_t)); }
                        Q[qt++] = nst;
                    }
                }
            }
        }
    }
#undef VIS
#undef SETVIS
    free(Q);
    free(vis);
    *orbH = oH;
    *orbHS = oHS;
}

/* ------------------------------------------------------------------------ */

typedef struct {
    u128 affine, proper;
    int improper_pg;    /* point group contains det -1 */
    int method;         /* 0 trivial H, 1 explicit, 2 Burnside */
    int badgens;        /* normalizer generators whose computed action on H
                           is not an automorphism (CARAT normalop failure) */
    long imgsize;
    double log2h;
} zresult;

static void analyze_zclass(bravais_TYP *Z, const char *fname, int zno, zresult *res)
{
    int i;
    memset(res, 0, sizeof(*res));

    int improper = 0;
    for (i = 0; i < Z->gen_no; i++) if (det_sign(Z->gen[i]) < 0) improper = 1;
    res->improper_pg = improper;

    /* presentation of the point group */
    matrix_TYP **base = get_base(Z);
    bahn **strong = strong_generators(base, Z, TRUE);
    matrix_TYP *relmat = pres(strong, Z, NULL);
    for (i = 0; i < Z->dim; i++) {
        free_mat(base[i]);
        free_bahn(strong[i]);
        free(strong[i]);
    }
    free(strong);
    free(base);

    word *relator = calloc(relmat->rows, sizeof(word));
    for (i = 0; i < relmat->rows; i++) matrix_2_word(relmat, relator + i, i);
    matrix_TYP **matinv = calloc(Z->gen_no, sizeof(matrix_TYP *));
    long cdim;
    matrix_TYP **X = cohomology(&cdim, Z->gen, matinv, relator, Z->gen_no, relmat->rows);

    /* the generators of N: the point group itself (acts trivially on H,
       but carries determinants), centralizer and normalizer generators */
    int nN = Z->cen_no + Z->normal_no;
    int anyneg = improper;
    int *nsign = malloc(sizeof(int) * (nN + 1));
    for (i = 0; i < nN; i++) {
        matrix_TYP *M = (i < Z->cen_no) ? Z->cen[i] : Z->normal[i - Z->cen_no];
        nsign[i] = det_sign(M);
        if (nsign[i] < 0) anyneg = 1;
    }

    if (X[0]->cols < 1) {
        if (COHOONLY) { printf("H\t%s\t%d\t1\n", fname, zno); fflush(stdout); }
        res->method = 0;
        res->affine = 1;
        res->proper = anyneg ? 1 : 2;
        res->log2h = 0;
    } else {
        hgrp H;
        matrix_TYP *D = X[1];
        int first, last;
        for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);
        for (last = 0; last < D->cols && D->array.SZ[last][last] != 0; last++);
        H.m = last - first;
        if (H.m > MAXM) die("H^1 rank too large", fname);
        H.log2h = 0;
        for (i = 0; i < H.m; i++) {
            H.d[i] = D->array.SZ[first + i][first + i];
            if (H.d[i] >= 65536) die("elementary divisor too large", fname);
            H.log2h += log2((double)H.d[i]);
        }
        for (i = 0; i + 1 < H.m; i++)
            if (H.d[i + 1] % H.d[i]) die("elementary divisors not in SNF order", fname);
        H.L = H.d[H.m - 1];
        res->log2h = H.log2h;
        if (COHOONLY) {
            printf("H\t%s\t%d\t", fname, zno);
            for (i = 0; i < H.m; i++) printf("%s%d", i ? "x" : "", H.d[i]);
            printf("\n");
            fflush(stdout);
            res->affine = res->proper = 0;
            goto cleanup;
        }

        genset G = {0};
        int m = H.m;
        for (i = 0; i < nN; i++) {
            matrix_TYP *M = (i < Z->cen_no) ? Z->cen[i] : Z->normal[i - Z->cen_no];
            matrix_TYP *Aop = normalop(X[0], X[1], X[2], Z, M, i == nN - 1);
            int *A = malloc(sizeof(int) * m * m);
            int isid = 1;
            for (int r = 0; r < m; r++)
                for (int c = 0; c < m; c++) {
                    int v = Aop->array.SZ[r][c] % H.d[r];
                    if (v < 0) v += H.d[r];
                    A[r * m + c] = v;
                    if (v != ((r == c) ? 1 % H.d[r] : 0)) isid = 0;
                }
            free_mat(Aop);
            {
                char why[128];
                if (!check_aut(&H, A, why)) {
                    res->badgens++;
                    fprintf(stderr, "BADGEN\t%s\tZ%d\tgen %d (det %+d)\t%s\n",
                            fname, zno, i, nsign[i], why);
                }
            }
            if ((isid && nsign[i] > 0) || gs_has(&G, A, nsign[i], m * m)) { free(A); continue; }
            gs_add(&G, A, nsign[i]);
        }
        if (improper) {
            /* some point-group element has det -1 and acts trivially */
            int *A = calloc(m * m, sizeof(int));
            for (int r = 0; r < m; r++) A[r * m + r] = 1 % H.d[r];
            gs_add(&G, A, -1);
        }

        if (DUMP) {
            fprintf(stderr, "Z-class %s #%d: H^1 =", fname, zno);
            for (i = 0; i < m; i++) fprintf(stderr, " Z/%d", H.d[i]);
            fprintf(stderr, "   (%d normalizer gens, %d nontrivial on H)\n", nN, G.n);
            for (int j = 0; j < G.n; j++) {
                fprintf(stderr, "  gen %d det %+d:", j, G.s[j]);
                for (int r = 0; r < m; r++) {
                    fprintf(stderr, " [");
                    for (int c = 0; c < m; c++) fprintf(stderr, "%s%d", c ? " " : "", G.A[j][r * m + c]);
                    fprintf(stderr, "]");
                }
                fprintf(stderr, "\n");
            }
        }

        u128 oH = 0, oHS = 0;
        int done = 0;
        if (H.log2h <= 16.0) {
            explicit_orbits(&H, &G, &oH, &oHS);
            res->method = 1;
            done = 1;
            if (CROSSCHECK) {
                u128 bH, bHS;
                long isz;
                if (burnside(&H, &G, &bH, &bHS, &isz)) {
                    if (bH != oH || bHS != oHS) {
                        char b1[64], b2[64], b3[64], b4[64];
                        fprintf(stderr, "MISMATCH %s Z%d explicit %s/%s burnside %s/%s\n",
                                fname, zno, u128str(oH, b1), u128str(oHS, b2),
                                u128str(bH, b3), u128str(bHS, b4));
                        die("method mismatch", fname);
                    }
                    res->imgsize = isz;
                }
            }
        }
        if (!done) {
            long isz;
            if (burnside(&H, &G, &oH, &oHS, &isz)) {
                res->method = 2;
                res->imgsize = isz;
                done = 1;
                if (CROSSCHECK && H.log2h <= EXPLICIT_LOG2) {
                    u128 eH, eHS;
                    explicit_orbits(&H, &G, &eH, &eHS);
                    if (eH != oH || eHS != oHS) die("method mismatch (burnside vs explicit)", fname);
                }
            }
        }
        if (!done) {
            if (H.log2h <= EXPLICIT_LOG2) {
                explicit_orbits(&H, &G, &oH, &oHS);
                res->method = 1;
                done = 1;
            } else {
                fprintf(stderr, "enantio: %s Z%d: |H| = 2^%.1f and image > %ld; giving up\n",
                        fname, zno, H.log2h, MAX_IMAGE);
                exit(3);
            }
        }
        res->affine = oH;
        res->proper = oHS;
        gs_free(&G);
    }

cleanup:   /* (COHOONLY jumps here before any normalop() call) */
    free(nsign);
    for (i = 0; i < 3; i++) free_mat(X[i]);
    free(X);
    for (i = 0; i < Z->gen_no; i++) if (matinv[i]) free_mat(matinv[i]);
    free(matinv);
    for (i = 0; i < relmat->rows; i++) wordfree(relator + i);
    free(relator);
    free_mat(relmat);
}

static void process_qclass(const char *fname)
{
    double t0 = now();
    bravais_TYP *G = get_bravais((char *)fname);
    int i;

    int improper = 0;
    for (i = 0; i < G->gen_no; i++) if (det_sign(G->gen[i]) < 0) improper = 1;
    if (MINIMAL && improper) {
        printf("%s\t-1\t-1\t-1\t0\t1\t0\t0\n", fname);
        fflush(stdout);
        free_bravais(G);
        return;
    }

    if (G->form == NULL || G->form_no == 0)
        G->form = formspace(G->gen, G->gen_no, 1, &G->form_no);
    if (G->order == 0) {
        matrix_TYP **basis = get_base(G);
        bahn **ST = strong_generators(basis, G, FALSE);
        G->order = 1;
        for (i = 0; i < G->dim; i++) {
            G->order *= ST[i]->length;
            free_mat(basis[i]);
            free_bahn(ST[i]);
            free(ST[i]);
        }
        factorize_new(G->order, G->divisors);
        free(ST);
        free(basis);
    }

    int no;
    bravais_TYP **Classes;
    if (ZFILES) {
        no = 1;
        Classes = malloc(2 * sizeof(bravais_TYP *));
        Classes[0] = G;
        G = NULL;
    } else {
        Classes = q2z(G, &no, FALSE, NULL, TRUE, FALSE);
    }

    /* optionally write every Z-class representative (with its normalizer)
       to $ENANTIO_DUMP_DIR/<family>__<name>__Z<k> for independent checks */
    const char *dumpdir = getenv("ENANTIO_DUMP_DIR");
    char qid[512] = "";
    if (dumpdir) {
        /* .../dir.<family>/ordnung.<o>/<type>/<name> */
        const char *name = strrchr(fname, '/');
        name = name ? name + 1 : fname;
        const char *fam = strstr(fname, "/dir.");
        char famb[256] = "x";
        if (fam) {
            fam += 5;
            int j = 0;
            for (; fam[j] && fam[j] != '/' && j < 255; j++) {
                char ch = fam[j];
                famb[j] = (ch == ';') ? 's' : (ch == ',') ? 'c' : (ch == '\'') ? 'p' : ch;
            }
            famb[j] = 0;
        }
        snprintf(qid, sizeof(qid), "%s__%s", famb, name);
    }

    u128 aff = 0, prop = 0;
    int bad = 0;
    for (i = 0; i < no; i++) {
        zresult r;
        if (dumpdir) {
            char out[1024];
            snprintf(out, sizeof(out), "%s/%s__Z%d", dumpdir, qid, i + 1);
            put_bravais(Classes[i], out, "Z-class representative with normalizer");
        }
        if (DUMPONLY) continue;
        analyze_zclass(Classes[i], fname, i + 1, &r);
        if (r.improper_pg != improper) die("Z-class/Q-class determinant mismatch", fname);
        aff += r.affine;
        prop += r.proper;
        if (r.badgens) bad++;
        if (ZLINES) {
            char b1[64], b2[64];
            fprintf(stderr, "Z\t%s\t%d\t%s\t%s\tmethod=%d\tlog2H=%.2f\timage=%ld\tbad=%d\n",
                    fname, i + 1, u128str(r.affine, b1), u128str(r.proper, b2),
                    r.method, r.log2h, r.imgsize, r.badgens);
        }
    }
    char b1[64], b2[64], b3[64];
    printf("%s\t%d\t%s\t%s\t%s\t%d\t%d\t%.2f\n", fname, no, u128str(aff, b1),
           u128str(prop, b2), u128str(prop - aff, b3), improper, bad, now() - t0);
    fflush(stdout);

    for (i = 0; i < no; i++) free_bravais(Classes[i]);
    free(Classes);
    if (G) free_bravais(G);
}

int main(int argc, char **argv)
{
    read_header(argc, argv);
    if (FILEANZ < 1) {
        fprintf(stderr, "usage: %s [-x=k] [-y=n] [-z] [-m] [-w] qclass_file...\n", argv[0]);
        exit(1);
    }
    if (is_option('x')) EXPLICIT_LOG2 = optionnumber('x');
    if (is_option('y')) MAX_IMAGE = optionnumber('y');
    CROSSCHECK = is_option('z');
    MINIMAL = is_option('m');
    ZLINES = is_option('w');
    DUMP = is_option('e');
    ZFILES = is_option('j');
    DUMPONLY = is_option('k');
    COHOONLY = is_option('b');
    for (int f = 0; f < FILEANZ; f++) process_qclass(FILENAMES[f]);
    cleanup_prime();
    return 0;
}
