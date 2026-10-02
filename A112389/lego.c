/*
 * lego.c — count contiguous buildings of n 2x4 LEGO bricks (OEIS A112389)
 *
 * A building is a finite set of 2x4 bricks placed on the stud lattice, each
 * brick lying horizontally (4 along x, 2 along y: "H") or vertically ("V")
 * in some layer z.  Bricks in one layer may not overlap; two bricks in
 * adjacent layers are joined iff their footprints share a stud cell.  The
 * building must be connected.  a(n) counts buildings up to translation and
 * rotation about the vertical axis (no reflections).
 *
 * Method
 * ------
 * Burnside over C4:  a(n) = (F + S180 + 2*S90) / 4, where F is the number of
 * buildings up to translation ("fixed"), S180 / S90 those invariant under a
 * rotation by 180 / 90 degrees.  All three are split by "refinement" (layer
 * profile) <z0 z1 ... z_{h-1}> as in Deleuran's github.com/LasseD/A112389.
 *
 *  - A profile with an interior layer of size 1 (a "bottleneck") factors:
 *      F(<R1 1 R2>) = F(<R1 1>) * F(<1 R2>) / 2,
 *      S180(<R1 1 R2>) = S180(<R1 1>) * S180(<1 R2>) / 2,   S90 = 0.
 *    (Fix the bottleneck brick; the parts below and above are independent.)
 *  - Bottleneck-free ("fat") profiles are enumerated directly:
 *      F: Redelmeier enumeration of connected brick sets on the infinite
 *         adjacency graph, rooted at the lexicographically least brick
 *         (z, y, x, o).  Same-layer overlaps are a hereditary constraint and
 *         are simply skipped.  Subtrees that cannot end fat are pruned
 *         (#interior single layers > bricks left).  The last brick is not
 *         placed but counted per layer ("leaf counting").
 *      S180/S90: Redelmeier over brick *orbits* under the rotation about a
 *         fixed centre class (4 classes for 180, 2 for 90), rooted at the
 *         least orbit, keeping only unions that are really connected.
 *
 *  - For one fat profile, "-r" uses layer aggregation instead: pick a set A
 *    of pairwise non-adjacent layers, enumerate only the other layers R,
 *    and count the A-layers in closed form by Moebius inversion over the
 *    partitions of R's components (see "Aggregated counter" below).  With
 *    A = every other layer this turns, e.g., <32321> into a 4-brick
 *    enumeration.  R-sets are further reduced by the D4 symmetry of the
 *    plane (only the least image of each R-set is evaluated).
 *
 * Output per refinement uses Deleuran's convention: "<R> count (symmetric)"
 * with count = classes up to rotation and symmetric = classes fixed by a
 * 180-degree rotation, so tables can be diffed against sum-for-size.py.
 *
 * Usage
 * -----
 *   ./lego N [threads]         compute a(1..N), a180(1..N), refinement table
 *                              (plain enumeration; practical for N <= 8)
 *   ./lego -s [threads]        self test against OEIS A112389 / A123829, n<=7
 *   ./lego -r PROFILE [opts] [threads]
 *                              count one fat refinement by aggregation, e.g.
 *                              ./lego -r 353 -a 02
 *        -a LAYERS   aggregated layers (digits, pairwise non-adjacent);
 *                    default: most bricks among layers of size <= 3
 *        -p FRAC     process a random fraction of tasks and extrapolate
 *                    (count and running time estimate)
 *        -d DEPTH    task split depth (default 4)
 *        -n          disable the D4 reduction (cross-check)
 *        -q          skip the task-count pre-pass (no % progress)
 *        -x 1        check every I_4 formula value against plain DFS
 *   ./lego -b PROFILE [threads]
 *                              same refinement by plain enumeration
 */

#include <inttypes.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

typedef unsigned __int128 u128;

#define MAXN 12
#define NZ 16            /* layers representable                           */
#define GW 128           /* grid width/height (offset coordinates)         */
#define OFF 64           /* raw (0,0) -> grid (OFF,OFF)                     */
#define NPL (1u << 19)   /* placement ids: z(4) y(7) x(7) o(1)             */
#define MAXPROF (1 << (MAXN - 1))

/* placement id: (z<<15)|(y<<8)|(x<<1)|o, lexicographic in (z,y,x,o) */
#define PZ(p) ((int)((p) >> 15))
#define PY(p) ((int)(((p) >> 8) & 127))
#define PX(p) ((int)(((p) >> 1) & 127))
#define PO(p) ((int)((p) & 1))
#define MKP(z, x, y, o) (((uint32_t)(z) << 15) | ((uint32_t)(y) << 8) | ((uint32_t)(x) << 1) | (uint32_t)(o))

static int N;                     /* target number of bricks              */
static int NTHR = 1;

/* neighbour tables: for orientation o, the 46 (dx,dy,o2) with overlap */
static int32_t OVD[2][46];        /* pid delta within a layer (incl. self) */
static int NOV[2];
static int32_t CFD[2][45];        /* same-layer conflicts, excluding self  */

static double now(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

static void init_tables(void)
{
    for (int o = 0; o < 2; o++) {
        int w1 = o ? 2 : 4, h1 = o ? 4 : 2, k = 0, kc = 0;
        for (int o2 = 0; o2 < 2; o2++) {
            int w2 = o2 ? 2 : 4, h2 = o2 ? 4 : 2;
            for (int dy = -h2 + 1; dy < h1; dy++)
                for (int dx = -w2 + 1; dx < w1; dx++) {
                    int32_t d = dy * 256 + dx * 2 + (o2 - o);
                    OVD[o][k++] = d;
                    if (d != 0) CFD[o][kc++] = d;
                }
        }
        NOV[o] = k;
        if (k != 46 || kc != 45) { fprintf(stderr, "table error\n"); exit(1); }
    }
}

/* ------------------------------------------------------------------ */
/* profiles (compositions of n) <-> index                              */

static int prof_index(const int *cnt, int h)
{
    int idx = 0, s = 0;
    for (int i = 0; i < h - 1; i++) { s += cnt[i]; idx |= 1 << (s - 1); }
    return idx;
}

static int prof_decode(int idx, int n, int *cnt)   /* returns height */
{
    int h = 0, run = 0;
    for (int i = 0; i < n; i++) {
        run++;
        if (i == n - 1 || (idx >> i & 1)) { cnt[h++] = run; run = 0; }
    }
    return h;
}

static int prof_fat(const int *cnt, int h)
{
    for (int i = 1; i < h - 1; i++) if (cnt[i] == 1) return 0;
    return 1;
}

static void prof_str(const int *cnt, int h, char *buf)
{
    for (int i = 0; i < h; i++) buf[i] = (char)('0' + cnt[i]);
    buf[h] = 0;
}

/* ------------------------------------------------------------------ */
/* per-thread geometry state                                            */

typedef struct {
    u128 occ[NZ][GW];
    uint8_t mark[NPL];
    uint32_t G[8192];
    int cnt[NZ];
    int np;
    uint32_t P[MAXN + 4];
} Geo;

static inline u128 bmask(uint32_t p) { return (u128)(PO(p) ? 0x3 : 0xF) << PX(p); }

static inline int conflict(const Geo *g, uint32_t p)
{
    int z = PZ(p), y = PY(p);
    u128 m = bmask(p);
    const u128 *r = g->occ[z];
    if (PO(p)) return ((r[y] | r[y + 1] | r[y + 2] | r[y + 3]) & m) != 0;
    return ((r[y] | r[y + 1]) & m) != 0;
}

static inline void toggle(Geo *g, uint32_t p)
{
    int z = PZ(p), y = PY(p), hgt = PO(p) ? 4 : 2;
    u128 m = bmask(p);
    for (int i = 0; i < hgt; i++) g->occ[z][y + i] ^= m;
}

/* overlap test between two placements in adjacent layers (or same layer) */
static inline int overlap2d(uint32_t a, uint32_t b)
{
    int wa = PO(a) ? 2 : 4, ha = PO(a) ? 4 : 2, wb = PO(b) ? 2 : 4, hb = PO(b) ? 4 : 2;
    int ax = PX(a), ay = PY(a), bx = PX(b), by = PY(b);
    return bx < ax + wa && ax < bx + wb && by < ay + ha && ay < by + hb;
}

static inline int adjacent(uint32_t a, uint32_t b)
{
    int dz = PZ(a) - PZ(b);
    return (dz == 1 || dz == -1) && overlap2d(a, b);
}

/* ------------------------------------------------------------------ */
/* F: Redelmeier enumeration of fat buildings, rooted at least brick    */

static int CAP[NZ];               /* per-layer caps (refinement mode)    */
static int HMAX;                  /* max height                          */
static int REFMODE;               /* 1: exact profile CAP[0..HMAX-1]      */
static int SPLIT;                 /* task split depth                    */
static atomic_int next_task;

typedef struct {
    Geo g;
    uint32_t root;
    int task, claimed;
    uint64_t F[MAXPROF];
    uint64_t nodes;
} FT;

static inline void tally_leaf(FT *t, const int *lc, int top)
{
    int *cnt = t->g.cnt;
    for (int z = 0; z <= top + 1 && z < NZ; z++) {
        if (!lc[z]) continue;
        cnt[z]++;
        int h = (z > top) ? z + 1 : top + 1;
        if (prof_fat(cnt, h)) t->F[prof_index(cnt, h)] += (uint64_t)lc[z];
        cnt[z]--;
    }
}

static void frec(FT *t, int s, int e, int top)
{
    Geo *g = &t->g;
    for (int i = s; i < e; i++) {
        uint32_t c = g->G[i];
        int z = PZ(c);
        if (g->cnt[z] >= CAP[z] || conflict(g, c)) continue;
        int np1 = g->np + 1;
        if (np1 == SPLIT && SPLIT > 0) {
            if (t->task++ != t->claimed) continue;
            t->claimed = atomic_fetch_add(&next_task, 1);
        }
        toggle(g, c);
        g->cnt[z]++;
        g->P[g->np++] = c;
        int ntop = z > top ? z : top;
        /* prune: interior single layers need a brick each */
        int need = 0;
        for (int l = 1; l < ntop; l++) need += (g->cnt[l] == 1);
        if (REFMODE) need = 0;
        if (need <= N - np1) {
            int k = e;
            const int32_t *ov = OVD[PO(c)];
            for (int j = 0; j < 46; j++) {
                if (z + 1 < HMAX) {
                    uint32_t q = (uint32_t)((int32_t)c + (1 << 15) + ov[j]);
                    if (q > t->root && !g->mark[q]) { g->mark[q] = 1; g->G[k++] = q; }
                }
                if (z > 0) {
                    uint32_t q = (uint32_t)((int32_t)c - (1 << 15) + ov[j]);
                    if (q > t->root && !g->mark[q]) { g->mark[q] = 1; g->G[k++] = q; }
                }
            }
            t->nodes++;
            if (np1 == N - 1) {
                int lc[NZ + 1] = {0};
                for (int j = i + 1; j < k; j++) {
                    uint32_t q = g->G[j];
                    int zq = PZ(q);
                    if (g->cnt[zq] < CAP[zq] && !conflict(g, q)) lc[zq]++;
                }
                tally_leaf(t, lc, ntop);
            } else {
                frec(t, i + 1, k, ntop);
            }
            for (int j = e; j < k; j++) g->mark[g->G[j]] = 0;
        }
        g->np--;
        g->cnt[z]--;
        toggle(g, c);
    }
}

static void *fthread(void *arg)
{
    FT *t = arg;
    t->claimed = atomic_fetch_add(&next_task, 1);
    for (int o = 0; o < 2; o++) {
        Geo *g = &t->g;
        memset(g->mark, 0, sizeof g->mark);
        memset(g->occ, 0, sizeof g->occ);
        memset(g->cnt, 0, sizeof g->cnt);
        uint32_t r = MKP(0, OFF, OFF, o);
        t->root = r;
        g->np = 0;
        /* root as a one-element untried list */
        g->G[0] = r;
        g->mark[r] = 1;
        frec(t, 0, 1, 0);
        g->mark[r] = 0;
    }
    return NULL;
}

/* count fat buildings with N bricks; F[prof] summed over threads */
static void count_fixed(uint64_t *Fout)
{
    memset(Fout, 0, sizeof(uint64_t) * MAXPROF);
    if (N == 1) { Fout[0] = 2; return; }
    if (N == 2) { /* <11>: 2 orientations * 46 */
        Fout[1] = 92;
        return;
    }
    SPLIT = N - 1 > 5 ? 5 : 0;
    if (SPLIT && N - 1 <= SPLIT) SPLIT = 0;
    int nt = SPLIT ? NTHR : 1;
    atomic_store(&next_task, 0);
    FT *ts = calloc((size_t)nt, sizeof(FT));
    pthread_t th[64];
    for (int i = 0; i < nt; i++) pthread_create(&th[i], NULL, fthread, &ts[i]);
    for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
    for (int i = 0; i < nt; i++)
        for (int p = 0; p < MAXPROF; p++) Fout[p] += ts[i].F[p];
    free(ts);
}

/* ------------------------------------------------------------------ */
/* S180 / S90: Redelmeier over orbits under rotation about a centre     */

static int GK;          /* group order: 2 (180) or 4 (90)                */
static int CX2, CY2;    /* twice the centre, raw coords, in {0,1}        */

/* rotate a placement (grid coords) by the group generator */
static inline uint32_t rot(uint32_t p)
{
    int z = PZ(p), x = PX(p), y = PY(p), o = PO(p);
    int nx, ny, no;
    if (GK == 2) {
        if (!o) { nx = CX2 + 2 * OFF - x - 4; ny = CY2 + 2 * OFF - y - 2; }
        else    { nx = CX2 + 2 * OFF - x - 2; ny = CY2 + 2 * OFF - y - 4; }
        no = o;
    } else {
        int s = (CX2 + CY2) / 2, d = (CY2 - CX2) / 2;
        if (!o) { nx = s + 2 * OFF - y - 2; ny = d + x; }
        else    { nx = s + 2 * OFF - y - 4; ny = d + x; }
        no = 1 - o;
    }
    return MKP(z, nx, ny, no);
}

static inline int orbit(uint32_t p, uint32_t *el)
{
    int m = 0;
    uint32_t q = p;
    do { el[m++] = q; q = rot(q); } while (q != p && m < 4);
    return m;
}

static inline uint32_t orep(uint32_t p)
{
    uint32_t el[4], r = p;
    int m = orbit(p, el);
    for (int i = 1; i < m; i++) if (el[i] < r) r = el[i];
    return r;
}

static int orbit_ok(const uint32_t *el, int m)
{
    for (int i = 0; i < m; i++)
        for (int j = i + 1; j < m; j++)
            if (overlap2d(el[i], el[j])) return 0;
    return 1;
}

/* lower bound on bricks strictly between a and b on a joining path */
static inline int gapbound(uint32_t a, uint32_t b)
{
    int dz = abs(PZ(a) - PZ(b)), dx = abs(PX(a) - PX(b)), dy = abs(PY(a) - PY(b));
    int L = dz;
    int lx = (dx + 2) / 3, ly = (dy + 2) / 3;
    if (lx > L) L = lx;
    if (ly > L) L = ly;
    if ((L - dz) & 1) L++;
    if (L == 0 && a != b) L = 2;
    return L > 0 ? L - 1 : 0;
}

typedef struct {
    Geo g;
    uint32_t root;
    int nb;
    uint32_t B[MAXN + 4];
    uint64_t S[MAXPROF];
    uint64_t nodes;
} ST;

static uint32_t *roots;
static int nroots;
static atomic_int next_root;

/* components of the placed bricks; returns number, fills comp[] */
static int components(const ST *t, int *comp)
{
    int nc = 0, st[MAXN + 4];
    for (int i = 0; i < t->nb; i++) comp[i] = -1;
    for (int i = 0; i < t->nb; i++) {
        if (comp[i] >= 0) continue;
        int sp = 0;
        st[sp++] = i;
        comp[i] = nc;
        while (sp) {
            int u = st[--sp];
            for (int v = 0; v < t->nb; v++)
                if (comp[v] < 0 && adjacent(t->B[u], t->B[v])) { comp[v] = nc; st[sp++] = v; }
        }
        nc++;
    }
    return nc;
}

static void srec(ST *t, int s, int e)
{
    Geo *g = &t->g;
    for (int i = s; i < e; i++) {
        uint32_t q = g->G[i];
        uint32_t el[4];
        int m = orbit(q, el);
        if (t->nb + m > N) continue;
        int z = PZ(q);
        if (g->cnt[z] + m > CAP[z]) continue;
        if (!orbit_ok(el, m)) continue;
        int bad = 0;
        for (int j = 0; j < m && !bad; j++) bad = conflict(g, el[j]);
        if (bad) continue;
        for (int j = 0; j < m; j++) { toggle(g, el[j]); t->B[t->nb++] = el[j]; }
        g->cnt[z] += m;
        t->nodes++;
        int comp[MAXN + 4];
        int nc = components(t, comp);
        if (t->nb == N) {
            if (nc == 1) {
                int h = 0;
                while (h < NZ && g->cnt[h]) h++;
                t->S[prof_index(g->cnt, h)]++;
            }
        } else {
            int ok = 1;
            if (nc > 1) {
                int lb = 1 << 20;
                for (int a = 0; a < t->nb; a++)
                    for (int b = a + 1; b < t->nb; b++)
                        if (comp[a] != comp[b]) {
                            int gb = gapbound(t->B[a], t->B[b]);
                            if (gb < lb) lb = gb;
                        }
                if (t->nb + lb > N) ok = 0;
            }
            if (ok) {
                int k = e;
                for (int j = 0; j < m; j++) {
                    uint32_t c = el[j];
                    const int32_t *ov = OVD[PO(c)];
                    int zc = PZ(c);
                    for (int d = 0; d < 46; d++) {
                        if (zc + 1 < HMAX) {
                            uint32_t f = (uint32_t)((int32_t)c + (1 << 15) + ov[d]);
                            uint32_t r = orep(f);
                            if (r > t->root && !g->mark[r]) { g->mark[r] = 1; g->G[k++] = r; }
                        }
                        if (zc > 0) {
                            uint32_t f = (uint32_t)((int32_t)c - (1 << 15) + ov[d]);
                            uint32_t r = orep(f);
                            if (r > t->root && !g->mark[r]) { g->mark[r] = 1; g->G[k++] = r; }
                        }
                    }
                }
                srec(t, i + 1, k);
                for (int j = e; j < k; j++) g->mark[g->G[j]] = 0;
            }
        }
        g->cnt[z] -= m;
        for (int j = 0; j < m; j++) { toggle(g, el[j]); t->nb--; }
    }
}

static void *sthread(void *arg)
{
    ST *t = arg;
    Geo *g = &t->g;
    memset(g->mark, 0, sizeof g->mark);
    memset(g->occ, 0, sizeof g->occ);
    memset(g->cnt, 0, sizeof g->cnt);
    for (;;) {
        int ri = atomic_fetch_add(&next_root, 1);
        if (ri >= nroots) break;
        uint32_t r = roots[ri];
        t->root = r;
        t->nb = 0;
        g->G[0] = r;
        g->mark[r] = 1;
        srec(t, 0, 1);
        g->mark[r] = 0;
    }
    return NULL;
}

/* symmetric counts under group of order gk, summed over centre classes */
static void count_sym(int gk, uint64_t *Sout)
{
    memset(Sout, 0, sizeof(uint64_t) * MAXPROF);
    GK = gk;
    int R = 3 * (N - 1) / 2 + 6;
    ST *ts = calloc((size_t)NTHR, sizeof(ST));
    roots = malloc(sizeof(uint32_t) * (size_t)(2 * (2 * R + 1) * (2 * R + 1)));
    for (int cls = 0; cls < 4; cls++) {
        CX2 = cls & 1;
        CY2 = cls >> 1;
        if (gk == 4 && CX2 != CY2) continue;
        nroots = 0;
        for (int y = -R; y <= R; y++)
            for (int x = -R; x <= R; x++)
                for (int o = 0; o < 2; o++) {
                    uint32_t p = MKP(0, OFF + x, OFF + y, o), el[4];
                    if (orep(p) != p) continue;
                    int m = orbit(p, el);
                    if (m != gk && !(gk == 2 && m == 1)) continue;
                    if (m > N || !orbit_ok(el, m)) continue;
                    if (m > 1) {  /* all orbit elements must be joinable */
                        int lb = gapbound(el[0], el[1]);
                        if (m + lb > N) continue;
                    }
                    roots[nroots++] = p;
                }
        atomic_store(&next_root, 0);
        pthread_t th[64];
        for (int i = 0; i < NTHR; i++) pthread_create(&th[i], NULL, sthread, &ts[i]);
        for (int i = 0; i < NTHR; i++) pthread_join(th[i], NULL);
    }
    for (int i = 0; i < NTHR; i++)
        for (int p = 0; p < MAXPROF; p++) Sout[p] += ts[i].S[p];
    free(roots);
    free(ts);
}

/* ------------------------------------------------------------------ */
/* Aggregated counter for one fat refinement.                           */
/*                                                                      */
/* Choose a set A of pairwise non-adjacent layers ("aggregated"); the   */
/* other layers R are enumerated explicitly.  A brick of an aggregated  */
/* layer touches only R-bricks, so given the R-bricks (possibly         */
/* disconnected, real components K_1..K_r) the aggregated layers are    */
/* independent of each other except through connectivity, and Moebius   */
/* inversion on the partition lattice of {K_i} gives                    */
/*                                                                      */
/*   #buildings with these R-bricks                                     */
/*      = sum_sigma mu(sigma, 1) * prod_{l in A} I_{z_l}(X^l_sigma)     */
/*                                                                      */
/* where X^l = placements in layer l touching an R-brick, X^l_sigma the */
/* ones whose touched components lie inside one block of sigma, and     */
/* I_k(X) = number of k-sets of pairwise non-overlapping bricks of X.   */
/* R is enumerated by Redelmeier on an augmented graph: real overlaps   */
/* between R-layers, plus "virtual" edges between two R-bricks that a   */
/* single brick of an aggregated layer could overlap both (same layer   */
/* or two layers apart).  A building's R-part is connected in this      */
/* graph, so each is generated exactly once, rooted at its least brick. */

static void u128str(u128 v, char *buf);

typedef struct { int8_t dx, dy, o2; } Off;
static Off OVT[2][46];
static int32_t VJD[2][512];      /* same layer, via one brick above/below */
static int NVJ[2];
static int32_t VJX[2][512];      /* two layers apart, via one brick       */
static int NVX[2];

static void init_vj(void)
{
    for (int o = 0; o < 2; o++) {
        int w1 = o ? 2 : 4, h1 = o ? 4 : 2, k = 0;
        for (int o2 = 0; o2 < 2; o2++) {
            int w2 = o2 ? 2 : 4, h2 = o2 ? 4 : 2;
            for (int dy = -h2 + 1; dy < h1; dy++)
                for (int dx = -w2 + 1; dx < w1; dx++) OVT[o][k++] = (Off){(int8_t)dx, (int8_t)dy, (int8_t)o2};
        }
    }
    for (int o = 0; o < 2; o++) {
        static uint8_t seen[64][64][2];
        memset(seen, 0, sizeof seen);
        int n = 0, nx = 0;
        for (int i = 0; i < 46; i++) {
            Off a = OVT[o][i];
            for (int j = 0; j < 46; j++) {
                Off b = OVT[a.o2][j];
                int dx = a.dx + b.dx, dy = a.dy + b.dy, o2 = b.o2;
                if (seen[dx + 32][dy + 32][o2]) continue;
                seen[dx + 32][dy + 32][o2] = 1;
                int32_t d = dy * 256 + dx * 2 + (o2 - o);
                VJX[o][nx++] = d;
                /* same layer: exclude bricks overlapping the origin brick */
                int ov = 0;
                for (int q = 0; q < 46; q++)
                    if (OVT[o][q].dx == dx && OVT[o][q].dy == dy && OVT[o][q].o2 == o2) ov = 1;
                if (!ov) VJD[o][n++] = d;
            }
        }
        NVJ[o] = n;
        NVX[o] = nx;
    }
}

#define MAXC 512
#define CW (MAXC / 64)
#define MAXAGG 6
#define CTQ 20                   /* max mask classes for class tables     */

typedef struct {
    int n, W;
    uint32_t p[MAXC];
    uint16_t mask[MAXC];
    uint64_t adj[MAXC][CW];
    int nm;                  /* distinct masks present */
    uint16_t mv[256];
    uint64_t mbits[256][CW];
    uint8_t cls[MAXC];       /* mask class of each candidate */
    int ct_ok;               /* class tables below are valid */
    int32_t nci[CTQ];
    int32_t eij[CTQ][CTQ];
    int64_t Q[CTQ][CTQ][CTQ];
    int32_t tri[CTQ][CTQ][CTQ];
} Cand;

#define MEMO 1024
typedef struct { uint64_t key; int64_t val; uint32_t gen; } MemoE;

static int AMASK, RMASK, ZR, MR, NAGG;
static int AGL[MAXAGG];          /* aggregated layers                    */

typedef struct {
    Geo g;
    uint32_t root;
    int task, claimed;
    u128 F;
    uint64_t cores;
    uint16_t cmark[NPL];
    Cand cd[MAXAGG];
    uint32_t memogen;
    MemoE memo[MEMO];
    uint8_t fg[4][GW][GW];
} AT;

static inline int inR(int z) { return z >= 0 && z < NZ && (RMASK >> z & 1); }
static inline int inA(int z) { return z >= 0 && z < NZ && (AMASK >> z & 1); }

static void collect(AT *t, Cand *c, int layer, int k, const int *comp)
{
    Geo *g = &t->g;
    c->n = 0;
    for (int i = 0; i < g->np; i++) {
        uint32_t b = g->P[i];
        int dz = layer - PZ(b);
        if (dz != 1 && dz != -1) continue;
        const int32_t *ov = OVD[PO(b)];
        int32_t base = (int32_t)b + dz * (1 << 15);
        for (int j = 0; j < 46; j++) {
            uint32_t q = (uint32_t)(base + ov[j]);
            int kk = t->cmark[q];
            if (!kk) {
                kk = ++c->n;
                t->cmark[q] = (uint16_t)kk;
                c->p[kk - 1] = q;
                c->mask[kk - 1] = 0;
            }
            c->mask[kk - 1] |= (uint16_t)(1u << comp[i]);
        }
    }
    c->W = (c->n + 63) / 64;
    if (k >= 2) {
        for (int i = 0; i < c->n; i++) {
            memset(c->adj[i], 0, sizeof(uint64_t) * (size_t)c->W);
            uint32_t q = c->p[i];
            const int32_t *cf = CFD[PO(q)];
            for (int j = 0; j < 45; j++) {
                int kk = t->cmark[(uint32_t)((int32_t)q + cf[j])];
                if (kk) c->adj[i][(kk - 1) >> 6] |= 1ull << ((kk - 1) & 63);
            }
        }
    }
    for (int i = 0; i < c->n; i++) t->cmark[c->p[i]] = 0;
    c->nm = 0;
    for (int i = 0; i < c->n; i++) {
        int m;
        for (m = 0; m < c->nm; m++) if (c->mv[m] == c->mask[i]) break;
        if (m == c->nm) {
            if (m == 256) { fprintf(stderr, "too many masks\n"); exit(1); }
            c->mv[m] = c->mask[i];
            memset(c->mbits[m], 0, sizeof(uint64_t) * CW);
            c->nm++;
        }
        c->mbits[m][i >> 6] |= 1ull << (i & 63);
        c->cls[i] = (uint8_t)m;
    }
    c->ct_ok = 0;
}

/* Per-class tallies so that I_1..I_3 of any union of classes costs O(q^3): */
/* n_i, edges e_ij, Q_ijl = sum_{v in class i} deg_j(v) deg_l(v), and        */
/* triangles by sorted class triple.  Then for S = union of classes          */
/*   E_S = sum e_ij, Ch_S = (sum_{i,j,l in S} Q_ijl - 2E_S)/2, Tr_S = sum tri */
static void class_tables(Cand *c, int k)
{
    int q = c->nm, n = c->n, W = c->W;
    c->ct_ok = 0;
    if (q > CTQ || k > 3) return;
    memset(c->nci, 0, sizeof c->nci);
    for (int i = 0; i < n; i++) c->nci[c->cls[i]]++;
    if (k >= 2) {
        for (int a = 0; a < q; a++) memset(c->eij[a], 0, sizeof(int32_t) * (size_t)q);
        static __thread int16_t degc[MAXC][CTQ];
        if (k == 3) {
            for (int a = 0; a < q; a++)
                for (int b = 0; b < q; b++) {
                    memset(c->Q[a][b], 0, sizeof(int64_t) * (size_t)q);
                    memset(c->tri[a][b], 0, sizeof(int32_t) * (size_t)q);
                }
            for (int i = 0; i < n; i++) memset(degc[i], 0, sizeof(int16_t) * (size_t)q);
        }
        for (int u = 0; u < n; u++) {
            int cu = c->cls[u];
            for (int w = u >> 6; w < W; w++) {
                uint64_t y = c->adj[u][w];
                if (w == (u >> 6)) y &= ~((2ull << (u & 63)) - 1);
                while (y) {
                    int v = w * 64 + __builtin_ctzll(y);
                    y &= y - 1;
                    int cv = c->cls[v];
                    if (cu <= cv) c->eij[cu][cv]++; else c->eij[cv][cu]++;
                    if (k < 3) continue;
                    degc[u][cv]++;
                    degc[v][cu]++;
                    for (int x = v >> 6; x < W; x++) {
                        uint64_t z = c->adj[u][x] & c->adj[v][x];
                        if (x == (v >> 6)) z &= ~((2ull << (v & 63)) - 1);
                        while (z) {
                            int tw = x * 64 + __builtin_ctzll(z);
                            z &= z - 1;
                            int a = cu, b = cv, d = c->cls[tw], tmp;
                            if (a > b) { tmp = a; a = b; b = tmp; }
                            if (b > d) { tmp = b; b = d; d = tmp; }
                            if (a > b) { tmp = a; a = b; b = tmp; }
                            c->tri[a][b][d]++;
                        }
                    }
                }
            }
        }
        if (k == 3) {
            for (int v = 0; v < n; v++) {
                int nz[CTQ], nn = 0, cv = c->cls[v];
                for (int j = 0; j < q; j++) if (degc[v][j]) nz[nn++] = j;
                for (int a = 0; a < nn; a++)
                    for (int b = 0; b < nn; b++)
                        c->Q[cv][nz[a]][nz[b]] += (int64_t)degc[v][nz[a]] * degc[v][nz[b]];
            }
        }
    }
    c->ct_ok = k;
}

static int64_t ind_class(const Cand *c, uint64_t key, int k)
{
    int s[CTQ], ns = 0;
    for (int i = 0; i < c->nm; i++) if (key >> i & 1) s[ns++] = i;
    int64_t N0 = 0;
    for (int a = 0; a < ns; a++) N0 += c->nci[s[a]];
    if (k == 1) return N0;
    int64_t E = 0;
    for (int a = 0; a < ns; a++)
        for (int b = a; b < ns; b++) E += c->eij[s[a]][s[b]];
    if (k == 2) return N0 * (N0 - 1) / 2 - E;
    int64_t SQ = 0, Tr = 0;
    for (int a = 0; a < ns; a++)
        for (int b = 0; b < ns; b++)
            for (int d = 0; d < ns; d++) SQ += c->Q[s[a]][s[b]][s[d]];
    for (int a = 0; a < ns; a++)
        for (int b = a; b < ns; b++)
            for (int d = b; d < ns; d++) Tr += c->tri[s[a]][s[b]][s[d]];
    int64_t Ch = (SQ - 2 * E) / 2;
    return N0 * (N0 - 1) * (N0 - 2) / 6 - E * (N0 - 2) + Ch - Tr;
}

/* I_4 by inclusion-exclusion over edge sets spanning <= 4 vertices:      */
/*   I4 = C(n,4) - E C(n-2,2) + P3 (n-3) + 2K2 - K3 (n-3) - P4 - K13       */
/*        + C4 + Paw - Diamond + K4                                        */
/* (subgraph, not induced, counts).  Bricks are axis-parallel boxes, so    */
/* pairwise-overlapping bricks share a cell (Helly); K4 is counted over    */
/* grid features (cells +, adjacent cell pairs -, 2x2 blocks +).           */
static int DEBUG4;

static inline int64_t C2(int64_t x) { return x * (x - 1) / 2; }
static inline int64_t C3(int64_t x) { return x * (x - 1) * (x - 2) / 6; }
static inline int64_t C4n(int64_t x) { return x * (x - 1) * (x - 2) * (x - 3) / 24; }

static int64_t k4_features(const Cand *c, const int *vs, int n, uint8_t (*fg)[GW][GW])
{
    int64_t K4 = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (int i = 0; i < n; i++) {
            uint32_t p = c->p[vs[i]];
            int x = PX(p), y = PY(p), w = PO(p) ? 2 : 4, h = PO(p) ? 4 : 2;
            for (int yy = y; yy < y + h; yy++)
                for (int xx = x; xx < x + w; xx++) {
                    if (pass) { fg[0][yy][xx] = 0; if (xx + 1 < x + w) fg[1][yy][xx] = 0;
                                if (yy + 1 < y + h) fg[2][yy][xx] = 0;
                                if (xx + 1 < x + w && yy + 1 < y + h) fg[3][yy][xx] = 0; continue; }
                    K4 += C3(fg[0][yy][xx]++);
                    if (xx + 1 < x + w) K4 -= C3(fg[1][yy][xx]++);
                    if (yy + 1 < y + h) K4 -= C3(fg[2][yy][xx]++);
                    if (xx + 1 < x + w && yy + 1 < y + h) K4 += C3(fg[3][yy][xx]++);
                }
        }
    }
    return K4;
}

static int64_t ind4_formula(const Cand *c, const uint64_t *X, uint8_t (*fg)[GW][GW])
{
    int W = c->W, vs[MAXC], n = 0;
    static __thread int d[MAXC];
    static __thread int64_t tv[MAXC];
    for (int w = 0; w < W; w++) {
        uint64_t x = X[w];
        while (x) { vs[n++] = w * 64 + __builtin_ctzll(x); x &= x - 1; }
    }
    if (n < 4) return 0;
    for (int i = 0; i < n; i++) {
        int u = vs[i], dd = 0;
        for (int w = 0; w < W; w++) dd += __builtin_popcountll(c->adj[u][w] & X[w]);
        d[u] = dd;
        tv[u] = 0;
    }
    int64_t E = 0, S3 = 0, P4a = 0, Dia = 0, C4p = 0;
    uint64_t H[CW];
    for (int i = 0; i < n; i++) {
        int u = vs[i];
        for (int w = 0; w < W; w++) H[w] = 0;
        for (int w = 0; w < W; w++) {
            uint64_t y = c->adj[u][w] & X[w];
            while (y) {
                int v = w * 64 + __builtin_ctzll(y);
                y &= y - 1;
                for (int q = 0; q < W; q++) H[q] |= c->adj[v][q];
                if (v < u) continue;
                int cc = 0;
                for (int q = 0; q < W; q++) cc += __builtin_popcountll(c->adj[u][q] & c->adj[v][q] & X[q]);
                E++;
                S3 += cc;
                tv[u] += cc;
                tv[v] += cc;
                Dia += C2(cc);
                C4p += C2(cc);
                P4a += (int64_t)(d[u] - 1) * (d[v] - 1);
            }
        }
        /* non-adjacent pairs u < w at distance 2 */
        for (int w = 0; w < W; w++) {
            uint64_t above = (w < u / 64) ? 0 : (w > u / 64) ? ~0ull : ~((2ull << (u & 63)) - 1);
            H[w] &= X[w] & ~c->adj[u][w] & above;
        }
        for (int w = 0; w < W; w++) {
            uint64_t y = H[w];
            while (y) {
                int v = w * 64 + __builtin_ctzll(y);
                y &= y - 1;
                int cc = 0;
                for (int q = 0; q < W; q++) cc += __builtin_popcountll(c->adj[u][q] & c->adj[v][q] & X[q]);
                C4p += C2(cc);
            }
        }
    }
    int64_t K3 = S3 / 3, P3 = 0, K13 = 0, Paw = 0;
    for (int i = 0; i < n; i++) {
        int u = vs[i];
        P3 += C2(d[u]);
        K13 += C3(d[u]);
        Paw += (tv[u] / 2) * (d[u] - 2);
    }
    int64_t D2 = C2(E) - P3, P4 = P4a - 3 * K3, C4 = C4p / 2;
    int64_t K4 = k4_features(c, vs, n, fg);
    return C4n(n) - E * C2(n - 2) + P3 * (n - 3) + D2 - K3 * (n - 3) - P4 - K13 + C4 + Paw - Dia + K4;
}

static int64_t ind_dfs(const Cand *c, const uint64_t *X, int k, int W);

/* I_k by one step of the recursion I_k = (1/k) sum_u I_{k-1}(X - N[u])   */
/* is not used; k >= 5 falls back to DFS over independent (k-4)-sets with */
/* the I_4 formula at the bottom.                                          */
static int64_t ind_big(const Cand *c, const uint64_t *X, int k, uint8_t (*fg)[GW][GW])
{
    if (k == 4) return ind4_formula(c, X, fg);
    int W = c->W;
    int64_t s = 0;
    uint64_t Y[CW];
    for (int w = 0; w < W; w++) {
        uint64_t x = X[w];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int u = w * 64 + b;
            for (int v = 0; v < W; v++) Y[v] = (v < w ? 0 : v == w ? x : X[v]) & ~c->adj[u][v];
            s += ind_big(c, Y, k - 1, fg);
        }
    }
    return s;
}

static int64_t ind_dfs(const Cand *c, const uint64_t *X, int k, int W)
{
    if (k == 1) {
        int64_t s = 0;
        for (int w = 0; w < W; w++) s += __builtin_popcountll(X[w]);
        return s;
    }
    int64_t s = 0;
    uint64_t Y[CW];
    for (int w = 0; w < W; w++) {
        uint64_t x = X[w];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int u = w * 64 + b;
            for (int v = 0; v < W; v++) Y[v] = (v < w ? 0 : v == w ? x : X[v]) & ~c->adj[u][v];
            s += ind_dfs(c, Y, k - 1, W);
        }
    }
    return s;
}

/* number of independent k-subsets of X */
static int64_t ind_count(const Cand *c, const uint64_t *X, int k, uint8_t (*fg)[GW][GW])
{
    int W = c->W;
    int64_t N0 = 0;
    for (int w = 0; w < W; w++) N0 += __builtin_popcountll(X[w]);
    if (k == 0) return 1;
    if (k == 1 || N0 < 2) return k == 1 ? N0 : 0;
    if (k > 3) {
        if (!fg || DEBUG4 == 2) return ind_dfs(c, X, k, W);
        int64_t v = ind_big(c, X, k, fg);
        if (DEBUG4) {
            int64_t v2 = ind_dfs(c, X, k, W);
            if (v != v2) { fprintf(stderr, "I_%d mismatch: formula %" PRId64 " dfs %" PRId64 "\n", k, v, v2); exit(2); }
        }
        return v;
    }
    int64_t E2 = 0, Ch = 0, Tr = 0;
    for (int w = 0; w < W; w++) {
        uint64_t x = X[w];
        while (x) {
            int b = __builtin_ctzll(x);
            x &= x - 1;
            int u = w * 64 + b;
            int64_t d = 0;
            for (int v = 0; v < W; v++) d += __builtin_popcountll(c->adj[u][v] & X[v]);
            E2 += d;
            Ch += d * (d - 1) / 2;
            if (k == 3) {
                /* triangles u < v < t */
                for (int v = w; v < W; v++) {
                    uint64_t y = c->adj[u][v] & X[v];
                    if (v == w) y &= ~((2ull << b) - 1);
                    while (y) {
                        int bb = __builtin_ctzll(y);
                        y &= y - 1;
                        int vv = v * 64 + bb;
                        for (int q = v; q < W; q++) {
                            uint64_t zz = c->adj[u][q] & c->adj[vv][q] & X[q];
                            if (q == v) zz &= ~((2ull << bb) - 1);
                            Tr += __builtin_popcountll(zz);
                        }
                    }
                }
            }
        }
    }
    int64_t E = E2 / 2;
    if (k == 2) return N0 * (N0 - 1) / 2 - E;
    return N0 * (N0 - 1) * (N0 - 2) / 6 - E * (N0 - 2) + Ch - Tr;
}

/* restricted-growth enumeration of partitions of r elements */
#define MAXR 10
static int NPART[MAXR + 1];
static uint8_t (*PARTS[MAXR + 1])[MAXR];
static uint8_t *PNB[MAXR + 1];   /* number of blocks                    */
static uint16_t (*PCM[MAXR + 1])[MAXR];  /* mask of the block holding i */

static void init_parts(void)
{
    static const int bell[] = {1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975};
    for (int r = 1; r <= MAXR; r++) {
        int cap = bell[r];
        PARTS[r] = malloc(sizeof(uint8_t[MAXR]) * (size_t)cap);
        PNB[r] = malloc((size_t)cap);
        PCM[r] = malloc(sizeof(uint16_t[MAXR]) * (size_t)cap);
        NPART[r] = 0;
        uint8_t a[MAXR] = {0};
        for (;;) {
            int k = 0;
            for (int i = 0; i < r; i++) if (a[i] + 1 > k) k = a[i] + 1;
            PNB[r][NPART[r]] = (uint8_t)k;
            for (int i = 0; i < r; i++) {
                uint16_t bm = 0;
                for (int j = 0; j < r; j++) if (a[j] == a[i]) bm |= (uint16_t)(1u << j);
                PCM[r][NPART[r]][i] = bm;
            }
            memcpy(PARTS[r][NPART[r]++], a, MAXR);
            int i = r - 1;
            for (; i > 0; i--) {
                int mx = 0;
                for (int j = 0; j < i; j++) if (a[j] > mx) mx = a[j];
                if (a[i] <= mx) { a[i]++; for (int j = i + 1; j < r; j++) a[j] = 0; break; }
            }
            if (i == 0) break;
        }
        if (NPART[r] != cap) { fprintf(stderr, "partition error\n"); exit(1); }
    }
}

static void xsigma(const Cand *c, const uint8_t *blk, uint64_t *X)
{
    memset(X, 0, sizeof(uint64_t) * (size_t)c->W);
    for (int m = 0; m < c->nm; m++) {
        unsigned mk = c->mv[m];
        int b0 = -1, ok = 1;
        for (int i = 0; mk; i++, mk >>= 1)
            if (mk & 1) { if (b0 < 0) b0 = blk[i]; else if (blk[i] != b0) { ok = 0; break; } }
        if (ok) for (int w = 0; w < c->W; w++) X[w] |= c->mbits[m][w];
    }
}

/* bitmask over c's mask classes that are allowed by partition s of r */
static inline uint64_t allowed_classes(const Cand *c, int r, int s)
{
    const uint16_t *cm = PCM[r][s];
    uint64_t key = 0;
    for (int m = 0; m < c->nm; m++) {
        unsigned mk = c->mv[m];
        if ((mk & ~(unsigned)cm[__builtin_ctz(mk)]) == 0) key |= 1ull << m;
    }
    return key;
}

static int SAME01;               /* h=3, A={0,2}: both layers see layer 1 */
static int USECT = 5;            /* class tables when >= this many sigmas  */


static int64_t memo_ind(AT *t, MemoE *tab, uint32_t gen, const Cand *c, uint64_t key, int k, int tag)
{
    uint64_t hk = key * 0x9E3779B97F4A7C15ull ^ (uint64_t)(tag * 131 + k);
    unsigned h = (unsigned)(hk >> 54) & (MEMO - 1);
    uint64_t full = key ^ ((uint64_t)(tag * 131 + k) << 58);
    for (;;) {
        MemoE *e = &tab[h];
        if (e->gen != gen) break;
        if (e->key == full) return e->val;
        h = (h + 1) & (MEMO - 1);
    }
    int64_t v;
    if (c->ct_ok >= k && !DEBUG4) {
        v = ind_class(c, key, k);
    } else {
        uint64_t X[CW] = {0};
        for (int m = 0; m < c->nm; m++)
            if (key >> m & 1) for (int w = 0; w < c->W; w++) X[w] |= c->mbits[m][w];
        v = ind_count(c, X, k, t->fg);
        if (DEBUG4 && c->ct_ok >= k && ind_class(c, key, k) != v) {
            fprintf(stderr, "class-table mismatch k=%d\n", k);
            exit(2);
        }
    }
    tab[h] = (MemoE){full, v, gen};
    (void)t;
    return v;
}

static __int128 core_count(AT *t)
{
    Geo *g = &t->g;
    int m = g->np, comp[MAXN + 4], r = 0, st[MAXN + 4];
    for (int i = 0; i < m; i++) comp[i] = -1;
    for (int i = 0; i < m; i++) {
        if (comp[i] >= 0) continue;
        int sp = 0;
        st[sp++] = i;
        comp[i] = r;
        while (sp) {
            int u = st[--sp];
            for (int v = 0; v < m; v++)
                if (comp[v] < 0 && adjacent(g->P[u], g->P[v])) { comp[v] = r; st[sp++] = v; }
        }
        r++;
    }
    if (r > MAXR) { fprintf(stderr, "too many components\n"); exit(1); }
    for (int a = 0; a < NAGG; a++) {
        if (a == 1 && SAME01) continue;
        collect(t, &t->cd[a], AGL[a], CAP[AGL[a]] > CAP[AGL[SAME01 ? 1 : a]] ? CAP[AGL[a]] : CAP[AGL[SAME01 ? 1 : a]], comp);
        if (t->cd[a].n < CAP[AGL[a]]) return 0;
    }
    if (SAME01 && t->cd[0].n < CAP[AGL[1]]) return 0;
    if (r == 1) {
        __int128 prod = 1;
        for (int a = 0; a < NAGG && prod; a++) {
            const Cand *c = &t->cd[(a == 1 && SAME01) ? 0 : a];
            uint64_t X[CW];
            for (int w = 0; w < c->W; w++) X[w] = ~0ull;
            if (c->n & 63) X[c->W - 1] = (1ull << (c->n & 63)) - 1;
            prod *= ind_count(c, X, CAP[AGL[a]], t->fg);
        }
        return prod;
    }
    for (int a = 0; a < NAGG; a++) {
        if (a == 1 && SAME01) continue;
        int kk = CAP[AGL[a]];
        if (a == 0 && SAME01 && CAP[AGL[1]] > kk) kk = CAP[AGL[1]];
        if (kk <= 3 && (!USECT || NPART[r] >= USECT)) class_tables(&t->cd[a], kk);
    }
    uint32_t gen = ++t->memogen;
    if (!gen) { memset(t->memo, 0, sizeof t->memo); gen = t->memogen = 1; }
    __int128 tot = 0;
    for (int s = 0; s < NPART[r]; s++) {
        int k = PNB[r][s];
        __int128 prod = 1;
        for (int a = 0; a < NAGG && prod; a++) {
            int ca = (a == 1 && SAME01) ? 0 : a;
            const Cand *c = &t->cd[ca];
            int kk = CAP[AGL[a]];
            if (c->nm <= 64) {
                uint64_t key = allowed_classes(c, r, s);
                prod *= memo_ind(t, t->memo, gen, c, key, kk, ca);
            } else {
                uint64_t X[CW];
                xsigma(c, PARTS[r][s], X);
                prod *= ind_count(c, X, kk, t->fg);
            }
        }
        if (!prod) continue;
        for (int i = 2; i < k; i++) prod *= i;
        tot += (k & 1) ? prod : -prod;
    }
    return tot;
}

/* D4 symmetry of the plane: Count(R) is invariant under rotations and   */
/* reflections of the R-bricks, so only the least image is evaluated and */
/* weighted by the orbit size.  Returns 0 if P is not the least image.   */
static const int D4M[8][4] = {
    {1, 0, 0, 1}, {0, -1, 1, 0}, {-1, 0, 0, -1}, {0, 1, -1, 0},
    {-1, 0, 0, 1}, {1, 0, 0, -1}, {0, 1, 1, 0}, {0, -1, -1, 0}};

static inline void isort(uint32_t *a, int m)
{
    for (int i = 1; i < m; i++) {
        uint32_t v = a[i];
        int j = i - 1;
        while (j >= 0 && a[j] > v) { a[j + 1] = a[j]; j--; }
        a[j + 1] = v;
    }
}

static int canon_weight(const uint32_t *P, int m)
{
    uint32_t base[MAXN + 4], img[MAXN + 4];
    int bx[MAXN + 4], by[MAXN + 4], bo[MAXN + 4], bz[MAXN + 4];
    memcpy(base, P, sizeof(uint32_t) * (size_t)m);
    isort(base, m);
    int stab = 0;
    for (int gi = 0; gi < 8; gi++) {
        const int *M = D4M[gi];
        int zmin = 99, ymin = 1 << 20, xmin = 1 << 20, omin = 9;
        for (int i = 0; i < m; i++) {
            uint32_t p = base[i];
            int x = PX(p) - OFF, y = PY(p) - OFF, w = PO(p) ? 2 : 4, h = PO(p) ? 4 : 2;
            int ax = M[0] * x + M[1] * y, ay = M[2] * x + M[3] * y;
            int cx = M[0] * (x + w) + M[1] * (y + h), cy = M[2] * (x + w) + M[3] * (y + h);
            int nx = ax < cx ? ax : cx, ny = ay < cy ? ay : cy, nw = ax < cx ? cx - ax : ax - cx;
            bx[i] = nx; by[i] = ny; bo[i] = (nw == 2); bz[i] = PZ(p);
            if (bz[i] < zmin || (bz[i] == zmin && (ny < ymin || (ny == ymin && (nx < xmin || (nx == xmin && bo[i] < omin)))))) {
                zmin = bz[i]; ymin = ny; xmin = nx; omin = bo[i];
            }
        }
        for (int i = 0; i < m; i++) img[i] = MKP(bz[i], bx[i] - xmin + OFF, by[i] - ymin + OFF, bo[i]);
        isort(img, m);
        int c = 0;
        for (int i = 0; i < m && !c; i++) c = (img[i] < base[i]) ? -1 : (img[i] > base[i]) ? 1 : 0;
        if (c < 0) return 0;
        if (c == 0) stab++;
    }
    return 8 / stab;
}

static int USE_D4 = 1;

static inline void push(AT *t, uint32_t q, int *k)
{
    Geo *g = &t->g;
    if (q > t->root && !g->mark[q]) { g->mark[q] = 1; g->G[(*k)++] = q; }
}

static int COUNTTASKS;            /* pre-pass: only count split nodes    */
static double SAMPLE = 1.0;       /* fraction of tasks processed          */
static uint64_t SAMPLE_SEED = 12345;

static inline int sampled(uint64_t task)
{
    if (SAMPLE >= 1.0) return 1;
    uint64_t x = (task + SAMPLE_SEED) * 0x9E3779B97F4A7C15ull;
    x ^= x >> 31; x *= 0xBF58476D1CE4E5B9ull; x ^= x >> 29;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0) < SAMPLE;
}

static void arec(AT *t, int s, int e)
{
    Geo *g = &t->g;
    for (int i = s; i < e; i++) {
        uint32_t c = g->G[i];
        int z = PZ(c);
        if (g->cnt[z] >= CAP[z] || conflict(g, c)) continue;
        int np1 = g->np + 1;
        if (np1 == SPLIT && SPLIT > 0) {
            if (COUNTTASKS) { t->task++; continue; }
            int my = t->task++;
            if (my != t->claimed) continue;
            t->claimed = atomic_fetch_add(&next_task, 1);
            if (!sampled((uint64_t)my)) continue;
        }
        toggle(g, c);
        g->cnt[z]++;
        g->P[g->np++] = c;
        if (np1 == MR) {
            int wgt = USE_D4 ? canon_weight(g->P, g->np) : 1;
            if (wgt) {
                t->F += (u128)(core_count(t) * wgt);
                t->cores++;
            }
        } else {
            int k = e, o = PO(c);
            const int32_t *ov = OVD[o];
            if (inR(z + 1)) for (int j = 0; j < 46; j++) push(t, (uint32_t)((int32_t)c + (1 << 15) + ov[j]), &k);
            if (inR(z - 1)) for (int j = 0; j < 46; j++) push(t, (uint32_t)((int32_t)c - (1 << 15) + ov[j]), &k);
            if (inA(z + 1) || inA(z - 1))
                for (int j = 0; j < NVJ[o]; j++) push(t, (uint32_t)((int32_t)c + VJD[o][j]), &k);
            if (inA(z + 1) && inR(z + 2))
                for (int j = 0; j < NVX[o]; j++) push(t, (uint32_t)((int32_t)c + (2 << 15) + VJX[o][j]), &k);
            if (inA(z - 1) && inR(z - 2))
                for (int j = 0; j < NVX[o]; j++) push(t, (uint32_t)((int32_t)c - (2 << 15) + VJX[o][j]), &k);
            arec(t, i + 1, k);
            for (int j = e; j < k; j++) g->mark[g->G[j]] = 0;
        }
        g->np--;
        g->cnt[z]--;
        toggle(g, c);
    }
}

static void *athread(void *arg)
{
    AT *t = arg;
    if (!COUNTTASKS) t->claimed = atomic_fetch_add(&next_task, 1);
    for (int o = 0; o < 2; o++) {
        Geo *g = &t->g;
        memset(g->mark, 0, sizeof g->mark);
        memset(g->occ, 0, sizeof g->occ);
        memset(g->cnt, 0, sizeof g->cnt);
        uint32_t r = MKP(ZR, OFF, OFF, o);
        t->root = r;
        g->np = 0;
        g->G[0] = r;
        g->mark[r] = 1;
        arec(t, 0, 1);
        g->mark[r] = 0;
    }
    return NULL;
}

static int SPLITMAX = 4;
static int NOPREPASS;
static int64_t TOTTASKS;
static volatile int progress_stop;
static double progress_t0;

static void *progress_thread(void *arg)
{
    (void)arg;
    int tty = isatty(2);
    double last = now();
    while (!progress_stop) {
        usleep(200000);
        double tn = now();
        if (tn - last < (tty ? 2.0 : 60.0)) continue;
        last = tn;
        int64_t done = atomic_load(&next_task);
        if (!TOTTASKS) {
            fprintf(stderr, "%s  tasks claimed %" PRId64 "  elapsed %.0fs   %s", tty ? "\r" : "", done,
                    tn - progress_t0, tty ? "" : "\n");
            fflush(stderr);
            continue;
        }
        if (done > TOTTASKS) done = TOTTASKS;
        double el = tn - progress_t0, fr = TOTTASKS ? (double)done / (double)TOTTASKS : 0;
        double eta = fr > 0 ? el / fr - el : 0;
        fprintf(stderr, "%s  tasks %" PRId64 "/%" PRId64 " (%.2f%%)  elapsed %.0fs  eta %.0fs   %s",
                tty ? "\r" : "", done, TOTTASKS, 100 * fr, el, eta, tty ? "" : "\n");
        fflush(stderr);
    }
    if (tty) fprintf(stderr, "\r%80s\r", "");
    return NULL;
}

/* default aggregated set: max bricks over independent layer sets with z<=3 */
static int choose_agg(const int *prof, int h)
{
    int best = -1, bestm = 0;
    for (int m = 1; m < (1 << h); m++) {
        if (m & (m >> 1)) continue;
        int s = 0, ok = 1;
        for (int l = 0; l < h; l++) if (m >> l & 1) { if (prof[l] > 3) ok = 0; s += prof[l]; }
        if (ok && s > best) { best = s; bestm = m; }
    }
    return bestm;
}

/* F for one fat profile prof[0..h-1] by the aggregated method */
static u128 count_ref_agg(const int *prof, int h, int amask, uint64_t *ncores)
{
    for (int z = 0; z < NZ; z++) CAP[z] = z < h ? prof[z] : 0;
    HMAX = h;
    AMASK = amask;
    RMASK = ((1 << h) - 1) & ~amask;
    SAME01 = (h == 3 && amask == 5);
    NAGG = 0;
    MR = 0;
    ZR = -1;
    for (int l = 0; l < h; l++) {
        if (amask >> l & 1) AGL[NAGG++] = l;
        else { MR += prof[l]; if (ZR < 0) ZR = l; }
    }
    SPLIT = MR - 1 > SPLITMAX ? SPLITMAX : MR - 1;
    if (SPLIT < 2) SPLIT = 0;
    int nt = SPLIT ? NTHR : 1;
    AT *ts = calloc((size_t)nt, sizeof(AT));
    TOTTASKS = 0;
    if (SPLIT && SAMPLE >= 1.0 && !NOPREPASS) {  /* pre-pass: number of tasks, for progress */
        COUNTTASKS = 1;
        athread(&ts[0]);
        TOTTASKS = ts[0].task;
        ts[0].task = 0;
        COUNTTASKS = 0;
    }
    atomic_store(&next_task, 0);
    pthread_t th[64], pth;
    progress_stop = 0;
    progress_t0 = now();
    if (SPLIT) pthread_create(&pth, NULL, progress_thread, NULL);
    for (int i = 0; i < nt; i++) pthread_create(&th[i], NULL, athread, &ts[i]);
    for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
    if (SPLIT) { progress_stop = 1; pthread_join(pth, NULL); }
    u128 F = 0;
    uint64_t nc = 0;
    for (int i = 0; i < nt; i++) { F += ts[i].F; nc += ts[i].cores; }
    free(ts);
    if (ncores) *ncores = nc;
    return F;
}

/* F for one profile by the baseline enumerator */
static u128 count_ref_base(const int *prof, int h)
{
    int n = 0;
    for (int z = 0; z < h; z++) n += prof[z];
    N = n;
    for (int z = 0; z < NZ; z++) CAP[z] = z < h ? prof[z] : 0;
    HMAX = h;
    REFMODE = 1;
    static uint64_t F[MAXPROF];
    count_fixed(F);
    REFMODE = 0;
    return F[prof_index(prof, h)];
}

static int parse_prof(const char *s, int *prof)
{
    int h = 0;
    for (; *s; s++) {
        if (*s < '1' || *s > '9' || h >= NZ) return 0;
        prof[h++] = *s - '0';
    }
    return h;
}

/* base: 1 = baseline enumerator; aspec: NULL (auto) or layer digits */
static void ref_mode(const char *ps, int base, const char *aspec)
{
    int prof[NZ], h = parse_prof(ps, prof), n = 0;
    if (h < 2) { fprintf(stderr, "bad profile\n"); exit(1); }
    for (int i = 0; i < h; i++) n += prof[i];
    if (!prof_fat(prof, h)) { fprintf(stderr, "profile has a bottleneck; use composition\n"); exit(1); }
    int amask = choose_agg(prof, h);
    if (aspec) {
        amask = 0;
        for (const char *p = aspec; *p; p++) if (*p >= '0' && *p <= '9') amask |= 1 << (*p - '0');
        if ((amask & (amask >> 1)) || amask >= (1 << h) || amask == (1 << h) - 1) {
            fprintf(stderr, "aggregated layers must be non-adjacent layers of the profile\n");
            exit(1);
        }
    }
    double t0 = now();
    uint64_t ncores = 0;
    u128 F = base ? count_ref_base(prof, h) : count_ref_agg(prof, h, amask, &ncores);
    double t1 = now();
    if (SAMPLE < 1.0) {
        double Fe = (double)F / SAMPLE;
        char fr[48];
        u128str(F, fr);
        printf("<%s> SAMPLE p=%g: F~%.4e  classes~%.4e  Rsets=%" PRIu64 "  time %.1fs -> full ~%.0fs (%.2f h) on %d threads  [raw %s]\n",
               ps, SAMPLE, Fe, Fe / 4, ncores, t1 - t0, (t1 - t0) / SAMPLE, (t1 - t0) / SAMPLE / 3600, NTHR, fr);
        return;
    }
    N = n;
    for (int z = 0; z < NZ; z++) CAP[z] = z < h ? prof[z] : 0;
    HMAX = h;
    static uint64_t S2[MAXPROF], S4[MAXPROF];
    count_sym(2, S2);
    count_sym(4, S4);
    double t2 = now();
    int idx = prof_index(prof, h);
    u128 s180 = S2[idx], s90 = S4[idx];
    char b1[48], b2[48], b3[48], as[NZ + 1];
    int na = 0;
    for (int l = 0; l < h; l++) if (amask >> l & 1) as[na++] = (char)('0' + l);
    as[na] = 0;
    u128str((F + s180 + 2 * s90) / 4, b1);
    u128str((s180 + s90) / 2, b2);
    u128str(F, b3);
    printf("<%s> %s (%s)   F=%s S180=%" PRIu64 " S90=%" PRIu64 "%s  agg={%s} Rsets=%" PRIu64 "  [F %.1fs, S %.1fs]\n",
           ps, b1, b2, b3, (uint64_t)s180, (uint64_t)s90,
           ((F + s180 + 2 * s90) % 4) ? " NON-INTEGRAL" : "", base ? "-" : as, ncores, t1 - t0, t2 - t1);
    fflush(stdout);
}

/* ------------------------------------------------------------------ */
/* composition over profiles                                            */

typedef struct { u128 F, S180, S90; int known; } PV;
static PV TAB[MAXN + 1][MAXPROF];

static void u128str(u128 v, char *buf)
{
    char tmp[48];
    int k = 0;
    if (!v) { strcpy(buf, "0"); return; }
    while (v) { tmp[k++] = (char)('0' + (int)(v % 10)); v /= 10; }
    for (int i = 0; i < k; i++) buf[i] = tmp[k - 1 - i];
    buf[k] = 0;
}

static PV *lookup(const int *cnt, int h)
{
    int n = 0;
    for (int i = 0; i < h; i++) n += cnt[i];
    return &TAB[n][prof_index(cnt, h)];
}

/* fill bottleneck profiles of size n from smaller ones */
static void compose(int n)
{
    int cnt[MAXN + 1];
    for (int idx = 0; idx < (1 << (n - 1)); idx++) {
        int h = prof_decode(idx, n, cnt);
        if (prof_fat(cnt, h)) continue;
        int b = 1;
        while (cnt[b] != 1) b++;
        PV *lo = lookup(cnt, b + 1), *hi = lookup(cnt + b, h - b);
        PV *v = &TAB[n][idx];
        if (!lo->known || !hi->known) { v->known = 0; continue; }
        v->F = lo->F * hi->F / 2;
        v->S180 = lo->S180 * hi->S180 / 2;
        v->S90 = 0;
        v->known = 1;
    }
}

static void run_all(int nmax, int verbose)
{
    for (int n = 1; n <= nmax; n++) {
        N = n;
        for (int z = 0; z < NZ; z++) CAP[z] = n;
        HMAX = n;
        REFMODE = 0;
        static uint64_t F[MAXPROF], S2[MAXPROF], S4[MAXPROF];
        double t0 = now();
        count_fixed(F);
        double t1 = now();
        count_sym(2, S2);
        count_sym(4, S4);
        double t2 = now();
        int cnt[MAXN + 1];
        for (int idx = 0; idx < (1 << (n - 1)); idx++) {
            int h = prof_decode(idx, n, cnt);
            PV *v = &TAB[n][idx];
            if (prof_fat(cnt, h)) {
                v->F = F[idx]; v->S180 = S2[idx]; v->S90 = S4[idx]; v->known = 1;
            }
        }
        compose(n);
        /* cross-check symmetric counts of bottleneck profiles */
        for (int idx = 0; idx < (1 << (n - 1)); idx++) {
            int h = prof_decode(idx, n, cnt);
            if (prof_fat(cnt, h)) continue;
            if (TAB[n][idx].S180 != S2[idx] || S4[idx]) {
                char ps[16];
                prof_str(cnt, h, ps);
                fprintf(stderr, "WARNING: S180 composition mismatch at <%s>\n", ps);
            }
        }
        u128 tot = 0, tots = 0;
        int ok = 1;
        for (int idx = 0; idx < (1 << (n - 1)); idx++) {
            PV *v = &TAB[n][idx];
            u128 num = v->F + v->S180 + 2 * v->S90;
            if (num % 4 || (v->S180 + v->S90) % 2) ok = 0;
            tot += num / 4;
            tots += (v->S180 + v->S90) / 2;
        }
        char b1[48], b2[48];
        u128str(tot, b1);
        u128str(tots, b2);
        printf("a(%d) = %s   a180(%d) = %s%s   [F %.1fs, S %.1fs]\n", n, b1, n, b2,
               ok ? "" : "  (NON-INTEGRAL!)", t1 - t0, t2 - t1);
        if (verbose && n == nmax) {
            for (int idx = 0; idx < (1 << (n - 1)); idx++) {
                int h = prof_decode(idx, n, cnt);
                PV *v = &TAB[n][idx];
                char ps[16];
                prof_str(cnt, h, ps);
                u128str((v->F + v->S180 + 2 * v->S90) / 4, b1);
                u128str((v->S180 + v->S90) / 2, b2);
                printf("  <%s> %s (%s)%s\n", ps, b1, b2, prof_fat(cnt, h) ? "" : " L1");
            }
        }
        fflush(stdout);
    }
}

static int selftest(void)
{
    static const char *A[] = {"", "1", "24", "1560", "119580", "10166403", "915103765", "85747377755"};
    static const char *S[] = {"", "1", "2", "44", "185", "3276", "15682", "282377"};
    int bad = 0, cnt[MAXN + 1];
    run_all(7, 0);
    for (int n = 1; n <= 7; n++) {
        u128 tot = 0, tots = 0;
        for (int idx = 0; idx < (1 << (n - 1)); idx++) {
            PV *v = &TAB[n][idx];
            tot += (v->F + v->S180 + 2 * v->S90) / 4;
            tots += (v->S180 + v->S90) / 2;
        }
        char b1[48], b2[48];
        u128str(tot, b1);
        u128str(tots, b2);
        int ok = !strcmp(b1, A[n]) && !strcmp(b2, S[n]);
        printf("n=%d  a=%s  a180=%s  %s\n", n, b1, b2, ok ? "OK" : "FAIL");
        bad += !ok;
    }
    (void)cnt;
    printf(bad ? "SELFTEST FAILED\n" : "SELFTEST PASSED\n");
    return bad != 0;
}

int main(int argc, char **argv)
{
    init_tables();
    init_vj();
    init_parts();
    long nc = sysconf(_SC_NPROCESSORS_ONLN);
    NTHR = nc > 0 ? (int)nc : 1;
    if (argc >= 2 && !strcmp(argv[1], "-s")) {
        if (argc >= 3) NTHR = atoi(argv[2]);
        return selftest();
    }
    if (argc >= 3 && (!strcmp(argv[1], "-r") || !strcmp(argv[1], "-b"))) {
        const char *aspec = NULL;
        int ai = 3;
        for (;;) {
            if (argc >= ai + 2 && !strcmp(argv[ai], "-a")) { aspec = argv[ai + 1]; ai += 2; continue; }
            if (argc >= ai + 2 && !strcmp(argv[ai], "-p")) { SAMPLE = atof(argv[ai + 1]); ai += 2; continue; }
            if (argc >= ai + 2 && !strcmp(argv[ai], "-d")) { SPLITMAX = atoi(argv[ai + 1]); ai += 2; continue; }
            if (argc >= ai + 1 && !strcmp(argv[ai], "-n")) { USE_D4 = 0; ai += 1; continue; }
            if (argc >= ai + 1 && !strcmp(argv[ai], "-q")) { NOPREPASS = 1; ai += 1; continue; }
            if (argc >= ai + 2 && !strcmp(argv[ai], "-x")) { DEBUG4 = atoi(argv[ai + 1]); ai += 2; continue; }
            if (argc >= ai + 2 && !strcmp(argv[ai], "-c")) { USECT = atoi(argv[ai + 1]); ai += 2; continue; }
            break;
        }
        if (argc > ai) NTHR = atoi(argv[ai]);
        ref_mode(argv[2], argv[1][1] == 'b', aspec);
        return 0;
    }
    if (argc < 2) {
        fprintf(stderr, "usage: %s N [threads] | -s [threads] | -r PROFILE [threads] | -b PROFILE [threads]\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (argc >= 3) NTHR = atoi(argv[2]);
    if (n < 1 || n > MAXN - 1) { fprintf(stderr, "N out of range\n"); return 1; }
    run_all(n, 1);
    return 0;
}
