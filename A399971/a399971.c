/*
 * a399971.c - terms of OEIS A399971
 *
 *   a(n) = number of compact Numberlink solutions on an n X n grid, counted up
 *          to the 8 symmetries of the square (link numbers are interchangeable).
 *
 * Definitions (https://oeis.org/A399971/a399971.txt):
 *   solution : partition of all n^2 cells into vertex-disjoint grid paths, each
 *              with >= 2 cells.  Its path endpoints form a puzzle (starting
 *              position).
 *   linkage  : any set of vertex-disjoint paths joining every endpoint pair of
 *              a puzzle; cells may be left empty.
 *   compact  : no linkage of the solution's puzzle leaves a cell empty.
 *
 * Method
 *   1. DFS over the cells (outer ring first, orbit by orbit under the square's
 *      symmetries, then the inside row by row), each cell choosing its edges to
 *      the cells after it.  A compact solution only has induced paths (a chord
 *      is a shortcut that frees cells), so two paths are merged only when they
 *      touch along exactly one grid edge; this also rules out cycles.
 *   2. Only the cover that is smallest among its 8 images (cell by cell in
 *      processing order) is kept; each symmetric image is compared as soon as
 *      the cells it needs are decided, which cuts most other branches early.
 *   3. Compactness: search for a "witness", a linkage of the endpoint pairs that
 *      leaves some cell empty.  Paths are grown one cell at a time, forced moves
 *      first; every open pair must stay connected through free cells; only
 *      induced paths are grown (shortcutting a witness gives an induced
 *      witness); the last open pair is finished by a BFS shortest path.  With
 *      no witness the search has visited every linkage of the puzzle (all of
 *      them full and induced), which also gives the starting positions count.
 *   The DFS is split into prefixes that a pool of threads works through.
 *
 * Build: cc -O3 -march=native -pthread -o a399971 a399971.c
 * Usage: a399971 [-t threads] [-s split] [-o order] [-u first:last] [-x K]
 *                [-e secs [-m subs] [-S seed]] [-p] [-P] [-q] n [n2]
 */
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef uint64_t u64;
#define BIT(i) (1ULL << (i))
#define MAXN 8

enum { UPB = 1, RTB = 2, DNB = 4, LTB = 8 };

/* -DXCHECK builds a variant for cross-checking: canonical forms compare the
   cells in reverse processing order (so symmetry is only tested at the leaves),
   and the original witness search is used with its heuristics off.  It must
   give the same counts. */
#ifdef XCHECK
enum { XCHK = 1 };
#else
enum { XCHK = 0 };
#endif

static int N, NN;
static u64 FULL, COL0, COLL;
static u64 NB[64];
static unsigned char SRC[8][64];  /* cell mapped onto cell j by symmetry g */
static unsigned char DMAP[8][16]; /* image of a set of directions under g */

static inline int ctz(u64 x) { return __builtin_ctzll(x); }
static inline int popc(u64 x) { return __builtin_popcountll(x); }

static inline u64 nbrs(u64 m)
{
    return (((m >> 1) & ~COLL) | ((m << 1) & ~COL0) | (m >> N) | (m << N)) & FULL;
}

/* number of grid-adjacent pairs (x, y) with x in a, y in b */
static inline int adjcount(u64 a, u64 b)
{
    return popc((a >> 1) & ~COLL & b) + popc((a << 1) & ~COL0 & b) +
           popc((a >> N) & b) + popc((a << N) & b);
}

static void sym_point(int g, int r, int c, int *R, int *C)
{
    const int m = N - 1;
    switch (g) {
    case 0: *R = r;     *C = c;     break;
    case 1: *R = c;     *C = m - r; break;
    case 2: *R = m - r; *C = m - c; break;
    case 3: *R = m - c; *C = r;     break;
    case 4: *R = r;     *C = m - c; break; /* left-right mirror */
    case 5: *R = m - r; *C = c;     break;
    case 6: *R = c;     *C = r;     break;
    default: *R = m - c; *C = m - r; break;
    }
}

static void init(int n)
{
    static const int dr[4] = {-1, 0, 1, 0}, dc[4] = {0, 1, 0, -1}; /* UP RT DN LT */
    N = n;
    NN = n * n;
    FULL = NN == 64 ? ~0ULL : BIT(NN) - 1;
    COL0 = COLL = 0;
    for (int r = 0; r < n; r++) {
        COL0 |= BIT(r * n);
        COLL |= BIT(r * n + n - 1);
    }
    for (int i = 0; i < NN; i++) NB[i] = nbrs(BIT(i));
    for (int g = 0; g < 8; g++) {
        for (int i = 0; i < NN; i++) {
            int R, C;
            sym_point(g, i / n, i % n, &R, &C);
            SRC[g][R * n + C] = i;
        }
        int img[4];
        for (int d = 0; d < 4; d++) {
            int a, b, x, y;
            sym_point(g, 0, 0, &a, &b);
            sym_point(g, dr[d], dc[d], &x, &y);
            for (int e = 0; e < 4; e++)
                if (x - a == dr[e] && y - b == dc[e]) img[d] = e;
        }
        for (int m = 0; m < 16; m++) {
            DMAP[g][m] = 0;
            for (int d = 0; d < 4; d++)
                if (m >> d & 1) DMAP[g][m] |= 1 << img[d];
        }
    }
}

/* Cells are processed in the order ORDER; processing a cell decides its edges
   to the cells processed after it.  Covers are compared cell by cell in the
   order CMPORD, each cell as its set of edge directions, and a cover is kept
   only if no symmetric image is smaller.  READY[g][k] is the depth from which
   positions 0..k of the image under g can be compared, so non-canonical
   branches are cut as early as the order allows.  Order kinds (-o):
     0: row by row
     1: ring by ring from the outside, each ring orbit by orbit, so that all
        8 symmetries are tested as early as possible
     2: the outer ring orbit by orbit, then the inside row by row */
static int ORDER[64], POS[64], CMPORD[64], DELTA[4];
static u64 LATER[65]; /* cells processed at depth >= k */
static unsigned char READY[8][64], CHECKAT[65];
static int order_kind = 2;

static void make_order(int kind)
{
    unsigned char used[64] = {0};
    int k = 0;
    if (kind > 0)
        for (int d = 0; d < (kind == 1 ? (N + 1) / 2 : 1); d++)
            for (int j = d; j <= (N - 1) / 2; j++)
                for (int g = 0; g < 8; g++) {
                    int R, C;
                    sym_point(g, d, j, &R, &C);
                    if (!used[R * N + C]) {
                        used[R * N + C] = 1;
                        ORDER[k++] = R * N + C;
                    }
                }
    for (int i = 0; i < NN; i++)
        if (!used[i]) ORDER[k++] = i;
    for (k = 0; k < NN; k++) POS[ORDER[k]] = k;
    for (k = NN, LATER[NN] = 0; k > 0; k--) LATER[k - 1] = LATER[k] | BIT(ORDER[k - 1]);
    for (k = 0; k < NN; k++) CMPORD[k] = XCHK ? ORDER[NN - 1 - k] : ORDER[k];
    memset(CHECKAT, 0, sizeof CHECKAT);
    for (int g = 1; g < 8; g++)
        for (int m = 0, k = 0; k < NN; k++) {
            const int c = CMPORD[k];
            if (POS[c] > m) m = POS[c];
            if (POS[SRC[g][c]] > m) m = POS[SRC[g][c]];
            READY[g][k] = m + 1;
            CHECKAT[m + 1] = 1;
        }
    DELTA[0] = -N; /* UPB RTB DNB LTB */
    DELTA[1] = 1;
    DELTA[2] = N;
    DELTA[3] = -1;
}

/* 1 if some image of cover t is smaller than cover s */
static int has_smaller_image(const unsigned char *t, const unsigned char *s)
{
    for (int g = 0; g < 8; g++) {
        const unsigned char *src = SRC[g], *dm = DMAP[g];
        for (int j = 0; j < NN; j++) {
            const int a = dm[t[src[CMPORD[j]]]], b = s[CMPORD[j]];
            if (a != b) {
                if (a < b) return 1;
                break;
            }
        }
    }
    return 0;
}

/* ------------------------------------------------------------------------ */
/* witness search                                                           */

typedef struct {
    u64 occ;               /* terminals and cells used by partial paths */
    int k, nopen;
    int head[32], tgt[32]; /* each open pair grows from head towards tgt */
    u64 pth[32];           /* cells of each partial path, head included */
    unsigned char open[32];
    /* the cover being tested, and what the search found out about its puzzle */
    u64 H, V;
    const unsigned char *dir;
    u64 ppath[32]; /* cells of the cover's own path for each pair */
    int nsol, self, designated, branched;
    u64 nodes;
} Search;

/* a full linkage was found; extra = cells taken by the last pair p (or p < 0) */
static void found_full_paths(Search *S, const u64 *pth, int p, u64 extra)
{
    u64 H = 0, V = 0;
    for (int q = 0; q < S->k; q++) {
        u64 m = pth[q] | BIT(S->tgt[q]);
        if (q == p) m |= extra;
        H |= m & (m >> 1) & ~COLL; /* full linkages are induced */
        V |= m & (m >> N);
    }
    S->nsol++;
    if (H == S->H && V == S->V) {
        S->self = 1;
        return;
    }
    unsigned char t[64];
    for (int i = 0; i < NN; i++)
        t[i] = (int)(H >> i & 1) * RTB | (int)(V >> i & 1) * DNB |
               (i % N ? (int)(H >> (i - 1) & 1) * LTB : 0) |
               (i >= N ? (int)(V >> (i - N) & 1) * UPB : 0);
    if (has_smaller_image(t, S->dir)) S->designated = 0;
}

static void found_full(Search *S, int p, u64 extra) { found_full_paths(S, S->pth, p, extra); }

/* 1 if the open pairs can be completed so that some cell stays empty */
static int witness(Search *S)
{
    S->nodes++;
    const u64 fr = FULL & ~S->occ;
    if (!S->nopen) {
        if (fr) return 1;
        found_full(S, -1, 0);
        return 0;
    }
    if (S->nopen == 1) {
        int p = 0;
        while (!S->open[p]) p++;
        const u64 nt = NB[S->tgt[p]];
        u64 lay = NB[S->head[p]] & fr, vis = lay;
        int d = 1;
        while (!(lay & nt)) {
            if (!lay) return 0;
            lay = nbrs(lay) & fr & ~vis;
            vis |= lay;
            d++;
        }
        if (d < popc(fr)) return 1;
        found_full(S, p, fr); /* the free cells form the only completion */
        return 0;
    }

    /* moves that keep each path induced; a pair with one move is forced */
    u64 mv[32], bestmv = 0;
    int best = -1, bestc = 5;
    for (int p = 0; p < S->k; p++) {
        if (!S->open[p]) continue;
        const u64 hb = BIT(S->head[p]);
        u64 m = 0;
        for (u64 cand = NB[S->head[p]] & fr; cand; cand &= cand - 1) {
            const int x = ctz(cand);
            if ((NB[x] & S->pth[p]) == hb) m |= BIT(x);
        }
        const int c = popc(m);
        if (!c) return 0;
        mv[p] = m;
        if (c < bestc) {
            bestc = c;
            best = p;
            bestmv = m;
        }
    }

    if (XCHK || bestc > 1) {
        /* the interior of a path lies in one free component touching both ends */
        u64 comp[64];
        int nc = 0;
        for (u64 rem = fr; rem;) {
            u64 x = rem & -rem, q;
            do {
                q = x;
                x = (x | nbrs(x)) & fr;
            } while (x != q);
            comp[nc++] = x;
            rem &= ~x;
        }
        bestc = 5;
        for (int p = 0; p < S->k; p++) {
            if (!S->open[p]) continue;
            const u64 nh = NB[S->head[p]], nt = NB[S->tgt[p]];
            u64 u = 0;
            for (int c = 0; c < nc; c++)
                if ((comp[c] & nh) && (comp[c] & nt)) u |= comp[c];
            const u64 m = mv[p] & u;
            const int c = popc(m);
            if (!c) return 0;
            if (c < bestc) {
                bestc = c;
                best = p;
                bestmv = m;
            }
        }
    }

    const int p = best, oh = S->head[p];
    for (u64 m = bestmv; m; m &= m - 1) {
        const int x = ctz(m);
        const int fin = (NB[x] >> S->tgt[p]) & 1;
        S->occ |= BIT(x);
        S->pth[p] |= BIT(x);
        S->head[p] = x;
        if (fin) {
            S->open[p] = 0;
            S->nopen--;
        }
        const int r = witness(S);
        if (fin) {
            S->open[p] = 1;
            S->nopen++;
        }
        S->head[p] = oh;
        S->pth[p] &= ~BIT(x);
        S->occ &= ~BIT(x);
        if (r) return 1;
    }
    return 0;
}

/* Fast witness search (the default; -DXCHECK uses witness() above).  Nearly
   every step of a compact cover's search is forced, so forced moves are
   applied in a loop that only updates the pairs a move touches, and the small
   state is copied at the rare real branch points instead of being undone. */
typedef struct {
    u64 occ;
    uint32_t openm, fm; /* open pairs; open pairs with exactly one move */
    u64 hm, tm;     /* heads and targets of the open pairs */
    unsigned char head[32], hp[64]; /* hp: pair whose head is at a cell */
    u64 pth[32];
    u64 fb[32]; /* cells next to a path cell before the head: moving there makes a chord */
    u64 mv[32]; /* cells the head can move to keeping the path induced */
} WS;

/* pair p moves its head to x; 0 if some open pair is left without a route */
static inline int ws_move(Search *S, WS *w, int p, int x)
{
    const u64 xb = BIT(x);
    const int h = w->head[p];
    S->nodes++;
    w->occ |= xb;
    w->pth[p] |= xb;
    w->fb[p] |= NB[h];
    w->hm &= ~BIT(h);
    w->head[p] = x;
    const u64 fr = FULL & ~w->occ;
    if (NB[x] >> S->tgt[p] & 1) {
        w->openm &= ~(1u << p);
        w->fm &= ~(1u << p);
        w->tm &= ~BIT(S->tgt[p]);
    } else {
        const u64 m = NB[x] & fr & ~w->fb[p];
        if (!m) return 0;
        w->mv[p] = m;
        w->fm = (w->fm & ~(1u << p)) | (uint32_t)!(m & (m - 1)) << p;
        w->hm |= xb;
        w->hp[x] = p;
    }
    /* open pairs whose head is next to x lose x as a move */
    for (u64 h = NB[x] & w->hm; h; h &= h - 1) {
        const int q = w->hp[ctz(h)];
        if (w->mv[q] & xb) {
            const u64 m = w->mv[q] &= ~xb;
            if (!m) return 0;
            if (!(m & (m - 1))) w->fm |= 1u << q;
        }
    }
    /* open pairs whose target is next to x may have lost their last way in
       (the head of an open pair never touches its target) */
    for (u64 t = NB[x] & w->tm; t; t &= t - 1)
        if (!(NB[ctz(t)] & fr)) return 0;
    return 1;
}

/* Finish pair p along the rest of the cover's own path in one go, if every
   remaining step would be forced: the path walked so far is a prefix of the
   cover's path, its remaining cells are free, and every free cell next to a
   remaining head (the cells before the one touching the target) lies off the
   path only where it is already next to an earlier cell of the walk.
   Returns 1 if done, 0 if not applicable, -1 if done but some pair is stuck. */
static inline int ws_batch(Search *S, WS *w, int p)
{
    const u64 pp = S->ppath[p], tb = BIT(S->tgt[p]);
    if (w->pth[p] & ~pp) return 0;
    const u64 rem = pp & ~w->pth[p] & ~tb;
    if (rem & w->occ) return 0;
    const u64 heads = (BIT(w->head[p]) | rem) & ~NB[S->tgt[p]];
    if (nbrs(heads) & ~w->occ & ~pp & ~w->fb[p] & FULL) return 0;
    S->nodes += popc(rem);
    w->occ |= rem;
    w->pth[p] |= rem;
    w->hm &= ~BIT(w->head[p]);
    w->openm &= ~(1u << p);
    w->fm &= ~(1u << p);
    w->tm &= ~tb;
    const u64 fr = FULL & ~w->occ, nn = nbrs(rem);
    for (u64 h = nn & w->hm; h; h &= h - 1) {
        const int q = w->hp[ctz(h)];
        if (w->mv[q] & rem) {
            const u64 m = w->mv[q] &= ~rem;
            if (!m) return -1;
            if (!(m & (m - 1))) w->fm |= 1u << q;
        }
    }
    for (u64 t = nn & w->tm; t; t &= t - 1)
        if (!(NB[ctz(t)] & fr)) return -1;
    return 1;
}

/* 1 if the open pairs can be completed so that some cell stays empty */
static int ws_search(Search *S, WS *w)
{
    for (;;) {
        int f = -1; /* forced moves while >= 2 pairs are open, same pair first */
        while (w->openm & (w->openm - 1) && w->fm) {
            if (f < 0 || !(w->fm >> f & 1)) {
                f = __builtin_ctz(w->fm);
                const int b = ws_batch(S, w, f);
                if (b < 0) return 0;
                if (b) {
                    f = -1;
                    continue;
                }
            }
            if (!ws_move(S, w, f, ctz(w->mv[f]))) return 0;
        }
        const u64 fr = FULL & ~w->occ;
        if (!w->openm) {
            if (fr) return 1;
            if (S->branched) found_full_paths(S, w->pth, -1, 0);
            else S->nsol++, S->self = 1; /* an unbranched search can only meet the cover */
            return 0;
        }
        if (!(w->openm & (w->openm - 1))) { /* last pair: a shortest path decides */
            const int p = __builtin_ctz(w->openm);
            if (fr == (S->ppath[p] & ~w->pth[p] & ~BIT(S->tgt[p]))) {
                /* the free cells are the rest of the pair's own induced path,
                   which is then its only route */
                if (S->branched) found_full_paths(S, w->pth, p, fr);
                else S->nsol++, S->self = 1;
                return 0;
            }
            const u64 nt = NB[S->tgt[p]];
            u64 lay = NB[w->head[p]] & fr, vis = lay;
            int d = 1;
            while (!(lay & nt)) {
                if (!lay) return 0;
                lay = nbrs(lay) & fr & ~vis;
                vis |= lay;
                d++;
            }
            if (d < popc(fr)) return 1;
            found_full_paths(S, w->pth, p, fr);
            return 0;
        }
        /* the interior of a path lies in one free component touching both ends */
        u64 comp[64];
        int nc = 0;
        for (u64 rem = fr; rem;) {
            u64 x = rem & -rem, q;
            do {
                q = x;
                x = (x | nbrs(x)) & fr;
            } while (x != q);
            comp[nc++] = x;
            rem &= ~x;
        }
        int best = -1, bestc = 99;
        for (uint32_t o = w->openm; o; o &= o - 1) {
            const int q = __builtin_ctz(o);
            const u64 nh = NB[w->head[q]], nt = NB[S->tgt[q]];
            u64 u = 0;
            for (int c = 0; c < nc; c++)
                if ((comp[c] & nh) && (comp[c] & nt)) u |= comp[c];
            const u64 m = w->mv[q] & u;
            if (!m) return 0;
            w->mv[q] = m;
            if (!(m & (m - 1))) w->fm |= 1u << q;
            const int c = popc(m);
            if (c < bestc) {
                bestc = c;
                best = q;
            }
        }
        if (bestc == 1) {
            if (!ws_move(S, w, best, ctz(w->mv[best]))) return 0;
            continue;
        }
        S->branched = 1;
        for (u64 m = w->mv[best]; m; m &= m - 1) {
            WS w2 = *w;
            if (ws_move(S, &w2, best, ctz(m)) && ws_search(S, &w2)) return 1;
        }
        return 0;
    }
}

/* ------------------------------------------------------------------------ */
/* enumeration of induced path covers                                       */

typedef struct {
    int deg[64], oth[64]; /* oth: other endpoint of the path (endpoints only) */
    u64 pm[64];           /* cells of the path (endpoints only) */
    u64 H, V;             /* chosen edges right of / below each cell */
    u64 ends;             /* cells of degree 1 */
    unsigned char dir[64];
    unsigned tied;          /* symmetries whose image equals the cover so far */
    unsigned char cpos[8];  /* first position not yet compared, per symmetry */
    int stop;               /* dfs depth at which a branch is complete */
    unsigned char opt[64];
    u64 covers, canon, compact, labeled, puzzles, multi, errors;
    u64 links[33]; /* compact solutions (reduced) by number of links */
    u64 xseen[2], xchecked, xmismatch;
    /* -e: reservoir of sampled sub-unit prefixes, sums of per-unit estimates */
    u64 rng, seen2;
    unsigned char res[8][64];
    double ek, ea, ea2, ev, ec, et, et2;
    Search S;
} Ctx;

typedef struct {
    int a, b, oa, ob;
    u64 ma, mb;
} Undo;

static inline int join(Ctx *C, int u, int v, Undo *U)
{
    if (C->oth[u] == v) return 0;
    const u64 A = C->pm[u], B = C->pm[v];
    if (adjcount(A, B) != 1) return 0;
    const int a = C->oth[u], b = C->oth[v];
    U->a = a;
    U->b = b;
    U->oa = C->oth[a];
    U->ob = C->oth[b];
    U->ma = C->pm[a];
    U->mb = C->pm[b];
    C->oth[a] = b;
    C->oth[b] = a;
    C->pm[a] = C->pm[b] = A | B;
    C->deg[u]++;
    C->deg[v]++;
    C->ends ^= BIT(u) | BIT(v); /* degrees 0->1 and 1->2 both toggle */
    return 1;
}

static inline void unjoin(Ctx *C, int u, int v, const Undo *U)
{
    C->deg[u]--;
    C->deg[v]--;
    C->ends ^= BIT(u) | BIT(v);
    C->oth[U->a] = U->oa;
    C->oth[U->b] = U->ob;
    C->pm[U->a] = U->ma;
    C->pm[U->b] = U->mb;
}

static void reset(Ctx *C)
{
    for (int i = 0; i < NN; i++) {
        C->deg[i] = 0;
        C->oth[i] = i;
        C->pm[i] = BIT(i);
        C->dir[i] = 0;
    }
    C->H = C->V = 0;
    C->ends = 0;
    C->tied = 0xFE;
    memset(C->cpos, 0, sizeof C->cpos);
}

static inline void toggle_edge(Ctx *C, int x, int y)
{
    const int a = x < y ? x : y, b = x < y ? y : x;
    if (b - a == 1) {
        C->H ^= BIT(a);
        C->dir[a] ^= RTB;
        C->dir[b] ^= LTB;
    } else {
        C->V ^= BIT(a);
        C->dir[a] ^= DNB;
        C->dir[b] ^= UPB;
    }
}

static inline int dirbit(int x, int y)
{
    return y == x - N ? UPB : y == x + 1 ? RTB : y == x + N ? DNB : LTB;
}

/* compare the images on the positions that became comparable at depth D;
   0 if some image is smaller (the branch is not canonical) */
static int symcheck(Ctx *C, int D)
{
    for (unsigned t = C->tied; t; t &= t - 1) {
        const int g = __builtin_ctz(t);
        int k = C->cpos[g];
        while (k < NN && READY[g][k] <= D) {
            const int c = CMPORD[k];
            const int a = DMAP[g][C->dir[SRC[g][c]]], b = C->dir[c];
            if (a != b) {
                if (a < b) return 0;
                C->tied &= ~(1u << g);
                break;
            }
            k++;
        }
        C->cpos[g] = (unsigned char)k;
    }
    return 1;
}

static int is_compact(Ctx *C)
{
    Search *S = &C->S;
    if (!XCHK) { /* one pass: each pair once, each grown from its less free end */
        const u64 occ = C->ends, fr = FULL & ~occ;
        WS w;
        w.occ = occ;
        w.openm = w.fm = 0;
        w.hm = w.tm = 0;
        int k = 0;
        for (u64 m = occ; m; k++) {
            const int e = ctz(m), o = C->oth[e];
            m &= ~(BIT(e) | BIT(o));
            const u64 ne = NB[e] & fr, no = NB[o] & fr;
            const int sw = popc(no) < popc(ne);
            const int a = sw ? o : e, b = sw ? e : o;
            const u64 mv = sw ? no : ne;
            const uint32_t open = !(NB[a] >> b & 1);
            S->tgt[k] = b;
            S->ppath[k] = C->pm[e];
            w.head[k] = a;
            w.hp[a] = k;
            w.pth[k] = BIT(a);
            w.fb[k] = 0;
            w.mv[k] = mv;
            w.openm |= open << k;
            w.fm |= (open & !(mv & (mv - 1))) << k;
            w.hm |= open ? BIT(a) : 0;
            w.tm |= open ? BIT(b) : 0;
        }
        S->k = k;
        S->H = C->H;
        S->V = C->V;
        S->dir = C->dir;
        S->nsol = S->self = S->branched = 0;
        S->designated = 1;
        return !ws_search(S, &w);
    }
    S->occ = 0;
    S->k = S->nopen = 0;
    for (u64 m = C->ends; m; m &= m - 1) {
        const int e = ctz(m);
        if (C->oth[e] > e) {
            S->head[S->k] = e;
            S->tgt[S->k] = C->oth[e];
            S->ppath[S->k] = C->pm[e];
            S->occ |= BIT(e) | BIT(C->oth[e]);
            S->k++;
        }
    }
    for (int p = 0; p < S->k; p++) {
        const int a = S->head[p], b = S->tgt[p];
        S->open[p] = 0;
        if (!(NB[a] & BIT(b))) { /* adjacent ends form a forced domino */
            if (!XCHK && popc(NB[b] & ~S->occ) < popc(NB[a] & ~S->occ)) {
                S->head[p] = b;
                S->tgt[p] = a;
            }
            S->open[p] = 1;
            S->nopen++;
        }
        S->pth[p] = BIT(S->head[p]);
    }
    S->H = C->H;
    S->V = C->V;
    S->dir = C->dir;
    S->nsol = S->self = S->branched = 0;
    S->designated = 1;
    if (XCHK) return !witness(S);
    WS w;
    w.occ = S->occ;
    w.openm = w.fm = 0;
    w.hm = w.tm = 0;
    for (int p = 0; p < S->k; p++) {
        w.head[p] = S->head[p];
        w.pth[p] = S->pth[p];
        if (S->open[p]) {
            w.openm |= 1u << p;
            w.mv[p] = NB[S->head[p]] & ~S->occ;
            if (!w.mv[p]) return 1; /* cannot happen for a cover: its own path leaves */
            if (!(w.mv[p] & (w.mv[p] - 1))) w.fm |= 1u << p;
            w.hm |= BIT(S->head[p]);
            w.hp[S->head[p]] = p;
            w.tm |= BIT(S->tgt[p]);
        }
    }
    return !ws_search(S, &w);
}

static int print_solutions; /* bit 0: draw the solution, bit 1: its starting position */
static pthread_mutex_t print_lock = PTHREAD_MUTEX_INITIALIZER;

static void print_cover(const Ctx *C)
{
    int lab[64], k = 0;
    for (int i = 0; i < NN; i++) lab[i] = 0;
    for (int i = 0; i < NN; i++) {
        if (lab[i] || C->deg[i] != 1) continue;
        k++;
        for (int x = i, prev = -1;;) {
            lab[x] = k;
            int nx = -1;
            if ((C->dir[x] & RTB) && x + 1 != prev) nx = x + 1;
            else if ((C->dir[x] & LTB) && x - 1 != prev) nx = x - 1;
            else if ((C->dir[x] & DNB) && x + N != prev) nx = x + N;
            else if ((C->dir[x] & UPB) && x - N != prev) nx = x - N;
            if (nx < 0) break;
            prev = x;
            x = nx;
        }
    }
    if (print_solutions & 2) /* end-point grid, laid out as in a399971.txt */
        for (int r = 0; r < N; r++) {
            for (int c = 0; c < N; c++) printf("%4d", C->deg[r * N + c] == 1 ? lab[r * N + c] : 0);
            printf(r + 1 < N ? "\n.\n" : "\n\n");
        }
    if (!(print_solutions & 1)) return;
    const int w = k < 10 ? 1 : 2;
    for (int r = 0; r < N; r++) {
        for (int c = 0; c < N; c++) {
            const int i = r * N + c;
            printf("%*d%s", w, lab[i], c + 1 == N ? "" : (C->dir[i] & RTB) ? "---" : "   ");
        }
        printf("\n");
        if (r + 1 < N) {
            for (int c = 0; c < N; c++)
                printf("%*s%s", w, (C->dir[r * N + c] & DNB) ? "|" : "", c + 1 < N ? "   " : "");
            printf("\n");
        }
    }
    printf("\n");
}

/* -x: independent exhaustive routing (every simple path, pairs in index
   order, plain reachability pruning); 1 if a routing leaves a cell empty */
typedef struct {
    int np, s[32], t[32];
    unsigned char used[64];
} Naive;

static int nv_reach(const Naive *V, int s, int t)
{
    int q[64], h = 0, e = 0;
    unsigned char seen[64] = {0};
    q[e++] = s;
    seen[s] = 1;
    while (h < e) {
        const int x = q[h++], r = x / N, c = x % N;
        const int nb[4] = {r ? x - N : -1, c + 1 < N ? x + 1 : -1, r + 1 < N ? x + N : -1,
                           c ? x - 1 : -1};
        for (int d = 0; d < 4; d++) {
            const int y = nb[d];
            if (y < 0) continue;
            if (y == t) return 1;
            if (!seen[y] && !V->used[y]) {
                seen[y] = 1;
                q[e++] = y;
            }
        }
    }
    return 0;
}

static int nv_route(Naive *V, int j);

static int nv_extend(Naive *V, int j, int x)
{
    const int r = x / N, c = x % N;
    const int nb[4] = {r ? x - N : -1, c + 1 < N ? x + 1 : -1, r + 1 < N ? x + N : -1,
                       c ? x - 1 : -1};
    for (int d = 0; d < 4; d++) {
        const int y = nb[d];
        if (y < 0) continue;
        if (y == V->t[j]) {
            if (nv_route(V, j + 1)) return 1;
        } else if (!V->used[y]) {
            V->used[y] = 1;
            const int f = nv_extend(V, j, y);
            V->used[y] = 0;
            if (f) return 1;
        }
    }
    return 0;
}

static int nv_route(Naive *V, int j)
{
    if (j == V->np) {
        for (int i = 0; i < NN; i++)
            if (!V->used[i]) return 1;
        return 0;
    }
    for (int i = j; i < V->np; i++)
        if (!nv_reach(V, V->s[i], V->t[i])) return 0;
    return nv_extend(V, j, V->s[j]);
}

static u64 xcheck_every;

static void xcheck(Ctx *C, int compact)
{
    /* sample non-compact covers (rare) 100 times more often */
    const u64 every = compact ? xcheck_every : (xcheck_every + 99) / 100;
    if (++C->xseen[compact] % every) return;
    Naive V;
    V.np = 0;
    memset(V.used, 0, sizeof V.used);
    for (int e = 0; e < NN; e++)
        if (C->deg[e] == 1 && C->oth[e] > e) {
            V.s[V.np] = e;
            V.t[V.np] = C->oth[e];
            V.used[e] = V.used[C->oth[e]] = 1;
            V.np++;
        }
    C->xchecked++;
    if ((!nv_route(&V, 0)) != compact) {
        C->xmismatch++;
        pthread_mutex_lock(&print_lock);
        printf("MISMATCH: fast search says %scompact, exhaustive routing disagrees\n",
               compact ? "" : "not ");
        print_cover(C);
        pthread_mutex_unlock(&print_lock);
    }
}

static void leaf(Ctx *C)
{
    C->covers++; /* every leaf is canonical: the others were cut on the way */
    const int st = 1 + __builtin_popcount(C->tied);
    C->canon++;
    const int compact = is_compact(C);
    if (xcheck_every) xcheck(C, compact);
    if (!compact) return;
    C->compact++;
    C->links[popc(C->ends) / 2]++;
    C->labeled += 8 / st;
    C->puzzles += C->S.designated;
    C->multi += C->S.nsol > 1;
    C->errors += !C->S.self; /* the search must meet the cover itself */
    if (print_solutions) {
        pthread_mutex_lock(&print_lock);
        print_cover(C);
        pthread_mutex_unlock(&print_lock);
    }
}

/* task prefixes: SPLIT option codes (directions of the edges each cell adds) */
static int SPLIT;
/* -e estimate mode: sub-unit depth, sub-units per unit, time budget, seed */
static int SPLIT2, est_sub = 4;
static int max_links = 64; /* -L: only solutions with at most this many links */
static double est_budget;
static u64 est_seed = 1, est_units; /* -E: fixed number of units instead of a time budget */
static unsigned char *tasks;
static size_t ntasks, captasks;

static void record(Ctx *C)
{
    if (ntasks == captasks) {
        captasks = captasks ? 2 * captasks : 1024;
        tasks = realloc(tasks, captasks * SPLIT);
        if (!tasks) {
            perror("realloc");
            exit(1);
        }
    }
    memcpy(tasks + ntasks * SPLIT, C->opt, SPLIT);
    ntasks++;
}

static inline u64 rnd(u64 *x)
{
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    return *x;
}

/* -e: keep a uniform sample of est_sub sub-unit prefixes (reservoir sampling) */
static void reservoir(Ctx *C)
{
    const u64 seen = ++C->seen2;
    const u64 slot = seen <= (u64)est_sub ? seen - 1 : rnd(&C->rng) % seen;
    if (slot < (u64)est_sub) memcpy(C->res[slot], C->opt, SPLIT2);
}

static void dfs(Ctx *C, int k, int gen)
{
    const unsigned tied0 = C->tied;
    unsigned char cpos0[8];
    const int chk = CHECKAT[k] && tied0;
    if (chk) {
        memcpy(cpos0, C->cpos, sizeof cpos0);
        if (!symcheck(C, k)) goto out;
    }
    if (max_links < 32 && popc(C->ends & ~LATER[k]) > 2 * max_links) goto out; /* terminals */
    if (k == C->stop) {
        if (gen == 1) record(C);
        else if (gen == 2) reservoir(C);
        else leaf(C);
        goto out;
    }
    const int x = ORDER[k], d0 = C->deg[x];
    const u64 later = LATER[k + 1];
    int nb[4], nn = 0;
    for (u64 m = NB[x] & later; m; m &= m - 1) nb[nn++] = ctz(m);
    for (int sub = 0; sub < (1 << nn); sub++) {
        const int sz = __builtin_popcount(sub);
        if (d0 + sz < 1 || d0 + sz > 2) continue;
        Undo u[2];
        int y[2], ny = 0, ok = 1;
        for (int j = 0; j < nn && ok; j++)
            if (sub >> j & 1) {
                if (join(C, x, nb[j], &u[ny])) y[ny++] = nb[j];
                else ok = 0;
            }
        /* a later cell left with no edge and no later neighbour is dead */
        for (int j = 0; j < nn && ok; j++)
            if (!C->deg[nb[j]] && !(NB[nb[j]] & later)) ok = 0;
        if (ok) {
            unsigned char o = 0;
            for (int j = 0; j < ny; j++) {
                toggle_edge(C, x, y[j]);
                o |= dirbit(x, y[j]);
            }
            C->opt[k] = o;
            dfs(C, k + 1, gen);
            for (int j = 0; j < ny; j++) toggle_edge(C, x, y[j]);
        }
        while (ny) {
            ny--;
            unjoin(C, x, y[ny], &u[ny]);
        }
    }
out:
    if (chk) {
        C->tied = tied0;
        memcpy(C->cpos, cpos0, sizeof cpos0);
    }
}

static void replay(Ctx *C, const unsigned char *opt, int len)
{
    for (int k = 0; k < len; k++) {
        const int x = ORDER[k];
        C->opt[k] = opt[k];
        for (int d = 0; d < 4; d++)
            if (opt[k] >> d & 1) {
                Undo u;
                join(C, x, x + DELTA[d], &u);
                toggle_edge(C, x, x + DELTA[d]);
            }
        if (k + 1 < len && CHECKAT[k + 1] && C->tied) symcheck(C, k + 1);
    }
}

static atomic_size_t next_task, done_tasks;
static size_t unit_first, unit_last = (size_t)-1; /* -u: work units [first, last) */
static int quiet;
static double t_start;

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void *worker(void *arg)
{
    Ctx *C = arg;
    C->stop = NN;
    for (;;) {
        const size_t t = atomic_fetch_add(&next_task, 1);
        if (t >= unit_last) break;
        reset(C);
        replay(C, tasks + t * SPLIT, SPLIT);
        dfs(C, SPLIT, 0);
        const size_t d = atomic_fetch_add(&done_tasks, 1) + 1, all = unit_last - unit_first;
        if (!quiet && (d % (all / 200 + 1) == 0 || d == all)) {
            const double el = now() - t_start;
            fprintf(stderr, "\r  n=%d: %zu/%zu tasks, %.0fs elapsed, ~%.0fs left   ", N, d, all,
                    el, el * (all - d) / d);
        }
    }
    return NULL;
}

/* -e: take work units in random order until the time budget is spent; the
   count of each unit is estimated from est_sub random sub-units of depth
   SPLIT2 (two-stage sampling, unbiased) */
static void *est_worker(void *arg)
{
    Ctx *C = arg;
    for (;;) {
        if (!est_units && now() - t_start >= est_budget) break;
        const size_t t = atomic_fetch_add(&next_task, 1);
        if (t >= ntasks || (est_units && t >= est_units)) break;
        C->rng = (est_seed + 1) * 0x9E3779B97F4A7C15ULL ^ (t * 0xBF58476D1CE4E5B9ULL + 1);
        if (!C->rng) C->rng = 1;
        reset(C);
        replay(C, tasks + t * SPLIT, SPLIT);
        C->seen2 = 0;
        C->stop = SPLIT2;
        dfs(C, SPLIT, 2);
        C->stop = NN;
        const int m = C->seen2 < (u64)est_sub ? (int)C->seen2 : est_sub;
        double sa = 0, sv = 0, sc = 0, st = 0;
        for (int j = 0; j < m; j++) {
            const u64 a0 = C->compact, v0 = C->covers, c0 = C->canon;
            const double t0 = now();
            reset(C);
            replay(C, C->res[j], SPLIT2);
            dfs(C, SPLIT2, 0);
            st += now() - t0;
            sa += C->compact - a0;
            sv += C->covers - v0;
            sc += C->canon - c0;
        }
        const double w = m ? (double)C->seen2 / m : 0;
        C->ek += 1;
        C->ea += w * sa;
        C->ea2 += w * sa * w * sa;
        C->ev += w * sv;
        C->ec += w * sc;
        C->et += w * st;
        C->et2 += w * st * w * st;
    }
    return NULL;
}

static void estimate(int nthreads)
{
    SPLIT2 = SPLIT + 2 * N < NN ? SPLIT + 2 * N : NN;
    u64 x = est_seed ? est_seed : 1;
    unsigned char tmp[64];
    for (size_t i = ntasks; i > 1; i--) { /* shuffle the work units */
        const size_t j = rnd(&x) % i;
        memcpy(tmp, tasks + (i - 1) * SPLIT, SPLIT);
        memcpy(tasks + (i - 1) * SPLIT, tasks + j * SPLIT, SPLIT);
        memcpy(tasks + j * SPLIT, tmp, SPLIT);
    }
    Ctx *ctx = calloc(nthreads, sizeof(Ctx));
    pthread_t *th = malloc(nthreads * sizeof(pthread_t));
    atomic_store(&next_task, 0);
    for (int t = 0; t < nthreads; t++) {
        pthread_create(&th[t], NULL, est_worker, &ctx[t]);
    }
    double k = 0, a = 0, a2 = 0, v = 0, c = 0, et = 0, et2 = 0;
    for (int t = 0; t < nthreads; t++) {
        pthread_join(th[t], NULL);
        k += ctx[t].ek;
        a += ctx[t].ea;
        a2 += ctx[t].ea2;
        v += ctx[t].ev;
        c += ctx[t].ec;
        et += ctx[t].et;
        et2 += ctx[t].et2;
    }
    const double T = (double)ntasks, fpc = 1; /* conservative: ignores the sampled fraction */
    const double ma = a / k, mt = et / k;
    const double sa = k > 1 ? sqrt((a2 - k * ma * ma) / (k - 1) / k) * fpc : 0;
    const double stt = k > 1 ? sqrt((et2 - k * mt * mt) / (k - 1) / k) * fpc : 0;
    printf("estimate for n=%d from %.0f of %zu work units in random order (%d sub-units "
           "each, seed %llu, %.0f s):\n", N, k, ntasks, est_sub, (unsigned long long)est_seed,
           now() - t_start);
    printf("  a(%d) ~ %.4g  (standard error %.1f%%)\n", N, T * ma, 100 * sa / ma);
    printf("  induced covers visited ~ %.3g (canonical ~ %.3g)\n", T * v / k, T * c / k);
    printf("  full run ~ %.3g thread-seconds = %.3g days with %d threads (standard error %.1f%%)\n",
           T * mt, T * mt / nthreads / 86400, nthreads, 100 * stt / mt);
    fflush(stdout);
    free(th);
    free(ctx);
}

static void run(int n, int nthreads, int split)
{
    init(n);
    make_order(order_kind);
    t_start = now();
    SPLIT = split > 0 ? split : n >= 6 ? 2 * n + 2 : n;
    if (SPLIT > NN) SPLIT = NN;
    ntasks = 0;

    Ctx *gen = calloc(1, sizeof(Ctx));
    reset(gen);
    gen->stop = SPLIT;
    dfs(gen, 0, 1);
    free(gen);

    if (nthreads < 1) nthreads = 1;
    if (est_budget > 0) {
        estimate(nthreads);
        return;
    }
    const size_t first = unit_first < ntasks ? unit_first : ntasks;
    const size_t last = unit_last < ntasks ? unit_last : ntasks;
    const int partial = first > 0 || last < ntasks || unit_first == unit_last;
    if (unit_first == unit_last) nthreads = 0; /* -u 0:0 only reports the unit count */
    unit_first = first;
    unit_last = last;
    Ctx *ctx = calloc(nthreads, sizeof(Ctx));
    pthread_t *th = malloc(nthreads * sizeof(pthread_t));
    atomic_store(&next_task, first);
    atomic_store(&done_tasks, 0);
    for (int t = 0; t < nthreads; t++) pthread_create(&th[t], NULL, worker, &ctx[t]);
    u64 covers = 0, canonc = 0, compact = 0, labeled = 0, puzzles = 0, multi = 0, errors = 0,
        nodes = 0, xchecked = 0, xmismatch = 0, links[33] = {0};
    for (int t = 0; t < nthreads; t++) {
        pthread_join(th[t], NULL);
        covers += ctx[t].covers;
        canonc += ctx[t].canon;
        compact += ctx[t].compact;
        labeled += ctx[t].labeled;
        puzzles += ctx[t].puzzles;
        multi += ctx[t].multi;
        errors += ctx[t].errors;
        xchecked += ctx[t].xchecked;
        xmismatch += ctx[t].xmismatch;
        nodes += ctx[t].S.nodes;
        for (int k = 0; k < 33; k++) links[k] += ctx[t].links[k];
    }
    if (!quiet && ntasks) fprintf(stderr, "\n");
    unit_first = 0;
    unit_last = (size_t)-1;
    if (partial) {
        /* one line, summed over all ranges by merge_ranges.py */
        printf("RANGE n=%d split=%d order=%d units=%zu first=%zu last=%zu compact=%llu labeled=%llu "
               "puzzles=%llu multi=%llu covers=%llu errors=%llu seconds=%.0f\n",
               n, SPLIT, order_kind, ntasks, first, last, (unsigned long long)compact,
               (unsigned long long)labeled, (unsigned long long)puzzles, (unsigned long long)multi,
               (unsigned long long)covers, (unsigned long long)errors, now() - t_start);
        if (xcheck_every)
            printf("  exhaustive-routing cross-check: %llu covers checked, %llu mismatches\n",
                   (unsigned long long)xchecked, (unsigned long long)xmismatch);
        fflush(stdout);
        free(th);
        free(ctx);
        return;
    }
    printf("a(%d) = %llu\n", n, (unsigned long long)compact);
    printf("  compact solutions, not reduced for symmetry: %llu\n", (unsigned long long)labeled);
    printf("  compact starting positions, reduced: %llu (%llu solutions share a puzzle with another)\n",
           (unsigned long long)puzzles, (unsigned long long)multi);
    printf("  induced path covers visited: %llu (%llu canonical), witness-search nodes: %llu\n",
           (unsigned long long)covers, (unsigned long long)canonc, (unsigned long long)nodes);
    printf("  %zu tasks, %d threads, %.2f s\n", ntasks, nthreads, now() - t_start);
    if (xcheck_every)
        printf("  exhaustive-routing cross-check: %llu covers checked, %llu mismatches\n",
               (unsigned long long)xchecked, (unsigned long long)xmismatch);
    printf("  by number of links:");
    for (int k = 0; k < 33; k++)
        if (links[k]) printf(" %d:%llu", k, (unsigned long long)links[k]);
    printf("\n");
    if (XCHK) printf("  (XCHECK variant build)\n");
    if (errors) printf("  INTERNAL ERROR: %llu searches missed their own cover\n",
                       (unsigned long long)errors);
    fflush(stdout);
    free(th);
    free(ctx);
}

int main(int argc, char **argv)
{
    int nthreads = (int)sysconf(_SC_NPROCESSORS_ONLN), split = 0, opt;
    while ((opt = getopt(argc, argv, "t:s:x:e:E:m:S:o:u:L:pPq")) != -1) {
        switch (opt) {
        case 't': nthreads = atoi(optarg); break;
        case 's': split = atoi(optarg); break;
        case 'p': print_solutions |= 1; break;
        case 'P': print_solutions |= 2; break;
        case 'q': quiet = 1; break;
        case 'x': xcheck_every = strtoull(optarg, NULL, 10); break;
        case 'e': est_budget = atof(optarg); break;
        case 'm': est_sub = atoi(optarg) < 1 ? 1 : atoi(optarg) > 8 ? 8 : atoi(optarg); break;
        case 'S': est_seed = strtoull(optarg, NULL, 10); break;
        case 'o': order_kind = atoi(optarg); break;
        case 'L': max_links = atoi(optarg); break;
        case 'u': {
            unsigned long long a = 0, b = 0;
            if (sscanf(optarg, "%llu:%llu", &a, &b) != 2 || a > b) goto usage;
            unit_first = a;
            unit_last = b;
            break;
        }
        case 'E': est_units = strtoull(optarg, NULL, 10); est_budget = 1; break;
        default: goto usage;
        }
    }
    if (optind >= argc) goto usage;
    const int n1 = atoi(argv[optind]);
    const int n2 = optind + 1 < argc ? atoi(argv[optind + 1]) : n1;
    if (n1 < 1 || n2 > MAXN || n1 > n2) {
        fprintf(stderr, "n must be in 1..%d\n", MAXN);
        return 1;
    }
    for (int n = n1; n <= n2; n++) run(n, nthreads, split);
    return 0;
usage:
    fprintf(stderr, "usage: %s [-t threads] [-s split] [-o order] [-u first:last] [-x K]\n"
                    "          [-e secs [-m subs] [-S seed]] [-p] [-P] [-q] n [n2]\n"
                    "  -t threads  worker threads (default: all cores)\n"
                    "  -s split    cells fixed per work unit (default: n, or 2n+2 for n >= 6)\n"
                    "  -o order    cell order: 0 rows, 1 rings by orbits, 2 outer ring by\n"
                    "              orbits then rows (default 2)\n"
                    "  -u a:b      only work units a..b-1, printed as one RANGE line; sum the\n"
                    "              lines of disjoint ranges with merge_ranges.py (-u 0:0: count)\n"
                    "  -x K        re-decide every K-th canonical cover (every K/100-th\n"
                    "              non-compact one) by exhaustive routing and compare\n"
                    "  -e secs     estimate a(n) and the full run time by sampling random\n"
                    "              work units for about secs seconds (-m sub-units per\n"
                    "              unit, default 4; -S random seed; -E K: exactly K units)\n"
                    "  -L max      only solutions with at most max links\n"
                    "  -p          draw every counted solution (use -t 1 for a fixed order)\n"
                    "  -P          print the starting position of every counted solution\n"
                    "  -q          no progress output\n", argv[0]);
    return 1;
}
