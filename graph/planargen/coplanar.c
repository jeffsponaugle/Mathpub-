/* coplanar.c -- PRUNE/PREPRUNE hooks that turn nauty's geng into a
 * generator of *complements* of connected planar graphs.
 *
 * Why complements?  geng grows a graph by adding a vertex that has
 * maximum degree in the new graph (canonical augmentation, ties broken by
 * canonical labelling).  For a hereditary class of sparse graphs such as
 * planar graphs this is the wrong way round: every planar graph has a
 * vertex of degree <= 5, so the *last* vertex would ideally be one of
 * minimum degree, and only neighbourhood subsets of size <= 6 would need
 * to be tried.  Generating the complement class instead achieves exactly
 * that: "maximum degree in the complement" == "minimum degree in the
 * planar graph".  Complements of planar graphs form a hereditary class,
 * so pruning by "complement is planar" is valid at every level.
 *
 * Compile geng with
 *    -DPREPRUNE=coplanar_preprune -DPRUNE=coplanar_prune
 * and link this file together with planarity.c.
 *
 * Usage (n vertices, planar graph with e edges <-> complement with
 * C(n,2)-e edges; a connected planar graph has n-1 <= e <= 3n-6):
 *    geng_coplanar -u -v n  [C(n,2)-(3n-6)]:[C(n,2)-(n-1)]  [res/mod]
 * The per-edge counts printed by -v are for the complements; edge count
 * k of the complement corresponds to C(n,2)-k edges of the planar graph.
 *
 * Statistics (calls, rejections, planarity tests) are reported on stderr
 * at exit through geng's SUMMARY hook when compiled with
 *    -DSUMMARY=coplanar_summary
 */

#include "gtools.h"
#include "planarity.h"
#ifdef USE_LR
#include "lrplanar.h"      /* -DUSE_LR: left-right test instead of nauty's */
#endif
#ifdef USE_LRSG
#include "lrplanar_sg.h"   /* -DUSE_LRSG: sparsegraph left-right test (regression testing) */
#endif

#if MAXN > 32
#error "coplanar.c is written for geng compiled with WORDSIZE=32, MAXN=32"
#endif

static nauty_counter np_calls, np_edge_rej, np_conn_rej, np_planar_tests,
                     np_nonplanar, np_trivial, np_deferred, np_bydeg[8],
                     pr_calls, pr_rej, pr_final;

/* Complement g (n vertices, nauty m=1 format) into gc; return edge count. */
static int
compl_bits(graph *g, int n, graph *gc)
{
    int i, ne = 0;
    setword all = ALLMASK(n);

    for (i = 0; i < n; ++i)
    {
        gc[i] = (g[i] ^ all) & ~bit[i];
        ne += POPCOUNT(gc[i]);
    }
    return ne / 2;
}

/* Number of connected components of gc (bitmask BFS). */
static int
ncomponents(graph *gc, int n)
{
    setword notvisited, queue;
    int nc = 0, i;

    notvisited = ALLMASK(n);
    while (notvisited)
    {
        ++nc;
        queue = SWHIBIT(notvisited);
        notvisited &= ~queue;
        while (queue)
        {
            TAKEBIT(i, queue);
            notvisited &= ~bit[i];
            queue |= gc[i] & notvisited;
        }
    }
    return nc;
}

#ifdef NBCACHE
/* Per-level cache of recent results, exploiting monotonicity in the new
 * vertex's neighbourhood N (McKay's suggestion, specialised): with the
 * parent fixed, parent + N' non-planar and N' subset of N implies parent +
 * N non-planar; parent + N'' planar and N subset of N'' implies parent + N
 * planar.  geng presents the children of one parent consecutively, so a
 * few recent neighbourhoods catch most siblings.  One cache per level
 * because the recursion interleaves levels.  With NBCACHE_MIN a non-planar
 * N is first shrunk greedily to a minimal non-planar neighbourhood (a few
 * extra tests) so that it catches more siblings. */
#ifndef NBCACHE_K
#define NBCACHE_K 16
#endif
typedef struct {
    graph parent[MAXN];        /* rows 0..n-2 of the complement, masked */
    boolean valid;
    setword np[NBCACHE_K]; int nnp, inp;     /* non-planar neighbourhoods */
    setword pl[NBCACHE_K]; int npl, ipl;     /* planar neighbourhoods */
} nbcache;
static nbcache cache[MAXN+1];
static nauty_counter cc_lookups, cc_hit_np, cc_hit_pl, cc_miss, cc_reset, cc_mintests;

static boolean isplanar_bitmask_raw(graph *gc, int n, int ne);

static boolean
isplanar_bitmask(graph *gc, int n, int ne)
{
    nbcache *ca = &cache[n];
    setword N = gc[n-1], mask = ALLMASK(n-1);
    int i;
    boolean ans;

    if (n <= 4) return TRUE;
    if (ne > 3*n - 6) return FALSE;
    ++cc_lookups;

    /* same parent as last time at this level? */
    for (i = 0; i < n-1; ++i)
        if ((gc[i] & mask) != ca->parent[i]) break;
    if (i < n-1 || !ca->valid)
    {
        for (i = 0; i < n-1; ++i) ca->parent[i] = gc[i] & mask;
        ca->valid = TRUE; ca->nnp = ca->npl = ca->inp = ca->ipl = 0;
        ++cc_reset;
    }
    else
    {
        for (i = 0; i < ca->nnp; ++i)
            if ((ca->np[i] & ~N) == 0) { ++cc_hit_np; return FALSE; }
        for (i = 0; i < ca->npl; ++i)
            if ((N & ~ca->pl[i]) == 0) { ++cc_hit_pl; return TRUE; }
    }
    ++cc_miss;
    ans = isplanar_bitmask_raw(gc, n, ne);
    if (ans)
    {
        ca->pl[ca->ipl] = N; ca->ipl = (ca->ipl + 1) % NBCACHE_K;
        if (ca->npl < NBCACHE_K) ++ca->npl;
    }
    else
    {
#ifdef NBCACHE_MIN
        /* Shrink N greedily while parent + N stays non-planar, on a private
           copy of the graph.  A single neighbour never causes non-planarity,
           so stop at two. */
        graph h[MAXN];
        setword rest = N;
        int j, m = ne;
        for (i = 0; i < n; ++i) h[i] = gc[i];
        while (rest && POPCOUNT(N) > 2)
        {
            TAKEBIT(j, rest);
            h[n-1] = N & ~bit[j]; h[j] &= ~bit[n-1];
            ++cc_mintests;
            if (!isplanar_bitmask_raw(h, n, m - 1)) { N &= ~bit[j]; --m; }
            else h[j] |= bit[n-1];
        }
#endif
        ca->np[ca->inp] = N; ca->inp = (ca->inp + 1) % NBCACHE_K;
        if (ca->nnp < NBCACHE_K) ++ca->nnp;
    }
    return ans;
}

static boolean
isplanar_bitmask_raw(graph *gc, int n, int ne)
{
#else
/* Planarity test of gc via nauty's Boyer-Myrvold style tester. */
static boolean
isplanar_bitmask(graph *gc, int n, int ne)
{
#endif
    t_ver_sparse_rep V[MAXN];
    t_adjl_sparse_rep A[MAXN*(MAXN-1) + 1];
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    int edge_pos, v, w, c;
    int i, j, k;
    setword gv;
    boolean ans;

    if (n <= 4) return TRUE;             /* K4 and smaller are planar */
    if (ne > 3*n - 6) return FALSE;      /* Euler bound */
#ifdef USE_LR
    ++np_planar_tests;
    ans = lr_isplanar(gc, n);
    if (!ans) ++np_nonplanar;
    return ans;
#endif
#ifdef USE_LRSG
    ++np_planar_tests;
    ans = lrplanar_dense(gc, 1, n, FALSE);
    if (!ans) ++np_nonplanar;
    return ans;
#endif

    k = 0;
    for (i = 0; i < n; ++i)
    {
        gv = gc[i];
        if (gv == 0)
            V[i].first_edge = NIL;
        else
        {
            V[i].first_edge = k;
            while (gv)
            {
                TAKEBIT(j, gv);
                A[k].end_vertex = j;
                A[k].next = k + 1;
                ++k;
            }
            A[k-1].next = NIL;
        }
    }

    ++np_planar_tests;
    arena_reset();      /* all memory of the previous call is released */
    ans = sparseg_adjl_is_planar(V, n, A, &c, &dfs_tree, &back_edges,
                                 &mult_edges, &embed_graph, &edge_pos, &v, &w);
    /* no explicit frees: everything lives in the arena */

    if (!ans) ++np_nonplanar;
    return ans;
}

/* PREPRUNE: called for every candidate child (n vertices) *before* the
 * canonicity test.  Reject (return 1) unless the complement is planar,
 * and at the final level also connected. */
int
coplanar_preprune(graph *g, int n, int maxn)
{
    graph gc[MAXN];
    int ne;

    ++np_calls;
    ne = compl_bits(g, n, gc);

    if (n >= 3 && ne > 3*n - 6) { ++np_edge_rej; return 1; }

    if (n == maxn && ncomponents(gc, n) != 1) { ++np_conn_rej; return 1; }

    /* The parent (vertices 0..n-2) already passed this test, so its
       complement is planar.  Adding a vertex of degree <= 1, or of degree 2
       whose two neighbours are adjacent, keeps it planar. */
    {
        setword nb = gc[n-1];
        int d = POPCOUNT(nb);
        ++np_bydeg[d < 8 ? d : 7];
        if (d <= 1) { ++np_trivial; return 0; }
        if (d == 2)
        {
            int u, w;
            TAKEBIT(u, nb);
            TAKEBIT(w, nb);
            if (gc[u] & bit[w]) { ++np_trivial; return 0; }
        }
    }

    /* At the final level the canonicity test (geng's accept2) is cheaper
       than a planarity test and rejects about half of the candidates, so
       defer the planarity test to PRUNE, which geng calls after accept2. */
    if (n == maxn) { ++np_deferred; return 0; }

    return !isplanar_bitmask(gc, n, ne);
}

/* PRUNE: called on each accepted graph before it is extended, and on final
 * graphs after the canonicity test.  At the final level this is where the
 * planarity test happens (see coplanar_preprune).  At level maxn-1 a parent
 * whose complement has an isolated vertex, or more components than
 * (minimum degree + 1), cannot be completed to a connected graph by one
 * more minimum-degree vertex. */
int
coplanar_prune(graph *g, int n, int maxn)
{
    graph gc[MAXN];
    int i, d, mindeg, nc, ne;

    if (n == maxn)
    {
        setword nb;

        ++pr_final;
        ne = compl_bits(g, n, gc);
        nb = gc[n-1];
        d = POPCOUNT(nb);
        if (d <= 1) return 0;
        if (d == 2)
        {
            int u, w;
            TAKEBIT(u, nb);
            TAKEBIT(w, nb);
            if (gc[u] & bit[w]) return 0;
        }
        return !isplanar_bitmask(gc, n, ne);
    }

    if (n != maxn - 1 || n < 2) return 0;

    ++pr_calls;
    compl_bits(g, n, gc);
    mindeg = n;
    for (i = 0; i < n; ++i)
    {
        d = POPCOUNT(gc[i]);
        if (d < mindeg) mindeg = d;
    }
    if (mindeg == 0) { ++pr_rej; return 1; }
    nc = ncomponents(gc, n);
    if (nc > mindeg + 1) { ++pr_rej; return 1; }
    return 0;
}

void
coplanar_summary(nauty_counter nout, double cpu)
{
    fprintf(stderr, ">C preprune calls=" COUNTER_FMT
            " edge_rej=" COUNTER_FMT " conn_rej=" COUNTER_FMT
            " planarity_tests=" COUNTER_FMT " nonplanar=" COUNTER_FMT
            " | prune(maxn-1) calls=" COUNTER_FMT " rej=" COUNTER_FMT "\n",
            np_calls, np_edge_rej, np_conn_rej, np_planar_tests,
            np_nonplanar, pr_calls, pr_rej);
    fprintf(stderr, ">C final-level candidates deferred to PRUNE=" COUNTER_FMT
            ", of which canonical (planarity-tested)=" COUNTER_FMT "\n",
            np_deferred, pr_final);
    fprintf(stderr, ">C trivially planar (skipped test)=" COUNTER_FMT
            "  new-vertex degree histogram 0..6,7+:", np_trivial);
    { int i; for (i = 0; i < 8; ++i) fprintf(stderr, " " COUNTER_FMT, np_bydeg[i]); }
    fprintf(stderr, "\n");
#ifdef NBCACHE
    fprintf(stderr, ">C nbcache(K=%d): lookups=" COUNTER_FMT " hits nonplanar=" COUNTER_FMT
            " planar=" COUNTER_FMT " misses(full tests)=" COUNTER_FMT " parent changes=" COUNTER_FMT
            " minimisation tests=" COUNTER_FMT "\n", NBCACHE_K, cc_lookups, cc_hit_np, cc_hit_pl,
            cc_miss, cc_reset, cc_mintests);
#endif
    if (nout > 0)
        fprintf(stderr, ">C %.3f us cpu per output graph, %.3f us per preprune call\n",
                1e6*cpu/(double)nout, 1e6*cpu/(double)(np_calls ? np_calls : 1));
}
