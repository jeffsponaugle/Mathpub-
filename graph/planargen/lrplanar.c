/* lrplanar.c -- Left-Right planarity test for small simple graphs.
 *
 * Algorithm: the left-right (LR) planarity test of de Fraysseix and
 * Rosenstiehl, in the formulation of U. Brandes, "The Left-Right Planarity
 * Test" (2009).  The networkx function check_planarity() is another
 * implementation of the same pseudo-code and was used as a cross-reference
 * while writing this one.  Only the decision "planar or not" is computed;
 * the embedding phase of the algorithm is omitted.
 *
 * Input: an undirected simple graph with n <= WORDSIZE vertices in nauty's
 * bitmask format (g[v] is the set of neighbours of v, one bit per vertex).
 * Output: 1 if planar, 0 if not.  Everything lives in static arrays, no
 * memory is allocated, and the function is not thread-safe.  It was
 * validated against nauty's own planarity tester (planarity.c) with zero
 * disagreements on all graphs with <= 9 vertices, all 11.7 million
 * connected graphs with 10 vertices, and 12 million random graphs with
 * 11..16 vertices (see lrtest.c).
 *
 * ---------------------------------------------------------------------
 * How the test works
 * ---------------------------------------------------------------------
 *
 * 1. Orient the graph by a depth-first search.  Every edge becomes either
 *    a TREE edge (parent -> child in the DFS tree) or a BACK edge (from a
 *    vertex up to one of its proper ancestors).  A DFS tree of an
 *    undirected graph has no other kind of edge.  height[v] is the depth
 *    of v in the tree (0 for a root).
 *
 * 2. The LR criterion.  Imagine a planar drawing in which the DFS tree is
 *    drawn.  A back edge (u -> a) closes a cycle with the tree path from a
 *    down to u, and it runs around either the LEFT or the RIGHT of that
 *    path.  Constraints between back edges arise only at "forks": a vertex
 *    v with two outgoing edges e1, e2 whose subtrees both contain back
 *    edges returning to proper ancestors of v ("return edges" of e1 and
 *    e2).  Return edges of the two branches whose return points interlace
 *    along the tree path must be drawn on different sides (they would
 *    cross otherwise); return edges of the same branch whose return points
 *    enclose a return point of the other branch must be on the same side.
 *    de Fraysseix and Rosenstiehl proved that a graph is planar if and
 *    only if the back edges can be split into a left class and a right
 *    class satisfying all fork constraints.  Brandes shows how to check
 *    this in linear time, which is what the code below does.
 *
 * 3. Quantities computed in the first DFS (dfs1):
 *      lowpt(e)   the lowest height (closest to the root) among the
 *                 source of e and the return points of all back edges in
 *                 the subtree below e.  For a back edge it is simply the
 *                 height of its target.  A tree edge e out of v "has a
 *                 return edge" iff lowpt(e) < height(v).
 *      lowpt2(e)  the second-lowest such height.
 *      nesting(e) = 2*lowpt(e) + [lowpt2(e) < height(source(e))].
 *    The outgoing edges of every vertex are then sorted by increasing
 *    nesting depth: edges whose subtrees return lowest come first, and
 *    among equal lowpt the "chordal" ones (with a second return point
 *    above the source) come last.  Processing edges in this order is what
 *    allows the constraints to be handled with a simple stack.
 *
 * 4. The second DFS (dfs2) walks the tree again in that order and keeps a
 *    stack S of CONFLICT PAIRS.  A conflict pair holds two INTERVALS of
 *    back edges, L and R, that must end up on opposite sides of the tree
 *    path.  An interval {low, high} is a run of return edges sorted by
 *    return point: `low` returns lowest (nearest the root), `high`
 *    returns highest (nearest the fork).  When a back edge is met it is
 *    pushed as a pair with an empty L and a one-edge R.  When the walk
 *    returns to a vertex v after finishing an outgoing edge e_i that has
 *    return edges, add_constraints() folds everything pushed while inside
 *    e_i into ONE interval (all return edges of e_i must be on a common
 *    side relative to the earlier edges of v), and moves the return edges
 *    of the earlier edges e_1..e_{i-1} that conflict with e_i, i.e. that
 *    return above lowpt(e_i), into the OTHER interval of a new pair.  A
 *    contradiction -- a pair whose two intervals would both have to be on
 *    the same side -- means the graph is not planar; there are exactly two
 *    places where this is detected, marked NOT PLANAR below.
 *
 * 5. When the walk returns from a vertex to its parent u, back edges that
 *    return exactly to u are finished: they cannot conflict with anything
 *    decided higher up, so trim_back_edges() removes them from the top of
 *    the intervals and drops pairs that have become empty.  The `ref`
 *    links (edges that must lie on the same side as another edge) make it
 *    possible to step from an interval's `high` edge to the next one below
 *    it without storing the interval explicitly.
 *
 * Arrays (all indexed by directed edge id unless noted):
 *   eid[v][w]        id of the undirected edge {v,w}
 *   esrc/etgt        orientation after dfs1 (source, target)
 *   oriented         edge already oriented (seen from one endpoint)
 *   height[v]        DFS depth of v, NONE if not yet visited
 *   parent_edge[v]   tree edge entering v, NONE for roots
 *   lowpt, lowpt2, nesting   see item 3
 *   adj_out[v][]     outgoing edges of v (nout[v] of them), later sorted
 *   lowpt_edge[e]    a return edge of e realising lowpt(e) (used to
 *                    align sides of the lowest return edges)
 *   ref_[e]          edge that e must be on the same side as (NONE if none)
 *   stack_bottom[e]  height of S when the processing of e started; the
 *                    pairs above it were pushed while inside e's subtree
 *   S, sp            the conflict-pair stack and its size
 * The `side` array of the original algorithm is only needed to construct
 * an embedding and is left out.
 */

#include "lrplanar.h"

#define LRN   WORDSIZE
#define LRE   (3*LRN)          /* we only test graphs with <= 3n-6 edges */
#define NONE  (-1)

typedef struct { int low, high; } interval;       /* edge ids, NONE = empty */
typedef struct { interval L, R; } cpair;

static int nv;
static int eid[LRN][LRN];                  /* undirected edge -> id */
static int esrc[LRE], etgt[LRE];           /* orientation (source, target) */
static int oriented[LRE];
static int height[LRN], parent_edge[LRN];
static int lowpt[LRE], lowpt2[LRE], nesting[LRE];
static int adj_out[LRN][LRN], nout[LRN];   /* outgoing edges per vertex */
static int lowpt_edge[LRE], ref_[LRE], stack_bottom[LRE];
static cpair S[LRE + 2];
static int sp;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define EMPTY(I) ((I).low == NONE && (I).high == NONE)

/* Phase 1: orientation DFS from v.  Orients every edge at v that has not
 * been oriented from the other end yet, recurses into unvisited
 * neighbours (tree edges), records back edges, and computes lowpt, lowpt2
 * and the nesting depth of each edge.  Neighbours are visited in
 * increasing vertex number (TAKEBIT extracts the lowest-numbered vertex of
 * a nauty set first), so the DFS is deterministic. */
static void
dfs1(const graph *g, int v)
{
    int e = parent_edge[v];
    setword nb = g[v];
    int w, ei;

    while (nb)
    {
        TAKEBIT(w, nb);
        ei = eid[v][w];
        /* Already oriented from w: it is either the tree edge into v or a
           back edge from a descendant of v.  Nothing to do here. */
        if (oriented[ei]) continue;
        oriented[ei] = 1;
        esrc[ei] = v; etgt[ei] = w;
        /* Before exploring below, the lowest point reachable through ei
           is v itself. */
        lowpt[ei] = lowpt2[ei] = height[v];
        adj_out[v][nout[v]++] = ei;
        if (height[w] == NONE)                 /* tree edge */
        {
            parent_edge[w] = ei;
            height[w] = height[v] + 1;
            dfs1(g, w);
            /* lowpt/lowpt2 of ei were updated by the recursive calls for
               w's outgoing edges (see the "update lowpoints" block). */
        }
        else                                   /* back edge */
            lowpt[ei] = height[w];             /* it returns to w */

        /* Nesting depth: order primarily by lowpt; among edges with the
           same lowpt, those that also return somewhere above the source
           ("chordal") are placed after those that do not. */
        nesting[ei] = 2*lowpt[ei] + (lowpt2[ei] < height[v] ? 1 : 0);

        if (e != NONE)                         /* update lowpoints of e */
        {
            /* Fold the two lowest return heights of ei into those of the
               tree edge e that leads into v, keeping lowpt[e] <= lowpt2[e]. */
            if (lowpt[ei] < lowpt[e])
            {
                lowpt2[e] = MIN(lowpt[e], lowpt2[ei]);
                lowpt[e] = lowpt[ei];
            }
            else if (lowpt[ei] > lowpt[e])
                lowpt2[e] = MIN(lowpt2[e], lowpt[ei]);
            else
                lowpt2[e] = MIN(lowpt2[e], lowpt2[ei]);
        }
    }
}

/* Does interval I conflict with back edge b?  It does when I is not empty
 * and its highest return edge returns ABOVE (further from the root than)
 * the return point of b: those edges would have to cross b if drawn on the
 * same side. */
static int
conflicting(interval I, int b)
{
    return !EMPTY(I) && I.high != NONE && lowpt[I.high] > lowpt[b];
}

/* Lowest return point of a conflict pair (over both intervals). */
static int
lowest(cpair P)
{
    if (EMPTY(P.L)) return lowpt[P.R.low];
    if (EMPTY(P.R)) return lowpt[P.L.low];
    return MIN(lowpt[P.L.low], lowpt[P.R.low]);
}

/* Called when the walk has finished outgoing edge ei of vertex v, where
 * ei is not the first outgoing edge of v (that one, e_1, sets the
 * reference side) and ei has at least one return edge.  e is the tree
 * edge entering v.  Builds a new conflict pair P: P.R collects all return
 * edges of ei, P.L collects the return edges of the earlier edges of v
 * that conflict with ei.  Returns 0 if a contradiction is found. */
static int
add_constraints(int ei, int e)
{
    cpair P, Q;
    interval t;

    P.L.low = P.L.high = P.R.low = P.R.high = NONE;

    /* merge return edges of e_i into P.R */
    /* Everything above stack_bottom[ei] was pushed while inside ei's
       subtree.  All of it has to end up on one side (relative to the
       other edges of v), so each popped pair Q must have one empty
       interval; swap so that the non-empty one is Q.R. */
    do
    {
        Q = S[--sp];
        if (!EMPTY(Q.L)) { t = Q.L; Q.L = Q.R; Q.R = t; }
        /* NOT PLANAR: both intervals of Q are non-empty, i.e. two groups of
           return edges of ei must be on different sides of each other, yet
           all return edges of ei must be on a common side. */
        if (!EMPTY(Q.L)) return 0;                      /* not planar */
        if (lowpt[Q.R.low] > lowpt[e])                  /* merge intervals */
        {
            /* Q.R returns above lowpt(e): it stays relevant higher up, so
               append it below what P.R already holds (intervals are
               chained through ref_ from high to low). */
            if (EMPTY(P.R)) P.R.high = Q.R.high;
            else if (P.R.low != NONE) ref_[P.R.low] = Q.R.high;
            P.R.low = Q.R.low;
        }
        else                                            /* align */
            /* Q.R returns as low as e itself does: no conflict is possible
               any more, only the side must agree with e's lowest return
               edge.  Record that and drop it from the stack. */
            ref_[Q.R.low] = lowpt_edge[e];
    } while (sp > stack_bottom[ei]);

    /* merge conflicting return edges of e_1..e_{i-1} into P.L */
    /* Pairs left on the stack belong to earlier outgoing edges of v.  As
       long as the top pair has an interval returning above lowpt(ei), that
       interval must go on the side opposite to ei's return edges: move it
       into P.L (swapping so the conflicting interval is Q.L, i.e. Q.R is
       the non-conflicting one). */
    while (sp > 0 && (conflicting(S[sp-1].L, ei) || conflicting(S[sp-1].R, ei)))
    {
        Q = S[--sp];
        if (conflicting(Q.R, ei)) { t = Q.L; Q.L = Q.R; Q.R = t; }
        /* NOT PLANAR: both intervals of Q conflict with ei, so both would
           have to be opposite to ei, i.e. on the same side as each other,
           but they form a conflict pair and must differ. */
        if (conflicting(Q.R, ei)) return 0;             /* not planar */
        /* merge interval below lowpt(e_i) into P.R */
        /* Q.R does not conflict with ei, so it joins ei's side. */
        if (P.R.low != NONE) ref_[P.R.low] = Q.R.high;
        if (Q.R.low != NONE) P.R.low = Q.R.low;
        /* Q.L conflicts with ei: it goes to the opposite side, P.L. */
        if (EMPTY(P.L)) P.L.high = Q.L.high;
        else if (P.L.low != NONE) ref_[P.L.low] = Q.L.high;
        P.L.low = Q.L.low;
    }
    /* The guards `P.R.low != NONE` / `P.L.low != NONE` above correspond to
       places where the reference pseudo-code writes ref of a possibly
       empty interval; there they are harmless no-ops and here they must be
       skipped explicitly. */
    if (!(EMPTY(P.L) && EMPTY(P.R))) S[sp++] = P;
    return 1;
}

/* Called when the walk returns from a child to its parent u.  Back edges
 * returning exactly to u are complete: remove them.  Pairs whose lowest
 * return point is u disappear entirely; in the (at most one) remaining
 * pair that still contains such edges, they sit at the `high` end of its
 * intervals and are peeled off by following the ref_ chain. */
static void
trim_back_edges(int u)
{
    cpair P;

    /* drop entire conflict pairs */
    while (sp > 0 && lowest(S[sp-1]) == height[u]) --sp;

    if (sp > 0)                       /* one more conflict pair to consider */
    {
        P = S[--sp];
        /* trim left interval */
        while (P.L.high != NONE && etgt[P.L.high] == u) P.L.high = ref_[P.L.high];
        if (P.L.high == NONE && P.L.low != NONE)        /* just emptied */
        {
            /* The whole interval returned to u.  Its lowest edge is now
               only tied to the other side's lowest edge (needed for the
               embedding, harmless for the decision). */
            ref_[P.L.low] = P.R.low;
            P.L.low = NONE;
        }
        /* trim right interval */
        while (P.R.high != NONE && etgt[P.R.high] == u) P.R.high = ref_[P.R.high];
        if (P.R.high == NONE && P.R.low != NONE)        /* just emptied */
        {
            ref_[P.R.low] = P.L.low;
            P.R.low = NONE;
        }
        S[sp++] = P;
    }
}

/* Phase 2: testing DFS from v, visiting the outgoing edges in nesting
 * order.  Returns 0 as soon as a contradiction is found. */
static int
dfs2(int v)
{
    int e = parent_edge[v];
    int i, ei, w, u, hl, hr;

    for (i = 0; i < nout[v]; ++i)
    {
        ei = adj_out[v][i];
        w = etgt[ei];
        /* Remember how deep the stack is before ei: everything pushed from
           now until ei is finished belongs to ei's subtree. */
        stack_bottom[ei] = sp;
        if (ei == parent_edge[w])                       /* tree edge */
        {
            if (!dfs2(w)) return 0;
        }
        else                                            /* back edge */
        {
            /* A back edge is its own (and only) return edge; push it as a
               pair with an empty left interval. */
            lowpt_edge[ei] = ei;
            S[sp].L.low = S[sp].L.high = NONE;
            S[sp].R.low = S[sp].R.high = ei;
            ++sp;
        }
        if (lowpt[ei] < height[v])                      /* e_i has return edge */
        {
            /* The first outgoing edge (lowest nesting depth) just defines
               the lowest return edge of e; every later edge is constrained
               against the earlier ones. */
            if (i == 0)
                lowpt_edge[e] = lowpt_edge[ei];
            else if (!add_constraints(ei, e))
                return 0;
        }
    }

    if (e != NONE)                                      /* v is not a root */
    {
        u = esrc[e];
        /* Back edges returning to the parent u are finished now. */
        trim_back_edges(u);
        if (lowpt[e] < height[u] && sp > 0)             /* e has return edge */
        {
            /* Side of e is the side of its highest remaining return edge;
               recorded through ref_ for the embedding phase, and used by
               later trims to walk the interval chain. */
            hl = S[sp-1].L.high;
            hr = S[sp-1].R.high;
            if (hl != NONE && (hr == NONE || lowpt[hl] > lowpt[hr])) ref_[e] = hl;
            else ref_[e] = hr;
        }
    }
    return 1;
}

/* Entry point: 1 if the graph g on n vertices is planar, else 0. */
int
lr_isplanar(const graph *g, int n)
{
    int v, w, m, i, j, k, x;
    setword nb;

    /* Exact shortcuts.  K5 has 10 edges and K3,3 has 9, so a graph with
       fewer than 9 edges contains no Kuratowski subdivision and is planar;
       a planar graph has at most 3n-6 edges (Euler's formula). */
    if (n <= 4) return 1;
    m = 0;
    for (v = 0; v < n; ++v) m += POPCOUNT(g[v]);
    m /= 2;
    if (m > 3*n - 6) return 0;
    if (m < 9) return 1;                 /* K5 has 10 edges, K3,3 has 9 */

    /* Number the undirected edges and reset per-vertex state. */
    nv = n;
    k = 0;
    for (v = 0; v < n; ++v)
    {
        height[v] = NONE; parent_edge[v] = NONE; nout[v] = 0;
        nb = g[v];
        while (nb)
        {
            TAKEBIT(w, nb);
            if (w > v) { eid[v][w] = eid[w][v] = k; ++k; }
        }
    }
    for (i = 0; i < k; ++i) { oriented[i] = 0; ref_[i] = NONE; lowpt_edge[i] = NONE; }

    /* orientation */
    /* One DFS per connected component; each component gets its own root at
       height 0. */
    for (v = 0; v < n; ++v)
        if (height[v] == NONE) { height[v] = 0; dfs1(g, v); }

    /* sort outgoing edges by nesting depth (insertion sort; lists are short) */
    for (v = 0; v < n; ++v)
        for (i = 1; i < nout[v]; ++i)
        {
            x = adj_out[v][i];
            for (j = i - 1; j >= 0 && nesting[adj_out[v][j]] > nesting[x]; --j)
                adj_out[v][j+1] = adj_out[v][j];
            adj_out[v][j+1] = x;
        }

    /* testing */
    /* The stack is shared by all components; it is empty between them. */
    sp = 0;
    for (v = 0; v < n; ++v)
        if (parent_edge[v] == NONE && !dfs2(v)) return 0;
    return 1;
}
