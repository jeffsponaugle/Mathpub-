/* lrplanar_sg.c -- Left-Right planarity test for graphs of any size.
 *
 * This is the sparsegraph version of lrplanar.c (which is restricted to
 * n <= WORDSIZE vertices in nauty's bitmask format and uses recursion).
 * The algorithm is unchanged: the left-right planarity test of de Fraysseix
 * and Rosenstiehl in the formulation of U. Brandes, "The Left-Right
 * Planarity Test" (2009); see the long comment in lrplanar.c for how it
 * works.  The differences are:
 *
 *   - Input is a nauty sparsegraph (nausparse.h): v[] indexes into e[],
 *     d[] holds degrees, e[] the neighbour lists.  For an undirected graph
 *     every edge must be present from both endpoints.  Loops and parallel
 *     edges are ignored, so multigraphs are accepted and answered for the
 *     underlying simple graph.  With digraph = TRUE the lists are taken to
 *     hold arcs and the question is asked of the underlying undirected
 *     graph (u-v is an edge iff u->v or v->u is an arc): the arc lists are
 *     first symmetrized into work arrays in linear time, so undirected input
 *     pays nothing for the feature.
 *   - No size limit.  Following nauty's conventions, vertex numbers,
 *     heights and degrees are int (n < 2^31) while every edge count and
 *     edge number is size_t, so graphs with more than 2^31 edges are
 *     handled on 64-bit systems.  Work space is allocated with nauty's
 *     DYNALLSTAT / DYNALLOC1 (thread-local when nauty is built with TLS),
 *     grows as needed and is kept between calls; lrplanar_freedyn()
 *     releases it.
 *   - Both depth-first searches are iterative, with explicit stacks, so a
 *     path on ten million vertices is as safe as a triangle.
 *   - Outgoing edges are sorted by nesting depth with one global counting
 *     sort (linear time) instead of an insertion sort per vertex, so a
 *     vertex of degree 10^6 does not cost 10^12 steps.
 *   - Edges are identified without a v x v table.  When a vertex v is
 *     entered, its neighbour list is copied without loops and duplicates
 *     into v's slice of the work array out[] (one pass, using a mark per
 *     vertex).  When the DFS at v then meets a neighbour w that is already
 *     visited, w is either v's parent (the tree edge into v: skip), a
 *     descendant (height[w] > height[v]: the edge was already oriented from
 *     w as a back edge: skip) or a proper ancestor (orient a new back edge
 *     v -> w).  Undirected DFS produces no other kind of edge.  The edge
 *     numbers oriented from v overwrite the copied list from the front; they
 *     never overtake the position being read.
 *
 * Memory: per vertex 5 ints + 3 size_t-or-int stack entries (about 44
 * bytes); per undirected edge 2 ints (etgt, lowpt), 3 size_t (esrc/
 * lowpt_edge, lowpt2/ref, nesting/stack_bottom: arrays needed only in the
 * first phase are reused in the second), the sort buffer (size_t), out[]
 * with 2 size_t per edge, and the conflict-pair stack (4 size_t per back
 * edge): about 100 bytes per edge, roughly 350 bytes per vertex of a
 * planar graph with 3n edges.
 *
 * Early exits (all exact): n <= 4 is planar; more than 3n-6 distinct edges
 * (found while orienting, so a dense input is rejected after 3n-5 edges)
 * is not; fewer than 9 distinct edges is planar (K3,3 has 9).
 *
 * Written 2026 for the enumeration of planar graphs in OEIS A003094 (see
 * README.md); validated against nauty's planarity.c (test suite: run_tests.sh).
 */

#include "lrplanar_sg.h"

#define NONE   (-1)                 /* no vertex / not visited (int) */
#define NOEDGE ((size_t)-1)         /* no edge (size_t) */
#define MIN(a,b) ((a) < (b) ? (a) : (b))

typedef struct { size_t low, high; } interval;   /* edge numbers, NOEDGE = empty */
typedef struct { interval L, R; } cpair;          /* conflict pair */
#define EMPTY(I) ((I).low == NOEDGE && (I).high == NOEDGE)

/* --- work space, grown on demand, kept between calls ------------------- */

/* per vertex */
DYNALLSTAT(int,height,height_sz);           /* DFS depth, NONE = unvisited */
DYNALLSTAT(size_t,parent_edge,parent_edge_sz); /* tree edge into it, NOEDGE for roots */
DYNALLSTAT(int,parent_v,parent_v_sz);       /* DFS parent vertex, NONE for roots */
DYNALLSTAT(int,nout,nout_sz);               /* number of outgoing edges (<= degree) */
DYNALLSTAT(int,mark,mark_sz);               /* vertex whose list is being deduplicated */
DYNALLSTAT(int,frame_v,frame_v_sz);         /* DFS stack: vertex */
DYNALLSTAT(int,frame_i,frame_i_sz);         /* DFS stack: position in its list */
DYNALLSTAT(int,frame_n,frame_n_sz);         /* DFS stack: length of its list */
DYNALLSTAT(size_t,count,count_sz);          /* counting sort buckets (2n+2) */

/* per edge */
DYNALLSTAT(int,etgt,etgt_sz);               /* target (head) after orientation */
DYNALLSTAT(size_t,esrc,esrc_sz);            /* source vertex; phase 2: lowpt_edge */
DYNALLSTAT(int,lowpt,lowpt_sz);             /* lowest return height */
DYNALLSTAT(size_t,lowpt2,lowpt2_sz);        /* second lowest height; phase 2: ref */
DYNALLSTAT(size_t,nesting,nesting_sz);      /* nesting depth; phase 2: stack_bottom */
DYNALLSTAT(size_t,out,out_sz);              /* per vertex slice (like sg->e): neighbours, then edge numbers */
DYNALLSTAT(size_t,sorted,sorted_sz);        /* sort buffer (one entry per edge) */
DYNALLSTAT(cpair,S,S_sz);                   /* conflict-pair stack */

/* symmetrized copy of a digraph (only used when digraph = TRUE) */
DYNALLSTAT(size_t,sym_v,sym_v_sz);
DYNALLSTAT(int,sym_d,sym_d_sz);
DYNALLSTAT(int,sym_e,sym_e_sz);
DYNALLSTAT(int,sym_fill,sym_fill_sz);
static TLS_ATTR sparsegraph symsg;

static TLS_ATTR size_t sp;                  /* conflict-pair stack size */
static TLS_ATTR size_t *lowpt_edge, *ref_, *stack_bottom;   /* phase-2 aliases */

void
lrplanar_freedyn(void)
{
    DYNFREE(height,height_sz); DYNFREE(parent_edge,parent_edge_sz);
    DYNFREE(parent_v,parent_v_sz); DYNFREE(nout,nout_sz); DYNFREE(mark,mark_sz);
    DYNFREE(frame_v,frame_v_sz); DYNFREE(frame_i,frame_i_sz);
    DYNFREE(frame_n,frame_n_sz); DYNFREE(count,count_sz);
    DYNFREE(etgt,etgt_sz); DYNFREE(esrc,esrc_sz); DYNFREE(lowpt,lowpt_sz);
    DYNFREE(lowpt2,lowpt2_sz); DYNFREE(nesting,nesting_sz); DYNFREE(out,out_sz);
    DYNFREE(sorted,sorted_sz); DYNFREE(S,S_sz);
    DYNFREE(sym_v,sym_v_sz); DYNFREE(sym_d,sym_d_sz);
    DYNFREE(sym_e,sym_e_sz); DYNFREE(sym_fill,sym_fill_sz);
}

/* --- digraph input: symmetrize the arc lists ---------------------------- */

/* Return a sparsegraph (static work space) in which every arc u->w of sg
 * appears as w in u's list and u in w's list.  Duplicates and loops are
 * left in; the DFS removes them as for any multigraph. */
static sparsegraph *
symmetrize(sparsegraph *sg)
{
    int n = sg->nv, u, w;
    size_t j, nde2 = 2*sg->nde, pos;

    DYNALLOC1(size_t,sym_v,sym_v_sz,n > 0 ? n : 1,"lrplanar_sg");
    DYNALLOC1(int,sym_d,sym_d_sz,n > 0 ? n : 1,"lrplanar_sg");
    DYNALLOC1(int,sym_fill,sym_fill_sz,n > 0 ? n : 1,"lrplanar_sg");
    DYNALLOC1(int,sym_e,sym_e_sz,nde2 > 0 ? nde2 : 1,"lrplanar_sg");

    for (u = 0; u < n; ++u) sym_d[u] = 0;
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            ++sym_d[u];
            ++sym_d[sg->e[j]];
        }
    pos = 0;
    for (u = 0; u < n; ++u) { sym_v[u] = pos; pos += (size_t)sym_d[u]; sym_fill[u] = 0; }
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            w = sg->e[j];
            sym_e[sym_v[u] + (size_t)sym_fill[u]++] = w;
            sym_e[sym_v[w] + (size_t)sym_fill[w]++] = u;
        }

    symsg.nv = n; symsg.nde = nde2;
    symsg.v = sym_v; symsg.d = sym_d; symsg.e = sym_e; symsg.w = NULL;
    symsg.vlen = sym_v_sz; symsg.dlen = sym_d_sz; symsg.elen = sym_e_sz; symsg.wlen = 0;
    return &symsg;
}

/* --- phase 1: orientation ---------------------------------------------- */

/* Copy v's neighbours without loops and duplicates into out[v0..]; return
 * how many.  Called once per vertex, when it is entered. */
static int
dedup_list(sparsegraph *sg, int v)
{
    int *e = sg->e + sg->v[v];
    size_t *o = out + sg->v[v];
    int j, w, cnt = 0, d = sg->d[v];

    for (j = 0; j < d; ++j)
    {
        w = e[j];
        if (w == v || mark[w] == v) continue;
        mark[w] = v;
        o[cnt++] = (size_t)w;
    }
    return cnt;
}

/* Bookkeeping after outgoing edge ei of v has been fully explored (for a
 * tree edge: after its subtree; for a back edge: immediately): nesting
 * depth of ei, and lowpt/lowpt2 of the tree edge e into v.  lowpt2 holds
 * heights (small ints) although its array is size_t for later reuse. */
static void
poststep1(int v, size_t ei)
{
    size_t e = parent_edge[v];

    nesting[ei] = 2*(size_t)lowpt[ei] + ((int)lowpt2[ei] < height[v] ? 1 : 0);

    if (e != NOEDGE)
    {
        if (lowpt[ei] < lowpt[e])
        {
            lowpt2[e] = (size_t)MIN(lowpt[e], (int)lowpt2[ei]);
            lowpt[e] = lowpt[ei];
        }
        else if (lowpt[ei] > lowpt[e])
            lowpt2[e] = (size_t)MIN((int)lowpt2[e], lowpt[ei]);
        else
            lowpt2[e] = (size_t)MIN((int)lowpt2[e], (int)lowpt2[ei]);
    }
}

/* Iterative DFS from root.  Assigns edge numbers *pk, *pk+1, ... to the
 * distinct non-loop edges met, orients them, fills the outgoing lists,
 * computes height, lowpt, lowpt2 and nesting.  Returns FALSE as soon as
 * more than maxk edges have been found (the graph is then not planar). */
static boolean
dfs1(sparsegraph *sg, int root, size_t *pk, size_t maxk)
{
    size_t *vv = sg->v;
    size_t k = *pk, ei;
    int fp, v, w, j, u;

    height[root] = 0;
    fp = 0;
    frame_v[0] = root;
    frame_i[0] = 0;
    frame_n[0] = dedup_list(sg, root);

    while (fp >= 0)
    {
        v = frame_v[fp];
        j = frame_i[fp];
        if (j < frame_n[fp])
        {
            frame_i[fp] = j + 1;
            w = (int)out[vv[v] + (size_t)j];
            if (w == parent_v[v]) continue;          /* the tree edge into v */
            if (height[w] != NONE && height[w] > height[v]) continue;   /* descendant: oriented from w */

            if (k > maxk) { *pk = k; return FALSE; } /* too many edges for a planar graph */
            ei = k++;
            esrc[ei] = (size_t)v;
            etgt[ei] = w;
            lowpt[ei] = height[v];
            lowpt2[ei] = (size_t)height[v];
            out[vv[v] + (size_t)nout[v]++] = ei;     /* nout[v] <= j: slot already read */

            if (height[w] == NONE)                   /* tree edge */
            {
                parent_edge[w] = ei;
                parent_v[w] = v;
                height[w] = height[v] + 1;
                ++fp;
                frame_v[fp] = w;
                frame_i[fp] = 0;
                frame_n[fp] = dedup_list(sg, w);
            }
            else                                     /* back edge to a proper ancestor */
            {
                lowpt[ei] = height[w];
                poststep1(v, ei);
            }
        }
        else                                         /* v is finished */
        {
            --fp;
            if (fp >= 0)                             /* its tree edge, seen from the parent u */
            {
                u = parent_v[v];
                poststep1(u, parent_edge[v]);
            }
        }
    }
    *pk = k;
    return TRUE;
}

/* --- phase 2: testing --------------------------------------------------- */

static boolean
conflicting(interval I, size_t b)
{
    return !EMPTY(I) && I.high != NOEDGE && lowpt[I.high] > lowpt[b];
}

static int
lowest(cpair P)
{
    if (EMPTY(P.L)) return lowpt[P.R.low];
    if (EMPTY(P.R)) return lowpt[P.L.low];
    return MIN(lowpt[P.L.low], lowpt[P.R.low]);
}

/* Outgoing edge ei of v (not the first one) has a return edge; e is the
 * tree edge into v.  Fold the constraints of ei into a new conflict pair.
 * Returns FALSE if the graph is not planar. */
static boolean
add_constraints(size_t ei, size_t e)
{
    cpair P, Q;
    interval t;

    P.L.low = P.L.high = P.R.low = P.R.high = NOEDGE;

    /* merge the return edges of ei (everything above stack_bottom[ei]) into P.R */
    do
    {
        Q = S[--sp];
        if (!EMPTY(Q.L)) { t = Q.L; Q.L = Q.R; Q.R = t; }
        if (!EMPTY(Q.L)) return FALSE;                  /* NOT PLANAR */
        if (lowpt[Q.R.low] > lowpt[e])                  /* merge intervals */
        {
            if (EMPTY(P.R)) P.R.high = Q.R.high;
            else if (P.R.low != NOEDGE) ref_[P.R.low] = Q.R.high;
            P.R.low = Q.R.low;
        }
        else                                            /* align */
            ref_[Q.R.low] = lowpt_edge[e];
    } while (sp > stack_bottom[ei]);

    /* merge the return edges of the earlier edges of v that conflict with ei into P.L */
    while (sp > 0 && (conflicting(S[sp-1].L, ei) || conflicting(S[sp-1].R, ei)))
    {
        Q = S[--sp];
        if (conflicting(Q.R, ei)) { t = Q.L; Q.L = Q.R; Q.R = t; }
        if (conflicting(Q.R, ei)) return FALSE;         /* NOT PLANAR */
        if (P.R.low != NOEDGE) ref_[P.R.low] = Q.R.high;
        if (Q.R.low != NOEDGE) P.R.low = Q.R.low;
        if (EMPTY(P.L)) P.L.high = Q.L.high;
        else if (P.L.low != NOEDGE) ref_[P.L.low] = Q.L.high;
        P.L.low = Q.L.low;
    }
    if (!(EMPTY(P.L) && EMPTY(P.R))) S[sp++] = P;
    return TRUE;
}

/* Back edges returning to u are finished when the walk returns to u. */
static void
trim_back_edges(int u)
{
    cpair P;

    while (sp > 0 && lowest(S[sp-1]) == height[u]) --sp;

    if (sp > 0)
    {
        P = S[--sp];
        while (P.L.high != NOEDGE && etgt[P.L.high] == u) P.L.high = ref_[P.L.high];
        if (P.L.high == NOEDGE && P.L.low != NOEDGE)
        {
            ref_[P.L.low] = P.R.low;
            P.L.low = NOEDGE;
        }
        while (P.R.high != NOEDGE && etgt[P.R.high] == u) P.R.high = ref_[P.R.high];
        if (P.R.high == NOEDGE && P.R.low != NOEDGE)
        {
            ref_[P.R.low] = P.L.low;
            P.R.low = NOEDGE;
        }
        S[sp++] = P;
    }
}

/* After outgoing edge number i (edge ei) of v has been processed. */
static boolean
poststep2(int v, int i, size_t ei)
{
    size_t e = parent_edge[v];

    if (e != NOEDGE && lowpt[ei] < height[v])           /* ei has a return edge (never at a root) */
    {
        if (i == 0)
            lowpt_edge[e] = lowpt_edge[ei];
        else if (!add_constraints(ei, e))
            return FALSE;
    }
    return TRUE;
}

/* Iterative testing DFS from root over the sorted outgoing lists.
 * Returns FALSE if the graph is not planar. */
static boolean
dfs2(sparsegraph *sg, int root)
{
    size_t *vv = sg->v;
    size_t ei, e, hl, hr;
    int fp, v, i, w, u;

    fp = 0;
    frame_v[0] = root;
    frame_i[0] = 0;

    while (fp >= 0)
    {
        v = frame_v[fp];
        i = frame_i[fp];
        if (i < nout[v])
        {
            ei = out[vv[v] + (size_t)i];
            w = etgt[ei];
            stack_bottom[ei] = sp;
            if (ei == parent_edge[w])                   /* tree edge: descend, poststep2 on return */
            {
                ++fp;
                frame_v[fp] = w;
                frame_i[fp] = 0;
            }
            else                                        /* back edge */
            {
                lowpt_edge[ei] = ei;
                S[sp].L.low = S[sp].L.high = NOEDGE;
                S[sp].R.low = S[sp].R.high = ei;
                ++sp;
                if (!poststep2(v, i, ei)) return FALSE;
                frame_i[fp] = i + 1;
            }
        }
        else                                            /* v is finished */
        {
            e = parent_edge[v];
            if (e != NOEDGE)
            {
                u = parent_v[v];
                trim_back_edges(u);
                if (lowpt[e] < height[u] && sp > 0)     /* e has a return edge: record its side */
                {
                    hl = S[sp-1].L.high;
                    hr = S[sp-1].R.high;
                    if (hl != NOEDGE && (hr == NOEDGE || lowpt[hl] > lowpt[hr])) ref_[e] = hl;
                    else ref_[e] = hr;
                }
            }
            --fp;
            if (fp >= 0)                                /* finish tree edge e in the parent's frame */
            {
                u = frame_v[fp];
                i = frame_i[fp];
                if (!poststep2(u, i, e)) return FALSE;
                frame_i[fp] = i + 1;
            }
        }
    }
    return TRUE;
}

/* --- driver -------------------------------------------------------------- */

boolean
lrplanar_sg(sparsegraph *sg, boolean digraph)
{
    int n = sg->nv;
    size_t nde, extent, x, k, maxk, nb, i, ei, maxkey;
    int v;

    if (n <= 4) return TRUE;                         /* every graph on <= 4 vertices is planar */
    if (digraph) sg = symmetrize(sg);                /* work on the underlying undirected graph */
    nde = sg->nde;

    maxk = 3*(size_t)n - 6;                /* a planar graph has at most this many edges */
    nb = MIN(nde, maxk + 1);               /* at most maxk+1 edges are ever numbered */
    if (nb < 1) nb = 1;

    /* out[] mirrors the layout of sg->e, whose lists may leave gaps */
    extent = 1;
    for (v = 0; v < n; ++v)
    {
        x = sg->v[v] + (size_t)sg->d[v];
        if (x > extent) extent = x;
    }

    DYNALLOC1(int,height,height_sz,n,"lrplanar_sg");
    DYNALLOC1(size_t,parent_edge,parent_edge_sz,n,"lrplanar_sg");
    DYNALLOC1(int,parent_v,parent_v_sz,n,"lrplanar_sg");
    DYNALLOC1(int,nout,nout_sz,n,"lrplanar_sg");
    DYNALLOC1(int,mark,mark_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_v,frame_v_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_i,frame_i_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_n,frame_n_sz,n,"lrplanar_sg");
    DYNALLOC1(int,etgt,etgt_sz,nb,"lrplanar_sg");
    DYNALLOC1(size_t,esrc,esrc_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,lowpt,lowpt_sz,nb,"lrplanar_sg");
    DYNALLOC1(size_t,lowpt2,lowpt2_sz,nb,"lrplanar_sg");
    DYNALLOC1(size_t,nesting,nesting_sz,nb,"lrplanar_sg");
    DYNALLOC1(size_t,out,out_sz,extent,"lrplanar_sg");

    for (v = 0; v < n; ++v)
    {
        height[v] = NONE; parent_edge[v] = NOEDGE; parent_v[v] = NONE;
        nout[v] = 0; mark[v] = NONE;
    }

    /* phase 1: orientation, one DFS per component */
    k = 0;
    for (v = 0; v < n; ++v)
        if (height[v] == NONE && !dfs1(sg, v, &k, maxk))
            return FALSE;                            /* more than 3n-6 distinct edges */
    if (k < 9) return TRUE;                          /* K5 has 10 edges, K3,3 has 9 */

    /* sort every vertex's outgoing edges by nesting depth: one stable
       counting sort over all edges, then redistribute into out[] */
    maxkey = 2*(size_t)n + 1;
    DYNALLOC1(size_t,count,count_sz,maxkey + 1,"lrplanar_sg");
    DYNALLOC1(size_t,sorted,sorted_sz,k,"lrplanar_sg");
    for (i = 0; i <= maxkey; ++i) count[i] = 0;
    for (i = 0; i < k; ++i) ++count[nesting[i]];
    for (i = 1; i <= maxkey; ++i) count[i] += count[i-1];
    for (i = k; i-- > 0; ) sorted[--count[nesting[i]]] = i;
    for (v = 0; v < n; ++v) nout[v] = 0;
    for (i = 0; i < k; ++i)
    {
        ei = sorted[i];
        v = (int)esrc[ei];
        out[sg->v[v] + (size_t)nout[v]++] = ei;
    }

    /* phase 2: reuse phase-1 arrays that are no longer needed */
    lowpt_edge = esrc;          /* esrc was only needed for the redistribution */
    ref_ = lowpt2;              /* lowpt2 was only needed for nesting depths */
    stack_bottom = nesting;     /* nesting was only needed for the sort */
    for (i = 0; i < k; ++i) { ref_[i] = NOEDGE; lowpt_edge[i] = NOEDGE; }
    DYNALLOC1(cpair,S,S_sz,k + 2,"lrplanar_sg");

    sp = 0;
    for (v = 0; v < n; ++v)
        if (parent_edge[v] == NOEDGE && !dfs2(sg, v)) return FALSE;
    return TRUE;
}

/* Dense-format wrapper: converts to a sparsegraph kept between calls. */
boolean
lrplanar_dense(graph *g, int m, int n, boolean digraph)
{
    static TLS_ATTR sparsegraph sg;
    static TLS_ATTR boolean init = FALSE;

    if (!init) { SG_INIT(sg); init = TRUE; }
    nauty_to_sg(g, &sg, m, n);
    return lrplanar_sg(&sg, digraph);
}
