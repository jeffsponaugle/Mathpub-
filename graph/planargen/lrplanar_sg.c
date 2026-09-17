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
 *     d[] holds degrees, e[] the neighbour lists; every edge must be present
 *     from both endpoints.  Loops and parallel edges are ignored, so
 *     multigraphs are accepted and answered for the underlying simple graph.
 *   - No size limit.  Work space is allocated with nauty's DYNALLSTAT /
 *     DYNALLOC1 (thread-local when nauty is built with TLS), grows as needed
 *     and is kept between calls; lrplanar_freedyn() releases it.  Vertex and
 *     edge numbers must fit in an int (up to 2^31-1 undirected edges).
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
 *     v -> w).  Undirected DFS produces no other kind of edge.  The edge ids
 *     oriented from v overwrite the copied list from the front; they never
 *     overtake the position being read.
 *
 * Memory: about 10 ints per vertex plus 11 ints per undirected edge (5 edge
 * arrays, out[] with 2 entries per edge, the sort buffer and the
 * conflict-pair stack), i.e. roughly 170 bytes per vertex of a planar
 * graph with 3n edges.  Arrays only needed in the first phase (esrc,
 * lowpt2, nesting) are reused in the second (lowpt_edge, ref, stack_bottom).
 *
 * Early exits (all exact): n <= 4 is planar; more than 3n-6 distinct edges
 * (found while orienting, so a dense input is rejected after 3n-5 edges)
 * is not; fewer than 9 distinct edges is planar (K3,3 has 9).
 *
 * Written 2026 for the enumeration of planar graphs in OEIS A003094 (see
 * README.md); validated against nauty's planarity.c (test suite: run_tests.sh).
 */

#include "lrplanar_sg.h"

#define NONE (-1)
#define MIN(a,b) ((a) < (b) ? (a) : (b))

typedef struct { int low, high; } interval;      /* edge ids, NONE = empty */
typedef struct { interval L, R; } cpair;          /* conflict pair */
#define EMPTY(I) ((I).low == NONE && (I).high == NONE)

/* --- work space, grown on demand, kept between calls ------------------- */

DYNALLSTAT(int,height,height_sz);        /* per vertex: DFS depth, NONE = unvisited */
DYNALLSTAT(int,parent_edge,parent_edge_sz); /* per vertex: tree edge into it, NONE for roots */
DYNALLSTAT(int,parent_v,parent_v_sz);    /* per vertex: DFS parent vertex, NONE for roots */
DYNALLSTAT(int,nout,nout_sz);            /* per vertex: number of outgoing edges */
DYNALLSTAT(int,mark,mark_sz);            /* per vertex: vertex whose list is being deduplicated */
DYNALLSTAT(int,frame_v,frame_v_sz);      /* DFS stack: vertex */
DYNALLSTAT(int,frame_i,frame_i_sz);      /* DFS stack: position in its list */
DYNALLSTAT(int,frame_n,frame_n_sz);      /* DFS stack: length of its list */
DYNALLSTAT(int,count,count_sz);          /* counting sort buckets (2n+2) */

DYNALLSTAT(int,etgt,etgt_sz);            /* per edge: target (head) after orientation */
DYNALLSTAT(int,esrc,esrc_sz);            /* per edge: source; phase 2: lowpt_edge */
DYNALLSTAT(int,lowpt,lowpt_sz);          /* per edge: lowest return height */
DYNALLSTAT(int,lowpt2,lowpt2_sz);        /* per edge: second lowest; phase 2: ref */
DYNALLSTAT(int,nesting,nesting_sz);      /* per edge: nesting depth; phase 2: stack_bottom */
DYNALLSTAT(int,out,out_sz);              /* per vertex slice (like sg->e): neighbours, then outgoing edge ids */
DYNALLSTAT(int,sorted,sorted_sz);        /* sort buffer (one entry per edge) */
DYNALLSTAT(cpair,S,S_sz);                /* conflict-pair stack */

static TLS_ATTR int sp;                  /* conflict-pair stack size */
static TLS_ATTR int *lowpt_edge, *ref_, *stack_bottom;   /* phase-2 aliases */

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
}

/* --- phase 1: orientation ---------------------------------------------- */

/* Copy v's neighbours without loops and duplicates into out[v0..]; return
 * how many.  Called once per vertex, when it is entered. */
static int
dedup_list(sparsegraph *sg, int v)
{
    int *e = sg->e + sg->v[v], *o = out + sg->v[v];
    int j, w, cnt = 0, d = sg->d[v];

    for (j = 0; j < d; ++j)
    {
        w = e[j];
        if (w == v || mark[w] == v) continue;
        mark[w] = v;
        o[cnt++] = w;
    }
    return cnt;
}

/* Bookkeeping after outgoing edge ei of v has been fully explored (for a
 * tree edge: after its subtree; for a back edge: immediately): nesting
 * depth of ei, and lowpt/lowpt2 of the tree edge e into v. */
static void
poststep1(int v, int ei)
{
    int e = parent_edge[v];

    nesting[ei] = 2*lowpt[ei] + (lowpt2[ei] < height[v] ? 1 : 0);

    if (e != NONE)
    {
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

/* Iterative DFS from root.  Assigns edge ids k, k+1, ... to the distinct
 * non-loop edges met, orients them, fills the outgoing lists, computes
 * height, lowpt, lowpt2 and nesting.  Returns the new k, or -1 as soon as
 * more than maxk edges have been found (the graph is then not planar). */
static int
dfs1(sparsegraph *sg, int root, int k, int maxk)
{
    size_t *vv = sg->v;
    int fp, v, w, j, ei, u;

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
            w = out[vv[v] + j];
            if (w == parent_v[v]) continue;          /* the tree edge into v */
            if (height[w] != NONE && height[w] > height[v]) continue;   /* descendant: oriented from w */

            if (k > maxk) return -1;                 /* too many edges for a planar graph */
            ei = k++;
            esrc[ei] = v;
            etgt[ei] = w;
            lowpt[ei] = lowpt2[ei] = height[v];
            out[vv[v] + nout[v]++] = ei;             /* nout[v] <= j: slot already read */

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
    return k;
}

/* --- phase 2: testing --------------------------------------------------- */

static int
conflicting(interval I, int b)
{
    return !EMPTY(I) && I.high != NONE && lowpt[I.high] > lowpt[b];
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
 * Returns 0 if the graph is not planar. */
static int
add_constraints(int ei, int e)
{
    cpair P, Q;
    interval t;

    P.L.low = P.L.high = P.R.low = P.R.high = NONE;

    /* merge the return edges of ei (everything above stack_bottom[ei]) into P.R */
    do
    {
        Q = S[--sp];
        if (!EMPTY(Q.L)) { t = Q.L; Q.L = Q.R; Q.R = t; }
        if (!EMPTY(Q.L)) return 0;                      /* NOT PLANAR */
        if (lowpt[Q.R.low] > lowpt[e])                  /* merge intervals */
        {
            if (EMPTY(P.R)) P.R.high = Q.R.high;
            else if (P.R.low != NONE) ref_[P.R.low] = Q.R.high;
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
        if (conflicting(Q.R, ei)) return 0;             /* NOT PLANAR */
        if (P.R.low != NONE) ref_[P.R.low] = Q.R.high;
        if (Q.R.low != NONE) P.R.low = Q.R.low;
        if (EMPTY(P.L)) P.L.high = Q.L.high;
        else if (P.L.low != NONE) ref_[P.L.low] = Q.L.high;
        P.L.low = Q.L.low;
    }
    if (!(EMPTY(P.L) && EMPTY(P.R))) S[sp++] = P;
    return 1;
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
        while (P.L.high != NONE && etgt[P.L.high] == u) P.L.high = ref_[P.L.high];
        if (P.L.high == NONE && P.L.low != NONE)
        {
            ref_[P.L.low] = P.R.low;
            P.L.low = NONE;
        }
        while (P.R.high != NONE && etgt[P.R.high] == u) P.R.high = ref_[P.R.high];
        if (P.R.high == NONE && P.R.low != NONE)
        {
            ref_[P.R.low] = P.L.low;
            P.R.low = NONE;
        }
        S[sp++] = P;
    }
}

/* After outgoing edge number i (edge ei) of v has been processed. */
static int
poststep2(int v, int i, int ei)
{
    int e = parent_edge[v];

    if (e != NONE && lowpt[ei] < height[v])            /* ei has a return edge (never at a root) */
    {
        if (i == 0)
            lowpt_edge[e] = lowpt_edge[ei];
        else if (!add_constraints(ei, e))
            return 0;
    }
    return 1;
}

/* Iterative testing DFS from root over the sorted outgoing lists.
 * Returns 0 if the graph is not planar. */
static int
dfs2(sparsegraph *sg, int root)
{
    size_t *vv = sg->v;
    int fp, v, i, ei, w, e, u, hl, hr;

    fp = 0;
    frame_v[0] = root;
    frame_i[0] = 0;

    while (fp >= 0)
    {
        v = frame_v[fp];
        i = frame_i[fp];
        if (i < nout[v])
        {
            ei = out[vv[v] + i];
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
                S[sp].L.low = S[sp].L.high = NONE;
                S[sp].R.low = S[sp].R.high = ei;
                ++sp;
                if (!poststep2(v, i, ei)) return 0;
                frame_i[fp] = i + 1;
            }
        }
        else                                            /* v is finished */
        {
            e = parent_edge[v];
            if (e != NONE)
            {
                u = parent_v[v];
                trim_back_edges(u);
                if (lowpt[e] < height[u] && sp > 0)     /* e has a return edge: record its side */
                {
                    hl = S[sp-1].L.high;
                    hr = S[sp-1].R.high;
                    if (hl != NONE && (hr == NONE || lowpt[hl] > lowpt[hr])) ref_[e] = hl;
                    else ref_[e] = hr;
                }
            }
            --fp;
            if (fp >= 0)                                /* finish tree edge e in the parent's frame */
            {
                u = frame_v[fp];
                i = frame_i[fp];
                if (!poststep2(u, i, e)) return 0;
                frame_i[fp] = i + 1;
            }
        }
    }
    return 1;
}

/* --- driver -------------------------------------------------------------- */

boolean
lrplanar_sg(sparsegraph *sg)
{
    int n = sg->nv;
    size_t nde = sg->nde, extent, x;
    int v, i, k, maxk, maxkey, nb;

    if (n <= 4) return TRUE;                         /* every graph on <= 4 vertices is planar */
    if (nde > (size_t)0x7FFFFFF0)
    {
        fprintf(stderr, ">E lrplanar_sg: more than 2^31 edges not supported\n");
        exit(1);
    }

    maxk = 3*n - 6;                        /* a planar graph has at most this many edges */
    nb = (int)MIN(nde, (size_t)maxk + 1);  /* at most maxk+1 edges are ever numbered */
    if (nb < 1) nb = 1;

    /* out[] mirrors the layout of sg->e, whose lists may leave gaps */
    extent = 1;
    for (v = 0; v < n; ++v)
    {
        x = sg->v[v] + (size_t)sg->d[v];
        if (x > extent) extent = x;
    }

    DYNALLOC1(int,height,height_sz,n,"lrplanar_sg");
    DYNALLOC1(int,parent_edge,parent_edge_sz,n,"lrplanar_sg");
    DYNALLOC1(int,parent_v,parent_v_sz,n,"lrplanar_sg");
    DYNALLOC1(int,nout,nout_sz,n,"lrplanar_sg");
    DYNALLOC1(int,mark,mark_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_v,frame_v_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_i,frame_i_sz,n,"lrplanar_sg");
    DYNALLOC1(int,frame_n,frame_n_sz,n,"lrplanar_sg");
    DYNALLOC1(int,etgt,etgt_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,esrc,esrc_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,lowpt,lowpt_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,lowpt2,lowpt2_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,nesting,nesting_sz,nb,"lrplanar_sg");
    DYNALLOC1(int,out,out_sz,extent,"lrplanar_sg");

    for (v = 0; v < n; ++v)
    {
        height[v] = NONE; parent_edge[v] = NONE; parent_v[v] = NONE;
        nout[v] = 0; mark[v] = NONE;
    }

    /* phase 1: orientation, one DFS per component */
    k = 0;
    for (v = 0; v < n; ++v)
        if (height[v] == NONE)
        {
            k = dfs1(sg, v, k, maxk);
            if (k < 0) return FALSE;                /* more than 3n-6 distinct edges */
        }
    if (k < 9) return TRUE;                          /* K5 has 10 edges, K3,3 has 9 */

    /* sort every vertex's outgoing edges by nesting depth: one stable
       counting sort over all edges, then redistribute into out[] */
    maxkey = 2*n + 1;
    DYNALLOC1(int,count,count_sz,maxkey + 1,"lrplanar_sg");
    DYNALLOC1(int,sorted,sorted_sz,k,"lrplanar_sg");
    for (i = 0; i <= maxkey; ++i) count[i] = 0;
    for (i = 0; i < k; ++i) ++count[nesting[i]];
    for (i = 1; i <= maxkey; ++i) count[i] += count[i-1];
    for (i = k - 1; i >= 0; --i) sorted[--count[nesting[i]]] = i;
    for (v = 0; v < n; ++v) nout[v] = 0;
    for (i = 0; i < k; ++i)
    {
        int ei = sorted[i];
        v = esrc[ei];
        out[sg->v[v] + nout[v]++] = ei;
    }

    /* phase 2: reuse phase-1 arrays that are no longer needed */
    lowpt_edge = esrc;          /* esrc was only needed for the redistribution */
    ref_ = lowpt2;              /* lowpt2 was only needed for nesting depths */
    stack_bottom = nesting;     /* nesting was only needed for the sort */
    for (i = 0; i < k; ++i) { ref_[i] = NONE; lowpt_edge[i] = NONE; }
    DYNALLOC1(cpair,S,S_sz,k + 2,"lrplanar_sg");

    sp = 0;
    for (v = 0; v < n; ++v)
        if (parent_edge[v] == NONE && !dfs2(sg, v)) return FALSE;
    return TRUE;
}

/* Dense-format wrapper: converts to a sparsegraph kept between calls. */
boolean
lrplanar_dense(graph *g, int m, int n)
{
    static TLS_ATTR sparsegraph sg;
    static TLS_ATTR boolean init = FALSE;

    if (!init) { SG_INIT(sg); init = TRUE; }
    nauty_to_sg(g, &sg, m, n);
    return lrplanar_sg(&sg);
}
