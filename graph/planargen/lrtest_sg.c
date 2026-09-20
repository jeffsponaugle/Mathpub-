/* lrtest_sg.c -- test harness for lrplanar_sg.c.
 *
 * Reads graphs of any size (graph6 / sparse6 / digraph6, loops and parallel
 * edges allowed) from stdin and runs up to three planarity testers on each:
 *
 *   N  nauty's own tester (planarity.c, Lieby & McKay; the reference)
 *   S  lrplanar_sg (the new sparsegraph left-right test under test)
 *   B  lrplanar (the original bitmask left-right test), only for simple
 *      graphs with n <= WORDSIZE
 *
 * A digraph6 input is passed to lrplanar_sg with digraph = TRUE; the
 * references are run on the underlying undirected graph, which the harness
 * builds itself (independently of the library's symmetrization), and
 * lrplanar_sg is additionally run on that underlying graph with
 * digraph = FALSE, so both code paths must agree.
 *
 * Every disagreement is reported (with the graph in graph6/sparse6 when it
 * is small) and counted.  Options:
 *
 *   -e p | -e n   every input is expected to be planar / non-planar; a
 *                 different answer counts as a failure
 *   -D            orient each undirected input graph inside the harness
 *                 (every edge becomes one arc in a pseudo-random direction,
 *                 every 7th edge is kept in both directions) and give the
 *                 result to lrplanar_sg with digraph = TRUE; the references
 *                 run on the original.  digraph6 is a dense format, so this
 *                 is the only way to test digraphs with millions of vertices.
 *   -N            skip nauty's tester (for timing the new code alone)
 *   -B            skip the bitmask tester
 *   -q            only print the summary line
 *   -v            print one line per graph (n, edges, answers, times)
 *
 * Exit status 0 iff there were no disagreements and no unexpected answers.
 *
 * usage:  geng -q 8 | ./lrtest_sg
 *         ./biggraphs grid 1000 1000 | ./lrtest_sg -e p
 */

#include "gtools.h"
#include "nausparse.h"
#include NAUTY_PLANARITY_H          /* the original nauty planarity.h, set by the Makefile */
#include "lrplanar_sg.h"
#include "lrplanar.h"

/* nauty's tester, on the same adjacency-list input planarg builds */
static boolean
nauty_isplanar_sg(sparsegraph *sg, int loops)
{
    DYNALLSTAT(t_ver_sparse_rep,V,V_sz);
    DYNALLSTAT(t_adjl_sparse_rep,A,A_sz);
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    int edge_pos, v, w, c, i, k, n = sg->nv;
    size_t j, ne = (sg->nde + loops) / 2;
    boolean ans;

    DYNALLOC1(t_ver_sparse_rep,V,V_sz,n > 0 ? n : 1,"lrtest_sg");
    DYNALLOC1(t_adjl_sparse_rep,A,A_sz,2*ne + 2,"lrtest_sg");
    k = 0;
    for (i = 0; i < n; ++i)
    {
        if (sg->d[i] == 0) { V[i].first_edge = NIL; continue; }
        V[i].first_edge = k;
        for (j = sg->v[i]; j < sg->v[i] + (size_t)sg->d[i]; ++j)
        {
            A[k].end_vertex = sg->e[j]; A[k].next = k + 1; ++k;
            if (A[k-1].end_vertex == i)            /* loops go in twice */
            { A[k].end_vertex = i; A[k].next = k + 1; ++k; }
        }
        A[k-1].next = NIL;
    }
    ans = sparseg_adjl_is_planar(V, n, A, &c, &dfs_tree, &back_edges,
                                 &mult_edges, &embed_graph, &edge_pos, &v, &w);
    sparseg_dlcl_delete(dfs_tree, n);
    sparseg_dlcl_delete(back_edges, n);
    sparseg_dlcl_delete(mult_edges, n);
    embedg_VES_delete(embed_graph, n);
    return ans;
}

/* Underlying undirected graph of a digraph, as a simple graph: for every
 * arc u->w with u != w, an edge u-w (once).  Written for independence from
 * symmetrize() in lrplanar_sg.c: uses sorted lists and a merge instead of
 * marks. */
static sparsegraph *
underlying(sparsegraph *sg)
{
    static SG_DECL(ug);
    DYNALLSTAT(int,mark,mark_sz);
    int n = sg->nv, u, w;
    size_t j, nde = 0, pos;

    DYNALLOC1(int,mark,mark_sz,n > 0 ? n : 1,"lrtest_sg");
    SG_ALLOC(ug, n > 0 ? n : 1, 2*sg->nde + 1, "lrtest_sg");
    for (u = 0; u < n; ++u) { ug.d[u] = 0; mark[u] = -1; }
    /* pass 1: degrees of the simple underlying graph (u-w counted once) */
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            w = sg->e[j];
            if (w == u) continue;
            if (w < u)      /* count only from the smaller endpoint, if arc w->u absent or first */
            {
                size_t jj; boolean seen = FALSE;
                for (jj = sg->v[w]; jj < sg->v[w] + (size_t)sg->d[w]; ++jj)
                    if (sg->e[jj] == u) { seen = TRUE; break; }
                if (seen) continue;             /* counted from w's side */
            }
            /* first occurrence of the pair from this side? (dedupe parallel arcs) */
            if (mark[w] == u) continue;
            mark[w] = u;
            ++ug.d[u]; ++ug.d[w]; ++nde; ++nde;
        }
    pos = 0;
    for (u = 0; u < n; ++u) { ug.v[u] = pos; pos += (size_t)ug.d[u]; ug.d[u] = 0; mark[u] = -1; }
    /* pass 2: fill, same rule */
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            w = sg->e[j];
            if (w == u) continue;
            if (w < u)
            {
                size_t jj; boolean seen = FALSE;
                for (jj = sg->v[w]; jj < sg->v[w] + (size_t)sg->d[w]; ++jj)
                    if (sg->e[jj] == u) { seen = TRUE; break; }
                if (seen) continue;
            }
            if (mark[w] == u) continue;
            mark[w] = u;
            ug.e[ug.v[u] + (size_t)ug.d[u]++] = w;
            ug.e[ug.v[w] + (size_t)ug.d[w]++] = u;
        }
    ug.nv = n; ug.nde = nde;
    return &ug;
}

/* -D: orient the undirected graph sg into a digraph (static work space). */
static sparsegraph *
orient_randomly(sparsegraph *sg)
{
    static SG_DECL(dg);
    static unsigned long long x = 0x9E3779B97F4A7C15ULL;
    unsigned long long x0 = x;               /* pass 2 must replay pass 1's decisions */
    int n = sg->nv, u, w;
    size_t j, k = 0, pos;

    SG_ALLOC(dg, n > 0 ? n : 1, sg->nde + 1, "lrtest_sg");
    for (u = 0; u < n; ++u) dg.d[u] = 0;
    /* each undirected edge is met twice (u<w and w<u); decide once, at u<w */
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            w = sg->e[j];
            if (w < u) continue;
            ++k;
            if (w == u) { ++dg.d[u]; continue; }          /* loop: one arc */
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            if (k % 7 == 0) { ++dg.d[u]; ++dg.d[w]; }
            else if (x & 1) ++dg.d[u];
            else ++dg.d[w];
        }
    pos = 0;
    for (u = 0; u < n; ++u) { dg.v[u] = pos; pos += (size_t)dg.d[u]; dg.d[u] = 0; }
    dg.nv = n; dg.nde = pos;
    x = x0; k = 0;                                        /* replay the same decisions */
    for (u = 0; u < n; ++u)
        for (j = sg->v[u]; j < sg->v[u] + (size_t)sg->d[u]; ++j)
        {
            w = sg->e[j];
            if (w < u) continue;
            ++k;
            if (w == u) { dg.e[dg.v[u] + (size_t)dg.d[u]++] = u; continue; }
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            if (k % 7 == 0) { dg.e[dg.v[u] + (size_t)dg.d[u]++] = w; dg.e[dg.v[w] + (size_t)dg.d[w]++] = u; }
            else if (x & 1) dg.e[dg.v[u] + (size_t)dg.d[u]++] = w;
            else dg.e[dg.v[w] + (size_t)dg.d[w]++] = u;
        }
    return &dg;
}

/* TRUE if sg has no loops or parallel edges */
static boolean
is_simple(sparsegraph *sg, int loops)
{
    DYNALLSTAT(int,mark,mark_sz);
    int i, w, n = sg->nv;
    size_t j;

    if (loops) return FALSE;
    DYNALLOC1(int,mark,mark_sz,n > 0 ? n : 1,"lrtest_sg");
    for (i = 0; i < n; ++i) mark[i] = -1;
    for (i = 0; i < n; ++i)
        for (j = sg->v[i]; j < sg->v[i] + (size_t)sg->d[i]; ++j)
        {
            w = sg->e[j];
            if (w == i || mark[w] == i) return FALSE;
            mark[w] = i;
        }
    return TRUE;
}

int
main(int argc, char *argv[])
{
    SG_DECL(sg);
    sparsegraph *ug;
    int loops, n, m, expected = -1, i;
    boolean doN = TRUE, doB = TRUE, quiet = FALSE, verbose = FALSE, digraph, orientD = FALSE;
    boolean a, b, c, d, simple;
    graph *g;
    nauty_counter nin = 0, nplanar = 0, nbad = 0, nB = 0, ndig = 0;
    double t0, tN = 0, tS = 0, tB = 0;
    size_t maxn = 0, maxe = 0;

    for (i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-N") == 0) doN = FALSE;
        else if (strcmp(argv[i], "-B") == 0) doB = FALSE;
        else if (strcmp(argv[i], "-q") == 0) quiet = TRUE;
        else if (strcmp(argv[i], "-v") == 0) verbose = TRUE;
        else if (strcmp(argv[i], "-D") == 0) orientD = TRUE;
        else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc)
        {
            ++i;
            expected = (argv[i][0] == 'p' || argv[i][0] == 'P') ? 1 : 0;
        }
        else { fprintf(stderr, "usage: lrtest_sg [-e p|n] [-D] [-N] [-B] [-q] [-v] < graphs\n"); return 2; }
    }

    while (read_sgg_loops(stdin, &sg, &loops, &digraph) != NULL)
    {
        ++nin;
        n = sg.nv;
        if ((size_t)n > maxn) maxn = n;
        if (sg.nde / 2 > maxe) maxe = sg.nde / 2;

        if (orientD && !digraph)             /* -D on an undirected input: orient it here */
        {
            sparsegraph *dg = orient_randomly(&sg);
            ++ndig;
            t0 = CPUTIME; b = lrplanar_sg(dg, TRUE); tS += CPUTIME - t0;
            ug = &sg;                        /* the exact underlying graph is the input itself */
            d = lrplanar_sg(ug, FALSE);
            a = b; c = b;
        }
        else
        {
            t0 = CPUTIME; b = lrplanar_sg(&sg, digraph); tS += CPUTIME - t0;
            a = b; c = b; d = b;
            if (digraph)
            {
                ++ndig;
                ug = underlying(&sg);        /* references see the underlying simple graph */
                loops = 0;
                d = lrplanar_sg(ug, FALSE);  /* the library's two code paths must agree */
            }
            else ug = &sg;
        }
        if (doN) { t0 = CPUTIME; a = nauty_isplanar_sg(ug, loops); tN += CPUTIME - t0; }
        simple = FALSE;
        if (doB && n <= WORDSIZE && (simple = is_simple(ug, loops)))
        {
            g = sg_to_nauty(ug, NULL, 1, &m);
            t0 = CPUTIME; c = lr_isplanar(g, n); tB += CPUTIME - t0;
            ++nB;
        }
        if (b) ++nplanar;

        if (a != b || b != c || b != d || (expected >= 0 && b != expected))
        {
            ++nbad;
            if (nbad <= 20)
            {
                fprintf(stderr, "DISAGREEMENT graph %llu: n=%d edges=%lu%s nauty=%d sg=%d%s%s%s",
                        (unsigned long long)nin, n, (unsigned long)(sg.nde/2),
                        (digraph || orientD) ? " (digraph)" : "", a, b,
                        (digraph || orientD) ? (d ? " sg_underlying=1" : " sg_underlying=0") : "",
                        simple ? " bitmask=" : "", simple ? (c ? "1" : "0") : "");
                if (expected >= 0) fprintf(stderr, " expected=%d", expected);
                if (n <= 100) { fprintf(stderr, " : "); if (digraph) writed6_sg(stderr, &sg); else writes6_sg(stderr, &sg); }
                else fprintf(stderr, "\n");
            }
        }
        if (verbose)
            fprintf(stderr, "graph %llu: n=%d edges=%lu%s sg=%d%s\n", (unsigned long long)nin, n,
                    (unsigned long)(sg.nde/2), loops ? " (loops)" : "", b,
                    doN ? (a ? " nauty=1" : " nauty=0") : "");
    }

    if (!quiet || nbad)
    {
        fprintf(stderr, "lrtest_sg: %llu graphs%s (max n=%lu, max edges=%lu), %llu planar, %llu failures",
                (unsigned long long)nin, ndig ? " (digraphs)" : "", (unsigned long)maxn, (unsigned long)maxe,
                (unsigned long long)nplanar, (unsigned long long)nbad);
        if (nin)
        {
            fprintf(stderr, "; cpu/graph: sg %.3f us", 1e6*tS/(double)nin);
            if (doN) fprintf(stderr, ", nauty %.3f us", 1e6*tN/(double)nin);
            if (nB) fprintf(stderr, ", bitmask %.3f us (%llu graphs)", 1e6*tB/(double)nB, (unsigned long long)nB);
        }
        fprintf(stderr, "\n");
    }
    lrplanar_freedyn();
    return nbad != 0;
}
