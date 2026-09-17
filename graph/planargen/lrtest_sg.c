/* lrtest_sg.c -- test harness for lrplanar_sg.c.
 *
 * Reads graphs of any size (graph6 / sparse6, loops and parallel edges
 * allowed) from stdin and runs up to three planarity testers on each:
 *
 *   N  nauty's own tester (planarity.c, Lieby & McKay; the reference)
 *   S  lrplanar_sg (the new sparsegraph left-right test under test)
 *   B  lrplanar (the original bitmask left-right test), only for simple
 *      graphs with n <= WORDSIZE
 *
 * Every disagreement is reported (with the graph in graph6/sparse6 when it
 * is small) and counted.  Options:
 *
 *   -e p | -e n   every input is expected to be planar / non-planar; a
 *                 different answer counts as a failure
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
    int loops, n, m, expected = -1, i;
    boolean doN = TRUE, doB = TRUE, quiet = FALSE, verbose = FALSE;
    boolean a, b, c, simple;
    graph *g;
    nauty_counter nin = 0, nplanar = 0, nbad = 0, nB = 0;
    double t0, tN = 0, tS = 0, tB = 0;
    size_t maxn = 0, maxe = 0;

    for (i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-N") == 0) doN = FALSE;
        else if (strcmp(argv[i], "-B") == 0) doB = FALSE;
        else if (strcmp(argv[i], "-q") == 0) quiet = TRUE;
        else if (strcmp(argv[i], "-v") == 0) verbose = TRUE;
        else if (strcmp(argv[i], "-e") == 0 && i + 1 < argc)
        {
            ++i;
            expected = (argv[i][0] == 'p' || argv[i][0] == 'P') ? 1 : 0;
        }
        else { fprintf(stderr, "usage: lrtest_sg [-e p|n] [-N] [-B] [-q] [-v] < graphs\n"); return 2; }
    }

    while (read_sg_loops(stdin, &sg, &loops) != NULL)
    {
        ++nin;
        n = sg.nv;
        if ((size_t)n > maxn) maxn = n;
        if (sg.nde / 2 > maxe) maxe = sg.nde / 2;

        t0 = CPUTIME; b = lrplanar_sg(&sg); tS += CPUTIME - t0;
        a = b; c = b;
        if (doN) { t0 = CPUTIME; a = nauty_isplanar_sg(&sg, loops); tN += CPUTIME - t0; }
        simple = FALSE;
        if (doB && n <= WORDSIZE && (simple = is_simple(&sg, loops)))
        {
            g = sg_to_nauty(&sg, NULL, 1, &m);
            t0 = CPUTIME; c = lr_isplanar(g, n); tB += CPUTIME - t0;
            ++nB;
        }
        if (b) ++nplanar;

        if (a != b || b != c || (expected >= 0 && b != expected))
        {
            ++nbad;
            if (nbad <= 20)
            {
                fprintf(stderr, "DISAGREEMENT graph %llu: n=%d edges=%lu nauty=%d sg=%d%s%s",
                        (unsigned long long)nin, n, (unsigned long)(sg.nde/2), a, b,
                        simple ? " bitmask=" : "", simple ? (c ? "1" : "0") : "");
                if (expected >= 0) fprintf(stderr, " expected=%d", expected);
                if (n <= 100) { fprintf(stderr, " : "); writes6_sg(stderr, &sg); }
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
        fprintf(stderr, "lrtest_sg: %llu graphs (max n=%lu, max edges=%lu), %llu planar, %llu failures",
                (unsigned long long)nin, (unsigned long)maxn, (unsigned long)maxe,
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
