/* lrtest.c -- compare lr_isplanar against nauty's planarity tester on
 * every graph read from stdin (graph6/sparse6).  Reports mismatches and
 * timing.  usage: geng ... | ./lrtest */
#include "gtools.h"
#include "planarity.h"
#include "lrplanar.h"

static boolean
nauty_isplanar(graph *g, int n)
{
    t_ver_sparse_rep V[MAXN];
    t_adjl_sparse_rep A[MAXN*(MAXN-1) + 1];
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    int edge_pos, v, w, c, i, j, k;
    setword gv;
    boolean ans;

    k = 0;
    for (i = 0; i < n; ++i)
    {
        gv = g[i];
        if (gv == 0) V[i].first_edge = NIL;
        else
        {
            V[i].first_edge = k;
            while (gv) { TAKEBIT(j, gv); A[k].end_vertex = j; A[k].next = k+1; ++k; }
            A[k-1].next = NIL;
        }
    }
    arena_reset();
    ans = sparseg_adjl_is_planar(V, n, A, &c, &dfs_tree, &back_edges,
                                 &mult_edges, &embed_graph, &edge_pos, &v, &w);
    return ans;
}

int
main(int argc, char *argv[])
{
    graph *g;
    int m, n, i;
    nauty_counter nin = 0, nplanar = 0, nmismatch = 0;
    double t0, tn = 0, tl = 0;
    int a, b;
    graph gg[MAXN];

    while ((g = readg(stdin, NULL, 0, &m, &n)) != NULL)
    {
        if (m != 1 || n > MAXN) { fprintf(stderr, "n too large\n"); exit(1); }
        for (i = 0; i < n; ++i) gg[i] = g[i];
        ++nin;
        t0 = CPUTIME; a = nauty_isplanar(gg, n); tn += CPUTIME - t0;
        t0 = CPUTIME; b = lr_isplanar(gg, n);    tl += CPUTIME - t0;
        if (a) ++nplanar;
        if (a != b)
        {
            ++nmismatch;
            if (nmismatch <= 10)
            {
                fprintf(stderr, "MISMATCH nauty=%d lr=%d : ", a, b);
                writeg6(stderr, gg, 1, n);
            }
        }
    }
    fprintf(stderr, "lrtest: " COUNTER_FMT " graphs, " COUNTER_FMT " planar, "
            COUNTER_FMT " mismatches; nauty %.3f us/graph, lr %.3f us/graph\n",
            nin, nplanar, nmismatch, 1e6*tn/(double)nin, 1e6*tl/(double)nin);
    return nmismatch != 0;
}
