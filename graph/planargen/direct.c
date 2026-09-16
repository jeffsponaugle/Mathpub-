/* direct.c -- the "obvious" planarity prune for geng, kept only as a
 * baseline for comparison with coplanar.c.  PREPRUNE rejects every
 * candidate that is itself non-planar; connectivity is left to geng -c.
 *
 * Compile geng with -DPREPRUNE=direct_preprune and link with planarity.c.
 * Usage:  geng_direct -c -u n 0:3n-6
 */

#include "gtools.h"
#include "planarity.h"

static nauty_counter np_calls, np_planar_tests, np_nonplanar;

int
direct_preprune(graph *g, int n, int maxn)
{
    t_ver_sparse_rep V[MAXN];
    t_adjl_sparse_rep A[MAXN*(MAXN-1) + 1];
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    int edge_pos, v, w, c;
    int i, j, k, ne;
    setword gv;
    boolean ans;

    ++np_calls;
    if (n <= 4) return 0;

    k = 0;
    for (i = 0; i < n; ++i)
    {
        gv = g[i];
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
    ne = k / 2;
    if (ne > 3*n - 6) { ++np_nonplanar; return 1; }

    ++np_planar_tests;
    arena_reset();
    ans = sparseg_adjl_is_planar(V, n, A, &c, &dfs_tree, &back_edges,
                                 &mult_edges, &embed_graph, &edge_pos, &v, &w);

    if (!ans) ++np_nonplanar;
    return !ans;
}

void
direct_summary(nauty_counter nout, double cpu)
{
    fprintf(stderr, ">C preprune calls=" COUNTER_FMT " planarity_tests="
            COUNTER_FMT " nonplanar=" COUNTER_FMT "\n",
            np_calls, np_planar_tests, np_nonplanar);
    if (nout > 0)
        fprintf(stderr, ">C %.3f us cpu per output graph, %.3f us per preprune call\n",
                1e6*cpu/(double)nout, 1e6*cpu/(double)(np_calls ? np_calls : 1));
}
