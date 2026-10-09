/* planarprune.c -- geng plugin: reject non-planar graphs during generation.

   Planarity is closed under vertex deletion, so an intermediate graph that is
   non-planar can never grow into a planar one; geng's PRUNE hook is designed
   for exactly this kind of property.  The output of geng then consists of
   planar graphs only, so no planarg pass is needed:

       geng_planar -c 14 0:36 $k/$m -u        # a(14) is the sum over chunks

   Build from the nauty source directory (tested with nauty 2.8.9):

   gcc -o geng_planar -O3 -march=native -DMAXN=WORDSIZE -DWORDSIZE=32 \
       -DPRUNE=planarprune -DSUMMARY=planarsummary \
       geng.c planarprune.c planarity.c \
       gtoolsW.o nautyW1.o nautilW1.o naugraphW1.o schreier.o naurng.o

   -DPREPRUNE=planarprune also works (tests candidates before the canonicity
   check) but measured slower at n = 10.  The >P line at exit reports how many
   planarity tests were run.
*/
#include "gtools.h"
#include "planarity.h"

static nauty_counter ntests = 0, nrejected = 0, nshortcut = 0;

int
planarprune(graph *g, int n, int maxn)
/* Return nonzero to reject. Vertex n-1 is the most recently added vertex and
   the graph without it has already passed. */
{
    static t_ver_sparse_rep V[WORDSIZE];
    static t_adjl_sparse_rep A[WORDSIZE*WORDSIZE];
    t_dlcl **dfs_tree, **back_edges, **mult_edges;
    t_ver_edge *embed_graph;
    int c, edge_pos, v, w, i, j, k, ne;
    setword x;
    boolean planar;

    if (n <= 4) return 0;
    if (POPCOUNT(g[n-1]) <= 1) { ++nshortcut; return 0; }   /* pendant or isolated vertex keeps planarity */

    ne = 0;
    for (i = 0; i < n; ++i) ne += POPCOUNT(g[i]);
    ne /= 2;
    if (ne > 3*n - 6) { ++nshortcut; ++nrejected; return 1; }

    k = 0;
    for (i = 0; i < n; ++i)
    {
        x = g[i];
        if (x == 0) { V[i].first_edge = NIL; continue; }
        V[i].first_edge = k;
        while (x)
        {
            TAKEBIT(j, x);
            A[k].end_vertex = j;
            A[k].next = k + 1;
            ++k;
        }
        A[k-1].next = NIL;
    }

    ++ntests;
    planar = sparseg_adjl_is_planar(V, n, A, &c, &dfs_tree, &back_edges,
                                    &mult_edges, &embed_graph, &edge_pos, &v, &w);
    sparseg_dlcl_delete(dfs_tree, n);
    sparseg_dlcl_delete(back_edges, n);
    sparseg_dlcl_delete(mult_edges, n);
    embedg_VES_delete(embed_graph, n);

    if (!planar) ++nrejected;
    return planar ? 0 : 1;
}

void
planarsummary(nauty_counter nout, double cpu)
{
    fprintf(stderr, ">P planarity tests=" COUNTER_FMT " rejected=" COUNTER_FMT
            " shortcuts=" COUNTER_FMT "\n", ntests, nrejected, nshortcut);
}
