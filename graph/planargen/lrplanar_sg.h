/* lrplanar_sg.h -- left-right planarity test for graphs of any size in
 * nauty's sparsegraph representation.  See lrplanar_sg.c. */

#ifndef LRPLANAR_SG_H
#define LRPLANAR_SG_H

#include "gtools.h"
#include "nausparse.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TRUE if the graph sg is planar.  Loops and parallel edges are ignored
 * (they do not affect planarity).  With digraph = FALSE, sg must be
 * undirected, i.e. every edge must appear in the lists of both endpoints.
 * With digraph = TRUE the lists hold arcs and the underlying undirected
 * graph is tested (u-v is an edge iff u->v or v->u is an arc).  No limit
 * on size other than memory: vertex numbers are int, edge counts size_t. */
boolean lrplanar_sg(sparsegraph *sg, boolean digraph);

/* The same for a graph in nauty's dense format (m setwords per row). */
boolean lrplanar_dense(graph *g, int m, int n, boolean digraph);

/* Release the work space kept between calls. */
void lrplanar_freedyn(void);

#ifdef __cplusplus
}
#endif

#endif
