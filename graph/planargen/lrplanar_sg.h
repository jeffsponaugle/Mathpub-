/* lrplanar_sg.h -- left-right planarity test for graphs of any size in
 * nauty's sparsegraph representation.  See lrplanar_sg.c. */

#ifndef LRPLANAR_SG_H
#define LRPLANAR_SG_H

#include "gtools.h"
#include "nausparse.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TRUE if the undirected graph sg is planar.  Loops and parallel edges are
 * ignored (they do not affect planarity); sg must be undirected, i.e. every
 * edge must appear in the lists of both endpoints.  No limit on size other
 * than memory; vertex and edge counts must fit in an int. */
boolean lrplanar_sg(sparsegraph *sg);

/* The same for a graph in nauty's dense format (m setwords per row). */
boolean lrplanar_dense(graph *g, int m, int n);

/* Release the work space kept between calls. */
void lrplanar_freedyn(void);

#ifdef __cplusplus
}
#endif

#endif
