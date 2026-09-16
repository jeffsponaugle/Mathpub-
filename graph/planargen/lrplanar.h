/* lrplanar.h -- left-right planarity test for small simple graphs given as
 * nauty bitmask adjacency (m = 1, n <= WORDSIZE).  Returns 1 if planar. */
#ifndef LRPLANAR_H
#define LRPLANAR_H
#include "gtools.h"
int lr_isplanar(const graph *g, int n);
#endif
