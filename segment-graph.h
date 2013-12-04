/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_GRAPH
#define SEGMENT_GRAPH

#include <algorithm>
#include <parallel/algorithm>
#include <cmath>
#include <omp.h>
#include "disjoint-set.h"

// threshold function
#define THRESHOLD(size, c) (c/size)

typedef struct {
  float w;
  int a, b;
} edge;

bool operator<(const edge &a, const edge &b) {
  return a.w < b.w;
}

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

void merge(edge *v1, int len1, edge *v2, int len2) {
  /* merge v1 and v2 into v1 -- they are in the same vector */
  edge *v = (edge *) malloc(len1 * sizeof(edge));
  memcpy(v, v1, len1 * sizeof(edge));
  int i, j, k;
  i = j = k = 0;

  while (i < len1 && j < len2) {
    if (v[i] < v2[j]) {
      v1[k++] = v[i++];
    } else {
      v1[k++] = v2[j++];
    }
  }
  while (i < len1) {
    v1[k++] = v[i++];
  }
  while (j < len2) {
    v1[k++] = v2[j++];
  }
  free(v);
}

void my_parallel_sort(edge *values, int len) {
  int nthreads = omp_thread_count();
  int each = len / nthreads;

  // sort chunks
#pragma omp parallel
{
  int rank = omp_get_thread_num();
  edge *end = values + (rank + 1) * each;
  if (rank == nthreads - 1) {
    end = values + len;
  }
  std::sort(values + rank * each, end);
}

  // merge chunks together
  int i;
  int iterations = log2(nthreads);
  int jobs = nthreads / 2;
  for (i = 0; i < iterations; ++i) {
#pragma omp parallel
    {
      int rank = omp_get_thread_num();
      int each = len / jobs / 2;
      if (rank % 2 == 0 && rank < jobs * 2) {
	merge(values + rank * each, each, values + (rank + 1) * each, each);
      }
    }
      jobs /= 2;
  }

}

/*
 * Segment a graph
 *
 * Returns a disjoint-set forest representing the segmentation.
 *
 * num_vertices: number of vertices in graph.
 * num_edges: number of edges in graph
 * edges: array of edges.
 * c: constant for treshold function.
 */
universe *segment_graph(int num_vertices, int num_edges, edge *edges,
			float c) {
  // sort edges by weight
  //my_parallel_sort(edges, num_edges);
  __gnu_parallel::sort(edges, edges + num_edges);

  // make a disjoint-set forest
  universe *u = new universe(num_vertices);

  // init thresholds
  float *threshold = new float[num_vertices];
  for (int i = 0; i < num_vertices; i++)
    threshold[i] = THRESHOLD(1,c);

  // for each edge, in non-decreasing weight order...
  for (int i = 0; i < num_edges; ++i) {
    edge *pedge = &edges[i];

    // components conected by this edge
    int a = u->find(pedge->a);
    int b = u->find(pedge->b);
    if (a != b) {
      if ((pedge->w <= threshold[a]) &&
          (pedge->w <= threshold[b])) {
        u->join(a, b);
        a = u->find(a);
        threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
      }
    }
  }

  // free up
  delete threshold;
  return u;
}

#endif
