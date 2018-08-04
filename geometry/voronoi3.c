// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/voronoi3.h"

voronoi3_t* voronoi3_new(point_t* generators,
                         size_t num_generators,
                         polyhedron_t** cells,
                         adj_graph_t* graph)
{
  ASSERT(num_generators > 0);
  voronoi3_t* diagram = polymec_malloc(sizeof(voronoi3_t));
  diagram->generators = generators;
  diagram->num_generators = num_generators;
  diagram->cells = cells;
  diagram->graph = graph;
  return diagram;
}

void voronoi3_free(voronoi3_t* diagram)
{
  adj_graph_free(diagram->graph);
  for (size_t i = 0; i < diagram->num_generators; ++i)
    polymec_release(diagram->cells[i]);
  polymec_free(diagram->cells);
  polymec_free(diagram->generators);
  polymec_free(diagram);
}

