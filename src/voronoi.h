#ifndef ARBI_VORONOI_H
#define ARBI_VORONOI_H

#include "mesh.h"
#include "point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This type represents a planar Voronoi graph.
typedef struct planar_voronoi_t planar_voronoi_t;

// This function creates a planar Voronoi graph using the Triangle library.
planar_voronoi_t* planar_voronoi_new(point_t* points, int num_points, bbox_t* bounding_box);

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space, returning a fully-featured mesh.
mesh_t* voronoi_tessellation(point_t* points, int num_points, bbox_t* bounding_box);

// This function produces a 3D mesh by linearly extruding the given planar graph 
// in the z direction.
mesh_t* extrusion(planar_voronoi_t* planar_graph, int num_extruded_cells, double length);

#ifdef __cplusplus
}
#endif

#endif

