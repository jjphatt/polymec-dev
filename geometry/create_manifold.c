#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/point.h"
#include "geometry/create_manifold.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function chases down a merged node.
static inline int map_node(int_int_unordered_map_t* node_map, int n)
{
  int* n_p = &n;
  while (n_p != NULL)
  {
    n = *n_p;
    n_p = int_int_unordered_map_get(node_map, n);
  }
  return *n_p;
}

voronoi_tessellation_t* create_manifold(voronoi_tessellation_t* tessellation)
{
  int num_cells = tessellation->num_cells;

  // Find the nodes in the manifold. We do this by finding those nodes that 
  // appear in the faces between every pair of adjacent cells. The tessellation
  // we are given has faces that have two distinct sides that contain 
  // (in general) different sets of nodes. So we create a node mapping from 
  // the original tessellation to the manifold.
  int_int_unordered_map_t* node_map = int_int_unordered_map_new();
  int_int_unordered_map_t* node_pairs = int_int_unordered_map_new();
  int_unordered_set_t* deleted_edges = int_unordered_set_new();
  int_unordered_set_t* deleted_faces = int_unordered_set_new();
  int_unordered_set_t* f1_nodes = int_unordered_set_new();
  int_unordered_set_t* f2_nodes = int_unordered_set_new();
  for (int c1 = 0; c1 < num_cells; ++c1)
  {
    // Go over the faces of this cell.
    int nf1 = tessellation->cells[c1].num_faces;
    for (int fc1 = 0; fc1 < nf1; ++fc1)
    {
      int f1 = tessellation->cells[c1].faces[fc1];
      ASSERT(c1 == tessellation->faces[f1].cell1);

      // Retrieve the face in the opposite cell corresponding to this
      // face.
      int c2 = tessellation->faces[f1].cell2;
      if (c2 < c1) continue; // We've already processed this interaction.

      int nf2 = tessellation->cells[c2].num_faces;
      int f2 = -1;
      for (int fc2 = 0; fc2 < nf2; ++fc2)
      {
        int f = tessellation->cells[c2].faces[fc2];
        if (tessellation->faces[f].cell2 == c1)
        {
          f2 = f;
          break;
        }
      }
      if (f2 == -1)
      {
        // If f2 doesn't exist, that means that the nodes of f1 have 
        // collapsed to a single node in f2's view. That means we must 
        // delete f1 and all of its edges from the tessellation, and 
        // merge all of its nodes with the appropriate node within 
        // c2.
        int_unordered_set_insert(deleted_faces, f1);
        for (int e = 0; e < tessellation->faces[f1].num_edges; ++e)
        {
          int edge = tessellation->faces[f1].edges[e];
          int_unordered_set_insert(deleted_edges, edge);
        }

        // FIXME: This doesn't work yet.
        polymec_error("create_manifold: Face degeneracies are not yet implemented.");
      }

      // Now we have f1 and f2, the two sides of the face between 
      // the cells c1 and c2. We must resolve the nodes of f1 and f2 into 
      // a single set.

      // Create sets of the nodes in the faces, and associate the nodes
      // with their edges.
      int num_edges1 = tessellation->faces[f1].num_edges;
      for (int e = 0; e < num_edges1; ++e)
      {
        int edge = tessellation->faces[f1].edges[e];
        int n1 = tessellation->edges[edge].node1;
        int n2 = tessellation->edges[edge].node2;
        int_unordered_set_insert(f1_nodes, n1);
        int_unordered_set_insert(f1_nodes, n2);
      }
      ASSERT(f1_nodes->size == num_edges1);
      int num_edges2 = tessellation->faces[f2].num_edges;
      for (int e = 0; e < num_edges2; ++e)
      {
        int edge = tessellation->faces[f2].edges[e];
        int n1 = tessellation->edges[edge].node1;
        int n2 = tessellation->edges[edge].node2;
        int_unordered_set_insert(f2_nodes, n1);
        int_unordered_set_insert(f2_nodes, n2);
      }
      ASSERT(f2_nodes->size == num_edges2);

      // Pair up the nodes so that each node is paired with the node closest 
      // to it.
      int pos1 = 0, n1;
      while (int_unordered_set_next(f1_nodes, &pos1, &n1))
      {
        double min_dist = FLT_MAX;
        double* xn1 = &tessellation->nodes[3*n1];
        point_t x1 = {.x = xn1[0], .y = xn1[1], .z = xn1[2]};
        int pos2 = 0, n2;
        while (int_unordered_set_next(f2_nodes, &pos2, &n2))
        {
          double* xn2 = &tessellation->nodes[3*n2];
          point_t x2 = {.x = xn2[0], .y = xn2[1], .z = xn2[2]};
          double dist = point_distance(&x1, &x2);
          if (dist < min_dist)
          {
            // Pair up these nodes.
            int_int_unordered_map_insert(node_pairs, n1, n2);
            int_int_unordered_map_insert(node_pairs, n2, n1);
          }
        }
      }

      // We must get rid of any unpaired nodes, merging them with their 
      // nearest node. A node is unpaired if its closest opposite node is 
      // not paired with it.
      {
        int pos = 0, n1, n2;
        while (int_int_unordered_map_next(node_pairs, &pos, &n1, &n2))
        {
          // Check to see whether n1 and n2 are actually a pair.
          int n2_nearest = *int_int_unordered_map_get(node_pairs, n2);
          if (n2_nearest != n1) 
          {
            // n1 and n2 are evidently not paired. There are two 
            // possibilities:
            // 1. n2 is paired with another node and n1 is unpaired.
            // 2. Both n1 and n2 are unpaired.
            // In either case, n1 must be mapped to (merged with) n2, and 
            // any edge containing n1 must be deleted.
            int_int_unordered_map_insert(node_map, n1, n2);
            for (int e = 0; e < num_edges1; ++e)
            {
              int edge = tessellation->faces[f1].edges[e];
              if ((tessellation->edges[edge].node1 == n1) || 
                  (tessellation->edges[edge].node2 == n1))
              {
                int_unordered_set_insert(deleted_edges, edge);
              }
            }

            int n2_nearest_nearest = *int_int_unordered_map_get(node_pairs, n2_nearest);
            if (n2_nearest_nearest != n2)
            {
              // n2 must be mapped to its nearest node.
              int_int_unordered_map_insert(node_map, n2, n2_nearest);

              for (int e = 0; e < num_edges2; ++e)
              {
                int edge = tessellation->faces[f2].edges[e];
                if ((tessellation->edges[edge].node1 == n2) || 
                    (tessellation->edges[edge].node2 == n2))
                {
                  int_unordered_set_insert(deleted_edges, edge);
                }
              }
            }

          }
        }
      }

      // Clear out the book-keeping data structures.
      int_unordered_set_clear(f1_nodes);
      int_unordered_set_clear(f2_nodes);
      int_int_unordered_map_clear(node_pairs);
    }
  }

  // Clean up what we can so far.
  int_unordered_set_free(f1_nodes);
  int_unordered_set_free(f2_nodes);
  int_int_unordered_map_free(node_pairs);

  // At this point, we have a mapping from the nodes in the original 
  // tessellation to those in the manifold, using the original index
  // space.
  int num_faces = tessellation->num_faces - deleted_faces->size;
  ASSERT(num_faces > 0);
  int num_edges = tessellation->num_edges - deleted_edges->size;
  ASSERT(num_edges > 0);
  int num_nodes = tessellation->num_nodes - node_map->size;
  ASSERT(num_nodes > 0);

  // Create the manifold.
  voronoi_tessellation_t* m = voronoi_tessellation_new(num_cells, num_faces, num_edges, num_nodes);

  // Initialize some data fields that we'll use.
  for (int f = 0; f < m->num_faces; ++f)
    m->faces[f].cell1 = m->faces[f].cell2 = -1;
  for (int e = 0; e < m->num_edges; ++e)
    m->edges[e].node1 = m->edges[e].node2 = -1;

  // Set up the nodes.
  int node_index = 0;
  for (int n = 0; n < tessellation->num_nodes; ++n)
  {
    // Skip the mapped nodes.
    if (int_int_unordered_map_contains(node_map, n))
      continue;

    m->nodes[3*node_index] = tessellation->nodes[3*n];
    m->nodes[3*node_index+1] = tessellation->nodes[3*n+1];
    m->nodes[3*node_index+2] = tessellation->nodes[3*n+2];

    // Map the node in the given tessellation to that in the manifold.
    ASSERT(!int_int_unordered_map_contains(node_map, n));
    int_int_unordered_map_insert(node_map, n, node_index);

    // Move along.
    ++node_index;
  }
  ASSERT(node_index == m->num_nodes);

  // Now wire up the cells and faces.
  int face_index = 0;
  for (int c = 0; c < m->num_cells; ++c)
  {
    // Count up the faces that weren't deleted.
    int num_faces = 0;
    int nf = tessellation->cells[c].num_faces;
    for (int f = 0; f < nf; ++f)
    {
      int face = tessellation->cells[c].faces[f];
      if (!int_unordered_set_contains(deleted_faces, face))
        ++num_faces;
    }
    m->cells[c].num_faces = num_faces;
    m->cells[c].faces = malloc(sizeof(int)*num_faces);

    // Now create the faces.
    int cell_face_index = 0;
    for (int f = 0; f < nf; ++f)
    {
      int face = tessellation->cells[c].faces[f];
      if (!int_unordered_set_contains(deleted_faces, face))
      {
        // Hook the face up to the cell.
        m->cells[c].faces[cell_face_index] = face_index;

        // Hook the cell up to the face.
        if (m->faces[face_index].cell1 == -1)
          m->faces[face_index].cell1 = c;
        else
          m->faces[face_index].cell2 = c;

        // Move everything along.
        ++cell_face_index;
        ++face_index;
      }
    }
  }
  ASSERT(face_index == m->num_faces);

  // Wire up the faces and edges.
  face_index = 0;
  int edge_index = 0;
  for (int f = 0; f < tessellation->num_faces; ++f)
  {
    // Skip deleted faces.
    if (int_unordered_set_contains(deleted_faces, f))
      continue;

    // Count up the edges that weren't deleted for this face.
    int num_edges = 0;
    int ne = tessellation->faces[f].num_edges;
    for (int e = 0; f < ne; ++e)
    {
      int edge = tessellation->faces[f].edges[e];
      if (!int_unordered_set_contains(deleted_edges, edge))
        ++num_edges;
    }
    m->faces[face_index].num_edges = num_edges;
    m->faces[face_index].edges = malloc(sizeof(int)*num_edges);

    // Now create the edges.
    int face_edge_index = 0;
    for (int e = 0; e < ne; ++e)
    {
      int edge = tessellation->faces[f].edges[e];
      if (!int_unordered_set_contains(deleted_edges, edge))
      {
        // Hook the edge up to the face.
        m->faces[face_index].edges[face_edge_index] = edge_index;

        // Move things along.
        ++face_edge_index;
        ++edge_index;
      }
    }
    ++face_index;
  }
  ASSERT(edge_index == m->num_edges);

  // Finally, wire up the edges and nodes.
  edge_index = 0;
  for (int e = 0; e < tessellation->num_edges; ++e)
  {
    // Skip deleted edges.
    if (int_unordered_set_contains(deleted_edges, e))
      continue;

    // Get the nodes for this edge within the manifold.
    int n1 = tessellation->edges[e].node1;
    m->edges[edge_index].node1 = map_node(node_map, n1);
    int n2 = tessellation->edges[e].node1;
    m->edges[edge_index].node2 = map_node(node_map, n2);

    // Move along.
    ++edge_index;
  }

  // Clean up.
  int_int_unordered_map_free(node_map);
  int_unordered_set_free(deleted_edges);
  int_unordered_set_free(deleted_faces);

  return m;
}

#ifdef __cplusplus
}
#endif

