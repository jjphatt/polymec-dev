// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/octree.h"

static double square(double x) 
{
  return x*x;
}

static double square_dist(double* x, double* y)
{
  return square(x[0]-y[0]) + square(x[1]-y[1]) + square(x[2]-y[2]);
}

// A bucket PR (point region) octree.

typedef enum
{
  OCTREE_BRANCH_NODE,
  OCTREE_LEAF_NODE
} octree_node_type_t;

typedef struct octree_node_t octree_node_t;

// Branch node.
typedef struct 
{
  // Children.
  octree_node_t* children[8];
} octree_branch_node_t;

// Leaf node.
typedef struct
{
  int index;
  point_t point;
} octree_leaf_node_t;

// Basic octree node -- can be either a branch node or a leaf node.
// We cast between this type and the actual node types below.
struct octree_node_t
{
  octree_node_type_t type;
  union 
  {
    octree_branch_node_t branch_node;
    octree_leaf_node_t leaf_node;
  };
};


static octree_node_t* branch_new()
{
  octree_node_t* node = malloc(sizeof(octree_node_t));
  node->type = OCTREE_BRANCH_NODE;
  memset(node->branch_node.children, 0, 8*sizeof(octree_node_t*));
  return node;
}

// This creates a new leaf node with a single occupant.
static octree_node_t* leaf_new(point_t* point, int index)
{
  octree_node_t* node = malloc(sizeof(octree_node_t));
  node->type = OCTREE_LEAF_NODE;
  node->leaf_node.index = index;
  node->leaf_node.point = *point;
  return node;
}

// This returns the slot within a branch into which the given point fits given 
// a rectangular domain centered at the given center point.
static int find_slot(point_t* center, point_t* point)
{
  if (point->x < center->x)
  {
    if (point->y < center->y)
    {
      if (point->z < center->z)
        return 0;
      else 
        return 1;
    }
    else
    {
      if (point->z < center->z)
        return 2;
      else 
        return 3;
    }
  }
  else
  {
    if (point->y < center->y)
    {
      if (point->z < center->z)
        return 4;
      else 
        return 5;
    }
    else
    {
      if (point->z < center->z)
        return 6;
      else 
        return 7;
    }
  }
}

struct octree_t 
{
  octree_node_t* root;  // The root node.
  bbox_t bbox; // The bounding box for the tree.
  int num_points; // Total number of points in the tree.
};

octree_t* octree_new(bbox_t* bounding_box)
{
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);

  octree_t* tree = malloc(sizeof(octree_t));
  tree->root = NULL;
  tree->bbox = *bounding_box;
  tree->num_points = 0;
  return tree;
}

void octree_free(octree_t* tree)
{
  octree_clear(tree);
  free(tree);
}

void octree_insert(octree_t* tree, point_t* point, int index)
{
  if (tree->root == NULL) // Empty tree
  {
    octree_node_t* node = leaf_new(point, index);
    tree->root = node;
    ++tree->num_points;
  }
  else if (tree->root->type == OCTREE_LEAF_NODE)
  {
    point_t center = {.x = 0.5 * (tree->bbox.x1 + tree->bbox.x2),
                      .y = 0.5 * (tree->bbox.y1 + tree->bbox.y2),
                      .z = 0.5 * (tree->bbox.z1 + tree->bbox.z2)};

    // The tree consists of a single node.
    octree_node_t* root = tree->root;

    // Does the given point already exist here?
    if (point_distance(&root->leaf_node.point, point) == 0.0)
      return;

    // We need to create a branch node here.
    octree_node_t* node = root;
    tree->root = branch_new();
    int slot = find_slot(&center, point);
    tree->root->branch_node.children[slot] = node;
  }
  
  // Now we proceed with the normal logic, given that the root node 
  // is a branch node.
  ASSERT(tree->root->type == OCTREE_BRANCH_NODE);
  octree_node_t* node = tree->root;
  point_t center = {.x = 0.5 * (tree->bbox.x1 + tree->bbox.x2),
                    .y = 0.5 * (tree->bbox.y1 + tree->bbox.y2),
                    .z = 0.5 * (tree->bbox.z1 + tree->bbox.z2)};
  double lx = tree->bbox.x2 - tree->bbox.x1;
  double ly = tree->bbox.y2 - tree->bbox.y1;
  double lz = tree->bbox.z2 - tree->bbox.z1;
  int slot = find_slot(&center, point);
  static double xf[] = {-0.25, -0.25, -0.25, -0.25, +0.25, +0.25, +0.25, +0.25};
  static double yf[] = {-0.25, -0.25, +0.25, +0.25, -0.25, -0.25, +0.25, +0.25};
  static double zf[] = {-0.25, +0.25, -0.25, +0.25, -0.25, +0.25, -0.25, +0.25};
  while ((node->branch_node.children[slot] != NULL) && 
         (node->branch_node.children[slot]->type == OCTREE_BRANCH_NODE))
  {
    node = node->branch_node.children[slot];
    center.x += xf[slot]*lx;
    lx *= 0.5;
    center.y += yf[slot]*ly;
    ly *= 0.5;
    center.z += zf[slot]*lz;
    lz *= 0.5;
    slot = find_slot(&center, point);
  }
  octree_node_t* leaf = node->branch_node.children[slot];
  if (leaf == NULL)
  {
    // No leaf here, so we create a new one!
    leaf = leaf_new(point, index);
    node->branch_node.children[slot] = leaf;
    ++tree->num_points;
  }
  else
  {
    // Is the point already in this node?
    if (point_distance(&leaf->leaf_node.point, point) == 0.0)
        return;
    else
    {
      // We have to make a new branch.
      int old_slot, new_slot; 
      do
      {
        node->branch_node.children[slot] = branch_new();
        node = node->branch_node.children[slot];
        center.x += xf[slot]*lx;
        lx *= 0.5;
        center.y += yf[slot]*ly;
        ly *= 0.5;
        center.z += zf[slot]*lz;
        lz *= 0.5;
        new_slot = find_slot(&center, point);
        old_slot = find_slot(&center, &leaf->leaf_node.point);
      }
      while (new_slot == old_slot);
      node->branch_node.children[old_slot] = leaf;
      octree_node_t* new_leaf = leaf_new(point, index);
      node->branch_node.children[new_slot] = new_leaf;
      ++tree->num_points;
    }
  }
}

void octree_delete(octree_t* tree, point_t* point, int index)
{
  // FIXME
}

int octree_size(octree_t* tree)
{
  return tree->num_points;
}

static void node_clear(octree_node_t* node)
{
  if (node == NULL) 
  {
    return;
  }
  else if (node->type == OCTREE_LEAF_NODE)
  {
    free(node);
  }
  else 
  {
    ASSERT(node->type == OCTREE_BRANCH_NODE);
    for (int i = 0; i < 8; ++i)
      node_clear(node->branch_node.children[i]);
  }
}

void octree_clear(octree_t* tree)
{
  node_clear(tree->root);
  tree->root = NULL;
  tree->num_points = 0;
}

static octree_node_t* find_leaf(octree_node_t* root, 
                                bbox_t* bounding_box, 
                                point_t* x)
{
  if (!bbox_contains(bounding_box, x)) 
    return NULL;

  return NULL;
}

int octree_nearest(octree_t* tree, point_t* point)
{
  if (tree->root == NULL)
    return -1;

  octree_node_t* leaf = find_leaf(tree->root, &tree->bbox, point);
  if (leaf == 0) 
    return -1;
  else return leaf->leaf_node.index;
}

int_slist_t* octree_within_radius(octree_t* tree, 
                                  point_t* point, 
                                  double radius)
{
  return NULL;
#if 0
  if (tree->root == NULL)
    return NULL;

  // Start with the root.
  octree_node_t* node = tree->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = square_dist(pos, node->pos);
  octree_rect_t rect;
  int_slist_t* results = int_slist_new();
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest node.
  find_within_radius(tree->root, pos, radius, &rect, results);
  return results;
#endif
}


