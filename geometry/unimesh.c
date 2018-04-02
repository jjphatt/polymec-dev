// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "geometry/unimesh.h"
#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"
#include "geometry/unimesh_field.h"

#if POLYMEC_HAVE_MPI
#include "core/adj_graph.h"
#include "ptscotch.h"
#endif

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

struct unimesh_observer_t
{
  void* context;
  unimesh_observer_vtable vtable;
};

DEFINE_ARRAY(unimesh_observer_array, unimesh_observer_t*)

// This maps patch indices to patch boundary conditions for the unimesh.
DEFINE_UNORDERED_MAP(patch_bc_map, int, unimesh_patch_bc_t**, int_hash, int_equals)

// This stuff allows us to perform local patch boundary updates.
typedef struct boundary_buffer_pool_t boundary_buffer_pool_t;
static boundary_buffer_pool_t* boundary_buffer_pool_new(unimesh_t* mesh);
static void boundary_buffer_pool_free(boundary_buffer_pool_t* pool);

struct unimesh_t
{
  // Bounding box and cell spacings.
  bbox_t bbox;
  real_t dx, dy, dz;

  // Intrinsic metadata.
  int npx, npy, npz, nx, ny, nz;
  bool periodic_in_x, periodic_in_y, periodic_in_z;
  
  // Information about which patches are present.
  int_unordered_set_t* patches;
  int* patch_indices;

  // Patch boundary conditions.
  patch_bc_map_t* patch_bcs;
  boundary_buffer_pool_t* boundary_buffers;
  int boundary_update_token; // current transaction
  int_ptr_unordered_map_t* boundary_updates; // maps tokens to patch arrays

  // Special, custom-made boundary conditions, owned and managed by the mesh.
  unimesh_patch_bc_t* copy_bc;
  unimesh_patch_bc_t* periodic_bc;
  unimesh_patch_bc_t* remote_bc;

  // Parallel metadata.
  MPI_Comm comm;
  int nproc, rank;
  int_int_unordered_map_t* owner_procs; // maps (patch index, boundary) pairs
                                        // to processes that own them.

  // Observers.
  unimesh_observer_array_t* observers;

  // This flag is set by unimesh_finalize() after a mesh has been assembled.
  bool finalized;
};

unimesh_t* create_empty_unimesh(MPI_Comm comm, bbox_t* bbox,
                                int npx, int npy, int npz, 
                                int nx, int ny, int nz,
                                bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
{
  ASSERT(!bbox_is_empty_set(bbox));
  ASSERT(npx > 0);
  ASSERT(npy > 0);
  ASSERT(npz > 0);
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);

  unimesh_t* mesh = polymec_malloc(sizeof(unimesh_t));
  mesh->bbox = *bbox;
  real_t Lx = bbox->x2 - bbox->x1,
         Ly = bbox->y2 - bbox->y1,
         Lz = bbox->z2 - bbox->z1;
  mesh->dx = Lx / (npx*nx);
  mesh->dy = Ly / (npy*ny);
  mesh->dz = Lz / (npz*nz);
  mesh->npx = npx;
  mesh->npy = npy;
  mesh->npz = npz;
  mesh->nx = nx;
  mesh->ny = ny;
  mesh->nz = nz;
  mesh->periodic_in_x = periodic_in_x;
  mesh->periodic_in_y = periodic_in_y;
  mesh->periodic_in_z = periodic_in_z;
  mesh->patches = int_unordered_set_new();
  mesh->patch_indices = NULL;
  mesh->patch_bcs = patch_bc_map_new();
  mesh->boundary_buffers = NULL;
  mesh->boundary_updates = int_ptr_unordered_map_new();
  mesh->boundary_update_token = -1;
  mesh->copy_bc = NULL;
  mesh->periodic_bc = NULL;
  mesh->remote_bc = NULL;
  mesh->comm = comm;
  MPI_Comm_rank(comm, &mesh->rank);
  MPI_Comm_size(comm, &mesh->nproc);
  mesh->owner_procs = int_int_unordered_map_new();
  mesh->observers = unimesh_observer_array_new();
  mesh->finalized = false;
  return mesh;
}

static inline int patch_index(unimesh_t* mesh, int i, int j, int k)
{
  return mesh->npy*mesh->npz*i + mesh->npz*j + k;
}

void unimesh_insert_patch(unimesh_t* mesh, int i, int j, int k)
{
  ASSERT(!mesh->finalized);
  ASSERT(i >= 0);
  ASSERT(i < mesh->npx);
  ASSERT(j >= 0);
  ASSERT(j < mesh->npy);
  ASSERT(k >= 0);
  ASSERT(k < mesh->npz);
  int index = patch_index(mesh, i, j, k);
  ASSERT(!int_unordered_set_contains(mesh->patches, index));
  int_unordered_set_insert(mesh->patches, index);
}

static void patch_bcs_free(unimesh_patch_bc_t** bcs)
{
  for (int b = 0; b < 6; ++b)
  {
    if (bcs[b] != NULL)
      polymec_release(bcs[b]);
  }
  polymec_free(bcs);
}

static void unimesh_set_patch_bc(unimesh_t* mesh,
                                 int i, int j, int k,
                                 unimesh_boundary_t patch_boundary,
                                 unimesh_patch_bc_t* patch_bc)
{
  ASSERT(unimesh_has_patch(mesh, i, j, k));
  int index = patch_index(mesh, i, j, k);
  unimesh_patch_bc_t*** bcs_p = patch_bc_map_get(mesh->patch_bcs, index);
  unimesh_patch_bc_t** bcs;
  if (bcs_p == NULL)
  {
    bcs = polymec_calloc(6*sizeof(unimesh_patch_bc_t*));
    patch_bc_map_insert_with_v_dtor(mesh->patch_bcs, index, bcs, patch_bcs_free);
  }
  else
    bcs = *bcs_p;
  polymec_retain(patch_bc);
  int b = (int)patch_boundary;
  bcs[b] = patch_bc;
}

static void set_up_patch_bcs(unimesh_t* mesh);
void unimesh_finalize(unimesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(!mesh->finalized);

  // Make an array of indices for locally-present patches.
  mesh->patch_indices = polymec_malloc(sizeof(int) * 3 * mesh->patches->size);
  int l = 0;
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        int index = patch_index(mesh, i, j, k);
        if (int_unordered_set_contains(mesh->patches, index))
        {
          mesh->patch_indices[3*l]   = i;
          mesh->patch_indices[3*l+1] = j;
          mesh->patch_indices[3*l+2] = k;
          ++l;
        }
      }
    }
  }
  ASSERT(l == mesh->patches->size);
  mesh->finalized = true;

  // Set up boundary buffers.
  mesh->boundary_buffers = boundary_buffer_pool_new(mesh);

  // Now make sure every patch has a set of boundary conditions.
  set_up_patch_bcs(mesh);

  STOP_FUNCTION_TIMER();
}

static int naive_rank_for_patch(unimesh_t* mesh, 
                                int start_patch_for_proc[mesh->nproc+1], 
                                int i, int j, int k)
{
  int rank = mesh->rank;
  if ((i < 0) && mesh->periodic_in_x)
    i = mesh->npx - 1;
  else if ((i >= mesh->npx) && mesh->periodic_in_x)
    i = 0;
  if ((j < 0) && mesh->periodic_in_y)
    j = mesh->npy - 1;
  else if ((j >= mesh->npy) && mesh->periodic_in_y)
    j = 0;
  if ((k < 0) && mesh->periodic_in_z)
    k = mesh->npz - 1;
  else if ((k >= mesh->npz) && mesh->periodic_in_z)
    k = 0;
  if ((i >= 0) && (i < mesh->npx) && 
      (j >= 0) && (j < mesh->npy) &&
      (k >= 0) && (k < mesh->npz))
  {
    int index = patch_index(mesh, i, j, k);
    rank = (int)(int_lower_bound(start_patch_for_proc, mesh->nproc+1, index));
    ASSERT(rank <= mesh->nproc);
    if (index != start_patch_for_proc[rank]) --rank;
  }
  return rank;
}

static void do_naive_partitioning(unimesh_t* mesh)
{
  START_FUNCTION_TIMER();
  // Total up the number of patches and allocate them to available 
  // processes.
  int npx = mesh->npx, npy = mesh->npy, npz = mesh->npz;
  int num_patches = npx * npy * npz;
  int num_local_patches = num_patches / mesh->nproc;

  // This is the index of the first patch stored on the local process.
  // The patches are alloted to ranks 0 thru nproc-1, with the first 
  // num_local_patches assigned to rank 0, the next num_local_patches 
  // to rank 1, and so on till all the patches are assigned.
  int start_patch_for_proc[mesh->nproc+1];
  start_patch_for_proc[0] = 0;
  for (int p = 0; p < mesh->nproc; ++p)
    start_patch_for_proc[p+1] = start_patch_for_proc[p] + num_local_patches;
  start_patch_for_proc[mesh->nproc] = num_patches;

  // We allocate patches to our own process, and track the processes that
  // own neighboring patches.
  for (int i = 0; i < npx; ++i)
  {
    for (int j = 0; j < npy; ++j)
    {
      for (int k = 0; k < npz; ++k)
      {
        // Which processes own this patch and its neighbors?
        int my_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i, j, k);
        if (my_rank == mesh->rank)
        {
          int my_index = patch_index(mesh, i, j, k);

          // Insert this patch locally.
          unimesh_insert_patch(mesh, i, j, k);

          // x1 boundary
          int x1_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i-1, j, k);
          if (x1_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index, x1_rank);

          // x2 boundary
          int x2_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i+1, j, k);
          if (x2_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index+1, x2_rank);

          // y1 boundary
          int y1_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i, j-1, k);
          if (y1_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index+2, y1_rank);

          // y2 boundary
          int y2_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i, j+1, k);
          if (y2_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index+3, y2_rank);

          // z1 boundary
          int z1_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i, j, k-1);
          if (z1_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index+4, z1_rank);

          // z2 boundary
          int z2_rank = naive_rank_for_patch(mesh, start_patch_for_proc, i, j, k+1);
          if (z2_rank != mesh->rank)
            int_int_unordered_map_insert(mesh->owner_procs, 6*my_index+5, z2_rank);
        }
      }
    }
  }
  STOP_FUNCTION_TIMER();
}

unimesh_t* unimesh_new(MPI_Comm comm, bbox_t* bbox,
                       int npx, int npy, int npz, 
                       int nx, int ny, int nz,
                       bool periodic_in_x, bool periodic_in_y, bool periodic_in_z)
{
  START_FUNCTION_TIMER();
  unimesh_t* mesh = create_empty_unimesh(comm, bbox, 
                                         npx, npy, npz,
                                         nx, ny, nz,
                                         periodic_in_x, periodic_in_y, periodic_in_z);

  if (comm == MPI_COMM_SELF) // every proc gets all patches!
  {
    for (int i = 0; i < npx; ++i)
      for (int j = 0; j < npy; ++j)
        for (int k = 0; k < npz; ++k)
          unimesh_insert_patch(mesh, i, j, k);
  }
  else // we do a naive allotment
    do_naive_partitioning(mesh);
  
  // Finalize and send 'er off.
  unimesh_finalize(mesh);
  STOP_FUNCTION_TIMER();
  return mesh;
}

void unimesh_free(unimesh_t* mesh)
{
  unimesh_observer_array_free(mesh->observers);
  int_int_unordered_map_free(mesh->owner_procs);
  int_ptr_unordered_map_free(mesh->boundary_updates);
  if (mesh->boundary_buffers != NULL)
    boundary_buffer_pool_free(mesh->boundary_buffers);
  patch_bc_map_free(mesh->patch_bcs);
  if (mesh->copy_bc != NULL)
    polymec_release(mesh->copy_bc);
  if (mesh->periodic_bc != NULL)
    polymec_release(mesh->periodic_bc);
  if (mesh->remote_bc != NULL)
    polymec_release(mesh->remote_bc);
  int_unordered_set_free(mesh->patches);
  if (mesh->patch_indices != NULL)
    polymec_free(mesh->patch_indices);
  polymec_free(mesh);
}

MPI_Comm unimesh_comm(unimesh_t* mesh)
{
  return mesh->comm;
}

bbox_t* unimesh_bbox(unimesh_t* mesh)
{
  return &(mesh->bbox);
}

void unimesh_get_spacings(unimesh_t* mesh, 
                          real_t* dx, real_t* dy, real_t* dz)
{
  *dx = mesh->dx;
  *dy = mesh->dy;
  *dz = mesh->dz;
}

bool unimesh_next_patch(unimesh_t* mesh, int* pos, 
                        int* i, int* j, int* k,
                        bbox_t* bbox)
{
  ASSERT(mesh->finalized);
  ASSERT(*pos >= 0);
#if POLYMEC_HAVE_OPENMP
  int num_threads = omp_get_num_threads();
  int tid = omp_get_thread_num();
#else
  int num_threads = 1;
  int tid = 0;
#endif
  if (*pos == 0) 
    *pos = tid;
  bool result = (*pos < mesh->patches->size);
  if (result)
  {
    int l = *pos;
    *i = mesh->patch_indices[3*l];
    *j = mesh->patch_indices[3*l+1];
    *k = mesh->patch_indices[3*l+2];
    *pos += num_threads;
    if (bbox != NULL)
    {
      real_t Lx = mesh->nx * mesh->dx,
             Ly = mesh->ny * mesh->dy,
             Lz = mesh->nz * mesh->dz;
      bbox->x1 = mesh->bbox.x1 + (*i) * Lx;
      bbox->x2 = bbox->x1 + Lx;
      bbox->y1 = mesh->bbox.y1 + (*j) * Ly;
      bbox->y2 = bbox->y1 + Ly;
      bbox->z1 = mesh->bbox.z1 + (*k) * Lz;
      bbox->z2 = bbox->z1 + Lz;
    }
  }
  return result;
}

void unimesh_get_extents(unimesh_t* mesh, int* npx, int* npy, int* npz)
{
  *npx = mesh->npx;
  *npy = mesh->npy;
  *npz = mesh->npz;
}

void unimesh_get_patch_size(unimesh_t* mesh, int* nx, int* ny, int* nz)
{
  *nx = mesh->nx;
  *ny = mesh->ny;
  *nz = mesh->nz;
}

int unimesh_num_patches(unimesh_t* mesh)
{
  return mesh->patches->size;
}

void unimesh_get_periodicity(unimesh_t* mesh, 
                             bool* periodic_in_x,
                             bool* periodic_in_y,
                             bool* periodic_in_z)
{
  *periodic_in_x = mesh->periodic_in_x;
  *periodic_in_y = mesh->periodic_in_y;
  *periodic_in_z = mesh->periodic_in_z;
}

bool unimesh_has_patch(unimesh_t* mesh, int i, int j, int k)
{
  if ((i < 0) || (i >= mesh->npx) ||
      (j < 0) || (j >= mesh->npy) ||
      (k < 0) || (k >= mesh->npz))
    return false;
  int index = patch_index(mesh, i, j, k);
  return int_unordered_set_contains(mesh->patches, index);
}

bool unimesh_has_patch_bc(unimesh_t* mesh, int i, int j, int k, 
                          unimesh_boundary_t patch_boundary)
{
  int index = patch_index(mesh, i, j, k);
  unimesh_patch_bc_t*** bcs_p = patch_bc_map_get(mesh->patch_bcs, index);
  if (bcs_p != NULL)
  {
    unimesh_patch_bc_t** bcs = *bcs_p;
    int b = (int)patch_boundary;
    return (bcs[b] != NULL);
  }
  else
    return false;
}

//------------------------------------------------------------------------ 
//                   Patch boundary condition machinery
//------------------------------------------------------------------------ 
// The following functions are not part of the proper API for the unimesh, 
// but need to be exposed to the unimesh_field class to enable patch 
// boundary updates.
//------------------------------------------------------------------------ 

// The boundary buffer class is an annotated blob of memory that stores 
// patch boundary data for patches in a unimesh with data of a given centering
// and number of components.
typedef struct
{
  unimesh_t* mesh;
  unimesh_centering_t centering;
  int nx, ny, nz, nc;
  bool in_use;
  int_int_unordered_map_t* patch_offsets;
  size_t boundary_offsets[6];
  real_t* storage;
} boundary_buffer_t;

static void boundary_buffer_reset(boundary_buffer_t* buffer, 
                                  unimesh_centering_t centering,
                                  int num_components)
{
  ASSERT(num_components > 0);
  ASSERT(!buffer->in_use);

  // Do we need to do anything?
  if ((buffer->centering == centering) && 
      (buffer->nc == num_components))
    return;

  START_FUNCTION_TIMER();
  // Compute buffer offsets based on centering and boundary.
  buffer->centering = centering;
  buffer->nc = num_components;
  int nx = buffer->nx, ny = buffer->ny, nz = buffer->nz, nc = buffer->nc;
  size_t patch_sizes[8] = {2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2) + (nx+2)*(ny+2)), // cells
                           2*nc*(ny*nz + (nx+1)*nz + (nx+1)*ny), // x faces
                           2*nc*((ny+1)*nz + nx*nz + nx*(ny+1)), // y faces
                           2*nc*(ny*(nz+1) + nx*(nz+1) + nx*ny), // z faces
                           2*nc*((ny+1)*(nz+1) + nx*(nz+1) + nx*(ny+1)),     // x edges
                           2*nc*(ny*(nz+1) + (nx+1)*(nz+1) + (nx+1)*ny), // y edges
                           2*nc*((ny+1)*nz + (nx+1)*nz + (nx+1)*(ny+1)), // z edges
                           2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1) + (nx+1)*(ny+1))}; // nodes

  size_t offsets[8][6] =  { // cells (including ghosts for simplicity)
                           {0, nc*(ny+2)*(nz+2), 
                            2*nc*(ny+2)*(nz+2), 2*nc*(ny+2)*(nz+2) + nc*(nx+2)*(nz+2),
                            2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2)), 2*nc*((ny+2)*(nz+2) + (nx+2)*(nz+2)) + nc*(nx+2)*(ny+2)},
                            // x faces
                           {0, nc*ny*nz,
                            2*nc*ny*nz, 2*nc*ny*nz + nc*(nx+1)*nz,
                            2*nc*(ny*nz + (nx+1)*nz), 2*nc*(ny*nz + (nx+1)*nz) + nc*(nx+1)*ny},
                            // y faces
                           {0, nc*(ny+1)*nz,
                            2*nc*(ny+1)*nz, 2*nc*(ny+1)*nz + nc*nx*nz,
                            2*nc*((ny+1)*nz + nx*nz), 2*nc*((ny+1)*nz + nx*nz) + nc*nx*(ny+1)},
                            // z faces
                           {0, nc*ny*(nz+1),
                            2*nc*ny*(nz+1), 2*nc*ny*(nz+1) + nc*nx*(nz+1),
                            2*nc*(ny*(nz+1) + nx*(nz+1)), 2*nc*(ny*(nz+1) + nx*(nz+1)) + nc*nx*ny},
                            // x edges
                           {0, nc*(ny+1)*(nz+1),
                            2*nc*(ny+1)*(nz+1), 2*nc*(ny+1)*(nz+1) + nc*nx*(nz+1),
                            2*nc*((ny+1)*(nz+1) + nx*(nz+1)), 2*nc*((ny+1)*(nz+1) + nx*(nz+1)) + nc*nx*(ny+1)},
                            // y edges
                           {0, nc*ny*(nz+1),
                            2*nc*ny*(nz+1), 2*nc*ny*(nz+1) + nc*(nx+1)*(nz+1),
                            2*nc*(ny*(nz+1) + (nx+1)*(nz+1)), 2*nc*(ny*(nz+1) + (nx+1)*(nz+1)) + nc*(nx+1)*ny},
                            // z edges
                           {0, nc*(ny+1)*nz,
                            2*nc*(ny+1)*nz, 2*nc*(ny+1)*nz + nc*(nx+1)*nz,
                            2*nc*((ny+1)*nz + (nx+1)*nz), 2*nc*((ny+1)*nz + (nx+1)*nz) + nc*(nx+1)*(ny+1)},
                            // nodes
                           {0, nc*(ny+1)*(nz+1),
                            2*nc*(ny+1)*(nz+1), 2*nc*(ny+1)*(nz+1) + nc*(nx+1)*(nz+1),
                            2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1)), 2*nc*((ny+1)*(nz+1) + (nx+1)*(nz+1)) + nc*(nx+1)*(ny+1)}};

  // Now compute offsets.
  int cent = (int)centering;
  int pos = 0, i, j, k;
  size_t last_offset = 0;
  while (unimesh_next_patch(buffer->mesh, &pos, &i, &j, &k, NULL))
  {
    int index = patch_index(buffer->mesh, i, j, k);
    int_int_unordered_map_insert(buffer->patch_offsets, index, (int)last_offset);
    last_offset += patch_sizes[cent];
  }
  memcpy(buffer->boundary_offsets, offsets[cent], 6*sizeof(size_t));

  // Allocate storage if needed.
  buffer->storage = polymec_realloc(buffer->storage, sizeof(real_t) * last_offset);
  STOP_FUNCTION_TIMER();
}

static boundary_buffer_t* boundary_buffer_new(unimesh_t* mesh, 
                                              unimesh_centering_t centering,
                                              int num_components)
{
  boundary_buffer_t* buffer = polymec_malloc(sizeof(boundary_buffer_t));
  buffer->mesh = mesh;
  buffer->centering = centering;
  unimesh_get_patch_size(mesh, &buffer->nx, &buffer->ny, &buffer->nz);
  buffer->nc = -1;
  buffer->in_use = false;
  buffer->patch_offsets = int_int_unordered_map_new();
  buffer->storage = NULL;
  boundary_buffer_reset(buffer, centering, num_components);
  return buffer;
}

static void boundary_buffer_free(boundary_buffer_t* buffer)
{
  if (buffer->storage != NULL)
    polymec_free(buffer->storage);
  int_int_unordered_map_free(buffer->patch_offsets);
  polymec_free(buffer);
}

static inline void* boundary_buffer_data(boundary_buffer_t* buffer,
                                         int i, int j, int k,
                                         unimesh_boundary_t boundary)
{
  int index = patch_index(buffer->mesh, i, j, k);
  int b = (int)boundary;
  size_t offset = *int_int_unordered_map_get(buffer->patch_offsets, index) + 
                  buffer->boundary_offsets[b];
  return &(buffer->storage[offset]);
}

DEFINE_ARRAY(boundary_buffer_array, boundary_buffer_t*)

// The boundary_buffer_pool class maintains a set of resources for supporting
// local patch boundary updates on the mesh. Specifically, the pool allows 
// several concurrent patch boundary updates for asynchronous communication 
// over several unimesh_fields.
struct boundary_buffer_pool_t
{
  unimesh_t* mesh;
  boundary_buffer_array_t* buffers;
};

// Creates a new boundary update pool.
static boundary_buffer_pool_t* boundary_buffer_pool_new(unimesh_t* mesh)
{
  boundary_buffer_pool_t* pool = polymec_malloc(sizeof(boundary_buffer_pool_t));
  pool->mesh = mesh;
  pool->buffers = boundary_buffer_array_new();

  // Start off with a handful of single-component, cell-centered buffers.
  for (int i = 0; i < 4; ++i)
  {
    boundary_buffer_t* buffer = boundary_buffer_new(mesh, UNIMESH_CELL, 1);
    boundary_buffer_array_append_with_dtor(pool->buffers, buffer, boundary_buffer_free);
  }
  
  return pool;
}

// Destroys the given boundary update pool.
static void boundary_buffer_pool_free(boundary_buffer_pool_t* pool)
{
  boundary_buffer_array_free(pool->buffers);
  polymec_free(pool);
}

// Returns an integer token that uniquely identifies a set of resources
// that can be used for patch boundary updates for data with the given 
// centering.
static int boundary_buffer_pool_acquire(boundary_buffer_pool_t* pool,
                                        unimesh_centering_t centering,
                                        int num_components)
{
  START_FUNCTION_TIMER();
  ASSERT(num_components > 0);

  size_t token = 0; 
  while (token < pool->buffers->size)
  {
    boundary_buffer_t* buffer = pool->buffers->data[token];
    if (!buffer->in_use)
    {
      // Repurpose this buffer if needed.
      boundary_buffer_reset(buffer, centering, num_components);
      buffer->in_use = true;
      break;
    }
    else
      ++token;
  }

  if (token == pool->buffers->size) // We're out of buffers!
  {
    // Add another one.
    boundary_buffer_t* buffer = boundary_buffer_new(pool->mesh, centering, num_components);
    boundary_buffer_array_append_with_dtor(pool->buffers, buffer, boundary_buffer_free);
  }

  STOP_FUNCTION_TIMER();
  return (int)token;
}

// Releases the resources associated with the given integer token, returning 
// them to the pool for later use.
static inline void boundary_buffer_pool_release(boundary_buffer_pool_t* pool,
                                                int token)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  pool->buffers->data[token]->in_use = false;
}

// Retrieves a buffer for the given token that stores patch boundary data 
// for the given boundary on patch (i, j, k). 
static inline void* boundary_buffer_pool_buffer(boundary_buffer_pool_t* pool,
                                                int token,
                                                int i, int j, int k,
                                                unimesh_boundary_t boundary)
{
  ASSERT(token >= 0);
  ASSERT((size_t)token < pool->buffers->size);
  return boundary_buffer_data(pool->buffers->data[token], i, j, k, boundary);
}

// Returns a unique token that can be used to identify a patch boundary 
// update operation on patches with the given centering, so that boundary 
// conditions can be enforced asynchronously.
int unimesh_patch_boundary_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components);
int unimesh_patch_boundary_buffer_token(unimesh_t* mesh, 
                                        unimesh_centering_t centering,
                                        int num_components)
{
  // Acquire a buffer and a token.
  return boundary_buffer_pool_acquire(mesh->boundary_buffers, 
                                      centering, num_components); 
}

// This allows access to the buffer that stores data for the specific boundary 
// of the (i, j, k)th patch in the current transaction.
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary);
void* unimesh_patch_boundary_buffer(unimesh_t* mesh, 
                                    int i, int j, int k, 
                                    unimesh_boundary_t boundary)
{
  ASSERT(mesh->boundary_update_token != -1);
  return boundary_buffer_pool_buffer(mesh->boundary_buffers, 
                                     mesh->boundary_update_token,
                                     i, j, k, boundary);
}

// A boundary update is a set of data that allows us to finish processing
// patch boundary updates in progress.
typedef struct
{
  int i, j, k;
  real_t t;
  unimesh_boundary_t boundary;
  unimesh_patch_t* patch;
} boundary_update_t;

static boundary_update_t* boundary_update_new(int i, int j, int k, real_t t,
                                              unimesh_boundary_t boundary,
                                              unimesh_patch_t* patch)
{
  ASSERT(patch != NULL);
  boundary_update_t* update = polymec_malloc(sizeof(boundary_update_t));
  update->i = i;
  update->j = j;
  update->k = k;
  update->t = t;
  update->boundary = boundary;
  update->patch = patch;
  return update;
}

static void boundary_update_free(boundary_update_t* update)
{
  polymec_free(update);
}

DEFINE_ARRAY(boundary_update_array, boundary_update_t*)

// This provides access to the mesh's current boundary update token.
int unimesh_boundary_update_token(unimesh_t* mesh);
int unimesh_boundary_update_token(unimesh_t* mesh)
{
  return mesh->boundary_update_token;
}

// This starts updating the given patch using boundary conditions in the mesh
// at the given time, tracking the transaction with the given token.
void unimesh_start_updating_patch_boundary(unimesh_t* mesh, int token,
                                           int i, int j, int k, real_t t,
                                           unimesh_boundary_t boundary,
                                           unimesh_patch_t* patch);
void unimesh_start_updating_patch_boundary(unimesh_t* mesh, int token,
                                           int i, int j, int k, real_t t,
                                           unimesh_boundary_t boundary,
                                           unimesh_patch_t* patch)
{
  START_FUNCTION_TIMER();
  ASSERT(mesh->finalized);
  ASSERT(unimesh_has_patch(mesh, i, j, k));
  int index = patch_index(mesh, i, j, k);

  // Fetch the boundary condition.
  int b = (int)boundary;
  unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];
  ASSERT(bc != NULL);

  // Set the boundary update token.
  mesh->boundary_update_token = token;

  // Start the update.
  unimesh_patch_bc_start_update(bc, i, j, k, t, boundary, patch);

  // Stash information for this patch in our boundary updates.
  boundary_update_array_t** updates_p = (boundary_update_array_t**)int_ptr_unordered_map_get(mesh->boundary_updates, token);
  boundary_update_array_t* updates;
  if (updates_p == NULL)
  {
    updates = boundary_update_array_new();
    int_ptr_unordered_map_insert_with_v_dtor(mesh->boundary_updates, token, 
                                             updates, DTOR(boundary_update_array_free));
  }
  else
    updates = *updates_p;
  boundary_update_t* update = boundary_update_new(i, j, k, t, boundary, patch);
  boundary_update_array_append_with_dtor(updates, update, boundary_update_free);

  // Inform our observers that we've started this boundary update. 
  for (size_t o = 0; o < mesh->observers->size; ++o)
  {
    unimesh_observer_t* obs = mesh->observers->data[o];
    if (obs->vtable.started_boundary_update != NULL)
    {
      obs->vtable.started_boundary_update(obs->context, mesh, token, 
                                          update->i, update->j, update->k,
                                          update->boundary, update->t, update->patch);
    }
  }

  mesh->boundary_update_token = -1;
  STOP_FUNCTION_TIMER();
}

void unimesh_start_updating_patch_boundaries(unimesh_t* mesh, int token);
void unimesh_start_updating_patch_boundaries(unimesh_t* mesh, int token)
{
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT((size_t)token < mesh->boundary_buffers->buffers->size);

  // Inform our observers that we've started these boundary updates. 
  boundary_buffer_t* buffer = mesh->boundary_buffers->buffers->data[token];
  for (size_t i = 0; i < mesh->observers->size; ++i)
  {
    unimesh_observer_t* obs = mesh->observers->data[i];
    if (obs->vtable.started_boundary_update != NULL)
    {
      obs->vtable.started_boundary_updates(obs->context, mesh, token, 
                                           buffer->centering, buffer->nc);
    }
  }
  STOP_FUNCTION_TIMER();
}

void unimesh_finish_updating_patch_boundaries(unimesh_t* mesh, int token);
void unimesh_finish_updating_patch_boundaries(unimesh_t* mesh, int token)
{
  START_FUNCTION_TIMER();
  ASSERT(token >= 0);
  ASSERT((size_t)token < mesh->boundary_buffers->buffers->size);

  // We're working on this transaction now.
  mesh->boundary_update_token = token;

  // Inform our observers that we're about to finish updating all patch
  // boundaries.
  boundary_buffer_t* buffer = mesh->boundary_buffers->buffers->data[token];
  for (size_t i = 0; i < mesh->observers->size; ++i)
  {
    unimesh_observer_t* obs = mesh->observers->data[i];
    if (obs->vtable.about_to_finish_boundary_updates != NULL)
    {
      obs->vtable.about_to_finish_boundary_updates(obs->context, mesh, token, 
                                                   buffer->centering, buffer->nc);
    }
  }

  // Go over the patches that correspond to this token.
  boundary_update_array_t* updates = *((boundary_update_array_t**)int_ptr_unordered_map_get(mesh->boundary_updates, token));
  for (size_t i = 0; i < updates->size; ++i)
  {
    boundary_update_t* update = updates->data[i];
    int index = patch_index(mesh, update->i, update->j, update->k);
    int b = (int)update->boundary;

    // Inform our observers that we're about to finish updating this particular
    // patch boundary.
    for (size_t o = 0; o < mesh->observers->size; ++o)
    {
      unimesh_observer_t* obs = mesh->observers->data[o];
      if (obs->vtable.about_to_finish_boundary_update != NULL)
      {
        obs->vtable.about_to_finish_boundary_update(obs->context, mesh, token, 
                                                    update->i, update->j, update->k,
                                                    update->boundary, update->t, update->patch);
      }
    }

    unimesh_patch_bc_t* bc = (*patch_bc_map_get(mesh->patch_bcs, index))[b];
    unimesh_patch_bc_finish_update(bc, update->i, update->j, update->k, 
                                   update->t, update->boundary, update->patch);

    // Inform our observers that we've just finished updating this particular
    // patch boundary.
    for (size_t o = 0; o < mesh->observers->size; ++o)
    {
      unimesh_observer_t* obs = mesh->observers->data[o];
      if (obs->vtable.finished_boundary_update != NULL)
      {
        obs->vtable.finished_boundary_update(obs->context, mesh, token, 
                                             update->i, update->j, update->k, 
                                             update->boundary, update->t, update->patch);
      }
    }
  }

  // Inform our observers that we're finished with this boundary update. 
  for (size_t i = 0; i < mesh->observers->size; ++i)
  {
    unimesh_observer_t* obs = mesh->observers->data[i];
    if (obs->vtable.finished_boundary_update != NULL)
    {
      obs->vtable.finished_boundary_updates(obs->context, mesh, token, 
                                            buffer->centering, buffer->nc);
    }
  }

  // Clear the updates array.
  boundary_update_array_clear(updates);

  // Release the boundary update corresponding to this token.
  boundary_buffer_pool_release(mesh->boundary_buffers, token);
  mesh->boundary_update_token = -1;
  STOP_FUNCTION_TIMER();
}

extern unimesh_patch_bc_t* unimesh_copy_bc_new(unimesh_t* mesh);
extern unimesh_patch_bc_t* unimesh_periodic_bc_new(unimesh_t* mesh);
extern unimesh_patch_bc_t* unimesh_remote_bc_new(unimesh_t* mesh);
static void set_up_patch_bcs(unimesh_t* mesh)
{
  // Here's some ready-made boundary conditions.
  mesh->copy_bc = unimesh_copy_bc_new(mesh);
  if (mesh->periodic_in_x || mesh->periodic_in_y || mesh->periodic_in_z)
    mesh->periodic_bc = unimesh_periodic_bc_new(mesh);
  if (mesh->nproc > 1)
    mesh->remote_bc = unimesh_remote_bc_new(mesh);

  // Assign these to each patch in the mesh.
  int pos = 0, i, j, k;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, NULL))
  {
    unimesh_patch_bc_t *x1_bc = NULL, 
                       *x2_bc = NULL, 
                       *y1_bc = NULL, 
                       *y2_bc = NULL, 
                       *z1_bc = NULL, 
                       *z2_bc = NULL;

    // x boundaries
    if (unimesh_has_patch(mesh, i-1, j, k))
      x1_bc = mesh->copy_bc;
    else if (mesh->periodic_in_x && (i == 0))
    {
      if (unimesh_has_patch(mesh, mesh->npx-1, j, k))
        x1_bc = mesh->periodic_bc;
      else
        x1_bc = mesh->remote_bc;
    }
    else if (i > 0)
      x1_bc = mesh->remote_bc;

    if (unimesh_has_patch(mesh, i+1, j, k))
      x2_bc = mesh->copy_bc;
    else if (mesh->periodic_in_x && (i == mesh->npx-1))
    {
      if (unimesh_has_patch(mesh, 0, j, k))
        x2_bc = mesh->periodic_bc;
      else
        x2_bc = mesh->remote_bc;
    }
    else if (i < mesh->npx-1)
      x2_bc = mesh->remote_bc;

    // y boundaries
    if (unimesh_has_patch(mesh, i, j-1, k))
      y1_bc = mesh->copy_bc;
    else if (mesh->periodic_in_y && (j == 0))
    {
      if (unimesh_has_patch(mesh, i, mesh->npy-1, k))
        y1_bc = mesh->periodic_bc;
      else 
        y1_bc = mesh->remote_bc;
    }
    else if (j > 0)
      y1_bc = mesh->remote_bc;

    if (unimesh_has_patch(mesh, i, j+1, k))
      y2_bc = mesh->copy_bc;
    else if (mesh->periodic_in_y && (j == mesh->npy-1))
    {
      if (unimesh_has_patch(mesh, i, 0, k))
        y2_bc = mesh->periodic_bc;
      else
        y2_bc = mesh->remote_bc;
    }
    else if (j < mesh->npy-1)
      y2_bc = mesh->remote_bc;

    // z boundaries
    if (unimesh_has_patch(mesh, i, j, k-1))
      z1_bc = mesh->copy_bc;
    else if (mesh->periodic_in_z && (k == 0))
    {
      if (unimesh_has_patch(mesh, i, j, mesh->npz-1))
        z1_bc = mesh->periodic_bc;
      else
        z1_bc = mesh->remote_bc;
    }
    else if (k > 0)
      z1_bc = mesh->remote_bc;

    if (unimesh_has_patch(mesh, i, j, k+1))
      z2_bc = mesh->copy_bc;
    else if (mesh->periodic_in_z && (k == mesh->npz-1))
    {
      if (unimesh_has_patch(mesh, i, j, 0))
        z2_bc = mesh->periodic_bc;
      else
        z2_bc = mesh->remote_bc;
    }
    else if (k < mesh->npz-1)
      z2_bc = mesh->remote_bc;

    if (x1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X1_BOUNDARY, x1_bc);
    if (x2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_X2_BOUNDARY, x2_bc);
    if (y1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y1_BOUNDARY, y1_bc);
    if (y2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Y2_BOUNDARY, y2_bc);
    if (z1_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z1_BOUNDARY, z1_bc);
    if (z2_bc != NULL)
      unimesh_set_patch_bc(mesh, i, j, k, UNIMESH_Z2_BOUNDARY, z2_bc);
  }
}

// This returns the process that owns the patch attached to the given 
// boundary of the given local patch (i, j, k).
int unimesh_owner_proc(unimesh_t* mesh, 
                       int i, int j, int k,
                       unimesh_boundary_t boundary);
int unimesh_owner_proc(unimesh_t* mesh, 
                       int i, int j, int k,
                       unimesh_boundary_t boundary)
{
  int index = patch_index(mesh, i, j, k);
  int b = (int)boundary;
  int key = 6*index + b;
  int* proc_p = int_int_unordered_map_get(mesh->owner_procs, key);
  if (proc_p == NULL)
    return mesh->rank;
  else
    return *proc_p;
}

// This provides direct access to our remote BC, in order to simplify access
// to send/receive buffers. Tight coupling much? :-)
unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh);
unimesh_patch_bc_t* unimesh_remote_bc(unimesh_t* mesh)
{
  return mesh->remote_bc;
}

unimesh_observer_t* unimesh_observer_new(void* context,
                                         unimesh_observer_vtable vtable)
{
  unimesh_observer_t* observer = polymec_malloc(sizeof(unimesh_observer_t));
  observer->context = context;
  observer->vtable = vtable;
  return observer;
}

void unimesh_observer_free(unimesh_observer_t* observer)
{
  if ((observer->vtable.dtor != NULL) && (observer->context != NULL))
    observer->vtable.dtor(observer->context);
  polymec_free(observer);
}

void unimesh_add_observer(unimesh_t* mesh,
                          unimesh_observer_t* observer)
{
  for (size_t i = 0; i < mesh->observers->size; ++i)
  {
    if (mesh->observers->data[i] == observer) // already there!
      return;
  }

  unimesh_observer_array_append_with_dtor(mesh->observers, observer,
                                          unimesh_observer_free);
}

void unimesh_remove_observer(unimesh_t* mesh,
                             unimesh_observer_t* observer)
{
  for (size_t i = 0; i < mesh->observers->size; ++i)
  {
    if (mesh->observers->data[i] == observer) 
    {
      unimesh_observer_array_remove(mesh->observers, i);
      unimesh_observer_free(observer);
      return;
    }
  }
}

#if POLYMEC_HAVE_MPI
extern int64_t* partition_graph(adj_graph_t* global_graph, 
                                MPI_Comm comm,
                                int* weights,
                                real_t imbalance_tol);

static adj_graph_t* graph_from_unimesh_patches(unimesh_t* mesh)
{
  return NULL;
}

static void redistribute_unimesh(unimesh_t** mesh, int64_t* partition)
{
}

static void redistribute_unimesh_field(unimesh_field_t* field, int64_t* partition)
{
}
#endif

void repartition_unimesh(unimesh_t** mesh, 
                         int* weights,
                         real_t imbalance_tol,
                         unimesh_field_t** fields,
                         size_t num_fields)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");

  // On a single process, repartitioning has no meaning.
  if ((*mesh)->nproc == 1) 
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // If meshes on rank != 0 are not NULL, we delete them.
  unimesh_t* m = *mesh;
  if ((m->rank != 0) && (m != NULL))
  {
    unimesh_free(m);
    *mesh = m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* graph = graph_from_unimesh_patches(m);

  log_debug("repartition_unimesh: Repartitioning mesh on %d subdomains.", m->nproc);

  // Map the graph to the different domains, producing a partition vector.
  int64_t* partition = (m->rank == 0) ? partition_graph(graph, m->comm, weights, imbalance_tol): NULL;

  // Redistribute the mesh. 
  redistribute_unimesh(mesh, partition);

  // Redistribute the fields.
  for (size_t f = 0; f < num_fields; ++f)
    redistribute_unimesh_field(fields[f], partition);

  // Clean up.
  adj_graph_free(graph);
  polymec_free(partition);

  STOP_FUNCTION_TIMER();
#endif
}
