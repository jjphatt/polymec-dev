// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/polynomial_fit.h"

struct polynomial_fit_t 
{
  char* name;
  void* context;
  polynomial_fit_vtable vtable;
  int num_comps;

  int degree;
  int point_index;
  int num_points;

  // Polynomials for fit.
  polynomial_t** polys;
};

polynomial_fit_t* polynomial_fit_new(const char* name,
                                     void* context,
                                     polynomial_fit_vtable vtable,
                                     int num_comps)
{
  ASSERT(num_comps > 0);
  ASSERT(vtable.num_interior_neighbors != NULL);
  ASSERT(vtable.get_interior_neighbors != NULL);
  ASSERT(vtable.get_interior_data != NULL);
  ASSERT(vtable.num_boundary_neighbors != NULL);
  ASSERT(vtable.get_boundary_neighbors != NULL);
  ASSERT(vtable.get_boundary_data != NULL);
  ASSERT(vtable.targeted_degree != NULL);

  polynomial_fit_t* fit = malloc(sizeof(polynomial_fit_t));
  fit->name = string_dup(name);
  fit->context = context;
  fit->vtable = vtable;
  fit->num_comps = num_comps;

  fit->degree = -1;
  fit->point_index = -1;
  fit->num_points = 0;
  fit->polys = NULL;

  return fit;
}

void polynomial_fit_free(polynomial_fit_t* fit)
{
  if (fit->polys != NULL)
  {
    for (int i = 0; i < fit->num_comps; ++i)
      fit->polys[i] = NULL;
  }
  if ((fit->context != NULL) && (fit->vtable.dtor != NULL))
    fit->vtable.dtor(fit->context);
  if (fit->name != NULL)
    free(fit->name);
  free(fit);
}

int polynomial_fit_num_comps(polynomial_fit_t* fit)
{
  return fit->num_comps;
}

void polynomial_fit_compute(polynomial_fit_t* fit, real_t* data, int point_index)
{
  ASSERT(data != NULL);
  ASSERT(point_index >= 0);
  fit->point_index = point_index;

  // Find the number of neighbors near the given point.
  int num_int_points = 1 + fit->vtable.num_interior_neighbors(fit->context, fit->point_index);
  int num_bnd_points = fit->vtable.num_boundary_neighbors(fit->context, fit->point_index);
  fit->num_points = num_int_points + num_bnd_points;

  // Determine the targeted degree for the number of point
  fit->degree = fit->vtable.targeted_degree(fit->context, fit->num_points);
  int dim = polynomial_basis_dim(fit->degree);
  ASSERT(dim >= fit->num_points);

  // Get the indices of all the points.
  int int_indices[num_int_points];
  int_indices[0] = fit->point_index;
  fit->vtable.get_interior_neighbors(fit->context, fit->point_index, &int_indices[1]);
  int bnd_indices[num_bnd_points];
  fit->vtable.get_boundary_neighbors(fit->context, fit->point_index, bnd_indices);

  // Now fit each component.
  point_t int_points[num_int_points], bnd_points[num_bnd_points];
  vector_t bnd_normals[num_bnd_points];
  void* bnd_conditions[num_bnd_points];

  // Fetch the interior points and values.
  real_t int_values[num_int_points][fit->num_comps];
  real_t* int_value_ptrs[num_int_points];
  for (int i = 0; i < num_int_points; ++i)
    int_value_ptrs[i] = &int_values[i][0];
  fit->vtable.get_interior_data(fit->context, data, fit->num_comps, int_indices, num_int_points,
                                int_points, int_value_ptrs);

  // Fetch the boundary points, normal vectors, and values.
  fit->vtable.get_boundary_data(fit->context, data, fit->num_comps, bnd_indices, num_bnd_points,
                                bnd_points, bnd_normals, bnd_conditions);

  // Now fit the component data to a polynomial.
  real_t poly_coeffs[dim][fit->num_comps];
  real_t* poly_coeff_ptrs[dim];
  for (int i = 0; i < dim; ++i)
    poly_coeff_ptrs[i] = &poly_coeffs[i][0];
  fit->vtable.fit_data(fit->context, fit->degree, 
                       int_points, int_value_ptrs, num_int_points,
                       bnd_points, bnd_normals, bnd_conditions, num_bnd_points,
                       poly_coeff_ptrs);

  for (int c = 0; c < fit->num_comps; ++c)
  {
    // Construct the polynomial for this component.
    real_t coeffs[dim];
    for (int i = 0; i < dim; ++i)
      coeffs[i] = poly_coeffs[i][c];
    if (fit->polys[c] == NULL)
      fit->polys[c] = polynomial_new(fit->degree, coeffs, &int_points[0]);
    else
    {
      memcpy(polynomial_coeffs(fit->polys[c]), coeffs, sizeof(real_t) * dim);
      *polynomial_x0(fit->polys[c]) = int_points[0];
    }
  }
}

int polynomial_fit_degree(polynomial_fit_t* fit)
{
  return fit->degree;
}

void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value)
{
  for (int i = 0; i < fit->num_comps; ++i)
    value[i] = polynomial_value(fit->polys[i], x);
}

void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv)
{
  for (int i = 0; i < fit->num_comps; ++i)
    deriv[i] = polynomial_deriv_value(fit->polys[i], x_deriv, y_deriv, z_deriv, x);
}

//------------------------------------------------------------------------
//                Freebie polynomial fit constructors
//------------------------------------------------------------------------

// Cell-centered fit.
typedef struct 
{
  mesh_t* mesh;
  int degree;

  // Mapping of face indices to boundary conditions.
  int_ptr_unordered_map_t* face_bc_map;

  // Machinery for fit_data().
  polynomial_fit_fit_data_func fit_data;
  void* fit_data_context;
  polynomial_fit_dtor dtor;
} cc_fit_t;

static int cc_num_interior_neighbors(void* context, int point_index)
{
  cc_fit_t* fit = context;
  int cell = point_index;

  // The number of neighbors is equal to the number of cells attached to 
  // this one via faces.
  int num_neighbors = 0, pos = 0, face;
  while (mesh_cell_next_face(fit->mesh, cell, &pos, &face))
  {
    if (mesh_face_opp_cell(fit->mesh, face, cell) != -1)
      ++num_neighbors;
  }

  return num_neighbors;
}

static void cc_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  cc_fit_t* fit = context;
  int cell = point_index;

  int offset = 0, pos = 0, face;
  while (mesh_cell_next_face(fit->mesh, cell, &pos, &face))
  {
    int opp_cell = mesh_face_opp_cell(fit->mesh, face, cell);
    if (opp_cell != -1)
    {
      neighbor_indices[offset] = opp_cell;
      ++offset;
    }
  }
}

static int cc_num_boundary_neighbors(void* context, int point_index)
{
  cc_fit_t* fit = context;
  int cell = point_index;

  // The number of neighbors is equal to the number of faces without 
  // opposite faces.
  int num_neighbors = 0, pos = 0, face;
  while (mesh_cell_next_face(fit->mesh, cell, &pos, &face))
  {
    if (mesh_face_opp_cell(fit->mesh, face, cell) == -1)
      ++num_neighbors;
  }
  return num_neighbors;
}

static void cc_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  cc_fit_t* fit = context;
  int cell = point_index;

  int offset = 0, pos = 0, face;
  while (mesh_cell_next_face(fit->mesh, cell, &pos, &face))
  {
    if (mesh_face_opp_cell(fit->mesh, face, cell) == -1)
    {
      neighbor_indices[offset] = face;
      ++offset;
    }
  }
}

static void cc_get_interior_data(void* context, real_t* data, int num_comps,
                                 int* point_indices, int num_points,
                                 point_t* point_coords, real_t** point_values)
{
  cc_fit_t* fit = context;

  // We assume the data is in component-minor form in the array, 
  // and that values are at cell centers.
  for (int i = 0; i < num_points; ++i)
  {
    int j = point_indices[i];
    point_coords[i] = fit->mesh->cell_centers[j];
    for (int c = 0; c < num_comps; ++c)
      point_values[i][c] = data[num_comps*j+c];
  }
}

static void cc_get_boundary_data(void* context, real_t* data, int num_comps,
                                 int* point_indices, int num_points,
                                 point_t* point_coords, vector_t* boundary_normals, void** boundary_conditions)
{
  cc_fit_t* fit = context;

  // We assume the data is in component-minor form in the array, 
  // and that boundary values are at face centers.
  // FIXME: For higher-order quadrature, probably have to be more clever.
  for (int i = 0; i < num_points; ++i)
  {
    int j = point_indices[i];
    point_coords[i] = fit->mesh->face_centers[j];
    boundary_normals[i] = fit->mesh->face_normals[j];

    void** bc_ptr = int_ptr_unordered_map_get(fit->face_bc_map, j);
    ASSERT(bc_ptr != NULL);
    boundary_conditions[i] = *bc_ptr;
  }
}

static int cc_targeted_degree(void* context, int num_points)
{
  cc_fit_t* fit = context;
  return fit->degree;
}

static void cc_fit_data(void* context, int degree,
                        point_t* interior_points, real_t** interior_values, int num_interior_points,
                        point_t* boundary_points, vector_t* boundary_normals, void** boundary_conditions, int num_boundary_points,
                        real_t** poly_coeffs)
{
  cc_fit_t* fit = context;
  fit->fit_data(fit->fit_data_context, degree, interior_points, interior_values, num_interior_points,
                boundary_points, boundary_normals, boundary_conditions, num_boundary_points, poly_coeffs);
}

static void cc_dtor(void* context)
{
  cc_fit_t* fit = context;
  if ((fit->fit_data_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_data_context);
  int_ptr_unordered_map_free(fit->face_bc_map);
}

polynomial_fit_t* cc_polynomial_fit_new(int num_comps,
                                        mesh_t* mesh,
                                        int degree,
                                        boundary_cell_map_t* boundary_cells,
                                        polynomial_fit_fit_data_func fit_data,
                                        void* fit_data_context,
                                        polynomial_fit_dtor dtor)
{
  ASSERT(mesh != NULL);
  ASSERT(degree >= 0);

  char name[1024];
  snprintf(name, 1024, "cell-centered (degree %d)", degree);
  cc_fit_t* context = malloc(sizeof(cc_fit_t));
  context->mesh = mesh;
  context->degree = degree;

  // Create a mapping from face indices to boundary conditions.
  context->face_bc_map = int_ptr_unordered_map_new();
  {
    int pos = 0, cell;
    boundary_cell_t* bcell;
    while (boundary_cell_map_next(boundary_cells, &pos, &cell, &bcell))
    {
      for (int i = 0; i < bcell->num_boundary_faces; ++i)
      {
        int f = bcell->boundary_faces[i];
        void* bc = bcell->bc_for_face[i];
        int_ptr_unordered_map_insert(context->face_bc_map, f, bc);
      }
    }
  }

  context->fit_data = fit_data;
  context->fit_data_context = fit_data_context;
  context->dtor = dtor;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cc_num_interior_neighbors,
                                  .get_interior_neighbors = cc_get_interior_neighbors,
                                  .get_interior_data = cc_get_interior_data,
                                  .num_boundary_neighbors = cc_num_boundary_neighbors,
                                  .get_boundary_neighbors = cc_get_boundary_neighbors,
                                  .get_boundary_data = cc_get_boundary_data,
                                  .targeted_degree = cc_targeted_degree,
                                  .fit_data = cc_fit_data,
                                  .dtor = cc_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

// Point cloud fit.
typedef struct 
{
  point_cloud_t* points;
  int degree;

  // Mapping of boundary points to boundary conditions.
  int_ptr_unordered_map_t* bc_for_point;

  // Machinery for fit_data().
  polynomial_fit_fit_data_func fit_data;
  void* fit_data_context;
  polynomial_fit_dtor dtor;
} cloud_fit_t;

static int cloud_num_interior_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static int cloud_num_boundary_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static void cloud_get_interior_data(void* context, real_t* data, int num_comps,
                                    int* point_indices, int num_points,
                                    point_t* point_coords, real_t** point_values)
{
  cloud_fit_t* fit = context;

  // We assume the data is in component-minor form in the array, 
  // and that values are at cell centers.
  for (int i = 0; i < num_points; ++i)
  {
    int j = point_indices[i];
    point_coords[i] = fit->points->point_coords[j];
    for (int c = 0; c < num_comps; ++c)
      point_values[i][c] = data[num_comps*j+c];
  }
}

static void cloud_get_boundary_data(void* context, real_t* data, int num_comps,
                                    int* point_indices, int num_points,
                                    point_t* point_coords, vector_t* boundary_normals, void** boundary_conditions)
{
  // FIXME
}

static int cloud_fit_targeted_degree(void* context, int num_points)
{
  cloud_fit_t* fit = context;
  return fit->degree;
}

static void cloud_fit_data(void* context, int degree,
                           point_t* interior_points, real_t** interior_values, int num_interior_points,
                           point_t* boundary_points, vector_t* boundary_normals, void** boundary_conditions, int num_boundary_points,
                           real_t** poly_coeffs)
{
  cloud_fit_t* fit = context;
  fit->fit_data(fit->fit_data_context, degree, interior_points, interior_values, num_interior_points,
                boundary_points, boundary_normals, boundary_conditions, num_boundary_points, poly_coeffs);
}

static void cloud_dtor(void* context)
{
  cloud_fit_t* fit = context;
  if ((fit->fit_data_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_data_context);
}

polynomial_fit_t* point_cloud_polynomial_fit_new(int num_comps,
                                                 point_cloud_t* points,
                                                 point_cloud_neighbor_search_t* search,
                                                 int degree,
                                                 polynomial_fit_fit_data_func fit_data,
                                                 void* fit_data_context,
                                                 polynomial_fit_dtor dtor)
{
  ASSERT(points != NULL);
  ASSERT(search != NULL);
  ASSERT(degree >= 0);

  char name[1024];
  snprintf(name, 1024, "point cloud fixed degree (%d)", degree);
  cloud_fit_t* context = malloc(sizeof(cloud_fit_t));
  context->points = points;
  context->degree = degree;
  context->fit_data = fit_data;
  context->fit_data_context = fit_data_context;
  context->dtor = dtor;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cloud_num_interior_neighbors,
                                  .get_interior_neighbors = cloud_get_interior_neighbors,
                                  .get_interior_data = cloud_get_interior_data,
                                  .num_boundary_neighbors = cloud_num_boundary_neighbors,
                                  .get_boundary_neighbors = cloud_get_boundary_neighbors,
                                  .get_boundary_data = cloud_get_boundary_data,
                                  .targeted_degree = cloud_fit_targeted_degree,
                                  .fit_data = cloud_fit_data,
                                  .dtor = cloud_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

