// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "model/gmls_matrix.h"
#include "make_mlpg_lattice.h"
#include "poisson_gmls_functional.h"

void test_gmls_matrix_ctor(void** state)
{
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  make_mlpg_lattice(10, 10, 10, 2.0, &points, &extents, &stencil);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(W, points, extents, stencil);

  // Clean up.
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

static void franke(void* context, point_t* x, real_t* u)
{
  real_t X = x->x, Y = x->y;
  u[0] = 0.75 * exp(-0.25 * (pow(9.0*X - 2.0, 2.0) + pow(9.0*Y - 2.0, 2.0))) + 
         0.75 * exp(-1.0/49.0 * pow(9.0*X + 1.0, 2.0) - 1.0/10.0 * pow(9.0*Y + 1.0, 2.0)) + 
         0.50 * exp(-0.25 * (pow(9.0*X - 7.0, 2.0) + pow(9.0*Y - 3.0, 2.0))) - 
         0.20 * exp(-pow(9.0*X - 4.0, 2.0) - pow(9.0*Y - 7.0, 2.0));
}

void test_gmls_matrix_with_frankes_function(void** state)
{
  point_cloud_t* points;
  real_t* extents;
  stencil_t* stencil;
  make_mlpg_lattice(10, 10, 1, 2.0, &points, &extents, &stencil);
  point_weight_function_t* W = gaussian_point_weight_function_new(4.0);
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(W, points, extents, stencil);
  gmls_functional_t* lambda = poisson_gmls_functional_new(2, points, extents);
  sp_func_t* F = sp_func_from_func("Franke's function", franke, SP_INHOMOGENEOUS, 1);

  // Set up our beloved linear system using a dense matrix.
  // FIXME: point lattice should have boundary points labeled.
  real_t A[100*100], B[100];

  // Treat boundary nodes first.
  int_unordered_set_t* boundary_nodes = int_unordered_set_new();
  int num_bnodes; 
  int* bnodes = point_cloud_tag(points, "boundary", &num_bnodes);
  ASSERT(bnodes != NULL);
  for (int b = 0; b < num_bnodes; ++b)
  {
    int bnode = bnodes[b];
    int num_cols = gmls_matrix_num_columns(matrix, bnode);
    int cols[num_cols];
    real_t coeffs[num_cols];
    gmls_matrix_compute_dirichlet_row(matrix, bnode, lambda, cols, coeffs);
    for (int c = 0; c < num_cols; ++c)
    {
      // FIXME
    }
    int_unordered_set_insert(boundary_nodes, bnode);
  }

  for (int r = 0; r < 100; ++r)
  {
    if (!int_unordered_set_contains(boundary_nodes, r))
    {
      int num_cols = gmls_matrix_num_columns(matrix, r);
      int cols[num_cols];
      real_t coeffs[num_cols];
      gmls_matrix_compute_row(matrix, r, lambda, 0.0, cols, coeffs);
      for (int c = 0; c < num_cols; ++c)
      {
        // FIXME
      }
    }
  }

  // Fill in the RHS vector.
  for (int b = 0; b < num_bnodes; ++b)
  {
    int bnode = bnodes[b];
    // FIXME
  }

  // Solve the linear system.

  // Clean up.
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  polymec_free(extents);
  stencil_free(stencil);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmls_matrix_ctor)
  };
  return run_tests(tests);
}
