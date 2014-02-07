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

#include "core/linear_algebra.h"
#include "core/sundials_helpers.h"
#include "integrators/block_jacobi_preconditioner.h"

typedef struct 
{
  adj_graph_t* sparsity;
  adj_graph_coloring_t* coloring;
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  void (*communicate)(void* context, real_t t, real_t* x);
  void* context;

  int num_block_rows;
  int block_size; 

  // Work vectors.
  int num_work_vectors;
  real_t** work;

} block_jacobi_preconditioner_t;

// Block-diagonal matrix.
typedef struct
{
  int num_block_rows;
  int block_size;
  real_t* coeffs;
} bd_mat_t;

static void bd_scale_and_shift(void* context, real_t gamma)
{
  bd_mat_t* mat = context;
  for (int i = 0; i < mat->num_block_rows*mat->block_size; ++i)
    mat->coeffs[i] = 1.0 + gamma * mat->coeffs[i];
  int n = mat->num_block_rows;
  int bs = mat->block_size;
  for (int i = 0; i < n; ++i)
  {
    real_t* A = &mat->coeffs[i*bs*bs];

    // Scale.
    for (int j = 0; j < bs*bs; ++j)
      A[j] *= gamma;

    // Shift.
    for (int j = 0; j < bs; ++j)
      A[j*bs+j] += 1.0;
  }
}

static real_t bd_coeff(void* context, int i, int j)
{
  bd_mat_t* mat = context;
  int bs = mat->block_size;
  if (abs(j - i) >= bs)
    return 0.0;
  int block_row = i/bs;
  int r = i % bs;
  int c = j - block_row*bs;
  if ((r < 0) || (r >= bs) || (c < 0) || (c >= bs))
    return 0.0;
  real_t* A = &mat->coeffs[block_row * bs * bs];
  return A[bs * c + r];
}

static void bd_fprintf(void* context, FILE* stream)
{
  bd_mat_t* mat = context;
  int n = mat->num_block_rows;
  int bs = mat->block_size;
  fprintf(stream, "Block diagonal matrix: (%d rows):", n * bs);
  for (int i = 0; i < n; ++i)
  {
    fprintf(stream, "\nRows %6d - %6d: ", i*bs, (i+1)*bs - 1);
    matrix_fprintf(&mat->coeffs[i * bs * bs], bs, bs, stream);
  }
  fprintf(stream, "\n");
}

static void bd_dtor(void* context)
{
  bd_mat_t* mat = context;
  free(mat->coeffs);
  free(mat);
}

static preconditioner_matrix_t* block_jacobi_preconditioner_matrix(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  preconditioner_matrix_vtable vtable = {.scale_and_shift = bd_scale_and_shift,
                                         .coeff = bd_coeff,
                                         .fprintf = bd_fprintf,
                                         .dtor = bd_dtor};
  int bs = precond->block_size;
  int num_rows = precond->num_block_rows * bs;
  bd_mat_t* mat = malloc(sizeof(bd_mat_t));
  mat->num_block_rows = precond->num_block_rows;
  mat->block_size = bs;
  mat->coeffs = malloc(sizeof(real_t) * mat->num_block_rows * bs * bs);
  memset(mat->coeffs, 0, sizeof(real_t) * mat->num_block_rows * bs * bs);
  return preconditioner_matrix_new("Block-diagonal", mat, vtable, num_rows);
}

// Here's our finite difference implementation of the Jacobian matrix-vector 
// product. 
static void finite_diff_Jv(int (*F)(void* context, real_t t, real_t* x, real_t* F), 
                           void (*communicate)(void* context, real_t t, real_t* x),
                           void* context, 
                           real_t* x, 
                           real_t t, 
                           int num_rows,
                           real_t* v, 
                           real_t** work, 
                           real_t* Jv)
{
  real_t eps = rsqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(x).
  // work[2] == u + eps*v
  // work[3] == F(x + eps*v)

  // u + eps*v -> work[2].
  for (int i = 0; i < num_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v) -> work[3].
  if (communicate != NULL)
    communicate(context, t, work[2]);
  F(context, t, work[2], work[3]);

  // (F(x + eps*v) - F(x)) / eps -> Jv
  for (int i = 0; i < num_rows; ++i)
    Jv[i] = (work[3][i] - work[1][i]) / eps;
}

static void insert_Jv_into_bd_mat(adj_graph_t* graph,
                                  adj_graph_coloring_t* coloring, 
                                  int color, 
                                  real_t* Jv, 
                                  bd_mat_t* J)
{
  int bs = J->block_size;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    int block_row = i / bs;
    for (int j = block_row*bs; j < (block_row+1)*bs; ++j)
    {
      int r = i % bs;
      int c = j % bs;
      J->coeffs[block_row*bs*bs + r*bs + c] = Jv[j];
    }
  }
}

static void block_jacobi_preconditioner_compute_jacobian(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat)
{
  block_jacobi_preconditioner_t* precond = context;
  adj_graph_t* graph = precond->sparsity;
  adj_graph_coloring_t* coloring = precond->coloring;
  real_t** work = precond->work;

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int num_rows = adj_graph_num_vertices(graph);
  ASSERT(num_rows == precond->block_size*precond->num_block_rows);
  real_t* Jv = malloc(sizeof(real_t) * num_rows);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * num_rows);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(x) and place it into work[1].
    if (precond->communicate != NULL)
      precond->communicate(precond->context, t, x);
    precond->F(precond->context, t, x, work[1]);

    // Now evaluate the matrix-vector product.
    memset(Jv, 0, sizeof(real_t) * num_rows);
    finite_diff_Jv(precond->F, precond->communicate, precond->context, x, t, num_rows, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    bd_mat_t* J = preconditioner_matrix_context(mat);
    insert_Jv_into_bd_mat(graph, coloring, c, Jv, J);
  }
  free(Jv);
}


static bool block_jacobi_preconditioner_solve(void* context, preconditioner_matrix_t* A, real_t* B)
{
  block_jacobi_preconditioner_t* precond = context;
  int bs = precond->block_size;
  bd_mat_t* mat = preconditioner_matrix_context(A);

  for (int i = 0; i < precond->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &mat->coeffs[i*bs*bs], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[i*bs], sizeof(real_t)*bs);

    // Solve the linear system.
    int one = 1, ipiv[bs], info;
    dgesv(&bs, &one, Aij, &bs, ipiv, bi, &bs, &info);
    ASSERT(info == 0);

    // Copy the solution into place.
    memcpy(&B[i*bs], bi, sizeof(real_t)*bs);
  }

  return true;
}

static void block_jacobi_preconditioner_dtor(void* context)
{
  block_jacobi_preconditioner_t* precond = context;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    free(precond->work[i]);
  free(precond->work);
  adj_graph_coloring_free(precond->coloring);
  adj_graph_free(precond->sparsity);
  free(precond);
}

preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  int (*residual_func)(void* context, real_t t, real_t* x, real_t* F),
                                                  void (*communication_func)(void* context, real_t t, real_t* x),
                                                  adj_graph_t* sparsity,
                                                  int num_block_rows,
                                                  int block_size)
{
  ASSERT(num_block_rows > 0);
  ASSERT(block_size > 0);

  block_jacobi_preconditioner_t* precond = malloc(sizeof(block_jacobi_preconditioner_t));

  // Do we have a block graph?
  int num_rows = adj_graph_num_vertices(sparsity);
  ASSERT((num_rows == num_block_rows) || (num_rows = block_size*num_block_rows));
  if (num_rows == num_block_rows)
    precond->sparsity = adj_graph_new_with_block_size(block_size, sparsity);
  else
    precond->sparsity = adj_graph_clone(sparsity);

  precond->coloring = adj_graph_coloring_new(precond->sparsity, SMALLEST_LAST);
  precond->F = residual_func;
  precond->communicate = communication_func;
  precond->context = context;
  precond->num_block_rows = num_block_rows;
  precond->block_size = block_size;

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = malloc(sizeof(real_t*) * precond->num_work_vectors);
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = malloc(sizeof(real_t) * precond->num_block_rows * precond->block_size);

  preconditioner_vtable vtable = {.matrix = block_jacobi_preconditioner_matrix,
                                  .compute_jacobian = block_jacobi_preconditioner_compute_jacobian,
                                  .solve = block_jacobi_preconditioner_solve,
                                  .dtor = block_jacobi_preconditioner_dtor};
  return preconditioner_new("Block Jacobi preconditioner", precond, vtable);
}
