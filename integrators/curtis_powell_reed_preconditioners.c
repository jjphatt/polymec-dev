// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/linear_algebra.h"
#include "core/array_utils.h"
#include "core/sundials_helpers.h"
#include "integrators/curtis_powell_reed_preconditioners.h"
#include "integrators/newton_preconditioner.h"

typedef struct
{
  int (*F)(void* context, real_t t, real_t* x, real_t* Fval);
  int (*F_dae)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;
  void (*set_identity_matrix)(void* context, real_t alpha);
  void (*add_Jv_into_matrix)(void* context, adj_graph_t* graph,
                             adj_graph_coloring_t* coloring, 
                             int color, 
                             real_t factor,
                             real_t* Jv);
  bool (*solve)(void* context, real_t* z);
  void (*fprintf)(void* context, FILE* stream);
  void (*dtor)(void* context);
} cpr_vtable;

typedef struct
{
  MPI_Comm comm;
  void* context;
  cpr_vtable vtable;

  // Graph coloring stuff.
  adj_graph_t* sparsity;
  adj_graph_coloring_t* coloring;
  int max_colors; // Maximum number of colors on all processes.

  // Numbers of local and remote rows.
  int num_local_rows, num_remote_rows;

  // Work vectors.
  int num_work_vectors;
  real_t** work;
} cpr_t;

// This function adapts non-DAE functions F(t, x) to DAE ones F(t, x, xdot).
static int F_adaptor(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval)
{
  ASSERT(xdot == NULL);

  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  cpr_t* pc = context;
  return pc->vtable.F(pc->vtable.F_context, t, x, Fval);
}

// Here's our finite difference implementation of the dF/dx matrix-vector 
// product. 
static void finite_diff_dFdx_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                               void* context, 
                               real_t t, 
                               real_t* x, 
                               real_t* xdot, 
                               int num_local_rows,
                               int num_remote_rows,
                               real_t* v, 
                               real_t** work, 
                               real_t* dFdx_v)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == x + eps*v
  // work[3] == F(t, x + eps*v)

  // x + eps*v -> work[2].
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = x[i] + eps*v[i];

  // F(t, x + eps*v, xdot) -> work[3].
  F(context, t, work[2], xdot, work[3]);

  // (F(t, x + eps*v, xdot) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdx_v[i] = (work[3][i] - work[1][i]) / eps;
}

// Here's the same finite difference calculation for dF/d(xdot).
static void finite_diff_dFdxdot_v(int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F), 
                                  void* context, 
                                  real_t t, 
                                  real_t* x, 
                                  real_t* xdot, 
                                  int num_local_rows,
                                  int num_remote_rows,
                                  real_t* v, 
                                  real_t** work, 
                                  real_t* dFdxdot_v)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(t, x, xdot).
  // work[2] == xdot + eps*v
  // work[3] == F(t, x, xdot + eps*v)

  // xdot + eps*v -> work[2].
  for (int i = 0; i < num_local_rows + num_remote_rows; ++i)
    work[2][i] = xdot[i] + eps*v[i];

  // F(t, x, xdot + eps*v) -> work[3].
  F(context, t, x, work[2], work[3]);

  // (F(t, x, xdot + eps*v) - F(t, x, xdot)) / eps -> (dF/dx) * v
  for (int i = 0; i < num_local_rows; ++i)
    dFdxdot_v[i] = (work[3][i] - work[1][i]) / eps;
}

static void cpr_compute_P(void* context, 
                          real_t alpha, real_t beta, real_t gamma,
                          real_t t, real_t* x, real_t* xdot)
{
  cpr_t* precond = context;
  adj_graph_t* graph = precond->sparsity;
  adj_graph_coloring_t* coloring = precond->coloring;
  real_t** work = precond->work;

  // Normalize F.
  int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval);
  void* F_context;
  if (precond->vtable.F_dae != NULL)
  {
    ASSERT(xdot != NULL);
    F = precond->vtable.F_dae;
    F_context = precond->vtable.F_context;
  }
  else
  {
    ASSERT(xdot == NULL);
    F = F_adaptor;
    F_context = precond;
  }

  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  log_debug("curtis_powell_reed_preconditioner: approximating Jacobian...");

  // First, set up an identity matrix.
  precond->vtable.set_identity_matrix(precond->context, alpha);

  // Now iterate over all of the colors in our coloring. 
  real_t* Jv = polymec_malloc(sizeof(real_t) * precond->num_local_rows);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(real_t) * (precond->num_local_rows + precond->num_remote_rows));
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      work[0][i] = 1.0;

    // We evaluate F(t, x, xdot) and place it into work[1].
    F(F_context, t, x, xdot, work[1]);

    // Evaluate dF/dx * d and stash it in P.
    memset(Jv, 0, sizeof(real_t) * precond->num_local_rows);
    finite_diff_dFdx_v(F, F_context, t, x, xdot, precond->num_local_rows, 
                       precond->num_remote_rows, work[0], work, Jv);
    precond->vtable.add_Jv_into_matrix(precond->context, graph, coloring, c, beta, Jv);

    if ((gamma != 0.0) && (xdot != NULL))
    {
      // Now evaluate dF/d(xdot) * d and do the same.
      memset(Jv, 0, sizeof(real_t) * precond->num_local_rows);
      finite_diff_dFdxdot_v(F, F_context, t, x, xdot, precond->num_local_rows, 
                            precond->num_remote_rows, work[0], work, Jv);
      precond->vtable.add_Jv_into_matrix(precond->context, graph, coloring, c, gamma, Jv);
    }
  }

  // Now call the RHS functions in the same way as we would for all the colors
  // we don't have, up through the maximum number, so our neighboring 
  // processes can get exchanged data from us if they need it.
  int num_neighbor_colors = precond->max_colors - num_colors;
  for (int c = 0; c < num_neighbor_colors; ++c)
  {
    F(F_context, t, x, xdot, work[1]);
    finite_diff_dFdx_v(F, F_context, t, x, xdot, precond->num_local_rows, 
                       precond->num_remote_rows, work[0], work, Jv);
    if ((gamma != 0.0) && (xdot != NULL))
    {
      finite_diff_dFdxdot_v(F, F_context, t, x, xdot, precond->num_local_rows, 
                            precond->num_remote_rows, work[0], work, Jv);
    }
  }

  polymec_free(Jv);
}

static bool cpr_solve(void* context, real_t* z)
{
  cpr_t* precond = context;
  return precond->vtable.solve(precond->context, z);
}

static void cpr_fprintf(void* context, FILE* stream)
{
  cpr_t* precond = context;
  precond->vtable.fprintf(precond->context, stream);
}

static void cpr_dtor(void* context)
{
  cpr_t* precond = context;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    polymec_free(precond->work[i]);
  polymec_free(precond->work);
  adj_graph_coloring_free(precond->coloring);
  adj_graph_free(precond->sparsity);
  if ((precond->vtable.dtor != NULL) && (precond->context != NULL))
    precond->vtable.dtor(precond->context);
  polymec_free(precond);
}

static preconditioner_t* curtis_powell_reed_preconditioner_new(const char* name,
                                                               MPI_Comm comm,
                                                               void* context,
                                                               cpr_vtable vtable,
                                                               adj_graph_t* sparsity,
                                                               int num_local_block_rows,
                                                               int num_remote_block_rows,
                                                               int block_size)
{
  ASSERT(num_remote_block_rows >= 0);
  ASSERT(vtable.set_identity_matrix != NULL);
  ASSERT(vtable.add_Jv_into_matrix != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(block_size >= 1);

  // Exactly one of F and F_dae must be given.
  ASSERT((vtable.F != NULL) || (vtable.F_dae != NULL));
  ASSERT((vtable.F == NULL) || (vtable.F_dae == NULL));

  cpr_t* precond = polymec_malloc(sizeof(cpr_t));
  precond->comm = comm;
  precond->context = context;
  precond->vtable = vtable;

  // Do we have a block graph?
  int num_local_rows = adj_graph_num_vertices(sparsity);
  ASSERT((num_local_rows == num_local_block_rows) || (num_local_rows == block_size*num_local_block_rows));
  if (num_local_rows == num_local_block_rows)
  {
    // We were given the number of vertices in the graph as the number of 
    // block rows, so we create a graph with a block size of 1.
    precond->sparsity = adj_graph_new_with_block_size(block_size, sparsity);
  }
  else
  {
    // The number of vertices in the graph is the number of degrees of freedom
    // in the solution, so we don't need to create
    precond->sparsity = adj_graph_clone(sparsity);
  }

  // In either case, we multiply by the block size to get the number of 
  // degrees of freedom.
  precond->num_local_rows = num_local_block_rows * block_size;
  precond->num_remote_rows = num_remote_block_rows * block_size;

  // Assemble a graph coloring for the preconditioner matrix.
  precond->coloring = adj_graph_coloring_new(precond->sparsity, SMALLEST_LAST);

  // Get the maximum number of colors on all MPI processes so that we can 
  // still compute preconditioners in lockstep.
  int num_colors = adj_graph_coloring_num_colors(precond->coloring);
  MPI_Allreduce(&num_colors, &precond->max_colors, 1, MPI_INT, MPI_MAX, precond->comm);
  log_debug("curtis_powell_reed_preconditioner: graph coloring produced %d colors.", 
            precond->max_colors);

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = polymec_malloc(sizeof(real_t*) * precond->num_work_vectors);
  int N = (num_local_block_rows + num_remote_block_rows) * block_size;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = polymec_malloc(sizeof(real_t) * N);

  newton_preconditioner_vtable nlpc_vtable = {.compute_P = cpr_compute_P,
                                              .solve = cpr_solve,
                                              .fprintf = cpr_fprintf,
                                              .dtor = cpr_dtor};
  return newton_preconditioner_new(name, precond, nlpc_vtable);
}

// Variable block-size version.
static preconditioner_t* vbs_curtis_powell_reed_preconditioner_new(const char* name,
                                                                   MPI_Comm comm,
                                                                   void* context,
                                                                   cpr_vtable vtable,
                                                                   adj_graph_t* sparsity,
                                                                   int num_local_block_rows,
                                                                   int num_remote_block_rows,
                                                                   int* block_sizes)
{
  ASSERT(num_remote_block_rows >= 0);
  ASSERT(vtable.set_identity_matrix != NULL);
  ASSERT(vtable.add_Jv_into_matrix != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(block_sizes != NULL);

  // Exactly one of F and F_dae must be given.
  ASSERT((vtable.F != NULL) || (vtable.F_dae != NULL));
  ASSERT((vtable.F == NULL) || (vtable.F_dae == NULL));

  cpr_t* precond = polymec_malloc(sizeof(cpr_t));
  precond->comm = comm;
  precond->context = context;
  precond->vtable = vtable;

  // Do we have a block graph?
  int num_local_rows = adj_graph_num_vertices(sparsity);
  int alleged_num_local_rows = 0;
  int max_block_size = 1;
  for (int r = 0; r < num_local_block_rows; ++r)
  {
    ASSERT(block_sizes[r] >= 1);
    max_block_size = MAX(block_sizes[r], max_block_size);
    alleged_num_local_rows += block_sizes[r];
  }
  ASSERT((num_local_rows == num_local_block_rows) || (num_local_rows == alleged_num_local_rows)); 
  if (num_local_rows == num_local_block_rows)
  {
    // We were given the number of vertices in the graph as the number of 
    // block rows, so we create a graph with a block size of 1.
    precond->sparsity = adj_graph_new_with_block_sizes(block_sizes, sparsity);
    ASSERT(adj_graph_num_vertices(precond->sparsity) == alleged_num_local_rows);
  }
  else
  {
    // The number of vertices in the graph is the number of degrees of freedom
    // in the solution, so we don't need to create
    precond->sparsity = adj_graph_clone(sparsity);
  }

  // The number of local rows is known at this point.
  precond->num_local_rows = alleged_num_local_rows;

  // However, we can't know the actual number of remote block rows without 
  // knowing the specific communication pattern! However, it is safe to just 
  // multiply the given number of remote block rows by the maximum block size, 
  // since only the underlying RHS function will actually use the data.
  precond->num_remote_rows = num_remote_block_rows * max_block_size;

  // Assemble a graph coloring for the preconditioner matrix.
  precond->coloring = adj_graph_coloring_new(precond->sparsity, SMALLEST_LAST);

  // Get the maximum number of colors on all MPI processes so that we can 
  // still compute preconditioners in lockstep.
  int num_colors = adj_graph_coloring_num_colors(precond->coloring);
  MPI_Allreduce(&num_colors, &precond->max_colors, 1, MPI_INT, MPI_MAX, precond->comm);
  log_debug("curtis_powell_reed_preconditioner: graph coloring produced %d colors.", 
            precond->max_colors);

  // Make work vectors.
  precond->num_work_vectors = 4;
  precond->work = polymec_malloc(sizeof(real_t*) * precond->num_work_vectors);
  int N = precond->num_local_rows + precond->num_remote_rows;
  for (int i = 0; i < precond->num_work_vectors; ++i)
    precond->work[i] = polymec_malloc(sizeof(real_t) * N);

  newton_preconditioner_vtable nlpc_vtable = {.compute_P = cpr_compute_P,
                                              .solve = cpr_solve,
                                              .fprintf = cpr_fprintf,
                                              .dtor = cpr_dtor};
  return newton_preconditioner_new(name, precond, nlpc_vtable);
}

//------------------------------------------------------------------------
//              Block-Jacobi Newton preconditioner.
//------------------------------------------------------------------------

typedef struct 
{
  void* context;
  void (*dtor)(void* context);

  // Preconditioner matrix data -- block diagonal.
  int num_block_rows, block_size;
  int *D_offsets, *B_offsets; // For variable block sizes.
  real_t* D;
} bjpc_t;

static bjpc_t* bjpc_new(void* context, 
                        void (*dtor)(void* context),
                        int num_block_rows,
                        int block_size, 
                        int* block_sizes)
{
  bjpc_t* pc = polymec_malloc(sizeof(bjpc_t));
  pc->context = context;
  pc->dtor = dtor;
  pc->num_block_rows = num_block_rows;
  pc->block_size = block_size;
  pc->D_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  pc->B_offsets = polymec_malloc(sizeof(int) * (num_block_rows+1));
  pc->D_offsets[0] = pc->B_offsets[0] = 0;
  for (int i = 0; i < num_block_rows; ++i)
  {
    int bs = (block_size == -1) ? block_sizes[i] : block_size;
    ASSERT(bs >= 1);
    pc->D_offsets[i+1] = pc->D_offsets[i] + bs*bs;
    pc->B_offsets[i+1] = pc->B_offsets[i] + bs;
  }
  int N = pc->D_offsets[pc->num_block_rows];
  pc->D = polymec_malloc(sizeof(real_t) * N);
  return pc;
}

static void bjpc_set_identity_matrix(void* context, real_t diag_val)
{
  bjpc_t* pc = context;

  // Zero the matrix coefficients.
  memset(pc->D, 0, sizeof(real_t) * pc->D_offsets[pc->num_block_rows]);

  // Set the diagonal values.
  for (int i = 0; i < pc->num_block_rows; ++i)
  {
    int bs = pc->B_offsets[i+1] - pc->B_offsets[i];
    int offset = pc->D_offsets[i];
    for (int j = 0; j < bs; ++j)
      pc->D[offset+bs*j+j] = diag_val;
  }
}

static void bjpc_add_Jv_into_matrix(void* context,
                                    adj_graph_t* graph,
                                    adj_graph_coloring_t* coloring, 
                                    int color, 
                                    real_t factor,
                                    real_t* Jv)
{
  bjpc_t* pc = context;
  real_t* D = pc->D;

  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    if (i >= pc->B_offsets[pc->num_block_rows]) continue;

    int bs = (pc->block_size == -1) ? (pc->B_offsets[i+1] - pc->B_offsets[i]) : pc->block_size;
    int block_col = i / bs;
    int c = i % bs;
    for (int j = block_col*bs; j < (block_col+1)*bs; ++j)
    {
      int r = j % bs;
      D[pc->D_offsets[block_col] + c*bs + r] += factor * Jv[j];
    }
  }
}

static bool bjpc_solve(void* context, real_t* B)
{
  bjpc_t* pc = context;
  real_t* D = pc->D;

  bool success = false;
  for (int i = 0; i < pc->num_block_rows; ++i)
  {
    // Copy the block for this row into place.
    int bs = pc->B_offsets[i+1] - pc->B_offsets[i];
    int D_offset = pc->D_offsets[i], B_offset = pc->B_offsets[i];
    real_t Aij[bs*bs], bi[bs];
    memcpy(Aij, &D[D_offset], sizeof(real_t)*bs*bs);
    memcpy(bi, &B[B_offset], sizeof(real_t)*bs);

    // Replace each zero on the diagonal of Aij with a small number.
    static const real_t epsilon = 1e-25;
    for (int j = 0; j < bs; ++j)
    {
      if (Aij[bs*j+j] == 0.0)
        Aij[bs*j+j] = epsilon;
    }

    // Solve the linear system.
    int one = 1, ipiv[bs], info;
    dgesv(&bs, &one, Aij, &bs, ipiv, bi, &bs, &info);
    success = (info == 0);

    if (success)
    {
      // Copy the solution into place.
      memcpy(&B[B_offset], bi, sizeof(real_t)*bs);
    }
    else
    {
      ASSERT(info > 0);
      log_debug("block_jacobi_preconditioner_solve: call to dgesv failed for block row %d.", i);
      log_debug("(U is singular).", i);
      break;
    }

    D_offset += bs*bs;
    B_offset += bs;
  }

  return success;
}

static void bjpc_fprintf(void* context, FILE* stream)
{
  bjpc_t* pc = context;
  int N = pc->num_block_rows;
  real_t* D = pc->D;
  if (pc->block_size != -1)
    fprintf(stream, "\nBlock diagonal matrix P (block size = %d):\n", pc->block_size);
  else
    fprintf(stream, "\nBlock diagonal matrix P (variable block size):\n");
  for (int i = 0; i < N; ++i)
  {
    fprintf(stream, "%d: [", i);
    int bs = pc->B_offsets[i+1] - pc->B_offsets[i];
    int offset = pc->D_offsets[i];
    for (int ii = 0; ii < bs; ++ii)
    {
      for (int jj = 0; jj < bs; ++jj)
        fprintf(stream, "%g ", D[offset+bs*ii+jj]);
      if (ii < (bs - 1))
        fprintf(stream, "; ");
    }
    fprintf(stream, "]\n");
  }
}

static void bjpc_dtor(void* context)
{
  bjpc_t* precond = context;
  polymec_free(precond->D);
  if ((precond->dtor != NULL) && (precond->context != NULL))
    precond->dtor(precond->context);
  polymec_free(precond);
}

preconditioner_t* block_jacobi_preconditioner_from_function(const char* name, 
                                                            MPI_Comm comm,
                                                            void* context,
                                                            int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                            void (*dtor)(void* context),
                                                            adj_graph_t* sparsity,
                                                            int num_local_block_rows,
                                                            int num_remote_block_rows,
                                                            int block_size)
{
  ASSERT(F != NULL);
  bjpc_t* pc = bjpc_new(context, dtor, num_local_block_rows, block_size, NULL);
  cpr_vtable vtable = {.F = F,
                       .F_context = context,
                       .set_identity_matrix = bjpc_set_identity_matrix,
                       .add_Jv_into_matrix = bjpc_add_Jv_into_matrix,
                       .solve = bjpc_solve,
                       .fprintf = bjpc_fprintf,
                       .dtor = bjpc_dtor};
  return curtis_powell_reed_preconditioner_new(name, comm, pc, vtable, sparsity, 
                                               num_local_block_rows, num_remote_block_rows, 
                                               block_size);
}

preconditioner_t* block_jacobi_preconditioner_from_dae_function(const char* name, 
                                                                MPI_Comm comm,
                                                                void* context,
                                                                int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                                void (*dtor)(void* context),
                                                                adj_graph_t* sparsity,
                                                                int num_local_block_rows,
                                                                int num_remote_block_rows,
                                                                int block_size)
{
  ASSERT(F != NULL);
  bjpc_t* pc = bjpc_new(context, dtor, num_local_block_rows, block_size, NULL);
  cpr_vtable vtable = {.F_dae = F,
                       .F_context = context,
                       .set_identity_matrix = bjpc_set_identity_matrix,
                       .add_Jv_into_matrix = bjpc_add_Jv_into_matrix,
                       .solve = bjpc_solve,
                       .dtor = bjpc_dtor};
  return curtis_powell_reed_preconditioner_new(name, comm, pc, vtable, sparsity, 
                                               num_local_block_rows, num_remote_block_rows,
                                               block_size);
}

//------------------------------------------------------------------------
//          Variable block-size Block-Jacobi Newton preconditioner.
//------------------------------------------------------------------------

preconditioner_t* var_block_jacobi_preconditioner_from_function(const char* name, 
                                                                MPI_Comm comm,
                                                                void* context,
                                                                int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                                void (*dtor)(void* context),
                                                                adj_graph_t* sparsity,
                                                                int num_local_block_rows,
                                                                int num_remote_block_rows,
                                                                int* block_sizes)
{
  ASSERT(F != NULL);
  ASSERT(block_sizes != NULL);
  bjpc_t* pc = bjpc_new(context, dtor, num_local_block_rows, -1, block_sizes);
  cpr_vtable vtable = {.F = F,
                       .F_context = context,
                       .set_identity_matrix = bjpc_set_identity_matrix,
                       .add_Jv_into_matrix = bjpc_add_Jv_into_matrix,
                       .solve = bjpc_solve,
                       .fprintf = bjpc_fprintf,
                       .dtor = bjpc_dtor};
  return vbs_curtis_powell_reed_preconditioner_new(name, comm, pc, vtable, sparsity, 
                                                   num_local_block_rows, num_remote_block_rows, 
                                                   block_sizes);
}

preconditioner_t* var_block_jacobi_preconditioner_from_dae_function(const char* name, 
                                                                    MPI_Comm comm,
                                                                    void* context,
                                                                    int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                                    void (*dtor)(void* context),
                                                                    adj_graph_t* sparsity,
                                                                    int num_local_block_rows,
                                                                    int num_remote_block_rows,
                                                                    int* block_sizes)
{
  ASSERT(F != NULL);
  ASSERT(block_sizes != NULL);
  bjpc_t* pc = bjpc_new(context, dtor, num_local_block_rows, -1, block_sizes);
  cpr_vtable vtable = {.F_dae = F,
                       .F_context = context,
                       .set_identity_matrix = bjpc_set_identity_matrix,
                       .add_Jv_into_matrix = bjpc_add_Jv_into_matrix,
                       .solve = bjpc_solve,
                       .dtor = bjpc_dtor};
  return vbs_curtis_powell_reed_preconditioner_new(name, comm, pc, vtable, sparsity, 
                                                   num_local_block_rows, num_remote_block_rows,
                                                   block_sizes);
}

//------------------------------------------------------------------------
//              LU Newton preconditioner.
//------------------------------------------------------------------------

typedef struct 
{
  void* context;
  void (*dtor)(void* context);

  int N; // Number of rows in matrix.

  // Matrix data.
  int* etree;
  SuperMatrix rhs, X, L, U;
  real_t *rhs_data, *X_data;
  int *rperm, *cperm;
  superlu_options_t options;
  SuperLUStat_t stat;
  SuperMatrix* P; // <-- preconditioner matrix.

  // ILU parameters (if any).
  ilu_params_t* ilu_params;
} lupc_t;

static void supermatrix_free(SuperMatrix* matrix)
{
  switch(matrix->Stype) {
  case SLU_DN:
    Destroy_Dense_Matrix(matrix);
    break;
  case SLU_NR:
    Destroy_CompRow_Matrix(matrix);
    break;
  case SLU_NC:
    Destroy_CompCol_Matrix(matrix);
    break;
  default:
    printf("ERROR: unknown matrix type passed to supermatrix_free!");
    break;
  }
  polymec_free(matrix);
}

static SuperMatrix* supermatrix_new(adj_graph_t* graph)
{
  SuperMatrix* A = polymec_malloc(sizeof(SuperMatrix));
  
  // Fetch sparsity information from the graph.
  int* edge_offsets = adj_graph_edge_offsets(graph);
  int num_rows = adj_graph_num_vertices(graph);
  int num_edges = edge_offsets[num_rows];
  int num_nz = num_rows + num_edges;

  // Translate the sparsity information in the graph to column indices and 
  // row pointers. Since the graph is effectively stored in compressed 
  // row format, we will use SuperLU's compressed row matrix format.
  int* row_ptrs = intMalloc(num_rows + 1);
  int* col_indices = intMalloc(num_nz);
  int* edges = adj_graph_adjacency(graph);
  int offset = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    row_ptrs[i] = offset;
    col_indices[offset++] = i; // diagonal column

    // Off-diagonal columns.
    for (int j = edge_offsets[i]; j < edge_offsets[i+1]; ++j)
      col_indices[offset++] = edges[j];
    int_qsort(&col_indices[row_ptrs[i]+1], edge_offsets[i+1]-edge_offsets[i]);
  }
  row_ptrs[num_rows] = offset;
  ASSERT(offset == num_nz);

  // Create zeros for the matrix.
  real_t* mat_zeros = SUPERLU_MALLOC(sizeof(real_t) * num_nz);
  memset(mat_zeros, 0, sizeof(real_t) * num_nz);

  // Hand over these resources to create the Supermatrix.
  dCreate_CompCol_Matrix(A, num_rows, num_rows, num_nz, 
                         mat_zeros, col_indices, row_ptrs, 
                         SLU_NC, SLU_D, SLU_GE);
  return A;
}

static void lupc_set_identity_matrix(void* context, real_t diag_val)
{
  lupc_t* pc = context;
  SuperMatrix* mat = pc->P;
  int num_cols = mat->ncol;
  NCformat* data = mat->Store;
  real_t* Aij = data->nzval;

  // Zero the matrix coefficients.
  memset(Aij, 0, data->colptr[num_cols]);

  // Set the diagonal values.
  for (int j = 0; j < num_cols; ++j)
  {
    int col_index = data->colptr[j];
    Aij[col_index] = diag_val;
  }
}

static void lupc_add_Jv_into_matrix(void* context, 
                                    adj_graph_t* graph, 
                                    adj_graph_coloring_t* coloring, 
                                    int color, 
                                    real_t factor, 
                                    real_t* Jv)
{
  lupc_t* precond = context;
  SuperMatrix* mat = precond->P;

  NCformat* data = mat->Store;
  real_t* Jij = data->nzval;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    if (i >= precond->N) continue;

    // Add in the diagonal element.
    Jij[data->colptr[i]] += factor * Jv[i];
      
    // Add in off-diagonal column values.
    int pos = 0, j;
    while (adj_graph_next_edge(graph, i, &pos, &j))
    {
      int col_index = data->colptr[i];
      size_t num_rows = data->colptr[i+1] - col_index;
      int* entry = int_bsearch(&data->rowind[col_index+1], num_rows - 1, j);
      ASSERT(entry != NULL);
      size_t offset = entry - &data->rowind[col_index];
      Jij[data->colptr[i] + offset] += factor * Jv[j];
    }
  }
}

static bool lupc_solve(void* context, real_t* B)
{
  lupc_t* precond = context;
  SuperMatrix* mat = precond->P;

  // Copy B to the rhs vector.
  DNformat* rhs = precond->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * precond->N);

  if (precond->cperm == NULL)
  {
    precond->cperm = intMalloc(precond->N);
    precond->rperm = intMalloc(precond->N);
  }
  else if (precond->options.Fact != SamePattern)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }

  // Do the solve.
  int info;
  dgssv(&precond->options, mat, precond->cperm, precond->rperm, 
        &precond->L, &precond->U, &precond->rhs, &precond->stat, &info);
  // SuperLU claims that it can re-use factorization info with the same 
  // sparsity pattern, but this appears not to be the case at the moment.
  // precond->options.Fact = SamePattern;

  bool success = (info == 0);

  if (success)
  {
    // Copy the rhs vector to B.
    memcpy(B, rhs->nzval, sizeof(real_t) * precond->N);
  }

  return success;
}

static void lupc_fprintf(void* context, FILE* stream)
{
  lupc_t* precond = context;
  SuperMatrix* A = precond->P;
  fprintf(stream, "\nCompCol preconditioner matrix P:\n");
  fprintf(stream, "Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
  int n = A->ncol;
  NCformat* Astore = (NCformat *) A->Store;
  real_t* dp = Astore->nzval;
  fprintf(stream, "nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
  fprintf(stream, "nzval: ");
  for (int i = 0; i < Astore->colptr[n]; ++i) 
    fprintf(stream, "%f  ", (double)dp[i]);
  fprintf(stream, "\nrowind: ");
  for (int i = 0; i < Astore->colptr[n]; ++i) 
    fprintf(stream, "%d  ", Astore->rowind[i]);
  fprintf(stream, "\ncolptr: ");
  for (int i = 0; i <= n; ++i) 
    fprintf(stream, "%d  ", Astore->colptr[i]);
  fprintf(stream, "\n");
}

static void lupc_dtor(void* context)
{
  lupc_t* precond = context;
  if (precond->cperm != NULL)
  {
    SUPERLU_FREE(precond->cperm);
    SUPERLU_FREE(precond->rperm);
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }
  Destroy_SuperMatrix_Store(&precond->rhs);
  polymec_free(precond->rhs_data);
  Destroy_SuperMatrix_Store(&precond->X);
  supermatrix_free(precond->P);
  polymec_free(precond->X_data);
  StatFree(&precond->stat);
  if (precond->etree != NULL)
    polymec_free(precond->etree);
  precond->ilu_params = NULL;
  if ((precond->dtor != NULL) && (precond->context != NULL))
    precond->dtor(precond->context);
  polymec_free(precond);
}

static void lupc_initialize_matrix_data(void* context, adj_graph_t* sparsity)
{
  lupc_t* precond = context;
  precond->ilu_params = NULL;
  precond->P = supermatrix_new(sparsity);

  // Preconditioner data.
  precond->N = adj_graph_num_vertices(sparsity);
  precond->rhs_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->rhs, precond->N, 1, precond->rhs_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  precond->X_data = polymec_malloc(sizeof(real_t) * precond->N);
  dCreate_Dense_Matrix(&precond->X, precond->N, 1, precond->X_data, precond->N, SLU_DN, SLU_D, SLU_GE);
  StatInit(&precond->stat);
  precond->cperm = NULL;
  precond->rperm = NULL;
  set_default_options(&precond->options);
  precond->options.ColPerm = NATURAL;
  precond->options.Fact = DOFACT;
  precond->etree = NULL;
}
   
preconditioner_t* lu_preconditioner_from_function(const char* name, 
                                                  MPI_Comm comm,
                                                  void* context,
                                                  int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                  void (*dtor)(void* context),
                                                  adj_graph_t* sparsity,
                                                  int num_local_block_rows, 
                                                  int num_remote_block_rows, 
                                                  int block_size)
{
  return ilu_preconditioner_from_function(name, comm, context, F, dtor, sparsity, 
                                          num_local_block_rows, num_remote_block_rows, 
                                          block_size, NULL);
}

preconditioner_t* lu_preconditioner_from_dae_function(const char* name, 
                                                      MPI_Comm comm, 
                                                      void* context,
                                                      int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                      void (*dtor)(void* context),
                                                      adj_graph_t* sparsity,
                                                      int num_local_block_rows, 
                                                      int num_remote_block_rows, 
                                                      int block_size)
{
  return ilu_preconditioner_from_dae_function(name, comm, context, F, dtor, sparsity, 
                                              num_local_block_rows, num_remote_block_rows, 
                                              block_size, NULL);
}

preconditioner_t* var_lu_preconditioner_from_function(const char* name, 
                                                      MPI_Comm comm,
                                                      void* context,
                                                      int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                      void (*dtor)(void* context),
                                                      adj_graph_t* sparsity,
                                                      int num_local_block_rows, 
                                                      int num_remote_block_rows, 
                                                      int* block_sizes)
{
  return var_ilu_preconditioner_from_function(name, comm, context, F, dtor, sparsity, 
                                              num_local_block_rows, num_remote_block_rows, 
                                              block_sizes, NULL);
}

preconditioner_t* var_lu_preconditioner_from_dae_function(const char* name, 
                                                          MPI_Comm comm, 
                                                          void* context,
                                                          int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                          void (*dtor)(void* context),
                                                          adj_graph_t* sparsity,
                                                          int num_local_block_rows, 
                                                          int num_remote_block_rows, 
                                                          int* block_sizes)
{
  return var_ilu_preconditioner_from_dae_function(name, comm, context, F, dtor, sparsity, 
                                                  num_local_block_rows, num_remote_block_rows, 
                                                  block_sizes, NULL);
}

// Globals borrowed from SuperLU.
const int ILU_DROP_BASIC = DROP_BASIC;
const int ILU_DROP_PROWS = DROP_PROWS;
const int ILU_DROP_COLUMN = DROP_COLUMN;
const int ILU_DROP_AREA = DROP_AREA;
const int ILU_DROP_DYNAMIC = DROP_DYNAMIC;
const int ILU_DROP_INTERP = DROP_INTERP;

ilu_params_t* ilu_params_new()
{
  ilu_params_t* params = GC_MALLOC(sizeof(ilu_params_t));
  params->diag_pivot_threshold = 0.1;
  params->row_perm = ILU_LARGE_DIAG_PERM;
  params->drop_rule = ILU_DROP_BASIC | ILU_DROP_AREA;
  params->drop_tolerance = 1e-4;
  params->fill_factor = 10.0;
  params->milu_variant = ILU_SILU;
  params->fill_tolerance = 0.01;
  params->norm = ILU_LINF;
  return params;
}

static bool ilupc_solve(void* context, real_t* B)
{
  lupc_t* precond = context;
  SuperMatrix* mat = precond->P;

  // Copy B to the rhs vector.
  DNformat* rhs = precond->rhs.Store;
  memcpy(rhs->nzval, B, sizeof(real_t) * precond->N);

  if (precond->cperm == NULL)
  {
    precond->cperm = intMalloc(precond->N);
    precond->rperm = intMalloc(precond->N);
  }
  else if (precond->options.Fact != SamePattern)
  {
    // Jettison the existing factorization.
    Destroy_SuperNode_Matrix(&precond->L);
    Destroy_CompCol_Matrix(&precond->U);
  }

  // Do the (approximate) solve.
  int info;
  char equed;
  int lwork = 0; // Indicates SuperLU internal memory allocation.
  real_t R[precond->N], C[precond->N], recip_pivot_growth, cond_number;
  mem_usage_t mem_usage;
  dgsisx(&precond->options, mat, precond->cperm, precond->rperm, 
         precond->etree, &equed, R, C, &precond->L, &precond->U, 
         NULL, lwork, &precond->rhs, &precond->X, &recip_pivot_growth, &cond_number,
         &mem_usage, &precond->stat, &info);

  bool success = (info == 0);

  // Copy the rhs vector to B.
  memcpy(B, rhs->nzval, sizeof(real_t) * precond->N);

  return success;
}

static void ilupc_set_ilu_params(lupc_t* precond, ilu_params_t* ilu_params)
{
  if (ilu_params != NULL)
  {
    // Copy options into place.
    ilu_set_default_options(&precond->options);
    precond->ilu_params = ilu_params;
    precond->options.ILU_DropRule = ilu_params->drop_rule;
    precond->options.ILU_DropTol = ilu_params->drop_tolerance;
    precond->options.ILU_FillFactor = ilu_params->fill_factor;
    precond->options.ILU_FillTol = ilu_params->fill_tolerance;
    if (ilu_params->norm == ILU_L1)
      precond->options.ILU_Norm = ONE_NORM;
    else if (ilu_params->norm == ILU_L2)
      precond->options.ILU_Norm = TWO_NORM;
    else
      precond->options.ILU_Norm = INF_NORM;
    if (ilu_params->milu_variant == ILU_SILU)
      precond->options.ILU_MILU = SILU;
    else if (ilu_params->milu_variant == ILU_MILU1)
      precond->options.ILU_MILU = SMILU_1;
    else if (ilu_params->milu_variant == ILU_MILU2)
      precond->options.ILU_MILU = SMILU_2;
    else 
      precond->options.ILU_MILU = SMILU_3;
    precond->options.ILU_MILU_Dim = 3;
    if (ilu_params->row_perm == ILU_NO_ROW_PERM)
      precond->options.RowPerm = NOROWPERM;
    else
      precond->options.RowPerm = LargeDiag;
  }
  else
    precond->ilu_params = NULL;
}

preconditioner_t* ilu_preconditioner_from_function(const char* name, 
                                                   MPI_Comm comm,
                                                   void* context,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity,
                                                   int num_local_block_rows, 
                                                   int num_remote_block_rows, 
                                                   int block_size,
                                                   ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  lupc_t* precond = polymec_malloc(sizeof(lupc_t));
  precond->context = context;
  precond->dtor = dtor;
  lupc_initialize_matrix_data(precond, sparsity);
  ilupc_set_ilu_params(precond, ilu_params);
  cpr_vtable vtable = {.F = F,
                       .F_context = context,
                       .set_identity_matrix = lupc_set_identity_matrix,
                       .add_Jv_into_matrix = lupc_add_Jv_into_matrix,
                       .solve = lupc_solve,
                       .fprintf = lupc_fprintf,
                       .dtor = lupc_dtor};
  if (ilu_params != NULL)
    vtable.solve = ilupc_solve;
  return curtis_powell_reed_preconditioner_new(name, comm, precond, vtable, sparsity, 
                                               num_local_block_rows, num_remote_block_rows,
                                               block_size);
}

preconditioner_t* ilu_preconditioner_from_dae_function(const char* name, 
                                                       MPI_Comm comm,
                                                       void* context,
                                                       int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                       void (*dtor)(void* context),
                                                       adj_graph_t* sparsity,
                                                       int num_local_block_rows, 
                                                       int num_remote_block_rows, 
                                                       int block_size,
                                                       ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  lupc_t* precond = polymec_malloc(sizeof(lupc_t));
  precond->context = context;
  precond->dtor = dtor;
  lupc_initialize_matrix_data(precond, sparsity);
  ilupc_set_ilu_params(precond, ilu_params);
  cpr_vtable vtable = {.F_dae = F,
                       .F_context = context,
                       .set_identity_matrix = lupc_set_identity_matrix,
                       .add_Jv_into_matrix = lupc_add_Jv_into_matrix,
                       .solve = lupc_solve,
                       .fprintf = lupc_fprintf,
                       .dtor = lupc_dtor};
  if (ilu_params != NULL)
    vtable.solve = ilupc_solve;
  return curtis_powell_reed_preconditioner_new(name, comm, precond, vtable, sparsity, 
                                               num_local_block_rows, num_remote_block_rows,
                                               block_size);
}

preconditioner_t* var_ilu_preconditioner_from_function(const char* name, 
                                                       MPI_Comm comm,
                                                       void* context,
                                                       int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                       void (*dtor)(void* context),
                                                       adj_graph_t* sparsity,
                                                       int num_local_block_rows, 
                                                       int num_remote_block_rows, 
                                                       int* block_sizes,
                                                       ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  lupc_t* precond = polymec_malloc(sizeof(lupc_t));
  precond->context = context;
  precond->dtor = dtor;
  lupc_initialize_matrix_data(precond, sparsity);
  ilupc_set_ilu_params(precond, ilu_params);
  cpr_vtable vtable = {.F = F,
                       .F_context = context,
                       .set_identity_matrix = lupc_set_identity_matrix,
                       .add_Jv_into_matrix = lupc_add_Jv_into_matrix,
                       .solve = lupc_solve,
                       .fprintf = lupc_fprintf,
                       .dtor = lupc_dtor};
  if (ilu_params != NULL)
    vtable.solve = ilupc_solve;
  return vbs_curtis_powell_reed_preconditioner_new(name, comm, precond, vtable, sparsity, 
                                                   num_local_block_rows, num_remote_block_rows,
                                                   block_sizes);
}

preconditioner_t* var_ilu_preconditioner_from_dae_function(const char* name, 
                                                           MPI_Comm comm,
                                                           void* context,
                                                           int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                           void (*dtor)(void* context),
                                                           adj_graph_t* sparsity,
                                                           int num_local_block_rows, 
                                                           int num_remote_block_rows, 
                                                           int* block_sizes,
                                                           ilu_params_t* ilu_params)
{
  ASSERT(F != NULL);
  lupc_t* precond = polymec_malloc(sizeof(lupc_t));
  precond->context = context;
  precond->dtor = dtor;
  lupc_initialize_matrix_data(precond, sparsity);
  ilupc_set_ilu_params(precond, ilu_params);
  cpr_vtable vtable = {.F_dae = F,
                       .F_context = context,
                       .set_identity_matrix = lupc_set_identity_matrix,
                       .add_Jv_into_matrix = lupc_add_Jv_into_matrix,
                       .solve = lupc_solve,
                       .fprintf = lupc_fprintf,
                       .dtor = lupc_dtor};
  if (ilu_params != NULL)
    vtable.solve = ilupc_solve;
  return vbs_curtis_powell_reed_preconditioner_new(name, comm, precond, vtable, sparsity, 
                                                   num_local_block_rows, num_remote_block_rows,
                                                   block_sizes);
}

