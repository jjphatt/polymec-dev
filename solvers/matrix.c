// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_matrix.h"

#include "core/timer.h"
#include "solvers/matrix.h"
#include "solvers/nvector.h"

struct matrix_t 
{
  MPI_Comm comm;
  int num_local_rows;
  index_t num_global_rows;
  SUNMatrix A;
  matrix_vtable vtable;
  bool owns;
};

static SUNMatrix_ID getid(SUNMatrix A)
{
  return SUNMATRIX_CUSTOM;
}

matrix_t* matrix_new(MPI_Comm comm, 
                     void* context, 
                     matrix_vtable vtable,
                     int num_local_rows,
                     index_t num_global_rows)
{
  ASSERT(num_local_rows > 0);
  ASSERT(num_global_rows >= (index_t)num_local_rows);
  ASSERT(vtable.clone != NULL);
  ASSERT(vtable.destroy != NULL);
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.copy != NULL);
  ASSERT(vtable.scaleadd != NULL);
  ASSERT(vtable.scaleaddi != NULL);
  ASSERT(vtable.matvec != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);

  matrix_t* A = polymec_malloc(sizeof(matrix_t));
  A->comm = comm;
  A->num_local_rows = num_local_rows;
  A->num_global_rows = num_global_rows;
  A->vtable = vtable;
  A->owns = true;

  // Set up the SUNMatrix.
  A->A = polymec_malloc(sizeof(struct _generic_SUNMatrix));
  A->A->content = context;
  A->A->ops = polymec_malloc(sizeof(struct _generic_SUNMatrix_Ops));
  A->A->ops->getid = getid;
  A->A->ops->clone = vtable.clone;
  A->A->ops->destroy = vtable.destroy;
  A->A->ops->zero = vtable.zero;
  A->A->ops->copy = vtable.copy;
  A->A->ops->scaleadd = vtable.scaleadd;
  A->A->ops->scaleaddi = vtable.scaleaddi;
  A->A->ops->matvec = vtable.matvec;
  A->A->ops->space = vtable.space;

  return A;
}

matrix_t* matrix_from_SUNMatrix(SUNMatrix A, bool assume_ownership)
{
  ASSERT(SUNMatGetID(A) != SUNMATRIX_CUSTOM); 

  matrix_t* mat = polymec_malloc(sizeof(matrix_t));
  mat->comm = MPI_COMM_SELF;
  mat->owns = true;
  mat->A = A;

  // Set up the vtable.
  matrix_vtable vtable = {
    .clone = A->ops->clone,
    .destroy = A->ops->destroy,
    .zero = A->ops->zero,
    .copy = A->ops->copy,
    .scaleadd = A->ops->scaleadd,
    .scaleaddi = A->ops->scaleaddi,
    .matvec = A->ops->matvec,
    .space = A->ops->space
  };
  mat->vtable = vtable;

  // FIXME: Fill out rest of vtable based on ID.
  mat->num_local_rows = -1;
  mat->num_global_rows = 0;

  return mat;
}

void matrix_free(matrix_t* A)
{
  if (A->owns)
  {
    SUNMatrix_ID id = SUNMatGetID(A->A);
    if (id == SUNMATRIX_CUSTOM)
      polymec_free(A->A->ops);
    SUNMatDestroy(A->A);
  }
  polymec_free(A);
}

SUNMatrix matrix_as_SUNMatrix(matrix_t* A)
{
  return A->A;
}

MPI_Comm matrix_comm(matrix_t* A)
{
  return A->comm;
}

int matrix_get_block_sizes(matrix_t* A, 
                           int num_block_rows, 
                           index_t* block_rows, 
                           int* block_sizes)
{
  if (A->vtable.get_block_sizes != NULL)
    return A->vtable.get_block_sizes(A->A, num_block_rows, block_rows, block_sizes);
  else
  {
    for (int r = 0; r < num_block_rows; ++r)
      block_sizes[r] = 1;
    return 0;
  }
}

matrix_t* matrix_clone(matrix_t* A)
{
  matrix_t* clone = polymec_malloc(sizeof(matrix_t));
  clone->comm = A->comm;
  clone->num_local_rows = A->num_local_rows;
  clone->num_global_rows = A->num_global_rows;
  clone->vtable = A->vtable;
  clone->A = A->vtable.clone(A->A);
  return clone;
}

int matrix_copy(matrix_t* A, matrix_t* copy)
{
  return A->vtable.copy(A->A, copy->A);
}

int matrix_num_local_rows(matrix_t* A)
{
  return A->num_local_rows;
}

index_t matrix_num_global_rows(matrix_t* A)
{
  return A->num_global_rows;
}

int matrix_zero(matrix_t* A)
{
  return A->vtable.zero(A->A);
}

int matrix_scale_add(matrix_t* A, real_t c, matrix_t* B)
{
  return A->vtable.scaleadd(c, A->A, B->A);
}

int matrix_scale_add_I(matrix_t* A, real_t c)
{
  return A->vtable.scaleaddi(c, A->A);
}

int matrix_mat_vec(matrix_t* A, nvector_t* x, nvector_t* y)
{
  return A->vtable.matvec(A->A, nvector_as_N_Vector(x), nvector_as_N_Vector(y));
}

int matrix_set_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values)
{
  START_FUNCTION_TIMER();
  int result = A->vtable.set_values(A->A, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
  return result;
}
                              
int matrix_add_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values)
{
  START_FUNCTION_TIMER();
  int result = A->vtable.add_values(A->A, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
  return result;
}
                              
int matrix_get_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values)
{
  START_FUNCTION_TIMER();
  int result = A->vtable.get_values(A->A, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
  return result;
}

int matrix_set_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values)
{
  START_FUNCTION_TIMER();
  int result;
  if (A->vtable.set_blocks != NULL)
    result = A->vtable.set_blocks(A->A, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
  return result;
}
                              
int matrix_add_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values)
{
  START_FUNCTION_TIMER();
  int result;
  if (A->vtable.add_blocks != NULL)
    result = A->vtable.add_blocks(A->A, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
  return result;
}
                              
int matrix_get_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values)
{
  START_FUNCTION_TIMER();
  int result;
  if (A->vtable.get_blocks != NULL)
    result = A->vtable.get_blocks(A->A, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
  return result;
}

int matrix_assemble(matrix_t* A)
{
  START_FUNCTION_TIMER();
  int result;
  if (A->vtable.assemble != NULL)
    result = A->vtable.assemble(A->A);
  STOP_FUNCTION_TIMER();
  return result;
}

void matrix_fprintf(matrix_t* A, FILE* stream)
{
  if ((A->vtable.fprintf != NULL) && (stream != NULL))
    A->vtable.fprintf(A->A, stream);
}

