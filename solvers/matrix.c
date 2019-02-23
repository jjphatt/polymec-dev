// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_matrix.h"
#include "solvers/matrix.h"

struct matrix_t 
{
  SUNMatrix A;
};

matrix_t* matrix_from_SUNMatrix(SUNMatrix sundials_matrix)
{
  matrix_t* A = polymec_malloc(sizeof(matrix_t));
  A->A = sundials_matrix;
  return A;
}

void matrix_free(matrix_t* A)
{
  SUNMatDestroy(A->A);
  polymec_free(A);
}

SUNMatrix matrix_as_SUNMatrix(matrix_t* A)
{
  return A->A;
}

