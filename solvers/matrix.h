// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MATRIX_H
#define POLYMEC_MATRIX_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class matrix
/// This class represents a matrix that appears in linear equations, like A
/// in Ax = b.
///
/// Matrices have several implementations corresponding to different linear 
/// solvers, so they're typically created with factory methods that belong to 
/// those solvers.
///
/// Polymec's matrix class is actually a thin wrapper around the Sundials 
/// SUNMatrix type. We use this wrapper to separate concerns.
typedef struct matrix_t matrix_t;

// Here's the Sundials SUNMatrix type.
typedef struct _generic_SUNMatrix* SUNMatrix;

/// Creates an instance of a matrix from the given Sundials SUNMatrix. The
/// matrix assumes ownership of the SUNMatrix.
/// \param [in] sundials_matrix A Sundials SUNMatrix from which to create this 
///                             matrix.
/// \memberof matrix
matrix_t* matrix_from_SUNMatrix(SUNMatrix sundials_matrix);

/// Frees a matrix.
/// \memberof matrix
void matrix_free(matrix_t* A);

/// Returns the underlying Sundials SUNMatrix.
/// \memberof matrix
SUNMatrix matrix_as_SUNMatrix(matrix_t* A);

///@}

#endif

