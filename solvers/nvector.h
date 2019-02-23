// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NVECTOR_H
#define POLYMEC_NVECTOR_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class nvector
/// An n-vector represents an n-dimensional vector that appears in linear 
/// equations like the x and b in Ax = b.
///
/// N-vectors have several implementations corresponding to different linear 
/// solvers, so they're typically created with factory methods that belong to 
/// those solvers.
///
/// Polymec's n-vector class is actually a thin wrapper around the Sundials 
/// N_Vector type. We use this wrapper to separate concerns.
typedef struct nvector_t nvector_t;

// Here's the Sundials N_Vector type.
typedef struct _generic_N_Vector* N_Vector;

/// Creates an instance of an n-vector from the given Sundials N_Vector. The
/// n-vector assumes ownership of the N_Vector.
/// \param [in] sundials_nvector A Sundials N_Vector from which to create this 
///                              n-vector.
/// \memberof nvector
nvector_t* nvector_from_N_Vector(N_Vector sundials_nvector);

/// Frees an n-vector.
/// \memberof nvector
void nvector_free(nvector_t* v);

/// Returns the underlying Sundials N_Vector.
/// \memberof nvector
N_Vector nvector_as_N_Vector(nvector_t* v);

///@}

#endif

