// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRECONDITIONER_H
#define POLYMEC_PRECONDITIONER_H

#include "solvers/nvector.h"

/// \addtogroup solvers solvers
///@{

/// \class preconditioner
/// This class represents a preconditioner for an iterative linear solver.
/// Preconditioners are strongly coupled to their related linear solvers.
typedef struct preconditioner_t preconditioner_t;

/// \enum preconditioner_side_t
/// Preconditioner sided-ness.
typedef enum
{
  PRECONDITIONER_LEFT,
  PRECONDITIONER_RIGHT,
  PRECONDITIONER_BOTH
} preconditioner_side_t;

/// This virtual table defines the behavior for the preconditioner.
typedef struct 
{
  /// Set up the preconditioner system Pz = r for solving.
  /// Returns 0 on success and nonzero on failure.
  int (*set_up)(void* context);

  /// Solve the system Pz = r for the vector z, given its data and sidedness. 
  /// If the preconditioner is iterative, the solution z should satisfy
  /// || Pz - r || < tolerance for the norm relevant to the linear solver 
  /// (often a weighted RMS norm).
  /// Returns 0 on success, a negative number for an unrecoverable condition,
  /// and a positive number for a recoverable condition.
  int (*solve)(void* context, preconditioner_side_t sidedness,
               nvector_t* r, real_t tolerance, nvector_t* z);

  /// This is a destructor function for the context pointer.
  void (*dtor)(void* context);
} preconditioner_vtable;

/// Creates an instance of a preconditioner, defining its name, behavior, 
/// and sided-ness.
/// \memberof preconditioner
preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable,
                                     preconditioner_side_t sidedness);

/// Frees a preconditioner.
/// \memberof preconditioner
void preconditioner_free(preconditioner_t* pc);

/// Returns the "sided-ness" of the preconditioner.
preconditioner_side_t preconditioner_sidedness(preconditioner_t* pc);

///@}

#endif

