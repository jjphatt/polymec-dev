// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LINEAR_SOLVER_H
#define POLYMEC_LINEAR_SOLVER_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class linear_solver
/// This class represents a type that solves a linear system Ax = b.
///
/// Polymec's linear_solver class is actually a thin wrapper around the Sundials 
/// SUNLinearSolver type. We use this wrapper to separate concerns.
typedef struct linear_solver_t linear_solver_t;

// Here's the Sundials linear solver type.
typedef struct _generic_SUNLinearSolver* SUNLinearSolver;

/// Creates an instance of a linear solver from the given Sundials solver. This
/// solver assumes ownership of the Sundials one.
/// \param [in] sundials_solver A Sundials linear solver from which to create 
///                             this one.
/// \memberof linear_solver
linear_solver_t* linear_solver_from_SUNLinearSolver(SUNLinearSolver sundials_solver);

/// Frees a linear solver.
/// \memberof linear_solver
void linear_solver_free(linear_solver_t* solver);

/// Returns the underlying Sundials linear solver.
/// \memberof linear_solver
SUNLinearSolver linear_solver_as_SUNLinearSolver(linear_solver_t* solver);

///@}

#endif

