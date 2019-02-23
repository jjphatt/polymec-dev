// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NONLINEAR_SOLVER_H
#define POLYMEC_NONLINEAR_SOLVER_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class nonlinear_solver
/// This class represents a type that solves a nonlinear system F(X) = 0.
///
/// Polymec's nonlinear_solver class is actually a thin wrapper around the 
/// Sundials SUNNonLinearSolver type. We use this wrapper to separate concerns.
typedef struct nonlinear_solver_t nonlinear_solver_t;

// Here's the Sundials nonlinear solver type.
typedef struct _generic_SUNNonlinearSolver* SUNNonlinearSolver;

/// Creates an instance of a nonlinear solver from the given Sundials solver. 
/// This solver assumes ownership of the Sundials one.
/// \param [in] sundials_solver A Sundials nonlinear solver from which to create 
///                             this one.
/// \memberof nonlinear_solver
nonlinear_solver_t* nonlinear_solver_from_SUNNonlinearSolver(SUNNonlinearSolver sundials_solver);

/// Frees a nonlinear solver.
/// \memberof nonlinear_solver
void nonlinear_solver_free(nonlinear_solver_t* solver);

/// Returns the underlying Sundials nonlinear solver.
/// \memberof nonlinear_solver
SUNNonlinearSolver nonlinear_solver_as_SUNNonlinearSolver(nonlinear_solver_t* solver);

///@}

#endif

