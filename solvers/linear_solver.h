// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LINEAR_SOLVER_H
#define POLYMEC_LINEAR_SOLVER_H

#include "solvers/preconditioner.h"

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

typedef struct _generic_SUNMatrix* SUNMatrix;

/// This virtual table defines the behavior for a linear solver. Because the 
/// solver class is implemented using the SUNLinearSolver type, all methods use 
/// SUNLinearSolver objects instead of raw context pointers. All methods
/// with int return values return 0 on success and nonzero on failure.
typedef struct
{
  /// Sets left and right scaling vectors s1 and s2 for the solve.
  int (*setscalingvectors)(SUNLinearSolver, N_Vector s1, N_Vector s2);
  /// Initializes the solver.
  int (*initialize)(SUNLinearSolver);
  /// Performs pre-solve setup. Matrix argument can be ignored for matrix-free solvers.
  int (*setup)(SUNLinearSolver, SUNMatrix A);
  /// Performs the linear solve, terminating when the residual norm falls below the tolerance.
  int (*solve)(SUNLinearSolver, SUNMatrix A, N_Vector x, N_Vector b, real_t tolerance);
  /// Estimates the space needed for the solver in real words (lrw) and integer words (liw). (Optional)
  int (*space)(SUNLinearSolver, long int* lrw, long int* liw);
  /// Destroys the solver.
  int (*free)(SUNLinearSolver);
} linear_solver_vtable;

/// \enum linear_solver_type_t;
/// This enumerated type describes types of linear solvers.
typedef enum
{
  DIRECT_LINEAR_SOLVER,
  ITERATIVE_LINEAR_SOLVER,
  MATRIX_FREE_LINEAR_SOLVER
} linear_solver_type_t;

/// Creates an instance of a direct linear solver on the given communicator. 
/// \param [in] comm The MPI communicator on which the solver lives.
/// \param [in] context A pointer to data for the solver.
/// \param [in] vtable A virtual table defining the behavior of the solver.
/// \memberof linear_solver
linear_solver_t* direct_linear_solver_new(MPI_Comm comm, 
                                          void* context, 
                                          linear_solver_vtable vtable);

/// Creates an instance of an iterative linear solver on the given communicator. 
/// This solver requires a matrix, but solves the linear system iteratively.
/// \param [in] comm The MPI communicator on which the solver lives.
/// \param [in] context A pointer to data for the solver.
/// \param [in] vtable A virtual table defining the behavior of the solver.
/// \param [in] p The preconditioner to use for the iterative solution.
/// \param [in] max_iterations The maximum number of iterations for the solve.
/// \param [in] tolerance The tolerance for the residual norm.
/// \memberof linear_solver
linear_solver_t* iterative_linear_solver_new(MPI_Comm comm, 
                                             void* context, 
                                             linear_solver_vtable vtable,
                                             preconditioner_t* p,
                                             int max_iterations,
                                             real_t tolerance);

/// Creates an instance of an iterative linear solver on the given communicator. 
/// This solver does not require a matrix. Instead, one provides a function that 
/// computes the product of a linear operator with a vector.
/// \param [in] comm The MPI communicator on which the solver lives.
/// \param [in] context A pointer to data for the solver.
/// \param [in] vtable A virtual table defining the behavior of the solver.
/// \param [in] A_times A function that computes the product of a linear operator 
///                     with a vector.
/// \param [in] p The preconditioner to use for the iterative solution.
/// \param [in] max_iterations The maximum number of iterations for the solve.
/// \param [in] tolerance The tolerance for the residual norm.
/// \memberof linear_solver
linear_solver_t* matrix_free_linear_solver_new(MPI_Comm comm, 
                                               void* context, 
                                               linear_solver_vtable vtable,
                                               int (*A_times)(void* context, N_Vector v, N_Vector Av),
                                               preconditioner_t* p,
                                               int max_iterations,
                                               real_t tolerance);

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

/// Returns the type of this linear solver.
linear_solver_type_t linear_solver_type(linear_solver_t* solver);

///@}

#endif

