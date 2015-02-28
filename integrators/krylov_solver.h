// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_SOLVER_H
#define POLYMEC_KRYLOV_SOLVER_H

#include "core/polymec.h"

// Types of Krylov solver.
typedef enum
{
  KRYLOV_GMRES,
  KRYLOV_BICGSTAB,
  KRYLOV_TFQMR,
} krylov_t;

// Objects of this type solve linear systems A*x = b using preconditioned 
// Krylov subspace methods (built from matrix-vector products A*x, A*Ax, etc).
typedef struct krylov_solver_t krylov_solver_t;

// Creates a Krylov linear solver with a given maximum subspace dimension of 
// max_krylov_dim. matrix_vector_product is a function that computes the product 
// A*y for the linear operator A and the vector y, storing the result in the 
// array Ay. If the solver_type is KRYLOV_GMRES, the maximum number of restarts 
// is given by max_restarts--otherwise that parameter is ignored. The system has 
// num_local_values local equations, and num_remote_values remote ones.
krylov_solver_t* krylov_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*matrix_vector_product)(void* context, real_t* y, real_t* Ay),
                                   void (*dtor)(void* context),
                                   krylov_t solver_type,
                                   int max_krylov_dim,
                                   int max_restarts);

// Frees a solver.
void krylov_solver_free(krylov_solver_t* solver);

// Returns the context pointer for the solver.
void* krylov_solver_context(krylov_solver_t* solver);

// Returns the number of (local) equations in the linear system.
int krylov_solver_num_equations(krylov_solver_t* solver);

// Sets the tolerances for the residual norm (norm_tolerance) = |R|, where 
// R = A*x - b, and the stop tolerance (fractional decrease in the residual
// norm).
void krylov_solver_set_tolerances(krylov_solver_t* solver, 
                                  real_t norm_tolerance, 
                                  real_t stop_tolerance);

// Sets the maximum number of iterations for the solver.
void krylov_solver_set_max_iterations(krylov_solver_t* solver, int max_iterations);

// Solves the linear system of equations A * x = b in place, storing the solution
// in the array b. Returns true if the solution was obtained, false if not. The 
// number of linear iterations will be stored in num_iterations upon success.
bool krylov_solver_solve(krylov_solver_t* solver, real_t* b, int* num_iterations);

// Diagnostics for the linear solver.
typedef struct
{
  char* status_message; // borrowed pointer from solver: do not free.
  long int num_function_evaluations;
  real_t function_norm;
  long int num_solve_iterations;
  long int num_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
  long int num_matrix_vector_product_evaluations;
} krylov_solver_diagnostics_t;

// Retrieve diagnostics for the nonlinear solver.
void krylov_solver_get_diagnostics(krylov_solver_t* solver, 
                                   krylov_solver_diagnostics_t* diagnostics);

// Writes nonlinear solver diagnostics to the given file.
void krylov_solver_diagnostics_fprintf(krylov_solver_diagnostics_t* diagnostics, 
                                       FILE* stream);

#endif
