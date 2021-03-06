// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_AM_ODE_SOLVER_H
#define POLYMEC_AM_ODE_SOLVER_H

#include "solvers/ode_solver.h"
#include "solvers/newton_pc.h"

// This type of ODE solver integrates a non-stiff set of ordinary 
// differential equations using the Adams-Moulton formulae. They can provide 
// accurate time integration at an order from 1 to 12. There are two kinds of 
// Adams-Moulton ODE solvers: those based on functional iteration, which 
// requires only function evaluations, and those based on Newton iteration, 
// which requires the representation of a Jacobian matrix. Both of these 
// types of solvers are implemented using CVODE from the Sundials suite 
// of nonlinear solvers.

/// \addtogroup solvers solvers
///@{

/// The functional iteration variant of the Adams-Moulton solver is particularly 
/// simple and can be constructed on the communicator comm for a number of 
/// non-stiff ordinary differential equations (with given numbers of local and 
/// remote values, for parallel capability), using a context pointer, a 
/// right-hand side function, and a destructor.
/// \relates ode_solver
ode_solver_t* functional_am_ode_solver_new(int order, 
                                           MPI_Comm comm, 
                                           int num_local_values,
                                           int num_remote_values,
                                           int num_accel_vectors,
                                           void* context, 
                                           int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                           void (*dtor)(void* context));

/// \enum jfnk_am_krylov_t
/// Describes different types of Jacobian-Free Newton-Krylov iteration 
/// available for the Adams-Moulton ode_solvers.
typedef enum 
{
  JFNK_AM_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_AM_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_AM_TFQMR     // Transpose-Free QMR Krylov solver
} jfnk_am_krylov_t;

/// Creates an Adams-Moulton solver that uses a Jacobian-Free
/// Newton-Krylov method to solve the underlying linearized equations. This 
/// method requires a preconditioner that is a coarse approximation of the 
/// Jacobian matrix, but captures its essential behavior. The right-hand side 
/// function and destructor are required, and a function to compute Jy, the 
/// product of the Jacobian with the vector y, can be optionally provided. Its 
/// signature is: 
/// Jy_func(context, t, x, rhs, y, temp, Jy), where context is the context pointer, 
/// t is the time at which Jy is evaluated, y is the vector the Jacobian operator
/// J is applied to, rhs is the right hand side function evaluated at t and y,
/// temp is a work vector the same size as y, and Jy is a vector that stores the 
/// product Jy.
/// If Jy_func is not given, a finite difference approximation of Jy will be used.
/// Additionally, the type of Krylov solver (JFNK_AM_GMRES, JFNK_AM_BICGSTAB, 
/// or JFNK_AM_TFQMR) must be given, along with the maximum dimension of the 
/// Krylov subspace. 
/// \relates ode_solver
ode_solver_t* jfnk_am_ode_solver_new(int order, 
                                     MPI_Comm comm,
                                     int num_local_equations, 
                                     int num_remote_values, 
                                     void* context, 
                                     int (*rhs_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                     int (*Jy_func)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
                                     void (*dtor)(void* context),
                                     newton_pc_t* precond,
                                     jfnk_am_krylov_t solver_type,
                                     int max_krylov_dim);

/// This returns the context pointer passed to the am_ode_solver 
/// constructor. In general, this will NOT return the same pointer as 
/// ode_solver_context!
/// \relates ode_solver
void* am_ode_solver_context(ode_solver_t* solver);

/// Sets the relative and absolute tolerances for integrated quantities.
/// \relates ode_solver
void am_ode_solver_set_tolerances(ode_solver_t* solver,
                                  real_t relative_tol, real_t absolute_tol);

/// Sets the error weights for evaluating the WRMS norm that is used 
/// as a proxy for the quality of the solution. This may be used in lieu of 
/// relative and absolute tolerances. Weights are copied into the solver.
/// \relates ode_solver
void am_ode_solver_set_error_weights(ode_solver_t* solver, real_t* weights);

/// Sets the error weight function for evaluating the WRMS norm that is used 
/// as a proxy for the quality of the solution. This may be used in lieu of 
/// relative and absolute tolerances.
/// \relates ode_solver
void am_ode_solver_set_error_weight_function(ode_solver_t* solver,
                                             void (*compute_weights)(void* context, real_t* y, real_t* weights));

/// Sets the maximum number of error test failures permitted in attempting 
/// a single time step. By default, this value is 7.
/// \relates ode_solver
void am_ode_solver_set_max_err_test_failures(ode_solver_t* solver,
                                             int max_failures);

/// Sets the maximum number of nonlinear solver iterations per time step.
/// By default, this value is 3.
/// \relates ode_solver
void am_ode_solver_set_max_nonlinear_iterations(ode_solver_t* solver,
                                                int max_iterations);

/// Sets the safety factor (coefficient) used in the nonlinear convergence test.
/// By default, this value is 0.1.
/// \relates ode_solver
void am_ode_solver_set_nonlinear_convergence_coeff(ode_solver_t* solver,
                                                   real_t coefficient);

/// Evaluates the right-hand side of the system at the given time and with the 
/// given solution X, placing the results in rhs.
/// \relates ode_solver
void am_ode_solver_eval_rhs(ode_solver_t* integ, real_t t, real_t* X, real_t* rhs);

/// Returns an internal pointer to the preconditioner passed to this 
/// solver during construction time. (This returns NULL for the functional 
/// solver, which doesn't need a preconditioner.)
/// \relates ode_solver
newton_pc_t* am_ode_solver_preconditioner(ode_solver_t* solver);

/// \class am_ode_solver_diagnostics
/// Diagnostics for the time solver.
typedef struct
{
  char* status_message; // borrowed pointer from solver: do not free.
  long int num_steps;
  int order_of_last_step, order_of_next_step;
  real_t last_step_size, next_step_size;
  long int num_rhs_evaluations;
  long int num_linear_solve_setups;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_error_test_failures;
  long int num_nonlinear_solve_iterations;
  long int num_nonlinear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
} am_ode_solver_diagnostics_t;

/// Retrieve diagnostics for the time solver.
/// \relates ode_solver
void am_ode_solver_get_diagnostics(ode_solver_t* solver, 
                                   am_ode_solver_diagnostics_t* diagnostics);

/// Writes out time solver diagnostics to the given file.
/// \memberof am_ode_solver_diagnostics
void am_ode_solver_diagnostics_fprintf(am_ode_solver_diagnostics_t* diagnostics, 
                                       FILE* stream);

/// \class am_ode_observer
/// This observer type can be used to define objects that respond to actions
/// taken by the am_ode_solver.
typedef struct am_ode_observer_t am_ode_observer_t;

/// Creates and returns a newly-allocated observer that observes an 
/// am_ode_solver. Responses are given by the arguments.
/// \param [in] rhs_computed This function is called when the right hand side of the ODE
///             system is computed by the solver. Accepts the context.
/// \param [in] Jy_computed This function is called when the Jacobian-vector product J*y 
///             is computed by the solver. Accepts the context.
/// \memberof am_ode_observer
am_ode_observer_t* am_ode_observer_new(void* context,
                                       void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                       void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
                                       void (*dtor)(void* context));

/// Adds the given observer to the given am_ode_solver. The observer 
/// is consumed, so no destructor is needed.
/// \relates ode_solver
void am_ode_solver_add_observer(ode_solver_t* solver,
                                am_ode_observer_t* observer);

///@}

#endif

