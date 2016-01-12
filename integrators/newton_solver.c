// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/sundials_helpers.h"
#include "core/timer.h"
#include "integrators/newton_solver.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spfgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

struct newton_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  void* context;
  newton_solver_vtable vtable;
  newton_krylov_t solver_type;
  int max_krylov_dim, max_restarts;

  int num_local_values, num_remote_values;

  // KINSol data structures.
  void* kinsol;
  int strategy; // Global strategy.
  N_Vector x, x_scale, F_scale; // Stores solution vector and scaling vectors.
  real_t* x_with_ghosts;
  char* status_message; // status of most recent integration.

  // Jacobian-vector product function.
  newton_solver_Jv_func Jv_func;

  // Preconditioning stuff.
  newton_pc_t* precond;

  // Current simulation time.
  real_t t;
};

// This function wraps around the user-supplied evaluation function.
static int evaluate_F(N_Vector x, N_Vector F, void* context)
{
  newton_solver_t* solver = context;
  real_t* xx = NV_DATA(x);
  real_t* FF = NV_DATA(F);

  // Evaluate the residual using a solution vector with ghosts.
  memcpy(solver->x_with_ghosts, xx, sizeof(real_t) * solver->num_local_values);
  return solver->vtable.eval(solver->context, solver->t, solver->x_with_ghosts, FF);
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  newton_solver_t* solver = context;
  real_t t = solver->t;
  newton_pc_setup(solver->precond, 0.0, 1.0, 0.0, t, NV_DATA(x), NULL);

  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(N_Vector x, N_Vector x_scale,
                                       N_Vector F, N_Vector F_scale,
                                       N_Vector r, void* context,
                                       N_Vector work)
{
  newton_solver_t* solver = context;

  // FIXME: Apply scaling if needed.

  if (newton_pc_solve(solver->precond, solver->t, NV_DATA(x), NULL,
                      NV_DATA(r), NV_DATA(work)))
  {
    // Copy the solution to r.
    memcpy(NV_DATA(r), NV_DATA(work), sizeof(real_t) * NV_LOCLENGTH(r));
    return 0;
  }
  else 
  {
    // Recoverable error.
    log_debug("newton_solver: preconditioner solve failed.");
    return 1; 
  }
}

// Generic constructor.
newton_solver_t* newton_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   newton_solver_vtable vtable,
                                   newton_pc_t* precond,
                                   newton_krylov_t solver_type,
                                   int max_krylov_dim, 
                                   int max_restarts)
{
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(vtable.eval != NULL);
  ASSERT(max_krylov_dim >= 3);
  ASSERT(((solver_type != NEWTON_GMRES) && (solver_type != NEWTON_FGMRES)) || 
         (max_restarts >= 0));
  ASSERT(precond != NULL);
  ASSERT(newton_pc_side(precond) == NEWTON_PC_LEFT);

  newton_solver_t* solver = polymec_malloc(sizeof(newton_solver_t));
  solver->context = context;
  solver->comm = comm;
  solver->vtable = vtable;
  solver->precond = precond;
  solver->solver_type = solver_type;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;

  // By default, we take a full Newton step.
  solver->strategy = KIN_NONE;

  // Set up KINSol and accessories.
  solver->kinsol = KINCreate();
  KINSetUserData(solver->kinsol, solver);
  solver->x = N_VNew(solver->comm, num_local_values);
  solver->x_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->x_scale = N_VNew(solver->comm, num_local_values);
  solver->F_scale = N_VNew(solver->comm, num_local_values);
  solver->status_message = NULL;

  KINInit(solver->kinsol, evaluate_F, solver->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver->solver_type == NEWTON_GMRES)
  {
    KINSpgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == NEWTON_FGMRES)
  {
    KINSpfgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == NEWTON_BICGSTAB)
    KINSpbcg(solver->kinsol, solver->max_krylov_dim);
  else
    KINSptfqmr(solver->kinsol, solver->max_krylov_dim);

  // Enable debugging diagnostics if logging permits.
  FILE* info_stream = log_stream(LOG_DEBUG);
  if (info_stream != NULL)
  {
    KINSetPrintLevel(solver->kinsol, 3);
    KINSetInfoFile(solver->kinsol, info_stream);
  }
  else
  {
    KINSetPrintLevel(solver->kinsol, 0);
    KINSetInfoFile(solver->kinsol, NULL);
  }

  // Set the constraints (if any) for the solution.
  if (solver->vtable.set_constraints != NULL)
  {
    N_Vector constraints = N_VNew(solver->comm, num_local_values);
    solver->vtable.set_constraints(solver->context, NV_DATA(constraints));
    KINSetConstraints(solver->kinsol, constraints);
    N_VDestroy(constraints);
  }

  // Set up the preconditioner.
  if (solver->precond != NULL)
  {
    KINSpilsSetPreconditioner(solver->kinsol, set_up_preconditioner,
                              solve_preconditioner_system);
  }

  // By default, we use a difference-quotient approximation to J*v.
  solver->Jv_func = NULL;

  solver->t = 0.0;

  return solver;
}

void newton_solver_free(newton_solver_t* solver)
{
  // Kill the preconditioner stuff.
  if (solver->precond != NULL)
    newton_pc_free(solver->precond);

  // Kill the KINSol stuff.
  N_VDestroy(solver->x);
  N_VDestroy(solver->x_scale);
  N_VDestroy(solver->F_scale);
  KINFree(&solver->kinsol);

  // Kill the rest.
  if ((solver->vtable.dtor != NULL) && (solver->context != NULL))
    solver->vtable.dtor(solver->context);
  // Kill the rest.
  if (solver->status_message != NULL)
    polymec_free(solver->status_message);
  polymec_free(solver->x_with_ghosts);
  polymec_free(solver);
}

// This is a KINSOL Jacobian-vector product function that wraps our own.
static int Jv_func_wrapper(N_Vector v, N_Vector Jv, N_Vector x,
                           booleantype* new_x, void* context)
{
  newton_solver_t* solver = context;
  ASSERT(solver->Jv_func != NULL);
  int result = solver->Jv_func(solver->context, (*new_x != 0), solver->t, 
                               NV_DATA(x), NV_DATA(v), NV_DATA(Jv));
  if ((result == 0) && (*new_x == 0))
    *new_x = 1;
  return result;
}

void newton_solver_set_jacobian_vector_product(newton_solver_t* solver,
                                               newton_solver_Jv_func Jv_func)
{
  solver->Jv_func = Jv_func;
  if (Jv_func != NULL)
    KINSpilsSetJacTimesVecFn(solver->kinsol, Jv_func_wrapper);
  else
    KINSpilsSetJacTimesVecFn(solver->kinsol, NULL);
}

void newton_solver_use_full_step(newton_solver_t* solver)
{
  log_debug("newton_solver: using full Newton step.");
  solver->strategy = KIN_NONE;
}

void newton_solver_use_line_search(newton_solver_t* solver)
{
  log_debug("newton_solver: using line search.");
  solver->strategy = KIN_LINESEARCH;
}

void newton_solver_use_fixed_point(newton_solver_t* solver, 
                                   int num_residuals)
{
  ASSERT(num_residuals >= 0);
  log_debug("newton_solver: using fixed point Anderson acceleration (%d residuals).", num_residuals);
  solver->strategy = KIN_FP;
  KINSetMAA(solver->kinsol, num_residuals);
}

void newton_solver_use_picard(newton_solver_t* solver, 
                              newton_solver_Jv_func Jv_func,
                              int num_residuals)
{
  ASSERT(num_residuals >= 0);
  log_debug("newton_solver: using Picard Anderson acceleration (%d residuals).", num_residuals);
  solver->strategy = KIN_PICARD;
  newton_solver_set_jacobian_vector_product(solver, Jv_func);
  KINSetMAA(solver->kinsol, num_residuals);
}

void* newton_solver_context(newton_solver_t* solver)
{
  return solver->context;
}

int newton_solver_num_equations(newton_solver_t* solver)
{
  return solver->num_local_values;
}

void newton_solver_set_tolerances(newton_solver_t* solver, real_t norm_tolerance, real_t step_tolerance)
{
  ASSERT(norm_tolerance > 0.0);
  ASSERT(step_tolerance > 0.0);
  KINSetFuncNormTol(solver->kinsol, norm_tolerance);
  KINSetScaledStepTol(solver->kinsol, step_tolerance);
}

void newton_solver_set_max_iterations(newton_solver_t* solver, int max_iterations)
{
  ASSERT(max_iterations > 0);
  KINSetNumMaxIters(solver->kinsol, max_iterations);
}

void newton_solver_set_linear_solver_stopping_criteria(newton_solver_t* solver,
                                                       newton_solver_stopping_criteria_t criteria)
{
  if (criteria == NEWTON_EISENSTAT_WALKER1)
    KINSetEtaForm(solver->kinsol, KIN_ETACHOICE1);
  else if (criteria == NEWTON_EISENSTAT_WALKER2)
    KINSetEtaForm(solver->kinsol, KIN_ETACHOICE2);
  else // (criteria == NEWTON_CONSTANT_ETA)
    KINSetEtaForm(solver->kinsol, KIN_ETACONSTANT);
}

void newton_solver_set_constant_eta(newton_solver_t* solver, real_t eta)
{
  ASSERT(eta > 0.0);
  KINSetEtaConstValue(solver->kinsol, eta);
}

newton_pc_t* newton_solver_preconditioner(newton_solver_t* solver)
{
  return solver->precond;
}

void newton_solver_eval_residual(newton_solver_t* solver, real_t t, real_t* X, real_t* F)
{
  solver->vtable.eval(solver->context, t, X, F);
}

bool newton_solver_solve(newton_solver_t* solver,
                         real_t t,
                         real_t* X,
                         int* num_iterations)
{
  START_FUNCTION_TIMER();
  ASSERT(X != NULL);

  // Set the current time in the state.
  solver->t = t;

  // Set the x_scale and F_scale vectors. If we don't have methods for doing 
  // this, the scaling vectors are set to 1.
  int N = solver->num_local_values;
  if (solver->vtable.set_x_scale != NULL)
    solver->vtable.set_x_scale(solver->context, NV_DATA(solver->x_scale));
  else
  {
    for (int i = 0; i < N; ++i)
      NV_Ith(solver->x_scale, i) = 1.0;
  }

  if (solver->vtable.set_F_scale != NULL)
    solver->vtable.set_F_scale(solver->context, NV_DATA(solver->F_scale));
  else
  {
    for (int i = 0; i < N; ++i)
      NV_Ith(solver->F_scale, i) = 1.0;
  }

  // Copy the values in X to the internal solution vector.
  memcpy(NV_DATA(solver->x), X, sizeof(real_t) * N);

  // Suspend the currently active floating point exceptions for now.
//  polymec_suspend_fpe_exceptions();

  // Solve.
  log_debug("newton_solver: solving...");
  int status = KINSol(solver->kinsol, solver->x, solver->strategy, 
                      solver->x_scale, solver->F_scale);

  // Clear the present status.
  if (solver->status_message != NULL)
  {
    polymec_free(solver->status_message);
    solver->status_message = NULL;
  }

  // Reinstate the floating point exceptions.
//  polymec_restore_fpe_exceptions();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(solver->kinsol, &num_iters);
    *num_iterations = (int)num_iters;
    log_debug("newton_solver: solved after %d iterations.", *num_iterations);

    // Copy the data back into X.
    memcpy(X, NV_DATA(solver->x), sizeof(real_t) * N);
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    if (status == KIN_STEP_LT_STPTOL)
      solver->status_message = string_dup("Nonlinear solve stalled because scaled Newton step is too small.");
    else if (status == KIN_LINESEARCH_NONCONV)
      solver->status_message = string_dup("Line search could not sufficiently decrease the error of the iterate.");
    else if (status == KIN_MAXITER_REACHED)
      solver->status_message = string_dup("Maximum number of nonlinear iterations was reached.");
    else if (status == KIN_MXNEWT_5X_EXCEEDED)
      solver->status_message = string_dup("Maximum Newton step size was exceeded 5 times.");
    else if (status == KIN_LINESEARCH_BCFAIL)
      solver->status_message = string_dup("Line search could not satisfy beta condition.");
    else if (status == KIN_LINSOLV_NO_RECOVERY)
      solver->status_message = string_dup("Preconditioner solve encountered a recoverable error after update.");
    else if (status == KIN_LINIT_FAIL)
      solver->status_message = string_dup("Linear solve setup failed.");
    else if (status == KIN_LSETUP_FAIL)
      solver->status_message = string_dup("Preconditioner setup failed unrecoverably.");
    else if (status == KIN_LSOLVE_FAIL)
      solver->status_message = string_dup("Linear solve failed (or preconditioner solve failed unrecoverably).");
    else if (status == KIN_SYSFUNC_FAIL)
      solver->status_message = string_dup("Nonlinear function evaluation failed unrecoverably.");
    else if (status == KIN_FIRST_SYSFUNC_ERR)
      solver->status_message = string_dup("First nonlinear function evaluation failed recoverably.");
    else if (status == KIN_REPTD_SYSFUNC_ERR)
      solver->status_message = string_dup("Nonlinear function evaluation repeatedly failed (no recovery possible).");
  }

  // Failed!
  STOP_FUNCTION_TIMER();
  return false;
}
                                  
void newton_solver_reset(newton_solver_t* solver, real_t t)
{
  // Reset the preconditioner.
  newton_pc_reset(solver->precond, t);
}

void newton_solver_get_diagnostics(newton_solver_t* solver, 
                                   newton_solver_diagnostics_t* diagnostics)
{
  diagnostics->status_message = solver->status_message; // borrowed!
  KINGetNumFuncEvals(solver->kinsol, &diagnostics->num_function_evaluations);
  KINGetNumBetaCondFails(solver->kinsol, &diagnostics->num_beta_condition_failures);
  KINGetNumNonlinSolvIters(solver->kinsol, &diagnostics->num_nonlinear_iterations);
  KINGetNumBacktrackOps(solver->kinsol, &diagnostics->num_backtrack_operations);
  KINGetFuncNorm(solver->kinsol, &diagnostics->scaled_function_norm);
  KINGetStepLength(solver->kinsol, &diagnostics->scaled_newton_step_length);
  KINSpilsGetNumLinIters(solver->kinsol, &diagnostics->num_linear_solve_iterations);
  KINSpilsGetNumConvFails(solver->kinsol, &diagnostics->num_linear_solve_convergence_failures);
  KINSpilsGetNumPrecEvals(solver->kinsol, &diagnostics->num_preconditioner_evaluations);
  KINSpilsGetNumPrecSolves(solver->kinsol, &diagnostics->num_preconditioner_solves);
  KINSpilsGetNumJtimesEvals(solver->kinsol, &diagnostics->num_jacobian_vector_product_evaluations);
  KINSpilsGetNumFuncEvals(solver->kinsol, &diagnostics->num_difference_quotient_function_evaluations);
}

void newton_solver_diagnostics_fprintf(newton_solver_diagnostics_t* diagnostics, 
                                       FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "Nonlinear solver diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num function evaluations: %d\n", (int)diagnostics->num_function_evaluations);
  fprintf(stream, "  Num beta condition failures: %d\n", (int)diagnostics->num_beta_condition_failures);
  fprintf(stream, "  Num backtrack operations: %d\n", (int)diagnostics->num_backtrack_operations);
  fprintf(stream, "  Num nonlinear iterations: %d\n", (int)diagnostics->num_nonlinear_iterations);
  fprintf(stream, "  Scaled function norm: %g\n", (double)diagnostics->scaled_function_norm);
  fprintf(stream, "  Scaled Newton step length: %g\n", (double)diagnostics->scaled_newton_step_length);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
  fprintf(stream, "  Num Jacobian-vector product evaluations: %d\n", (int)diagnostics->num_jacobian_vector_product_evaluations);
  fprintf(stream, "  Num difference quotient function evaluations: %d\n", (int)diagnostics->num_difference_quotient_function_evaluations);
}

