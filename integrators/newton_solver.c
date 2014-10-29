// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <float.h>
#include "core/sundials_helpers.h"
#include "integrators/newton_solver.h"
#include "integrators/newton_preconditioner.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

// Types of linear solvers.
typedef enum 
{
  GMRES,
  BICGSTAB,
  TFQMR
} solver_type_t;

struct newton_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  char* name;
  void* context;
  newton_solver_vtable vtable;
  solver_type_t solver_type;
  int max_krylov_dim, max_restarts;

  int N; // Number of degrees of freedom.

  // KINSol data structures.
  void* kinsol;
  int strategy; // Global strategy.
  N_Vector x, x_scale, F_scale; // Stores solution vector and scaling vectors.
  char* status_message; // status of most recent integration.

  // Preconditioning stuff.
  preconditioner_t* precond;

  // Null space information.
  bool homogeneous_functions_in_null_space;
  real_t** null_space_vectors;
  int null_dim;

  // Current simulation time.
  real_t current_time;
};

static void project_out_of_null_space(newton_solver_t* solver,
                                      real_t* R)
{
  // If homogeneous functions are in the null space, subtract the spatial 
  // mean from R.
  if (solver->homogeneous_functions_in_null_space)
  {
    real_t mean = 0.0;
    for (int i = 0; i < solver->N; ++i)
      mean += R[i];
    mean *= 1.0/solver->N;
    for (int i = 0; i < solver->N; ++i)
      R[i] -= mean;
  }

  // FIXME: Do the other stuff here.
}

// This function wraps around the user-supplied evaluation function.
static int evaluate_F(N_Vector x, N_Vector F, void* context)
{
  newton_solver_t* solver = context;
  real_t* xx = NV_DATA(x);
  real_t* FF = NV_DATA(F);

  // Evaluate the residual.
  int retval = solver->vtable.eval(solver->context, solver->current_time, xx, FF);
  project_out_of_null_space(solver, FF);
  return retval;
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(N_Vector x, N_Vector x_scale, 
                                 N_Vector F, N_Vector f_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  newton_solver_t* solver = context;
  real_t t = solver->current_time;
  newton_preconditioner_setup(solver->precond, 0.0, 1.0, 0.0, t, NV_DATA(x), NULL);
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

  // Project r out of the null space.
  project_out_of_null_space(solver, NV_DATA(r));

  if (preconditioner_solve(solver->precond, NV_DATA(r)))
    return 0;
  else 
  {
    // Recoverable error.
    log_debug("newton_solver: preconditioner solve failed.");
    return 1; 
  }
}

// Generic constructor.
static newton_solver_t* newton_solver_new(const char* name, 
                                          void* context,
                                          MPI_Comm comm,
                                          int N,
                                          newton_solver_vtable vtable,
                                          newton_solver_strategy_t global_strategy,
                                          solver_type_t solver_type,
                                          int max_krylov_dim, 
                                          int max_restarts)
{
  ASSERT(N > 0);
  ASSERT(vtable.eval != NULL);
  ASSERT(max_krylov_dim >= 3);
  ASSERT(max_restarts >= 0);

  newton_solver_t* solver = polymec_malloc(sizeof(newton_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->comm = comm;
  solver->vtable = vtable;
  solver->solver_type = solver_type;
  solver->strategy = (global_strategy == LINE_SEARCH) ? KIN_LINESEARCH : KIN_NONE;
  solver->N = N;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;

  // Set up KINSol and accessories.
  solver->kinsol = KINCreate();
  KINSetUserData(solver->kinsol, solver);
  solver->x = N_VNew(solver->comm, N);
  solver->x_scale = N_VNew(solver->comm, N);
  solver->F_scale = N_VNew(solver->comm, N);
  solver->status_message = NULL;

  KINInit(solver->kinsol, evaluate_F, solver->x);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver->solver_type == GMRES)
  {
    KINSpgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == BICGSTAB)
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
    N_Vector constraints = N_VNew(solver->comm, N);
    solver->vtable.set_constraints(solver->context, NV_DATA(constraints));
    KINSetConstraints(solver->kinsol, constraints);
    N_VDestroy(constraints);
  }

  solver->precond = NULL;
  solver->current_time = 0.0;

  // Set up the null space.
  solver->homogeneous_functions_in_null_space = false;
  solver->null_space_vectors = NULL;
  solver->null_dim = 0;

  return solver;
}

newton_solver_t* gmres_newton_solver_new(const char* name,
                                         void* context,
                                         MPI_Comm comm,
                                         int N,
                                         newton_solver_vtable vtable,
                                         newton_solver_strategy_t global_strategy,
                                         int max_krylov_dim,
                                         int max_restarts)
{
  return newton_solver_new(name, context, comm, N, vtable, global_strategy,
                                  GMRES, max_krylov_dim, max_restarts);
}

newton_solver_t* bicgstab_newton_solver_new(const char* name,
                                            void* context,
                                            MPI_Comm comm,
                                            int N,
                                            newton_solver_vtable vtable,
                                            newton_solver_strategy_t global_strategy,
                                            int max_krylov_dim)
{
  return newton_solver_new(name, context, comm, N, vtable, global_strategy,
                                  BICGSTAB, max_krylov_dim, 0);
}

newton_solver_t* tfqmr_newton_solver_new(const char* name,
                                         void* context,
                                         MPI_Comm comm,
                                         int N,
                                         newton_solver_vtable vtable,
                                         newton_solver_strategy_t global_strategy,
                                         int max_krylov_dim)
{
  return newton_solver_new(name, context, comm, N, vtable, global_strategy,
                           TFQMR, max_krylov_dim, 0);
}

void newton_solver_free(newton_solver_t* solver)
{
  // Kill the null space.
  newton_solver_set_null_space(solver, false, NULL, 0);

  // Kill the preconditioner stuff.
  if (solver->precond != NULL)
    preconditioner_free(solver->precond);

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
  polymec_free(solver->name);
  polymec_free(solver);
}

char* newton_solver_name(newton_solver_t* solver)
{
  return solver->name;
}

void* newton_solver_context(newton_solver_t* solver)
{
  return solver->context;
}

int newton_solver_num_equations(newton_solver_t* solver)
{
  return solver->N;
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

void newton_solver_set_preconditioner(newton_solver_t* solver,
                                      preconditioner_t* precond)
{
  solver->precond = precond;

  // Set up the preconditioner.
  if (solver->precond != NULL)
  {
    KINSpilsSetPreconditioner(solver->kinsol, set_up_preconditioner,
                              solve_preconditioner_system);
  }
}

preconditioner_t* newton_solver_preconditioner(newton_solver_t* solver)
{
  return solver->precond;
}

void newton_solver_set_null_space(newton_solver_t* solver,
                                  bool homogeneous_functions,
                                  real_t** null_space_vectors,
                                  int null_dim)
{
  ASSERT(((null_dim == 0) && (null_space_vectors == NULL)) ||
         ((null_dim > 0) && (null_space_vectors != NULL)));

  solver->homogeneous_functions_in_null_space = homogeneous_functions;
  if (solver->null_space_vectors != NULL)
  {
    for (int i = 0; i < solver->null_dim; ++i)
      polymec_free(solver->null_space_vectors[i]);
    polymec_free(solver->null_space_vectors);
    solver->null_dim = 0;
    solver->null_space_vectors = NULL;
  }
  if (null_space_vectors != NULL)
  {
    solver->null_dim = null_dim;
    solver->null_space_vectors = polymec_malloc(sizeof(real_t*));
    for (int i = 0; i < null_dim; ++i)
    {
      solver->null_space_vectors[i] = polymec_malloc(sizeof(real_t) * solver->N);
      memcpy(solver->null_space_vectors[i], null_space_vectors[i], solver->N * sizeof(real_t));
    }
  }
}

void newton_solver_eval_residual(newton_solver_t* solver, real_t t, real_t* X, real_t* F)
{
  solver->vtable.eval(solver->context, t, X, F);
  project_out_of_null_space(solver, F);
}

bool newton_solver_solve(newton_solver_t* solver,
                            real_t t,
                            real_t* X,
                            int* num_iterations)
{
  ASSERT(X != NULL);

  // Set the current time in the state.
  solver->current_time = t;

  // Set the x_scale and F_scale vectors. If we don't have methods for doing 
  // this, the scaling vectors are set to 1.
  int N = solver->N;
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

  if (solver->vtable.initial_guess != NULL)
  {
    // Form the initial guess magically.
    log_debug("newton_solver: forming initial guess...");
    solver->vtable.initial_guess(solver->context, t, NV_DATA(solver->x));
  }
  else
  {
    // Copy the values in X to the internal solution vector.
    memcpy(NV_DATA(solver->x), X, sizeof(real_t) * N);
  }

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
  return false;
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
