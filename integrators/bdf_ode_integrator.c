// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "integrators/bdf_ode_integrator.h"
#include "integrators/newton_preconditioner.h"
#include "cvode/cvode.h"

// JFNK stuff.
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"
#include "cvode/cvode_sptfqmr.h"

typedef struct
{
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot);
  void (*dtor)(void* context);

  // CVODE data structures.
  void* cvode;
  N_Vector x; 
  real_t* x_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.

  // JFNK stuff.
  int max_krylov_dim;
  preconditioner_t* precond;
  int (*Jy)(void* context, real_t t, real_t* x, real_t* rhstx, real_t* y, real_t* temp, real_t* Jy);

  // Error weight function.
  bdf_ode_integrator_error_weight_func compute_weights;
} bdf_ode_t;

static char* get_status_message(int status, real_t current_time)
{
  char* status_message = NULL;
  if (status == CV_TOO_CLOSE)
    status_message = string_dup("t1 and t2 are too close to each other.");
  else if (status == CV_TOO_MUCH_WORK)
  {
    char err[1024];
    snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == CV_TOO_MUCH_ACC)
    status_message = string_dup("Integrator could not achieve desired level of accuracy.");
  else if (status == CV_ERR_FAILURE)
    status_message = string_dup("Integrator encountered too many error test failures.");
  else if (status == CV_CONV_FAILURE)
    status_message = string_dup("Integrator encountered too many convergence test failures.");
  else if (status == CV_LINIT_FAIL)
    status_message = string_dup("Integrator's linear solver failed to initialize.");
  else if (status == CV_LSETUP_FAIL)
    status_message = string_dup("Integrator's linear solver setup failed.");
  else if (status == CV_LSOLVE_FAIL)
    status_message = string_dup("Integrator's linear solver failed.");
  else if (status == CV_RHSFUNC_FAIL)
    status_message = string_dup("Integrator's RHS function failed unrecoverably.");
  else if (status == CV_FIRST_RHSFUNC_ERR)
    status_message = string_dup("Integrator's first call to RHS function failed.");
  else if (status == CV_REPTD_RHSFUNC_ERR)
    status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
  else if (status == CV_UNREC_RHSFUNC_ERR)
    status_message = string_dup("Integrator failed to recover from a recoverable RHS failure.");
  else if (status == CV_RTFUNC_FAIL)
    status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied right hand side.
static int bdf_evaluate_rhs(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  bdf_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->x_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->rhs(integ->context, t, integ->x_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->x_with_ghosts, sizeof(real_t) * integ->num_local_values);
  }

  return status;
}

static bool bdf_step(void* context, real_t max_dt, real_t* t, real_t* x)
{
  bdf_ode_t* integ = context;
  int status;

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = CVode(integ->cvode, t2, integ->x, &integ->t, CV_ONE_STEP);
    if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN))
    {
      integ->status_message = get_status_message(status, integ->t);
      return false;
    }

    if ((t2 - *t) < (integ->t - *t))
      log_detail("bdf_ode_integrator: took internal step dt = %g", integ->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (integ->t > t2)
  {
    status = CVodeGetDky(integ->cvode, t2, 0, integ->x);
    *t = t2;
  }
  else
    *t = integ->t;
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->num_local_values); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static bool bdf_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  bdf_ode_t* integ = context;
  integ->t = t1;
  CVodeReInit(integ->cvode, t1, integ->x);
  CVodeSetStopTime(integ->cvode, t2);
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->num_local_values); 

  // Integrate.
  int status = CVode(integ->cvode, t2, integ->x, &integ->t, CV_NORMAL);
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == CV_SUCCESS) || (status == CV_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->num_local_values); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static void bdf_reset(void* context, real_t t, real_t* x)
{
  bdf_ode_t* integ = context;

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->num_local_values); 
  CVodeReInit(integ->cvode, t, integ->x);
  integ->t = t;
}

static void bdf_dtor(void* context)
{
  bdf_ode_t* integ = context;

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    preconditioner_free(integ->precond);

  // Kill the CVode stuff.
  polymec_free(integ->x_with_ghosts);
  N_VDestroy(integ->x);
  CVodeFree(&integ->cvode);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  polymec_free(integ);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  bdf_ode_t* integ = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    memcpy(integ->x_with_ghosts, NV_DATA(x), sizeof(real_t) * integ->num_local_values);
    newton_preconditioner_setup(integ->precond, 1.0, -gamma, 0.0, t, integ->x_with_ghosts, NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector x, N_Vector F, 
                                       N_Vector r, N_Vector z, 
                                       real_t gamma, real_t delta, 
                                       int lr, void* context, 
                                       N_Vector work)
{
  ASSERT(lr == 1); // Left preconditioning only.

  bdf_ode_t* integ = context;
  
  // FIXME: Apply scaling if needed.
  // Copy the contents of the RHS to the output vector.
  memcpy(NV_DATA(z), NV_DATA(r), sizeof(real_t) * integ->num_local_values);

  // Solve it.
  if (preconditioner_solve(integ->precond, NV_DATA(z)))
    return 0;
  else 
    return 1; // recoverable error.
}

// Adaptor for J*y function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector x, N_Vector rhs, void* context, N_Vector tmp)
{
  bdf_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(rhs);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);
  return integ->Jy(integ->context, t, xx, my_rhs, yy, temp, jjy);
}

ode_integrator_t* jfnk_bdf_ode_integrator_new(int order,
                                             MPI_Comm comm,
                                             int num_local_values, 
                                             int num_remote_values, 
                                             void* context, 
                                             int (*rhs_func)(void* context, real_t t, real_t* x, real_t* xdot),
                                             int (*Jy_func)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
                                             void (*dtor)(void* context),
                                             preconditioner_t* precond,
                                             jfnk_bdf_krylov_t solver_type,
                                             int max_krylov_dim)
{
  ASSERT(order >= 1);
  ASSERT(order <= 12);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim > 3);

  bdf_ode_t* integ = polymec_malloc(sizeof(bdf_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = max_krylov_dim;
  integ->Jy = Jy_func;
  integ->t = 0.0;

  // Set up KINSol and accessories.
  integ->x = N_VNew(integ->comm, integ->num_local_values);
  integ->x_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetMaxOrd(integ->cvode, order);
  CVodeSetUserData(integ->cvode, integ);
  CVodeInit(integ->cvode, bdf_evaluate_rhs, 0.0, integ->x);

  // Set up the solver type.
  if (solver_type == JFNK_BDF_GMRES)
  {
    CVSpgmr(integ->cvode, PREC_LEFT, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    CVSpilsSetGSType(integ->cvode, MODIFIED_GS);
  }
  else if (solver_type == JFNK_BDF_BICGSTAB)
    CVSpbcg(integ->cvode, PREC_LEFT, max_krylov_dim);
  else
    CVSptfqmr(integ->cvode, PREC_LEFT, max_krylov_dim);

  // Set up the Jacobian function and preconditioner.
  if (Jy_func != NULL)
    CVSpilsSetJacTimesVecFn(integ->cvode, eval_Jy);
  integ->precond = precond;
  CVSpilsSetPreconditioner(integ->cvode, set_up_preconditioner,
                           solve_preconditioner_system);

  ode_integrator_vtable vtable = {.step = bdf_step, .advance = bdf_advance, .reset = bdf_reset, .dtor = bdf_dtor};
  char name[1024];
  snprintf(name, 1024, "JFNK Backwards-Difference-Formulae (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  bdf_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* bdf_ode_integrator_context(ode_integrator_t* integrator)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  return integ->context;
}

void bdf_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures)
{
  ASSERT(max_failures > 0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetMaxErrTestFails(integ->cvode, max_failures);
}

void bdf_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations)
{
  ASSERT(max_iterations > 0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetMaxNonlinIters(integ->cvode, max_iterations);
}

void bdf_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient)
{
  ASSERT(coefficient > 0.0);
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetNonlinConvCoef(integ->cvode, (double)coefficient);
}

void bdf_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                      real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  bdf_ode_t* integ = ode_integrator_context(integrator);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;

  // Set the tolerances.
  CVodeSStolerances(integ->cvode, relative_tol, absolute_tol);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  bdf_ode_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void bdf_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                 bdf_ode_integrator_error_weight_func compute_weights)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  CVodeWFtolerances(integ->cvode, compute_error_weights);
}

void bdf_integrator_set_stability_limit_detection(ode_integrator_t* integrator,
                                                  bool use_detection)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  CVodeSetStabLimDet(integ->cvode, use_detection);
}

void bdf_ode_integrator_eval_rhs(ode_integrator_t* integrator, real_t t, real_t* X, real_t* rhs)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  integ->rhs(integ->context, t, X, rhs);
}

preconditioner_t* bdf_ode_integrator_preconditioner(ode_integrator_t* integrator)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  return integ->precond;
}

void bdf_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        bdf_ode_integrator_diagnostics_t* diagnostics)
{
  bdf_ode_t* integ = ode_integrator_context(integrator);
  diagnostics->status_message = integ->status_message; // borrowed!
  CVodeGetNumSteps(integ->cvode, &diagnostics->num_steps);
  CVodeGetLastOrder(integ->cvode, &diagnostics->order_of_last_step);
  CVodeGetCurrentOrder(integ->cvode, &diagnostics->order_of_next_step);
  CVodeGetLastStep(integ->cvode, &diagnostics->last_step_size);
  CVodeGetCurrentStep(integ->cvode, &diagnostics->next_step_size);
  CVodeGetNumRhsEvals(integ->cvode, &diagnostics->num_rhs_evaluations);
  CVodeGetNumLinSolvSetups(integ->cvode, &diagnostics->num_linear_solve_setups);
  CVodeGetNumErrTestFails(integ->cvode, &diagnostics->num_error_test_failures);
  CVodeGetNumNonlinSolvIters(integ->cvode, &diagnostics->num_nonlinear_solve_iterations);
  CVodeGetNumNonlinSolvConvFails(integ->cvode, &diagnostics->num_nonlinear_solve_convergence_failures);
  CVSpilsGetNumLinIters(integ->cvode, &diagnostics->num_linear_solve_iterations);
  CVSpilsGetNumPrecEvals(integ->cvode, &diagnostics->num_preconditioner_evaluations);
  CVSpilsGetNumPrecSolves(integ->cvode, &diagnostics->num_preconditioner_solves);
  CVSpilsGetNumConvFails(integ->cvode, &diagnostics->num_linear_solve_convergence_failures);
}

void bdf_ode_integrator_diagnostics_fprintf(bdf_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "ODE integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Order of next step: %d\n", diagnostics->order_of_next_step);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Next step size: %g\n", diagnostics->next_step_size);
  fprintf(stream, "  Num RHS evaluations: %d\n", (int)diagnostics->num_rhs_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}

