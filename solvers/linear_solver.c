// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_linearsolver.h"
#include "solvers/linear_solver.h"

struct linear_solver_t 
{
  linear_solver_type_t type;
  SUNLinearSolver s;
  linear_solver_vtable vtable;
  preconditioner_t* pc;
  int num_iters;
  real_t res_norm;
  long int last_flag;
};

static SUNLinearSolver_Type gettype(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  switch (solver->type)
  {
    case DIRECT_LINEAR_SOLVER: return SUNLINEARSOLVER_DIRECT;
    case ITERATIVE_LINEAR_SOLVER: return SUNLINEARSOLVER_MATRIX_ITERATIVE;
    case MATRIX_FREE_LINEAR_SOLVER: return SUNLINEARSOLVER_ITERATIVE;
  }
}

static int setpreconditioner(SUNLinearSolver s, 
                             void* context, 
                             PSetupFn p_setup, 
                             PSolveFn p_solve)
{
  // FIXME
}

static int numiters(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return s->num_iters;
}

static real_t resnorm(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return s->res_norm;
}

static long int lastflag(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return s->last_flag;
}

linear_solver_t* direct_linear_solver_new(MPI_Comm comm, 
                                          void* context, 
                                          linear_solver_vtable vtable)
{
  ASSERT(vtable.initialize != NULL);
  ASSERT(vtable.setup != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.space != NULL);
  ASSERT(vtable.free != NULL);
}

linear_solver_t* iterative_linear_solver_new(MPI_Comm comm, 
                                             void* context, 
                                             linear_solver_vtable vtable,
                                             preconditioner_t* p,
                                             int max_iterations,
                                             real_t tolerance)
{
  ASSERT(p != NULL);
  ASSERT(max_iterations > 0);
  ASSERT(tolerance > 0);
  ASSERT(vtable.setscalingvectors != NULL);
  ASSERT(vtable.initialize != NULL);
  ASSERT(vtable.setup != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.space != NULL);
  ASSERT(vtable.free != NULL);
}

linear_solver_t* matrix_free_linear_solver_new(MPI_Comm comm, 
                                               void* context, 
                                               linear_solver_vtable vtable,
                                               int (*A_times)(void* context, N_vector v, N_Vector Av),
                                               preconditioner_t* p,
                                               int max_iterations,
                                               real_t tolerance)
{
  ASSERT(A_times != NULL);
  ASSERT(p != NULL);
  ASSERT(max_iterations > 0);
  ASSERT(tolerance > 0);
  ASSERT(vtable.initialize != NULL);
  ASSERT(vtable.setup != NULL);
  ASSERT(vtable.solve != NULL);
  ASSERT(vtable.space != NULL);
  ASSERT(vtable.free != NULL);
}

linear_solver_t* linear_solver_from_SUNLinearSolver(SUNLinearSolver sundials_solver)
{
  linear_solver_t* solver = polymec_malloc(sizeof(linear_solver_t));
  solver->s = sundials_solver;
  return solver;
}

void linear_solver_free(linear_solver_t* solver)
{
  SUNLinSolFree(solver->s);
  polymec_free(solver);
}

void linear_solver_set_preconditioner(linear_solver_t* solver,
                                      preconditioner_t* pc)
{
  // FIXME
}

SUNLinearSolver linear_solver_as_SUNLinearSolver(linear_solver_t* solver)
{
  return solver->s;
}

linear_solver_type_t linear_solver_type(linear_solver_t* solver)
{
  return solver->type;
}

