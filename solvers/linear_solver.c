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
  void* context;

  MPI_Comm comm;
  SUNLinearSolver s;
  linear_solver_vtable vtable;
  preconditioner_t* pc;

  int (*A_times)(void* context, N_Vector v, N_Vector Av);

  int num_iters;
  real_t res_norm;
  long int last_flag;
  N_Vector resid;
  int max_iters;
  real_t tol;
};

static SUNLinearSolver_Type gettype(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  SUNLinearSolver_Type type;
  switch (solver->type)
  {
    case DIRECT_LINEAR_SOLVER: type = SUNLINEARSOLVER_DIRECT; break;
    case ITERATIVE_LINEAR_SOLVER: type = SUNLINEARSOLVER_MATRIX_ITERATIVE; break;
    case MATRIX_FREE_LINEAR_SOLVER: type = SUNLINEARSOLVER_ITERATIVE; 
  }
  return type;
}

static int numiters(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return solver->num_iters;
}

static real_t resnorm(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return solver->res_norm;
}

static long int lastflag(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return solver->last_flag;
}

static N_Vector resid(SUNLinearSolver s)
{
  linear_solver_t* solver = s->content;
  return solver->resid;
}

static int space(SUNLinearSolver s, long int* lrw, long int* liw)
{
  *lrw = *liw = 0;
  return 0;
}

static int setatimes(SUNLinearSolver s, void* context, ATimesFn func)
{
  // We don't do anything here, since we've already got an A*v function.
  return 0;
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

  linear_solver_t* solver = polymec_malloc(sizeof(linear_solver_t));
  solver->s = polymec_malloc(sizeof(struct _generic_SUNLinearSolver));
  solver->context = context;
  solver->type = DIRECT_LINEAR_SOLVER;
  solver->comm = comm;
  solver->s->content = solver;
  solver->pc = NULL;
  solver->A_times = NULL;
  solver->s->ops = polymec_malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  solver->s->ops->gettype = gettype;
  solver->s->ops->setatimes = NULL;
  solver->s->ops->setpreconditioner = NULL;
  solver->s->ops->setscalingvectors = vtable.setscalingvectors;
  solver->s->ops->initialize = vtable.initialize;
  solver->s->ops->setup = vtable.setup;
  solver->s->ops->solve = vtable.solve;
  solver->s->ops->numiters = numiters;
  solver->s->ops->resnorm = resnorm;
  solver->s->ops->lastflag = lastflag;
  solver->s->ops->space = (vtable.space != NULL) ? vtable.space : space;
  solver->s->ops->resid = resid;
  solver->s->ops->free = NULL;
  return solver;
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

  linear_solver_t* solver = polymec_malloc(sizeof(linear_solver_t));
  solver->s = polymec_malloc(sizeof(struct _generic_SUNLinearSolver));
  solver->context = context;
  solver->type = DIRECT_LINEAR_SOLVER;
  solver->comm = comm;
  solver->pc = p;
  solver->max_iters = max_iterations;
  solver->tol = tolerance;
  solver->A_times = NULL;
  solver->s->content = solver;
  solver->s->ops = polymec_malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  solver->s->ops->gettype = gettype;
  solver->s->ops->setatimes = NULL;
  solver->s->ops->setpreconditioner = NULL;
  solver->s->ops->setscalingvectors = vtable.setscalingvectors;
  solver->s->ops->initialize = vtable.initialize;
  solver->s->ops->setup = vtable.setup;
  solver->s->ops->solve = vtable.solve;
  solver->s->ops->numiters = numiters;
  solver->s->ops->resnorm = resnorm;
  solver->s->ops->lastflag = lastflag;
  solver->s->ops->space = (vtable.space != NULL) ? vtable.space : space;
  solver->s->ops->resid = resid;
  solver->s->ops->free = vtable.free;
  return solver;
}

linear_solver_t* matrix_free_linear_solver_new(MPI_Comm comm, 
                                               void* context, 
                                               linear_solver_vtable vtable,
                                               int (*A_times)(void* context, N_Vector v, N_Vector Av),
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

  linear_solver_t* solver = polymec_malloc(sizeof(linear_solver_t));
  solver->s = polymec_malloc(sizeof(struct _generic_SUNLinearSolver));
  solver->context = context;
  solver->type = DIRECT_LINEAR_SOLVER;
  solver->comm = comm;
  solver->pc = p;
  solver->max_iters = max_iterations;
  solver->tol = tolerance;
  solver->A_times = A_times;
  solver->s->content = solver;
  solver->s->ops = polymec_malloc(sizeof(struct _generic_SUNLinearSolver_Ops));
  solver->s->ops->gettype = gettype;
  solver->s->ops->setatimes = setatimes;
  solver->s->ops->setpreconditioner = NULL;
  solver->s->ops->setscalingvectors = vtable.setscalingvectors;
  solver->s->ops->initialize = vtable.initialize;
  solver->s->ops->setup = vtable.setup;
  solver->s->ops->solve = vtable.solve;
  solver->s->ops->numiters = numiters;
  solver->s->ops->resnorm = resnorm;
  solver->s->ops->lastflag = lastflag;
  solver->s->ops->space = (vtable.space != NULL) ? vtable.space : space;
  solver->s->ops->resid = resid;
  solver->s->ops->free = vtable.free;
  return solver;
}

linear_solver_t* linear_solver_from_SUNLinearSolver(SUNLinearSolver sundials_solver)
{
  linear_solver_t* solver = polymec_malloc(sizeof(linear_solver_t));
  solver->s = sundials_solver;
  return solver;
}

void linear_solver_free(linear_solver_t* solver)
{
  N_VDestroy(solver->resid);
  if (solver->pc != NULL)
    preconditioner_free(solver->pc);
  SUNLinSolFree(solver->s);
  if ((solver->vtable.free != NULL) && (solver->context != NULL))
    solver->vtable.free(solver->context);
  polymec_free(solver);

}

SUNLinearSolver linear_solver_as_SUNLinearSolver(linear_solver_t* solver)
{
  return solver->s;
}

linear_solver_type_t linear_solver_type(linear_solver_t* solver)
{
  return solver->type;
}

