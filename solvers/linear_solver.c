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
  SUNLinearSolver s;
};

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

