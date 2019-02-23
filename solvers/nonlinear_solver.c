// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_nonlinearsolver.h"
#include "solvers/nonlinear_solver.h"

struct nonlinear_solver_t 
{
  SUNNonlinearSolver s;
};

nonlinear_solver_t* nonlinear_solver_from_SUNNonlinearSolver(SUNNonlinearSolver sundials_solver)
{
  nonlinear_solver_t* solver = polymec_malloc(sizeof(nonlinear_solver_t));
  solver->s = sundials_solver;
  return solver;
}

void nonlinear_solver_free(nonlinear_solver_t* solver)
{
  SUNNonlinSolFree(solver->s);
  polymec_free(solver);
}

SUNNonlinearSolver nonlinear_solver_as_SUNNonlinearSolver(nonlinear_solver_t* solver)
{
  return solver->s;
}

