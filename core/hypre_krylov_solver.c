// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include <dlfcn.h>
#include "core/krylov_solver.h"
#include "core/timer.h"
#include "core/text_buffer.h"
#include "core/string_utils.h"

//------------------------------------------------------------------------
// This file implements the dynamically-loadable HYPRE Krylov solver.
//------------------------------------------------------------------------

// HYPRE types and definitions.
#define HYPRE_PARCSR  5555
typedef double HYPRE_Real;
typedef double _Complex HYPRE_Complex;
typedef long long HYPRE_Int;
typedef void* HYPRE_Solver;
typedef void* HYPRE_Matrix;
typedef void* HYPRE_Vector;
typedef void* HYPRE_IJMatrix;
typedef void* HYPRE_IJVector;
typedef void* HYPRE_ParCSRMatrix;
typedef void* HYPRE_ParVector;
typedef HYPRE_Int (*HYPRE_PtrToSolverFcn)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector);
typedef HYPRE_Int (*HYPRE_PtrToParSolverFcn)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
#if POLYMEC_HAVE_MPI
typedef MPI_Comm HYPRE_MPI_Comm;
#else
typedef HYPRE_Int HYPRE_MPI_Comm;
#endif

// Types of supported solver methods.
typedef enum
{
  HYPRE_PCG,
  HYPRE_GMRES,
  HYPRE_BICGSTAB,
  HYPRE_BOOMERANG,
} hypre_solver_method_t;

// Here's a table of function pointers for the HYPRE library.
typedef struct
{
  HYPRE_Int (*HYPRE_GetError)();
  HYPRE_Int (*HYPRE_ClearAllErrors)();

  HYPRE_Int (*HYPRE_ParCSRPCGCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRPCGDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRPCGSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRPCGSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRPCGSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRPCGSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRPCGSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRPCGSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRPCGSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRPCGGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRPCGGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRPCGGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);

  HYPRE_Int (*HYPRE_ParCSRGMRESCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRGMRESDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRGMRESSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetKDim)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);

  HYPRE_Int (*HYPRE_ParCSRBiCGSTABCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);

  HYPRE_Int (*HYPRE_BoomerAMGCreate)(HYPRE_Solver*);
  HYPRE_Int (*HYPRE_BoomerAMGDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_BoomerAMGSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_BoomerAMGSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_BoomerAMGGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_BoomerAMGGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);
  HYPRE_Int (*HYPRE_BoomerAMGSetNumFunctions)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetDofFunc)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_BoomerAMGSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxCoarseSize)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMinCoarseSize)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxLevels)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetStrongThreshold)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxRowSum)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetCoarsenType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMeasureType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetAggNumLevels)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNumPaths)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetCGCIts)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNodal)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNodalDiag)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetInterpType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetTruncFactor)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetPMaxElmts)(HYPRE_Solver, HYPRE_Int);

  HYPRE_Int (*HYPRE_EuclidCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_EuclidDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_EuclidSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_EuclidSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_EuclidSetLevel)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetBJ)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetSparseA)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_EuclidSetRowScale)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetILUT)(HYPRE_Solver, HYPRE_Real);

  HYPRE_Int (*HYPRE_ParaSailsCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParaSailsDestroy)(HYPRE_Solver); 
  HYPRE_Int (*HYPRE_ParaSailsSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParaSailsSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParaSailsSetParams)(HYPRE_Solver, HYPRE_Real, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParaSailsSetFilter)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParaSailsSetSym)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParaSailsSetLoadbal)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParaSailsSetReuse)(HYPRE_Solver, HYPRE_Int);

  HYPRE_Int (*HYPRE_IJMatrixCreate)(HYPRE_MPI_Comm, HYPRE_Int, HYPRE_Int, HYPRE_Int, HYPRE_Int, HYPRE_IJMatrix*);
  HYPRE_Int (*HYPRE_IJMatrixDestroy)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixInitialize)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixSetValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, const HYPRE_Int*, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixAddToValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, const HYPRE_Int*, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixAssemble)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixGetValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, HYPRE_Int*, HYPRE_Int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixSetObjectType)(HYPRE_IJMatrix, HYPRE_Int);
  HYPRE_Int (*HYPRE_IJMatrixGetObject)(HYPRE_IJMatrix, void**);
  HYPRE_Int (*HYPRE_IJMatrixSetRowSizes)(HYPRE_IJMatrix, const HYPRE_Int*);
  HYPRE_Int (*HYPRE_IJMatrixSetDiagOffdSizes)(HYPRE_IJMatrix, const HYPRE_Int*, const HYPRE_Int*);
  HYPRE_Int (*HYPRE_IJMatrixGetRowCounts)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, HYPRE_Int*);
  HYPRE_Int (*HYPRE_IJMatrixPrint)(HYPRE_IJMatrix, const char*);

  HYPRE_Int (*HYPRE_IJVectorCreate)(HYPRE_MPI_Comm, HYPRE_Int, HYPRE_Int, HYPRE_IJVector*);
  HYPRE_Int (*HYPRE_IJVectorDestroy)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorInitialize)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorSetValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAddToValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAssemble)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorGetValues)(HYPRE_IJVector, HYPRE_Int, HYPRE_Int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorSetObjectType)(HYPRE_IJVector, HYPRE_Int);
  HYPRE_Int (*HYPRE_IJVectorGetObject)(HYPRE_IJVector, void**);
  HYPRE_Int (*HYPRE_IJVectorPrint)(HYPRE_IJVector, const char*);

  // Boneheaded built-in preconditioner.
  HYPRE_Int (*HYPRE_ParCSRDiagScaleSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRDiagScale)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

  // Internals needed for scaling and zeroing matrices.
  HYPRE_Int (*HYPRE_ParCSRMatrixGetRow)(void*, HYPRE_Int, HYPRE_Int*, HYPRE_Int**, HYPRE_Real**);
  HYPRE_Int (*HYPRE_ParCSRMatrixRestoreRow)(void*, HYPRE_Int, HYPRE_Int*, HYPRE_Int**, HYPRE_Real**);
} hypre_methods_table;

typedef struct
{
  void* hypre;
  hypre_methods_table methods;
} hypre_factory_t;

typedef struct
{
  hypre_factory_t* factory;
  hypre_solver_method_t type;
  HYPRE_Solver solver;
  HYPRE_ParCSRMatrix op;
} hypre_solver_t;

typedef struct
{
  HYPRE_Solver pc;
  HYPRE_Int (*setup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*solve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*dtor)(HYPRE_Solver);
} hypre_pc_t;

typedef struct
{
  MPI_Comm comm;
  index_t ilow, ihigh;
  hypre_factory_t* factory;
  HYPRE_IJMatrix A;
  bool initialized;
} hypre_matrix_t;

typedef struct
{
  MPI_Comm comm;
  index_t ilow, ihigh;
  hypre_factory_t* factory;
  HYPRE_IJVector v;
  bool initialized;
} hypre_vector_t;

static void hypre_solver_set_tolerances(void* context,
                                        real_t rel_tol,
                                        real_t abs_tol,
                                        real_t div_tol)

{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*set_tol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*set_abs_tol)(HYPRE_Solver, HYPRE_Real);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_tol = solver->factory->methods.HYPRE_ParCSRPCGSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRPCGSetAbsoluteTol; 
      break;
    case HYPRE_GMRES: 
      set_tol = solver->factory->methods.HYPRE_ParCSRGMRESSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRGMRESSetAbsoluteTol; 
      break;
    case HYPRE_BICGSTAB:
      set_tol = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetAbsoluteTol; 
      break;
    case HYPRE_BOOMERANG:
      set_tol = solver->factory->methods.HYPRE_BoomerAMGSetTol; 
      set_abs_tol = NULL;
  }
  set_tol(solver->solver, (HYPRE_Real)rel_tol);
  if (set_abs_tol != NULL)
    set_abs_tol(solver->solver, (HYPRE_Real)abs_tol);
}

static void hypre_solver_set_max_iterations(void* context,
                                            int max_iters)
{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*set_max_iters)(HYPRE_Solver, HYPRE_Int);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_max_iters = solver->factory->methods.HYPRE_ParCSRPCGSetMaxIter; 
      break;
    case HYPRE_GMRES: 
      set_max_iters = solver->factory->methods.HYPRE_ParCSRGMRESSetMaxIter; 
      break;
    case HYPRE_BICGSTAB:
      set_max_iters = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetMaxIter; 
      break;
    case HYPRE_BOOMERANG:
      set_max_iters = solver->factory->methods.HYPRE_BoomerAMGSetMaxIter; 
  }
  set_max_iters(solver->solver, (HYPRE_Int)max_iters);
}

static void hypre_solver_set_operator(void* context,
                                      void* op)
{
  hypre_solver_t* solver = context;
  HYPRE_IJMatrix A = ((hypre_matrix_t*)op)->A;
  solver->factory->methods.HYPRE_IJMatrixGetObject(A, &solver->op);

  // NOTE: We can't do setup here, since we don't have x or b.
}

static void hypre_solver_set_pc(void* context,
                                void* pc)
{
  hypre_solver_t* solver = context;
  hypre_pc_t* p = pc;

  // Dispatch.
  HYPRE_Int (*set_pc)(HYPRE_Solver, HYPRE_PtrToSolverFcn, HYPRE_PtrToSolverFcn, HYPRE_Solver);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_pc = solver->factory->methods.HYPRE_ParCSRPCGSetPrecond; 
      break;
    case HYPRE_GMRES: 
      set_pc = solver->factory->methods.HYPRE_ParCSRGMRESSetPrecond; 
      break;
    case HYPRE_BICGSTAB:
      set_pc = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetPrecond; 
      break;
    default:
      set_pc = NULL;
  }
  if (set_pc != NULL)
    set_pc(solver, (HYPRE_PtrToSolverFcn)p->solve, (HYPRE_PtrToSolverFcn)p->setup, p->pc);
}

static bool hypre_solver_solve(void* context,
                               void* b,
                               void* x,
                               real_t* res_norm,
                               int* num_iters)
{
  hypre_solver_t* solver = context;

  // Get our objects.
  HYPRE_IJVector B = ((hypre_vector_t*)b)->v;
  HYPRE_IJVector X = ((hypre_vector_t*)x)->v;
  HYPRE_ParVector par_X, par_B;
  solver->factory->methods.HYPRE_IJVectorGetObject(X, &par_X);
  solver->factory->methods.HYPRE_IJVectorGetObject(B, &par_B);

  // Dispatch.
  HYPRE_Int (*setup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*solve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*get_num_iters)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*get_norm)(HYPRE_Solver, HYPRE_Real*);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      setup = solver->factory->methods.HYPRE_ParCSRPCGSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRPCGSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRPCGGetNumIterations;
      get_norm = solver->factory->methods.HYPRE_ParCSRPCGGetFinalRelativeResidualNorm;
      break;
    case HYPRE_GMRES: 
      setup = solver->factory->methods.HYPRE_ParCSRGMRESSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRGMRESSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRGMRESGetNumIterations;
      get_norm = solver->factory->methods.HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;
      break;
    case HYPRE_BICGSTAB:
      setup = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRBiCGSTABSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRBiCGSTABGetNumIterations; 
      get_norm = solver->factory->methods.HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm; 
      break;
    case HYPRE_BOOMERANG:
      setup = solver->factory->methods.HYPRE_BoomerAMGSetup; 
      solve = solver->factory->methods.HYPRE_BoomerAMGSolve; 
      get_num_iters = solver->factory->methods.HYPRE_BoomerAMGGetNumIterations; 
      get_norm = solver->factory->methods.HYPRE_BoomerAMGGetFinalRelativeResidualNorm; 
  }
  log_debug("hypre_solver_solve: Setting up linear system.");
  HYPRE_Int result = setup(solver->solver, solver->op, par_B, par_X);
  if (result != 0)
    return false;

  log_debug("hypre_solver_solve: Solving linear system.");

  // HYPRE's Krylov methods can produce some pretty small numbers, which can be 
  // interpreted as denormalized garbage by polymec, so we tell polymec to chill.
  polymec_suspend_fpe();
  result = solve(solver->solver, solver->op, par_B, par_X);
  polymec_restore_fpe();

  bool success = (result == 0);
  if (success)
    log_debug("hypre_solver_solve: Success!");
  else
    log_debug("hypre_solver_solve: Failure (error code = %d)", result);

  HYPRE_Int iters;
  get_num_iters(solver->solver, &iters);
  *num_iters = iters;
  HYPRE_Real norm;
  get_norm(solver->solver, &norm);
  *res_norm = norm;

  return success;
}

static void hypre_solver_dtor(void* context)
{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*destroy)(HYPRE_Solver);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      destroy = solver->factory->methods.HYPRE_ParCSRPCGDestroy; 
      break;
    case HYPRE_GMRES: 
      destroy = solver->factory->methods.HYPRE_ParCSRGMRESDestroy; 
      break;
    case HYPRE_BICGSTAB:
      destroy = solver->factory->methods.HYPRE_ParCSRBiCGSTABDestroy; 
      break;
    case HYPRE_BOOMERANG:
      destroy = solver->factory->methods.HYPRE_BoomerAMGDestroy; 
  }
  destroy(solver->solver);
  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* hypre_factory_pcg_solver(void* context,
                                                 MPI_Comm comm)
     
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_PCG;
  HYPRE_MPI_Comm hypre_comm = comm;
  HYPRE_Int result;
  result = solver->factory->methods.HYPRE_ParCSRPCGCreate(hypre_comm, &solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_pcg_solver: Failed to create PCG solver (error = %d).", result);
    polymec_free(solver);
    return NULL;
  }
  log_debug("hypre_factory_pcg_solver: Created PCG solver");

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRPCGSetPrecond(solver->solver, 
                                                     (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                     (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                     solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_pcg_solver: Failed to set up diagonal scaling PC (error = %d).", result);
    solver->factory->methods.HYPRE_ParCSRPCGDestroy(solver->solver);
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE PCG", solver, vtable);
}

static krylov_solver_t* hypre_factory_gmres_solver(void* context,
                                                   MPI_Comm comm,
                                                   int krylov_dimension)
     
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_GMRES;
  HYPRE_MPI_Comm hypre_comm = comm;
  HYPRE_Int result;
  result = solver->factory->methods.HYPRE_ParCSRGMRESCreate(hypre_comm, &solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_gmres_solver: Failed to create GMRES solver (error = %d).", result);
    polymec_free(solver);
    return NULL;
  }
  result = solver->factory->methods.HYPRE_ParCSRGMRESSetKDim(solver->solver, (HYPRE_Int)krylov_dimension);
  if (result != 0)
  {
    log_urgent("hypre_factory_gmres_solver: Failed to set Krylov dimension (error = %d).", result);
    polymec_free(solver);
    return NULL;
  }
  log_debug("hypre_factory_gmres_solver: Created solver with Krylov dim = %d", krylov_dimension);

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRGMRESSetPrecond(solver->solver, 
                                                       (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                       (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                       solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_gmres_solver: Failed to set up diagonal scaling PC (error = %d).", result);
    solver->factory->methods.HYPRE_ParCSRGMRESDestroy(solver->solver);
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE GMRES", solver, vtable);
}

static krylov_solver_t* hypre_factory_bicgstab_solver(void* context,
                                                      MPI_Comm comm)
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_BICGSTAB;
  HYPRE_MPI_Comm hypre_comm = comm;
  HYPRE_Int result = solver->factory->methods.HYPRE_ParCSRBiCGSTABCreate(hypre_comm, &solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_bicgstab_solver: Failed to create BiCGSTAB solver (error = %d).", result);
    polymec_free(solver);
    return NULL;
  }

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRBiCGSTABSetPrecond(solver->solver, 
                                                          (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                          (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                          solver->solver);
  if (result != 0)
  {
    log_urgent("hypre_factory_bicgstab_solver: Failed to set up diagonal scaling PC (error = %d).", result);
    solver->factory->methods.HYPRE_ParCSRBiCGSTABDestroy(solver->solver);
    polymec_free(solver);
    return NULL;
  }


  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE Bi-CGSTAB", solver, vtable);
}

static krylov_solver_t* hypre_factory_special_solver(void* context,
                                                     MPI_Comm comm,
                                                     const char* solver_name,
                                                     string_string_unordered_map_t* options)
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;

  HYPRE_Int result;
  if ((string_casecmp(solver_name, "boomerang") == 0) ||
      (string_casecmp(solver_name, "boomeramg") == 0))
  {
    solver->type = HYPRE_BOOMERANG;
    result = solver->factory->methods.HYPRE_BoomerAMGCreate(&solver->solver);
  }
  else
    return NULL;

  if (result != 0)
  {
    log_urgent("hypre_factory_special_solver: Failed to create %s solver (error = %d).", solver_name, result);
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new(solver_name, solver, vtable);
}

static void hypre_pc_free(void* pc)
{
  hypre_pc_t* p = pc;
  if ((p->dtor != NULL) && (p->pc != NULL))
    p->dtor(p->pc);
  polymec_free(p);
}

static krylov_pc_t* hypre_factory_pc(void* context,
                                     MPI_Comm comm,
                                     const char* pc_name,
                                     string_string_unordered_map_t* options)
{
  hypre_factory_t* factory = context;
  hypre_pc_t* pc = polymec_malloc(sizeof(hypre_pc_t));

  if ((string_casecmp(pc_name, "boomerang") == 0) ||
      (string_casecmp(pc_name, "boomeramg") == 0) || 
      (string_casecmp(pc_name, "amg") == 0))
  {
    factory->methods.HYPRE_BoomerAMGCreate(&pc->pc);
    pc->setup = factory->methods.HYPRE_BoomerAMGSetup;
    pc->solve = factory->methods.HYPRE_BoomerAMGSolve;
    pc->dtor = factory->methods.HYPRE_BoomerAMGDestroy;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "num_functions") == 0) && string_is_number(value))
          factory->methods.HYPRE_BoomerAMGSetNumFunctions(pc->pc, atoi(value));
        // FIXME: etc etc
      }
    }
  }
  else if (string_casecmp(pc_name, "euclid") == 0)
  {
    HYPRE_MPI_Comm hypre_comm = comm;
    factory->methods.HYPRE_EuclidCreate(hypre_comm, &pc->pc);
    pc->setup = factory->methods.HYPRE_EuclidSetup;
    pc->solve = factory->methods.HYPRE_EuclidSolve;
    pc->dtor = factory->methods.HYPRE_EuclidDestroy;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "level") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetLevel(pc->pc, atoi(value));
        else if (strcmp(key, "bj") == 0)
          factory->methods.HYPRE_EuclidSetBJ(pc->pc, string_as_boolean(value));
        else if ((strcmp(key, "sparseA") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetSparseA(pc->pc, atof(value));
        else if (strcmp(key, "rowScale") == 0)
          factory->methods.HYPRE_EuclidSetRowScale(pc->pc, string_as_boolean(value));
        else if ((strcmp(key, "ilut") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetILUT(pc->pc, atof(value));
      }
    }
  }
  else if (string_casecmp(pc_name, "parasails") == 0)
  {
    HYPRE_MPI_Comm hypre_comm = comm;
    factory->methods.HYPRE_ParaSailsCreate(hypre_comm, &pc->pc);
    pc->setup = factory->methods.HYPRE_ParaSailsSetup;
    pc->solve = factory->methods.HYPRE_ParaSailsSolve;
    pc->dtor = factory->methods.HYPRE_ParaSailsDestroy;
    HYPRE_Int nlevel = 1;
    HYPRE_Real thresh = 0.1;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "symmetry") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetSym(pc->pc, atoi(value));
        else if ((strcmp(key, "thresh") == 0) && string_is_number(value))
          thresh = atof(value);
        else if ((strcmp(key, "nlevel") == 0) && string_is_number(value))
          nlevel = atoi(value);
        else if ((strcmp(key, "filter") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetFilter(pc->pc, atof(value));
        else if ((strcmp(key, "loadbal") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetLoadbal(pc->pc, atof(value));
        else if (strcmp(key, "reuse") == 0)
          factory->methods.HYPRE_ParaSailsSetReuse(pc->pc, string_as_boolean(value));
      }
    }
    factory->methods.HYPRE_ParaSailsSetParams(pc->pc, thresh, nlevel);
  }
  else
    return NULL;

  // Set up the virtual table.
  krylov_pc_vtable vtable = {.dtor = hypre_pc_free};
  return krylov_pc_new(pc_name, pc, vtable);
}

static void hypre_matrix_set_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  if (!A->initialized)
  {
    A->factory->methods.HYPRE_IJMatrixInitialize(A->A);
    A->initialized = true;
  }
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    A->factory->methods.HYPRE_IJMatrixSetValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, (HYPRE_Real*)values);
  else
  {
    // We have to convert to HYPRE_Real.
    int tot_num_values = 0;
    for (int r = 0; r < num_rows; ++r)
      tot_num_values += num_columns[r];
    HYPRE_Real dvals[tot_num_values];
    for (int i = 0; i < tot_num_values; ++i)
      dvals[i] = (HYPRE_Real)values[i];
    A->factory->methods.HYPRE_IJMatrixSetValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, dvals);
  }
  HYPRE_Int error = A->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_matrix_set_values: error occurred (%d)", error);
}

static void hypre_matrix_add_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  if (!A->initialized)
  {
    A->factory->methods.HYPRE_IJMatrixInitialize(A->A);
    A->initialized = true;
  }
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    A->factory->methods.HYPRE_IJMatrixAddToValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, values);
  else
  {
    // We have to convert to HYPRE_Real.
    int tot_num_values = 0;
    for (int r = 0; r < num_rows; ++r)
      tot_num_values += num_columns[r];
    HYPRE_Real dvals[tot_num_values];
    for (int i = 0; i < tot_num_values; ++i)
      dvals[i] = (HYPRE_Real)values[i];
    A->factory->methods.HYPRE_IJMatrixAddToValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, dvals);
  }
  HYPRE_Int error = A->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_matrix_add_values: error occurred (%d)", error);
}

static void hypre_matrix_start_assembly(void* context)
{
  // Nothing to do!
}

static void hypre_matrix_finish_assembly(void* context)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixAssemble(A->A);
  A->initialized = false;
  HYPRE_Int error = A->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_matrix_finish_assembly: error occurred (%d)", error);
}

static void* hypre_matrix_clone(void* context)
{
  hypre_matrix_t* A = context;
  hypre_matrix_t* clone = polymec_malloc(sizeof(hypre_matrix_t));
  clone->comm = A->comm;
  clone->factory = A->factory;
  clone->ilow = A->ilow;
  clone->ihigh = A->ihigh;
  HYPRE_MPI_Comm hypre_comm = A->comm;
  clone->factory->methods.HYPRE_IJMatrixCreate(hypre_comm, 
                                               clone->ilow, clone->ihigh, 
                                               clone->ilow, clone->ihigh,
                                               &clone->A);
  clone->factory->methods.HYPRE_IJMatrixSetObjectType(clone->A, HYPRE_PARCSR);
  clone->factory->methods.HYPRE_IJMatrixInitialize(clone->A);

  // Copy over the nonzero structure and values of A.
  index_t num_rows = clone->ihigh - clone->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows];
  index_t tot_num_columns = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int row = clone->ilow + r, ncols;
    clone->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
    num_columns[r] = ncols;
    rows[r] = row;
    tot_num_columns += ncols;
  }
  index_t columns[tot_num_columns];
  real_t values[tot_num_columns];
  void* par_mat;
  clone->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
  index_t col_offset = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int ncols, *cols;
    HYPRE_Real* vals;
    HYPRE_Int row = clone->ilow + r;
    clone->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
    for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
    {
      columns[col_offset] = cols[c];
      values[col_offset] = (real_t)vals[c];
    }
    clone->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
  }
  ASSERT(col_offset == tot_num_columns);

  // Stuff the values in.
  clone->initialized = false;
  hypre_matrix_set_values(clone, num_rows, num_columns, rows, columns, values);
  hypre_matrix_finish_assembly(clone);

  return clone;
}

static void hypre_matrix_zero(void* context)
{
  hypre_matrix_t* A = context;

  // Get the information about the matrix's nonzero structure.
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows];
  index_t tot_num_columns = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int row = A->ilow + r, ncols;
    A->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
    tot_num_columns += ncols; 
  }
  index_t columns[tot_num_columns];
  void* par_mat;
  A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
  index_t col_offset = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int ncols, *cols;
    HYPRE_Real* vals;
    index_t row = A->ilow + r;
    A->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
    for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
      columns[col_offset] = cols[c];
    A->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
  }

  // Stuff a bunch of zeros in.
  real_t zeros[tot_num_columns];
  memset(zeros, 0, sizeof(real_t) * tot_num_columns);
  hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, zeros);
  hypre_matrix_finish_assembly(context);
}

static void hypre_matrix_scale(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;

  // Get the information about the matrix's nonzero structure.
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows];
  index_t tot_num_columns = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int row = A->ilow + r, ncols;
    A->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
    tot_num_columns += ncols; 
  }
  index_t columns[tot_num_columns];
  real_t values[tot_num_columns];
  void* par_mat;
  A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
  index_t col_offset = 0;
  for (index_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int ncols, *cols;
    HYPRE_Real* vals;
    index_t row = A->ilow + r;
    A->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
    for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
    {
      columns[col_offset] = cols[c];
      values[col_offset] = (real_t)(scale_factor * vals[c]);
    }
    A->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
  }

  // Scale the values in the matrix.
  hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, values);
  hypre_matrix_finish_assembly(context);
}

static void hypre_matrix_add_identity(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
    values[r] = scale_factor;
  }
  hypre_matrix_add_values(context, num_rows, num_columns, rows, columns, values); 
  hypre_matrix_finish_assembly(context);
}

// This is used for the diagonal methods below, and defined later in the file.
static void hypre_vector_get_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values);

static void hypre_matrix_set_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
  }
  hypre_vector_get_values(context, num_rows, rows, values); // get values from D
  hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, values); // set diag A.
  hypre_matrix_finish_assembly(context);
}

static void hypre_matrix_add_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
  }
  hypre_vector_get_values(context, num_rows, rows, values); // get values from D
  hypre_matrix_add_values(context, num_rows, num_columns, rows, columns, values); // add to diag A.
  hypre_matrix_finish_assembly(context);
}

static void hypre_matrix_get_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    A->factory->methods.HYPRE_IJMatrixGetValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, values);
  else
  {
    // We have to convert from HYPRE_Real.
    int tot_num_values = 0;
    for (int r = 0; r < num_rows; ++r)
      tot_num_values += num_columns[r];
    HYPRE_Real dvals[tot_num_values];
    A->factory->methods.HYPRE_IJMatrixGetValues(A->A, num_rows, (HYPRE_Int*)num_columns, (HYPRE_Int*)rows, (HYPRE_Int*)columns, dvals);
    for (int i = 0; i < tot_num_values; ++i)
      values[i] = (real_t)dvals[i];
  }
}

static void hypre_matrix_dtor(void* context)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixDestroy(A->A);
  A->factory = NULL;
  polymec_free(A);
}

static krylov_matrix_t* hypre_factory_block_matrix(void* context,
                                                   adj_graph_t* sparsity,
                                                   int block_size)
{
  ASSERT(block_size >= 1);

  hypre_matrix_t* A = polymec_malloc(sizeof(hypre_matrix_t));
  A->factory = context;
  A->comm = adj_graph_comm(sparsity);
  int rank, nprocs;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &nprocs);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  A->ilow = vtx_dist[rank];
  A->ihigh = vtx_dist[rank+1] - 1;
  HYPRE_MPI_Comm hypre_comm = A->comm;
  A->factory->methods.HYPRE_IJMatrixCreate(hypre_comm, 
                                           A->ilow, A->ihigh, A->ilow, A->ihigh,
                                           &A->A);
  A->factory->methods.HYPRE_IJMatrixSetObjectType(A->A, HYPRE_PARCSR);

  // Preallocate non-zero storage.
  HYPRE_Int N_local = block_size * (vtx_dist[rank+1] - vtx_dist[rank]);
  if (nprocs == 1)
  {
    HYPRE_Int nnz[N_local];
    for (int v = 0; v < N_local; ++v)
      nnz[v] = block_size * (HYPRE_Int)(1 + adj_graph_num_edges(sparsity, v/block_size));
    A->factory->methods.HYPRE_IJMatrixSetRowSizes(A->A, nnz);
  }
  else
  {
    HYPRE_Int d_nnz[N_local], o_nnz[N_local];
    for (int v = 0; v < N_local; ++v)
    {
      d_nnz[v] = block_size; // Diagonal entry.
      o_nnz[v] = 0; // Diagonal entry.
      int num_edges = adj_graph_num_edges(sparsity, v);
      int* edges = adj_graph_edges(sparsity, v);
      for (int e = 0; e < num_edges; ++e)
      {
        int edge = edges[e];
        if (edge >= N_local)
          o_nnz[v] += block_size;
        else
          d_nnz[v] += block_size;
      }
    }
    A->factory->methods.HYPRE_IJMatrixSetDiagOffdSizes(A->A, d_nnz, o_nnz);
  }

  // We're ready to set values.
  A->factory->methods.HYPRE_IJMatrixInitialize(A->A);
  A->initialized = true;

  // Set the "non-zero" values to zero initially. This constructs the specific non-zero structure.
  index_t num_rows = A->ihigh - A->ilow + 1;
  index_t num_columns[num_rows], rows[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    rows[r] = A->ilow + r;
    num_columns[r] = block_size * (1 + adj_graph_num_edges(sparsity, r));
  }
  int tot_num_values = 0;
  for (int r = 0; r < num_rows; ++r)
    tot_num_values += num_columns[r];

  index_t col_offset = 0;
  index_t columns[tot_num_values];
  for (int r = 0; r < num_rows; ++r)
  {
    columns[col_offset++] = A->ilow + r;
    int* edges = adj_graph_edges(sparsity, r);
    for (int c = 0; c < num_columns[r]-1; ++c, ++col_offset)
      columns[col_offset] = A->ilow + edges[c];
  }
  ASSERT(col_offset == tot_num_values);

  real_t zeros[tot_num_values];
  memset(zeros, 0, sizeof(real_t) * tot_num_values);
  hypre_matrix_set_values(A, num_rows, num_columns, rows, columns, zeros);

  // Assemble the matrix.
  hypre_matrix_finish_assembly(A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = hypre_matrix_clone,
                                 .zero = hypre_matrix_zero,
                                 .scale = hypre_matrix_scale,
                                 .add_identity = hypre_matrix_add_identity,
                                 .add_diagonal = hypre_matrix_add_diagonal,
                                 .set_diagonal = hypre_matrix_set_diagonal,
                                 .set_values = hypre_matrix_set_values,
                                 .add_values = hypre_matrix_add_values,
                                 .start_assembly = hypre_matrix_start_assembly,
                                 .finish_assembly = hypre_matrix_finish_assembly,
                                 .get_values = hypre_matrix_get_values,
                                 .dtor = hypre_matrix_dtor};
  HYPRE_Int N_global = block_size * vtx_dist[nprocs];
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static krylov_matrix_t* hypre_factory_matrix(void* context,
                                             adj_graph_t* sparsity)
{
  return hypre_factory_block_matrix(context, sparsity, 1);
}

static void hypre_vector_set_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    v->factory->methods.HYPRE_IJVectorSetValues(v->v, num_values, (HYPRE_Int*)indices, values);
  else
  {
    // We need to convert to HYPRE_Real.
    HYPRE_Real dvals[num_values];
    for (int i = 0; i < num_values; ++i)
      dvals[i] = (HYPRE_Real)values[i];
    v->factory->methods.HYPRE_IJVectorSetValues(v->v, num_values, (HYPRE_Int*)indices, dvals);
  }
  HYPRE_Int error = v->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_vector_set_values: error occurred (%d)", error);
}

static void hypre_vector_add_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    v->factory->methods.HYPRE_IJVectorAddToValues(v->v, num_values, (HYPRE_Int*)indices, values);
  else
  {
    // We need to convert to HYPRE_Real.
    HYPRE_Real dvals[num_values];
    for (int i = 0; i < num_values; ++i)
      dvals[i] = (HYPRE_Real)values[i];
    v->factory->methods.HYPRE_IJVectorAddToValues(v->v, num_values, (HYPRE_Int*)indices, dvals);
  }
  HYPRE_Int error = v->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_vector_add_values: error occurred (%d)", error);
}

static void hypre_vector_start_assembly(void* context)
{
  // Nothing to do!
}

static void hypre_vector_finish_assembly(void* context)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorAssemble(v->v);
  v->initialized = false;
  HYPRE_Int error = v->factory->methods.HYPRE_GetError();
  if (error != 0)
    log_urgent("hypre_vector_finish_assembly: error occurred (%d)", error);
}

static void hypre_vector_get_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;
  if (sizeof(real_t) == sizeof(HYPRE_Real))
    v->factory->methods.HYPRE_IJVectorGetValues(v->v, num_values, (HYPRE_Int*)indices, values);
  else
  {
    // We have to convert from HYPRE_Real.
    HYPRE_Real dvals[num_values];
    v->factory->methods.HYPRE_IJVectorGetValues(v->v, num_values, (HYPRE_Int*)indices, dvals);
    for (int i = 0; i < num_values; ++i)
      values[i] = (real_t)dvals[i];
  }
}

static void* hypre_vector_clone(void* context)
{
  hypre_vector_t* v = context;

  // Get data from the original vector and use it to create a new one.
  hypre_vector_t* clone = polymec_malloc(sizeof(hypre_vector_t));
  clone->comm = v->comm;
  clone->factory = v->factory;
  clone->ilow = v->ilow;
  clone->ihigh = v->ihigh;
  clone->initialized = v->initialized;
  HYPRE_MPI_Comm hypre_comm = clone->comm;
  clone->factory->methods.HYPRE_IJVectorCreate(hypre_comm, clone->ilow, clone->ihigh, &clone->v);
  clone->factory->methods.HYPRE_IJVectorSetObjectType(clone->v, HYPRE_PARCSR);
  clone->factory->methods.HYPRE_IJVectorInitialize(clone->v);

  index_t num_rows = clone->ihigh - clone->ilow + 1;
  index_t rows[num_rows];
  for (index_t r = 0; r < num_rows; ++r)
    rows[r] = clone->ilow + r;
  real_t vals[num_rows];
  hypre_vector_get_values(v, num_rows, rows, vals);
  hypre_vector_set_values(clone, num_rows, rows, vals);
  return clone;
}


static void hypre_vector_set_value(void* context, real_t value)
{
  hypre_vector_t* v = context;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
  {
    rows[r] = v->ilow + r;
    values[r] = value;
  }
  hypre_vector_set_values(context, num_rows, rows, values);
  hypre_vector_finish_assembly(context);
}

static void hypre_vector_zero(void* context)
{
  hypre_vector_set_value(context, 0.0);
}

static void hypre_vector_scale(void* context, real_t scale_factor)
{
  hypre_vector_t* v = context;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, values);
  for (index_t r = 0; r < num_rows; ++r) 
    values[r] *= scale_factor;
  hypre_vector_set_values(context, num_rows, rows, values);
  hypre_vector_finish_assembly(context);
}

static real_t hypre_vector_norm(void* context, int p)
{
  // SIGH. HYPRE doesn't do vector norms, so we have to do this manually.
  hypre_vector_t* v = context;
  
  // Accumulate the local part of the norm.
  real_t local_norm = 0.0;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, values);
  if (p == 0)
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm = MAX(local_norm, values[r]);
  }
  else if (p == 1)
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm += ABS(values[r]);
  }
  else
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm += values[r] * values[r];
  }

  // Now mash together all the parallel portions.
  real_t global_norm = 0.0;
  if (p == 0)
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_MAX, v->comm);
  else
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_SUM, v->comm);

  if (p == 2)
    global_norm = sqrt(global_norm);
  return global_norm;
}

static void hypre_vector_dtor(void* context)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorDestroy(v->v);
  v->factory = NULL;
  polymec_free(v);
}

static krylov_vector_t* hypre_factory_vector(void* context,
                                             adj_graph_t* dist_graph)
{
  hypre_vector_t* v = polymec_malloc(sizeof(hypre_vector_t));
  v->factory = context;
  v->comm = adj_graph_comm(dist_graph);
  int rank, nprocs;
  MPI_Comm_rank(v->comm, &rank);
  MPI_Comm_size(v->comm, &nprocs);
  index_t* vtx_dist = adj_graph_vertex_dist(dist_graph);
  v->ilow = vtx_dist[rank];
  v->ihigh = vtx_dist[rank+1]-1;
  HYPRE_MPI_Comm hypre_comm = v->comm;
  v->factory->methods.HYPRE_IJVectorCreate(hypre_comm, v->ilow, v->ihigh, &v->v);
  v->factory->methods.HYPRE_IJVectorSetObjectType(v->v, HYPRE_PARCSR);
  v->factory->methods.HYPRE_IJVectorInitialize(v->v);
  v->initialized = true;

  // Set all the entries of the vector to zero.
  hypre_vector_zero(v);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = hypre_vector_clone,
                                 .zero = hypre_vector_zero,
                                 .set_value = hypre_vector_set_value,
                                 .scale = hypre_vector_scale,
                                 .set_values = hypre_vector_set_values,
                                 .add_values = hypre_vector_add_values,
                                 .start_assembly = hypre_vector_start_assembly,
                                 .finish_assembly = hypre_vector_finish_assembly,
                                 .get_values = hypre_vector_get_values,
                                 .norm = hypre_vector_norm,
                                 .dtor = hypre_vector_dtor};
  return krylov_vector_new(v, vtable, vtx_dist[rank+1]-vtx_dist[rank], vtx_dist[nprocs]);
}

static void hypre_factory_dtor(void* context)
{
  hypre_factory_t* factory = context;
  log_debug("hypre_krylov_factory: Closing HYPRE library.");
  dlclose(factory->hypre);
  polymec_free(factory);
}

// Use this to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr, fail_label) \
  { \
    void* ptr = dlsym(dylib, symbol_name); \
    if (ptr == NULL) \
    { \
      log_urgent("%s: unable to find %s in dynamic library.", __func__, symbol_name); \
      goto fail_label; \
    } \
    *((void**)&(function_ptr)) = ptr; \
  } 

krylov_factory_t* hypre_krylov_factory(const char* hypre_dir)
{
  hypre_factory_t* factory = polymec_malloc(sizeof(hypre_factory_t));

  // Try to find HYPRE.
  char hypre_path[FILENAME_MAX+1];
  snprintf(hypre_path, FILENAME_MAX, "%s/libHYPRE%s", hypre_dir, SHARED_LIBRARY_SUFFIX);

  // Try to open libPETSc and mine it for symbols.
  log_debug("hypre_krylov_factory: Opening HYPRE library at %s.", hypre_path);
  void* hypre = dlopen(hypre_path, RTLD_NOW);
  if (hypre == NULL)
  {
    char* msg = dlerror();
    polymec_error("hypre_krylov_factory: %s.", msg);
  }

#if APPLE
  // Mac-specific trick: 
  // If the HYPRE library is parallel, it will contain a reference to MPI symbols, 
  // even if those symbols are not defined within the library itself.
  bool hypre_is_parallel = (dlsym(hypre, "MPI_Abort") != NULL);
#if POLYMEC_HAVE_MPI
  if (!hypre_is_parallel)
  {
    log_urgent("hypre_krylov_factory: Polymec is configured to use MPI, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=OFF."); 
    goto failure;
  }
#else
  if (hypre_is_parallel)
  {
    log_urgent("hypre_krylov_factory: Polymec is configured serially, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=ON."); 
    goto failure;
  }
#endif
#endif
  log_debug("hypre_krylov_factory: Succeeded.");

  // HYPRE defines a function called hypre_printf if it is configured for "big" (64-bit) integers.
  bool hypre_uses_64bit_indices = (dlsym(hypre, "hypre_printf") != NULL);
  if (sizeof(index_t) == sizeof(int64_t))
  {
    if (!hypre_uses_64bit_indices)
    {
      log_urgent("hypre_krylov_factory: Since polymec is configured for 64-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=ON.");
      goto failure;
    }
  }
  else
  {
    if (hypre_uses_64bit_indices)
    {
      log_urgent("hypre_krylov_factory: Since polymec is configured for 32-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=OFF.");
      goto failure;
    }
  }

  // Get the symbols.
#define FETCH_HYPRE_SYMBOL(symbol_name) \
  FETCH_SYMBOL(hypre, #symbol_name, factory->methods.symbol_name, failure);

  FETCH_HYPRE_SYMBOL(HYPRE_GetError);
  FETCH_HYPRE_SYMBOL(HYPRE_ClearAllErrors);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetFinalRelativeResidualNorm);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetKDim);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm);

  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGGetFinalRelativeResidualNorm);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNumFunctions);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetDofFunc);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxCoarseSize);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMinCoarseSize);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxLevels);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetStrongThreshold);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxRowSum);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetCoarsenType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMeasureType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetAggNumLevels);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNumPaths);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetCGCIts);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNodal);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNodalDiag);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetInterpType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetTruncFactor);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetPMaxElmts);

  FETCH_HYPRE_SYMBOL(HYPRE_EuclidCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetLevel);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetBJ);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetSparseA);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetRowScale);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetILUT);

  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetParams);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetFilter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetSym);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetLoadbal);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetReuse);

  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetObjectType);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetObject);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetRowSizes);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetDiagOffdSizes);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetRowCounts);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixPrint);

  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetObjectType);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetObject);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorPrint);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRDiagScaleSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRDiagScale);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixGetRow);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixRestoreRow);
#undef FETCH_HYPRE_SYMBOL

  log_debug("hypre_krylov_factory: Got HYPRE symbols.");

  // Stash the library.
  factory->hypre = hypre; 

  // Construct the factory.
  krylov_factory_vtable vtable = {.pcg_solver = hypre_factory_pcg_solver,
                                  .gmres_solver = hypre_factory_gmres_solver,
                                  .bicgstab_solver = hypre_factory_bicgstab_solver,
                                  .special_solver = hypre_factory_special_solver,
                                  .preconditioner = hypre_factory_pc,
                                  .matrix = hypre_factory_matrix,
                                  .block_matrix = hypre_factory_block_matrix,
                                  .vector = hypre_factory_vector,
                                  .dtor = hypre_factory_dtor};
  return krylov_factory_new("2.10", factory, vtable);

failure:
  dlclose(hypre);
  polymec_free(factory);
  return NULL;
}

