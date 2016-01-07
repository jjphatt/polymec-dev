// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/krylov_solver.h"
#include "core/timer.h"
#include "core/string_utils.h"

//------------------------------------------------------------------------
//                Krylov data types and virtual tables
//------------------------------------------------------------------------

struct krylov_solver_t
{
  char* name;
  void* context;
  krylov_solver_vtable vtable;
  krylov_pc_t* pc;
};

struct krylov_pc_t
{
  char* name;
  void* context;
  krylov_pc_vtable vtable;
};

struct krylov_matrix_t
{
  void* context;
  krylov_matrix_vtable vtable;
  int num_local_rows;
  int num_global_rows;
};

struct krylov_vector_t
{
  void* context;
  krylov_vector_vtable vtable;
  int local_size;
  index_t global_size;
};

struct krylov_factory_t
{
  char* name;
  void* context;
  krylov_factory_vtable vtable;
};

//------------------------------------------------------------------------
//                          Krylov solver
//------------------------------------------------------------------------

krylov_solver_t* krylov_solver_new(const char* name,
                                   void* context,
                                   krylov_solver_vtable vtable)
{
  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->pc = NULL;
  return solver;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  if (solver->pc != NULL)
    krylov_pc_free(solver->pc);
  string_free(solver->name);
  polymec_free(solver);
}

char* krylov_solver_name(krylov_solver_t* solver)
{
  return solver->name;
}

void* krylov_solver_impl(krylov_solver_t* solver)
{
  return solver->context;
}

void krylov_solver_set_tolerances(krylov_solver_t* solver, 
                                  real_t relative_tolerance,
                                  real_t absolute_tolerance,
                                  real_t divergence_tolerance)
{
  ASSERT(relative_tolerance > 0.0);
  ASSERT(absolute_tolerance > 0.0);
  ASSERT(divergence_tolerance > 0.0);
  solver->vtable.set_tolerances(solver->context, 
                                relative_tolerance,
                                absolute_tolerance,
                                divergence_tolerance);
}

void krylov_solver_set_max_iterations(krylov_solver_t* solver, 
                                      int max_iterations)
{
  ASSERT(max_iterations > 0);
  solver->vtable.set_max_iterations(solver->context, max_iterations);
}

void krylov_solver_set_operator(krylov_solver_t* solver, 
                                krylov_matrix_t* op)
{
  solver->vtable.set_operator(solver->context, op->context);
}

void krylov_solver_set_preconditioner(krylov_solver_t* solver,
                                      krylov_pc_t* preconditioner)
{
  if (solver->pc != NULL)
    krylov_pc_free(solver->pc);
  solver->pc = preconditioner;
}

krylov_pc_t* krylov_solver_preconditioner(krylov_solver_t* solver)
{
  return solver->pc;
}

bool krylov_solver_solve(krylov_solver_t* solver, 
                         krylov_vector_t* x, 
                         krylov_vector_t* b, 
                         real_t* residual_norm, 
                         int* num_iterations)
{
  return solver->vtable.solve(solver->context, x->context, b->context,
                              residual_norm, num_iterations);
}

//------------------------------------------------------------------------
//                          Krylov preconditioner
//------------------------------------------------------------------------

krylov_pc_t* krylov_pc_new(const char* name,
                           void* context,
                           krylov_pc_vtable vtable)
{
  krylov_pc_t* pc = polymec_malloc(sizeof(krylov_pc_t));
  pc->name = string_dup(name);
  pc->context = context;
  pc->vtable = vtable;
  return pc;
}

void krylov_pc_free(krylov_pc_t* preconditioner)
{
  if ((preconditioner->context != NULL) && (preconditioner->vtable.dtor != NULL))
    preconditioner->vtable.dtor(preconditioner->context);
  polymec_free(preconditioner);
}

char* krylov_pc_name(krylov_solver_t* preconditioner)
{
  return preconditioner->name;
}

//------------------------------------------------------------------------
//                          Krylov matrix
//------------------------------------------------------------------------

krylov_matrix_t* krylov_matrix_new(void* context,
                                   krylov_matrix_vtable vtable,
                                   int num_local_rows,
                                   index_t num_global_rows)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_diagonal != NULL);
  ASSERT(vtable.set_diagonal != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(num_local_rows > 0);
  ASSERT(num_global_rows >= num_local_rows);
  krylov_matrix_t* A = polymec_malloc(sizeof(krylov_matrix_t));
  A->context = context;
  A->vtable = vtable;
  A->num_local_rows = num_local_rows;
  A->num_global_rows = num_global_rows;
  return A;
}

void krylov_matrix_free(krylov_matrix_t* A)
{
  if ((A->context != NULL) && (A->vtable.dtor != NULL))
    A->vtable.dtor(A->context);
  polymec_free(A);
}

krylov_matrix_t* krylov_matrix_clone(krylov_matrix_t* A)
{
  krylov_matrix_t* B = polymec_malloc(sizeof(krylov_matrix_t));
  B->context = A->vtable.clone(A->context);
  B->vtable = A->vtable;
  B->num_local_rows = A->num_local_rows;
  B->num_global_rows = A->num_global_rows;
  return B;
}

void* krylov_matrix_impl(krylov_matrix_t* A)
{
  return A->context;
}

int krylov_matrix_num_local_rows(krylov_matrix_t* A)
{
  return A->num_local_rows;
}

int krylov_matrix_num_global_rows(krylov_matrix_t* A)
{
  return A->num_global_rows;
}

void krylov_matrix_zero(krylov_matrix_t* A)
{
  A->vtable.zero(A->context);
}

void krylov_matrix_add_identity(krylov_matrix_t* A,
                                real_t scale_factor)
{
  A->vtable.add_identity(A->context, scale_factor);
}

void krylov_matrix_scale(krylov_matrix_t* A,
                         real_t scale_factor)
{
  A->vtable.scale(A->context, scale_factor);
}

void krylov_matrix_add_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D)
{
  A->vtable.add_diagonal(A->context, D->context);
}

void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D)
{
  A->vtable.set_diagonal(A->context, D->context);
}

void krylov_matrix_set_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.set_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_add_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.add_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_assemble(krylov_matrix_t* A)
{
  krylov_matrix_start_assembly(A);
  krylov_matrix_finish_assembly(A);
}

void krylov_matrix_start_assembly(krylov_matrix_t* A)
{
  if (A->vtable.start_assembly != NULL)
    A->vtable.start_assembly(A->context);
}

void krylov_matrix_finish_assembly(krylov_matrix_t* A)
{
  if (A->vtable.finish_assembly != NULL)
    A->vtable.finish_assembly(A->context);
}

void krylov_matrix_get_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.get_values(A->context, num_rows, num_columns, rows, columns, values);
}

//------------------------------------------------------------------------
//                          Krylov vector
//------------------------------------------------------------------------

krylov_vector_t* krylov_vector_new(void* context,
                                   krylov_vector_vtable vtable,
                                   int local_size,
                                   index_t global_size)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(local_size > 0);
  ASSERT(global_size >= local_size);
  krylov_vector_t* v = polymec_malloc(sizeof(krylov_vector_t));
  v->context = context;
  v->vtable = vtable;
  v->local_size = local_size;
  v->global_size = global_size;
  return v;
}

void krylov_vector_free(krylov_vector_t* v)
{
  if ((v->context != NULL) && (v->vtable.dtor != NULL))
    v->vtable.dtor(v->context);
  polymec_free(v);
}

krylov_vector_t* krylov_vector_clone(krylov_vector_t* v)
{
  krylov_vector_t* u = polymec_malloc(sizeof(krylov_vector_t));
  u->context = v->vtable.clone(v->context);
  u->vtable = v->vtable;
  u->local_size = v->local_size;
  u->global_size = v->global_size;
  return u;
}

void* krylov_vector_impl(krylov_vector_t* v)
{
  return v->context;
}

int krylov_vector_local_size(krylov_vector_t* v)
{
  return v->local_size;
}

int krylov_vector_global_size(krylov_vector_t* v)
{
  return v->global_size;
}

void krylov_vector_zero(krylov_vector_t* v)
{
  v->vtable.zero(v->context);
}

void krylov_vector_set_value(krylov_vector_t* v,
                             real_t value)
{
  v->vtable.set_value(v->context, value);
}

void krylov_vector_scale(krylov_vector_t* v,
                         real_t scale_factor)
{
  v->vtable.scale(v->context, scale_factor);
}

void krylov_vector_set_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.set_values(v->context, num_values, indices, values);
}
                              
void krylov_vector_add_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.add_values(v->context, num_values, indices, values);
}

void krylov_vector_assemble(krylov_vector_t* v)
{
  krylov_vector_start_assembly(v);
  krylov_vector_finish_assembly(v);
}

void krylov_vector_start_assembly(krylov_vector_t* v)
{
  if (v->vtable.start_assembly != NULL)
    v->vtable.start_assembly(v->context);
}

void krylov_vector_finish_assembly(krylov_vector_t* v)
{
  if (v->vtable.finish_assembly != NULL)
    v->vtable.finish_assembly(v->context);
}

void krylov_vector_get_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.get_values(v->context, num_values, indices, values);
}

real_t krylov_vector_norm(krylov_vector_t* v, int p)
{
  ASSERT((p == 0) || (p == 1) || (p == 2));
  return v->vtable.norm(v->context, p);
}

//------------------------------------------------------------------------
//                  Factories for creating Krylov solvers
//------------------------------------------------------------------------

krylov_factory_t* krylov_factory_new(const char* name,
                                     void* context,
                                     krylov_factory_vtable vtable)
{
  ASSERT(vtable.pcg_solver != NULL);
  ASSERT(vtable.gmres_solver != NULL);
  ASSERT(vtable.bicgstab_solver != NULL);
  ASSERT(vtable.preconditioner != NULL);
  ASSERT(vtable.matrix != NULL);
  ASSERT(vtable.block_matrix != NULL);
  ASSERT(vtable.vector != NULL);
  krylov_factory_t* factory = polymec_malloc(sizeof(krylov_factory_t));
  factory->name = string_dup(name);
  factory->context = context;
  factory->vtable = vtable;
  return factory;
}

void krylov_factory_free(krylov_factory_t* factory)
{
  if ((factory->context != NULL) && (factory->vtable.dtor != NULL))
    factory->vtable.dtor(factory->context);
  string_free(factory->name);
  polymec_free(factory);
}

char* krylov_factory_name(krylov_factory_t* factory)
{
  return factory->name;
}

krylov_matrix_t* krylov_factory_matrix(krylov_factory_t* factory, 
                                       adj_graph_t* sparsity)
{
  return factory->vtable.matrix(factory->context, sparsity);
}

krylov_matrix_t* krylov_factory_block_matrix(krylov_factory_t* factory, 
                                             adj_graph_t* sparsity,
                                             int block_size)
{
  ASSERT(block_size > 0);
  return factory->vtable.block_matrix(factory->context, sparsity, block_size);
}

krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       adj_graph_t* dist_graph)
{
  return factory->vtable.vector(factory->context, dist_graph);
}

krylov_solver_t* krylov_factory_pcg_solver(krylov_factory_t* factory,
                                           MPI_Comm comm)
{
  return factory->vtable.pcg_solver(factory->context, comm);
}

krylov_solver_t* krylov_factory_gmres_solver(krylov_factory_t* factory,
                                             MPI_Comm comm,
                                             int krylov_dimension)
{
  ASSERT(krylov_dimension >= 3);
  return factory->vtable.gmres_solver(factory->context, comm, krylov_dimension);
}

krylov_solver_t* krylov_factory_bicgstab_solver(krylov_factory_t* factory,
                                                MPI_Comm comm)
{
  return factory->vtable.bicgstab_solver(factory->context, comm);
}

krylov_solver_t* krylov_factory_special_solver(krylov_factory_t* factory,
                                               MPI_Comm comm,
                                               const char* solver_name,
                                               string_string_unordered_map_t* options)
{
  if (factory->vtable.special_solver != NULL)
    return factory->vtable.special_solver(factory->context, comm, solver_name, options);
  else
    return NULL;
}

krylov_pc_t* krylov_factory_preconditioner(krylov_factory_t* factory,
                                           MPI_Comm comm,
                                           const char* pc_name,
                                           string_string_unordered_map_t* options)
{
  return factory->vtable.preconditioner(factory->context, comm, pc_name, options);
}

