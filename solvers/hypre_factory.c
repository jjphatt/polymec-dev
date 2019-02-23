// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/hypre_factory.h"

typedef struct hypre_factory_32_t hypre_factory_32_t;
typedef struct hypre_factory_64_t hypre_factory_64_t;

struct hypre_factory_t 
{
  hypre_factory_32_t* f32;
  hypre_factory_64_t* f64;
};

bool hypre_factory_available(bool use_64_bit_indices)
{
  hypre_factory_t* factory = hypre_factory_new(use_64_bit_indices);
  bool result = (factory != NULL);
  hypre_factory_free(factory);
  return result;
}

extern hypre_factory_32_t* hypre_factory_32_new(void);
extern void hypre_factory_32_free(hypre_factory_32_t* f);
extern hypre_factory_64_t* hypre_factory_64_new(void);
extern void hypre_factory_64_free(hypre_factory_64_t* f);
hypre_factory_t* hypre_factory_new(bool use_64_bit_indices)
{
  if (use_64_bit_indices)
  {
    hypre_factory_64_t* f = hypre_factory_64_new();
    if (f != NULL)
    {
      hypre_factory_t* factory = polymec_malloc(sizeof(hypre_factory_t));
      factory->f64 = f;
      return factory;
    }
    else
      return NULL;
  else
  {
    hypre_factory_32_t* f = hypre_factory_32_new();
    if (f != NULL)
    {
      hypre_factory_t* factory = polymec_malloc(sizeof(hypre_factory_t));
      factory->f32 = f;
      return factory;
    }
    else
      return NULL;
  }
}

void hypre_factory_free(hypre_factory_t* factory)
{
  if (factory->f64 != NULL)
    hypre_factory_64_free(factory->f64);
  else
    hypre_factory_32_free(factory->f32);
  polymec_free(hypre_factory_t);
}

bool hypre_factory_has_64_bit_indices(hypre_factory_t* factory)
{
  return (factory->f64 != NULL);
}

extern matrix_t* hypre_factory_32_matrix(hypre_factory_32_t*, matrix_sparsity_t*);
extern matrix_t* hypre_factory_64_matrix(hypre_factory_64_t*, matrix_sparsity_t*);
matrix_t* hypre_factory_matrix(hypre_factory_t* factory, 
                               matrix_sparsity_t* sparsity)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_matrix(factory->f64, sparsity);
  else
    return hypre_factory_32_matrix(factory->f32, sparsity);
}

extern matrix_t* hypre_factory_32_block_matrix(hypre_factory_32_t*, matrix_sparsity_t*, int);
extern matrix_t* hypre_factory_64_block_matrix(hypre_factory_64_t*, matrix_sparsity_t*, int);
matrix_t* hypre_factory_block_matrix(hypre_factory_t* factory, 
                                     matrix_sparsity_t* sparsity,
                                     int block_size)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_block_matrix(factory->f64, sparsity, block_size);
  else
    return hypre_factory_32_block_matrix(factory->f32, sparsity, block_size);
}

extern matrix_t* hypre_factory_32_var_block_matrix(hypre_factory_32_t*, matrix_sparsity_t*, int*);
extern matrix_t* hypre_factory_64_var_block_matrix(hypre_factory_64_t*, matrix_sparsity_t*, int*);
matrix_t* hypre_factory_var_block_matrix(hypre_factory_t* factory, 
                                         matrix_sparsity_t* sparsity,
                                         int* block_sizes)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_var_block_matrix(factory->f64, sparsity, block_sizes);
  else
    return hypre_factory_32_var_block_matrix(factory->f32, sparsity, block_sizes);
}

extern nvector_t* hypre_factory_32_nvector(hypre_factory_32_t*, MPI_Comm, index_t*);
extern nvector_t* hypre_factory_64_nvector(hypre_factory_64_t*, MPI_Comm, index_t*);
nvector_t* hypre_factory_nvector(hypre_factory_t* factory,
                                 MPI_Comm comm,
                                 index_t* row_distribution)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_nvector(factory->f64, comm, row_distribution);
  else
    return hypre_factory_32_nvector(factory->f32, comm, row_distribution);
}

extern linear_solver_t* hypre_factory_32_pcg(hypre_factory_32_t*, MPI_Comm);
extern linear_solver_t* hypre_factory_64_pcg(hypre_factory_64_t*, MPI_Comm);
linear_solver_t* hypre_factory_pcg_linear_solver(hypre_factory_t* factory,
                                                 MPI_Comm comm)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_pcg(factory->f64, comm);
  else
    return hypre_factory_32_pcg(factory->f32, comm);
}
                                             
extern linear_solver_t* hypre_factory_32_gmres(hypre_factory_32_t*, MPI_Comm, int);
extern linear_solver_t* hypre_factory_64_gmres(hypre_factory_64_t*, MPI_Comm, int);
linear_solver_t* hypre_factory_gmres_linear_solver(hypre_factory_t* factory,
                                                   MPI_Comm comm,
                                                   int dimension)
{
  ASSERT(dimension > 0);
  if (factory->f64 != NULL)
    return hypre_factory_64_gmres(factory->f64, comm, dimension);
  else
    return hypre_factory_32_gmres(factory->f32, comm, dimension);
}
                                             
extern linear_solver_t* hypre_factory_32_bicgstab(hypre_factory_32_t*, MPI_Comm);
extern linear_solver_t* hypre_factory_64_bicgstab(hypre_factory_64_t*, MPI_Comm);
linear_solver_t* hypre_factory_bicgstab_solver(hypre_factory_t* factory,
                                               MPI_Comm comm)
{
  if (factory->f64 != NULL)
    return hypre_factory_64_bicgstab(factory->f64, comm);
  else
    return hypre_factory_32_bicgstab(factory->f32, comm);
}

preconditioner_t* hypre_factory_diag_scaling_preconditioner(hypre_factory_t* factory)
{
  return NULL;
}

