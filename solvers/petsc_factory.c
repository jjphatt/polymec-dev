// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/petsc_factory.h"

typedef struct petsc_factory_32_t petsc_factory_32_t;
typedef struct petsc_factory_64_t petsc_factory_64_t;

struct petsc_factory_t 
{
  petsc_factory_32_t* f32;
  petsc_factory_64_t* f64;
};

bool petsc_factory_available(bool use_64_bit_indices)
{
  petsc_factory_t* factory = petsc_factory_new(use_64_bit_indices);
  bool result = (factory != NULL);
  petsc_factory_free(factory);
  return result;
}

extern petsc_factory_32_t* petsc_factory_32_new(void);
extern void petsc_factory_32_free(petsc_factory_32_t* f);
extern petsc_factory_64_t* petsc_factory_64_new(void);
extern void petsc_factory_64_free(petsc_factory_64_t* f);
petsc_factory_t* petsc_factory_new(bool use_64_bit_indices)
{
  if (use_64_bit_indices)
  {
    petsc_factory_64_t* f = petsc_factory_64_new();
    if (f != NULL)
    {
      petsc_factory_t* factory = polymec_malloc(sizeof(petsc_factory_t));
      factory->f64 = f;
      return factory;
    }
    else
      return NULL;
  else
  {
    petsc_factory_32_t* f = petsc_factory_32_new();
    if (f != NULL)
    {
      petsc_factory_t* factory = polymec_malloc(sizeof(petsc_factory_t));
      factory->f32 = f;
      return factory;
    }
    else
      return NULL;
  }
}

void petsc_factory_free(petsc_factory_t* factory)
{
  if (factory->f64 != NULL)
    petsc_factory_64_free(factory->f64);
  else
    petsc_factory_32_free(factory->f32);
  polymec_free(petsc_factory_t);
}

bool petsc_factory_has_64_bit_indices(petsc_factory_t* factory)
{
  return (factory->f64 != NULL);
}

extern matrix_t* petsc_factory_32_matrix(petsc_factory_32_t*, matrix_sparsity_t*);
extern matrix_t* petsc_factory_64_matrix(petsc_factory_64_t*, matrix_sparsity_t*);
matrix_t* petsc_factory_matrix(petsc_factory_t* factory, 
                               matrix_sparsity_t* sparsity)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_matrix(factory->f64, sparsity);
  else
    return petsc_factory_32_matrix(factory->f32, sparsity);
}

extern matrix_t* petsc_factory_32_block_matrix(petsc_factory_32_t*, matrix_sparsity_t*, int);
extern matrix_t* petsc_factory_64_block_matrix(petsc_factory_64_t*, matrix_sparsity_t*, int);
matrix_t* petsc_factory_block_matrix(petsc_factory_t* factory, 
                                     matrix_sparsity_t* sparsity,
                                     int block_size)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_block_matrix(factory->f64, sparsity, block_size);
  else
    return petsc_factory_32_block_matrix(factory->f32, sparsity, block_size);
}

extern matrix_t* petsc_factory_32_var_block_matrix(petsc_factory_32_t*, matrix_sparsity_t*, int*);
extern matrix_t* petsc_factory_64_var_block_matrix(petsc_factory_64_t*, matrix_sparsity_t*, int*);
matrix_t* petsc_factory_var_block_matrix(petsc_factory_t* factory, 
                                         matrix_sparsity_t* sparsity,
                                         int* block_sizes)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_var_block_matrix(factory->f64, sparsity, block_sizes);
  else
    return petsc_factory_32_var_block_matrix(factory->f32, sparsity, block_sizes);
}

extern nvector_t* petsc_factory_32_nvector(petsc_factory_32_t*, MPI_Comm, index_t*);
extern nvector_t* petsc_factory_64_nvector(petsc_factory_64_t*, MPI_Comm, index_t*);
nvector_t* petsc_factory_nvector(petsc_factory_t* factory,
                                 MPI_Comm comm,
                                 index_t* row_distribution)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_nvector(factory->f64, comm, row_distribution);
  else
    return petsc_factory_32_nvector(factory->f32, comm, row_distribution);
}

extern linear_solver_t* petsc_factory_32_pcg(petsc_factory_32_t*, MPI_Comm);
extern linear_solver_t* petsc_factory_64_pcg(petsc_factory_64_t*, MPI_Comm);
linear_solver_t* petsc_factory_pcg_linear_solver(petsc_factory_t* factory,
                                                 MPI_Comm comm)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_pcg(factory->f64, comm);
  else
    return petsc_factory_32_pcg(factory->f32, comm);
}
                                             
extern linear_solver_t* petsc_factory_32_gmres(petsc_factory_32_t*, MPI_Comm, int);
extern linear_solver_t* petsc_factory_64_gmres(petsc_factory_64_t*, MPI_Comm, int);
linear_solver_t* petsc_factory_gmres_linear_solver(petsc_factory_t* factory,
                                                   MPI_Comm comm,
                                                   int dimension)
{
  ASSERT(dimension > 0);
  if (factory->f64 != NULL)
    return petsc_factory_64_gmres(factory->f64, comm, dimension);
  else
    return petsc_factory_32_gmres(factory->f32, comm, dimension);
}
                                             
extern linear_solver_t* petsc_factory_32_bicgstab(petsc_factory_32_t*, MPI_Comm);
extern linear_solver_t* petsc_factory_64_bicgstab(petsc_factory_64_t*, MPI_Comm);
linear_solver_t* petsc_factory_bicgstab_solver(petsc_factory_t* factory,
                                               MPI_Comm comm)
{
  if (factory->f64 != NULL)
    return petsc_factory_64_bicgstab(factory->f64, comm);
  else
    return petsc_factory_32_bicgstab(factory->f32, comm);
}

preconditioner_t* petsc_factory_diag_scaling_preconditioner(petsc_factory_t* factory)
{
  return NULL;
}

