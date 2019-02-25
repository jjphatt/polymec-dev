// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_HYPRE_FACTORY_H
#define POLYMEC_HYPRE_FACTORY_H

#include "solvers/nvector.h"
#include "solvers/matrix.h"
#include "solvers/linear_solver.h"
#include "solvers/nonlinear_solver.h"
#include "solvers/matrix_sparsity.h"

/// \addtogroup solvers solvers
///@{

/// \class hypre_factory
/// This factory provides access to the solvers, matrices, and vectors in 
/// the HYPRE library.
/// \refcounted
typedef struct hypre_factory_t hypre_factory_t;

/// Returns true if the HYPRE factory is available, depending on whether 
/// 64-bit indices are expected. False otherwise. HYPRE libraries are 
/// sought using polymec's dynamic loading search path.
/// \param [in] use_64_bit_indices If true, this function returns true only if 
///                                a hypre library is available with support 
///                                for 64-bit indices.
/// \related hypre_factory
bool hypre_factory_available(bool use_64_bit_indices);

/// Creates a HYPRE factory that can be used to construct solvers, matrices, 
/// and vectors.
/// \param [in] use_64_bit_indices If true, this function attempts to construct
///                                a hypre factory for a HYPRE installation that
///                                supports 64-bit indices.
/// \returns A newly created factory object, or NULL if no factory can be 
///          constructed.
/// \memberof hypre_factory
hypre_factory_t* hypre_factory_new(bool use_64_bit_indices);

/// Returns true if the factory supports 64-bit indices, false if not.
/// \memberof hypre_factory
bool hypre_factory_has_64_bit_indices(hypre_factory_t* factory);

/// Constructs a (square) sparse matrix with the given sparsity pattern.
/// \param [in] sparsity The desired matrix sparsity pattern.
/// \memberof hypre_factory
matrix_t* hypre_factory_matrix(hypre_factory_t* factory, 
                               matrix_sparsity_t* sparsity);

/// Constructs a (square) sparse block matrix with the given sparsity pattern 
/// and block size.
/// \param [in] sparsity The desired matrix sparsity pattern.
/// \param [in] block_size The desired block size for the matrix.
/// \memberof hypre_factory
matrix_t* hypre_factory_block_matrix(hypre_factory_t* factory, 
                                     matrix_sparsity_t* sparsity,
                                     int block_size);

/// Constructs a (square) sparse block matrix with the given sparsity pattern 
/// and variable block sizes (one per matrix row).
/// \param [in] sparsity The desired matrix sparsity pattern.
/// \param [in] block_sizes The desired block sizes for the matrix rows.
/// \memberof hypre_factory
matrix_t* hypre_factory_var_block_matrix(hypre_factory_t* factory, 
                                         matrix_sparsity_t* sparsity,
                                         int* block_sizes);

/// Constructs a vector on the given communicator with its local and global 
/// dimensions defined by the given row distribution array (which is 
/// nprocs + 1 in length). 
/// \param [in] comm The MPI communicator for the newly-created n-vector.
/// \param [in] row_distribution An array whose pth element gives the global 
///                              index of the first row stored on process p. 
///                              The number of rows local to that process is 
///                              given by row_distribution[p+1] - row_distribution[p].
/// \memberof hypre_factory
nvector_t* hypre_factory_nvector(hypre_factory_t* factory,
                                 MPI_Comm comm,
                                 index_t* row_distribution);

/// Constructs a linear solver that uses the preconditioned conjugate gradient 
/// (PCG) method. This method can only be used for systems having symmetric, 
/// positive-definite matrices.
/// \param [in] comm The communicator on which the solver operates.
/// \memberof hypre_factory
linear_solver_t* hypre_factory_pcg_linear_solver(hypre_factory_t* factory,
                                                 MPI_Comm comm);
                                             
/// Constructs a linear solver that uses the Generalized Minimum Residual 
/// (GMRES) method with the given Krylov subspace dimension.
/// \param [in] comm The communicator on which the solver operates.
/// \param [in] dimension The dimension of the associated Krylov subspace.
/// \memberof hypre_factory
linear_solver_t* hypre_factory_gmres_linear_solver(hypre_factory_t* factory,
                                                   MPI_Comm comm,
                                                   int dimension);
                                             
/// Constructs a linear solver that uses the stabilized bi-conjugate gradient
/// method.
/// \param [in] comm The communicator on which the solver operates.
/// \memberof hypre_factory
linear_solver_t* hypre_factory_bicgstab_solver(hypre_factory_t* factory,
                                               MPI_Comm comm);

/// Constructs a simple diagonal scaling preconditioner for a linear solver.
/// function returns NULL.
/// \memberof hypre_factory
preconditioner_t* hypre_factory_diag_scaling_preconditioner(hypre_factory_t* factory);

///@}

#endif

