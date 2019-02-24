// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MATRIX_H
#define POLYMEC_MATRIX_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class matrix
/// This class represents a real-valued matrix that appears in linear equations, 
/// like A in Ax = b.
///
/// Matrices have several implementations corresponding to different linear 
/// solvers, so they're typically created with factory methods that belong to 
/// those solvers.
///
/// Polymec's matrix class is actually a thin wrapper around the Sundials 
/// SUNMatrix type. We use this wrapper to separate concerns.
typedef struct matrix_t matrix_t;

typedef struct nvector_t nvector_t;

/// Here's the Sundials SUNMatrix type.
typedef struct _generic_SUNMatrix* SUNMatrix;

/// This virtual table defines the behavior for a matrix. Because the matrix
/// class is implemented using the SUNMatrix type, all matrix methods use 
/// SUNMatrix objects instead of raw context pointers. All methods
/// with int return values return 0 on success and nonzero on failure.
typedef struct
{
  /// Clones the matrix's data.
  SUNMatrix (*clone)(SUNMatrix);
  /// Destroys the matrix's data.
  void (*destroy)(SUNMatrix);
  /// Fills the matrix with zeros.
  int (*zero)(SUNMatrix A);
  /// Copies data from A to B.
  int (*copy)(SUNMatrix A, SUNMatrix B);
  /// A <- s * A + B
  int (*scaleadd)(real_t s, SUNMatrix A, SUNMatrix B);
  /// A <- I + s*A
  int (*scaleaddi)(real_t s, SUNMatrix A);
  /// y <- A * x
  int (*matvec)(SUNMatrix A, N_Vector x, N_Vector y);
  /// Computes storage requirements for A in real words (lrw) and integer words (liw) (optional).
  int (*space)(SUNMatrix A, long int* lrw, long int* liw);
  /// Returns the block sizes of the given row in A (optional).
  int (*get_block_sizes)(SUNMatrix A, int num_block_rows, index_t* block_rows, int* block_sizes);
  /// Sets values of A in the specified rows and columns.
  int (*set_values)(SUNMatrix A, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  /// Adds to values of A in the specified rows and columns.
  int (*add_values)(SUNMatrix A, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  /// Retrieves values of A in the specified rows and columns.
  int (*get_values)(SUNMatrix A, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  /// Sets values of specific blocks in A in the specified rows and columns (optional).
  int (*set_blocks)(SUNMatrix A, int num_block, index_t* block_rows, index_t* block_columns, real_t* block_values);
  /// Adds to values of specific blocks in A in the specified rows and columns (optional).
  int (*add_blocks)(SUNMatrix A, int num_block, index_t* block_rows, index_t* block_columns, real_t* block_values);
  /// Retrieves values of specific blocks in A in the specified rows and columns (optional).
  int (*get_blocks)(SUNMatrix A, int num_block, index_t* block_rows, index_t* block_columns, real_t* block_values);
  /// Performs matrix assembly, committing any values to the matrix (optional).
  int (*assemble)(SUNMatrix A);
  /// Writes a text representation of A to the given file (optional).
  void (*fprintf)(SUNMatrix A, FILE* stream);
} matrix_vtable;

/// Creates an instance of a custom matrix type on the given communicator. 
/// \param [in] comm The MPI communicator on which the matrix lives.
/// \param [in] context A pointer to data for the matrix.
/// \param [in] vtable A virtual table defining the behavior of the matrix.
/// \param [in] num_local_rows The number of locally-stored matrix rows.
/// \param [in] num_global_rows The number of globally-stored matrix rows.
/// \memberof matrix
matrix_t* matrix_new(MPI_Comm comm, 
                     void* context, 
                     matrix_vtable vtable,
                     int num_local_rows,
                     index_t num_global_rows);

/// Creates a matrix from a pre-defined (serial) SUNMatrix and 
/// assumes ownership of this object.
/// \param [in] A The SUNMatrix to be represented by this matrix. Must not be 
///               a custom matrix type.
/// \param [in] assume_ownership If true, the matrix assumes ownership of the 
///                              Sundials matrix. Otherwise it doesn't.
/// \memberof matrix
matrix_t* matrix_from_SUNMatrix(SUNMatrix A, bool assume_ownership);

/// Frees a matrix.
/// \memberof matrix
void matrix_free(matrix_t* A);

/// Returns the underlying Sundials SUNMatrix.
/// \memberof matrix
SUNMatrix matrix_as_SUNMatrix(matrix_t* A);

/// Returns the communicator on which the matrix is defined.
/// \memberof matrix
MPI_Comm matrix_comm(matrix_t* A);

/// Retrieves the sizes of the given block rows. If the matrix doesn't 
/// use block rows, the block size is 1.
/// \param [in] num_block_rows The number of block rows for which block sizes are 
///                            retrieved.
/// \param [in] block_rows An array of length num_block_rows storing the indices of
///                        the block rows for which sizes are retrieved.
/// \param [out] block_sizes An array of length num_block_rows that stores the block
///                        sizes for the given rows.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_get_block_sizes(matrix_t* A, 
                           int num_block_rows, 
                           index_t* block_rows,
                           int* block_sizes);

/// Creates and returns a deep copy of a matrix.
/// \memberof matrix
matrix_t* matrix_clone(matrix_t* A);

/// Copies the contents of the matrix A to those of copy.
/// \param [out] copy The (allocated) matrix to which A's values are copied.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_copy(matrix_t* A, matrix_t* copy);

/// Returns the number of locally stored rows in the matrix.
/// \memberof matrix
int matrix_num_local_rows(matrix_t* A);

/// Returns the number of globally stored rows in the matrix.
/// \memberof matrix
index_t matrix_num_global_rows(matrix_t* A);

/// Zeros all of the entries in the given matrix.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
/// \collective
int matrix_zero(matrix_t* A);

/// Performs the operation \f$\mathbf{A} \leftarrow c \mathbf{A} + \mathbf{B}\f$.
/// \param [in] c A scale factor for the matrix.
/// \param [in] B The matrix to be added.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
/// \collective
int matrix_scale_add(matrix_t* A, real_t c, matrix_t* B);

/// Performs the operation \f$\mathbf{A} \leftarrow c \mathbf{A} + \mathbf{I}\f$.
/// \param [in] c A scale factor for the matrix.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
/// \collective
int matrix_scale_add_I(matrix_t* A, real_t c);

/// Computes the matrix-vector product \f$\mathbf{A}*\mathbf{x}\f$, storing the 
/// result in the vector \f$\mathbf{y}\f$. 
/// \param [in] x The vector to be multiplied by the matrix.
/// \param [out] y The vector that stores the matrix-vector product.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_mat_vec(matrix_t* A, nvector_t* x, nvector_t* y);

/// Sets the values of the elements in the matrix identified by the given 
/// _globally-indexed_ rows and columns. The rows and columns being set must 
/// exist on the local process. 
/// \param num_rows [in] The number of rows in which values are set.
/// \param num_columns [in] An array of length num_rows containing the number of columns
///                         whose values are set for each row.
/// \param rows [in] An array of length num_rows indentifying the indices of the rows
///                  whose values are to be set. These rows must be locally stored.
/// \param columns [in] An array of length num_rows identifying indices of columns
///                     whose values are to be set, stored in row-major order.
/// \param values [in] An array of values to be stored, indexed the same way as the columns
///                    array.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_set_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values);
                              
/// Adds the given values of the elements to those in the matrix.
/// The values are identified by the given _globally-indexed_ rows and columns.
/// \param num_rows [in] The number of rows in which values are added to.
/// \param num_columns [in] An array of length num_rows containing the number of columns
///                         whose values are added to for each row.
/// \param rows [in] An array of length num_rows indentifying the indices of the rows
///                  whose values are to be added to. These rows must be locally stored.
/// \param columns [in] An array of length num_rows identifying indices of columns
///                     whose values are to be added to, stored in row-major order.
/// \param values [in] An array of values to add to existing columns, indexed the same 
///                    way as the columns array.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_add_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values);
                              
/// Retrieves the values of the elements in the matrix identified by the 
/// given *globally-indexed* rows and columns, storing them in the values array.
/// \param num_rows [in] The number of rows in which values are retrieved.
/// \param num_columns [in] An array of length num_rows containing the number of columns
///                         whose values are retrieved for each row.
/// \param rows [in] An array of length num_rows indentifying the indices of the rows
///                  whose values are to be retrieved. These rows must be locally stored.
/// \param columns [in] An array of length num_rows identifying indices of columns
///                     whose values are to be retrieved, stored in row-major order.
/// \param values [in] An array that stores retrieved column values, indexed the same 
///                    way as the columns array.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_get_values(matrix_t* A,
                      int num_rows,
                      int* num_columns,
                      index_t* rows, 
                      index_t* columns,
                      real_t* values);

/// Sets the values of the blocks in the matrix identified by the given 
/// *globally-indexed* block rows and columns. The block rows and columns being 
/// set must be locally stored.
/// \param num_blocks [in] The number of blocks in which values are set.
/// \param block_rows [in] An array of length num_blocks indentifying the indices 
///                        of the block rows whose values are set.
/// \param block_columns [in] An array of length num_blocks identifying indices of 
///                           of the block columns whose values are set, stored in 
///                           row-major order.
/// \param block_values [in] An array of length block_size*block_size*num_blocks 
///                          containing values to set within the blocks. This array
///                          stores blocks consecutively, each in row-major order.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_set_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values);
                              
/// Adds in the values of the blocks in the matrix identified by the given 
/// *globally-indexed* block rows and columns. The block rows and columns being 
/// set must be locally stored.
/// \param num_blocks [in] The number of blocks in which values are added to.
/// \param block_rows [in] An array of length num_blocks indentifying the indices 
///                        of the block rows whose values are added to.
/// \param block_columns [in] An array of length num_blocks identifying indices of 
///                           of the block columns whose values are added to, stored in 
///                           row-major order.
/// \param block_values [in] An array of length block_size*block_size*num_blocks 
///                          containing values to add to within the blocks. This array
///                          stores blocks consecutively, each in row-major order.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_add_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values);
                              
/// Retrieves the values of the blocks in the matrix identified by the given 
/// *globally-indexed* block rows and columns. The block rows and columns being 
/// set must be locally stored.
/// \param num_blocks [in] The number of blocks in which values are retrieved.
/// \param block_rows [in] An array of length num_blocks indentifying the indices 
///                        of the block rows whose values are retrieved.
/// \param block_columns [in] An array of length num_blocks identifying indices of 
///                           of the block columns whose values are retrieved, stored in 
///                           row-major order.
/// \param block_values [in] An array of length block_size*block_size*num_blocks 
///                          containing values to retrieve within the blocks. This array
///                          stores blocks consecutively, each in row-major order.
/// \returns 0 on success, nonzero on failure.
/// \memberof matrix
int matrix_get_blocks(matrix_t* A,
                      int num_blocks,
                      index_t* block_rows, 
                      index_t* block_columns,
                      real_t* block_values);

/// Assembles added/inserted values into the matrix, allowing all the processes
/// a consistent representation of the matrix. This should be called after calls
/// to \ref matrix_set_values, \ref matrix_add_values, \ref matrix_set_blocks, 
/// and/or \ref matrix_add_blocks, and should be placed in between sets and adds.
/// \returns 0 on success, or nonzero on failure.
/// \memberof matrix
/// \collective
int matrix_assemble(matrix_t* A);

/// Writes a text representation of the matrix (or portion stored on the local
/// MPI process) to the given stream.
/// \param [out] stream The file to which the representation is written.
/// \memberof matrix
void matrix_fprintf(matrix_t* A, FILE* stream);

///@}

#endif

