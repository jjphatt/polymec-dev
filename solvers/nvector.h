// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NVECTOR_H
#define POLYMEC_NVECTOR_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class nvector
/// An n-vector represents an n-dimensional vector that appears in linear 
/// equations like the x and b in Ax = b.
///
/// N-vectors have several implementations corresponding to different linear 
/// solvers, so they're typically created with factory methods that belong to 
/// those solvers.
///
/// Polymec's n-vector class is actually a thin wrapper around the Sundials 
/// N_Vector type. We use this wrapper to separate concerns.
typedef struct nvector_t nvector_t;

// Here's the Sundials N_Vector type.
typedef struct _generic_N_Vector* N_Vector;

/// This virtual table defines the behavior for an n-vector. Because the nvector
/// class is implemented using the Sundials N_Vector type, all vector methods use 
/// N_Vector objects instead of raw context pointers. All methods with int 
/// return values return 0 on success and nonzero on failure.
typedef struct
{
  /// Clones the n-vector's data.
  N_Vector (*nvclone)(N_Vector);
  /// Creates an empty clone of the vector's data.
  N_Vector (*nvcloneempty)(N_Vector);
  /// Destroy's the n-vector's data.
  void (*nvdestroy)(N_Vector);
  /// Computes storage requirements for v in real words (lrw) and integer words (liw) (optional).
  void (*nvspace)(N_Vector v, int64_t* lrw, int64_t* liw);
  /// Performs z <- a*x + b*y
  void (*nvlinearsum)(real_t a, N_Vector x, real_t b, N_Vector y, N_Vector z);
  /// Sets all components of v to the constant c.
  void (*nvconst)(real_t c, N_Vector v);
  /// Performs z <- x .* y (componentwise product).
  void (*nvprod)(N_Vector x , N_Vector y, N_Vector z);
  /// Performs z <- x ./ y (componentwise product).
  void (*nvdiv)(N_Vector x, N_Vector y, N_Vector z);
  /// Performs z <- c * x
  void (*nvscale)(real_t c, N_Vector x, N_Vector z);
  /// Performs z[i] <- abs(x[i]) for all i.
  void (*nvabs)(N_Vector x, N_Vector z);
  /// Performs z[i] <- 1.0 / x[i] for all i.
  void (*nvinv)(N_Vector x, N_Vector z);
  /// Performs z[i] <- x[i] + b for all i.
  void (*nvaddconst)(N_Vector x, real_t b, N_Vector z);
  /// Computes the dot product of x and y.
  real_t (*nvdotprod)(N_Vector x, N_Vector y);
  /// Computes the max norm of x.
  real_t (*nvmaxnorm)(N_Vector x);
  /// Computes the weighted root-mean squared (WRMS) norm of x with weights in w.
  real_t (*nvwrmsnorm)(N_Vector x, N_Vector w);
  /// Computes the WRMS norm summed over all i for which mask[i] > 0.
  real_t (*nvwrmsnormmask)(N_Vector x, N_Vector w, N_Vector mask);
  /// Returns the smallest component of x.
  real_t (*nvmin)(N_Vector x);
  /// Returns the weighted l2 norm of x with weights in w.
  real_t (*nvwl2norm)(N_Vector x, N_Vector w);
  /// Returns the l1 norm of x.
  real_t (*nvl1norm)(N_Vector);
  /// Performs z[i] = (|x[i]| >= c) ? 1 : 0.
  void (*nvcompare)(real_t c, N_Vector x, N_Vector z);
  /// nvinv with testing for zeros in x. Returns true if all x values are nonzero, false otherwise.
  int (*nvinvtest)(N_Vector x, N_Vector z);
  /// constraint masking. See Sundials documentation.
  int (*nvconstrmask)(N_Vector, N_Vector, N_Vector);
  /// Returns the minimum componentwise quotient of num and denom.
  real_t (*nvminquotient)(N_Vector num, N_Vector denom);
  /// Computes the linear combination of n scalar-vector products {c[i]*x[i]}, storing the result in z.
  int (*nvlinearcombination)(int n, real_t* c, N_Vector* x, N_Vector z);
  /// Scales and adds y to n vectors {x[i]}, storing the results in z. See Sundial documentation.
  int (*nvscaleaddmulti)(int, real_t*, N_Vector, N_Vector*, N_Vector*);
  /// Computes the dot product of x with n vectors y, storing the products in d.
  int (*nvdotprodmulti)(int n, N_Vector x, N_Vector* y, real_t* d);
  /// Computes the linear sum of vector arrays x and y, storing the results in the array z.
  int (*nvlinearsumvectorarray)(int n, real_t a, N_Vector* x, real_t b, N_Vector* y, N_Vector* z);
  /// Scales n vectors x with different constants c, storing the results in z.
  int (*nvscalevectorarray)(int n, real_t* c, N_Vector* x, N_Vector* z);
  /// Sets all components in the array of n vectors x to c.
  int (*nvconstvectorarray)(int n, real_t c, N_Vector* x);
  /// Computes the WRMS norm of n vectors x using weights w, storing the result in m.
  int (*nvwrmsnormvectorarray)(int n, N_Vector* x, N_Vector* w, real_t* m);
  /// Computes WRMS norms of n nvectors x with weights w, subject to masking.
  int (*nvwrmsnormmaskvectorarray)(int n, N_Vector* x, N_Vector* w, N_Vector mask, real_t* m);
  /// See Sundials documentation for this one.
  int (*nvscaleaddmultivectorarray)(int, int, real_t*, N_Vector*, N_Vector**, N_Vector**);
  /// See Sundials documentation for this one, too.
  int (*nvlinearcombinationvectorarray)(int, int, real_t*, N_Vector**, N_Vector*);

  /// Sets values within a vector.
  void (*set_values)(N_Vector v, int num_values, index_t* indices, real_t* values);
  /// Adds values into a vector.
  void (*add_values)(N_Vector v, int num_values, index_t* indices, real_t* values);
  /// Retrieves values from a vector.
  void (*get_values)(N_Vector v, int num_values, index_t* indices, real_t* values);
  /// Assembles a vector, committing all set/added values (optional).
  void (*assemble)(N_Vector v);
  /// Prints a vector to a file (optional.
  void (*fprintf)(N_Vector v, FILE* stream);
} nvector_vtable;

/// Creates an instance of a custom n-vector type on the given communicator. 
/// \param [in] comm The MPI communicator on which the n-vector lives.
/// \param [in] context A pointer to data for the n-vector.
/// \param [in] vtable A virtual table defining the behavior of the n-vector.
/// \param [in] local_size The number of locally-stored n-vector rows.
/// \param [in] global_size The number of globally-stored n-vector rows.
/// \memberof nvector
nvector_t* nvector_new(MPI_Comm comm, 
                       void* context, 
                       nvector_vtable vtable,
                       int local_size,
                       index_t global_size);

/// Creates an instance of an n-vector from the given (serial) pre-fab 
/// Sundials N_Vector. 
/// \param [in] sundials_nvector A Sundials N_Vector from which to create this 
///                              n-vector. Must not be a custom n-vector type.
/// \param [in] assume_ownership If true, the nvector assumes ownership of the Sundials 
///                              vector. Otherwise it doesn't.
/// \memberof nvector
nvector_t* nvector_from_N_Vector(N_Vector sundials_nvector,
                                 bool assume_ownership);

/// Frees an n-vector.
/// \memberof nvector
void nvector_free(nvector_t* v);

/// Returns the underlying Sundials N_Vector.
/// \memberof nvector
N_Vector nvector_as_N_Vector(nvector_t* v);

/// Creates and returns a deep copy of a vector.
/// \memberof nvector
nvector_t* nvector_clone(nvector_t* v);

/// Copies the contents of the vector v to those of copy.
/// \memberof nvector
void nvector_copy(nvector_t* v, nvector_t* copy);

/// Returns the locally-stored size (dimension) of the vector.
/// \memberof nvector
int nvector_local_size(nvector_t* v);

/// Returns the global size (dimension) of the vector.
/// \memberof nvector
index_t nvector_global_size(nvector_t* v);

/// Sets all of the components of v to the given constant c.
/// \memberof nvector
/// \collective
void nvector_set(nvector_t* v, real_t c);

/// Scales the vector by the given factor c, storing the result in w.
/// \param [out] w The scaled vector.
/// \memberof nvector
/// \collective
void nvector_scale(nvector_t* v, real_t c, nvector_t* w);

/// Sets the values of the elements in the vector identified by the given 
/// _globally-indexed_ indices.
/// \param [in] num_values The number of (row) values to set in the vector.
/// \param [in] indices An array of length num_values storing row indices of v to be set.
/// \param [in] values An array of length num_values storing the values.
/// \memberof nvector
void nvector_set_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values);
                              
/// Adds the given values of the elements to those in the vector.
/// These values are identified by the given *globally-indexed* indices.
/// \param [in] num_values The number of (row) values to add into the vector.
/// \param [in] indices An array of length num_values storing row indices of v to be added.
/// \param [in] values An array of length num_values storing the values.
/// \memberof nvector
void nvector_add_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values);

/// Retrieves the values of the elements in the vector identified by the 
/// given _global_ indices, storing them in the values array. The values 
/// must exist on the local process.
/// \param [in] num_values The number of (row) values to retrieve from the vector.
/// \param [in] indices An array of length num_values storing row indices of v to be retrieved.
/// \param [in] values An array of length num_values that stores the retrieved values.
/// \memberof nvector
void nvector_get_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values);

/// Assembles added/inserted values into the vector, allowing all the processes
/// a consistent representation of the vector. This should be called after calls
/// to \ref nvector_set_values and \ref nvector_add_values, and should be placed in between 
/// sets and adds. 
/// \memberof nvector
/// \collective
void nvector_assemble(nvector_t* A);

/// Computes and returns the dot product of the vector v with the vector w.
/// This is collective, and must be called by all MPI processes.
/// \param [in] w The vector dotted by v.
/// \memberof nvector
real_t nvector_dot(nvector_t* v, nvector_t* w);

/// Computes the weighted l2-norm for this vector, using the weights
/// in the vector w. If w is NULL, the 2-norm is returned.
/// The weighted l2-norm is \f$\sqrt{\sum_i (w_i\cdot v_i)^2}\f$.
/// \param [in] w A vector of weights.
/// \memberof nvector
/// \collective
real_t nvector_wl2_norm(nvector_t* v, nvector_t* w);

/// Computes and returns a weighted root-mean-squared (WRMS) norm for 
/// this vector, using the weights in the given vector w. w must be non-NULL.
/// The WRMS norm is \f$\frac{\sqrt{\sum_i, (w_i\cdoti v_i)^2}}{N}.
/// \param [in] w A vector of weights.
/// \memberof nvector
/// \collective
real_t nvector_wrms_norm(nvector_t* v, nvector_t* w);
 
/// Writes a text representation of the vector (or portion stored on the local
/// MPI process) to the given stream.
/// \param [out] stream The file to which the vector is written.
/// \memberof nvector
void nvector_fprintf(nvector_t* v, FILE* stream);

///@}

#endif

