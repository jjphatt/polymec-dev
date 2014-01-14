// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_POLYNOMIAL_FIT_H
#define POLYMEC_POLYNOMIAL_FIT_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/least_squares.h"
#include "core/mesh.h"
#include "core/point_cloud.h"

// This type represents a mechanism for generating least-squares polynomial 
// fits for multi-component quantities on discrete domains.
typedef struct polynomial_fit_t polynomial_fit_t;

// This virtual table defines the behavior of a polynomial fit.
typedef struct 
{
  // Gets the number of neighbors associated with the point at the given index.
  int (*num_neighbors)(void* context, int point_index);

  // Retrieves the indices and values of the polynomial components at points 
  // in the vicinity of the point with the given index. The component values 
  // at the desired point are placed into the array point_values. The indices 
  // of the neighbors are stored in the array neighbor_indices, and the 
  // values at the neighbors are stored in the neighbor_values array in 
  // neighbor-major, component-minor order.
  void (*get_data)(void* context, int point_index, 
                   point_t* point, real_t* point_value,
                   point_t* neighbor_points, real_t* neighbor_values); 

  // Returns the targeted degree of accuracy for the polynomial fit, given 
  // a number of neighbors. This method effectively determines how aggressively
  // a polynomial fit will pursue higher-order approximations to the solution.
  int (*targeted_degree)(void* context, int num_neighbors);

  // Destroys the context pointer.
  void (*dtor)(void* context);
} polynomial_fit_vtable;

// Constructs a generic polynomial fit that uses the provided functions to 
// compute fits on discrete domains.
polynomial_fit_t* polynomial_fit_new(const char* name,
                                     void* context,
                                     polynomial_fit_vtable vtable,
                                     int num_comps,
                                     ls_weight_func_t* weight_function);

// Creates a polynomial fit that fits cell-centered data on a mesh at the 
// given degree.
polynomial_fit_t* cc_fixed_degree_polynomial_fit_new(int num_comps,
                                                     ls_weight_func_t* weight_function,
                                                     mesh_t* mesh,
                                                     real_t* data,
                                                     int degree);

// Creates a polynomial fit that fits cell-centered data on a mesh at a 
// variable degree based on the number of neighbors it finds within a given 
// "depth."
polynomial_fit_t* cc_variable_degree_polynomial_fit_new(int num_comps,
                                                        ls_weight_func_t* weight_function,
                                                        mesh_t* mesh,
                                                        real_t* data,
                                                        int neighbor_search_depth);

// Creates a polynomial fit that fits point cloud data at the given degree.
polynomial_fit_t* point_cloud_fixed_degree_polynomial_fit_new(int num_comps,
                                                              ls_weight_func_t* weight_function,
                                                              point_cloud_t* points,
                                                              point_cloud_neighbor_search_t* search,
                                                              real_t* data,
                                                              int degree);

// Creates a polynomial fit that fits point cloud data at a variable degree 
// based on the number of neighbors it finds within a given "depth."
polynomial_fit_t* point_cloud_variable_degree_polynomial_fit_new(int num_comps,
                                                                 ls_weight_func_t* weight_function,
                                                                 point_cloud_t* points,
                                                                 point_cloud_neighbor_search_t* search,
                                                                 real_t* data,
                                                                 int neighbor_search_depth);

// Destroys the given polynomial fit.
void polynomial_fit_free(polynomial_fit_t* fit);

// Returns the number of components in the quantity represented by the 
// polynomial fit.
int polynomial_fit_num_comps(polynomial_fit_t* fit);

// Computes the polynomial fit at the given indexed point in the discrete 
// domain. This discrete domain can be a mesh, a point_cloud, or another 
// abstract discrete domain, as long as neighbor relationships are defined 
// for points.
void polynomial_fit_compute(polynomial_fit_t* fit, int point_index);

// Returns the degree of the present polynomial fit. Returns -1 if no fit 
// has yet been computed.
int polynomial_fit_degree(polynomial_fit_t* fit);

// Evaluate the computed polynomial fit at the given point in space.
void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value);

// Evaluate the given partial derivative of the computed polynomial fit at 
// the given point in space.
void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv);

#endif
