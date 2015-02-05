// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CPR_DIFFERENCER_H
#define POLYMEC_CPR_DIFFERENCER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/local_matrix.h"

// The cpr_differencer is an object that computes a Jacobian matrix for a 
// function F. The sparsity of the matrix is given by a graph, and the 
// differencer computes the matrix using the method of Curtis, Powell, and 
// Reed, in which the graph is colored so that the Jacobian can be computed 
// using finite differencing, with a minimum number of calls to F.
typedef struct cpr_differencer_t cpr_differencer_t;

// Creates a new Curtis-Powell-Reed differencer, associated with the given 
// function F *OR* F_dae and a graph, which represents the sparsity of the 
// Jacobian matrix. The numbers of rows local to this process and "remote"
// (stored on other processes) are given to allow F to exchange data.
// The sparsity graph is consumed by this constructor.
cpr_differencer_t* cpr_differencer_new(MPI_Comm comm,
                                       void* F_context,
                                       int (*F)(void* context, real_t, real_t* x, real_t* Fval),
                                       int (*F_dae)(void* context, real_t, real_t* x, real_t* xdot, real_t* Fval),
                                       void (*F_dtor)(void* context),
                                       adj_graph_t* sparsity,
                                       int num_local_rows,
                                       int num_remote_rows);

// Frees the differencer.
void cpr_differencer_free(cpr_differencer_t* diff);

// Uses the differencer to compute a generalized Jacobian matrix of the form
// J = alpha * I + beta * dF/dx + gamma * dF/d(xdot). NULL can be passed for 
// xdot if gamma == 0.0.
void cpr_differencer_compute(cpr_differencer_t* diff, 
                             real_t alpha, real_t beta, real_t gamma,  
                             real_t t, real_t* x, real_t* xdot,
                             local_matrix_t* matrix);

#endif