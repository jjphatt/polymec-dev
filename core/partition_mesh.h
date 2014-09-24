// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#ifndef POLYMEC_PARTITION_MESH_H
#define POLYMEC_PARTITION_MESH_H

#include "core/mesh.h"

// This function partitions the given mesh on rank 0 with the given cell-centered 
// load weights, distributing the cells to parallel domains on the given 
// communicator to balance the load. It creates and returns an exchanger 
// object that can be used to distribute data from the rank 0 to the partition. 
// The mesh on rank 0 (as well as any non-NULL mesh on rank != 0) is consumed. 
// In each case, the mesh is replaced with a partitioned mesh.
exchanger_t* partition_mesh(mesh_t** mesh, MPI_Comm comm, int* weights, real_t imbalance_tol);

// This function repartitions the given mesh with the given load weights, 
// alloting the cells to parallel domains to balance their load.
// It creates and returns an exchanger object that can be used to migrate 
// data from the old partition to the new. The mesh is consumed and replaced
// with a repartitioned mesh.
exchanger_t* repartition_mesh(mesh_t** mesh, int* weights, real_t imbalance_tol);

#endif
