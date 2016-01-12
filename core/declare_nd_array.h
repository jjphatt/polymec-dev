// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DECLARE_ND_ARRAY_H
#define POLYMEC_DECLARE_ND_ARRAY_H

#include "core/polymec.h"

// These macros are used to interpret existing 1D "storage" arrays as 
// multi-dimensional arrays using C99's multidimensional array machinery.
// They take the following arguments.
//   type      - The (primitive) data type of the array.
//   array_var - A multidimensional array that refers to data in storage.
//   storage   - A "flat" array with enough storage to store all the data in 
//               the array.
//   dimX      - The extent of the Xth dimension of the array.

// DECLARE_2D_ARRAY(type, array_var, storage, dim1, dim2)
// DECLARE_3D_ARRAY(type, array_var, storage, dim1, dim2, dim3)
// DECLARE_4D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4)
// DECLARE_5D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5)
// DECLARE_6D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6)
// DECLARE_7D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6, dim7)
// 
// We go to 7D just to make Fortran feel less special. :-)
//
#define DECLARE_2D_ARRAY(type, array_var, storage, dim1, dim2) \
type (*array_var)[dim2] = (void*)storage
#define DECLARE_3D_ARRAY(type, array_var, storage, dim1, dim2, dim3) \
type (*array_var)[dim2][dim3] = (void*)storage
#define DECLARE_4D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4) \
type (*array_var)[dim2][dim3][dim4] = (void*)storage
#define DECLARE_5D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5) \
type (*array_var)[dim2][dim3][dim4][dim5] = (void*)storage
#define DECLARE_6D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6] = (void*)storage
#define DECLARE_7D_ARRAY(type, array_var, storage, dim1, dim2, dim3, dim4, dim5, dim6, dim7) \
type (*array_var)[dim2][dim3][dim4][dim5][dim6][dim7] = (void*)storage

#endif
