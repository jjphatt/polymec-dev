// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_UNION_H
#define POLYMEC_UNION_H

#include "core/sp_func.h"

// This signed distance function represents the union of the set of 
// given surfaces represented by signed distance functions.
sp_func_t* union_new(sp_func_t** surfaces, int num_surfaces);

// This is a shorthand function that creates the union of two 
// surfaces.
sp_func_t* union_new2(sp_func_t* surface1, sp_func_t* surface2);

#endif

