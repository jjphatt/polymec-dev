// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DIFFERENCE_SD_FUNC_H
#define POLYMEC_DIFFERENCE_SD_FUNC_H

#include "geometry/sd_func.h"

// This signed distance function represents the difference of two 
// given surfaces represented by signed distance functions.
sd_func_t* difference_sd_func_new(sd_func_t* surface1, sd_func_t* surface2);

#endif

