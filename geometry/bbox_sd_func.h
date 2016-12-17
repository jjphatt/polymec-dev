// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BBOX_SD_FUNC_H
#define POLYMEC_BBOX_SD_FUNC_H

#include "geometry/sd_func.h"

// This creates a signed distance function equivalent to the given 
// bounding box.
sd_func_t* bbox_sd_func_new(bbox_t* bounding_box);

#endif
