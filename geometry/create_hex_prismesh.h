// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_HEX_PRISMESH_H
#define POLYMEC_CREATE_HEX_PRISMESH_H

#include "geometry/prismesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a prismesh consisting of nx x ny x nz
/// hexagonal prism cells spanning the given bounding box. The box can be 
/// periodic in x, y, and/or z. The cells on the xy boundary are cropped to 
/// fit the bounding box unless the box is periodic on that boundary.
prismesh_t* create_hex_prismesh(int nx, int ny, int nz,
                                bbox_t* bbox,
                                bool periodic_in_x,
                                bool periodic_in_y,
                                bool periodic_in_z);

///@}

#endif

