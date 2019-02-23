// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_nvector.h"
#include "solvers/nvector.h"

struct nvector_t 
{
  N_Vector v;
  bool owns;
};

nvector_t* nvector_from_N_Vector(N_Vector sundials_nvector,
                                 bool assume_ownership)
{
  nvector_t* v = polymec_malloc(sizeof(nvector_t));
  v->v = sundials_nvector;
  v->owns = assume_ownership;
  return v;
}

void nvector_free(nvector_t* v)
{
  if (v->owns)
    N_VDestroy(v->v);
  polymec_free(v);
}

N_Vector nvector_as_N_Vector(nvector_t* v)
{
  return v->v;
}

