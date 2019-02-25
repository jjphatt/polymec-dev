// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/preconditioner.h"

struct preconditioner_t 
{
  char* name;
  void* context;
  preconditioner_side_t side;
  preconditioner_vtable vtable;
};

preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable,
                                     preconditioner_side_t sidedness)
{
  preconditioner_t* pc = polymec_malloc(sizeof(preconditioner_t));
  pc->name = string_dup(name);
  pc->context = context;
  pc->vtable = vtable;
  pc->side = sidedness;
  return pc;
}

void preconditioner_free(preconditioner_t* pc)
{
  if ((pc->vtable.dtor != NULL) && (pc->context != NULL))
    pc->vtable.dtor(pc->context);
  polymec_free(pc);
}

preconditioner_side_t preconditioner_sidedness(preconditioner_t* pc)
{
  return pc->side;
}

int preconditioner_set_up(preconditioner_t* pc)
{
  // FIXME
  return 0;
}

int preconditioner_solve(preconditioner_t* pc, nvector_t* z, nvector_t* r)
{
  // FIXME
  return 0;
}

