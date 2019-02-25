// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_nvector.h"

#include "core/timer.h"
#include "solvers/nvector.h"

struct nvector_t 
{
  MPI_Comm comm;
  int local_size;
  index_t global_size;
  N_Vector v;
  nvector_vtable vtable;
  bool owns;
};

static N_Vector_ID getvectorid(N_Vector v)
{
  return SUNDIALS_NVEC_CUSTOM;
}

nvector_t* nvector_new(MPI_Comm comm, 
                       void* context, 
                       nvector_vtable vtable,
                       int local_size,
                       index_t global_size)
{
  ASSERT(local_size > 0);
  ASSERT(global_size >= (index_t)local_size);
  ASSERT(vtable.nvclone != NULL);
  ASSERT(vtable.nvcloneempty != NULL);
  ASSERT(vtable.nvdestroy != NULL);
  ASSERT(vtable.nvlinearsum != NULL);
  ASSERT(vtable.nvconst != NULL);
  ASSERT(vtable.nvprod != NULL);
  ASSERT(vtable.nvdiv != NULL);
  ASSERT(vtable.nvscale != NULL);
  ASSERT(vtable.nvabs != NULL);
  ASSERT(vtable.nvinv != NULL);
  ASSERT(vtable.nvaddconst != NULL);
  ASSERT(vtable.nvdotprod != NULL);
  ASSERT(vtable.nvmaxnorm != NULL);
  ASSERT(vtable.nvwrmsnorm != NULL);
  ASSERT(vtable.nvwrmsnormmask != NULL);
  ASSERT(vtable.nvmin != NULL);
  ASSERT(vtable.nvwl2norm != NULL);
  ASSERT(vtable.nvl1norm != NULL);
  ASSERT(vtable.nvcompare != NULL);
  ASSERT(vtable.nvinvtest != NULL);
  ASSERT(vtable.nvconstrmask != NULL);
  ASSERT(vtable.nvminquotient != NULL);
  ASSERT(vtable.nvlinearcombination != NULL);
  ASSERT(vtable.nvscaleaddmulti != NULL);
  ASSERT(vtable.nvdotprodmulti != NULL);
  ASSERT(vtable.nvlinearsumvectorarray != NULL);
  ASSERT(vtable.nvconstvectorarray != NULL);
  ASSERT(vtable.nvwrmsnormvectorarray != NULL);
  ASSERT(vtable.nvwrmsnormmaskvectorarray != NULL);
  ASSERT(vtable.nvscaleaddmultivectorarray != NULL);
  ASSERT(vtable.nvlinearcombinationvectorarray != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);

  nvector_t* v = polymec_malloc(sizeof(nvector_t));
  v->comm = comm;
  v->local_size = local_size;
  v->global_size = global_size;
  v->vtable = vtable;
  v->owns = true;

  // Set up the N_Vector.
  v->v = polymec_malloc(sizeof(struct _generic_N_Vector));
  v->v->content = context;
  v->v->ops = polymec_malloc(sizeof(struct _generic_N_Vector_Ops));
  v->v->ops->nvgetvectorid = getvectorid;
  v->v->ops->nvclone = vtable.nvclone;
  v->v->ops->nvcloneempty = vtable.nvcloneempty;
  v->v->ops->nvdestroy = vtable.nvdestroy;
  v->v->ops->nvspace = vtable.nvspace;
  v->v->ops->nvlinearsum = vtable.nvlinearsum;
  v->v->ops->nvconst = vtable.nvconst;
  v->v->ops->nvprod = vtable.nvprod;
  v->v->ops->nvdiv = vtable.nvdiv;
  v->v->ops->nvscale = vtable.nvscale;
  v->v->ops->nvabs = vtable.nvabs;
  v->v->ops->nvinv = vtable.nvinv;
  v->v->ops->nvaddconst = vtable.nvaddconst;
  v->v->ops->nvdotprod = vtable.nvdotprod;
  v->v->ops->nvmaxnorm = vtable.nvmaxnorm;
  v->v->ops->nvwrmsnorm = vtable.nvwrmsnorm;
  v->v->ops->nvwrmsnormmask = vtable.nvwrmsnormmask;
  v->v->ops->nvmin = vtable.nvmin;
  v->v->ops->nvwl2norm = vtable.nvwl2norm;
  v->v->ops->nvl1norm = vtable.nvl1norm;
  v->v->ops->nvcompare = vtable.nvcompare;
  v->v->ops->nvinvtest = vtable.nvinvtest;
  v->v->ops->nvconstrmask = vtable.nvconstrmask;
  v->v->ops->nvminquotient = vtable.nvminquotient;
  v->v->ops->nvlinearcombination = vtable.nvlinearcombination;
  v->v->ops->nvscaleaddmulti = vtable.nvscaleaddmulti;
  v->v->ops->nvdotprodmulti = vtable.nvdotprodmulti;
  v->v->ops->nvlinearsumvectorarray = vtable.nvlinearsumvectorarray;
  v->v->ops->nvconstvectorarray = vtable.nvconstvectorarray;
  v->v->ops->nvwrmsnormvectorarray = vtable.nvwrmsnormvectorarray;
  v->v->ops->nvwrmsnormmaskvectorarray = vtable.nvwrmsnormmaskvectorarray;
  v->v->ops->nvscaleaddmultivectorarray = vtable.nvscaleaddmultivectorarray;
  v->v->ops->nvlinearcombinationvectorarray = vtable.nvlinearcombinationvectorarray;

  return v;
}

nvector_t* nvector_from_N_Vector(N_Vector sundials_nvector,
                                 bool assume_ownership)
{
  ASSERT(N_VGetVectorID(sundials_nvector) != SUNDIALS_NVEC_CUSTOM); 

  nvector_t* v = polymec_malloc(sizeof(nvector_t));
  v->comm = MPI_COMM_SELF;
  // FIXME: Fill in other fields.
  v->v = sundials_nvector;
  v->owns = assume_ownership;
  return v;
}

void nvector_free(nvector_t* v)
{
  if (v->owns)
  {
    N_Vector_ID id = N_VGetVectorID(v->v);
    if (id == SUNDIALS_NVEC_CUSTOM)
      polymec_free(v->v->ops);
    N_VDestroy(v->v);
  }
  polymec_free(v);
}

N_Vector nvector_as_N_Vector(nvector_t* v)
{
  return v->v;
}

nvector_t* nvector_clone(nvector_t* v)
{
  nvector_t* clone = polymec_malloc(sizeof(nvector_t));
  clone->comm = v->comm;
  clone->local_size = v->local_size;
  clone->global_size = v->global_size;
  clone->vtable = v->vtable;
  clone->v = v->vtable.nvclone(v->v);
  return clone;
}

void nvector_copy(nvector_t* v, nvector_t* copy)
{
  v->vtable.nvscale(1.0, v->v, copy->v);
}

int nvector_local_size(nvector_t* v)
{
  return v->local_size;
}

index_t nvector_global_size(nvector_t* v)
{
  return v->global_size;
}

void nvector_set(nvector_t* v, real_t c)
{
  v->vtable.nvconst(c, v->v);
}

void nvector_scale(nvector_t* v, real_t c, nvector_t* w)
{
  v->vtable.nvscale(c, v->v, w->v);
}

void nvector_set_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.set_values(v->v, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}
                              
void nvector_add_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.add_values(v->v, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}

void nvector_get_values(nvector_t* v,
                        int num_values,
                        index_t* indices,
                        real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.get_values(v->v, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}

void nvector_assemble(nvector_t* v)
{
  START_FUNCTION_TIMER();
  if (v->vtable.assemble != NULL)
    v->vtable.assemble(v->v);
  STOP_FUNCTION_TIMER();
}

real_t nvector_dot(nvector_t* v, nvector_t* w)
{
  return v->vtable.nvdotprod(v->v, w->v);
}

real_t nvector_wl2_norm(nvector_t* v, nvector_t* w)
{
  return v->vtable.nvwl2norm(v->v, w->v);
}

real_t nvector_wrms_norm(nvector_t* v, nvector_t* w)
{
  return v->vtable.nvwrmsnorm(v->v, w->v);
}
 
void nvector_fprintf(nvector_t* v, FILE* stream)
{
  if ((v->vtable.fprintf != NULL) && (stream != NULL))
    v->vtable.fprintf(v->v, stream);
}

