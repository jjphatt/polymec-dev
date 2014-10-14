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

#include <gc/gc.h>
#include "core/rng.h"

struct rng_t 
{
  char* name;
  void* context;
  rng_vtable vtable;
  uint32_t min, max;
  bool has_global_state, is_thread_safe;
};

static void rng_free(void* ctx, void* dummy)
{
  rng_t* rng = ctx;
  string_free(rng->name);
  if ((rng->vtable.dtor != NULL) && (rng->context != NULL))
    rng->vtable.dtor(rng->context);
}

rng_t* rng_new(const char* name, void* context, 
               uint32_t min, uint32_t max, rng_vtable vtable,
               bool has_global_state, bool is_thread_safe)
{
  ASSERT(min < max);
  ASSERT(vtable.get != NULL);
  rng_t* rng = GC_MALLOC(sizeof(rng_t));
  rng->name = string_dup(name);
  rng->context = context;
  rng->min = min;
  rng->max = max;
  rng->vtable = vtable;
  rng->has_global_state = has_global_state;
  rng->is_thread_safe = is_thread_safe;
  GC_register_finalizer(rng, rng_free, rng, NULL, NULL);
  return rng;
}

const char* rng_name(rng_t* rng)
{
  return (const char*)rng->name;
}

void rng_set_seed(rng_t* rng, uint32_t seed)
{
  if (rng->vtable.set_seed != NULL)
    rng->vtable.set_seed(rng->context, seed);
}

uint32_t rng_min(rng_t* rng)
{
  return rng->min;
}

uint32_t rng_max(rng_t* rng)
{
  return rng->max;
}

uint32_t rng_get(rng_t* rng)
{
  return rng->vtable.get(rng->context);
}

real_t rng_uniform(rng_t* rng)
{
  if (rng->vtable.uniform != NULL)
    return rng->vtable.uniform(rng->context);
  else
    return (real_t)(1.0 * rng_get(rng) / (rng->max + 1.0));
}

real_t rng_uniform_positive(rng_t* rng)
{
  if (rng->vtable.uniform_positive != NULL)
    return rng->vtable.uniform_positive(rng->context);
  else
  {
    while (true)
    {
      real_t val = rng_uniform(rng);
      if (val != 0.0)
        return val;
    }
  }
}

uint32_t rng_uniform_int(rng_t* rng, uint32_t n)
{
  if (rng->vtable.uniform_int != NULL)
    return rng->vtable.uniform_int(rng->context, n);
  else
  {
    uint32_t reduction_factor = rng->max / n;
    return rng_get(rng) / reduction_factor;
  }
}

#ifdef _BSD_SOURCE

// POSIX random() generator.

static const uint32_t posix_max = (uint32_t)(-1);

static void posix_set_seed(void* context, uint32_t seed)
{
  char* state = context;
  setstate((const char*)state);
  srandom(seed);
}

static uint32_t posix_get(void* context)
{
  return (uint32_t)random();
}

static void posix_free(void* context)
{
  polymec_free(context);
}

rng_t* posix_rng_new(size_t state_size)
{
  static const int num_bytes = 256;
  char* state = polymec_malloc(sizeof(char) * num_bytes);
  rng_vtable vtable = {.set_seed = posix_set_seed,
                       .get = posix_get,
                       .dtor = posix_free};
  initstate(random(), state, num_bytes);
  return rng_new("posix RNG", state, 0, posix_max, vtable, true, false);
}
#endif

#if APPLE 

// ARC4 random generator.

static const uint32_t arc4_max = (uint32_t)(-1);

static uint32_t arc4_get(void* context)
{
  return arc4random();
}

static uint32_t arc4_uniform_int(void* context, uint32_t n)
{
  return arc4random_uniform(n);
}

static rng_t* arc4_rng_new()
{
  rng_vtable vtable = {.get = arc4_get,
                       .uniform_int = arc4_uniform_int};
  return rng_new("arc4 RNG", NULL, 0, arc4_max, vtable, true, true);
}
#endif

// Standard C rand() generator -- always available.

static const uint32_t rand_max = (uint32_t)(-1);

static void rand_set_seed(void* context, uint32_t seed)
{
  srand(seed);
}

static uint32_t rand_get(void* context)
{
  return (uint32_t)rand();
}

rng_t* rand_rng_new()
{
  rng_vtable vtable = {.set_seed = rand_set_seed,
                       .get = rand_get};
  return rng_new("rand (standard C) RNG", NULL, 0, rand_max, vtable,
                 true, false);
}

rng_t* host_rng_new()
{
#if APPLE 
  return arc4_rng_new();
#elif defined(_BSD_SOURCE)
  return posix_rng_new();
#else
  return rand_rng_new();
#endif
}
