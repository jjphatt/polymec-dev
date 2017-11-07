// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_RNG_H
#define POLYMEC_RNG_H

#include "core/polymec.h"

// This type represents a random number generator with a specific min, max, 
// and underlying algorithm. Objects of this type are garbage-collected.
typedef struct rng_t rng_t;

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  // Optional method to set the seed to the random number generator. 
  // rng_set_seed() is a no-op if not provided.
  void (*set_seed)(void* context, uint32_t seed);

  // Required method to return a random integer in [min, max].
  uint32_t (*get)(void* context);

  // Optional method to return a random real_t in [0, 1].
  // Defaults to (rng_get(rng) / (max + 1.0)).
  real_t (*uniform)(void* context);

  // Optional method to return a random real in (0, 1].
  // Defaults to rng_uniform(rng), which is called until it returns 
  // a positive value.
  real_t (*uniform_positive)(void* context);

  // Optional method to return a random integer in [0, n-1].
  // Defaults to a scaled implementation using rng_get(rng).
  uint32_t (*uniform_int)(void* context, uint32_t n);

  // Optional method to provide a destructor for the context.
  void (*dtor)(void* context);
} rng_vtable;

// Constructs a random number generator from the given name, context, and 
// virtual table.
rng_t* rng_new(const char* name, void* context, 
               uint32_t min, uint32_t max, rng_vtable vtable,
               bool has_global_state, bool is_thread_safe);

// Returns the name of the random number generator.
const char* rng_name(rng_t* rng);

// Returns true if the random number generator has a global state, false if 
// not. Most "bundled" generators have a global state in the sense that two 
// different rng objects of a given type will generate random numbers using 
// common machinery, and using one will affect the state of the other.
bool rng_has_global_state(rng_t* rng);

// Returns true if the random number generator is thread-safe (re-entrant), 
// false if not. Note that the question of thread-safety is separate from 
// that of a global state (see above).
bool rng_is_thread_safe(rng_t* rng);

// Sets the seed used by the random number generator.
void rng_set_seed(rng_t* rng, uint32_t seed);

// Returns the minimum (positive) random integer that can be generated 
// by this generator.
uint32_t rng_min(rng_t* rng);

// Returns the maximum (positive) random integer that can be generated 
// by this generator.
uint32_t rng_max(rng_t* rng);

// Generates and returns a random integer within the range [min, max], where 
// these bounds are, respectively, the minimum and maximum numbers that can 
// be generated by the random number generator. All integers in this range 
// are equally likely.
uint32_t rng_get(rng_t* rng);

// Generates a random floating point number within the range [0, 1), with 
// uniform probability.
real_t rng_uniform(rng_t* rng);

// Generates a random floating point number within the range (0, 1), with 
// uniform probability.
real_t rng_uniform_positive(rng_t* rng);

// Generates an integer within the range [0, n-1], with 
// uniform probability.
uint32_t rng_uniform_int(rng_t* rng, uint32_t n);

#ifdef _BSD_SOURCE
// On (BSD-like) systems with random() available, this creates a random 
// number with the given number of bytes in the internal state. Optimal sizes 
// are 8, 32, 64, 128, and 256.
rng_t* posix_rng_new(size_t state_size);
#endif

// Creates the standard C random number generator implemented using rand().
// This is available on every platform.
rng_t* rand_rng_new(void);

// Creates the best vanilla random number generator available on this system.
rng_t* host_rng_new(void);

#endif

