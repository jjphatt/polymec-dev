// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/polymec.h"
#include "solvers/async_solver.h"

//------------------------------------------------------------------------
//                    Simple 1D harmonic oscillator
//------------------------------------------------------------------------
typedef struct 
{
  int N;
  real_t *k, *x, *v, *a;
} sho_t;

static int sho_num_entities(void* context)
{
  sho_t* sho = context;
  return sho->N;
}

static real_t sho_entity_dt(void* context, int i, real_t t)
{
  sho_t* sho = context;
  return MIN(1.0/ABS(sho->v[i] + 1e-8), sqrt(1e-5 / ABS(sho->k[i] * sho->x[i])));
}

static void sho_active_step(void* context, 
                            int* entities, size_t num_entities,
                            real_t t, real_t dt)
{
  sho_t* sho = context;
  for (int i = 0; i < sho->N; ++i)
  {
    sho->a[i] = -sho->k[i] * sho->x[i];
    sho->v[i] += sho->a[i] * dt;
    sho->x[i] += sho->v[i] * dt;
  }
}

static void sho_inactive_step(void* context, 
                              int* entities, size_t num_entities,
                              real_t t, real_t dt)
{
  sho_t* sho = context;
  for (int i = 0; i < sho->N; ++i)
  {
    sho->v[i] += sho->a[i] * dt;
    sho->x[i] += sho->v[i] * dt;
  }
}

static void sho_dtor(void* context)
{
  sho_t* sho = context;
  polymec_free(sho->k);
  polymec_free(sho->x);
  polymec_free(sho->v);
  polymec_free(sho->a);
  polymec_free(sho);
}

static void test_sho(void** state, int N)
{
  // Create an async_solver for a set of harmonic oscillators.
  sho_t* sho = polymec_malloc(sizeof(sho_t));
  sho->N = N;
  sho->k = polymec_malloc(sizeof(real_t) * N);
  sho->x = polymec_malloc(sizeof(real_t) * N);
  sho->v = polymec_malloc(sizeof(real_t) * N);
  sho->a = polymec_malloc(sizeof(real_t) * N);
  async_solver_vtable vtable = {.num_entities = sho_num_entities,
                                .entity_dt = sho_entity_dt,
                                .active_step = sho_active_step,
                                .inactive_step = sho_inactive_step,
                                .dtor = sho_dtor};
  async_solver_t* solver = async_solver_new(sho, vtable);

  // Set initial conditions and spring constants.
  real_t k_min = 1.0;
  for (int i = 0; i < N; ++i)
  {
    sho->k[i] = (i+1)*k_min;
    sho->x[i] = 0.5;
    sho->v[i] = 0.0;
    sho->a[i] = 0.0;
  }

  // Integrate over the maximum period of all oscillators. Since all the 
  // periods are multiples of each other, they should all end up in their 
  // original positions (+ error).
  real_t omega = sqrt(k_min);
  real_t T = 2.0 * M_PI / omega;
  log_debug("Integrating from t = 0 to %g", T);
  real_t t = 0.0;
  while (t < T)
  {
    real_t max_dt = async_solver_max_dt(solver, t);
    real_t dt = async_solver_step(solver, t, max_dt);
    t += dt;
  }
  for (int i = 0; i < N; ++i)
  {
    printf("%d: x = %g, v = %g\n", i, sho->x[i], sho->v[i]);
    assert_true(reals_nearly_equal(sho->x[i], 0.5, 0.05));
    assert_true(reals_nearly_equal(sho->v[i], 0.0, 1e-3));
  }

  // Clean up.
  async_solver_free(solver);
}

static void test_1_sho(void** state)
{
  test_sho(state, 1);
}

static void test_2_shos(void** state)
{
  test_sho(state, 2);
}

static void test_10_shos(void** state)
{
  test_sho(state, 10);
}

static void test_100_shos(void** state)
{
  test_sho(state, 100);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_1_sho),
    cmocka_unit_test(test_2_shos),
    cmocka_unit_test(test_10_shos),
    cmocka_unit_test(test_100_shos)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
