// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "solvers/async_solver.h"

struct async_solver_t 
{
  void* context;
  async_solver_vtable vtable;
};

async_solver_t* async_solver_new(void* context,
                                 async_solver_vtable vtable)
{
  // Validate the vtable.
  ASSERT(vtable.num_entities != NULL);
  ASSERT(vtable.entity_dt != NULL);
  ASSERT(vtable.active_step != NULL);
  ASSERT(vtable.inactive_step != NULL);

  // Put everything together.
  async_solver_t* solver = polymec_malloc(sizeof(async_solver_t));
  solver->context = context;
  solver->vtable = vtable;
  return solver;
}

void async_solver_free(async_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  polymec_free(solver);
}

void async_solver_init(async_solver_t* solver, real_t t)
{
  if (solver->vtable.init != NULL)
  {
    START_FUNCTION_TIMER();
    solver->vtable.init(solver->context, t);
    STOP_FUNCTION_TIMER();
  }
}

real_t async_solver_max_dt(async_solver_t* solver, real_t t)
{
  START_FUNCTION_TIMER();
  real_t max_dt = -REAL_MAX;
  int N = solver->vtable.num_entities(solver->context);
  ASSERT(N > 0);
  for (int i = 0; i < N; ++i)
    max_dt = MAX(max_dt, solver->vtable.entity_dt(solver->context, i, t));
  STOP_FUNCTION_TIMER();
  return max_dt;
}

static inline bool bin_is_active(int bin, int substep, int num_bins)
{
  return ((substep % (1 << bin)) == 0);
}

real_t async_solver_step(async_solver_t* solver, real_t t, real_t dt)
{
  START_FUNCTION_TIMER();
  int N = solver->vtable.num_entities(solver->context);
  ASSERT(N > 0);

  // Compute steps for each entity, fitting each one into a bin for the nearest 
  // power of two down from the given time step. Also keep track of the number 
  // of bins we are creating.
  int num_bins = 1;
  int* dt_bins = polymec_malloc(sizeof(int) * N);
  real_t max_dt = -REAL_MAX;
  for (int i = 0; i < N; ++i)
  {
    real_t dt_i = solver->vtable.entity_dt(solver->context, i, t);
    max_dt = MAX(max_dt, dt_i);
    int bin = 0;
    while (dt_i < dt)
    {
      dt_i *= 2.0;
      ++bin;
    }
    dt_bins[i] = bin;
    num_bins = MAX(num_bins, bin+1);
  }

  // Make sure we're not exceeding the given step size.
  max_dt = MIN(max_dt, dt);

  // Now bin 0 holds all entities with the largest time step size. Bin i+1
  // has a time step that is half the size of bin i, and so on.

  // Call our pre-step thingy if we have one.
  if (solver->vtable.pre_step != NULL)
    solver->vtable.pre_step(solver->context, t);

  // If there's only one bin, we do the obvious thing.
  if (num_bins == 1)
  {
    // All entities are active, and are integrated over the entire step.
    int_array_t* all_entities = int_array_new();
    for (int i = 0; i < N; ++i)
      int_array_append(all_entities, i);
    solver->vtable.active_step(solver->context, 
                               all_entities->data, all_entities->size, 
                               t, max_dt);
    log_debug("async_solver_step: Took 1 step.");
    int_array_free(all_entities);
  }

  // Otherwise, the fun begins!
  else
  {
    // The number of substeps in a time step is 2^(num_bins-1). 
    int num_substeps = 1 << (num_bins-1);
    real_t min_dt = max_dt;

    // Make lists of which entities are in which bins.
    int_array_t** entities_for_bin = polymec_malloc(sizeof(int_array_t*) * num_bins);
    for (int bin = 0; bin < num_bins; ++bin)
      entities_for_bin[bin] = int_array_new();
    for (int i = 0; i < N; ++i)
      int_array_append(entities_for_bin[dt_bins[i]], i);

    // Now perform each of the substeps.
    int actual_substeps = 0;
    real_t t_cur = t;
    for (int substep = 0; substep < num_substeps; ++substep, ++actual_substeps)
    {
      real_t substep_dt = max_dt * pow(0.5, num_bins-1);
      min_dt = MIN(min_dt, substep_dt);
      for (int bin = num_bins-1; bin >= 0; --bin) // backwards, from smallest to largest dt
      {
        int_array_t* entities = entities_for_bin[bin];
        if (bin_is_active(bin, substep, num_bins))
        {
          // Integrate active entities.
//          real_t dt_for_bin = max_dt * pow(0.5, bin);
          solver->vtable.active_step(solver->context, 
                                     entities->data, entities->size, 
                                     t_cur, substep_dt);
        }
        else
        {
          // Inactive entities are just along for the ride.
          solver->vtable.inactive_step(solver->context, 
                                       entities->data, entities->size, 
                                       t_cur, substep_dt);
        }
      }

      // Now recompute the allowed time step for each active entity.
      for (int bin = num_bins-1; bin >= 0; --bin) // backwards, from smallest to largest dt
      {
        int_array_t* entities = entities_for_bin[bin];
        real_t dt_for_bin = max_dt * pow(0.5, bin);

        if (bin_is_active(bin, substep, num_bins))
        {
          for (size_t i = 0; i < entities->size; ++i) 
          {
            int e = entities->data[i];
            real_t new_dt = solver->vtable.entity_dt(solver->context, e, t_cur);

            // Did the timestep decrease enough for the entity to move to a 
            // different bin? If so, we do a little surgery.
            if (new_dt < 0.5 * dt_for_bin)
            {
              if (bin == num_bins-1)
              {
                // We're in the bin with the smallest time step already, so we
                // have to add a new bin with an even smaller time step.
                ++num_bins;
                entities_for_bin = polymec_realloc(entities_for_bin, sizeof(int_array_t*) * num_bins);
                entities_for_bin[num_bins-1] = int_array_new();

                // Oh, hey, we also have to keep track of more substeps now.
                num_substeps *= 2;
                // We're changing a loop invariant, so we have to change the loop 
                // variable too. Gross.
                substep *= 2;
              }

              // Remove the entity from this bin.
              int_array_swap(entities, i, entities->size-1);
              int_array_resize(entities, entities->size-1);

              // Add it to the bin with the smaller time step.
              int smaller_bin = bin + 1;
              int_array_append(entities_for_bin[smaller_bin], e);
            }

            // Otherwise, it's possible that the new time step is bigger, and 
            // that the entity could be moved to another bin because of it. 
            // However, this is only possible if we're on a substep that aligns
            // this bin with the next biggest one.
            else if (new_dt > dt_for_bin)
            {
              if ((bin > 0) && (substep % 2 == 0))
              {
                // Remove the entity from this bin.
                int_array_swap(entities, i, entities->size-1);
                int_array_resize(entities, entities->size-1);

                // Add it to the bin with the bigger time step.
                int bigger_bin = bin - 1;
                int_array_append(entities_for_bin[bigger_bin], e);
              }
            }
          }
        }
      }

      // Update the current time.
      t_cur += substep_dt;
    }

    // Clean up.
    for (int bin = 0; bin < num_bins; ++bin)
      int_array_free(entities_for_bin[bin]);
    polymec_free(entities_for_bin);

    // Report!
    log_debug("async_solver_step: Took %d steps (max dt %g, min dt %g).", 
              actual_substeps, max_dt, min_dt);
  }

  // Call our post-step thingy if we have one.
  if (solver->vtable.pre_step != NULL)
    solver->vtable.post_step(solver->context, t + max_dt);

  // Clean up.
  polymec_free(dt_bins);

  STOP_FUNCTION_TIMER();
  return max_dt;
}

