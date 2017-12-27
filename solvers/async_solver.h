// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ASYNC_SOLVER_H
#define POLYMEC_ASYNC_SOLVER_H

#include "core/polymec.h"

// The async_solver class provides an interface for solving a system of 
// differential equations asyncronously in time. A system consists of a set 
// of entities which can each be evolved using different time steps. An 
// entity can be a degree of freedom or a set of related degrees of freedom--
// the async_solver doesn't care, as long as you tell it how to integrate 
// each entity.
//
// To asynchronously integrate its entities, the solver determines valid 
// time step sizes for each entity, and then creates a nested hierarchy of 
// time steps in powers of two, starting with the largest step size. The 
// solver begins a step by integrating all entities using the smallest step, 
// and then proceeds with a subset of "active" entities that share that 
// smallest step. The other entities, which can be integrate with a larger step, 
// are "inactive," and are evolved according to some sort of "drift" approximation.
// 
// Unlike most other solvers, the async_solver doesn't impose a structure 
// or representation on the solution to the system of equations. All solution 
// information is encapsulated within the context manipulated by the solver.
// In this sense, the async_solver is really more of a "meta-solver" than 
// a solver.
typedef struct async_solver_t async_solver_t;

// This virtual table determines the implementation of the async_solver.
typedef struct
{
  // This function returns the number of entities in the system of differential
  // equations.
  int (*num_entities)(void* context);

  // This function initializes the entities in the system of differential 
  // equations at time t. Optional.
  void (*init)(void* context, real_t t);

  // This function returns the maximum allowed time step for the ith entity
  // at time t.
  real_t (*entity_dt)(void* context, int i, real_t t);

  // This function is called once at the beginning of an asynchronous step.
  // Optional.
  void (*pre_step)(void* context, real_t t);

  // This function integrates the given set of active entities over the 
  // interval [t, t + dt]. //active_dt is the timestep over which the 
//  // entities compute their time derivatives, and is greater than or 
//  // equal to dt.
  void (*active_step)(void* context, 
                      int* entities, size_t num_entities, 
                      real_t t, real_t dt);//, real_t active_dt);

  // This function integrates the given set of inactive entities over the 
  // interval [t, t + dt] using a "drift" approximation in which the 
  // derivatives remain the same over a larger step.
  void (*inactive_step)(void* context, 
                        int* entities, size_t num_entities, 
                        real_t t, real_t dt);

  // This function is called once at the end of an asynchronous step. Optional.
  void (*post_step)(void* context, real_t t);

  // This function destroys the state (context) when the solver 
  // is destroyed. Optional.
  void (*dtor)(void* context);

} async_solver_vtable;

// Creates an async_solver that uses the given context and virtual 
// table to integrate a system of differential equations. 
async_solver_t* async_solver_new(void* context,
                                 async_solver_vtable vtable);

// Frees the given async_solver.
void async_solver_free(async_solver_t* solver);

// Initializes the solution to the system of equations at time t.
void async_solver_init(async_solver_t* solver, real_t t);

// Returns the maximum time step size for all entities within the solver
// at time t.
real_t async_solver_max_dt(async_solver_t* solver, real_t t);

// Integrates the solution to the equations from time t to t + dt.
// Returns the size of the actual time step taken, which is the largest step
// for the entities in the system.
real_t async_solver_step(async_solver_t* solver, real_t t, real_t dt);

#endif

