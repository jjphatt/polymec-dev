// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"
#include "core/timer.h"
#include "geometry/polymesh_field.h"

// Constructs a new polymesh field with the given number of components
// on the given mesh.
polymesh_field_t* polymesh_field_new(polymesh_t* mesh,
                                     polymesh_centering_t centering,
                                     size_t num_components)
{
  ASSERT(num_components > 0);
  polymesh_field_t* field = polymec_malloc(sizeof(polymesh_field_t));
  field->mesh = mesh;
  field->centering = centering;
  field->num_components = num_components;
  switch (centering)
  {
    case POLYMESH_CELL: 
      field->num_local_values = mesh->num_cells; 
      field->ex = polymesh_cell_exchanger(mesh);
      break;
    case POLYMESH_FACE: 
      field->num_local_values = mesh->num_faces; 
      field->ex = polymesh_1v_face_exchanger_new(mesh);
      break;
    case POLYMESH_EDGE: 
      field->num_local_values = mesh->num_edges; 
      field->ex = NULL;
      break;
    case POLYMESH_NODE: 
      field->num_local_values = mesh->num_nodes;
      field->ex = polymesh_1v_node_exchanger_new(mesh);
  }
  if (field->ex != NULL)
    polymec_retain(field->ex);
  field->ex_token = -1;
  field->num_ghost_values = (centering == POLYMESH_CELL) ? mesh->num_ghost_cells : 0;
  field->capacity = field->num_local_values + field->num_ghost_values;
  field->data = polymec_calloc(sizeof(real_t) * num_components * field->capacity);

  return field;
}

void polymesh_field_free(polymesh_field_t* field)
{
  if (field->ex != NULL)
    polymec_release(field->ex);
  polymec_free(field->data);
  polymec_free(field);
}

void polymesh_field_exchange(polymesh_field_t* field)
{
  polymesh_field_start_exchange(field);
  polymesh_field_finish_exchange(field);
}

void polymesh_field_start_exchange(polymesh_field_t* field)
{
  ASSERT(field->centering != POLYMESH_EDGE);
  ASSERT(!polymesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();

  // Start the xy exchange.
  int stride = field->num_components;
  field->ex_token = exchanger_start_exchange(field->ex, field->data, stride, 0, MPI_REAL_T);
  STOP_FUNCTION_TIMER();
}

void polymesh_field_finish_exchange(polymesh_field_t* field)
{
  ASSERT(polymesh_field_is_exchanging(field));
  START_FUNCTION_TIMER();
  if (field->ex_token != -1)
    exchanger_finish_exchange(field->ex, field->ex_token);
  field->ex_token = -1;
  STOP_FUNCTION_TIMER();
}

bool polymesh_field_is_exchanging(polymesh_field_t* field)
{
  return (field->ex_token != -1);
}

real_enumerable_generator_t* polymesh_field_enumerate(polymesh_field_t* field)
{
  size_t num_values = field->num_components * field->capacity;
  return real_enumerable_generator_from_array(field->data, num_values, NULL);
}

