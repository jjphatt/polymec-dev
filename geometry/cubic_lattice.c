// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include "core/polymec.h"
#include "core/mesh.h"
#include "geometry/cubic_lattice.h"
#include <gc/gc.h>

cubic_lattice_t* cubic_lattice_new(uint64_t nx, uint64_t ny, uint64_t nz)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  cubic_lattice_t* l = GC_MALLOC(sizeof(cubic_lattice_t));
  l->nx = nx;
  l->ny = ny;
  l->nz = nz;
  return l;
}

static int_int_unordered_map_t* cubic_lattice_generate_x_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the x-faces and map everything.
  for (int j = 0; j < lattice->ny; ++j)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_x_face(lattice, 0, j, k);
      int high = cubic_lattice_x_face(lattice, lattice->nx, j, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

static int_int_unordered_map_t* cubic_lattice_generate_y_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the y-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int k = 0; k < lattice->nz; ++k)
    {
      int low = cubic_lattice_y_face(lattice, i, 0, k);
      int high = cubic_lattice_y_face(lattice, i, lattice->ny, k);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

static int_int_unordered_map_t* cubic_lattice_generate_z_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int_int_unordered_map_t* map = int_int_unordered_map_new();
 
  // Go over the z-faces and map everything.
  for (int i = 0; i < lattice->nx; ++i)
  {
    for (int j = 0; j < lattice->ny; ++j)
    {
      int low = cubic_lattice_z_face(lattice, i, j, 0);
      int high = cubic_lattice_z_face(lattice, i, j, lattice->nz);
      int_int_unordered_map_insert(map, low, high);
      int_int_unordered_map_insert(map, high, low);
    }
  }

  return map;
}

static size_t cl_byte_size(void* obj)
{
  return 3 * sizeof(uint64_t);
}

static void* cl_byte_read(byte_array_t* bytes, size_t* offset)
{
  uint64_t data[3];
  byte_array_read_uint64s(bytes, 3, data, offset);
  return cubic_lattice_new(data[0], data[1], data[2]);
}

static void cl_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  uint64_t* l = obj;
  byte_array_write_uint64s(bytes, 3, l, offset);
}

serializer_t* cubic_lattice_serializer()
{
  return serializer_new("cubic_lattice", cl_byte_size,
                        cl_byte_read, cl_byte_write, NULL);
}

#if 0
periodic_bc_t* cubic_lattice_x_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_x_periodic_map, NULL);
}

periodic_bc_t* cubic_lattice_y_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_y_periodic_map, NULL);
}

periodic_bc_t* cubic_lattice_z_periodic_bc_new(const char* tag1, const char* tag2)
{
  return periodic_bc_new_with_map_func(tag1, tag2, cubic_lattice_generate_z_periodic_map, NULL);
}
#endif

