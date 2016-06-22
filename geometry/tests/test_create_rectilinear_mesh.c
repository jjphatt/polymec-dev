// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "core/silo_file.h"
#include "geometry/create_rectilinear_mesh.h"

void test_create_rectilinear_mesh(void** state)
{
  // Create a 10x10x10 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 11, ys, 11, zs, 11);
  assert_true(mesh_verify_topology(mesh, polymec_error));

  int nproc;
  MPI_Comm_size(mesh->comm, &nproc);
  int num_cells, num_faces, num_edges, num_nodes;
  MPI_Allreduce(&mesh->num_cells, &num_cells, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_faces, &num_faces, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_edges, &num_edges, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_nodes, &num_nodes, 1, MPI_INT, MPI_SUM, mesh->comm);
  assert_int_equal(10*10*10, num_cells);
  if (nproc > 1)
  {
    assert_true(mesh->num_ghost_cells > 0);
    assert_true(num_faces > 3*10*10*11);
    assert_true(num_edges > 3*10*11*11);
    assert_true(num_nodes > 11*11*11);
  }
  else
  {
    assert_int_equal(0, mesh->num_ghost_cells);
    assert_int_equal(3*10*10*11, num_faces);
    assert_int_equal(3*10*11*11, num_edges);
    assert_int_equal(11*11*11, num_nodes);
  }

  mesh_free(mesh);
}

void test_plot_rectilinear_mesh(void** state)
{
  // Create a 4x4x4 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 5, ys, 5, zs, 5);

  // Plot it.
  double ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "rectilinear_4x4x4", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "solution", "mesh", ones, NULL);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_4x4x1_rectilinear_mesh(void** state)
{
  // Create a 4x4x1 rectilinear mesh, which isn't a big deal in and of itself, 
  // but does seem to exhibit problems on certain numbers of processes.
  double xs[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double ys[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double zs[] = {0.0, 1.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 5, ys, 5, zs, 2);
  assert_true(mesh_verify_topology(mesh, polymec_error));

  int nproc, rank;
  MPI_Comm_size(mesh->comm, &nproc);
  MPI_Comm_rank(mesh->comm, &rank);

  // Check for correctness on the problematic 3-process case.
  if (nproc == 3)
  {
    // Create global cell indices.
    index_t vtx_dist[4] = {0, 5, 10, 16};
    index_t G[mesh->num_cells + mesh->num_ghost_cells];
    for (int i = 0; i < mesh->num_cells; ++i)
      G[i] = vtx_dist[rank] + i;
    exchanger_exchange(mesh_exchanger(mesh), G, 1, 0, MPI_INDEX_T);

    if (rank == 0)
    {
      assert_int_equal(5, mesh->num_cells);
      assert_int_equal(5, mesh->num_ghost_cells);

      // Check the layouts of local and ghost cells.
      // Note that ghost cells are created by walking through the 
      // local cells and trying to add ghosts in each of the -/+ x/y/z 
      // directions (in that order).
      index_t G0[] = {0, 1, 2, 3, 4, 5, 5, 6, 7, 8};
      for (int i = 0; i < 10; ++i)
        assert_int_equal(G0[i], G[i]);
    }
    else if (rank == 1)
    {
      assert_int_equal(5, mesh->num_cells);
      assert_int_equal(10, mesh->num_ghost_cells);

      // Check the layouts of local and ghost cells.
      // Note that ghost cells are created by walking through the 
      // local cells and trying to add ghosts in each of the -/+ x/y/z 
      // directions (in that order).
      index_t G1[] = {5, 6, 7, 8, 9, 1, 2, 3, 10, 4, 10, 4, 11, 12, 13};
      for (int i = 0; i < 15; ++i)
        assert_int_equal(G1[i], G[i]);
    }
    else if (rank == 2)
    {
      assert_int_equal(6, mesh->num_cells);
      assert_int_equal(5, mesh->num_ghost_cells);

      // Check the layouts of local and ghost cells.
      // Note that ghost cells are created by walking through the 
      // local cells and trying to add ghosts in each of the -/+ x/y/z 
      // directions (in that order).
      index_t G2[] = {10, 11, 12, 13, 14, 15, 6, 7, 8, 9, 9};
      for (int i = 0; i < 11; ++i)
        assert_int_equal(G2[i], G[i]);
    }
  }

  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_rectilinear_mesh),
    cmocka_unit_test(test_plot_rectilinear_mesh),
    cmocka_unit_test(test_create_4x4x1_rectilinear_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
