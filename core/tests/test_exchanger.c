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
#include "cmockery.h"
#include "core/exchanger.h"

void test_exchanger_new(void** state)
{
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  exchanger_free(exchanger);
}

void test_exchanger_construct(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
  }
  exchanger_free(exchanger);
}

void test_exchanger_construct_and_delete(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);

    exchanger_delete_send(exchanger, 1);
    exchanger_delete_receive(exchanger, 1);
  }
  exchanger_free(exchanger);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exchanger_new),
    unit_test(test_exchanger_construct),
    unit_test(test_exchanger_construct_and_delete)
  };
  return run_tests(tests);
}
