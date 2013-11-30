// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string.h>
#include <strings.h>
#include "core/unordered_map.h"
#include "core/constant_st_func.h"
#include "core/create_uniform_mesh.h"
#include "core/least_squares.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"
#include "poisson/interpreter_register_poisson_functions.h"
#include "poisson/register_poisson_benchmarks.h"

typedef enum
{
  FV, FVPM
} poisson_sim_type;

static void run_analytic_problem_on_mesh(mesh_t* mesh, 
                                         st_func_t* lambda, 
                                         st_func_t* rhs, 
                                         string_ptr_unordered_map_t* bcs, 
                                         options_t* options, 
                                         double t1, 
                                         double t2, 
                                         st_func_t* solution, 
                                         double* lp_norms)
{
  // Create the model.
  model_t* model = create_fv_poisson(mesh, lambda, rhs, bcs, solution, options);

  // Run the thing.
  model_run(model, t1, t2, INT_MAX);

  // Calculate the Lp norms of the error and write it to Lp_norms.
  model_compute_error_norms(model, solution, lp_norms);

  // Clean up.
  model_free(model);
}

static void run_analytic_problem_on_points(point_cloud_t* point_cloud, 
                                           st_func_t* lambda, 
                                           st_func_t* rhs, 
                                           string_ptr_unordered_map_t* bcs, 
                                           options_t* options, 
                                           double t1, 
                                           double t2, 
                                           st_func_t* solution, 
                                           double* lp_norms)
{
  // Create the model.
  model_t* model = create_fvpm_poisson(point_cloud, lambda, rhs, bcs, solution, options);

  // Run the thing.
  model_run(model, t1, t2, INT_MAX);

  // Calculate the Lp norms of the error and write it to Lp_norms.
  model_compute_error_norms(model, solution, lp_norms);

  // Clean up.
  model_free(model);
}

static void laplace_1d_solution(void* context, point_t* x, double t, double* phi)
{
  phi[0] = 1.0 + 2.0*x->x;
}

static void laplace_1d_solution_grad(void* context, point_t* x, double t, double* grad_phi)
{
  grad_phi[0] = 2.0;
  grad_phi[1] = 0.0;
  grad_phi[2] = 0.0;
}

static void poisson_run_laplace_1d(options_t* options, 
                                   poisson_sim_type sim_type, 
                                   int dim)
{
  // Parse any benchmark-specific options.
  bool all_dirichlet = false;
  bool reversed_bcs = false;
  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
    else if (!strcmp(bcs_opt, "reversed"))
      reversed_bcs = true;
  }

  // Lambda is one for Laplace's equation.
  double o = 1.0;
  st_func_t* one = constant_st_func_new(1, &o);

  // RHS function is zero for Laplace's equation.
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);

  // Analytic solution and gradient.
  st_func_t* sol = st_func_from_func("laplace_1d_sol", laplace_1d_solution,
                                     ST_INHOMOGENEOUS, ST_CONSTANT, 1);
  st_func_t* sol_grad = st_func_from_func("laplace_1d_sol_grad", laplace_1d_solution_grad,
                                          ST_INHOMOGENEOUS, ST_CONSTANT, 3);

  // Boundary conditions: Dirichlet on -x/+x (unless they've been reversed).
  string_ptr_unordered_map_t* bcs = string_ptr_unordered_map_new();
  if (!reversed_bcs)
  {
    string_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    string_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(0.0, 1.0, sol_grad));
    string_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(0.0, 1.0, sol_grad));
  }

  // Transverse faces.
  if (all_dirichlet || reversed_bcs)
  {
    // Dirichlet BCs.
    string_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    // Homogeneous Neumann BCs.
    string_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Run time.
  double t = 0.0;

  // Base resolution, number of runs.
  int N0 = 1;
  int num_runs = 1;
  switch(dim)
  {
    case 1: 
      N0 = 32;
      num_runs = 4;
      break;
    case 2:
      N0 = 16;
      num_runs = 3;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }

  // Do a convergence study.
  double Lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    double dx = 1.0/Nx;
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    if (dim == 1)
      bbox.y2 = bbox.z2 = dx;
    if (dim == 2)
    {
      Ny = Nx;
      bbox.z2 = dx;
    }
    if (dim == 3)
      Nz = Nx;

    if (sim_type == FV)
    {
      mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);
      tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
      string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
      run_analytic_problem_on_mesh(mesh, one, zero, bcs_copy, options, t, t, sol, Lp_norms[iter]);
    }
    else
    {
      point_t* points = malloc(sizeof(point_t) * Nx*Ny*Nz);
//      tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
      point_cloud_t* cloud = point_cloud_new(MPI_COMM_WORLD, points, Nx*Ny*Nz);
      free(points);
      string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
      run_analytic_problem_on_points(cloud, one, zero, bcs_copy, options, t, t, sol, Lp_norms[iter]);
    }

    // If we run in 1D or 2D, we need to adjust the norms.
    if (dim == 1)
    {
      Lp_norms[iter][1] *= Nx*Nx;
      Lp_norms[iter][2] *= Nx*Nx;
    }
    else if (dim == 2)
    {
      Lp_norms[iter][1] *= Nx;
      Lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, Lp_norms[iter][1], Lp_norms[iter][2], Lp_norms[iter][0]);
  }

  if ((num_runs > 2) && (Lp_norms[0][2] > 0.0))
  {
    // Fit the log of the L2 norms to a line.
    double log_N_ratios[num_runs-1], log_L2_ratios[num_runs-1];
    for (int i = 0; i < num_runs-1; ++i)
    {
      log_N_ratios[i] = log(pow(2.0, i));
      log_L2_ratios[i] = log(Lp_norms[i+1][2] / Lp_norms[0][2]);
    }
    double A, B, sigma;
    linear_regression(log_N_ratios, log_L2_ratios, num_runs-1, &A, &B, &sigma);
    double rate = -A;
    model_report_conv_rate(options, rate, sigma);
  }

  // Clean up.
  string_ptr_unordered_map_free(bcs);
}

static void run_fv_laplace_1d(options_t* options)
{
  poisson_run_laplace_1d(options, FV, 1);
}

static void run_fv_laplace_1d_2(options_t* options)
{
  poisson_run_laplace_1d(options, FV, 2);
}

static void run_fv_laplace_1d_3(options_t* options)
{
  poisson_run_laplace_1d(options, FV, 3);
}

static void run_fvpm_laplace_1d(options_t* options)
{
  poisson_run_laplace_1d(options, FVPM, 1);
}

static void run_fvpm_laplace_1d_2(options_t* options)
{
  poisson_run_laplace_1d(options, FVPM, 2);
}

static void run_fvpm_laplace_1d_3(options_t* options)
{
  poisson_run_laplace_1d(options, FVPM, 3);
}

static void paraboloid_solution(void* context, point_t* x, double t, double* phi)
{
  double r2 = x->x*x->x + x->y*x->y; // Distance from center axis.
  phi[0] = 1.0 + r2;
//  printf("phi(%g, %g) = %g\n", x->x, x->y, phi[0]);
}

static void poisson_run_paraboloid(options_t* options, poisson_sim_type sim_type, int dim)
{
  ASSERT((dim == 2) || (dim == 3));

  // Parse model-specific options.
  bool offcenter = false;
  bool all_dirichlet = false;

  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
  }
  char *geom = options_value(options, "geometry");
  if (geom != NULL)
  {
    // The mesh can be generated off-center so that the origin is
    // at the lower left.
    if (!strcasecmp(geom, "offcenter"))
      offcenter = true;
  }

  // Lambda is one for Laplace's equation.
  double o = 1.0;
  st_func_t* one = constant_st_func_new(1, &o);

  // RHS function.
  double four = 4.0;
  st_func_t* rhs = constant_st_func_new(1, &four);

  // Analytic solution.
  st_func_t* sol = st_func_from_func("paraboloid", paraboloid_solution,
                                     ST_INHOMOGENEOUS, ST_CONSTANT, 1);

  // Set up a Dirichlet boundary condition along each of the outside faces.
  string_ptr_unordered_map_t* bcs = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(1.0, 0.0, sol));
  string_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(1.0, 0.0, sol));
  
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  if (all_dirichlet)
  {
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(1.0, 0.0, sol));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    // Set up a homogeneous Neumann boundary condition on +/- z.
    string_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
    string_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Start/end time.
  double t = 0.0;

  // Base resolution and number of runs.
  int N0 = 1;
  int num_runs = 1;
  switch (dim)
  {
    case 2:
      N0 = 16;
      num_runs = 4;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }
 
  // Do a convergence study.
  double Lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    if (dim > 1) Ny = Nx;
    if (dim > 2) Nz = Nx;
    bbox_t bbox;
    if (offcenter)
    {
      bbox.x1 = 0.0, bbox.x2 = 0.5, 
      bbox.y1 = 0.0, bbox.y2 = 0.5;
    }
    else
    {
      bbox.x1 = -0.5, bbox.x2 = 0.5, 
      bbox.y1 = -0.5, bbox.y2 = 0.5; 
    }
    bbox.z1 = 0.0, bbox.z2 = 1.0;
    if (dim == 2)
      bbox.z2 = 1.0/Nx;
    if (sim_type == FV)
    {
      mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);
      tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
      string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
      run_analytic_problem_on_mesh(mesh, one, rhs, bcs_copy, options, t, t, sol, Lp_norms[iter]);
    }
    else
    {
      point_t* points = malloc(sizeof(point_t) * Nx*Ny*Nz);
      point_cloud_t* cloud = point_cloud_new(MPI_COMM_WORLD, points, Nx*Ny*Nz);
      free(points);
//      tag_rectilinear_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
      string_ptr_unordered_map_t* bcs_copy = string_ptr_unordered_map_copy(bcs);
      run_analytic_problem_on_points(cloud, one, rhs, bcs_copy, options, t, t, sol, Lp_norms[iter]);
    }

    // If we run in 2D, we need to adjust the norms.
    if (dim == 2)
    {
      Lp_norms[iter][1] *= Nx;
      Lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = Ny = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, Lp_norms[iter][1], Lp_norms[iter][2], Lp_norms[iter][0]);
  }

  if ((num_runs > 2) && (Lp_norms[0][2] > 0.0))
  {
    // Fit the log of the L2 norms to a line.
    double log_N_ratios[num_runs-1], log_L2_ratios[num_runs-1];
    for (int i = 0; i < num_runs-1; ++i)
    {
      log_N_ratios[i] = log(pow(2.0, i));
      log_L2_ratios[i] = log(Lp_norms[i+1][2] / Lp_norms[0][2]);
    }
    double A, B, sigma;
    linear_regression(log_N_ratios, log_L2_ratios, num_runs-1, &A, &B, &sigma);
    double rate = -A;
    model_report_conv_rate(options, rate, sigma);
  }

  // Clean up.
  string_ptr_unordered_map_free(bcs);
}

static void run_fv_paraboloid(options_t* options)
{
  poisson_run_paraboloid(options, FV, 2);
}

static void run_fv_paraboloid_3(options_t* options)
{
  poisson_run_paraboloid(options, FV, 3);
}

static void run_fvpm_paraboloid(options_t* options)
{
  poisson_run_paraboloid(options, FVPM, 2);
}

static void run_fvpm_paraboloid_3(options_t* options)
{
  poisson_run_paraboloid(options, FVPM, 3);
}

void register_poisson_benchmarks(model_t* model)
{
  model_register_benchmark(model, "fv_laplace_1d", run_fv_laplace_1d, "Laplace's equation in 1D Cartesian coordinates (FV).");
  model_register_benchmark(model, "fv_laplace_1d_2", run_fv_laplace_1d_2, "Laplace's equation in 1D Cartesian coordinates (FV, run in 2D).");
  model_register_benchmark(model, "fv_laplace_1d_3", run_fv_laplace_1d_3, "Laplace's equation in 1D Cartesian coordinates (FV, run in 3D).");
  model_register_benchmark(model, "fvpm_laplace_1d", run_fvpm_laplace_1d, "Laplace's equation in 1D Cartesian coordinates (FVPM).");
  model_register_benchmark(model, "fvpm_laplace_1d_2", run_fvpm_laplace_1d_2, "Laplace's equation in 1D Cartesian coordinates (FVPM, run in 2D).");
  model_register_benchmark(model, "fvpm_laplace_1d_3", run_fvpm_laplace_1d_3, "Laplace's equation in 1D Cartesian coordinates (FVPM, run in 3D).");
  model_register_benchmark(model, "fv_paraboloid", run_fv_paraboloid, "A paraboloid solution to Poisson's equation (FV, 2D).");
  model_register_benchmark(model, "fv_paraboloid_3", run_fv_paraboloid_3, "A paraboloid solution to Poisson's equation (FV, 3D).");
  model_register_benchmark(model, "fvpm_paraboloid", run_fvpm_paraboloid, "A paraboloid solution to Poisson's equation (FVPM, 2D).");
  model_register_benchmark(model, "fvpm_paraboloid_3", run_fvpm_paraboloid_3, "A paraboloid solution to Poisson's equation (FVPM, 3D).");
}

