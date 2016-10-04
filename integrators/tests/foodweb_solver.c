// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/declare_nd_array.h"
#include "integrators/newton_solver.h"
#include "integrators/bj_newton_pc.h"

// We use this for some of the underlying data structures.
#include "sundials/sundials_direct.h"

//------------------------------------------------------------------------
//                          FOOD WEB PROBLEM
//------------------------------------------------------------------------
// This was adapted from KINSol's food web example problem, programmed by 
// Allan Taylor, Alan Hindmarsh, and Radu Serban @ LLNL.
// This test solves a nonlinear system that arises from a system
// of partial differential equations. The PDE system is a food web
// population model, with predator-prey interaction and diffusion
// on the unit square in two dimensions. The dependent variable
// vector is the following:
// 
//       1   2         ns
// c = (c , c ,  ..., c  )     (denoted by the variable cc)
// 
// and the PDE's are as follows:
//
//                    i       i
//         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
//                    xx      yy         i
//
//   where
//
//                   i             ns         j
//   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
//    i                           j=1
//
// The number of species is ns = 2 * np, with the first np being
// prey and the last np being predators. The number np is both the
// number of prey and predator species. The coefficients a(i,j),
// b(i), d(i) are:
//
//   a(i,i) = -AA   (all i)
//   a(i,j) = -GG   (i <= np , j >  np)
//   a(i,j) =  EE   (i >  np,  j <= np)
//   b(i) = BB * (1 + alpha * x * y)   (i <= np)
//   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
//   d(i) = DPREY   (i <= np)
//   d(i) = DPRED   ( i > np)
//
// The various scalar parameters are set using define's or in
// routine InitUserData.
//
// The boundary conditions are: normal derivative = 0, and the
// initial guess is constant in x and y, but the final solution
// is not.
//
// The PDEs are discretized by central differencing on an MX by
// MY mesh.
// 
// The nonlinear system is solved by KINSOL using the method
// specified in local variable globalstrat.
//
// Constraints are imposed to make all components of the solution
// positive.
// -----------------------------------------------------------------
// References:
//
// 1. Peter N. Brown and Youcef Saad,
//    Hybrid Krylov Methods for Nonlinear Systems of Equations
//    LLNL report UCRL-97645, November 1987.
//
// 2. Peter N. Brown and Alan C. Hindmarsh,
//    Reduced Storage Matrix Methods in Stiff ODE systems,
//    Lawrence Livermore National Laboratory Report  UCRL-95088,
//    Rev. 1, June 1987, and  Journal of Applied Mathematics and
//    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
//    description of the time-dependent version of this test
//    problem.)
// -----------------------------------------------------------------

/* Problem Constants */

#define NUM_SPECIES     6  /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */ 

#define PI       RCONST(3.1415926535898)   /* pi */ 

#define MX          8              // MX = number of x mesh points */
#define MY          8              // MY = number of y mesh points */
#define NSMX        (NUM_SPECIES * MX)
#define NEQ         (NSMX * MY)    // number of equations in the system 
#define AA          RCONST(1.0)    // value of coefficient AA in above eqns 
#define EE          RCONST(10000.) // value of coefficient EE in above eqns 
#define GG          RCONST(0.5e-6) // value of coefficient GG in above eqns 
#define BB          RCONST(1.0)    // value of coefficient BB in above eqns 
#define DPREY       RCONST(1.0)    // value of coefficient dprey above 
#define DPRED       RCONST(0.5)    // value of coefficient dpred above 
#define ALPHA       RCONST(1.0)    // value of coefficient alpha above 
#define AX          RCONST(1.0)    // total range of x variable 
#define AY          RCONST(1.0)    // total range of y variable 
#define FTOL        RCONST(1.e-7)  // ftol tolerance 
#define STOL        RCONST(1.e-13) // stol tolerance 
#define THOUSAND    RCONST(1000.0) // one thousand 
#define ZERO        RCONST(0.)     // 0. 
#define ONE         RCONST(1.0)    // 1. 
#define TWO         RCONST(2.0)    // 2. 
#define PREYIN      RCONST(1.0)    // initial guess for prey concentrations. 
#define PREDIN      RCONST(30000.0)// initial guess for predator concs.  

// User-defined vector access macro: IJ_Vptr 

// IJ_Vptr is defined in order to translate from the underlying 3D structure
// of the dependent variable vector to the 1D storage scheme for an N-vector.
// IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
// indices is = 0, jx = i, jy = j.    

#define IJ_Vptr(vv,i,j)   (&vv[i*NUM_SPECIES + j*NSMX])

typedef struct 
{
  real_t acoef[NUM_SPECIES][NUM_SPECIES], bcoef[NUM_SPECIES];
  real_t *rates;
  real_t cox[NUM_SPECIES], coy[NUM_SPECIES];
  real_t ax, ay, dx, dy;
  real_t uround, sqruround;
  long int mx, my, ns, np;

  // Sparsity graph for preconditioner matrix.
  adj_graph_t* sparsity;

} foodweb_t;

// Readability definitions used in other routines below.
// Note that these may cause warnings.
#pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
#define acoef  (data->acoef)
#define bcoef  (data->bcoef)
#define cox    (data->cox)
#define coy    (data->coy)

// Newly initialized food web data context.
static foodweb_t* foodweb_new()
{
  real_t *a1,*a2, *a3, *a4;

  foodweb_t* data = polymec_malloc(sizeof(foodweb_t));
  
  data->mx = MX;
  data->my = MY;
  data->ns = NUM_SPECIES;
  data->np = NUM_SPECIES/2;
  data->ax = AX;
  data->ay = AY;
  data->dx = (data->ax)/(MX-1);
  data->dy = (data->ay)/(MY-1);
  data->uround = UNIT_ROUNDOFF;
  data->sqruround = sqrt(data->uround);
  data->rates = polymec_malloc(sizeof(real_t) * NEQ);

  // Set up the coefficients a and b plus others found in the equations.
  long np = data->np;

  real_t dx2=(data->dx)*(data->dx), dy2=(data->dy)*(data->dy);

  for (long i = 0; i < np; i++) 
  {
    a1= &(acoef[i][np]);
    a2= &(acoef[i+np][0]);
    a3= &(acoef[i][0]);
    a4= &(acoef[i+np][np]);

    // Fill in the portion of acoef in the four quadrants, row by row...
    for (long j = 0; j < np; j++) 
    {
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    // ...and then change the diagonal elements of acoef to -AA.
    acoef[i][i]=-AA;
    acoef[i+np][i+np] = -AA;

    bcoef[i] = BB;
    bcoef[i+np] = -BB;

    cox[i]=DPREY/dx2;
    cox[i+np]=DPRED/dx2;

    coy[i]=DPREY/dy2;
    coy[i+np]=DPRED/dy2;
  }  

  // Construct a sparsity graph.
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_SELF, MX*MY);
  for (int i = 0; i < MX; ++i) 
  {
    // Set left/right index shifts, special at boundaries. 
    int i_left  = (i !=  0  ) ? -1 :  1;
    int i_right = (i != MX-1) ?  1 : -1;

    for (int j = 0; j < MY; ++j) 
    {
      // Set lower/upper index shifts, special at boundaries. 
      int j_down = (j != 0   ) ? -1 :  1;
      int j_up   = (j != MY-1) ?  1 : -1;
    
      // Find the edges for the vertex corresponding to (jx, jy).
      int num_edges = 0;
      int edges[4];
      if (j > 0)
        edges[num_edges++] = ARRAY_INDEX_2D(MX, MY, i, j+j_down); // lower
      if (j < (MY-1))
        edges[num_edges++] = ARRAY_INDEX_2D(MX, MY, i, j+j_up); // upper
      if (i > 0)
        edges[num_edges++] = ARRAY_INDEX_2D(MX, MY, i+i_left, j); // left
      if (i < (MX-1))
        edges[num_edges++] = ARRAY_INDEX_2D(MX, MY, i+i_right, j); // right

      // Set the edges within the sparsity graph.
      int idx = ARRAY_INDEX_2D(MX, MY, i, j);
      adj_graph_set_num_edges(sparsity, idx, num_edges);
      memcpy(adj_graph_edges(sparsity, idx), edges, sizeof(int) * num_edges);
    }
  }

  data->sparsity = adj_graph_new_with_block_size(sparsity, NUM_SPECIES);
  adj_graph_free(sparsity);

  return data;
}

// Food web data context destructor.
static void foodweb_dtor(void* context)
{
  foodweb_t* data = context;
  
  polymec_free(data->rates);
  adj_graph_free(data->sparsity);
  polymec_free(data);
}

// Dot product routine for real_t arrays 
static real_t dot_prod(long int size, real_t* x1, real_t* x2)
{
  long int i;
  real_t *xx1, *xx2, temp = ZERO;
  
  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) 
    temp += (*xx1++) * (*xx2++);

  return temp;
}

// Interaction rate function routine 
static void web_rate(void* context, real_t xx, real_t yy, real_t *cxy, real_t *ratesxy)
{
  foodweb_t* data = context;
  
  for (int i = 0; i < NUM_SPECIES; ++i)
    ratesxy[i] = dot_prod(NUM_SPECIES, cxy, acoef[i]);
  
  real_t fac = ONE + ALPHA * xx * yy;
  
  for (int i = 0; i < NUM_SPECIES; i++)
    ratesxy[i] = cxy[i] * ( bcoef[i] * fac + ratesxy[i] );  
}

// System function
static int foodweb_func(void* context, real_t t, real_t* U, real_t* F)
{
  foodweb_t* data = context;
  
  DECLARE_3D_ARRAY(real_t, U_ijk, U, MX, MY, NUM_SPECIES);
  DECLARE_3D_ARRAY(real_t, F_ijk, F, MX, MY, NUM_SPECIES);
  DECLARE_3D_ARRAY(real_t, r_ijk, data->rates, MX, MY, NUM_SPECIES);

  real_t delx = data->dx;
  real_t dely = data->dy;
  
  // Loop over all mesh points, evaluating rate array at each point
  for (int i = 0; i < MX; ++i) 
  {
    real_t xi = delx*i;

    // Set left/right index shifts, special at boundaries. 
    int i_left  = (i !=  0  ) ? -1 :  1;
    int i_right = (i != MX-1) ?  1 : -1;

    for (int j = 0; j < MY; ++j) 
    {
      real_t yj = dely*j;

      // Set lower/upper index shifts, special at boundaries. 
      int j_down = (j != 0   ) ? -1 :  1;
      int j_up   = (j != MY-1) ?  1 : -1;

      // Get species interaction rate array at (xi,yj) 
      web_rate(data, xi, yj, U_ijk[i][j], r_ijk[i][j]);
      
      for(int s = 0; s < NUM_SPECIES; s++) 
      {
        // Differencing in x direction 
        real_t dcyli = U_ijk[i][j][s] - U_ijk[i][j+j_down][s];
        real_t dcyui = U_ijk[i][j+j_up][s] - U_ijk[i][j][s];
        
        // Differencing in y direction 
        real_t dcxli = U_ijk[i][j][s] - U_ijk[i+i_left][j][s];
        real_t dcxri = U_ijk[i+i_right][j][s] - U_ijk[i][j][s];

        // Compute F at (xi,yj) 
        F_ijk[i][j][s] = (coy)[s] * (dcyui - dcyli) +
                         (cox)[s] * (dcxri - dcxli) + 
                         r_ijk[i][j][s];
      }
    }
  }

  return 0;
}

// Function for accumulating a Jacobian column value into a column map.
static void accumulate_J_value(index_real_unordered_map_t* col_map, 
                               index_t column, 
                               real_t value)
{
  real_t* val_ptr = index_real_unordered_map_get(col_map, column);
  if (val_ptr == NULL) // New value is inserted.
    index_real_unordered_map_insert(col_map, column, value);
  else // accumulate the given value into the existing one.
    index_real_unordered_map_insert(col_map, column, value + *val_ptr);
}

// Function for inserting a row of values into the Jacobian matrix.
static void insert_J_values(index_t row, 
                            index_real_unordered_map_t* col_map, 
                            krylov_matrix_t* J)
{
  index_t num_cols = (index_t)col_map->size;
  index_t indices[num_cols];
  real_t values[num_cols];
  int pos = 0, k = 0;
  while (index_real_unordered_map_next(col_map, &pos, &indices[k], &values[k])) ++k;
  krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
  index_real_unordered_map_clear(col_map);
}

// Jacobian function
static int foodweb_J(void* context, 
                     real_t t, real_t* U, real_t* F, 
                     krylov_matrix_t* J)
{
  foodweb_t* data = context;
  
  DECLARE_3D_ARRAY(real_t, U_ijk, U, MX, MY, NUM_SPECIES);
  DECLARE_3D_ARRAY(real_t, r_ijk, data->rates, MX, MY, NUM_SPECIES);

  real_t delx = data->dx;
  real_t dely = data->dy;
  
  // Maps from indices to their values in the Jacobian.
  index_real_unordered_map_t* I_map = index_real_unordered_map_new();

  // Loop over all mesh points, evaluating rate array at each point
  for (int i = 0; i < MX; ++i) 
  {
    real_t xi = delx*i;

    // Set left/right index shifts, special at boundaries. 
    int i_left  = (i !=  0  ) ? -1 :  1;
    int i_right = (i != MX-1) ?  1 : -1;

    for (int j = 0; j < MY; ++j) 
    {
      real_t yj = dely*j;

      // Set lower/upper index shifts, special at boundaries. 
      int j_down = (j != 0   ) ? -1 :  1;
      int j_up   = (j != MY-1) ?  1 : -1;

      // Get species interaction rate array at (xi,yj) 
      web_rate(data, xi, yj, U_ijk[i][j], r_ijk[i][j]);
      for(int s = 0; s < NUM_SPECIES; s++) 
      {
        // There are up to 4 + NUM_SPECIES Jacobian contributions for each 
        // species: 5 stencil points, and NUM_SPECIES-1 reaction terms in 
        // the same location. We lump the diagonal in with the reactions 
        // for simplicity.
        real_t J_left = 0.0, J_right = 0.0, 
               J_up = 0.0, J_down = 0.0,
               J_rxn[NUM_SPECIES];
        index_t I_left  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_left, j, 0),
                I_right = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_right, j, 0),
                I_up    = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+j_up, 0),
                I_down  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+j_down, 0),
                I_rxn[NUM_SPECIES];

        // Compute reaction rate derivatives.
        for (int s1 = 0; s1 < NUM_SPECIES; ++s1)
        {
          I_rxn[s1] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j, s1);
          J_rxn[s1] = 2.0 * r_ijk[i][j][s]/U_ijk[i][j][s1];
        }

        // Differencing in x direction.
        J_rxn[s] += -2.0 * (cox)[s];
        J_left   = (cox)[s];
        J_right  = (cox)[s];

        // Differencing in y direction.
        J_rxn[s] += -2.0 * (coy)[s];
        J_up     = (coy)[s];
        J_down   = (coy)[s];

        // Accumulate the column values.
        accumulate_J_value(I_map, I_left, J_left);
        accumulate_J_value(I_map, I_right, J_right);
        accumulate_J_value(I_map, I_up, J_up);
        accumulate_J_value(I_map, I_down, J_down);
        for (int s1 = 0; s1 < NUM_SPECIES; ++s1)
          accumulate_J_value(I_map, I_rxn[s1], J_rxn[s1]);

        // Stick the data into the matrix.
        index_t row = I_rxn[s];
        insert_J_values(row, I_map, J);
      }
    }
  }

  // Make sure the matrix is assembled.
  krylov_matrix_assemble(J);
//static bool first = true;
//if (first)
//{
//FILE* f = fopen("J.txt", "w");
//krylov_matrix_fprintf(J, f);
//fclose(f);
//first = false;
//}

  // Clean up.
  index_real_unordered_map_free(I_map);

  return 0;
}

// Constructor for Jacobian-Free Newton-Krylov food web solver with 
// no preconditioner.
static newton_solver_t* jfnk_foodweb_solver_new(foodweb_t* data, newton_pc_t* precond)
{
  // Set up a nonlinear solver using GMRES with a full Newton step.
  newton_solver_t* solver = jfnk_newton_solver_new(MPI_COMM_SELF, NEQ, 0, data,
                                                   foodweb_func, NULL, foodweb_dtor, 
                                                   precond, NEWTON_GMRES, 15, 2);

  // Enforce positivity on all components.
  real_t constraints[NEQ];
  for (int i = 0; i < NEQ; ++i)
    constraints[i] = 2.0;
  newton_solver_set_constraints(solver, constraints);

  // Scale the U and F vectors.
  real_t species_scale[NUM_SPECIES];
  for (int s = 0; s < NUM_SPECIES/2; s++) 
    species_scale[s] = ONE;
  for (int s = NUM_SPECIES/2; s < NUM_SPECIES; s++)
    species_scale[s] = RCONST(0.00001);

  real_t scale[NEQ];
  for (int jy = 0; jy < MY; jy++) 
  {
    for (int jx = 0; jx < MX; jx++) 
    {
      real_t* sloc = IJ_Vptr(scale,jx,jy);
      for (int s = 0; s < NUM_SPECIES; s++) 
        sloc[s] = species_scale[s];
    }
  }
  newton_solver_set_U_scale(solver, scale);
  newton_solver_set_F_scale(solver, scale);

  return solver;
}

// Constructor for block-Jacobi-preconditioned food web solver.
newton_solver_t* block_jacobi_precond_foodweb_solver_new(void);
newton_solver_t* block_jacobi_precond_foodweb_solver_new()
{
  foodweb_t* data = foodweb_new();
  int block_size = NUM_SPECIES;
  newton_pc_t* precond = cpr_bj_newton_pc_new(MPI_COMM_WORLD, data, foodweb_func, NULL, NEWTON_PC_LEFT, data->sparsity, NEQ/block_size, 0, block_size);
  return jfnk_foodweb_solver_new(data, precond);
}

// Returns initial conditions.
real_t* foodweb_initial_conditions(void);
real_t* foodweb_initial_conditions()
{
  real_t* cc = polymec_malloc(sizeof(real_t) * NEQ);
  int i, jx, jy;
  real_t *cloc;
  real_t  ctemp[NUM_SPECIES];
  
  for (i = 0; i < NUM_SPECIES/2; i++)
    ctemp[i] = PREYIN;
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++) 
    ctemp[i] = PREDIN;

  for (jy = 0; jy < MY; jy++) 
  {
    for (jx = 0; jx < MX; jx++) 
    {
      cloc = IJ_Vptr(cc,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++) 
        cloc[i] = ctemp[i];
    }
  }
  return cc;
}

// Constructor for Inexact Newton-Krylov food web solver.
newton_solver_t* ink_foodweb_solver_new(krylov_factory_t* factory);
newton_solver_t* ink_foodweb_solver_new(krylov_factory_t* factory)
{
  foodweb_t* data = foodweb_new();
  matrix_sparsity_t* J_sparsity = matrix_sparsity_from_graph(data->sparsity, NULL);

  // Set up a nonlinear solver using GMRES with a full Newton step.
  newton_solver_t* solver = ink_newton_solver_new(MPI_COMM_SELF, factory, 
                                                  J_sparsity, data,
                                                  foodweb_func, foodweb_J, 
                                                  foodweb_dtor);

  // Enforce positivity on all components.
  real_t constraints[NEQ];
  for (int i = 0; i < NEQ; ++i)
    constraints[i] = 2.0;
  newton_solver_set_constraints(solver, constraints);

  // Scale the U and F vectors.
  real_t species_scale[NUM_SPECIES];
  for (int s = 0; s < NUM_SPECIES/2; s++) 
    species_scale[s] = ONE;
  for (int s = NUM_SPECIES/2; s < NUM_SPECIES; s++)
    species_scale[s] = RCONST(0.00001);

  real_t scale[NEQ];
  for (int jy = 0; jy < MY; jy++) 
  {
    for (int jx = 0; jx < MX; jx++) 
    {
      real_t* sloc = IJ_Vptr(scale,jx,jy);
      for (int s = 0; s < NUM_SPECIES; s++) 
        sloc[s] = species_scale[s];
    }
  }
  newton_solver_set_U_scale(solver, scale);
  newton_solver_set_F_scale(solver, scale);

  return solver;
}

