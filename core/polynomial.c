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

#include <gc/gc.h>
#include "core/polynomial.h"
#include "core/slist.h"

struct polynomial_t 
{
  int degree, num_terms;
  double* coeffs;
  int *x_pow, *y_pow, *z_pow;
  point_t x0;
};

static void polynomial_free(void* ctx, void* dummy)
{
  polynomial_t* p = ctx;
  free(p->coeffs);
  free(p->x_pow);
  free(p->y_pow);
  free(p->z_pow);
}

static const int N_coeffs[5] = {1, 4, 10, 20, 35};

int polynomial_basis_dim(int degree)
{
  return N_coeffs[degree];
}

// Powers of x, y, and z in the standard basis.
static int std_x_pow[35] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 
                            3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 
                            4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 
                            0, 0, 0, 0, 0};
static int std_y_pow[35] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        4, 3, 2, 1, 0};
static int std_z_pow[35] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 
                            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                            0, 1, 2, 3, 4};

polynomial_t* polynomial_new(int degree, double* coeffs, point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(degree <= 4);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = malloc(sizeof(double) * N_coeffs[degree]);
  memcpy(p->coeffs, coeffs, sizeof(double) * N_coeffs[degree]);
  p->x_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->x_pow, std_x_pow, sizeof(int) * N_coeffs[degree]);
  p->y_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->y_pow, std_y_pow, sizeof(int) * N_coeffs[degree]);
  p->z_pow = malloc(sizeof(int) * N_coeffs[degree]);
  memcpy(p->z_pow, std_z_pow, sizeof(int) * N_coeffs[degree]);
  p->num_terms = N_coeffs[degree];
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

polynomial_t* polynomial_from_monomials(int degree, int num_monomials, double* coeffs, 
                                        int* x_powers, int* y_powers, int* z_powers, 
                                        point_t* x0)
{
  ASSERT(degree >= 0);
  ASSERT(num_monomials > 0);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->degree = degree;
  p->coeffs = malloc(sizeof(double) * num_monomials);
  memcpy(p->coeffs, coeffs, sizeof(double) * num_monomials);
  p->x_pow = malloc(sizeof(int) * num_monomials);
  memcpy(p->x_pow, x_powers, sizeof(int) * num_monomials);
  p->y_pow = malloc(sizeof(int) * num_monomials);
  memcpy(p->y_pow, y_powers, sizeof(int) * num_monomials);
  p->z_pow = malloc(sizeof(int) * num_monomials);
  memcpy(p->z_pow, z_powers, sizeof(int) * num_monomials);
  if (x0 != NULL)
    p->x0 = *x0;
  else
  {
    p->x0.x = 0.0, p->x0.y = 0.0, p->x0.z = 0.0;
  }
  p->num_terms = num_monomials;
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

polynomial_t* polynomial_clone(polynomial_t* p)
{
  polynomial_t* q = GC_MALLOC(sizeof(polynomial_t));
  q->degree = p->degree;
  q->num_terms = p->num_terms;
  q->coeffs = malloc(sizeof(double) * q->num_terms);
  memcpy(q->coeffs, p->coeffs, sizeof(double) * q->num_terms);
  q->x_pow = malloc(sizeof(int) * q->num_terms);
  memcpy(q->x_pow, p->x_pow, sizeof(int) * q->num_terms);
  q->y_pow = malloc(sizeof(int) * q->num_terms);
  memcpy(q->y_pow, p->y_pow, sizeof(int) * q->num_terms);
  q->z_pow = malloc(sizeof(int) * q->num_terms);
  memcpy(q->z_pow, p->z_pow, sizeof(int) * q->num_terms);
  q->x0 = p->x0;
  GC_register_finalizer(q, polynomial_free, q, NULL, NULL);
  return q;
}

int polynomial_degree(polynomial_t* p)
{
  return p->degree;
}

int polynomial_num_terms(polynomial_t* p)
{
  return p->num_terms;
}

double* polynomial_coeffs(polynomial_t* p)
{
  return p->coeffs;
}

point_t* polynomial_x0(polynomial_t* p)
{
  return &p->x0;
}

double polynomial_value(polynomial_t* p, point_t* x)
{
  int pos = 0, x_pow, y_pow, z_pow;
  double coeff, val = 0.0;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    val += coeff * pow(x->x - p->x0.x, x_pow) * 
                   pow(x->y - p->x0.y, y_pow) * 
                   pow(x->z - p->x0.z, z_pow);
  }
  return val;
}

static int fact(int x)
{
  if ((x == 0) || (x == 1))
    return 1;
  else return fact(x-1);
}

double polynomial_deriv(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x)
{
  ASSERT(x_deriv >= 0);
  ASSERT(y_deriv >= 0);
  ASSERT(z_deriv >= 0);

  if (x_deriv + y_deriv + z_deriv > p->degree)
    return 0.0;

  int pos = 0, x_pow, y_pow, z_pow;
  double coeff, val = 0.0;
  while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    double x_term = pow(x->x - p->x0.x, x_pow - x_deriv) * fact(x_pow)/fact(x_deriv);
    double y_term = pow(x->y - p->x0.y, y_pow - y_deriv) * fact(y_pow)/fact(y_deriv);
    double z_term = pow(x->z - p->x0.z, z_pow - z_deriv) * fact(z_pow)/fact(z_deriv);
    val += coeff * x_term * y_term * z_term;
  }
  return val;
}

bool polynomial_next(polynomial_t* p, int* pos, double* coeff, int* x_power, int* y_power, int* z_power)
{
  if (*pos >= N_coeffs[p->degree])
    return false;
  *coeff = p->coeffs[*pos];
  *x_power = p->x_pow[*pos];
  *y_power = p->y_pow[*pos];
  *z_power = p->z_pow[*pos];
  ++(*pos);
  return true;
}

// Returns the index of the term within p for which the x, y, and z powers 
// match those given. If such a term is not found, returns -1. This function 
// uses a linear search.
static int matching_term_index(polynomial_t* p, int x_pow, int y_pow, int z_pow)
{
  for (int i = 0; i < p->num_terms; ++i)
  {
    if ((p->x_pow[i] == x_pow) && (p->y_pow[i] == y_pow) && (p->z_pow[i] == z_pow))
      return i;
  }
  return -1;
}

void polynomial_add(polynomial_t* p, double factor, polynomial_t* q)
{
  // There are two cases we need to consider for each term in q:
  // 1. We have a term matching this term in p, and the coefficient should 
  //    simply be added in.
  // 2. We have no such matching term, and a new term should be appended 
  //    to p.
  int_slist_t* terms_to_append = int_slist_new();
  for (int i = 0; i < q->num_terms; ++i)
  {
    int index = matching_term_index(p, q->x_pow[i], q->y_pow[i], q->z_pow[i]);
    if (index != -1)
      p->coeffs[index] += factor * q->coeffs[i];
    else
      int_slist_append(terms_to_append, i);
  }

  if (!int_slist_empty(terms_to_append))
  {
    int old_size = p->num_terms;
    p->num_terms = old_size + terms_to_append->size;
    p->coeffs = realloc(p->coeffs, sizeof(double) * p->num_terms);
    p->x_pow = realloc(p->x_pow, sizeof(int) * p->num_terms);
    p->y_pow = realloc(p->y_pow, sizeof(int) * p->num_terms);
    p->z_pow = realloc(p->z_pow, sizeof(int) * p->num_terms);

    for (int i = old_size; i < p->num_terms; ++i)
    {
      int j = int_slist_pop(terms_to_append, NULL);
      p->coeffs[i] = factor * q->coeffs[j];
      p->x_pow[i] = q->x_pow[j];
      p->y_pow[i] = q->y_pow[j];
      p->z_pow[i] = q->z_pow[j];
    }
  }

  int_slist_free(terms_to_append);
}

static void wrap_eval(void* context, point_t* x, double* result)
{
  polynomial_t* p = context;
  *result = polynomial_value(p, x);
}

static void wrap_eval_deriv(void* context, int deriv, point_t* x, double* result)
{
  polynomial_t* p = context;
  int result_size = pow(3, deriv);
  if (deriv > p->degree)
    memset(result, 0, sizeof(double) * result_size);
  else if (deriv == 0)
    *result = polynomial_value(p, x);
  else 
  {
    double coeff;
    int pos = 0, x_pow, y_pow, z_pow;
    while (polynomial_next(p, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    {
      // FIXME
      POLYMEC_NOT_IMPLEMENTED
    }
  }
}

static bool wrap_has_deriv(void* context, int deriv)
{
  // Polynomials are analytic.
  return true;
}

sp_func_t* polynomial_sp_func(polynomial_t* p)
{
  sp_vtable vtable = {.eval = wrap_eval,
                      .eval_deriv = wrap_eval_deriv,
                      .has_deriv = wrap_has_deriv,
                      .dtor = NULL};
  char name[128];
  snprintf(name, 128, "polynomial (p = %d)", p->degree);
  return sp_func_new(name, p, vtable, SP_INHOMOGENEOUS, 1);
}

