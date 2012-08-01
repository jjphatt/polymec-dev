// poisson.c - The Poisson engine for Arbi.

#include <stdlib.h>
#include "poisson/poisson.h"
//#include "core/lin_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
} poisson_t;

// Benchmarks
static void run_paraboloid(poisson_t* p, options_t* opts)
{
  printf("rat = %s\n", options_value(opts, "rat"));
  printf("poo\n");
}

// Vtable stuff
static void* poisson_ctor(options_t* opts)
{
  poisson_t* p = malloc(sizeof(poisson_t));
  return p;
}

static void poisson_run_benchmark(void* p, const char* benchmark, options_t* opts)
{
  if (!strcmp(benchmark, "paraboloid"))
  {
    run_paraboloid(p, opts);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "poisson: unknown benchmark: '%s'", benchmark);
    arbi_error(err);
  }
}

static void poisson_init(void* p, double t)
{
}

static void poisson_advance(void* p, double t, double dt)
{
}

static void poisson_plot(void* p, plot_interface_t* plot, double t, int step)
{
}

static void poisson_dtor(void* p)
{
  free(p);
}

static model_vtable poisson_vtable = 
{
  .ctor          = &poisson_ctor,
  .run_benchmark = &poisson_run_benchmark,
  .init          = &poisson_init,
  .advance       = &poisson_advance,
  .plot          = &poisson_plot,
  .dtor          = &poisson_dtor
};

// FIXME
static const char* poisson_usage = "poisson model\n\n";

void register_poisson()
{
  register_model("poisson", poisson_usage, poisson_vtable);
}

#ifdef __cplusplus
}
#endif

