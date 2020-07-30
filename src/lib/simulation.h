#ifndef _SIMULATION_
#define _SIMULATION_

#include <stdint.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "infindividual.h"

struct sim_pars
{
  double tbar;
  double p;
  double lambda;
  double kappa;
  double q;
  double mbar;
  double kappaq;
  double tmax;
  uint32_t nsim;
  uint32_t nstart;
};

struct sim_vars
{
  struct sim_pars const* pars;
  gsl_rng const* r;
  struct infindividual* ii;
  int layer;
};

int var_sim_init(struct sim_pars* sim, const struct sim_pars sim_in);

#define sim_init(handle, ...) var_sim_init(handle, (const struct sim_pars){__VA_ARGS__});

int simulate(struct sim_pars const* sim, const gsl_rng* r);

int infect_attendees(struct sim_vars const* sv, struct iillist* list);
int infect_attendees_isolation(struct sim_vars const* sv, struct iillist* list);

#endif
