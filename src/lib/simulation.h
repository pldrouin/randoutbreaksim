#ifndef _SIMULATION_
#define _SIMULATION_

#include <stdint.h>

#include <gsl/gsl_rng.h>

#include "infindividual.h"

struct simulation
{
  double tbar;
  double p;
  double lambda;
  double kappa;
  double q;
  double mbar;
  double kappaq;
  uint32_t nsim;
  uint32_t nstart;
  uint32_t tmax;
};

int var_sim_init(struct simulation* sim, const struct simulation sim_in);

#define sim_init(handle, ...) var_sim_init(handle, (const struct simulation){__VA_ARGS__});

int simulate(struct simulation* sim, const gsl_rng* r);

#endif
