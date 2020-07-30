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
  struct infindividual* iis;
  struct infindividual* ii;
  int nlayers;
};

int var_sim_init(struct sim_pars* sim, const struct sim_pars sim_in);

#define sim_init(handle, ...) var_sim_init(handle, (const struct sim_pars){__VA_ARGS__});

int simulate(struct sim_pars const* sim, const gsl_rng* r);

static inline void gen_comm_period(struct sim_vars* sv)
{
  double m;

  sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars->kappa*sv->pars->tbar, 1./sv->pars->kappa);

  if(gsl_rng_uniform(sv->r) >= sv->pars->q) {

    if(isinf(sv->pars->kappaq)) m=sv->pars->mbar;
    else m=gsl_ran_gamma(sv->r, sv->pars->kappaq*sv->pars->mbar, 1./sv->pars->kappaq);

    if(m<sv->ii->comm_period) sv->ii->comm_period=m;
  }
  printf("Comm period is %f\n",sv->ii->comm_period);
}

static inline void gen_comm_period_isolation(struct sim_vars* sv)
{
  sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars->kappa*sv->pars->tbar, 1./sv->pars->kappa);
  printf("Comm period is %f\n",sv->ii->comm_period);
}

#endif
