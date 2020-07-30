#ifndef _SIMULATION_
#define _SIMULATION_

#include <stdint.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "infindividual.h"

#define INIT_N_LAYERS (16)

#define DEBUG_PRINTF(fmt, ...)
//#define DEBUG_PRINTF(fmt, ...) printf(fmt,__VA_ARGS__)

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
  struct sim_pars pars;
  gsl_rng const* r;
  struct infindividual* iis;
  struct infindividual* ii;
  uint32_t nlayers;
  void (*gen_trunc_comm_period_func)(struct sim_vars*);
};

void var_sim_init(struct sim_pars* sim, const struct sim_pars sim_in);
void sim_vars_init(struct sim_vars* sv, const gsl_rng* r);

#define sim_init(handle, r, ...) {var_sim_init(&(handle)->pars, (const struct sim_pars){__VA_ARGS__}); sim_vars_init(handle,r);}
void sim_free(struct sim_vars* sim){free(sim->iis);}

int simulate(struct sim_vars* sv);

static inline void gen_trunc_comm_period(struct sim_vars* sv)
{
  sv->ii->trunc_comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa);
  double time_left=sv->pars.tmax-(sv->ii-1)->event_time;

  if(sv->ii->trunc_comm_period > time_left) {
    sv->ii->trunc_comm_period=time_left;
    sv->ii->infectious_at_tmax=true;

  } else sv->ii->infectious_at_tmax=false;
  DEBUG_PRINTF("Comm period is %f%s\n",sv->ii->trunc_comm_period,(sv->ii->infectious_at_tmax?" (reached end)":""));
}

static inline void gen_trunc_comm_period_isolation(struct sim_vars* sv)
{
  double m;

  sv->ii->trunc_comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa);

  if(gsl_rng_uniform(sv->r) >= sv->pars.q) {

    if(isinf(sv->pars.kappaq)) m=sv->pars.mbar;
    else m=gsl_ran_gamma(sv->r, sv->pars.kappaq*sv->pars.mbar, 1./sv->pars.kappaq);

    if(m<sv->ii->trunc_comm_period) sv->ii->trunc_comm_period=m;
  }
  double time_left=sv->pars.tmax-(sv->ii-1)->event_time;

  if(sv->ii->trunc_comm_period > time_left) {
    sv->ii->trunc_comm_period=time_left;
    sv->ii->infectious_at_tmax=true;

  } else sv->ii->infectious_at_tmax=false;
  DEBUG_PRINTF("Comm period is %f%s\n",sv->ii->trunc_comm_period,(sv->ii->infectious_at_tmax?" (reached end)":""));
}

#endif
