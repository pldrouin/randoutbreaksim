#include "simulation.h"

int var_sim_init(struct simulation* sim, const struct simulation sim_in)
{
  sim->tbar=(sim_in.tbar?sim_in.tbar:10);
  sim->p=(sim_in.p?sim_in.p:0.5);
  sim->lambda=(sim_in.lambda?sim_in.lambda:0.11);
  sim->kappa=(sim_in.kappa?sim_in.kappa:1);
  sim->q=(sim_in.q?sim_in.q:0.6);
  sim->mbar=(sim_in.mbar?sim_in.mbar:5);
  sim->kappaq=(sim_in.kappaq?sim_in.kappaq:3);
  sim->nsim=(sim_in.nsim?sim_in.nsim:10);
  sim->nstart=(sim_in.nstart?sim_in.nstart:1);
  sim->tmax=(sim_in.tmax?sim_in.tmax:100);
  return 0;
}

int simulate(struct simulation* sim, const gsl_rng* r)
{
  struct iillist list;
  struct infindividual* ii;
  int i;

  iillist_init(&list);

  for(i=sim->nstart-1; i>=0; --i) {
    ii=iillist_newelement(&list);
  }

  iillist_clear(&list);

  return 0;
}
