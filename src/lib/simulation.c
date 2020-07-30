#include "simulation.h"

#define INITIAL_LAYER_NUMBER (1)

int var_sim_init(struct sim_pars* sim, const struct sim_pars sim_in)
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

int simulate(struct sim_pars const* sim, const gsl_rng* r)
{
  struct iillist list;
  int i,j;
  double t, m;
  struct sim_vars sv={.pars=sim, .r=r, .layer=INITIAL_LAYER_NUMBER};
  double t_infect=0;

  iillist_init(&list);

  for(i=sim->nstart-1; i>=0; --i) {
    printf("initial individual %i\n",i);
    sv.ii=iillist_new_element(&list);
    printf("New individual %p\n",sv.ii);
    t=gsl_ran_gamma(r, sim->kappa*sim->tbar, 1./sim->kappa);

    if(sim->q && gsl_rng_uniform(r) >= sim->q) {

      if(isinf(sim->kappaq)) m=sim->mbar;
      else m=gsl_ran_gamma(r, sim->kappaq*sim->mbar, 1./sim->kappaq);

      if(m<t) t=m;
    }
    printf("Time is %f\n",t);
    sv.ii->nevents=gsl_ran_poisson(r, sim->lambda*t);
    printf("Nevents is %i\n",sv.ii->nevents);

    if(sv.ii->nevents) {
      sv.ii->etimes=(float*)realloc(sv.ii->etimes,sv.ii->nevents*sizeof(float));

      for(j=sv.ii->nevents-1; j>=0; --j) sv.ii->etimes[j]=t_infect+t*gsl_rng_uniform(r); //Note: unsorted
      sv.ii->etimesptr=sv.ii->etimes-INITIAL_LAYER_NUMBER-1;
      iillist_add_element(&list, sv.ii);

    } else iillist_recycle_detached_element(&list, sv.ii);
  }

  int (*infect_func)(struct sim_vars const*, struct iillist*)=(sim->q?infect_attendees_isolation:infect_attendees);
  struct infindividual* prev_ii=NULL;
  struct infindividual* next_ii=NULL;

  do {
    ++(sv.layer);
    printf("Layer is %i\n",sv.layer);
    sv.ii=list.head;
    printf("Head is %p\n",sv.ii);

    do {
      infect_func(&sv, &list);

      if(sv.ii->etimesptr+sv.layer-sv.ii->etimes==sv.ii->nevents-1) {
	next_ii=sv.ii->next;
	iillist_recycle_element(&list, prev_ii, sv.ii);
        prev_ii=sv.ii;
	sv.ii=next_ii;

      } else {
        prev_ii=sv.ii;
        sv.ii=sv.ii->next;
      }
      printf("Next individual is %p\n",sv.ii);

    } while(sv.ii);

  } while(list.head);

  iillist_clear(&list);

  return 0;
}

int infect_attendees(struct sim_vars const* sv, struct iillist* list)
{
  int ninf=gsl_ran_logarithmic(sv->r, sv->pars->p);
  int i,j;
  struct infindividual* ii;
  double t,m;

  for(i=ninf-1; i>=0; --i) {
    ii=iillist_new_element(list);
    printf("New individual %p\n",ii);
    t=gsl_ran_gamma(sv->r, sv->pars->kappa*sv->pars->tbar, 1./sv->pars->kappa);

    if(gsl_rng_uniform(sv->r) >= sv->pars->q) {

      if(isinf(sv->pars->kappaq)) m=sv->pars->mbar;
      else m=gsl_ran_gamma(sv->r, sv->pars->kappaq*sv->pars->mbar, 1./sv->pars->kappaq);

      if(m<t) t=m;
    }
    printf("Time is %f\n",t);
    ii->nevents=gsl_ran_poisson(sv->r, sv->pars->lambda*t);
    printf("Nevents is %i\n",ii->nevents);

    if(ii->nevents) {
      ii->etimes=(float*)realloc(ii->etimes,ii->nevents*sizeof(float));

      for(j=ii->nevents-1; j>=0; --j) ii->etimes[j]=sv->ii->etimesptr[sv->layer]+t*gsl_rng_uniform(sv->r); //Note: unsorted
      ii->etimesptr=ii->etimes-sv->layer-1;
      iillist_add_element(list, ii);

    } else iillist_recycle_detached_element(list, ii);
  }
  return 0;
}

int infect_attendees_isolation(struct sim_vars const* sv, struct iillist* list)
{
  int ninf=gsl_ran_logarithmic(sv->r, sv->pars->p);
  int i,j;
  struct infindividual* ii;
  double t;

  for(i=ninf-1; i>=0; --i) {
    ii=iillist_new_element(list);
    printf("New individual %p\n",ii);
    t=gsl_ran_gamma(sv->r, sv->pars->kappa*sv->pars->tbar, 1./sv->pars->kappa);
    printf("Time is %f\n",t);
    ii->nevents=gsl_ran_poisson(sv->r, sv->pars->lambda*t);
    printf("Nevents is %i\n",ii->nevents);

    if(ii->nevents) {
      ii->etimes=(float*)realloc(ii->etimes,ii->nevents*sizeof(float));

      for(j=ii->nevents-1; j>=0; --j) ii->etimes[j]=sv->ii->etimesptr[sv->layer]+t*gsl_rng_uniform(sv->r); //Note: unsorted
      ii->etimesptr=ii->etimes-sv->layer-1;
      iillist_add_element(list, ii);

    } else iillist_recycle_detached_element(list, ii);
  }
  return 0;
}
