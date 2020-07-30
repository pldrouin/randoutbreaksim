#include "simulation.h"

#define INIT_N_LAYERS (16)
#define II_ARRAY_GROW_FACT (1.5)

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
  int i;
  double m;
  struct sim_vars sv={.pars=sim, .r=r, .iis=(struct infindividual*)malloc(INIT_N_LAYERS*sizeof(struct infindividual)), .nlayers=INIT_N_LAYERS};
  void (*gen_comm_period_func)(struct sim_vars*)=(sim->q?gen_comm_period_isolation:gen_comm_period);
  sv.iis[0].event_time=0;

  for(i=sim->nstart-1; i>=0; --i) {
    printf("initial individual %i\n",i);
    sv.iis[1].comm_period=gsl_ran_gamma(r, sim->kappa*sim->tbar, 1./sim->kappa);

    if(sim->q && gsl_rng_uniform(r) >= sim->q) {

      if(isinf(sim->kappaq)) m=sim->mbar;
      else m=gsl_ran_gamma(r, sim->kappaq*sim->mbar, 1./sim->kappaq);

      if(m<sv.iis[1].comm_period) sv.iis[1].comm_period=m;
    }
    printf("Comm period is %f\n",sv.iis[1].comm_period);
    sv.ii=sv.iis+1;
    sv.ii->nevents=gsl_ran_poisson(r, sim->lambda*sv.ii->comm_period);
    printf("Nevents is %i\n",sv.ii->nevents);

    //If events for the current individual in the primary layer
    if(sv.ii->nevents) {
      sv.ii->curevent=0;
      sv.ii->event_time=sv.ii->comm_period*gsl_rng_uniform(r);
      printf("Event %i/%i at time %f\n",sv.ii->curevent,sv.ii->nevents,sv.ii->event_time);

      //Skip the events exceeding the limit
      while(sv.ii->event_time > sv.pars->tmax) {
	++(sv.ii->curevent);

	//If all events have been exhausted, continue with the next individual
	//in the primary layer
	if(sv.ii->curevent == sv.ii->nevents) goto next_pri_ii;
	sv.ii->event_time=sv.ii->comm_period*gsl_rng_uniform(r);
        printf("Event %i/%i at time %f\n",sv.ii->curevent,sv.ii->nevents,sv.ii->event_time);
      }
      sv.ii->ninfections=gsl_ran_logarithmic(r, sv.pars->p);
      sv.ii->curinfection=0;
      printf("Infection %i/%i\n",sv.ii->curinfection,sv.ii->ninfections);

      //Create a new infected individual
      for(;;) {
next_layer:
	printf("Move to next layer\n");
	++sv.ii;

	//If reaching the end of the allocated array, increase its size
	if(sv.ii==sv.iis+sv.nlayers) {
	  sv.nlayers*=II_ARRAY_GROW_FACT;
	  sv.iis=(struct infindividual*)realloc(sv.iis,sv.nlayers*sizeof(struct infindividual));
	}
	//Generate the communicable period appropriately
	gen_comm_period_func(&sv);

	//Generate the number of events
	sv.ii->nevents=gsl_ran_poisson(r, sim->lambda*sv.ii->comm_period);
	printf("Nevents is %i\n",sv.ii->nevents);

	//If the number of events is non-zero
	if(sv.ii->nevents) {
	  sv.ii->curevent=0;
	  //Generate the event time
gen_event:
	  sv.ii->event_time=(sv.ii-1)->event_time+sv.ii->comm_period*gsl_rng_uniform(r);
          printf("Event %i/%i at time %f\n",sv.ii->curevent,sv.ii->nevents,sv.ii->event_time);

	  //Loop over the events
	  for(;;) {

	    //If the event time does not exceed the limit
	    if(sv.ii->event_time <= sv.pars->tmax) {
	      //Generate the number of infections and the associated index for
	      //the current event
	      //Move to the next layer
	      sv.ii->ninfections=gsl_ran_logarithmic(r, sv.pars->p);
	      sv.ii->curinfection=0;
              printf("Infection %i/%i\n",sv.ii->curinfection,sv.ii->ninfections);
	      goto next_layer;
	    }
	    //Otherwise if the limit was exceeded, look at the next event
	    ++(sv.ii->curevent);

	    //If all events have been exhausted, exit the loop
	    if(sv.ii->curevent == sv.ii->nevents) break;
	    //Otherwise, generate the event time for the next event
	    sv.ii->event_time=(sv.ii-1)->event_time+sv.ii->comm_period*gsl_rng_uniform(r);
            printf("Event %i/%i at time %f\n",sv.ii->curevent,sv.ii->nevents,sv.ii->event_time);
	  }
	}

	//All events for the current individual have been exhausted
	for(;;) {
	  //Move down one layer
	  printf("Move to previous layer\n");
	  --sv.ii;

	  //If done processing, exit
	  if(sv.ii == sv.iis) goto done_parsing;

	  //Look at the next infected individual in the current event
	  ++(sv.ii->curinfection);

	  //If the infections have been exhausted
	  if(sv.ii->curinfection == sv.ii->ninfections) {
	    //Move to the next event for the individual
	    ++(sv.ii->curevent);

	    //If the events have been exhausted, go down another layer
	    if(sv.ii->curevent == sv.ii->nevents) continue;

	    //Else
	    //Generate the number of infections for the event and the associated
	    //index
	    goto gen_event;
	  }
          printf("Infection %i/%i\n",sv.ii->curinfection,sv.ii->ninfections);
	  break;
	}
      }
    }
next_pri_ii:
    ;
  }
done_parsing:

  return 0;
}
