#include "simulation.h"

#define II_ARRAY_GROW_FACT (1.5)

void var_sim_init(struct sim_pars* sim, const struct sim_pars sim_in)
{
  sim->tbar=(sim_in.tbar?sim_in.tbar:10);
  sim->p=(sim_in.p?sim_in.p:0.5);
  sim->lambda=(sim_in.lambda?sim_in.lambda:0.11);
  sim->kappa=(sim_in.kappa?sim_in.kappa:1);
  sim->q=(sim_in.q?sim_in.q:0);
  sim->mbar=(sim_in.mbar?sim_in.mbar:5);
  sim->kappaq=(sim_in.kappaq?sim_in.kappaq:3);
  sim->nsim=(sim_in.nsim?sim_in.nsim:10);
  sim->nstart=(sim_in.nstart?sim_in.nstart:1);
  sim->tmax=(sim_in.tmax?sim_in.tmax:100);
}

void sim_vars_init(struct sim_vars* sv, const gsl_rng* r)
{
  sv->r=r;
  sv->iis=(struct infindividual*)malloc(INIT_N_LAYERS*sizeof(struct infindividual));
  sv->nlayers=INIT_N_LAYERS;
  sv->gen_comm_period_func=(sv->pars.q?gen_comm_period_isolation:gen_comm_period);
  sv->iis[0].event_time=0;
  sv->dataptr=NULL;
  sv->increase_layers_proc_func=default_increase_layers_proc_func;
  sv->new_event_proc_func=default_event_proc_func;
  sv->new_inf_proc_func=dummy_proc_func_one_par;
  sv->end_inf_proc_func=dummy_proc_func_two_pars;
  sv->inf_proc_func_noevent=dummy_proc_func_two_pars;
}

int simulate(struct sim_vars* sv)
{
  int i;
  struct sim_pars const* sim=&(sv->pars);

  for(i=sim->nstart-1; i>=0; --i) {
    DEBUG_PRINTF("initial individual %i\n",i);
    sv->ii=sv->iis+1;
    //Generate the communicable period appropriately
    sv->gen_comm_period_func(sv);
    sv->ii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->ii->comm_period);
    DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->ii->comm_period, sv->ii->nevents);

    //If no event for the current individual in the primary layer
    if(!sv->ii->nevents) {
      sv->inf_proc_func_noevent(sv->ii, sv->dataptr);
      continue;
    }

    sv->new_inf_proc_func(sv->ii);
    sv->ii->curevent=0;

    for(;;) {
      sv->ii->event_time=sv->ii->comm_period*(1-gsl_rng_uniform(sv->r));
      DEBUG_PRINTF("Event %i/%i at time %f\n",sv->ii->curevent,sv->ii->nevents,sv->ii->event_time);

      sv->ii->ninfections=gsl_ran_logarithmic(sv->r, sv->pars.p);

      if(!sv->new_event_proc_func(sv)) {

	//If the events have been exhausted, go down another layer
	if(sv->ii->curevent == sv->ii->nevents-1) {
	  sv->end_inf_proc_func(sv->ii, sv->dataptr);
	  goto done_parsing;
	}

	//Else
	//Move to the next event for the individual
	++(sv->ii->curevent);

      } else break;
    }
    sv->ii->curinfection=0;
    DEBUG_PRINTF("Infection %i/%i\n",sv->ii->curinfection,sv->ii->ninfections);

    //Create a new infected individual
    for(;;) {
      ++(sv->ii);
      DEBUG_PRINTF("Move to next layer (%li)\n",sv->ii-sv->iis);

      //If reaching the end of the allocated array, increase its size
      if(sv->ii==sv->iis+sv->nlayers) {
	sv->nlayers*=II_ARRAY_GROW_FACT;
	DEBUG_PRINTF("Growing layers to %i\n",sv->nlayers);
	uint64_t layer=sv->ii-sv->iis;
	sv->iis=(struct infindividual*)realloc(sv->iis,sv->nlayers*sizeof(struct infindividual));
	sv->increase_layers_proc_func(sv->iis+layer,sv->nlayers-layer);
	sv->ii=sv->iis+layer;
      }
      //Generate the communicable period appropriately
      sv->gen_comm_period_func(sv);

      //Generate the number of events
      sv->ii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->ii->comm_period);
      DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->ii->comm_period, sv->ii->nevents);

      //If the number of events is non-zero
      if(sv->ii->nevents) {
	sv->ii->curevent=0;
	sv->new_inf_proc_func(sv->ii);
	//Generate the event time
gen_event:
	sv->ii->event_time=(sv->ii-1)->event_time+sv->ii->comm_period*(1-gsl_rng_uniform(sv->r));
	DEBUG_PRINTF("Event %i/%i at time %f\n",sv->ii->curevent,sv->ii->nevents,sv->ii->event_time);

	//Generate the number of infections and the associated index for
	//the current event
	//Move to the next layer
	sv->ii->ninfections=gsl_ran_logarithmic(sv->r, sv->pars.p);

	if(sv->new_event_proc_func(sv)) {
	  sv->ii->curinfection=0;
	  DEBUG_PRINTF("Infection %i/%i\n",sv->ii->curinfection,sv->ii->ninfections);
	  continue;

	} else {

	  //If the events have not been exhausted
	  if(sv->ii->curevent < sv->ii->nevents-1) {
	    //Move to the next event for the individual
	    ++(sv->ii->curevent);
	    goto gen_event;
	  }
	  sv->end_inf_proc_func(sv->ii, sv->dataptr);
	}

      } else sv->inf_proc_func_noevent(sv->ii, sv->dataptr);

      //All events for the current individual have been exhausted
      for(;;) {

	if(sv->ii == sv->iis+1) goto done_parsing;
	//Move down one layer
	--(sv->ii);
	DEBUG_PRINTF("Move to previous layer (%li)\n",sv->ii-sv->iis);

	//If the infections have been exhausted
	if(sv->ii->curinfection == sv->ii->ninfections-1) {

	  //If the events have been exhausted, go down another layer
	  if(sv->ii->curevent == sv->ii->nevents-1) {
	    sv->end_inf_proc_func(sv->ii, sv->dataptr);
	    continue;
	  }

	  //Else
	  //Move to the next event for the individual
	  ++(sv->ii->curevent);
	  goto gen_event;
	}

	//Look at the next infected individual in the current event
	++(sv->ii->curinfection);
	DEBUG_PRINTF("Infection %i/%i\n",sv->ii->curinfection,sv->ii->ninfections);
	break;
      }
    }

    //Else if no event for the current infected individual
done_parsing:
    ;
  }

  return 0;
}
