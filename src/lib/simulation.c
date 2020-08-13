/**
 * @file simulation.c
 * @brief Simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
 * <david.maybury@tpsgc-pwgsc.gc.ca>
 */

#include "simulation.h"

#define II_ARRAY_GROW_FACT (1.5)  //!< Growing factor for the array of current infectious individuals across all layers.

int sim_pars_check(sim_pars const* pars)
{
  int ret=0;

  if(pars->tbar<=0) {
    fprintf(stderr,"%s: Error: tbar must be greater than 0\n",__func__);
    ret-=1;
  }

  if(pars->p<0) {
    fprintf(stderr,"%s: Error: p must be non-negative\n",__func__);
    ret-=2;
  }

  if(pars->lambda<=0) {
    fprintf(stderr,"%s: Error: lambda must be greater than 0\n",__func__);
    ret-=4;
  }

  if(pars->kappa<=0) {
    fprintf(stderr,"%s: Error: kappa must be greater than 0\n",__func__);
    ret-=8;
  }

  if(pars->lbar<0) {
    fprintf(stderr,"%s: Error: lbar must be non-negative\n",__func__);
    ret-=16;

  } else if(pars->lbar>0 && pars->kappal<=0) {
      fprintf(stderr,"%s: Error: kappal must be greater than 0 if lbar>0\n",__func__);
      ret-=32;
  }

  if(pars->q<0) {
    fprintf(stderr,"%s: Error: q must be non-negative\n",__func__);
    ret-=64;

  } else if(pars->q>0) {

    if(pars->mbar<0) {
      fprintf(stderr,"%s: Error: mbar must be non-negative if q>0\n",__func__);
      ret-=128;
    }

    if(pars->kappaq<=0) {
      fprintf(stderr,"%s: Error: kappaq must be greater than 0 if q>0\n",__func__);
      ret-=256;
    }
  }

  if(pars->tmax<=0) {
    fprintf(stderr,"%s: Error: tmax must be greater than 0\n",__func__);
    ret-=512;
  }

  if(pars->nstart<=0) {
    fprintf(stderr,"%s: Error: nstart must be greater than 0\n",__func__);
    ret-=1024;
  }
  return ret;
}

void sim_pars_init(sim_pars* pars)
{
  pars->tbar=NAN;
  pars->p=NAN;
  pars->lambda=NAN;
  pars->kappa=NAN;
  pars->lbar=0;
  pars->kappal=0;
  pars->q=0;
  pars->mbar=0;
  pars->kappaq=0;
  pars->R0=NAN;
  pars->mu=NAN;
  pars->t95=NAN;
  pars->m95=NAN;
  pars->l95=NAN;
  pars->tmax=INFINITY;
  pars->nstart=1;
}

int sim_init(sim_vars* sv, sim_pars* pars, const gsl_rng* r)
{
  int ret=sim_pars_check(pars);

  if(ret) return ret;

  sv->pars=*pars;
  sv->r=r;
  sv->iis=(infindividual*)malloc(INIT_N_LAYERS*sizeof(infindividual));
  sv->nlayers=INIT_N_LAYERS;

  if(sv->pars.lbar || sv->pars.kappal) sv->gen_time_periods_func=(sv->pars.q?gen_time_periods_isolation:gen_time_periods);

  else sv->gen_time_periods_func=(sv->pars.q?gen_comm_period_isolation:gen_comm_period);
  sv->iis[0].event_time=0;
  sv->dataptr=NULL;
  sv->increase_layers_proc_func=default_increase_layers_proc_func;
  sv->new_event_proc_func=default_event_proc_func;
  sv->new_inf_proc_func=dummy_proc_func_one_par;
  sv->end_inf_proc_func=dummy_proc_func_two_pars;
  sv->inf_proc_func_noevent=dummy_proc_func_two_pars;
  return 0;
}

int simulate(sim_vars* sv)
{
  int i;
  sim_pars const* sim=&(sv->pars);

  sv->ii=sv->iis;
  sv->ii->event_time=0;
  sv->ii->nevents=1;
  sv->ii->curevent=0;
  sv->ii->ninfections=sim->nstart;
  sv->new_event_proc_func(sv);

  for(i=sim->nstart-1; i>=0; --i) {
    DEBUG_PRINTF("initial individual %i\n",i);
    sv->ii=sv->iis+1;
    //Generate the communicable period appropriately
    sv->gen_time_periods_func(sv);
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
      sv->ii->event_time=sv->ii->latent_period+sv->ii->comm_period*(1-gsl_rng_uniform(sv->r));
      DEBUG_PRINTF("Event %i/%i at time %f\n",sv->ii->curevent,sv->ii->nevents,sv->ii->event_time);

      sv->ii->ninfections=gsl_ran_logarithmic(sv->r, sv->pars.p);

      if(!sv->new_event_proc_func(sv)) {
	DEBUG_PRINTF("New event returned false\n");

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
	sv->iis=(infindividual*)realloc(sv->iis,sv->nlayers*sizeof(infindividual));
	sv->increase_layers_proc_func(sv->iis+layer,sv->nlayers-layer);
	sv->ii=sv->iis+layer;
      }
      //Generate the communicable period appropriately
      sv->gen_time_periods_func(sv);

      //Generate the number of events
      sv->ii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->ii->comm_period);
      DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->ii->comm_period, sv->ii->nevents);

      //If the number of events is non-zero
      if(sv->ii->nevents) {
	sv->ii->curevent=0;
	sv->new_inf_proc_func(sv->ii);
	//Generate the event time
gen_event:
	sv->ii->event_time=(sv->ii-1)->event_time+sv->ii->latent_period+sv->ii->comm_period*(1-gsl_rng_uniform(sv->r));
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
