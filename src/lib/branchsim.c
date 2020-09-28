/**
 * @file branchsim.c
 * @brief Branching simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Original model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
 * <david.maybury@tpsgc-pwgsc.gc.ca>
 */

#include "branchsim.h"

void branchsim_init(sim_vars* sv)
{
  sv->brsim.iis=(infindividual*)malloc(INIT_N_LAYERS*sizeof(infindividual));

  for(uint32_t i=0; i<INIT_N_LAYERS; ++i) sv->ii_alloc_proc_func(sv->brsim.iis+i);
  sv->brsim.nlayers=INIT_N_LAYERS;
  ran_log_init(&sv->rl, (rng_stream*)sv->r->state, sv->pars.p);

  BR_GENINF_COND;
}

int branchsim(sim_vars* sv)
{
  int i;
  model_pars const* sim=&(sv->pars);

  sv->brsim.iis[0].nevents=1;
  sv->brsim.iis[0].curevent=0;
  sv->brsim.iis[0].nattendees=1;
  sv->brsim.iis[0].ninfections=1;
  sv->brsim.iis[0].event_time=sv->brsim.iis[1].event_time=0;

  for(i=sim->nstart-1; i>=0; --i) {
    DEBUG_PRINTF("initial individual %i\n",i);
    //Generate the communicable period appropriately
    sv->gen_pri_time_periods_func(sv, sv->brsim.iis+1, sv->brsim.iis, 0);

    sv->gen_time_origin_func(sv);
    sv->brsim.iis[1].commpertype|=ro_commper_tmax*(sv->brsim.iis[1].end_comm_period > sv->pars.tmax);

    DEBUG_PRINTF("Comm period is %f%s\n",sv->brsim.iis[1].comm_period,(sv->brsim.iis[1].commpertype&ro_commper_tmax?" (reached end)":"")); \

    sv->pri_init_proc_func(sv, sv->brsim.iis+1);

    sv->curii=sv->brsim.iis;
    sv->new_event_proc_func(sv);
    sv->curii=sv->brsim.iis+1;
    DEBUG_PRINTF("Move to primary layer (%li)\n",sv->curii-sv->brsim.iis);

    sv->curii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->curii->comm_period);
    DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->curii->comm_period, sv->curii->nevents);

    //If no event for the current individual in the primary layer
    if(!sv->curii->nevents) {
      sv->new_inf_proc_func_noevent(sv, sv->curii, sv->brsim.iis);
      continue;
    }
    sv->new_inf_proc_func(sv, sv->curii, sv->brsim.iis);

    sv->curii->curevent=0;

    for(;;) {
      sv->curii->event_time=sv->curii->end_comm_period-sv->curii->comm_period*gsl_rng_uniform(sv->r);
      //sv->curii->event_time=sv->curii->latent_period+sv->curii->comm_period*rng_rand_pu01d((rng_stream*)sv->r->state);
      DEBUG_PRINTF("Event %i/%i at time %f\n",sv->curii->curevent,sv->curii->nevents,sv->curii->event_time);

      //sv->curii->ninfections=gsl_ran_logarithmic(sv->r, sv->pars.p);
      sv->gen_att_inf_func(sv);
      DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curii->nattendees,sv->curii->ninfections);

      if(sv->curii->nattendees<2 || !sv->new_event_proc_func(sv)) {
	DEBUG_PRINTF("New event returned false\n");

	//If the events have been exhausted, go down another layer
	if(sv->curii->curevent == sv->curii->nevents-1) {
	  sv->end_inf_proc_func(sv->curii, sv->dataptr);
	  goto done_parsing;
	}

	//Else
	//Move to the next event for the individual
	++(sv->curii->curevent);

      } else break;
    }
    sv->curii->curinfection=0;
    DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfection,sv->curii->ninfections);

    //Create a new infected individual
    for(;;) {
      ++(sv->curii);
      DEBUG_PRINTF("Move to next layer (%li)\n",sv->curii-sv->brsim.iis);

      //If reaching the end of the allocated array, increase its size
      if(sv->curii==sv->brsim.iis+sv->brsim.nlayers) {
	sv->brsim.nlayers*=II_ARRAY_GROW_FACT;
	DEBUG_PRINTF("Growing layers to %i\n",sv->brsim.nlayers);
	uint32_t layer=sv->curii-sv->brsim.iis;
	sv->brsim.iis=(infindividual*)realloc(sv->brsim.iis,sv->brsim.nlayers*sizeof(infindividual));

	for(uint32_t i=sv->brsim.nlayers-1; i>=layer; --i) sv->ii_alloc_proc_func(sv->brsim.iis+i);
	sv->curii=sv->brsim.iis+layer;
      }
      //Generate the communicable period appropriately
      sv->gen_time_periods_func(sv, sv->curii, sv->curii-1, (sv->curii-1)->event_time);
      sv->curii->commpertype|=ro_commper_tmax*(sv->curii->end_comm_period > sv->pars.tmax);
      DEBUG_PRINTF("Comm period is %f, type is 0x%08x, end comm is %f%s\n",sv->curii->comm_period,sv->curii->commpertype,sv->curii->end_comm_period,(sv->curii->commpertype&ro_commper_tmax?" (reached end)":"")); \

      //Generate the number of events
      sv->curii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->curii->comm_period);
      DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->curii->comm_period, sv->curii->nevents);

      //If the number of events is non-zero
      if(sv->curii->nevents) {
	sv->curii->curevent=0;
	sv->new_inf_proc_func(sv, sv->curii, sv->curii-1);
	//Generate the event time
gen_event:
        sv->curii->event_time=sv->curii->end_comm_period-sv->curii->comm_period*gsl_rng_uniform(sv->r);
	//sv->curii->event_time=(sv->curii-1)->event_time+sv->curii->latent_period+sv->curii->comm_period*rng_rand_pu01d((rng_stream*)sv->r->state);
	DEBUG_PRINTF("Event %i/%i at time %f\n",sv->curii->curevent,sv->curii->nevents,sv->curii->event_time);

	//Generate the number of infections and the associated index for
	//the current event
	//Move to the next layer
	//sv->curii->ninfections=gsl_ran_logarithmic(sv->r, sv->pars.p);
        sv->gen_att_inf_func(sv);
        DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curii->nattendees,sv->curii->ninfections);

	if(sv->curii->nattendees>1 && sv->new_event_proc_func(sv)) {
	  sv->curii->curinfection=0;
	  DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfection,sv->curii->ninfections);
	  continue;

	} else {

	  //If the events have not been exhausted
	  if(sv->curii->curevent < sv->curii->nevents-1) {
	    //Move to the next event for the individual
	    ++(sv->curii->curevent);
	    goto gen_event;
	  }
	  sv->end_inf_proc_func(sv->curii, sv->dataptr);
	}

      } else sv->new_inf_proc_func_noevent(sv, sv->curii, sv->curii-1);

      //All events for the current individual have been exhausted
      for(;;) {

	if(sv->curii == sv->brsim.iis+1) goto done_parsing;
	//Move down one layer
	--(sv->curii);
	DEBUG_PRINTF("Move to previous layer (%li)\n",sv->curii-sv->brsim.iis);

	//If the infections have been exhausted
	if(sv->curii->curinfection == sv->curii->ninfections-1) {

	  //If the events have been exhausted, go down another layer
	  if(sv->curii->curevent == sv->curii->nevents-1) {
	    sv->end_inf_proc_func(sv->curii, sv->dataptr);
	    continue;
	  }

	  //Else
	  //Move to the next event for the individual
	  ++(sv->curii->curevent);
	  goto gen_event;
	}

	//Look at the next infected individual in the current event
	++(sv->curii->curinfection);
	DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfection,sv->curii->ninfections);
	break;
      }
    }

    //Else if no event for the current primary infected individual
done_parsing:
    ;
  }

  return 0;
}

void branchsim_free(sim_vars* sv)
{
  for(uint32_t i=0; i<sv->brsim.nlayers; ++i) if(sv->brsim.iis[i].dataptr) free(sv->brsim.iis[i].dataptr);
  free(sv->brsim.iis);
}
