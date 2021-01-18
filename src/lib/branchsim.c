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
  sv->brsim.iis[0].commpertype=0;
  sv->brsim.iis[0].nevents=1;
  sv->brsim.iis[0].cureventi=0;
  sv->brsim.iis[0].nattendees=1;
  sv->brsim.iis[0].ninfections=1;

  for(uint32_t i=0; i<INIT_N_LAYERS; ++i) sv->ii_alloc_proc_func(sv->brsim.iis+i);
  sv->brsim.nlayers=INIT_N_LAYERS;
  ran_log_init(&sv->rl, (rng_stream*)sv->r->state, sv->pars.p);

#ifdef DUAL_PINF
  BR_GENINF_COND(&& (sv->pars.ppip==0 || sv->pars.rpinfp==1));
#else
  BR_GENINF_COND();
#endif
}

int branchsim(sim_vars* sv)
{
  int i;
#ifdef CT_OUTPUT
  uint32_t npevents;
  int e;
  double end_latent_per;
  double ct_latent_overlap;
#endif
  model_pars const* sim=&(sv->pars);

  do {
    sv->path_init_proc_func(sv);
    sv->brsim.iis[0].event_time=sv->brsim.iis[1].event_time=0;

    for(i=sim->nstart-1; i>=0; --i) {
      DEBUG_PRINTF("initial individual %i\n",i);
      //Generate the communicable period appropriately
      sv->gen_pri_time_periods_func(sv, sv->brsim.iis+1, sv->brsim.iis, 0);

      sv->gen_time_origin_func(sv);
      DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",sv->brsim.iis[1].latent_period,sv->brsim.iis[1].comm_period,sv->brsim.iis[1].commpertype,sv->brsim.iis[1].end_comm_period);

      sv->curii=sv->brsim.iis;
      sv->pri_init_proc_func(sv, sv->brsim.iis+1);
      sv->curii=sv->brsim.iis+1;
      DEBUG_PRINTF("Move to primary layer (%li)\n",sv->curii-sv->brsim.iis);

#ifdef CT_OUTPUT
#define GEN_LATENT_CONTACTS \
      npevents=0; \
      \
      /*If the CT window starts before the communicable period, generate */ \
      /*pseudo-events and calculate the number of traced contacts */ \
      if(sv->curii->comm_period<sim->ctwindow) { \
	ct_latent_overlap=sim->ctwindow-sv->curii->comm_period; \
	npevents=gsl_ran_poisson(sv->r, sim->lambda*ct_latent_overlap); \
	DEBUG_PRINTF("Number of pre-events is %u during %f\n",npevents,ct_latent_overlap); \
	\
	if(npevents) { \
	  sv->new_inf_proc_func(sv, sv->curii, sv->curii-1); \
	  sv->curii->ninfections=0; \
	  end_latent_per=sv->curii->end_comm_period-sv->curii->comm_period; \
	  \
	  for(e=npevents-1; e>=0; --e) { \
	    sv->curii->event_time=end_latent_per-ct_latent_overlap*gsl_rng_uniform(sv->r); \
	    sv->curii->nattendees=sv->gen_att_func(sv); \
	    sv->curii->ntracednicts=gsl_ran_binomial(sv->r, sim->pt, sv->curii->nattendees-1); \
	    DEBUG_PRINTF("%u attendees, %u successfully traced contacts were generated\n",sv->curii->nattendees,sv->curii->ntracednicts); \
	    sv->new_event_proc_func(sv); \
	  } \
	} \
      }
      GEN_LATENT_CONTACTS;
#endif

      sv->curii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->curii->comm_period);
      DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->curii->comm_period, sv->curii->nevents);

      //If no event for the current individual in the primary layer
      if(!sv->curii->nevents) {
#ifdef CT_OUTPUT
	if(!npevents) sv->new_inf_proc_func_noevent(sv, sv->curii, sv->brsim.iis);
	else sv->end_inf_proc_func(sv, sv->curii, sv->curii-1);
#else
	sv->new_inf_proc_func_noevent(sv, sv->curii, sv->brsim.iis);
#endif
	continue;
      }

#ifdef CT_OUTPUT
      if(!npevents) sv->new_inf_proc_func(sv, sv->curii, sv->brsim.iis);
#else
      sv->new_inf_proc_func(sv, sv->curii, sv->brsim.iis);
#endif
      sv->curii->cureventi=0;

      for(;;) {
	sv->curii->event_time=sv->curii->end_comm_period-sv->curii->comm_period*gsl_rng_uniform(sv->r);
	//sv->curii->event_time=sv->curii->latent_period+sv->curii->comm_period*rng_rand_pu01d((rng_stream*)sv->r->state);
	DEBUG_PRINTF("Event %i/%i at time %f\n",sv->curii->cureventi,sv->curii->nevents,sv->curii->event_time);

	//sv->curii->ninfections=gsl_ran_logarithmic(sv->r, sim->p);
	sv->gen_att_inf_func(sv);
#ifdef CT_OUTPUT
#define GEN_CONTACTS_AND_TRACE \
	if((sv->curii->commpertype&ro_commper_true_positive_test) && sv->curii->event_time>=sv->curii->end_comm_period-sim->ctwindow) { \
	  sv->curii->ntracednicts=gsl_ran_binomial(sv->r, sim->pt, sv->curii->nattendees-1-sv->curii->ninfections); \
	  \
	  if(sv->curii->ninfections) sv->curii->ntracedicts=gsl_ran_binomial(sv->r, sim->pt, sv->curii->ninfections); \
	  else sv->curii->ntracedicts=0; \
	  sv->curii->gen_ct_time_periods_func=sv->gen_time_periods_func; \
	  \
	} else { \
	  sv->curii->ntracednicts=sv->curii->ntracedicts=0; \
	  sv->curii->gen_ct_time_periods_func=sv->gen_time_periods_func_no_int; \
	}
	GEN_CONTACTS_AND_TRACE;
	DEBUG_PRINTF("%u attendees, %u infections, %u / %u non-infected/infected successfully traced contacts were generated (%f)\n",sv->curii->nattendees,sv->curii->ninfections,sv->curii->ntracednicts,sv->curii->ntracedicts,sv->curii->event_time-(sv->curii->end_comm_period-sim->ctwindow));
#else
	DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curii->nattendees,sv->curii->ninfections);
#endif

	if(!sv->new_event_proc_func(sv)) {
	  DEBUG_PRINTF("New event returned false\n");

	  //If the events have been exhausted, go down another layer
	  if(sv->curii->cureventi == sv->curii->nevents-1) {
	    sv->end_inf_proc_func(sv, sv->curii, sv->curii-1);
	    goto done_parsing;
	  }

	  //Else
	  //Move to the next event for the individual
	  ++(sv->curii->cureventi);

	} else break;
      }
      sv->curii->curinfectioni=0;
      DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfectioni,sv->curii->ninfections);

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
#ifdef CT_OUTPUT
	(sv->curii-1)->gen_ct_time_periods_func(sv, sv->curii, sv->curii-1, (sv->curii-1)->event_time);
#else
	sv->gen_time_periods_func(sv, sv->curii, sv->curii-1, (sv->curii-1)->event_time);
#endif
	DEBUG_PRINTF("Event time: %f, latent period is %f, comm period is %f, type is %u, end comm is %f\n",(sv->curii-1)->event_time,sv->curii->latent_period,sv->curii->comm_period,sv->curii->commpertype,sv->curii->end_comm_period);

#ifdef CT_OUTPUT
	GEN_LATENT_CONTACTS;
#endif

	//Generate the number of events
	sv->curii->nevents=gsl_ran_poisson(sv->r, sim->lambda*sv->curii->comm_period);
	DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, sv->curii->comm_period, sv->curii->nevents);

	//If the number of events is non-zero
	if(sv->curii->nevents) {
	  sv->curii->cureventi=0;

#ifdef CT_OUTPUT
	  if(!npevents) sv->new_inf_proc_func(sv, sv->curii, sv->curii-1);
#else
	  sv->new_inf_proc_func(sv, sv->curii, sv->curii-1);
#endif
	  //Generate the event time
gen_event:
	  sv->curii->event_time=sv->curii->end_comm_period-sv->curii->comm_period*gsl_rng_uniform(sv->r);
	  //sv->curii->event_time=(sv->curii-1)->event_time+sv->curii->latent_period+sv->curii->comm_period*rng_rand_pu01d((rng_stream*)sv->r->state);
	  DEBUG_PRINTF("Event %i/%i at time %f\n",sv->curii->cureventi,sv->curii->nevents,sv->curii->event_time);

	  //Generate the number of infections and the associated index for
	  //the current event
	  //Move to the next layer
	  //sv->curii->ninfections=gsl_ran_logarithmic(sv->r, sim->p);
	  sv->gen_att_inf_func(sv);
#ifdef CT_OUTPUT
	  GEN_CONTACTS_AND_TRACE;
	  DEBUG_PRINTF("%u attendees, %u infections, %u / %u non-infected/infected successfully traced contacts were generated (%f)\n",sv->curii->nattendees,sv->curii->ninfections,sv->curii->ntracednicts,sv->curii->ntracedicts,sv->curii->event_time-(sv->curii->end_comm_period-sim->ctwindow));
#else
	  DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curii->nattendees,sv->curii->ninfections);
#endif

	  if(sv->new_event_proc_func(sv)) {
	    sv->curii->curinfectioni=0;
	    DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfectioni,sv->curii->ninfections);
	    continue;

	  } else {

	    //If the events have not been exhausted
	    if(sv->curii->cureventi < sv->curii->nevents-1) {
	      //Move to the next event for the individual
	      ++(sv->curii->cureventi);
	      goto gen_event;
	    }
	    sv->end_inf_proc_func(sv, sv->curii, sv->curii-1);
	  }

	} else {
#ifdef CT_OUTPUT
	  if(!npevents) sv->new_inf_proc_func_noevent(sv, sv->curii, sv->brsim.iis);
	  else sv->end_inf_proc_func(sv, sv->curii, sv->curii-1);
#else
	  sv->new_inf_proc_func_noevent(sv, sv->curii, sv->curii-1);
#endif
	}

	//All events for the current individual have been exhausted
	for(;;) {

	  if(sv->curii == sv->brsim.iis+1) goto done_parsing;
	  //Move down one layer
	  --(sv->curii);
	  DEBUG_PRINTF("Move to previous layer (%li)\n",sv->curii-sv->brsim.iis);

	  //If the infections have been exhausted
	  if(sv->curii->curinfectioni == sv->curii->ninfections-1) {

	    //If the events have been exhausted, go down another layer
	    if(sv->curii->cureventi == sv->curii->nevents-1) {
	      sv->end_inf_proc_func(sv, sv->curii, sv->curii-1);
	      continue;
	    }

	    //Else
	    //Move to the next event for the individual
	    ++(sv->curii->cureventi);
	    goto gen_event;
	  }

	  //Look at the next infected individual in the current event
	  ++(sv->curii->curinfectioni);
	  DEBUG_PRINTF("Infection %i/%i\n",sv->curii->curinfectioni,sv->curii->ninfections);
	  break;
	}
      }

      //Else if no event for the current primary infected individual
done_parsing:
      ;
    }

  } while(!sv->path_end_proc_func(sv));

  return 0;
}

void branchsim_free(sim_vars* sv)
{
  for(uint32_t i=0; i<sv->brsim.nlayers; ++i) if(sv->brsim.iis[i].dataptr) free(sv->brsim.iis[i].dataptr);
  free(sv->brsim.iis);
}
