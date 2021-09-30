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
  sv->brsim.layers=(inflayer*)malloc(INIT_N_LAYERS*sizeof(inflayer));
  sv->brsim.layers[0].ii.commpertype=0;
  sv->brsim.layers[0].nevents=1;
  sv->brsim.layers[0].cureventi=0;
  sv->brsim.layers[0].ii.nattendees=1;
  sv->brsim.layers[0].ii.ninfections=1;

  for(uint32_t i=0; i<INIT_N_LAYERS; ++i) {
    sv->brsim.layers[i].ii.generation=i;
    sv->ii_alloc_proc_func(&sv->brsim.layers[i].ii);
  }
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
  brsim_vars* brsim=&sv->brsim;
  inflayer* curlayer;

#ifdef DUAL_PINF
  const double pinfpinf=sv->pars.ppip*sv->pars.rpinfp/(1+sv->pars.ppip*(sv->pars.rpinfp-1));
#endif

  do {
    sv->path_init_proc_func(sv);
    sv->event_time=0;

    for(i=sim->nstart-1; i>=0; --i) {
      DEBUG_PRINTF("initial individual %i\n",i);

      #ifdef DUAL_PINF
      if(gsl_rng_uniform(sv->r) < pinfpinf) {
        #ifdef SEC_INF_TIMELINES
        brsim->layers[0].ii.ninfectionsf=0;
        brsim->layers[0].ii.ninfectionsp=1;
        #endif
	brsim->layers[1].ii.inftypep=true;
	brsim->layers[1].ii.q=sv->pars.qp;
	brsim->layers[1].ii.pinf=sv->pars.pinf*sv->pars.rpshedp;

      } else {
        #ifdef SEC_INF_TIMELINES
        brsim->layers[0].ii.ninfectionsf=1;
        brsim->layers[0].ii.ninfectionsp=0;
        #endif
	brsim->layers[1].ii.inftypep=false;
	brsim->layers[1].ii.q=sv->pars.q;
	brsim->layers[1].ii.pinf=sv->pars.pinf;
      }
      #endif
      //Generate the communicable period appropriately
      sv->gen_pri_time_periods_func(sv, &brsim->layers[1].ii, &brsim->layers[0].ii, 0);

      sv->gen_time_origin_func(sv, &brsim->layers[1].ii);
      DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",brsim->layers[1].ii.latent_period,brsim->layers[1].ii.comm_period,brsim->layers[1].ii.commpertype,brsim->layers[1].ii.end_comm_period);

      sv->pri_init_proc_func(sv, &brsim->layers[0].ii, &brsim->layers[1].ii);
      curlayer=brsim->layers+1;
      DEBUG_PRINTF("Move to primary layer (%li)\n",curlayer->ii.generation);

#ifdef CT_OUTPUT
#define GEN_LATENT_CONTACTS \
      npevents=0; \
      \
      /*If the CT window starts before the communicable period, generate */ \
      /*pseudo-events and calculate the number of traced contacts */ \
      if(curlayer->ii.comm_period<sim->ctwindow) { \
	ct_latent_overlap=sim->ctwindow-curlayer->ii.comm_period; \
	npevents=gsl_ran_poisson(sv->r, sim->lambda*ct_latent_overlap); \
	DEBUG_PRINTF("Number of pre-events is %u during %f\n",npevents,ct_latent_overlap); \
	\
	if(npevents) { \
	  sv->new_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii); \
	  curlayer->ii.ninfections=0; \
	  end_latent_per=curlayer->ii.end_comm_period-curlayer->ii.comm_period; \
	  \
	  for(e=npevents-1; e>=0; --e) { \
	    sv->event_time=end_latent_per-ct_latent_overlap*gsl_rng_uniform(sv->r); \
	    curlayer->ii.nattendees=sv->gen_att_func(sv); \
	    curlayer->ii.ntracednicts=gsl_ran_binomial(sv->r, sim->pt, curlayer->ii.nattendees-1); \
	    DEBUG_PRINTF("%u attendees, %u successfully traced contacts were generated\n",curlayer->ii.nattendees,curlayer->ntracednicts); \
	    sv->new_event_proc_func(sv, &curlayer->ii); \
	  } \
	} \
      }
      GEN_LATENT_CONTACTS;
#endif

      curlayer->nevents=gsl_ran_poisson(sv->r, sim->lambda*curlayer->ii.comm_period);
      DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, curlayer->ii.comm_period, curlayer->nevents);

      //If no event for the current individual in the primary layer
      if(!curlayer->nevents) {
#ifdef CT_OUTPUT
	if(!npevents) sv->new_inf_proc_func_noevent(sv, &curlayer->ii, &brsim->layers[0].ii);
	else sv->end_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
#else
	sv->new_inf_proc_func_noevent(sv, &curlayer->ii, &brsim->layers[0].ii);
#endif
	continue;
      }

#ifdef CT_OUTPUT
      if(!npevents) sv->new_inf_proc_func(sv, &curlayer->ii, &brsim->layers[0].ii);
#else
      sv->new_inf_proc_func(sv, &curlayer->ii, &brsim->layers[0].ii);
#endif
      curlayer->cureventi=0;

      for(;;) {
	sv->event_time=curlayer->ii.end_comm_period-curlayer->ii.comm_period*gsl_rng_uniform(sv->r);
	//sv->event_time=curlayer->ii.latent_period+curlayer->ii.comm_period*rng_rand_pu01d((rng_stream*)sv->r->state);
	DEBUG_PRINTF("Event %i/%i at time %f\n",curlayer->cureventi,curlayer->nevents,sv->event_time);

	//curlayer->ninfections=gsl_ran_logarithmic(sv->r, sim->p);
	brsim->gen_att_inf_func(sv, &curlayer->ii);
#ifdef CT_OUTPUT
#define GEN_CONTACTS_AND_TRACE \
	if((curlayer->ii.commpertype&ro_commper_true_positive_test) && sv->event_time>=curlayer->ii.end_comm_period-sim->ctwindow) { \
	  curlayer->ii.ntracednicts=gsl_ran_binomial(sv->r, sim->pt, curlayer->ii.nattendees-1-curlayer->ii.ninfections); \
	  \
	  if(curlayer->ii.ninfections) curlayer->ii.ntracedicts=gsl_ran_binomial(sv->r, sim->pt, curlayer->ii.ninfections); \
	  else curlayer->ii.ntracedicts=0; \
	  curlayer->ii.gen_ct_time_periods_func=sv->gen_time_periods_func; \
	  \
	} else { \
	  curlayer->ii.ntracednicts=curlayer->ii.ntracedicts=0; \
	  curlayer->ii.gen_ct_time_periods_func=sv->gen_time_periods_func_no_int; \
	}
	GEN_CONTACTS_AND_TRACE;
	DEBUG_PRINTF("%u attendees, %u infections, %u / %u non-infected/infected successfully traced contacts were generated (%f)\n",curlayer->ii.nattendees,curlayer->i.ninfections,curlayer->i.ntracednicts,curlayer->i.ntracedicts,sv->event_time-(curlayer->ii.end_comm_period-sim->ctwindow));
#else
	DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curevent.nattendees,curlayer->ninfections);
#endif

	if(!sv->new_event_proc_func(sv, &curlayer->ii)) {
	  DEBUG_PRINTF("New event returned false\n");

	  //If the events have been exhausted, go down another layer
	  if(curlayer->cureventi == curlayer->nevents-1) {
	    sv->end_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
	    goto done_parsing;
	  }

	  //Else
	  //Move to the next event for the individual
	  ++(curlayer->cureventi);

	} else break;
      }
      curlayer->curinfectioni=0;
      DEBUG_PRINTF("Infection %i/%i\n",curlayer->curinfectioni,curlayer->ninfections);

      //Create a new infected individual
      for(;;) {
	++(curlayer);
	DEBUG_PRINTF("Move to next layer (%li)\n",curlayer->ii.layer);

	//If reaching the end of the allocated array, increase its size
	if(curlayer->ii.generation==brsim->nlayers-1) {
	  brsim->nlayers*=II_ARRAY_GROW_FACT;
	  DEBUG_PRINTF("Growing layers to %i\n",brsim->nlayers);
	  uint32_t layer=curlayer->ii.generation;
	  brsim->layers=(inflayer*)realloc(brsim->layers,brsim->nlayers*sizeof(inflayer));

	  for(uint32_t i=brsim->nlayers-1; i>=layer; --i) {
	    brsim->layers[i].ii.generation=i;
	    sv->ii_alloc_proc_func(&brsim->layers[i].ii);
	  }
	  curlayer=brsim->layers+layer;
	}

      #ifdef DUAL_PINF
	//The number of infections is known, so the infection category must be
	//randomly assigned based on the counts of remaining individuals for
	//each category.
	if(gsl_rng_uniform(sv->r) < ((double)(curlayer-1)->ii.ninfectionsp) / ((curlayer-1)->ii.ninfectionsf + (curlayer-1)->ii.ninfectionsp)) {
	  --((curlayer-1)->ii.ninfectionsp);
	  curlayer->ii.inftypep=true;
	  curlayer->ii.q=sv->pars.qp;
	  curlayer->ii.pinf=sv->pars.pinf*sv->pars.rpshedp;

	} else {
	  --((curlayer-1)->ii.ninfectionsf);
	  curlayer->ii.inftypep=false;
	  curlayer->ii.q=sv->pars.q;
	  curlayer->ii.pinf=sv->pars.pinf;
	}
        #endif
	//Generate the communicable period appropriately
#ifdef CT_OUTPUT
        //We don't need to draw a random number to find which infection indices can be traced since all infections are drawn independently. It is thus possible to compare the infection index to the number of successfully traced infection contacts
	if((curlayer-1)->curinfectioni < (curlayer-1)->ii.ntracedicts) curlayer->ii.traced=true;

	else curlayer->ii.traced=false;

	(curlayer-1)->ii.gen_ct_time_periods_func(sv, &curlayer->ii, &(curlayer-1)->ii, sv->event_time);
#else
	sv->gen_time_periods_func(sv, &curlayer->ii, &(curlayer-1)->ii, sv->event_time);
#endif
	DEBUG_PRINTF("Event time: %f, latent period is %f, comm period is %f, type is %u, end comm is %f\n",sv->event_time,curlayer->ii.latent_period,curlayer->ii.comm_period,curlayer->ii.commpertype,curlayer->ii.end_comm_period);

#ifdef CT_OUTPUT
	GEN_LATENT_CONTACTS;
#endif

	//Generate the number of events
	curlayer->nevents=gsl_ran_poisson(sv->r, sim->lambda*curlayer->ii.comm_period);
	DEBUG_PRINTF("Nevents (%f*%f) is %i\n", sim->lambda, curlayer->ii.comm_period, curlayer->nevents);

	//If the number of events is non-zero
	if(curlayer->nevents) {
	  curlayer->cureventi=0;

#ifdef CT_OUTPUT
	  if(!npevents) sv->new_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
#else
	  sv->new_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
#endif
	  //Generate the event time
gen_event:
	  sv->event_time=curlayer->ii.end_comm_period-curlayer->ii.comm_period*gsl_rng_uniform(sv->r);
	  DEBUG_PRINTF("Event %i/%i at time %f\n",curlayer->cureventi,curlayer->nevents,sv->event_time);

	  //Generate the number of infections and the associated index for
	  //the current event
	  //Move to the next layer
	  //curlayer->ninfections=gsl_ran_logarithmic(sv->r, sim->p);
	  brsim->gen_att_inf_func(sv, &curlayer->ii);
#ifdef CT_OUTPUT
	  GEN_CONTACTS_AND_TRACE;
	  DEBUG_PRINTF("%u attendees, %u infections, %u / %u non-infected/infected successfully traced contacts were generated (%f)\n",sv->curevent.nattendees,curlayer->ninfections,curlayer->ntracednicts,curlayer->ntracedicts,sv->event_time-(curlayer->ii.end_comm_period-sim->ctwindow));
#else
	  DEBUG_PRINTF("%u attendees and %u infections were generated\n",sv->curevent.nattendees,curlayer->ninfections);
#endif

	  if(sv->new_event_proc_func(sv, &curlayer->ii)) {
	    curlayer->curinfectioni=0;
	    DEBUG_PRINTF("Infection %i/%i\n",curlayer->curinfectioni,curlayer->ninfections);
	    continue;

	  } else {

	    //If the events have not been exhausted
	    if(curlayer->cureventi < curlayer->nevents-1) {
	      //Move to the next event for the individual
	      ++(curlayer->cureventi);
	      goto gen_event;
	    }
	    sv->end_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
	  }

	} else {
#ifdef CT_OUTPUT
	  if(!npevents) sv->new_inf_proc_func_noevent(sv, &curlayer->ii, &brsim->layers[0].ii);
	  else sv->end_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
#else
	  sv->new_inf_proc_func_noevent(sv, &curlayer->ii, &(curlayer-1)->ii);
#endif
	}

	//All events for the current individual have been exhausted
	for(;;) {

	  if(curlayer->ii.generation == 1) goto done_parsing;
	  //Move down one layer
	  --(curlayer);
	  DEBUG_PRINTF("Move to previous layer (%li)\n",curlayer->ii.generation);

	  //If the infections have been exhausted
	  if(curlayer->curinfectioni == curlayer->ii.ninfections-1) {

	    //If the events have been exhausted, go down another layer
	    if(curlayer->cureventi == curlayer->nevents-1) {
	      sv->end_inf_proc_func(sv, &curlayer->ii, &(curlayer-1)->ii);
	      continue;
	    }

	    //Else
	    //Move to the next event for the individual
	    ++(curlayer->cureventi);
	    goto gen_event;
	  }

	  //Look at the next infected individual in the current event
	  ++(curlayer->curinfectioni);
	  DEBUG_PRINTF("Infection %i/%i\n",curlayer->curinfectioni,curlayer->ninfections);
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
  for(uint32_t i=0; i<sv->brsim.nlayers; ++i) if(sv->brsim.layers[i].ii.dataptr) free(sv->brsim.layers[i].ii.dataptr);
  free(sv->brsim.layers);
}
