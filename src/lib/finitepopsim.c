/**
 * @file finitepopsim.c
 * @brief Finite population simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "finitepopsim.h"

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers
#define II_ARRAY_GROW_FACT (1.5)  //!< Growing factor for the array of current infectious individuals across all layers.

enum ro_fp_sus_state{ro_fp_sus_state_uninitialized=0, ro_fp_sus_state_high_level=1, ro_fp_sus_state_detailed=2};

void finitepopsim_init(sim_vars* sv)
{
  sv->fpsim.rooti.ii.commpertype=0;
  sv->fpsim.rooti.ii.nevents=1;
  sv->fpsim.rooti.ii.generation=0;

  sv->fpsim.is=(individual*)malloc(sv->pars.popsize*sizeof(individual));

  sv->ii_alloc_proc_func(&sv->fpsim.rooti.ii);
  for(int32_t i=sv->pars.popsize-1; i>=0; --i) sv->ii_alloc_proc_func(&sv->fpsim.is[i].ii);
  ran_log_init(&sv->rl, (rng_stream*)sv->r->state, sv->pars.p);

  sv->fpsim.ctx_activated=sarray_init((void*)&sv->fpsim.activated, &sv->fpsim.nactivated, sizeof(individual*), 0);
  sarray_alloc(sv->fpsim.ctx_activated, sv->pars.popsize);
  sv->fpsim.ctx_einfectious=sarray_init((void*)&sv->fpsim.einfectious, &sv->fpsim.neinfectious, sizeof(individual*), 0);
  sarray_alloc(sv->fpsim.ctx_einfectious, sv->pars.popsize);

  FP_GENINF_COND();
}

int finitepopsim(sim_vars* sv)
{
  int i,j,k;
  model_pars const* sim=&(sv->pars);
  fpsim_vars* fpsim=&sv->fpsim;
  void* ptr;
  individual* ind;
  uint32_t neinvitees;   //Number of invitees at the event
  uint32_t nepainvitees; //Number of potentially activated invitees at the event
  uint32_t nsusceptibles; //Total number of remaining susceptible individuals within the whole population
  uint32_t nesusceptibles; //Total number of susceptible individuals at the event
#ifdef DUAL_PINF
  uint32_t nsusceptiblesf; //Number of remaining susceptible individuals of the first category within the whole population
  uint32_t nsusceptiblesp; //Number of remaining susceptible individuals of the second category within the whole population
  uint32_t nesusceptiblesf; //Number of susceptible individuals of the first category at the event
  uint32_t nesusceptiblesp; //Number of susceptible individuals of the second category at the event
  uint32_t neinfectionsf;   //Number of infected individuals of the first category at the event
  uint32_t neinfectionsp;   //Number of infected individuals of the second category at the event
  double epninf;            //Probability of non-infection at the event for the first category of infected individual
  double etpinf;            //Sum of the probabilities of infections at the event

  const double ppinf=sim->pinf*sim->rpshedp;
  const double pinfpinf=sim->ppip*sim->rpinfp/(1+sim->ppip*(sim->rpinfp-1));
#else
  uint32_t neinfections;  //Total number of infections at the event
#endif
  bool wasinfectious;
  bool positiveevents;
  uint8_t susstate;       //State for the susceptible individual information

  //Nstart individuals are assumed to be infected when introduced to the
  //simulation

  do {
    ptr=fpsim->rooti.ii.dataptr;
    memset(&fpsim->rooti, 0, sizeof(individual));
    fpsim->rooti.ii.dataptr=ptr;
    sarray_empty(fpsim->ctx_activated);
    susstate=ro_fp_sus_state_uninitialized;

    for(i=sim->popsize-1; i>=0; --i) {
      ptr=fpsim->is[i].ii.dataptr;
      memset(fpsim->is+i, 0, sizeof(individual));
      fpsim->is[i].ii.dataptr=ptr;
    }

    sv->path_init_proc_func(sv);
    sv->event_time=0;

    for(i=sim->nstart-1; i>=0; --i) {
      DEBUG_PRINTF("initial individual %i\n",i);

#ifdef DUAL_PINF
      if(gsl_rng_uniform(sv->r) < pinfpinf) {
#ifdef SEC_INF_TIMELINES
	++fpsim->rooti.ii.ninfectionsp;
#endif
	fpsim->is[i].ii.inftypep=true;
	fpsim->is[i].ii.q=sim->qp;
	fpsim->is[i].ii.pinf=ppinf;

      } else {
#ifdef SEC_INF_TIMELINES
	++fpsim->rooti.ii.ninfectionsf;
#endif
	fpsim->is[i].ii.inftypep=false;
	fpsim->is[i].ii.q=sim->q;
	fpsim->is[i].ii.pinf=sim->pinf;
      }
#endif
      fpsim->is[i].parent=&fpsim->rooti;
      fpsim->is[i].ii.generation=1;

      //Generate the communicable period appropriately
      sv->gen_pri_time_periods_func(sv, &fpsim->is[i].ii, &fpsim->rooti.ii, 0);

      sv->gen_time_origin_func(sv, &fpsim->is[i].ii);
      DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",fpsim->is[i].ii.latent_period,fpsim->is[i].ii.comm_period,fpsim->is[i].ii.commpertype,fpsim->is[i].ii.end_comm_period);

      ind_init_next_change_time(fpsim->is+i);
      fpsim->activated[fpsim->nactivated++]=fpsim->is+i;

      sv->pri_init_proc_func(sv, &fpsim->rooti.ii, &fpsim->is[i].ii);
    }

    //Susceptible array is not up to date

    //Event loop for the current path
    for(;;) {
      sv->event_time+=gsl_ran_exponential(sv->r, 1./sim->lambdap);
      DEBUG_PRINTF("Event at time %f\n",sv->event_time);
      //Activated array is not up to date

      //Generate the number of invitees for the current event
      neinvitees=sv->gen_att_func(sv);

      //Generate a number of previously known activated individuals amongst
      //these invitees. The drawn is without replacement here, so it is a
      //hypergeometric distribution and not a binomial
      nepainvitees=gsl_ran_hypergeometric(sv->r, fpsim->nactivated, sim->popsize, neinvitees);

      //If the potential number of activated individuals is non-zero
      if(nepainvitees>0) {
	//Identify the infectious individuals amongst the previously known
	//activated individuals, updating the activated array at the same time
	fpsim->neinfectious=0;
#ifdef DUAL_PINF
	epninf=1;
	etpinf=0;
#endif

	//Loop over the selected previously activated individuals
	for(i=nepainvitees-1; i>=0; --i) {
	  //Pick a random previously activated individual
	  j=gsl_rng_uniform(sv->r)*fpsim->nactivated;
          DEBUG_PRINTF("Previously activated invitee %i/%i at index %i\n",i,nepainvitees,j);
	  wasinfectious=(fpsim->activated[j]->indinfstatus==ro_ind_inf_status_infectious);

	  //If the individual is no longer activated
	  if(!ind_update_next_change_time(fpsim->activated[j], sv->event_time)) {
	    DEBUG_PRINTF("Individual at index %i is no longer activated\n",j);

	    //If the individual has not participated to any event
	    if(fpsim->activated[j]->ii.nevents==0) sv->new_inf_proc_func_noevent(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);

	    //Else if the individual participated to at least one event
	    else sv->end_inf_proc_func(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);

	    //Loop over the next previously activated individuals to see if they
	    //are no longer activated as well
	    for(k=j+1; k<fpsim->nactivated; ++k) {

	      if(ind_update_next_change_time(fpsim->activated[k], sv->event_time)) break;
	      DEBUG_PRINTF("Individual at index %i is no longer activated as well\n",k);
	    }
	    sarray_remove_many_at(fpsim->ctx_activated, k-j, j);

	    if(fpsim->nactivated==0) {
	      DEBUG_PRINTF("No previously activated individual left for this event\n");
	      break;
	    }

	    //Else if the selected individual is still activated and is infectious,
	    //add it to the list of infectious individuals for the event
	  } else if(fpsim->activated[j]->indinfstatus==ro_ind_inf_status_infectious) {
	    ++fpsim->activated[j]->ii.nevents;
	    fpsim->einfectious[fpsim->neinfectious++]=fpsim->activated[j];
#ifdef DUAL_PINF
	    epninf*=1-fpsim->activated[j]->ii.pinf;
	    etpinf+=fpsim->activated[j]->ii.pinf;
#endif

	    //If the infectious individuals was not previously infectious, call
	    //the appropriate function when there is at least one event
	    if(!wasinfectious) sv->new_inf_proc_func(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);
	  }
	}

	//If there is at least one infectious individual at the event
	if(fpsim->neinfectious>0) {

#ifdef DUAL_PINF
	  //If the susceptible individuals have not been initialised yet, it is
	  //time to generate the number of susceptible individuals in each category within
	  //the whole population
	  if(susstate==ro_fp_sus_state_uninitialized) {
	    nsusceptibles=sim->popsize-sim->nstart;
	    nsusceptiblesp=gsl_ran_binomial(sv->r, sim->ppip, nsusceptibles);
	    nsusceptiblesf=nsusceptibles-nsusceptiblesp;
	    DEBUG_PRINTF("Initialising the population to %u susceptible individuals (%u first category, %u second category)\n",nsusceptibles,nsusceptiblesf,nsusceptiblesp);
	  }

	  //Generate the number of susceptible individuals of each category for
	  //the current event
	  nesusceptibles=neinvitees-nepainvitees;
	  nesusceptiblesf=gsl_ran_hypergeometric(sv->r, nsusceptiblesf, nsusceptibles, nesusceptibles);
	  nesusceptiblesp=nesusceptibles-nesusceptiblesf;

#else
	  if(susstate==ro_fp_sus_state_uninitialized) {
	    nsusceptibles=sim->popsize-sim->nstart;
	    DEBUG_PRINTF("Initialising the population to %u susceptible individuals\n",nsusceptibles);
	  }

	  //Number of susceptible individuals for the current event
	  nesusceptibles=neinvitees-nepainvitees;
#endif

	  //If there is only one infectious individual at the event
	  if(fpsim->neinfectious==1) {

#ifdef DUAL_PINF
	    //Generate the number of infections occurring during the event for
	    //each category of individual, based on the probability of infection when
	    //the individuals are exposed to the single infectious individual
	    fpsim->einfectious[0]->ii.ninfectionsf=gsl_ran_binomial(sv->r, 1-epninf, nesusceptiblesf); 
	    fpsim->einfectious[0]->ii.ninfectionsp=gsl_ran_binomial(sv->r, 1-epninf*sim->rpinfp, nesusceptiblesp); 
	    fpsim->einfectious[0]->ii.ninfections=fpsim->einfectious[0]->ii.ninfectionsf+fpsim->einfectious[0]->ii.ninfectionsp;
            DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible and %u infections (%u first category, %u second category) were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,fpsim->einfectious[0]->ii.ninfections,fpsim->einfectious[0]->ii.ninfectionsf,fpsim->einfectious[0]->ii.ninfectionsp);

#else

	    //Generate the number of infections occurring during the event, based on the probability of
	    //infection when the individuals are exposed to the single infectious individual
	    fpsim->einfectious[0]->ii.ninfections=gsl_ran_binomial(sv->r, sim->pinf, nesusceptibles); 
            DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible and %u infections were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,fpsim->einfectious[0]->ii.ninfections);
#endif

	    //Create the event for the single infectious individual.
	    //If additional generations of infection should occur
	    if(sv->new_event_proc_func(sv, &fpsim->einfectious[0]->ii)) {

#ifdef DUAL_PINF
	      //Create new activated individuals out of susceptible individuals 
	      //First category
	      for(i=fpsim->einfectious[0]->ii.ninfectionsf-1; i>=0; --i) {
                DEBUG_PRINTF("Infection first category %i/%i\n",i,fpsim->einfectious[0]->ii.ninfectionsf);
		ind=fpsim->is+sim->nstart+nsusceptibles-i;
		ind->parent=fpsim->einfectious[0];
		ind->ii.generation=fpsim->einfectious[0]->ii.generation+1;
		ind->ii.q=sim->q;
		ind->ii.pinf=sim->pinf;

	        sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[0]->ii, sv->event_time);
                DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		ind_init_next_change_time(ind);
		fpsim->activated[fpsim->nactivated++]=ind;
	      }
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfectionsf;
	      nsusceptiblesf-=fpsim->einfectious[0]->ii.ninfectionsf;

	      //Second category
	      for(i=fpsim->einfectious[0]->ii.ninfectionsp-1; i>=0; --i) {
                DEBUG_PRINTF("Infection second category %i/%i\n",i,fpsim->einfectious[0]->ii.ninfectionsp);
		ind=fpsim->is+sim->nstart+nsusceptibles-i;
		ind->parent=fpsim->einfectious[0];
		ind->ii.generation=fpsim->einfectious[0]->ii.generation+1;
		ind->ii.inftypep=true;
		ind->ii.q=sim->qp;
		ind->ii.pinf=ppinf;

	        sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[0]->ii, sv->event_time);
                DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		ind_init_next_change_time(ind);
		fpsim->activated[fpsim->nactivated++]=ind;
	      }
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfectionsp;
	      nsusceptiblesp-=fpsim->einfectious[0]->ii.ninfectionsp;
#else
	      //Create new activated individuals out of susceptible individuals 
	      for(i=fpsim->einfectious[0]->ii.ninfections-1; i>=0; --i) {
                DEBUG_PRINTF("Infection %i/%i\n",i,fpsim->einfectious[0]->ii.ninfections);
		ind=fpsim->is+fpsim->nstart+nsusceptibles-i;
		ind->parent=fpsim->einfectious[0];
		ind->ii.generation==fpsim->einfectious[0]->ii.generation+1;
	        sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[0]->ii, sv->event_time);
                DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		ind_init_next_change_time(ind);
		fpsim->activated[fpsim->nactivated++]=ind;
	      }
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfections;


#endif

	    //Else we are done with this path
	    } else {
	      DEBUG_PRINTF("New event returned false\n");
	      break;
	    }

	    //Else if there is more than one infectious individual at the event
	  } else {

#ifdef DUAL_PINF
	    //Generate the number of infections occurring during the event for
	    //each category of individual, based on the combined probability of non-infection  when
	    //the individuals are exposed to the combination of all infectious individuals
	    neinfectionsf=gsl_ran_binomial(sv->r, 1-epninf, nesusceptiblesf); 
	    neinfectionsp=gsl_ran_binomial(sv->r, 1-epninf*pow(sim->rpinfp,fpsim->neinfectious), nesusceptiblesp); 
            DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible and %u infections (%u first category, %u second category) were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,neinfectionsf+neinfectionsp,neinfectionsf,neinfectionsp);

#else

	    //Generate the number of infections occurring during the event, based on the combined probability of
	    //non-infection when the individuals are exposed to the combination of all infectious individuals
	    neinfections=gsl_ran_binomial(sv->r, 1-pow(1-sim->pinf,fpsim->neinfectious), nesusceptibles); 
            DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible and %u infections were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,neinfections);
#endif

	    positiveevents=true;

	    //Loop over each infectious individual at the event
	    for(j=fpsim->neinfectious-1; j>=0; --j) {
              DEBUG_PRINTF("Infectious individual %i/%lu\n",j,fpsim->neinfectious);

	      //Create the event for the current infectious individual
	      //If additional generations of infection should occur
	      if(sv->new_event_proc_func(sv, &fpsim->einfectious[j]->ii)) {
#ifdef DUAL_PINF
		//Generate the number of new infections associated to the current infectious
		//individual, based on the overall weight of his probability of
		//infecting others amongst all remaining infectious individuals.
		//This should be equivalent to assign a number of infections to
		//each infectious individual using a multinomial distribution.
		fpsim->einfectious[j]->ii.ninfectionsf=gsl_ran_binomial(sv->r, fpsim->activated[j]->ii.pinf/etpinf, neinfectionsf); 
		fpsim->einfectious[j]->ii.ninfectionsp=gsl_ran_binomial(sv->r, fpsim->activated[j]->ii.pinf/etpinf, neinfectionsp); 
		fpsim->einfectious[j]->ii.ninfections=fpsim->einfectious[j]->ii.ninfectionsf+fpsim->einfectious[j]->ii.ninfectionsp;
		neinfectionsf-=fpsim->einfectious[j]->ii.ninfectionsf;
		neinfectionsp-=fpsim->einfectious[j]->ii.ninfections;
		etpinf-=fpsim->activated[j]->ii.pinf;

		//Create new activated individuals out of susceptible individuals 
		//First category
		for(i=fpsim->einfectious[j]->ii.ninfectionsf-1; i>=0; --i) {
                  DEBUG_PRINTF("Infection first category %i/%i\n",i,fpsim->einfectious[j]->ii.ninfectionsf);
		  ind=fpsim->is+sim->nstart+nsusceptibles-i;
		  ind->parent=fpsim->einfectious[j];
		  ind->ii.generation=fpsim->einfectious[j]->ii.generation+1;
		  ind->ii.q=sim->q;
		  ind->ii.pinf=sim->pinf;

		  sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time);
                  DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		  ind_init_next_change_time(ind);
		  fpsim->activated[fpsim->nactivated++]=ind;
		}
		nsusceptibles-=fpsim->einfectious[j]->ii.ninfectionsf;
		nsusceptiblesf-=fpsim->einfectious[j]->ii.ninfectionsf;

		//Second category
		for(i=fpsim->einfectious[j]->ii.ninfectionsp-1; i>=0; --i) {
                  DEBUG_PRINTF("Infection second category %i/%i\n",i,fpsim->einfectious[j]->ii.ninfectionsp);
		  ind=fpsim->is+sim->nstart+nsusceptibles-i;
		  ind->parent=fpsim->einfectious[j];
		  ind->ii.generation=fpsim->einfectious[j]->ii.generation+1;
		  ind->ii.inftypep=true;
		  ind->ii.q=sim->qp;
		  ind->ii.pinf=ppinf;

		  sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time);
                  DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		  ind_init_next_change_time(ind);
		  fpsim->activated[fpsim->nactivated++]=ind;
		}
		nsusceptibles-=fpsim->einfectious[j]->ii.ninfectionsp;
		nsusceptiblesp-=fpsim->einfectious[j]->ii.ninfectionsp;
#else
		//Generate the number of new infections associated to the
		//current infectious individual
		fpsim->einfectious[j]->ii.ninfections=gsl_ran_binomial(sv->r, 1./fpsim->neinfectious, neinfections);

		//Create new activated individuals out of susceptible individuals 
		for(i=fpsim->einfectious[j]->ii.ninfections-1; i>=0; --i) {
                  DEBUG_PRINTF("Infection %i/%i\n",i,fpsim->einfectious[j]->ii.ninfections);
		  ind=fpsim->is+fpsim->nstart+nsusceptibles-i;
		  ind->parent=fpsim->einfectious[j];
		  ind->ii.generation==fpsim->einfectious[j]->ii.generation+1;
		  sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time);
                  DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		  ind_init_next_change_time(ind);
		  fpsim->activated[fpsim->nactivated++]=ind;
		}
	        nsusceptibles-=fpsim->einfectious[j]->ii.ninfections;
#endif
	      } else {
		DEBUG_PRINTF("New event returned false\n");
		positiveevents=false;
	      }
	    }

	    //If event creation returns false for at least one infectious
	    //individual, then we are done with this path
	    if(!positiveevents) break;
	  }
	}

      } else {
	DEBUG_PRINTF("No previously activated invtee\n");
      }
    }
    
    //Loop over the remaining activated individuals
    for(i=fpsim->nactivated-1; i>=0; --i) {

      //If the individual has not participated to any event
      if(fpsim->activated[j]->ii.nevents==0) sv->new_inf_proc_func_noevent(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);

      //Else if the individual participated to at least one event
      else sv->end_inf_proc_func(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);
    }

  } while(!sv->path_end_proc_func(sv));

  return 0;
}

void finitepopsim_free(sim_vars* sv)
{
  if(sv->fpsim.rooti.ii.dataptr) free(sv->fpsim.rooti.ii.dataptr);

  for(int32_t i=sv->pars.popsize-1; i>=0; --i) if(sv->fpsim.is[i].ii.dataptr) free(sv->fpsim.is[i].ii.dataptr);
  free(sv->fpsim.is);

  sarray_free(sv->fpsim.ctx_activated);
  sarray_free(sv->fpsim.ctx_einfectious);
}
