/**
 * @file finitepopsim.c
 * @brief Finite population simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 *
 * Note: Code does not currently supports contacts. Also it only supports
 * group_invitees
 */


#include <assert.h>
#include "finitepopsim.h"

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers
#define II_ARRAY_GROW_FACT (1.5)  //!< Growing factor for the array of current infectious individuals across all layers.

void finitepopsim_init(sim_vars* sv)
{
  memset(&sv->fpsim.rooti, 0, sizeof(individual));
  sv->fpsim.rooti.ii.commpertype=0;
  sv->fpsim.rooti.ii.nevents=1;
  sv->fpsim.rooti.ii.nattendees=1;
  sv->fpsim.rooti.ii.ninfections=1;
  sv->fpsim.rooti.ii.generation=0;

  sv->fpsim.is=(individual*)malloc(sv->pars.popsize*sizeof(individual));

  sv->ii_alloc_proc_func(&sv->fpsim.rooti.ii);
  for(int32_t i=sv->pars.popsize-1; i>=0; --i) sv->ii_alloc_proc_func(&sv->fpsim.is[i].ii);
  ran_log_init(&sv->rl, (rng_stream*)sv->r->state, sv->pars.p);

  sv->fpsim.nactivated=0;
  sv->fpsim.activated=malloc(sv->pars.popsize*sizeof(individual*));
  sv->fpsim.neinfectious=0;
  sv->fpsim.einfectious=malloc(sv->pars.popsize*sizeof(individual*));

  FP_GENINF_COND();
}

int finitepopsim(sim_vars* sv)
{
  int i,j;
  model_pars const* sim=&(sv->pars);
  fpsim_vars* fpsim=&sv->fpsim;
  void* ptr;
  individual* ind;
  uint32_t neinvitees;   //Number of invitees at the event
  uint32_t nepainvitees; //Number of potentially activated invitees at the event
  uint32_t nsusceptibles=UINT32_MAX; //Total number of remaining susceptible individuals within the whole population
  uint32_t nesusceptibles; //Total number of susceptible individuals at the event
#ifdef DUAL_PINF
  uint32_t nsusceptiblesf=UINT32_MAX; //Number of remaining susceptible individuals of the first category within the whole population
  uint32_t nsusceptiblesp=UINT32_MAX; //Number of remaining susceptible individuals of the second category within the whole population
  uint32_t nesusceptiblesf; //Number of susceptible individuals of the first category at the event
  uint32_t nesusceptiblesp; //Number of susceptible individuals of the second category at the event
  uint32_t neinfectionsf;   //Number of infected individuals of the first category at the event
  uint32_t neinfectionsp;   //Number of infected individuals of the second category at the event
  double epninff;           //Probability of non-infection at the event for the first category of infected individual
  double epninfp;           //Probability of non-infection at the event for the second category of infected individual
  double etpinf;            //Sum of the probabilities of infections at the event

  const double ppinf=sim->pinf*sim->rpshedp;
  const double pinfpinf=sim->ppip*sim->rpinfp/(1+sim->ppip*(sim->rpinfp-1));
  bool initsus;
#else
  uint32_t neinfections;  //Total number of infections at the event
#endif
  bool wasinfectious;

  //Nstart individuals are assumed to be infected when introduced to the
  //simulation

  do {
    fpsim->nactivated=0;
    nsusceptibles=sim->popsize-sim->nstart;
#ifdef DUAL_PINF
    initsus=false;
#endif

    for(i=sim->popsize-1; i>=0; --i) {
      DEBUG_PRINTF("Individual %i is %p\n",i,fpsim->is+i);
      ptr=fpsim->is[i].ii.dataptr;
      memset(fpsim->is+i, 0, sizeof(individual));
      fpsim->is[i].ii.dataptr=ptr;
    }

    sv->path_init_proc_func(sv);
    sv->event_time=0;

    for(i=sim->nstart-1; i>=0; --i) {
      DEBUG_PRINTF("initial individual %i is %p\n",i,fpsim->is+i);

#ifdef DUAL_PINF
      if(gsl_rng_uniform(sv->r) < pinfpinf) {
#ifdef SEC_INF_TIMELINES
	fpsim->rooti.ii.ninfectionsf=0;
	fpsim->rooti.ii.ninfectionsp=1;
#endif
	fpsim->is[i].ii.inftypep=true;
	fpsim->is[i].ii.q=sim->qp;
	fpsim->is[i].ii.pinf=ppinf;

      } else {
#ifdef SEC_INF_TIMELINES
	fpsim->rooti.ii.ninfectionsf=1;
	fpsim->rooti.ii.ninfectionsp=0;
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
      //Activated array is not up to date

      //Generate the number of invitees for the current event
      neinvitees=sv->gen_att_func(sv);
      DEBUG_PRINTF("Event at time %f with %u invitees\n",sv->event_time,neinvitees);

      //Generate a number of previously known activated individuals amongst
      //these invitees. The drawn is without replacement here, so it is a
      //hypergeometric distribution and not a binomial
      nepainvitees=gsl_ran_hypergeometric(sv->r, fpsim->nactivated, sim->popsize-fpsim->nactivated, neinvitees);
      //DEBUG_PRINTF("# event previously activated is %li/%u * %u = %u\n",fpsim->nactivated,sim->popsize,neinvitees,nepainvitees);

      //If there is no previously known activated individual selected for this
      //event, go to the next event
      //if(nepainvitees==0 || nepainvitees==neinvitees) continue;

      //Generate the number of susceptible individuals for the current event
      //Generate a number of susceptible individuals amongst
      //the invitees who are not previously activated invididuals.
      //The drawn is without replacement here, so it is a
      //hypergeometric distribution and not a binomial
      //This works because previous activated individuals are not returned
      //as susceptible individuals
      nesusceptibles=gsl_ran_hypergeometric(sv->r, nsusceptibles, sim->popsize-fpsim->nactivated-nsusceptibles, neinvitees-nepainvitees);
      //DEBUG_PRINTF("# event susceptibles is %u/%li * %u = %u\n",nsusceptibles,sim->popsize-fpsim->nactivated,neinvitees-nepainvitees,nesusceptibles);

      //If there is no susceptible individual selected for this
      //event, go to the next event
      if(nesusceptibles==0) continue;

      //DEBUG_PRINTF("Event at time %f\n",sv->event_time);
      //Identify the infectious individuals amongst the previously known
      //activated individuals, updating the activated array at the same time
      fpsim->neinfectious=0;
#ifdef DUAL_PINF
      epninff=epninfp=1;
      etpinf=0;
#endif

      //Loop over the selected previously activated individuals
      for(i=nepainvitees-1; i>=0; --i) {
	//Pick a random previously activated individual amongst the previously
	//activated individuals who have not been picked yet
	assert(fpsim->nactivated>fpsim->neinfectious);
	double dbuf=fpsim->neinfectious+gsl_rng_uniform(sv->r)*(fpsim->nactivated-fpsim->neinfectious);
	//DEBUG_PRINTF("Uniform %li -> %li = %f\n",fpsim->neinfectious,fpsim->nactivated,dbuf);
	j=dbuf;
	DEBUG_PRINTF("Previously activated invitee %i/%i %p\n",i,nepainvitees,fpsim->activated[j]);
	wasinfectious=(fpsim->activated[j]->indinfstatus==ro_ind_inf_status_infectious);

	//If the individual is no longer activated
	if(!ind_update_next_change_time(fpsim->activated[j], sv->event_time)) {

	  //If the individual has not participated to any event
	  if(fpsim->activated[j]->ii.nevents==0) sv->new_inf_proc_func_noevent(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);

	  //Else if the individual participated to at least one event
	  else sv->end_inf_proc_func(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);
	  DEBUG_PRINTF("Individual %p is no longer activated\n",fpsim->activated[j]);

	  //Replace the individual with the last previously activated one from
	  //the end
	  fpsim->activated[j]=fpsim->activated[fpsim->nactivated-1];
	  --(fpsim->nactivated);
	  DEBUG_PRINTF("Number of activated individuals down to %li\n",fpsim->nactivated);

	  if(fpsim->nactivated==fpsim->neinfectious) {

	    if(fpsim->nactivated==0) {
	      DEBUG_PRINTF("No previously activated individual left in the current path\n");
	      goto done_parsing;
	    }
	    DEBUG_PRINTF("No previously activated individual left for this event\n");
	    break;
	  }

	  //Else if the selected individual is still activated and is infectious,
	  //add it to the list of infectious individuals for the event
	} else if(fpsim->activated[j]->indinfstatus==ro_ind_inf_status_infectious) {
	  DEBUG_PRINTF("Individual %p is currently infectious\n",fpsim->activated[j]);
	  ++(fpsim->activated[j]->ii.nevents);
	  fpsim->einfectious[fpsim->neinfectious++]=fpsim->activated[j];
#ifdef DUAL_PINF
	  epninff*=1-fpsim->activated[j]->ii.pinf;
	  epninfp*=1-fpsim->activated[j]->ii.pinf*sim->rpinfp;
	  etpinf+=fpsim->activated[j]->ii.pinf;
#endif

	  //If the infectious individuals was not previously infectious, call
	  //the appropriate function when there is at least one event
	  if(!wasinfectious) sv->new_inf_proc_func(sv, &fpsim->activated[j]->ii, &fpsim->activated[j]->parent->ii);

	  //Move the infectious individual at the beginning of the array.
	  ind=fpsim->activated[fpsim->neinfectious-1];
	  fpsim->activated[fpsim->neinfectious-1]=fpsim->activated[j];
	  fpsim->activated[j]=ind;

	  if(fpsim->nactivated==fpsim->neinfectious) {
	    DEBUG_PRINTF("No previously activated individual left for this event\n");
	    break;
	  }

	} else {
	  DEBUG_PRINTF("Individual %p is not infectious yet\n",fpsim->activated[j]);
	}
      }

      //If there is at least one infectious individual at the event
      if(fpsim->neinfectious>0) {

#ifdef DUAL_PINF
	//If the susceptible individuals have not been initialised yet, it is
	//time to generate the number of susceptible individuals in each category within
	//the whole population
	if(!initsus) {
	  nsusceptiblesp=gsl_ran_binomial(sv->r, sim->ppip, nsusceptibles);
	  nsusceptiblesf=nsusceptibles-nsusceptiblesp;
	  DEBUG_PRINTF("Initialising the population to %u susceptible individuals (%u first category, %u second category)\n",nsusceptibles,nsusceptiblesf,nsusceptiblesp);
	  initsus=true;
	}
	nesusceptiblesf=gsl_ran_hypergeometric(sv->r, nsusceptiblesf, nsusceptiblesp, nesusceptibles);
	nesusceptiblesp=nesusceptibles-nesusceptiblesf;
#endif

	//If there is only one infectious individual at the event
	if(fpsim->neinfectious==1) {

#ifdef DUAL_PINF
	  //Generate the number of infections occurring during the event for
	  //each category of individual, based on the probability of infection when
	  //the individuals are exposed to the single infectious individual
	  assert(nesusceptibles==nesusceptiblesf+nesusceptiblesp);
	  fpsim->einfectious[0]->ii.ninfectionsf=gsl_ran_binomial(sv->r, 1-epninff, nesusceptiblesf); 
	  fpsim->einfectious[0]->ii.ninfectionsp=gsl_ran_binomial(sv->r, 1-epninfp, nesusceptiblesp); 
	  fpsim->einfectious[0]->ii.ninfections=fpsim->einfectious[0]->ii.ninfectionsf+fpsim->einfectious[0]->ii.ninfectionsp;
	  DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible individuals and %u infections (%u first category, %u second category) were generated for this event with infection probabilities %f%% and %f%%\n",neinvitees,fpsim->neinfectious,nesusceptibles,fpsim->einfectious[0]->ii.ninfections,fpsim->einfectious[0]->ii.ninfectionsf,fpsim->einfectious[0]->ii.ninfectionsp,100*(1-epninff),100*(1-epninfp));

#else

	  //Generate the number of infections occurring during the event, based on the probability of
	  //infection when the individuals are exposed to the single infectious individual
	  fpsim->einfectious[0]->ii.ninfections=gsl_ran_binomial(sv->r, sim->pinf, nesusceptibles); 
	  DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible individuals and %u infections were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,fpsim->einfectious[0]->ii.ninfections);
#endif

	  //Create the event for the single infectious individual.
	  //If additional generations of infection should occur
	  if(fpsim->einfectious[0]->ii.ninfections>0) {

	    if(sv->new_event_proc_func(sv, &fpsim->einfectious[0]->ii)) {

#ifdef DUAL_PINF
	      //Create new activated individuals out of susceptible individuals 
	      //First category
	      for(i=fpsim->einfectious[0]->ii.ninfectionsf-1; i>=0; --i) {
		ind=fpsim->is+sim->nstart+nsusceptibles-1-i;
		DEBUG_PRINTF("Infection first category %i/%i for individual %p\n",i,fpsim->einfectious[0]->ii.ninfectionsf,ind);
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
		ind=fpsim->is+sim->nstart+nsusceptibles-1-i;
		DEBUG_PRINTF("Infection second category %i/%i for individual %p\n",i,fpsim->einfectious[0]->ii.ninfectionsp,ind);
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
	      DEBUG_PRINTF("Total number of susceptible individuals down to %u (%u first category, %u second category)\n",nsusceptibles,nsusceptiblesf,nsusceptiblesp);
#else
	      //Create new activated individuals out of susceptible individuals 
	      for(i=fpsim->einfectious[0]->ii.ninfections-1; i>=0; --i) {
		ind=fpsim->is+sim->nstart+nsusceptibles-1-i;
		DEBUG_PRINTF("Infection %i/%i for individual %p\n",i,fpsim->einfectious[0]->ii.ninfections,ind);
		ind->parent=fpsim->einfectious[0];
		ind->ii.generation=fpsim->einfectious[0]->ii.generation+1;
		sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[0]->ii, sv->event_time);
		DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period);

		ind_init_next_change_time(ind);
		fpsim->activated[fpsim->nactivated++]=ind;
	      }
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfections;
	      DEBUG_PRINTF("Total number of susceptible individuals down to %u\n",nsusceptibles);
#endif

	      //Else we are done with this event
	    } else {
	      DEBUG_PRINTF("New event returned false\n");
#ifdef DUAL_PINF
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfectionsf+fpsim->einfectious[0]->ii.ninfectionsp;
	      nsusceptiblesf-=fpsim->einfectious[0]->ii.ninfectionsf;
	      nsusceptiblesp-=fpsim->einfectious[0]->ii.ninfectionsp;
	      DEBUG_PRINTF("Total number of susceptible individuals down to %u (%u first category, %u second category)\n",nsusceptibles,nsusceptiblesf,nsusceptiblesp);
#else
	      nsusceptibles-=fpsim->einfectious[0]->ii.ninfections;
	      DEBUG_PRINTF("Total number of susceptible individuals down to %u\n",nsusceptibles);
#endif
	    }
	  }

	  //Else if there is more than one infectious individual at the event
	} else {

#ifdef DUAL_PINF
	  //Generate the number of infections occurring during the event for
	  //each category of individual, based on the combined probability of non-infection when
	  //the individuals are exposed to the combination of all infectious individuals
	  assert(nesusceptibles==nesusceptiblesf+nesusceptiblesp);
	  neinfectionsf=gsl_ran_binomial(sv->r, 1-epninff, nesusceptiblesf); 
	  neinfectionsp=gsl_ran_binomial(sv->r, 1-epninfp, nesusceptiblesp); 
	  DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible individuals and %u infections (%u first category, %u second category) were generated for this event with infection probabilities %f%% and %f%%\n",neinvitees,fpsim->neinfectious,nesusceptibles,neinfectionsf+neinfectionsp,neinfectionsf,neinfectionsp,100*(1-epninff),100*(1-epninfp));
	  if(neinfectionsf>0 || neinfectionsp>0) {

#else

	    //Generate the number of infections occurring during the event, based on the combined probability of
	    //non-infection when the individuals are exposed to the combination of all infectious individuals
	    neinfections=gsl_ran_binomial(sv->r, 1-pow(1-sim->pinf,fpsim->neinfectious), nesusceptibles); 
	    DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible individuals and %u infections were generated for this event\n",neinvitees,fpsim->neinfectious,nesusceptibles,neinfections);
	    if(neinfections>0) {
#endif

	      //Loop over each infectious individual at the event
	      for(j=fpsim->neinfectious-1; j>0; --j) {
		DEBUG_PRINTF("Infectious individual %i/%lu\n",j,fpsim->neinfectious);
#ifdef DUAL_PINF
		//Generate the number of new infections associated to the current infectious
		//individual, based on the overall weight of his probability of
		//infecting others amongst all remaining infectious individuals.
		//This should be equivalent to assign a number of infections to
		//each infectious individual using a multinomial distribution.
		DEBUG_PRINTF("Infection assignment probability to this infectious individual is %f/%f (%f%%)\n",fpsim->einfectious[j]->ii.pinf,etpinf,100*fpsim->einfectious[j]->ii.pinf/etpinf);
		fpsim->einfectious[j]->ii.ninfectionsf=gsl_ran_binomial(sv->r, fpsim->einfectious[j]->ii.pinf/etpinf, neinfectionsf); 
		fpsim->einfectious[j]->ii.ninfectionsp=gsl_ran_binomial(sv->r, fpsim->einfectious[j]->ii.pinf/etpinf, neinfectionsp); 
		fpsim->einfectious[j]->ii.ninfections=fpsim->einfectious[j]->ii.ninfectionsf+fpsim->einfectious[j]->ii.ninfectionsp;
		DEBUG_PRINTF("Number of infections assigned to this infectious individual are %u/%u first category and %u/%u second category\n",fpsim->einfectious[j]->ii.ninfectionsf,neinfectionsf,fpsim->einfectious[j]->ii.ninfectionsp,neinfectionsp);
		neinfectionsf-=fpsim->einfectious[j]->ii.ninfectionsf;
		neinfectionsp-=fpsim->einfectious[j]->ii.ninfectionsp;
		etpinf-=fpsim->einfectious[j]->ii.pinf;
#else
		//Generate the number of new infections associated to the
		//current infectious individual
		fpsim->einfectious[j]->ii.ninfections=gsl_ran_binomial(sv->r, 1./(j+1), neinfections);
		neinfections-=fpsim->einfectious[j]->ii.ninfections;
#endif

#ifdef DUAL_PINF
#define FP_CREATE_EVENT_MULTI_INFECTIOUS(j) \
		/*Create the event for the current infectious individual*/ \
		/*If additional generations of infection should occur*/ \
		if(sv->new_event_proc_func(sv, &fpsim->einfectious[j]->ii)) { \
		  \
		  /*Create new activated individuals out of susceptible individuals*/ \
		  /*First category*/ \
		  for(i=fpsim->einfectious[j]->ii.ninfectionsf-1; i>=0; --i) { \
		    ind=fpsim->is+sim->nstart+nsusceptibles-1-i; \
		    DEBUG_PRINTF("Infection first category %i/%i for individual %p\n",i,fpsim->einfectious[j]->ii.ninfectionsf,ind); \
		    ind->parent=fpsim->einfectious[j]; \
		    ind->ii.generation=fpsim->einfectious[j]->ii.generation+1; \
		    ind->ii.q=sim->q; \
		    ind->ii.pinf=sim->pinf; \
		    \
		    sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time); \
		    DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period); \
		    \
		    ind_init_next_change_time(ind); \
		    fpsim->activated[fpsim->nactivated++]=ind; \
		  } \
		  nsusceptibles-=fpsim->einfectious[j]->ii.ninfectionsf; \
		  nsusceptiblesf-=fpsim->einfectious[j]->ii.ninfectionsf; \
		  \
		  /*Second category*/ \
		  for(i=fpsim->einfectious[j]->ii.ninfectionsp-1; i>=0; --i) { \
		    ind=fpsim->is+sim->nstart+nsusceptibles-1-i; \
		    DEBUG_PRINTF("Infection second category %i/%i for individual %p\n",i,fpsim->einfectious[j]->ii.ninfectionsp,ind); \
		    ind->parent=fpsim->einfectious[j]; \
		    ind->ii.generation=fpsim->einfectious[j]->ii.generation+1; \
		    ind->ii.inftypep=true; \
		    ind->ii.q=sim->qp; \
		    ind->ii.pinf=ppinf; \
		    \
		    sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time); \
		    DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period); \
		    \
		    ind_init_next_change_time(ind); \
		    fpsim->activated[fpsim->nactivated++]=ind; \
		  } \
		  nsusceptibles-=fpsim->einfectious[j]->ii.ninfectionsp; \
		  nsusceptiblesp-=fpsim->einfectious[j]->ii.ninfectionsp; \
		  \
		} else { \
		  DEBUG_PRINTF("New event returned false\n"); \
		  nsusceptibles-=fpsim->einfectious[j]->ii.ninfectionsf+fpsim->einfectious[j]->ii.ninfectionsp; \
		  nsusceptiblesf-=fpsim->einfectious[j]->ii.ninfectionsf; \
		  nsusceptiblesp-=fpsim->einfectious[j]->ii.ninfectionsp; \
		} \
		DEBUG_PRINTF("Total number of susceptible individuals down to %u (%u first category, %u second category)\n",nsusceptibles,nsusceptiblesf,nsusceptiblesp);
#else
#define FP_CREATE_EVENT_MULTI_INFECTIOUS(j) \
		/*Create the event for the current infectious individual*/ \
		/*If additional generations of infection should occur*/ \
		if(sv->new_event_proc_func(sv, &fpsim->einfectious[j]->ii)) { \
		  \
		  /*Create new activated individuals out of susceptible individuals*/ \
		  for(i=fpsim->einfectious[j]->ii.ninfections-1; i>=0; --i) { \
		    ind=fpsim->is+sim->nstart+nsusceptibles-1-i; \
		    DEBUG_PRINTF("Infection %i/%i for individual %p\n",i,fpsim->einfectious[j]->ii.ninfections,ind); \
		    ind->parent=fpsim->einfectious[j]; \
		    ind->ii.generation=fpsim->einfectious[j]->ii.generation+1; \
		    sv->gen_time_periods_func(sv, &ind->ii, &fpsim->einfectious[j]->ii, sv->event_time); \
		    DEBUG_PRINTF("Latent period is %f, comm period is %f, type is %u, end comm is %f\n",ind->ii.latent_period,ind->ii.comm_period,ind->ii.commpertype,ind->ii.end_comm_period); \
		    \
		    ind_init_next_change_time(ind); \
		    fpsim->activated[fpsim->nactivated++]=ind; \
		  } \
		  nsusceptibles-=fpsim->einfectious[j]->ii.ninfections; \
		  \
		} else { \
		  nsusceptibles-=fpsim->einfectious[j]->ii.ninfections; \
		} \
		DEBUG_PRINTF("Total number of susceptible individuals down to %u\n",nsusceptibles);
#endif
		if(fpsim->einfectious[j]->ii.ninfections>0) {
		  FP_CREATE_EVENT_MULTI_INFECTIOUS(j);
		}
#ifdef DUAL_PINF
		if(neinfectionsp==0 && neinfectionsf==0) goto done_multi_infections;
#else
		if(neinfections==0) goto done_multi_infections;
#endif
	      }

	      DEBUG_PRINTF("Infectious individual %i/%lu\n",0,fpsim->neinfectious);
#ifdef DUAL_PINF
	      assert(neinfectionsf+neinfectionsp>0);
	      //Generate the number of new infections associated to the current infectious
	      //individual, based on the overall weight of his probability of
	      //infecting others amongst all remaining infectious individuals.
	      //This should be equivalent to assign a number of infections to
	      //each infectious individual using a multinomial distribution.
	      DEBUG_PRINTF("Infection assignment probability to this infectious individual is %f/%f (%f%%)\n",fpsim->einfectious[0]->ii.pinf,etpinf,100*fpsim->einfectious[0]->ii.pinf/etpinf);
	      fpsim->einfectious[0]->ii.ninfectionsf=neinfectionsf; 
	      fpsim->einfectious[0]->ii.ninfectionsp=neinfectionsp; 
	      fpsim->einfectious[0]->ii.ninfections=neinfectionsf+neinfectionsp;
	      DEBUG_PRINTF("Number of infections assigned to this infectious individual are %u/%u first category and %u/%u second category\n",fpsim->einfectious[0]->ii.ninfectionsf,neinfectionsf,fpsim->einfectious[0]->ii.ninfectionsp,neinfectionsp);
#else
	      //Generate the number of new infections associated to the
	      //current infectious individual
	      fpsim->einfectious[0]->ii.ninfections=neinfections;
#endif
	      FP_CREATE_EVENT_MULTI_INFECTIOUS(0);
done_multi_infections:
	      ;
	    }
	  }

	  //If there is no susceptible individual left, we are done for this path
	  if(nsusceptibles==0) break;

        } else {
	  DEBUG_PRINTF("%u invitees, %lu infectious, %u susceptible individuals and %u infections (%u first category, %u second category) were generated for this event with infection probabilities %f%% and %f%%\n",neinvitees,fpsim->neinfectious,0,0,0,0,0.,0.);
        }
      }
done_parsing:

    //Loop over the remaining activated individuals
    for(i=fpsim->nactivated-1; i>=0; --i) {

      //If the individual has not participated to any event
      if(fpsim->activated[i]->ii.nevents==0) sv->new_inf_proc_func_noevent(sv, &fpsim->activated[i]->ii, &fpsim->activated[i]->parent->ii);

      //Else if the individual participated to at least one event
      else sv->end_inf_proc_func(sv, &fpsim->activated[i]->ii, &fpsim->activated[i]->parent->ii);
    }

  } while(!sv->path_end_proc_func(sv));

  return 0;
}

void finitepopsim_free(sim_vars* sv)
{
  if(sv->fpsim.rooti.ii.dataptr) free(sv->fpsim.rooti.ii.dataptr);

  for(int32_t i=sv->pars.popsize-1; i>=0; --i) if(sv->fpsim.is[i].ii.dataptr) free(sv->fpsim.is[i].ii.dataptr);
  free(sv->fpsim.is);

  free(sv->fpsim.activated);
  free(sv->fpsim.einfectious);
}
