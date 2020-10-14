/**
 * @file standard_summary_stats.h
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _STANDARD_SUMMARY_STATS_
#define _STANDARD_SUMMARY_STATS_

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include <assert.h>

#include "infindividual.h"
#include "simulation.h"

#define INIT_NINF_ALLOC (16) //!< Initial number of allocated bins for the ninf histogram
#define INIT_NACTENTRIES (16)
#define CTENTRIES_GROWFACT (1.5)

extern int __ro_debug;
#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) {if(__ro_debug) printf(__VA_ARGS__);} //!< Debug print function

#ifdef CT_OUTPUT
typedef struct {
  int32_t postesttime;
  int32_t presymtime;
  uint32_t id;
  uint32_t pid;
  uint32_t ntracedcts;
} ctposinf;
#endif

typedef struct
{
  uint32_t   n;	               //!< Number of infectious creating infections in the timeline bin.
  uint32_t rsum;	       //!< Sum of the number of individuals that get infected in the timeline bin.
  double commpersum;   	       //!< Sum of communicable periods for all infectious individuals whose communicable period start before tmax in the timeline bin.
#ifdef NUMEVENTSSTATS
  uint32_t neventssum;	       //!< Sum of the number of transmission events for all infectious individuals whose communicable period occurs before tmax in the timeline bin.
#endif
  uint64_t* ngeninfs;	       //!< Number of generated infections from each infectious individual in the timeline bin.
} ext_timeline_info;

/**
 * Simulation-level standard summary statistics data struct.
 */
typedef struct
{
  double extinction_time;	//!< *Path extinction time, if any. 
  uint32_t* inf_timeline;	//!< For each integer interval between 0 and floor(nbinsperunit*tmax)-1, the number of individuals that are infected, but not isolated, at some point in this interval.
  uint32_t* newinf_timeline;	//!< For each integer interval between 0 and floor(nbinsperunit*tmax)-1, the number of individuals that get infected at some point in this interval. For the last interval,
  uint32_t* postest_timeline;	//!< For each integer interval between 0 and floor(nbinsperunit*tmax)-1, the number of individuals that have a recent positive test result at the end of this interval.
  uint32_t* newpostest_timeline;//!< For each integer interval between 0 and floor(nbinsperunit*tmax), the number of individuals that receive a positive test result at some point in this interval.
  ext_timeline_info* ext_timeline;     //!< Extended timeline for parameters that cannot correctly be calculated when a dynamic time cur is used.
  uint32_t nainfbins;		//!< Number of allocated infectious individual bins.
  uint32_t ninfbins;		//!< Number of used infectious individual bins.
#ifdef CT_OUTPUT
  ctposinf** ctentries;
  uint32_t nctentries;
  uint32_t nactentries;
  int32_t curctid;
#endif
  int32_t npers;		//!< Number of positive integer intervals
  int32_t nbinsperunit;		//!< Number of timeline bins per unit of time.
  int32_t  tlshifta;            //!< *Allocated integral shift of the timeline origin (to allow for negative time bins). Corresponds also to the number of negative integer intervals
  int32_t  tlshift;             //!< *Integral shift of the timeline origin (to allow for negative time bins). Corresponds also to the number of negative integer intervals
  uint32_t tnpersa;             //!< *Total number of allocated integer intervals (negative+positive)
  int32_t  tlppt0idx;           //!< Timeline post-processing t=0 index
  int32_t tlppnnpers;		//!< Timeline post-processing number of negative periods
  uint32_t tlpptnvpers;		//!< Timeline post-processing total number of valid periods
  uint32_t lmax;                //!< Maximum number of layers for the simulation. lmax=1 means only primary infectious individuals.
  uint32_t nimax;               //!< Maximum number of infectious individuals for a given integer interval between 0 and floor(nbinsperunit*tmax)-1. Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  uint32_t npostestmax;         //!< Maximum number of positive test results during an interval of duration npostestmaxnpers for each individual that starts when the test results are received. Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  uint32_t npostestmaxnpers;    //!< Interval duration for the maximum number of positive test results
  int32_t maxedoutmintimeindex; //!< *Minimum time index which maxed out the allowed number of infected individuals or positive test results.
  //uint32_t n_ended_infections;
  bool extinction;		//!< *Set to true if extinction does not occur before tmax.
  bool timerelfirstpostestresults; //!< Indicate that the output time values are relative to the time of the first positive test results. This require the times to be recomputed in post processing.
} std_summary_stats;

/**
 * @brief Initialises the standard summary statistics.
 *
 * This function must be called once to initialise the standard summary
 * statistics, if used.
 *
 * @param sv: Pointer to the simulation variables.
 */
void std_stats_init(sim_vars *sv, const uint32_t nbinsperunit, bool ngeninfs, bool timerelfirstpostestresults);

/**
 * @brief Initialises elements of the standard summary statistics
 * before a path simulation.
 *
 * This function must be called before the simulation of each path to initialise
 * elements of the standard summary statistics, if used
 *
 * @param stats: Pointer to the standard summary statistics whose
 * elements will be initialised.
 */
inline static void std_stats_path_init(std_summary_stats* stats)
{
  stats->extinction_time=-INFINITY;
  const uint32_t nerase=stats->tlshift+stats->npers;
  memset(stats->inf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newinf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->postest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newpostest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  int32_t i;
  ext_timeline_info* const set=stats->ext_timeline-stats->tlshift;

  if(stats->nainfbins) {
    stats->ninfbins=1;

    for(i=nerase-1; i>=0; --i) {
      set[i].n=set[i].rsum=set[i].commpersum=0;
#ifdef NUMEVENTSSTATS
      set[i].neventssum=0;
#endif
      memset(set[i].ngeninfs,0,stats->nainfbins*sizeof(uint64_t));
    }

  } else memset(stats->ext_timeline-stats->tlshift,0,nerase*sizeof(ext_timeline_info));

  stats->extinction=true;
  stats->maxedoutmintimeindex=INT32_MAX;
  stats->tlshift=0;

#ifdef CT_OUTPUT
  stats->nctentries=0;
  stats->curctid=0;
#endif
}

/**
 * @brief Performs task at the end of a path simulation.
 *
 * This function must be called after the simulation of each path.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void std_stats_path_end(sim_vars* sv)
{
  const int32_t maxedoutmintimeindex=((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex;

  int32_t i,j;
  std_summary_stats* sss=(std_summary_stats*)sv->dataptr;
  uint32_t tnvpers;

  if(maxedoutmintimeindex<INT32_MAX) {
    tnvpers=(maxedoutmintimeindex<sss->npers?maxedoutmintimeindex+1:sss->npers)+sss->tlshift;

  } else {
    tnvpers=sss->npers+sss->tlshift;
  }

  //printf("Shift %i/%i, tnpers %u/%u\n",sss->tlshift,sss->tlshifta,tnvpers,sss->tnpersa);

  ext_timeline_info* const et=sss->ext_timeline-sss->tlshift;

  //for(i=-sss->tlshift; i<(int32_t)tnvpers; ++i) printf("rsum[%i]=%u %u\n",i,et[i].rsum,et[i].n);

  //for(i=tnvpers-sss->tlshift-1; i>=-sss->tlshift; --i) printf("newinf_timeline[%i]=%i\n",i,sss->newinf_timeline[i]);
  //for(i=sss->npers-1; i>=0; --i) printf("inf_timeline[%i]=%i\n",i,sss->inf_timeline[i]);

  if(sss->ninfbins) {

    for(i=tnvpers-2; i>=0; --i) {
      //printf("et[%i].n (%u) += %u\n",i,et[i].n,et[i+1].n);
      et[i].n+=et[i+1].n;
      et[i].rsum+=et[i+1].rsum;
      et[i].commpersum+=et[i+1].commpersum;
#ifdef NUMEVENTSSTATS
      et[i].neventssum+=et[i+1].neventssum;
#endif
      //printf("Address et[%i].ngeninfs=%p\n",i-sss->tlshift,et[i].ngeninfs);

      for(j=sss->ninfbins-1; j>=0; --j) {
	//printf("et[%i].ngeninfs[%i] (%lu) += %lu\n",i-sss->tlshift,j,et[i].ngeninfs[j],et[i+1].ngeninfs[j]);
	et[i].ngeninfs[j]+=et[i+1].ngeninfs[j];
	//if(i==0) printf("et[%i].ngeninfs[%i] = %lu\n",i-sss->tlshift,j,et[i].ngeninfs[j]);
      }
    }

  } else {

    for(i=tnvpers-2; i>=0; --i) {
      et[i].n+=et[i+1].n;
      et[i].rsum+=et[i+1].rsum;
      et[i].commpersum+=et[i+1].commpersum;
#ifdef NUMEVENTSSTATS
      et[i].neventssum+=et[i+1].neventssum;
#endif
    }
  }

  if(sss->timerelfirstpostestresults) {
    uint32_t const* const nptt=sss->newpostest_timeline-sss->tlshift;
    int32_t k,l;
    int32_t z;
    bool single;

    sss->tlpptnvpers=sss->tlppt0idx=sss->tlppnnpers=0;

    for(z=0; z<tnvpers; ++z) if(nptt[z]) {
    //z=sss->tlshift;
    //while(true) {
      sss->tlppt0idx=z-sss->tlshift;
      sss->tlppnnpers=(uint32_t)ceil(z*0.5);
      sss->tlpptnvpers=(uint32_t)ceil((tnvpers-z)*0.5)+sss->tlppnnpers;

      if(sss->maxedoutmintimeindex<INT32_MAX) sss->maxedoutmintimeindex=(int32_t)floor((sss->maxedoutmintimeindex-sss->tlppt0idx)*0.5);
      sss->extinction_time-=sss->tlppt0idx/(double)sss->nbinsperunit;
      //printf("Before merging: %i negatives and %u total. First positive test found at %i\n",sss->tlshift,tnvpers,z);
      //printf("First positive test found at %i => %i negative periods and %u total valid periods. maxedoutmintimeindex set to %i\n",sss->tlppt0idx,sss->tlppnnpers,sss->tlpptnvpers,sss->maxedoutmintimeindex);
      break;
    }

    if(sss->tlpptnvpers) {
      //Check if the last index is even. If it is, there is a single fine bin in
      //the last merged bin
      single=!((tnvpers-1-z)%2);
      j=-sss->tlppnnpers+sss->tlpptnvpers-1-single;

      for(k=0; k<=j; ++k) {
	i=sss->tlppt0idx+2*k;
	l=sss->tlppt0idx+k;
	//printf("[%i] = [%i] (%u) + [%i] (%u)\n",l,i,sss->newpostest_timeline[i],i+1,sss->newpostest_timeline[i+1]);
	sss->inf_timeline[l]=(sss->inf_timeline[i]>sss->inf_timeline[i+1]?sss->inf_timeline[i]:sss->inf_timeline[i+1]);
	sss->newinf_timeline[l]=sss->newinf_timeline[i]+sss->newinf_timeline[i+1];
	sss->newpostest_timeline[l]=sss->newpostest_timeline[i]+sss->newpostest_timeline[i+1];
      }

      if(single) {
	i=-sss->tlshift+tnvpers;
	l=sss->tlppt0idx+j+1;
	//printf("[%i] = [%i] (%u)\n",l,i,sss->newpostest_timeline[i]);
	sss->inf_timeline[l]=sss->inf_timeline[i];
	sss->newinf_timeline[l]=sss->newinf_timeline[i];
	sss->newpostest_timeline[l]=sss->newpostest_timeline[i];
      } 

      //Check if the first (negative) index is odd. If it is, there is a single fine bin in
      //the first merged bin
      single=z%2;
      j=-sss->tlppnnpers+single;

      for(k=-1; k>=j; --k) {
	i=sss->tlppt0idx+2*k;
	l=sss->tlppt0idx+k;
	//printf("[%i] = [%i] (%u) + [%i] (%u)\n",l,i,sss->newpostest_timeline[i],i+1,sss->newpostest_timeline[i+1]);
	sss->inf_timeline[l]=(sss->inf_timeline[i]>sss->inf_timeline[i+1]?sss->inf_timeline[i]:sss->inf_timeline[i+1]);
	sss->newinf_timeline[l]=sss->newinf_timeline[i]+sss->newinf_timeline[i+1];
	sss->newpostest_timeline[l]=sss->newpostest_timeline[i]+sss->newpostest_timeline[i+1];
      }

      if(single) {
	i=-sss->tlshift;
	l=sss->tlppt0idx+j-1;
	//printf("[%i] = [%i] (%u)\n",l,i,sss->newpostest_timeline[i]);
	sss->inf_timeline[l]=sss->inf_timeline[i];
	sss->newinf_timeline[l]=sss->newinf_timeline[i];
	sss->newpostest_timeline[l]=sss->newpostest_timeline[i];
      } 
    }

  } else {
    sss->tlppt0idx=0;
    sss->tlppnnpers=sss->tlshift;
    sss->tlpptnvpers=tnvpers;
  }
}

/**
 * @brief Frees memory used by the standard summary statistics.
 *
 * This function must be called once to free the memory that was allocated by
 * std_stats_init.
 *
 * @param stats: Pointer to the standard summary statistics.
 **/
void std_stats_free(std_summary_stats* stats);

inline static void std_stats_pri_init(sim_vars* sv, infindividual* ii) {
  //We have to use curii here!
  if((int)(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time) < ((std_summary_stats*)sv->dataptr)->npers && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
    DEBUG_PRINTF("Pri inf at %i\n",(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time));
    ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
  }
}

/**
 * @brief Allocates memory for a new primary individual.
 *
 * This function must be assigned to the simulation engine through a call of
 * sim_set_pri_init_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Infectious individuals.
 **/
inline static void std_stats_pri_init_rel(sim_vars* sv, infindividual* ii)
{
  const int32_t newshift=ceil(((std_summary_stats*)sv->dataptr)->nbinsperunit*(-ii->end_comm_period+(ii->comm_period+ii->latent_period)));

  if(newshift > ((std_summary_stats*)sv->dataptr)->tlshift) {
    std_summary_stats* stats=(std_summary_stats*)sv->dataptr;
    stats->tlshift=newshift;

    if(newshift > stats->tlshifta) {
      const uint32_t dshift=newshift-stats->tlshifta;
      const uint32_t newsize=newshift+stats->npers;

      uint32_t* newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->inf_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->inf_timeline-stats->tlshifta);
      stats->inf_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->newinf_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->newinf_timeline-stats->tlshifta);
      stats->newinf_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->postest_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->postest_timeline-stats->tlshifta);
      stats->postest_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->newpostest_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->newpostest_timeline-stats->tlshifta);
      stats->newpostest_timeline=newarray+newshift;

      ext_timeline_info* newarray_ext=malloc(newsize*sizeof(ext_timeline_info));
      memset(newarray_ext,0,dshift*sizeof(ext_timeline_info));
      memcpy(newarray_ext+dshift,stats->ext_timeline-stats->tlshifta,stats->tnpersa*sizeof(ext_timeline_info));
      free(stats->ext_timeline-stats->tlshifta);

      int32_t i;

      if(stats->nainfbins) {
	uint64_t* newarray64;

	for(i=dshift-1; i>=0; --i) {
	  newarray64=(uint64_t*)malloc(stats->nainfbins*sizeof(uint64_t));
	  memset(newarray64,0,stats->nainfbins*sizeof(uint64_t));
	  newarray_ext[i].ngeninfs=newarray64;
	  //printf("ext_timeline[%i].ngeninfs=%p\n",i-newshift,newarray64);
	}
      }
      stats->ext_timeline=newarray_ext+newshift;

      stats->tlshifta=newshift;
      stats->tnpersa=newsize;
    }
  }

  //We have to use curii here!
  if((int)(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time) < ((std_summary_stats*)sv->dataptr)->npers && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
    DEBUG_PRINTF("Pri inf at %i\n",(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time));
    ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
  }
}

/**
 * @brief Allocates memory for the user-defined function data.
 *
 * This function must be assigned to the simulation engine through a call of
 * sim_set_ii_alloc_proc_func.
 *
 * @param ii: Infectious individuals.
 **/
inline static void std_stats_ii_alloc(infindividual* ii){
#ifdef CT_OUTPUT
  ii->dataptr=malloc(3*sizeof(uint32_t));
#else
  ii->dataptr=malloc(sizeof(uint32_t));
#endif
}

#ifdef CT_OUTPUT
inline static void std_stats_add_ct_entry(std_summary_stats* stats, const double postesttime, const double presymtime, const uint32_t id, const uint32_t pid, const uint32_t ntracedcts){
  ++(stats->nctentries);
  DEBUG_PRINTF("%s: %u %22.15e %22.15e %u %u %u\n",__func__, stats->nctentries, postesttime, presymtime, id, pid, ntracedcts);

  if(stats->nctentries==stats->nactentries) {
    stats->nactentries*=CTENTRIES_GROWFACT;
    stats->ctentries=(ctposinf**)realloc(stats->ctentries,stats->nactentries*sizeof(ctposinf*));

    int32_t i;
    for(i=stats->nctentries; i<stats->nactentries; ++i) stats->ctentries[i]=(ctposinf*)malloc(sizeof(ctposinf));
  }

  stats->ctentries[stats->nctentries-1]->postesttime=postesttime*1440; //time in minutes
  stats->ctentries[stats->nctentries-1]->presymtime=(isinf(presymtime)?INT32_MAX:presymtime*1440); //time in minutes
  stats->ctentries[stats->nctentries-1]->id=id;
  stats->ctentries[stats->nctentries-1]->pid=pid;
  stats->ctentries[stats->nctentries-1]->ntracedcts=ntracedcts;
  DEBUG_PRINTF("%s: Encoded: %i %i %u %u %u\n",__func__,  stats->ctentries[stats->nctentries-1]->postesttime, stats->ctentries[stats->nctentries-1]->presymtime, stats->ctentries[stats->nctentries-1]->id=id, stats->ctentries[stats->nctentries-1]->pid=pid, stats->ctentries[stats->nctentries-1]->ntracedcts);
}
#endif

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before tmax.  This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((uint32_t*)sv->curii->dataptr)[2]+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",(((uint32_t*)sv->curii->dataptr)[2]));
#endif

  if(sv->curii->ninfections) {
    ((uint32_t*)sv->curii->dataptr)[0]+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((uint32_t*)sv->curii->dataptr)[0]);

    if(sv->curii->event_time < sv->pars.tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
      ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
      return true;
    }
  }
  return false;
}

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before tmax. If the incremented bin for the total infection timeline
 * exceeds nimax, then set the extinction for the current path to false and
 * update the value for the minimum time index where the maximum number of
 * infectious individuals was exceeded if required. This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event_nimax(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((uint32_t*)sv->curii->dataptr)[2]+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",(((uint32_t*)sv->curii->dataptr)[2]));
#endif

  if(sv->curii->ninfections) {
    ((uint32_t*)sv->curii->dataptr)[0]+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((uint32_t*)sv->curii->dataptr)[0]);

    if(sv->curii->event_time < sv->pars.tmax) {
      const int eti=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time);

      if(sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {

	if(((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+sv->curii->ninfections < ((std_summary_stats*)sv->dataptr)->nimax) {
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;

	} else if(((std_summary_stats*)sv->dataptr)->newinf_timeline[eti] < ((std_summary_stats*)sv->dataptr)->nimax) {
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->newinf_timeline[eti],((std_summary_stats*)sv->dataptr)->nimax);

	} else {
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->newinf_timeline[eti],((std_summary_stats*)sv->dataptr)->nimax);
	  return false;
	}
	return true;
      }
    }
  }
  return false;
}

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before tmax. If the bin for the recent positive test result timeline
 * exceeds npostestmax, then set the extinction for the current path to false and
 * update the value for the minimum time index where the maximum number of
 * recent positive test results was exceeded if required. This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event_npostestmax(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((uint32_t*)sv->curii->dataptr)[2]+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",(((uint32_t*)sv->curii->dataptr)[2]));
#endif

  if(sv->curii->ninfections) {
    ((uint32_t*)sv->curii->dataptr)[0]+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((uint32_t*)sv->curii->dataptr)[0]);

    if(sv->curii->event_time < sv->pars.tmax) {
      const int eti=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time);

      if(sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {

	if(((std_summary_stats*)sv->dataptr)->postest_timeline[eti] < ((std_summary_stats*)sv->dataptr)->npostestmax)
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;

	else {
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  {
	    DEBUG_PRINTF("maxedoutmintimeindex reduced from %i to %i\n",((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex,eti);
	    ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  }
	  DEBUG_PRINTF("npostestmax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->postest_timeline[eti],((std_summary_stats*)sv->dataptr)->npostestmax);
	  return false;
	}
	return true;
      }
    }
  }
  return false;
}

inline static void std_stats_fill_newpostest(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  if(ii->commpertype&ro_commper_true_positive_test) {
    const int trt=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*(ii->end_comm_period+sv->pars.tdeltat));

    DEBUG_PRINTF("New pos test at %i\n",trt);

    if(trt<((std_summary_stats*)sv->dataptr)->npers) ++(((std_summary_stats*)sv->dataptr)->newpostest_timeline[trt]);

#ifdef CT_OUTPUT
    ((uint32_t*)ii->dataptr)[1]=++(((std_summary_stats*)sv->dataptr)->curctid);
    DEBUG_PRINTF("True positive test with ID %u\n",((uint32_t*)ii->dataptr)[1]);
    ((uint32_t*)ii->dataptr)[2]=0;
#endif
    int32_t i;
    const int32_t end_comm_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*(ii->end_comm_period+sv->pars.tdeltat+((std_summary_stats*)sv->dataptr)->npostestmaxnpers));

    for(i=(end_comm_per_i>=((std_summary_stats*)sv->dataptr)->npers?((std_summary_stats*)sv->dataptr)->npers-1:end_comm_per_i); i>=trt; --i) ++(((std_summary_stats*)sv->dataptr)->postest_timeline[i]);
  }
}

/**
 * @brief Initialises the processing variable for a new infectious individual.
 *
 * This function initialises the number of infections from a new infectious
 * individual to zero. It is only called for infectious individuals that
 * participate to a non-zero number of transmission events, except for a primary
 * individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Infectious individual.
 * @param parent: Infectious individual's parent.
 **/
inline static void std_stats_new_inf(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  ((uint32_t*)ii->dataptr)[0]=0;
  DEBUG_PRINTF("%s: Number of infections initialized to 0.\n",__func__);

  //++(*(uint32_t*)(ii-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(ii-1)->dataptr);
  DEBUG_PRINTF("%s\n",__func__);
  std_stats_fill_newpostest(sv, ii, parent);
}

inline static void std_stats_fill_inf_ext_n(sim_vars* sv, infindividual* ii)
{
  const double start_comm_per=ii->end_comm_period-ii->comm_period;
  const int32_t start_latent_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*(ii->end_comm_period-(ii->comm_period+ii->latent_period)));
  const int32_t end_comm_per_i=(((std_summary_stats*)sv->dataptr)->nbinsperunit*ii->end_comm_period >= ((std_summary_stats*)sv->dataptr)->npers ? ((std_summary_stats*)sv->dataptr)->npers-1 : floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*ii->end_comm_period));
  int32_t i;

  if(start_comm_per<sv->pars.tmax) {
    const int32_t start_comm_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*start_comm_per);
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].rsum+=((uint32_t*)ii->dataptr)[0];
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].n+=1;
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].commpersum+=ii->comm_period;
#ifdef NUMEVENTSSTATS
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].neventssum+=ii->nevents;
#endif
  }

  for(i=start_latent_per_i; i<=end_comm_per_i; ++i) {
    ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
  }
}

/**
 * @brief Process the statistics after the last transmission event for an
 * infectious individual that participates to some transmission events.
 *
 * After an infectious infectious individual participated to its last
 * transmission event, this function adds the total number of infections from
 * the individual to the R sum, the individual's communicable period to the sum
 * of communicable periods and the number of transmission events this individual
 * participated to to the sum of events. If the individual was still infectious
 * at tmax, then the path is set to not go extinct. Otherwise, the
 * extinction time for the generating path is updated if its current value was
 * exceeded, using the current event time for the individual's parent and the
 * individual's communicable period. The timeline for the number of infectious
 * individuals is also updated.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Infectious individual.
 * @param parent: Infectious individual's parent.
 **/
inline static void std_stats_end_inf(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  DEBUG_PRINTF("Number of infections was %u\n",((uint32_t*)ii->dataptr)[0]);

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY), ((uint32_t*)ii->dataptr)[1], ((ii->commpertype&ro_commper_int)?((uint32_t*)parent->dataptr)[1]:0), ((uint32_t*)ii->dataptr)[2]);
#endif

  //If truncated by tmax
  if(ii->commpertype&ro_commper_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

  else {
    //++((std_summary_stats*)sv->dataptr)->n_ended_infections;

    if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->extinction_time) ((std_summary_stats*)sv->dataptr)->extinction_time=ii->end_comm_period;
  }

  std_stats_fill_inf_ext_n(sv, ii);
}

/**
 * @brief Process the statistics after the last transmission event for an
 * infectious individual that participates to some transmission events and
 * record the number of infections generated by each infectious individual.
 *
 * In addition to calling the function std_stats_end_inf, this function records
 * the number of infections generated by each infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Infectious individual.
 * @param parent: Infectious individual's parent.
 **/
inline static void std_stats_end_inf_rec_ninfs(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  const double start_comm_per=ii->end_comm_period-ii->comm_period;

  if(start_comm_per<sv->pars.tmax) {
    std_summary_stats* const stats=(std_summary_stats*)sv->dataptr;

    if(((uint32_t*)ii->dataptr)[0] >= stats->ninfbins) {
      stats->ninfbins=((uint32_t*)ii->dataptr)[0]+1;

      if(stats->ninfbins>stats->nainfbins) {
	ext_timeline_info* const set=stats->ext_timeline-stats->tlshifta;
	int32_t i;
	const uint32_t dnbins=stats->ninfbins-stats->nainfbins;

	for(i=stats->tnpersa-1; i>=0; --i) {
	  set[i].ngeninfs=(uint64_t*)realloc(set[i].ngeninfs,stats->ninfbins*sizeof(uint64_t));
	  //printf("ext_timeline[%i].ngeninfs=%p\n",i-stats->tlshifta,set[i].ngeninfs);
	  memset(set[i].ngeninfs+stats->nainfbins,0,dnbins*sizeof(uint64_t));
	}
	stats->nainfbins=stats->ninfbins;
      }
    }
    const int32_t start_comm_per_i=floor(stats->nbinsperunit*start_comm_per);
    ++(stats->ext_timeline[start_comm_per_i].ngeninfs[((uint32_t*)ii->dataptr)[0]]);
  }

  std_stats_end_inf(sv, ii, parent);
}

/**
 * @brief Process the statistics for an infectious individual that does not
 * participate to any transmission event.
 *
 * This function adds the individual's communicable period to the sum
 * of communicable periods. If the individual was still infectious
 * at tmax, then the path is set to not go extinct. Otherwise, the
 * extinction time for the generating path is updated if its current value was
 * exceeded, using the current event time for the individual's parent and the
 * individual's communicable period. The timeline for the number of infectious
 * individuals is also updated.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Pointer to the infectious individual.
 * @param parent: Pointer to the infectious individual's parent.
 **/
inline static void std_stats_noevent_new_inf(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  DEBUG_PRINTF("%s\n",__func__);
  DEBUG_PRINTF("Number of infections was 0\n");

  ((uint32_t*)ii->dataptr)[0]=0;

  std_stats_fill_newpostest(sv, ii, parent);

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY), ((uint32_t*)ii->dataptr)[1], ((ii->commpertype&ro_commper_int)?((uint32_t*)parent->dataptr)[1]:0), ((uint32_t*)ii->dataptr)[2]);
#endif

  //If truncated by tmax
  if(ii->commpertype&ro_commper_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

  else {
    //++((std_summary_stats*)sv->dataptr)->n_ended_infections;

    if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->extinction_time) ((std_summary_stats*)sv->dataptr)->extinction_time=ii->end_comm_period;
  }

  std_stats_fill_inf_ext_n(sv, ii);
}

/**
 * @brief Process the statistics for an infectious individual that does not
 * participate to any transmission event and record the number of infections
 * generated by each infectious individual.
 *
 * In addition to calling the function std_stats_noevent_new_inf, this function records
 * the number of infections generated by each infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Pointer to the infectious individual.
 * @param parent: Pointer to the infectious individual's parent.
 **/
inline static void std_stats_noevent_new_inf_rec_ninfs(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  const double start_comm_per=ii->end_comm_period-ii->comm_period;

  if(start_comm_per<sv->pars.tmax) {
    const int32_t start_comm_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*start_comm_per);
    ++(((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].ngeninfs[0]);
  }

  std_stats_noevent_new_inf(sv, ii, parent);
}

#endif
