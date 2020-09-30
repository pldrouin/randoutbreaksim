/**
 * @file standard_summary_stats.h
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _STANDARD_SUMMARY_STATS_
#define _STANDARD_SUMMARY_STATS_

#include <stdlib.h>
#include <stdbool.h>

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

/**
 * Simulation-level standard summary statistics data struct.
 */
typedef struct
{
  double extinction_time;	//!< Path extinction time, if any.
  double commpersum;		//!< Sum of communicable periods for all infectious individuals whose communicable period does not occur after tmax.
  uint32_t* inf_timeline;	//!< For each integer interval between 0 and floor(tmax), the number of individuals that are infected at some point in this interval (lower bound included, upper bound excluded).
  uint32_t* newinf_timeline;	//!< For each integer interval between 0 and floor(tmax), the number of individuals that get infected at some point in this interval (lower bound included, upper bound excluded). For the last interval, it includes the infectious that occur between floor(tmax) and tmax.
  uint32_t* newpostest_timeline;	//!< For each integer interval between 0 and floor(tmax), the number of individuals that receive a positive test result at some point in this interval with an individual whose communicable period is the alternate period (lower bound included, upper bound excluded). For the last interval, it includes positive test results that occur between floor(tmax) and tmax.
  uint64_t** ngeninfs;	        //!< Number of generated infections from each infectious individual.
  uint32_t* ninfbins;		//!< Number of allocated infectious individuals.
#ifdef CT_OUTPUT
  ctposinf** ctentries;
  uint32_t nctentries;
  uint32_t nactentries;
  int32_t curctid;
#endif
  uint32_t npers;		//!< Number of positive integer intervals
  int32_t timelineshift;       //!< Integral shift of the timeline origin (to allow for negative time bins). Corresponds also to the number of negative integer intervals
  uint32_t tnpersa;              //!< Total number of allocated integer intervals (negative+positive)
  uint32_t rsum;		//!< Sum of the number of individuals that get infected during the simulation by the infectious individuals whose last transmission event does not occur after tmax.
#ifdef NUMEVENTSSTATS
  uint32_t neventssum;		//!< Sum of the number of transmission events for all infectious individuals whose communicable period does not occur after tmax.
#endif
  uint32_t lmax;                //!< Maximum number of layers for the simulation. lmax=1 means only primary infectious individuals.
  uint32_t nimax;               //!< Maximum number of infectious individuals for a given integer interval between 0 and floor(tmax). Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  int32_t nimaxedoutmintimeindex; //!< Minimum time index which maxed out the allowed number of infected individuals.
  //uint32_t n_ended_infections;
  bool extinction;		//!< Set to true if extinction does not occur before or at tmax.
} std_summary_stats;

/**
 * @brief Initialises the standard summary statistics.
 *
 * This function must be called once to initialise the standard summary
 * statistics, if used.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ngeninfs: Pointer to an array of number of generated infections to be
 * filled. Ignored if NULL.
 * @param ninfbins: Pointer to the length of the ngeninfs array. Ignored if
 * NULL.
 */
void std_stats_init(sim_vars *sv, uint64_t** ngeninfs, uint32_t* ninfbins);

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
  stats->commpersum=0;
  memset(stats->inf_timeline-stats->timelineshift,0,stats->tnpersa*sizeof(uint32_t));
  memset(stats->newinf_timeline-stats->timelineshift,0,stats->tnpersa*sizeof(uint32_t));
  memset(stats->newpostest_timeline-stats->timelineshift,0,stats->tnpersa*sizeof(uint32_t));
  stats->rsum=0;
#ifdef NUMEVENTSSTATS
  stats->neventssum=0;
#endif
  stats->extinction=true;
  stats->nimaxedoutmintimeindex=UINT32_MAX;

#ifdef CT_OUTPUT
  stats->nctentries=0;
  stats->curctid=0;
#endif
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
inline static void std_stats_pri_init(sim_vars* sv, infindividual* ii)
{
  const int32_t newshift=ceil(-ii->end_comm_period+(ii->comm_period+ii->latent_period));

  if(newshift>((std_summary_stats*)sv->dataptr)->timelineshift) {
    std_summary_stats* stats=(std_summary_stats*)sv->dataptr;
    const uint32_t dshift=newshift-stats->timelineshift;
    const uint32_t newsize=newshift+stats->npers;
    uint32_t* newarray;

    newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
    memset(newarray,0,dshift*sizeof(uint32_t));
    memcpy(newarray+dshift,stats->inf_timeline-stats->timelineshift,stats->tnpersa*sizeof(uint32_t));
    free(stats->inf_timeline-stats->timelineshift);
    stats->inf_timeline=newarray+newshift;

    newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
    memset(newarray,0,dshift*sizeof(uint32_t));
    memcpy(newarray+dshift,stats->newinf_timeline-stats->timelineshift,stats->tnpersa*sizeof(uint32_t));
    free(stats->newinf_timeline-stats->timelineshift);
    stats->newinf_timeline=newarray+newshift;

    newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
    memset(newarray,0,dshift*sizeof(uint32_t));
    memcpy(newarray+dshift,stats->newpostest_timeline-stats->timelineshift,stats->tnpersa*sizeof(uint32_t));
    free(stats->newpostest_timeline-stats->timelineshift);
    stats->newpostest_timeline=newarray+newshift;

    stats->timelineshift=newshift;
    stats->tnpersa=newsize;
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
  stats->ctentries[stats->nctentries-1]->presymtime=presymtime*1440; //time in minutes
  stats->ctentries[stats->nctentries-1]->id=id;
  stats->ctentries[stats->nctentries-1]->pid=pid;
  stats->ctentries[stats->nctentries-1]->ntracedcts=ntracedcts;
  DEBUG_PRINTF("%s: Encoded: %i %i %u %u %u\n",__func__,  stats->ctentries[stats->nctentries-1]->postesttime, stats->ctentries[stats->nctentries-1]->presymtime, stats->ctentries[stats->nctentries-1]->id=id, stats->ctentries[stats->nctentries-1]->pid=pid, stats->ctentries[stats->nctentries-1]->ntracedcts);
}
#endif

/**
 * @brief Frees memory used by the standard summary statistics.
 *
 * This function must be called once to free the memory that was allocated by
 * std_stats_init.
 *
 * @param stats: Pointer to the standard summary statistics.
 **/
void std_stats_free(std_summary_stats* stats);

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time does not exceed tmax.  This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time does not exceed tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((uint32_t*)sv->curii->dataptr)[2]+=sv->curii->ntracednicts;
  DEBUG_PRINTF("ID %u: Successfully traced contacts incremented to %u from non-infected contacts\n",(((uint32_t*)sv->curii->dataptr)[1]),(((uint32_t*)sv->curii->dataptr)[2]));
#endif

  if(sv->curii->ninfections) {
    ((uint32_t*)sv->curii->dataptr)[0]+=sv->curii->ninfections;
    DEBUG_PRINTF("Number of infections incremented to %u\n",((uint32_t*)sv->curii->dataptr)[0]);

    //The following condition on event_time is looser than the one in the
    //returned value
    if((int)sv->curii->event_time <= (int)sv->pars.tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
      ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(sv->curii->event_time)]+=sv->curii->ninfections;
      return (sv->curii->event_time <= sv->pars.tmax);
    }
  }
  return false;
}

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time does not exceed tmax. If the incremented bin for the total infection timeline
 * exceeds nimax, then set the extinction for the current path to false and
 * update the value for the minimum time index where the maximum number of
 * infectious individuals was exceeded if required. This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time does not exceed tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event_nimax(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((uint32_t*)sv->curii->dataptr)[2]+=sv->curii->ntracednicts;
  DEBUG_PRINTF("ID %u: Successfully traced contacts incremented to %u from non-infected contacts\n",(((uint32_t*)sv->curii->dataptr)[1]),(((uint32_t*)sv->curii->dataptr)[2]));
#endif

  if(sv->curii->ninfections) {
    ((uint32_t*)sv->curii->dataptr)[0]+=sv->curii->ninfections;
    DEBUG_PRINTF("Number of infections incremented to %u\n",((uint32_t*)sv->curii->dataptr)[0]);

    if((int)sv->curii->event_time <= (int)sv->pars.tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
      const int eti=floor(sv->curii->event_time);

      if(((std_summary_stats*)sv->dataptr)->newinf_timeline[eti] <= ((std_summary_stats*)sv->dataptr)->nimax)
	((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;

      else {
	((std_summary_stats*)sv->dataptr)->extinction=false;

	if(eti < ((std_summary_stats*)sv->dataptr)->nimaxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->nimaxedoutmintimeindex=eti;
	DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->newinf_timeline[eti],((std_summary_stats*)sv->dataptr)->nimax);
	return false;
      }
      return (sv->curii->event_time <= sv->pars.tmax);
    }
  }
  return false;
}

inline static void std_stats_fill_newpostest(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  if(ii->commpertype&ro_commper_true_positive_test) {
    const int trt=floor(ii->end_comm_period+sv->pars.tdeltat);

    DEBUG_PRINTF("New pos test at %i\n",trt);

    if(trt<(int32_t)((std_summary_stats*)sv->dataptr)->npers) ++(((std_summary_stats*)sv->dataptr)->newpostest_timeline[trt]);

#ifdef CT_OUTPUT
    ((uint32_t*)ii->dataptr)[1]=++(((std_summary_stats*)sv->dataptr)->curctid);
    ((uint32_t*)ii->dataptr)[2]=0;
    DEBUG_PRINTF("Infectious individual ID is %u\n",((uint32_t*)ii->dataptr)[1]);
    DEBUG_PRINTF("ID %u: Successfully traced contacts initialized to 0\n",(((uint32_t*)ii->dataptr)[1]));
#endif
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

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_int) {
    ++(((uint32_t*)parent->dataptr)[2]);
    DEBUG_PRINTF("ID %u: Successfully traced contacts incremented to %u\n",(((uint32_t*)parent->dataptr)[1]),(((uint32_t*)parent->dataptr)[2]));
  }
#endif
  //++(*(uint32_t*)(ii-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(ii-1)->dataptr);
  DEBUG_PRINTF("%s\n",__func__);
  std_stats_fill_newpostest(sv, ii, parent);
}

/**
 * @brief Process the statistics after the last transmission event for an
 * infectious individual that participates to some transmission events.
 *
 * After an infectious infectious individual participated to its last
 * transmission event, this function adds the total number of infections from
 * the individual to the R sum, the individual's communicable period to the sum
 * of communicable periods and the number of transmission events t+=((uint32_t*)ii->dataptr)[1]his individual
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
  ((std_summary_stats*)sv->dataptr)->rsum+=((uint32_t*)ii->dataptr)[0];
  ((std_summary_stats*)sv->dataptr)->commpersum+=ii->comm_period;
#ifdef NUMEVENTSSTATS
  ((std_summary_stats*)sv->dataptr)->neventssum+=ii->nevents;
#endif

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY), ((uint32_t*)ii->dataptr)[1], ((ii->commpertype&ro_commper_int)?((uint32_t*)parent->dataptr)[1]:0), ((uint32_t*)ii->dataptr)[2]);
#endif

  //If truncated by tmax
  if(ii->commpertype&ro_commper_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

  else {
    //++((std_summary_stats*)sv->dataptr)->n_ended_infections;

    if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->extinction_time) ((std_summary_stats*)sv->dataptr)->extinction_time=ii->end_comm_period;
  }

  const int end_comm_per=(ii->end_comm_period >= (int32_t)((std_summary_stats*)sv->dataptr)->npers ? ((std_summary_stats*)sv->dataptr)->npers-1 : floor(ii->end_comm_period));
  int i=floor(ii->end_comm_period-(ii->comm_period+ii->latent_period));

  for(; i<=end_comm_per; ++i) ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
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
  if(*(uint32_t*)ii->dataptr >= *((std_summary_stats*)sv->dataptr)->ninfbins) {
     *((std_summary_stats*)sv->dataptr)->ngeninfs=(uint64_t*)realloc(*((std_summary_stats*)sv->dataptr)->ngeninfs,(*(uint32_t*)ii->dataptr+1)*sizeof(uint64_t));
     memset(*((std_summary_stats*)sv->dataptr)->ngeninfs+*((std_summary_stats*)sv->dataptr)->ninfbins,0,(*(uint32_t*)ii->dataptr+1-*((std_summary_stats*)sv->dataptr)->ninfbins)*sizeof(uint64_t));
    *((std_summary_stats*)sv->dataptr)->ninfbins=*(uint64_t*)ii->dataptr+1;
  }
  ++((*((std_summary_stats*)sv->dataptr)->ngeninfs)[*(uint32_t*)ii->dataptr]);

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

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_int) {
    ++(((uint32_t*)parent->dataptr)[2]);
    DEBUG_PRINTF("ID %u: Successfully traced contacts incremented to %u\n",(((uint32_t*)parent->dataptr)[1]),(((uint32_t*)parent->dataptr)[2]));
  }
#endif

  std_stats_fill_newpostest(sv, ii, parent);

  ((std_summary_stats*)sv->dataptr)->commpersum+=ii->comm_period;

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY), ((uint32_t*)ii->dataptr)[1], ((ii->commpertype&ro_commper_int)?((uint32_t*)parent->dataptr)[1]:0), ((uint32_t*)ii->dataptr)[2]);
#endif

  //If truncated by tmax
  if(ii->commpertype&ro_commper_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

  else {
    //++((std_summary_stats*)sv->dataptr)->n_ended_infections;

    if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->extinction_time) ((std_summary_stats*)sv->dataptr)->extinction_time=ii->end_comm_period;
  }

  const int end_comm_per=(ii->end_comm_period >= (int32_t)((std_summary_stats*)sv->dataptr)->npers ? ((std_summary_stats*)sv->dataptr)->npers-1 : floor(ii->end_comm_period));
  int i=floor(ii->end_comm_period-(ii->comm_period+ii->latent_period));

  for(; i<=end_comm_per; ++i) ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
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
  ++((*((std_summary_stats*)sv->dataptr)->ngeninfs)[0]);

  std_stats_noevent_new_inf(sv, ii, parent);
}

#endif
