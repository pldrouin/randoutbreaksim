/**
 * @file standard_summary_stats.h
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _STANDARD_SUMMARY_STATS_
#define _STANDARD_SUMMARY_STATS_

#include <stdlib.h>
#include <stdbool.h>

#include "infindividual.h"
#include "simulation.h"

#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__) //!< Debug print function

/**
 * Simulation-level standard summary statistics data struct.
 */
typedef struct
{
  double extinction_time;	//!< Path extinction time, if any.
  double commpersum;		//!< Sum of communicable periods for all infectious individuals whose communicable period does not occur after tmax.
  uint32_t* inf_timeline;	//!< For each integer interval between 0 and floor(tmax), the number of individuals that are infectious at some point in this interval (lower bound included, upper bound excluded).
  uint32_t* totinf_timeline;	//!< For each integer interval between 0 and floor(tmax), the number of individuals that get infected at some point in this interval (lower bound included, upper bound excluded). For the last interval, it includes the infectious that occur between floor(tmax) and tmax.
  uint32_t npers;		//!< Number of integer intervals (floor(tmax)+1)
  uint32_t rsum;		//!< Sum of the number of individuals that get infected during the simulation by the infectious individuals whose last transmission event does not occur after tmax.
  uint32_t neventssum;		//!< Sum of the number of transmission events for all infectious individuals whose communicable period does not occur after tmax.
  uint32_t nimax;               //!< Maximum number of infectious individuals for a given integer inerval between 0 and floor(tmax). Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  uint32_t nimaxedoutmintimeindex; //!< Minimum time index which maxed out the allowed number of infected individuals.
  //uint32_t n_ended_infections;
  bool extinction;		//!< Set to true if extinction does not occur before or at tmax.
} std_summary_stats;

/**
 * @brief Initialises the standard summary statistics.
 *
 * This function must be called once to initialise the standard summary
 * statistics, if used, as well as to allocate required memory in the
 * simulation variables.
 *
 * @param sv: Pointer to the simulation variables were memory for the
 * user-defined function data will be allocated.
 * @param stats: Pointer to the standard summary statistics that will be
 * initialised.
 */
void std_stats_init(sim_vars* sv, std_summary_stats* stats);

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
  stats->extinction_time=0;
  stats->commpersum=0;
  memset(stats->inf_timeline,0,stats->npers*sizeof(uint32_t));
  memset(stats->totinf_timeline,0,stats->npers*sizeof(uint32_t));
  stats->rsum=0;
  stats->neventssum=0;
  stats->extinction=true;
  stats->nimaxedoutmintimeindex=UINT32_MAX;
}

/**
 * @brief Allocates additional memory for the user-defined function data.
 *
 * This function must be assigned to the simulation engine through a call of
 * sim_set_increase_layers_proc_func.
 *
 * @param iis: First newly allocated layer element of the infectious individuals array.
 * @param n: Number of layers that have been added
 * */
void std_stats_increase_layers(infindividual* iis, uint32_t n);

/**
 * @brief Frees memory used by the standard summary statistics.
 *
 * This function must be called once to free the memory that was allocated by
 * std_stats_init.
 *
 * @param sv: Pointer to the simulation variables were memory for the
 * user-defined function data will be freed.
 * @param stats: Pointer to the standard summary statistics.
 */
void std_stats_free(sim_vars* sv, std_summary_stats* stats);

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
 * @return true if the event time does not exceed tmax and false otherwise.
 * */
inline static bool std_stats_new_event(sim_vars* sv)
{
  (*(uint32_t*)sv->ii->dataptr)+=sv->ii->ninfections;
  DEBUG_PRINTF("Number of infections incremented to %u\n",*(uint32_t*)sv->ii->dataptr);

  if((int)sv->ii->event_time <= (int)sv->pars.tmax) {

    if(((std_summary_stats*)sv->dataptr)->totinf_timeline[(int)(sv->ii->event_time)] <= ((std_summary_stats*)sv->dataptr)->nimax)
      ((std_summary_stats*)sv->dataptr)->totinf_timeline[(int)(sv->ii->event_time)]+=sv->ii->ninfections;

    else {
      ((std_summary_stats*)sv->dataptr)->extinction=false;

      if((int)(sv->ii->event_time) < ((std_summary_stats*)sv->dataptr)->nimaxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->nimaxedoutmintimeindex=(int)(sv->ii->event_time);
      DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",(int)(sv->ii->event_time),((std_summary_stats*)sv->dataptr)->totinf_timeline[(int)(sv->ii->event_time)],((std_summary_stats*)sv->dataptr)->nimax);
      return false;
    }
  }

  return (sv->ii->event_time <= sv->pars.tmax);
}

/**
 * @brief Initialises the processing variable for a new infectious individual.
 *
 * This function initialises the number of infectins from a new infectious
 * individual to zero. It is only called for infectious individuals that
 * participate to a non-zero number of transmission events.
 *
 * @param inf: Pointer to the new infectious individual.
 * */
inline static void std_stats_new_inf(infindividual* inf)
{
  *(uint32_t*)inf->dataptr=0;
  //++(*(uint32_t*)(inf-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(inf-1)->dataptr);
  DEBUG_PRINTF("%s\n",__func__);
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
 * @param inf: Pointer to the infectious individual.
 * @param ptr: Pointer to the summary statistics.
 * */
inline static void std_stats_end_inf(infindividual* inf, void* ptr)
{
  DEBUG_PRINTF("Number of infections was %u\n",*(uint32_t*)inf->dataptr);
  ((std_summary_stats*)ptr)->rsum+=*(uint32_t*)inf->dataptr;
  ((std_summary_stats*)ptr)->commpersum+=inf->comm_period;
  ((std_summary_stats*)ptr)->neventssum+=inf->nevents;

  const double inf_end=(inf-1)->event_time+inf->comm_period;

  //If truncated by tmax
  if(inf->infectious_at_tmax) ((std_summary_stats*)ptr)->extinction=false;

  else {
    //++((std_summary_stats*)ptr)->n_ended_infections;

    if(inf_end > ((std_summary_stats*)ptr)->extinction_time) ((std_summary_stats*)ptr)->extinction_time=inf_end;
  }

  int i;
  const int end_comm_per=(inf_end >= ((std_summary_stats*)ptr)->npers ? ((std_summary_stats*)ptr)->npers-1 : (int)inf_end);

  for(i=(int)((inf-1)->event_time); i<=end_comm_per; ++i) ++(((std_summary_stats*)ptr)->inf_timeline[i]);
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
 * @param inf: Pointer to the infectious individual.
 * @param ptr: Pointer to the summary statistics.
 * */
inline static void std_stats_noevent_inf(infindividual* inf, void* ptr)
{
  DEBUG_PRINTF("%s\n",__func__);
  //++(*(uint32_t*)(inf-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(inf-1)->dataptr);
  DEBUG_PRINTF("Number of infections was 0\n");
  ((std_summary_stats*)ptr)->commpersum+=inf->comm_period;

  const double inf_end=(inf-1)->event_time+inf->comm_period;

  //If truncated by tmax
  if(inf->infectious_at_tmax) ((std_summary_stats*)ptr)->extinction=false;

  else {
    //++((std_summary_stats*)ptr)->n_ended_infections;

    if(inf_end > ((std_summary_stats*)ptr)->extinction_time) ((std_summary_stats*)ptr)->extinction_time=inf_end;
  }

  int i;
  const int end_comm_per=(inf_end >= ((std_summary_stats*)ptr)->npers ? ((std_summary_stats*)ptr)->npers-1 : (int)inf_end);

  for(i=(int)((inf-1)->event_time); i<=end_comm_per; ++i) ++(((std_summary_stats*)ptr)->inf_timeline[i]);
}

#endif
