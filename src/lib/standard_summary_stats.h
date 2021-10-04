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
  int32_t id;
  int32_t pid;
  uint32_t ntracedcts;
} ctposinf;
#endif

typedef struct
{
  uint32_t   n;	               //!< Number of infectious creating infections in the timeline bin.
  uint32_t rsum;	       //!< Sum of the number of individuals that get infected in the timeline bin.
  uint64_t r2sum;	       //!< Sum of the square of the number of individuals that get infected in the timeline bin.
  double commpersum;   	       //!< Sum of communicable periods for all infectious individuals whose communicable period start before tmax in the timeline bin.
#ifdef NUMEVENTSSTATS
  uint32_t neventssum;	       //!< Sum of the number of transmission events for all infectious individuals whose communicable period occurs before tmax in the timeline bin.
#endif
#ifdef OBSREFF_OUTPUT
  uint32_t nobs;               //!< Number of observable (who test positive) infectious creating infections in the timeline bin
  uint32_t robssum;            //!< Sum of the number of observable (who test positive) individuals that get infected in the timeline bin.
  uint64_t robs2sum;           //!< Sum of the square of the number of observable (who test positive) individuals that get infected in the timeline bin.
#endif
  uint64_t* ngeninfs;	       //!< Number of generated infections from each infectious individual in the timeline bin.
} ext_timeline_info;

/**
 * Simulation-level standard summary statistics data struct.
 */
typedef struct
{
  double extinction_time;	//!< *Path extinction time, if any. 
  double abs_tmax;		//!< Absolute maximum simulation time, including period before post-processing origin.
  double first_pos_test_results_time; //!< Absolute time for the first positive test results. Only populated with time_rel_first_pos_test_results.
#ifdef OBSREFF_OUTPUT
  infindividual iibuf;
#endif
  uint32_t* inf_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that are infected, but not isolated, at some point in this interval. First index is -tlshift.
  uint32_t* newinf_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that get infected at some point in this interval. First index is -tlshift.
  uint32_t* postest_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that have a recent positive test result at the end of this interval. First index is -tlshift.
  uint32_t* newpostest_timeline;//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that receive a positive test result at some point in this interval. First index is -tlshift.
  uint32_t* pp_inf_timeline;	//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that are infected, but not isolated, at some point in this interval. First index is -tlshift.
  uint32_t* pp_newinf_timeline;	//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that get infected at some point in this interval. For the last interval, First index is -tlshift.
  uint32_t* pp_newpostest_timeline;//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that receive a positive test result at some point in this interval. First index is -tlshift.
#ifdef SEC_INF_TIMELINES
  uint32_t* secinf_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that are infected, but not isolated, at some point in this interval. First index is -tlshift.
  uint32_t* newsecinf_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that get infected at some point in this interval. First index is -tlshift.
  uint32_t* secpostest_timeline;	//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that have a recent positive test result at the end of this interval. First index is -tlshift.
  uint32_t* newsecpostest_timeline;//!< For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that receive a positive test result at some point in this interval. First index is -tlshift.
  uint32_t* pp_secinf_timeline;	//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that are infected, but not isolated, at some point in this interval. First index is -tlshift.
  uint32_t* pp_newsecinf_timeline;	//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that get infected at some point in this interval. For the last interval, First index is -tlshift.
  uint32_t* pp_newsecpostest_timeline;//!< Post-processing timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the number of individuals that receive a positive test result at some point in this interval. First index is -tlshift.
#endif
  ext_timeline_info* pp_ext_timeline;//!< Post-processing extended timeline. For each integer interval between 0 and nbinsperunit*abs_tmax-1, the extended information. First index is -tlshift.
  ext_timeline_info* ext_timeline;     //!< Extended timeline for parameters that cannot correctly be calculated when a dynamic time cur is used. First index is -tlshift.
  uint32_t nainfbins;		//!< Number of allocated infectious individual bins.
  uint32_t ninfbins;		//!< Number of used infectious individual bins.
#ifdef CT_OUTPUT
  ctposinf** ctentries;
  uint32_t nctentries;
  uint32_t nactentries;
  int32_t curctid;
#endif
  int32_t abs_maxnpers;		//!< Absolute maximum number of potential positive integer intervals for the simulation. This number can decrease during the simulation of a path (=nbinsperunit*abs_tmax).
  int32_t abs_npers;		//!< Actual number of absolute simulated positive integer intervals. This number can only increase during the simulation of a path (=floor(nbinsperunit*(end_comm_period+tdeltat+npostestmaxnunits)) for positive tests and floor(nbinsperunit*end_comm_period) for negative tests when using time_rel_first_pos_test).
  int32_t npers;		//!< Maximum number of desired positive integer intervals (=nbinsperunit*tmax).
  int32_t nbinsperunit;		//!< Number of timeline bins per unit of time.
  int32_t  tlshifta;            //!< *Allocated integral shift of the timeline origin (to allow for negative time bins). Corresponds also to the number of negative integer intervals
  int32_t  tlshift;             //!< *Integral shift of the timeline origin (to allow for negative time bins). Corresponds also to the number of negative integer intervals
  uint32_t tnpersa;             //!< *Total number of allocated integer intervals (negative+positive)
  int32_t tlppnnpers;		//!< Timeline post-processing number of negative periods
  uint32_t tlpptnvpers;		//!< Timeline post-processing total number of valid periods
  uint32_t lmax;                //!< Maximum number of layers for the simulation. lmax=1 means only primary infectious individuals.
  uint32_t nimax;               //!< Maximum number of infectious individuals for a given integer interval between 0 and nbinsperunit*abs_tmax-1. Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  uint32_t npostestmax;         //!< Maximum number of positive test results during an interval of duration npostestmaxnunits for each individual that starts when the test results are received. Extinction is set to false and the simulation does not proceed further if this maximum is exceeded.
  uint32_t npostestmaxnunits;    //!< Interval duration for the maximum number of positive test results
  int32_t maxedoutmintimeindex; //!< *Minimum time index which maxed out the allowed number of infected individuals or positive test results.
  //uint32_t n_ended_infections;
  bool extinction;		//!< *Set to true if extinction does not occur before abs_tmax.
} std_summary_stats;

typedef struct {
  uint32_t ninf;
#ifdef OBSREFF_OUTPUT
  uint32_t nobsinf;
#endif
#ifdef CT_OUTPUT
  int32_t id;
  uint32_t ntracedcts;
#endif
} std_stats_inf_data;

/**
 * @brief Initialises the standard summary statistics.
 *
 * This function must be called once to initialise the standard summary
 * statistics, if used.
 *
 * @param sv: Pointer to the simulation variables.
 */
void std_stats_init(sim_vars *sv, const uint32_t nbinsperunit, bool ngeninfs);

/**
 * @brief Initialises elements of the standard summary statistics
 * before a path simulation.
 *
 * This function must be called before the simulation of each path to initialise
 * elements of the standard summary statistics, if used
 *
 * @param sv: Pointer to the simulation variables.
 * elements will be initialised.
 */
inline static void std_stats_path_init(sim_vars* sv)
{
  std_summary_stats* stats=(std_summary_stats*)sv->dataptr;
  stats->extinction_time=-INFINITY;
  const uint32_t nerase=stats->tlshift+stats->abs_npers;
  memset(stats->inf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newinf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->postest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newpostest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
#ifdef SEC_INF_TIMELINES
  memset(stats->secinf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newsecinf_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->secpostest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
  memset(stats->newsecpostest_timeline-stats->tlshift,0,nerase*sizeof(uint32_t));
#endif
  int32_t i;
  ext_timeline_info* const set=stats->ext_timeline-stats->tlshift;

  if(stats->nainfbins) {
    stats->ninfbins=1;

    for(i=nerase-1; i>=0; --i) {
      set[i].n=set[i].rsum=set[i].r2sum=set[i].commpersum=0;
#ifdef NUMEVENTSSTATS
      set[i].neventssum=0;
#endif
#ifdef OBSREFF_OUTPUT
      set[i].nobs=set[i].robssum=set[i].robs2sum=0;
#endif
      memset(set[i].ngeninfs,0,stats->nainfbins*sizeof(uint64_t));
    }

  } else memset(stats->ext_timeline-stats->tlshift,0,nerase*sizeof(ext_timeline_info));

  stats->extinction=true;
  stats->maxedoutmintimeindex=INT32_MAX;

  if(sv->pars.timetype==ro_time_first_pos_test_results) {
    stats->abs_maxnpers=INT32_MAX;
    stats->abs_tmax=((double)stats->abs_maxnpers)/stats->nbinsperunit;
    stats->first_pos_test_results_time=INFINITY;
    stats->abs_npers=0;

  } else stats->tlshift=0;

#ifdef CT_OUTPUT
  stats->nctentries=0;
  stats->curctid=0;
#endif
  if(stats->ext_timeline[0].n!=0) {
    printf("ext: %u, %u, %u, %u\n",stats->ext_timeline[0].n,nerase,stats->tlshift,stats->abs_npers);
  }
}

/**
 * @brief Performs task at the end of a path simulation.
 *
 * This function must be called after the simulation of each path.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the path is valid and false otherwise.
 */
inline static bool std_stats_path_end(sim_vars* sv)
{
  const int32_t maxedoutmintimeindex=((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex;

  int32_t i,j;
  std_summary_stats* stats=(std_summary_stats*)sv->dataptr;
  uint32_t tnvpers;
  bool includepath;

  if(maxedoutmintimeindex<INT32_MAX) {
    tnvpers=(maxedoutmintimeindex<stats->abs_npers?maxedoutmintimeindex+1:stats->abs_npers)+stats->tlshift;

  } else tnvpers=stats->abs_npers+stats->tlshift;

  ext_timeline_info* const et=stats->ext_timeline-stats->tlshift;

  if(sv->pars.timetype==ro_time_first_pos_test_results) {

    if(isinf(stats->first_pos_test_results_time)) return false;
    includepath=true;

    int32_t k;
    bool single;
    int32_t tlppt0idx=floor(stats->nbinsperunit*stats->first_pos_test_results_time);

    stats->tlppnnpers=(uint32_t)ceil(tlppt0idx*0.5);

    //if(tnvpers-tlppt0idx<stats->npers && stats->extinction==false && maxedoutmintimeindex!=tnvpers) {
    //  printf("tnvpers=%u, tlppt0idx=%i, diff=%i, npers=%i, ext=%i, maxoutidx=%i\n",tnvpers,tlppt0idx,tnvpers-tlppt0idx,stats->npers,stats->extinction,maxedoutmintimeindex);
    //}
    //assert(tnvpers-tlppt0idx>=stats->npers || stats->extinction==true || maxedoutmintimeindex==tnvpers);

    if(tnvpers-tlppt0idx>stats->npers) {
      tnvpers=tlppt0idx+stats->npers;

      if(stats->inf_timeline[tnvpers]) stats->extinction=false;
    }
    stats->tlpptnvpers=(uint32_t)ceil((tnvpers-tlppt0idx)*0.5)+stats->tlppnnpers;
    //if(stats->extinction && stats->tlpptnvpers-stats->tlppnnpers>0.5*stats->npers) __ro_debug=1;
    //if(stats->extinction && tnvpers-tlppt0idx>stats->npers && stats->inf_timeline[tlppt0idx+stats->npers]) __ro_debug=1;
    //else __ro_debug=0;

    if(stats->maxedoutmintimeindex<INT32_MAX) stats->maxedoutmintimeindex=(int32_t)floor((stats->maxedoutmintimeindex-tlppt0idx)*0.5);
    stats->extinction_time-=stats->first_pos_test_results_time;
    DEBUG_PRINTF("Before merging: %i negatives and %u total. First positive test found at %i\n",stats->tlshift,tnvpers,tlppt0idx);
    DEBUG_PRINTF("First positive test found at %i => %i negative periods and %u total valid periods. maxedoutmintimeindex set to %i\n",tlppt0idx,stats->tlppnnpers,stats->tlpptnvpers,stats->maxedoutmintimeindex);

    stats->pp_inf_timeline=stats->inf_timeline+tlppt0idx;
    stats->pp_newinf_timeline=stats->newinf_timeline+tlppt0idx;
    stats->pp_newpostest_timeline=stats->newpostest_timeline+tlppt0idx;
#ifdef SEC_INF_TIMELINES
    stats->pp_secinf_timeline=stats->secinf_timeline+tlppt0idx;
    stats->pp_newsecinf_timeline=stats->newsecinf_timeline+tlppt0idx;
    stats->pp_newsecpostest_timeline=stats->newsecpostest_timeline+tlppt0idx;
#endif
    stats->pp_ext_timeline=stats->ext_timeline+tlppt0idx;
    //Check if the last index is even. If it is, there is a single fine bin in
    //the last merged bin
    j=(tnvpers-stats->tlshift-tlppt0idx)/2;
    single=(j!=stats->tlpptnvpers-stats->tlppnnpers);;

    for(k=0; k<j; ++k) {
      i=2*k;
      DEBUG_PRINTF("[%i] = [%i] (%u) + [%i] (%u)\n",k,i+tlppt0idx,stats->pp_inf_timeline[i],i+1+tlppt0idx,stats->pp_inf_timeline[i+1]);
      stats->pp_inf_timeline[k]=(stats->pp_inf_timeline[i]>stats->pp_inf_timeline[i+1]?stats->pp_inf_timeline[i]:stats->pp_inf_timeline[i+1]);
      stats->pp_newinf_timeline[k]=stats->pp_newinf_timeline[i]+stats->pp_newinf_timeline[i+1];
      stats->pp_newpostest_timeline[k]=stats->pp_newpostest_timeline[i]+stats->pp_newpostest_timeline[i+1];
#ifdef SEC_INF_TIMELINES
      stats->pp_secinf_timeline[k]=(stats->pp_secinf_timeline[i]>stats->pp_secinf_timeline[i+1]?stats->pp_secinf_timeline[i]:stats->pp_secinf_timeline[i+1]);
      stats->pp_newsecinf_timeline[k]=stats->pp_newsecinf_timeline[i]+stats->pp_newsecinf_timeline[i+1];
      stats->pp_newsecpostest_timeline[k]=stats->pp_newsecpostest_timeline[i]+stats->pp_newsecpostest_timeline[i+1];
#endif

      stats->pp_ext_timeline[k].n=stats->pp_ext_timeline[i].n+stats->pp_ext_timeline[i+1].n;
      stats->pp_ext_timeline[k].rsum=stats->pp_ext_timeline[i].rsum+stats->pp_ext_timeline[i+1].rsum;
      stats->pp_ext_timeline[k].r2sum=stats->pp_ext_timeline[i].r2sum+stats->pp_ext_timeline[i+1].r2sum;
      #ifdef OBSREFF_OUTPUT
      stats->pp_ext_timeline[k].nobs=stats->pp_ext_timeline[i].nobs+stats->pp_ext_timeline[i+1].nobs;
      stats->pp_ext_timeline[k].robssum=stats->pp_ext_timeline[i].robssum+stats->pp_ext_timeline[i+1].robssum;
      stats->pp_ext_timeline[k].robs2sum=stats->pp_ext_timeline[i].robs2sum+stats->pp_ext_timeline[i+1].robs2sum;
      #endif
    }

    if(single) {
      i=-stats->tlshift+tnvpers-1;
      DEBUG_PRINTF("[%i] = [%i] (%u)\n",j,i,stats->inf_timeline[i]);
      stats->pp_inf_timeline[j]=stats->inf_timeline[i];
      stats->pp_newinf_timeline[j]=stats->newinf_timeline[i];
      stats->pp_newpostest_timeline[j]=stats->newpostest_timeline[i];
#ifdef SEC_INF_TIMELINES
      stats->pp_secinf_timeline[j]=stats->secinf_timeline[i];
      stats->pp_newsecinf_timeline[j]=stats->newsecinf_timeline[i];
      stats->pp_newsecpostest_timeline[j]=stats->newsecpostest_timeline[i];
#endif

      stats->pp_ext_timeline[j].n=stats->ext_timeline[i].n;
      stats->pp_ext_timeline[j].rsum=stats->ext_timeline[i].rsum;
      stats->pp_ext_timeline[j].r2sum=stats->ext_timeline[i].r2sum;
      #ifdef OBSREFF_OUTPUT
      stats->pp_ext_timeline[j].nobs=stats->ext_timeline[i].nobs;
      stats->pp_ext_timeline[j].robssum=stats->ext_timeline[i].robssum;
      stats->pp_ext_timeline[j].robs2sum=stats->ext_timeline[i].robs2sum;
      #endif
    } 

    //Check if the first (negative) index is odd. If it is, there is a single fine bin in
    //the first merged bin
    single=tlppt0idx%2;
    j=-stats->tlppnnpers+single;

    for(k=-1; k>=j; --k) {
      i=2*k;
      DEBUG_PRINTF("[%i] = [%i] (%u) + [%i] (%u)\n",k,i+tlppt0idx,stats->pp_inf_timeline[i],i+1+tlppt0idx,stats->pp_inf_timeline[i+1]);
      stats->pp_inf_timeline[k]=(stats->pp_inf_timeline[i]>stats->pp_inf_timeline[i+1]?stats->pp_inf_timeline[i]:stats->pp_inf_timeline[i+1]);
      stats->pp_newinf_timeline[k]=stats->pp_newinf_timeline[i]+stats->pp_newinf_timeline[i+1];
      #ifdef SEC_INF_TIMELINES
      stats->pp_secinf_timeline[k]=(stats->pp_secinf_timeline[i]>stats->pp_secinf_timeline[i+1]?stats->pp_secinf_timeline[i]:stats->pp_secinf_timeline[i+1]);
      stats->pp_newsecinf_timeline[k]=stats->pp_newsecinf_timeline[i]+stats->pp_newsecinf_timeline[i+1];
      #endif

      stats->pp_ext_timeline[k].n=stats->pp_ext_timeline[i].n+stats->pp_ext_timeline[i+1].n;
      stats->pp_ext_timeline[k].rsum=stats->pp_ext_timeline[i].rsum+stats->pp_ext_timeline[i+1].rsum;
      stats->pp_ext_timeline[k].r2sum=stats->pp_ext_timeline[i].rsum+stats->pp_ext_timeline[i+1].r2sum;
      #ifdef OBSREFF_OUTPUT
      stats->pp_ext_timeline[k].nobs=stats->pp_ext_timeline[i].nobs+stats->pp_ext_timeline[i+1].nobs;
      stats->pp_ext_timeline[k].robssum=stats->pp_ext_timeline[i].robssum+stats->pp_ext_timeline[i+1].robssum;
      stats->pp_ext_timeline[k].robs2sum=stats->pp_ext_timeline[i].robs2sum+stats->pp_ext_timeline[i+1].robs2sum;
      #endif
    }

    if(single) {
      i=-stats->tlshift;
      DEBUG_PRINTF("[%i] = [%i] (%u)\n",j-1,i,stats->inf_timeline[i]);
      stats->pp_inf_timeline[j-1]=stats->inf_timeline[i];
      stats->pp_newinf_timeline[j-1]=stats->newinf_timeline[i];
      #ifdef SEC_INF_TIMELINES
      stats->pp_secinf_timeline[j-1]=stats->secinf_timeline[i];
      stats->pp_newsecinf_timeline[j-1]=stats->newsecinf_timeline[i];
      #endif

      stats->pp_ext_timeline[j-1].n=stats->ext_timeline[i].n;
      stats->pp_ext_timeline[j-1].rsum=stats->ext_timeline[i].rsum;
      stats->pp_ext_timeline[j-1].r2sum=stats->ext_timeline[i].r2sum;
      #ifdef OBSREFF_OUTPUT
      stats->pp_ext_timeline[j-1].nobs=stats->ext_timeline[i].nobs;
      stats->pp_ext_timeline[j-1].robssum=stats->ext_timeline[i].robssum;
      stats->pp_ext_timeline[j-1].robs2sum=stats->ext_timeline[i].robs2sum;
      #endif
    } 

    memset(stats->pp_newpostest_timeline-stats->tlppnnpers,0,stats->tlppnnpers*sizeof(uint32_t));
    #ifdef SEC_INF_TIMELINES
    memset(stats->pp_newsecpostest_timeline-stats->tlppnnpers,0,stats->tlppnnpers*sizeof(uint32_t));
    #endif

  } else {
    stats->tlppnnpers=stats->tlshift;
    stats->tlpptnvpers=tnvpers;
    stats->pp_inf_timeline=stats->inf_timeline;
    stats->pp_newinf_timeline=stats->newinf_timeline;
    stats->pp_newpostest_timeline=stats->newpostest_timeline;
    #ifdef SEC_INF_TIMELINES
    stats->pp_secinf_timeline=stats->secinf_timeline;
    stats->pp_newsecinf_timeline=stats->newsecinf_timeline;
    stats->pp_newsecpostest_timeline=stats->newsecpostest_timeline;
    #endif
    stats->pp_ext_timeline=stats->ext_timeline;

    if(sv->pars.pathtype==ro_all_paths) includepath=true;

    else if(sv->pars.pathtype==ro_observable_paths_only) {
      uint32_t* const ptt=stats->postest_timeline-stats->tlshift;
      includepath=false;

      for(i=tnvpers-1; i>=0; --i) if(ptt[i]) {
	includepath=true;
	break;
      }

    } else { //Non-observable paths only. Current implementation is not efficient at all
      uint32_t* const ptt=stats->postest_timeline-stats->tlshift;
      includepath=true;

      for(i=tnvpers-1; i>=0; --i) if(ptt[i]) {
	includepath=false;
	break;
      }
    }

  }

  if(stats->ninfbins) {

    for(i=tnvpers-2; i>=0; --i) {
      //printf("et[%i].n (%u) += %u\n",i,et[i].n,et[i+1].n);
      //et[i].n+=et[i+1].n;
      //et[i].rsum+=et[i+1].rsum;
      et[i].commpersum+=et[i+1].commpersum;
#ifdef NUMEVENTSSTATS
      et[i].neventssum+=et[i+1].neventssum;
#endif
#ifdef OBSREFF_OUTPUT
      //et[i].nobs+=et[i+1].nobs;
      //et[i].robssum+=et[i+1].robssum;
#endif
      //printf("Address et[%i].ngeninfs=%p\n",i-stats->tlshift,et[i].ngeninfs);

      for(j=stats->ninfbins-1; j>=0; --j) {
	//printf("et[%i].ngeninfs[%i] (%lu) += %lu\n",i-stats->tlshift,j,et[i].ngeninfs[j],et[i+1].ngeninfs[j]);
	et[i].ngeninfs[j]+=et[i+1].ngeninfs[j];
	//if(i==0) printf("et[%i].ngeninfs[%i] = %lu\n",i-stats->tlshift,j,et[i].ngeninfs[j]);
      }
    }

  } else {

    for(i=tnvpers-2; i>=0; --i) {
      //printf("et[%i].n (%u) += %u\n",i,et[i].n,et[i+1].n);
      //et[i].n+=et[i+1].n;
      //et[i].rsum+=et[i+1].rsum;
      et[i].commpersum+=et[i+1].commpersum;
#ifdef NUMEVENTSSTATS
      et[i].neventssum+=et[i+1].neventssum;
#endif
#ifdef OBSREFF_OUTPUT
      //et[i].nobs+=et[i+1].nobs;
      //et[i].robssum+=et[i+1].robssum;
#endif
    }
  }
  return includepath;
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
  if(sv->curii->event_time < ((std_summary_stats*)sv->dataptr)->abs_tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
    DEBUG_PRINTF("Pri inf at %i\n",(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time));
    ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
    ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfectionsp;
#endif
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

      #ifdef SEC_INF_TIMELINES
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->secinf_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->secinf_timeline-stats->tlshifta);
      stats->secinf_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->newsecinf_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->newsecinf_timeline-stats->tlshifta);
      stats->newsecinf_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->secpostest_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->secpostest_timeline-stats->tlshifta);
      stats->secpostest_timeline=newarray+newshift;

      newarray=(uint32_t*)malloc(newsize*sizeof(uint32_t));
      memset(newarray,0,dshift*sizeof(uint32_t));
      memcpy(newarray+dshift,stats->newsecpostest_timeline-stats->tlshifta,stats->tnpersa*sizeof(uint32_t));
      free(stats->newsecpostest_timeline-stats->tlshifta);
      stats->newsecpostest_timeline=newarray+newshift;
      #endif

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
  if(sv->curii->event_time < ((std_summary_stats*)sv->dataptr)->abs_tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
    DEBUG_PRINTF("Pri inf at %i\n",(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time));
    ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
    ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfectionsp;
#endif
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
  ii->dataptr=malloc(sizeof(std_stats_inf_data));
}

#ifdef CT_OUTPUT
inline static void std_stats_add_ct_entry(std_summary_stats* stats, const double postesttime, const double presymtime, const int32_t id, const int32_t pid, const uint32_t ntracedcts){
  ++(stats->nctentries);
  DEBUG_PRINTF("%s: %u %22.15e %22.15e %i %i %u\n",__func__, stats->nctentries, postesttime, presymtime, id, pid, ntracedcts);

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
 * @brief Compute the number of observed child infections for events occurring
 * after a time cut.
 *
 * Adds the number of observed infections for the current event to the number of observed infections from
 * the current infectious individual, for events that return false.
 *
 * @param sv: Pointer to the simulation variables.
 */
#ifdef OBSREFF_OUTPUT
inline static void std_stats_calc_obs_child_inf_after_time_cut(sim_vars* sv)
{
    if(sv->curii->commpertype&ro_commper_true_positive_test) {
      infindividual* ii=sv->curii;
      std_summary_stats* stats=(std_summary_stats*)sv->dataptr;

      for(ii->curinfectioni=0; ii->curinfectioni<ii->ninfections; ++(ii->curinfectioni)) {
#ifdef CT_OUTPUT
	ii->gen_ct_time_periods_func(sv, &stats->iibuf, ii, ii->event_time);
#else
	sv->gen_time_periods_func(sv, &stats->iibuf, ii, ii->event_time);
#endif
	((std_stats_inf_data*)sv->curii->dataptr)->nobsinf+=((stats->iibuf.commpertype&ro_commper_int_true_positive_test)==ro_commper_int_true_positive_test);
	DEBUG_PRINTF("Adding %u to the number of observed infections, for a total of %u\n",((stats->iibuf.commpertype&ro_commper_true_positive_test)!=0),((std_stats_inf_data*)sv->curii->dataptr)->nobsinf);
      }
    }
}
#endif

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before abs_tmax.  This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before abs_tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts);
#endif

  if(sv->curii->ninfections) {
    ((std_stats_inf_data*)sv->curii->dataptr)->ninf+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((std_stats_inf_data*)sv->curii->dataptr)->ninf);

    if(sv->curii->event_time < ((std_summary_stats*)sv->dataptr)->abs_tmax && sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {
      ((std_summary_stats*)sv->dataptr)->newinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
      ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[(int)floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time)]+=sv->curii->ninfectionsp;
#endif
      return true;
    }
#ifdef OBSREFF_OUTPUT
    std_stats_calc_obs_child_inf_after_time_cut(sv);
#endif
  }
  return false;
}

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before abs_tmax. If the incremented bin for the total infection timeline
 * exceeds nimax, then set the extinction for the current path to false and
 * update the value for the minimum time index where the maximum number of
 * infectious individuals was exceeded if required. This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before abs_tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event_nimax(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts);
#endif

  if(sv->curii->ninfections) {
    ((std_stats_inf_data*)sv->curii->dataptr)->ninf+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((std_stats_inf_data*)sv->curii->dataptr)->ninf);

    if(sv->curii->event_time < ((std_summary_stats*)sv->dataptr)->abs_tmax) {
      const int eti=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time);

      if(sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {

	if(((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+sv->curii->ninfections < ((std_summary_stats*)sv->dataptr)->nimax) {
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
	  ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[eti]+=sv->curii->ninfectionsp;
#endif

	} else if(((std_summary_stats*)sv->dataptr)->newinf_timeline[eti] < ((std_summary_stats*)sv->dataptr)->nimax) {
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
	  ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[eti]+=sv->curii->ninfectionsp;
#endif
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->newinf_timeline[eti],((std_summary_stats*)sv->dataptr)->nimax);

	} else {
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  DEBUG_PRINTF("nimax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->newinf_timeline[eti],((std_summary_stats*)sv->dataptr)->nimax);
	  goto nimax_event_false;
	}
	return true;
      }
    }
nimax_event_false:
    ;
#ifdef OBSREFF_OUTPUT
    std_stats_calc_obs_child_inf_after_time_cut(sv);
#endif
  }
  return false;
}

/**
 * @brief Processes the number of infections for this new event.
 *
 * Adds the number of infections for this new event to the number of infections from
 * the current infectious individual, and to the total infection timeline if the event
 * time is before abs_tmax. If the bin for the recent positive test result timeline
 * exceeds npostestmax, then set the extinction for the current path to false and
 * update the value for the minimum time index where the maximum number of
 * recent positive test results was exceeded if required. This function must be assigned to
 * the simulation engine through a call of sim_set_new_event_proc_func.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event time is before abs_tmax and new infections were
 * generated, and false otherwise.
 **/
inline static bool std_stats_new_event_npostestmax(sim_vars* sv)
{
#ifdef CT_OUTPUT
  ((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts+=sv->curii->ntracednicts+sv->curii->ntracedicts;
  DEBUG_PRINTF("Successfully traced contacts incremented to %u\n",((std_stats_inf_data*)sv->curii->dataptr)->ntracedcts);
#endif

  if(sv->curii->ninfections) {
    ((std_stats_inf_data*)sv->curii->dataptr)->ninf+=sv->curii->ninfections;
    DEBUG_PRINTF("%s: Number of infections incremented to %u\n",__func__,((std_stats_inf_data*)sv->curii->dataptr)->ninf);

    if(sv->curii->event_time < ((std_summary_stats*)sv->dataptr)->abs_tmax) {
      const int eti=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*sv->curii->event_time);

      if(sv->curii <= sv->brsim.iis+((std_summary_stats*)sv->dataptr)->lmax) {

	if(((std_summary_stats*)sv->dataptr)->postest_timeline[eti] < ((std_summary_stats*)sv->dataptr)->npostestmax) {
	  ((std_summary_stats*)sv->dataptr)->newinf_timeline[eti]+=sv->curii->ninfections;
#ifdef SEC_INF_TIMELINES
	  ((std_summary_stats*)sv->dataptr)->newsecinf_timeline[eti]+=sv->curii->ninfectionsp;
#endif

	} else {
	  ((std_summary_stats*)sv->dataptr)->extinction=false;

	  if(eti < ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex)  {
	    DEBUG_PRINTF("maxedoutmintimeindex reduced from %i to %i\n",((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex,eti);
	    ((std_summary_stats*)sv->dataptr)->maxedoutmintimeindex=eti;
	  }
	  DEBUG_PRINTF("npostestmax exceeded for time index %i (%u vs %u)\n",eti,((std_summary_stats*)sv->dataptr)->postest_timeline[eti],((std_summary_stats*)sv->dataptr)->npostestmax);
	  goto npostestmax_event_false;
	}
	return true;
      }
    }
npostestmax_event_false:
    ;
#ifdef OBSREFF_OUTPUT
    std_stats_calc_obs_child_inf_after_time_cut(sv);
#endif
  }
  return false;
}

inline static void std_stats_fill_newpostest(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  if(ii->commpertype&ro_commper_true_positive_test) {
    const int trt=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*(ii->end_comm_period+sv->pars.tdeltat));

    DEBUG_PRINTF("New pos test at %i\n",trt);

    if(trt<((std_summary_stats*)sv->dataptr)->abs_maxnpers) {
      ++(((std_summary_stats*)sv->dataptr)->newpostest_timeline[trt]);

      #ifdef SEC_INF_TIMELINES
      if(ii->inftypep) ++(((std_summary_stats*)sv->dataptr)->newsecpostest_timeline[trt]);
      #endif
    }

#ifdef CT_OUTPUT
    ((std_stats_inf_data*)ii->dataptr)->id=++(((std_summary_stats*)sv->dataptr)->curctid);
    DEBUG_PRINTF("True positive test with ID %u\n",((std_stats_inf_data*)ii->dataptr)->id);
    ((std_stats_inf_data*)ii->dataptr)->ntracedcts=0;
#endif
#ifdef OBSREFF_OUTPUT
    if(ii->commpertype&ro_commper_int) ++(((std_stats_inf_data*)parent->dataptr)->nobsinf);
    DEBUG_PRINTF("Number of observed infections for the parent increased by 1, for a total of %u\n",((std_stats_inf_data*)parent->dataptr)->nobsinf);
#endif
    int32_t i;
    const int32_t end_npostestmax_per_i=trt+((std_summary_stats*)sv->dataptr)->nbinsperunit*((std_summary_stats*)sv->dataptr)->npostestmaxnunits-1;

    #ifdef SEC_INF_TIMELINES
    if(ii->inftypep) for(i=(end_npostestmax_per_i>=((std_summary_stats*)sv->dataptr)->abs_maxnpers?((std_summary_stats*)sv->dataptr)->abs_maxnpers-1:end_npostestmax_per_i); i>=trt; --i) {
      ++(((std_summary_stats*)sv->dataptr)->postest_timeline[i]);
      ++(((std_summary_stats*)sv->dataptr)->secpostest_timeline[i]);
    }

    else for(i=(end_npostestmax_per_i>=((std_summary_stats*)sv->dataptr)->abs_maxnpers?((std_summary_stats*)sv->dataptr)->abs_maxnpers-1:end_npostestmax_per_i); i>=trt; --i) ++(((std_summary_stats*)sv->dataptr)->postest_timeline[i]);
    #else
    for(i=(end_npostestmax_per_i>=((std_summary_stats*)sv->dataptr)->abs_maxnpers?((std_summary_stats*)sv->dataptr)->abs_maxnpers-1:end_npostestmax_per_i); i>=trt; --i) ++(((std_summary_stats*)sv->dataptr)->postest_timeline[i]);
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
  ((std_stats_inf_data*)ii->dataptr)->ninf=0;
  DEBUG_PRINTF("%s: Number of infections initialized to 0.\n",__func__);
#ifdef OBSREFF_OUTPUT
  ((std_stats_inf_data*)ii->dataptr)->nobsinf=0;
  DEBUG_PRINTF("%s: Number of observed infections initialized to 0.\n",__func__);
#endif

  //++(*(uint32_t*)(ii-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(ii-1)->dataptr);
  DEBUG_PRINTF("%s\n",__func__);
  std_stats_fill_newpostest(sv, ii, parent);
}

/**
 * @brief Update the parameters and timeline memory allocation for a simulation
 * using a time relative to the first positive test results.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Infectious individual.
 */
inline static void first_pos_test_results_update(sim_vars* sv, infindividual* ii)
{
  std_summary_stats* const stats=(std_summary_stats*)sv->dataptr;

  int32_t newsize;

  if(ii->commpertype&ro_commper_true_positive_test && ii->end_comm_period+sv->pars.tdeltat < stats->first_pos_test_results_time) {
    std_summary_stats* const stats=(std_summary_stats*)sv->dataptr;
    stats->first_pos_test_results_time=ii->end_comm_period+sv->pars.tdeltat;

    stats->abs_maxnpers=floor(stats->nbinsperunit*stats->first_pos_test_results_time)+stats->npers;
    stats->abs_tmax=((double)stats->abs_maxnpers)/stats->nbinsperunit;

    newsize=(int32_t)(stats->nbinsperunit*(ii->end_comm_period+sv->pars.tdeltat+stats->npostestmaxnunits))+1;

    if(newsize > stats->abs_npers) stats->abs_npers=(newsize<=stats->abs_maxnpers?newsize:stats->abs_maxnpers);

    newsize=stats->abs_maxnpers;

  } else {
    newsize=(int32_t)(stats->nbinsperunit*ii->end_comm_period)+1;

    if(newsize > stats->abs_npers) stats->abs_npers=(newsize<=stats->abs_maxnpers?newsize:stats->abs_maxnpers);
  }

  if(newsize > stats->tnpersa) {
    const uint32_t dsize=newsize-stats->tnpersa;

    stats->inf_timeline=(uint32_t*)realloc(stats->inf_timeline,newsize*sizeof(uint32_t));
    memset(stats->inf_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->newinf_timeline=(uint32_t*)realloc(stats->newinf_timeline,newsize*sizeof(uint32_t));
    memset(stats->newinf_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->postest_timeline=(uint32_t*)realloc(stats->postest_timeline,newsize*sizeof(uint32_t));
    memset(stats->postest_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->newpostest_timeline=(uint32_t*)realloc(stats->newpostest_timeline,newsize*sizeof(uint32_t));
    memset(stats->newpostest_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    #ifdef SEC_INF_TIMELINES
    stats->secinf_timeline=(uint32_t*)realloc(stats->secinf_timeline,newsize*sizeof(uint32_t));
    memset(stats->secinf_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->newsecinf_timeline=(uint32_t*)realloc(stats->newsecinf_timeline,newsize*sizeof(uint32_t));
    memset(stats->newsecinf_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->secpostest_timeline=(uint32_t*)realloc(stats->secpostest_timeline,newsize*sizeof(uint32_t));
    memset(stats->secpostest_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));

    stats->newsecpostest_timeline=(uint32_t*)realloc(stats->newsecpostest_timeline,newsize*sizeof(uint32_t));
    memset(stats->newsecpostest_timeline+stats->tnpersa,0,dsize*sizeof(uint32_t));
    #endif

    stats->ext_timeline=(ext_timeline_info*)realloc(stats->ext_timeline,newsize*sizeof(ext_timeline_info));
    memset(stats->ext_timeline+stats->tnpersa,0,dsize*sizeof(ext_timeline_info));

    if(stats->nainfbins) {
      uint64_t* newarray64;
      int32_t i;

      for(i=dsize-1; i>=0; --i) {
	newarray64=(uint64_t*)malloc(stats->nainfbins*sizeof(uint64_t));
	memset(newarray64,0,stats->nainfbins*sizeof(uint64_t));
	stats->ext_timeline[stats->tnpersa+i].ngeninfs=newarray64;
      }
    }
    stats->tnpersa=newsize;
  }
}

/**
 * @brief Initialises the processing variable for a new infectious individual
 * for a simulation where the time is relative to the time of the first positive
 * test results.
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
inline static void std_stats_new_inf_first_pos_test_results(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  first_pos_test_results_update(sv, ii);
  std_stats_new_inf(sv, ii, parent);
}

inline static void std_stats_fill_inf_ext_n(sim_vars* sv, infindividual* ii)
{
  const double start_comm_per=ii->end_comm_period-ii->comm_period;
  const int32_t start_latent_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*(ii->end_comm_period-(ii->comm_period+ii->latent_period)));
  const int32_t end_comm_per_i=(((std_summary_stats*)sv->dataptr)->nbinsperunit*ii->end_comm_period >= ((std_summary_stats*)sv->dataptr)->abs_maxnpers ? ((std_summary_stats*)sv->dataptr)->abs_maxnpers-1 : floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*ii->end_comm_period));
  int32_t i;

  if(start_comm_per < ((std_summary_stats*)sv->dataptr)->abs_tmax) {
    const int32_t start_comm_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*start_comm_per);
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].rsum+=((std_stats_inf_data*)ii->dataptr)->ninf;
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].r2sum+=((std_stats_inf_data*)ii->dataptr)->ninf*((std_stats_inf_data*)ii->dataptr)->ninf;
    ++(((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].n);
#ifdef OBSREFF_OUTPUT
    if(ii->commpertype&ro_commper_true_positive_test) {
      ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].robssum+=((std_stats_inf_data*)ii->dataptr)->nobsinf;
      ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].robs2sum+=((std_stats_inf_data*)ii->dataptr)->nobsinf*=((std_stats_inf_data*)ii->dataptr)->nobsinf;
      ++(((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].nobs);
    }
#endif
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].commpersum+=ii->comm_period;
#ifdef NUMEVENTSSTATS
    ((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].neventssum+=ii->nevents;
#endif
  }

  #ifdef SEC_INF_TIMELINES
  if(ii->inftypep) for(i=start_latent_per_i; i<=end_comm_per_i; ++i) {
    ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
    ++(((std_summary_stats*)sv->dataptr)->secinf_timeline[i]);
  }

  else for(i=start_latent_per_i; i<=end_comm_per_i; ++i) ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
  #else
  for(i=start_latent_per_i; i<=end_comm_per_i; ++i) ++(((std_summary_stats*)sv->dataptr)->inf_timeline[i]);
  #endif
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
 * at abs_tmax, then the path is set to not go extinct. Otherwise, the
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
  DEBUG_PRINTF("Number of infections was %u\n",((std_stats_inf_data*)ii->dataptr)->ninf);

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY),  ((std_stats_inf_data*)ii->dataptr)->id, ((ii->commpertype&ro_commper_int)? ((std_stats_inf_data*)parent->dataptr)->id:-((std_stats_inf_data*)parent->dataptr)->id), ((std_stats_inf_data*)ii->dataptr)->ntracedcts);
#endif

  //If truncated by abs_tmax
  if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->abs_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

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

  if(start_comm_per < ((std_summary_stats*)sv->dataptr)->abs_tmax) {
    std_summary_stats* const stats=(std_summary_stats*)sv->dataptr;

    if(((std_stats_inf_data*)ii->dataptr)->ninf >= stats->ninfbins) {
      stats->ninfbins=((std_stats_inf_data*)ii->dataptr)->ninf+1;

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
    ++(stats->ext_timeline[start_comm_per_i].ngeninfs[((std_stats_inf_data*)ii->dataptr)->ninf]);
  }

  std_stats_end_inf(sv, ii, parent);
}

/**
 * @brief Process the statistics for an infectious individual that does not
 * participate to any transmission event.
 *
 * This function adds the individual's communicable period to the sum
 * of communicable periods. If the individual was still infectious
 * at abs_tmax, then the path is set to not go extinct. Otherwise, the
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

  ((std_stats_inf_data*)ii->dataptr)->ninf=0;
#ifdef OBSREFF_OUTPUT
  ((std_stats_inf_data*)ii->dataptr)->nobsinf=0;
#endif

  std_stats_fill_newpostest(sv, ii, parent);

#ifdef CT_OUTPUT
  if(ii->commpertype&ro_commper_true_positive_test) std_stats_add_ct_entry((std_summary_stats*)sv->dataptr, ii->end_comm_period+sv->pars.tdeltat, ((ii->commpertype&ro_commper_alt)?ii->end_comm_period-ii->comm_period+ii->presym_comm_period:INFINITY), ((std_stats_inf_data*)ii->dataptr)->id, ((ii->commpertype&ro_commper_int)?((std_stats_inf_data*)parent->dataptr)->id:-((std_stats_inf_data*)parent->dataptr)->id), ((std_stats_inf_data*)ii->dataptr)->ntracedcts);
#endif

  //If truncated by tmax
  if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->abs_tmax) ((std_summary_stats*)sv->dataptr)->extinction=false;

  else {
    //++((std_summary_stats*)sv->dataptr)->n_ended_infections;

    if(ii->end_comm_period > ((std_summary_stats*)sv->dataptr)->extinction_time) ((std_summary_stats*)sv->dataptr)->extinction_time=ii->end_comm_period;
  }

  std_stats_fill_inf_ext_n(sv, ii);
}

/**
 * @brief Process the statistics for an infectious individual that does not
 * participate to any transmission event, for a simulation where the time is
 * relative to the time of the first positive test results.
 *
 * This function adds the individual's communicable period to the sum
 * of communicable periods. If the individual was still infectious
 * at abs_tmax, then the path is set to not go extinct. Otherwise, the
 * extinction time for the generating path is updated if its current value was
 * exceeded, using the current event time for the individual's parent and the
 * individual's communicable period. The timeline for the number of infectious
 * individuals is also updated.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Pointer to the infectious individual.
 * @param parent: Pointer to the infectious individual's parent.
 **/
inline static void std_stats_noevent_new_inf_first_pos_test_results(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  first_pos_test_results_update(sv, ii);
  std_stats_noevent_new_inf(sv, ii, parent);
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

  if(start_comm_per < ((std_summary_stats*)sv->dataptr)->abs_tmax) {
    const int32_t start_comm_per_i=floor(((std_summary_stats*)sv->dataptr)->nbinsperunit*start_comm_per);
    ++(((std_summary_stats*)sv->dataptr)->ext_timeline[start_comm_per_i].ngeninfs[0]);
  }

  std_stats_noevent_new_inf(sv, ii, parent);
}

/**
 * @brief Process the statistics for an infectious individual that does not
 * participate to any transmission event and record the number of infections
 * generated by each infectious individual, for a simulation where the time is
 * relative to the time of the first positive test results.
 *
 * In addition to calling the function std_stats_noevent_new_inf, this function records
 * the number of infections generated by each infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii: Pointer to the infectious individual.
 * @param parent: Pointer to the infectious individual's parent.
 **/
inline static void std_stats_noevent_new_inf_rec_ninfs_first_pos_test_results(sim_vars* sv, infindividual* ii, infindividual* parent)
{
  first_pos_test_results_update(sv, ii);
  std_stats_noevent_new_inf_rec_ninfs(sv, ii, parent);
}

#endif
