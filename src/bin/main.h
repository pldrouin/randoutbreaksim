/**
 * @file main.h
 * @brief Main function for the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include <endian.h>

#include <assert.h>

#include <pthread.h>

#include <gsl/gsl_rng.h>

#include "rngstream_gsl.h"

#include "config.h"
#include "model_parameters.h"
#include "infindividual.h"
#include "branchsim.h"
#include "standard_summary_stats.h"

/**
 * @brief Main function.
 */
int main(const int nargs, const char* args[]);

typedef struct {
  config_pars const* cp;
  double npathsperset;
  uint32_t nsets;
  uint32_t npers;
  uint32_t tnpersa;
  uint32_t volatile* set;
  double commper_mean;
#ifdef NUMEVENTSSTATS
  double nevents_mean;
#endif
  double r_mean;
  double pe;
  double te_mean;
  double te_std;
  double* inf_timeline_mean_ext;
  double* inf_timeline_std_ext;
  double* inf_timeline_mean_noext;
  double* inf_timeline_std_noext;
  double* newinf_timeline_mean_ext;
  double* newinf_timeline_std_ext;
  double* newinf_timeline_mean_noext;
  double* newinf_timeline_std_noext;
  double* totaltctc_timeline_mean_ext;
  double* totaltctc_timeline_std_ext;
  double* totaltctc_timeline_mean_noext;
  double* totaltctc_timeline_std_noext;
  uint64_t* ngeninfs;
  uint32_t ninfbins;
  uint32_t nimaxedoutmintimeindex;
  gsl_rng* r;
  pthread_mutex_t* tlflock;
} thread_data;

void* simthread(void* arg);

inline static ssize_t tlo_write_reg_path(std_summary_stats const* stats, char* buf)
{
  int32_t b;

  for(b=stats->npers-1; b>0; --b) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->newinf_timeline[b],b,stats->inf_timeline[b]);

    if(stats->inf_timeline[b]) break;
    assert(stats->newinf_timeline[b]==0);
  }

  //If timelines are all zero
  if(b==0 && stats->inf_timeline[b]==0) {
    assert(stats->newinf_timeline[b]==0);
    *(uint32_t*)buf=0;
    buf[4]=1;
    return 5;
  }

  const uint32_t nbins=b+1;
  *(uint32_t*)buf=htole32(nbins);
  buf[4]=(char)stats->extinction;;
  buf+=5;

  for(b=0; b<nbins; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    ((uint32_t*)buf)[2*b]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[2*b+1]=htole32(stats->newinf_timeline[b]);
  }

  return 5+nbins*8*4;
}

inline static ssize_t tlo_write_reltime_path(std_summary_stats const* stats, char* buf)
{
  int32_t bmin=-(int32_t)stats->timelineshift;
  int32_t bmax;

  for(bmax=stats->npers-1; bmax>bmin; --bmax) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",bmax,stats->newinf_timeline[bmax],bmax,stats->inf_timeline[bmax]);

    if(stats->inf_timeline[bmax]) break;
    assert(stats->newinf_timeline[bmax]==0);
  }

  //If timelines are all zero
  if(bmax==-bmin && stats->inf_timeline[bmax]==0) {
    assert(stats->newinf_timeline[bmax]==0);
    *(uint64_t*)buf=0;
    buf[8]=1;
    return 9;
  }


  for(;; ++bmin) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",bmin,stats->newinf_timeline[bmin],bmin,stats->inf_timeline[bmin]);

    if(stats->inf_timeline[bmin]) break;
    assert(stats->newinf_timeline[bmin]==0);
  }
  assert(bmin<=0);

  const uint32_t nbins=bmax-bmin+1;

  ((uint32_t*)buf)[0]=htole32(nbins);
  ((uint32_t*)buf)[1]=htole32(-bmin);
  buf[8]=(char)stats->extinction;
  buf+=9;

  int32_t b;
  uint32_t bp;

  for(b=bmin; b<=bmax; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    bp=2*(b-bmin);
    ((uint32_t*)buf)[bp]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[bp+1]=htole32(stats->newinf_timeline[b]);
  }

  return 9+nbins*8;
}
