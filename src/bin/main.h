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
  uint32_t id;
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
  double* newpostest_timeline_mean_ext;
  double* newpostest_timeline_std_ext;
  double* newpostest_timeline_mean_noext;
  double* newpostest_timeline_std_noext;
  uint64_t* ngeninfs;
  uint32_t ninfbins;
  int32_t maxedoutmintimeindex;
  gsl_rng* r;
  pthread_mutex_t* tlflock;
  pthread_mutex_t* ctflock;
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
  /*
  if(b==0 && stats->inf_timeline[b]==0) {
    assert(stats->newinf_timeline[b]==0);
    *(uint32_t*)buf=0;
    buf[4]=1;
    return 5;
  }
  */
  assert(b>0 || (stats->inf_timeline[b] && stats->newinf_timeline[b]));

  const uint32_t nbins=b+1;
  *(uint32_t*)buf=htole32(nbins);
  *(uint32_t*)(buf+4)=htole32(stats->maxedoutmintimeindex);
  *(uint32_t*)(buf+8)=htole64((int32_t)(stats->extinction?floor(stats->extinction_time):-INT32_MAX));
  buf+=12;

  for(b=0; b<nbins; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[b+nbins]=htole32(stats->newinf_timeline[b]);
  }

  return 12+nbins*8;
}

inline static ssize_t tlo_write_reg_postest_path(std_summary_stats const* stats, char* buf)
{
  int32_t b;

  for(b=stats->npers-1; b>0; --b) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->newinf_timeline[b],b,stats->inf_timeline[b]);

    if(stats->inf_timeline[b] || stats->newpostest_timeline[b]) break;
    assert(stats->newinf_timeline[b]==0);
  }

  //If timelines are all zero
  /*
  if(b==0 && stats->inf_timeline[b]==0) {
    assert(stats->newinf_timeline[b]==0);
    assert(stats->newpostest_timeline[b]==0);
    *(uint32_t*)buf=0;
    buf[4]=1;
    return 5;
  }
  */
  assert(b>0 || (stats->inf_timeline[b] && stats->newinf_timeline[b] && stats->newpostest_timeline[b]));

  const uint32_t nbins=b+1;
  *(uint32_t*)buf=htole32(nbins);
  *(uint32_t*)(buf+4)=htole32(stats->maxedoutmintimeindex);
  *(uint32_t*)(buf+8)=htole64((int32_t)(stats->extinction?floor(stats->extinction_time):-INT32_MAX));
  buf+=12;
  const uint32_t tnbins=2*nbins;

  for(b=0; b<nbins; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[b+nbins]=htole32(stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b+tnbins]=htole32(stats->newpostest_timeline[b]);
  }

  return 12+nbins*12;
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
  /*
  if(bmax==bmin && stats->inf_timeline[bmax]==0) {
    assert(stats->newinf_timeline[bmax]==0);
    *(uint64_t*)buf=0;
    buf[8]=1;
    return 9;
  }
  */
  assert(bmax>bmin || (stats->inf_timeline[bmax] && stats->newinf_timeline[bmax]));

  for(;; ++bmin) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",bmin,stats->newinf_timeline[bmin],bmin,stats->inf_timeline[bmin]);

    if(stats->inf_timeline[bmin]) break;
    assert(stats->newinf_timeline[bmin]==0);
  }
  assert(bmin<=0);

  const uint32_t nbins=bmax-bmin+1;

  ((uint32_t*)buf)[0]=htole32(nbins);
  ((uint32_t*)buf)[1]=htole32(-bmin);
  *(uint32_t*)(buf+8)=htole32(stats->maxedoutmintimeindex);
  *(uint32_t*)(buf+12)=htole64((int32_t)(stats->extinction?floor(stats->extinction_time):-INT32_MAX));
  buf+=16;

  int32_t b;
  const int32_t nbmbm=-bmin+nbins;

  for(b=bmin; b<=bmax; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b-bmin]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[b+nbmbm]=htole32(stats->newinf_timeline[b]);
  }

  return 16+nbins*8;
}

inline static ssize_t tlo_write_reltime_postest_path(std_summary_stats const* stats, char* buf)
{
  int32_t bmin=-(int32_t)stats->timelineshift;
  int32_t bmax;

  for(bmax=stats->npers-1; bmax>bmin; --bmax) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",bmax,stats->newinf_timeline[bmax],bmax,stats->inf_timeline[bmax]);

    if(stats->inf_timeline[bmax] || stats->newpostest_timeline[bmax]) break;
    assert(stats->newinf_timeline[bmax]==0);
  }

  //If timelines are all zero
  /*
  if(bmax==bmin && stats->inf_timeline[bmax]==0) {
    
    assert(stats->newinf_timeline[bmax]==0);
    assert(stats->newpostest_timeline[bmax]==0);
    *(uint64_t*)buf=0;
    buf[8]=1;
    return 9;
  }
  */
  assert(bmax>bmin || (stats->inf_timeline[bmax] && stats->newinf_timeline[bmax] && stats->newpostest_timeline[bmax]));

  for(;; ++bmin) {
    //printf("TotInf[%" PRIi32 "]=%" PRIu32 ",\tInf[%" PRIi32 "]=%" PRIu32 "\n",bmin,stats->newinf_timeline[bmin],bmin,stats->inf_timeline[bmin]);

    if(stats->inf_timeline[bmin]) break;
    assert(stats->newinf_timeline[bmin]==0);
    assert(stats->newpostest_timeline[bmin]==0);
  }
  assert(bmin<=0);

  const uint32_t nbins=bmax-bmin+1;

  ((uint32_t*)buf)[0]=htole32(nbins);
  ((uint32_t*)buf)[1]=htole32(-bmin);
  *(uint32_t*)(buf+8)=htole32(stats->maxedoutmintimeindex);
  *(uint32_t*)(buf+12)=htole64((int32_t)(stats->extinction?floor(stats->extinction_time):-INT32_MAX));
  buf+=16;

  int32_t b;
  const int32_t nbmbm=-bmin+nbins;
  const int32_t tnbmbm=-bmin+2*nbins;

  for(b=bmin; b<=bmax; ++b) {
    //printf("Inf[%" PRIi32 "]=%" PRIu32 ",\tTotInf[%" PRIi32 "]=%" PRIu32 "\n",b,stats->inf_timeline[b],b,stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b-bmin]=htole32(stats->inf_timeline[b]);
    ((uint32_t*)buf)[b+nbmbm]=htole32(stats->newinf_timeline[b]);
    ((uint32_t*)buf)[b+tnbmbm]=htole32(stats->newpostest_timeline[b]);
  }

  return 16+nbins*12;
}

#ifdef CT_OUTPUT
inline static int ctcompar(const void* first, const void* second){return ((*(ctposinf**)first)->postesttime<(*(ctposinf**)second)->postesttime?-1:((*(ctposinf**)first)->postesttime>(*(ctposinf**)second)->postesttime?1:0));}

inline static ssize_t ct_write_func(std_summary_stats const* stats, char* buf)
{
  uint32_t i;
  char* buf0=buf;

  //printf("Extinction: %u, nctentries: %u\n",stats->extinction,stats->nctentries);

  if(stats->maxedoutmintimeindex==INT32_MAX) {

    for(i=0; i<stats->nctentries; ++i) {
      *((uint32_t*)(buf))=htole32(stats->ctentries[i]->postesttime);
      *((uint32_t*)(buf+4))=htole32(stats->ctentries[i]->presymtime);
      *((uint32_t*)(buf+8))=htole32(stats->ctentries[i]->id);
      *((uint32_t*)(buf+12))=htole32(stats->ctentries[i]->pid);
      *((uint32_t*)(buf+16))=htole32(stats->ctentries[i]->ntracedcts);
      buf+=20;
    }

  } else {

    for(i=0; i<stats->nctentries; ++i) {

      if(floor(stats->ctentries[i]->postesttime/1440.) <= stats->maxedoutmintimeindex) {
	*((uint32_t*)(buf))=htole32(stats->ctentries[i]->postesttime);
	*((uint32_t*)(buf+4))=htole32(stats->ctentries[i]->presymtime);
	*((uint32_t*)(buf+8))=htole32(stats->ctentries[i]->id);
	*((uint32_t*)(buf+12))=htole32(stats->ctentries[i]->pid);
	*((uint32_t*)(buf+16))=htole32(stats->ctentries[i]->ntracedcts);
	buf+=20;
      }
    }
  }
  return buf-buf0;
}
#endif

