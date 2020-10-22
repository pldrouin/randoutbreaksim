/**
 * @file standard_summary_stats.c
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "standard_summary_stats.h"

void std_stats_init(sim_vars* sv, const uint32_t nbinsperunit, bool ngeninfs)
{
  std_summary_stats* stats=(std_summary_stats*)sv->dataptr;

  if(sv->pars.timetype==ro_time_first_pos_test_results) {
    stats->nbinsperunit=2*nbinsperunit;
    stats->tnpersa=stats->npers=(uint32_t)(stats->nbinsperunit*sv->pars.tmax);
    stats->abs_maxnpers=INT32_MAX;
    stats->abs_tmax=stats->abs_maxnpers/stats->nbinsperunit;
    stats->first_pos_test_results_time=INFINITY;
    stats->abs_npers=0;

  } else {
    stats->nbinsperunit=nbinsperunit;
    stats->tnpersa=stats->npers=stats->abs_maxnpers=stats->abs_npers=(uint32_t)(stats->nbinsperunit*sv->pars.tmax);
    stats->abs_tmax=sv->pars.tmax;
  }

  stats->tlshift=stats->tlshifta=0;
  stats->inf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newinf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->postest_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newpostest_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->ext_timeline=(ext_timeline_info*)malloc(stats->tnpersa*sizeof(ext_timeline_info));

  int32_t i;

  if(ngeninfs) {
    stats->nainfbins=INIT_NINF_ALLOC;

    for(i=stats->tnpersa-1; i>=0; --i) {
      stats->ext_timeline[i].ngeninfs=(uint64_t*)malloc(INIT_NINF_ALLOC*sizeof(uint64_t));
      //printf("ext_timeline[%i].ngeninfs=%p\n",i,stats->ext_timeline[i].ngeninfs);
      memset(stats->ext_timeline[i].ngeninfs,0,INIT_NINF_ALLOC*sizeof(uint64_t));
    }

  } else stats->nainfbins=stats->ninfbins=0;

  stats->lmax=UINT32_MAX;
  stats->nimax=UINT32_MAX;

#ifdef CT_OUTPUT
  stats->nactentries=INIT_NACTENTRIES;
  stats->ctentries=(ctposinf**)malloc(INIT_NACTENTRIES*sizeof(ctposinf*));

  for(i=stats->nactentries-1; i>=0; --i) stats->ctentries[i]=(ctposinf*)malloc(sizeof(ctposinf));
#endif
}

void std_stats_free(std_summary_stats* stats)
{
  free(stats->inf_timeline-stats->tlshifta);
  free(stats->newinf_timeline-stats->tlshifta);
  free(stats->postest_timeline-stats->tlshifta);
  free(stats->newpostest_timeline-stats->tlshifta);

  ext_timeline_info* const set=stats->ext_timeline-stats->tlshifta;
  int32_t i;

  if(stats->nainfbins) {

    for(i=stats->tnpersa-1; i>=0; --i) {
      //printf("Free ext_timeline[%i] (%p)\n",i-stats->tlshifta,set[i].ngeninfs);
      free(set[i].ngeninfs);
    }
  }
  free(stats->ext_timeline-stats->tlshifta);

#ifdef CT_OUTPUT

  for(i=stats->nactentries-1; i>=0; --i) free(stats->ctentries[i]);
  free(stats->ctentries);
#endif
}
