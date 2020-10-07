/**
 * @file standard_summary_stats.c
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "standard_summary_stats.h"

void std_stats_init(sim_vars* sv, uint64_t** ngeninfs, uint32_t* ninfbins)
{
  std_summary_stats* stats=(std_summary_stats*)sv->dataptr;

  stats->timelineshift=0;
  stats->tnpersa=stats->npers=(int)sv->pars.tmax+1;
  //stats->inf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newinf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->postest_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newpostest_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->ext_timeline=(ext_timeline_info*)malloc(stats->tnpersa*sizeof(ext_timeline_info));

  int32_t i;

  for(i=stats->tnpersa-1; i>=0; --i) {
    stats->ext_timeline[i].inf=((uint32_t*)malloc((stats->tnpersa-i)*sizeof(uint32_t)))-i;
    //printf("%i: Allocating at %p->%p over %lu\n",i,stats->ext_timeline+i,stats->ext_timeline[i].inf+i,(stats->tnpersa-i)*sizeof(uint32_t));
  }

  if(ngeninfs && ninfbins) {
    *ngeninfs=(uint64_t*)malloc(INIT_NINF_ALLOC*sizeof(uint64_t));
    stats->ngeninfs=ngeninfs;
    *ninfbins=INIT_NINF_ALLOC;
    stats->ninfbins=ninfbins;
    memset(*stats->ngeninfs,0,*stats->ninfbins*sizeof(uint64_t));
  }
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
  //free(stats->inf_timeline-stats->timelineshift);
  free(stats->newinf_timeline-stats->timelineshift);
  free(stats->postest_timeline-stats->timelineshift);
  free(stats->newpostest_timeline-stats->timelineshift);

  ext_timeline_info* const set=stats->ext_timeline-stats->timelineshift;
  int32_t i;

  for(i=stats->tnpersa-1; i>=0; --i) {
    //printf("%i: Free at %p\n",i,set[i].inf+i-stats->timelineshift);
    free(set[i].inf+i-stats->timelineshift);
  }
  free(stats->ext_timeline-stats->timelineshift);

#ifdef CT_OUTPUT

  for(i=stats->nactentries-1; i>=0; --i) free(stats->ctentries[i]);
  free(stats->ctentries);
#endif
}
