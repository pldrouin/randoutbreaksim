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
  stats->inf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newinf_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));
  stats->newpostest_timeline=(uint32_t*)malloc(stats->tnpersa*sizeof(uint32_t));

  if(ngeninfs && ninfbins) {
    *ngeninfs=(uint64_t*)malloc(INIT_NINF_ALLOC*sizeof(uint64_t));
    stats->ngeninfs=ngeninfs;
    *ninfbins=INIT_NINF_ALLOC;
    stats->ninfbins=ninfbins;
    memset(*stats->ngeninfs,0,*stats->ninfbins*sizeof(uint64_t));
  }
  stats->lmax=UINT32_MAX;
  stats->nimax=UINT32_MAX;
}
