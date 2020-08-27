/**
 * @file standard_summary_stats.c
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "standard_summary_stats.h"

void std_stats_init(sim_vars* sv)
{
  std_summary_stats* stats=(std_summary_stats*)sv->dataptr;
  stats->npers=(int)sv->pars.tmax+1;
  stats->inf_timeline=(uint32_t*)malloc(stats->npers*sizeof(uint32_t));
  stats->totinf_timeline=(uint32_t*)malloc(stats->npers*sizeof(uint32_t));
  stats->nimax=UINT32_MAX;
}
