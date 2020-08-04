/**
 * @file standard_summary_stats.c
 * @brief User-defined functions to compute standard summary statistics.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "standard_summary_stats.h"

void std_stats_init(struct sim_vars* sv, struct std_summary_stats* stats)
{
  int i;
  for(i=sv->nlayers-1; i>=0; --i) {
    sv->iis[i].dataptr=malloc(sizeof(double));
  }
  stats->npers=(int)sv->pars.tmax+1;
  stats->inf_timeline=(uint32_t*)malloc(stats->npers*sizeof(uint32_t));
  stats->totinf_timeline=(uint32_t*)malloc(stats->npers*sizeof(uint32_t));
}

void std_stats_increase_layers(struct infindividual* iis, uint32_t n)
{
  int i;
  for(i=n-1; i>=0; --i) {
    iis[i].dataptr=malloc(sizeof(double));
  }
}

void std_stats_free(struct sim_vars* sv, struct std_summary_stats* stats)
{
  int i;
  for(i=sv->nlayers-1; i>=0; --i) {
    free(sv->iis[i].dataptr);
  }
  free(stats->inf_timeline);
  free(stats->totinf_timeline);
}
