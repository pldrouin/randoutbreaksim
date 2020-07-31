#include "standard_summary_stats.h"

void std_stats_init(struct sim_vars* sv)
{
  int i;
  for(i=sv->nlayers-1; i>0; --i) {
    sv->iis[i].dataptr=malloc(sizeof(double));
  }
}

void std_stats_increase_layers(struct infindividual* iis, uint32_t n)
{
  int i;
  for(i=n-1; i>=0; --i) {
    iis[i].dataptr=malloc(sizeof(double));
  }
}

void std_stats_free(struct sim_vars* sv)
{
  int i;
  for(i=sv->nlayers-1; i>0; --i) {
    free(sv->iis[i].dataptr);
  }
}
