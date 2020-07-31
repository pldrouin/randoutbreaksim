#ifndef _STANDARD_SUMMARY_STATS_
#define _STANDARD_SUMMARY_STATS_

#include <stdlib.h>
#include <stdbool.h>

#include "infindividual.h"
#include "simulation.h"

#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
#define DEBUG_PRINTF(...)
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__)

struct std_summary_stats
{
  double extinction_time;
  double commpersum;
  uint32_t* inf_timeline;
  uint32_t* totinf_timeline;
  uint32_t npers;
  uint32_t rsum;
  uint32_t neventssum;
  //uint32_t n_ended_infections;
  uint32_t total_n_infections;
  bool extinction;
};

void std_stats_init(struct sim_vars* sv, struct std_summary_stats* stats);

inline static void std_stats_path_init(struct std_summary_stats* stats)
{
  stats->extinction_time=0;
  stats->commpersum=0;
  memset(stats->inf_timeline,0,stats->npers*sizeof(uint32_t));
  memset(stats->totinf_timeline,0,stats->npers*sizeof(uint32_t));
  stats->rsum=0;
  stats->neventssum=0;
  stats->total_n_infections=0;
  stats->extinction=true;
}

void std_stats_increase_layers(struct infindividual* iis, uint32_t n);
void std_stats_free(struct sim_vars* sv, struct std_summary_stats* stats);

inline static bool std_stats_new_event(struct sim_vars* sv)
{
  (*(uint32_t*)sv->ii->dataptr)+=sv->ii->ninfections;
  DEBUG_PRINTF("Number of infections incremented to %u\n",*(uint32_t*)sv->ii->dataptr);
  return (sv->ii->event_time <= sv->pars.tmax);
}

inline static void std_stats_new_inf(struct infindividual* inf)
{
  *(uint32_t*)inf->dataptr=0;
  //++(*(uint32_t*)(inf-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(inf-1)->dataptr);
}

inline static void std_stats_end_inf(struct infindividual* inf, void* ptr)
{
  ++((struct std_summary_stats*)ptr)->total_n_infections;
  DEBUG_PRINTF("Number of infections was %u\n",*(uint32_t*)inf->dataptr);
  ((struct std_summary_stats*)ptr)->rsum+=*(uint32_t*)inf->dataptr;
  ((struct std_summary_stats*)ptr)->commpersum+=inf->comm_period;
  ((struct std_summary_stats*)ptr)->neventssum+=inf->nevents;

  //If truncated by tmax
  if(inf->infectious_at_tmax) ((struct std_summary_stats*)ptr)->extinction=false;

  else {
    //++((struct std_summary_stats*)ptr)->n_ended_infections;
    const double inf_end=(inf-1)->event_time+inf->comm_period;

    if(inf_end > ((struct std_summary_stats*)ptr)->extinction_time) ((struct std_summary_stats*)ptr)->extinction_time=inf_end;
  }

  int i;
  int end_comm_per=(int)((inf-1)->event_time+inf->comm_period);

  if(end_comm_per>=((struct std_summary_stats*)ptr)->npers) end_comm_per=((struct std_summary_stats*)ptr)->npers-1;

  for(i=(int)((inf-1)->event_time); i<=end_comm_per; ++i) ++(((struct std_summary_stats*)ptr)->inf_timeline[i]);
  ++(((struct std_summary_stats*)ptr)->totinf_timeline[(int)((inf-1)->event_time)]);
}

inline static void std_stats_noevent_inf(struct infindividual* inf, void* ptr)
{
  //++(*(uint32_t*)(inf-1)->dataptr);
  //DEBUG_PRINTF("Number of parent infections incremented to %u\n",*(uint32_t*)(inf-1)->dataptr);
  ++((struct std_summary_stats*)ptr)->total_n_infections;
  DEBUG_PRINTF("Number of infections was 0\n");
  ((struct std_summary_stats*)ptr)->commpersum+=inf->comm_period;
  ((struct std_summary_stats*)ptr)->neventssum+=inf->nevents;

  //If truncated by tmax
  if(inf->infectious_at_tmax) ((struct std_summary_stats*)ptr)->extinction=false;

  else {
    //++((struct std_summary_stats*)ptr)->n_ended_infections;
    const double inf_end=(inf-1)->event_time+inf->comm_period;

    if(inf_end > ((struct std_summary_stats*)ptr)->extinction_time) ((struct std_summary_stats*)ptr)->extinction_time=inf_end;
  }

  int i;
  int end_comm_per=(int)((inf-1)->event_time+inf->comm_period);

  if(end_comm_per>=((struct std_summary_stats*)ptr)->npers) end_comm_per=((struct std_summary_stats*)ptr)->npers-1;

  for(i=(int)((inf-1)->event_time); i<=end_comm_per; ++i) ++(((struct std_summary_stats*)ptr)->inf_timeline[i]);
  ++(((struct std_summary_stats*)ptr)->totinf_timeline[(int)((inf-1)->event_time)]);
}

#endif
