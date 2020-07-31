#ifndef _STANDARD_SUMMARY_STATS_
#define _STANDARD_SUMMARY_STATS_

#include <stdlib.h>
#include <stdbool.h>

#include "infindividual.h"

#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
#define DEBUG_PRINTF(...)
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__)

struct std_summary_stats
{
  double extinction_time;
  double commpersum;
  uint32_t rsum;
  uint32_t neventssum;
  //uint32_t n_ended_infections;
  uint32_t total_n_infections;
  bool extinction;
};

inline static void std_stats_pri_inf(struct infindividual* inf)
{
  DEBUG_PRINTF("Alloc for primary %p\n",inf);
  inf->dataptr=malloc(sizeof(uint32_t));
  *(uint32_t*)inf->dataptr=0;
}

inline static bool std_stats_new_event(struct sim_vars* sv)
{
  (*(uint32_t*)sv->ii->dataptr)+=sv->ii->ninfections;
  DEBUG_PRINTF("Number of infections incremented to %u\n",*(uint32_t*)sv->ii->dataptr);
  return (sv->ii->event_time <= sv->pars.tmax);
}

inline static void std_stats_new_inf(struct infindividual* inf)
{
  DEBUG_PRINTF("Alloc for %p\n",inf);
  inf->dataptr=malloc(sizeof(uint32_t));
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
  DEBUG_PRINTF("Free for %p\n",inf);
  free(inf->dataptr);
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
}

#endif
