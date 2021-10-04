/**
 * @file finitepopsim.c
 * @brief Finite population simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "finitepopsim.h"

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers
#define II_ARRAY_GROW_FACT (1.5)  //!< Growing factor for the array of current infectious individuals across all layers.

void finitepopsim_init(sim_vars* sv)
{
  sv->fpsim.rooti.ii.commpertype=0;

  sv->fpsim.is=(individual*)malloc(sv->pars.popsize*sizeof(individual));

  sv->ii_alloc_proc_func(&sv->fpsim.rooti.ii);
  for(int32_t i=sv->pars.popsize-1; i>=0; --i) sv->ii_alloc_proc_func(&sv->fpsim.is[i].ii);
  ran_log_init(&sv->rl, (rng_stream*)sv->r->state, sv->pars.p);

  /*
#ifdef DUAL_PINF
  FP_GENINF_COND(&& (sv->pars.ppip==0 || sv->pars.rpinfp==1));
#else
  FP_GENINF_COND();
#endif
*/
}

int finitepopsim(sim_vars* sv)
{
}

void finitepopsim_free(sim_vars* sv)
{
  if(sv->fpsim.rooti.ii.dataptr) free(sv->fpsim.rooti.ii.dataptr);

  for(int32_t i=sv->pars.popsize; i>=0; --i) if(sv->fpsim.is[i].ii.dataptr) free(sv->fpsim.is[i].ii.dataptr);
  free(sv->fpsim.is);
}
