/**
 * @file finitepopsim.h
 * @brief Finite population simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _FINITEPOPSIM_
#define _FINITEPOPSIM_

#include <stdlib.h>

#include "simulation.h"

#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
//#define DEBUG_PRINTF(...) //!< Debug print function
#define DEBUG_PRINTF(...) printf(__VA_ARGS__) //!< Debug print function

/**
 * @brief Initialises the branching simulation.
 *
 * This function must be called to initialise the branching simulation.
 *
 * @param sv: Pointer to the simulation handle.
 */
void finitepopsim_init(sim_vars* sv);

/**
 * @brief Performs the branching simulation.
 *
 * Performs the branching simulation, as configured through the simulation variables. This
 * function can be called multiple times in a row.
 *
 * @param sv: Pointer to the simulation variables.
 * @return 0 if there is no error.
 */
int finitepopsim(sim_vars* sv);

/**
 * @brief Frees the dynamic memory used in the simulation handle.
 *
 * Does not free the memory related to the simulation-level data pointer
 * for the user-defined functions, but it does free the user-allocated memory
 * for the infectious individuals.
 *
 * @param sv: Pointer to the simulation variables.
 */
void finitepopsim_free(sim_vars* sv);

inline static void fp_sort_changetimes(fpsim_vars* fpsim)
{
  qsort(fpsim->activated, fpsim->nactivated, sizeof(individual*), ind_ct_comp);
}

inline static uint32_t gen_att_finpop_log_plus_1(sim_vars* sv){return (uint32_t)ran_log_capped(&sv->rl, sv->pars.popsize-1)+1;}
inline static uint32_t gen_att_finpop_log(sim_vars* sv){return (uint32_t)ran_log_capped_gt1(&sv->rl, sv->pars.popsize);}
inline static uint32_t gen_att_finpop_log_p0(sim_vars* sv){return 2;}
inline static uint32_t gen_att_finpop_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5 && ret>=sv->pars.popsize+0.5); return (uint32_t)(ret+0.5);} //Not very efficient implementation

#define FP_GENINF_COND() { \
  if(sv->pars.grouptype&ro_group_gauss) {sv->gen_att_func=gen_att_finpop_gauss;} \
  else if(sv->pars.p == 0)  {sv->gen_att_func=gen_att_finpop_log_p0;} \
  else if(sv->pars.grouptype&ro_group_log_plus_1) {sv->gen_att_func=gen_att_finpop_log_plus_1;} \
  else {sv->gen_att_func=gen_att_finpop_log;} \
}

#endif
