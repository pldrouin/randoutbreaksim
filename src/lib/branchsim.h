/**
 * @file branchsim.h
 * @brief Branching simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Original model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
 * <david.maybury@tpsgc-pwgsc.gc.ca>
 */

#ifndef _BRANCHSIM_
#define _BRANCHSIM_

#include "simulation.h"

#include "ran_log.h"

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers
#define II_ARRAY_GROW_FACT (1.5)  //!< Growing factor for the array of current infectious individuals across all layers.

extern int __ro_debug;
#ifdef DEBUG_PRINTF
#undef DEBUG_PRINTF
#endif
#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) {if(__ro_debug) printf(__VA_ARGS__);} //!< Debug print function

/**
 * @brief Initialises the branching simulation.
 *
 * This function must be called to initialise the branching simulation.
 *
 * @param sv: Pointer to the simulation handle.
 */
void branchsim_init(sim_vars* sv);

/**
 * @brief Performs the branching simulation.
 *
 * Performs the branching simulation, as configured through the simulation variables. This
 * function can be called multiple times in a row.
 *
 * @param sv: Pointer to the simulation variables.
 * @return 0 if there is no error.
 */
int branchsim(sim_vars* sv);

/**
 * @brief Frees the dynamic memory used in the simulation handle.
 *
 * Does not free the memory related to the simulation-level data pointer
 * for the user-defined functions, but it does free the user-allocated memory
 * for the infectious individuals.
 *
 * @param sv: Pointer to the simulation variables.
 */
void branchsim_free(sim_vars* sv);

inline static uint32_t gen_att_infpop_log_plus_1(sim_vars* sv){return (uint32_t)ran_log_finite(&sv->rl)+1;}
inline static uint32_t gen_att_infpop_log(sim_vars* sv){return (uint32_t)ran_log_finite_gt1(&sv->rl);}
inline static uint32_t gen_att_infpop_log_p0(sim_vars* sv){return 2;}
inline static uint32_t gen_att_infpop_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); return (uint32_t)(ret+0.5);} //Not very efficient implementation

inline static void gen_att_inf_infpop_pinf1_log_plus_1(sim_vars* sv){sv->curii->ninfections=(uint32_t)ran_log_finite(&sv->rl); sv->curii->nattendees=sv->curii->ninfections+1;}
inline static void gen_att_inf_infpop_pinf1_log(sim_vars* sv){sv->curii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl); sv->curii->ninfections=sv->curii->nattendees-1;}
inline static void gen_att_inf_infpop_pinf1_log_p0(sim_vars* sv){sv->curii->ninfections=1; sv->curii->nattendees=2;}
inline static void gen_att_inf_infpop_pinf1_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); sv->curii->nattendees=(uint32_t)(ret+0.5); sv->curii->ninfections=sv->curii->nattendees-1;} //Not very efficient implementation

#ifdef DUAL_PINF
inline static void gen_att_inf_infpop_log_plus_1(sim_vars* sv){sv->curii->nattendees=(uint32_t)ran_log_finite(&sv->rl)+1; uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, sv->curii->nattendees-1); sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1-ninfpatt)+gsl_ran_binomial(sv->r, sv->pars.rpinfp*sv->pars.pinf, ninfpatt);}
inline static void gen_att_inf_infpop_log(sim_vars* sv){sv->curii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl); uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, sv->curii->nattendees-1); sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1-ninfpatt)+gsl_ran_binomial(sv->r, sv->pars.rpinfp*sv->pars.pinf, ninfpatt);}
inline static void gen_att_inf_infpop_log_p0(sim_vars* sv){sv->curii->ninfections=(gsl_rng_uniform(sv->r) < sv->pars.pinf*(1+sv->pars.ppip*(sv->pars.rpinfp-1))); sv->curii->nattendees=2;}
inline static void gen_att_inf_infpop_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); sv->curii->nattendees=(uint32_t)(ret+0.5); uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, sv->curii->nattendees-1); sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1-ninfpatt)+gsl_ran_binomial(sv->r, sv->pars.rpinfp*sv->pars.pinf, ninfpatt);} //Not very efficient implementation
#else
inline static void gen_att_inf_infpop_log_plus_1(sim_vars* sv){sv->curii->nattendees=(uint32_t)ran_log_finite(&sv->rl)+1; sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1);}
inline static void gen_att_inf_infpop_log(sim_vars* sv){sv->curii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl);  sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1);}
inline static void gen_att_inf_infpop_log_p0(sim_vars* sv){sv->curii->ninfections=(gsl_rng_uniform(sv->r) < sv->pars.pinf); sv->curii->nattendees=2;}
inline static void gen_att_inf_infpop_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); sv->curii->nattendees=(uint32_t)(ret+0.5); sv->curii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, sv->curii->nattendees-1);} //Not very efficient implementation
#endif

#define BR_GENINF_COND(EXTRA_COND) if(sv->pars.pinf==1 EXTRA_COND) { \
  if(sv->pars.grouptype&ro_group_gauss) {sv->gen_att_func=gen_att_infpop_gauss; sv->gen_att_inf_func=gen_att_inf_infpop_pinf1_gauss;} \
  else if(sv->pars.p == 0)  {sv->gen_att_func=gen_att_infpop_log_p0; sv->gen_att_inf_func=gen_att_inf_infpop_pinf1_log_p0;} \
  else if(sv->pars.grouptype&ro_group_log_plus_1) {sv->gen_att_func=gen_att_infpop_log_plus_1; sv->gen_att_inf_func=gen_att_inf_infpop_pinf1_log_plus_1;} \
  else {sv->gen_att_func=gen_att_infpop_log; sv->gen_att_inf_func=gen_att_inf_infpop_pinf1_log;} \
} else { \
  if(sv->pars.grouptype&ro_group_gauss) {sv->gen_att_func=gen_att_infpop_gauss; sv->gen_att_inf_func=gen_att_inf_infpop_gauss;} \
  else if(sv->pars.p == 0)  {sv->gen_att_func=gen_att_infpop_log_p0; sv->gen_att_inf_func=gen_att_inf_infpop_log_p0;} \
  else if(sv->pars.grouptype&ro_group_log_plus_1) {sv->gen_att_func=gen_att_infpop_log_plus_1; sv->gen_att_inf_func=gen_att_inf_infpop_log_plus_1;} \
  else {sv->gen_att_func=gen_att_infpop_log; sv->gen_att_inf_func=gen_att_inf_infpop_log;} \
}

#endif
