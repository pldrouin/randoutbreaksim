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

typedef struct inflayer_ {
  infindividual ii;
  double last_event_time;     //!< Last event time used for this layer, as recorded before moving to the next layer, and as retrieved when going back to the previous layer
  uint32_t cureventi;	      //!< Index of the current iteration event
  uint32_t curinfectioni;     //!< Index of the current iteration infection
} inflayer;

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
inline static uint32_t gen_att_infpop_geom(sim_vars* sv){return 1+gsl_ran_geometric(sv->r, 1-sv->pars.p);}
inline static uint32_t gen_att_infpop_gauss(sim_vars* sv){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); return (uint32_t)(ret+0.5);} //Not very efficient implementation

inline static void gen_att_inf_infpop_pinf1_log_plus_1(sim_vars* sv, infindividual* ii){ii->ninfections=(uint32_t)ran_log_finite(&sv->rl); ii->nattendees=ii->ninfections+1;}
inline static void gen_att_inf_infpop_pinf1_log(sim_vars* sv, infindividual* ii){ii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl); ii->ninfections=ii->nattendees-1;}
inline static void gen_att_inf_infpop_pinf1_log_p0(sim_vars* sv, infindividual* ii){ii->ninfections=1; ii->nattendees=2;}
inline static void gen_att_inf_infpop_pinf1_geom(sim_vars* sv, infindividual* ii){ii->ninfections=gsl_ran_geometric(sv->r, 1-sv->pars.p); ii->nattendees=ii->ninfections+1;}
inline static void gen_att_inf_infpop_pinf1_gauss(sim_vars* sv, infindividual* ii){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); ii->nattendees=(uint32_t)(ret+0.5); ii->ninfections=ii->nattendees-1;} //Not very efficient implementation

#ifdef DUAL_PINF
inline static void gen_att_inf_infpop_log_plus_1(sim_vars* sv, infindividual* ii){ii->nattendees=(uint32_t)ran_log_finite(&sv->rl)+1; uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, ii->nattendees-1); ii->ninfectionsf=gsl_ran_binomial(sv->r, ii->pinf, ii->nattendees-1-ninfpatt); ii->ninfectionsp=gsl_ran_binomial(sv->r, sv->pars.rpinfp*ii->pinf, ninfpatt); ii->ninfections=ii->ninfectionsf+ii->ninfectionsp;}
inline static void gen_att_inf_infpop_log(sim_vars* sv, infindividual* ii){ii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl); uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, ii->nattendees-1); ii->ninfectionsf=gsl_ran_binomial(sv->r, ii->pinf, ii->nattendees-1-ninfpatt); ii->ninfectionsp=gsl_ran_binomial(sv->r, sv->pars.rpinfp*ii->pinf, ninfpatt); ii->ninfections=ii->ninfectionsf+ii->ninfectionsp;}
inline static void gen_att_inf_infpop_log_p0(sim_vars* sv, infindividual* ii){if(gsl_rng_uniform(sv->r) < sv->pars.ppip) {ii->ninfections=ii->ninfectionsp=(gsl_rng_uniform(sv->r) < ii->pinf*sv->pars.rpinfp); ii->ninfectionsf=0;} else {ii->ninfections=ii->ninfectionsf=(gsl_rng_uniform(sv->r) < ii->pinf); ii->ninfectionsp=0;} ii->nattendees=2;}
inline static void gen_att_inf_infpop_geom(sim_vars* sv, infindividual* ii){ii->nattendees=1+gsl_ran_geometric(sv->r, 1-sv->pars.p); uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, ii->nattendees-1); ii->ninfectionsf=gsl_ran_binomial(sv->r, ii->pinf, ii->nattendees-1-ninfpatt); ii->ninfectionsp=gsl_ran_binomial(sv->r, sv->pars.rpinfp*ii->pinf, ninfpatt); ii->ninfections=ii->ninfectionsf+ii->ninfectionsp;}
inline static void gen_att_inf_infpop_gauss(sim_vars* sv, infindividual* ii){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); ii->nattendees=(uint32_t)(ret+0.5); uint32_t ninfpatt=gsl_ran_binomial(sv->r, sv->pars.ppip, ii->nattendees-1); ii->ninfectionsf=gsl_ran_binomial(sv->r, ii->pinf, ii->nattendees-1-ninfpatt); ii->ninfectionsp=gsl_ran_binomial(sv->r, sv->pars.rpinfp*ii->pinf, ninfpatt); ii->ninfections=ii->ninfectionsf+ii->ninfectionsp;} //Not very efficient implementation
#else
inline static void gen_att_inf_infpop_log_plus_1(sim_vars* sv, infindividual* ii){ii->nattendees=(uint32_t)ran_log_finite(&sv->rl)+1; ii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, ii->nattendees-1);}
inline static void gen_att_inf_infpop_log(sim_vars* sv, infindividual* ii){ii->nattendees=(uint32_t)ran_log_finite_gt1(&sv->rl);  ii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, ii->nattendees-1);}
inline static void gen_att_inf_infpop_log_p0(sim_vars* sv, infindividual* ii){ii->ninfections=(gsl_rng_uniform(sv->r) < sv->pars.pinf); ii->nattendees=2;}
inline static void gen_att_inf_infpop_geom(sim_vars* sv, infindividual* ii){ii->nattendees=1+gsl_ran_geometric(sv->r, 1-sv->pars.p);  ii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, ii->nattendees-1);}
inline static void gen_att_inf_infpop_gauss(sim_vars* sv, infindividual* ii){double ret; do {ret=sv->pars.mu+gsl_ran_gaussian_ziggurat(sv->r,sv->pars.sigma);} while(ret<1.5); ii->nattendees=(uint32_t)(ret+0.5); ii->ninfections=gsl_ran_binomial(sv->r, sv->pars.pinf, ii->nattendees-1);} //Not very efficient implementation
#endif

#define BR_GENINF_COND(EXTRA_COND) if(sv->pars.pinf==1 EXTRA_COND) { \
  if(sv->pars.grouptype&ro_group_gauss) {sv->gen_att_func=gen_att_infpop_gauss; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_pinf1_gauss;} \
  else if(sv->pars.grouptype&ro_group_geom) {sv->gen_att_func=gen_att_infpop_geom; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_pinf1_geom;} \
  else if(sv->pars.p == 0)  {sv->gen_att_func=gen_att_infpop_log_p0; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_pinf1_log_p0;} \
  else if(sv->pars.grouptype&ro_group_log_plus_1) {sv->gen_att_func=gen_att_infpop_log_plus_1; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_pinf1_log_plus_1;} \
  else {sv->gen_att_func=gen_att_infpop_log; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_pinf1_log;} \
} else { \
  if(sv->pars.grouptype&ro_group_gauss) {sv->gen_att_func=gen_att_infpop_gauss; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_gauss;} \
  else if(sv->pars.grouptype&ro_group_geom) {sv->gen_att_func=gen_att_infpop_geom; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_geom;} \
  else if(sv->pars.p == 0)  {sv->gen_att_func=gen_att_infpop_log_p0; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_log_p0;} \
  else if(sv->pars.grouptype&ro_group_log_plus_1) {sv->gen_att_func=gen_att_infpop_log_plus_1; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_log_plus_1;} \
  else {sv->gen_att_func=gen_att_infpop_log; sv->brsim.gen_att_inf_func=gen_att_inf_infpop_log;} \
}

#endif
