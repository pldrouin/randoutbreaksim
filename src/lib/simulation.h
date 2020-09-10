/**
 * @file simulation.h
 * @brief Common simulation data structures and functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _SIMULATION_
#define _SIMULATION_

#include <stdint.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "infindividual.h"
#include "model_parameters.h"

#include "ran_log.h"

#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__) //!< Debug print function

typedef struct {
  infindividual* iis;	//!< Array of current infectious individuals across all layers
  uint32_t nlayers;	//!< Current maximum number of layers that has been used so far 
} brsim_vars;

typedef struct {
} fpsim_vars;

/**
 * Simulation variables
 */
typedef struct sim_vars_
{
  model_pars pars;		//!< Simulation input parameters
  gsl_rng const* r;		//!< Pointer to GSL random number generator
  infindividual* curii;		//!< Pointer to current iteration infectious individual
  void* dataptr;		//!< Simulation-level data pointer for user-defined functions
  void (*gen_time_periods_func)(struct sim_vars_*);				//!< Pointer to the function used to generate time periods for a given infectious individual
  uint32_t (*gen_infections_func)(struct sim_vars_*);				//!< Pointer to the function used to generate new infections during one event
  void (*ii_alloc_proc_func)(infindividual* ii);	//!< Pointer to the user-defined processing function that is called when memory for a new infectious individual is allocated.
  bool (*new_event_proc_func)(struct sim_vars_* sv);				//!< Pointer to the user-defined processing function that is called when a new transmission event is created, after an event time and the number of new infections have been assigned. The function is also called at the beginning of the simulation to account for the initial infectious individuals. The returned value from this function determines if new infectious individuals are instantiated for this event.
  void (*new_inf_proc_func)(infindividual* newinf);			//!< Pointer to the user-defined processing function that is called when a new infected individual is created, after the communicable period and the number of transmission events have been assigned. The function is only called if the number of transmission events is non-zero. 
  void (*end_inf_proc_func)(infindividual* inf, void* dataptr); 		//!< Pointer to the user-defined processing function that is called once all transmission events for a given infectious individual have been generated.
  void (*inf_proc_func_noevent)(infindividual* inf, void* dataptr);	//!< Pointer to the user-defined processing function that is called for an infectious individual that does not generate any transmission event.
  ran_log rl;	//!< Handle for the logarithmica random variate generator.

  union{
    brsim_vars brsim;
    fpsim_vars fpsim;
  };
} sim_vars;

/**
 * @brief Initialises simulation parameters.
 *
 * This function must be called to initialise simulation
 * parameters.
 *
 * @param pars: Pointer to the simulation parameters
 */
void sim_pars_init(model_pars* pars);

/**
 * @brief Initialises the simulation variables.
 *
 * Simulation variables  initialisation function.
 *
 * @param sv: Pointer to the simulation handle.
 * @param pars: Pointer to the simulation parameters
 * @param r: Pointer to the GSL random number generator.
 */
void sim_init(sim_vars* sv, model_pars const* pars, const gsl_rng* r);

/**
 * @brief Initialises the simulation-level data pointer for user-defined functions.
 *
 * This function can be used to set the simulation-level data pointer for user-defined functions.
 *
 * @param sv: Pointer to the simulation variables.
 * @param dataptr: Data pointer.
 */
inline static void sim_set_proc_data(sim_vars* sv, void* dataptr){sv->dataptr=dataptr;}

/**
 * @brief Sets the user-defined processing function that is called when a new transmission event is created.
 *
 * This function is called after an event time and the number of new infections
 * have been assigned. The function is also called at the beginning of the simulation
 * to account for the initial infectious individuals.
 *
 * @param sv: Pointer to the simulation variables.
 * @param new_event_proc_func: Pointer to the user-defined function.
 * @return true if new infectious individuals are to be instantiated from this
 * event, false otherwise.
 */
inline static void sim_set_new_event_proc_func(sim_vars* sv, bool (*new_event_proc_func)(sim_vars* sv)){sv->new_event_proc_func=new_event_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called when a new infected individual is created.
 *
 * This function is called after the communicable period and the number of
 * transmission events have been assigned. The function is only called if the number
 * of transmission events is non-zero.
 *
 * @param sv: Pointer to the simulation variables.
 * @param new_inf_proc_func: Pointer to the user-defined function.
 */
inline static void sim_set_new_inf_proc_func(sim_vars* sv, void (*new_inf_proc_func)(infindividual* newinf)){sv->new_inf_proc_func=new_inf_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called when memory
 * for a new infectious individual is allocated.
 *
 * @param sv: Pointer to the simulation variables.
 * @param ii_alloc_proc_func: Pointer to the user-defined function. The
 * argument of this function points to the new infectious individual,
 */
inline static void sim_set_ii_alloc_proc_func(sim_vars* sv, void (*ii_alloc_proc_func)(infindividual* ii)){sv->ii_alloc_proc_func=ii_alloc_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called once all
 * transmission events for a given infectious individual have been generated.
 *
 * @param sv: Pointer to the simulation variables.
 * @param end_inf_proc_func: Pointer to the user-defined function. The second
 * argument for this function is the simulation-level data pointer.
 */
inline static void sim_set_end_inf_proc_func(sim_vars* sv, void (*end_inf_proc_func)(infindividual* inf, void* dataptr)){sv->end_inf_proc_func=end_inf_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called for an infectious individual that does not generate any transmission event.
 *
 * This function is called after the communicable period has been assigned for
 * the current infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param inf_proc_func_noevent: Pointer to the user-defined function. The second
 * argument for this function is the simulation-level data pointer.
 */
inline static void sim_set_inf_proc_noevent_func(sim_vars* sv, void (*inf_proc_func_noevent)(infindividual* inf, void* dataptr)){sv->inf_proc_func_noevent=inf_proc_func_noevent;}

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main GEN_PER macro.
 */
#define GEN_PER_LATENT_0 sv->curii->latent_period=0;
#define GEN_PER_LATENT_1 sv->curii->latent_period=sv->pars.lbar;
#define GEN_PER_LATENT_2 sv->curii->latent_period=gsl_ran_gamma(sv->r, sv->pars.kappal*sv->pars.lbar, 1./sv->pars.kappal);

#define GEN_PER_INTERRUPTED_MAIN_0
#define GEN_PER_INTERRUPTED_MAIN_1 if(gsl_rng_uniform(sv->r) < sv->pars.pit && sv->pars.itbar < sv->curii->comm_period) sv->curii->comm_period=sv->pars.itbar;
#define GEN_PER_INTERRUPTED_MAIN_2 if(gsl_rng_uniform(sv->r) < sv->pars.pit) {const double time=gsl_ran_gamma(sv->r, sv->pars.kappait*sv->pars.itbar, 1./sv->pars.kappait); if(time < sv->curii->comm_period) sv->curii->comm_period=time;}

#define GEN_PER_MAIN_1(IT) sv->curii->comm_period=sv->pars.tbar; GEN_PER_INTERRUPTED_MAIN_ ## IT;
#define GEN_PER_MAIN_2(IT) sv->curii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa); GEN_PER_INTERRUPTED_MAIN_ ## IT;

#define GEN_PER_INTERRUPTED_ALT_0
#define GEN_PER_INTERRUPTED_ALT_1 if(gsl_rng_uniform(sv->r) < sv->pars.pim && sv->pars.imbar < sv->curii->comm_period) sv->curii->comm_period=sv->pars.imbar;
#define GEN_PER_INTERRUPTED_ALT_2 if(gsl_rng_uniform(sv->r) < sv->pars.pim) {const double time=gsl_ran_gamma(sv->r, sv->pars.kappaim*sv->pars.imbar, 1./sv->pars.kappaim); if(time < sv->curii->comm_period) sv->curii->comm_period=time;}

#define GEN_PER_ALTERNATE_0(MAIN,IT,IM) GEN_PER_MAIN_ ## MAIN(IT);
#define GEN_PER_ALTERNATE_1(MAIN,IT,IM) if(gsl_rng_uniform(sv->r) < sv->pars.q) {sv->curii->comm_period=sv->pars.mbar; GEN_PER_INTERRUPTED_ALT_ ## IM} else GEN_PER_MAIN_ ## MAIN(IT);
#define GEN_PER_ALTERNATE_2(MAIN,IT,IM) if(gsl_rng_uniform(sv->r) < sv->pars.q) {sv->curii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappaq*sv->pars.mbar, 1./sv->pars.kappaq); GEN_PER_INTERRUPTED_ALT_ ## IM} else GEN_PER_MAIN_ ## MAIN(IT);
//! @endcond

/**
 * @brief Generates one function that generates a specific combination of
 * communicable and latent periods.
 *
 * This preprocessing macro generates a function that
 * generates the communicable period for the current
 * infectious individual and a latent period, as determined by the
 * simulation variables and parameters. It considers the alternate
 * communicable period based on the q simulation parameter, as well
 * as the two types of interrupted periods. It also sets the
 * infectious_at_tmax parameter for the current infectious
 * individual.
 *
 * @param LATENT: Type of latent period (0=none, 1=fixed, 2=variable)
 * @param MAIN: Type of main period (1=fixed, 2=variable)
 * @param IT: Type of interrupted main period (0=none, 1=fixed, 2=variable)
 * @param ALTERNATE: Type of alternate period (0=none, 1=fixed, 2=variable)
 * @param IM: Type of interrupted alternate period (0=none, 1=fixed, 2=variable)
 * latent)
 */
#define GEN_PER(LATENT,MAIN,IT,ALTERNATE,IM) static inline void gen_comm_ ## LATENT ## _ ## MAIN ## _ ## IT ## _ ## ALTERNATE ## _ ## IM ## _periods(sim_vars* sv) \
{ \
  GEN_PER_LATENT_ ## LATENT \
  GEN_PER_ALTERNATE_ ## ALTERNATE(MAIN,IT,IM) \
}

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main GEN_PERS_MAIN macro.
 */
#define GEN_PERS_IM(LATENT,MAIN,IT,ALTERNATE) GEN_PER(LATENT,MAIN,IT,ALTERNATE,0) GEN_PER(LATENT,MAIN,IT,ALTERNATE,1) GEN_PER(LATENT,MAIN,IT,ALTERNATE,2)
#define GEN_PERS_ALT(LATENT,MAIN,IT) GEN_PER(LATENT,MAIN,IT,0,0) GEN_PERS_IM(LATENT,MAIN,IT,1) GEN_PERS_IM(LATENT,MAIN,IT,2)
#define GEN_PERS_IT(LATENT,MAIN) GEN_PERS_ALT(LATENT,MAIN,0) GEN_PERS_ALT(LATENT,MAIN,1) GEN_PERS_ALT(LATENT,MAIN,2)
//! @endcond

/**
 * @brief For a given type of latent period, generates functions that each
 * generates a specific combination of communicable and latent periods.
 *
 * @param LATENT: Type of latent period (0=none, 1=fixed, 2=variable)
 */
#define GEN_PERS_MAIN(LATENT) GEN_PERS_IT(LATENT,1) GEN_PERS_IT(LATENT,2)
GEN_PERS_MAIN(0)
GEN_PERS_MAIN(1)
GEN_PERS_MAIN(2)

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main PER_COND macro.
 */
#define PER_COND_IM(IMPREFIX) if(!(sv->pars.pim>0)) sv->gen_time_periods_func=gen_comm_ ## IMPREFIX ## _0_periods; else if(isinf(sv->pars.kappaim)) sv->gen_time_periods_func=gen_comm_ ## IMPREFIX ## _1_periods; else sv->gen_time_periods_func=gen_comm_ ## IMPREFIX ## _2_periods;
#define PER_COND_ALT(ALTPREFIX) if(!(sv->pars.q>0)) sv->gen_time_periods_func=gen_comm_ ## ALTPREFIX ## _0_0_periods; else if(isinf(sv->pars.kappaq)) {PER_COND_IM(ALTPREFIX ## _1)} else {PER_COND_IM(ALTPREFIX ## _2)};
#define PER_COND_IT(ITPREFIX) if(!(sv->pars.pit>0)) {PER_COND_ALT(ITPREFIX ## _0)} else if(isinf(sv->pars.kappait)) {PER_COND_ALT(ITPREFIX ## _1)} else {PER_COND_ALT(ITPREFIX ## _2)};
#define PER_COND_MAIN(MAINPREFIX) if(isinf(sv->pars.kappa)) {PER_COND_IT(MAINPREFIX ## _1)} else {PER_COND_IT(MAINPREFIX ## _2)};
//! @endcond

/**
 * @brief Generates the conditional code that assigns the function that
 * generates a specific combination of communicable and latent periods.
 *
 * This preprocessing macro generates the conditional code that assigns the
 * function that generates a specific combination of communicable and latent
 * periods to the gen_time_periods_func pointer of the simulation, given the
 * configuration of the model parameters.
 */
#define PER_COND if(isnan(sv->pars.kappal)) {PER_COND_MAIN(0)} else if(isinf(sv->pars.kappal)) {PER_COND_MAIN(1)} else {PER_COND_MAIN(2)};

/**
 * @brief Default processing function that is called when a new transmission event is created.
 *
 * This function is called by default if a user-defined function has not been
 * set.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event does not occur after tmax, false otherwise.
 */
inline static bool default_event_proc_func(sim_vars* sv){return (sv->curii->event_time <= sv->pars.tmax);}

/**
 * @brief Default processing function that is called memory for a new infectious
 * individual is allocated.
 *
 * @param ii: New infectious individual.
 */
inline static void default_ii_alloc_proc_func(infindividual* ii){ii->dataptr=NULL;}

/**
 * @brief Default processing function.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 */
inline static void dummy_proc_func_sv(sim_vars* sv){}

/**
 * @brief Default processing function.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 */
inline static void dummy_proc_func_one_par(infindividual* inf){}

/**
 * @brief Default processing function.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 */
inline static void dummy_proc_func_two_pars(infindividual* inf, void* ptr){}

#endif
