/**
 * @file simulation.h
 * @brief Simulation data structures and functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Original model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
 * <david.maybury@tpsgc-pwgsc.gc.ca>
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

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers

#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__) //!< Debug print function

/**
 * Simulation variables
 */
typedef struct sim_vars_
{
  model_pars pars;		//!< Simulation input parameters
  gsl_rng const* r;		//!< Pointer to GSL random number generator
  infindividual* iis;	//!< Array of current infectious individuals across all layers
  infindividual* ii;	//!< Pointer to current iteration infectious individual
  uint32_t nlayers;		//!< Current maximum number of layers that has been used so far 
  void* dataptr;		//!< Simulation-level data pointer for user-defined functions
  void (*gen_time_periods_func)(struct sim_vars_*);				//!< Pointer to the function used to generate time periods for a given infectious individual
  void (*increase_layers_proc_func)(infindividual* iis, uint32_t n);	//!< Pointer to the user-defined processing function that is called when the maximum number of layers is increased.
  bool (*new_event_proc_func)(struct sim_vars_* sv);				//!< Pointer to the user-defined processing function that is called when a new transmission event is created, after an event time and the number of new infections have been assigned. The function is also called at the beginning of the simulation to account for the initial infectious individuals. The returned value from this function determines if new infectious individuals are instantiated for this event.
  void (*new_inf_proc_func)(infindividual* newinf);			//!< Pointer to the user-defined processing function that is called when a new infected individual is created, after the communicable period and the number of transmission events have been assigned. The function is only called if the number of transmission events is non-zero. 
  void (*end_inf_proc_func)(infindividual* inf, void* dataptr); 		//!< Pointer to the user-defined processing function that is called once all transmission events for a given infectious individual have been generated.
  void (*inf_proc_func_noevent)(infindividual* inf, void* dataptr);	//!< Pointer to the user-defined processing function that is called for an infectious individual that does not generate any transmission event.
} sim_vars;

/**
 * @brief Verifies the validity of simulation parameters.
 *
 * This internal function verifies if the set of provided simulation parameter values are
 * valid.
 *
 * @param pars: Pointer to the simulation parameters.
 * @return 0 if the parameters are valid. If the parameters are
 * invalid, an error is printed on stderr.
 */
int sim_pars_check(model_pars const* pars);

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
 * @brief Initialises the simulation.
 *
 * This function must be called to initialise the simulation.
 *
 * @param sv: Pointer to the simulation handle.
 * @param pars: Pointer to the simulation parameters
 * @param r: Pointer to the GSL random number generator.
 * @return 0 if the simulation was initialised successfully.
 */
int sim_init(sim_vars* sv, model_pars* pars, const gsl_rng* r);

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
 * @brief Sets the user-defined processing function that is called when the maximum number of layers is increased.
 *
 * This function is called after the array of infectious individuals has been
 * resized.
 *
 * @param sv: Pointer to the simulation variables.
 * @param increase_layers_proc_func: Pointer to the user-defined function. The
 * first argument of this function points to the infectious individuals array,
 * starting with the first newly allocated layer, and the second argument is the
 * number of layers that have been added.
 */
inline static void sim_set_increase_layers_proc_func(sim_vars* sv, void (*increase_layers_proc_func)(infindividual* iis, uint32_t n)){sv->increase_layers_proc_func=increase_layers_proc_func;}

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

/**
 * @brief Frees the dynamic memory used in the simulation handle.
 *
 * Does not free
 * any memory related to the simulation-level data pointer for the user-defined
 * functions.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void sim_free(sim_vars* sv){free(sv->iis);}

/**
 * @brief Performs the simulation.
 *
 * Performs the simulation, as configured through the simulation variables. This
 * function can be called multiple times in a row.
 *
 * @param sv: Pointer to the simulation variables.
 * @return 0 if there is no error.
 */
int simulate(sim_vars* sv);

//! @cond Doxygen_Suppress
#define GEN_PER_LATENT_0 sv->ii->latent_period=0;
#define GEN_PER_LATENT_1 sv->ii->latent_period=sv->pars.lbar;
#define GEN_PER_LATENT_2 sv->ii->latent_period=gsl_ran_gamma(sv->r, sv->pars.kappal*sv->pars.lbar, 1./sv->pars.kappal);

#define GEN_PER_INTERRUPTED_MAIN_0
#define GEN_PER_INTERRUPTED_MAIN_1 if(gsl_rng_uniform(sv->r) < sv->pars.pit && sv->pars.itbar < sv->ii->comm_period) sv->ii->comm_period=sv->pars.itbar;
#define GEN_PER_INTERRUPTED_MAIN_2 if(gsl_rng_uniform(sv->r) < sv->pars.pit) {const double time=gsl_ran_gamma(sv->r, sv->pars.kappait*sv->pars.itbar, 1./sv->pars.kappait); if(time < sv->ii->comm_period) sv->ii->comm_period=time;}

#define GEN_PER_MAIN_1(IT) sv->ii->comm_period=sv->pars.tbar; GEN_PER_INTERRUPTED_MAIN_##IT;
#define GEN_PER_MAIN_2(IT) sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa); GEN_PER_INTERRUPTED_MAIN_##IT;

#define GEN_PER_INTERRUPTED_ALT_0
#define GEN_PER_INTERRUPTED_ALT_1 if(gsl_rng_uniform(sv->r) < sv->pars.pim && sv->pars.imbar < sv->ii->comm_period) sv->ii->comm_period=sv->pars.imbar;
#define GEN_PER_INTERRUPTED_ALT_2 if(gsl_rng_uniform(sv->r) < sv->pars.pim) {const double time=gsl_ran_gamma(sv->r, sv->pars.kappaim*sv->pars.imbar, 1./sv->pars.kappaim); if(time < sv->ii->comm_period) sv->ii->comm_period=time;}

#define GEN_PER_ALTERNATE_0(MAIN,IT,IM) GEN_PER_MAIN_##MAIN(IT);
#define GEN_PER_ALTERNATE_1(MAIN,IT,IM) if(gsl_rng_uniform(sv->r) < sv->pars.q) {sv->ii->comm_period=sv->pars.mbar; GEN_PER_INTERRUPTED_ALT_##IM} else GEN_PER_MAIN_##MAIN(IT);
#define GEN_PER_ALTERNATE_2(MAIN,IT,IM) if(gsl_rng_uniform(sv->r) < sv->pars.q) {sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappaq*sv->pars.mbar, 1./sv->pars.kappaq); GEN_PER_INTERRUPTED_ALT_##IM} else GEN_PER_MAIN_##MAIN(IT);
//! @endcond

/**
 * @brief Generates a communicable period and a latent period.
 *
 * This function generates the communicable period for the current
 * infectious individual and a latent period, as determined by the
 * simulation variables and parameters. It considers the alternate
 * communicable period based on the q simulation parameter. It also
 * sets the infectious_at_tmax parameter for the current infectious
 * individual.
 */
#define GEN_PER(NAME,LATENT,MAIN,ALTERNATE,IT,IM) static inline void NAME(sim_vars* sv) \
{ \
  GEN_PER_LATENT_##LATENT \
  GEN_PER_ALTERNATE_##ALTERNATE(MAIN,IT,IM) \
  double time_left=sv->pars.tmax-(sv->ii-1)->event_time; \
 \
  if(sv->ii->comm_period > time_left) { \
    sv->ii->infectious_at_tmax=true; \
 \
  } else sv->ii->infectious_at_tmax=false; \
  DEBUG_PRINTF("Comm period is %f%s\n",sv->ii->comm_period,(sv->ii->infectious_at_tmax?" (reached end)":"")); \
}

//! @cond Doxygen_Suppress
GEN_PER(gen_fixed_comm_period,0,1,0,0,0);
GEN_PER(gen_fixed_comm_fixed_int_periods,0,1,0,1,0);
GEN_PER(gen_fixed_comm_int_periods,0,1,0,2,0);
GEN_PER(gen_comm_period,0,2,0,0,0);
GEN_PER(gen_comm_fixed_int_periods,0,2,0,1,0);
GEN_PER(gen_comm_int_periods,0,2,0,2,0);


GEN_PER(gen_fixed_comm_fixed_alternate_periods,0,1,1,0,0);
GEN_PER(gen_fixed_comm_fixed_alternate_fixed_int_periods,0,1,1,0,1);
GEN_PER(gen_fixed_comm_fixed_alternate_int_periods,0,1,1,0,2);
GEN_PER(gen_fixed_comm_fixed_int_fixed_alternate_periods,0,1,1,1,0);
GEN_PER(gen_fixed_comm_fixed_int_fixed_alternate_fixed_int_periods,0,1,1,1,1);
GEN_PER(gen_fixed_comm_fixed_int_fixed_alternate_int_periods,0,1,1,1,2);
GEN_PER(gen_fixed_comm_int_fixed_alternate_periods,0,1,1,2,0);
GEN_PER(gen_fixed_comm_int_fixed_alternate_fixed_int_periods,0,1,1,2,1);
GEN_PER(gen_fixed_comm_int_fixed_alternate_int_periods,0,1,1,2,2);

GEN_PER(gen_fixed_comm_alternate_periods,0,1,2,0,0);
GEN_PER(gen_fixed_comm_alternate_fixed_int_periods,0,1,2,0,1);
GEN_PER(gen_fixed_comm_alternate_int_periods,0,1,2,0,2);
GEN_PER(gen_fixed_comm_fixed_int_alternate_periods,0,1,2,1,0);
GEN_PER(gen_fixed_comm_fixed_int_alternate_fixed_int_periods,0,1,2,1,1);
GEN_PER(gen_fixed_comm_fixed_int_alternate_int_periods,0,1,2,1,2);
GEN_PER(gen_fixed_comm_int_alternate_periods,0,1,2,2,0);
GEN_PER(gen_fixed_comm_int_alternate_fixed_int_periods,0,1,2,2,1);
GEN_PER(gen_fixed_comm_int_alternate_int_periods,0,1,2,2,2);

GEN_PER(gen_comm_fixed_alternate_periods,0,2,1,0,0);
GEN_PER(gen_comm_fixed_alternate_fixed_int_periods,0,2,1,0,1);
GEN_PER(gen_comm_fixed_alternate_int_periods,0,2,1,0,2);
GEN_PER(gen_comm_fixed_int_fixed_alternate_periods,0,2,1,1,0);
GEN_PER(gen_comm_fixed_int_fixed_alternate_fixed_int_periods,0,2,1,1,1);
GEN_PER(gen_comm_fixed_int_fixed_alternate_int_periods,0,2,1,1,2);
GEN_PER(gen_comm_int_fixed_alternate_periods,0,2,1,2,0);
GEN_PER(gen_comm_int_fixed_alternate_fixed_int_periods,0,2,1,2,1);
GEN_PER(gen_comm_int_fixed_alternate_int_periods,0,2,1,2,2);

GEN_PER(gen_comm_alternate_periods,0,2,2,0,0);
GEN_PER(gen_comm_alternate_fixed_int_periods,0,2,2,0,1);
GEN_PER(gen_comm_alternate_int_periods,0,2,2,0,2);
GEN_PER(gen_comm_fixed_int_alternate_periods,0,2,2,1,0);
GEN_PER(gen_comm_fixed_int_alternate_fixed_int_periods,0,2,2,1,1);
GEN_PER(gen_comm_fixed_int_alternate_int_periods,0,2,2,1,2);
GEN_PER(gen_comm_int_alternate_periods,0,2,2,2,0);
GEN_PER(gen_comm_int_alternate_fixed_int_periods,0,2,2,2,1);
GEN_PER(gen_comm_int_alternate_int_periods,0,2,2,2,2);


GEN_PER(gen_fixed_latent_fixed_comm_periods,1,1,0,0,0);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_periods,1,1,0,1,0);
GEN_PER(gen_fixed_latent_fixed_comm_int_periods,1,1,0,2,0);

GEN_PER(gen_fixed_latent_comm_periods,1,2,0,0,0);
GEN_PER(gen_fixed_latent_comm_fixed_int_periods,1,2,0,1,0);
GEN_PER(gen_fixed_latent_comm_int_periods,1,2,0,2,0);

GEN_PER(gen_latent_fixed_comm_periods,2,1,0,0,0);
GEN_PER(gen_latent_fixed_comm_fixed_int_periods,2,1,0,1,0);
GEN_PER(gen_latent_fixed_comm_int_periods,2,1,0,2,0);

GEN_PER(gen_latent_comm_periods,2,2,0,0,0);
GEN_PER(gen_latent_comm_fixed_int_periods,2,2,0,1,0);
GEN_PER(gen_latent_comm_int_periods,2,2,0,2,0);


GEN_PER(gen_fixed_latent_fixed_comm_fixed_alternate_periods,1,1,1,0,0);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_alternate_fixed_int_periods,1,1,1,0,1);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_alternate_int_periods,1,1,1,0,2);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_fixed_alternate_periods,1,1,1,1,0);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_fixed_alternate_fixed_int_periods,1,1,1,1,1);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_fixed_alternate_int_periods,1,1,1,1,2);
GEN_PER(gen_fixed_latent_fixed_comm_int_fixed_alternate_periods,1,1,1,2,0);
GEN_PER(gen_fixed_latent_fixed_comm_int_fixed_alternate_fixed_int_periods,1,1,1,2,1);
GEN_PER(gen_fixed_latent_fixed_comm_int_fixed_alternate_int_periods,1,1,1,2,2);

GEN_PER(gen_fixed_latent_fixed_comm_alternate_periods,1,1,2,0,0);
GEN_PER(gen_fixed_latent_fixed_comm_alternate_fixed_int_periods,1,1,2,0,1);
GEN_PER(gen_fixed_latent_fixed_comm_alternate_int_periods,1,1,2,0,2);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_alternate_periods,1,1,2,1,0);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_alternate_fixed_int_periods,1,1,2,1,1);
GEN_PER(gen_fixed_latent_fixed_comm_fixed_int_alternate_int_periods,1,1,2,1,2);
GEN_PER(gen_fixed_latent_fixed_comm_int_alternate_periods,1,1,2,2,0);
GEN_PER(gen_fixed_latent_fixed_comm_int_alternate_fixed_int_periods,1,1,2,2,1);
GEN_PER(gen_fixed_latent_fixed_comm_int_alternate_int_periods,1,1,2,2,2);

GEN_PER(gen_fixed_latent_comm_fixed_alternate_periods,1,2,1,0,0);
GEN_PER(gen_fixed_latent_comm_fixed_alternate_fixed_int_periods,1,2,1,0,1);
GEN_PER(gen_fixed_latent_comm_fixed_alternate_int_periods,1,2,1,0,2);
GEN_PER(gen_fixed_latent_comm_fixed_int_fixed_alternate_periods,1,2,1,1,0);
GEN_PER(gen_fixed_latent_comm_fixed_int_fixed_alternate_fixed_int_periods,1,2,1,1,1);
GEN_PER(gen_fixed_latent_comm_fixed_int_fixed_alternate_int_periods,1,2,1,1,2);
GEN_PER(gen_fixed_latent_comm_int_fixed_alternate_periods,1,2,1,2,0);
GEN_PER(gen_fixed_latent_comm_int_fixed_alternate_fixed_int_periods,1,2,1,2,1);
GEN_PER(gen_fixed_latent_comm_int_fixed_alternate_int_periods,1,2,1,2,2);

GEN_PER(gen_fixed_latent_comm_alternate_periods,1,2,2,0,0);
GEN_PER(gen_fixed_latent_comm_alternate_fixed_int_periods,1,2,2,0,1);
GEN_PER(gen_fixed_latent_comm_alternate_int_periods,1,2,2,0,2);
GEN_PER(gen_fixed_latent_comm_fixed_int_alternate_periods,1,2,2,1,0);
GEN_PER(gen_fixed_latent_comm_fixed_int_alternate_fixed_int_periods,1,2,2,1,1);
GEN_PER(gen_fixed_latent_comm_fixed_int_alternate_int_periods,1,2,2,1,2);
GEN_PER(gen_fixed_latent_comm_int_alternate_periods,1,2,2,2,0);
GEN_PER(gen_fixed_latent_comm_int_alternate_fixed_int_periods,1,2,2,2,1);
GEN_PER(gen_fixed_latent_comm_int_alternate_int_periods,1,2,2,2,2);

GEN_PER(gen_latent_fixed_comm_fixed_alternate_periods,2,1,1,0,0);
GEN_PER(gen_latent_fixed_comm_fixed_alternate_fixed_int_periods,2,1,1,0,1);
GEN_PER(gen_latent_fixed_comm_fixed_alternate_int_periods,2,1,1,0,2);
GEN_PER(gen_latent_fixed_comm_fixed_int_fixed_alternate_periods,2,1,1,1,0);
GEN_PER(gen_latent_fixed_comm_fixed_int_fixed_alternate_fixed_int_periods,2,1,1,1,1);
GEN_PER(gen_latent_fixed_comm_fixed_int_fixed_alternate_int_periods,2,1,1,1,2);
GEN_PER(gen_latent_fixed_comm_int_fixed_alternate_periods,2,1,1,2,0);
GEN_PER(gen_latent_fixed_comm_int_fixed_alternate_fixed_int_periods,2,1,1,2,1);
GEN_PER(gen_latent_fixed_comm_int_fixed_alternate_int_periods,2,1,1,2,2);

GEN_PER(gen_latent_fixed_comm_alternate_periods,2,1,2,0,0);
GEN_PER(gen_latent_fixed_comm_alternate_fixed_int_periods,2,1,2,0,1);
GEN_PER(gen_latent_fixed_comm_alternate_int_periods,2,1,2,0,2);
GEN_PER(gen_latent_fixed_comm_fixed_int_alternate_periods,2,1,2,1,0);
GEN_PER(gen_latent_fixed_comm_fixed_int_alternate_fixed_int_periods,2,1,2,1,1);
GEN_PER(gen_latent_fixed_comm_fixed_int_alternate_int_periods,2,1,2,1,2);
GEN_PER(gen_latent_fixed_comm_int_alternate_periods,2,1,2,2,0);
GEN_PER(gen_latent_fixed_comm_int_alternate_fixed_int_periods,2,1,2,2,1);
GEN_PER(gen_latent_fixed_comm_int_alternate_int_periods,2,1,2,2,2);

GEN_PER(gen_latent_comm_fixed_alternate_periods,2,2,1,0,0);
GEN_PER(gen_latent_comm_fixed_alternate_fixed_int_periods,2,2,1,0,1);
GEN_PER(gen_latent_comm_fixed_alternate_int_periods,2,2,1,0,2);
GEN_PER(gen_latent_comm_fixed_int_fixed_alternate_periods,2,2,1,1,0);
GEN_PER(gen_latent_comm_fixed_int_fixed_alternate_fixed_int_periods,2,2,1,1,1);
GEN_PER(gen_latent_comm_fixed_int_fixed_alternate_int_periods,2,2,1,1,2);
GEN_PER(gen_latent_comm_int_fixed_alternate_periods,2,2,1,2,0);
GEN_PER(gen_latent_comm_int_fixed_alternate_fixed_int_periods,2,2,1,2,1);
GEN_PER(gen_latent_comm_int_fixed_alternate_int_periods,2,2,1,2,2);

GEN_PER(gen_latent_comm_alternate_periods,2,2,2,0,0);
GEN_PER(gen_latent_comm_alternate_fixed_int_periods,2,2,2,0,1);
GEN_PER(gen_latent_comm_alternate_int_periods,2,2,2,0,2);
GEN_PER(gen_latent_comm_fixed_int_alternate_periods,2,2,2,1,0);
GEN_PER(gen_latent_comm_fixed_int_alternate_fixed_int_periods,2,2,2,1,1);
GEN_PER(gen_latent_comm_fixed_int_alternate_int_periods,2,2,2,1,2);
GEN_PER(gen_latent_comm_int_alternate_periods,2,2,2,2,0);
GEN_PER(gen_latent_comm_int_alternate_fixed_int_periods,2,2,2,2,1);
GEN_PER(gen_latent_comm_int_alternate_int_periods,2,2,2,2,2);
//! @endcond

/*static inline void gen_time_periods_isolation(sim_vars* sv)
{
  double m;

  sv->ii->latent_period=gsl_ran_gamma(sv->r, sv->pars.kappal*sv->pars.lbar, 1./sv->pars.kappal);
  sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa);

  if(gsl_rng_uniform(sv->r) >= sv->pars.q) {

    if(isinf(sv->pars.kappaq)) m=sv->pars.mbar;
    else m=gsl_ran_gamma(sv->r, sv->pars.kappaq*sv->pars.mbar, 1./sv->pars.kappaq);

    if(m<sv->ii->comm_period) sv->ii->comm_period=m;
  }
  double time_left=sv->pars.tmax-(sv->ii-1)->event_time;

  if(sv->ii->comm_period > time_left) {
    //sv->ii->comm_period=time_left;
    sv->ii->infectious_at_tmax=true;

  } else sv->ii->infectious_at_tmax=false;
  DEBUG_PRINTF("Comm period is %f%s\n",sv->ii->comm_period,(sv->ii->infectious_at_tmax?" (reached end)":""));
}
*/

/**
 * @brief Default processing function that is called when a new transmission event is created.
 *
 * This function is called by default if a user-defined function has not been
 * set.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event does not occur after tmax, false otherwise.
 */
inline static bool default_event_proc_func(sim_vars* sv){return (sv->ii->event_time <= sv->pars.tmax);}

/**
 * @brief Default processing function that is called when the maximum number of
 * layers is increased.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 *
 * @param iis: First newly allocated layer element of the infectious individuals array.
 * @param n: Number of layers that have been added
 */
inline static void default_increase_layers_proc_func(infindividual* iis, uint32_t n){}

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
