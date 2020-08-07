/**
 * @file simulation.h
 * @brief Simulation data structures and functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
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

#define INIT_N_LAYERS (16) //!< Initial number of simulation layers

#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) printf(__VA_ARGS__) //!< Debug print function

/**
 * Simulation input parameters.
 */
typedef struct 
{
  double tbar;		//!< Mean uninterrupted communicable period
  double p;		//!< Parameter for the logarithmic distribution used to draw number of new infections for a given transmission event
  double lambda;	//!< Rate of transmission events
  double kappa;		//!< branchim's kappa parameter for the gamma distribution used to generate the uninterrupted communicable period
  double q;		//!< Probability of alternate communicable period
  double mbar;		//!< Mean period for the alternate communicable period
  double kappaq;	//!< branchim's kappa parameter for the gamma distribution used to generate the alternate communicable period
  double tmax;		//!< Maximum simulation period used to instantiate new infectious individuals.
  uint32_t nstart;	//!< Initial number of infectious individuals
} sim_pars;

/**
 * Simulation variables
 */
typedef struct sim_vars_
{
  sim_pars pars;		//!< Simulation input parameters
  gsl_rng const* r;		//!< Pointer to GSL random number generator
  infindividual* iis;	//!< Array of current infectious individuals across all layers
  infindividual* ii;	//!< Pointer to current iteration infectious individual
  uint32_t nlayers;		//!< Current maximum number of layers that has been used so far 
  void* dataptr;		//!< Simulation-level data pointer for user-defined functions
  void (*gen_comm_period_func)(struct sim_vars_*);				//!< Pointer to the function used to generate a communicable period for a given infectious individual
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
int sim_pars_check(sim_pars const* pars);

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
int sim_init(sim_vars* sv, sim_pars* pars, const gsl_rng* r);

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

/**
 * @brief Generates an uninterrupted communicable period.
 *
 * This function generates the uninterrupted communicable period for the current
 * infectious individual, as determined by the simulation variables and
 * parameters. It also sets the infectious_at_tmax parameter for the current
 * infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 */
static inline void gen_comm_period(sim_vars* sv)
{
  sv->ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.kappa*sv->pars.tbar, 1./sv->pars.kappa);
  double time_left=sv->pars.tmax-(sv->ii-1)->event_time;

  if(sv->ii->comm_period > time_left) {
    //sv->ii->comm_period=time_left;
    sv->ii->infectious_at_tmax=true;

  } else sv->ii->infectious_at_tmax=false;
  DEBUG_PRINTF("Comm period is %f%s\n",sv->ii->comm_period,(sv->ii->infectious_at_tmax?" (reached end)":""));
}

/**
 * @brief Generates a communicable period.
 *
 * This function generates the communicable period for the current
 * infectious individual, as determined by the simulation variables and
 * parameters. It considers the alternate communicable period based on the q
 * simulation parameter, and selects the shortest period between the
 * uninterrupted period and the alternate period. It also sets the
 * infectious_at_tmax parameter for the current infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 */
static inline void gen_comm_period_isolation(sim_vars* sv)
{
  double m;

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
