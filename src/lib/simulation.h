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
#include "individual.h"
#include "model_parameters.h"
#include "simple_array.h"

#include "ran_log.h"

extern int __ro_debug;
#define DEBUG_PRINTF(...) //!< Debug print function
//#define DEBUG_PRINTF(...) {if(__ro_debug) printf(__VA_ARGS__);} //!< Debug print function

typedef struct {
  infindividual* iis;	//!< Array of current infectious individuals across all layers
  uint32_t nlayers;	//!< Current maximum number of layers that has been used so far 
  uint32_t naevents;    //!< Number of allocated events for each layer
} brsim_vars;

typedef struct {
  individual rooti;     //!< Root individual for the simulation
  individual* is;       //!< All the individuals in the simulation
} fpsim_vars;

/**
 * Simulation variables
 */
typedef struct sim_vars_
{
  model_pars pars;		//!< Simulation input parameters
  gsl_rng const* r;		//!< Pointer to GSL random number generator
  infindividual* curii;		//!< Pointer to current iteration infectious individual
  double event_time;	        //!< Start time for the current iteration event
  void* dataptr;		//!< Simulation-level data pointer for user-defined functions
  void (*gen_time_origin_func)(struct sim_vars_*, infindividual* ii);	//!<Pointer to the function used to apply a time shift
  void (*gen_pri_time_periods_func)(struct sim_vars_*, infindividual* ii, infindividual* iiparent, const double inf_start);	//!< Pointer to the function used to generate time periods for a given primary infectious individual
  void (*gen_time_periods_func)(struct sim_vars_*, infindividual* ii, infindividual* iiparent, const double inf_start);	//!< Pointer to the function used to generate time periods for a given infectious individual
  void (*gen_time_periods_func_no_int)(struct sim_vars_*, infindividual* ii, infindividual* iiparent, const double inf_start);	//!< Pointer to the function used to generate time periods without interruption for a given infectious individual
  uint32_t (*gen_att_func)(struct sim_vars_*);				        //!< Pointer to the function used to generate attendees during one event
  void (*gen_att_inf_func)(struct sim_vars_*);              	        //!< Pointer to the function used to generate attendees and new infections during one event
  void (*path_init_proc_func)(struct sim_vars_*);	                //!< Pointer to the user-defined path initialisation function.
  bool (*path_end_proc_func)(struct sim_vars_*);	                //!< Pointer to the user-defined path termination function. The returned value from this function determines if the simulated path is to be included in the simulation.
  void (*pri_init_proc_func)(struct sim_vars_*, infindividual* ii);	//!< Pointer to the user-defined initialisation function for a given primary infectious individual
  void (*ii_alloc_proc_func)(infindividual* ii);	//!< Pointer to the user-defined processing function that is called when memory for a new infectious individual is allocated.
  bool (*new_event_proc_func)(struct sim_vars_* sv);				//!< Pointer to the user-defined processing function that is called when a new transmission event is created, after an event time and the number of new infections have been assigned. The returned value from this function determines if new infectious individuals are instantiated for this event. The function can also be called in CT_OUTPUT mode for contact events during the latent phase of the individual.
  void (*new_inf_proc_func)(struct sim_vars_* sv, infindividual* ii, infindividual* parent);			//!< Pointer to the user-defined processing function that is called when a new infected individual is created, after the communicable period and the number of transmission events have been assigned. The function is only called if the number of transmission events is non-zero. 
  void (*new_inf_proc_func_noevent)(struct sim_vars_* sv, infindividual* ii, infindividual* parent);	//!< Pointer to the user-defined processing function that is called for a new infected individual that does not generate any transmission event.
  void (*end_inf_proc_func)(struct sim_vars_* sv, infindividual* ii, infindividual* parent); 		//!< Pointer to the user-defined processing function that is called once all transmission events for a given infectious individual have been generated.
  ran_log rl;	//!< Handle for the logarithmic random variate generator.

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
 * have been assigned.
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
inline static void sim_set_new_inf_proc_func(sim_vars* sv, void (*new_inf_proc_func)(sim_vars* sv, infindividual* ii, infindividual* parent)){sv->new_inf_proc_func=new_inf_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called when a new
 * path is simulated.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void sim_set_path_init_proc_func(sim_vars* sv, void (*path_init_proc_func)(sim_vars* sv)){sv->path_init_proc_func=path_init_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called when a simulated path has ended. The returned value
 * from this function determines if the simulated path is to be included in the simulation.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void sim_set_path_end_proc_func(sim_vars* sv, bool (*path_end_proc_func)(sim_vars* sv)){sv->path_end_proc_func=path_end_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called when a new
 * primary individual is created.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void sim_set_pri_init_proc_func(sim_vars* sv, void (*pri_init_proc_func)(sim_vars* sv, infindividual* ii)){sv->pri_init_proc_func=pri_init_proc_func;}

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
inline static void sim_set_end_inf_proc_func(sim_vars* sv, void (*end_inf_proc_func)(sim_vars* sv, infindividual* ii, infindividual* parent)){sv->end_inf_proc_func=end_inf_proc_func;}

/**
 * @brief Sets the user-defined processing function that is called for an infectious individual that does not generate any transmission event.
 *
 * This function is called after the communicable period has been assigned for
 * the current infectious individual.
 *
 * @param sv: Pointer to the simulation variables.
 * @param new_inf_proc_func_noevent: Pointer to the user-defined function. The second
 * argument for this function is the simulation-level data pointer.
 */
inline static void sim_set_new_inf_proc_noevent_func(sim_vars* sv, void (*new_inf_proc_func_noevent)(struct sim_vars_* sv, infindividual* ii, infindividual* parent)){sv->new_inf_proc_func_noevent=new_inf_proc_func_noevent;}

/**
 * @brief Function that modifies the simulation time to use the creation time of
 * a primary individual as the time origin.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void gen_time_origin_pri_created(sim_vars* sv){}

/**
 * @brief Function that modifies the simulation time to use the time of
 * a primary individual becomes infectious as the time origin.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void gen_time_origin_pri_infectious(sim_vars* sv, infindividual* ii){sv->event_time=-ii->end_comm_period+ii->comm_period; ii->end_comm_period=ii->comm_period;}

/**
 * @brief Function that modifies the simulation time to use the end of the
 * communicable period for a primary individual as the time origin.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void gen_time_origin_pri_end_comm(sim_vars* sv, infindividual* ii){sv->event_time=-ii->end_comm_period; ii->end_comm_period=0;}

/**
 * @brief Function that modifies the simulation time to time of test results for a
 * primary individual as the time origin.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void gen_time_origin_pri_test_results(sim_vars* sv, infindividual* ii){sv->event_time=-ii->end_comm_period-sv->pars.tdeltat; ii->end_comm_period=-sv->pars.tdeltat;}

/**
 * @brief Function that modifies the simulation such that t=0 is uniformly
 * distributed within the communicable period and the primary individual is
 * assumed to enter the simulation at t=0, so no infection occurs prior to that.
 *
 * @param sv: Pointer to the simulation variables.
 */
inline static void gen_time_origin_pri_flat_comm(sim_vars* sv, infindividual* ii){
  ii->latent_period=0;
  ii->end_comm_period=ii->comm_period*=gsl_rng_uniform(sv->r);
#ifdef CT_OUTPUT
  if(ii->commpertype&&ro_commper_alt) ii->presym_comm_period=ii->comm_period;
#endif
}

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main GEN_PER macro.
 */
#define GEN_PER_LATENT_0 ii->latent_period=0;
#define GEN_PER_LATENT_1 ii->latent_period=sv->pars.lbar;
#define GEN_PER_LATENT_2 ii->latent_period=gsl_ran_gamma(sv->r, sv->pars.la, sv->pars.lb);

#define GEN_PER_INTERRUPTED_MAIN_0

#ifdef CT_OUTPUT
  #define GEN_PER_INTERRUPTED_MAIN_1_PRE if(ii->traced && gsl_rng_uniform(sv->r) < sv->pars.pitnet) { const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + sv->pars.itbar;
#else
  #define GEN_PER_INTERRUPTED_MAIN_1_PRE if(iiparent->commpertype&ro_commper_true_positive_test && gsl_rng_uniform(sv->r) < sv->pars.pit) { const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + sv->pars.itbar;
#endif
#define GEN_PER_INTERRUPTED_MAIN_POST if (ecp < ii->end_comm_period) { \
    ii->comm_period=ecp-(inf_start+ii->latent_period); \
    \
    if(ii->comm_period<0) { \
      ii->latent_period+=ii->comm_period; \
      ii->comm_period=0; \
      ii->end_comm_period=inf_start+ii->latent_period; \
      \
    } else ii->end_comm_period=ecp; \
    \
    if(gsl_rng_uniform(sv->r) < sv->pars.ttpr) ii->commpertype|=ro_commper_int|ro_commper_true_positive_test; \
    \
    else ii->commpertype|=ro_commper_int; \
  } else ii->commpertype|=ro_commper_int; \
}
#define GEN_PER_INTERRUPTED_MAIN_1 GEN_PER_INTERRUPTED_MAIN_1_PRE GEN_PER_INTERRUPTED_MAIN_POST

#ifdef CT_OUTPUT
#define GEN_PER_INTERRUPTED_MAIN_2_PRE if(ii->traced && gsl_rng_uniform(sv->r) < sv->pars.pitnet) { \
  const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + gsl_ran_gamma(sv->r, sv->pars.ita, sv->pars.itb);
#else
#define GEN_PER_INTERRUPTED_MAIN_2_PRE if(iiparent->commpertype&ro_commper_true_positive_test && gsl_rng_uniform(sv->r) < sv->pars.pit) { \
  const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + gsl_ran_gamma(sv->r, sv->pars.ita, sv->pars.itb);
#endif
#define GEN_PER_INTERRUPTED_MAIN_2 GEN_PER_INTERRUPTED_MAIN_2_PRE GEN_PER_INTERRUPTED_MAIN_POST

#define GEN_TESTED_ALT_0 ii->commpertype=ro_commper_alt;
#define GEN_TESTED_ALT_1 ii->commpertype=ro_commper_alt|ro_commper_true_positive_test;
#define GEN_TESTED_ALT_2 if(gsl_rng_uniform(sv->r) < sv->pars.mtpr) ii->commpertype=ro_commper_alt|ro_commper_true_positive_test; else ii->commpertype=ro_commper_alt;

#define GEN_PER_MAIN_1(IT) {ii->comm_period=sv->pars.tbar; ii->end_comm_period=inf_start+ii->latent_period+ii->comm_period; ii->commpertype=ro_commper_main; GEN_PER_INTERRUPTED_MAIN_ ## IT;}
#define GEN_PER_MAIN_2(IT) {ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.ta, sv->pars.tb); ii->end_comm_period=inf_start+ii->latent_period+ii->comm_period; ii->commpertype=ro_commper_main; GEN_PER_INTERRUPTED_MAIN_ ## IT;}

#define GEN_PER_INTERRUPTED_ALT_0

#ifdef CT_OUTPUT
#define GEN_PER_INTERRUPTED_ALT_1_PRE if(ii->traced && gsl_rng_uniform(sv->r) < sv->pars.pimnet) { const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + sv->pars.imbar;
#else
#define GEN_PER_INTERRUPTED_ALT_1_PRE if(iiparent->commpertype&ro_commper_true_positive_test && gsl_rng_uniform(sv->r) < sv->pars.pim) { const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + sv->pars.imbar;
#endif
#define GEN_PER_INTERRUPTED_ALT_POST if (ecp < ii->end_comm_period) { \
    ii->comm_period=ecp-(inf_start+ii->latent_period); \
    \
    if(ii->comm_period<0) { \
      ii->latent_period+=ii->comm_period; \
      ii->comm_period=0; \
      ii->end_comm_period=inf_start+ii->latent_period; \
      \
    } else ii->end_comm_period=ecp; \
    \
    if(gsl_rng_uniform(sv->r) < sv->pars.mtpr) ii->commpertype|=ro_commper_int|ro_commper_true_positive_test; \
    \
    else ii->commpertype|=ro_commper_int; \
  } else ii->commpertype|=ro_commper_int; \
}
#define GEN_PER_INTERRUPTED_ALT_1 GEN_PER_INTERRUPTED_ALT_1_PRE GEN_PER_INTERRUPTED_ALT_POST

#ifdef CT_OUTPUT
#define GEN_PER_INTERRUPTED_ALT_2_PRE if(ii->traced && gsl_rng_uniform(sv->r) < sv->pars.pimnet) { \
  const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + gsl_ran_gamma(sv->r, sv->pars.ima, sv->pars.imb);
#else
#define GEN_PER_INTERRUPTED_ALT_2_PRE if(iiparent->commpertype&ro_commper_true_positive_test && gsl_rng_uniform(sv->r) < sv->pars.pim) { \
  const double ecp=iiparent->end_comm_period + sv->pars.tdeltat + gsl_ran_gamma(sv->r, sv->pars.ima, sv->pars.imb);
#endif
#define GEN_PER_INTERRUPTED_ALT_2 GEN_PER_INTERRUPTED_ALT_2_PRE GEN_PER_INTERRUPTED_ALT_POST

#define GEN_PER_ALTERNATE_ONLY_1_PRE ii->comm_period=sv->pars.mbar;
#define GEN_PER_ALTERNATE_ONLY_2_PRE ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.ma, sv->pars.mb);
#define GEN_PER_ALTERNATE_ONLY_POST(IM,TESTING) ii->end_comm_period=inf_start+ii->latent_period+ii->comm_period; GEN_TESTED_ALT_ ## TESTING GEN_PER_INTERRUPTED_ALT_ ## IM;

#ifdef DUAL_PINF
#define GEN_PER_ALTERNATE_PRE if(gsl_rng_uniform(sv->r) < ii->q) {
#else
#define GEN_PER_ALTERNATE_PRE if(gsl_rng_uniform(sv->r) < sv->pars.q) {
#endif
#define GEN_PER_ALTERNATE_1_POST(MAIN,IT,IM,TESTING) ii->comm_period=sv->pars.mbar; ii->end_comm_period=inf_start+ii->latent_period+ii->comm_period; GEN_TESTED_ALT_ ## TESTING GEN_PER_INTERRUPTED_ALT_ ## IM} else {GEN_PER_MAIN_ ## MAIN(IT);}
#define GEN_PER_ALTERNATE_2_POST(MAIN,IT,IM,TESTING) ii->comm_period=gsl_ran_gamma(sv->r, sv->pars.ma, sv->pars.mb); ii->end_comm_period=inf_start+ii->latent_period+ii->comm_period; GEN_TESTED_ALT_ ## TESTING GEN_PER_INTERRUPTED_ALT_ ## IM} else {GEN_PER_MAIN_ ## MAIN(IT);}

#ifdef CT_OUTPUT
#define GEN_PER_ALTERNATE_ONLY_1(IM,TESTING) ii->presym_comm_period=GEN_PER_ALTERNATE_ONLY_1_PRE GEN_PER_ALTERNATE_ONLY_POST(IM,TESTING)
#define GEN_PER_ALTERNATE_ONLY_2(IM,TESTING) ii->presym_comm_period=GEN_PER_ALTERNATE_ONLY_2_PRE GEN_PER_ALTERNATE_ONLY_POST(IM,TESTING)
#define GEN_PER_ALTERNATE_1(MAIN,IT,IM,TESTING) GEN_PER_ALTERNATE_PRE ii->presym_comm_period=GEN_PER_ALTERNATE_1_POST(MAIN,IT,IM,TESTING)
#define GEN_PER_ALTERNATE_2(MAIN,IT,IM,TESTING) GEN_PER_ALTERNATE_PRE ii->presym_comm_period=GEN_PER_ALTERNATE_2_POST(MAIN,IT,IM,TESTING)
#else
#define GEN_PER_ALTERNATE_ONLY_1(IM,TESTING) GEN_PER_ALTERNATE_ONLY_1_PRE GEN_PER_ALTERNATE_ONLY_POST(IM,TESTING)
#define GEN_PER_ALTERNATE_ONLY_2(IM,TESTING) GEN_PER_ALTERNATE_ONLY_2_PRE GEN_PER_ALTERNATE_ONLY_POST(IM,TESTING)
#define GEN_PER_ALTERNATE_1(MAIN,IT,IM,TESTING) GEN_PER_ALTERNATE_PRE GEN_PER_ALTERNATE_1_POST(MAIN,IT,IM,TESTING)
#define GEN_PER_ALTERNATE_2(MAIN,IT,IM,TESTING) GEN_PER_ALTERNATE_PRE GEN_PER_ALTERNATE_2_POST(MAIN,IT,IM,TESTING)
#endif

#define GEN_PER_MAIN_ALTERNATE_0_1(IT,IM,TESTING) GEN_PER_ALTERNATE_ONLY_1(IM,TESTING);
#define GEN_PER_MAIN_ALTERNATE_0_2(IT,IM,TESTING) GEN_PER_ALTERNATE_ONLY_2(IM,TESTING);
#define GEN_PER_MAIN_ALTERNATE_1_0(IT,IM,TESTING) GEN_PER_MAIN_1(IT);
#define GEN_PER_MAIN_ALTERNATE_2_0(IT,IM,TESTING) GEN_PER_MAIN_2(IT);
#define GEN_PER_MAIN_ALTERNATE_1_1(IT,IM,TESTING)  GEN_PER_ALTERNATE_1(1,IT,IM,TESTING);
#define GEN_PER_MAIN_ALTERNATE_2_1(IT,IM,TESTING)  GEN_PER_ALTERNATE_1(2,IT,IM,TESTING);
#define GEN_PER_MAIN_ALTERNATE_1_2(IT,IM,TESTING)  GEN_PER_ALTERNATE_2(1,IT,IM,TESTING);
#define GEN_PER_MAIN_ALTERNATE_2_2(IT,IM,TESTING)  GEN_PER_ALTERNATE_2(2,IT,IM,TESTING);
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
 * @param TESTING: Indicates if testing is performed (0=no, 1=yes)
 * latent)
 */
#define GEN_PER(LATENT,MAIN,IT,ALTERNATE,IM,TESTING) static inline void gen_comm_ ## LATENT ## _ ## MAIN ## _ ## IT ## _ ## ALTERNATE ## _ ## IM ## _ ## TESTING ## _periods(sim_vars* sv, infindividual* ii, infindividual* iiparent, const double inf_start) \
{ \
  DEBUG_PRINTF("Function is %s\n",__func__); \
  GEN_PER_LATENT_ ## LATENT \
  GEN_PER_MAIN_ALTERNATE_ ## MAIN ## _ ## ALTERNATE(IT,IM,TESTING) \
}

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main GEN_PERS_MAIN macro.
 */
#define GEN_PERS_IM(LATENT,MAIN,IT,ALTERNATE,TESTING) GEN_PER(LATENT,MAIN,IT,ALTERNATE,0,TESTING) GEN_PER(LATENT,MAIN,IT,ALTERNATE,1,TESTING) GEN_PER(LATENT,MAIN,IT,ALTERNATE,2,TESTING)
#define GEN_PERS_ALT(LATENT,MAIN,IT,TESTING) GEN_PER(LATENT,MAIN,IT,0,0,TESTING) GEN_PERS_IM(LATENT,MAIN,IT,1,TESTING) GEN_PERS_IM(LATENT,MAIN,IT,2,TESTING)
#define GEN_PERS_IT_0(LATENT,TESTING) GEN_PERS_IM(LATENT,0,0,1,TESTING) GEN_PERS_IM(LATENT,0,0,2,TESTING)
#define GEN_PERS_IT(LATENT,MAIN,TESTING) GEN_PERS_ALT(LATENT,MAIN,0,TESTING) GEN_PERS_ALT(LATENT,MAIN,1,TESTING) GEN_PERS_ALT(LATENT,MAIN,2,TESTING)
//! @endcond

/**
 * @brief For a given type of latent period, generates functions that each
 * generates a specific combination of communicable and latent periods.
 *
 * @param LATENT: Type of latent period (0=none, 1=fixed, 2=variable)
 * @param TESTING: Indicates if testing is performed (0=no, 1=yes, true positive
 * rate of 1, 2=yes, true positive rate < 1)
 */
#define GEN_PERS_MAIN(LATENT,TESTING) GEN_PERS_IT_0(LATENT,TESTING) GEN_PERS_IT(LATENT,1,TESTING) GEN_PERS_IT(LATENT,2,TESTING)
GEN_PERS_MAIN(0,0)
GEN_PERS_MAIN(1,0)
GEN_PERS_MAIN(2,0)
GEN_PERS_MAIN(0,1)
GEN_PERS_MAIN(1,1)
GEN_PERS_MAIN(2,1)
GEN_PERS_MAIN(0,2)
GEN_PERS_MAIN(1,2)
GEN_PERS_MAIN(2,2)

//! @cond Doxygen_Suppress
/**
 * The preprocessing macros below are used by the main PER_COND macro.
 */
#define PRI_PER_COND_FULL(TESTING,LATENT,MAIN,ALT) if(sv->pars.pricommpertype&ro_pricommper_alt_use_tpr) sv->gen_pri_time_periods_func=gen_comm_ ## LATENT ## _ ## MAIN ## _0_ ## ALT ## _0_ ## TESTING ## _periods; else sv->gen_pri_time_periods_func=gen_comm_ ## LATENT ## _ ## MAIN ## _0_ ## ALT ## _0_1_periods; 
#define PER_COND_WITH_ALT(TESTING,LATENT,MAIN,IT,ALT,IM)  sv->gen_time_periods_func=gen_comm_ ##  LATENT ## _ ## MAIN ## _ ## IT ## _ ## ALT ## _ ## IM ## _ ## TESTING ## _periods; sv->gen_time_periods_func_no_int=gen_comm_ ##  LATENT ## _ ## MAIN ## _0_ ## ALT ## _0_ ## TESTING ## _periods; if(sv->pars.pricommpertype&ro_pricommper_main) {if(sv->pars.pricommpertype&ro_pricommper_alt) {PRI_PER_COND_FULL(TESTING,LATENT,MAIN,ALT)} else {PRI_PER_COND_FULL(TESTING,LATENT,MAIN,0)}} else {PRI_PER_COND_FULL(TESTING,LATENT,0,ALT)}
#define PER_COND_IM(TESTING,LATENT,MAIN,IT,ALT) if(!(sv->pars.pim>0)) {PER_COND_WITH_ALT(TESTING,LATENT,MAIN,IT,ALT,0)} else if(isinf(sv->pars.kappaim)) {PER_COND_WITH_ALT(TESTING,LATENT,MAIN,IT,ALT,1)} else {PER_COND_WITH_ALT(TESTING,LATENT,MAIN,IT,ALT,2)}
#define PER_COND_ALT(TESTING,LATENT,MAIN,IT) if(!(sv->pars.q>0)) {sv->gen_time_periods_func=gen_comm_ ## LATENT ## _ ## MAIN ## _ ## IT ## _0_0_ ## TESTING ## _periods; sv->gen_time_periods_func_no_int=gen_comm_ ## LATENT ## _ ## MAIN ## _0_0_0_ ## TESTING ## _periods; sv->gen_pri_time_periods_func=gen_comm_ ## LATENT ## _ ## MAIN ## _0_0_0_ ## TESTING ## _periods;} else if(isinf(sv->pars.kappaq)) {PER_COND_IM(TESTING,LATENT,MAIN,IT,1)} else {PER_COND_IM(TESTING,LATENT,MAIN,IT,2)};
#define PER_COND_IT(TESTING,LATENT,MAIN) if(!(sv->pars.pit>0)) {PER_COND_ALT(TESTING,LATENT,MAIN,0)} else if(isinf(sv->pars.kappait)) {PER_COND_ALT(TESTING,LATENT,MAIN,1)} else {PER_COND_ALT(TESTING,LATENT,MAIN,2)};
#define PER_COND_MAIN(TESTING,LATENT) if(isinf(sv->pars.kappa)) {PER_COND_IT(TESTING,LATENT,1)} else {PER_COND_IT(TESTING,LATENT,2)};
#define PER_COND_LATENT(TESTING) if(isnan(sv->pars.kappal)) {PER_COND_MAIN(TESTING,0)} else if(isinf(sv->pars.kappal)) {PER_COND_MAIN(TESTING,1)} else {PER_COND_MAIN(TESTING,2)};
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
#define PER_COND if(isnan(sv->pars.tdeltat)) {PER_COND_LATENT(0)} else if(sv->pars.mtpr==1) {PER_COND_LATENT(1)} else {PER_COND_LATENT(2)};

/**
 * @brief Default processing function that is called when a new transmission event is created.
 *
 * This function is called by default if a user-defined function has not been
 * set.
 *
 * @param sv: Pointer to the simulation variables.
 * @return true if the event does not occur after tmax, false otherwise.
 */
inline static bool default_event_proc_func(sim_vars* sv){return (sv->event_time <= sv->pars.tmax);}

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
inline static bool dummy_proc_bool_func_sv(sim_vars* sv){return true;}

/**
 * @brief Default processing function.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 */
inline static void dummy_proc_func_sv_ii(sim_vars* sv, infindividual* ii){}

/**
 * @brief Default processing function.
 *
 * This function is called by default if a user-defined function has not been
 * set. The function does not do anything.
 */
inline static void dummy_proc_func_sv_ii2(sim_vars* sv, infindividual* ii, infindividual* ii2){}

#endif
