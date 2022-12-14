/**
 * @file infindividual.h
 * @brief Data structure containing infected individual information.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
*/

#ifndef _INFINDIVIDUAL_
#define _INFINDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

struct sim_vars_;

/**
 * Type of communicable period type for a given individual.
 **/
enum ro_commper_type {ro_commper_main=1, ro_commper_alt=2, ro_commper_int=4, ro_commper_true_positive_test=8, ro_commper_main_int=5, ro_commper_alt_int=6, ro_commper_int_true_positive_test=12};

/**
 * Infected individual
 */
typedef struct infindividual_
{
    void* dataptr;	      //!< Data pointer for user-defined functions.
    double latent_period;     //!< Latent period
    double comm_period;       //!< Communicable period
#ifdef CT_OUTPUT
    double presym_comm_period; //!< Pre-symptomatic communicable period for symptomatic individuals
    void (*gen_ct_time_periods_func)(struct sim_vars_*, struct infindividual_* ii, struct infindividual_* iiparent, const double inf_start);	//!< Effective time period function used for the current event
#endif
    double end_comm_period;   //!< End of communicable period
#ifdef DUAL_PINF
    double q;                 //!< Probability of alternate communicable period for the infectious individual. Depends on the infection category
    double pinf;              //!< Probability of infection from the infectious individual. Depends on the infection category
#endif
    uint32_t generation;      //!< Infection generation
    uint8_t commpertype;      //!< Flags for the type of communicable period (filled using ro_commper_type flags)
#ifdef CT_OUTPUT
    bool traced;	      //!> Flag that indicates if this individual has been successfully traced (traced!=interrupted)
#endif
#ifdef DUAL_PINF
    bool inftypep;            //!> Flag that indicates if the infectious individual is in the second category
#endif
  uint32_t nevents;	      //!< Number of events
  uint32_t nattendees;        //!< Number of attendees for the current iteration event
#ifdef CT_OUTPUT
  uint32_t ntracednicts;    //!< Number of successfully traced non-infected contacts for the current iteration event
  uint32_t ntracedicts;     //!< Number of successfully traced infected contacts for the current iteration event
#endif
  uint32_t ninfections;     //!< Number of infections for the current iteration event
#ifdef DUAL_PINF
  uint32_t ninfectionsf;    //<! Number of remaining infections in the first category for the current iteration event. 
  uint32_t ninfectionsp;    //<! Number of remaining infections in the second category for the current iteration event. 
#endif
} infindividual;

#endif
