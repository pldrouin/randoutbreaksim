/**
 * @file infindividual.h
 * @brief Data structure containing infected individual information.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 * Model from <jerome.levesque@tpsgc-pwgsc.gc.ca> and
 * <david.maybury@tpsgc-pwgsc.gc.ca>
*/

#ifndef _INFINDIVIDUAL_
#define _INFINDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

/**
 * Type of communicable period type for a given individual.
 **/
enum ro_commper_type {ro_commper_main=1, ro_commper_alt=2, ro_commper_int=4, ro_commper_true_positive_test=8, ro_commper_main_int=5, ro_commper_alt_int=6, ro_commper_tmax=16};

/**
 * Infected individual
 */
typedef struct
{
    void* dataptr;	      //!< Data pointer for user-defined functions
    double latent_period;     //!< Latent period
    double comm_period;       //!< Communicable period
    double end_comm_period;   //!< End of communicable period
    double event_time;	      //!< Start time for the current iteration event
    uint32_t nevents;	      //!< Number of events
    uint32_t curevent;	      //!< Index of the current iteration event
    uint32_t nattendees;      //!< Number of attendees for the current iteration event
    uint32_t ninfections;     //!< Number of infections for the current iteration event
    uint32_t curinfection;    //!< Index of the current iteration infection
    uint8_t commpertype;      //!< Flags for the type of communicable period (filled using ro_commpertype flags)
} infindividual;

#endif
