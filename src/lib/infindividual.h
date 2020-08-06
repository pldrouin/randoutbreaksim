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
 * Infected individual
 */
typedef struct
{
    void* dataptr;	      //!< Data pointer for user-defined functions
    double comm_period;       //!< Communicable period
    double event_time;	      //!< Start time for the current iteration event
    uint32_t nevents;	      //!< Number of events
    uint32_t curevent;	      //!< Index of the current iteration event
    uint32_t ninfections;     //!< Number of infections for the current iteration event
    uint32_t curinfection;    //!< Index of the current iteration infection
    bool infectious_at_tmax;  //!< Is the infected individual still infectious at tmax?
} infindividual;

#endif
