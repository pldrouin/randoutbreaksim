/**
 * @file individual.h
 * @brief Data structure containing individual information for finitepopsim.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
*/

#ifndef _INDIVIDUAL_
#define _INDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "infindividual.h"

/**
 * Type of infection status for a given individual.
 **/
enum ro_ind_inf_status{ro_ind_inf_status_susceptible=0, ro_ind_inf_status_latent, ro_ind_inf_status_infectious, ro_ind_inf_status_hospitalised, ro_ind_inf_status_recovered, ro_ind_inf_status_dead};

/**
 * Type of infection status for a given individual.
 **/
enum ro_ind_interact_status{ro_ind_interact_status_interacting=0, ro_ind_interact_status_isolated};

/**
 * Infected individual
 */
typedef struct individual_
{
  struct individual_* parent;   //!< Last who has infected this individual => *********** Careful with this!
  infindividual ii;		//!< Infectious individual properties for this individual
  uint8_t indinfstatus;         //!< Flags for the infection status of this individual
  uint8_t indinteractstatus;    //!< Flags for the interaction status of this individual
} individual;

#endif
