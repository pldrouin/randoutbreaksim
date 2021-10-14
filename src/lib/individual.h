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
enum ro_ind_inf_status{ro_ind_inf_status_latent=0, ro_ind_inf_status_infectious, ro_ind_inf_status_hospitalised, ro_ind_inf_status_recovered, ro_ind_inf_status_dead};

/**
 * Infected individual
 */
typedef struct individual_
{
  infindividual ii;		//!< Infectious individual properties for this individual
  struct individual_* parent;   //!< Last who has infected this individual => *********** Careful with this!
  double nextchangetime;        //!< Time where the infection status of the individual will change next (for an active individual)
  uint8_t indinfstatus;         //!< Flags for the infection status of this individual
} individual;

inline static void ind_init_next_change_time(individual* ind) {
  ind->nextchangetime=-INFINITY;
}

inline static bool ind_update_next_change_time(individual* ind, const double time) {
  //This function updates the individual's infection status from
  //ro_ind_inf_status_latent to ro_ind_inf_status_infectious.

  if(time < ind->nextchangetime) return true;
  
  if(time >= ind->ii.end_comm_period) return false;

  const double timebuf=ind->ii.end_comm_period-ind->ii.comm_period;

  if(time < timebuf) {
    ind->nextchangetime=timebuf;
    return true;
  }

  ind->indinfstatus=ro_ind_inf_status_infectious;
  ind->nextchangetime=ind->ii.end_comm_period;
  return true;
}

inline static int ind_ct_comp(const void* ind1, const void* ind2)
{
  individual* left=(individual*)ind1;
  individual* right=(individual*)ind2;

  if(left->nextchangetime > right->nextchangetime) return 1;
  return -1;
}
#endif
