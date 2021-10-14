/**
 * @file inflayer.h
 * @brief Data structure containing a layer of infectious individual for
 * branchsim.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
*/

#ifndef _INFLAYER_
#define _INFLAYER_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#include "infindividual.h"

/**
 * Infectious individual layer
 */
typedef struct inflayer_
{
  infindividual ii;
  double event_time;	      //!< Start time for the current iteration event
} inflayer;

#endif
