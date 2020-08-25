/**
 * @file rnstream_gsl.h
 * @brief GSL adapter interface for rngstream.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _RNGSTREAM_GSL_
#define _RNGSTREAM_GSL_

#include <gsl/gsl_rng.h>

#include "rngstream.h"

inline static unsigned long int rngstream_get(void *vstate){return rng_rand_m1((rng_stream*)vstate);}

inline static double rngstream_get_double(void *vstate){return rng_rand_u01dm((rng_stream*)vstate);}
//inline static double rngstream_get_double(void *vstate){return rng_rand_u01d((rng_stream*)vstate);}
//inline static double rngstream_get_double(void *vstate){return rng_rand_u01((rng_stream*)vstate);}

inline static void rngstream_set(void *state, unsigned long int s) {rng_init((rng_stream*)state);}

static const gsl_rng_type rngstream_type = 
{"rngstream",                       /* name */
 __rngstream_m1-1,              /* RAND_MAX */
 0,                             /* RAND_MIN */
 sizeof (rng_stream),
 &rngstream_set,
 &rngstream_get,
 &rngstream_get_double
};

const gsl_rng_type *rngstream_gsl = &rngstream_type;

#endif
