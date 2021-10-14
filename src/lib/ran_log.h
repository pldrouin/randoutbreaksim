/**
 * @file ran_log.h
 * @brief Common simulation data structures and functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _RAN_LOG_
#define _RAN_LOG_

#include <stdint.h>
#include <math.h>

#include "rngstream.h"

typedef struct {
  rng_stream* s;
  double p;
  double r;
} ran_log;

inline static void ran_log_init(ran_log* rl, rng_stream* s, const double p){rl->s=s; rl->p=p; rl->r=log(1-p);}

/**
 * @brief Finite logarithmic deviate.
 *
 * Returns a logarithmic deviate. Returned value is finite, and can be cast to
 * the desired integer data type. Modified algorithm from "Non-Uniform Random Variate Generation
", by Luc Devroye.
 *
 * @param rl handle
 * @return finite logarithmic deviate
 */
#define ran_log_finite_gen(rl) \
{\
  double v=rng_rand_pu01(rl->s); /*We want a finite returned value;*/ \
 \
  if (v>=rl->p) return 1; \
 \
  else { \
    double q=1-exp(rl->r*rng_rand_pu01(rl->s)); \
 \
    if (v<=q*q) return 1+log(v)/log(q); \
 \
    else if (v <= q) return 2; \
 \
    else return 1 ; \
  } \
}
inline static uint32_t ran_log_finite(ran_log* rl){ran_log_finite_gen(rl);}
inline static uint64_t ran_log_finite_l(ran_log* rl){ran_log_finite_gen(rl);}

/**
 * @brief Capped logarithmic deviate.
 *
 * Returns a logarithmic deviate. Returned value is capped, and can be cast to
 * the desired integer data type. Modified algorithm from "Non-Uniform Random Variate Generation
", by Luc Devroye.
 *
 * @param rl handle
 * @return finite logarithmic deviate
 */
#define ran_log_capped_gen(rl, max) \
{\
  for(;;) {\
    double v=rng_rand_pu01(rl->s); /*We want a finite returned value;*/ \
    \
    if (v>=rl->p) return 1; \
    \
    else { \
      double q=1-exp(rl->r*rng_rand_pu01(rl->s)); \
      \
      if (v<=q*q) {\
	v=1+log(v)/log(q);\
	\
	if(v>max) continue; \
	return v; \
      \
      } else if (v <= q) { \
	\
	if(2>max) continue; \
	return 2; \
      } \
      \
      else return 1 ; \
    } \
  }\
}
inline static uint32_t ran_log_capped(ran_log* rl, const uint32_t max){ran_log_capped_gen(rl, max);}
inline static uint64_t ran_log_capped_l(ran_log* rl, const uint64_t max){ran_log_capped_gen(rl, max);}

/**
 * @brief Finite logarithmic deviate with a lower bound of 2.
 *
 * Returns a logarithmic deviate greater than 1. Returned value is finite, and can be cast to
 * the desired integer data type. Modified algorithm from "Non-Uniform Random Variate Generation"
 * by Luc Devroye. Loop seems necessary. Note: There is no support for p=0 by
 * this function due to optimisation.
 *
 * @param rl handle
 * @return finite logarithmic deviate greater than 1
 */
#define ran_log_finite_gt1_gen(rl) \
{\
  for(;;) { \
    double v=rl->p*rng_rand_pu01(rl->s); /*We want a finite returned value;*/ \
    \
    while(v==rl->p) v=rl->p*rng_rand_pu01(rl->s); \
    \
    double q=1-exp(rl->r*rng_rand_pu01(rl->s)); \
    \
    if (v<=q*q) return 1+log(v)/log(q); \
    \
    else if (v <= q) return 2; \
    \
    else continue ; \
  } \
}
inline static uint32_t ran_log_finite_gt1(ran_log* rl){ran_log_finite_gt1_gen(rl);}
inline static uint64_t ran_log_finite_gt1_l(ran_log* rl){ran_log_finite_gt1_gen(rl);}

/**
 * @brief Capped logarithmic deviate with a lower bound of 2. Cap value must be
 * at least 2
 *
 * Returns a logarithmic deviate greater than 1. Returned value is capped, and can be cast to
 * the desired integer data type. Modified algorithm from "Non-Uniform Random Variate Generation"
 * by Luc Devroye. Loop seems necessary. Note: There is no support for p=0 by
 * this function due to optimisation.
 *
 * @param rl handle
 * @return finite logarithmic deviate greater than 1
 */
#define ran_log_capped_gt1_gen(rl, max) \
{\
  for(;;) { \
    double v=rl->p*rng_rand_pu01(rl->s); /*We want a finite returned value;*/ \
    \
    while(v==rl->p) v=rl->p*rng_rand_pu01(rl->s); \
    \
    double q=1-exp(rl->r*rng_rand_pu01(rl->s)); \
    \
    if (v<=q*q) {\
      v=1+log(v)/log(q); \
      \
      if(v>max) continue; \
      return v; \
      \
    } else if (v <= q) return 2; \
    \
    else continue ; \
  } \
}
inline static uint32_t ran_log_capped_gt1(ran_log* rl, const uint32_t max){ran_log_capped_gt1_gen(rl, max);}
inline static uint64_t ran_log_capped_gt1_l(ran_log* rl, const uint64_t max){ran_log_capped_gt1_gen(rl, max);}

#endif
