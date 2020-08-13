/**
 * @file root_finder.c
 * @brief Finds the root of a function depending on a single parameter.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "root_finder.h"

int root_finder_find(root_finder* rf, const double eps, const uint32_t maxiter, const double xmin, const double xmax, double* x)
{
  double diff;
  uint32_t iter=0;
  //printf("x=%22.15e\n",*x);

  do {
    rf->func(x, &diff, rf->params);
    //printf("x=%22.15e, diff=%22.15e\n",*x,diff);

    if(*x>xmax) *x=xmax;
    else if(*x<xmin) *x=xmin;

  } while(!(fabs(diff) < eps) && ++iter < maxiter);
  return -(iter==maxiter);
}
