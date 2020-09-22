/**
 * @file root_finder.c
 * @brief Finds the root of a function depending on a single parameter.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "root_finder.h"

int root_finder_find(root_finder* rf, const double eps, const uint32_t maxiter, const double xmin, const double xmax, double* x, double* calcdiff)
{
  double diff;
  uint32_t iter=0;
  double oldx=NAN;
  bool samex;
  //printf("x=%22.15e (xmin=%22.15e, xmax=%22.15e)\n",*x,xmin,xmax);

  do {
    rf->func(x, &diff, rf->params);
    //printf("x=%22.15e, diff=%22.15e",*x,diff);

    if(*x==oldx) samex=true;

    else {
      samex=false;
      oldx=*x;
    }

    if(*x>xmax) *x=xmax;
    else if(*x<xmin) *x=xmin;

    //printf(" (corrected x=%22.15e)\n",*x);

  } while(!(fabs(diff) < eps) && ++iter < maxiter);

  if(calcdiff) *calcdiff=diff;
  return -samex-2*(iter==maxiter);
}
