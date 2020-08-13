#include "root_finder.h"

int root_finder_find(root_finder* rf, const double epsx, const double epsf, const uint32_t maxiter, const double xmin, const double xmax, double* x)
{
  double f;
  double lastx;
  uint32_t iter=0;
  printf("x=%22.15e\n",*x);

  do {
    lastx=*x;
    f=rf->func(*x, rf->params);
    printf("f=%f\n",f);
    *x -= f/rf->funcderiv(*x, rf->params);
    printf("x=%22.15e\n",*x);

    if(*x>xmax) *x=xmax;
    else if(*x<xmin) *x=xmin;

  } while((fabs(*x-lastx) > epsx || fabs(f) > epsf) && ++iter < maxiter);
  return -(iter==maxiter);
}
