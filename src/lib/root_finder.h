#ifndef _ROOT_FINDER_
#define _ROOT_FINDER_

#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct
{
  void (*func)(double* x, double* diff, void* params);
  void* params;
} root_finder;

inline static root_finder* root_finder_init(void (*func)(double* x, double* diff, void* params), void* params){
  root_finder* ret=(root_finder*)malloc(sizeof(root_finder));

  if(ret) {
    ret->func=func;
    ret->params=params;
  }
  return ret;
}

inline static void root_finder_free(root_finder* rf){free(rf);}

int root_finder_find(root_finder* rf, const double eps, const uint32_t maxiter, const double xmin, const double xmax, double* x);

#endif
