#include <stdio.h>

#include "main.h"

int main(const int nargs, const char* args[])
{
    gsl_rng_env_setup();

    gsl_rng* r = gsl_rng_alloc(gsl_rng_taus2);
    struct sim_vars sv;
    //sim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);
    sim_init(&sv,r,.nstart=1,.tmax=25,.tbar=3,.lambda=0.1,.p=0.91,.kappa=0.47,.q=0); //"B5"

    //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

    for(int i=999999; i>=0; --i) simulate(&sv);

    sim_free(&sv);
    gsl_rng_free(r);
    return 0;
}
