#include <stdio.h>

#include "main.h"

int main(const int nargs, const char* args[])
{
    struct sim_pars sim;
    sim_init(&sim);

    gsl_rng_env_setup();

    gsl_rng* r = gsl_rng_alloc(gsl_rng_taus2);

    //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

    simulate(&sim, r);

    gsl_rng_free(r);
    return 0;
}
