#include "main.h"

int main(const int nargs, const char* args[])
{
    printf("Hello World!\n");
    struct simulation sim;
    sim_init(&sim);

    gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);

    simulate(&sim, r);

    gsl_rng_free(r);
    return 0;
}
