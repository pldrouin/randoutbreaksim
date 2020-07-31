#ifndef _INFINDIVIDUAL_
#define _INFINDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

struct infindividual
{
    void* dataptr;
    double comm_period;       //Communicable period
    double event_time;
    uint32_t nevents;
    uint32_t curevent;
    uint32_t ninfections;
    uint32_t curinfection;
    bool infectious_at_tmax;
};

#endif
