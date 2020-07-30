#ifndef _INFINDIVIDUAL_
#define _INFINDIVIDUAL_

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

struct infindividual
{
    void* dataptr;
    double comm_period;
    double event_time;
    int nevents;
    int curevent;
    int ninfections;
    int curinfection;
};

#endif
