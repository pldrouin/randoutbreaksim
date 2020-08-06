#ifndef _CONFIG_
#define _CONFIG_

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>

#include "args.h"
#include "simulation.h"

int config(sim_pars* pars, uint32_t* npaths, int* oout, int* eout, const int nargs, const char* args[]);
void printusage(const char* name);

#endif
