/**
 * @file config.h
 * @brief Configuration functions the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _CONFIG_
#define _CONFIG_

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>

#include "args.h"
#include "simulation.h"

/**
 * @brief Configures the input parameters for the executable.
 *
 * @param pars: Simulation parameters
 * @param npaths: Number of simulation paths.
 * @param oout: File descriptor for the standard output.
 * @param eout: File descriptor for the standard error.
 * @param nargs: Number of arguments from the main function, shifted by 1.
 * @param args: Arguments from the main function, shifted by 1.
 */
int config(sim_pars* pars, uint32_t* npaths, int* oout, int* eout, const int nargs, const char* args[]);

/**
 * @brief Prints usage information for the executable.
 */
void printusage(const char* name);

#endif
