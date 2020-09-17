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
#include <math.h>

#include "args.h"
#include "simulation.h"
#include "model_parameters.h"

/**
 * @brief Configures the input parameters for the executable.
 *
 * @param pars: Simulation parameters
 * @param ninfhist: Request or the generation of a histogram of number of
 * infections.
 * @param npaths: Number of simulation paths.
 * @param nthreads: Number of threads used to perform the simulation.
 * @param nsetsperthread: Number of path sets used for each thread.
 * @param lmax: Maximum number of layers for the simulation. lmax=1 means only primary infectious individuals.
 * @param nimax: Maximum number of infectious individuals for a given time
 * integet interval.
 * @param oout: File descriptor for the standard output.
 * @param eout: File descriptor for the standard error.
 * @param nargs: Number of arguments from the main function, shifted by 1.
 * @param args: Arguments from the main function, shifted by 1.
 */
int config(model_pars* pars, bool* ninfhist, uint32_t* npaths, uint32_t* nthreads, uint32_t* nsetsperthread, uint32_t* lmax, uint32_t* nimax, int* oout, int* eout, const int nargs, const char* args[]);

/**
 * @brief Prints usage information for the executable.
 */
void printusage(const char* name);
#endif
