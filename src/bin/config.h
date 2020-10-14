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
 * @brief Struct used to store configuration parmeters.
 */
typedef struct
{
model_pars pars; 		//!< Simulation parameters
bool ninfhist;			//!< Request or the generation of a histogram of number of infections.
bool timerelfirstpostestresults;  //!< Indicate that the output time values are relative to the time of the first positive test results.
uint32_t npaths;		//!< Number of simulation paths.
uint32_t lmax;			//!< Maximum number of layers for the simulation. lmax=1 means only primary infectious individuals.
int32_t nbinsperunit;		//!< Number of timeline bins per unit of time.
uint32_t nimax;			//!< Maximum number of infectious individuals for a given time integet interval.
uint32_t npostestmax;		//!< Maximum number of positive test results during an interval of duration npostestmaxnpers for each individual that starts when the test results are received.
uint32_t npostestmaxnpers;      //!< Interval duration for the maximum number of positive test results
uint32_t nthreads;		//!< Number of threads used to perform the simulation.
uint32_t nsetsperthread;	//!< Number of path sets used for each thread.
uint32_t stream;		//!< RNG stream index.
uint32_t tloutbufsize;		//!< Per-thread memory buffer size (in MB) used to accumulate data for timeline output before writing them to disk.
int tlout;			//!< File descriptor used to record timeline data for each simulated path.
#ifdef CT_OUTPUT
uint32_t ctoutbufsize;		//!< Per-thread memory buffer size (in MB) used to accumulate data for contact tracing output before writing them to disk.
int ctout;			//!< File descriptor used to CT data for each simulated path.
#endif
int oout;			//!< File descriptor for the standard output.
int eout;			//!< File descriptor for the standard error.
} config_pars;

/**
 * @brief Configures the input parameters for the executable.
 *
 * @param cp: Configuration parameters
 * @param nargs: Number of arguments from the main function, shifted by 1.
 * @param args: Arguments from the main function, shifted by 1.
 */
int config(config_pars *cp, const int nargs, const char* args[]);

/**
 * @brief Prints usage information for the executable.
 */
void printusage(const char* name);
#endif
