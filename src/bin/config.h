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

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "args.h"
#include "simulation.h"
#include "root_finder.h"

#define RF_P_EPS (1e-15)
#define RF_GPERC_EPSF (1e-15)
#define RF_GKAPPA_EPSF (1e-15)

/**
 * @brief Configures the input parameters for the executable.
 *
 * @param pars: Simulation parameters
 * @param npaths: Number of simulation paths.
 * @param nimax: Maximum number of infectious individuals for a given time
 * integet interval.
 * @param oout: File descriptor for the standard output.
 * @param eout: File descriptor for the standard error.
 * @param nargs: Number of arguments from the main function, shifted by 1.
 * @param args: Arguments from the main function, shifted by 1.
 */
int config(sim_pars* pars, uint32_t* npaths, uint32_t* nimax, int* oout, int* eout, const int nargs, const char* args[]);

/**
 * @brief Prints usage information for the executable.
 */
void printusage(const char* name);

int config_solve_pars(sim_pars* pars);
int config_solve_R0_group(sim_pars* pars);
int config_solve_gamma_group(double* ave, double* kappa, double* p95);

//Finding p in logarithmic distribution using Newton's method
inline static void logroot(double* x, double* diff, void* params){const double x0=*x; const double l=log(1-x0); *x-=(*(double*)params*l*(1-x0)+x0)*(1-x0)/(x0/l+1); *diff=*x-x0;} 

inline static double gpercrootfunc(const double a, const double t){return gsl_sf_gamma_inc_P(a,t)-0.95;}
//Finding p95 for the gamma cumulative distribution function using Newton's method
inline static void gpercroot(double* x, double* diff, void* params){const double t=*(((double*)params)+1)* *x; *diff=gpercrootfunc(*(double*)params,t); *x-=*diff * tgamma(*(double*)params) / (pow(t,*(double*)params-1)*exp(-t) * *(((double*)params)+1));} 

//Finding kappa for the gamma cumulative distribution function using secant method
inline static void gkapparoot(double* x, double* diff, void* params){const double oldx=*x; *diff=gpercrootfunc(*(double*)params * *x,*(((double*)params)+1) * *x); *x-=*diff * (*x - *(((double*)params)+2)) / (*diff - *(((double*)params)+3)); *(((double*)params)+2)=oldx; *(((double*)params)+3)=*diff;} 

#endif
