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

inline static double logroot(double x, void* params){return *(double*)params+x/((1-x)*log(1-x));} 
inline static double logrootderiv(double x, void* params){const double l=log(1-x); return (x/l+1)/(l*(1-x)*(1-x));} 
inline static double gpercroot(double x, void* params){return gsl_sf_gamma_inc_P(*(double*)params,*(((double*)params)+1)*x)-0.95;} 
inline static double gpercrootderiv(double x, void* params){const double t=*(((double*)params)+1)*x; return pow(t,*(double*)params-1)*exp(-t)/tgamma(*(double*)params) * *(((double*)params)+1);} 
inline static double gkapparoot(double x, void* params){return gsl_sf_gamma_inc_P(*(double*)params*x,*(((double*)params)+1)*x)-0.95;} 
inline static double gkapparootderiv(double x, void* params){const double a=*(double*)params * x; const double t=*(((double*)params)+1)*x; return -gsl_sf_psi(a) * *(double*)params * gsl_sf_gamma_inc_P(a,t) +  pow(t,a-1) * exp(-t)/tgamma(a) * *(((double*)params)+1);} 

#endif
