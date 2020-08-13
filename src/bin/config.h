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

#define RF_P_EPS (1e-15)	//!< EPS for the change in p parameter
#define RF_GPERC_EPSF (1e-15)	//!< EPS for the x95 CDF discrepancy
#define RF_GKAPPA_EPSF (1e-15)  //!< EPS for the kappa CDF discrepancy

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

/**
 * @brief Solve for all simulation parameters.
 *
 * This function solves for the simulation parameters that have not been
 * provided as an input.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int config_solve_pars(sim_pars* pars);

/**
 * @brief Solve for the R0-related simulation parameters.
 *
 * This function solves for the R0-related simulation parameters
 * that have not been provided as an input.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int config_solve_R0_group(sim_pars* pars);

/**
 * @brief Solve for gamma distribution related simulation parameters.
 *
 * This function solves for gamma distribution related simulation parameters
 * that have not been provided as an input.
 *
 * @param ave: Gamma distribution average parameter.
 * @param kappa: Gamma distribution kappa parameter.
 * @param x95: Gamma distribution 95th percentile parameter.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int config_solve_gamma_group(double* ave, double* kappa, double* x95);

/**
 * @brief Computes p (logarithmic distribution) using Newton's method.
 *
 * This function estimates the value of the p parameter from a logarithmic
 * distribution, given the my parameter, using Newton's method
 *
 * @param x: Current value for p (input/output).
 * @param diff: Change in value for p (output).
 * @param params: Pointer to mu parameter (input)
 */
inline static void logroot(double* x, double* diff, void* params){const double x0=*x; const double l=log(1-x0); *x-=(*(double*)params*l*(1-x0)+x0)*(1-x0)/(x0/l+1); *diff=1 - x0 / *x;} 

/**
 * @brief Computes the discrepancy between an evaluation of gamma cumulative
 * distribution function and the 95th percentile.
 *
 * This function computes the discrepancy between an evaluation of gamma cumulative
 * distribution function and the 95th percentile.
 *
 * @param a: a parameter of the gamma distribution.
 * @param t: Evaluation point for the CDF.
 * @return the discrepancy between the evaluated CDF and the 95th percentile.
 */
inline static double gpercrootfunc(const double a, const double t){return gsl_sf_gamma_inc_P(a,t)-0.95;}

/**
 * @brief Computes x95 for a gamma cumulative distribution using Newton's
 * method.
 *
 * this function estimates the value of x95 (the x value corresponding to the
 * 95th percentile) for a gamma distribution, using Newton's method.
 *
 * @param x: Current value for x95 (input/output)
 * @param diff: Current value for the discrepancy between the CDF evaluated at
 * the current x95 parameter value and the 95th percentile. (output)
 * @param params: Pointer to an array of parameters. First parameter is the a parameter
 * for the gamma distribution and the second parameters is kappa. (input)
 */
inline static void gpercroot(double* x, double* diff, void* params){const double t=*(((double*)params)+1)* *x; *diff=gpercrootfunc(*(double*)params,t); *x-=*diff * tgamma(*(double*)params) / (pow(t,*(double*)params-1)*exp(-t) * *(((double*)params)+1));} 

/**
 * @brief Computes kappa for a gamma cumulative distribution using the secant
 * method.
 *
 * this function estimates the value of kappa for a gamma distribution,
 * using the secant method.
 *
 * @param x: Current value for kappa (input/output)
 * @param diff: Current value for the discrepancy between the CDF evaluated
 * using the current kappa parameter value and the 95th percentile. (output)
 * @param params: Pointer to an array of parameters. First parameter is average parameter
 * for the gamma distribution, the second parameters is x95, the third parameter
 * is the former kappa value and the fourth parameter is the former discrepancy
 * of the evaluated CDF value with the 95th percentile. (input/output)
 */
inline static void gkapparoot(double* x, double* diff, void* params){const double oldx=*x; *diff=gpercrootfunc(*(double*)params * *x,*(((double*)params)+1) * *x); *x-=*diff * (*x - *(((double*)params)+2)) / (*diff - *(((double*)params)+3)); *(((double*)params)+2)=oldx; *(((double*)params)+3)=*diff;} 

#endif
