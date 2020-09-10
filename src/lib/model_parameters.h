/**
 * @file model_parameters.h
 * @brief Model parameter functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _MODEL_PARAMETERS_
#define _MODEL_PARAMETERS_

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "root_finder.h"

#define RF_P_EPSF (1e-15)	//!< EPS for the mu discrepancy
#define RF_GPERC_EPSF (1e-15)	//!< EPS for the x95 CDF discrepancy
#define RF_GKAPPA_EPSF (1e-15)  //!< EPS for the kappa CDF discrepancy

/**
 * Model type. Flags used to specify the model. ro_log_group_plus_1, ro_log_attendees and
 * ro_log_invitees are mutually exclusive options that indicate how the number
 * of individuals at events are distributed.  
 **/
enum ro_model_flags {ro_group_log_plus_1=0, ro_group_log_attendees=1, ro_group_log_invitees=2};

#define RO_GROUP_DIST_MASK (0x3)

/**
 * Model parameters.
 */
typedef struct 
{
  double tbar;		//!< Mean main communicable period
  double p;		//!< Parameter for the logarithmic distribution used to draw then number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (0 <= p < 1).
  double lambda;	//!< Rate of events for a given individual. Events are defined to have more than one invitees
  double lambdap;	//!< Total rate of events for a finite population. Events are defined to have more than one invitees
  double kappa;		//!< branchim's kappa parameter for the gamma distribution used to generate the main communicable period
  double lbar;		//!< Mean latent period
  double kappal;   	//!< kappa parameter for the gamma distribution used to generate the latent period
  double q;		//!< Probability of alternate communicable period
  double mbar;		//!< Mean period for the alternate communicable period
  double kappaq;	//!< branchim's kappa parameter for the gamma distribution used to generate the alternate communicable period
  double pit;		//!< Probability of main communicable period interruption
  double itbar;		//!< Mean period for the interrupted main communicable period
  double kappait;	//!< kappa parameter for the gamma distribution used to generate the interrupted main communicable period
  double pim;		//!< Probability of alternate communicable period interruption
  double imbar;		//!< Mean period for the interrupted alternate communicable period
  double kappaim;	//!< kappa parameter for the gamma distribution used to generate the interrupted alternate communicable period
  double R0;		//!< Basic reproduction number
  double mu;		//!< Parameter for the mean of an unbounded logarithmic distribution used to draw number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (mu=-1/log(1-p)*p/(1-p), mu >= 1)
  double t95;	        //!< Parameter for the 95th percentile of the main communicable period (depends on tbar and kappa)
  double m95;	        //!< Parameter for the 95th percentile of the alternate communicable period (depends on mbar and kappaq)
  double l95;	        //!< Parameter for the 95th percentile of the latent period (depends on lbar and kappal)
  double it95;	        //!< Parameter for the 95th percentile of the interrupted main communicable period (depends on itbar and kappait)
  double im95;	        //!< Parameter for the 95th percentile of the interrupted alternate communicable period (depends on imbar and kappaim)
  double tmax;		//!< Maximum simulation period used to instantiate new infectious individuals
  double pinf;		//!< Probability that a given susceptible individual gets infected when exposed to one infectious individual during one event.
  uint32_t nstart;	//!< Initial number of infectious individuals
  uint32_t popsize;     //!< Population size (finite population simulation)
  uint8_t grouptype;   	//!< Simulation type bit field (composed using ro_model_type)
} model_pars;

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
int model_solve_pars(model_pars* pars);

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
int model_solve_R0_group(model_pars* pars);

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
int model_solve_gamma_group(double* ave, double* kappa, double* x95);

/**
 * @brief Verifies the validity of the model parameters.
 *
 * This internal function verifies if the set of provided model parameter values are
 * valid.
 *
 * @param pars: Pointer to the model parameters.
 * @return 0 if the parameters are valid. If the parameters are
 * invalid, an error is printed on stderr.
 */
int model_pars_check(model_pars const* pars);

/**
 * @brief Computes p (logarithmic distribution) using Newton's method.
 *
 * This function estimates the value of the p parameter from a logarithmic
 * distribution, given the my parameter, using Newton's method
 *
 * @param x: Current value for p (input/output).
 * @param diff: mu value discrepancy (output).
 * @param params: Pointer to mu parameter (input)
 */
inline static void logroot(double* x, double* diff, void* params){const double omx=1-*x; const double l=log(omx); *diff=*(double*)params+ *x/(omx*l); *x-=*diff*l*omx*omx/(*x/l+1); *diff/=*(double*)params;} 

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
