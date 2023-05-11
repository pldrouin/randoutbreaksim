/**
 * @file model_parameters.h
 * @brief Model parameter functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#ifndef _MODEL_PARAMETERS_
#define _MODEL_PARAMETERS_

#include <stdbool.h>
#include <inttypes.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

#include "root_finder.h"

#define RF_P_EPSF (1e-15)	//!< EPS for the mu and g_ave discrepancy
#define RF_GPERC_EPSF (1e-15)	//!< EPS for the x95 CDF discrepancy
#define RF_GKAPPA_EPSF (1e-15)  //!< EPS for the kappa CDF discrepancy
#define RF_GAUSSMU_EPSF (1e-15)  //!< EPS for the Gaussian mu mean discrepancy

/**
 * Primary individual communicable period model type. Flags used to specify the model.
 **/
enum ro_pricommper_model_flags {ro_pricommper_main=1, ro_pricommper_alt=2, ro_pricommper_alt_use_tpr=4, ro_pricommper_first_cat=8, ro_pricommper_second_cat=16};

/**
 * Time model type.
 **/
enum ro_time_model{ro_time_pri_created=1, ro_time_pri_flat_comm=2, ro_time_pri_infectious=3, ro_time_pri_end_comm=4, ro_time_pri_test_results=5, ro_time_first_pos_test_results=6};

/**
 * Group model type. Flags used to specify the model. ro_log_group_attendees_plus_1,
 * ro_log_group_invitees_plus_1, ro_log_attendees and ro_log_invitees are mutually
 * exclusive options that indicate how the number of individuals at events are distributed.  
 **/
enum ro_group_model_flags {ro_group_invitees=1, ro_group_log_plus_1=2, ro_group_log=4, ro_group_gauss=8, ro_group_dist_mask=ro_group_log_plus_1|ro_group_log|ro_group_gauss};

/**
 * Path model.
 */
enum ro_path_model {ro_all_paths, ro_observable_paths_only, ro_non_observable_paths_only};

/**
 * Model parameters.
 */
typedef struct 
{
  double tbar;		//!< Mean main communicable period
  double p;		//!< Parameter for the logarithmic distribution used to draw then number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (0 <= p < 1).
  double mu;		//!< Parameter for the mean of an unbounded logarithmic distribution (mu=-1/log(1-p)*p/(1-p), mu >= 1) or of an unbounded Gaussian used to draw number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.
  double sigma;         //!< Parameter for the standard deviation of an unbounded Gaussian used to draw the number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.
  double rsigma;        //!< Parameter for the relative standard deviation of an unbounded Gaussian used to draw the number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.
  double g_ave;		//!< Parameter for the average group size for one event. These individuals can correspond to invitees or attendees depending on the choice of group type. Events are defined to include at least two invitees, and their definition also depends on the usage of group_interactions or group_transmission
  double g_ave_transm;	//!< Parameter for the average group size for one transmission event. These individuals can correspond to invitees or attendees depending on the choice of group type. Events are defined to include at least two invitees, and with one infectious individual. Equal to g_ave_transm if group_transmission is used
  double lambda;	//!< Rate of events for a given individual. Events are defined to include at least two invitees.
  double lambda_uncut;  //!< Rate of events for a given individual, including events of one invitee.
  double lambdap;	//!< Total rate of events for a finite population. Events are defined to include at least two invitees
  double pinfpri;       //!< Probability that a primary individual be infections.
  double pinf;		//!< Probability that a given susceptible individual gets infected when exposed to one infectious individual during one event.
#ifdef DUAL_PINF
  double ppip;          //!< Probability that a given susceptible individual be in the second infection probability category.
  double rpinfp;        //!< Relative probability that a given susceptible individual of the second category gets infected when exposed to one infectious individual during one event (value relative to pinf, 0 < rpinfp * pinf <= 1, default value of 1).
  double rpshedp;       //!> Relative strength of infectiousness from an infectious individual of the second category vs the fist category (value relative to pinf, 0 < rpshedp * pinf <=1, default value of 1).
  double qp;            //!> Probability of alternate communicable period for an infectious individual in the second category.
#endif
  double R0;		//!< Basic reproduction number
  double kappa;		//!< branchim's kappa parameter for the gamma distribution used to generate the main communicable period
  double lbar;		//!< Mean latent period
  double kappal;   	//!< kappa parameter for the gamma distribution used to generate the latent period
  double q;		//!< Probability of alternate communicable period
  double mbar;		//!< Mean period for the alternate communicable period
  double kappaq;	//!< branchim's kappa parameter for the gamma distribution used to generate the alternate communicable period
#ifdef CT_OUTPUT
  double ctwindow;      //!< Period prior to individual isolation during which contacts are considered
  double pt;		//!< Probability of successful contact tracing. Probability must be larger than pit and pim, as it is considered to be applicable to all contacts.
  double pitnet;        //!< Probability of main communicable period interruption, net of the probability of tracing the contact.
  double pimnet;        //!< Probability of alternate communicable period interruption, net of the probability of tracing the contact.
#endif
  double pit;		//!< Overall probability of main communicable period interruption
  double itbar;		//!< Mean period for the interrupted main communicable period
  double kappait;	//!< kappa parameter for the gamma distribution used to generate the interrupted main communicable period
  double pim;		//!< Overall probability of alternate communicable period interruption
  double imbar;		//!< Mean period for the interrupted alternate communicable period
  double kappaim;	//!< kappa parameter for the gamma distribution used to generate the interrupted alternate communicable period
  double t95;	        //!< Parameter for the 95th percentile of the main communicable period (depends on tbar and kappa)
  double ta;		//!< a parameter for the main communicable period (kappa * tbar)
  double tb;		//!< b parameter for the main communicable period (1 / kappa)
  double m95;	        //!< Parameter for the 95th percentile of the alternate communicable period (depends on mbar and kappaq)
  double ma;		//!< a parameter for the alternate communicable period (kappaq * mbar)
  double mb;		//!< b parameter for the alternate communicable period (1 / kappaq)
  double l95;	        //!< Parameter for the 95th percentile of the latent period (depends on lbar and kappal)
  double la;		//!< a parameter for the latent period (kappal * lbar)
  double lb;		//!< b parameter for the latent period (1 / kappal)
  double it95;	        //!< Parameter for the 95th percentile of the interrupted main communicable period (depends on itbar and kappait)
  double ita;		//!< a parameter for the interrupted main communicable period (kappait * itbar)
  double itb;		//!< b parameter for the interrupted main communicable period (1 / kappait)
  double im95;	        //!< Parameter for the 95th percentile of the interrupted alternate communicable period (depends on imbar and kappaim)
  double ima;		//!< a parameter for the interrupted alternate communicable period (kappaim * imbar)
  double imb;		//!< b parameter for the interrupted alternate communicable period (1 / kappaim)
  double ttpr;		//!< True positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the main period, for which a positive test would allow for the interruption of a child's communicable period
  double mtpr;		//!< True positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the alternate period, for which a positive test would allow for the interruption of a child's communicable period
  double tdeltat;	//!< Time delay between the end of the applicable communicable period and test results.
  int32_t tmax;		//!< Maximum simulation time units used to instantiate new infectious individuals
  uint32_t nstart;	//!< Initial number of individuals
  uint32_t popsize;     //!< Population size (finite population simulation)
  uint8_t pricommpertype;	//!< Primary individual communicable period type bit field (composed using ro_pricommper_model_flags) 
  uint8_t grouptype;   	//!< Group type bit field (composed using ro_group_model_flags)
  uint8_t timetype;     //!< Recorded time origin (value set using ro_time_model)
  uint8_t pathtype;        //!< Indicate which type of path, observable or not, should be included in the simulation results (value set using ro_path_model)
  bool groupinteractions;
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
 * @brief Solve for the value of the p, mu and g_ave parameters of the
 * logarithmic plus 1 distribution, given one of the parameters.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_log_plus_1_group(model_pars* pars);

/**
 * @brief Solve for the value of the p, mu and g_ave parameters of the
 * logarithmic distribution, given one of the parameters.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_log_group(model_pars* pars);

/**
 * @brief Solve for the value of the mu, sigma and g_ave parameters of the
 * Gaussian distribution, given two of the parameters.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_gauss_group(model_pars* pars);

/**
 * @brief Solve for the value of lambda uncut, given lambda, for the log
 * plus 1 distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_log_plus_1_lambda_uncut_from_lambda(model_pars* pars){pars->lambda_uncut=pars->lambda; return 0;}

/**
 * @brief Solve for the value of lambda uncut, given lambda, for the log
 * distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_log_lambda_uncut_from_lambda(model_pars* pars){const double l1mp=log(1-pars->p); pars->lambda_uncut=(pars->p==0?INFINITY:l1mp/(l1mp+pars->p)*pars->lambda); return 0;}

/**
 * @brief Solve for the value of lambda uncut, given lambda, for the Gaussian
 * distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_gauss_lambda_uncut_from_lambda(model_pars* pars){pars->lambda_uncut=pars->lambda/gsl_cdf_ugaussian_Q((1.5-pars->mu)/pars->sigma); return 0;}

/**
 * @brief Solve for the value of lambda, given lambda uncut, for the log
 * plus 1 distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_log_plus_1_lambda_from_lambda_uncut(model_pars* pars){pars->lambda=pars->lambda_uncut; return 0;}

/**
 * @brief Solve for the value of lambda, given lambda uncut, for the log
 * distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_log_lambda_from_lambda_uncut(model_pars* pars){

  if(pars->p==0) {
    fprintf(stderr,"%s: Error: lambda cannot be computed from lambda_uncut for group_log_attendees distribution if p=0\n",__func__);
    return -1;
  }
  const double l1mp=log(1-pars->p);
  pars->lambda=(l1mp+pars->p)/l1mp*pars->lambda_uncut;
  return 0;
}

/**
 * @brief Solve for the value of lambda, given lambda uncut, for the Gaussian
 * distribution.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
inline static int model_solve_gauss_lambda_from_lambda_uncut(model_pars* pars){pars->lambda=pars->lambda_uncut*gsl_cdf_ugaussian_Q((1.5-pars->mu)/pars->sigma); return 0;}

/**
 * @brief Solve for the value of the p parameter of the logarithmic
 * distribution, given mu.
 *
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_log_p_from_mu(model_pars* pars);

/**
 * @brief Solve for the value of the p parameter of the truncated logarithmic
 * distribution below 2, given a mean.
 *
 * @param mean: Mean of the truncated logarithmic distribution below 2.
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_trunc_log_p_from_mean(const double mean, model_pars* pars);

/**
 * @brief Solve for the value of the p parameter of the logarithmic plus 1
 * distribution, given the mean for the transmission events.
 *
 * @param mean: Mean of invitees for the transmission events.
 * @param pars: Simulation parameters.
 * @return 0 if the parameters could be determined, and a non-zero value
 * otherwise.
 */
int model_solve_log_p_plus_1_from_transm_mean(model_pars* pars);

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
 * @brief Computes p (logarithmic distribution) from mu using Newton's method.
 *
 * This function estimates the value of the p parameter from a logarithmic
 * distribution, given the mu parameter, using Newton's method
 *
 * @param x: Current value for p (input/output).
 * @param diff: mu value discrepancy (output).
 * @param params: Pointer to mu parameter (input)
 */
inline static void logroot(double* x, double* diff, void* params){const double omx=1-*x; const double l=log(omx); *diff=*(double*)params+ *x/(omx*l); *x-=*diff*l*omx*omx/(*x/l+1); *diff/=*(double*)params;} 

/**
 * @brief Computes p (logarithmic distribution) from the mean of the
 * transmission events when a logarithmic distribution plus one is used for
 * all group interaction events, using Newton's method.
 *
 * @param x: Current value for p (input/output).
 * @param diff: transmission event mean value discrepancy (output).
 * @param params: Pointer to transmission mean parameter (input)
 */
inline static void gilogp1root(double* x, double* diff, void* params){
  const double omx=1-*x;
  const double l=log(omx);
  const double xpl=*x+l;
  const double omxl=omx*l;
  const double omxlmx=omxl-*x;
  const double xpldomxlmx=xpl/omxlmx;
  const double opxpldomxlmx = 1 + xpldomxlmx;
  const double mean=-*x/omxl * opxpldomxlmx;
  *diff=mean - *(double*)params + 1;
  *x-=*diff * omxl / (-opxpldomxlmx + mean*(1+l) + *x / omxlmx * (*x / omx + xpldomxlmx * (l+2)));
  *diff/=*(double*)params - 1;
} 

/**
 * @brief Computes p (truncated logarithmic distribution) from the mean using Newton's method.
 *
 * This function estimates the value of the p parameter from a truncated logarithmic
 * distribution below 2, given the mean parameter, using Newton's method
 *
 * @param x: Current value for p (input/output).
 * @param diff: mean value discrepancy (output).
 * @param params: Pointer to mean parameter (input)
 */
inline static void loggt1root(double* x, double* diff, void* params){const double omx=1-*x; const double l=log(omx); const double lpx=l+*x; *diff=*(double*)params+ *x * *x/(omx*lpx); *x-=*diff*lpx*lpx*omx*omx/(*x * (2 * lpx - *x * l)); *diff/=*(double*)params;} 

/**
 * @brief Computes g_ave (truncated Gaussian distribution) from the mean and the
 * standard deviation.
 */
inline static double gauss_trunc_g_ave(const double mu, const double sigma)
{
  int32_t i; //Discrete random value relative to mui
  int32_t mui=floor(mu+0.5);
  const double dmu=mu-mui;
  const double psmu=dmu-0.5;
  const double nsmu=dmu+0.5;
  int32_t nbins=mui-2;
  int32_t lasti=0;
  double lastrangeint;
  double newrangeint;
  double lastint=gsl_cdf_ugaussian_P((-nbins-1-psmu)/sigma);;
  double newint;
  double mean=0;
  double fint=0;
  double dint;
  double meanbuf;
  double fintbuf;

  if(mu != mui) {

    for(i=-nbins; i<=0; ++i) {
      //printf("i: %i, lastint: %22.15e, fint: %22.15e, mean: %22.15e\n", i, lastint, fint, mean);

      newint=gsl_cdf_ugaussian_P((i-psmu)/sigma);
      dint=newint-lastint;
      //printf("dint: %22.15e\n",dint);
      fint+=dint;
      mean+=dint*i;
      lastint=newint;
    }
    lastrangeint=lastint;

  } else {
    mean=0;
    lastrangeint=gsl_cdf_ugaussian_P((nbins-psmu)/sigma);
    fint=lastrangeint-lastint;
    lasti=nbins;
    nbins*=2;
    printf("lasti: %i, nbins: %i\n",lasti,nbins);
  }

  for(;;) {
    fintbuf=0;
    meanbuf=0;
    newrangeint=lastint=gsl_cdf_ugaussian_P((lasti+nbins-psmu)/sigma);

    for(i=lasti+nbins; i>lasti+1; --i) {
      printf("i: %i, lastint: %22.15e, fintbuf: %22.15e, meanbuf: %22.15e\n", i, lastint, fintbuf, meanbuf);
      newint=gsl_cdf_ugaussian_P((i-nsmu)/sigma);
      dint=lastint-newint;
      printf("dint: %22.15e\n",dint);
      fintbuf+=dint;
      meanbuf+=dint*i;
      lastint=newint;
    }
    printf("i: %i, lastint: %22.15e, fintbuf: %22.15e, meanbuf: %22.15e\n", lasti+1, lastint, fintbuf, meanbuf);
    dint=lastint-lastrangeint;
    printf("dint: %22.15e\n",dint);
    fintbuf+=dint;
    meanbuf+=dint*(lasti+1);

    fint+=fintbuf;
    mean+=meanbuf;
    if(newrangeint==1) break;
    lastrangeint=newrangeint;
    lasti+=nbins;
    nbins*=2;
    printf("nbins to %i\n",nbins);
  }
  printf("mean=%22.15e, fint=%22.15e, mui=%i\n",mean,fint,mui);

  return mean/fint+mui;
}

/**
 * @brief Computes mu for a discretized truncated Gaussian distribution, given g_ave and sigma,
 * using the secant method.
 *
 * this function estimates the value of mu for a discretized truncaed Gaussian distribution,
 * using the secant method.
 *
 * @param x: Current value for mu (input/output)
 * @param diff: Current value for the discrepancy between the mean evaluated
 * using the current mu parameter value and g_ave. (output)
 * @param params: Pointer to an array of parameters. First parameter is the sigma parameter
 * for the Gaussian distribution, the second parameters is g_ave, the third parameter
 * is the former mu value and the fourth parameter is the former discrepancy
 * of the evaluated mean value with g_ave. (input/output)
 */
inline static void gaussmuroot(double* x, double* diff, void* params){const double oldx=*x; *diff=gauss_trunc_g_ave(*x, (*(double*)params))-*(((double*)params)+1); /*printf("%22.15e %22.15e %22.15e %22.15e\n",*diff,*x,*(((double*)params)+2),*(((double*)params)+3));*/ *x-=*diff * (*x - *(((double*)params)+2)) / (*diff - *(((double*)params)+3)); *(((double*)params)+2)=oldx; *(((double*)params)+3)=*diff;} 

/**
 * @brief Computes mu for a discretized truncated Gaussian distribution, given g_ave and rsigma,
 * using the secant method.
 *
 * this function estimates the value of mu for a discretized truncaed Gaussian distribution,
 * using the secant method.
 *
 * @param x: Current value for mu (input/output)
 * @param diff: Current value for the discrepancy between the mean evaluated
 * using the current mu parameter value and g_ave. (output)
 * @param params: Pointer to an array of parameters. First parameter is the rsigma parameter
 * for the Gaussian distribution, the second parameters is g_ave, the third parameter
 * is the former mu value and the fourth parameter is the former discrepancy
 * of the evaluated mean value with g_ave. (input/output)
 */
inline static void gaussrmuroot(double* x, double* diff, void* params){const double oldx=*x; *diff=gauss_trunc_g_ave(*x, *x * (*(double*)params))-*(((double*)params)+1); /*printf("%22.15e %22.15e %22.15e %22.15e\n",*diff,*x,*(((double*)params)+2),*(((double*)params)+3));*/ *x-=*diff * (*x - *(((double*)params)+2)) / (*diff - *(((double*)params)+3)); *(((double*)params)+2)=oldx; *(((double*)params)+3)=*diff;} 

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
inline static void gkapparoot(double* x, double* diff, void* params){const double oldx=*x; *diff=gpercrootfunc(*(double*)params * *x,*(((double*)params)+1) * *x); /*printf("%22.15e %22.15e %22.15e %22.15e\n",*diff,*x,*(((double*)params)+2),*(((double*)params)+3));*/ *x-=*diff * (*x - *(((double*)params)+2)) / (*diff - *(((double*)params)+3)); *(((double*)params)+2)=oldx; *(((double*)params)+3)=*diff;} 

#endif
