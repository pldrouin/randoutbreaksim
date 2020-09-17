/**
 * @file main.h
 * @brief Main function for the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include <stdio.h>
#include <unistd.h>

#include <math.h>

#include <pthread.h>

#include <gsl/gsl_rng.h>

#include "rngstream_gsl.h"

#include "config.h"
#include "model_parameters.h"
#include "infindividual.h"
#include "branchsim.h"
#include "standard_summary_stats.h"

/**
 * @brief Main function.
 */
int main(const int nargs, const char* args[]);

typedef struct {
  double npathsperset;
  uint32_t nsets;
  uint32_t npers;
  uint32_t tnpersa;
  uint32_t lmax;
  uint32_t nimax;
  model_pars const* pars;
  uint32_t volatile* set;
  double commper_mean;
#ifdef NUMEVENTSSTATS
  double nevents_mean;
#endif
  double r_mean;
  double pe;
  double te_mean;
  double te_std;
  double* inf_timeline_mean_ext;
  double* inf_timeline_std_ext;
  double* inf_timeline_mean_noext;
  double* inf_timeline_std_noext;
  double* totinf_timeline_mean_ext;
  double* totinf_timeline_std_ext;
  double* totinf_timeline_mean_noext;
  double* totinf_timeline_std_noext;
  double* totmainatt_timeline_mean_ext;
  double* totmainatt_timeline_std_ext;
  double* totmainatt_timeline_mean_noext;
  double* totmainatt_timeline_std_noext;
  uint64_t* ngeninfs;
  uint32_t ninfbins;
  uint32_t nimaxedoutmintimeindex;
  gsl_rng* r;
  bool rec_ninfs;
} thread_data;

void* simthread(void* arg);
