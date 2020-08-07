/**
 * @file main.c
 * @brief Main function for the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include "main.h"

int main(const int nargs, const char* args[])
{
  sim_pars pars={.tbar=0, .p=0, .lambda=0, .kappa=0, .q=0, .mbar=0, .kappaq=0, .tmax=INFINITY, .nstart=1};
  uint32_t npaths=10000;
  int oout=STDOUT_FILENO;
  int eout=STDERR_FILENO;

  if(config(&pars, &npaths, &oout, &eout, nargs-1, args+1)) return 1;
  gsl_rng_env_setup();

  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus2);
  sim_vars sv;
  //sim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);

  if(sim_init(&sv,&pars,r)) {
    fprintf(stderr,"%s: Error: While attempting to initialise the simulation\n",args[0]);
    return 1;
  }

  std_summary_stats stats;

  sim_set_proc_data(&sv, &stats);
  sim_set_increase_layers_proc_func(&sv, std_stats_increase_layers);
  sim_set_new_event_proc_func(&sv, std_stats_new_event);
  sim_set_new_inf_proc_func(&sv, std_stats_new_inf);
  sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
  sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf);

  //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

  double commper_mean=0;
  double nevents_mean=0;
  double r_mean=0;
  //uint32_t nr=0;
  double pe=0;
  double te_mean=0, te_std=0;
  double inf_timeline_mean_ext[(int)sv.pars.tmax+1];
  double inf_timeline_std_ext[(int)sv.pars.tmax+1];
  double inf_timeline_mean_noext[(int)sv.pars.tmax+1];
  double inf_timeline_std_noext[(int)sv.pars.tmax+1];
  double totinf_timeline_mean_ext[(int)sv.pars.tmax+1];
  double totinf_timeline_std_ext[(int)sv.pars.tmax+1];
  double totinf_timeline_mean_noext[(int)sv.pars.tmax+1];
  double totinf_timeline_std_noext[(int)sv.pars.tmax+1];

  std_stats_init(&sv, &stats);

  memset(inf_timeline_mean_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_mean_noext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_noext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_mean_ext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_std_ext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_mean_noext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_std_noext, 0, stats.npers*sizeof(double));
  int j;

  for(int i=npaths-1; i>=0; --i) {
    std_stats_path_init(&stats);
    simulate(&sv);
    r_mean+=stats.rsum;
    commper_mean+=stats.commpersum;
    nevents_mean+=stats.neventssum;
    //nr+=stats.n_ended_infections;

    if(stats.extinction) {
      pe+=stats.extinction;
      te_mean+=stats.extinction_time;
      te_std+=stats.extinction_time*stats.extinction_time;

      inf_timeline_mean_ext[0]+=stats.inf_timeline[0];
      inf_timeline_std_ext[0]+=stats.inf_timeline[0]*stats.inf_timeline[0];
      totinf_timeline_mean_ext[0]+=stats.totinf_timeline[0];
      totinf_timeline_std_ext[0]+=stats.totinf_timeline[0]*stats.totinf_timeline[0];

      for(j=1; j<stats.npers; ++j) {
	inf_timeline_mean_ext[j]+=stats.inf_timeline[j];
	inf_timeline_std_ext[j]+=stats.inf_timeline[j]*stats.inf_timeline[j];
	stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	totinf_timeline_mean_ext[j]+=stats.totinf_timeline[j];
	totinf_timeline_std_ext[j]+=stats.totinf_timeline[j]*stats.totinf_timeline[j];
      }

    } else {
      inf_timeline_mean_noext[0]+=stats.inf_timeline[0];
      inf_timeline_std_noext[0]+=stats.inf_timeline[0]*stats.inf_timeline[0];
      totinf_timeline_mean_noext[0]+=stats.totinf_timeline[0];
      totinf_timeline_std_noext[0]+=stats.totinf_timeline[0]*stats.totinf_timeline[0];

      for(j=1; j<stats.npers; ++j) {
	inf_timeline_mean_noext[j]+=stats.inf_timeline[j];
	inf_timeline_std_noext[j]+=stats.inf_timeline[j]*stats.inf_timeline[j];
	stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	totinf_timeline_mean_noext[j]+=stats.totinf_timeline[j];
	totinf_timeline_std_noext[j]+=stats.totinf_timeline[j]*stats.totinf_timeline[j];
      }
    }

  }
  std_stats_free(&sv, &stats);
  const double ninf=totinf_timeline_mean_ext[stats.npers-2]+totinf_timeline_mean_noext[stats.npers-2];
  const double ninf_per_event_mean=r_mean/nevents_mean;
  const double nnoe=npaths-pe;
  r_mean/=ninf;
  commper_mean/=ninf;
  nevents_mean/=ninf;
  te_mean/=pe;
  te_std=sqrt(pe/(pe-1)*(te_std/pe-te_mean*te_mean));

  for(j=stats.npers-1; j>=0; --j) {
    inf_timeline_mean_ext[j]/=pe;
    inf_timeline_std_ext[j]=sqrt(pe/(pe-1)*(inf_timeline_std_ext[j]/pe-inf_timeline_mean_ext[j]*inf_timeline_mean_ext[j]));
    totinf_timeline_mean_ext[j]/=pe;
    totinf_timeline_std_ext[j]=sqrt(pe/(pe-1)*(totinf_timeline_std_ext[j]/pe-totinf_timeline_mean_ext[j]*totinf_timeline_mean_ext[j]));

    inf_timeline_mean_noext[j]/=nnoe;
    inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1)*(inf_timeline_std_noext[j]/nnoe-inf_timeline_mean_noext[j]*inf_timeline_mean_noext[j]));
    totinf_timeline_mean_noext[j]/=nnoe;
    totinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1)*(totinf_timeline_std_noext[j]/nnoe-totinf_timeline_mean_noext[j]*totinf_timeline_mean_noext[j]));
  }
  pe/=npaths;

  printf("Mean R is %f\n",r_mean);
  printf("Uninterrupted communication period is %f\n",commper_mean);
  printf("Number of events per infectious individual is %f\n",nevents_mean);
  printf("Number of infections per event is %f\n",ninf_per_event_mean);
  printf("Probability of extinction is %f\n",pe);
  printf("Extinction time, if it occurs is %f +/- %f\n",te_mean,te_std);

  printf("Current infection timeline, for paths with extinction vs no extinction is:\n");
  for(j=0; j<stats.npers; ++j) printf("%3i: %9.4f +/- %9.4f\t%9.4f +/- %9.4f\n",j,inf_timeline_mean_ext[j],inf_timeline_std_ext[j],inf_timeline_mean_noext[j],inf_timeline_std_noext[j]);

  printf("Total infections timeline, for paths with extinction vs no extinction is:\n");
  for(j=0; j<stats.npers; ++j) printf("%3i: %9.4f +/- %9.4f\t%9.4f +/- %9.4f\n",j,totinf_timeline_mean_ext[j],totinf_timeline_std_ext[j],totinf_timeline_mean_noext[j],totinf_timeline_std_noext[j]);

  sim_free(&sv);
  gsl_rng_free(r);
  close(oout);
  close(eout);
  return 0;
}
