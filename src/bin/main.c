/**
 * @file main.c
 * @brief Main function for the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "main.h"

int main(const int nargs, const char* args[])
{
  model_pars pars;
  uint32_t npaths=10000;
  uint32_t nimax=UINT32_MAX;
  int oout=STDOUT_FILENO;
  int eout=STDERR_FILENO;

  sim_pars_init(&pars);

  if(config(&pars, &npaths, &nimax, &oout, &eout, nargs-1, args+1)) return 1;

  if(model_solve_pars(&pars)) return 1;
  gsl_rng_env_setup();

  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus2);
  sim_vars sv;
  //sim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);

  if(sim_init(&sv,&pars,r)) {
    fprintf(stderr,"%s: Error: While attempting to initialise the simulation\n",args[0]);
    return 1;
  }

  //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

  double commper_mean=0;
#ifdef NUMEVENTSSTATS
  double nevents_mean=0;
#endif
  double r_mean=0;
  //uint32_t nr=0;
  double pe=0;
  double te_mean=0, te_std=0;
  double inf_timeline_mean[(int)sv.pars.tmax+1];
  double inf_timeline_std[(int)sv.pars.tmax+1];
  double inf_timeline_mean_ext[(int)sv.pars.tmax+1];
  double inf_timeline_std_ext[(int)sv.pars.tmax+1];
  double inf_timeline_mean_noext[(int)sv.pars.tmax+1];
  double inf_timeline_std_noext[(int)sv.pars.tmax+1];
  double totinf_timeline_mean[(int)sv.pars.tmax+1];
  double totinf_timeline_std[(int)sv.pars.tmax+1];
  double totinf_timeline_mean_ext[(int)sv.pars.tmax+1];
  double totinf_timeline_std_ext[(int)sv.pars.tmax+1];
  double totinf_timeline_mean_noext[(int)sv.pars.tmax+1];
  double totinf_timeline_std_noext[(int)sv.pars.tmax+1];

  std_summary_stats stats;

  sim_set_proc_data(&sv, &stats);
  sim_set_increase_layers_proc_func(&sv, std_stats_increase_layers);

  if(nimax == UINT32_MAX) sim_set_new_event_proc_func(&sv, std_stats_new_event);
  else sim_set_new_event_proc_func(&sv, std_stats_new_event_nimax);
  sim_set_new_inf_proc_func(&sv, std_stats_new_inf);
  sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
  sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf);

  std_stats_init(&sv);
  stats.nimax=nimax;

  memset(inf_timeline_mean_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_mean_noext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_noext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_mean_ext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_std_ext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_mean_noext, 0, stats.npers*sizeof(double));
  memset(totinf_timeline_std_noext, 0, stats.npers*sizeof(double));
  int j;
  uint32_t nimaxedoutmintimeindex=UINT32_MAX;

  for(int i=npaths-1; i>=0; --i) {
    std_stats_path_init(&stats);
    simulate(&sv);
    r_mean+=stats.rsum;
    commper_mean+=stats.commpersum;
#ifdef NUMEVENTSSTATS
    nevents_mean+=stats.neventssum;
#endif
    //nr+=stats.n_ended_infections;

    if(stats.extinction) {
      pe+=stats.extinction;
      te_mean+=stats.extinction_time;
      te_std+=stats.extinction_time*stats.extinction_time;

      inf_timeline_mean_ext[0]+=stats.inf_timeline[0];
      inf_timeline_std_ext[0]+=(double)stats.inf_timeline[0]*stats.inf_timeline[0];
      totinf_timeline_mean_ext[0]+=stats.totinf_timeline[0];
      totinf_timeline_std_ext[0]+=(double)stats.totinf_timeline[0]*stats.totinf_timeline[0];

      for(j=1; j<stats.npers; ++j) {
	inf_timeline_mean_ext[j]+=stats.inf_timeline[j];
	inf_timeline_std_ext[j]+=(double)stats.inf_timeline[j]*stats.inf_timeline[j];
	stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	totinf_timeline_mean_ext[j]+=stats.totinf_timeline[j];
	totinf_timeline_std_ext[j]+=(double)stats.totinf_timeline[j]*stats.totinf_timeline[j];
      }

    } else {

      if(stats.nimaxedoutmintimeindex < nimaxedoutmintimeindex) nimaxedoutmintimeindex=stats.nimaxedoutmintimeindex;
      inf_timeline_mean_noext[0]+=stats.inf_timeline[0];
      inf_timeline_std_noext[0]+=(double)stats.inf_timeline[0]*stats.inf_timeline[0];
      totinf_timeline_mean_noext[0]+=stats.totinf_timeline[0];
      totinf_timeline_std_noext[0]+=(double)stats.totinf_timeline[0]*stats.totinf_timeline[0];

      for(j=1; j<stats.npers; ++j) {
	inf_timeline_mean_noext[j]+=stats.inf_timeline[j];
	inf_timeline_std_noext[j]+=(double)stats.inf_timeline[j]*stats.inf_timeline[j];
	stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	totinf_timeline_mean_noext[j]+=stats.totinf_timeline[j];
	totinf_timeline_std_noext[j]+=(double)stats.totinf_timeline[j]*stats.totinf_timeline[j];
      }
    }
  }
  std_stats_free(&sv, &stats);
  const double ninf=totinf_timeline_mean_ext[stats.npers-2]+totinf_timeline_mean_noext[stats.npers-2];
#ifdef NUMEVENTSSTATS
  const double ninf_per_event_mean=r_mean/nevents_mean;
#endif
  const double nnoe=npaths-pe;
  r_mean/=ninf;
  commper_mean/=ninf;
#ifdef NUMEVENTSSTATS
  nevents_mean/=ninf;
#endif
  te_mean/=pe;
  te_std=sqrt(pe/(pe-1.)*(te_std/pe-te_mean*te_mean));

  for(j=stats.npers-1; j>=0; --j) {
    inf_timeline_mean[j]=inf_timeline_mean_ext[j]+inf_timeline_mean_noext[j];
    inf_timeline_std[j]=inf_timeline_std_ext[j]+inf_timeline_std_noext[j];
    totinf_timeline_mean[j]=totinf_timeline_mean_ext[j]+totinf_timeline_mean_noext[j];
    totinf_timeline_std[j]=totinf_timeline_std_ext[j]+totinf_timeline_std_noext[j];

    inf_timeline_mean[j]/=npaths;
    inf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(inf_timeline_std[j]/npaths-inf_timeline_mean[j]*inf_timeline_mean[j]));
    totinf_timeline_mean[j]/=npaths;
    totinf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(totinf_timeline_std[j]/npaths-totinf_timeline_mean[j]*totinf_timeline_mean[j]));

    inf_timeline_mean_ext[j]/=pe;
    inf_timeline_std_ext[j]=sqrt(pe/(pe-1.)*(inf_timeline_std_ext[j]/pe-inf_timeline_mean_ext[j]*inf_timeline_mean_ext[j]));
    totinf_timeline_mean_ext[j]/=pe;
    totinf_timeline_std_ext[j]=sqrt(pe/(pe-1.)*(totinf_timeline_std_ext[j]/pe-totinf_timeline_mean_ext[j]*totinf_timeline_mean_ext[j]));

    inf_timeline_mean_noext[j]/=nnoe;
    inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(inf_timeline_std_noext[j]/nnoe-inf_timeline_mean_noext[j]*inf_timeline_mean_noext[j]));
    totinf_timeline_mean_noext[j]/=nnoe;
    totinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(totinf_timeline_std_noext[j]/nnoe-totinf_timeline_mean_noext[j]*totinf_timeline_mean_noext[j]));
  }
  pe/=npaths;

  printf("Mean R is %f\n",r_mean);
  printf("Communicable period is %f\n",commper_mean);
#ifdef NUMEVENTSSTATS
  printf("Number of events per infectious individual is %f\n",nevents_mean);
  printf("Number of infections per event is %f\n",ninf_per_event_mean);
#endif
  printf("Probability of extinction and its statistical uncertainty: %f +/- %f%s\n",pe,sqrt(pe*(1.-pe)/(npaths-1.)),(nimaxedoutmintimeindex<UINT32_MAX?" (nimax reached, could be biased)":""));
  printf("Extinction time, if it occurs is %f +/- %f%s\n",te_mean,te_std,(nimaxedoutmintimeindex<UINT32_MAX?" (nimax reached, could be biased)":""));

  printf("Current infection timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<stats.npers; ++j) printf("%3i: %11.4f +/- %11.4f\t%11.4f +/- %11.4f\t%11.4f +/- %11.4f%s\n",j,inf_timeline_mean_ext[j],inf_timeline_std_ext[j],inf_timeline_mean_noext[j],inf_timeline_std_noext[j],inf_timeline_mean[j],inf_timeline_std[j],(j<nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  printf("Total infections timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<stats.npers; ++j) printf("%3i: %11.4f +/- %11.4f\t%11.4f +/- %11.4f\t%11.4f +/- %11.4f%s\n",j,totinf_timeline_mean_ext[j],totinf_timeline_std_ext[j],totinf_timeline_mean_noext[j],totinf_timeline_std_noext[j],totinf_timeline_mean[j],totinf_timeline_std[j],(j<nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  sim_free(&sv);
  gsl_rng_free(r);
  fflush(stdout);
  fflush(stderr);
  close(oout);
  close(eout);
  return 0;
}
