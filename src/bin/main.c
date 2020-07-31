#include <stdio.h>

#include "main.h"

int main(const int nargs, const char* args[])
{
  gsl_rng_env_setup();

  gsl_rng* r = gsl_rng_alloc(gsl_rng_taus2);
  struct sim_vars sv;
  //sim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);
  const int npaths=1000000;
  sim_init(&sv,r,.nstart=1,.tmax=25,.tbar=3,.lambda=0.1,.p=0.913068,.kappa=0.466367,.q=0); //"B5"

  struct std_summary_stats stats;

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
  double ne_mean=0, ne_std=0;
  double te_mean=0, te_std=0;
  double ng_mean=0, ng_std=0;
  double inf_timeline_mean_ext[(int)sv.pars.tmax+1];
  double inf_timeline_std_ext[(int)sv.pars.tmax+1];
  double inf_timeline_mean_noext[(int)sv.pars.tmax+1];
  double inf_timeline_std_noext[(int)sv.pars.tmax+1];

  std_stats_init(&sv, &stats);

  memset(inf_timeline_mean_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_ext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_mean_noext, 0, stats.npers*sizeof(double));
  memset(inf_timeline_std_noext, 0, stats.npers*sizeof(double));
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
      ne_mean+=stats.total_n_infections;
      ne_std+=stats.total_n_infections*stats.total_n_infections;
      te_mean+=stats.extinction_time;
      te_std+=stats.extinction_time*stats.extinction_time;

      for(j=stats.npers-1; j>=0; --j) {
	inf_timeline_mean_ext[j]+=stats.inf_timeline[j];
	inf_timeline_std_ext[j]+=stats.inf_timeline[j]*stats.inf_timeline[j];
      }

    } else {
      ng_mean+=stats.total_n_infections;
      ng_std+=stats.total_n_infections*=stats.total_n_infections;

      for(j=stats.npers-1; j>=0; --j) {
	inf_timeline_mean_noext[j]+=stats.inf_timeline[j];
	inf_timeline_std_noext[j]+=stats.inf_timeline[j]*stats.inf_timeline[j];
      }
    }

  }
  std_stats_free(&sv, &stats);
  const double ninf=ne_mean+ng_mean;
  const double ninf_per_event_mean=r_mean/nevents_mean;
  const double nnoe=npaths-pe;
  r_mean/=ninf;
  commper_mean/=ninf;
  nevents_mean/=ninf;
  ne_mean/=pe;
  ne_std=sqrt(pe/(pe-1)*(ne_std/pe-ne_mean*ne_mean));
  te_mean/=pe;
  te_std=sqrt(pe/(pe-1)*(te_std/pe-te_mean*te_mean));
  ng_mean/=nnoe;
  ng_std=sqrt(nnoe/(nnoe-1)*(ng_std/nnoe-ng_mean*ng_mean));

  for(j=stats.npers-1; j>=0; --j) {
    inf_timeline_mean_ext[j]/=pe;
    inf_timeline_std_ext[j]=sqrt(pe/(pe-1)*(inf_timeline_std_ext[j]/pe-inf_timeline_mean_ext[j]*inf_timeline_mean_ext[j]));
  }

  for(j=stats.npers-1; j>=0; --j) {
    inf_timeline_mean_noext[j]/=nnoe;
    inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1)*(inf_timeline_std_noext[j]/nnoe-inf_timeline_mean_noext[j]*inf_timeline_mean_noext[j]));
  }
  pe/=npaths;


  printf("Mean R is %f\n",r_mean);
  printf("Uninterrupted communication period is %f\n",commper_mean);
  printf("Number of events per infectious individual is %f\n",nevents_mean);
  printf("Number of infections per event is %f\n",ninf_per_event_mean);
  printf("Probability of extinction is %f\n",pe);
  printf("Total number of infected individuals at extinction is %f +/- %f\n",ne_mean,ne_std);
  printf("Extinction time, if it occurs is %f +/- %f\n",te_mean,te_std);
  printf("Total number of infected individuals if no extinction is %f +/- %f\n",ng_mean,ng_std);

  printf("Average infection timeline, for paths with extinction vs no extinction is:\n");
  for(j=0; j<stats.npers; ++j) printf("%3i: %8.4f +/- %8.4f\t%8.4f +/- %8.4f\n",j,inf_timeline_mean_ext[j],inf_timeline_std_ext[j],inf_timeline_mean_noext[j],inf_timeline_std_noext[j]);

  sim_free(&sv);
  gsl_rng_free(r);
  return 0;
}
