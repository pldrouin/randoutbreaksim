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
  int ninfrec=0;
  uint32_t npaths=10000;
  uint32_t nthreads=1;
  uint32_t nsetsperthread=(nthreads>1?100:1);
  uint32_t nimax=UINT32_MAX;
  int oout=STDOUT_FILENO;
  int eout=STDERR_FILENO;

  sim_pars_init(&pars);

  if(config(&pars, &ninfrec, &npaths, &nthreads, &nsetsperthread, &nimax, &oout, &eout, nargs-1, args+1)) return 1;

  if(model_solve_pars(&pars)) return 1;

  gsl_rng_env_setup();

  if(model_pars_check(&pars)) {
    fprintf(stderr,"%s: Error: While verifying the validity of the simulation parameters.\n",args[0]);
    return 1;
  }

  thread_data* tdata=(thread_data*)malloc(nthreads*sizeof(thread_data));
  int j;
  const uint32_t nsets=nthreads*nsetsperthread;
  const double npathsperset=((double)npaths)/nsets;
  const uint32_t npers=pars.tmax+1;
  volatile uint32_t set=0;
  int t;
  uint32_t** ngeninfs=NULL;
  uint32_t* ninfs=NULL;

  if(ninfrec>0) {
    ngeninfs=(uint32_t**)malloc(npaths*sizeof(uint32_t*));
    ninfs=(uint32_t*)malloc(npaths*sizeof(uint32_t));
  }

  if(nthreads>1) {
    pthread_t* threads=(pthread_t*)malloc(nthreads*sizeof(pthread_t));

    for(t=nthreads-1; t>=0; --t){
      tdata[t].npathsperset=npathsperset;
      tdata[t].nsets=nsets;
      tdata[t].npers=npers;
      tdata[t].nimax=nimax;
      tdata[t].pars=&pars;
      tdata[t].set=&set;
      tdata[t].ngeninfs=ngeninfs;
      tdata[t].ninfs=ninfs;
      //tdata[t].r = gsl_rng_alloc(gsl_rng_taus2);
      tdata[t].r = gsl_rng_alloc(rngstream_gsl);
      tdata[t].rec_ninfs=(ninfrec>0);
      pthread_create(threads+t,NULL,simthread,tdata+t);
    }

    pthread_join(threads[0],NULL);
    gsl_rng_free(tdata[0].r);

    for(t=1; t<nthreads; ++t) {
      pthread_join(threads[t],NULL);
      gsl_rng_free(tdata[t].r);
      tdata[0].r_mean+=tdata[t].r_mean;
      tdata[0].commper_mean+=tdata[t].commper_mean;
#ifdef NUMEVENTSSTATS
      tdata[0].nevents_mean+=tdata[t].nevents_means;
#endif
      tdata[0].pe+=tdata[t].pe;
      tdata[0].te_mean+=tdata[t].te_mean;
      tdata[0].te_std+=tdata[t].te_std;

      for(j=0; j<npers; ++j) {
        tdata[0].inf_timeline_mean_ext[j]+=tdata[t].inf_timeline_mean_ext[j];
        tdata[0].inf_timeline_std_ext[j]+=tdata[t].inf_timeline_std_ext[j];
        tdata[0].totinf_timeline_mean_ext[j]+=tdata[t].totinf_timeline_mean_ext[j];
        tdata[0].totinf_timeline_std_ext[j]+=tdata[t].totinf_timeline_std_ext[j];

        tdata[0].inf_timeline_mean_noext[j]+=tdata[t].inf_timeline_mean_noext[j];
        tdata[0].inf_timeline_std_noext[j]+=tdata[t].inf_timeline_std_noext[j];
        tdata[0].totinf_timeline_mean_noext[j]+=tdata[t].totinf_timeline_mean_noext[j];
        tdata[0].totinf_timeline_std_noext[j]+=tdata[t].totinf_timeline_std_noext[j];
      }

      if(tdata[t].nimaxedoutmintimeindex < tdata[0].nimaxedoutmintimeindex) tdata[0].nimaxedoutmintimeindex = tdata[t].nimaxedoutmintimeindex;
    }
    free(threads);

  } else {
    tdata[0].npathsperset=npathsperset;
    tdata[0].nsets=nsets;
    tdata[0].npers=npers;
    tdata[0].nimax=nimax;
    tdata[0].pars=&pars;
    tdata[0].set=&set;
    tdata[0].ngeninfs=ngeninfs;
    tdata[0].ninfs=ninfs;
    //tdata[0].r = gsl_rng_alloc(gsl_rng_taus2);
    tdata[0].r = gsl_rng_alloc(rngstream_gsl);
    tdata[0].rec_ninfs=(ninfrec>0);
    simthread(tdata);
    gsl_rng_free(tdata[0].r);
  }

  const double ninf=tdata[0].totinf_timeline_mean_ext[npers-2]+tdata[0].totinf_timeline_mean_noext[npers-2];
#ifdef NUMEVENTSSTATS
  const double ninf_per_event_mean=tdata[0].r_mean/tdata[0].nevents_mean;
#endif
  const double nnoe=npaths-tdata[0].pe;
  tdata[0].r_mean/=ninf;
  tdata[0].commper_mean/=ninf;
#ifdef NUMEVENTSSTATS
  tdata[0].nevents_mean/=ninf;
#endif
  tdata[0].te_mean/=tdata[0].pe;
  tdata[0].te_std=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[0].te_std/tdata[0].pe-tdata[0].te_mean*tdata[0].te_mean));

  double inf_timeline_mean[npers];
  double inf_timeline_std[npers];
  double totinf_timeline_mean[npers];
  double totinf_timeline_std[npers];

  for(j=npers-1; j>=0; --j) {
    inf_timeline_mean[j]=tdata[0].inf_timeline_mean_ext[j]+tdata[0].inf_timeline_mean_noext[j];
    inf_timeline_std[j]=tdata[0].inf_timeline_std_ext[j]+tdata[0].inf_timeline_std_noext[j];
    totinf_timeline_mean[j]=tdata[0].totinf_timeline_mean_ext[j]+tdata[0].totinf_timeline_mean_noext[j];
    totinf_timeline_std[j]=tdata[0].totinf_timeline_std_ext[j]+tdata[0].totinf_timeline_std_noext[j];

    inf_timeline_mean[j]/=npaths;
    inf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(inf_timeline_std[j]/npaths-inf_timeline_mean[j]*inf_timeline_mean[j]));
    totinf_timeline_mean[j]/=npaths;
    totinf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(totinf_timeline_std[j]/npaths-totinf_timeline_mean[j]*totinf_timeline_mean[j]));

    tdata[0].inf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[0].inf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[0].inf_timeline_std_ext[j]/tdata[0].pe-tdata[0].inf_timeline_mean_ext[j]*tdata[0].inf_timeline_mean_ext[j]));
    tdata[0].totinf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[0].totinf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[0].totinf_timeline_std_ext[j]/tdata[0].pe-tdata[0].totinf_timeline_mean_ext[j]*tdata[0].totinf_timeline_mean_ext[j]));

    tdata[0].inf_timeline_mean_noext[j]/=nnoe;
    tdata[0].inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[0].inf_timeline_std_noext[j]/nnoe-tdata[0].inf_timeline_mean_noext[j]*tdata[0].inf_timeline_mean_noext[j]));
    tdata[0].totinf_timeline_mean_noext[j]/=nnoe;
    tdata[0].totinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[0].totinf_timeline_std_noext[j]/nnoe-tdata[0].totinf_timeline_mean_noext[j]*tdata[0].totinf_timeline_mean_noext[j]));
  }
  tdata[0].pe/=npaths;

  printf("Mean R is %f\n",tdata[0].r_mean);
  printf("Communicable period is %f\n",tdata[0].commper_mean);
#ifdef NUMEVENTSSTATS
  printf("Number of events per infectious individual is %f\n",tdata[0].nevents_mean);
  printf("Number of infections per event is %f\n",ninf_per_event_mean);
#endif
  printf("Probability of extinction and its statistical uncertainty: %f +/- %f%s\n",tdata[0].pe,sqrt(tdata[0].pe*(1.-tdata[0].pe)/(npaths-1.)),(tdata[0].nimaxedoutmintimeindex<UINT32_MAX?" (nimax reached, could be biased)":""));
  printf("Extinction time, if it occurs is %f +/- %f%s\n",tdata[0].te_mean,tdata[0].te_std,(tdata[0].nimaxedoutmintimeindex<UINT32_MAX?" (nimax reached, could be biased)":""));

  printf("Current infection timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<npers; ++j) printf("%3i: %11.4f +/- %11.4f\t%11.4f +/- %11.4f\t%11.4f +/- %11.4f%s\n",j,tdata[0].inf_timeline_mean_ext[j],tdata[0].inf_timeline_std_ext[j],tdata[0].inf_timeline_mean_noext[j],tdata[0].inf_timeline_std_noext[j],inf_timeline_mean[j],inf_timeline_std[j],(j<tdata[0].nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  printf("Total infections timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<npers; ++j) printf("%3i: %11.4f +/- %11.4f\t%11.4f +/- %11.4f\t%11.4f +/- %11.4f%s\n",j,tdata[0].totinf_timeline_mean_ext[j],tdata[0].totinf_timeline_std_ext[j],tdata[0].totinf_timeline_mean_noext[j],tdata[0].totinf_timeline_std_noext[j],totinf_timeline_mean[j],totinf_timeline_std[j],(j<tdata[0].nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  for(t=nthreads-1; t>=0; --t) {
    free(tdata[t].inf_timeline_mean_ext);
    free(tdata[t].inf_timeline_std_ext);
    free(tdata[t].inf_timeline_mean_noext);
    free(tdata[t].inf_timeline_std_noext);
    free(tdata[t].totinf_timeline_mean_ext);
    free(tdata[t].totinf_timeline_std_ext);
    free(tdata[t].totinf_timeline_mean_noext);
    free(tdata[t].totinf_timeline_std_noext);
  }
  free(tdata);

  if(ninfrec>0) {
    uint32_t p;

    for(p=0; p<npaths; ++p) {

      if(write(ninfrec, ngeninfs[p], ninfs[p]*sizeof(uint32_t)) !=  ninfs[p]*sizeof(uint32_t)) {
	perror(args[0]);
	return 1;
      }
      free(ngeninfs[p]);
    }
    free(ngeninfs);
    free(ninfs);
  }

  fflush(stdout);
  fflush(stderr);
  close(oout);
  close(eout);

  if(ninfrec>0) close(ninfrec);
  return 0;
}

void* simthread(void* arg)
{
  thread_data* data=(thread_data*)arg;
  data->commper_mean=0;
#ifdef NUMEVENTSSTATS
  data->nevents_mean=0;
#endif
  data->r_mean=0;
  data->pe=0;
  data->te_mean=0;
  data->te_std=0;
  data->nimaxedoutmintimeindex=UINT32_MAX;

  data->inf_timeline_mean_ext=(double*)malloc(data->npers*sizeof(double));
  data->inf_timeline_std_ext=(double*)malloc(data->npers*sizeof(double));
  data->inf_timeline_mean_noext=(double*)malloc(data->npers*sizeof(double));
  data->inf_timeline_std_noext=(double*)malloc(data->npers*sizeof(double));
  data->totinf_timeline_mean_ext=(double*)malloc(data->npers*sizeof(double));
  data->totinf_timeline_std_ext=(double*)malloc(data->npers*sizeof(double));
  data->totinf_timeline_mean_noext=(double*)malloc(data->npers*sizeof(double));
  data->totinf_timeline_std_noext=(double*)malloc(data->npers*sizeof(double));

  memset(data->inf_timeline_mean_ext, 0, data->npers*sizeof(double));
  memset(data->inf_timeline_std_ext, 0, data->npers*sizeof(double));
  memset(data->inf_timeline_mean_noext, 0, data->npers*sizeof(double));
  memset(data->inf_timeline_std_noext, 0, data->npers*sizeof(double));
  memset(data->totinf_timeline_mean_ext, 0, data->npers*sizeof(double));
  memset(data->totinf_timeline_std_ext, 0, data->npers*sizeof(double));
  memset(data->totinf_timeline_mean_noext, 0, data->npers*sizeof(double));
  memset(data->totinf_timeline_std_noext, 0, data->npers*sizeof(double));

  sim_vars sv;
  //branchsim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);

  sim_init(&sv,data->pars,data->r);

  //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

  //uint32_t nr=0;

  std_summary_stats stats;

  sim_set_proc_data(&sv, &stats);
  sim_set_ii_alloc_proc_func(&sv, std_stats_ii_alloc);

  if(data->nimax == UINT32_MAX) sim_set_new_event_proc_func(&sv, std_stats_new_event);
  else sim_set_new_event_proc_func(&sv, std_stats_new_event_nimax);
  sim_set_new_inf_proc_func(&sv, std_stats_new_inf);

  if(data->rec_ninfs) {
    sim_set_end_inf_proc_func(&sv, std_stats_end_inf_rec_ninfs);
    sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf_rec_ninfs);

  } else {
    sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
    sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf);
  }

  branchsim_init(&sv);
  std_stats_init(&sv);
  stats.nimax=data->nimax;
  int j;
  uint32_t curset;
  uint32_t initpath;
  uint32_t npaths;

  for(;;) {
    curset=__sync_fetch_and_add(data->set,1);

    if(curset>=data->nsets) break;

    //printf("%22.15e\t%22.15e\n",curset*data->npathsperset,(curset+1)*data->npathsperset);
    initpath=round(curset*data->npathsperset);
    npaths=round((curset+1)*data->npathsperset)-initpath;
    //printf("npaths %u\n",npaths);

    for(int i=npaths-1; i>=0; --i) {
      std_stats_path_init(&stats);
      branchsim(&sv);
      data->r_mean+=stats.rsum;
      data->commper_mean+=stats.commpersum;
#ifdef NUMEVENTSSTATS
      data->nevents_mean+=stats.neventssum;
#endif
      //nr+=stats.n_ended_infections;

      if(stats.extinction) {
	data->pe+=stats.extinction;
	data->te_mean+=stats.extinction_time;
	data->te_std+=stats.extinction_time*stats.extinction_time;

	data->inf_timeline_mean_ext[0]+=stats.inf_timeline[0];
	data->inf_timeline_std_ext[0]+=(double)stats.inf_timeline[0]*stats.inf_timeline[0];
	data->totinf_timeline_mean_ext[0]+=stats.totinf_timeline[0];
	data->totinf_timeline_std_ext[0]+=(double)stats.totinf_timeline[0]*stats.totinf_timeline[0];

	for(j=1; j<stats.npers; ++j) {
	  data->inf_timeline_mean_ext[j]+=stats.inf_timeline[j];
	  data->inf_timeline_std_ext[j]+=(double)stats.inf_timeline[j]*stats.inf_timeline[j];
	  stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totinf_timeline_mean_ext[j]+=stats.totinf_timeline[j];
	  data->totinf_timeline_std_ext[j]+=(double)stats.totinf_timeline[j]*stats.totinf_timeline[j];
	}

      } else {

	if(stats.nimaxedoutmintimeindex < data->nimaxedoutmintimeindex) data->nimaxedoutmintimeindex=stats.nimaxedoutmintimeindex;
	data->inf_timeline_mean_noext[0]+=stats.inf_timeline[0];
	data->inf_timeline_std_noext[0]+=(double)stats.inf_timeline[0]*stats.inf_timeline[0];
	data->totinf_timeline_mean_noext[0]+=stats.totinf_timeline[0];
	data->totinf_timeline_std_noext[0]+=(double)stats.totinf_timeline[0]*stats.totinf_timeline[0];

	for(j=1; j<stats.npers; ++j) {
	  data->inf_timeline_mean_noext[j]+=stats.inf_timeline[j];
	  data->inf_timeline_std_noext[j]+=(double)stats.inf_timeline[j]*stats.inf_timeline[j];
	  stats.totinf_timeline[j]+=stats.totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totinf_timeline_mean_noext[j]+=stats.totinf_timeline[j];
	  data->totinf_timeline_std_noext[j]+=(double)stats.totinf_timeline[j]*stats.totinf_timeline[j];
	}
      }

      if(data->rec_ninfs) {
	stats.ngeninfs=(uint32_t*)realloc(stats.ngeninfs,stats.ninfs*sizeof(uint32_t));
	data->ngeninfs[initpath+i]=stats.ngeninfs;
	data->ninfs[initpath+i]=stats.ninfs;
	stats.ngeninfs=(uint32_t*)malloc(stats.nainfs*sizeof(uint32_t));
      }
    }
  }
  std_stats_free(&stats);
  branchsim_free(&sv);

  return NULL;
}
