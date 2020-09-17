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
  bool ninfhist=false;
  uint32_t npaths=10000;
  uint32_t nthreads=1;
  uint32_t nsetsperthread=(nthreads>1?100:1);
  uint32_t lmax=UINT32_MAX;
  uint32_t nimax=UINT32_MAX;
  int oout=STDOUT_FILENO;
  int eout=STDERR_FILENO;

  sim_pars_init(&pars);

  if(config(&pars, &ninfhist, &npaths, &nthreads, &nsetsperthread, &lmax, &nimax, &oout, &eout, nargs-1, args+1)) return 1;

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
  int tmaxnpersa=0;

  if(nthreads>1) {
    pthread_t* threads=(pthread_t*)malloc(nthreads*sizeof(pthread_t));

    for(t=nthreads-1; t>=0; --t){
      tdata[t].npathsperset=npathsperset;
      tdata[t].nsets=nsets;
      tdata[t].tnpersa=tdata[t].npers=npers;
      tdata[t].lmax=lmax;
      tdata[t].nimax=nimax;
      tdata[t].pars=&pars;
      tdata[t].set=&set;
      tdata[t].ngeninfs=NULL;
      tdata[t].ninfbins=0;
      //tdata[t].r = gsl_rng_alloc(gsl_rng_taus2);
      tdata[t].r = gsl_rng_alloc(rngstream_gsl);
      tdata[t].rec_ninfs=ninfhist;
      pthread_create(threads+t,NULL,simthread,tdata+t);
    }

    pthread_join(threads[0],NULL);
    gsl_rng_free(tdata[0].r);
    int tp;
    uint32_t shift;

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

      if(tdata[t].tnpersa > tdata[tmaxnpersa].tnpersa) {
	tp=tmaxnpersa;
	tmaxnpersa=t;

      } else tp=t;
      shift=tdata[tmaxnpersa].tnpersa-tdata[tp].tnpersa;

      for(j=shift; j<tdata[tmaxnpersa].tnpersa; ++j) {
        tdata[tmaxnpersa].inf_timeline_mean_ext[j]+=tdata[tp].inf_timeline_mean_ext[j-shift];
        tdata[tmaxnpersa].inf_timeline_std_ext[j]+=tdata[tp].inf_timeline_std_ext[j-shift];
        tdata[tmaxnpersa].totinf_timeline_mean_ext[j]+=tdata[tp].totinf_timeline_mean_ext[j-shift];
        tdata[tmaxnpersa].totinf_timeline_std_ext[j]+=tdata[tp].totinf_timeline_std_ext[j-shift];
        tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j]+=tdata[tp].totmainatt_timeline_mean_ext[j-shift];
        tdata[tmaxnpersa].totmainatt_timeline_std_ext[j]+=tdata[tp].totmainatt_timeline_std_ext[j-shift];

        tdata[tmaxnpersa].inf_timeline_mean_noext[j]+=tdata[tp].inf_timeline_mean_noext[j-shift];
        tdata[tmaxnpersa].inf_timeline_std_noext[j]+=tdata[tp].inf_timeline_std_noext[j-shift];
        tdata[tmaxnpersa].totinf_timeline_mean_noext[j]+=tdata[tp].totinf_timeline_mean_noext[j-shift];
        tdata[tmaxnpersa].totinf_timeline_std_noext[j]+=tdata[tp].totinf_timeline_std_noext[j-shift];
        tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j]+=tdata[tp].totmainatt_timeline_mean_noext[j-shift];
        tdata[tmaxnpersa].totmainatt_timeline_std_noext[j]+=tdata[tp].totmainatt_timeline_std_noext[j-shift];
      }

      if(tdata[t].nimaxedoutmintimeindex < tdata[0].nimaxedoutmintimeindex) tdata[0].nimaxedoutmintimeindex = tdata[t].nimaxedoutmintimeindex;
    }
    free(threads);

  } else {
    tdata[0].npathsperset=npathsperset;
    tdata[0].nsets=nsets;
    tdata[0].tnpersa=tdata[0].npers=npers;
    tdata[0].lmax=lmax;
    tdata[0].nimax=nimax;
    tdata[0].pars=&pars;
    tdata[0].set=&set;
    tdata[0].ngeninfs=NULL;
    tdata[0].ninfbins=0;
    //tdata[0].r = gsl_rng_alloc(gsl_rng_taus2);
    tdata[0].r = gsl_rng_alloc(rngstream_gsl);
    tdata[0].rec_ninfs=ninfhist;
    simthread(tdata);
    gsl_rng_free(tdata[0].r);
  }

  const double ninf=tdata[tmaxnpersa].totinf_timeline_mean_ext[tdata[tmaxnpersa].tnpersa-2]+tdata[tmaxnpersa].totinf_timeline_mean_noext[tdata[tmaxnpersa].tnpersa-2];
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

  double inf_timeline_mean[tdata[tmaxnpersa].tnpersa];
  double inf_timeline_std[tdata[tmaxnpersa].tnpersa];
  double totinf_timeline_mean[tdata[tmaxnpersa].tnpersa];
  double totinf_timeline_std[tdata[tmaxnpersa].tnpersa];
  double totmainatt_timeline_mean[tdata[tmaxnpersa].tnpersa];
  double totmainatt_timeline_std[tdata[tmaxnpersa].tnpersa];

  for(j=tdata[tmaxnpersa].tnpersa-1; j>=0; --j) {
    inf_timeline_mean[j]=tdata[tmaxnpersa].inf_timeline_mean_ext[j]+tdata[tmaxnpersa].inf_timeline_mean_noext[j];
    inf_timeline_std[j]=tdata[tmaxnpersa].inf_timeline_std_ext[j]+tdata[tmaxnpersa].inf_timeline_std_noext[j];
    totinf_timeline_mean[j]=tdata[tmaxnpersa].totinf_timeline_mean_ext[j]+tdata[tmaxnpersa].totinf_timeline_mean_noext[j];
    totinf_timeline_std[j]=tdata[tmaxnpersa].totinf_timeline_std_ext[j]+tdata[tmaxnpersa].totinf_timeline_std_noext[j];
    totmainatt_timeline_mean[j]=tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j]+tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j];
    totmainatt_timeline_std[j]=tdata[tmaxnpersa].totmainatt_timeline_std_ext[j]+tdata[tmaxnpersa].totmainatt_timeline_std_noext[j];

    inf_timeline_mean[j]/=npaths;
    inf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(inf_timeline_std[j]/npaths-inf_timeline_mean[j]*inf_timeline_mean[j]));
    totinf_timeline_mean[j]/=npaths;
    totinf_timeline_std[j]=sqrt(npaths/(npaths-1.)*(totinf_timeline_std[j]/npaths-totinf_timeline_mean[j]*totinf_timeline_mean[j]));
    totmainatt_timeline_mean[j]/=npaths;
    totmainatt_timeline_std[j]=sqrt(npaths/(npaths-1.)*(totmainatt_timeline_std[j]/npaths-totmainatt_timeline_mean[j]*totmainatt_timeline_mean[j]));

    tdata[tmaxnpersa].inf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpersa].inf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpersa].inf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpersa].inf_timeline_mean_ext[j]*tdata[tmaxnpersa].inf_timeline_mean_ext[j]));
    tdata[tmaxnpersa].totinf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpersa].totinf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpersa].totinf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpersa].totinf_timeline_mean_ext[j]*tdata[tmaxnpersa].totinf_timeline_mean_ext[j]));
    tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpersa].totmainatt_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpersa].totmainatt_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j]*tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j]));

    tdata[tmaxnpersa].inf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpersa].inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpersa].inf_timeline_std_noext[j]/nnoe-tdata[tmaxnpersa].inf_timeline_mean_noext[j]*tdata[tmaxnpersa].inf_timeline_mean_noext[j]));
    tdata[tmaxnpersa].totinf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpersa].totinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpersa].totinf_timeline_std_noext[j]/nnoe-tdata[tmaxnpersa].totinf_timeline_mean_noext[j]*tdata[tmaxnpersa].totinf_timeline_mean_noext[j]));
    tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpersa].totmainatt_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpersa].totmainatt_timeline_std_noext[j]/nnoe-tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j]*tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j]));
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

  int shift=tdata[tmaxnpersa].tnpersa-npers;
  printf("\nCurrent infection timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpersa].tnpersa; ++j) printf("%3i: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",j-shift,tdata[tmaxnpersa].inf_timeline_mean_ext[j],tdata[tmaxnpersa].inf_timeline_std_ext[j],tdata[tmaxnpersa].inf_timeline_mean_noext[j],tdata[tmaxnpersa].inf_timeline_std_noext[j],inf_timeline_mean[j],inf_timeline_std[j],(j<tdata[0].nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  printf("\nTotal infections timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpersa].tnpersa; ++j) printf("%3i: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",j-shift,tdata[tmaxnpersa].totinf_timeline_mean_ext[j],tdata[tmaxnpersa].totinf_timeline_std_ext[j],tdata[tmaxnpersa].totinf_timeline_mean_noext[j],tdata[tmaxnpersa].totinf_timeline_std_noext[j],totinf_timeline_mean[j],totinf_timeline_std[j],(j<tdata[0].nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  printf("\nTotal number of direct attendees timeline from individuals whose ended infectious communicable period was the main period, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpersa].tnpersa; ++j) printf("%3i: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",j-shift,tdata[tmaxnpersa].totmainatt_timeline_mean_ext[j],tdata[tmaxnpersa].totmainatt_timeline_std_ext[j],tdata[tmaxnpersa].totmainatt_timeline_mean_noext[j],tdata[tmaxnpersa].totmainatt_timeline_std_noext[j],totmainatt_timeline_mean[j],totmainatt_timeline_std[j],(j<tdata[0].nimaxedoutmintimeindex?"":" (nimax reached, biased)"));

  if(ninfhist) {
    uint32_t maxninfnbins=tdata[nthreads-1].ninfbins;
    uint32_t b;
    int tmaxninfnbins=nthreads-1;
    int tp;

    for(t=nthreads-2; t>=0; --t) {

      if(tdata[t].ninfbins > maxninfnbins) {
	maxninfnbins=tdata[t].ninfbins;
	tp=tmaxninfnbins;
	tmaxninfnbins=t;

      } else tp=t;

      for(b=0; b<tdata[tp].ninfbins; ++b) tdata[tmaxninfnbins].ngeninfs[b]+=tdata[tp].ngeninfs[b];
      free(tdata[tp].ngeninfs);
    }

    printf("\nDistribution of number of generated infections per infectious individual:\n");
    printf(" n inf\t               count\n");
    
    for(b=0; b<maxninfnbins; ++b) if(tdata[tmaxninfnbins].ngeninfs[b] > 0) printf("%6" PRIu32 "\t%20" PRIu64 "\n",b,tdata[tmaxninfnbins].ngeninfs[b]);
    free(tdata[tmaxninfnbins].ngeninfs);
  }

  for(t=nthreads-1; t>=0; --t) {
    free(tdata[t].inf_timeline_mean_ext);
    free(tdata[t].inf_timeline_std_ext);
    free(tdata[t].inf_timeline_mean_noext);
    free(tdata[t].inf_timeline_std_noext);
    free(tdata[t].totinf_timeline_mean_ext);
    free(tdata[t].totinf_timeline_std_ext);
    free(tdata[t].totinf_timeline_mean_noext);
    free(tdata[t].totinf_timeline_std_noext);
    free(tdata[t].totmainatt_timeline_mean_ext);
    free(tdata[t].totmainatt_timeline_std_ext);
    free(tdata[t].totmainatt_timeline_mean_noext);
    free(tdata[t].totmainatt_timeline_std_noext);
  }
  free(tdata);

  fflush(stdout);
  fflush(stderr);
  close(oout);
  close(eout);

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

  data->inf_timeline_mean_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->inf_timeline_std_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->inf_timeline_mean_noext=(double*)malloc(data->tnpersa*sizeof(double));
  data->inf_timeline_std_noext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totinf_timeline_mean_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totinf_timeline_std_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totinf_timeline_mean_noext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totinf_timeline_std_noext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totmainatt_timeline_mean_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totmainatt_timeline_std_ext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totmainatt_timeline_mean_noext=(double*)malloc(data->tnpersa*sizeof(double));
  data->totmainatt_timeline_std_noext=(double*)malloc(data->tnpersa*sizeof(double));

  memset(data->inf_timeline_mean_ext, 0, data->tnpersa*sizeof(double));
  memset(data->inf_timeline_std_ext, 0, data->tnpersa*sizeof(double));
  memset(data->inf_timeline_mean_noext, 0, data->tnpersa*sizeof(double));
  memset(data->inf_timeline_std_noext, 0, data->tnpersa*sizeof(double));
  memset(data->totinf_timeline_mean_ext, 0, data->tnpersa*sizeof(double));
  memset(data->totinf_timeline_std_ext, 0, data->tnpersa*sizeof(double));
  memset(data->totinf_timeline_mean_noext, 0, data->tnpersa*sizeof(double));
  memset(data->totinf_timeline_std_noext, 0, data->tnpersa*sizeof(double));
  memset(data->totmainatt_timeline_mean_ext, 0, data->tnpersa*sizeof(double));
  memset(data->totmainatt_timeline_std_ext, 0, data->tnpersa*sizeof(double));
  memset(data->totmainatt_timeline_mean_noext, 0, data->tnpersa*sizeof(double));
  memset(data->totmainatt_timeline_std_noext, 0, data->tnpersa*sizeof(double));

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
  sim_set_new_pri_inf_proc_func(&sv, std_stats_new_pri_inf);
  sim_set_new_inf_proc_func(&sv, std_stats_new_inf);

  if(data->rec_ninfs) {
    sim_set_end_inf_proc_func(&sv, std_stats_end_inf_rec_ninfs);
    sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf_rec_ninfs);

  } else {
    sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
    sim_set_inf_proc_noevent_func(&sv, std_stats_noevent_inf);
  }

  branchsim_init(&sv);

  if(data->rec_ninfs) std_stats_init(&sv, &data->ngeninfs, &data->ninfbins);

  else std_stats_init(&sv, NULL, NULL);
  
  stats.lmax=data->lmax;
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
      uint32_t* abs_inf_timeline=stats.inf_timeline-stats.timelineshift;
      uint32_t* abs_totinf_timeline=stats.totinf_timeline-stats.timelineshift;
      uint32_t* abs_totmainatt_timeline=stats.totmainatt_timeline-stats.timelineshift;

      if(stats.tnpersa > data->tnpersa) {
	uint32_t diff=stats.tnpersa-data->tnpersa;
	double* newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->inf_timeline_mean_ext,data->tnpersa*sizeof(double));
	free(data->inf_timeline_mean_ext);
	data->inf_timeline_mean_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->inf_timeline_std_ext,data->tnpersa*sizeof(double));
	free(data->inf_timeline_std_ext);
	data->inf_timeline_std_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totinf_timeline_mean_ext,data->tnpersa*sizeof(double));
	free(data->totinf_timeline_mean_ext);
	data->totinf_timeline_mean_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totinf_timeline_std_ext,data->tnpersa*sizeof(double));
	free(data->totinf_timeline_std_ext);
	data->totinf_timeline_std_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totmainatt_timeline_mean_ext,data->tnpersa*sizeof(double));
	free(data->totmainatt_timeline_mean_ext);
	data->totmainatt_timeline_mean_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totmainatt_timeline_std_ext,data->tnpersa*sizeof(double));
	free(data->totmainatt_timeline_std_ext);
	data->totmainatt_timeline_std_ext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->inf_timeline_mean_noext,data->tnpersa*sizeof(double));
	free(data->inf_timeline_mean_noext);
	data->inf_timeline_mean_noext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->inf_timeline_std_noext,data->tnpersa*sizeof(double));
	free(data->inf_timeline_std_noext);
	data->inf_timeline_std_noext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totinf_timeline_mean_noext,data->tnpersa*sizeof(double));
	free(data->totinf_timeline_mean_noext);
	data->totinf_timeline_mean_noext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totinf_timeline_std_noext,data->tnpersa*sizeof(double));
	free(data->totinf_timeline_std_noext);
	data->totinf_timeline_std_noext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totmainatt_timeline_mean_noext,data->tnpersa*sizeof(double));
	free(data->totmainatt_timeline_mean_noext);
	data->totmainatt_timeline_mean_noext=newarray;

	newarray=(double*)malloc(stats.tnpersa*sizeof(double));
	memset(newarray,0,diff*sizeof(double));
	memcpy(newarray+diff,data->totmainatt_timeline_std_noext,data->tnpersa*sizeof(double));
	free(data->totmainatt_timeline_std_noext);
	data->totmainatt_timeline_std_noext=newarray;

	data->tnpersa=stats.tnpersa;
      }

      if(stats.extinction) {
	data->pe+=stats.extinction;
	data->te_mean+=stats.extinction_time;
	data->te_std+=stats.extinction_time*stats.extinction_time;

	data->inf_timeline_mean_ext[0]+=abs_inf_timeline[0];
	data->inf_timeline_std_ext[0]+=(double)abs_inf_timeline[0]*abs_inf_timeline[0];
	data->totinf_timeline_mean_ext[0]+=abs_totinf_timeline[0];
	data->totinf_timeline_std_ext[0]+=(double)abs_totinf_timeline[0]*abs_totinf_timeline[0];
	data->totmainatt_timeline_mean_ext[0]+=abs_totmainatt_timeline[0];
	data->totmainatt_timeline_std_ext[0]+=(double)abs_totmainatt_timeline[0]*abs_totmainatt_timeline[0];

	for(j=1; j<stats.tnpersa; ++j) {
	  data->inf_timeline_mean_ext[j]+=abs_inf_timeline[j];
	  data->inf_timeline_std_ext[j]+=(double)abs_inf_timeline[j]*abs_inf_timeline[j];
	  abs_totinf_timeline[j]+=abs_totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totinf_timeline_mean_ext[j]+=abs_totinf_timeline[j];
	  data->totinf_timeline_std_ext[j]+=(double)abs_totinf_timeline[j]*abs_totinf_timeline[j];
	  abs_totmainatt_timeline[j]+=abs_totmainatt_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totmainatt_timeline_mean_ext[j]+=abs_totmainatt_timeline[j];
	  data->totmainatt_timeline_std_ext[j]+=(double)abs_totmainatt_timeline[j]*abs_totmainatt_timeline[j];
	}

      } else {

	if(stats.nimaxedoutmintimeindex < data->nimaxedoutmintimeindex) data->nimaxedoutmintimeindex=stats.nimaxedoutmintimeindex;
	data->inf_timeline_mean_noext[0]+=abs_inf_timeline[0];
	data->inf_timeline_std_noext[0]+=(double)abs_inf_timeline[0]*abs_inf_timeline[0];
	data->totinf_timeline_mean_noext[0]+=abs_totinf_timeline[0];
	data->totinf_timeline_std_noext[0]+=(double)abs_totinf_timeline[0]*abs_totinf_timeline[0];
	data->totmainatt_timeline_mean_noext[0]+=abs_totmainatt_timeline[0];
	data->totmainatt_timeline_std_noext[0]+=(double)abs_totmainatt_timeline[0]*abs_totmainatt_timeline[0];

	for(j=1; j<stats.tnpersa; ++j) {
	  data->inf_timeline_mean_noext[j]+=abs_inf_timeline[j];
	  data->inf_timeline_std_noext[j]+=(double)abs_inf_timeline[j]*abs_inf_timeline[j];
	  abs_totinf_timeline[j]+=abs_totinf_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totinf_timeline_mean_noext[j]+=abs_totinf_timeline[j];
	  data->totinf_timeline_std_noext[j]+=(double)abs_totinf_timeline[j]*abs_totinf_timeline[j];
	  abs_totmainatt_timeline[j]+=abs_totmainatt_timeline[j-1]; //Warning: Operation on simulation output variable
	  data->totmainatt_timeline_mean_noext[j]+=abs_totmainatt_timeline[j];
	  data->totmainatt_timeline_std_noext[j]+=(double)abs_totmainatt_timeline[j]*abs_totmainatt_timeline[j];
	}
      }
    }
  }
  std_stats_free(&stats);
  branchsim_free(&sv);

  return NULL;
}
