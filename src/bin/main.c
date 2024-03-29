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
  config_pars cp={.ninfhist=false, .npaths=10000, .lmax=UINT32_MAX, .nbinsperunit=1, .nimax=UINT32_MAX, .npostestmax=UINT32_MAX, .npostestmaxnunits=1, .nthreads=1, .stream=0, .tloutbufsize=10, .tlout=0, 
#ifdef CT_OUTPUT
    .ctoutbufsize=10, 
    .ctout=0, 
#endif
    .oout=STDOUT_FILENO, .eout=STDERR_FILENO};
  cp.nsetsperthread=(cp.nthreads>1?100:1);
  pthread_mutex_t tlflock;
#ifdef CT_OUTPUT
  pthread_mutex_t ctflock;
#endif

  sim_pars_init(&cp.pars);

  if(config(&cp, nargs-1, args+1)) return 1;

  if(model_solve_pars(&cp.pars)) return 1;

  gsl_rng_env_setup();

  rng_skipstreams(cp.nthreads*cp.stream);

  if(model_pars_check(&cp.pars)) {
    fprintf(stderr,"%s: Error: While verifying the validity of the simulation parameters.\n",args[0]);
    return 1;
  }

  if(cp.tlout) {
    pthread_mutex_init(&tlflock, NULL);
    uint32_t ubuf=htole32(cp.nbinsperunit*cp.pars.tmax);

    if(write(cp.tlout,&ubuf,4)!=4) {
      perror("tlout");
      return 1;
    }
#ifdef SEC_INF_TIMELINES
    ubuf=cp.pars.timetype|((!isnan(cp.pars.tdeltat))<<3)|16;
#else
    ubuf=cp.pars.timetype|((!isnan(cp.pars.tdeltat))<<3);
#endif

    if(write(cp.tlout,&ubuf,1)!=1) {
      perror("tlout");
      return 1;
    }
  }

#ifdef CT_OUTPUT
  if(cp.ctout) {
    pthread_mutex_init(&ctflock, NULL);
  }
#endif

  thread_data* tdata=(thread_data*)malloc(cp.nthreads*sizeof(thread_data));
  int j;
  const uint32_t nsets=cp.nthreads*cp.nsetsperthread;
  const double npathsperset=((double)cp.npaths)/nsets;
  const uint32_t npers=cp.nbinsperunit*cp.pars.tmax;
  volatile uint32_t set=cp.nthreads;
  int t;
  int tmaxnpers=0;

  if(cp.nthreads>1) {
    pthread_t* threads=(pthread_t*)malloc(cp.nthreads*sizeof(pthread_t));

    for(t=cp.nthreads-1; t>=0; --t){
      tdata[t].cp=&cp;
      tdata[t].npathsperset=npathsperset;
      tdata[t].id=t;
      tdata[t].nsets=nsets;
      tdata[t].nbinsperunit=cp.nbinsperunit;
      tdata[t].npers=npers;
      tdata[t].set=&set;
      //tdata[t].r = gsl_rng_alloc(gsl_rng_taus2);
      tdata[t].r = gsl_rng_alloc(rngstream_gsl);
      //rng_writestatefull((rng_stream*)tdata[t].r->state);
      tdata[t].tlflock = &tlflock;
#ifdef CT_OUTPUT
      tdata[t].ctflock = &ctflock;
#endif
      pthread_create(threads+t,NULL,simthread,tdata+t);
    }

    pthread_join(threads[0],NULL);
    gsl_rng_free(tdata[0].r);
    int tp;
    uint32_t maxper;
    int32_t ndiff, pdiff;

    for(t=1; t<cp.nthreads; ++t) {
      pthread_join(threads[t],NULL);
      gsl_rng_free(tdata[t].r);
      tdata[0].commper_mean+=tdata[t].commper_mean;
#ifdef NUMEVENTSSTATS
      tdata[0].nevents_mean+=tdata[t].nevents_mean;
#endif
      tdata[0].nnzpaths+=tdata[t].nnzpaths;
      tdata[0].pe+=tdata[t].pe;
      tdata[0].penz+=tdata[t].penz;
      tdata[0].pm+=tdata[t].pm;
      tdata[0].tenz_mean+=tdata[t].tenz_mean;
      tdata[0].tenz_std+=tdata[t].tenz_std;

      ndiff=tdata[t].tlppnnpers-tdata[tmaxnpers].tlppnnpers;
      pdiff=(int32_t)tdata[t].tlpptnvpers-tdata[tmaxnpers].tlpptnvpers-ndiff;
      //printf("Thread %i: nnpers: %i, tnvpers: %u\n",t,tdata[t].tlppnnpers,tdata[t].tlpptnvpers);
      //printf("Thread %i: nnpers: %i, tnvpers: %u\n",tmaxnpers,tdata[tmaxnpers].tlppnnpers,tdata[tmaxnpers].tlpptnvpers);
      //printf("ndiff: %i, pdiff: %i\n",ndiff,pdiff);

      if(ndiff>0) {
	tp=tmaxnpers;
	//printf("tmaxnpers: %i -> %i\n",tmaxnpers,t);
	tmaxnpers=t;

	if(pdiff<0) {
	  realloc_thread_timelines(tdata+t,ndiff,0);
	}
	ndiff*=-1;

      } else {
	tp=t;

	if(pdiff>0) realloc_thread_timelines(tdata+tmaxnpers,0,pdiff);
      }

      maxper=tdata[tp].tlpptnvpers-ndiff;
      //printf("Storing from %i to %i\n",-ndiff,maxper-1);

      for(j=-ndiff; j<maxper; ++j) {
        tdata[tmaxnpers].inf_timeline_mean_ext[j]+=tdata[tp].inf_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].inf_timeline_std_ext[j]+=tdata[tp].inf_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].newinf_timeline_mean_ext[j]+=tdata[tp].newinf_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].newinf_timeline_std_ext[j]+=tdata[tp].newinf_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].newpostest_timeline_mean_ext[j]+=tdata[tp].newpostest_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].newpostest_timeline_std_ext[j]+=tdata[tp].newpostest_timeline_std_ext[j+ndiff];
        #ifdef SEC_INF_TIMELINES
        tdata[tmaxnpers].secinf_timeline_mean_ext[j]+=tdata[tp].secinf_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].secinf_timeline_std_ext[j]+=tdata[tp].secinf_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].newsecinf_timeline_mean_ext[j]+=tdata[tp].newsecinf_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].newsecinf_timeline_std_ext[j]+=tdata[tp].newsecinf_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j]+=tdata[tp].newsecpostest_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].newsecpostest_timeline_std_ext[j]+=tdata[tp].newsecpostest_timeline_std_ext[j+ndiff];
        #endif
        tdata[tmaxnpers].reff_timeline_mean_ext[j]+=tdata[tp].reff_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].reff_timeline_std_ext[j]+=tdata[tp].reff_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].reff_timeline_n_ext[j]+=tdata[tp].reff_timeline_n_ext[j+ndiff];
	#ifdef OBSREFF_OUTPUT
        tdata[tmaxnpers].reffobs_timeline_mean_ext[j]+=tdata[tp].reffobs_timeline_mean_ext[j+ndiff];
        tdata[tmaxnpers].reffobs_timeline_std_ext[j]+=tdata[tp].reffobs_timeline_std_ext[j+ndiff];
        tdata[tmaxnpers].reffobs_timeline_n_ext[j]+=tdata[tp].reffobs_timeline_n_ext[j+ndiff];
        #endif

        tdata[tmaxnpers].inf_timeline_mean_noext[j]+=tdata[tp].inf_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].inf_timeline_std_noext[j]+=tdata[tp].inf_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].newinf_timeline_mean_noext[j]+=tdata[tp].newinf_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].newinf_timeline_std_noext[j]+=tdata[tp].newinf_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].newpostest_timeline_mean_noext[j]+=tdata[tp].newpostest_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].newpostest_timeline_std_noext[j]+=tdata[tp].newpostest_timeline_std_noext[j+ndiff];
        #ifdef SEC_INF_TIMELINES
        tdata[tmaxnpers].secinf_timeline_mean_noext[j]+=tdata[tp].secinf_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].secinf_timeline_std_noext[j]+=tdata[tp].secinf_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].newsecinf_timeline_mean_noext[j]+=tdata[tp].newsecinf_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].newsecinf_timeline_std_noext[j]+=tdata[tp].newsecinf_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j]+=tdata[tp].newsecpostest_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].newsecpostest_timeline_std_noext[j]+=tdata[tp].newsecpostest_timeline_std_noext[j+ndiff];
        #endif
        tdata[tmaxnpers].reff_timeline_mean_noext[j]+=tdata[tp].reff_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].reff_timeline_std_noext[j]+=tdata[tp].reff_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].reff_timeline_n_noext[j]+=tdata[tp].reff_timeline_n_noext[j+ndiff];
	#ifdef OBSREFF_OUTPUT
        tdata[tmaxnpers].reffobs_timeline_mean_noext[j]+=tdata[tp].reffobs_timeline_mean_noext[j+ndiff];
        tdata[tmaxnpers].reffobs_timeline_std_noext[j]+=tdata[tp].reffobs_timeline_std_noext[j+ndiff];
        tdata[tmaxnpers].reffobs_timeline_n_noext[j]+=tdata[tp].reffobs_timeline_n_noext[j+ndiff];
        #endif
      }
      free(tdata[tp].inf_timeline_mean_ext);
      free(tdata[tp].inf_timeline_std_ext);
      free(tdata[tp].inf_timeline_mean_noext);
      free(tdata[tp].inf_timeline_std_noext);
      free(tdata[tp].newinf_timeline_mean_ext);
      free(tdata[tp].newinf_timeline_std_ext);
      free(tdata[tp].newinf_timeline_mean_noext);
      free(tdata[tp].newinf_timeline_std_noext);
      free(tdata[tp].newpostest_timeline_mean_ext);
      free(tdata[tp].newpostest_timeline_std_ext);
      free(tdata[tp].newpostest_timeline_mean_noext);
      free(tdata[tp].newpostest_timeline_std_noext);
      #ifdef SEC_INF_TIMELINES
      free(tdata[tp].secinf_timeline_mean_ext);
      free(tdata[tp].secinf_timeline_std_ext);
      free(tdata[tp].secinf_timeline_mean_noext);
      free(tdata[tp].secinf_timeline_std_noext);
      free(tdata[tp].newsecinf_timeline_mean_ext);
      free(tdata[tp].newsecinf_timeline_std_ext);
      free(tdata[tp].newsecinf_timeline_mean_noext);
      free(tdata[tp].newsecinf_timeline_std_noext);
      free(tdata[tp].newsecpostest_timeline_mean_ext);
      free(tdata[tp].newsecpostest_timeline_std_ext);
      free(tdata[tp].newsecpostest_timeline_mean_noext);
      free(tdata[tp].newsecpostest_timeline_std_noext);
      #endif
      free(tdata[tp].reff_timeline_mean_ext);
      free(tdata[tp].reff_timeline_std_ext);
      free(tdata[tp].reff_timeline_n_ext);
      free(tdata[tp].reff_timeline_mean_noext);
      free(tdata[tp].reff_timeline_std_noext);
      free(tdata[tp].reff_timeline_n_noext);
      #ifdef OBSREFF_OUTPUT
      free(tdata[tp].reffobs_timeline_mean_ext);
      free(tdata[tp].reffobs_timeline_std_ext);
      free(tdata[tp].reffobs_timeline_n_ext);
      free(tdata[tp].reffobs_timeline_mean_noext);
      free(tdata[tp].reffobs_timeline_std_noext);
      free(tdata[tp].reffobs_timeline_n_noext);
      #endif

      if(tdata[t].maxedoutmintimeindex < tdata[0].maxedoutmintimeindex) tdata[0].maxedoutmintimeindex = tdata[t].maxedoutmintimeindex;
    }
    free(threads);

  } else {
    tdata[0].cp=&cp;
    tdata[0].npathsperset=npathsperset;
    tdata[0].id=0;
    tdata[0].nsets=nsets;
    tdata[0].nbinsperunit=cp.nbinsperunit;
    tdata[0].npers=npers;
    tdata[0].set=&set;
    //tdata[0].r = gsl_rng_alloc(gsl_rng_taus2);
    tdata[0].r = gsl_rng_alloc(rngstream_gsl);
    //rng_writestatefull((rng_stream*)tdata[0].r->state);
    tdata[0].tlflock = &tlflock;
#ifdef CT_OUTPUT
    tdata[0].ctflock = &ctflock;
#endif
    simthread(tdata);
    gsl_rng_free(tdata[0].r);
  }

  if(cp.tlout) {
    close(cp.tlout);
    pthread_mutex_destroy(&tlflock);
  }

#ifdef CT_OUTPUT
  if(cp.ctout) {
    close(cp.ctout);
    pthread_mutex_destroy(&ctflock);
  }
#endif

  //printf("R sum is %f\n",tdata[0].r_mean);
  //printf("Total number of infectious is %f\n",ninf);
//#ifdef NUMEVENTSSTATS
//  const double ninf_per_event_mean=tdata[0].r_mean/tdata[0].nevents_mean;
//#endif
  const double nnoe=cp.npaths-tdata[0].pe;

  double inf_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double inf_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  double newinf_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double newinf_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  double newpostest_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double newpostest_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  #ifdef SEC_INF_TIMELINES
  double secinf_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double secinf_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  double newsecinf_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double newsecinf_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  double newsecpostest_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double newsecpostest_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  #endif
  double reff_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double reff_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  uint64_t reff_timeline_n[tdata[tmaxnpers].tlpptnvpers];
  double reff_mean_ext=0, reff_mean_noext=0, reff_mean;
  double reff_std_ext=0, reff_std_noext=0, reff_std;
  uint64_t reff_ext_n=0, reff_noext_n=0, reff_n;
  #ifdef OBSREFF_OUTPUT
  double reffobs_timeline_mean[tdata[tmaxnpers].tlpptnvpers];
  double reffobs_timeline_std[tdata[tmaxnpers].tlpptnvpers];
  uint64_t reffobs_timeline_n[tdata[tmaxnpers].tlpptnvpers];
  double reffobs_mean_ext=0, reffobs_mean_noext=0, reffobs_mean=0;
  double reffobs_std_ext=0, reffobs_std_noext=0, reffobs_std=0;
  uint64_t reffobs_ext_n=0, reffobs_noext_n=0, reffobs_n=0;
  #endif

  for(j=tdata[tmaxnpers].tlpptnvpers-1; j>=0; --j) {
    inf_timeline_mean[j]=tdata[tmaxnpers].inf_timeline_mean_ext[j]+tdata[tmaxnpers].inf_timeline_mean_noext[j];
    inf_timeline_std[j]=tdata[tmaxnpers].inf_timeline_std_ext[j]+tdata[tmaxnpers].inf_timeline_std_noext[j];
    newinf_timeline_mean[j]=tdata[tmaxnpers].newinf_timeline_mean_ext[j]+tdata[tmaxnpers].newinf_timeline_mean_noext[j];
    newinf_timeline_std[j]=tdata[tmaxnpers].newinf_timeline_std_ext[j]+tdata[tmaxnpers].newinf_timeline_std_noext[j];
    newpostest_timeline_mean[j]=tdata[tmaxnpers].newpostest_timeline_mean_ext[j]+tdata[tmaxnpers].newpostest_timeline_mean_noext[j];
    newpostest_timeline_std[j]=tdata[tmaxnpers].newpostest_timeline_std_ext[j]+tdata[tmaxnpers].newpostest_timeline_std_noext[j];
    #ifdef SEC_INF_TIMELINES
    secinf_timeline_mean[j]=tdata[tmaxnpers].secinf_timeline_mean_ext[j]+tdata[tmaxnpers].secinf_timeline_mean_noext[j];
    secinf_timeline_std[j]=tdata[tmaxnpers].secinf_timeline_std_ext[j]+tdata[tmaxnpers].secinf_timeline_std_noext[j];
    newsecinf_timeline_mean[j]=tdata[tmaxnpers].newsecinf_timeline_mean_ext[j]+tdata[tmaxnpers].newsecinf_timeline_mean_noext[j];
    newsecinf_timeline_std[j]=tdata[tmaxnpers].newsecinf_timeline_std_ext[j]+tdata[tmaxnpers].newsecinf_timeline_std_noext[j];
    newsecpostest_timeline_mean[j]=tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j]+tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j];
    newsecpostest_timeline_std[j]=tdata[tmaxnpers].newsecpostest_timeline_std_ext[j]+tdata[tmaxnpers].newsecpostest_timeline_std_noext[j];
    #endif
    reff_timeline_mean[j]=tdata[tmaxnpers].reff_timeline_mean_ext[j]+tdata[tmaxnpers].reff_timeline_mean_noext[j];
    reff_timeline_std[j]=tdata[tmaxnpers].reff_timeline_std_ext[j]+tdata[tmaxnpers].reff_timeline_std_noext[j];
    reff_timeline_n[j]=tdata[tmaxnpers].reff_timeline_n_ext[j]+tdata[tmaxnpers].reff_timeline_n_noext[j];
    reff_ext_n+=tdata[tmaxnpers].reff_timeline_n_ext[j];
    reff_noext_n+=tdata[tmaxnpers].reff_timeline_n_noext[j];
    reff_mean_ext+=tdata[tmaxnpers].reff_timeline_mean_ext[j];
    reff_mean_noext+=tdata[tmaxnpers].reff_timeline_mean_noext[j];
    reff_std_ext+=tdata[tmaxnpers].reff_timeline_std_ext[j];
    reff_std_noext+=tdata[tmaxnpers].reff_timeline_std_noext[j];
    #ifdef OBSREFF_OUTPUT
    reffobs_timeline_mean[j]=tdata[tmaxnpers].reffobs_timeline_mean_ext[j]+tdata[tmaxnpers].reffobs_timeline_mean_noext[j];
    reffobs_timeline_std[j]=tdata[tmaxnpers].reffobs_timeline_std_ext[j]+tdata[tmaxnpers].reffobs_timeline_std_noext[j];
    reffobs_timeline_n[j]=tdata[tmaxnpers].reffobs_timeline_n_ext[j]+tdata[tmaxnpers].reffobs_timeline_n_noext[j];
    reffobs_ext_n+=tdata[tmaxnpers].reffobs_timeline_n_ext[j];
    reffobs_noext_n+=tdata[tmaxnpers].reffobs_timeline_n_noext[j];
    reffobs_mean_ext+=tdata[tmaxnpers].reffobs_timeline_mean_ext[j];
    reffobs_mean_noext+=tdata[tmaxnpers].reffobs_timeline_mean_noext[j];
    reffobs_std_ext+=tdata[tmaxnpers].reffobs_timeline_std_ext[j];
    reffobs_std_noext+=tdata[tmaxnpers].reffobs_timeline_std_noext[j];
    #endif

    inf_timeline_mean[j]/=cp.npaths;
    inf_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(inf_timeline_std[j]/cp.npaths-inf_timeline_mean[j]*inf_timeline_mean[j]));
    newinf_timeline_mean[j]/=cp.npaths;
    newinf_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(newinf_timeline_std[j]/cp.npaths-newinf_timeline_mean[j]*newinf_timeline_mean[j]));
    newpostest_timeline_mean[j]/=cp.npaths;
    newpostest_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(newpostest_timeline_std[j]/cp.npaths-newpostest_timeline_mean[j]*newpostest_timeline_mean[j]));
    #ifdef SEC_INF_TIMELINES
    secinf_timeline_mean[j]/=cp.npaths;
    secinf_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(secinf_timeline_std[j]/cp.npaths-secinf_timeline_mean[j]*secinf_timeline_mean[j]));
    newsecinf_timeline_mean[j]/=cp.npaths;
    newsecinf_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(newsecinf_timeline_std[j]/cp.npaths-newsecinf_timeline_mean[j]*newsecinf_timeline_mean[j]));
    newsecpostest_timeline_mean[j]/=cp.npaths;
    newsecpostest_timeline_std[j]=sqrt(cp.npaths/(cp.npaths-1.)*(newsecpostest_timeline_std[j]/cp.npaths-newsecpostest_timeline_mean[j]*newsecpostest_timeline_mean[j]));
    #endif
    reff_timeline_mean[j]=(reff_timeline_n[j]?reff_timeline_mean[j]/reff_timeline_n[j]:NAN);
    reff_timeline_std[j]=(reff_timeline_n[j]>1?sqrt(reff_timeline_n[j]/(reff_timeline_n[j]-1.)*(reff_timeline_std[j]/reff_timeline_n[j]-reff_timeline_mean[j]*reff_timeline_mean[j])):(reff_timeline_n[j]?INFINITY:NAN));
    #ifdef OBSREFF_OUTPUT
    reffobs_timeline_mean[j]=(reffobs_timeline_n[j]?reffobs_timeline_mean[j]/reffobs_timeline_n[j]:NAN);
    reffobs_timeline_std[j]=(reffobs_timeline_n[j]>1?sqrt(reffobs_timeline_n[j]/(reffobs_timeline_n[j]-1.)*(reffobs_timeline_std[j]/reffobs_timeline_n[j]-reffobs_timeline_mean[j]*reffobs_timeline_mean[j])):(reffobs_timeline_n[j]?INFINITY:NAN));
    #endif

    tdata[tmaxnpers].inf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].inf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].inf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].inf_timeline_mean_ext[j]*tdata[tmaxnpers].inf_timeline_mean_ext[j]));
    tdata[tmaxnpers].newinf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].newinf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].newinf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].newinf_timeline_mean_ext[j]*tdata[tmaxnpers].newinf_timeline_mean_ext[j]));
    tdata[tmaxnpers].newpostest_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].newpostest_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].newpostest_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].newpostest_timeline_mean_ext[j]*tdata[tmaxnpers].newpostest_timeline_mean_ext[j]));
    #ifdef SEC_INF_TIMELINES
    tdata[tmaxnpers].secinf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].secinf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].secinf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].secinf_timeline_mean_ext[j]*tdata[tmaxnpers].secinf_timeline_mean_ext[j]));
    tdata[tmaxnpers].newsecinf_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].newsecinf_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].newsecinf_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].newsecinf_timeline_mean_ext[j]*tdata[tmaxnpers].newsecinf_timeline_mean_ext[j]));
    tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j]/=tdata[0].pe;
    tdata[tmaxnpers].newsecpostest_timeline_std_ext[j]=sqrt(tdata[0].pe/(tdata[0].pe-1.)*(tdata[tmaxnpers].newsecpostest_timeline_std_ext[j]/tdata[0].pe-tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j]*tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j]));
    #endif
    tdata[tmaxnpers].reff_timeline_mean_ext[j]=(tdata[tmaxnpers].reff_timeline_n_ext[j]?tdata[tmaxnpers].reff_timeline_mean_ext[j]/tdata[tmaxnpers].reff_timeline_n_ext[j]:NAN);
    tdata[tmaxnpers].reff_timeline_std_ext[j]=(tdata[tmaxnpers].reff_timeline_n_ext[j]>1?sqrt(tdata[tmaxnpers].reff_timeline_n_ext[j]/(tdata[tmaxnpers].reff_timeline_n_ext[j]-1.)*(tdata[tmaxnpers].reff_timeline_std_ext[j]/tdata[tmaxnpers].reff_timeline_n_ext[j]-tdata[tmaxnpers].reff_timeline_mean_ext[j]*tdata[tmaxnpers].reff_timeline_mean_ext[j])):(tdata[tmaxnpers].reff_timeline_n_ext[j]?INFINITY:NAN));
    #ifdef OBSREFF_OUTPUT
    tdata[tmaxnpers].reffobs_timeline_mean_ext[j]=(tdata[tmaxnpers].reffobs_timeline_n_ext[j]?tdata[tmaxnpers].reffobs_timeline_mean_ext[j]/tdata[tmaxnpers].reffobs_timeline_n_ext[j]:NAN);
    tdata[tmaxnpers].reffobs_timeline_std_ext[j]=(tdata[tmaxnpers].reffobs_timeline_n_ext[j]>1?sqrt(tdata[tmaxnpers].reffobs_timeline_n_ext[j]/(tdata[tmaxnpers].reffobs_timeline_n_ext[j]-1.)*(tdata[tmaxnpers].reffobs_timeline_std_ext[j]/tdata[tmaxnpers].reffobs_timeline_n_ext[j]-tdata[tmaxnpers].reffobs_timeline_mean_ext[j]*tdata[tmaxnpers].reffobs_timeline_mean_ext[j])):(tdata[tmaxnpers].reffobs_timeline_n_ext[j]?INFINITY:NAN));
    #endif

    tdata[tmaxnpers].inf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].inf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].inf_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].inf_timeline_mean_noext[j]*tdata[tmaxnpers].inf_timeline_mean_noext[j]));
    tdata[tmaxnpers].newinf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].newinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].newinf_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].newinf_timeline_mean_noext[j]*tdata[tmaxnpers].newinf_timeline_mean_noext[j]));
    tdata[tmaxnpers].newpostest_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].newpostest_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].newpostest_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].newpostest_timeline_mean_noext[j]*tdata[tmaxnpers].newpostest_timeline_mean_noext[j]));
    #ifdef SEC_INF_TIMELINES
    tdata[tmaxnpers].secinf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].secinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].secinf_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].secinf_timeline_mean_noext[j]*tdata[tmaxnpers].secinf_timeline_mean_noext[j]));
    tdata[tmaxnpers].newsecinf_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].newsecinf_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].newsecinf_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].newsecinf_timeline_mean_noext[j]*tdata[tmaxnpers].newsecinf_timeline_mean_noext[j]));
    tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j]/=nnoe;
    tdata[tmaxnpers].newsecpostest_timeline_std_noext[j]=sqrt(nnoe/(nnoe-1.)*(tdata[tmaxnpers].newsecpostest_timeline_std_noext[j]/nnoe-tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j]*tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j]));
    #endif
    tdata[tmaxnpers].reff_timeline_mean_noext[j]=(tdata[tmaxnpers].reff_timeline_n_noext[j]?tdata[tmaxnpers].reff_timeline_mean_noext[j]/tdata[tmaxnpers].reff_timeline_n_noext[j]:NAN);
    tdata[tmaxnpers].reff_timeline_std_noext[j]=(tdata[tmaxnpers].reff_timeline_n_noext[j]>1?sqrt(tdata[tmaxnpers].reff_timeline_n_noext[j]/(tdata[tmaxnpers].reff_timeline_n_noext[j]-1.)*(tdata[tmaxnpers].reff_timeline_std_noext[j]/tdata[tmaxnpers].reff_timeline_n_noext[j]-tdata[tmaxnpers].reff_timeline_mean_noext[j]*tdata[tmaxnpers].reff_timeline_mean_noext[j])):(tdata[tmaxnpers].reff_timeline_n_noext[j]?INFINITY:NAN));
    #ifdef OBSREFF_OUTPUT
    tdata[tmaxnpers].reffobs_timeline_mean_noext[j]=(tdata[tmaxnpers].reffobs_timeline_n_noext[j]?tdata[tmaxnpers].reffobs_timeline_mean_noext[j]/tdata[tmaxnpers].reffobs_timeline_n_noext[j]:NAN);
    tdata[tmaxnpers].reffobs_timeline_std_noext[j]=(tdata[tmaxnpers].reffobs_timeline_n_noext[j]>1?sqrt(tdata[tmaxnpers].reffobs_timeline_n_noext[j]/(tdata[tmaxnpers].reffobs_timeline_n_noext[j]-1.)*(tdata[tmaxnpers].reffobs_timeline_std_noext[j]/tdata[tmaxnpers].reffobs_timeline_n_noext[j]-tdata[tmaxnpers].reffobs_timeline_mean_noext[j]*tdata[tmaxnpers].reffobs_timeline_mean_noext[j])):(tdata[tmaxnpers].reffobs_timeline_n_noext[j]?INFINITY:NAN));
    #endif
  }

  reff_mean=reff_mean_ext+reff_mean_noext;
  reff_std=reff_std_ext+reff_std_noext;
  reff_n=reff_ext_n+reff_noext_n;

  printf("r_mean %22.15e %" PRIu64 "\n",reff_mean,reff_n);

  reff_mean_ext=(reff_ext_n?reff_mean_ext/reff_ext_n:NAN);
  reff_std_ext=(reff_ext_n>1?sqrt(reff_ext_n/(reff_ext_n-1.)*(reff_std_ext/reff_ext_n-reff_mean_ext*reff_mean_ext)):(reff_ext_n?INFINITY:NAN));
  reff_mean_noext=(reff_noext_n?reff_mean_noext/reff_noext_n:NAN);
  reff_std_noext=(reff_noext_n>1?sqrt(reff_noext_n/(reff_noext_n-1.)*(reff_std_noext/reff_noext_n-reff_mean_noext*reff_mean_noext)):(reff_noext_n?INFINITY:NAN));
  reff_mean=(reff_n?reff_mean/reff_n:NAN);
  reff_std=(reff_n>1?sqrt(reff_n/(reff_n-1.)*(reff_std/reff_n-reff_mean*reff_mean)):(reff_n?INFINITY:NAN));
#ifdef OBSREFF_OUTPUT
  if(!isnan(cp.pars.tdeltat)) {
    reffobs_mean=reffobs_mean_ext+reffobs_mean_noext;
    reffobs_std=reffobs_std_ext+reffobs_std_noext;
    reffobs_n=reffobs_ext_n+reffobs_noext_n;

    printf("robs_mean %22.15e %" PRIu64 "\n",reffobs_mean,reffobs_n);

    reffobs_mean_ext=(reffobs_ext_n?reffobs_mean_ext/reffobs_ext_n:NAN);
    reffobs_std_ext=(reffobs_ext_n>1?sqrt(reffobs_ext_n/(reffobs_ext_n-1.)*(reffobs_std_ext/reffobs_ext_n-reffobs_mean_ext*reffobs_mean_ext)):(reffobs_ext_n?INFINITY:NAN));
    reffobs_mean_noext=(reffobs_noext_n?reffobs_mean_noext/reffobs_noext_n:NAN);
    reffobs_std_noext=(reffobs_noext_n>1?sqrt(reffobs_noext_n/(reffobs_noext_n-1.)*(reffobs_std_noext/reffobs_noext_n-reffobs_mean_noext*reffobs_mean_noext)):(reffobs_noext_n?INFINITY:NAN));
    reffobs_mean=(reffobs_n?reffobs_mean/reffobs_n:NAN);
    reffobs_std=(reffobs_n>1?sqrt(reffobs_n/(reffobs_n-1.)*(reffobs_std/reffobs_n-reffobs_mean*reffobs_mean)):(reffobs_n?INFINITY:NAN));
  }
#endif

  tdata[0].commper_mean/=reff_n;
#ifdef NUMEVENTSSTATS
  tdata[0].nevents_mean/=reff_n;
#endif
  tdata[0].pe/=cp.npaths;
  tdata[0].tenz_mean/=tdata[0].penz;
  tdata[0].tenz_std=sqrt(tdata[0].penz/(tdata[0].penz-1.)*(tdata[0].tenz_std/tdata[0].penz-tdata[0].tenz_mean*tdata[0].tenz_mean));
  tdata[0].penz/=tdata[0].nnzpaths;
  tdata[0].pm/=cp.npaths;

  printf("\nComputed simulation results:\n");
  printf("Mean R:\n\t    Extinct: %22.15e +/- %22.15e\n\tNon-extinct: %22.15e +/- %22.15e\n\t      Total: %22.15e +/- %22.15e\n",reff_mean_ext,reff_std_ext/sqrt(reff_ext_n),reff_mean_noext,reff_std_noext/sqrt(reff_noext_n),reff_mean,reff_std/sqrt(reff_n));
#ifdef OBSREFF_OUTPUT
  if(!isnan(cp.pars.tdeltat)) printf("Mean observed R:\n\t    Extinct: %22.15e +/- %22.15e\n\tNon-extinct: %22.15e +/- %22.15e\n\t      Total: %22.15e +/- %22.15e\n",reffobs_mean_ext,reffobs_std_ext/sqrt(reffobs_ext_n),reffobs_mean_noext,reffobs_std_noext/sqrt(reffobs_noext_n),reffobs_mean,reffobs_std/sqrt(reffobs_n));
#endif
  printf("Communicable period is %22.15e\n",tdata[0].commper_mean);
#ifdef NUMEVENTSSTATS
  printf("Number of events per infectious individual is %22.15e\n",tdata[0].nevents_mean);
  //printf("Number of infections per event is %22.15e\n",ninf_per_event_mean);
#endif
  printf("Probability of extinction and its statistical uncertainty: %22.15e +/- %22.15e%s\n",tdata[0].penz,sqrt(tdata[0].penz*(1.-tdata[0].penz)/(tdata[0].nnzpaths-1.)),(tdata[0].maxedoutmintimeindex<INT32_MAX?" (max reached, could be biased if simulation cut)":""));
  printf("Probability of non outgoing outbreak and its statistical uncertainty: %22.15e +/- %22.15e%s\n",tdata[0].pe,sqrt(tdata[0].pe*(1.-tdata[0].pe)/(cp.npaths-1.)),(tdata[0].maxedoutmintimeindex<INT32_MAX?" (max reached, could be biased if simulation cut)":""));
  printf("Probability of reaching maximum as defined by nimax/npostestmax and its statistical uncertainty: %22.15e +/- %22.15e\n",tdata[0].pm,sqrt(tdata[0].pm*(1.-tdata[0].pm)/(cp.npaths-1.)));
  printf("Probability of reaching maximum as defined by nimax/npostestmax and its statistical uncertainty: %22.15e +/- %22.15e\n",tdata[0].pm,sqrt(tdata[0].pm*(1.-tdata[0].pm)/(cp.npaths-1.)));
  printf("Extinction time, if it occurs is %22.15e +/- %22.15e%s\n",tdata[0].tenz_mean,tdata[0].tenz_std,(tdata[0].maxedoutmintimeindex<INT32_MAX?" (max reached, could be biased if simulation cut)":""));

  int shift=tdata[tmaxnpers].tlppnnpers;
  printf("\nCurrent infection (non-isolated infected individuals) timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].inf_timeline_mean_ext[j],tdata[tmaxnpers].inf_timeline_std_ext[j],tdata[tmaxnpers].inf_timeline_mean_noext[j],tdata[tmaxnpers].inf_timeline_std_noext[j],inf_timeline_mean[j],inf_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));

  printf("\nNew infections (new infected individuals) timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].newinf_timeline_mean_ext[j],tdata[tmaxnpers].newinf_timeline_std_ext[j],tdata[tmaxnpers].newinf_timeline_mean_noext[j],tdata[tmaxnpers].newinf_timeline_std_noext[j],newinf_timeline_mean[j],newinf_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));

  if(!isnan(cp.pars.tdeltat)) {
    printf("\nNew positive test timeline, for paths with extinction vs no extinction vs overall is:\n");
    for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].newpostest_timeline_mean_ext[j],tdata[tmaxnpers].newpostest_timeline_std_ext[j],tdata[tmaxnpers].newpostest_timeline_mean_noext[j],tdata[tmaxnpers].newpostest_timeline_std_noext[j],newpostest_timeline_mean[j],newpostest_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));
  }

#ifdef SEC_INF_TIMELINES
  printf("\nCurrent infection (non-isolated infected individuals) timeline for the second infection category, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].secinf_timeline_mean_ext[j],tdata[tmaxnpers].secinf_timeline_std_ext[j],tdata[tmaxnpers].secinf_timeline_mean_noext[j],tdata[tmaxnpers].secinf_timeline_std_noext[j],secinf_timeline_mean[j],secinf_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));

  printf("\nNew infections (new infected individuals) timeline for the second infection category, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].newsecinf_timeline_mean_ext[j],tdata[tmaxnpers].newsecinf_timeline_std_ext[j],tdata[tmaxnpers].newsecinf_timeline_mean_noext[j],tdata[tmaxnpers].newsecinf_timeline_std_noext[j],newsecinf_timeline_mean[j],newsecinf_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));

  if(!isnan(cp.pars.tdeltat)) {
    printf("\nNew positive test timeline for the second infection category, for paths with extinction vs no extinction vs overall is:\n");
    for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].newsecpostest_timeline_mean_ext[j],tdata[tmaxnpers].newsecpostest_timeline_std_ext[j],tdata[tmaxnpers].newsecpostest_timeline_mean_noext[j],tdata[tmaxnpers].newsecpostest_timeline_std_noext[j],newsecpostest_timeline_mean[j],newsecpostest_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));
  }
#endif

  printf("\nReff timeline, for paths with extinction vs no extinction vs overall is:\n");
  for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].reff_timeline_mean_ext[j],tdata[tmaxnpers].reff_timeline_std_ext[j],tdata[tmaxnpers].reff_timeline_mean_noext[j],tdata[tmaxnpers].reff_timeline_std_noext[j],reff_timeline_mean[j],reff_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));

  #ifdef OBSREFF_OUTPUT
  if(!isnan(cp.pars.tdeltat)) {
    printf("\nObservable Reff timeline, for paths with extinction vs no extinction vs overall is:\n");
    for(j=0; j<tdata[tmaxnpers].tlpptnvpers; ++j) printf("%6.2f: %22.15e +/- %22.15e\t%22.15e +/- %22.15e\t%22.15e +/- %22.15e%s\n",(j-shift)/(double)cp.nbinsperunit,tdata[tmaxnpers].reffobs_timeline_mean_ext[j],tdata[tmaxnpers].reffobs_timeline_std_ext[j],tdata[tmaxnpers].reffobs_timeline_mean_noext[j],tdata[tmaxnpers].reffobs_timeline_std_noext[j],reffobs_timeline_mean[j],reffobs_timeline_std[j],(j-shift<tdata[0].maxedoutmintimeindex?"":" (max reached, biased if simulation cut)"));
  }
  #endif

  if(cp.ninfhist) {
    uint32_t maxninfnbins=tdata[cp.nthreads-1].ninfbins;
    uint32_t b;
    int tmaxninfnbins=cp.nthreads-1;
    int tp;

    for(t=cp.nthreads-2; t>=0; --t) {

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

  free(tdata[tmaxnpers].inf_timeline_mean_ext);
  free(tdata[tmaxnpers].inf_timeline_std_ext);
  free(tdata[tmaxnpers].inf_timeline_mean_noext);
  free(tdata[tmaxnpers].inf_timeline_std_noext);
  free(tdata[tmaxnpers].newinf_timeline_mean_ext);
  free(tdata[tmaxnpers].newinf_timeline_std_ext);
  free(tdata[tmaxnpers].newinf_timeline_mean_noext);
  free(tdata[tmaxnpers].newinf_timeline_std_noext);
  free(tdata[tmaxnpers].newpostest_timeline_mean_ext);
  free(tdata[tmaxnpers].newpostest_timeline_std_ext);
  free(tdata[tmaxnpers].newpostest_timeline_mean_noext);
  free(tdata[tmaxnpers].newpostest_timeline_std_noext);
  #ifdef SEC_INF_TIMELINES
  free(tdata[tmaxnpers].secinf_timeline_mean_ext);
  free(tdata[tmaxnpers].secinf_timeline_std_ext);
  free(tdata[tmaxnpers].secinf_timeline_mean_noext);
  free(tdata[tmaxnpers].secinf_timeline_std_noext);
  free(tdata[tmaxnpers].newsecinf_timeline_mean_ext);
  free(tdata[tmaxnpers].newsecinf_timeline_std_ext);
  free(tdata[tmaxnpers].newsecinf_timeline_mean_noext);
  free(tdata[tmaxnpers].newsecinf_timeline_std_noext);
  free(tdata[tmaxnpers].newsecpostest_timeline_mean_ext);
  free(tdata[tmaxnpers].newsecpostest_timeline_std_ext);
  free(tdata[tmaxnpers].newsecpostest_timeline_mean_noext);
  free(tdata[tmaxnpers].newsecpostest_timeline_std_noext);
  #endif
  free(tdata[tmaxnpers].reff_timeline_mean_ext);
  free(tdata[tmaxnpers].reff_timeline_std_ext);
  free(tdata[tmaxnpers].reff_timeline_n_ext);
  free(tdata[tmaxnpers].reff_timeline_mean_noext);
  free(tdata[tmaxnpers].reff_timeline_std_noext);
  free(tdata[tmaxnpers].reff_timeline_n_noext);
  #ifdef OBSREFF_OUTPUT
  free(tdata[tmaxnpers].reffobs_timeline_mean_ext);
  free(tdata[tmaxnpers].reffobs_timeline_std_ext);
  free(tdata[tmaxnpers].reffobs_timeline_n_ext);
  free(tdata[tmaxnpers].reffobs_timeline_mean_noext);
  free(tdata[tmaxnpers].reffobs_timeline_std_noext);
  free(tdata[tmaxnpers].reffobs_timeline_n_noext);
  #endif
  free(tdata);

  fflush(stdout);
  fflush(stderr);
  close(cp.oout);
  close(cp.eout);

  return 0;
}

void* simthread(void* arg)
{
  thread_data* data=(thread_data*)arg;
  config_pars const* cp=data->cp;
  data->nnzpaths=0;
  data->tlppnnpers=0;
  data->tlpptnvpers=data->npers;
  data->commper_mean=0;
#ifdef NUMEVENTSSTATS
  data->nevents_mean=0;
#endif
  data->pe=0;
  data->penz=0;
  data->pm=0;
  data->tenz_mean=0;
  data->tenz_std=0;
  data->maxedoutmintimeindex=INT32_MAX;

  data->inf_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->inf_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->inf_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->inf_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newinf_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newinf_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newinf_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newinf_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newpostest_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newpostest_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newpostest_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newpostest_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  #ifdef SEC_INF_TIMELINES
  data->secinf_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->secinf_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->secinf_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->secinf_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecinf_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecinf_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecinf_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecinf_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecpostest_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecpostest_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecpostest_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->newsecpostest_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  #endif
  data->reff_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reff_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reff_timeline_n_ext=(uint64_t*)malloc(data->tlpptnvpers*sizeof(uint64_t));
  data->reff_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reff_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reff_timeline_n_noext=(uint64_t*)malloc(data->tlpptnvpers*sizeof(uint64_t));
  #ifdef OBSREFF_OUTPUT
  data->reffobs_timeline_mean_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reffobs_timeline_std_ext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reffobs_timeline_n_ext=(uint64_t*)malloc(data->tlpptnvpers*sizeof(uint64_t));
  data->reffobs_timeline_mean_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reffobs_timeline_std_noext=(double*)malloc(data->tlpptnvpers*sizeof(double));
  data->reffobs_timeline_n_noext=(uint64_t*)malloc(data->tlpptnvpers*sizeof(uint64_t));
  #endif

  memset(data->inf_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->inf_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->inf_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->inf_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newinf_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newinf_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newinf_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newinf_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newpostest_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newpostest_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newpostest_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newpostest_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  #ifdef SEC_INF_TIMELINES
  memset(data->secinf_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->secinf_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->secinf_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->secinf_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecinf_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecinf_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecinf_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecinf_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecpostest_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecpostest_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecpostest_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->newsecpostest_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  #endif
  memset(data->reff_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reff_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reff_timeline_n_ext, 0, data->tlpptnvpers*sizeof(uint64_t));
  memset(data->reff_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reff_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reff_timeline_n_noext, 0, data->tlpptnvpers*sizeof(uint64_t));
  #ifdef OBSREFF_OUTPUT
  memset(data->reffobs_timeline_mean_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reffobs_timeline_std_ext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reffobs_timeline_n_ext, 0, data->tlpptnvpers*sizeof(uint64_t));
  memset(data->reffobs_timeline_mean_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reffobs_timeline_std_noext, 0, data->tlpptnvpers*sizeof(double));
  memset(data->reffobs_timeline_n_noext, 0, data->tlpptnvpers*sizeof(uint64_t));
  #endif

  data->ninfbins=0;
  data->ngeninfs=NULL;

  char* tloutbuf=NULL;
  const ssize_t tlobasize=cp->tloutbufsize*INT64_C(1024*1024);
  ssize_t tlobsize=0;
  ssize_t (*buf_write_func)(std_summary_stats const* stats, char* buf)=NULL;

  if(cp->tlout) {
    tloutbuf=(char*)malloc(tlobasize);

    if(!tloutbuf) {
      fprintf(stderr,"%s: Error: Cannot allocate memory buffer of size %" PRIu32 " MB for timeline file output buffer\n",__func__,cp->tloutbufsize);
      exit(1);
    }

    if(cp->pars.timetype!=ro_time_pri_created && cp->pars.timetype!=ro_time_pri_flat_comm) {

      if(isnan(cp->pars.tdeltat)) buf_write_func=tlo_write_reltime_path;
      else buf_write_func=tlo_write_reltime_postest_path;

    } else {

      if(isnan(cp->pars.tdeltat)) buf_write_func=tlo_write_reg_path;
      else buf_write_func=tlo_write_reg_postest_path;
    }
  }

#ifdef CT_OUTPUT
  const ssize_t ctobasize=cp->ctoutbufsize*INT64_C(1024*1024);
  ssize_t ctobsize=0;
  char* ctoutbuf=NULL;

  if(cp->ctout) {
    ctoutbuf=(char*)malloc(ctobasize);

    if(!ctoutbuf) {
      fprintf(stderr,"%s: Error: Cannot allocate memory buffer of size %" PRIu32 " MB for contact tracing output buffer\n",__func__,cp->ctoutbufsize);
      exit(1);
    }
  }
#endif

  sim_vars sv;
  //branchsim_init(&sv,.lambda=0.5,.p=0.8,.tmax=6);

  sim_init(&sv,&cp->pars,data->r);

  //printf("Min is %lu, max is %lu\n",gsl_rng_min(r),gsl_rng_max(r));

  //uint32_t nr=0;

  std_summary_stats stats;

  sim_set_proc_data(&sv, &stats);
  sim_set_path_init_proc_func(&sv, std_stats_path_init);
  sim_set_path_end_proc_func(&sv, std_stats_path_end);

  if(cp->pars.timetype==ro_time_pri_created || cp->pars.timetype==ro_time_pri_flat_comm || cp->pars.timetype==ro_time_first_pos_test_results) sim_set_pri_init_proc_func(&sv, std_stats_pri_init);
  else sim_set_pri_init_proc_func(&sv, std_stats_pri_init_rel);

  if( cp->pars.timetype==ro_time_first_pos_test_results) {
    sim_set_new_inf_proc_func(&sv, std_stats_new_inf_first_pos_test_results);

    if(cp->ninfhist) {
      sim_set_end_inf_proc_func(&sv, std_stats_end_inf_rec_ninfs);
      sim_set_new_inf_proc_noevent_func(&sv, std_stats_noevent_new_inf_rec_ninfs_first_pos_test_results);

    } else {
      sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
      sim_set_new_inf_proc_noevent_func(&sv, std_stats_noevent_new_inf_first_pos_test_results);
    }

  } else {
    sim_set_new_inf_proc_func(&sv, std_stats_new_inf);

    if(cp->ninfhist) {
      sim_set_end_inf_proc_func(&sv, std_stats_end_inf_rec_ninfs);
      sim_set_new_inf_proc_noevent_func(&sv, std_stats_noevent_new_inf_rec_ninfs);

    } else {
      sim_set_end_inf_proc_func(&sv, std_stats_end_inf);
      sim_set_new_inf_proc_noevent_func(&sv, std_stats_noevent_new_inf);
    }
  }

  sim_set_ii_alloc_proc_func(&sv, std_stats_ii_alloc);

  if(cp->nimax == UINT32_MAX) {

    if(cp->npostestmax == UINT32_MAX) sim_set_new_event_proc_func(&sv, std_stats_new_event);
    else sim_set_new_event_proc_func(&sv, std_stats_new_event_npostestmax);

  } else sim_set_new_event_proc_func(&sv, std_stats_new_event_nimax);

  int (*simfunc)(sim_vars*);

  if(cp->pars.popsize>0) {
    finitepopsim_init(&sv);
    simfunc=finitepopsim;

  } else {
    branchsim_init(&sv);
    simfunc=branchsim;
  }

  if(cp->ninfhist) std_stats_init(&sv, cp->nbinsperunit, true);

  else std_stats_init(&sv, cp->nbinsperunit, false);

  stats.lmax=cp->lmax;
  stats.nimax=cp->nimax;
  stats.npostestmax=cp->npostestmax;
  stats.npostestmaxnunits=cp->npostestmaxnunits;
  int i,j,k;
  uint32_t curset=data->id;
  uint32_t initpath;
  uint32_t npaths;
  uint32_t* abs_inf_timeline;
  uint32_t* abs_newinf_timeline;
  uint32_t* abs_newpostest_timeline;
  #ifdef SEC_INF_TIMELINES
  uint32_t* abs_secinf_timeline;
  uint32_t* abs_newsecinf_timeline;
  uint32_t* abs_newsecpostest_timeline;
  #endif
  ext_timeline_info* abs_ext_timeline;
  int32_t dshift;
  ssize_t maxwrite;
#ifdef SEC_INF_TIMELINES
  const ssize_t binsize=2*(2+(!isnan(cp->pars.tdeltat)))*sizeof(uint32_t);
#else
  const ssize_t binsize=(2+(!isnan(cp->pars.tdeltat)))*sizeof(uint32_t);
#endif
  ext_timeline_info* eti;

  do {

    //printf("%22.15e\t%22.15e\n",curset*data->npathsperset,(curset+1)*data->npathsperset);
    initpath=round(curset*data->npathsperset);
    npaths=round((curset+1)*data->npathsperset)-initpath;
    //printf("npaths %u\n",npaths);

    for(i=npaths-1; i>=0; --i) {
      simfunc(&sv);

      eti=stats.ext_timeline-stats.tlshift;
      data->commper_mean+=eti->commpersum;
#ifdef NUMEVENTSSTATS
      data->nevents_mean+=stats.ext_timeline->neventssum;
#endif
      //nr+=stats.n_ended_infections;
      abs_inf_timeline=stats.pp_inf_timeline-stats.tlppnnpers;
      abs_newinf_timeline=stats.pp_newinf_timeline-stats.tlppnnpers;
      abs_newpostest_timeline=stats.pp_newpostest_timeline-stats.tlppnnpers;
      #ifdef SEC_INF_TIMELINES
      abs_secinf_timeline=stats.pp_secinf_timeline-stats.tlppnnpers;
      abs_newsecinf_timeline=stats.pp_newsecinf_timeline-stats.tlppnnpers;
      abs_newsecpostest_timeline=stats.pp_newsecpostest_timeline-stats.tlppnnpers;
      #endif
      abs_ext_timeline=stats.pp_ext_timeline-stats.tlppnnpers;
      dshift=stats.tlshifta-stats.tlshift;

      if(cp->tlout) {
	maxwrite=16+binsize*stats.tlpptnvpers;

	if(tlobsize+maxwrite > tlobasize) {

	  if(maxwrite > tlobasize) {
	    fprintf(stderr,"%s: Error: Timeline output from a single path cannot exceed the allocated per-thread memory buffer size!\n",__func__);
	    exit(1);
	  }
	  pthread_mutex_lock(data->tlflock);
	  //printf("Writing %" PRIi64 " bytes\n",(int64_t)tlobsize);

	  if(write(cp->tlout, tloutbuf, tlobsize)!=tlobsize) {
	    perror(__func__);
	    exit(1);
	  }
	  pthread_mutex_unlock(data->tlflock);
	  tlobsize=0;
	}
	tlobsize+=buf_write_func(&stats, tloutbuf+tlobsize);
      }

#ifdef CT_OUTPUT
      if(cp->ctout) {
	qsort(stats.ctentries, stats.nctentries, sizeof(ctposinf*), ctcompar);
	//printf("Entries: %u, tnctevents: %li\n",stats.nctentries,stats.tnctevents);
	maxwrite=stats.nctentries*20;
	//printf("Maxwrite: %li\n",maxwrite);

	if(ctobsize+maxwrite > ctobasize) {

	  if(maxwrite > ctobasize) {
	    fprintf(stderr,"%s: Error: Contact tracing output from a single path cannot exceed the allocated per-thread memory buffer size!\n",__func__);
	    exit(1);
	  }
	  pthread_mutex_lock(data->ctflock);
	  //printf("Writing %" PRIi64 " bytes\n",(int64_t)ctobsize);

	  if(write(cp->ctout, ctoutbuf, ctobsize)!=ctobsize) {
	    perror(__func__);
	    exit(1);
	  }
	  pthread_mutex_unlock(data->ctflock);
	  ctobsize=0;
	}
	ctobsize+=ct_write_func(&stats, ctoutbuf+ctobsize);
      }
#endif

      int32_t ndiff=stats.tlppnnpers-data->tlppnnpers;
      int32_t pdiff=(int32_t)stats.tlpptnvpers-data->tlpptnvpers-ndiff;
      dshift=(ndiff<0?-ndiff:0);

      if(ndiff<0) ndiff=0;
      if(pdiff<0) pdiff=0;

      realloc_thread_timelines(data, ndiff, pdiff);

      if(stats.ninfbins > data->ninfbins) {
	data->ngeninfs=(uint64_t*)realloc(data->ngeninfs,stats.ninfbins*sizeof(uint64_t));
	memset(data->ngeninfs+data->ninfbins,0,(stats.ninfbins-data->ninfbins)*sizeof(uint64_t));
	data->ninfbins=stats.ninfbins;
      }

      data->pm+=(stats.maxedoutmintimeindex<INT32_MAX);

      if(stats.extinction) {
	data->pe+=stats.extinction;

	if(isinf(stats.extinction_time)!=-1) {
	  data->penz++;
	  data->nnzpaths++;
	  data->tenz_mean+=stats.extinction_time;
	  data->tenz_std+=stats.extinction_time*stats.extinction_time;
	}

	for(j=stats.tlpptnvpers-1; j>=0; --j) {
	  k=dshift+j;
	  //printf("data->inf_timeline_mean_ext[%i]+=abs_inf_timeline[%i] (%u)\n",k,j,abs_inf_timeline[j]);
	  data->inf_timeline_mean_ext[k]+=abs_inf_timeline[j];
	  data->inf_timeline_std_ext[k]+=(double)abs_inf_timeline[j]*abs_inf_timeline[j];
	  data->newinf_timeline_mean_ext[k]+=abs_newinf_timeline[j];
	  data->newinf_timeline_std_ext[k]+=(double)abs_newinf_timeline[j]*abs_newinf_timeline[j];
	  data->newpostest_timeline_mean_ext[k]+=abs_newpostest_timeline[j];
	  data->newpostest_timeline_std_ext[k]+=(double)abs_newpostest_timeline[j]*abs_newpostest_timeline[j];
          #ifdef SEC_INF_TIMELINES
       	  data->secinf_timeline_mean_ext[k]+=abs_secinf_timeline[j];
	  data->secinf_timeline_std_ext[k]+=(double)abs_secinf_timeline[j]*abs_secinf_timeline[j];
	  data->newsecinf_timeline_mean_ext[k]+=abs_newsecinf_timeline[j];
	  data->newsecinf_timeline_std_ext[k]+=(double)abs_newsecinf_timeline[j]*abs_newsecinf_timeline[j];
	  data->newsecpostest_timeline_mean_ext[k]+=abs_newsecpostest_timeline[j];
	  data->newsecpostest_timeline_std_ext[k]+=(double)abs_newsecpostest_timeline[j]*abs_newsecpostest_timeline[j];
          #endif

	  if(abs_ext_timeline[j].n) {
	    data->reff_timeline_mean_ext[k]+=abs_ext_timeline[j].rsum;
	    data->reff_timeline_std_ext[k]+=abs_ext_timeline[j].r2sum;
	    data->reff_timeline_n_ext[k]+=abs_ext_timeline[j].n;
	  }

	  #ifdef OBSREFF_OUTPUT
	  if(abs_ext_timeline[j].nobs) {
	    data->reffobs_timeline_mean_ext[k]+=abs_ext_timeline[j].robssum;
	    data->reffobs_timeline_std_ext[k]+=abs_ext_timeline[j].robs2sum;
	    data->reffobs_timeline_n_ext[k]+=abs_ext_timeline[j].nobs;
	  }
          #endif
	}

      } else {
	data->nnzpaths++;

	if(stats.maxedoutmintimeindex < data->maxedoutmintimeindex) data->maxedoutmintimeindex=stats.maxedoutmintimeindex;

	for(j=stats.tlpptnvpers-1; j>=0; --j) {
	  k=dshift+j;
	  //printf("data->inf_timeline_mean_ext[%i]+=abs_inf_timeline[%i] (%u)\n",k,j,abs_inf_timeline[j]);
	  data->inf_timeline_mean_noext[k]+=abs_inf_timeline[j];
	  data->inf_timeline_std_noext[k]+=(double)abs_inf_timeline[j]*abs_inf_timeline[j];
	  data->newinf_timeline_mean_noext[k]+=abs_newinf_timeline[j];
	  data->newinf_timeline_std_noext[k]+=(double)abs_newinf_timeline[j]*abs_newinf_timeline[j];
	  data->newpostest_timeline_mean_noext[k]+=abs_newpostest_timeline[j];
	  data->newpostest_timeline_std_noext[k]+=(double)abs_newpostest_timeline[j]*abs_newpostest_timeline[j];
          #ifdef SEC_INF_TIMELINES
	  data->secinf_timeline_mean_noext[k]+=abs_secinf_timeline[j];
	  data->secinf_timeline_std_noext[k]+=(double)abs_secinf_timeline[j]*abs_secinf_timeline[j];
	  data->newsecinf_timeline_mean_noext[k]+=abs_newsecinf_timeline[j];
	  data->newsecinf_timeline_std_noext[k]+=(double)abs_newsecinf_timeline[j]*abs_newsecinf_timeline[j];
	  data->newsecpostest_timeline_mean_noext[k]+=abs_newsecpostest_timeline[j];
	  data->newsecpostest_timeline_std_noext[k]+=(double)abs_newsecpostest_timeline[j]*abs_newsecpostest_timeline[j];
          #endif

	  if(abs_ext_timeline[j].n) {
	    data->reff_timeline_mean_noext[k]+=abs_ext_timeline[j].rsum;
	    data->reff_timeline_std_noext[k]+=abs_ext_timeline[j].r2sum;
	    data->reff_timeline_n_noext[k]+=abs_ext_timeline[j].n;
	  }

	  #ifdef OBSREFF_OUTPUT
	  if(abs_ext_timeline[j].nobs) {
	    data->reffobs_timeline_mean_noext[k]+=abs_ext_timeline[j].robssum;
	    data->reffobs_timeline_std_noext[k]+=abs_ext_timeline[j].robs2sum;
	    data->reffobs_timeline_n_noext[k]+=abs_ext_timeline[j].nobs;
	  }
          #endif
	}
      }

      for(j=stats.ninfbins-1; j>=0; --j) {
	data->ngeninfs[j]+=eti->ngeninfs[j];
      }
    }
    curset=__sync_fetch_and_add(data->set,1);

  } while(curset<data->nsets);

  if(cp->tlout) {
    pthread_mutex_lock(data->tlflock);
    //printf("Writing %" PRIi64 " bytes\n",(int64_t)tlobsize);

    if(write(cp->tlout, tloutbuf, tlobsize)!=tlobsize) {
      perror(__func__);
      exit(1);
    }
    pthread_mutex_unlock(data->tlflock);
    free(tloutbuf);
  }

#ifdef CT_OUTPUT
  if(cp->ctout) {
    pthread_mutex_lock(data->ctflock);
    //printf("Writing %" PRIi64 " bytes\n",(int64_t)ctobsize);

    if(write(cp->ctout, ctoutbuf, ctobsize)!=ctobsize) {
      perror(__func__);
      exit(1);
    }
    pthread_mutex_unlock(data->ctflock);
    free(ctoutbuf);
  }
#endif

  std_stats_free(&stats);
  if(cp->pars.popsize>0) finitepopsim_free(&sv);
  else branchsim_free(&sv);

  return NULL;
}
