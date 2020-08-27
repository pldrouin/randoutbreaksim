/**
 * @file simulation.c
 * @brief Simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "simulation.h"

int sim_pars_check(model_pars const* pars)
{
  int ret=0;

  if(pars->tbar<=0) {
    fprintf(stderr,"%s: Error: tbar must be greater than 0\n",__func__);
    ret-=1;
  }

  if(pars->p<0) {
    fprintf(stderr,"%s: Error: p must be non-negative\n",__func__);
    ret-=2;
  }

  if(pars->lambda<=0) {
    fprintf(stderr,"%s: Error: lambda must be greater than 0\n",__func__);
    ret-=4;
  }

  if(pars->kappa<=0) {
    fprintf(stderr,"%s: Error: kappa must be greater than 0\n",__func__);
    ret-=8;
  }

  if(pars->lbar<0) {
    fprintf(stderr,"%s: Error: lbar must be non-negative\n",__func__);
    ret-=16;

  } else if(pars->lbar>0 && pars->kappal<=0) {
      fprintf(stderr,"%s: Error: kappal must be greater than 0 if lbar>0\n",__func__);
      ret-=32;
  }

  if(pars->pit<0) {
    fprintf(stderr,"%s: Error: pit must be non-negative\n",__func__);
    ret-=64;

  } else if(pars->pit>0) {

    if(pars->itbar<0) {
      fprintf(stderr,"%s: Error: itbar must be non-negative if pit>0\n",__func__);
      ret-=128;
    }

    if(pars->kappait<=0) {
      fprintf(stderr,"%s: Error: kappait must be greater than 0 if pit>0\n",__func__);
      ret-=256;
    }
  }

  if(pars->q<0) {
    fprintf(stderr,"%s: Error: q must be non-negative\n",__func__);
    ret-=512;

  } else if(pars->q>0) {

    if(pars->mbar<0) {
      fprintf(stderr,"%s: Error: mbar must be non-negative if q>0\n",__func__);
      ret-=1024;
    }

    if(pars->kappaq<=0) {
      fprintf(stderr,"%s: Error: kappaq must be greater than 0 if q>0\n",__func__);
      ret-=2048;
    }

    if(pars->pim<0) {
      fprintf(stderr,"%s: Error: pim must be non-negative\n",__func__);
      ret-=4096;

    } else if(pars->pim>0) {

      if(pars->imbar<0) {
	fprintf(stderr,"%s: Error: imbar must be non-negative if pim>0\n",__func__);
	ret-=8192;
      }

      if(pars->kappaim<=0) {
	fprintf(stderr,"%s: Error: kappaim must be greater than 0 if pim>0\n",__func__);
	ret-=16384;
      }
    }
  }

  if(pars->tmax<=0) {
    fprintf(stderr,"%s: Error: tmax must be greater than 0\n",__func__);
    ret-=32768;
  }

  if(pars->nstart<=0) {
    fprintf(stderr,"%s: Error: nstart must be greater than 0\n",__func__);
    ret-=65536;
  }
  return ret;
}

void sim_pars_init(model_pars* pars)
{
  pars->tbar=NAN;
  pars->p=NAN;
  pars->lambda=NAN;
  pars->kappa=NAN;
  pars->lbar=NAN;
  pars->kappal=NAN;
  pars->q=0;
  pars->mbar=NAN;
  pars->kappaq=NAN;
  pars->pit=0;
  pars->itbar=NAN;
  pars->kappait=NAN;
  pars->pim=NAN;
  pars->imbar=NAN;
  pars->kappaim=NAN;
  pars->R0=NAN;
  pars->mu=NAN;
  pars->t95=NAN;
  pars->m95=NAN;
  pars->l95=NAN;
  pars->it95=NAN;
  pars->im95=NAN;
  pars->tmax=INFINITY;
  pars->nstart=1;
}

void sim_init(sim_vars* sv, model_pars const* pars, const gsl_rng* r)
{
  sv->pars=*pars;
  sv->r=r;

  PER_COND;

  sv->dataptr=NULL;
  sv->ii_alloc_proc_func=default_ii_alloc_proc_func;
  sv->new_event_proc_func=default_event_proc_func;
  sv->new_inf_proc_func=dummy_proc_func_one_par;
  sv->end_inf_proc_func=dummy_proc_func_two_pars;
  sv->inf_proc_func_noevent=dummy_proc_func_two_pars;
}

void sim_free(sim_vars* sv)
{
  for(uint32_t i=0; i<sv->brsim.nlayers; ++i) if(sv->brsim.iis[i].dataptr) free(sv->brsim.iis[i].dataptr);
  free(sv->brsim.iis);
}
