/**
 * @file simulation.c
 * @brief Common simulation functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "simulation.h"

#define DO_EXPAND(VAL)  VAL ## 1
#define EXPAND(VAL)     DO_EXPAND(VAL)

#if defined(DEBUG_PRINTF) || (EXPAND(DEBUG_PRINTF) != 1)
int __ro_debug=1;
#endif

void sim_pars_init(model_pars* pars)
{
  pars->pinfpri=1;
  pars->tbar=NAN;
  pars->p=NAN;
  pars->mu=NAN;
  pars->sigma=NAN;
  pars->rsigma=NAN;
  pars->g_ave=NAN;
  pars->g_ave_transm=NAN;
  pars->lambda=NAN;
  pars->lambda_uncut=NAN;
  pars->lambdap=NAN;
  pars->pinf=NAN;
#ifdef DUAL_PINF
  pars->ppip=0;
  pars->rpinfp=NAN;
  pars->rpshedp=NAN;
  pars->qp=NAN;
#endif
  pars->R0=NAN;
  pars->kappa=NAN;
  pars->lbar=NAN;
  pars->kappal=NAN;
  pars->q=NAN;
  pars->mbar=NAN;
  pars->kappaq=NAN;
#ifdef CT_OUTPUT
  pars->ctwindow=0;
  pars->pt=NAN;
#endif
  pars->pit=NAN;
  pars->itbar=NAN;
  pars->kappait=NAN;
  pars->pim=NAN;
  pars->imbar=NAN;
  pars->kappaim=NAN;
  pars->t95=NAN;
  pars->m95=NAN;
  pars->l95=NAN;
  pars->it95=NAN;
  pars->im95=NAN;
  pars->ttpr=NAN;
  pars->mtpr=NAN;
  pars->tdeltat=NAN;
  pars->tmax=INT32_MAX;
  pars->nstart=1;
  pars->popsize=0;
  pars->pricommpertype=ro_pricommper_main|ro_pricommper_alt|ro_pricommper_alt_use_tpr;
  pars->grouptype=ro_group_log_plus_1;
  pars->timetype=ro_time_pri_created;
  pars->pathtype=ro_all_paths;
  pars->groupinteractions=false;
}

void sim_init(sim_vars* sv, model_pars const* pars, const gsl_rng* r)
{
  sv->pars=*pars;
  sv->r=r;

  if(pars->pinfpri==1) sv->gen_n_pri_inf=gen_n_pri_inf_nstart;
  else sv->gen_n_pri_inf=gen_n_pri_inf_binom_pinfpri_nstart;

  switch(pars->timetype) {
    case ro_time_pri_created:
    case ro_time_first_pos_test_results:
      sv->gen_time_origin_func=gen_time_origin_pri_created;
      break;

    case ro_time_pri_infectious:
      sv->gen_time_origin_func=gen_time_origin_pri_infectious;
      break;

    case ro_time_pri_end_comm:
      sv->gen_time_origin_func=gen_time_origin_pri_end_comm;
      break;

    case ro_time_pri_test_results:
      sv->gen_time_origin_func=gen_time_origin_pri_test_results;
      break;

    case ro_time_pri_flat_comm:
      sv->gen_time_origin_func=gen_time_origin_pri_flat_comm;
      break;
  }

  PER_COND;

  sv->dataptr=NULL;
  sv->path_init_proc_func=dummy_proc_func_sv;
  sv->path_end_proc_func=dummy_proc_bool_func_sv;
  sv->pri_init_proc_func=dummy_proc_func_sv_ii2;
  sv->ii_alloc_proc_func=default_ii_alloc_proc_func;
  sv->new_event_proc_func=default_event_proc_func;
  sv->new_inf_proc_func=dummy_proc_func_sv_ii2;
  sv->new_inf_proc_func_noevent=dummy_proc_func_sv_ii2;
  sv->end_inf_proc_func=dummy_proc_func_sv_ii2;
}
