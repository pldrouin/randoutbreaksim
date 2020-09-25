/**
 * @file model_parameters.c
 * @brief Model parameter functions.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "model_parameters.h"

int model_solve_pars(model_pars* pars)
{
  if((isnan(pars->tbar)==0) + (isnan(pars->lambda)==0) + (isnan(pars->g_ave)==0 || isnan(pars->p)==0 || isnan(pars->mu)==0) + (isnan(pars->pinf)==0) + (isnan(pars->R0)==0) != 4) {
    fprintf(stderr,"%s: Error: An invalid combination of tbar, lambda, g_ave, p, mu, pinf and R0 parameters was provided.\n",__func__);
    return -1;
  }

  if(model_solve_R0_group(pars)) return -2;
  
  printf("\nBasic reproduction parameters are:\n");
  printf("lambda:\t%22.15e\n",pars->lambda);
  printf("tbar:\t%22.15e\n",pars->tbar);
  printf("g_ave:\t%22.15e\n",pars->g_ave);
  printf("mu:\t%22.15e\n",pars->mu);
  printf("p:\t%22.15e\n",pars->p);
  printf("pinf:\t%22.15e\n",pars->pinf);
  printf("R0:\t%22.15e\n",pars->R0);

  if((isnan(pars->kappa)==0) + (isnan(pars->t95)==0) != 1) {
    fprintf(stderr,"%s: Error: Either the kappa parameter or the t95 parameter must be provided.\n",__func__);
    return -3;
  }

  if(model_solve_gamma_group(&pars->tbar, &pars->kappa, &pars->t95)) {
    fprintf(stderr,"%s: Error: Cannot solve parameters for the main time gamma distribution\n",__func__);
    return -4;

  }
  pars->ta=pars->tbar*pars->kappa;
  pars->tb=1/pars->kappa;
  printf("\nParameters for the main time gamma distribution:\n");
  printf("tbar:\t%22.15e\n",pars->tbar);
  printf("kappa:\t%22.15e\n",pars->kappa);
  printf("t95:\t%22.15e\n",pars->t95);
  printf("ta:\t%22.15e\n",pars->ta);
  printf("tb:\t%22.15e\n",pars->tb);

  if(pars->pit>0) {

    if((isnan(pars->kappait)==0) + (isnan(pars->it95)==0) != 1) {
      fprintf(stderr,"%s: Error: Either the kappait parameter or the it95 parameter must be provided.\n",__func__);
      return -5;
    }

    if(model_solve_gamma_group(&pars->itbar, &pars->kappait, &pars->it95)) {
      fprintf(stderr,"%s: Error: Cannot solve parameters for the interrupted main time gamma distribution\n",__func__);
      return -6;

    } else {
      pars->ita=pars->itbar*pars->kappait;
      pars->itb=1/pars->kappait;
      printf("\nParameters for the interrupted main time gamma distribution:\n");
      printf("pit:\t%22.15e\n",pars->pit);
      printf("itbar:\t%22.15e\n",pars->itbar);
      printf("kappait:%22.15e\n",pars->kappait);
      printf("it95:\t%22.15e\n",pars->it95);
      printf("ita:\t%22.15e\n",pars->ita);
      printf("itb:\t%22.15e\n",pars->itb);
    }
  }

  if(pars->q>0) {

    if((isnan(pars->kappaq)==0) + (isnan(pars->m95)==0) != 1) {
      fprintf(stderr,"%s: Error: Either the kappaq parameter or the m95 parameter must be provided.\n",__func__);
      return -7;
    }

    if(model_solve_gamma_group(&pars->mbar, &pars->kappaq, &pars->m95)) {
      fprintf(stderr,"%s: Error: Cannot solve parameters for the alternate time gamma distribution\n",__func__);
      return -8;

    } else {
      pars->ma=pars->mbar*pars->kappaq;
      pars->mb=1/pars->kappaq;
      printf("\nParameters for the alternate time gamma distribution:\n");
      printf("q:\t%22.15e\n",pars->q);
      printf("mbar:\t%22.15e\n",pars->mbar);
      printf("kappaq:\t%22.15e\n",pars->kappaq);
      printf("m95:\t%22.15e\n",pars->m95);
      printf("ma:\t%22.15e\n",pars->ma);
      printf("mb:\t%22.15e\n",pars->mb);
    }

    if(isnan(pars->pim)) pars->pim=pars->pit;

    if(pars->pim>0) {

      if(isnan(pars->imbar) && isnan(pars->kappaim) && isnan(pars->im95)) {
	pars->imbar=pars->itbar;
	pars->kappaim=pars->kappait;
	pars->im95=pars->it95;
        pars->ima=pars->ita;
        pars->imb=pars->itb;

	printf("\nParameters for the interrupted alternate time gamma distribution:\n");
	printf("pim:\t%22.15e\n",pars->pim);
	printf("imbar:\t%22.15e\n",pars->imbar);
	printf("kappaim:%22.15e\n",pars->kappaim);
	printf("im95:\t%22.15e\n",pars->im95);
        printf("ima:\t%22.15e\n",pars->ima);
        printf("imb:\t%22.15e\n",pars->imb);

      } else {

	if(isnan(pars->imbar)) pars->imbar=pars->itbar;

	if((isnan(pars->kappaim)==0) + (isnan(pars->im95)==0) != 1) {
	  fprintf(stderr,"%s: Error: Either the kappaim parameter or the im95 parameter must be provided.\n",__func__);
	  return -9;
	}

	if(model_solve_gamma_group(&pars->imbar, &pars->kappaim, &pars->im95)) {
	  fprintf(stderr,"%s: Error: Cannot solve parameters for the interrupted alternate time gamma distribution\n",__func__);
	  return -10;

	} else {
          pars->ima=pars->imbar*pars->kappaim;
          pars->imb=1/pars->kappaim;
	  printf("\nParameters for the interrupted alternate time gamma distribution:\n");
	  printf("pim:\t%22.15e\n",pars->pim);
	  printf("imbar:\t%22.15e\n",pars->imbar);
	  printf("kappaim:%22.15e\n",pars->kappaim);
	  printf("im95:\t%22.15e\n",pars->im95);
          printf("ima:\t%22.15e\n",pars->ima);
          printf("imb:\t%22.15e\n",pars->imb);
	}
      }
    }
  }

  if(!isnan(pars->kappal) || !isnan(pars->l95)) {

    if((isnan(pars->kappal)==0) + (isnan(pars->l95)==0) != 1) {
      fprintf(stderr,"%s: Error: Either the kappal parameter or the l95 parameter must be provided.\n",__func__);
      return -11;
    }

    if(model_solve_gamma_group(&pars->lbar, &pars->kappal, &pars->l95)) {
      fprintf(stderr,"%s: Error: Cannot solve parameters for the latent time gamma distribution\n",__func__);
      return -12;

    } else {
      pars->la=pars->lbar*pars->kappal;
      pars->lb=1/pars->kappal;
      printf("\nParameters for the latent time gamma distribution:\n");
      printf("lbar:\t%22.15e\n",pars->lbar);
      printf("kappal:\t%22.15e\n",pars->kappal);
      printf("l95:\t%22.15e\n",pars->l95);
      printf("la:\t%22.15e\n",pars->la);
      printf("lb:\t%22.15e\n",pars->lb);
    }
  }

  printf("\nBranching process effective reproduction number:\n");
  printf("brReff:\t%22.15e\n",pars->R0*(1+(isnan(pars->q)?0:pars->q*(pars->mbar/pars->tbar-1))));

  return 0;
}

int model_solve_R0_group(model_pars* pars)
{
  double l1mp;
  int ret;

  if(!isnan(pars->pinf) && (!(pars->pinf>=0) || !(pars->pinf<=1))) {
    fprintf(stderr,"%s: Error: The pinf parameter must have a value in the interval (0,1]\n",__func__);
    return -2;
  }

  if(!isnan(pars->tbar) && pars->tbar<=0) {
    fprintf(stderr,"%s: Error: tbar must be greater than 0\n",__func__);
    return -3;
  }

  if(!isnan(pars->lambda) && pars->lambda<=0) {
    fprintf(stderr,"%s: Error: lambda must be greater than 0\n",__func__);
    return -4;
  }

  if(!isnan(pars->R0) && pars->R0<=0) {
    fprintf(stderr,"%s: Error: R0j must be greater than 0\n",__func__);
    return -5;
  }

  //If g_ave, p or mu is provided as an input
  if(!isnan(pars->g_ave) || !isnan(pars->p) || !isnan(pars->mu)) {

    //If g_ave is provided
    if(!isnan(pars->g_ave)) {

      if(!(pars->g_ave>=2)) {
	fprintf(stderr,"%s: Error: g_ave must be greater than or equal to 2\n",__func__);
	return -1;
      }

      if(pars->grouptype&ro_group_log_attendees_plus_1) {
	pars->mu=pars->g_ave-1;
	//Solve for p numerically from mu
	ret=model_solve_p_from_mu(pars);

	if(ret) return ret;

      } else {
	//Solve for p numerically from g_ave
	ret=model_solve_p_from_mean(pars->g_ave,pars);

	if(ret) return ret;
	pars->mu=-pars->p/((1-pars->p)*log(1-pars->p));
      }

    } else {

      //Else if p is provided
      if(!isnan(pars->p)) {

	if(!(pars->p>=0) || !(pars->p<1)) {
	  fprintf(stderr,"%s: Error: p must be non-negative and smaller than 1\n",__func__);
	  return -1;
	}
	pars->mu=(pars->p>0?-pars->p/((1-pars->p)*log(1-pars->p)):1);

	//Else if it is mu that is provided
      } else {

	if(!(pars->mu>=1)) {
	  fprintf(stderr,"%s: Error: mu must be greater than or equal to 1\n",__func__);
	  return -1;
	}
	//Solve for p numerically from mu
	ret=model_solve_p_from_mu(pars);

	if(ret) return ret;
      }

      if(pars->grouptype&ro_group_log_attendees_plus_1) pars->g_ave=pars->mu+1;

      else if(pars->p == 0) pars->g_ave=2;
      
      else {
	l1mp=log(1-pars->p);
	pars->g_ave=-pars->p*pars->p/((1-pars->p)*(l1mp+pars->p));
      }
    }

    //Solve for the missing parameter

    if(isnan(pars->R0)) pars->R0=pars->lambda*pars->tbar*(pars->g_ave-1)*pars->pinf;

    else if(isnan(pars->lambda)) pars->lambda=pars->R0/(pars->tbar*(pars->g_ave-1)*pars->pinf);

    else if(isnan(pars->tbar)) pars->tbar=pars->R0/(pars->lambda*(pars->g_ave-1)*pars->pinf);

    else pars->pinf=pars->R0/(pars->lambda*pars->tbar*(pars->g_ave-1));

  //Else if g_ave, p and mu are unknown
  } else {
    pars->g_ave=pars->R0/(pars->lambda*pars->tbar*pars->pinf)+1;

    if(!(pars->g_ave>=2)) {
      fprintf(stderr,"%s: Error: Computed g_ave is smaller than 2. The choice of the other input parameters that were provided is invalid.\n",__func__);
      return -1;
    }

    if(pars->grouptype&ro_group_log_attendees_plus_1) {
      pars->mu=pars->g_ave-1;

      if(pars->mu==1) pars->p=0;

      else {
	//Solve for p numerically from mu
	ret=model_solve_p_from_mu(pars);

	if(ret) return ret;
      }

    } else {
      //Solve for p numerically from g_ave
      ret=model_solve_p_from_mean(pars->g_ave,pars);

      if(ret) return ret;
      pars->mu=-pars->p/((1-pars->p)*log(1-pars->p));
    }
  }
  return 0;
}

int model_solve_p_from_mu(model_pars* pars)
{
  if(pars->mu==1) pars->p=0;

  else {
    root_finder* rf=root_finder_init(logroot, &pars->mu);
    pars->p=0.999;
    double diff;

    int ret=root_finder_find(rf, RF_P_EPSF, 100, RF_P_EPSF, 1-RF_P_EPSF, &pars->p, &diff);

    root_finder_free(rf);

    if(ret) {

      if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,RF_P_EPSF);

      else if(ret==-2) {
	fprintf(stderr,"%s: Error: Root could not be found!\n",__func__);
	return ret;
      }
    }
  }

  return 0;
}

int model_solve_p_from_mean(const double mean, model_pars* pars)
{
  if(mean==2) pars->p=0;

  else {
    root_finder* rf=root_finder_init(loggt1root, (void*)&mean);
    pars->p=0.999;
    double diff;

    int ret=root_finder_find(rf, RF_P_EPSF, 100, RF_P_EPSF, 1-RF_P_EPSF, &pars->p, &diff);

    root_finder_free(rf);

    if(ret) {

      if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,RF_P_EPSF);

      else if(ret==-2) {
	fprintf(stderr,"%s: Error: Root could not be found!\n",__func__);
	return ret;
      }
    }
  }

  return 0;
}

int model_solve_gamma_group(double* ave, double* kappa, double* x95)
{
  if(!(*ave>=0)) {
    fprintf(stderr,"%s: Error: The average of the distribution must be non-negative.\n",__func__);
    return -1;
  }

  if(isnan(*x95)) {

    if(!(*kappa>=0)) {
      fprintf(stderr,"%s: Error: The kappa parameter of the distribution must have a positive value.\n",__func__);
      return -1;

    } else if(!(*kappa>1/ *ave)) fprintf(stderr,"%s: Warning: The selected kappa value will generate a monotonically decreasing distribution!.\n",__func__);

    if(*kappa != INFINITY) {
      double pars[2]={*ave * *kappa, *kappa};
      root_finder* rf=root_finder_init(gpercroot, pars);
      *x95=*ave;
      double diff;

      int ret=root_finder_find(rf, RF_GPERC_EPSF, 100, *ave, 1e100, x95, &diff);

      root_finder_free(rf);

      if(ret) {

	if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,RF_GPERC_EPSF);

	else if(ret==-2) {
	  fprintf(stderr,"%s: Error: Root could not be found!\n",__func__);
	  return ret;
	}
      }

    } else *x95 = *ave;

  } else {

    if(!(*x95>=*ave)) {
      fprintf(stderr,"%s: Error: The 95th percentile of the distribution cannot be smaller than the average\n",__func__);
      return -1;
    }

    if(*x95 != *ave) {
      *kappa=1;
      double otherkappa=*kappa*0.9;
      double pars[4]={*ave, *x95, otherkappa, gpercrootfunc(*ave * otherkappa, *x95 * otherkappa)};
      root_finder* rf=root_finder_init(gkapparoot, pars);
      double diff;

      int ret=root_finder_find(rf, RF_GKAPPA_EPSF, 100, 1/ *ave, 1e100, kappa, &diff);


      if(ret) {

	if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,RF_GKAPPA_EPSF);

	else if(ret==-2) {
	  fprintf(stderr,"%s: Warning: Root could not be found with a mode of the gamma distribution above 0. Now searching for a monotonically decreasing solution!\n",__func__);
	  *kappa=1./ *ave;
	  otherkappa=*kappa*0.9;
	  pars[2]=otherkappa;
	  pars[3]=gpercrootfunc(*ave * otherkappa, *x95 * otherkappa);

          ret=root_finder_find(rf, RF_GKAPPA_EPSF, 100, 0, 1/ *ave, kappa, &diff);

	  if(ret==-3) fprintf(stderr,"%s: Warning: Convergence seems to have been reached, but the root discrepancy (%22.15e) is larger than required (%22.15e)!\n",__func__,diff,RF_GKAPPA_EPSF);
	  else if(ret==-2) {
	    fprintf(stderr,"%s: Error: Root could not be found with a mode of the gamma distribution at 0.!\n",__func__);
	    root_finder_free(rf);
	    return ret;
	  }
	}
      }
      root_finder_free(rf);

    } else *kappa = INFINITY;
  }
  return 0;
}

int model_pars_check(model_pars const* pars)
{
  int ret=0;

  if(pars->lambdap<=0) {
    fprintf(stderr,"%s: Error: If defined, lambdap must be greater than 0\n",__func__);
    ret-=1;

  } else if(pars->lambdap>0 && pars->popsize==0) {
    fprintf(stderr,"%s: Error: lambdap cannot be used with an infinite population\n",__func__);
    ret-=2;
  }

  if(pars->pit<0 || pars->pit>1) {
    fprintf(stderr,"%s: Error: pit must be in the interval [0.1]\n",__func__);
    ret-=4;
  }

  if(pars->q<0 || pars->q>1) {
    fprintf(stderr,"%s: Error: q must be in the [0,1] interval\n",__func__);
    ret-=8;

  } else if(pars->q>0) {

    if(pars->pim<0 || pars->pim>1) {
      fprintf(stderr,"%s: Error: pim must be in the interval [0.1]\n",__func__);
      ret-=16;
    }

    if(pars->q==1 && !(pars->pricommpertype&ro_pricommper_alt)) {
      fprintf(stderr,"%s: Error: Invalid configuration for the communicable period of the primary infectious individuals. The alternate communicable period distributions cannot be excluded if the probability for the alternate communicable period is 1\n",__func__);
      ret-=32;
    }

  } else if(!(pars->pricommpertype&ro_pricommper_main)) {
    fprintf(stderr,"%s: Error: Invalid configuration for the communicable period of the primary infectious individuals. The main communicable period distributions cannot be excluded if the probability for the alternate communicable period is 0\n",__func__);
    ret-=64;
  }

  if(!(pars->pricommpertype&(ro_pricommper_main|ro_pricommper_alt))) {
    fprintf(stderr,"%s: Error: Invalid configuration for the communicable period of the primary infectious individuals. Both communicable period distributions cannot be excluded\n",__func__);
    ret-=128;
  }

  if((pars->pricommpertype&ro_pricommper_main) && (pars->timetype==ro_time_pri_test_results)) {
    fprintf(stderr,"%s: Error: Invalid configuration for the combination of communicable period and time of the primary infectious individuals. Time relative to test results cannot be used if the main communicable period is allowed for primary infectious individuals\n",__func__);
    ret-=256;
  }

  //If testing is performed
  if(!isnan(pars->tdeltat)) {

    if(!(pars->tdeltat>=0)) {
      fprintf(stderr,"%s: Error: A tdelta value larger or equal to 0 must be defined\n",__func__);
      ret-=512;
    }

    //If alternate communicable period exists
    if(pars->q > 0) {

      //The true positive rate for the alternate communicable period must be
      //defined
      if(!(pars->mtpr>=0) || !(pars->mtpr<=1)) {
	fprintf(stderr,"%s: Error: An mtpr value in the interval [0,1] must be defined\n",__func__);
	ret-=1024;
      }
    }

    //If the main communicable period can be interrupted
    if(!isnan(pars->pit)) {

      //The true positive rate for the main communicable period must be
      //defined
      if(!(pars->ttpr>=0) || !(pars->ttpr<=1)) {
	fprintf(stderr,"%s: Error: A ttpr value in the interval [0,1] must be defined\n",__func__);
	ret-=2048;
      }
    }
  }

  if(pars->pit>0 || pars->pim>0 || !isnan(pars->tdeltat)) {

    if(!(pars->ttpr>=0) || !(pars->ttpr<=1)) {
      fprintf(stderr,"%s: Error: A ttpr value in the interval [0,1] must be defined\n",__func__);
      ret-=4096;
    }

    if(!(pars->tdeltat>=0)) {
      fprintf(stderr,"%s: Error: A tdelta value larger or equal to 0 must be defined\n",__func__);
      ret-=8192;
    }
  }

  if(isnan(pars->mtpr) && !(pars->pricommpertype&ro_pricommper_alt_use_tpr)) {
    fprintf(stderr,"%s: Error: pri_no_alt_test_fnr cannot be used if testing is not activated for the alternate communicable period\n",__func__);
    ret-=8192;
  }

  if(pars->tmax<=0) {
    fprintf(stderr,"%s: Error: tmax must be greater than 0\n",__func__);
    ret-=16384;
  }

  if(pars->nstart<=0) {
    fprintf(stderr,"%s: Error: nstart must be greater than 0\n",__func__);
    ret-=32768;
  }

  if(pars->popsize==0 && (pars->grouptype&ro_group_log_invitees)) {
    fprintf(stderr,"%s: Error: If modeling an infinite population, the groups of individuals cannot be generated using a logarithmically-distributed number of invitees\n",__func__);
    ret-=65536;
  }

  return ret;
}
