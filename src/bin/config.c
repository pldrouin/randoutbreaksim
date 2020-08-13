/**
 * @file config.c
 * @brief Configuration functions the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "config.h"

int config(sim_pars* pars, uint32_t* npaths, uint32_t* nimax, int* oout, int* eout, const int nargs, const char* args[])
{
  int plength=1;
  FILE **fptra=NULL;
  int fptri=-1;
  char pbuf[1024];
  int parc=0;

  if(!nargs) {
    printusage(args[-1]);
    return -1;
  }

  while(plength>0) {

    while((plength=getnextparam(fptra,&fptri,false,nargs,args,&parc,pbuf))>0) {

      if(!argsdiffer(pbuf, "config")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);

	++fptri;
	fptra=(FILE**)realloc(fptra,(fptri+1)*sizeof(FILE*));

	if(!(fptra[fptri]=fopen(pbuf,"r"))) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in read mode\n",__func__,pbuf);
	  return -1;
	}

      } else if(!argsdiffer(pbuf, "olog")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	fflush(stdout);

	if((*oout=open(pbuf,O_RDWR|O_CREAT|O_APPEND,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}
	dup2(*oout,STDOUT_FILENO);

      } else if(!argsdiffer(pbuf, "elog")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	fflush(stderr);

	if((*eout=open(pbuf,O_RDWR|O_CREAT|O_APPEND,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}
	dup2(*eout,STDERR_FILENO);

      } else if(!argsdiffer(pbuf, "tbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->tbar);

      } else if(!argsdiffer(pbuf, "p")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->p);

      } else if(!argsdiffer(pbuf, "lambda")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->lambda);

      } else if(!argsdiffer(pbuf, "kappa")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappa);

      } else if(!argsdiffer(pbuf, "t95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->t95);

      } else if(!argsdiffer(pbuf, "lbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->lbar);

      } else if(!argsdiffer(pbuf, "kappal")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappal);

      } else if(!argsdiffer(pbuf, "l95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->l95);

      } else if(!argsdiffer(pbuf, "q")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->q);

      } else if(!argsdiffer(pbuf, "mbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->mbar);

      } else if(!argsdiffer(pbuf, "kappaq")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappaq);

      } else if(!argsdiffer(pbuf, "m95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->m95);

      } else if(!argsdiffer(pbuf, "R0")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->R0);

      } else if(!argsdiffer(pbuf, "mu")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->mu);

      } else if(!argsdiffer(pbuf, "tmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->tmax);

      } else if(!argsdiffer(pbuf, "nstart")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&pars->nstart);

      } else if(!argsdiffer(pbuf, "npaths")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,npaths);

      } else if(!argsdiffer(pbuf, "nimax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,nimax);

      } else {

	if(argsdiffer(pbuf, "help")) fprintf(stderr,"%s: Error: Option '%s' is unknown\n",__func__,pbuf);
	printusage(args[-1]);
	return -1;
      }
    }
    fflush(stdout);
    fflush(stderr);
  }
  free(fptra);
  return 0;
}

void printusage(const char* name)
{
  fprintf(stderr,"Usage: %s [OPTION]\n\n",name);
  printf("Options\n");
  printf("\t--config FILENAME\tRead configuration options from FILENAME\n");
  printf("\t--olog FILENAME\t\tRedirect standard output to FILENAME\n");
  printf("\t--elog FILENAME\t\tRedirect standard error to FILENAME\n");
  printf("\t--tbar VALUE\t\tbranchsim's tbar parameter\n");
  printf("\t--p VALUE\t\tbranchsim's p parameter\n");
  printf("\t--lambda VALUE\t\tbranchsim's lambda parameter\n");
  printf("\t--kappa VALUE\t\tbranchsim's kappa parameter\n");
  printf("\t--t95 VALUE\t\tbranchsim's t95 parameter\n");
  printf("\t--lbar VALUE\t\tlbar parameter (default value of 0)\n");
  printf("\t--kappal VALUE\t\tkappal parameter (required if lbar>0)\n");
  printf("\t--l95 VALUE\t\tbranchsim's l95 parameter\n");
  printf("\t--q VALUE\t\tbranchsim's q parameter (default value of 0)\n");
  printf("\t--mbar VALUE\t\tbranchsim's mbar parameter (required if q>0)\n");
  printf("\t--kappaq VALUE\t\tbranchsim's kappaq parameter (required if q>0)\n");
  printf("\t--m95 VALUE\t\tbranchsim's m95 parameter\n");
  printf("\t--R0 VALUE\t\tbranchsim's R0 parameter\n");
  printf("\t--mu VALUE\t\tbranchsim's mu parameter\n");
  printf("\t--tmax VALUE\t\tbranchsim's tmax parameter (default value of INFINITY)\n");
  printf("\t--nstart VALUE\t\tbranchsim's nstart parameter (default value of 1)\n");
  printf("\t--npaths VALUE\t\tbranchsim's npaths parameter (default value of 10000)\n");
  printf("\t--nimax VALUE\t\tMaximum number of infectious individuals for a given time integer interval (default value of UINT32_MAX)\n");
  printf("\t--help\t\t\tPrint this usage information and exit\n");
  printf("\nEach option can be used as shown above from the command line. Dash(es) for option names are optional. For configuration files, '=', ':' or spaces as defined by isspace() can be used to separate option names from arguments. Characters following '#' on one line are considered to be comments.\nOptions can be used multiple times and configuration files can be read from configuration files.\n"); 
}

int config_solve_pars(sim_pars* pars)
{
  if((isnan(pars->tbar)==0) + (isnan(pars->lambda)==0) + (isnan(pars->p)==0 || isnan(pars->mu)==0) + (isnan(pars->R0)==0) != 3) {
    fprintf(stderr,"%s: Error: An invalid combination of tbar, lambda, p, mu and R0 parameters was provided.\n",__func__);
    return -1;
  }

  if((isnan(pars->kappa)==0) + (isnan(pars->t95)==0) != 1) {
    fprintf(stderr,"%s: Error: Either the kappa parameter or the t95 parameter must be provided.\n",__func__);
    return -2;
  }

  if((isnan(pars->kappaq)==0) + (isnan(pars->m95)==0) == 2) {
    fprintf(stderr,"%s: Error: Either the kappaq parameter or the m95 parameter must be provided.\n",__func__);
    return -3;
  }

  if((isnan(pars->kappal)==0) + (isnan(pars->l95)==0) == 2) {
    fprintf(stderr,"%s: Error: Either the kappal parameter or the l95 parameter must be provided.\n",__func__);
    return -4;
  }

  if(config_solve_R0_group(pars)) return -5;
  
  printf("Basic reproduction parameters are:\n");
  printf("lambda:\t%22.15e\n",pars->lambda);
  printf("tbar:\t%22.15e\n",pars->tbar);
  printf("mu:\t%22.15e\n",pars->mu);
  printf("p:\t%22.15e\n",pars->p);
  printf("R0:\t%22.15e\n",pars->R0);

  if(config_solve_gamma_group(&pars->tbar, &pars->kappa, &pars->t95)) {
    fprintf(stderr,"%s: Error: Cannot solve parameters for the main time gamma distribution\n",__func__);
    return -6;

  }
  printf("Parameters for the main time gamma distribution:\n");
  printf("tbar:\t%22.15e\n",pars->tbar);
  printf("kappa:\t%22.15e\n",pars->kappa);
  printf("t95:\t%22.15e\n",pars->t95);

  int ret;

  if((ret=config_solve_gamma_group(&pars->mbar, &pars->kappaq, &pars->m95))) {

    if(ret<0) {
      fprintf(stderr,"%s: Error: Cannot solve parameters for the alternate time gamma distribution\n",__func__);
      return -7;
    }

  } else {
    printf("Parameters for the alternate time gamma distribution:\n");
    printf("mbar:\t%22.15e\n",pars->mbar);
    printf("kappaq:\t%22.15e\n",pars->kappaq);
    printf("m95:\t%22.15e\n",pars->m95);
  }

  if((ret=config_solve_gamma_group(&pars->lbar, &pars->kappal, &pars->l95))) {

    if(ret<0) {
      fprintf(stderr,"%s: Error: Cannot solve parameters for the latent time gamma distribution\n",__func__);
      return -8;
    }
    
  } else {
    printf("Parameters for the latent time gamma distribution:\n");
    printf("lbar:\t%22.15e\n",pars->lbar);
    printf("kappal:\t%22.15e\n",pars->kappal);
    printf("l95:\t%22.15e\n",pars->l95);
  }

  return 0;
}

int config_solve_R0_group(sim_pars* pars)
{
  //If p is provided as an input
  if(!isnan(pars->p)) pars->mu=-pars->p/((1-pars->p)*log(1-pars->p));

  //Solve for the missing parameter
  if(isnan(pars->R0)) pars->R0=pars->lambda*pars->tbar*pars->mu;

  else if(isnan(pars->lambda)) pars->lambda=pars->R0/(pars->tbar*pars->mu);

  else if(isnan(pars->tbar)) pars->tbar=pars->R0/(pars->lambda*pars->mu);

  else pars->mu=pars->R0/(pars->lambda*pars->tbar);

  //If p is unknown solve for it numerically
  if(isnan(pars->p)) {

    if(pars->mu > 1) {
      root_finder* rf=root_finder_init(logroot, logrootderiv, &pars->mu);
      pars->p=0.999;

      int ret=root_finder_find(rf, RF_P_EPS, INFINITY, 100, RF_P_EPS, 1-RF_P_EPS, &pars->p);

      root_finder_free(rf);

      if(ret) return ret;

    } else pars->p=0;
  }
  return 0;
}

int config_solve_gamma_group(double* ave, double* kappa, double* p95)
{
  if(*ave) {

    if(isnan(*p95)) {

      if(*kappa != INFINITY) {
	double pars[2]={*ave * *kappa, *kappa};
	root_finder* rf=root_finder_init(gpercroot, gpercrootderiv, pars);
	*p95=*ave;

	int ret=root_finder_find(rf, INFINITY, RF_GPERC_EPSF, 100, *ave, 1e100, p95);

	root_finder_free(rf);

	if(ret) return ret;

      } else *p95 = *ave;

    } else {

      if(*p95 != *ave) {
	double pars[2]={*ave, *p95};
	root_finder* rf=root_finder_init(gkapparoot, gkapparootderiv, pars);
	*kappa=1;

	int ret=root_finder_find(rf, INFINITY, RF_GKAPPA_EPSF, 100, 1e-100, 1e100, kappa);

	root_finder_free(rf);

	if(ret) return ret;

      } else *kappa = INFINITY;
    }
    return 0;
  }
  return 1;
}
