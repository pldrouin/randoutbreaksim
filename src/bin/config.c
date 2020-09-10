/**
 * @file config.c
 * @brief Configuration functions the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "config.h"

int config(model_pars* pars, bool* ninfhist, uint32_t* npaths, uint32_t* nthreads, uint32_t* nsetsperthread, uint32_t* nimax, int* oout, int* eout, const int nargs, const char* args[])
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

      } else if(!argsdiffer(pbuf, "kappa")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappa);

      } else if(!argsdiffer(pbuf, "t95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->t95);

      } else if(!argsdiffer(pbuf, "lambda")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->lambda);

      } else if(!argsdiffer(pbuf, "lambdap")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->lambdap);

      } else if(!argsdiffer(pbuf, "group_log_plus_1")) {
	pars->grouptype=ro_group_log_plus_1;

      } else if(!argsdiffer(pbuf, "group_log_attendees")) {
	pars->grouptype=ro_group_log_attendees;

      } else if(!argsdiffer(pbuf, "group_log_invitees")) {
	pars->grouptype=ro_group_log_invitees;

      } else if(!argsdiffer(pbuf, "p")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->p);

      } else if(!argsdiffer(pbuf, "mu")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->mu);

      } else if(!argsdiffer(pbuf, "pinf")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->pinf);

      } else if(!argsdiffer(pbuf, "popsize")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&pars->popsize);

      } else if(!argsdiffer(pbuf, "R0")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->R0);

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

      } else if(!argsdiffer(pbuf, "pit")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->pit);

      } else if(!argsdiffer(pbuf, "itbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->itbar);

      } else if(!argsdiffer(pbuf, "kappait")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappait);

      } else if(!argsdiffer(pbuf, "it95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->it95);

      } else if(!argsdiffer(pbuf, "pim")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->pim);

      } else if(!argsdiffer(pbuf, "imbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->imbar);

      } else if(!argsdiffer(pbuf, "kappaim")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappaim);

      } else if(!argsdiffer(pbuf, "im95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->im95);

      } else if(!argsdiffer(pbuf, "tmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->tmax);

      } else if(!argsdiffer(pbuf, "nstart")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&pars->nstart);

      } else if(!argsdiffer(pbuf, "ninfhist")) {
	*ninfhist=true;

      } else if(!argsdiffer(pbuf, "npaths")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,npaths);

      } else if(!argsdiffer(pbuf, "nthreads")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,nthreads);

      } else if(!argsdiffer(pbuf, "nsetsperthread")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,nsetsperthread);

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
  printf("\t--tbar VALUE\t\tMean main communicable period\n");
  printf("\t--kappa VALUE\t\tkappa parameter for the gamma distribution used to generate the main communicable period\n");
  printf("\t--t95 VALUE\t\t95th percentile of the main communicable period\n");
  printf("\t--lambda VALUE\t\tRate of events for a given individual. Events are defined to have more than one invitees\n");
  //printf("\t--lambdap VALUE\t\tTotal rate of events for a finite population. Events are defined to have more than one invitees\n");
  printf("\t--group_log_plus_1\tNumber of individuals in an event to be distributed as a logarithmically-distributed variable plus 1 (default)\n");
  printf("\t--group_log_attendees\tNumber of attendees in an event to be distributed as a logarithmically-distributed variable\n");
  printf("\t--group_log_invitees\tNumber of invitees for an event to be distributed as a logarithmically-distributed variable\n");
  printf("\t--p VALUE\t\tParameter for the logarithmic distribution used to draw the number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (0 <= p < 1)\n");
  printf("\t--mu VALUE\t\tParameter for the mean of an unbounded logarithmic distribution used to draw number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (mu=-1/log(1-p)*p/(1-p), mu >= 1)\n");
  printf("\t--pinf VALUE\t\tProbability that a given susceptible individual gets infected when exposed to one infectious individual during one event (default value of 1)\n");
  //printf("\t--popsize VALUE\t\tPopulation size (default value of 0, for an infinite population)\n");
  printf("\t--R0 VALUE\t\tBasic reproduction number\n");
  printf("\t--lbar VALUE\t\tMean latent period (default value of 0)\n");
  printf("\t--kappal VALUE\t\tkappa parameter for the gamma distribution used to generate the latent period\n");
  printf("\t--l95 VALUE\t\t95th percentile of the latent period\n");
  printf("\t--q VALUE\t\tProbability of alternate communicable period (default value of 0)\n");
  printf("\t--mbar VALUE\t\tMean period for the alternate communicable period (required if q>0)\n");
  printf("\t--kappaq VALUE\t\tkappa parameter for the gamma distribution used to generate the alternate communicable period\n");
  printf("\t--m95 VALUE\t\t95th percentile of the alternate communicable period\n");
  printf("\t--pit VALUE\t\tProbability of main communicable period interruption (default value of 0)\n");
  printf("\t--itbar VALUE\t\tMean period for the interrupted main communicable period (required if pit>0)\n");
  printf("\t--kappait VALUE\t\tkappa parameter for the gamma period used to generate the interrupted main communicable period\n");
  printf("\t--it95 VALUE\t\t95th percentile of the interrupted main communicable period\n");
  printf("\t--pim VALUE\t\tProbability of alternate communicable period interruption (default value of pit)\n");
  printf("\t--imbar VALUE\t\tMean period for the interrupted alternate communicable period (default value of itbar)\n");
  printf("\t--kappaim VALUE\t\tkappa parameter for the gamma period used to generate the interrupted alternate communicable period (default value of kappait)\n");
  printf("\t--im95 VALUE\t\t95th percentile of the interrupted alternate communicable period (default value of it95)\n");
  printf("\t--tmax VALUE\t\tMaximum simulation period used to instantiate new infectious individuals (default value of INFINITY)\n");
  printf("\t--nstart VALUE\t\tInitial number of infectious individuals (default value of 1)\n");
  printf("\t--ninfhist\t\tCompute a histogram of the number of infected individuals for each infectious individual.\n");
  printf("\t--npaths VALUE\t\tNumber of generated simulation paths (default value of 10000)\n");
  printf("\t--nthreads VALUE\tNumber of threads used to perform the simulation (default value of 1)\n");
  printf("\t--nsetsperthread VALUE\tNumber of path sets used for each thread (default value of 100 when nthreads>1, and of 1 otherwise)\n");
  printf("\t--nimax VALUE\t\tMaximum number of infectious individuals for a given time integer interval (default value of UINT32_MAX)\n");
  printf("\t--help\t\t\tPrint this usage information and exit\n");
  printf("\nEach option can be used as shown above from the command line. Dash(es) for option names are optional. For configuration files, '=', ':' or spaces as defined by isspace() can be used to separate option names from arguments. Characters following '#' on one line are considered to be comments.\nOptions can be used multiple times and configuration files can be read from configuration files.\n"); 
}
