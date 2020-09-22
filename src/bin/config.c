/**
 * @file config.c
 * @brief Configuration functions the simulation executable.
 * @author <Pierre-Luc.Drouin@drdc-rddc.gc.ca>, Defence Research and Development Canada Ottawa Research Centre.
 */

#include "config.h"

int config(config_pars* cp, const int nargs, const char* args[])
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

	if((cp->oout=open(pbuf,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}
	dup2(cp->oout,STDOUT_FILENO);

      } else if(!argsdiffer(pbuf, "elog")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	fflush(stderr);

	if((cp->eout=open(pbuf,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}
	dup2(cp->eout,STDERR_FILENO);

      } else if(!argsdiffer(pbuf, "tbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.tbar);

      } else if(!argsdiffer(pbuf, "kappa")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.kappa);

      } else if(!argsdiffer(pbuf, "t95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.t95);

      } else if(!argsdiffer(pbuf, "lambda")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.lambda);

      } else if(!argsdiffer(pbuf, "lambdap")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.lambdap);

      } else if(!argsdiffer(pbuf, "group_log_attendees_plus_1")) {
	cp->pars.grouptype=ro_group_log_attendees_plus_1;

      } else if(!argsdiffer(pbuf, "group_log_attendees")) {
	cp->pars.grouptype=ro_group_log_attendees;

      } else if(!argsdiffer(pbuf, "group_log_invitees")) {
	cp->pars.grouptype=ro_group_log_invitees;

      } else if(!argsdiffer(pbuf, "p")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.p);

      } else if(!argsdiffer(pbuf, "mu")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.mu);

      } else if(!argsdiffer(pbuf, "pinf")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pinf);

      } else if(!argsdiffer(pbuf, "popsize")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->pars.popsize);

      } else if(!argsdiffer(pbuf, "R0")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.R0);

      } else if(!argsdiffer(pbuf, "lbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.lbar);

      } else if(!argsdiffer(pbuf, "kappal")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.kappal);

      } else if(!argsdiffer(pbuf, "l95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.l95);

      } else if(!argsdiffer(pbuf, "q")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.q);

      } else if(!argsdiffer(pbuf, "mbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.mbar);

      } else if(!argsdiffer(pbuf, "kappaq")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.kappaq);

      } else if(!argsdiffer(pbuf, "m95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.m95);

      } else if(!argsdiffer(pbuf, "pit")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pit);

      } else if(!argsdiffer(pbuf, "itbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.itbar);

      } else if(!argsdiffer(pbuf, "kappait")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.kappait);

      } else if(!argsdiffer(pbuf, "it95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.it95);

      } else if(!argsdiffer(pbuf, "pim")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pim);

      } else if(!argsdiffer(pbuf, "imbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.imbar);

      } else if(!argsdiffer(pbuf, "kappaim")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.kappaim);

      } else if(!argsdiffer(pbuf, "im95")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.im95);

      } else if(!argsdiffer(pbuf, "ttpr")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.ttpr);

      } else if(!argsdiffer(pbuf, "mtpr")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.mtpr);

      } else if(!argsdiffer(pbuf, "tdeltat")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.tdeltat);

      } else if(!argsdiffer(pbuf, "pri_no_main_period")) {
	cp->pars.pricommpertype&=~ro_pricommper_main;

      } else if(!argsdiffer(pbuf, "pri_no_alt_period")) {
	cp->pars.pricommpertype&=~ro_pricommper_alt;

      } else if(!argsdiffer(pbuf, "pri_no_alt_test_fnr")) {
	cp->pars.pricommpertype&=~ro_pricommper_alt_use_tpr;

      } else if(!argsdiffer(pbuf, "time_rel_pri_end")) {
	cp->pars.trelpriend=true;

      } else if(!argsdiffer(pbuf, "tmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.tmax);

      } else if(!argsdiffer(pbuf, "nstart")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->pars.nstart);

      } else if(!argsdiffer(pbuf, "tlout")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);

	if((cp->tlout=open(pbuf,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}

      } else if(!argsdiffer(pbuf, "tloutbufsize")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->tloutbufsize);

      } else if(!argsdiffer(pbuf, "ninfhist")) {
	cp->ninfhist=true;

      } else if(!argsdiffer(pbuf, "npaths")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->npaths);

      } else if(!argsdiffer(pbuf, "nthreads")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nthreads);

      } else if(!argsdiffer(pbuf, "nsetsperthread")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nsetsperthread);

      } else if(!argsdiffer(pbuf, "lmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->lmax);

      } else if(!argsdiffer(pbuf, "nimax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nimax);

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
  printf("Stochastic simulation of outbreaks, using gamma distributions for the different time periods and a Poisson distribution for the number of interaction events where transmission can occur.\n");
  printf("\nBasic reproduction parameters:\n\n");
  printf("\tThe basic reproduction number R0 is defined by the expression\n");
  printf("\t\tR0 = lambda * tbar * (g_ave - 1) * pinf.\n");
  printf("\n\tA sufficient number of input parameters must be provided to determine, without overdetermining, the above expression.\n");
  printf("\tThe mu and p parameters are alternate parameters that can be provided instead of g_ave.\n");
  printf("\tmu is the mean of an unbounded logarithmic distribution with parameter p (mu = -p / ((1 - p) * log(1 - p))).\n");
  printf("\tThe expression of g_ave as a function of p depends on the type of group distribution that is selected for the events.\n");
  printf("\tAn event is defined to include at least two invitees.\n");
  printf("\n\t--group_log_attendees_plus_1, the default distribution from branchsim, indicates that the number of attendees in an event is to be distributed according to a logarithmically-distributed variable plus 1. For an infinite population, a fixed communicable period, and when pinf=1, this results in the total number of infections from a given infectious individual to follow a negative binomial distribution. The expression for g_ave with this distribution is\n");
  printf("\t\tg_ave = mu + 1.\n");

  printf("\n\t--group_log_attendees indicates that the number of attendees in an event is to be distributed according to a logarithmically-distributed variable (truncated below 2). In this case, it is the distribution of the number of individuals in a group that is motivated from empirical evidence, instead of the distribution for the total number of infections from a given infectious individual. The expression for g_ave with this distribution is\n");
  printf("\t\tg_ave = -p * p / ((1 - p) * (log(1 - p) + p)).\n");
  printf("\n\n");

  printf("Options\n\n");
  printf("\t--config FILENAME\t\tRead configuration options from FILENAME.\n");
  printf("\t--olog FILENAME\t\t\tRedirect standard output to FILENAME.\n");
  printf("\t--elog FILENAME\t\t\tRedirect standard error to FILENAME.\n");
  printf("\t--tbar VALUE\t\t\tMean main communicable period.\n");
  printf("\t--kappa VALUE\t\t\tkappa parameter for the gamma distribution used to generate the main communicable period.\n");
  printf("\t--t95 VALUE\t\t\t95th percentile of the main communicable period.\n");
  printf("\t--lambda VALUE\t\t\tRate of events for a given individual. Events are defined to include at least two invitees.\n");
  //printf("\t--lambdap VALUE\t\t\tTotal rate of events for a finite population. Events are defined to include at least two invitees.\n");
  printf("\t--group_log_attendees_plus_1\tNumber of individuals in an event to be distributed as a logarithmically-distributed variable plus 1 (default).\n");
  printf("\t--group_log_attendees\t\tNumber of attendees in an event to be distributed as a logarithmically-distributed variable truncated below 2.\n");
  printf("\t--group_log_invitees\t\tNumber of invitees for an event to be distributed as a logarithmically-distributed variable truncated below 2.\n");
  printf("\t--g_ave VALUE\t\t\tParameter for the average group size for one event. These individuals can correspond to invitees or attendees depending on the choice of group type. Events are defined to include at least two invitees (g_ave>=2).\n");
  printf("\t--p VALUE\t\t\tParameter for the logarithmic distribution used to draw the number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (0 <= p < 1).\n");
  printf("\t--mu VALUE\t\t\tParameter for the mean of an unbounded logarithmic distribution used to draw number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (mu >= 1).\n");
  printf("\t--pinf VALUE\t\t\tProbability that a given susceptible individual gets infected when exposed to one infectious individual during one event.\n");
  //printf("\t--popsize VALUE\t\t\tPopulation size (default value of 0, for an infinite population).\n");
  printf("\t--R0 VALUE\t\t\tBasic reproduction number.\n");
  printf("\t--lbar VALUE\t\t\tMean latent period (default value of 0).\n");
  printf("\t--kappal VALUE\t\t\tkappa parameter for the gamma distribution used to generate the latent period.\n");
  printf("\t--l95 VALUE\t\t\t95th percentile of the latent period.\n");
  printf("\t--q VALUE\t\t\tProbability of alternate communicable period (default value of 0).\n");
  printf("\t--mbar VALUE\t\t\tMean period for the alternate communicable period (required if q>0).\n");
  printf("\t--kappaq VALUE\t\t\tkappa parameter for the gamma distribution used to generate the alternate communicable period.\n");
  printf("\t--m95 VALUE\t\t\t95th percentile of the alternate communicable period.\n");
  printf("\t--pit VALUE\t\t\tProbability of main communicable period interruption (default value of 0).\n");
  printf("\t--itbar VALUE\t\t\tMean period for the interrupted main communicable period (required if pit>0).\n");
  printf("\t--kappait VALUE\t\t\tkappa parameter for the gamma period used to generate the interrupted main communicable period.\n");
  printf("\t--it95 VALUE\t\t\t95th percentile of the interrupted main communicable period.\n");
  printf("\t--pim VALUE\t\t\tProbability of alternate communicable period interruption (default value of pit).\n");
  printf("\t--imbar VALUE\t\t\tMean period for the interrupted alternate communicable period (default value of itbar).\n");
  printf("\t--kappaim VALUE\t\t\tkappa parameter for the gamma period used to generate the interrupted alternate communicable period (default value of kappait).\n");
  printf("\t--im95 VALUE\t\t\t95th percentile of the interrupted alternate communicable period (default value of it95).\n");
  printf("\t--ttpr VALUE\t\t\tTrue positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the main period, for which a positive test would allow for the interruption of a child's communicable period.\n");
  printf("\t--mtpr VALUE\t\t\tTrue positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the alternate period, for which a positive test would allow for the interruption of a child's communicable period.\n");
  printf("\t--tdeltat VALUE\t\t\tTime delay between the end of the applicable communicable period and test results.\n");
  printf("\t--pri_no_main_period\t\tThe communicable period for a primary infectious individual cannot be the main period.\n");
  printf("\t--pri_no_alt_period\t\tThe communicable period for a primary infectious individual cannot be the alternate period.\n");
  printf("\t--pri_no_alt_test_fnr\t\tThe alternate communicable period for a primary infectious individual cannot result in a false negative test.\n");

  printf("\t--time_rel_pri_end\t\tRecorded event time is relative to the end of the communicable period for the generated primary infectious individuals.\n");
  printf("\t--tmax VALUE\t\t\tMaximum simulation period used to instantiate new infectious individuals (default value of INFINITY).\n");
  printf("\t--nstart VALUE\t\t\tInitial number of infectious individuals (default value of 1).\n");
  printf("\t--tlout VALUE\t\t\tOutput timeline information for each simulated path into the provided file in the binary format as described below.\n");
  printf("\t--tloutbufsize VALUE\t\tPer-thread memory buffer size (in MB) used to accumulate data for timeline output before writing them to disk (default value of 10 MB).\n");
  printf("\t--ninfhist\t\t\tCompute a histogram of the number of infected individuals for each infectious individual.\n");
  printf("\t--npaths VALUE\t\t\tNumber of generated simulation paths (default value of 10000).\n");
  printf("\t--nthreads VALUE\t\tNumber of threads used to perform the simulation (default value of 1).\n");
  printf("\t--nsetsperthread VALUE\t\tNumber of path sets used for each thread (default value of 100 when nthreads>1, and of 1 otherwise).\n");
  printf("\t--lmax VALUE\t\t\tMaximum number of layers (generations) for the simulation (value of 1 signifies only primary individuals, default value of UINT32_MAX).\n");
  printf("\t--nimax VALUE\t\t\tMaximum number of infectious individuals for a given time integer interval (default value of UINT32_MAX).\n");
  printf("\t--help\t\t\t\tPrint this usage information and exit.\n");
  printf("\nEach option can be used as shown above from the command line. Dash(es) for option names are optional. For configuration files, '=', ':' or spaces as defined by isspace() can be used to separate option names from arguments. Characters following '#' on one line are considered to be comments.\nOptions can be used multiple times and configuration files can be read from configuration files.\n"); 

  printf("\nBinary timeline output file:\n");
  printf("\n\tAll fields are stored in little endian.\n");
  printf("\n\tFile header:\n");
  printf("\t\t-Unsigned 32 bit value: floor(tmax)+1, the number of time bins starting from t=0.\n");
  printf("\t\t-Unsigned 8 bit field: value of 1 if time_rel_pri_end is used, and 0 otherwise.\n");

  printf("\n\tSimulation path records:\n");
  printf("\t\t-Unsigned 32 bit value: The number of written successive time bins.\n");
  printf("\t\t-Unsigned 32 bit value: Field is written only if time_rel_pri_end is used. Value is the number of time bins before t=0.\n");
  printf("\t\t-Unsigned 8 bit value: value of 1 if the path goes extinct, and 0 otherwise.\n");
  printf("\n\t\t-For each time bin, chronologically written:\n");
  printf("\t\t\t-Unsigned 32 bit value: Number of active infections.\n");
  printf("\t\t\t-Unsigned 32 bit value: Number of new infections.\n");
}
