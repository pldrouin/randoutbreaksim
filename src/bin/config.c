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

      } else if(!argsdiffer(pbuf, "pinfpri")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pinfpri);

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

      } else if(!argsdiffer(pbuf, "lambda_uncut")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.lambda_uncut);

      } else if(!argsdiffer(pbuf, "lambdap")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.lambdap);

      } else if(!argsdiffer(pbuf, "group_attendees")) {
	cp->pars.grouptype=(cp->pars.grouptype&ro_group_dist_mask);

      } else if(!argsdiffer(pbuf, "group_invitees")) {
	cp->pars.grouptype=(cp->pars.grouptype&ro_group_dist_mask)|ro_group_invitees;

      } else if(!argsdiffer(pbuf, "group_interactions")) {
	cp->pars.groupinteractions=true;

      } else if(!argsdiffer(pbuf, "group_transmissions")) {
	cp->pars.groupinteractions=false;

      } else if(!argsdiffer(pbuf, "group_log_plus_1")) {
	cp->pars.grouptype=(cp->pars.grouptype&~ro_group_dist_mask)|ro_group_log_plus_1;

      } else if(!argsdiffer(pbuf, "group_log")) {
	cp->pars.grouptype=(cp->pars.grouptype&~ro_group_dist_mask)|ro_group_log;

      } else if(!argsdiffer(pbuf, "group_gauss")) {
	cp->pars.grouptype=(cp->pars.grouptype&~ro_group_dist_mask)|ro_group_gauss;

      } else if(!argsdiffer(pbuf, "g_ave")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.g_ave);

      } else if(!argsdiffer(pbuf, "p")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.p);

      } else if(!argsdiffer(pbuf, "mu")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.mu);

      } else if(!argsdiffer(pbuf, "sigma")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.sigma);

      } else if(!argsdiffer(pbuf, "rsigma")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.rsigma);

      } else if(!argsdiffer(pbuf, "pinf")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pinf);

#ifdef DUAL_PINF
      } else if(!argsdiffer(pbuf, "ppip")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.ppip);

      } else if(!argsdiffer(pbuf, "rpinfp")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.rpinfp);

      } else if(!argsdiffer(pbuf, "rpshedp")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.rpshedp);

      } else if(!argsdiffer(pbuf, "qp")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.qp);
#endif

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

#ifdef CT_OUTPUT
      } else if(!argsdiffer(pbuf, "ctwindow")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.ctwindow);

      } else if(!argsdiffer(pbuf, "pt")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&cp->pars.pt);
#endif

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

#ifdef DUAL_PINF
      } else if(!argsdiffer(pbuf, "pri_first_category_only")) {
	cp->pars.pricommpertype=(cp->pars.pricommpertype|ro_pricommper_first_cat) - (cp->pars.pricommpertype&ro_pricommper_second_cat);

      } else if(!argsdiffer(pbuf, "pri_second_category_only")) {
	cp->pars.pricommpertype=(cp->pars.pricommpertype|ro_pricommper_second_cat) - (cp->pars.pricommpertype&ro_pricommper_first_cat);
#endif

      } else if(!argsdiffer(pbuf, "time_rel_pri_created")) {
	cp->pars.timetype=ro_time_pri_created;

      } else if(!argsdiffer(pbuf, "time_rel_pri_infectious")) {
	cp->pars.timetype=ro_time_pri_infectious;

      } else if(!argsdiffer(pbuf, "time_rel_pri_end_comm")) {
	cp->pars.timetype=ro_time_pri_end_comm;

      } else if(!argsdiffer(pbuf, "time_rel_pri_test_results")) {
	cp->pars.timetype=ro_time_pri_test_results;

      } else if(!argsdiffer(pbuf, "time_rel_first_pos_test_results")) {
	cp->pars.timetype=ro_time_first_pos_test_results;

      } else if(!argsdiffer(pbuf, "time_rel_pri_flat_comm")) {
	cp->pars.timetype=ro_time_pri_flat_comm;

      } else if(!argsdiffer(pbuf, "include_all_paths")) {
	cp->pars.pathtype=ro_all_paths;

      } else if(!argsdiffer(pbuf, "observable_paths_only")) {
	cp->pars.pathtype=ro_observable_paths_only;

      } else if(!argsdiffer(pbuf, "non-observable_paths_only")) {
	cp->pars.pathtype=ro_non_observable_paths_only;

      } else if(!argsdiffer(pbuf, "tmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%i",&cp->pars.tmax);

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

#ifdef CT_OUTPUT
      } else if(!argsdiffer(pbuf, "ctout")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);

	if((cp->ctout=open(pbuf,O_RDWR|O_CREAT|O_TRUNC,S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH))<0) {
	  fprintf(stderr,"%s: Error: Cannot open file '%s' in write mode\n",__func__,pbuf);
	  return -1;
	}

      } else if(!argsdiffer(pbuf, "ctoutbufsize")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->ctoutbufsize);
#endif

      } else if(!argsdiffer(pbuf, "ninfhist")) {
	cp->ninfhist=true;

      } else if(!argsdiffer(pbuf, "npaths")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->npaths);

      } else if(!argsdiffer(pbuf, "lmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->lmax);

      } else if(!argsdiffer(pbuf, "nbinsperunit")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIi32,&cp->nbinsperunit);

      } else if(!argsdiffer(pbuf, "nimax")) {

	if(cp->npostestmax!=UINT32_MAX) {
	  fprintf(stderr,"%s: Error: nimax and npostestmax cannot be both used at the same time\n",__func__);
	  return -1;
	}
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nimax);

      } else if(!argsdiffer(pbuf, "npostestmax")) {

	if(cp->nimax!=UINT32_MAX) {
	  fprintf(stderr,"%s: Error: nimax and npostestmax cannot be both used at the same time\n",__func__);
	  return -1;
	}
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->npostestmax);

      } else if(!argsdiffer(pbuf, "npostestmaxnunits")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->npostestmaxnunits);

      } else if(!argsdiffer(pbuf, "nthreads")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nthreads);

      } else if(!argsdiffer(pbuf, "nsetsperthread")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->nsetsperthread);

      } else if(!argsdiffer(pbuf, "stream")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&cp->stream);

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
  fprintf(stderr,"\nUSAGE: %s [OPTION]\n\n",name);
  printf("Stochastic simulation of outbreaks, using gamma distributions for the different time periods and a Poisson distribution for the number of interaction events where transmission can occur.\n");
  printf("\n\nBASIC REPRODUCTION PARAMETERS:\n\n");
  printf("\tThe basic reproduction number R0 is defined by the expression\n");
  printf("\t\tR0 = lambda * tbar * (g_ave - 1) * pinf,\n");
  printf("\tif group_transmissions is used, and\n");
  printf("\t\tR0 = lambda * tbar * (g_ave - 1 + g_sigma^2/g_ave) * pinf\n");
  printf("\tif group_interaction is used instead. R0 assumes an infinite population of susceptible individuals with a single infectious individual.\n");
  printf("\n\tA sufficient number of input parameters must be provided to determine, without overdetermining, the above expression.\n");
  printf("\tmu and p parameters are alternate parameters that can be provided instead of g_ave.\n");
  printf("\tmu is the mean of an unbounded logarithmic distribution with parameter p (mu = -p / ((1 - p) * log(1 - p))).\n");
  printf("\tThe expression of g_ave as a function of p depends on the type of group distribution that is selected for the events.\n");
  printf("\tAn event is defined to include at least two invitees.\n");
  printf("\n\t--group_log_plus_1, the default distribution from branchsim, indicates that the number of invitees/attendees in an event is to be distributed according to a logarithmically-distributed variable plus 1. For an infinite population, a fixed communicable period, and when pinf=1 and group_attendees are used, this results in the total number of infections from a given infectious individual to follow a negative binomial distribution. The expression for g_ave with this distribution is\n");
  printf("\t\tg_ave = mu + 1.\n");

  printf("\n\t--group_log indicates that the number of invitees/attendees in an event is to be distributed according to a logarithmically-distributed variable (truncated below 2). In this case, it is the distribution of the number of individuals in a group that is motivated from empirical evidence, instead of the distribution for the total number of infections from a given infectious individual. When using group_attendees, the expression for g_ave with this distribution is\n");
  printf("\t\tg_ave = -p * p / ((1 - p) * (log(1 - p) + p)).\n");

  printf("\n\nBRANCHING PROCESS EFFECTIVE REPRODUCTION NUMBER:\n\n");
  printf("\tIf an alternate communicable period of average duration mbar is defined, and if there is a probability q that an individual's communicable period be the alternate communicable instead of the main communicable period, then an effective reproduction number can be expressed as\n");
#ifdef DUAL_PINF
  printf("\t\tbrReff =  lambda * (g_ave - 1) * pinf * {(1 - ppip) * [(1 - q) * tbar + q * mbar] + rpshedp * ppip * rpinfp * [(1 - qp) * tbar + qp * mbar]}\n");
  printf("\t\t       =  R0 * {(1 - ppip) * [1 + q * ( mbar / tbar - 1)] + rpshedp * ppip * rpinfp * [1 + qp * ( mbar / tbar - 1)]}.\n\n");
#else
  printf("\t\tbrReff =  lambda * (g_ave - 1) * pinf * [(1 - q) * tbar + q * mbar]\n");
  printf("\t\t       =  R0 * [1 + q * ( mbar / tbar - 1)].\n\n");
#endif
  printf("\tThe expected effective reproduction number of the simulation will be given by the above expression if it consists of a branching process characterised by the model described above. For such a process, all generations of infections occur using the same static distributions. As identified below, some of the available options can make the simulation deviate from a branching process, in which case the effective reproduction number will deviate accordingly.\n");

  printf("\n\nOPTIONS\n\n");
  printf("\t--config FILENAME\t\tRead configuration options from FILENAME.\n");
  printf("\t--olog FILENAME\t\t\tRedirect standard output to FILENAME.\n");
  printf("\t--elog FILENAME\t\t\tRedirect standard error to FILENAME.\n");
  printf("\t--tbar VALUE\t\t\tMean main communicable period.\n");
  printf("\t--kappa VALUE\t\t\tkappa parameter for the gamma distribution used to generate the main communicable period.\n");
  printf("\t--t95 VALUE\t\t\t95th percentile of the main communicable period.\n");
  printf("\t--lambda VALUE\t\t\tRate of events for a given individual. Events are defined to include at least two invitees.\n");
  printf("\t--lambda_uncut VALUE\t\tRate of events for a given individual, including events of one invitee.\n");
  printf("\t--lambdap VALUE\t\t\tTotal rate of events for a finite population. Events are defined to include at least two invitees.\n");
  printf("\t--group_attendees\t\tThe group distributions are applicable to the number of attendees (default).\n");
  printf("\t--group_invitees\t\tThe group distributions are applicable to the number of invitees.\n");
  printf("\t--group_interactions\t\tThe group distribution is applicable to any interactions (no infectious individual required). This option is required for a finite population.\n");
  printf("\t--group_transmissions\t\tThe group distribution is applicable to interactions involving one infectious individual (default).\n");
  printf("\t--group_log_plus_1\t\tNumber of invitees/attendees in an event to be distributed as a logarithmically-distributed variable plus 1 (default).\n");
  printf("\t--group_log\t\t\tNumber of invitees/attendees in an event to be distributed as a logarithmically-distributed variable truncated below 2.\n");
  printf("\t--group_gauss\t\t\tNumber of invitees/attendees in an event to be distributed as a Gaussian-distributed variable truncated below 2.\n");
  printf("\t--g_ave VALUE\t\t\tParameter for the average group size for one event. These individuals can correspond to invitees or attendees depending on the choice of group type. The average group size for transmission events will be higher if group_interactions is used. Events are defined to include at least two invitees (g_ave>=2).\n");
  printf("\t--p VALUE\t\t\tParameter for the logarithmic distribution used to draw the number of individuals during one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type (0 <= p < 1).\n");
  printf("\t--mu VALUE\t\t\tParameter for the mean of an unbounded logarithmic distribution (mu >= 1) or of an unbounded Gaussian distribution used to draw number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.\n");
  printf("\t--sigma VALUE\t\t\tParameter for the standard deviation of an unbounded Gaussian used to draw the number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.\n");
  printf("\t--rsigma VALUE\t\t\tParameter for the relative standard deviation of an unbounded Gaussian used to draw the number of individuals for one event. These individuals can correspond to invitees, attendees or infected individuals depending on the choice of group type.\n");
  printf("\t--pinf VALUE\t\t\tProbability that a given susceptible individual gets infected when exposed to one infectious individual during one event.\n");
#ifdef DUAL_PINF
  printf("\t--ppip VALUE\t\t\tProbability that a given susceptible individual be in the second infection probability category (0 <= ppip <= 1, default value of 0).\n");
  printf("\t--rpinfp VALUE\t\t\tRelative probability that a given susceptible individual of the second category gets infected when exposed to one infectious individual during one event (value relative to pinf, 0 < rpinfp * pinf <= 1, default value of 1).\n");
  printf("\t--rpshedp VALUE\t\t\tRelative strength of infectiousness from an infectious individual of the second category vs the fist category (value relative to pinf, 0 < rpshedp * pinf <=1, default value of 1).\n");
  printf("\t--qp VALUE\t\t\tProbability of alternate communicable period for an infectious individual in the second category.\n");
#endif
  printf("\t--popsize VALUE\t\t\tPopulation size (default value of 0, for an infinite population).\n");
  printf("\t--R0 VALUE\t\t\tBasic reproduction number.\n");
  printf("\t--lbar VALUE\t\t\tMean latent period (default value of 0).\n");
  printf("\t--kappal VALUE\t\t\tkappa parameter for the gamma distribution used to generate the latent period.\n");
  printf("\t--l95 VALUE\t\t\t95th percentile of the latent period.\n");
  printf("\t--q VALUE\t\t\tProbability of alternate communicable period.\n");
  printf("\t--mbar VALUE\t\t\tMean period for the alternate communicable period (required if q>0).\n");
  printf("\t--kappaq VALUE\t\t\tkappa parameter for the gamma distribution used to generate the alternate communicable period.\n");
  printf("\t--m95 VALUE\t\t\t95th percentile of the alternate communicable period.\n");
#ifdef CT_OUTPUT
  printf("\t--ctwindow VALUE\t\tPeriod prior to individual isolation during which contacts are considered (default value of 0).\n");
  printf("\t--pt VALUE\t\t\tProbability of successful contact tracing. Probability must be larger than pit and pim, as it is considered to be applicable to all contacts.\n");
#endif
  printf("\t--pit VALUE\t\t\tProbability of main communicable period interruption This option makes a model diverge from a branching process.\n");
  printf("\t--itbar VALUE\t\t\tMean period for the interrupted main communicable period (required if pit>0).\n");
  printf("\t--kappait VALUE\t\t\tkappa parameter for the gamma period used to generate the interrupted main communicable period.\n");
  printf("\t--it95 VALUE\t\t\t95th percentile of the interrupted main communicable period.\n");
  printf("\t--pim VALUE\t\t\tProbability of alternate communicable period interruption (default value of pit). This option makes a model diverge from a branching process.\n");
  printf("\t--imbar VALUE\t\t\tMean period for the interrupted alternate communicable period (default value of itbar).\n");
  printf("\t--kappaim VALUE\t\t\tkappa parameter for the gamma period used to generate the interrupted alternate communicable period (default value of kappait).\n");
  printf("\t--im95 VALUE\t\t\t95th percentile of the interrupted alternate communicable period (default value of it95).\n");
  printf("\t--ttpr VALUE\t\t\tTrue positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the main period, for which a positive test would allow for the interruption of a child's communicable period.\n");
  printf("\t--mtpr VALUE\t\t\tTrue positive rate (= 1 - false negative rate) for the testing of a parent, whose communicable period is the alternate period, for which a positive test would allow for the interruption of a child's communicable period.\n");
  printf("\t--tdeltat VALUE\t\t\tTime delay between the end of the applicable communicable period and test results.\n");
  printf("\t--pri_no_main_period\t\tThe communicable period for a primary infectious individual cannot be the main period. This option makes a model diverge from a branching process.\n");
  printf("\t--pri_no_alt_period\t\tThe communicable period for a primary infectious individual cannot be the alternate period. This option makes a model diverge from a branching process.\n");
  printf("\t--pri_no_alt_test_fnr\t\tThe alternate communicable period for a primary infectious individual cannot result in a false negative test. This option makes a model diverge from a branching process.\n");
#ifdef DUAL_PINF
  printf("\t--pri_first_category_only\t\tA primary infectious individual can only be part of the first category (disables pri_second_category_only).\n");
  printf("\t--pri_second_category_only\t\tA primary infectious individual can only be part of the second category (disables pri_first_category_only).\n");
#endif

  printf("\t--time_rel_pri_created\t\tRecorded event time is relative to the time the primary individuals are generated.\n");
  printf("\t--time_rel_pri_infectious\tRecorded event time is relative to the time the generated primary individuals become infectious.\n");
  printf("\t--time_rel_pri_end_comm\t\tRecorded event time is relative to the end of the communicable period for the generated primary individuals.\n");
  printf("\t--time_rel_pri_flat_comm\tThe primary individuals are assumed to enter the simulation at a random time within their communicable period. There is thus no latent period for the these individuals and the duration of their communicable period within the simulation is truncated with a uniform probability. This option makes a model diverge from a branching process.\n");
  printf("\t--time_rel_pri_test_results\tRecorded event time is relative to the time the generated primary individuals receive test results.\n");
  printf("\t--time_rel_first_pos_test_results\tRecorded event time is relative to the time of the first positive test result. This operation is performed in post-processing.\n");
  printf("\t--include_all_paths\t\tIndicate that observable and non-observable paths should be included in the simulation results (default).\n");
  printf("\t--observable_paths_only\t\tIndicate that only observable paths should be included in the simulation results.\n");
  printf("\t--non-observable_paths_only\tIndicate that only non-observable paths should be included in the simulation results.\n");
  printf("\t--tmax VALUE\t\t\tMaximum simulation time used to instantiate new infectious individuals (default value of INFINITY).\n");
  printf("\t--nstart VALUE\t\t\tInitial number of individuals (default value of 1).\n");
  printf("\t--pinfpri VALUE\t\t\tProbability that an initial individual be infectious (default value of 1).\n");
  printf("\t--lmax VALUE\t\t\tMaximum number of layers (generations) for the simulation (value of 1 signifies only primary individuals, default value of UINT32_MAX).\n");
  printf("\t--nbinsperunit VALUE\t\tNumber of timeline bins per unit of time.\n");
  printf("\t--nimax VALUE\t\t\tMaximum number of infectious individuals for a given time integer interval (default value of UINT32_MAX). This option makes a model diverge from a branching process, but does not affect the expected effective reproduction number value.\n");
  printf("\t--npostestmax VALUE\t\tMaximum number of positive test results during an interval of duration npostestmaxunits that starts when the test results are received. (default value of UINT32_MAX). This option makes a model diverge from a branching process, but does not affect the expected effective reproduction number value.\n");
  printf("\t--npostestmaxnunits VALUE\tInterval duration for the maximum number of positive test results (default value of 1).\n");
  printf("\t--tlout FILENAME\t\tOutput timeline information for each simulated path into the provided file in the binary format as described below.\n");
  printf("\t--tloutbufsize VALUE\t\tPer-thread memory buffer size (in MB) used to accumulate data for timeline output before writing them to disk (default value of 10 MB).\n");
#ifdef CT_OUTPUT
  printf("\t--ctout FILENAME\t\tOutput contact tracing information for each simulated path into the provided file.\n");
  printf("\t--ctoutbufsize VALUE\t\tPer-thread memory buffer size (in MB) used to accumulate data for contact tracing output before writing them to disk (default value of 10 MB).\n");
#endif
  printf("\t--ninfhist\t\t\tCompute a histogram of the number of infected individuals for each infectious individual.\n");
  printf("\t--npaths VALUE\t\t\tNumber of generated simulation paths (default value of 10000).\n");
  printf("\t--nthreads VALUE\t\tNumber of threads used to perform the simulation (default value of 1).\n");
  printf("\t--nsetsperthread VALUE\t\tNumber of path sets used for each thread (default value of 100 when nthreads>1, and of 1 otherwise). Using a value of 1 guarantees the same stream of random numbers from one run to another, while using a larger value increases performance by assigning sets to available processing resources. In either case, the RNG stream algorithm is used to guarantee non-overlapping seed streams between threads.\n");
  printf("\t--stream VALUE\t\t\tSelect an RNG stream. Use to set the initial seed of the random number generator (default value of 0).\n");
  printf("\t--help\t\t\t\tPrint this usage information and exit.\n");
  printf("\n\tEach option can be used as shown above from the command line. Dash(es) for option names are optional. For configuration files, '=', ':' or spaces as defined by isspace() can be used to separate option names from arguments. Characters following '#' on one line are considered to be comments.\n");
  printf("\tOptions can be used multiple times and configuration files can be read from configuration files.\n"); 

  printf("\n\nBINARY TIMELINE OUTPUT FILE:\n");
  printf("\n\tAll fields are stored in little endian.\n");
  printf("\n\tFile header:\n");
  printf("\t\t-Unsigned 32 bit value: tmax, the number of time bins starting from t=0.\n");
  printf("\t\t-8 bit field:\n");
  printf("\t\t\tBits 0 to 2: A value from the three lower significant bits is used to indicate the model of the time origin. A value of 1 for primary individual creation time (time_rel_pri_created), 2 for primary individual entering the simulation at a random time during his communicable period (time_rel_pri_flat_comm), 3 for time primary individual becomes infectious (time_rel_pri_infectious), 4 for end of communicable period for primary individual (time_rel_pri_end_comm), 5 for test results for primary individual (time_rel_pri_test_results).\n");
  printf("\t\t\tBit 3: Indicates if a timeline is included for positive test results.\n");
#ifdef SEC_INF_TIMELINES
  printf("\t\t\tBit 4: Indicates that second series of timelines is included for the second category of infection.\n");
#endif

  printf("\n\tSimulation path records:\n");
  printf("\t\t-Unsigned 32 bit value: The number of written successive time bins.\n");
  printf("\t\t-Unsigned 32 bit value: Field is written only if the time mode is not the primary individual creation time. Value is the number of time bins before t=0.\n");
  printf("\t\t-Signed 32 bit value: Period (defined as floor(t)) where the path maxes out an nimax or npostestmax limit, if any. Otherwise, a value of INT32_MAX.\n");
  printf("\t\t-Signed 32 bit value: Period (defined as floor(t)) where the path goes extinct, if it does. Otherwise, a value of INT32_MAX. For paths without any initial infection, it is set to -INT32_MAX\n");
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written: Number of active infections.\n");
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written: Number of new infections.\n");
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written (written only if indicated in the file header): Number of new positive test results.\n");
#ifdef SEC_INF_TIMELINES
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written: Number of active infections for the second category of infection.\n");
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written: Number of new infections for the second category of infection.\n");
  printf("\t\t-Unsigned 32 bit value, for each time bin, chronologically written (written only if indicated in the file header): Number of new positive test results for the second category of infection.\n");
#endif
}
