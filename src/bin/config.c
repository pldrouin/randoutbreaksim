#include "config.h"

int config(sim_pars* pars, uint32_t* npaths, int* oout, int* eout, const int nargs, const char* args[])
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

      } else if(!argsdiffer(pbuf, "q")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->q);

      } else if(!argsdiffer(pbuf, "mbar")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->mbar);

      } else if(!argsdiffer(pbuf, "kappaq")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->kappaq);

      } else if(!argsdiffer(pbuf, "tmax")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%lf",&pars->tmax);

      } else if(!argsdiffer(pbuf, "nstart")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,&pars->nstart);

      } else if(!argsdiffer(pbuf, "npaths")) {
	safegetnextparam(fptra,&fptri,true,nargs,args,&parc,pbuf);
	sscanf(pbuf,"%"PRIu32,npaths);

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
  printf("\t--tbar VALUE\t\tbranchsim's tbar parameter (required)\n");
  printf("\t--p VALUE\t\tbranchsim's tbar parameter (required)\n");
  printf("\t--lambda VALUE\t\tbranchsim's tbar parameter (required)\n");
  printf("\t--kappa VALUE\t\tbranchsim's tbar parameter (required)\n");
  printf("\t--q VALUE\t\tbranchsim's tbar parameter (default value of 0)\n");
  printf("\t--mbar VALUE\t\tbranchsim's tbar parameter (required if q>0)\n");
  printf("\t--kappaq VALUE\t\tbranchsim's tbar parameter (required if q>0)\n");
  printf("\t--tmax VALUE\t\tbranchsim's tbar parameter (default value of INFINITY)\n");
  printf("\t--nstart VALUE\t\tbranchsim's tbar parameter (default value of 1)\n");
  printf("\t--npaths VALUE\t\tbranchsim's tbar parameter (default value of 10000)\n");
  printf("\t--help\t\t\tPrint this usage information and exit\n");
  printf("\nEach option can be used as shown above from the command line. Dash(es) for option names are optional. For configuration files, '=', ':' or spaces as defined by isspace() can be used to separate option names from arguments. Characters following '#' on one line are considered to be comments.\nOptions can be used multiple times and configuration files can be read from configuration files.\n"); 
}
