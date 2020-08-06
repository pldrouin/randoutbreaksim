#include "args.h"

int getnextparam(FILE **fptra, int *fptri, const bool isarg, const int nargs, char const* const args[], int *parc, char *param)
{
  int i=0;

  while(*fptri>-1 || nargs>*parc) {

    if(*fptri==-1) {
      if(!isarg) while(args[*parc][i] == '-') i++;
      strcpy(param,args[*parc]+i);
      (*parc)++;
      //printf("param is '%s'\n",param);
      return strlen(param);

    } else {
      signed char c;
      bool brs=false, brd=false, com=false;

      while(isspace((c=fgetc(fptra[*fptri]))) || ((c == '=' || c == ':') && isarg) || c == '#' || (com && c != 0 && c != EOF)) {if(com && c == '\n') com=false; if(c == '#') com=true;}
      if(!isarg && c == '-') while((c=fgetc(fptra[*fptri])) == '-') {}

      if(c != 0 && c != EOF) {

	if(c=='\'') brs=!brs;
	else if(c=='\"') brd=!brd;
	else param[i++]=c;

	while(((!isspace((c=fgetc(fptra[*fptri]))) && ((c != '=' && c != ':') || isarg) && c != '#') || brs || brd) && c != 0 && c != EOF) {if(c=='\'') brs=!brs; else if(c=='\"') brd=!brd; else param[i++]=c;}
	param[i++]=0;

	if(c == 0 || c == EOF) {
	  if(fptra[*fptri]!=stdin) fclose(fptra[*fptri]);
	  (*fptri)--;
	}
	//printf("param is '%s'\n",param);
	return i;

      } else {
	if(fptra[*fptri]!=stdin) fclose(fptra[*fptri]);
	(*fptri)--;
      }
    }
  }
  return -1;
}

void safegetnextparam(FILE **fptra, int *fptri, const bool isarg, const int nargs, char const* const args[], int *parc, char *param)
{
  if(getnextparam(fptra,fptri,isarg,nargs,args,parc,param)<0) {
    fprintf(stderr,"Error: Missing parameter\n");
    exit(1);
  }
}
