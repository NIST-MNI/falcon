/*  minc_add_history
 *
 *  Copyright 2014  Kunio Nakamura
 *
 *
 * compile example
 * gcc -o ../bin/minc_add_history minc_add_history.c -lm -lz -fopenmp -ltiff -lgsl -lgslcblas -L/Kunio/minc/lib -lminc2 -lhdf5 -lnetcdf -ldl /Kunio/minc/lib/*.so  -fopenmp -O3 -DNDEBUG -I. -I/usr/include/ -Wall -I/Kunio/minc/include
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif //HAVE_CONFIG_H

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <volume_io.h>
#include <time.h>

#include "basic.h"
#include "mincio.h"
#include <time_stamp.h>


int main(int argc, char  *argv[] ) {
  char
  path[256]="",
            old_history[65536],
            new_history[65536];
  char fcname[64]="minc_add_history";
  int i;
  int verbose=0;
  int mi_status;
  mihandle_t v;
  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  fprintf(stdout,"[%s] argc = %i\n",fcname,argc);
  fprintf(stdout,"[%s] argv = ",fcname);
  for(i=1; i<argc; i++) {
    fprintf(stdout,"%s ",argv[i]);
  }
  fprintf(stdout,"\n");

  fprintf(stdout,"[%s] reading input %s\n",fcname,argv[1]);
  if((mi_status = miopen_volume(argv[1], MI2_OPEN_RDWR, &v)) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miopen_volume\n",fcname);
    exit(1);
  }

  for(i=0; i<65536; i++) new_history[i]=old_history[i]=0;

  if(verbose>0) fprintf(stdout,"[%s]   miget_attr_values\n",fcname);
  if((mi_status = miget_attr_values ( v, MI_TYPE_STRING, path, "history", 65536, old_history)) == MI_ERROR) {
    fprintf(stdout,"[%s] no old history\n",fcname);
  } else {
    fprintf(stdout,"[%s] OLD HISTORY %i\n%s-------------------\n",fcname,(int)strlen(old_history),old_history);
  }

  strcat(old_history,tmstr);
  strcat(old_history,">>>");
  for(i=2; i<argc; i++) {
    strcat(old_history," ");
    strcat(old_history,argv[i]);
    if(verbose>0) fprintf(stdout,"\t%s\n",argv[i]);
  }
  sprintf(new_history,"%s\n",old_history);
  fprintf(stdout,"[%s] NEW HISTORY %i\n%s-------------------\n",fcname,(int)strlen(new_history),new_history);

  if(verbose>0) fprintf(stdout,"[%s]   miadd_history_attr\n",fcname);
  if((mi_status = miadd_history_attr(v,(int)strlen(new_history),new_history)) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miadd_history_attr\n",fcname);
    exit(1);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miflush_volume\n",fcname);
  if((mi_status = miflush_volume( v )) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miflush_volume\n",fcname);
    exit(1);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miclose_volume\n",fcname);
  if((mi_status = miclose_volume ( v )) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miclose_volume\n",fcname);
    exit(1);
  }

  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/