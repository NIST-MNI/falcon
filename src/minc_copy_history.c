/*  minc_add_history
 *
 *  Copyright 2014  Kunio Nakamura
 *
 *
 * compile example
 * gcc -o ../bin/minc_copy_history minc_copy_history.c -lm -lz -fopenmp -ltiff -lgsl -lgslcblas -L/Kunio/minc/lib -lminc2 -lhdf5 -lnetcdf -ldl /Kunio/minc/lib/*.so  -fopenmp -O3 -DNDEBUG -I. -I/usr/include/ -Wall -I/Kunio/minc/include
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
            history[65536];
  char fcname[64]="minc_copy_history";
  int verbose=0;
  int mi_status;
  mihandle_t s,t;
  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  if(argc!=3) {
    fprintf(stderr,"[%s] usage: <from.mnc> <to.mnc>\n",fcname);
    exit(1);
  }
  fprintf(stdout,"[%s] reading input %s\n",fcname,argv[1]);
  if((mi_status = miopen_volume(argv[1], MI2_OPEN_READ, &s)) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miopen_volume\n",fcname);
    exit(1);
  }
  fprintf(stdout,"[%s] reading input %s\n",fcname,argv[2]);
  if((mi_status = miopen_volume(argv[2], MI2_OPEN_RDWR, &t)) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miopen_volume\n",fcname);
    exit(1);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miget_attr_values\n",fcname);
  if((mi_status = miget_attr_values (s, MI_TYPE_STRING, path, "history", 65536, history)) == MI_ERROR) {
    fprintf(stdout,"[%s] no old history\n",fcname);
    exit(1);
  } else {
    fprintf(stdout,"[%s] HISTORY %i\n%s-------------------\n",fcname,(int)strlen(history),history);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miadd_history_attr\n",fcname);
  if((mi_status = miadd_history_attr(t,(int)strlen(history),history)) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miadd_history_attr\n",fcname);
    exit(1);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miflush_volume\n",fcname);
  if((mi_status = miflush_volume( t )) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miflush_volume\n",fcname);
    exit(1);
  }

  if(verbose>0) fprintf(stdout,"[%s]   miclose_volume\n",fcname);
  if((mi_status = miclose_volume ( t )) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miclose_volume\n",fcname);
    exit(1);
  }
  if((mi_status = miclose_volume ( s )) == MI_ERROR) {
    fprintf(stderr,"[%s] ERROR: miclose_volume\n",fcname);
    exit(1);
  }

  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/