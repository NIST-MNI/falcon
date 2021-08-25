/* FILENAME:     test_off.c
 * DESCRIPTION:  test basic object routines
 * AUTHOR:       Vladimir S. FONOV
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

void show_usage (const char *name) {
  fprintf(stdout,"Usage: %s <input.ply> <radius> [--tolerance <f>, default 0.001 --output <output.ply> write error map]\n", name);
}

int check_sphere(kobj *obj, double radius, double *rms_error, double *max_error, int *err_idx, double **err)
{
  kvert *v;
  int i ;

  *max_error=0.0;
  *err_idx=0;
  *rms_error=0.0;

  if(err)
    *err=calloc(obj->nvert,sizeof(double));

  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    double r = sqrt(v->v.x*v->v.x+v->v.y*v->v.y+v->v.z*v->v.z);
    double d = r-radius;
    if(err) (*err)[i] = d;
    *rms_error+=d*d;
    d=fabs(d);
    if(d>(*max_error) ) 
    {
      *max_error=d;
      *err_idx=i+1;
    }
  }
  *rms_error=sqrt(*rms_error/i);
  return 0;
}

int main(int argc, char **argv) {
  const char *fcname="test_sphere";
  int clobber=0;
  int verbose=0;
  int c;
  int i;
  int r=0;
  const char *in_off=NULL;
  double radius=0.0;
  double tolerance=1e-3;
  const char *out_err_off=NULL;

  int n;
  kobj *obj=NULL;

  struct option long_options[] = {
    {"verbose",  no_argument, &verbose, 1},
    {"tolerance",required_argument, 0, 't'},
    {"output",required_argument, 0, 'o'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "t:o:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 't':
      tolerance=atof(optarg);
      break;
    case 'o':
      out_err_off=optarg;
      break;
    case '?':
    default:
      show_usage (argv[0]);
      return 1;
    }
  }

  if((argc - optind)<2) {
    show_usage(argv[0]);
    return 1;
  }

  in_off =argv[optind];
  radius =atof(argv[optind+1]);

  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_EXIT( ((obj=off_kobj_read_offply(in_off))==NULL), fcname,"off_kobj_read_offply", 10 );

  double *err;
  double max_err,rms_err;
  int err_idx;
  
  check_sphere(obj, radius, &rms_err, &max_err,&err_idx,&err);

  fprintf(stderr,"RMS error:%f Max error:%f at %d\n",rms_err,max_err,err_idx);

  if(rms_err>tolerance)
    r=1;

  if(out_err_off)
  {
    const char *meas[]={"error"};
    NIIK_EXIT( ((off_kobj_write_ply_ex(out_err_off,obj, 0,1,1,0,1,meas, (const double**) &err))==0), fcname,"off_kobj_write_ply_ex", 20 );
  }

  off_kobj_free(obj);
  free(err);
  if(verbose>0) niik_fc_display(fcname,0);
  return r;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
