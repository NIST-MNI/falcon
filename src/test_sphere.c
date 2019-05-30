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
  fprintf(stdout,"Usage: %s <input.ply> <radius> [--tolerance <f>]\n", name);
}

int check_sphere(kobj *obj,double radius,double tolerance,int verbose)
{
  kvert *v;
  int i ;
  int ret=0;

  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    double r = sqrt(v->v.x*v->v.x+v->v.y*v->v.y+v->v.z*v->v.z);
    if(fabs(r-radius)>tolerance) {
      ret=1;
      if(verbose)
        fprintf(stdout,"Non-spherical coordinate found: %d (%f,%f,%f) r:%f expected %f\n",i,v->v.x,v->v.y,v->v.z,r,radius);
      break;
    }
  }
  return ret;
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

  int n;
  kobj *obj=NULL;

  struct option long_options[] = {
    {"verbose",  no_argument, &verbose, 1},
    {"tolerance",required_argument, 0, 't'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "t:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 't':
      tolerance=atof(optarg);
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

  r = check_sphere(obj,radius,tolerance,verbose);

  off_kobj_free(obj);
  if(verbose>0) niik_fc_display(fcname,0);
  return r;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
