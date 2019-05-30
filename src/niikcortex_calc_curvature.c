/* FILENAME: niikcortex_calc_thickness.c
 * DESCRIPTION: calculate cortical thickness using simple algorithm
 * AUTHOR:       Kunio Nakamura 2013
 *               Vladimir S. FONOV 2017
 */

#include "falcon.h"
#include "falcon_surfaces.h"

#include <unistd.h>
#include <getopt.h>

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"Surface curvature calculation\n");
  fprintf(stdout,"  usage: [options] <input.off> <out1.txt/ply> [out2.txt/ply]\n\n");
  fprintf(stdout,"Optional arguments:\n"
          "\t--old - use old-style curvature estimation\n"
          "\t--min - in case of one output output min abs curvature\n"
          "\t--max - in case of one output output max abs curvature\n"
          "\t--gauss - output gaussian curvature\n"
          "\t--abs   - output absolute curvature\n"
          "\t--mean  - output mean curvature\n"
          "\t--ply   - output .ply file with measurement\n"
          );

}


int main(int argc,char *argv[],char *envp[]) {
  kobj *obj;
  const char *fcname="niikcortex_calc_curvature.c";
  int clobber=0;
  int verbose=0;
  int output_sph=0;
  int output_ply=0;
  int old_style=0;
  int iter=3;
  int nval=10;
  int i;

  const char *header1;
  const char *header2;

  int output_min=0;
  int output_max=0;
  int output_mean=0;
  int output_gauss=0;
  int output_abs=0;

  const char *in_off,*out_curvature,*out_curvature2=NULL;
  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"sph",              no_argument, &output_sph, 1},
    {"old",              no_argument, &old_style, 1},
    {"min",              no_argument, &output_min, 1},
    {"max",              no_argument, &output_max, 1},
    {"avg",              no_argument, &output_mean, 1},
    {"mean",             no_argument, &output_mean, 1},
    {"gauss",            no_argument, &output_gauss, 1},
    {"abs",              no_argument, &output_abs, 1},
    {"ply",              no_argument, &output_ply, 1},
    {"help",             no_argument, 0, 'h'},
    {"version",          no_argument, 0, 'v'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "uhv", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
      return 1;
    case '?':
    case 'u':
    case 'h':
    default:
      usage();
      return 1;
    }
  }

  if((argc - optind)<2) {
    usage();
    return 1;
  }
  in_off=argv[optind];
  out_curvature=argv[optind+1];

  if((argc-optind)>2) {
    out_curvature2=argv[optind+2];
  }

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading surface       %s\n",fcname, in_off);
  NIIK_RETURN(((obj=off_kobj_read_offply(in_off))==NULL),"reading surface",1);

  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);

  if(iter<0) iter=3;

  if(old_style) {
    niikvec *vec;

    vec=niikvec_init(obj->nvert);
    if(!niik_off_curvature_map_update(obj,vec,iter)) {
      fprintf(stderr,"ERROR: niik_off_curvature_map_update\n");
      exit(1);
    }
    fprintf(stdout,"[%s] writing output text        %s\n",fcname,out_curvature);
    NIIK_RETURN((!niikvec_write(out_curvature,vec)),"niikvec_write",1);
    vec=niikvec_free(vec);
    header1="curvature";
  } else {

    double *cornerareas[3]= {
      (double*)calloc(obj->nface,sizeof(double)),
      (double*)calloc(obj->nface,sizeof(double)),
      (double*)calloc(obj->nface,sizeof(double))
    };

    double *pointareas=(double*)calloc(obj->nvert,sizeof(double));
    niikpt *pdir1=(niikpt*)calloc(obj->nvert, sizeof(niikpt));
    niikpt *pdir2=(niikpt*)calloc(obj->nvert, sizeof(niikpt));
    double *curv1=(double*)calloc(obj->nvert, sizeof(double));
    double *curv2=(double*)calloc(obj->nvert, sizeof(double));
    double *curv12=(double*)calloc(obj->nvert,sizeof(double));

    if(!niik_off_pointareas(obj, pointareas, cornerareas)) {
      fprintf(stderr,"ERROR: niik_off_pointareas\n");
      exit(1);
    }

    if(!niik_off_curvature_map_rusinkiewicz(obj, pointareas, cornerareas, pdir1, pdir2, curv1, curv2, curv12)) {
      fprintf(stderr,"ERROR: niik_off_curvature_map_rusinkiewicz\n");
      exit(1);
    }

    if( !out_curvature2 ) {
      int i;
      for(i=0; i<obj->nvert; i++) {
        if( output_min ) {
          header1="min_curvature";
          if(fabs(curv2[i])<fabs(curv1[i])) curv1[i]=curv2[i];
        } else if( output_max ) {
          header1="max_curvature";
          if(fabs(curv2[i])>fabs(curv1[i])) curv1[i]=curv2[i];
        } else if( output_mean ) {
          header1="mean_curvature";
          curv1[i]=(curv1[i]+curv2[i])/2.0;
        } else if( output_gauss ) {
          header1="gaussian_curvature";
          curv1[i]=curv1[i]*curv2[i];
        } else {
          header1="max_curvature";
          if(fabs(curv2[i])>fabs(curv1[i])) curv1[i]=curv2[i]; //min
        }
        if( output_abs )
          curv1[i]=fabs(curv1[i]);
      }
    } else {
      header1="curvature1";
      header2="curvature2";
      
      NIIK_RETURN((!niik_write_double_vector_ex(out_curvature2, curv2, obj->nvert,header2)),"niik_write_double_vector",1);
    }

    if(output_ply)
    {
      const char * meas[]={header1};
      double *val[]={curv1};
      NIIK_RETURN((!off_kobj_write_ply_ex(out_curvature, obj ,0, 1, 1,0,1,meas,val)),"off_kobj_write_ply_ex",1);
    }
    else 
      NIIK_RETURN((!niik_write_double_vector_ex(out_curvature, curv1, obj->nvert,header1)),"niik_write_double_vector",1);

    free(curv12);
    free(curv1);
    free(curv2);
    free(pdir1);
    free(pdir2);
    free(pointareas);
    free(cornerareas[0]);
    free(cornerareas[1]);
    free(cornerareas[2]);
  }

  obj=off_kobj_free(obj);

  niik_fc_display(fcname,0);
  exit(0);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/