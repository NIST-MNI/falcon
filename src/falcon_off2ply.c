/* FILENAME:     off2asc.c
 * DESCRIPTION:  Conversion from off to text asc file
 * AUTHOR:       Vladimir Fonov
 *
 *
 */

#include "falcon.h"
#include "falcon_surfaces.h"
#include "rply.h"

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

/*needed for ply file format, probably should be used everywhere else too*/
#include <locale.h>

static char *prog_version[] = {
  NIIK_VERSION "\n"
  "  program history\n"
  "  0.0.0   : August 16, 2016, Vladimir FONOV <vladimir.fonov@gmail.com>\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  off2ply:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help             : show this usage\n"
  "  --version                   : show version info\n"
};

void usage() {
  fprintf(stdout,"obj2ply\n");
  fprintf(stdout,"  usage: [options] <in.off> [measurements] <out.ply>\n");
  fprintf(stdout,"  options:\n");
  fprintf(stdout,"\t--sph     - output spherical coordinates\n");
  fprintf(stdout,"\t--normal  - output vertex normals\n");
  fprintf(stdout,"\t--edge    - output edges\n");
  fprintf(stdout,"\t--clobber - clobber output\n");
}

int main(int argc,char *argv[],char *envp[]) {
  kobj *obj=NULL;
  int output_sph=0;
  int output_normal=0;
  int output_edge=0;
  int clobber=0;
  const char  *fcname="off2ply";
  double *meas=NULL;
  int meas_num;

  const char *in_off=NULL;
  const char *in_meas=NULL;
  const char *out_fname=NULL;
  const char *old_locale;

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"help",no_argument, 0, 'h'},
    {"version",no_argument, 0, 'v'},
    {"sph", no_argument, &output_sph, 1},
    {"normal", no_argument, &output_normal, 1},
    {"edge", no_argument, &output_edge, 1},
    {0, 0, 0, 0}
  };
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "hv", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      fprintf(stdout,"%s\n",*prog_version);
      return 0;
    case '?':
    case 'h':
    default:
      usage ();
      return 1;
    }
  }

  if((argc - optind)<2) {
    usage();
    return 1;
  }

  in_off =argv[optind];
  if((argc - optind)>2) {
    in_meas=argv[optind+1];
    out_fname=argv[optind+2];
  } else {
    out_fname=argv[optind+1];
  }

  if (!clobber && !access (out_fname, F_OK)) {
    fprintf(stderr,"%s Exists!\n", out_fname);
    return 1;
  }

  niik_fc_display(fcname,1);

  old_locale = setlocale(LC_NUMERIC, NULL);
  setlocale(LC_NUMERIC, "C");
  NIIK_EXIT(((obj=off_kobj_read_offply(in_off))==NULL),fcname,"niik_kobj_read_off",1);
  if(in_meas) {
    meas=niik_read_double_vector(in_meas,&meas_num);

    if(meas_num!=obj->nvert) {
      fprintf(stderr,"[%s] Inconsistent number of measurement: %i , expected %i\n",fcname,meas_num,obj->nvert);
      return 1;
    }
  }
  fprintf(stdout,"[%s] #points = %i\n",fcname,obj->nvert);
  fprintf(stdout,"[%s] #triangles = %i\n",fcname,obj->nface);


  if(output_normal) {
    /*update all normals?*/
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
  }

  NIIK_EXIT((!niik_off_write_ply(out_fname, obj, meas, output_sph, output_normal, output_edge)),fcname,"niik_off_write_ply",1);

  setlocale(LC_NUMERIC, old_locale);
  obj=off_kobj_free(obj);
  if(meas!=NULL)
    free(meas);

  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
} /* off2mniasc */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
