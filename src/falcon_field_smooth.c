/* FILENAME: falcon_meas_smooth.c
 * DESCRIPTION: smooth measurements on the surface
 * AUTHOR:       Kunio Nakamura
 */

#include "falcon.h"
#include "falcon_cortex.h"

#include <unistd.h>
#include <getopt.h>

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

void usage() {
  fprintf(stdout,"Itarative Gaussian smoothing of signals on the surface, will try to smooth all columns in the input  \n");
  fprintf(stdout,"  usage: [options] <surface.ply> <in.csv/in.csv.gz/in.txt> <out.csv/out.csv.gz/out.txt> \n\n");
  fprintf(stdout,"Options:\n");
  fprintf(stdout,"\t--fwhm/--smooth <fwhm> - apply smoothing kernel (default: 1.0) - will try to apply it iteratively\n");
  /*fprintf(stdout,"\t--column  <name> - use column name (default 1st)\n");
  fprintf(stdout,"\t--column_n  <n> - use column no (default 1st)\n");*/
  fprintf(stdout,"\t--clobber     - overwrite output\n");
  fprintf(stdout,"\t--help        - this message\n");
}


int main(int argc,char *argv[],char *envp[])
{
  kobj *ctx[1];
  double *vectors[3];
  double smooth=0.0;
  const char *fcname="falcon_meas_smooth.c";
  int clobber=0;
  int output_sph=0;
  int output_ply=0;
  int output_white=0;
  int col_id=-1;

  niiktable *meas_in;

  const char *in_ics,*in_csv,*out_csv,*col_name=NULL;

  double it_smooth,it_sigma;
  double mean_elen;

  int iterations;
  int i,c;


  struct option long_options[] = {
    {"clobber",          no_argument, &clobber, 1},
    {"help",             no_argument, 0, 'h'},
    {"version",          no_argument, 0, 'v'},
    {"smooth",     required_argument, 0, 'f'},
    {"fwhm",       required_argument, 0, 'f'},
/*    {"column",     required_argument, 0, 'c'},
    {"column_n",   required_argument, 0, 'n'},*/
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long (argc, argv, "uhvf:c:n:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'n':
        col_id=atoi(optarg);
        col_name=NULL;
        break;
    case 'c':
        col_id=-1;
        col_name=optarg;
        break;
    case 'v':
      fprintf(stdout,"[%s] version %i.%i.%i\n",fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
      return 1;
    case 'f':
      smooth=atof(optarg);
      break;
    case '?':
    case 'u':
    case 'h':
    default:
      usage();
      return 1;
    }
  }

  if((argc - optind)<3) {
    usage();
    return 1;
  }
  in_ics=argv[optind];
  in_csv=argv[optind+1];
  out_csv=argv[optind+2];

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading surface       %s\n",fcname,in_ics);
  NIIK_RETURN(((ctx[0]=off_kobj_read_offply(in_ics))==NULL),"reading white surface",1);


  NIIK_EXIT(((meas_in=niiktable_read(in_csv))==NULL),fcname,"niiktable_read",1);
  if(meas_in->col[0]->num != ctx[0]->nvert) {
    fprintf(stderr,"[%s] Inconsistent number of measurement: %i , expected %i\n",fcname, meas_in->col[0]->num, ctx[0]->nvert);
    exit(1);
  }

  mean_elen=off_get_kobj_mean_edge_length(ctx[0]);

  iterations=ceil((smooth*smooth)/(mean_elen*mean_elen*4));
  /*apply small kernel several times, since each iteration only works with the local neighbourhood*/
  /*HACK: mean_elen*2  - is an arbitrary factor */

  if(iterations<1) iterations=1;
  it_smooth=sqrt(smooth*smooth/iterations);
  it_sigma=it_smooth/2.355;

  fprintf(stdout,"[%s] Smoothing: %f x %d\n",fcname, it_smooth, iterations);

  /*TODO: identify columns that shouldn't be smoothed!*/
  for(c=0;c<meas_in->ncol;c++) {
        NIIK_EXIT((!off_surface_gauss_smooth_using_vert(ctx[0], meas_in->col[c]->v, it_sigma, iterations)),fcname,"off_surface_gauss_smooth_using_vert",1);
  }

  fprintf(stdout,"[%s] writing output file        %s\n",fcname, out_csv);

  NIIK_EXIT((niiktable_write(out_csv, meas_in)==0), fcname, "niiktable_write",1);

  niiktable_free(meas_in);

  ctx[0]=off_kobj_free(ctx[0]);

  niik_fc_display(fcname,0);
  return 0;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/