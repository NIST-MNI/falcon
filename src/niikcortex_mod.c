/* FILENAME:     niikcortex_mod.c
 * DESCRIPTION:  Kunio's nifti1 for cortical deformation
 * AUTHOR:       Kunio Nakamura
 * DATE:         May 8, 2014
 * AUTHOR:       Vladimir FONOV
 * DATE:         Nov 22,2018
 */

#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include "falcon.h"
#include "falcon_cortex.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (1)


#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


void prog_version() {
  fprintf(stdout,"  niik_niikcortex_mod history\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  0.0.1  Sep 2, 2016, Vladimir S. FONOV - moved from niikmath\n");
  fprintf(stdout,"\n");
}

void prog_usage() {
  fprintf(stdout,"niikcortex_mod\n");
  fprintf(stdout,"  usage: [options] <input.off> <output.off> \n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  niikcortex_initocs_help:\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help                   : show this usage\n");
  fprintf(stdout,"  --iter <iter>                     : maximum number of outer iterations [default=20]\n");
  fprintf(stdout,"  --iter2 <iter>                    : maximum number of inner iterations [default=5]\n");
  fprintf(stdout,"  --step  <f>                       : Step size for outer iteration: positive - baloon, negative - shrink [default=0.1] \n");
  fprintf(stdout,"  --version                         : show version info\n");
}

int main(int argc,char *argv[],char *envp[]) {
  int maxiter=5;
  int maxiter2=20;
  double step=0.1;
  int clobber=0;

  kobj *obj=NULL;

  const char *in_ctx=NULL;
  const char *out_ocs=NULL;
  const char *fcname="niikcortex_mod";
  int i,n;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"version", no_argument, 0, 'v'},
    {"help",    no_argument, 0, 'h'},
    {"usage",   no_argument, 0, 'U'},
    {"iter", required_argument,  0,      'T'},
    {"iter2",required_argument,  0,      't'},
    {"step", required_argument, 0,       's'},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hUvT:t:s:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      prog_version();
      return 0;
    case 't':
      maxiter = atoi(optarg);
      break;
    case 'T':
      maxiter2 = atoi(optarg);
      break;
    case 's':
      step = atof(optarg);
      break;
    case 'h':
    case 'U':
    case '?':
    default:
      prog_usage ();
      return 1;
    }
  }

  if((argc - optind)<2) {
    prog_usage();
    return 1;
  }
  
#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  in_ctx  = argv[optind];
  out_ocs = argv[optind+1];

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);

  if((obj=off_kobj_read_offply(in_ctx))==NULL) {
    fprintf(stderr,"[%s] ERROR: off_kobj_read_off\n",fcname);
    exit(1);
  }

  fprintf(stdout,"[%s] test intersections\n",fcname);
  if((n=off_count_self_intersection_add_color(NULL,obj,1))<0) {
    fprintf(stderr,"[%s] ERROR: off_count_self_intersection_add_color(bb,obj,1)\n",fcname);
    exit(1);
  }

  if(n>0) {
    fprintf(stdout,"[%s] ERROR: ics has surface %i intersection(s)\n",fcname,n);
    exit(1);
  }

  for(i=0;i<maxiter2;i++) {
    fprintf(stdout,"[%s] iteration:%d\n",fcname,i);
    if(!off_kobj_balloon(obj,maxiter,step/maxiter,0)) {
      fprintf(stderr,"[%s] ERROR: off_kobj_balloon\n",fcname);
      exit(1);
    }
  }

  fprintf(stdout,"[%s] test intersections\n",fcname);
  if((n=off_count_self_intersection_add_color(NULL,obj,1))<0) {
    fprintf(stderr,"[%s] ERROR: off_count_self_intersection_add_color(bb,obj,1)\n",fcname);
  }

  if(n>0)
  {
    fprintf(stdout,"[%s] ERROR: ics has surface %i intersection(s)\n",fcname,n);
  } else if(!off_kobj_add_one_color(obj,1,0,0)) {
    fprintf(stderr,"[%s] ERROR: off_kobj_add_one_color(obj,1,0,0)\n",fcname);
    exit(1);
  }

  /*append metadata*/
  NIIK_RET1((!off_kobj_add_comment(obj,timestamp)),fcname,"off_kobj_add_comment");

  fprintf(stdout,"[%s] writing output  %s\n",fcname,out_ocs);
  if(!off_kobj_write_offply(out_ocs,obj,0)) {
    fprintf(stderr,"[%s] ERROR: off_kobj_write_offply(out_ocs,obj,0)\n",fcname);
    exit(1);
  }

  /*free memory*/
  obj=off_kobj_free(obj);
  niik_fc_display(fcname,0);
  free(timestamp);
  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/