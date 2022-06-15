/* FILENAME:     niikcortex_initocs.c
 * DESCRIPTION:  Kunio's nifti1 for cortical deformation
 * AUTHOR:       Kunio Nakamura
 * DATE:         May 8, 2014
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
  fprintf(stdout,"  niik_niikcortex_initocs_simple history\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  0.0.1  Jun 14, 2022, Vladimir S. FONOV - branched from niikcortex_initocs.c\n");
  fprintf(stdout,"\n");
}

void prog_usage() {
  fprintf(stdout,"niikcortex_initocs_simple\n");
  fprintf(stdout,"  usage: [options] <mask.mnc> <init_ics.ply> <out_init_ocs.ply> \n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  niikcortex_initocs_simple help:\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help            show this usage\n");
  fprintf(stdout,"  --vlist <step1>[,<step2>...] expansion step sizes (mm) [default=0.1,0.05,0.01]\n");
  fprintf(stdout,"  --max-thick <MT>           maximum allowed thickness (mm) [default=3]\n");
  fprintf(stdout,"  --smooth-iter <iter>       smooth iteration [default=3]\n");
  fprintf(stdout,"  --smooth-factor <f>        smooth factor [default=0.75]\n");
  fprintf(stdout,"  --iter <iter>              maximum number of iterations [default=20]\n");
  fprintf(stdout,"  --delta <init_cth>         initial cortical thickness [default=0.1]\n");
  fprintf(stdout,"  --version                  show version info\n");
}

int main(int argc,char *argv[],char *envp[]) {
  int clobber=0;

  //double niikcortex_tissue_val[12]= {-1.0,-1.0,-1.0,-1.0, -1.0,-1.0,-1.0,-1.0, -1.0,-1.0,-1.0,-1.0};
  double thresh=50.0;
  double omax=3.0;
  double initial=0.1;
  int    smooth_iter=3;
  double smooth_factor=0.75;
  double delta=0.1;
  double *vlist=NULL;
  int numvlist=0;
  int numimglist=3;
  int maxiter2=20;
  int alt_mode=0;
  nifti_image
    *img=NULL,
   **imglist=NULL,
    *border=NULL;


  kobj *obj=NULL;
  kobj *obj2=NULL;


  const char *in_ctx=NULL;
  const char *in_brain_mask=NULL;
  const char *in_vent=NULL;
  const char *in_gwi_mask=NULL;
  const char *in_mask=NULL;
  const char *in_border=NULL;
  const char *in_init_ics=NULL;
  const char *out_init_ocs=NULL;
  const char *fcname="niikcortex_initocs";
  int i,n;
  char* timestamp = niik_create_minc_timestamp(argc,argv);

  struct option long_options[] = {
    {"clobber", no_argument, &clobber, 1},
    {"version", no_argument, 0, 'v'},
    {"help",    no_argument, 0, 'h'},
    {"usage",   no_argument, 0, 'U'},

    {"vlist",required_argument,0,'V'},
    {"border",required_argument, 0, 'B'},

    {"max-thick",required_argument, 0,   'M'},
    {"smooth-iter",required_argument, 0, 'S'},
    {"smooth-factor",required_argument, 0,'F'},
    {"iter",required_argument,  0,       'T'},
    {"delta",required_argument, 0,       'd'},
    {"alt", no_argument,        &alt_mode, 1},
    {0, 0, 0, 0}
  };

  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vhUV:B:M:S:F:T:d:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      break;
    case 'v':
      prog_version();
      return 0;
    case 'V':
      if( (vlist = niik_csstring_to_double_list(optarg,&numvlist))==NULL ) {
        fprintf(stderr,"ERROR: niik_csstring_to_list: %s\n",optarg);
        exit(9);
      }
      break;
    case 'M':
      omax = atof(optarg);
      break;
    case 'S':
      smooth_iter = atoi(optarg);
      break;
    case 'F':
      smooth_factor = atoi(optarg);
      break;
    case 'd':
      delta = atof(optarg);
      break;
    case 'T':
      maxiter2 = atoi(optarg);
      break;
    case 'h':
    case 'U':
    case '?':
    default:
      prog_usage ();
      return 1;
    }
  }

  if((argc - optind)<3) {
    prog_usage();
    return 1;
  }
#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  in_mask      = argv[optind];
  in_init_ics  = argv[optind+1];
  out_init_ocs = argv[optind+2];

  if( (img=niik_image_read(in_mask))==NULL) {
    fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,in_ctx);
    exit(1);
  }
  NIIK_RET1((!niik_image_type_convert(img,NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert");

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);

  /*from niikmath*/

  if(vlist==NULL) {  /* default expansion step sizes */
    numvlist=3;
    vlist=(double *)calloc(numvlist,sizeof(double));
    vlist[0]=0.10; /*0.1,0.05,0.01*/
    vlist[1]=0.05;
    vlist[2]=0.01;
  }

  if(niik_check_double_problem(omax)) { /* maximum thickness */
    omax = 3.0;
  }
  if(niik_check_double_problem(delta)) { /* initial cortical thickness */
    delta = 0.1;
  }
  if(niik_check_double_problem(smooth_factor)) { /* smoothing */
    smooth_factor = 0.75;
  }
  if(smooth_iter<0) {
    smooth_iter=3;
  }

  if((obj=off_kobj_read_offply(in_init_ics))==NULL) {
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

  fprintf(stdout,"[%s] reading object    %s as initial pial surface\n",fcname,obj->fname);

  NIIK_RET1((obj2 = off_kobj_copy(obj))==NULL,fcname,"off_kobj_copy");

  if(!niikcortex_initocs_expand(img,
                                img,vlist,numvlist,
                                -1.0,-1.0,-1.0,
                                -1.0,-1.0,
                                0.0,
                                omax,delta,
                                3,smooth_iter,
                                maxiter2,
                                smooth_factor,
                                obj,obj2,
                                border,
                                alt_mode)) {
    fprintf(stderr,"[%s] ERROR: niikcortex_initocs_expand\n",fcname);
    exit(1);
  }

  if(!off_kobj_add_one_color(obj,1,0,0)) {
    fprintf(stderr,"[%s] ERROR: off_kobj_add_one_color(obj,1,0,0)\n",fcname);
    exit(1);
  }

  /*append metadata*/
  NIIK_RET1((!off_kobj_add_comment(obj2,timestamp)),fcname,"off_kobj_add_comment");

  fprintf(stdout,"[%s] writing output  %s\n",fcname,out_init_ocs);

  if(!off_kobj_write_offply(out_init_ocs,obj2,0)) {
    fprintf(stderr,"[%s] ERROR: off_kobj_write_offply(outname,obj2,0)\n",fcname);
    exit(1);
  }

  /*free memory*/

  obj=off_kobj_free(obj);
  obj2=off_kobj_free(obj2);

  free(vlist);
  img=niik_image_free(img);

  niik_fc_display(fcname,0);
  free(timestamp);
  exit(0);
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
