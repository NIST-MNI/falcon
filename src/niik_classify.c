/* Filename:     niik_classify.c
 * Description:  Kunio's nifti1 classification program
 * Author:       Kunio Nakamura
 * Date:         December 28, 2012
 */

#include "falcon.h"

#define MAJOR_VERSION (NIIK_MAJOR_VERSION)
#define MINOR_VERSION (NIIK_MINOR_VERSION)
#define MICRO_VERSION (NIIK_MICRO_VERSION)


static char *niik_classify_version[] = {
  "  niik_classify version history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *niik_classify_help[] = {
  "  simple classification program\n"
  "  usage: <in.nii> <mask.nii> <out.nii> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help                   : show this usage\n"
  "  --version                  : show version info\n"
  "  -kmean                     : use k-means classification (default)\n"
  "  -fcm                       : use fuzzy c-means classification (under developement)\n"
  "\n"
  "  <in.nii> can be 4-dimensional in u-dimension\n"
  "\n"
};


void usage() {
  fprintf(stdout,"%s",*niik_classify_help);
  return;
}

int main(int argc,char *argv[],char *envp[]) {
  niik_classify_kmeans *K=NULL;
  niik_classify_fcm *F=NULL;
  nifti_image
  *img=NULL,
   *brainmask=NULL;
  char
  *outcheckname=NULL,
   fcname[20]="niik_classify";
  int
  class_type=NIIK_CLASSIFY_KMEANS,
  maxiter=5,
  num_class=3,
  verbose=2,
  n,
  nc=1,sc=1;
  double
  FWHM=0;

  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  if(argc==1) {
    usage();
    exit(0);
  }

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);
  fprintf(stdout,"  niik_classifier: version %i.%i.%i\n",MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  fprintf(stdout,"  niik_classifier: executed at %s\n",tmstr);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niik_classify_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*niik_classify_help);
        exit(0);
      } else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      } /* outcheck */

      else if(!strncmp(argv[nc],"-kmeans",7)) {
        class_type=NIIK_CLASSIFY_KMEANS;
      } /* kmeans */

      else if(!strncmp(argv[nc],"-nclass",7)) {
        num_class=atoi(argv[++nc]);
      } /* class */

      else if(!strncmp(argv[nc],"-iter",5)) {
        maxiter=atoi(argv[++nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-FWHM",5)) {
        FWHM=atof(argv[++nc]);
      } /* FWHM */

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niik_classify_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-fcm",4)) {
        class_type=NIIK_CLASSIFY_FCM;
      } /* fcm */

      else if(!strncmp(argv[nc],"-u",2)) {
        usage();
        exit(0);
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  NIIK_EXIT((argc<4),fcname,"too few argments",9);
  NIIK_EXIT((argc>4),fcname,"too many argments",9);

  fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading brain    %s\n",fcname,argv[2]);
  NIIK_EXIT(((brainmask=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);

  if(fabs(FWHM)>=1e-4) {
    if(verbose>=2) fprintf(stdout,"[%s] gaussian filter %8.3f\n",fcname,FWHM);
    NIIK_RET0((!niik_image_filter_gaussian_update(img,17,FWHM)),fcname,"niik_image_filter_gaussian_update");
  }

  switch(class_type) {
  case NIIK_CLASSIFY_KMEANS:
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_kmeans_init\n",fcname);
    K = niik_classify_kmeans_init(img,brainmask,num_class,maxiter,0);
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_kmeans_initialize_means\n",fcname);
    NIIK_RET0((!niik_classify_kmeans_initialize_means(K,0)),fcname,"niik_classify_kmeans_initialize_means");
    NIIK_RET0((!niik_classify_kmeans_display(K)),fcname,"niik_classify_kmeans_display");
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_kmeans_iterate\n",fcname);
    NIIK_RET0((!niik_classify_kmeans_iterate(K,0)),fcname,"niik_classify_kmeans_iterate");
    niik_classify_kmeans_display(K);
    if(K->img->ndim==3) {
      for(n=0; n<=K->num_class; n++) {
        niik_image_label_display_stats(K->img,K->mask,n);
      }
    }
    fprintf(stdout,"[%s] writing output   %s\n",fcname,argv[3]);
    NIIK_EXIT((!niik_image_write(argv[3],K->mask)),fcname,"niik_image_write",9);
    K=niik_classify_kmeans_free(K);
    break;

  case NIIK_CLASSIFY_FCM:
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_fcm_init\n",fcname);
    F = niik_classify_fcm_init(img,brainmask,num_class,maxiter,0);
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_fcm_initialize_means\n",fcname);
    niik_classify_fcm_initialize_means(F,0);
    niik_classify_fcm_display(F);
    if(verbose>=2) fprintf(stdout,"[%s] niik_classify_kmeans_iterate\n",fcname);
    niik_classify_fcm_iterate(F,0);
    niik_classify_fcm_display(F);
    if(F->img->ndim==3) {
      for(n=0; n<=F->num_class; n++) {
        niik_image_label_display_stats(F->img,F->mask,n);
      }
    }
    fprintf(stdout,"[%s] writing output   %s\n",fcname,argv[3]);
    NIIK_EXIT((!niik_image_write(argv[3],F->primg)),fcname,"niik_image_write",9);
    niik_classify_fcm_display(F);
    F=niik_classify_fcm_free(F);
    break;

  default:
    fprintf(stderr,"[%s] ERROR: unknown classification option, %i\n",fcname,class_type);
    exit(0);
  }

  niik_fc_display(fcname,0);
  exit(0);
} /* niik_classify */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
