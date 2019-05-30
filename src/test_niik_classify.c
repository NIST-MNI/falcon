/* Filename:     niik_classify.c
 * Description:  Kunio's nifti1 classification program
 * Author:       Kunio Nakamura
 * Date:         December 28, 2012
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)



static char *niik_classify_version[] = {
  "  niik_classify version history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *niik_classify_help[] = {
  "  usage: <out.nii> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help                   : show this usage\n"
  "  --version                  : show version info\n"
  "  -t1w <t1w.nii>             : input T1-weighted image\n"
  "  -pdw <pdw.nii>             : input PD-weighted image\n"
  "  -t2w <t2w.nii>             : input T2-weighted image\n"
  "  -flr <flr.nii>             : input FLAIR image\n"
  "  -mask <brain.nii>          : input brain parenchyma mask (GM+WM, excluding CSF)\n"
  "\n"
  "  n.b.\n"
  "  -all images need to be in the same space, resampled to the same dimension\n"
};


void usage() {
  fprintf(stdout,"%s",*niik_classify_help);
  return;
}

int niik_classify_prep(niik_naive_bayesian_classifier *C,nifti_image *t1wimg,nifti_image *t1cimg,nifti_image *t2wimg,nifti_image *pdwimg,nifti_image *flrimg,nifti_image *brainmask,nifti_image *a_gm,nifti_image *a_wm,nifti_image *a_csf,nifti_image *a_les,int verbose) {
  char fcname[64]="niik_classify_func_prep";
  int n;

  double
  dval=0,
  imin=0,
  imax=0;
  niikmat
  *histomat=NULL;
  int histo_num=-1;

  if(verbose>=1) niik_fc_display(fcname,1);

  /* add images to the classifier */
  if(!niik_naive_bayesian_classifier_brain_classify_add_images(C,t1wimg,t1cimg,t2wimg,pdwimg,flrimg)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_images\n",fcname);
    return 0;
  }
  /* memory allocation for probability images */
  if(!niik_naive_bayesian_classifier_alloc_primg(C,C->mask)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_alloc_primg\n",fcname);
    return 0;
  }
  niik_naive_bayesian_classifier_display(C);

  /* write probability maps */
  if(!niik_naive_bayesian_classifier_brain_classify_calc_probability(C)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_calc_probability\n",fcname);
    return 0;
  }
  if(!niik_image_combine_and_write_as_type("tmp_class_pr0.nii.gz",C->primg,5,'t',0,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write_as_type\n",fcname);
    return 0;
  }

  /* brain histogram fitting from t1w */
  imin=niik_image_get_min(t1wimg,NULL);
  imax=0.8*niik_image_get_max(t1wimg,brainmask) + 0.2*niik_image_get_max(t1wimg,NULL);
  histomat=niikmat_init(3,5);
  if(histo_num<0) histo_num=101;
  fprintf(stdout,"[%s] min/max %9.5f %9.5f\n",fcname,imin,imax);
  if(!niik_image_fit_gaussian_tissues(C->img[0],C->mask,imin,imax,histo_num,0,histomat->m[0],histomat->m[1],histomat->m[2],&dval)) {
    fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian_tissues\n",fcname);
    exit(0);
  }
  fprintf(stdout,"[%s] GM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][0],histomat->m[0][0],histomat->m[1][0]);
  fprintf(stdout,"[%s] WM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][1],histomat->m[0][1],histomat->m[1][1]);
  fprintf(stdout,"[%s] CSF   %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][2],histomat->m[0][2],histomat->m[1][2]);

  C->mean->m[0][0] = histomat->m[0][2]/2.0;
  C->mean->m[0][1] = histomat->m[0][2];
  C->mean->m[0][2] = histomat->m[0][0];
  C->mean->m[0][3] = histomat->m[0][1];
  C->mean->m[0][4] = histomat->m[0][0];
  C->stdv->m[0][0] = histomat->m[0][2]/2.5;
  C->stdv->m[0][1] = histomat->m[1][2];
  C->stdv->m[0][2] = histomat->m[1][0];
  C->stdv->m[0][3] = histomat->m[1][1];
  C->stdv->m[0][4] = histomat->m[1][0] * 2.0;

  C-> num_img = 1;
  C-> num_class = 4;
  niik_naive_bayesian_classifier_display(C);

  niik_image_morph_close_brain(C->mask,5.5,3.5);

  if(!niik_naive_bayesian_classifier_brain_classify_calc_probability(C)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_calc_probability\n",fcname);
    return 0;
  }

  /* write probability maps */
  if(!niik_naive_bayesian_classifier_brain_classify_calc_probability(C)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_calc_probability\n",fcname);
    return 0;
  }
  if(!niik_image_combine_and_write_as_type("tmp_class_pr1.nii.gz",C->primg,5,'t',0,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write_as_type\n",fcname);
    return 0;
  }

  /*
  niik_naive_bayesian_classifier_brain_classify_calc_ML(C);
  if(!niik_naive_bayesian_classifier_brain_classify_calc_probability(C)){
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_calc_probability\n",fcname);
    return 0; }
  if(!niik_image_combine_and_write_as_type("tmp_class_pr1.nii.gz",C->primg,5,'t',0,NIFTI_TYPE_FLOAT32)){
    fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write_as_type\n",fcname);
    return 0; }
  */

  niik_image_write("tmp_pr_bg.nii.gz",C->primg[0]);
  niik_image_write("tmp_pr_csf.nii.gz",C->primg[1]);
  niik_image_write("tmp_pr_gm.nii.gz",C->primg[2]);
  niik_image_write("tmp_pr_wm.nii.gz",C->primg[3]);
  niik_image_write("tmp_pr_les.nii.gz",C->primg[4]);

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}

int niik_classify_test1(niik_naive_bayesian_classifier *C,nifti_image *t1wimg,nifti_image *t1cimg,nifti_image *t2wimg,nifti_image *pdwimg,nifti_image *flrimg,nifti_image *brainmask,nifti_image *a_gm,nifti_image *a_wm,nifti_image *a_csf,nifti_image *a_les,int verbose) {
  char fcname[64]="niik_classify_test1";

  nifti_image *brainroi=NULL;
  double
  dval=0,
  imin=0,
  imax=0;
  niikmat
  *histomat=NULL;
  int histo_num=-1;

  int i,n;
  double
  dlist[12];

  if(verbose>=1) niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] add images\n",fcname);
  /* add images to the classifier */
  if(!niik_naive_bayesian_classifier_brain_classify_add_images(C,t1wimg,t1cimg,t2wimg,pdwimg,flrimg)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_images\n",fcname);
    return 0;
  }
  /* memory allocation for probability images */
  if(!niik_naive_bayesian_classifier_alloc_primg(C,C->mask)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_alloc_primg\n",fcname);
    return 0;
  }
  niik_naive_bayesian_classifier_display(C);


  /* brain histogram fitting from t1w */
  imin=niik_image_get_min(t1wimg,NULL);
  imax=0.8*niik_image_get_max(t1wimg,brainmask) + 0.2*niik_image_get_max(t1wimg,NULL);
  histomat=niikmat_init(3,5);
  if(histo_num<0) histo_num=101;
  fprintf(stdout,"[%s] min/max %9.5f %9.5f\n",fcname,imin,imax);
  if(!niik_image_fit_gaussian_tissues(t1wimg,brainmask,imin,imax,histo_num,0,histomat->m[0],histomat->m[1],histomat->m[2],&dval)) {
    fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian_tissues\n",fcname);
    return 0;
  }
  fprintf(stdout,"[%s] GM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][0],histomat->m[0][0],histomat->m[1][0]);
  fprintf(stdout,"[%s] WM    %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][1],histomat->m[0][1],histomat->m[1][1]);
  fprintf(stdout,"[%s] CSF   %9.2f   %9.5f  +/-  %9.5f\n",fcname,histomat->m[2][2],histomat->m[0][2],histomat->m[1][2]);

  fprintf(stdout,"[%s] brainroi\n",fcname);
  brainroi = niik_image_copy_as_type(brainmask,NIFTI_TYPE_UINT8);
  niik_image_morph_close_brain(brainroi,5.5,3.5);

  fprintf(stdout,"[%s] process\n",fcname);
  for(i=0; i<brainmask->nvox; i++) {
    if(niik_image_get_voxel(brainroi,i)==0) continue;
    dval = niik_image_get_voxel(t1wimg,i);
    dlist[0] = NIIK_GaussPDF(dval - histomat->m[0][0],histomat->m[1][0]);
    dlist[1] = NIIK_GaussPDF(dval - histomat->m[0][1],histomat->m[1][1]);
    dlist[2] = NIIK_GaussPDF(dval - histomat->m[0][2],histomat->m[1][2]);
    if(dlist[0]>=dlist[1]) {
      if(dlist[0]>=dlist[2])
        niik_image_set_voxel(brainroi,i,2);
    } else if(dlist[1]>=dlist[0]) {
      if(dlist[1]>=dlist[2])
        niik_image_set_voxel(brainroi,i,3);
    } else if(dlist[2]>=dlist[0]) {
      if(dlist[2]>=dlist[1])
        niik_image_set_voxel(brainroi,i,1);
    }
  }

  fprintf(stdout,"[%s] write\n",fcname);
  niik_image_write("tmp_class_1.nii.gz",brainroi);

  fprintf(stdout,"[%s] calculation\n",fcname);
  niik_naive_bayesian_classifier_display(C);

  for(n=1; n<=3; n++) {
    C->mean->m[1][n] = niik_image_get_mode_label(pdwimg,brainroi,n,0,5000,1000,10);
    C->mean->m[2][n] = niik_image_get_mode_label(t2wimg,brainroi,n,0,5000,1000,10);
    C->stdv->m[1][n] = niik_image_get_stdv_label(pdwimg,brainroi,n);
    C->stdv->m[2][n] = niik_image_get_stdv_label(t2wimg,brainroi,n);
  }

  C->mean->m[0][2] = histomat->m[0][0];
  C->stdv->m[0][2] = histomat->m[1][0];
  C->mean->m[0][3] = histomat->m[0][1];
  C->stdv->m[0][3] = histomat->m[1][1];
  C->mean->m[0][1] = histomat->m[0][2];
  C->stdv->m[0][1] = histomat->m[1][2];

  C->num_class = 4;
  C->num_img = 2;

  for(n=0; n<C->num_img; n++) {
    C->stdv->m[n][0] = C->mean->m[n][1] / 2.0;
    C->mean->m[n][0] = 0.0;
  }
  C->stdv->m[0][0]=C->mean->m[0][1] / 3.0;
  C->stdv->m[1][0]=C->mean->m[1][3] / 3.0;
  C->stdv->m[2][0]=C->mean->m[2][3] / 3.0;

  niik_naive_bayesian_classifier_display(C);
  niik_image_morph_close_brain(C->mask,5.5,3.5);

  /* write probability maps */
  if(!niik_naive_bayesian_classifier_brain_classify_calc_probability(C)) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_calc_probability\n",fcname);
    return 0;
  }
  if(!niik_image_combine_and_write_as_type("tmp_class_pr1.nii.gz",C->primg,C->num_class,'t',0,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_combine_and_write_as_type\n",fcname);
    return 0;
  }

  histomat = niikmat_free(histomat);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}


int main(int argc,char *argv[],char *envp[]) {
  niik_naive_bayesian_classifier *C=NULL;
  nifti_image
  *a_gm=NULL,
   *a_wm=NULL,
    *a_csf=NULL,
     *a_les=NULL,
      *pdwimg=NULL,
       *t2wimg=NULL,
        *t1wimg=NULL,
         *t1cimg=NULL,
          *flrimg=NULL,
           *brainmask=NULL;
  char
  *outcheckname=NULL,
   fcname[20]="niik_classify";
  int
  num_img=0,
  num_class=5,
  verbose=2,
  nc=1,sc=1;
  niikmat
  *mnireg=NULL;

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

      else if(!strncmp(argv[nc],"-mask",5)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(brainmask!=NULL) {
          fprintf(stderr,"[%s] ERROR: brain mask already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading brain    %s\n",fcname,argv[nc]);
        if((brainmask=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* mask */
      else if(!strncmp(argv[nc],"-t1w",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(t1wimg!=NULL) {
          fprintf(stderr,"[%s] ERROR: t1w already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading T1w      %s\n",fcname,argv[nc]);
        if((t1wimg=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* t1w */
      else if(!strncmp(argv[nc],"-t2w",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(t2wimg!=NULL) {
          fprintf(stderr,"[%s] ERROR: t2w already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading T2w      %s\n",fcname,argv[nc]);
        if((t2wimg=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* t2w */
      else if(!strncmp(argv[nc],"-flr",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(flrimg!=NULL) {
          fprintf(stderr,"[%s] ERROR: FLAIR already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading FLAIR    %s\n",fcname,argv[nc]);
        if((flrimg=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* t1w */
      else if(!strncmp(argv[nc],"-t1c",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(t1cimg!=NULL) {
          fprintf(stderr,"[%s] ERROR: t1w+contrast already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading T1c     %s\n",fcname,argv[nc]);
        if((t1cimg=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* t1c */
      else if(!strncmp(argv[nc],"-pdw",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(0);
        }
        if(pdwimg!=NULL) {
          fprintf(stderr,"[%s] ERROR: PDw already used\n",fcname);
          exit(0);
        }
        fprintf(stdout,"[%s] reading PDw      %s\n",fcname,argv[nc]);
        if((pdwimg=niik_image_read(argv[nc]))==NULL) {
          fprintf(stderr,"[%s] ERROR: niik_image_read %s\n",fcname,argv[nc]);
          exit(0);
        }
      } /* pdw */

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niik_classify_help);
        exit(0);
      } else if(!strncmp(argv[nc],"-u",2)) {
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

  if(argc!=2) {
    fprintf(stderr,"[%s] ERROR: wrong usage\n",fcname);
    exit(0);
  }

  if(brainmask==NULL) {
    fprintf(stderr,"[%s] ERROR: missing brain mask (-mask)\n",fcname);
    exit(0) ;
  }

  num_img=0;
  num_img=(t1wimg!=NULL)?num_img+1:num_img;
  num_img=(t1cimg!=NULL)?num_img+1:num_img;
  num_img=(t2wimg!=NULL)?num_img+1:num_img;
  num_img=(pdwimg!=NULL)?num_img+1:num_img;
  num_img=(flrimg!=NULL)?num_img+1:num_img;
  fprintf(stdout,"[%s] using %i images\n",fcname,num_img);

  if((C = niik_naive_bayesian_classifier_init(num_img,num_class))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_init\n",fcname);
    exit(0);
  }
  C->mask = niik_image_copy_as_type(brainmask,NIFTI_TYPE_UINT8);

  fprintf(stdout,"[%s] start niik_classify_test1\n",fcname);
  if(!niik_classify_test1(C,t1wimg,t1cimg,t2wimg,pdwimg,flrimg,brainmask,a_gm,a_wm,a_csf,a_les,verbose)) {
    fprintf(stderr,"[%s] ERROR: niik_classify_test1\n",fcname);
    exit(0);
  }

  brainmask=niik_image_free(brainmask);
  t1wimg = niik_image_free(t1wimg);
  t1cimg = niik_image_free(t1cimg);
  t2wimg = niik_image_free(t2wimg);
  pdwimg = niik_image_free(pdwimg);
  flrimg = niik_image_free(flrimg);

  niik_fc_display(fcname,0);
  exit(0);
} /* niik_classify */



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/