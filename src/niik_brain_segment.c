/* Filename:     niik_brain_segment.c
 * Description:  Kunio's nifti1 rough brain segmentation program
 * Author:       Kunio Nakamura
 * Date:         January 18, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  niik_brain_segment history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_help[] = {
  "  optional usage:\n"
  "  -u -help --help              : show this usage\n"
  "  --version                    : show version info\n"
  "  -close_dilate_radius=<r>     : radius for brain close dilation\n"
  "  -close_erode_radius=<r>      : radius for brain close erosion\n"
  "  -uthresh2=<ut2>              : upper threshold (2)\n"
  "  -uthresh=<ut>                : upper threshold\n"
  "  -ithresh=<it>                : initial threshold\n"
  "  -thresh=<thresh>             : threshold\n"
  "  -radius=<r>                  : radius\n"
  "  -delta=<delta>               : delta\n"
  "  -FWHM=<FWHM>                 : FWHM\n"
  "  -imat=<M.mat>                : affine stx matrix\n"
  "  -imatxfm=<M.xfm>             : linear stx xfm file\n"
};

void usage() {
  fprintf(stdout,"niik_brain_segment\n");
  fprintf(stdout,"  usage: [options] <img.nii|img.mnc> <out.nii|out.mnc>\n");
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *stx_img=NULL,
   *stx_mask=NULL,
    *maskimg=NULL,
     *img=NULL;
  int
  iter=2,
  nc,sc;
  const char
  *NIIKDIR=NULL;
  char
  fcname[32]="niik_brain_segment";
  char
  stx_img_name[4096],
               stx_mask_name[4096];
  double
  thresh=NIIKMAX,
  ithresh=NIIKMAX,
  uthresh=NIIKMAX,
  uthresh2=NIIKMAX,
  radius=3.6,
  *affpar=NULL,
   delta=3.2,
   FWHM=6.0,
   close_dilate_radius=5.5,
   close_erode_radius=4.5;
  int afmat_xfm=0;
  niikmat
  *afmat=NULL;
  niikpt pt,ctr;

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;
  stx_img_name[0]=stx_mask_name[0]=0;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*prog_help);
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*prog_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"--stx-mask=",11)) {
        sprintf(stx_mask_name,"%s",argv[nc]+11);
        fprintf(stdout,"[%s] switched to %s for stx mask\n",
                fcname,stx_mask_name);
      } else if(!strncmp(argv[nc],"--stx-img=",10)) {
        sprintf(stx_img_name,"%s",argv[nc]+10);
        fprintf(stdout,"[%s] switched to %s for stx img\n",
                fcname,stx_img_name);
      }

      else if(!strncmp(argv[nc],"-close_dilate_radius=",20)) {
        close_dilate_radius=atof(argv[nc]+20);
      } else if(!strncmp(argv[nc],"-close_erode_radius=",19)) {
        close_erode_radius=atof(argv[nc]+19);
      } else if(!strncmp(argv[nc],"-uthresh2=",10)) {
        uthresh2=atof(argv[nc]+10);
      }

      else if(!strncmp(argv[nc],"-imatxfm=",9)) {
        afmat_xfm=1;
        NIIK_EXIT(((afmat = niikmat_read_xfm(argv[nc]+9))==NULL),
                  fcname,"niikmat_read_xfm",9);
      }

      else if(!strncmp(argv[nc],"-ithresh=",9)) {
        ithresh=atof(argv[nc]+9);
      } else if(!strncmp(argv[nc],"-uthresh=",9)) {
        uthresh=atof(argv[nc]+9);
      } else if(!strncmp(argv[nc],"-radius=",8)) {
        radius=atof(argv[nc]+8);
      } else if(!strncmp(argv[nc],"-delta=",7)) {
        delta=atof(argv[nc]+7);
      } else if(!strncmp(argv[nc],"-FWHM=",6)) {
        FWHM=atof(argv[nc]+6);
      }

      else if(!strncmp(argv[nc],"-imat=",6)) {
        fprintf(stdout,"[%s] reading matrix %s\n",fcname,argv[nc]+6);
        NIIK_EXIT(((afmat=niikmat_read(argv[nc]+6))==NULL),fcname,"niikmat_read",9);
      }

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

  if(argc!=3) {
    fprintf(stderr,"[%s] ERROR: wrong usage\n",fcname);
    exit(0);
  }

  fprintf(stdout,"[%s] reading input image:    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img,NIFTI_TYPE_FLOAT32,1)),
            fcname,"niik_image_type_convert",9);

  if(afmat==NULL) {
    fprintf(stdout,"[%s] stx 12param alignment:\n",fcname);
    fprintf(stdout,"[%s]   FWHM    : %7.2f\n",fcname,FWHM);
    fprintf(stdout,"[%s]   delta   : %7.2f\n",fcname,delta);
    if(stx_img_name[0]==0) {
      NIIK_EXIT(((NIIKDIR=get_NIIKDIR())==NULL),
                fcname,
                "please setenv | export  NIIKDIR",9);
      sprintf(stx_img_name,"%s/data/bseg/MNI152_T1_1mm.nii.gz",NIIKDIR);
    }
    if(stx_mask_name[0]==0) {
      NIIK_EXIT(((NIIKDIR=get_NIIKDIR())==NULL),
                fcname,
                "please setenv | export  NIIKDIR",9);
      sprintf(stx_mask_name,"%s/data/bseg/MNI152_T1_1mm_brain_mask.nii.gz",NIIKDIR);
    }

    fprintf(stdout,"[%s] reading stx img         %s\n",fcname,stx_img_name);
    NIIK_EXIT(((stx_img=niik_image_read(stx_img_name))==NULL),
              fcname,"niik_image_read",9);
    fprintf(stdout,"[%s] reading stx mask        %s\n",fcname,stx_mask_name);
    NIIK_EXIT(((stx_mask=niik_image_read(stx_mask_name))==NULL),
              fcname,"niik_image_read",9);

    affpar = (double *)calloc(20,sizeof(double));
    affpar[7]=affpar[8]=affpar[9]=affpar[10]=1;
    ctr=niikpt_image_get_centroid(img,NULL);
    affpar[14]=ctr.x;
    affpar[15]=ctr.y;
    affpar[16]=ctr.z;
    pt=niikpt_sub(niikpt_image_get_centroid(img,NULL),niikpt_image_get_centroid(stx_img,NULL));
    if(affpar[4]==0) affpar[4]=pt.x;
    if(affpar[5]==0) affpar[5]=pt.y;
    if(affpar[6]==0) affpar[6]=pt.z;
    niik_aregister_display_affine(affpar);

    NIIK_EXIT((!niik_aregister_align_mni(stx_img,stx_mask,img,NULL,affpar,NIIK_REGISTER_NMI,FWHM,delta)),
              fcname,"niik_aregister_align_mni",9);
    niik_aregister_display_affine(affpar);

    NIIK_EXIT(((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL),fcname,
              "niik_aregister_matrix_from_affpar",9);
  }

  else if(afmat_xfm) {
    fprintf(stdout,"[%s] convert xfm to niikmat\n",fcname);
    if(stx_img_name[0]==0) {
      NIIK_EXIT(((NIIKDIR=get_NIIKDIR())==NULL),
                fcname,"please setenv | export  NIIKDIR",9);
      sprintf(stx_img_name,"%s/data/bseg/MNI152_T1_1mm.nii.gz",NIIKDIR);
    }
    NIIK_EXIT(((stx_img=niik_image_read(stx_img_name))==NULL),
              fcname,"niik_image_read",9);
    NIIK_EXIT((!niikmat_convert_from_xfm(afmat,img,stx_img)),
              fcname,"niikmat_convert_from_xfm",9);
    stx_img=niik_image_free(stx_img);
  } /* afmat_xfm */

  fprintf(stdout,"[%s] segmentation:\n",fcname);
  if(niik_check_double_problem(  thresh)) fprintf(stdout,"[%s]   thresh    : undef\n",fcname);
  else fprintf(stdout,"[%s]   thresh    : %7.2f\n",fcname,thresh);
  if(niik_check_double_problem( ithresh)) fprintf(stdout,"[%s]   ithresh   : undef\n",fcname);
  else fprintf(stdout,"[%s]   ithresh   : %7.2f\n",fcname,ithresh);
  if(niik_check_double_problem( uthresh)) fprintf(stdout,"[%s]   uthresh   : undef\n",fcname);
  else fprintf(stdout,"[%s]   uthresh   : %7.2f\n",fcname,uthresh);
  if(niik_check_double_problem(uthresh2)) fprintf(stdout,"[%s]   uthresh2  : undef\n",fcname);
  else fprintf(stdout,"[%s]   uthresh2  : %7.2f\n",fcname,uthresh);
  fprintf(stdout,"[%s]   radius    : %7.2f\n",fcname,radius);

  NIIK_EXIT(((maskimg=niik_image_bseg_basic(img,iter,thresh,ithresh,uthresh,uthresh2,radius,afmat,0))==NULL),fcname,
            "niik_image_bseg_basic",9);
  maskimg->sto_xyz=img->sto_xyz;
  maskimg->qto_xyz=img->qto_xyz;

  fprintf(stdout,"[%s] brain closing: %7.3f %7.3f\n",fcname,close_dilate_radius,close_erode_radius);
  NIIK_EXIT((!niik_image_morph_close_brain(maskimg,close_dilate_radius,close_erode_radius)),
            fcname,"niik_image_morph_close_brain",9);
  fprintf(stdout,"[%s] vol %.4f\n",fcname,niik_image_get_mask_vol(maskimg));

  fprintf(stdout,"[%s] writing brain mask      %s\n",fcname,argv[2]);
  NIIK_EXIT((!niik_image_write(argv[2],maskimg)),fcname,"niik_image_write",9);

  free(affpar);
  img=niik_image_free(img);
  stx_img=niik_image_free(stx_img);
  stx_mask=niik_image_free(stx_mask);
  exit(0);
} /* niik_brain_segment */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/