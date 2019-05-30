/* FILENAME:     niik_stx_registration.c
 * DESCRIPTION:  Kunio's nifti1 stx registration program
 * AUTHOR:       Kunio Nakamura
 * DATE:         January 29, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  niik_stx_registration history\n"
  "  0.0.0 February 2, 2013; knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_describe[] = {
  "  [niik_stx_registration] description\n"
};

static char *prog_help[] = {
  "  niik_stx_registration:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
  "  -check=<out.nii>                  : output image to check registration\n"
  "  -nmi                              : use normalized mutual info [default]\n"
  "  -cc                               : use correlation coefficient\n"
};

void usage() {
  fprintf(stdout,"niik_stx_registration\n");
  fprintf(stdout,"  usage: [options] <mov.mnc> <target.mnc> <out.xfm>\n\n");
  /*(  fprintf(stdout,"\n  use nifti format for images and brain masks\n");*/
}


int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  **img=NULL;
  int
  verbose=1,
  n,nc,sc;
  char
  fcname[32]="niik_stx_registration";
  niikmat
  *afmat=NULL,
   *imat=NULL;
  double
  imax_percent=0.9,
  omin=0,
  omax=140,
  *affpar,
  dist=3.5,
  filFWHM=5.0;
  int
  areg_cost=NIIK_REGISTER_NMI;
  char *outcheckname=NULL;

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--describe",10)) {
        fprintf(stdout,"%s",*prog_describe);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-imax_percent=",14)) {
        imax_percent=atof(argv[nc]+14);
      }

      else if(!strncmp(argv[nc],"-check=",7)) {
        outcheckname=argv[nc]+7;
      }

      else if(!strncmp(argv[nc],"-verbose",8)) {
        verbose=1;
      } else if(!strncmp(argv[nc],"-quiet",6)) {
        verbose=0;
      }

      else if(!strncmp(argv[nc],"-init=",6)) {
        fprintf(stdout,"[%s] reading init xfm:       %s\n",fcname,argv[nc]+6);
        NIIK_EXIT(((imat = niikmat_read_xfm(argv[nc]+6))==NULL),fcname,"niikmat_read_xfm",9);
      }

      else if(!strncmp(argv[nc],"-omin=",6)) {
        omin=atof(argv[nc]+6);
      } else if(!strncmp(argv[nc],"-omax=",6)) {
        omax=atof(argv[nc]+6);
      }

      else if(!strncmp(argv[nc],"-nmi",4)) {
        areg_cost=NIIK_REGISTER_NMI;
      } else if(!strncmp(argv[nc],"-cc",3)) {
        areg_cost=NIIK_REGISTER_CC;
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
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

  if(argc>4) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(0);
  } else if(argc<4) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(0);
  }

  niik_fc_display(fcname,1);
  NIIK_EXIT(((    img=(nifti_image **)calloc(4,sizeof(nifti_image *)))==NULL),fcname,"calloc for img",9);

  fprintf(stdout,"[%s] reading moving image:   %s\n",fcname,argv[1]);
  NIIK_EXIT(((img[0]=niik_image_read(argv[1]))==NULL),fcname,"reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img[0],NIFTI_TYPE_FLOAT32,1)),fcname,"niik_image_type_convert",9);
  NIIK_EXIT((!niik_image_iscale(img[0],niik_image_get_min(img[0],NULL),niik_image_get_percentile(img[0],NULL,imax_percent),omin,omax)),
            fcname,"niik_image_iscale",9);

  fprintf(stdout,"[%s] reading target image:   %s\n",fcname,argv[2]);
  NIIK_EXIT(((img[1]=niik_image_read(argv[2]))==NULL),fcname,"reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img[1],NIFTI_TYPE_FLOAT32,1)),fcname,"niik_image_type_convert",9);

  /*mat44_display(img[0]->sto_xyz);
    mat44_display(img[1]->sto_xyz);
    exit(0);*/

  if(imat!=NULL) {
    fprintf(stdout,"[%s] using initial matrix\n",fcname);
    niikmat_display(imat);
    if(verbose>=1) fprintf(stdout,"[%s] convert niikmat\n",fcname);
    NIIK_EXIT((!niikmat_convert_from_xfm(imat,img[0],img[1])),fcname,"niikmat_convert_from_xfm",9);
    if(verbose>=1) niikmat_display(imat);
  }

  affpar=(double *)calloc(20,sizeof(double));
  affpar[7]=affpar[8]=affpar[9]=affpar[10]=1;

  NIIK_EXIT((!niik_aregister_align_mni(img[1],NULL,img[0],NULL,affpar,areg_cost,filFWHM,dist)),fcname,"niik_aregister_align_mni",9);
  NIIK_EXIT(((afmat=niik_aregister_matrix_from_affpar(affpar))==NULL),fcname,"niik_aregister_matrxi_from_affpar",9);
  NIIK_EXIT((!niikmat_convert_to_xfm(afmat,img[0],img[1])),fcname,"niikmat_convert_to_xfm",9);

  fprintf(stdout,"[%s] writing output :        %s\n",fcname,argv[3]);
  NIIK_EXIT((!niikmat_write_as_linear_xfm(argv[3],afmat)),fcname,"niikmat_write_as_linear_xfm",9);

  if(outcheckname!=NULL) {
    fprintf(stdout,"[%s] writing output check: %s\n",fcname,outcheckname);
    NIIK_EXIT((!niik_aregister_matrix_from_affpar_update(afmat,affpar)),fcname,"niik_aregister_matrxi_from_affpar_update",9);
    NIIK_EXIT(((img[2]=niik_image_affine_transform_3d(img[0],img[1],afmat,NIIK_INTERP_BSPLINE))==NULL),fcname,"niik_image_affine_transform_3d",9);
    img[3]=img[1];
    NIIK_EXIT((!niik_image_combine_and_write_as_type(outcheckname,img+2,2,'t',120,NIFTI_TYPE_UINT8)),fcname,"niik_image_combine_and_write_as_type",9);
    img[2]=niik_image_free(img[2]);
  }

  imat = niikmat_free(imat);
  afmat = niikmat_free(afmat);
  free(affpar);
  affpar=NULL;
  for(n=0; n<2; n++) {
    img[n]=niik_image_free(img[n]);
  }
  free(img);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_stx_registration */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/