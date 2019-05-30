/* Filename:     niik_classify2.c
 * Description:  Kunio's nifti1 classification program
 * Author:       Kunio Nakamura
 * Date:         April 9, 2012
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)


/*********************************************
 * object for classification
 *********************************************/
typedef struct {
  nifti_image **img;
  nifti_image **pr;
  nifti_image **pp;
  nifti_image *brainmask;
  int verbose;
} niik_classifier2;


/*********************************************
 * tissue classification
 *********************************************/
int niik_classifier2_execute(niik_classifier2 *c) {
  int verbose=0;
  char fcname[64]="niik_classifier2_execute";
  nifti_image *subimg=NULL;
  int n,i;
  double thresh;

  verbose=c->verbose;
  if(verbose>=1) niik_fc_display(fcname,1);

  for(n=0; n<3; n++) NIIK_RET0((!niik_image_type_convert(c->img[n],NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");

  fprintf(stdout,"[%s] stats\n",fcname);
  fprintf(stdout,"  t1w        %15.9f  %-15.5f\n",
          niik_image_get_mean(c->img[0],c->brainmask),
          niik_image_get_stdv(c->img[0],c->brainmask));
  fprintf(stdout,"  pdw        %15.9f  %-15.5f\n",
          niik_image_get_mean(c->img[1],c->brainmask),
          niik_image_get_stdv(c->img[1],c->brainmask));
  fprintf(stdout,"  t2w        %15.9f  %-15.5f\n",
          niik_image_get_mean(c->img[2],c->brainmask),
          niik_image_get_stdv(c->img[2],c->brainmask));

  NIIK_RET0(((subimg=niik_image_copy(c->img[1]))==NULL),fcname,"niik_image_copy");

  for(i=0; i<subimg->nvox; i++) {
    niik_image_mul_voxel(subimg,i,2);
    niik_image_add_voxel(subimg,i,-niik_image_get_voxel(c->img[2],i));
  }
  // NIIK_RET0((!niik_image_write("tmp_subimg.mnc",subimg)),fcname,"niik_image_write");

  if(verbose>1) fprintf(stdout,"[%s] threshold calculation\n",fcname);
  NIIK_RET0((!niik_image_thresh_otsu(subimg,c->brainmask,&thresh)),fcname,"niik_image_thresh_otsu");
  fprintf(stdout,"[%s] thresh %.4f\n",fcname,thresh);

  NIIK_RET0(((c->pp=(nifti_image **)calloc(4,sizeof(nifti_image *)))==NULL),fcname,"calloc pp");
  for(n=0; n<4; n++) {
    if(verbose>=2) fprintf(stdout,"[%s] create probability map memory %i\n",fcname,n);
    NIIK_RET0(((c->pp[n]=niik_image_copy_as_type(c->pr[0],NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
  }

  for(i=0; i<subimg->nvox; i++) {
    niik_image_set_voxel(c->pp[0],i,
                         NIIK_Heaviside(niik_image_get_voxel(subimg,i)-thresh,thresh/2.0));
  }

  NIIK_RET0((!niik_image_mask(c->pp[0],c->brainmask)),fcname,"niik_image_mask");
  c->pp[0]->fname=(char *)calloc(512,sizeof(char));
  c->pp[0]->iname=(char *)calloc(512,sizeof(char));
  //NIIK_RET0((!niik_image_write("tmp_pp_brain.nii.gz",c->pp[0])),fcname,"niik_image_write tmp_pp_brain.nii.gz");

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}


static char *niik_classify_version[] = {
  "  niik_classify version history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *niik_classify_help[] = {
  "  simple classification program\n"
  "  usage: <in_T1W.nii> <in_PDW.nii> <in_T2W.nii> <brainmask.nii> <out.nii> [options]\n"
  "\n"
  "  optional usage:\n"
  "  -u -help                   : show this usage\n"
  "  --version                  : show version info\n"
  "\n"
};

void usage() {
  fprintf(stdout,"%s",*niik_classify_help);
  return;
}

int main(int argc,char *argv[],char *envp[]) {
  niik_classifier2 *c=NULL;
  nifti_image
  **img=NULL,
    **pr=NULL,
      *brainmask=NULL;
  char file_pr_brain[512]="/lab2/Kunio/kproj/data/classify/avg152/avg152T1_brain.nii.gz";
  char file_pr_gm   [512]="/lab2/Kunio/kproj/data/classify/avg152/avg152T1_gray.nii.gz";
  char file_pr_wm   [512]="/lab2/Kunio/kproj/data/classify/avg152/avg152T1_white.nii.gz";
  char file_pr_csf  [512]="/lab2/Kunio/kproj/data/classify/avg152/avg152T1_csf.nii.gz";
  char
  fcname[24]="niik_classify2";
  int
  verbose=2,
  n,
  nc=1,sc=1;
  niikmat *afmat=NULL;

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

  NIIK_EXIT(((pr=(nifti_image **)calloc(4,sizeof(nifti_image *)))==NULL),fcname,"calloc for pr",9);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*niik_classify_version);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*niik_classify_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-pr-brain",9)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading prior (BR)   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((pr[0]=niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read pr-brain",9);
      } /* prior */

      else if(!strncmp(argv[nc],"-xfminv",7)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading xfm      %s\n",fcname,argv[nc]);
        NIIK_EXIT(((afmat = niikmat_read_xfm(argv[nc]))==NULL),
                  fcname,"niikmat_read_xfm",9);
        NIIK_EXIT((!niikmat_inverse_update(afmat)),fcname,"niikmat_inverse_update",9);
      } /* xfminv */

      else if(!strncmp(argv[nc],"-pr-CSF",7)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading prior (CSF)  %s\n",fcname,argv[nc]);
        NIIK_EXIT(((pr[3]=niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read pr-GM",9);
      } /* prior */

      else if(!strncmp(argv[nc],"-pr-WM",6)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading prior (WM)   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((pr[2]=niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read pr-WM",9);
      } /* prior */

      else if(!strncmp(argv[nc],"-pr-GM",6)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading prior (GM)   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((pr[1]=niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read pr-GM",9);
      } /* prior */

      else if(!strncmp(argv[nc],"-xfm",4)) {
        NIIK_EXIT((argc<=nc++),fcname,"missing arg",9);
        fprintf(stdout,"[%s] reading xfm      %s\n",fcname,argv[nc]);
        NIIK_EXIT(((afmat = niikmat_read_xfm(argv[nc]))==NULL),
                  fcname,"niikmat_read_xfm",9);
      } /* xfm */

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

  NIIK_EXIT((argc<6),fcname,"too few argments",9);
  NIIK_EXIT((argc>6),fcname,"too many argments",9);


  /*********************************************
   * reading priors
   *********************************************/
  NIIK_EXIT(((img=(nifti_image **)calloc(3,sizeof(nifti_image *)))==NULL),fcname,"calloc for img",9);
  if(pr[0]==NULL) NIIK_EXIT(((pr[0]=niik_image_read(file_pr_brain))==NULL),fcname,"niik_image_read brain",9);
  if(pr[1]==NULL) NIIK_EXIT(((pr[1]=niik_image_read(file_pr_gm   ))==NULL),fcname,"niik_image_read gray",9);
  if(pr[2]==NULL) NIIK_EXIT(((pr[2]=niik_image_read(file_pr_wm   ))==NULL),fcname,"niik_image_read white",9);
  if(pr[3]==NULL) NIIK_EXIT(((pr[3]=niik_image_read(file_pr_csf  ))==NULL),fcname,"niik_image_read csf",9);


  /*********************************************
   * reading images and brain mask
   *********************************************/
  for(n=0; n<3; n++) {
    if(!strncmp(argv[n+1],"NULL",4)) {
      img[n]=NULL;
    } else {
      fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[n+1]);
      NIIK_EXIT(((img[n]=niik_image_read(argv[n+1]))==NULL),fcname,"niik_image_read t1w",9);
    }
  }
  fprintf(stdout,"[%s] reading mask     %s\n",fcname,argv[4]);
  NIIK_EXIT(((brainmask=niik_image_read(argv[4]))==NULL),fcname,"niik_image_read t1w",9);


  /*********************************************
   * bring the priors to image space
   *********************************************/
  if(verbose>=1) fprintf(stdout,"[%s] convert xfm to niikmat\n",fcname);
  if(afmat==NULL) {
    afmat=niikmat_identity(4,4);
  }
  NIIK_EXIT((!niikmat_convert_from_xfm(afmat,pr[0],img[0])),fcname,"niikmat_convert_from_xfm",9);
  for(n=0; n<4; n++) {
    if(verbose>=3) fprintf(stdout,"[%s]   transformation %i\n",fcname,n+1);
    NIIK_EXIT((!niik_image_affine_transform_3d_update(pr[n],img[0],afmat,NIIK_INTERP_LINEAR)),fcname,
              "niik_image_affine_transform_3d_update",9);
  }
  if(verbose>2) {
    niik_image_write("tmp_pr_brain.nii.gz",pr[0]);
    niik_image_write("tmp_pr_gray.nii.gz",pr[1]);
    niik_image_write("tmp_pr_white.nii.gz",pr[2]);
    niik_image_write("tmp_pr_csf.nii.gz",pr[3]);
  }


  /*********************************************
   * prepare / run classification
   *********************************************/
  NIIK_EXIT(((c=(niik_classifier2 *)calloc(1,sizeof(niik_classifier2 *)))==NULL),fcname,"niik_classifier2 calloc",9);
  c->img=img;
  c->brainmask=brainmask;
  c->pr=pr;

  if(verbose>=1) fprintf(stdout,"[%s] run classification\n",fcname);
  NIIK_EXIT((!niik_classifier2_execute(c)),fcname,"niik_classifier2_update",9);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);


  /*********************************************
   * clean up
   *********************************************/
  /*  for(n=0;n<3;n++){
    img[n]=niik_image_free(img[n]); }
  free(img);
  brainmask=niik_image_free(brainmask);
  for(n=0;n<4;n++){
    pr[n]=niik_image_free(pr[n]); }
  free(pr);

  niik_fc_display(fcname,0);
  */
  return 0;
} /* niik_classify */


