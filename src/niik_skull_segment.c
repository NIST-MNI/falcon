/* FILENAME:     niik_skull_segment.c
 * DESCRIPTION:  Kunio's nifti1 skull-based segmentation program
 * AUTHOR:       Kunio Nakamura
 * DATE:         July 4, 2013
 */

#include "falcon.h"
#include "falcon_morph.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define MAJOR_VERSION (0)
#define MINOR_VERSION (1)
#define MICRO_VERSION (0)

static char *prog_skull_segment_version[] = {
  "  niik_skull_segment history\n"
  "  0.1.0  December 19, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -completely updated\n"
  "  -multispectral version\n"
  "  0.0.0  July 4, 2013, Kunio Nakamura, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
};

static char *prog_skull_segment_help[] = {
  "  niik_skull_registration:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
  "  -outimg-prefx <prefix>            : output image prefix for 3d probabilities\n"
};

typedef struct {          /* class for Kunio's segmentation */
  nifti_image *t1img;       /* input image */
  nifti_image *pdimg;       /* input image */
  nifti_image *t2img;       /* input image */
  nifti_image *gdimg;       /* input image */
  nifti_image *flimg;       /* input image */
  nifti_image *mask;        /* input mask */
  nifti_image *pcsf;        /* anatomic prior */
  nifti_image *pgm;         /* anatomic prior */
  nifti_image *pwm;         /* anatomic prior */
  nifti_image *ples;        /* anatomic prior */
  nifti_image *pskull;      /* anatomic prior */
  nifti_image *pves;        /* anatomic prior */
  nifti_image *pfat;        /* anatomic prior */
  nifti_image *picc;        /* anatomic prior */
} niik_segment; /* class for segmentation */

niik_segment *niik_segment_init() {
  niik_segment *ns=NULL;
  ns->t1img=NULL;
  ns->pdimg=NULL;
  ns->t2img=NULL;
  ns->gdimg=NULL;
  ns->flimg=NULL;
  ns->mask=NULL;
  ns->pcsf=NULL;
  ns->pgm=NULL;
  ns->pwm=NULL;
  ns->ples=NULL;
  ns->pskull=NULL;
  ns->pves=NULL;
  ns->pfat=NULL;
  ns->picc=NULL;
  return ns;
} /* constructor for niik_segment */

niik_segment *niik_segment_free(niik_segment *ns) {
  if(ns==NULL) return NULL;
  ns->t1img=niik_image_free(ns->t1img);
  ns->pdimg=niik_image_free(ns->pdimg);
  ns->t2img=niik_image_free(ns->t2img);
  ns->gdimg=niik_image_free(ns->gdimg);
  ns->flimg=niik_image_free(ns->flimg);
  ns->mask=niik_image_free(ns->mask);
  ns->pcsf=niik_image_free(ns->pcsf);
  ns->pgm=niik_image_free(ns->pgm);
  ns->pwm=niik_image_free(ns->pwm);
  ns->ples=niik_image_free(ns->ples);
  ns->pskull=niik_image_free(ns->pskull);
  ns->pves=niik_image_free(ns->pves);
  ns->pfat=niik_image_free(ns->pfat);
  ns->picc=niik_image_free(ns->picc);
  free(ns);
  return NULL;
} /* constructor for niik_segment */


void usage() {
  fprintf(stdout,"niik_skull_segment\n");
  fprintf(stdout,"  usage: [options] <img.nii> <anat.nii> <out.nii>\n\n");
  fprintf(stdout,"\n");
}

char *niik_skull_segment_tissue_string(int i);
nifti_image *niik_head_segment(nifti_image *img);
nifti_image *niik_skull_segment(nifti_image *img,nifti_image *anat);
void niik_skull_segment_disp_stats(niikmat *mu,niikmat *sd,int nimg);

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *tmpimg=NULL,
   *tmpimgs[9],
   *label=NULL,
    *out=NULL,
     *anat=NULL,
      *img=NULL;
  char
  fname[4096],
        *outimgprefix=NULL,
         *outcheckname=NULL,
          *outlabelname=NULL,
           fcname[20]="niik_skull_segment";
  int
  maxiter=5,
  num_class=4,
  verbose=2,
  m,n,k,n3,
  nc=1,sc=1;
  double
  v,
  dlist[12],
  aUniP=0.5,
  aFWHM=2.5;

  if(argc==1) {
    usage();
    exit(0);
  }

  niik_version_display(fcname,MAJOR_VERSION,MINOR_VERSION,MICRO_VERSION);
  niik_fc_display(fcname,1);



  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s",*prog_skull_segment_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s",*prog_skull_segment_help);
        exit(0);
      }

      else if(!strncmp(argv[nc],"-outimg-prefix",14)) {
        outimgprefix=argv[++nc];
        fprintf(stdout,"[%s] out image prefix %s\n",fcname,outimgprefix);
      } /* outimg prefix */

      else if(!strncmp(argv[nc],"-outcheck",9)) {
        outcheckname=argv[++nc];
        fprintf(stdout,"  outcheck %s\n",outcheckname);
      } /* outcheck */

      else if(!strncmp(argv[nc],"-outlabel",9)) {
        outlabelname=argv[++nc];
        fprintf(stdout,"[%s] out label name %s\n",fcname,outlabelname);
      } /* outlabel */

      else if(!strncmp(argv[nc],"-nclass",7)) {
        num_class=atoi(argv[++nc]);
      } /* class */

      else if(!strncmp(argv[nc],"-aFWHM",6)) {
        aFWHM=atof(argv[++nc]);
      } /* FWHM */

      else if(!strncmp(argv[nc],"-aUniP",6)) {
        aUniP=atof(argv[++nc]);
      } /* FWHM */

      else if(!strncmp(argv[nc],"-iter",5)) {
        maxiter=atoi(argv[++nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s",*prog_skull_segment_help);
        exit(0);
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

  NIIK_EXIT((argc<4),fcname,"too few argments",9);
  NIIK_EXIT((argc>4),fcname,"too many argments",9);

  fprintf(stdout,"[%s] reading image    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] img dim    [%3i,%3i,%3i,%3i]\n",fcname,img->nx,img->ny,img->nz,img->nu);
  fprintf(stdout,"[%s] reading anat img %s\n",fcname,argv[2]);
  NIIK_EXIT(((anat=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);


  fprintf(stdout,"[%s] anat processing   blur %0.3f  UniformMin Prob %0.3f\n",fcname,aFWHM,aUniP);
  n3=anat->nx*anat->ny*anat->nz;
  NIIK_EXIT((!niik_image_filter_gaussian_update(anat,15,aFWHM)),fcname,"niik_image_filter_gaussian_update",9);
  fprintf(stdout,"[%s] anat dim   [%3i,%3i,%3i,%3i]\n",fcname,anat->nx,anat->ny,anat->nz,anat->nu);
  for(n=0; n<n3*3; n++) { /* CSF,GM,WM */
    niik_image_set_voxel(anat,n,niik_image_get_voxel(anat,n)*(1.0-aUniP)+aUniP);
  }
  for(n=0; n<n3; n++) { /* Vessels */
    niik_image_set_voxel(anat,n+n3*4,niik_image_get_voxel(anat,n+n3*4)*(1.0-aUniP)+aUniP);
  }
  for(n=0; n<n3; n++) { /* apply icc on GM/WM */
    v=niik_image_get_voxel(anat,n+n3*5);
    niik_image_mul_voxel(anat,n+n3,v);
    niik_image_mul_voxel(anat,n+n3*2,v);
  }
  NIIK_EXIT(((tmpimg=niik_image_unmerge3d_extract_one(anat,5))==NULL),fcname,"niik_image_unmerge3d_extract_one",9);
  NIIK_EXIT((!niik_image_morph_3d_radius(tmpimg,NIIK_MORPH_DILATE,1.0)),fcname,"niik_image_morph_3d_radius",9);
  for(n=0; n<n3; n++) { /* add icc and skull */
    niik_image_add_voxel(tmpimg,n,niik_image_get_voxel(anat,n+n3*4));
  }
  NIIK_EXIT((!niik_image_filter_gaussian_update(tmpimg,15,aFWHM)),fcname,"niik_image_filter_gaussian_update",9);
  for(n=0; n<n3; n++) { /* apply icc on CSF and vessels */
    v=niik_image_get_voxel(tmpimg,n); /* v=niik_image_get_voxel(anat,n+n3*5); */
    niik_image_mul_voxel(anat,n,v);
    niik_image_mul_voxel(anat,n+n3*4,v);
  }
  for(n=0; n<n3; n++) {
    for(m=0,v=0; m<5; m++) {
      v+=niik_image_get_voxel(anat,n+m*n3);
    }
    if(v<=1e-5) continue;
    v=1.0/v;
    for(m=0; m<5; m++) {
      niik_image_mul_voxel(anat,n+m*n3,v);
    }
  } /* normalize */



  if(verbose>10) {
    fprintf(stdout,"[%s] writing output image test_anat.nii.gz\n",fcname);
    NIIK_EXIT((!niik_image_write("test_anat.nii.gz",anat)),fcname,"niik_image_write",9);
  }


  /************************************************************************
   *
   * input must be mask, t1p, pdw, t2w,flr,t1c in u-dimension
   *
   ************************************************************************/

  if(verbose>1) fprintf(stdout,"[%s] calling niik_skull_segment\n",fcname);
  NIIK_EXIT(((out=niik_skull_segment(img,anat))==NULL),fcname,"niik_skull_segment",9);

  if(outimgprefix!=NULL) { /* output each probability image */
    /* CSF, GM, WM, Les, Skull, etc, air/bg */
    for(n=0; n<7; n++) {
      NIIK_EXIT(((tmpimg=niik_image_unmerge3d_extract_one(out,n))==NULL),fcname,"niik_image_unmerge3d_extract_one",9);
      if((tmpimg->sform_code = img->sform_code)>0) {
        tmpimg->sto_xyz = img->sto_xyz;
        tmpimg->sto_ijk = img->sto_ijk;
      }
      sprintf(fname,"%s_%s.mnc",outimgprefix,niik_skull_segment_tissue_string(n));
      fprintf(stdout,"[%s] writing %s\n",fcname,fname);
      NIIK_EXIT((!niik_image_write(fname,tmpimg)),fcname,"niik_image_write",9);
      tmpimg=niik_image_free(tmpimg);
    } /* each tissue class */
  } /* out image */


  if(outlabelname!=NULL) {/* label */
    NIIK_RET0(((label=niik_image_init(img->nx,img->ny,img->nz,1,1,1,1,
                                      img->dx,img->dy,img->dz,1,1,1,1,NIFTI_TYPE_UINT8))==NULL),
              fcname,"niik_image_init");
    for(n=0; n<label->nvox; n++) {
      dlist[0]=niik_image_get_voxel(out,n);
      dlist[1]=niik_image_get_voxel(out,n+label->nvox);
      dlist[2]=niik_image_get_voxel(out,n+label->nvox*2);
      dlist[3]=niik_image_get_voxel(out,n+label->nvox*3);
      dlist[4]=niik_image_get_voxel(out,n+label->nvox*4);
      dlist[5]=niik_image_get_voxel(out,n+label->nvox*5);
      if(niik_get_sum_from_double_vector(dlist,6)<0.1) continue; /* low probability */
      m=niik_get_max_index_double_vector(dlist,6);
      if(m>=6) continue; /* background */
      k=1;
      niik_image_set_voxel(label,n,(k<<m));
    }
    if(verbose>0) fprintf(stdout,"[%s] writing output image %s\n",fcname,outlabelname);
    if(img->sform_code>0) {
      label->sform_code=img->sform_code;
      label->sto_xyz=img->sto_xyz;
      label->sto_ijk=img->sto_ijk;
    } else if(img->qform_code>0) {
      label->qform_code=img->qform_code;
      label->qto_xyz=img->qto_xyz;
      label->qto_ijk=img->qto_ijk;
    }
    NIIK_EXIT((!niik_image_write(outlabelname,label)),fcname,"niik_image_write",9);
    label=niik_image_free(label);
  } /* label */


  /* one type of output:
   * includes input images and output probability maps
   */
  for(n=0; n< out->nvox; n++) {
    niik_image_mul_voxel( out,n,100.0);
  }
  for(n=0; n<anat->nvox; n++) {
    niik_image_mul_voxel(anat,n,100.0);
  }
  tmpimgs[0]=img;
  tmpimgs[1]=out;
  tmpimgs[2]=anat;
  if(verbose>0) fprintf(stdout,"[%s] combining images\n",fcname);
  NIIK_EXIT(((tmpimg=niik_image_combine(tmpimgs,3,5,0))==NULL),fcname,"niik_image_combine",9);
  img=niik_image_free(img);
  anat=niik_image_free(anat);
  if(verbose>0) fprintf(stdout,"[%s] converting image\n",fcname);
  NIIK_EXIT((!niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert",9);
  if(verbose>0) fprintf(stdout,"[%s] writing output image %s\n",fcname,argv[3]);
  NIIK_EXIT((!niik_image_write(argv[3],tmpimg)),fcname,"niik_image_write",9);
  tmpimg=niik_image_free(tmpimg);

  niik_fc_display(fcname,0);
  exit(0);
} /* niik_skull_segment */


void niik_skull_segment_disp_stats(niikmat *mu,niikmat *sd,int nimg) {
  fprintf(stdout,"  mean / stdv matrix\n");
  fprintf(stdout,"          mask    t1p     pdw      t2w      flr      t1c\n");
  fprintf(stdout,"  CSF:  ");
  niik_display_double_vector(mu->m[0],nimg);
  fprintf(stdout,"   GM:  ");
  niik_display_double_vector(mu->m[1],nimg);
  fprintf(stdout,"   WM:  ");
  niik_display_double_vector(mu->m[2],nimg);
  fprintf(stdout,"  LES:  ");
  niik_display_double_vector(mu->m[3],nimg);
  fprintf(stdout,"   SK:  ");
  niik_display_double_vector(mu->m[4],nimg);
  fprintf(stdout,"   VS:  ");
  niik_display_double_vector(mu->m[5],nimg);
  fprintf(stdout,"   BG:  ");
  niik_display_double_vector(mu->m[6],nimg);
  fprintf(stdout,"  CSF:  ");
  niik_display_double_vector(sd->m[0],nimg);
  fprintf(stdout,"   GM:  ");
  niik_display_double_vector(sd->m[1],nimg);
  fprintf(stdout,"   WM:  ");
  niik_display_double_vector(sd->m[2],nimg);
  fprintf(stdout,"  LES:  ");
  niik_display_double_vector(sd->m[3],nimg);
  fprintf(stdout,"   SK:  ");
  niik_display_double_vector(sd->m[4],nimg);
  fprintf(stdout,"   VS:  ");
  niik_display_double_vector(sd->m[5],nimg);
  fprintf(stdout,"   BG:  ");
  niik_display_double_vector(sd->m[6],nimg);
  return;
}


nifti_image *niik_head_segment(nifti_image *img)
/* roughly segment head
 * -by thresholding all modalities (t1p,pd2,t2w,flr)
 *   -threshold is hard coded
 * -by morphologic operations
 *   -opening to remove noise
 *   -dilationg to diconnect outside to inside
 *   -seed fill from image edge
 *   -erode to get back to the similar size
 */
{
  char fcname[32]="niik_head_segment";
  int verbose=2;
  nifti_image
  *out=NULL;
  double
  thresh=25.0,
  radius=5.0;
  int
  i;
  unsigned char
  *mask=NULL;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0(((out=niik_image_init(img->nx,img->ny,img->nz,1,1,1,1,
                                  img->dx,img->dy,img->dz,1,1,1,1,NIFTI_TYPE_UINT8))==NULL),
            fcname,"niik_image_init");
  mask=(unsigned char *)out->data;
  for(i=0; i<out->nvox; i++) {
    mask[i]=0;
    if(niik_image_get_voxel(img,i+out->nvox  )>thresh) {
      mask[i]=1;
      continue;
    }
    if(niik_image_get_voxel(img,i+out->nvox*2)>thresh) {
      mask[i]=1;
      continue;
    }
    if(niik_image_get_voxel(img,i+out->nvox*3)>thresh) {
      mask[i]=1;
      continue;
    }
    if(niik_image_get_voxel(img,i+out->nvox*4)>thresh) {
      mask[i]=1;
      continue;
    }
  }
  NIIK_RET0((!niik_image_morph_3d_radius(out,NIIK_MORPH_OPEN,1.05)),fcname,"niik_image_morph_3d_radius <open>");
  NIIK_RET0((!niik_image_morph_3d_radius(out,NIIK_MORPH_DILATE,radius)),fcname,"niik_image_morph_3d_radius <dilate>");
  NIIK_RET0((!niik_image_flip_bgfg(out)),fcname,"niik_image_flip_bgfg");
  NIIK_RET0((!niik_image_seed_fill_edge(out,0)),fcname,"niik_image_seed_fill_edge");
  NIIK_RET0((!niik_image_flip_bgfg(out)),fcname,"niik_image_flip_bgfg");
  NIIK_RET0((!niik_image_morph_3d_radius(out,NIIK_MORPH_ERODE,radius)),fcname,"niik_image_morph_3d_radius <erode>");
  NIIK_RET0((!niik_image_seed_fill_from_middle(out,0)),fcname,"niik_image_seed_fill_from_middle");
  if(verbose>0) niik_fc_display(fcname,0);
  return out;
} /* niik_head_segment(nifti_image *img) */


char *niik_skull_segment_tissue_string(int i) {
  switch(i) {
  case 0:
    return "CSF";
  case 1:
    return "GM";
  case 2:
    return "WM";
  case 3:
    return "Lesion";
  case 4:
    return "Skull";
  case 5:
    return "Vessel";
  case 6:
    return "Bkgd";
  default:
    return "Unknown";
  }
  return "Unknown";
}

char *niik_skull_segment_image_string(int i) {
  switch(i) {
  case 0:
    return "Mask";
  case 1:
    return "T1";
  case 2:
    return "PD";
  case 3:
    return "T2";
  case 4:
    return "FLR";
  case 5:
    return "T1C";
  default:
    return "Unknown";
  }
  return "Unknown";
}


nifti_image *niik_skull_segment(nifti_image *img,nifti_image *anat)
/* -img needs dim in u
 *    mask,t1p,pdw,t2w,flr
 * -anat has dim in u
 *    csf,gm,wm,skull
 *
 */
{
  char fcname[32]="niik_skull_segment";
  int verbose=2;
  nifti_image
  *out=NULL;
  niikmat *mu=NULL,*sd=NULL,*pfunc=NULL;
  niikvec *p=NULL,*v=NULL;
  int
  nclass=7,  /* 0=CSF, 1=GM, 2=WM, 3=LES, 4=SKULL, 5=VESSELS, 6=BG */
  nimg=-1;
  int
  check_voxel=(42+55*193+77*193*229),
  i,   /*vox idx*/
  n,m, /*class,img idx*/
  n3,  /*#voxel in 3d*/
  ni,  /*img idx*/
  nmsk=0,nt1p=1,npdw=2,nt2w=3,nflr=4,nt1c=5,
  nc=0,
  ng=1,
  nw=2,
  nl=3,
  ns=4,
  nv=5,
  nb=6;
  unsigned char
  *mask=NULL;
  double
  dv,sv,st,sx;
  int
  iter,maxiter=1;

  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim<=4),fcname,"img is missing u-dim");

  nimg=img->nu;
  NIIK_RET0((nimg<5),fcname,"img needs (mask,t1p,pdw,t2w,flr(,t1c))");
  NIIK_RET0((nimg>6),fcname,"img needs (mask,t1p,pdw,t2w,flr(,t1c))");

  fprintf(stdout,"[%s] #class = %i\n",fcname,nclass);
  fprintf(stdout,"[%s] #img   = %i\n",fcname,nimg);
  fprintf(stdout,"[%s] initialization\n",fcname);
  NIIK_RET0(((out=niik_image_init(img->nx,img->ny,img->nz,1,nclass,1,1,
                                  img->dx,img->dy,img->dz,1,1,1,1,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_image_init");
  n3=img->nx*img->ny*img->nz;


  mu=niikmat_init(nclass,nimg);
  sd=niikmat_init(nclass,nimg);
  pfunc=niikmat_init(nclass,nimg);
  p=niikvec_init(nclass);
  v=niikvec_init(nclass);
  NIIK_RET0(((mask=(unsigned char *)calloc(n3,sizeof(char)))==NULL),fcname,"calloc for mask");
  for(i=0; i<n3; i++) {
    mask[i]=(niik_image_get_voxel(img,i)>0);
  }

  ni=nmsk;
  mu->m[nc][ni]=  1;
  mu->m[ng][ni]= 1;
  mu->m[nw][ni]= 1;
  mu->m[nl][ni]=  1;
  mu->m[nb][ni]= 0;
  mu->m[ns][ni]= 0;
  mu->m[nv][ni]= 1;  /*mask*/
  ni=nt1p;
  mu->m[nc][ni]= 30;
  mu->m[ng][ni]=70;
  mu->m[nw][ni]=90;
  mu->m[nl][ni]= 75;
  mu->m[nb][ni]=-9;
  mu->m[ns][ni]=25;
  mu->m[nv][ni]= 50; /*t1p*/
  ni=npdw;
  mu->m[nc][ni]= 70;
  mu->m[ng][ni]=85;
  mu->m[nw][ni]=65;
  mu->m[nl][ni]= 80;
  mu->m[nb][ni]=-9;
  mu->m[ns][ni]=10;
  mu->m[nv][ni]= 30; /*pdw*/
  ni=nt2w;
  mu->m[nc][ni]=120;
  mu->m[ng][ni]=75;
  mu->m[nw][ni]=70;
  mu->m[nl][ni]= 90;
  mu->m[nb][ni]=-9;
  mu->m[ns][ni]=15;
  mu->m[nv][ni]= 10; /*t2w*/
  ni=nflr;
  mu->m[nc][ni]= 40;
  mu->m[ng][ni]=85;
  mu->m[nw][ni]=75;
  mu->m[nl][ni]=100;
  mu->m[nb][ni]=-9;
  mu->m[ns][ni]= 5;
  mu->m[nv][ni]= 10; /*flr*/
  ni=nt1c;
  mu->m[nc][ni]= 20;
  mu->m[ng][ni]=60;
  mu->m[nw][ni]=85;
  mu->m[nl][ni]=65;
  mu->m[nb][ni]=-5;
  mu->m[ns][ni]= 10;
  mu->m[nv][ni]=100; /*t1c*/

  ni=nmsk;
  sd->m[nc][ni]=  0;
  sd->m[ng][ni]= 0;
  sd->m[nw][ni]= 0;
  sd->m[nl][ni]= 0;
  sd->m[nb][ni]=  0;
  sd->m[ns][ni]= 0;
  sd->m[nv][ni]= 0; /*mask*/
  ni=nt1p;
  sd->m[nc][ni]= 35;
  sd->m[ng][ni]=25;
  sd->m[nw][ni]=10;
  sd->m[nl][ni]=25;
  sd->m[nb][ni]= 10;
  sd->m[ns][ni]=25;
  sd->m[nv][ni]=20; /*t1p*/
  ni=npdw;
  sd->m[nc][ni]= 35;
  sd->m[ng][ni]=25;
  sd->m[nw][ni]=10;
  sd->m[nl][ni]=10;
  sd->m[nb][ni]= 10;
  sd->m[ns][ni]=15;
  sd->m[nv][ni]=20; /*pdw*/
  ni=nt2w;
  sd->m[nc][ni]= 30;
  sd->m[ng][ni]=25;
  sd->m[nw][ni]=15;
  sd->m[nl][ni]=10;
  sd->m[nb][ni]= 10;
  sd->m[ns][ni]=25;
  sd->m[nv][ni]=20; /*t2w*/
  ni=nflr;
  sd->m[nc][ni]= 40;
  sd->m[ng][ni]=20;
  sd->m[nw][ni]=20;
  sd->m[nl][ni]=20;
  sd->m[nb][ni]= 10;
  sd->m[ns][ni]=15;
  sd->m[nv][ni]=20; /*flr*/
  ni=nt1c;
  sd->m[nc][ni]= 35;
  sd->m[ng][ni]=20;
  sd->m[nw][ni]=20;
  sd->m[nl][ni]=20;
  sd->m[nb][ni]=10;
  sd->m[ns][ni]= 15;
  mu->m[nv][ni]= 35; /*t1c*/

  /* pfunc code: 0 = ignore, 1 = gaussian, 2 = heaviside, negative = (1.0-prob) */
  pfunc->m[nc][nmsk]=0;
  pfunc->m[ng][nmsk]=0;
  pfunc->m[nw][nmsk]=0;
  pfunc->m[nl][nmsk]=0;
  pfunc->m[ns][nmsk]=0;
  pfunc->m[nv][nmsk]=0;
  pfunc->m[nb][nmsk]=0;
  pfunc->m[nc][nt1p]=1;
  pfunc->m[ng][nt1p]=1;
  pfunc->m[nw][nt1p]=1;
  pfunc->m[nl][nt1p]=1;
  pfunc->m[ns][nt1p]=1;
  pfunc->m[nv][nt1p]=1;
  pfunc->m[nb][nt1p]=1;
  pfunc->m[nc][npdw]=1;
  pfunc->m[ng][npdw]=1;
  pfunc->m[nw][npdw]=1;
  pfunc->m[nl][npdw]=1;
  pfunc->m[ns][npdw]=1;
  pfunc->m[nv][npdw]=1;
  pfunc->m[nb][npdw]=1;
  pfunc->m[nc][nt2w]=1;
  pfunc->m[ng][nt2w]=1;
  pfunc->m[nw][nt2w]=1;
  pfunc->m[nl][nt2w]=2;
  pfunc->m[ns][nt2w]=1;
  pfunc->m[nv][nt2w]=1;
  pfunc->m[nb][nt2w]=1;
  pfunc->m[nc][nflr]=-2;
  pfunc->m[ng][nflr]=1;
  pfunc->m[nw][nflr]=1;
  pfunc->m[nl][nflr]=2;
  pfunc->m[ns][nflr]=1;
  pfunc->m[nv][nflr]=1;
  pfunc->m[nb][nflr]=1;
  pfunc->m[nc][nt1c]=1;
  pfunc->m[ng][nt1c]=1;
  pfunc->m[nw][nt1c]=1;
  pfunc->m[nl][nt1c]=1;
  pfunc->m[ns][nt1c]=1;
  pfunc->m[nv][nt1c]=2;
  pfunc->m[nb][nt1c]=1;


  p->v[nc]=0.05;
  p->v[ng]=0.50;
  p->v[nw]=0.40;
  p->v[nl]=0.05;
  p->v[nb]=0.001;
  p->v[ns]=0.001;
  niikvec_set_all(p,1);

  niik_skull_segment_disp_stats(mu,sd,nimg);

  for(iter=0; iter<maxiter; iter++) {

    /* update voxel-wise probabilities */
    fprintf(stdout,"[ x   y   z  cl im] prob           gray\n");
    for(i=0; i<n3; i++) { /*voxel*/
      NIIK_RET0((niikvec_set_all(v,1)),fcname,"niikvec_set_all");
      for(n=0; n<nclass; n++) { /*class*/
        for(m=1; m<nimg; m++) { /*img*/

          /* probability calculation */
          switch((int)fabs(pfunc->m[n][m])) {
          case 0:
            break;
          case 1:
            dv=NIIK_GaussPDF(niik_image_get_voxel(img,m*n3+i)-mu->m[n][m],sd->m[n][m]);
            if(pfunc->m[n][m]<0) {
              dv=1.0-dv;
            }
            break;
          case 2:
            dv=NIIK_Heaviside(niik_image_get_voxel(img,m*n3+i)-
                              mu->m[nl][m],sd->m[nl][m]);
            if(pfunc->m[n][m]<0) {
              dv=1.0-dv;
            }
            break;
          default:
            break;
          }
          v->v[n]*=dv;

          if(i==check_voxel) {
            fprintf(stdout,"[%3i %3i %3i   %i %i] %1.9f   %3i   %s %s\n",
                    i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,
                    n,m,dv,
                    (int)niik_image_get_voxel(img,m*n3+i),
                    niik_skull_segment_tissue_string(n),
                    niik_skull_segment_image_string(m));
          }

          /*anatomy : csf,gm,wm*/
          if(anat!=NULL) {
            if     (n < 3) v->v[n]*=niik_image_get_voxel(anat, n*n3+i); /* csf,gm,wm */
            else if(n==nl) v->v[n]*=niik_image_get_voxel(anat,nw*n3+i); /* lesion */
            else if(n==ns) v->v[n]*=niik_image_get_voxel(anat, 3*n3+i); /* skull */
            else if(n==nv) v->v[n]*=niik_image_get_voxel(anat, 4*n3+i); /* vessels */
            else if(n==nb) v->v[n]*=(1.0-niik_image_get_voxel(anat,3*n3+i)); /* not likely to be skull */
          }
        } /*img*/

        if(i==check_voxel) {
          if(anat!=NULL) {
            if     (n < 3) fprintf(stdout,"[%3i %3i %3i   %i a] %1.9f\n",
                                     i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,niik_image_get_voxel(anat,n *n3+i));
            else if(n==nl) fprintf(stdout,"[%3i %3i %3i   %i a] %1.9f\n",
                                     i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,niik_image_get_voxel(anat,nw*n3+i));
            else if(n==ns) fprintf(stdout,"[%3i %3i %3i   %i a] %1.9f\n",
                                     i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,niik_image_get_voxel(anat, 3*n3+i));
            else if(n==nv) fprintf(stdout,"[%3i %3i %3i   %i a] %1.9f\n",
                                     i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,niik_image_get_voxel(anat, 4*n3+i));
            else if(n==nb) fprintf(stdout,"[%3i %3i %3i   %i a] %1.9f\n",
                                     i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,1.0-niik_image_get_voxel(anat, 3*n3+i));
          }
          fprintf(stdout,"[%3i %3i %3i   %i  ] %1.9f\n",
                  i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,n,v->v[n]);
        }
      } /*class*/
      sv=niik_get_sum_from_double_vector(v->v,nclass);
      if(i==check_voxel) {
        fprintf(stdout,"\n[%3i %3i %3i      ] %1.9f (sum)\n",
                i%img->nx,(i/img->nx)%img->ny,(i%n3)/img->nx/img->ny,sv);
        fprintf(stdout,"      CSF          GRAY         WHITE        LES          SKULL        VESSEL       BG     \n");
        niikvec_display(v);
      }
      if(sv>1e-6) {
        for(n=0; n<nclass; n++) { /*class*/
          v->v[n]/=sv;
          niik_image_set_voxel(out,n*n3+i,v->v[n]);
        } /*class*/
      } else {
        for(n=0; n<nclass; n++) { /*class*/
          niik_image_set_voxel(out,n*n3+i,0);
        } /*class*/
      } /* can't tell for sure */

      if(i==check_voxel) {
        niikvec_display(v);
      }

    } /*voxel*/

    /* update group stats */
    for(n=0; n<nclass; n++) { /*class*/
      for(m=1; m<nimg; m++) { /*img*/
        mu->m[n][m]=sd->m[n][m]=0;
        st=sx=0;
        for(i=0; i<n3; i++) { /*voxel*/
          /*if(mask[i]==0) continue; */
          sv=niik_image_get_voxel(out,n*n3+i);
          sx+=sv*niik_image_get_voxel(img,m*n3+i);
          st+=sv;
        } /*vox*/
        mu->m[n][m]=sx/st;
        sx=0;
        for(i=0; i<n3; i++) { /*voxel*/
          /*if(mask[i]==0) continue;*/
          sv=niik_image_get_voxel(out,n*n3+i);
          sx+=sv*NIIK_SQ(niik_image_get_voxel(img,m*n3+i)-mu->m[n][m]);
        } /*vox*/
        sd->m[n][m]=sqrt(sx/st);
      } /*img*/

      for(i=0,p->v[n]=0; i<n3; i++) { /*voxel*/
        p->v[n]+=niik_image_get_voxel(out,n*n3+i);
      }
    } /*class*/
    for(n=0,sv=0; n<nclass; n++) {
      sv+=p->v[n];
    }
    for(n=0; n<nclass; n++) {
      p->v[n]/=sv;
    }

    niik_skull_segment_disp_stats(mu,sd,nimg);
  } /* iter */


  /* write skull probability? */
  /*headimg->sform_code=img->sform_code;
    headimg->sto_xyz=img->sto_xyz;
    headimg->sto_ijk=img->sto_ijk;
    fprintf(stdout,"[%s] writing output image test_skullp.nii.gz\n",fcname);
    NIIK_EXIT((!niik_image_write("test_skullp.nii.gz",headimg)),fcname,"niik_image_write",9);
  */

  mu=niikmat_free(mu);
  sd=niikmat_free(sd);
  pfunc=niikmat_free(pfunc);
  p=niikvec_free(p);
  v=niikvec_free(v);
  free(mask);
  mask=NULL;
  if(verbose>0) niik_fc_display(fcname,0);
  return out;
} /* niik_skull_segment */


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/