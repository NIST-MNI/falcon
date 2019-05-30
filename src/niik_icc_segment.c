/* FILENAME:
 * DESCRIPTION:
 * AUTHOR:       Kunio Nakamura
 * DATE:         August 29, 2014
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  niik_icc_segment history\n"
  "\n"
  "  0.0.0     August 29, 2014, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -initial version\n"
  "\n"
};

static char *prog_help[] = {
  "  niik_icc_segment:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};


typedef struct {
  int verbose;          /* verbose level */
  int debug_keep_tmp;   /* debug flag */
  int feature_method;   /* feature method, only one for now */
  niikmat *fw;          /* feature weighting */
  char *library_dir;    /* library dir */
  char *outname;        /* output filename */
  nifti_image *t1img;   /* t1w input image */
  nifti_image *pdimg;   /* pdw input image */
  nifti_image *t2img;   /* t2w input image */
  nifti_image *gdimg;   /* gdw input image */
  nifti_image *margin_mask;  /* margin mask */
  nifti_image **libmask;     /* libmask[model]  */
  nifti_image ***libfm;      /* libfm[contrast][model] where contrast=[t1w,t2w,pdw] */
  nifti_image *mask;         /* mask output image */
  int num_model;
  int nsearch;
  int ijk[4];
} niik_icc_segment;

niik_icc_segment *niik_icc_segment_init() {
  niik_icc_segment *icc=NULL;
  char fcname[32]="niik_icc_segment_init";
  if((icc=(niik_icc_segment *)malloc(sizeof(niik_icc_segment)))==NULL) {
    fprintf(stderr,"[%s] ERROR: malloc, icc\n",fcname);
    return NULL;
  }
  icc->num_model=0;
  icc->debug_keep_tmp=0;
  icc->library_dir=NULL;
  icc->outname=NULL;
  icc->t1img=icc->t2img=icc->pdimg=icc->gdimg=icc->mask=NULL;
  icc->verbose=1;
  icc->fw=niikmat_init(3,8);
  niikmat_set_all(icc->fw,1);
  icc->fw->m[0][0]=icc->fw->m[1][0]=icc->fw->m[2][0]=2;
  icc->fw->m[0][7]=icc->fw->m[1][7]=icc->fw->m[2][7]=10;
  icc->feature_method = 0;
  icc->margin_mask=NULL;
  icc->libfm=NULL;
  icc->ijk[0]=icc->ijk[1]=icc->ijk[2]=icc->ijk[3]=-1;
  icc->nsearch=2;
  return icc;
} /* constructor */

niik_icc_segment *niik_icc_segment_free(niik_icc_segment *icc) {
  char fcname[32]="niik_icc_segment_free";
  if(icc==NULL) return NULL;
  if(icc->verbose) niik_fc_display(fcname,1);
  if(icc->library_dir!=NULL) free(icc->library_dir);
  if(icc->outname!=NULL) free(icc->outname);
  if(icc->t1img!=NULL) icc->t1img=niik_image_free(icc->t1img);
  if(icc->t2img!=NULL) icc->t2img=niik_image_free(icc->t2img);
  if(icc->pdimg!=NULL) icc->pdimg=niik_image_free(icc->pdimg);
  if(icc->gdimg!=NULL) icc->gdimg=niik_image_free(icc->gdimg);
  if(icc->mask!=NULL) icc->mask=niik_image_free(icc->mask);
  if(icc->fw!=NULL) icc->fw=niikmat_free(icc->fw);
  return NULL;
} /* destructor */


int niik_icc_segment_feature_matching_voxel(niik_icc_segment *icc,int vox,nifti_image **imgfm,unsigned char ***libfm,unsigned char **libmask,nifti_image **img,int *min_model_idx) {
  int
  m,mm,j,
  n,nn,numc=3,
       dim2,dim3,
       ii[3],imin,imax,
       ijk[3];
  double dsum;
  float *fimg[3];
  niikvec *v=NULL;

  v=niikvec_init(icc->num_model);
  ijk[0]=vox%icc->t1img->nx;
  ijk[1]=(vox/icc->t1img->nx)%icc->t1img->ny;
  ijk[2]=(vox/icc->t1img->nx)/icc->t1img->ny;
  dim2=icc->t1img->nx*icc->t1img->ny;
  dim3=dim2*icc->t1img->nz;

  for(mm=0; mm<numc; mm++) {
    fimg[mm]=NULL;
    if(img[mm]==NULL) continue;
    fimg[mm]=(float *)img[mm]->data;
  }

  imin=imax=libmask[0][vox];
  for(m=1; m<icc->num_model; m++) {
    imin=(imin>libmask[m][vox])?libmask[m][vox]:imin;
    imax=(imax<libmask[m][vox])?libmask[m][vox]:imax;
  }
  if(imin==imax) {
    /*fprintf(stdout,"%3i %3i %3i   %i  min/max\n",ijk[0],ijk[1],ijk[2],imin);*/
    return imin;
  }

  for(m=0; m<icc->num_model; m++) {
    v->v[m]=0;
    for(mm=0; mm<numc; mm++) {
      if(img[mm]==NULL) continue;
      dsum=0;
      j=0;
      for(ii[0]=ijk[0]-icc->nsearch; ii[0]<=ijk[0]+icc->nsearch; ii[0]++) {
        if(ii[0]<0) continue;
        else if(ii[0]>=icc->t1img->nx) break;
        for(ii[1]=ijk[1]-icc->nsearch; ii[1]<=ijk[1]+icc->nsearch; ii[1]++) {
          if(ii[1]<0) continue;
          else if(ii[1]>=icc->t1img->ny) break;
          for(ii[2]=ijk[2]-icc->nsearch; ii[2]<=ijk[2]+icc->nsearch; ii[2]++) {
            if(ii[2]<0) continue;
            else if(ii[2]>=icc->t1img->nz) break;
            nn=ii[0]+ii[1]*icc->t1img->nx+ii[2]*dim2;
            j++;
            dsum+=fabs(fimg[mm][nn]-libfm[mm][m][nn]);
          }
        }
      }
      v->v[m]=icc->fw->m[mm][7]*dsum/j;
      for(n=0; n<7; n++) {
        v->v[m]+=fabs(niik_image_get_voxel(imgfm[mm],vox+n*dim3) -
                      niik_image_get_voxel(icc->libfm[mm][m],vox+n*dim3)) *
                 icc->fw->m[mm][n];
      }
    }
  } /* each image */
  m=niik_get_min_index_double_vector(v->v,v->num);

  if(vox==icc->ijk[0])  {
    for(n=0; n<v->num; n++) {
      if(m==n) {
        fprintf(stdout,"%3i *%3i* %12.6f   %s\n",n,libmask[n][vox],v->v[n],icc->libmask[n]->fname);
      } else {
        fprintf(stdout,"%3i  %3i  %12.6f   %s\n",n,libmask[n][vox],v->v[n],icc->libmask[n]->fname);
      }
    }
  }

  *min_model_idx=niik_get_min_index_double_vector(v->v,v->num);
  v=niikvec_free(v);
  return 2;
} /* niik_icc_segment_feature_matching_voxel */

int niik_icc_segment_main(niik_icc_segment *icc) {
  char fcname[32]="niik_icc_segment_main";
  char library_name[512];
  nifti_image
  **imgfm=NULL,
    **img=NULL;
  int
  i,
  m,numc=3,n,nimg=0;
  char *curr_dir=NULL;
  float *fimg[12];
  unsigned char ***libfm=NULL,**libmask=NULL;

  if(icc->verbose) {
    niik_fc_display(fcname,1);
    fprintf(stdout,"[%s] nsearch = %i\n",fcname,icc->nsearch);
  }


  NIIK_RET0((icc==NULL),fcname,"icc is null");
  NIIK_RET0((icc->margin_mask==NULL),fcname,"margin_mask is null");
  NIIK_RET0((icc->t1img==NULL),fcname,"missing t1w img");

  curr_dir = getcwd(curr_dir, 0);

  NIIK_RET0(((icc->libfm=(nifti_image ***)calloc(numc,sizeof(nifti_image **)))==NULL),fcname,"calloc for libfm");
  NIIK_RET0(((libfm=(unsigned char ***)calloc(numc,sizeof(char **)))==NULL),fcname,"calloc for char");
  for(n=0; n<3; n++) libfm[n]=NULL;
  NIIK_RET0((icc->library_dir==NULL),fcname,"missing library_dir");

  fprintf(stdout,"[%s] using lib, %s\n",fcname,icc->library_dir);
  if((chdir(icc->library_dir))!=0) {
    fprintf(stderr,"[%s] ERROR: chdir, %s\n",fcname,icc->library_dir);
    return 0;
  }

  /* reading masks */
  fprintf(stdout,"[%s] reading library masks\n",fcname);
  sprintf(library_name,"%s/library.masks.2mm",icc->library_dir);
  NIIK_RET0(((icc->libmask=niik_image_read_multiple_from_file(library_name,&icc->num_model))==NULL),fcname,"reading libmask");
  libmask=(unsigned char **)calloc(icc->num_model,sizeof(char *));
  for(n=0; n<icc->num_model; n++) {
    NIIK_RET0((!niik_image_type_convert(icc->libmask[n],NIFTI_TYPE_UINT8)),fcname,"niik_image_type_convert");
    fprintf(stdout,"  %5i %s\n",n,icc->libmask[n]->fname);
    libmask[n]=(unsigned char *)icc->libmask[n]->data;
  }

  /* reading library T1 images */
  fprintf(stdout,"[%s] reading library t1w images\n",fcname);
  sprintf(library_name,"%s/library.stx.2mm.fm4",icc->library_dir);
  icc->libfm[0]=niik_image_read_multiple_from_file(library_name,&nimg);
  NIIK_RET0((icc->num_model!=nimg),fcname,"lib num does not match");
  libfm[0]=(unsigned char **)calloc(icc->num_model,sizeof(char *));
  for(n=0; n<icc->num_model; n++) {
    libfm[0][n]=(unsigned char *)icc->libfm[0][n]->data;
  }

  /* reading library T2 images */
  if(icc->t2img!=NULL) {
    fprintf(stdout,"[%s] reading library t2w images\n",fcname);
    sprintf(library_name,"%s/library.stx.t2.2mm.fm4",icc->library_dir);
    icc->libfm[1]=niik_image_read_multiple_from_file(library_name,&nimg);
    NIIK_RET0((icc->num_model!=nimg),fcname,"lib num does not match");
    libfm[1]=(unsigned char **)calloc(icc->num_model,sizeof(char *));
    for(n=0; n<icc->num_model; n++) {
      libfm[1][n]=(unsigned char *)icc->libfm[1][n]->data;
    }
  }

  /* reading library PD images */
  if(icc->pdimg!=NULL) {
    fprintf(stdout,"[%s] reading library pdw images\n",fcname);
    sprintf(library_name,"%s/library.stx.pd.2mm.fm4",icc->library_dir);
    icc->libfm[2]=niik_image_read_multiple_from_file(library_name,&nimg);
    NIIK_RET0((icc->num_model!=nimg),fcname,"lib num does not match");
    libfm[2]=(unsigned char **)calloc(icc->num_model,sizeof(char *));
    for(n=0; n<icc->num_model; n++) {
      libfm[2][n]=(unsigned char *)icc->libfm[2][n]->data;
    }
  }


  NIIK_RET0(((chdir(curr_dir))!=0),fcname,"chdir curr_dir");

  fprintf(stdout,"[%s] image feature map\n",fcname);
  img=(nifti_image **)calloc(numc,sizeof(nifti_image **));
  imgfm=(nifti_image **)calloc(numc,sizeof(nifti_image **));
  img[0]=icc->t1img;
  img[1]=icc->t2img;
  img[2]=icc->pdimg;
  for(n=0; n<numc; n++) {
    if(img[n]==NULL) continue;
    fprintf(stdout,"[%s] feature %i image\n",fcname,n);
    NIIK_RET0((!niik_image_type_convert(img[n],NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
    NIIK_RET0(((imgfm[n]=niik_image_feature_map_type4(img[n]))==NULL),fcname,"niik_image_feature_map_type4");
    fimg[n]=(float *)img[n]->data;
  }
  icc->mask=niik_image_copy_as_type(img[0],NIFTI_TYPE_UINT8);
  niik_image_clear(icc->mask);

  fprintf(stdout,"[%s] margin mask\n",fcname);
  niik_image_resample_3d_update(icc->margin_mask,
                                img[0]->dx,img[0]->dy,img[0]->dz,
                                img[0]->nx,img[0]->ny,img[0]->nz,
                                NIIK_INTERP_NN);


  /* matching features
   * each voxel
   *   each modality
   *     each feature
   */
  fprintf(stdout,"[%s] feature matching\n",fcname);
  if(icc->ijk[0]>0) {
    icc->ijk[0]=icc->ijk[1]+icc->ijk[2]*img[0]->nx+icc->ijk[3]*img[0]->nx*img[0]->ny;
  }

  #pragma omp parallel for private(m,n)
  for(i=0; i<img[0]->nvox; i++) {
    //    if(niik_image_get_voxel(icc->margin_mask,i)<0.5) continue;
    n=niik_icc_segment_feature_matching_voxel(icc,i,imgfm,libfm,libmask,img,&m);
    if(n>1)
      niik_image_set_voxel(icc->mask,i,niik_image_get_voxel(icc->libmask[m],i));
    else
      niik_image_set_voxel(icc->mask,i,n);
  }

  if(icc->verbose) niik_fc_display(fcname,0);
  return 1;
} /* niik_icc_segment_main */


void usage() {
  fprintf(stdout,"niik_inpaint\n");
  fprintf(stdout,"  usage: [options] -lib=<library> -t1img=<in_T1.mnc> -out=<out.mnc>\n\n");
}


int main(int argc,char *argv[],char *envp[]) {
  int
  m,n,
  nc=1;
  char
  **CSlist,
  fcname[32]="niik_icc_segment";
  niik_icc_segment *icc=NULL;

  if(argc==1) {
    usage();
    exit(1);
  }

  niik_fc_display(fcname,1);
  NIIK_EXIT(((icc=niik_icc_segment_init())==NULL),fcname,"niik_icc_segment_init",1);

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(1);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(1);
      }


      else if(!strncmp(argv[nc],"-nsearch=",9)) {
        icc->nsearch=atoi(argv[nc]+9);
      }

      else if(!strncmp(argv[nc],"-margin=",8)) {
        NIIK_EXIT((icc->margin_mask!=NULL),fcname,"margin_mask is already used",1);
        NIIK_EXIT(((icc->margin_mask=niik_image_read(argv[nc]+8))==NULL),fcname,"niik_image_read margin mask",1);
      } else if(!strncmp(argv[nc],"-t2img=",7)) {
        NIIK_EXIT((icc->t2img!=NULL),fcname,"t2img is already used",1);
        NIIK_EXIT(((icc->t2img=niik_image_read(argv[nc]+7))==NULL),fcname,"niik_image_read T2 img",1);
      } else if(!strncmp(argv[nc],"-pdimg=",7)) {
        NIIK_EXIT((icc->pdimg!=NULL),fcname,"pdimg is already used",1);
        NIIK_EXIT(((icc->pdimg=niik_image_read(argv[nc]+7))==NULL),fcname,"niik_image_read PD img",1);
      } else if(!strncmp(argv[nc],"-gdimg=",7)) {
        NIIK_EXIT((icc->gdimg!=NULL),fcname,"gdimg is already used",1);
        NIIK_EXIT(((icc->gdimg=niik_image_read(argv[nc]+7))==NULL),fcname,"niik_image_read Gd img",1);
      } else if(!strncmp(argv[nc],"-t1img=",7)) {
        NIIK_EXIT((icc->t1img!=NULL),fcname,"t1img is already used",1);
        NIIK_EXIT(((icc->t1img=niik_image_read(argv[nc]+7))==NULL),fcname,"niik_image_read T1 img",1);
      }

      else if(!strncmp(argv[nc],"-img=",5)) {
        NIIK_EXIT((icc->t1img!=NULL),fcname,"t1img is already used",1);
        NIIK_EXIT(((icc->t1img=niik_image_read(argv[nc]+5))==NULL),fcname,"niik_image_read T1 img",1);
      } else if(!strncmp(argv[nc],"-out=",5)) {
        NIIK_EXIT((icc->outname!=NULL),fcname,"outname is already used",1);
        if((icc->outname=(char *)calloc(4096,sizeof(char)))==NULL) {
          fprintf(stderr,"[%s] ERROR: calloc, outname\n",fcname);
          free(icc);
          exit(1);
        }
        strcpy(icc->outname,argv[nc]+5);
      }

      else if(!strncmp(argv[nc],"-lib=",5)) {
        NIIK_EXIT((icc->library_dir!=NULL),fcname,"library_dir is already used",1);
        if((icc->library_dir=(char *)calloc(4096,sizeof(char)))==NULL) {
          fprintf(stderr,"[%s] ERROR: calloc, library_dir\n",fcname);
          free(icc);
          exit(1);
        }
        strcpy(icc->library_dir,argv[nc]+5);
      }

      else if(!strncmp(argv[nc],"-ijk=",5)) {
        NIIK_EXIT((icc->ijk[0]>0),fcname,"-ijk is already defined",1);
        icc->ijk[0]=1;
        NIIK_EXIT(((CSlist=niik_csstring_to_list(argv[nc]+5,&m))==NULL),fcname,"niik_csstring_to_list",1);
        NIIK_EXIT((m!=3),fcname,"ijk should have 3 numnbers",1);
        for(n=0; n<m; n++) {
          icc->ijk[n+1]=atoi(CSlist[n]);
        }
        free(CSlist);
      }

      else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(1);
      } else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(1);
      } else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(1);
      }
      nc++;
    } else {
      fprintf(stderr,"[%s] ERROR: unknown argument: %s\n",fcname,argv[nc]);
      exit(1);
    }
  } /* reading options (while) */


  NIIK_EXIT((!niik_icc_segment_main(icc)),fcname,"niik_icc_segment_main",1);

  fprintf(stdout,"[%s] writing output %s\n",fcname,icc->outname);
  NIIK_EXIT((!niik_image_write(icc->outname,icc->mask)),fcname,"niik_image_write",1);

  icc=niik_icc_segment_free(icc);
  niik_fc_display(fcname,0);
  exit(0);
} /* main */
