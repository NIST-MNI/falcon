/* FILENAME:     niik_inpaint.c
 * DESCRIPTION:  Kunio's nifti1 MS lesion inpainting program
 * AUTHOR:       Kunio Nakamura
 * DATE:         July 30, 2013
 */

#include "falcon.h"

#define MAJOR_VERSION (0)
#define MINOR_VERSION (1)
#define MICRO_VERSION (0)

#define INPAINT_METHOD_MEAN_ONLY 0
#define INPAINT_METHOD_MEAN_STD  1
#define INPAINT_METHOD_F1        2

int niik_image_inpaint(nifti_image *img,nifti_image *maskimg,nifti_image *wmimg,int inpaint_method);


static char *prog_version[] = {
  "  niik_inpaint history\n"
  "  0.0.0 July 30, 2013, Kunio Nakamura, knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
  "\n  0.1.0 July 30, 2013, Kunio Nakamura <knakamura@mrs.bic.mcgill.ca>\n"
  "  -fixed memory issue\n"
  "\n"
};

static char *prog_describe[] = {
  "  [niik_inpaint] description\n"
};

static char *prog_help[] = {
  "  niik_inpaint:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
};

void usage() {
  fprintf(stdout,"niik_inpaint\n");
  fprintf(stdout,"  usage: [options] <in.mnc> <white_matter_mask.mnc> <lesion_mask.mnc> <output.mnc>\n\n");
}




int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *wmimg=NULL,
   *maskimg=NULL,
    *tmpimg[3],
    *img=NULL;
  int
  nc,sc;
  char
  *outcheck=NULL,
   fcname[32]="niik_inpaint";
  int
  inpaint_method=INPAINT_METHOD_F1;
  double
  wmmask_thresh=0.5;


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

      else if(!strncmp(argv[nc],"-meanstd",8)) {
        inpaint_method=INPAINT_METHOD_MEAN_STD;
      }

      else if(!strncmp(argv[nc],"-wmthresh=",10)) {
        wmmask_thresh = atof(argv[nc]+10);
      }

      else if(!strncmp(argv[nc],"-outcheck=",10)) {
        outcheck = argv[nc]+10;
      }

      else if(!strncmp(argv[nc],"-mean",5)) {
        inpaint_method=INPAINT_METHOD_MEAN_ONLY;
      }

      else if(!strncmp(argv[nc],"-f1",2)) {
        inpaint_method=INPAINT_METHOD_F1;
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

  NIIK_EXIT((argc>5),fcname,"too many arguments",9);
  NIIK_EXIT((argc<5),fcname,"too few arguments",9);

  niik_fc_display(fcname,1);

  fprintf(stdout,"[%s] reading input image:    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img=niik_image_read(argv[1]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img,NIFTI_TYPE_FLOAT32,1)),
            fcname,
            "niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading input mask:     %s\n",fcname,argv[2]);
  NIIK_EXIT(((wmimg=niik_image_read(argv[2]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_threshold(wmimg,wmmask_thresh)),
            fcname,
            "niik_image_threshold",9);

  fprintf(stdout,"[%s] reading input image:    %s\n",fcname,argv[3]);
  NIIK_EXIT(((maskimg=niik_image_read(argv[3]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(maskimg,NIFTI_TYPE_UINT8,1)),
            fcname,
            "niik_image_type_convert",9);

  NIIK_EXIT(((tmpimg[0]=niik_image_copy(img))==NULL),fcname,"niik_image_copy",9);
  NIIK_EXIT((!niik_image_inpaint(img,maskimg,wmimg,inpaint_method)),fcname,"niik_image_inpaint",9);
  tmpimg[1]=img;

  fprintf(stdout,"[%s] writing output image    %s\n",fcname,argv[4]);
  NIIK_EXIT((!niik_image_write(argv[4],img)),fcname,"niik_image_write",9);

  if(outcheck!=NULL)  {
    fprintf(stdout,"[%s] writing output image    %s\n",fcname,outcheck);
    NIIK_EXIT((!niik_image_combine_and_write_as_type(outcheck,tmpimg,2,'t',140,NIFTI_TYPE_UINT8)),fcname,"niik_image_combine_and_write_as_type",9);
  }

  img=niik_image_free(img);
  wmimg=niik_image_free(wmimg);
  maskimg=niik_image_free(maskimg);
  tmpimg[0]=niik_image_free(tmpimg[0]);
  exit(0);
} /* main */


int niik_image_get_neighbors(nifti_image *img,nifti_image *mask,int *ijk,int *nei,int *nnei,double radius) {
  unsigned char
  *bimg=NULL;
  int
  mi,mj,mk,ni,nj,nk,
  i,j,k,n1,n2,n3;
  double r1,r2,r3;
  char fcname[32]="niik_image_get_neighbors";
  NIIK_RET0((img==NULL),fcname,"img is null");
  if(mask!=NULL) {
    NIIK_RET0((mask->datatype!=NIFTI_TYPE_UINT8),fcname,"mask is not uint8");
  }
  mi=ijk[0]-(radius/img->dx)-1;
  mj=ijk[1]-(radius/img->dy)-1;
  mk=ijk[2]-(radius/img->dz)-1;
  ni=ijk[0]+(radius/img->dx)+1;
  nj=ijk[1]+(radius/img->dy)+1;
  nk=ijk[2]+(radius/img->dz)+1;
  mi=(mi<0)?0:mi;
  mj=(mj<0)?0:mj;
  mk=(mk<0)?0:mk;
  ni=(ni>=img->nx)?img->nx-1:ni;
  nj=(nj>=img->ny)?img->ny-1:nj;
  nk=(nk>=img->nz)?img->nz-1:nk;
  radius*=radius;
  *nnei=0;
  if(mask!=NULL) {
    /*NIIK_RET0(((bimg=niik_image_get_voxels_as_uint8_vector(mask))==NULL),fcname,"niik_image_get_voxels_as_uint8_vector");*/
    bimg=mask->data;
    for(k=mk; k<=nk; k++) {
      n3=k*img->nx*img->ny;
      r3=NIIK_SQ((k-ijk[2])*img->dz);
      for(j=mj; j<=nj; j++) {
        n2=j*img->nx;
        r2=NIIK_SQ((j-ijk[1])*img->dy);
        for(i=mi; i<=ni; i++) {
          n1=i+n2+n3;
          if(bimg[n1]==0) continue;
          r1=NIIK_SQ((i-ijk[0])*img->dx);
          if(r1+r2+r3>radius) continue;
          nei[*nnei]=n1;
          *nnei = *nnei + 1;
        }
      }
    }
  }  /* with mask */
  else {
    for(k=mk; k<=nk; k++) {
      n3=k*img->nx*img->ny;
      r3=NIIK_SQ((k-ijk[2])*img->dz);
      for(j=mj; j<=nj; j++) {
        n2=j*img->nx;
        r2=NIIK_SQ((j-ijk[1])*img->dy);
        for(i=mi; i<=ni; i++) {
          n1=i+n2+n3;
          r1=NIIK_SQ((i-ijk[0])*img->dx);
          if(r1+r2+r3>radius) continue;
          nei[*nnei]=n1;
          *nnei = *nnei + 1;
        }
      }
    }
  } /* withing mask */
  return 1;
}



int niik_image_inpaint(nifti_image *img,nifti_image *maskimg,nifti_image *wmimg,int inpaint_method)
/* inpainting function
 * img  t1-weighted input image, to be changed, type = float32
 * maskimg   lesion mask of type = uint8
 * wmimg   white matter mask to estimate the lesions, type = uint8
 * inpaint_method    choose from the list
 */
{
  double
  *gnei,
  trim=0.1,
  radius=4.00,
  wmcount=0,
  wmstdv=0,
  wmmean=0;
  float
  *fimg=NULL;
  unsigned char
  *wm=NULL,
   *tmpmask=NULL,
    *mask=NULL;
  char fcname[32]="niik_image_inpaint";
  int
  verbose=0,
  iter=0,maxiter=20,
  mskip,
  area,
  i,ijk[3],*nei,n,nnei=0,ni,nf,
                  nneilim=4;

  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((maskimg==NULL),fcname,"lesion mask is null");
  NIIK_RET0((wmimg==NULL),fcname,"white matter mask is null");
  NIIK_RET0((img->nvox!=maskimg->nvox),fcname,"difference #voxel for img and mask");
  NIIK_RET0((img->nvox!=  wmimg->nvox),fcname,"difference #voxel for img and WM mask");

  fprintf(stdout,"[%s] mean / stdv only\n",fcname);
  fimg=(float         *)    img->data;
  wm  =(unsigned char *)  wmimg->data;
  mask=(unsigned char *)maskimg->data;
  area=img->nx*img->ny;
  NIIK_RET0(((tmpmask=niik_image_get_voxels_as_uint8_vector(maskimg))==NULL),fcname,"niik_image_get_voxels_as_uint8_vector");
  for(i=0; i<img->nvox; i++) {
    if(mask[i] >0) continue;
    if(  wm[i]==0) continue;
    wmmean+=fimg[i];
    wmstdv+=fimg[i]*fimg[i];
    wmcount+=1;
  } /* each voxel */
  wmmean/=wmcount;
  wmstdv=sqrt(wmstdv/wmcount-wmmean*wmmean);
  fprintf(stdout,"[%s] NAWM mean (SD) = %12.6f (%2.6f) %i\n",fcname,wmmean,wmstdv,(int)wmcount);

  switch(inpaint_method) {
  default:
    fprintf(stdout,"[%s] ERROR: unkown option %i\n",fcname,inpaint_method);
    exit(9);
  case INPAINT_METHOD_MEAN_ONLY:
    for(i=0; i<img->nvox; i++) {
      if(mask[i]==0) continue;
      fimg[i]=wmmean;
    } /* add mean +/- stdv */
    break;
  case INPAINT_METHOD_MEAN_STD:
    for(i=0; i<img->nvox; i++) {
      if(mask[i]==0) continue;
      fimg[i]=wmmean+niik_get_rand_normal()*wmstdv;
    } /* add mean +/- stdv */
    break;
  case INPAINT_METHOD_F1:  /* testing ... */
    NIIK_RET0((( nei=(int    *)calloc(img->nvox,sizeof(int)))==NULL),fcname,"calloc for nei");
    NIIK_RET0(((gnei=(double *)calloc(img->nvox,sizeof(double)))==NULL),fcname,"calloc for gnei");

    for(i=0; i<img->nvox; i++) {
      if(wm[i]>0 && mask[i]>0) wm[i]=0;
    }

    for(;;) {
      fprintf(stdout,"[%s] #voxel to fill = %i\n",fcname,niik_image_count_mask(maskimg));
      for(i=0; i<img->nvox; i++) {
        if(mask[i]==0) continue;
        ijk[0]=i%img->nx;
        ijk[1]=(i/img->nx) % img->ny;
        ijk[2]=(i/img->nx/img->ny);
        /* check if neighbours are connected to mask */
        if(iter<10) {
          mskip=1;
          if(ijk[0]>0) if(mask[i-1]>0) mskip=0;
          if(ijk[1]>0) if(mask[i-img->nx]>0) mskip=0;
          if(ijk[2]>0) if(mask[i-area]>0) mskip=0;
          if(ijk[0]<img->nx-1) if(mask[i+1]>0)       mskip=0;
          if(ijk[1]<img->ny-1) if(mask[i+img->nx]>0) mskip=0;
          if(ijk[2]<img->nz-1) if(mask[i+area]>0)    mskip=0;
        } else
          mskip=0;
        if(mskip) continue;
        NIIK_RET0((!niik_image_get_neighbors(img,wmimg,ijk,nei,&nnei,radius)),fcname,"niik_image_get_neighbors");
        if(iter+1>=maxiter) nneilim=1e9;
        if(nnei<nneilim) continue;
        tmpmask[i]=0;
        wm[i]=1;
        for(n=0; n<nnei; n++) {
          gnei[n]=(double)fimg[nei[n]];
        }
        NIIK_RET0((!niik_sort_double_vector(gnei,nnei)),fcname,"niik_sort_double_vector");
        wmmean=wmstdv=wmcount=0;
        ni=(int)(nnei*trim);
        nf=(int)(nnei*(1.0-trim));
        ni=(ni<0)?0:ni;
        nf=(nf>nnei)?nnei:nf;
        for(n=ni; n<nf; n++) {
          wmmean+=fimg[nei[n]];
          wmstdv+=fimg[nei[n]]*fimg[nei[n]];
          wmcount+=1;
        }
        if(wmcount>=3) {
          wmmean/=wmcount;
          wmstdv=sqrt(wmstdv/wmcount-wmmean*wmmean);
          fimg[i]=wmmean+niik_get_rand_normal()*wmstdv;
        } else
          fimg[i]=wmmean/wmcount;
        if(verbose>0) {
          fprintf(stdout,"[%s] %3i,%3i,%3i    %12.3f %6.2f %4i     %.6f\n",fcname,ijk[0],ijk[1],ijk[2],wmmean,wmstdv,nnei,fimg[i]);
          if(ijk[0]==86) {
            if(ijk[1]==100) {
              if(ijk[2]==61) {
                if(verbose>0) {
                  niik_display_int_vector(nei,nnei);
                  for(n=ni; n<nf; n++) {
                    fprintf(stdout,"%f ",fimg[nei[n]]);
                  }
                  fprintf(stdout,"\n");
                }
              }
            }
          }
        }
      } /* each voxel */
      for(i=0; i<img->nvox; i++) mask[i]=tmpmask[i];
      n=niik_image_count_mask(maskimg);
      if(verbose>0) fprintf(stdout,"[%s] iter %i, stop? %i\n",fcname,iter,n);
      if(iter++>maxiter) break;
      if(!n) break;
    } /* loop until we fill everything */
    free(nei);
    free(gnei);
    break;
  } /* inpaint methods */
  return 1;
} /* niik_image_inpaint */
