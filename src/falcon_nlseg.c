/* Filename:     nifti1_kunio_venseg.c
 * Description:
 * Author:       Kunio Nakamura
 * Date:
 */


#ifndef _FALCON_NLSEG_C_
#define _FALCON_NLSEG_C_

#include "falcon.h"

int niik_image_nlseg_extract_3dpatch(void *img,void *patch,int datatype,int x,int y,int z,int size,int sx,int sy,int sz) {
  NIIK_RET0((!niik_image_nlseg_extract_3dpatch_mean_var(img,patch,datatype,x,y,z,size,sx,sy,sz,NULL,NULL)),__func__,"niik_image_nlseg_extract_3dpatch_mean_var");
  return 1;
}

int niik_image_nlseg_extract_3dpatch_mean_var(void *img,void *patch,int datatype,int x,int y,int z,int size,int sx,int sy,int sz,float *M,float *V) {
  char fcname[64]="niik_image_nlseg_extract_3dpatch";
  unsigned char *bimg=NULL,*bout=NULL;
  float *fimg=NULL,*fout=NULL;
  double dm=0,ds=0,*dimg=NULL,*dout=NULL;
  int
  sxy,
  psize,
  ni,nj,nk,
  i,j,k;

  psize=2*size+1;
  sxy=sx*sy;
  switch(datatype) {
  default:
    fprintf(stderr,"[%s] ERROR: datatype is not allowed, %s\n",fcname,nifti_datatype_string(datatype));
    break;
  case NIFTI_TYPE_UINT8:
    bimg=(unsigned char *)img;
    bout=(unsigned char *)patch;
    for(k=-size; k<=size; k++) {
      for(j=-size; j<=size; j++) {
        for(i=-size; i<=size; i++) {
          ni=i+x;
          nj=j+y;
          nk=k+z;
          if(ni<0) ni=0;
          else if(ni>=sx) ni=sx-1;
          if(nj<0) nj=0;
          else if(nj>=sy) nj=sy-1;
          if(nk<0) nk=0;
          else if(nk>=sz) nk=sz-1;
          bout[(i+size) + (j+size)*psize + (k+size)*psize*psize] = bimg[ni+nj*sx+nk*sxy];
          dm+=bimg[ni+nj*sx+nk*sxy];
          dm+=(double)bimg[ni+nj*sx+nk*sxy]*bimg[ni+nj*sx+nk*sxy];
        }
      }
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fimg=(float *)img;
    fout=(float *)patch;
    for(k=-size; k<=size; k++) {
      for(j=-size; j<=size; j++) {
        for(i=-size; i<=size; i++) {
          ni=i+x;
          nj=j+y;
          nk=k+z;
          if(ni<0) ni=0;
          else if(ni>=sx) ni=sx-1;
          if(nj<0) nj=0;
          else if(nj>=sy) nj=sy-1;
          if(nk<0) nk=0;
          else if(nk>=sz) nk=sz-1;
          fout[(i+size) + (j+size)*psize + (k+size)*psize*psize] = fimg[ni+nj*sx+nk*sxy];
          dm+=fimg[ni+nj*sx+nk*sxy];
          dm+=(double)fimg[ni+nj*sx+nk*sxy]*fimg[ni+nj*sx+nk*sxy];
        }
      }
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dimg=(double *)img;
    dout=(double *)patch;
    for(k=-size; k<=size; k++) {
      for(j=-size; j<=size; j++) {
        for(i=-size; i<=size; i++) {
          ni=i+x;
          nj=j+y;
          nk=k+z;
          if(ni<0) ni=0;
          else if(ni>=sx) ni=sx-1;
          if(nj<0) nj=0;
          else if(nj>=sy) nj=sy-1;
          if(nk<0) nk=0;
          else if(nk>=sz) nk=sz-1;
          dout[(i+size) + (j+size)*psize + (k+size)*psize*psize] = dimg[ni+nj*sx+nk*sxy];
          dm+=dimg[ni+nj*sx+nk*sxy];
          dm+=dimg[ni+nj*sx+nk*sxy]*dimg[ni+nj*sx+nk*sxy];
        }
      }
    }
    break;
  } /* switch */
  dm/=pow(size*2+1,3);
  if(M!=NULL) {
    *M=dm;
    *V=ds/pow(size*2+1,3)-dm*dm;
  }
  return 1;
}


int niik_image_nlseg_patch_sad_float_images(nifti_image *img,nifti_image *ref,int *ijk,int nsearch,double *sad,int *num) {
  float *fimg,*fref;
  int ii[3],j=0,dim2,nn;
  double d=0;
  fimg=(float *)img->data;
  fref=(float *)ref->data;
  dim2=img->nx*img->ny;
  for(ii[0]=ijk[0]-nsearch; ii[0]<=ijk[0]+nsearch; ii[0]++) {
    if(ii[0]<0) continue;
    else if(ii[0]>=img->nx) break;
    for(ii[1]=ijk[1]-nsearch; ii[1]<=ijk[1]+nsearch; ii[1]++) {
      if(ii[1]<0) continue;
      else if(ii[1]>=img->ny) break;
      for(ii[2]=ijk[2]-nsearch; ii[2]<=ijk[2]+nsearch; ii[2]++) {
        if(ii[2]<0) continue;
        else if(ii[2]>=img->nz) break;
        nn=ii[0]+ii[1]*img->nx+ii[2]*dim2;
        j++;
        d+=fabs(fimg[nn]-fref[nn]);
      }
    }
  }
  *sad=d;
  *num=j;
  return 1;
}


#endif /* _FALCON_NLSEG_C_ */
