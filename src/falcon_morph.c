/* Filename:     niikmath.c
 * Description:  general nifti1 program
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 *
 * Revision:     Kunio Nakamura, May 25, 2012
 *               -removed parallel processing using 'pragma omp for'
 *               -this should have removed race condition
 *
 */


#include "falcon.h"
#include "falcon_morph.h"

int niik_image_morph_3d_radius(nifti_image *img,int morph_type,double radius) {
  if(!niik_image_morph_3d_radius_mask(img,NULL,morph_type,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return 0;
  }
  return 1;
}


int niik_image_morph_3d_radius_mask(nifti_image *img,nifti_image *maskimg,int morph_type,double radius) {
  int
  nx,ny,nz,
  area,
  kx,ky,kz,ksize,
  ii,jj,kk,
  p1,p2,p3,
  m,n,i,j,k;
  unsigned char
  *bmask=NULL,
   *bimg=NULL,
    *btmp=NULL;
  double
  rad2,
  d3,d2,d1,dd,
  dx,dy,dz;
  int verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  /* if radius is smaller than the smallest pixdim, then
   * it's the same image */
  dd=NIIK_DMIN(NIIK_DMIN(img->pixdim[1],img->pixdim[2]),img->pixdim[3]);
  if(radius<dd) return 1;
  if(radius==0) return 1;
  kx=floor(radius/img->pixdim[1]+1.0);
  ky=floor(radius/img->pixdim[2]+1.0);
  kz=floor(radius/img->pixdim[3]+1.0);
  ksize=(2*kx+1)*(2*ky+1)*(2*kz+1);
  rad2 = radius * radius;
  area = img->nx * img->ny;
  if(verbose) {
    fprintf(stdout,"\t  kernel %3i %3i %3i  |  ksize = %i \n",kx,ky,kz,ksize);
  }
  /* copy image for dilation and erosion */
  switch(morph_type) {
  case NIIK_MORPH_DILATE:
  case NIIK_MORPH_ERODE:
    if((bimg = niik_image_get_voxels_as_uint8_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector(img)\n");
      return 0;
    }
    if((btmp = niik_image_get_voxels_as_uint8_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector(img)\n");
      return 0;
    }
    if(maskimg!=NULL) {
      if((bmask = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
        fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
        return 0;
      }
    }
    break;
  default:
    break;
  }
  /* do processing */
  switch(morph_type) {

  case NIIK_MORPH_DILATE:
    if(maskimg==NULL) {
      /* #pragma omp parallel for private(j,i,m,kk,nz,dz,d3,p3,jj,ny,dy,d2,p2,ii,nx,dx,p1)*/
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = k * area + j * img->nx;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              dz=kk*img->pixdim[3];
              d3=dz*dz;
              if(d3>rad2) continue;
              p3 = nz*area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                dy=jj*img->pixdim[2];
                d2=d3+dy*dy;
                if(d2>rad2) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  dx=ii*img->pixdim[1];
                  d1=d2+dx*dx;
                  if(d1>rad2) continue;
                  p1 = nx + p2;
                  bimg[p1]=1;
                }
              }
            }
          }
        }
      }
    }

    else { /* DILATION with mask */
      /* #pragma omp parallel for private(j,i,m,kk,nz,dz,d3,p3,jj,ny,dy,d2,p2,ii,nx,dx,p1)*/
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = k * area + j * img->nx;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            if(bmask[m]==0) continue;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              dz=kk*img->pixdim[3];
              d3=dz*dz;
              if(d3>rad2) continue;
              p3 = nz*area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                dy=jj*img->pixdim[2];
                d2=d3+dy*dy;
                if(d2>rad2) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  dx=ii*img->pixdim[1];
                  d1=d2+dx*dx;
                  if(d1>rad2) continue;
                  p1 = nx + p2;
                  if(bmask[p1]==0) continue;
                  bimg[p1]=1;
                }
              }
            }
          }
        }
      }
      free(bmask);
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break; /* NIIK_MORPH_DILATE */

  case NIIK_MORPH_ERODE:
    if(maskimg==NULL) {
      /* #pragma omp parallel for private(j,i,m,n,kk,nz,dz,d3,p3,jj,ny,dy,d2,p2,ii,nx,dx,p1)*/
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = j * img->nx + k * area;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            n=1;
            for(kk=-kz; kk<=kz && n; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              dz=kk*img->pixdim[3];
              d3=dz*dz;
              if(d3>rad2) continue;
              p3 = nz * area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                dy=jj*img->pixdim[2];
                d2=d3+dy*dy;
                if(d2>rad2) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  dx=ii*img->pixdim[1];
                  d1=d2+dx*dx;
                  if(d1>rad2) continue;
                  p1 = p2 + nx;
                  if(btmp[p1]==0) {
                    n=0;
                  }
                }
              }
            }
            if(!n) {
              bimg[m]=0;
            }
          }
        }
      }
    }

    else { /* EROSION with mask */
      /* #pragma omp parallel for private(j,i,m,n,kk,nz,dz,d3,p3,jj,ny,dy,d2,p2,ii,nx,dx,p1)*/
      for(k=0; k<img->nz; k++) {
        m = k * area;
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; m++,i++) {
            if( btmp[m]==0) continue;
            if(bmask[m]==0) continue;
            n=1;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              dz=kk*img->pixdim[3];
              d3=dz*dz;
              if(d3>rad2) continue;
              p3 = nz * area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                dy=jj*img->pixdim[2];
                d2=d3+dy*dy;
                if(d2>rad2) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  dx=ii*img->pixdim[1];
                  d1=d2+dx*dx;
                  if(d1>rad2) continue;
                  p1 = p2 + nx;
                  if(bmask[p1]==0) continue;
                  if(btmp[p1]==0) {
                    n=0;
                    ii=jj=kk=99999;
                  }
                }
              }
            }
            if(!n) {
              bimg[m]=0;
            }
          }
        }
      }
      free(bmask);
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break; /* NIIK_MORPH_ERODE */

  case NIIK_MORPH_CLOSE:
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_DILATE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_ERODE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    break;
  case NIIK_MORPH_OPEN:
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_ERODE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    if(!niik_image_morph_3d_radius_mask(img,maskimg,NIIK_MORPH_DILATE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknonw morph type %i\n",morph_type);
    return 0;
  }
  return 1;
}




/******************************************************************
 * morphologic dilation/erosion/opening/closing function
 *
 * -Rmap is a map (float type) with voxel-wise radius for operations
 * -returns 1 for success
 * -img is updated
 ******************************************************************/

int niik_image_morph_3d_radius_map(nifti_image *img,nifti_image *Rmap,int morph_type) {
  float
  *fimg;
  double
  d1,d2,d3,
  R,R2;
  unsigned char
  *bimg,*btmp;
  int
  area,
  nx,ny,nz,
  i,j,k,m,n,
  ii,jj,kk,
  kx,ky,kz;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(Rmap==NULL) {
    fprintf(stderr,"ERROR: Rmap is null\n");
    return 0;
  }
  if(Rmap->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"ERROR: Rmap is not float32\n");
    return 0;
  }
  R=niik_image_get_max(Rmap,NULL);
  kx=R/img->dx+1.0;
  ky=R/img->dy+1.0;
  kz=R/img->dz+1.0;
  area = img->nx*img->ny;
  switch(morph_type) {
  case NIIK_MORPH_DILATE:
    bimg=niik_image_get_voxels_as_uint8_vector(img);
    btmp=niik_image_get_voxels_as_uint8_vector(img);
    fimg=(float *)Rmap->data;
    for(k=m=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; m++,i++) {
          if(bimg[m]==0) continue;
          R2=fimg[m]*fimg[m];
          if(R2==0) continue;
          for(kk=-kz; kk<=kz; kk++) {
            nz=kk+k;
            if(nz<0) continue;
            if(nz>=img->nz) continue;
            d3=NIIK_SQ(kk*img->dz);
            if(d3>R2) continue;
            for(jj=-ky; jj<=ky; jj++) {
              ny=jj+j;
              if(ny<0) continue;
              if(ny>=img->ny) continue;
              d2=d3+NIIK_SQ(jj*img->dy);
              if(d2>R2) continue;
              for(ii=-kx; ii<=kx; ii++) {
                nx=ii+i;
                if(nx<0) continue;
                if(nx>=img->nx) continue;
                d1=d2+NIIK_SQ(ii*img->dx);
                if(d1>R2) continue;
                n=nx+ny*img->nx+nz*area;
                btmp[n]=1;
              }
            }
          }
        }
      }
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,btmp)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector(img,btmp)\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break;

  case NIIK_MORPH_ERODE:
    bimg=niik_image_get_voxels_as_uint8_vector(img);
    btmp=niik_image_get_voxels_as_uint8_vector(img);
    fimg=(float *)Rmap->data;
    for(k=m=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; m++,i++) {
          if(bimg[m]==0) continue;
          R2=fimg[m]*fimg[m];
          for(kk=-kz; kk<=kz; kk++) {
            nz=kk+k;
            if(nz<0) continue;
            if(nz>=img->nz) continue;
            d3=NIIK_SQ(kk*img->dz);
            if(d3>R2) continue;
            for(jj=-ky; jj<=ky; jj++) {
              ny=jj+j;
              if(ny<0) continue;
              if(ny>=img->ny) continue;
              d2=d3+NIIK_SQ(jj*img->dy);
              if(d2>R2) continue;
              for(ii=-kx; ii<=kx; ii++) {
                nx=ii+i;
                if(nx<0) continue;
                if(nx>=img->nx) continue;
                d1=d2+NIIK_SQ(ii*img->dx);
                if(d1>R2) continue;
                n=nx+ny*img->nx+nz*area;
                if(!bimg[n]) {
                  btmp[m]=0;
                  ii=jj=kk=9999;
                }
              }
            }
          }
        }
      }
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,btmp)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector(img,btmp)\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break;

  case NIIK_MORPH_CLOSE:
    if(!niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)\n");
      return 0;
    }
    if(!niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_ERODE)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)\n");
      return 0;
    }
    break;

  case NIIK_MORPH_OPEN:
    if(!niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_ERODE)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)\n");
      return 0;
    }
    if(!niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius_map(img,Rmap,NIIK_MORPH_DILATE)\n");
      return 0;
    }
    break;

  default:
    fprintf(stderr,"ERROR: unknown morph_type %i\n",morph_type);
    return 0;
  }
  return 1;
}



/******************************************************************
 * seed fill program
 *
 * -flag_grad = 1 -> output image (initimg) will be a gradient image (UINT64)
 *            = 0 -> output image (initimg) will be a mask image (UINT8)
 ******************************************************************/

int niik_image_seed_fill(nifti_image *img,nifti_image *initimg,int flag_grad)
/****************************************************
 * niik_image_seed_fill
 * -img is replaced
 * -idx has the ijk coordinate
 * -flag_grad is the same as before
 * -returns 0 for failure / 1 for success
 * -initial version was seed_fill_slow
 * -memory intensive
 * -2012-09-11, Kunio
 * --corrected problem with (flag_grad>0) that 'iter' was not updated properly
 *
 *****************************************************/
{
  double
  *dimg;
  unsigned long
  *din1,*din2;
  int
  xdim,area,
       i,j,k,n,m,nn,mm,
       iter=1,num=1;
  unsigned char *bimg;
  int verbose=0;
  if( img     == NULL ) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if( initimg == NULL ) {
    fprintf(stderr,"ERROR: initimg is a null pointer\n");
    return 0;
  }
  if(verbose>0) fprintf(stdout,"niik_image_cmp_dim\n");
  if(niik_image_cmp_dim(img,initimg)!=0) {
    fprintf(stderr,"ERROR: niik_image_cmp_dim\n");
    return 0;
  }
  /* prepare for iteration */
  if((num = niik_image_count_mask(img))<0) {
    fprintf(stderr,"ERROR: niik_image_count_mask\n");
    return 0;
  }
  if(verbose>1) fprintf(stdout,"    mask count = %i \n",num);
  din1 = (unsigned long *)calloc(num,sizeof(unsigned long));
  din2 = (unsigned long *)calloc(num,sizeof(unsigned long));
  dimg = niik_image_get_voxels_as_double_vector(initimg);
  bimg = niik_image_get_voxels_as_uint8_vector(img);
  xdim = img->nx;
  area = img->nx * img->ny;
  if(  din1==NULL ||
       din2==NULL ||
       dimg==NULL ||
       bimg==NULL ) {
    fprintf(stderr,"ERROR: memory allocation\n");
    return 0;
  }
  /* get the voxel positions of the initial mask */
  for(i=k=0; i<initimg->nvox; i++) {
    if(dimg[i]>0 && bimg[i]>0) {
      din2[k++]=i;
    }
    NIIK_RET0((k>=num),__func__,"num is too small");
  }
  for(n=0; n<num; n++) {
    din1[n] = din2[n];
  }
  if(verbose>1) fprintf(stdout,"  num = %i \n",num);
  while(num>0) {
    iter++;
    /*if(verbose>1) fprintf(stdout,"  iter = %i \n",iter);*/
    for(n=m=0; n<num; n++) {
      nn = (int)din1[n];
      i =  nn % xdim;
      j = (nn % area) / xdim;
      k =  nn / area;
      /*fprintf(stdout,"%5i [%3i %3i %3i] %i\n",n,i,j,k,nn); */
      /* on the img
       * but not on the initimg */
      if(i>0) {
        mm = nn-1;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
      if(j>0) {
        mm = nn-xdim;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
      if(k>0) {
        mm = nn-area;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
      if(i<img->nx-1) {
        mm = nn+1;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
      if(j<img->ny-1) {
        mm = nn+xdim;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
      if(k<img->nz-1) {
        mm = nn+area;
        if(bimg[mm]>0) {
          if(dimg[mm]==0) {
            dimg[mm] = iter;
            din2[m++] = mm;
          }
        }
      }
    } /* for each frontier voxels */
    num=m;
    for(n=0; n<num; n++) {
      din1[n] = din2[n];
    }
    if(verbose) fprintf(stdout,"\tseed iter %4i %8i\n",iter-1,num);
  }
  if(verbose) {
    fprintf(stdout,"\tfreeing memory\n");
  }
  free(bimg);
  free(din1);
  free(din2);
  if(verbose) {
    fprintf(stdout,"\tfreed memory\n");
    for(i=k=0; i<img->nvox; i++) {
      if(dimg[i]>0) k++;
    }
    fprintf(stdout,"\t\tnew count = %i\n",k);
  }
  /* threshold or smaller image */
  if(!flag_grad) {
    if(verbose) fprintf(stdout,"\t  threshold \n");
    for(i=0; i<img->nvox; i++) {
      dimg[i] = (dimg[i]>0);
    }
    if(!niik_image_set_voxels_from_double_vector(initimg,dimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      return 0;
    }
  } else {
    if(!niik_image_set_voxels_from_double_vector(initimg,dimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      return 0;
    }
    if(!niik_image_auto_type_convert(initimg)) {
      fprintf(stderr,"ERROR: niik_image_auto_type_convert \n");
      return 0;
    }
  }
  free(dimg);
  if(verbose) {
    fprintf(stdout,"\t\tnew count = %i\n",niik_image_count_mask(initimg));
  }
  return 1;
}




/*
 * niik_image_seed_fill_slow
 *
 * -slower version of seed fill
 * -the original version
 * -less memory intensive
 *
 */

int niik_image_seed_fill_slow(nifti_image *img,nifti_image *initimg,int flag_grad) {
  nifti_image *tmpimg;
  int
  imin[8],imax[8],
       i,j,k,n,nch,dim12,
       iter=1;
  unsigned long *ulptr,*ultmp;
  int verbose=0;
  if( img == NULL ) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if( initimg == NULL ) {
    fprintf(stderr,"ERROR: initimg is a null pointer\n");
    return 0;
  }
  if(verbose>1) fprintf(stdout,"niik_image_cmp_dim(img,initimg)\n");
  if(niik_image_cmp_dim(img,initimg)!=0) {
    fprintf(stderr,"ERROR: niik_image_cmp_dim(%s,%s)\n",img->fname,initimg->fname);
    return 0;
  }
  /* prepare for iteration */
  if(verbose>1) fprintf(stdout,"niik_image_type_convert(initimg)\n");
  if(!niik_image_type_convert(initimg,NIFTI_TYPE_UINT64)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return 0;
  }
  if(verbose>1) fprintf(stdout,"niik_image_copy(img)\n");
  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return 0;
  }
  if(verbose>1) fprintf(stdout,"niik_image_type_convert\n");
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT64)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return 0;
  }
  ulptr = (unsigned long *)initimg->data;
  ultmp = (unsigned long *)tmpimg ->data;
  dim12 = initimg->nx*initimg->ny;
  for(i=0; i<initimg->nvox; i++) {
    if(ulptr[i]>0) ulptr[i]=(unsigned long)1;
    else           ulptr[i]=(unsigned long)0;
    if(ultmp[i]>0) ultmp[i]=(unsigned long)1;
    else           ultmp[i]=(unsigned long)0;
  }
  if(verbose>1) fprintf(stdout,"  sum = %i \n",niik_image_count_mask(initimg));
  /* determine the ROI */
  imin[1]=imin[2]=imin[3]=NIIK_VAL_IMAX;
  imax[1]=imax[2]=imax[3]=-NIIK_VAL_IMAX;
  for(k=nch=n=0; k<initimg->nz; k++) {
    for(j=0; j<initimg->ny; j++) {
      for(i=0; i<initimg->nx; n++,i++) {
        if(ulptr[n]==0) continue;
        imin[1]=NIIK_IMIN(i,imin[1]);
        imax[1]=NIIK_IMAX(i,imax[1]);
        imin[2]=NIIK_IMIN(j,imin[2]);
        imax[2]=NIIK_IMAX(j,imax[2]);
        imin[3]=NIIK_IMIN(k,imin[3]);
        imax[3]=NIIK_IMAX(k,imax[3]);
      }
    }
  }
  for(;;) {
    if(verbose>1) fprintf(stdout,"  iter = %i \n",iter);
    for(k=imin[3],nch=n=0; k<=imax[3]; k++) {
      for(j=imin[2]; j<=imax[2]; j++) {
        n = imin[1] + j*initimg->nx + k*dim12;
        for(i=imin[1]; i<=imax[1]; n++,i++) {
          if(ulptr[n]!=(unsigned long)iter) continue;
          if(verbose>3) fprintf(stdout,"[%3i %3i %3i] %lu\n",i,j,k,(unsigned long)iter);
          if(i>0) {
            if(ulptr[n-1]==0) {
              if(ultmp[n-1]>0) {
                nch++;
                imin[1]=NIIK_IMIN(imin[1],i-1);
                ulptr[n-1]=iter+1;
              }
            }
          }
          if(i<initimg->nx-1) {
            if(ulptr[n+1]==0) {
              if(ultmp[n+1]>0) {
                nch++;
                imax[1]=NIIK_IMAX(imax[1],i+1);
                ulptr[n+1]=iter+1;
              }
            }
          }
          if(j>0) {
            if(ulptr[n-(initimg->nx)]==0) {
              if(ultmp[n-(initimg->nx)]>0) {
                nch++;
                imin[2]=NIIK_IMIN(imin[2],j-1);
                ulptr[n-(initimg->nx)]=iter+1;
              }
            }
          }
          if(j<initimg->ny-1) {
            if(ulptr[n+(initimg->nx)]==0) {
              if(ultmp[n+(initimg->nx)]>0) {
                nch++;
                imax[2]=NIIK_IMAX(imax[2],j+1);
                ulptr[n+(initimg->nx)]=iter+1;
              }
            }
          }
          if(k>0) {
            if(ulptr[n-dim12]==0) {
              if(ultmp[n-dim12]>0) {
                nch++;
                imin[3]=NIIK_IMIN(imin[3],k-1);
                ulptr[n-dim12]=iter+1;
              }
            }
          }
          if(k<initimg->nz-1) {
            if(ulptr[n+dim12]==0) {
              if(ultmp[n+dim12]>0) {
                nch++;
                imax[3]=NIIK_IMAX(imax[3],k+1);
                ulptr[n+dim12]=iter+1;
              }
            }
          }
        }
      }
    }
    if(nch==0) {
      if(verbose) fprintf(stdout,"\tseed iter %4i %8i done\n",iter,nch);
      break;
    }
    if(verbose) fprintf(stdout,"\tseed iter %4i %8i\n",iter,nch);
    iter++;
  }
  /* threshold or smaller image */
  if(!flag_grad) {
    if(verbose) fprintf(stdout,"\t  threshold \n");
    if(!niik_image_threshold(initimg,1.0)) {
      fprintf(stderr,"ERROR: niik_image_threshold\n");
      return 0;
    }
  } else if(!niik_image_auto_type_convert(initimg)) {
    fprintf(stderr,"ERROR: niik_image_auto_type_convert \n");
    return 0;
  }
  return 1;
}





/****************************************************
 * niik_image_seed_fill_xyz
 * -img is replaced
 * -idx has the ijk coordinate
 * -flag_grad is the same as before
 * -returns 0 for failure / 1 for success
 *****************************************************/

int niik_image_seed_fill_xyz(nifti_image *img,int *idx,int flag_grad) {
  nifti_image *origimg;
  int verbose=0;
  char fcname[64]="niik_image_seed_fill_xyz";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(idx==NULL) {
    fprintf(stderr,"[%s] ERROR: idx is a null pointer\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  if((origimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] clear image\n",fcname);
  if(!niik_image_clear(img)) {
    fprintf(stderr,"ERROR: niik_image_clear(%s)\n",img->fname);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] set voxel %i %i %i to 1\n",fcname,idx[1],idx[2],idx[3]);
  if(!niik_image_set_voxel(img,idx[1]+idx[2]*img->nx+idx[3]*img->nx*img->ny,1.0)) {
    fprintf(stderr,"ERROR: nift1_k_set_voxel_value(%s, %i, %i, %i)  \n",img->fname,idx[1],idx[2],idx[3]);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] seed fill\n",fcname);
  if(!niik_image_seed_fill(origimg,img,flag_grad)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill(origimg,img,flag_grad) \n");
    return 0;
  }
  nifti_image_free(origimg);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}



/**********************************************************
 * niik_image_seed_fill_xyz
 * -img is replaced
 * -idx has the ijk coordinate
 * -flag_grad is the same as before
 * -returns 0 for failure / 1 for success
 **********************************************************/

int niik_image_seed_fill_edge(nifti_image *img,int flag_grad) {
  nifti_image *origimg;
  int verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if(verbose) fprintf(stdout,"\tniik_image_copy \n");
  if((origimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return 0;
  }
  if(verbose) fprintf(stdout,"\tniik_image_clear \n");
  if(!niik_image_clear(img)) {
    fprintf(stderr,"ERROR: niik_image_clear(%s)\n",img->fname);
    return 0;
  }
  /* 3d image has more edges */
  if(img->nx>2) {
    if(!niik_image_set_voxels_ROI(img,0,0,0,img->nx,img->ny,0,1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,0,0,0,%i,%i,0,1) \n",img->nx,img->ny);
      nifti_image_free(origimg);
      return 0;
    }
    if(!niik_image_set_voxels_ROI(img,0,0,img->nz-1,img->nx,img->ny,img->nz-1,1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,0,0,0,%i,%i,%i,1) \n",img->nx,img->ny,img->nz-1);
      nifti_image_free(origimg);
      return 0;
    }
  }
  if(verbose) fprintf(stdout,"  count %i\n",niik_image_count_mask(img));
  if(!niik_image_set_voxels_ROI(img,0,0,0,0,img->ny,img->nz,1)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,0,0,0,0,%i,%i,1) \n",img->ny,img->nz);
    nifti_image_free(origimg);
    return 0;
  }
  if(!niik_image_set_voxels_ROI(img,img->nx-1,0,0,img->nx-1,img->ny,img->nz,1)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,%i,0,0,%i,%i,%i,1) \n",img->nx-1,img->nx-1,img->ny,img->nz);
    nifti_image_free(origimg);
    return 0;
  }
  if(verbose) fprintf(stdout,"  count %i\n",niik_image_count_mask(img));
  if(!niik_image_set_voxels_ROI(img,0,0,0,img->nx,0,img->nz,1)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,0,0,0,0,%i,0,%i,1) \n",img->nx,img->nz);
    nifti_image_free(origimg);
    return 0;
  }
  if(!niik_image_set_voxels_ROI(img,0,img->ny-1,0,img->nx-1,img->ny-1,img->nz,1)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_ROI (img,0,%i,0,%i,%i,%i,1) \n",img->ny-1,img->nx-1,img->ny-1,img->nz);
    nifti_image_free(origimg);
    return 0;
  }
  if(verbose) fprintf(stdout,"  count %i\n",niik_image_count_mask(img));
  /* make sure that there is no object at the edge */
  if(!niik_image_mask(img,origimg)) {
    fprintf(stderr,"ERROR: niik_image_maskout\n");
    return 0;
  }
  /* run seed_fill function */
  if(verbose) fprintf(stdout,"\tniik_image_seed_fill\n");
  if(!niik_image_seed_fill(origimg,img,flag_grad)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill (origimg,img,flag_grad) \n");
    return 0;
  }
  if(verbose) fprintf(stdout,"\tnifti_image_free\n");
  nifti_image_free(origimg);
  return 1;
} /* niik_image_seed_fill_edge */


int niik_image_seed_fill_from_middle(nifti_image *img,int flag_grad) {
  int
  i,num,n,
  idx[9];
  unsigned char *bimg;
  bimg = niik_image_get_voxels_as_uint8_vector(img);
  for(i=num=0; i<img->nvox; i++) {
    num += (bimg[i]>0);
  }
  /* fprintf(stdout,"  n = %i\n",num);*/
  for(i=n=0; i<img->nvox; i++) {
    if(bimg[i]==0) continue;
    n++;
    if(n > num/2) break;
  }
  /*fprintf(stdout,"  n = %i %i\n",i,bimg[i]); */
  free(bimg);
  idx[1] = i % img->nx;
  idx[2] = floor((i % (img->nx*img->ny)) / img->nx);
  idx[3] = floor(i / (img->nx*img->ny));
  /*  fprintf(stdout,"%3i %3i %3i \n",idx[1],idx[2],idx[3]);*/
  if(!niik_image_seed_fill_xyz (img,idx,flag_grad)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill_xyz\n");
    return 0;
  }
  return 1;
}



int niik_image_close_holes(nifti_image *img) {
  unsigned char *bimg;
  int i,idx[4];
  int verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if(verbose) fprintf(stdout,"    flip obj/background \n");
  bimg = niik_image_get_voxels_as_uint8_vector(img);
  for(i=0; i<img->nvox; i++) {
    bimg[i] = (bimg[i]==0);
  }
  if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
    return 0;
  }
  free(bimg);
  if(0) {
    if(verbose) fprintf(stdout,"    seed fill from 0/0/0 \n");
    idx[1] = idx[2] = idx[3] = 0;
    if(!niik_image_seed_fill_xyz(img,idx,0)) {
      fprintf(stderr,"ERROR: niik_image_seed_fill_xyz\n");
      return 0;
    }
  } else {
    if(verbose) fprintf(stdout,"    seed fill from edge \n");
    if(!niik_image_seed_fill_edge (img,0) ) {
      fprintf(stderr,"ERROR: niik_image_seed_fill_edge\n");
      return 0;
    }
  }
  if(verbose) fprintf(stdout,"    flip back obj/background \n");
  bimg = niik_image_get_voxels_as_uint8_vector(img);
  for(i=0; i<img->nvox; i++) {
    bimg[i] = (bimg[i]==0);
  }
  if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
    return 0;
  }
  free(bimg);
  return 1;
}


int niik_image_morph_close_brain(nifti_image *segimg,double dilate_radius,double erode_radius) {
  if(dilate_radius>0) {
    if(!niik_image_morph_3d_radius(segimg,NIIK_MORPH_DILATE,dilate_radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
      return 0;
    }
  }
  if(!niik_image_close_holes(segimg)) {
    fprintf(stderr,"ERROR: niik_close_holes\n");
    return 0;
  }
  if(erode_radius>0) {
    if(!niik_image_morph_3d_radius(segimg,NIIK_MORPH_ERODE,erode_radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_raduis\n");
      return 0;
    }
  }
  return 1;
}




/****************************************************
 * niik_image_morph_3d_mask
 * -morphologic filter without pixel size consideration
 *****************************************************/


int niik_image_morph_3d_mask(nifti_image *img,nifti_image *maskimg,int morph_type,int shape, int kdim) {
  int
  nx,ny,nz,
  area,
  kx,ky,kz,ksize,
  ii,jj,kk,
  p1,p2,p3,
  m,n,i,j,k;
  unsigned char
  *bmask=NULL,
   *bimg=NULL,
    *btmp=NULL,
     *bker=NULL;
  int verbose=0;

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if((kdim%2)!=1) {
    fprintf(stderr,"ERROR: kernel size should have be an odd number\n");
    return 0;
  }
  kx=ky=kz=(kdim-1)/2;
  switch(shape) {
  case NIIK_MORPH_3D_KERNEL_DIAMOND:
    ksize = kdim*kdim*kdim;
    bker=(unsigned char *)calloc(ksize,sizeof(char));
    m=(kdim-1)/2;
    for(k=n=0; k<kdim; k++) {
      kk=fabs(k-kz);
      for(j=0; j<kdim; j++) {
        jj=fabs(j-ky);
        for(i=0; i<kdim; n++,i++) {
          ii=fabs(i-kx);
          if(ii+jj+kk<=m) bker[n]=1;
        }
      }
    }
    break;
  case NIIK_MORPH_3D_KERNEL_SQUARE:
    ksize = kdim*kdim*kdim;
    bker=(unsigned char *)calloc(ksize,sizeof(char));
    for(i=0; i<ksize; i++) bker[i]=1;
    break;
  case NIIK_MORPH_2D_KERNEL_DIAMOND:
  case NIIK_MORPH_2D_KERNEL_SQUARE:
    fprintf(stderr,"ERROR: 2d kernel\n");
    return 0;
  default:
    fprintf(stderr,"ERROR: kernel shape is unknown %i\n",shape);
    return 0;
  }
  /* if radius is smaller than the smallest pixdim, then
   * it's the same image */
  area = img->nx * img->ny;
  if(verbose) {
    fprintf(stdout,"\t  kernel %3i %3i %3i  |  ksize = %i \n",kx,ky,kz,kdim);
  }
  /* copy image for dilation and erosion */
  switch(morph_type) {
  case NIIK_MORPH_DILATE:
  case NIIK_MORPH_ERODE:
    if((bimg = niik_image_get_voxels_as_uint8_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector(img)\n");
      return 0;
    }
    if((btmp = niik_image_get_voxels_as_uint8_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector(img)\n");
      return 0;
    }
    if(maskimg!=NULL) {
      if((bmask = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
        fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
        return 0;
      }
    }
    break;
  default:
    break;
  }
  /* do processing */
  switch(morph_type) {
  case NIIK_MORPH_DILATE:
    if(maskimg==NULL) {
      /* #pragma omp parallel for private(j,i,m,kk,nz,p3,jj,ny,p2,ii,nx,p1)*/
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = k * area + j * img->nx;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              p3 = nz*area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  p1 = nx + p2;
                  if(bker[(ii+kx)+(jj+ky)*kdim+(kk+kz)*kdim*kdim]==0) continue;
                  bimg[p1]=1;
                }
              }
            }
          }
        }
      }
    }

    else { /* DILATION with mask */
      /* #pragma omp parallel for private(j,i,m,kk,nz,p3,jj,ny,p2,ii,nx,p1) */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = k * area + j * img->nx;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            if(bmask[m]==0) continue;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              p3 = nz*area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  p1 = nx + p2;
                  if(bmask[p1]==0) continue;
                  if(bker[(ii+kx)+(jj+ky)*kdim+(kk+kz)*kdim*kdim]==0) continue;
                  bimg[p1]=1;
                }
              }
            }
          }
        }
      }
      free(bmask);
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break; /* NIIK_MORPH_DILATE */

  case NIIK_MORPH_ERODE:
    if(maskimg==NULL) {
      /* #pragma omp parallel for private(j,i,m,n,kk,nz,p3,jj,ny,p2,ii,nx,p1) */
      for(k=0; k<img->nz; k++) {
        for(j=0; j<img->ny; j++) {
          m = j * img->nx + k * area;
          for(i=0; i<img->nx; m++,i++) {
            if(btmp[m]==0) continue;
            n=1;
            for(kk=-kz; kk<=kz && n; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              p3 = nz * area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  if(bker[(ii+kx)+(jj+ky)*kdim+(kk+kz)*kdim*kdim]==0) continue;
                  p1 = p2 + nx;
                  if(btmp[p1]==0) {
                    n=0;
                  }
                }
              }
            }
            if(!n) {
              bimg[m]=0;
            }
          }
        }
      }
    }

    else { /* EROSION with mask */
      /* #pragma omp parallel for private(j,i,m,n,kk,nz,p3,jj,ny,p2,ii,nx,p1)*/
      for(k=0; k<img->nz; k++) {
        m = k * area;
        for(j=0; j<img->ny; j++) {
          for(i=0; i<img->nx; m++,i++) {
            if( btmp[m]==0) continue;
            if(bmask[m]==0) continue;
            n=1;
            for(kk=-kz; kk<=kz; kk++) {
              nz=kk+k;
              if(nz<0) continue;
              if(nz>=img->nz) continue;
              p3 = nz * area;
              for(jj=-ky; jj<=ky; jj++) {
                ny=jj+j;
                if(ny<0) continue;
                if(ny>=img->ny) continue;
                p2 = p3 + ny * img->nx;
                for(ii=-kx; ii<=kx; ii++) {
                  nx=ii+i;
                  if(nx<0) continue;
                  if(nx>=img->nx) continue;
                  p1 = p2 + nx;
                  if(bmask[p1]==0) continue;
                  if(bker[(ii+kx)+(jj+ky)*kdim+(kk+kz)*kdim*kdim]==0) continue;
                  if(btmp[p1]==0) {
                    n=0;
                    ii=jj=kk=99999;
                  }
                }
              }
            }
            if(!n) {
              bimg[m]=0;
            }
          }
        }
      }
      free(bmask);
    }
    if(!niik_image_set_voxels_from_uint8_vector(img,bimg)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_uint8_vector\n");
      return 0;
    }
    free(bimg);
    free(btmp);
    break; /* NIIK_MORPH_ERODE */

  case NIIK_MORPH_CLOSE:
    if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_DILATE,shape,kdim)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_ERODE,shape,kdim)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    break;
  case NIIK_MORPH_OPEN:
    if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_ERODE,shape,kdim)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    if(!niik_image_morph_3d_mask(img,maskimg,NIIK_MORPH_DILATE,shape,kdim)) {
      fprintf(stderr,"ERROR: niik_image_morph dilation \n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown morph type %i\n",morph_type);
    return 0;
  }
  return 1;
}



/*
 kate: space-indent on; indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
 */