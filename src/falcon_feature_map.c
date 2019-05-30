/* Filename:     nifti1_kunio_feature_map.c
 * Description:  feature map functions
 * Author:       Kunio Nakamura
 * Date:         June 4, 2012
 *
 *
 * feature includes:
 *   mean, stdv, skew, kurtosis, median, adev (absolute deviation), mode (not included here)
 *
 * type1 = mean,stdv,skew,kurt,median,abs_dev,gradx,grady,gradz,grad_mag
 * type2 = type1 for 3x3x3,5x5x,7x7x7,9x9x9
 * type3 = mean,var
 */

#ifndef _FALCON_FEATURE_MAP_C_
#define _FALCON_FEATURE_MAP_C_

#include "falcon.h"





double niik_image_voxel_gabor_filter(nifti_image *img,int x,int y,int z,double lambdax,double lambday,double lambdaz,double sigmax,double sigmay,double sigmaz,double rx,double ry,double rz)
/* img = image
 * x,y,z = coordinate
 * lamda x,y,z = frequencies
 * sigma x,y,z = blurring sigma
 * rx,ry,rz = rotation in degrees
 */
{
  char fcname[32]="niik_image_voxel_gabor_filter";
  double const pi2_3div2=15.74960994572243;
  double const sigmamax=4.0;
  double
  fc,vv,
  v=0;
  int
  sx,sy,sz,
  cx,cy,cz,
  xsize,ysize,zsize;
  int const verbose=0;
  niikmat *R=NULL;
  niikpt d,p;
  niikvec *vec=NULL;
  int nv=0;
  if(verbose>0) niik_fc_display(fcname,1);
  // NIIK_EXIT((img==NULL),fcname,"img is null in niik_image_voxel_gabor_filter",1);
  R=niikmat_rotate_matrix(rx,ry,rz);
  fc = 1.0 / (pi2_3div2 * sigmax * sigmay * sigmaz);
  /*lambdax+=1e-9;
    lambday+=1e-9;
    lambdaz+=1e-9;
    sigmax+=1e-9;
    sigmay+=1e-9;
    sigmaz+=1e-9;
  */
  xsize = ceil(sigmax / img->dx * sigmamax);
  ysize = ceil(sigmay / img->dy * sigmamax);
  zsize = ceil(sigmaz / img->dz * sigmamax);
  if(verbose>0) {
    fprintf(stdout,"[%s] sigmas %f %f %f\n",__func__,sigmax,sigmay,sigmaz);
    fprintf(stdout,"[%s] lambda %f %f %f\n",__func__,lambdax,lambday,lambdaz);
    fprintf(stdout,"[%s] sizes %i %i %i\n",__func__,xsize,ysize,zsize);
    fprintf(stdout,"[%s] fc %f\n",__func__,fc);
    niikmat_display(R);
    vec=niikvec_init((2*zsize+1)*(2*ysize+1)*(2*xsize+1));
  }
  d.w = 0;
  for(cz=-zsize; cz<=zsize; cz++) {
    sz=NIIK_IMINMAX(cz+z,0,img->nz-1);
    for(cy=-ysize; cy<=ysize; cy++) {
      sy=NIIK_IMINMAX(cy+y,0,img->ny-1);
      for(cx=-xsize; cx<=xsize; cx++) {
        sx=NIIK_IMINMAX(cx+x,0,img->nx-1);
        d.x = cx*img->dx;
        d.y = cy*img->dy;
        d.z = cz*img->dz;
        p = niikpt_affine_transform(R,d);

        vv = niik_image_get_voxel(img,sx+sy*img->nx+sz*img->nx*img->ny);
        if(verbose>1) fprintf(stdout,"%3i %3i %3i | %6.2f %6.2f %6.2f -> %6.2f %6.2f %6.2f | %12.4f %12.6f %12.6f | %.5g\n",
                                x+cx,y+cy,z+cz,
                                d.x,d.y,d.z,p.x,p.y,p.z,
                                vv,
                                exp(-(pow(d.x/sigmax,2) + pow(d.y/sigmay,2) + pow(d.z/sigmaz,2))*0.5),
                                cos(-NIIK_PI2*(p.x/lambdax+p.y/lambday+p.z/lambdaz)),v);
        v +=
          vv *
          exp(-(pow(d.x/sigmax,2) + pow(d.y/sigmay,2) + pow(d.z/sigmaz,2))*0.5) *
          cos(-NIIK_PI2*(p.x/lambdax))*cos(-NIIK_PI2*(p.y/lambday))*cos(-NIIK_PI2*(p.z/lambdaz));
        if(verbose>0) {
          vec->v[nv++]=
            exp(-(pow(d.x/sigmax,2) + pow(d.y/sigmay,2) + pow(d.z/sigmaz,2))*0.5) *
            //cos(-NIIK_PI2*(p.x/lambdax))*cos(-NIIK_PI2*(p.y/lambday))*cos(-NIIK_PI2*(p.z/lambdaz)); }
            cos(-NIIK_PI2*(p.x/lambdax+p.y/lambday+p.z/lambdaz));
        }
      }
    }
  }
  if(verbose>0) {
    niikvec_write("tmp.txt",vec);
    vec=niikvec_free(vec);
    exit(0);
  }
  // v*=fc;
  R=niikmat_free(R);
  if(verbose>0) niik_fc_display(fcname,0);
  return v;
} /* niik_image_voxel_gabor_filter */


int niik_image_feature_voxel_type1(nifti_image *img,int voxel,int kernel,niikvec *fv)
/* -get the features
 * -img is the image of interest (3d)
 * -voxel is the voxel position in absolute index
 *  e.g., [x,y,z] -> x+y*(img->nx)+z*(img->nx)*(img->ny)
 * -kernel is the half kernel size
 *  e.g., kernel=1 is 3x3x3 square kernel
 * -fv is the output vector, containing
 *     0 = mean, 1 = stdv, 2 = skew, 3 = kurtosis, 4 = median, 5 = abs_dev,
 *     6 = gradx, 7 = grady, 8 = gradz, 9 = grad_mag
 */
{
  niikvec *v=NULL;
  int
  xlo,ylo,zlo,
      xhi,yhi,zhi,
      n,
      x,y,z;
  int const verbose=0;
  niikpt pt;
  char fcname[32]="niik_image_feature_voxel_type1";
  if(verbose) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((fv==NULL),fcname,"fv is null");
  NIIK_RET0((fv->num<10),fcname,"fv is too small");
  x = voxel % img->nx;
  y = (voxel / img->nx) % img->ny;
  z = voxel / img->nx / img->ny;
  if(verbose>1) {
    fprintf(stdout,"[%s]   index %3i %3i %3i\n",fcname,x,y,z);
  }
  NIIK_RET0((x<0),fcname,"x is too small");
  NIIK_RET0((y<0),fcname,"y is too small");
  NIIK_RET0((z<0),fcname,"z is too small");
  NIIK_RET0((x>=img->nx),fcname,"x is too large");
  NIIK_RET0((y>=img->ny),fcname,"y is too large");
  NIIK_RET0((z>=img->nz),fcname,"z is too large");
  xlo = NIIK_IMAX(x-kernel,0);
  ylo = NIIK_IMAX(y-kernel,0);
  zlo = NIIK_IMAX(z-kernel,0);
  xhi = NIIK_IMIN(x+kernel+1,img->nx);
  yhi = NIIK_IMIN(y+kernel+1,img->ny);
  zhi = NIIK_IMIN(z+kernel+1,img->nz);
  v=niikvec_init((xhi-xlo)*(yhi-ylo)*(zhi-zlo));
  for(z=zlo,n=0; z<zhi; z++) {
    for(y=ylo; y<yhi; y++) {
      for(x=xlo; x<xhi; n++,x++) {
        v->v[n]=niik_image_get_voxel(img,x+y*img->nx+z*img->nx*img->ny);
      }
    }
  }
  if(verbose>1) {
    fprintf(stdout,"  [%3i,%3i,%3i] %5i ",x,y,z,v->num);
    niikvec_display(v);
  }
  /* mean, stdv, skew, kurtosis, median and abs_dev */
  if(!niik_get_moments_from_double_vector(v->v,v->num,fv->v,fv->v+1,fv->v+2,fv->v+3,fv->v+5)) {
    fprintf(stderr,"[%s] ERROR: niik_get_moments_from_double_vector\n",fcname);
    return 0;
  }
  fv->v[4]=niik_median_quicksort_double(v->v,v->num);
  pt = niik_image_sobel_filter_voxel(img,voxel);
  fv->v[6]=pt.x;
  fv->v[7]=pt.y;
  fv->v[8]=pt.z;
  fv->v[9]=pt.w;
  /*niikvec_display(fv);*/
  v=niikvec_free(v);
  if(verbose) {
    niik_fc_display(fcname,0);
  }
  return 1;
} /* niik_image_feature_voxel_type1 */

nifti_image *niik_image_feature_map_type1(nifti_image *img,int kernel)
/* create a feature map from 3d image
 *     0 = mean, 1 = stdv, 2 = skew, 3 = kurtosis, 4 = median, 5 = abs_dev,
 *     6 = gradx, 7 = grady, 8 = gradz, 9 = grad_mag
 * output has u-dim
 */
{
  nifti_image
  *outimg=NULL;
  int
  m,n,i,j,k,dim3;
  int const verbose=2;
  niikvec
  **fv=NULL;
  char fcname[32]="niik_image_feature_map_type1";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim!=3),fcname,"img is not 3d");
  fv=(niikvec **)calloc(img->nz,sizeof(niikvec *));
  for(k=0; k<img->nz; k++) fv[k]=niikvec_init(10);
  NIIK_RET0(((outimg = niik_image_init(img->nx,img->ny,img->nz,1,10,0,0,
                                       img->dx,img->dy,img->dz,1,1,0,0,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_iamge_init");
  dim3=img->nx*img->ny*img->nz;
  #pragma omp parallel for private(i,j,n,m)
  for(k=0; k<img->nz; k++) {
    m=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        if(!niik_image_feature_voxel_type1(img,m,kernel,fv[k])) {
          fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type1\n",__func__);
          continue;
        }
        for(n=0; n<10; n++) {
          niik_image_set_voxel(outimg,m+dim3*n,fv[k]->v[n]);
        }  /* each freature dim */
      }
    }
    fv[k]=niikvec_free(fv[k]);
  } /* voxel */
  free(fv);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type1 */


int niik_image_feature_voxel_type2(nifti_image *img,int voxel,niikvec *fv)
/* -get the features
 * -img is the image of interest (3d)
 * -voxel is the voxel position in absolute index
 *  e.g., [x,y,z] -> x+y*(img->nx)+z*(img->nx)*(img->ny)
 * -fv is the output vector, containing
 *     0 = mean,  1 = stdv,  2 = skew,  3 = kurtosis,  4 = median,  5 = abs_dev, (3x3x3)
 *     6 = mean,  7 = stdv,  8 = skew,  9 = kurtosis, 10 = median, 11 = abs_dev, (5x5x5)
 *    12 = mean, 13 = stdv, 14 = skew, 15 = kurtosis, 16 = median, 17 = abs_dev, (7x7x7)
 *    18 = mean, 19 = stdv, 20 = skew, 21 = kurtosis, 22 = median, 23 = abs_dev, (9x9x9)
 */
{
  niikvec *v=NULL;
  double *dptr=NULL;
  int
  kernel,
  xlo,ylo,zlo,
  xhi,yhi,zhi,
  n,
  x0,y0,z0,
  x,y,z;
  int const verbose=0;
  char fcname[32]="niik_image_feature_voxel_type2";
  if(verbose) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((fv==NULL),fcname,"fv is null");
  NIIK_RET0((fv->num!=24),fcname,"fv size is wrong");
  x0 = voxel % img->nx;
  y0 = (voxel / img->nx) % img->ny;
  z0 = voxel / img->nx / img->ny;
  if(verbose>1) {
    fprintf(stdout,"[%s]   index %3i %3i %3i\n",fcname,x0,y0,z0);
  }
  NIIK_RET0((x0<0),fcname,"x is too small");
  NIIK_RET0((y0<0),fcname,"y is too small");
  NIIK_RET0((z0<0),fcname,"z is too small");
  NIIK_RET0((x0>=img->nx),fcname,"x is too large");
  NIIK_RET0((y0>=img->ny),fcname,"y is too large");
  NIIK_RET0((z0>=img->nz),fcname,"z is too large");
  dptr=fv->v;
  for(kernel=1; kernel<=4; kernel++,dptr+=6) {
    xlo = NIIK_IMAX(x0-kernel,0);
    ylo = NIIK_IMAX(y0-kernel,0);
    zlo = NIIK_IMAX(z0-kernel,0);
    xhi = NIIK_IMIN(x0+kernel+1,img->nx);
    yhi = NIIK_IMIN(y0+kernel+1,img->ny);
    zhi = NIIK_IMIN(z0+kernel+1,img->nz);
    v=niikvec_init((xhi-xlo)*(yhi-ylo)*(zhi-zlo));
    for(z=zlo,n=0; z<zhi; z++) {
      for(y=ylo; y<yhi; y++) {
        for(x=xlo; x<xhi; n++,x++) {
          v->v[n]=niik_image_get_voxel(img,x+y*img->nx+z*img->nx*img->ny);
        }
      }
    }
    if(verbose>1) {
      fprintf(stdout,"  [%3i,%3i,%3i] %5i ",x,y,z,v->num);
      niikvec_display(v);
    }
    /* mean, stdv, skew, kurtosis, median and abs_dev */
    if(!niik_get_moments_from_double_vector(v->v,v->num,dptr,dptr+1,dptr+2,dptr+3,dptr+5)) {
      fprintf(stderr,"ERROR: niik_get_moments_from_double_vector\n");
      return 0;
    }
    dptr[4]=niik_median_quicksort_double(v->v,v->num);
    v=niikvec_free(v);
  } /* each kernel size */
  if(verbose) {
    niik_fc_display(fcname,0);
  }
  return 1;
} /* niik_image_feature_voxel_type2 */

nifti_image *niik_image_feature_map_type2(nifti_image *img)
/* create a feature map from 3d image
 */
{
  nifti_image
  *outimg=NULL;
  int
  m,n,i,j,k,dim3;
  int const verbose=2;
  niikvec
  **fv=NULL;
  char fcname[32]="niik_image_feature_map_type2";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim!=3),fcname,"img is not 3d");
  fv=(niikvec **)calloc(img->nz,sizeof(niikvec *));
  for(k=0; k<img->nz; k++) fv[k]=niikvec_init(24);
  NIIK_RET0(((outimg = niik_image_init(img->nx,img->ny,img->nz,1,24,0,0,
                                       img->dx,img->dy,img->dz,1,1,0,0,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_image_init");
  dim3=img->nx*img->ny*img->nz;
  #pragma omp parallel for private(i,j,n,m)
  for(k=0; k<img->nz; k++) {
    m=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        if(!niik_image_feature_voxel_type2(img,m,fv[k])) {
          fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type2\n",fcname);
          continue;
        }
        for(n=0; n<24; n++) {
          niik_image_set_voxel(outimg,m+dim3*n,fv[k]->v[n]);
        }  /* each freature dim */
      }
    }
    fv[k]=niikvec_free(fv[k]);
    // fprintf(stdout,"z = %i\n",k);
  } /* voxel */
  free(fv);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type1 */


int niik_image_feature_voxel_type3(nifti_image *img,int voxel,int kernel,double *fv,int method)
/* -get the features
 * -img is the image of interest (3d)
 * -voxel is the voxel position in absolute index
 *  e.g., [x,y,z] -> x+y*(img->nx)+z*(img->nx)*(img->ny)
 * -kernel is the half kernel size
 *  e.g., kernel=1 is 3x3x3 square kernel
 */
{
  niikvec *v=NULL;
  int
  xlo,ylo,zlo,
      xhi,yhi,zhi,
      n,
      x,y,z;
  int const verbose=0;
  char fcname[32]="niik_image_feature_voxel_type3";
  if(verbose) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((fv==NULL),fcname,"fv is null");
  x = voxel % img->nx;
  y = (voxel / img->nx) % img->ny;
  z = voxel / img->nx / img->ny;
  if(verbose>1) {
    fprintf(stdout,"[niik_image_feature_voxel_type3]   index %3i %3i %3i\n",x,y,z);
  }
  NIIK_RET0((x<0),fcname,"x is too small");
  NIIK_RET0((y<0),fcname,"y is too small");
  NIIK_RET0((z<0),fcname,"z is too small");
  NIIK_RET0((x>=img->nx),fcname,"x is too large");
  NIIK_RET0((y>=img->ny),fcname,"y is too large");
  NIIK_RET0((z>=img->nz),fcname,"z is too large");
  xlo = NIIK_IMAX(x-kernel,0);
  ylo = NIIK_IMAX(y-kernel,0);
  zlo = NIIK_IMAX(z-kernel,0);
  xhi = NIIK_IMIN(x+kernel+1,img->nx);
  yhi = NIIK_IMIN(y+kernel+1,img->ny);
  zhi = NIIK_IMIN(z+kernel+1,img->nz);
  v=niikvec_init((xhi-xlo)*(yhi-ylo)*(zhi-zlo));
  for(z=zlo,n=0; z<zhi; z++) {
    for(y=ylo; y<yhi; y++) {
      for(x=xlo; x<xhi; n++,x++) {
        v->v[n]=niik_image_get_voxel(img,x+y*img->nx+z*img->nx*img->ny);
      }
    }
  }
  if(verbose>1) {
    fprintf(stdout,"  [%3i,%3i,%3i] %5i ",x,y,z,v->num);
    niikvec_display(v);
  }
  /* mean, stdv */
  switch(method) {
  case 1:
    *fv=niik_get_mean_from_double_vector(v->v,v->num);
    break;
  case 2:
    *fv=niik_get_var_from_double_vector(v->v,v->num);
    break;
  default:
    return 0;
  }
  v=niikvec_free(v);
  if(verbose) {
    niik_fc_display(fcname,0);
  }
  return 1;
} /* niik_image_feature_voxel_type3 */

nifti_image *niik_image_feature_map_type3(nifti_image *img,int kernel,int method)
/* create a feature map from 3d image */
{
  nifti_image
  *outimg=NULL;
  int
  m,i,j,k,nxy,
  ii,jj,kk,
  ni,nj,nk,nn,ns;
  int const verbose=0;
  double sum,ssq,*fv,*dimg=NULL;
  char fcname[32]="niik_image_feature_map_type3";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim!=3),fcname,"img is not 3d");
  fv=(double *)calloc(img->nz,sizeof(double));
  NIIK_RET0(((outimg = niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  nxy=img->nx*img->ny;
  NIIK_RET0(((dimg=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector_within_mask");
  #pragma omp parallel for private(i,j,m,ii,jj,kk,ni,nj,nk,nn,sum,ssq,ns)
  for(k=0; k<img->nz; k++) {
    m=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        sum=ssq=0;
        ns=0;
        for(kk=-kernel; kk<=kernel; kk++) {
          nk=kk+k;
          if(nk<0) continue;
          else if(nk>=img->nz) continue;
          for(jj=-kernel; jj<=kernel; jj++) {
            nj=jj+j;
            if(nj<0) continue;
            else if(nj>=img->ny) continue;
            for(ii=-kernel; ii<=kernel; ii++) {
              ni=ii+i;
              if(ni<0) continue;
              else if(ni>=img->nx) continue;
              nn=ni+nj*img->nx+nk*nxy;
              sum+=dimg[nn];
              ssq+=dimg[nn]*dimg[nn];
              ns+=1;
            }
          }
        }
        switch(method) {
        case 1:
          fv[k]=sum/ns;
          break;
        case 2:
          fv[k]=(ssq/ns)-NIIK_SQ(sum/ns);
          break;
        }
        niik_image_set_voxel(outimg,m,fv[k]);
      }
    }
  } /* voxel */
  free(fv);
  free(dimg);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type3 */

nifti_image **niik_image_feature_map_mean_var(nifti_image *img,int kernel)
/* create a feature map from 3d image */
{
  nifti_image
  **outimg=NULL;
  int
  m,i,j,k,nxy,
  ii,jj,kk,
  ni,nj,nk,nn,ns;
  int const verbose=0;
  double sum,ssq,*fv,*dimg=NULL;
  char fcname[32]="niik_image_feature_map_mean_var";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim!=3),fcname,"img is not 3d");
  fv=(double *)calloc(img->nz,sizeof(double));
  NIIK_RET0(((outimg = calloc(2,sizeof(nifti_image *)))==NULL),fcname,"calloc for outimg");
  NIIK_RET0(((outimg[0] = niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  NIIK_RET0(((outimg[1] = niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  nxy=img->nx*img->ny;
  NIIK_RET0(((dimg=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector_within_mask");
  #pragma omp parallel for private(i,j,m,ii,jj,kk,ni,nj,nk,nn,sum,ssq,ns)
  for(k=0; k<img->nz; k++) {
    m=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        sum=ssq=0;
        ns=0;
        for(kk=-kernel; kk<=kernel; kk++) {
          nk=kk+k;
          if(nk<0) continue;
          else if(nk>=img->nz) continue;
          for(jj=-kernel; jj<=kernel; jj++) {
            nj=jj+j;
            if(nj<0) continue;
            else if(nj>=img->ny) continue;
            for(ii=-kernel; ii<=kernel; ii++) {
              ni=ii+i;
              if(ni<0) continue;
              else if(ni>=img->nx) continue;
              nn=ni+nj*img->nx+nk*nxy;
              sum+=dimg[nn];
              ssq+=dimg[nn]*dimg[nn];
              ns+=1;
            }
          }
        }
        niik_image_set_voxel(outimg[0],m,sum/ns);
        niik_image_set_voxel(outimg[1],m,(ssq/ns)-NIIK_SQ(sum/ns));
      }
    }
  } /* voxel */
  free(fv);
  free(dimg);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type3 */

int niik_image_feature_voxel_type4(nifti_image *img,int voxel,niikvec *fv)
/* -get the features
 * -img is the image of interest (3d)
 * -voxel is the voxel position in absolute index
 *  e.g., [x,y,z] -> x+y*(img->nx)+z*(img->nx)*(img->ny)
 * -fv is the output vector, containing
 *     0 = mean, 1 = stdv,   from 3x3x3
 *     2 = mean, 3 = stdv,   from 5x5x5
 *     4 = mean, 5 = stdv,   from 7x7x7
 */
{
  niikvec *v=NULL;
  double *dptr=NULL;
  int
  kernel,
  xlo,ylo,zlo,
  xhi,yhi,zhi,
  n,
  x0,y0,z0,
  x,y,z;
  int const verbose=0;
  char fcname[32]="niik_image_feature_voxel_type4";
  if(verbose) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((fv==NULL),fcname,"fv is null");
  NIIK_RET0((fv->num!=6),fcname,"fv size is wrong");
  x0 = voxel % img->nx;
  y0 = (voxel / img->nx) % img->ny;
  z0 = voxel / img->nx / img->ny;
  if(verbose>1) {
    fprintf(stdout,"[niik_image_feature_voxel_type4]   index %3i %3i %3i\n",x0,y0,z0);
  }
  NIIK_RET0((x0<0),fcname,"x is too small");
  NIIK_RET0((y0<0),fcname,"y is too small");
  NIIK_RET0((z0<0),fcname,"z is too small");
  NIIK_RET0((x0>=img->nx),fcname,"x is too large");
  NIIK_RET0((y0>=img->ny),fcname,"y is too large");
  NIIK_RET0((z0>=img->nz),fcname,"z is too large");
  dptr=fv->v;
  for(kernel=1; kernel<=3; kernel++,dptr+=2) {
    xlo = NIIK_IMAX(x0-kernel,0);
    ylo = NIIK_IMAX(y0-kernel,0);
    zlo = NIIK_IMAX(z0-kernel,0);
    xhi = NIIK_IMIN(x0+kernel+1,img->nx);
    yhi = NIIK_IMIN(y0+kernel+1,img->ny);
    zhi = NIIK_IMIN(z0+kernel+1,img->nz);
    v=niikvec_init((xhi-xlo)*(yhi-ylo)*(zhi-zlo));
    for(z=zlo,n=0; z<zhi; z++) {
      for(y=ylo; y<yhi; y++) {
        for(x=xlo; x<xhi; n++,x++) {
          v->v[n]=niik_image_get_voxel(img,x+y*img->nx+z*img->nx*img->ny);
        }
      }
    }
    if(verbose>1) {
      fprintf(stdout,"  [%3i,%3i,%3i] %5i ",x,y,z,v->num);
      //niikvec_display(v);
    }
    /* mean, stdv */
    dptr[0]=niik_get_mean_from_double_vector(v->v,v->num);
    dptr[1]=niik_get_var_from_double_vector(v->v,v->num);
    v=niikvec_free(v);
  } /* each kernel size */
  if(verbose) {
    niik_fc_display(fcname,0);
  }
  return 1;
} /* niik_image_feature_voxel_type4 */

nifti_image *niik_image_feature_map_type4(nifti_image *img)
/* create a feature map from 3d image
 *     0 = image
 *     1 = mean, 2 = stdv,   from 3x3x3
 *     3 = mean, 4 = stdv,   from 5x5x5
 *     5 = mean, 6 = stdv,   from 7x7x7
 * output has u-dim
 */
{
  nifti_image
  *outimg=NULL;
  int
  m,n,i,j,k,dim3;
  int const verbose=2;
  niikvec
  **fv=NULL;
  char fcname[64]="niik_image_feature_map_type4";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((img->ndim!=3),fcname,"img is not 3d");
  fv=(niikvec **)calloc(img->nz,sizeof(niikvec *));
  for(k=0; k<img->nz; k++) fv[k]=niikvec_init(6);
  NIIK_RET0(((outimg = niik_image_init(img->nx,img->ny,img->nz,7,0,0,0,
                                       img->dx,img->dy,img->dz,1,1,0,0,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_iamge_init");
  dim3=img->nx*img->ny*img->nz;
  #pragma omp parallel for private(i,j,n,m)
  for(k=0; k<img->nz; k++) {
    m=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        if(!niik_image_feature_voxel_type4(img,m,fv[k])) {
          fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type4\n",fcname);
          continue;
        }
        niik_image_set_voxel(outimg,m,niik_image_get_voxel(img,m));
        for(n=1; n<outimg->nt; n++) {
          niik_image_set_voxel(outimg,m+dim3*n,fv[k]->v[n-1]);
        } /* each freature dim */
      }
    }
    fv[k]=niikvec_free(fv[k]);
  } /* voxel */
  free(fv);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type4 */

nifti_image *niik_image_feature_map_type4_multi(nifti_image **imglist,int nimg) {
  nifti_image
  *img=NULL,
   *outimg=NULL;
  int
  m,n,i,j,k,b,dim3,dim4;
  int const verbose=2;
  niikvec
  **fv=NULL;
  char fcname[64]="niik_image_feature_map_type4_multi";
  if(verbose>0) {
    niik_fc_display(fcname,1);
  }
  NIIK_RET0((imglist==NULL),fcname,"imglist is null");
  img=imglist[0];
  for(b=1; b<nimg; b++) {
    if(img->nx!=imglist[b]->nx) {
      fprintf(stderr,"[%s] img nx is different, %i %i at %i\n",fcname,img->nx,imglist[b]->nx,b);
      return NULL;
    }
    if(img->ny!=imglist[b]->ny) {
      fprintf(stderr,"[%s] img ny is different, %i %i at %i\n",fcname,img->ny,imglist[b]->ny,b);
      return NULL;
    }
    if(img->nz!=imglist[b]->nz) {
      fprintf(stderr,"[%s] img nz is different, %i %i at %i\n",fcname,img->nz,imglist[b]->nz,b);
      return NULL;
    }
  }
  fv=(niikvec **)calloc(img->nz,sizeof(niikvec *));
  for(k=0; k<img->nz; k++) fv[k]=niikvec_init(6);
  NIIK_RET0(((outimg = niik_image_init(img->nx,img->ny,img->nz,7,nimg,0,0,
                                       img->dx,img->dy,img->dz,1,1,0,0,NIFTI_TYPE_FLOAT32))==NULL),
            fcname,"niik_iamge_init");
  dim3=outimg->nx*outimg->ny*outimg->nz;
  dim4=dim3*outimg->nt;
  for(b=0; b<nimg; b++) {
    if(verbose>0) {
      fprintf(stdout,"[%s] image list [%i], %s\n",fcname,b,imglist[b]->fname);
    }
    #pragma omp parallel for private(i,j,n,m)
    for(k=0; k<img->nz; k++) {
      m=k*img->nx*img->ny;
      for(j=0; j<img->ny; j++) {
        for(i=0; i<img->nx; m++,i++) {
          if(!niik_image_feature_voxel_type4(imglist[b],m,fv[k])) {
            fprintf(stderr,"[%s] ERROR: niik_image_feature_voxel_type4\n",fcname);
            continue;
          }
          niik_image_set_voxel(outimg,m+dim4*b,niik_image_get_voxel(imglist[b],m));
          for(n=1; n<outimg->nt; n++) {
            niik_image_set_voxel(outimg,m+dim3*n+dim4*b,fv[k]->v[n-1]);
          }  /* each freature dim */
        }
      }
    } /* voxel */
  } /* each 3d image */
  for(k=0; k<img->nz; k++) fv[k]=niikvec_free(fv[k]);
  free(fv);
  if(verbose>0) {
    niik_fc_display(fcname,0);
  }
  return outimg;
} /* niik_image_feature_map_type4_multiple */


#endif /* _FALCON_FEATURE_MAP_C_ */
