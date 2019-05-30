/* Filename:     nifti1_kunio_median.c
 * Description:  median functions by Kunio
 * Author:       Kunio Nakamura
 * Date:         February 23, 2012
 *
 *
 *
 * -to calculate the median of an image within an optional mask,
 *  use this function, which is written in nifti1_kunio.c
 *
 *      double niik_image_get_median(nifti_image *img,nifti_image *mask);
 *
 * -a similar function is here,
 *
 *      int niik_image_get_median2(nifti_image *img,nifti_image *mask,double *out);
 *
 *
 */

#ifndef _FALCON_MEDIAN_C_
#define _FALCON_MEDIAN_C_

#include "falcon.h"
#include "falcon_morph.h"


int niik_image_filter_median_radius(nifti_image *img,nifti_image *maskimg,double radius) {
  nifti_image *tmpimg;
  double
  *dimg,
  dd;
  int
  kx,ky,kz,ksize,
  i;
  unsigned char *bimg=NULL;
  int verbose=0;
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) start analysis \n");
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) copying the image \n");
  }
  /* get the image intensities (x2) */
  if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) convert to double \n");
  }
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT64)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return 0;
  }
  /* mask image data pointer */
  if(maskimg!=NULL) {
    if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"ERROR: maskimg->datatype is not NIFTI_TYPE_UINT8\n");
      return 0;
    }
    bimg=maskimg->data;
  } else bimg=NULL;
  /* if radius is smaller than the smallest pixdim, then
   * it's the same image */
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) radius check %f \n",radius);
  }
  dd=NIIK_DMIN(NIIK_DMIN(img->pixdim[1],img->pixdim[2]),img->pixdim[3]);
  if(radius<dd) {
    return 1;
  }
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) pixdim %9.5f %9.5f %9.5f | R = %9.5f \n",
            img->pixdim[1],img->pixdim[2],img->pixdim[3],radius);
  }
  /* create a kernel */
  kx=floor(radius/img->pixdim[1]+0.5);
  ky=floor(radius/img->pixdim[2]+0.5);
  kz=floor(radius/img->pixdim[3]+0.5);
  ksize=(2*kx+1)*(2*ky+1)*(2*kz+1);
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) kernel %3i %3i %3i  |  ksize = %i \n",2*kx+1,2*ky+1,2*kz+1,ksize);
  }
  /* main loop */
  if(maskimg==NULL) {
    #pragma omp parallel for
    for(i=0; i<img->nvox; i++) {
      dimg[i] = niik_image_filter_median_radius_ijk(tmpimg,maskimg,i%img->nx,(i/img->nx)%img->ny,i/img->nx/img->ny,kx,ky,kz,radius);
    }
  } else {
    #pragma omp parallel for
    for(i=0; i<img->nvox; i++) {
      if(!bimg[i]) continue; /* don't change */
      dimg[i] = niik_image_filter_median_radius_ijk(tmpimg,maskimg,i%img->nx,(i/img->nx)%img->ny,i/img->nx/img->ny,kx,ky,kz,radius);
    }
  }
  /* finished main loop */
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) checking for errors \n");
  }
  for(i=0; i<img->nvox; i++) {
    if(niik_check_double_problem(dimg[i])) {
      fprintf(stderr,"ERROR: niik_image_filter_median_radius_ijk\n");
      return 0;
    }
  }
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) update the image \n");
  }
  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    free(dimg);
    nifti_image_free(tmpimg);
    return 0;
  }
  free(dimg);
  nifti_image_free(tmpimg);
  if(verbose) {
    fprintf(stdout,"-v1 (niik_image_filter_median_radius) finishing niik_image_filter_median_radius \n");
  }
  return 1;
}


int niik_image_filter_median_radius_dilated_area(nifti_image *maskimg,double radius) {
  nifti_image
  *roiimg=NULL;
  if((roiimg=niik_image_copy(maskimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy(maskimg)\n");
    return 0;
  }
  if(!niik_image_morph_3d_radius(roiimg,NIIK_MORPH_DILATE,radius*1.5)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius(maskimg,NIIK_MORPH_DILATE,radius*1.5)\n");
    return 0;
  }
  if(!niik_image_filter_median_radius(maskimg,roiimg,radius)) {
    fprintf(stderr,"ERROR: niik_image_filter_median_radius(maskimg,roiimg,radius)\n");
    roiimg=niik_image_free(roiimg);
    return 0;
  }
  roiimg=niik_image_free(roiimg);
  return 1;
}


double niik_image_filter_median_radius_ijk(nifti_image *img,nifti_image *maskimg,int i,int j,int k,int kx,int ky,int kz,double radius) {
  int
  ii,jj,kk,
  nx,ny,nz,
  qx,qy,qz,
  ksize,
  n;
  double
  rad2,
  voxval,
  *voxlist,
  dmin,dmax,
  d1,d2,d3,
  dx,dy,dz;
  unsigned char *bimg=NULL;

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NIIKMAX;
  }
  if(maskimg!=NULL) if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
      fprintf(stderr,"ERROR: maskimg is not uint8\n");
      return NIIKMAX;
    }

  n=0; /* list counter */

  if(maskimg!=NULL) bimg = maskimg->data;
  dmin = dmax = niik_image_get_voxel(img,i+j*img->nx+k*img->nx*img->ny);
  rad2 = radius * radius;

  ksize=(2*kx+1)*(2*ky+1)*(2*kz+1);
  voxlist = niik_calloc_double_vector(ksize);

  if(maskimg==NULL) {
    for(kk=-kz; kk<=kz; kk++) {
      nz = kk + k;
      if(nz<0) continue;
      if(nz>=img->nz) continue;
      dz = kk*img->dz;
      d3 = dz*dz;
      if(d3>rad2) continue;
      qz = nz * img->dim[1]*img->dim[2];
      for(jj=-ky; jj<=ky; jj++) {
        ny=jj+j;
        if(ny<0) continue;
        if(ny>=img->dim[2]) continue;
        dy=jj*img->pixdim[2];
        d2=d3+dy*dy;
        if(d2>rad2) continue;
        qy = qz + ny * img->dim[1];
        for(ii=-kx; ii<=kx; ii++) {
          nx=ii+i;
          if(nx<0) continue;
          if(nx>=img->dim[1]) continue;
          dx=ii*img->pixdim[1];
          d1=d2+dx*dx;
          if(d1>rad2) continue;
          qx = qy + nx;
          /* get the voxel value */
          voxval = niik_image_get_voxel(img,qx);
          /* put the voxel value in the list  */
          voxlist[n] = voxval;
          if     (dmin>voxval) dmin=voxval;
          else if(dmax<voxval) dmax=voxval;
          n++;
        }
      }
    } /* local area */
  } /* no mask */

  else {
    for(kk=-kz; kk<=kz; kk++) {
      nz = kk + k;
      if(nz<0) continue;
      if(nz>=img->nz) continue;
      dz = kk*img->dz;
      d3 = dz*dz;
      if(d3>rad2) continue;
      qz = nz * img->dim[1]*img->dim[2];
      for(jj=-ky; jj<=ky; jj++) {
        ny=jj+j;
        if(ny<0) continue;
        if(ny>=img->dim[2]) continue;
        dy=jj*img->pixdim[2];
        d2=d3+dy*dy;
        if(d2>rad2) continue;
        qy = qz + ny * img->dim[1];
        for(ii=-kx; ii<=kx; ii++) {
          nx=ii+i;
          if(nx<0) continue;
          if(nx>=img->dim[1]) continue;
          dx=ii*img->pixdim[1];
          d1=d2+dx*dx;
          if(d1>rad2) continue;
          qx = qy + nx;
          if(!bimg[qx]) continue;
          /* get the voxel value */
          voxval = niik_image_get_voxel(img,qx);
          /* put the voxel value in the list  */
          voxlist[n] = voxval;
          if     (dmin>voxval) dmin=voxval;
          else if(dmax<voxval) dmax=voxval;
          n++;
        }
      }
    } /* local area */
  } /* using mask */

  /* find the median value  */
  if(dmin==dmax)  {
    voxval=dmin;
  }

  else {

    voxval = niik_median_quicksort_double(voxlist,n);
    /*    else if(0){
      gsl_sort(gsl_voxlist->data,gsl_voxlist->stride,(size_t)n);
      if(n%2) dmin = gsl_voxlist->data[n/2];
      else    dmin = (gsl_voxlist->data[n/2]+gsl_voxlist->data[n/2-1])/2.0;
      if(dmin != voxval ){ fprintf(stdout,"%3i %3i %3i %5.1f %5.1f\n",i,j,k,dmin,voxval); } }
    else {
      gsl_sort(gsl_voxlist->data,gsl_voxlist->stride,(size_t)n);
      voxval=gsl_stats_median_from_sorted_data(gsl_voxlist->data,gsl_voxlist->stride,(size_t)n);
      }*/
  }

  free(voxlist);
  return voxval;
}

double niik_median_quicksort_double_untouch(double *v,int num) {
  double *w,d;
  int n;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return NIIKMAX;
  }
  w=(double *)calloc(num,sizeof(double));
  for(n=0; n<num; n++) w[n]=v[n];
  d = niik_median_quicksort_double(w,num);
  free(w);
  return d;
}

double niik_median_quicksort_double(double *v,int num)
/* median value using halfway quicksort
 * v will be re-arranged, at least up to half */
{
  int i;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return NIIKMAX;
  }
  niik_median_quicksort_double_func(v,0,num-1,num);
  if(num%2) return v[num/2];
  i=num/2;
  return (v[i] + v[i-1])/2.0;
}

void niik_median_quicksort_double_func(double *v,int nlo,int nhi,int num)
/* quicksort for median
 * -faster but only sorts up the some value above median */
{
  double x;
  int i,j;
  x = v[nhi];
  i=nlo;
  j=nhi;
  do {
    while(v[i]<x) {
      i++;
    }
    while(v[j]>x) {
      j--;
    }
    if(i<=j) {
      if(v[i]!=v[j]) {
        NIIK_DSWAP(&v[j],&v[i]);
      }
      i++;
      j--;
    }
  } while(i<=j);
  if(nlo<j) if(num/2<=j+1) if(nlo-1<=num/2) niik_median_quicksort_double_func(v,nlo,j,num);
  if(i<nhi) if(i-1<=num/2) if(num/2<=nhi+1) niik_median_quicksort_double_func(v,i,nhi,num);
}

int niik_image_get_median2(nifti_image *img,nifti_image *mask,double *out) {
  niikvec *v;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  /*if(mask==NULL) { fprintf(stderr,"ERROR: mask is null\n"); return 0; }*/
  if(mask==NULL) {
    v=niikvec_init(img->nvox);
    free(v->v);
    if((v->v=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
      return 0;
    }
  } else {
    if((v=niik_image_get_voxels_as_double_vector_within_mask(img,mask))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector_within_mask\n");
      return 0;
    }
  }
  if(!niikvec_sort(v)) {
    fprintf(stderr,"ERROR: niikvec_sort\n");
    return 0;
  }
  *out=niik_get_median_from_sorted_double_vector(v->v,v->num);
  v=niikvec_free(v);
  return 1;
} /* niik_image_get_median */


#endif /* _FALCON_MEDIAN_C_ */
