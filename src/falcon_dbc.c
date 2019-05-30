/* Filename:     nifti1_kunio_dbc.c
 * Description:  differential bias correction
 * Author:       Kunio Nakamura
 * Date:         April 15, 2012
 *
 * Reference:    Lewis and Fox. 2004. Correction of differential intensity inhomogeneity in longitudinal MR images. NeuroImage 23: 75-83.
 */

#ifndef _FALCON_DBC_C_
#define _FALCON_DBC_C_

#include "falcon.h"


nifti_image *niik_image_dbc(nifti_image *refimg, nifti_image *img, nifti_image *maskimg, niikmat *img2ref, double radius, double lim)
/* differential bias correction
 * -img is corrected
 * -returns the bias field image
 * -refimg is the reference image
 * -maskimg is an optional mask
 * -img2ref is the affine matrix trasnforming the img to refimg
 * -radius is the filter size, should be larege like 7.5
 * -lim is the limit values -should be close to 1, like 1.2
 *
 */
{
  nifti_image *bias_img;
  if((bias_img=niik_image_dbc_with_scaling(refimg,img,maskimg,img2ref,radius,lim,0))==NULL) {
    fprintf(stderr,"ERROR: niik_image_dbc_with_scaling\n");
    return NULL;
  }
  return bias_img;
} /* niik_image_dbc */

nifti_image *niik_image_dbc_with_scaling(nifti_image *refimg, nifti_image *img, nifti_image *maskimg,
    niikmat *img2ref, double radius, double lim, int flag_scale) {
  int verbose=0;
  nifti_image *bias_img=NULL;
  char fcname[32]="niik_image_dbc_with_scaling";
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0(((bias_img=niik_image_dbc_with_scaling_resample(refimg,img,maskimg,img2ref,radius,lim,flag_scale,NULL))==NULL),
            fcname,"niik_image_dbc_with_scaling_resample");
  if(verbose>0) niik_fc_display(fcname,0);
  return bias_img;
} /* niik_image_dbc_with_scaling */

nifti_image *niik_image_dbc_with_scaling_resample(nifti_image *refimg, nifti_image *img, nifti_image *maskimg,
    niikmat *img2ref, double radius, double lim, int flag_scale, double *xyz)
/*
 * differential bias correction
 * -in addition to niik_image_dbc, linear scaling is included
 * -also added optional resampling feature
 * -img is corrected
 * -returns the bias field image
 * -refimg is the reference image
 * -maskimg is an optional mask in img space
 * -img2ref is the affine matrix trasnforming the img to refimg
 * -radius is the filter size, should be large like 7.5
 * -lim is the limit values -should be close to 1, like 1.2
 * -flag_scale is helpful when there is a large offset between refimg and img
 */
{
  nifti_image
  *bias_img=NULL;
  unsigned char
  *bimg;
  niikpt
  pt,qt;
  niikmat
  *afmat=NULL;
  double
  sclpar,
  *sclx=NULL,*scly=NULL,
   *sclpp=NULL,
    hilim,lolim,
    d;
  int
  verbose=1,
  numscl,
  m,n,i,j,k;
  char fcname[64]="niik_image_dbc_with_scaling_resample";

  if(verbose) niik_fc_display(fcname,1);
  NIIK_RET0((refimg==NULL),fcname,"refimg is null");
  NIIK_RET0((   img==NULL),fcname,"img is null");

  NIIK_RET0(((bias_img = niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET0((!niik_image_one(bias_img)),fcname,"niik_image_one");

  if(img2ref!=NULL) {
    afmat = niikmat_copy(img2ref);
  } else {
    afmat = niikmat_identity(4,4);
  }
  if(verbose) {
    fprintf(stdout,"[%s]   afmat\n",fcname);
    fprintf(stdout,"   %8.3f %8.3f %8.3f %8.3f\n",afmat->m[0][0],afmat->m[0][1],afmat->m[0][2],afmat->m[0][3]);
    fprintf(stdout,"   %8.3f %8.3f %8.3f %8.3f\n",afmat->m[1][0],afmat->m[1][1],afmat->m[1][2],afmat->m[1][3]);
    fprintf(stdout,"   %8.3f %8.3f %8.3f %8.3f\n",afmat->m[2][0],afmat->m[2][1],afmat->m[2][2],afmat->m[2][3]);
    fprintf(stdout,"   %8.3f %8.3f %8.3f %8.3f\n",afmat->m[3][0],afmat->m[3][1],afmat->m[3][2],afmat->m[3][3]);
  }

  if(xyz!=NULL) {
    fprintf(stderr,"[%s] NOT IMPLEMENTED YET TO DO RESAMPLING\n",fcname);
    return 0;
  }

  /* linear scaling */
  if(!flag_scale) {
    if(verbose) fprintf(stdout,"[%s]   skip linear scaling\n",fcname);
  } else {
    if(verbose) fprintf(stdout,"[%s]   linear scaling\n",fcname);
    if(maskimg!=NULL) {
      numscl = niik_image_count_mask(maskimg);
    } else {
      numscl = img->nvox;
    }
    sclx=(double *)calloc(numscl,sizeof(double));
    scly=(double *)calloc(numscl,sizeof(double));
    if(maskimg!=NULL) {
      bimg = niik_image_get_voxels_as_uint8_vector(maskimg);
    } else {
      bimg = (unsigned char *)calloc(img->nvox,sizeof(char));
      for(i=0; i<img->nvox; i++) {
        bimg[i]=1;
      }
    }
    for(k=n=m=0; k<img->nz; k++) {
      pt.z = k *img->dz;
      for(j=0; j<img->ny; j++) {
        pt.y = j *img->dy;
        for(i=0; i<img->nx; n++,i++) {
          pt.x = i *img->dx;
          if(bimg[n]==0) continue;
          qt = niikpt_affine_transform(afmat,pt);
          sclx[m]=niik_image_get_voxel(img,n);
          scly[m]=niik_image_interpolate_3d_linear(refimg,qt);
          m++;
        }
      }
    }
    fprintf(stdout,"[%s] linear regression (slope)\n",fcname);
    if(!niik_slope_from_double_vector(sclx,scly,numscl,&sclpar)) {
      fprintf(stderr,"ERROR: niik_slope_from_double_vector\n");
      return NULL;
    }
    fprintf(stdout,"[%s]   slope = %.5f\n",fcname,sclpar);
    if(0) { /* What's this? -> [regression using absolute difference] So is this better? ...??? -Kunio 2012-06-29 */
      sclpp=(double *)calloc(5,sizeof(double));
      niik_slope_from_double_vector_least_absolute_difference(sclx,scly,numscl,sclpp,3);
      /*if(!niik_slope_from_double_vector_least_absolute_difference(sclx,scly,numscl,sclpp+1,sclpp+2)){
      fprintf(stderr,"ERROR: \n");
      exit(0); }*/
      for(n=0; n<3; n++) {
        fprintf(stdout,"\tsclpp[%i] %.6f\n",n,sclpp[n]);
      }
      exit(0);
    }
    free(sclx);
    free(scly);
    free(bimg);
    for(i=0; i<img->nvox; i++) {
      niik_image_mul_voxel(img,i,sclpar);
    }
  } /* flag_scale */

  if(lim>1) hilim=lim;
  else hilim=1.0/lim;
  lolim = 1.0/hilim;
  if(verbose) fprintf(stdout,"[%s] lim %8.3f %8.3f\n",fcname,lolim,hilim);

  if(verbose) fprintf(stdout,"[%s] main loop\n",fcname);
  for(k=n=0; k<img->nz; k++) {
    pt.z = k *img->dz;
    for(j=0; j<img->ny; j++) {
      pt.y = j *img->dy;
      for(i=0; i<img->nx; n++,i++) {
        pt.x = i *img->dx;
        qt = niikpt_affine_transform(afmat,pt);
        d = (1e-4 + niik_image_interpolate_3d(refimg,qt,NIIK_INTERP_LINEAR)) /
            (1e-4 + niik_image_get_voxel(img,n));
        if(d>hilim) d=hilim;
        else if(d<lolim) d=lolim;
        niik_image_set_voxel(bias_img,n,d);
      }
    }
  }
  /* niik_image_write("tmp_ratio.nii.gz",bias_img); */

  afmat = niikmat_free(afmat);

  if(verbose) fprintf(stdout,"[%s] median filter %6.2f\n",fcname,radius);
  NIIK_RET0((!niik_image_filter_median_radius(bias_img,maskimg,radius)),fcname,"niik_image_filter_median_radius");

  /* multiply the bias values */
  if(verbose) fprintf(stdout,"[%s] apply correction\n",fcname);
  if(maskimg!=NULL) {
    NIIK_RET0(((bimg = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL),fcname,"niik_image_get_voxels_as_uint8_vector");
    for(n=0; n<img->nvox; n++) {
      if(bimg[n]>0)
        niik_image_mul_voxel(img,n,niik_image_get_voxel(bias_img,n));
    }
    free(bimg);
  } /* using a mask */
  else {
    for(n=0; n<img->nvox; n++) {
      niik_image_mul_voxel(img,n,niik_image_get_voxel(bias_img,n));
    }
  } /* no mask */

  if(verbose) niik_fc_display(fcname,0);
  return bias_img;
} /* niik_image_dbc_with_scaling_resample */


#endif /* _FALCON_DBC_C_ */
