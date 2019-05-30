/* Filename:     nifti1_kunio_sobel.c
 * Description:  sobel filtering functions
 * Author:       Kunio Nakamura
 * Date:         March 4, 2012
 *
 * MAJOR REVISION: August 26, 2012, Kunio Nakamura
 * -added non-maximum suppression
 *  nifti_image *niik_image_non_maximum_suppression(nifti_image *img,nifti_image *maskimg,double FWHM,double mean, double stdv);
 * -
 *
 */

#ifndef _FALCON_SOBEL_C_
#define _FALCON_SOBEL_C_

#include "falcon.h"


nifti_image *niik_image_gauss_sobel_filter(nifti_image *img,double FWHM) {
  nifti_image *outimg=NULL;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(outimg,(int)(FWHM*2.0),FWHM)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return NULL;
  }
  if(!niik_image_sobel_filter_update(outimg,'m')) {
    fprintf(stderr,"ERROR: niik_image_sobel_filter_update\n");
    return NULL;
  }
  return outimg;
}

int niik_image_sobel_filter_update(nifti_image *img,char dir)
/* sobel filter with dir
 *   x,y,z,m for directions and magnitude
 *
 * 2013-05-26, Kunio Nakamura, knakamura@mrs.mni.mcgill.ca
 * -fixed error for j-ny in a loop
 *
 */
{
  const char *fcname="niik_image_sobel_filter_update";
  nifti_image *tmp[4];
  double
  *data1=NULL,*data2=NULL,*data3=NULL;
  int
  dim12,
  n,i,j,k,
  verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0(((data1=niik_image_get_voxels_as_double_vector(img))==NULL),fcname,"niik_image_get_voxels_as_double_vector");
  NIIK_RET0(((data2=(double *)calloc(img->nvox,sizeof(double)))==NULL), fcname,"calloc for data2");

  dim12=img->nx*img->ny;
  for(i=0; i<img->nvox; i++) data2[i]=data1[i];

  switch(dir) {
  case 'x':
    if(verbose>0) fprintf(stdout,"-v sobel filter xdir\n");
    /* x-dir: [-1 0 1] */
    for(k=n=0; k<img->nz; k++) {
      for(j=0; j<img->ny; j++)  {
        data2[n]=data1[n+1]-data1[n];
        n++;
        for(i=1; i<img->nx-1; n++,i++)  {
          data2[n]=data1[n+1]-data1[n-1];
        }
        data2[n]=data1[n]-data1[n-1];
        n++;
      }
    }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dx;
    /* y-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]+3*data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n-img->nx]+3*data1[n];
          } else {
            data2[n]=data1[n-img->nx]+2*data1[n]+data1[n+img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* z-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]+3*data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n-dim12]+3*data1[n];
          } else {
            data2[n]=data1[n-dim12]+2*data1[n]+data1[n+dim12];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    if(!niik_image_set_voxels_from_double_vector(img,data1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      free(data1);
      free(data2);
      return 0;
    }
    free(data1);
    free(data2);
    break;

  case 'y':
    if(verbose>0) fprintf(stdout,"-v sobel filter ydir\n");
    /* y-dir: [-1 0 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]-data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n]-data1[n-img->nx];
          } else {
            data2[n]=data1[n+img->nx]-data1[n-img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dy;
    /* x-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(i==0) {
            data2[n]=data1[n+1]+3*data1[n];
          } else if(i==img->nx-1) {
            data2[n]=data1[n-1]+3*data1[n];
          } else {
            data2[n]=data1[n-1]+2*data1[n]+data1[n+1];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* z-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]+3*data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n-dim12]+3*data1[n];
          } else {
            data2[n]=data1[n-dim12]+2*data1[n]+data1[n+dim12];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;

    if(!niik_image_set_voxels_from_double_vector(img,data1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      free(data1);
      free(data2);
      return 0;
    }
    free(data1);
    free(data2);
    break;

  case 'z':
    if(verbose>0) fprintf(stdout,"-v sobel filter zdir\n");
    /* x-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(i==0) {
            data2[n]=data1[n+1]+3*data1[n];
          } else if(i==img->nx-1) {
            data2[n]=data1[n-1]+3*data1[n];
          } else {
            data2[n]=data1[n-1]+2*data1[n]+data1[n+1];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* y-dir: [1 2 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]+3*data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n-img->nx]+3*data1[n];
          } else {
            data2[n]=data1[n-img->nx]+2*data1[n]+data1[n+img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* z-dir: [-1 0 1] */
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]-data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n]-data1[n-dim12];
          } else {
            data2[n]=data1[n+dim12]-data1[n-dim12];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dz;
    if(!niik_image_set_voxels_from_double_vector(img,data1)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector \n");
      free(data1);
      free(data2);
      return 0;
    }
    free(data1);
    free(data2);
    break;

  case 'm':
    /* previous code had an error, 2014-07-22 */
    tmp[0]=niik_image_sobel_filter(img,'x');
    for(i=0; i<img->nvox; i++)
      niik_image_set_voxel(tmp[0],i,NIIK_SQ(niik_image_get_voxel(tmp[0],i)));

    tmp[1]=niik_image_sobel_filter(img,'y');
    for(i=0; i<img->nvox; i++)
      niik_image_add_voxel(tmp[0],i,NIIK_SQ(niik_image_get_voxel(tmp[1],i)));
    tmp[1]=niik_image_free(tmp[1]);

    tmp[1]=niik_image_sobel_filter(img,'z');
    for(i=0; i<img->nvox; i++)
      niik_image_add_voxel(tmp[0],i,NIIK_SQ(niik_image_get_voxel(tmp[1],i)));
    tmp[1]=niik_image_free(tmp[1]);

    for(i=0; i<img->nvox; i++)
      niik_image_set_voxel(   img,i,sqrt(niik_image_get_voxel(tmp[0],i)));
    tmp[0]=niik_image_free(tmp[0]);

    free(data1);
    free(data2);

    return 1;

    if((data3=(double *)calloc(img->nvox,sizeof(double)))==NULL) {
      fprintf(stderr,"ERROR: calloc for data3\n");
      return 0;
    }
    if(verbose>0) fprintf(stdout,"-v sobel filter magnitude\n");
    /* X_DIR */
    /* x-dir: [-1 0 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter x-dir (x)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(i==0) {
            data2[n]=data1[n+1]-data1[n];
          } else if(i==img->nx-1) {
            data2[n]=data1[n]-data1[n-1];
          } else {
            data2[n]=data1[n+1]-data1[n-1];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dx;
    /* y-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter x-dir (y)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]+3*data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n-img->nx]+3*data1[n];
          } else {
            data2[n]=data1[n-img->nx]+2*data1[n]+data1[n+img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* z-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter x-dir (z)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]+3*data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n-dim12]+3*data1[n];
          } else {
            data2[n]=data1[n-dim12]+2*data1[n]+data1[n+dim12];
          }
    /* add to the image */
    for(i=0; i<img->nvox; i++) {
      data3[i] = NIIK_SQ(data2[i]/4.0);
    }
    /* Y_DIR */
    /* y-dir: [-1 0 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter y-dir (y)\n");
    if(!niik_image_get_voxels_as_double_vector_update(img,data1)) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector_update\n");
      free(data1);
      free(data2);
      free(data3);
      return 0;
    }
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]-data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n]-data1[n-img->nx];
          } else {
            data2[n]=data1[n+img->nx]-data1[n-img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dy;
    /* x-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter y-dir (x)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(i==0) {
            data2[n]=data1[n+1]+3*data1[n];
          } else if(i==img->nx-1) {
            data2[n]=data1[n-1]+3*data1[n];
          } else {
            data2[n]=data1[n-1]+2*data1[n]+data1[n+1];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* z-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter y-dir (z)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]+3*data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n-dim12]+3*data1[n];
          } else {
            data2[n]=data1[n-dim12]+2*data1[n]+data1[n+dim12];
          }
    /* add to the image */
    for(i=0; i<img->nvox; i++) {
      data3[i] += NIIK_SQ(data2[i]/4.0);
    }
    /* Z_DIR */
    if(!niik_image_get_voxels_as_double_vector_update(img,data1)) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector_update\n");
      return 0;
    }
    /* z-dir: [-1 0 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter z-dir (z)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(k==0) {
            data2[n]=data1[n+dim12]-data1[n];
          } else if(k==img->nz-1) {
            data2[n]=data1[n]-data1[n-dim12];
          } else {
            data2[n]=data1[n+dim12]-data1[n-dim12];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/img->dz;
    /* y-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter z-dir (y)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(j==0) {
            data2[n]=data1[n+img->nx]+3*data1[n];
          } else if(j==img->ny-1) {
            data2[n]=data1[n-img->nx]+3*data1[n];
          } else {
            data2[n]=data1[n-img->nx]+2*data1[n]+data1[n+img->nx];
          }
    for(i=0; i<img->nvox; i++) data1[i]=data2[i]/4.0;
    /* x-dir: [1 2 1] */
    if(verbose>0) fprintf(stdout,"-v sobel filter z-dir (x)\n");
    for(k=n=0; k<img->nz; k++)
      for(j=0; j<img->ny; j++)
        for(i=0; i<img->nx; n++,i++)
          if(i==0) {
            data2[n]=data1[n+1]+3*data1[n];
          } else if(i==img->nx-1) {
            data2[n]=data1[n-1]+3*data1[n];
          } else {
            data2[n]=data1[n-1]+2*data1[n]+data1[n+1];
          }
    free(data1);
    for(i=0; i<img->nvox; i++) {
      data3[i] = sqrt(data3[i] + NIIK_SQ(data2[i]/4.0));
    }
    free(data2);
    if(verbose>0) fprintf(stdout,"-v sobel filter update\n");
    if(!niik_image_set_voxels_from_double_vector(img,data3)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      free(data3);
      return 0;
    }
    free(data3);
    break;

  default:
    fprintf(stderr,"ERORR: unknown dir %c\n",dir);
    free(data1);
    free(data2);
    return 0;
  }
  if(dir!='m') {
    img->cal_min = niik_image_get_min(img,NULL);
    img->cal_max = niik_image_get_max(img,NULL);
  } else {
    img->cal_min = niik_image_get_min(img,NULL);
    img->cal_max = niik_image_get_percentile(img,NULL,0.90);
  }
  if(verbose>0) fprintf(stdout,"-v sobel finish\n");

  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_sobel_filter */


nifti_image *niik_image_sobel_filter(nifti_image *img,char dir) {
  nifti_image *gimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return NULL;
  }
  if((gimg = niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_sobel_filter_update(gimg,dir)) {
    fprintf(stderr,"ERROR: niik_image_sobel_filter_update\n");
    return NULL;
  }
  return gimg;
} /* niik_image_sobel_filter */


nifti_image *niik_image_gauss_sobel_filters_with_mag_single_output(nifti_image *img,double FWHM) {
  nifti_image *outimg=NULL,*tmpimg=NULL;
  char fcname[64]="niik_image_gauss_sobel_filters_with_mag_single_output";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  if(fabs(FWHM)>1e-5) {
    if((tmpimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
      return NULL;
    }
    if(!niik_image_filter_gaussian_update(tmpimg,(int)(FWHM*2.0),FWHM)) {
      fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian_update\n",fcname);
      return NULL;
    }
    if((outimg=niik_image_sobel_filters_with_mag_single_output(tmpimg))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag_single_output\n",fcname);
      return NULL;
    }
    tmpimg=niik_image_free(tmpimg);
  } else {
    if((outimg=niik_image_sobel_filters_with_mag_single_output(img))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag_single_output\n",fcname);
      return NULL;
    }
  }
  return outimg;
} /* niik_image_gauss_sobel_filters_with_mag_single_output */


nifti_image *niik_image_sobel_filters_with_mag_single_output(nifti_image *img) {
  nifti_image **g=NULL,*out=NULL;
  int n;
  if(img==NULL) {
    fprintf(stderr,"[niik_image_sobel_filters_with_mag_single_output] ERROR: img is null\n");
    return NULL;
  }
  if((g=niik_image_sobel_filters_with_mag(img))==NULL) {
    fprintf(stderr,"[niik_image_sobel_filters_with_mag_single_output] ERROR: niik_image_sobel_filters_with_mag\n");
    return NULL;
  }
  if((out=niik_image_combine(g,4,5,0))==NULL) {
    fprintf(stderr,"[niik_image_sobel_filters_with_mag_single_output] ERROR: niik_image_combine\n");
    for(n=0; n<4; n++) {
      g[n]=niik_image_free(g[n]);
    }
    free(g);
    return NULL;
  }
  for(n=0; n<4; n++) {
    g[n]=niik_image_free(g[n]);
  }
  free(g);
  return out;
} /* niik_image_sobel_filters_with_mag_single_output */


nifti_image **niik_image_sobel_filters_with_mag(nifti_image *img)
/* 2012-04-24 Kunio
 * -returns arrays of images for mag,x,y,z-sobel filters
 * 2012-08-19 Kunio
 * -parallel processing
 * 2013-05-26 Kunio
 * -parallel processing is back...
 */
{
  nifti_image **outimg=NULL;
  int n,nf=0,verbose=0;
  const char *fcname="niik_image_sobel_filters_with_mag";
  if(verbose>=1) niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NULL;
  }
  outimg=(nifti_image **)calloc(4,sizeof(nifti_image *));
  for(n=0; n<=3; n++) {
    NIIK_RET0(((outimg[n] = niik_image_copy(img))==NULL),fcname,"niik_image_copy");
  }

  /*#pragma omp parallel for reduction(+:nf)*/
  for(n=1; n<4; n++) {
    if(verbose>=1) fprintf(stdout,"[%s]   sobel dir %c\n",fcname,n+'x'-1);
    if(!niik_image_sobel_filter_update(outimg[n],'x'+n-1)) {
      fprintf(stderr,"[%s] ERROR: niik_image_sobel_filter_update, %i\n",fcname,n);
      nf++;
    }
    if(verbose>=1) fprintf(stdout,"[%s]   sobel dir %c\n",fcname,n+'x'-1);
  }

  if(verbose>=1) fprintf(stdout,"[%s]   sobel dir m\n",fcname);
  if(nf) {
    fprintf(stderr,"[%s] ERORR: sobel filter\n",fcname);
    outimg[0]=niik_image_free(outimg[0]);
    outimg[1]=niik_image_free(outimg[1]);
    outimg[2]=niik_image_free(outimg[2]);
    outimg[3]=niik_image_free(outimg[3]);
    free(outimg);
    return NULL;
  }

  if(verbose>=1) fprintf(stdout,"[%s]   sobel dir m\n",fcname);

  /*parallelize ?*/
  for(n=0; n<outimg[0]->nvox; n++) {
    niik_image_set_voxel(outimg[0],n,
                         sqrt(NIIK_SQ(niik_image_get_voxel(outimg[1],n)) +
                              NIIK_SQ(niik_image_get_voxel(outimg[2],n)) +
                              NIIK_SQ(niik_image_get_voxel(outimg[3],n))));
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_sobel_filters_with_mag */



nifti_image *niik_image_divergence(nifti_image **grad_with_mag,int norm)
/* 2018-09-27 VF
 * -returns divergence image based on grad array
 * norm - if True, normalize gradients by their norm
 */
{
  nifti_image *outimg=NULL;

  int n,nf=0,verbose=0;
  const char *fcname=__func__;
  if(verbose>=1) niik_fc_display(fcname,1);
  if(grad_with_mag==NULL) {
    fprintf(stderr,"[%s] ERROR: grad_with_mag is a null pointer\n",fcname);
    return NULL;
  }
  if(*grad_with_mag==NULL) {
    fprintf(stderr,"[%s] ERROR: *grad_with_mag is a null pointer\n",fcname);
    return NULL;
  }


  for(n=0; n<3; n++) {
    nifti_image *tmpimg = NULL;
    int k;

    NIIK_RET0(((tmpimg = niik_image_copy(grad_with_mag[n+1]))==NULL), fcname,"niik_image_copy");

    if(norm) {
      for(k=0; k<tmpimg->nvox; k++) {
        double grad_norm = niik_image_get_voxel(grad_with_mag[0], k)+1e-6;
        niik_image_set_voxel(tmpimg,k,
                            niik_image_get_voxel(tmpimg,k)/grad_norm
                            );
      }
    }

    NIIK_RET0(!niik_image_sobel_filter_update(tmpimg, 'x'+n), fcname,"niik_image_sobel_filter_update");

    if(n==0) {
        NIIK_RET0(((outimg = niik_image_copy(tmpimg))==NULL), fcname, "niik_image_copy");
    } else {
      for(k=0; k<outimg->nvox; k++)
        niik_image_set_voxel(outimg, k, 
          niik_image_get_voxel(outimg, k) + 
          niik_image_get_voxel(tmpimg, k) );
    }
    niik_image_free(tmpimg);
  }

  if(verbose>=1) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_sobel_filters_with_mag */



nifti_image **niik_image_sobel_filters(nifti_image *img)
/* 2012-04-09 Kunio
 * -returns arrays of images for x,y,z-sobel filters
 */
{
  nifti_image **outimg;
  int n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return NULL;
  }
  outimg=(nifti_image **)calloc(3,sizeof(nifti_image *));
  for(n=0; n<3; n++) {
    if((outimg[n] = niik_image_copy(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy, %i\n",n);
      return NULL;
    }
    if(!niik_image_sobel_filter_update(outimg[n],'x'+n)) {
      fprintf(stderr,"ERROR: niik_image_sobel_filter_update, %i\n",n);
      return NULL;
    }
  }
  return outimg;
} /* niik_image_sobel_filters */


niikpt niik_image_sobel_filter_niikpt(nifti_image *img,niikpt p) {
  int n;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return niikpt_problem();
  }
  if(p.x<0) return niikpt_problem();
  if(p.y<0) return niikpt_problem();
  if(p.z<0) return niikpt_problem();
  p.x/=img->dx;
  p.y/=img->dy;
  p.z/=img->dz;
  if(p.x>=img->nx) return niikpt_problem();
  if(p.y>=img->ny) return niikpt_problem();
  if(p.z>=img->nz) return niikpt_problem();
  n = floor(p.x) + floor(p.y) * img->nx + floor(p.z) * img->nx * img->ny;
  return niik_image_sobel_filter_voxel(img,n);
}


niikpt niik_image_sobel_filter_voxel_par(nifti_image *img,int p,double a1,double a2,double a3)
/* sobel filter with dir
 *   x,y,z,m for directions and magnitude
 * -a1 = farthest, a2 = intermediate, a3 = direct
 * -normally [a1,a2,a3] = [1,2,4]
 */
{
  double
  sum,
  vals[27];
  niikpt pt;
  int
  idx[27],
      m,i,j,k,
      verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return niikpt_problem();
  }
  for(m=0; m<27; m++) {
    idx[m]=p;
    idx[m]+=(m%3)-1;
    idx[m]+=(((m%9)/3)-1)*img->nx;
    idx[m]+=((m/9)-1)*img->nx*img->ny;
  }
  i=p%img->nx;
  j=(p/img->nx)%img->ny;
  k=p/img->nx/img->ny;
  if(verbose) fprintf(stdout,"-v (niik_image_sobel_filter_voxel) ijk = %i %i %i\n",i,j,k);
  if(i==0) {
    for(m=0; m<27; m+=3)
      idx[m]=p;
  } else if(i==img->nx-1) {
    for(m=0; m<27; m++)  {
      if((m%3)==2) {
        idx[m]=p;
      }
    }
  }
  if(j==0) {
    for(m=0; m<27; m++)
      if((m%9)<3)
        idx[m]=p;
  } else if(j==img->ny-1) {
    for(m=0; m<27; m++)
      if((m%9)>5)
        idx[m]=p;
  }
  if(k==0) {
    for(m=0; m<9; m++) {
      idx[m]=p;
    }
  } else if(k==img->nz-1) {
    for(m=18; m<27; m++) {
      idx[m]=p;
    }
  }
  for(m=0; m<27; m++) {
    vals[m]=niik_image_get_voxel(img,idx[m]);
  }
  if(verbose) {
    fprintf(stdout,"-v (niik_image_sobel_filter_voxel) voxel intensity\n");
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[0],vals[1],vals[2],vals[3],vals[4],vals[5],vals[6],vals[7],vals[8]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[ 9],vals[10],vals[11],vals[12],vals[13],vals[14],vals[15],vals[16],vals[17]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[18],vals[19],vals[20],vals[21],vals[22],vals[23],vals[24],vals[25],vals[26]);
    fprintf(stdout,"-v sobel filter \n");
  }
  sum=a1*4+a2*4+a3;
  pt.x =
    vals[ 2]*a1 - vals[ 0]*a1 +
    vals[ 5]*a2 - vals[ 3]*a2 +
    vals[ 8]*a1 - vals[ 6]*a1 +
    vals[20]*a1 - vals[18]*a1 +
    vals[23]*a2 - vals[21]*a2 +
    vals[26]*a1 - vals[24]*a1 +
    vals[11]*a2 - vals[ 9]*a2 +
    vals[14]*a3 - vals[12]*a3 +
    vals[17]*a2 - vals[15]*a2;
  pt.x /= img->dx * sum;
  pt.y =
    vals[ 6]*a1 - vals[ 0]*a1 +
    vals[ 7]*a2 - vals[ 1]*a2 +
    vals[ 8]*a1 - vals[ 2]*a1 +
    vals[24]*a1 - vals[18]*a1 +
    vals[25]*a2 - vals[19]*a2 +
    vals[26]*a1 - vals[20]*a1 +
    vals[15]*a2 - vals[ 9]*a2 +
    vals[16]*a3 - vals[10]*a3 +
    vals[17]*a2 - vals[11]*a2;
  pt.y /= img->dy * sum;
  pt.z =
    vals[18]*a1 - vals[0]*a1 +
    vals[19]*a2 - vals[1]*a2 +
    vals[20]*a1 - vals[2]*a1 +
    vals[21]*a2 - vals[3]*a2 +
    vals[22]*a3 - vals[4]*a3 +
    vals[23]*a2 - vals[5]*a2 +
    vals[24]*a1 - vals[6]*a1 +
    vals[25]*a2 - vals[7]*a2 +
    vals[26]*a1 - vals[8]*a1;
  pt.z /= img->dz * sum;
  pt.w = sqrt(NIIK_SQ(pt.x) + NIIK_SQ(pt.y) + NIIK_SQ(pt.z));
  return pt;
} /* niik_image_sobel_filter_voxel */


niikpt niik_image_sobel_filter_voxel(nifti_image *img,int p)
/* sobel filter with dir
 *   x,y,z,m for directions and magnitude */
{
  double
  vals[27];
  niikpt pt;
  int
  idx[27],
      m,i,j,k,
      verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return niikpt_problem();
  }
  for(m=0; m<27; m++) {
    idx[m]=p;
    idx[m]+=(m%3)-1;
    idx[m]+=(((m%9)/3)-1)*img->nx;
    idx[m]+=((m/9)-1)*img->nx*img->ny;
  }
  i=p%img->nx;
  j=(p/img->nx)%img->ny;
  k=p/img->nx/img->ny;
  if(verbose) fprintf(stdout,"-v (niik_image_sobel_filter_voxel) ijk = %i %i %i\n",i,j,k);
  if(i==0) {
    for(m=0; m<27; m+=3)
      idx[m]=p;
  } else if(i==img->nx-1) {
    for(m=0; m<27; m++)  {
      if((m%3)==2) {
        idx[m]=p;
      }
    }
  }
  if(j==0) {
    for(m=0; m<27; m++)
      if((m%9)<3)
        idx[m]=p;
  } else if(j==img->ny-1) {
    for(m=0; m<27; m++)
      if((m%9)>5)
        idx[m]=p;
  }
  if(k==0) {
    for(m=0; m<9; m++) {
      idx[m]=p;
    }
  } else if(k==img->nz-1) {
    for(m=18; m<27; m++) {
      idx[m]=p;
    }
  }
  for(m=0; m<27; m++) {
    vals[m]=niik_image_get_voxel(img,idx[m]);
  }
  if(verbose) {
    fprintf(stdout,"-v (niik_image_sobel_filter_voxel) voxel intensity\n");
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[0],vals[1],vals[2],vals[3],vals[4],vals[5],vals[6],vals[7],vals[8]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[ 9],vals[10],vals[11],vals[12],vals[13],vals[14],vals[15],vals[16],vals[17]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[18],vals[19],vals[20],vals[21],vals[22],vals[23],vals[24],vals[25],vals[26]);
    fprintf(stdout,"-v sobel filter \n");
  }
  pt.x =
    vals[ 2]   - vals[ 0]   +
    vals[ 5]*2 - vals[ 3]*2 +
    vals[ 8]   - vals[ 6]   +
    vals[20]   - vals[18]   +
    vals[23]*2 - vals[21]*2 +
    vals[26]   - vals[24]   +
    vals[11]*2 - vals[ 9]*2 +
    vals[14]*4 - vals[12]*4 +
    vals[17]*2 - vals[15]*2;
  pt.x /= img->dx * 16;
  pt.y =
    vals[ 6]   - vals[ 0]   +
    vals[ 7]*2 - vals[ 1]*2 +
    vals[ 8]   - vals[ 2]   +
    vals[24]   - vals[18]   +
    vals[25]*2 - vals[19]*2 +
    vals[26]   - vals[20]   +
    vals[15]*2 - vals[ 9]*2 +
    vals[16]*4 - vals[10]*4 +
    vals[17]*2 - vals[11]*2;
  pt.y /= img->dy * 16;
  pt.z =
    vals[18]   - vals[0]   +
    vals[19]*2 - vals[1]*2 +
    vals[20]   - vals[2]   +
    vals[21]*2 - vals[3]*2 +
    vals[22]*4 - vals[4]*4 +
    vals[23]*2 - vals[5]*2 +
    vals[24]   - vals[6]   +
    vals[25]*2 - vals[7]*2 +
    vals[26]   - vals[8];
  pt.z /= img->dz * 16;
  pt.w = sqrt(NIIK_SQ(pt.x) + NIIK_SQ(pt.y) + NIIK_SQ(pt.z));
  return pt;
} /* niik_image_sobel_filter_voxel */


niikpt niik_image_sobel_filter_voxel2(nifti_image *img,int p)
/* sobel filter with dir
 *   x,y,z,m for directions and magnitude
 * -similar to FSL/SIENA -more focal gradient */
{
  double
  vals[27];
  niikpt pt;
  int
  idx[27],
      m,i,j,k,
      verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return niikpt_problem();
  }
  for(m=0; m<27; m++) {
    idx[m]=p;
    idx[m]+=(m%3)-1;
    idx[m]+=(((m%9)/3)-1)*img->nx;
    idx[m]+=((m/9)-1)*img->nx*img->ny;
  }
  i=p%img->nx;
  j=(p/img->nx)%img->ny;
  k=p/img->nx/img->ny;
  if(verbose) fprintf(stdout,"-v (niik_image_sobel_filter_voxel) ijk = %i %i %i\n",i,j,k);
  if(i==0) {
    for(m=0; m<27; m+=3)
      idx[m]=p;
  } else if(i==img->nx-1) {
    for(m=0; m<27; m++)  {
      if((m%3)==2) {
        idx[m]=p;
      }
    }
  }
  if(j==0) {
    for(m=0; m<27; m++)
      if((m%9)<3)
        idx[m]=p;
  } else if(j==img->ny-1) {
    for(m=0; m<27; m++)
      if((m%9)>5)
        idx[m]=p;
  }
  if(k==0) {
    for(m=0; m<9; m++) {
      idx[m]=p;
    }
  } else if(k==img->nz-1) {
    for(m=18; m<27; m++) {
      idx[m]=p;
    }
  }
  for(m=0; m<27; m++) {
    vals[m]=niik_image_get_voxel(img,idx[m]);
  }
  if(verbose) {
    fprintf(stdout,"-v (niik_image_sobel_filter_voxel) voxel intensity\n");
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[0],vals[1],vals[2],vals[3],vals[4],vals[5],vals[6],vals[7],vals[8]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[ 9],vals[10],vals[11],vals[12],vals[13],vals[14],vals[15],vals[16],vals[17]);
    fprintf(stdout,"   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            vals[18],vals[19],vals[20],vals[21],vals[22],vals[23],vals[24],vals[25],vals[26]);
    fprintf(stdout,"-v sobel filter \n");
  }
  pt.x =
    (vals[ 2] - vals[ 0]   +
     vals[ 8] - vals[ 6]   +
     vals[20] - vals[18]   +
     vals[26] - vals[24] ) * 2.0  +
    (vals[23] - vals[21] +
     vals[ 5] - vals[ 3] +
     vals[11] - vals[ 9] +
     vals[17] - vals[15] ) * 5.0 +
    (vals[14] - vals[12] ) *10;
  pt.x /= img->dx * 38;
  pt.y =
    ( vals[ 6] - vals[ 0] +
      vals[ 8] - vals[ 2] +
      vals[24] - vals[18] +
      vals[26] - vals[20] ) * 2.0 +
    ( vals[ 7] - vals[ 1] +
      vals[25] - vals[19] +
      vals[15] - vals[ 9] +
      vals[17] - vals[11] ) * 5.0 +
    ( vals[16] - vals[10] ) * 10.0;
  pt.y /= img->dy * 38;
  pt.z =
    ( vals[22] - vals[4] ) * 10.0 +
    ( vals[19] - vals[1] +
      vals[21] - vals[3] +
      vals[23] - vals[5] +
      vals[25] - vals[7] ) * 5.0 +
    ( vals[18] - vals[0]   +
      vals[20] - vals[2]   +
      vals[24] - vals[6]   +
      vals[26] - vals[8] ) * 2.0;
  pt.z /= img->dz * 38;
  pt.w = sqrt(NIIK_SQ(pt.x) + NIIK_SQ(pt.y) + NIIK_SQ(pt.z));
  return pt;
} /* niik_image_sobel_filter_voxel */


nifti_image *niik_image_non_maximum_suppression(nifti_image *img,nifti_image *maskimg,double FWHM,double mean, double stdv,double thresh,double pthresh)
/* non-maximum suppression
 */
{
  nifti_image
  **gradimg=NULL,
    *outimg=NULL;
  int
  verbose=1,
  nx,xy,
  i,j,k,m,n,nn;
  char
  fname[512],
        fcname[64]="niik_image_non_maximum_suppression";
  niikpt g,nv[27];
  double dv[27];
  unsigned char *bimg=NULL;
  float *fgrad;

  if(verbose>=1) niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return NULL;
  }

  if(verbose>=1) fprintf(stdout,"[%s] gaussian filtering\n",fcname);
  if((outimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  if(!niik_image_filter_gaussian_update(outimg,(int)(FWHM*2.5),FWHM)) {
    fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian_update\n",fcname);
    return NULL;
  }
  if((gradimg=niik_image_sobel_filters_with_mag(outimg))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag\n",fcname);
    return NULL;
  }

  if(verbose>=2) {
    sprintf(fname,"tmp_thinedge_img.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,outimg);
    sprintf(fname,"tmp_thinedge_gradimgs.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_combine_and_write(fname,gradimg,4,5,0);
  }

  if(niik_check_double_problem(pthresh)) {
    pthresh=0.5;
    if(verbose>=1) fprintf(stdout,"[%s] automatic threshold percentile for gradient image %4.2f\n",fcname,pthresh);
  } else {
    if(verbose>=1) fprintf(stdout,"[%s] threshold percentile for gradient image  %4.2f\n",fcname,pthresh);
  }
  if(niik_check_double_problem(thresh)) {
    thresh = niik_image_get_percentile(gradimg[0],maskimg,pthresh);
    if(verbose>=1) fprintf(stdout,"[%s] threshold for gradient image = %7.2f using %2.0f%% percentile\n",fcname,thresh,100.0*pthresh);
  } else  {
    if(verbose>=1) fprintf(stdout,"[%s] threshold for gradient image = %7.2f\n",fcname,thresh);
  }

  /* mean / stdv */
  if(mean < 0 && stdv < 0) {
    if(verbose>=1) fprintf(stdout,"[%s] no gradient modulation using mean/stdv\n",fcname);
  } else {
    if(verbose>=1) fprintf(stdout,"[%s] gradient modulation using mean/stdv = %8.3f %8.3f\n",fcname,mean,stdv);
    for(i=0; i<img->nvox; i++) {
      niik_image_mul_voxel(gradimg[0],i,NIIK_GaussPDF(niik_image_get_voxel(outimg,i)-mean,stdv));
    }
  }

  if(verbose>=2) {
    sprintf(fname,"tmp_thinedge_pgrad.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,gradimg[0]);
  }

  if(verbose>=1) fprintf(stdout,"[%s] threshold for gradient image\n",fcname);
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  if(!niik_image_copy_data(maskimg,outimg)) {
    fprintf(stderr,"ERROR: niik_image_copy_data\n");
    return NULL;
  }
  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(gradimg[0],i)<thresh)
      niik_image_set_voxel(outimg,i,0);
  }


  /* 3d-directions
   *   x y z xy x-y x-z
   */
  if(verbose>=1) fprintf(stdout,"[%s] prepare directional vectors\n",fcname);
  nv[0]=niikpt_val(1,0,0,0);
  nv[1]=niikpt_val(0,1,0,0);
  nv[2]=niikpt_val(0,0,1,0);
  nv[3]=niikpt_val(1,1,0,0);
  nv[4]=niikpt_val(1,0,1,0);
  nv[5]=niikpt_val(0,1,1,0);
  nv[6]=niikpt_val(1,-1,0,0);
  nv[7]=niikpt_val(1,0,-1,0);
  nv[8]=niikpt_val(0,1,-1,0);
  nv[9]=niikpt_val(1,1,1,0);
  nv[10]=niikpt_val(-1,1,1,0);
  nv[11]=niikpt_val(1,-1,1,0);
  nv[12]=niikpt_val(1,1,-1,0);
  if(verbose>=2) {
    for(n=0; n<13; n++) {
      nv[n]=niikpt_unit(nv[n]);
      fprintf(stdout,"vector %2i   %8.3f %8.3f %8.3f\n",n,nv[n].x,nv[n].y,nv[n].z);
    }
  }

  nx = img->nx;
  xy = img->nx*img->ny;
  bimg=outimg->data;
  fgrad=gradimg[0]->data;
  if(verbose>=1) fprintf(stdout,"[%s] main loop\n",fcname);
  for(k=n=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; n++,i++) {
        if(i==71 && j==146 && k==115 && verbose>=2) {
          fprintf(stdout,"[%3i %3i %3i] %i\n",i,j,k,bimg[n]);
        }
        if(bimg[n]==0) {
          continue;
        }
        g.x=niik_image_get_voxel(gradimg[1],n);
        g.y=niik_image_get_voxel(gradimg[2],n);
        g.z=niik_image_get_voxel(gradimg[3],n);
        for(m=nn=0; m<13; m++) {
          dv[m]=fabs(niikpt_dot(nv[m],g));
          if(dv[m]>dv[nn]) nn=m;
        }

        switch(nn) {
        case 0:
          if     (fgrad[n-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 1:
          if     (fgrad[n-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 2:
          if     (fgrad[n-xy]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 3:
          if     (fgrad[n-nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 4:
          if     (fgrad[n-xy-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 5:
          if     (fgrad[n-xy-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 6:
          if     (fgrad[n+nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 7:
          if     (fgrad[n+xy-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 8:
          if     (fgrad[n+xy-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 9:
          if     (fgrad[n+xy+nx+1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy-nx-1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 10:
          if     (fgrad[n+xy+nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy-nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 11:
          if     (fgrad[n+xy-nx+1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx-1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 12:
          if     (fgrad[n+xy-nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        default:
          break;
        }

        dv[nn]=0;
        for(m=nn=0; m<13; m++) {
          if(dv[m]>dv[nn]) nn=m;
        }
        switch(nn) {
        case 0:
          if     (fgrad[n-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 1:
          if     (fgrad[n-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 2:
          if     (fgrad[n-xy]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 3:
          if     (fgrad[n-nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 4:
          if     (fgrad[n-xy-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 5:
          if     (fgrad[n-xy-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n+xy+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 6:
          if     (fgrad[n+nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 7:
          if     (fgrad[n+xy-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 8:
          if     (fgrad[n+xy-nx]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 9:
          if     (fgrad[n+xy+nx+1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy-nx-1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 10:
          if     (fgrad[n+xy+nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy-nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 11:
          if     (fgrad[n+xy-nx+1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx-1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        case 12:
          if     (fgrad[n+xy-nx-1]>fgrad[n]) bimg[n]=0;
          else if(fgrad[n-xy+nx+1]>fgrad[n]) bimg[n]=0;
          else bimg[n]=1;
          break;
        default:
          break;
        }

        if(i==71 && j==146 && k==115 && verbose>=2) {
          niikpt_disp(g);
          for(m=0; m<13; m++) {
            if(m==nn) fprintf(stdout,"%3i %9.3f   %7.4f %7.4f %7.4f *\n",m,dv[m],nv[m].x,nv[m].y,nv[m].z);
            else      fprintf(stdout,"%3i %9.3f   %7.4f %7.4f %7.4f\n",m,dv[m],nv[m].x,nv[m].y,nv[m].z);
          }
          fprintf(stdout,"  new mask %i\n",bimg[n]);
        } /* verbose */

      }
    }
  }  /* x,y,z */

  if(verbose>=1) niik_fc_display(fcname,0);
  return outimg;
} /* niik_image_non_maximum_suppression */


#endif /* _FALCON_SOBEL_FILTER_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/