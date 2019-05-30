/* Filename:     nifti1_kunio_bounds.c
 * Description:  functions for creating boundary image
 * Author:       Kunio Nakamura
 * Date:         March 4, 2012
 */

#ifndef _FALCON_BOUNDS_C_
#define _FALCON_BOUNDS_C_

#include "falcon.h"

int niik_image_boundary(nifti_image *img,nifti_image *maskimg,double *dval, int flag_dim_check, int flag_color, int flag_inclusive)
/* -puts boundary color or intensity
 * -img is the overlay image (3D)
 * -maskimg is the mask image
 * -dval is the intensity or color
 * -flag_dim_check is dimension check flag
 * --if slice thickness is large, then it becomes 2D
 * -flag_color is the flag for color
 * -flag_inclusive puts the boundary inside the mask
 */
{
  int
  xdim,area,size,
       n,i,j,k;
  double
  *dimg;
  unsigned char *bimg;
  char fcname[32]="niik_image_boundary";
  int verbose=0;

  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((maskimg==NULL),fcname,"maskimg is null");
  if(img->nx != maskimg->nx) {
    fprintf(stderr,"ERROR: xdim is different %i %i\n",img->nx,maskimg->nx);
    return 0;
  }
  if(img->ny != maskimg->ny) {
    fprintf(stderr,"ERROR: ydim is different %i %i\n",img->ny,maskimg->ny);
    return 0;
  }
  if(img->nz != maskimg->nz) {
    fprintf(stderr,"ERROR: zdim is different %i %i\n",img->nz,maskimg->nz);
    return 0;
  }

  dimg = niik_image_get_voxels_as_double_vector(img);
  bimg = niik_image_get_voxels_as_uint8_vector(maskimg);

  xdim = img->nx;
  area = xdim * img->ny;
  size = area * img->nz;

  if(flag_dim_check) {
    if(img->dz < img->dx*2 && img->dz < img->dy*2) {
      flag_dim_check=0;
    }
  }

  if(!flag_color) { /* intensity values */

    if(niik_image_cmp_dim(img,maskimg)) {
      fprintf(stderr,"ERROR: img and maskimg have different dim\n");
      fprintf(stderr,"       img  %i  %3i %3i %3i %3i %3i  %s\n",img->ndim,img->nx,img->ny,img->nz,img->nt,img->nu,img->fname);
      fprintf(stderr,"       mask %i  %3i %3i %3i %3i %3i  %s\n",maskimg->ndim,maskimg->nx,maskimg->ny,maskimg->nz,maskimg->nt,maskimg->nu,maskimg->fname);
      return 0;
    }

    if(img->nu>1) {
      fprintf(stderr,"ERROR: image can't have multiple bands (u %i)\n",img->nu);
      return 0;
    }
    if(maskimg->nu>1) {
      fprintf(stderr,"ERROR: mask image can't have multiple bands (u %i)\n",maskimg->nu);
      return 0;
    }
    if(img->nt>1) {
      fprintf(stderr,"ERROR: image can't have time component (t %i)\n",img->nt);
      return 0;
    }
    if(maskimg->nt>1) {
      fprintf(stderr,"ERROR: mask image can't have time components (t %i)\n",maskimg->nt);
      return 0;
    }

    if(flag_inclusive) {
      for(k=n=0; k<img->dim[3]; k++) {
        for(j=0; j<img->dim[2]; j++) {
          for(i=0; i<img->dim[1]; n++,i++) {
            if(bimg[n]==0) continue;
            if(i>0) {
              if(!bimg[n-   1]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(j>0) {
              if(!bimg[n-xdim]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k>0) {
                if(!bimg[n-area]) {
                  dimg[n] = dval[0];
                  continue;
                }
              }
            }
            if(i<img->dim[1]-1) {
              if(!bimg[n+   1]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(j<img->dim[2]-1) {
              if(!bimg[n+xdim]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k<img->dim[3]-1) {
                if(!bimg[n+area]) {
                  dimg[n] = dval[0];
                  continue;
                }
              }
            }
          }
        }
      }
    } /* inclusive = change occurs within the mask */

    else {
      for(k=n=0; k<img->dim[3]; k++) {
        for(j=0; j<img->dim[2]; j++) {
          for(i=0; i<img->dim[1]; n++,i++) {
            if(bimg[n]>0) continue;
            if(i>0) {
              if(bimg[n-   1]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(j>0) {
              if(bimg[n-xdim]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k>0) {
                if(bimg[n-area]) {
                  dimg[n] = dval[0];
                  continue;
                }
              }
            }
            if(i<img->dim[1]-1) {
              if(bimg[n+   1]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(j<img->dim[2]-1) {
              if(bimg[n+xdim]) {
                dimg[n] = dval[0];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k<img->dim[3]-1) {
                if(bimg[n+area]) {
                  dimg[n] = dval[0];
                  continue;
                }
              }
            }
          }
        }
      }
    } /* exclusive = change occurs outside the mask */
  } /* gray scale */


  else { /* with color */

    if(img->nv!=3) {
      free(dimg);
      img->nv=img->dim[5]=3;
      img->dv=img->pixdim[5]=1;
      img->ndim=6;
      img->nvox=img->nx*img->ny*img->nz*img->nu*img->nv;
      dimg = (double *)calloc(img->nx*img->ny*img->nz*img->nu*img->nv,sizeof(double));
      size = img->nvox / 3;
      for(i=0; i<img->nvox; i++) {
        dimg[i] = niik_image_get_voxel(img,i%size);
      }
      free(img->data);
      img->data = (void *)calloc(img->nvox,img->nbyper);
    }

    if(flag_inclusive) {
      fprintf(stderr,"ERROR: not implemented yet\n");
      return 0;
    }  /* color and inclusive */

    else {
      for(k=n=0; k<img->dim[3]; k++) {
        for(j=0; j<img->dim[2]; j++) {
          for(i=0; i<img->dim[1]; n++,i++) {
            if(bimg[n]>0) continue;
            if(i>0) {
              if(bimg[n-   1]) {
                dimg[n] = dval[0];
                dimg[n+size] = dval[1];
                dimg[n+size*2] = dval[2];
                continue;
              }
            }
            if(j>0) {
              if(bimg[n-xdim]) {
                dimg[n] = dval[0];
                dimg[n+size] = dval[1];
                dimg[n+size*2] = dval[2];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k>0) {
                if(bimg[n-area]) {
                  dimg[n] = dval[0];
                  dimg[n+size] = dval[1];
                  dimg[n+size*2] = dval[2];
                  continue;
                }
              }
            }
            if(i<img->dim[1]-1) {
              if(bimg[n+   1]) {
                dimg[n] = dval[0];
                dimg[n+size] = dval[1];
                dimg[n+size*2] = dval[2];
                continue;
              }
            }
            if(j<img->dim[2]-1) {
              if(bimg[n+xdim]) {
                dimg[n] = dval[0];
                dimg[n+size] = dval[1];
                dimg[n+size*2] = dval[2];
                continue;
              }
            }
            if(!flag_dim_check) {
              if(k<img->dim[3]-1) {
                if(bimg[n+area]) {
                  dimg[n] = dval[0];
                  dimg[n+size] = dval[1];
                  dimg[n+size*2] = dval[2];
                  continue;
                }
              }
            }
          }
        }
      }
    } /* color and not inclusive */

  } /* color */

  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }

  free(dimg);
  free(bimg);

  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* int niik_image_boundary(nifti_image *img,nifti_image *maskimg,double *dval, int flag_color,int flag_inclusive); */


nifti_image *niik_image_mask_add_red_color_uint8(nifti_image *img,nifti_image *maskimg) {
  nifti_image *outimg=NULL;
  if((outimg=niik_image_mask_add_red_color(img,maskimg,255))==NULL) {
    fprintf(stderr,"[niik_image_mask_add_red_color_uint8] ERROR: niik_image_mask_add_red_color\n");
    return NULL;
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"[niik_image_mask_add_red_color_uint8] ERROR: niik_image_type_convert\n");
    exit(0);
  }
  return outimg;
}

nifti_image *niik_image_mask_add_red_color(nifti_image *img,nifti_image *maskimg,double scale_max) {
  nifti_image *outimg=NULL;
  int i,size;
  char fcname[64]="niik_mask_add_red_color";
  fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return NULL;
  }
  fprintf(stdout,"[%s]   convert to color image\n",fcname);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
    return NULL;
  }
  size=img->nvox;
  if(!niik_check_double_problem(scale_max)) {
    fprintf(stdout,"[%s]   scaling max %f\n",fcname,scale_max);
    if(!niik_image_iscale(outimg,niik_image_get_percentile(img,NULL,0.1),niik_image_get_percentile(img,maskimg,0.98),0,scale_max)) {
      fprintf(stderr,"[%s] ERROR: niik_image_iscale\n",fcname);
      return NULL;
    }
  }
  if(!niik_image_convert_to_color_image(outimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_convert_to_color_image\n",fcname);
    exit(0);
  }
  /*fprintf(stdout,"[%s]   mask is red\n",fcname);*/
  for(i=0; i<size; i++) {
    if(niik_image_get_voxel(maskimg,i)>0)  {
      niik_image_mul_voxel(outimg,i+size,0.5);
      niik_image_mul_voxel(outimg,i+size*2,0.5);
    }
  }
  return outimg;
} /* nifti_image *niik_image_mask_add_red_color(nifti_image *img,nifti_image *maskimg,double scale_max); */


#endif /* _FALCON_BOUNDS_C_ */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/