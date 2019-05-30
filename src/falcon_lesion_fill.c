/* Filename:     nifti1_kunio_lesion_fill.c
 * Description:  lesion filling functions
 * Author:       Kunio Nakamura
 * Date:         June 4, 2012
 *
 */

#ifndef _FALCON_LESION_FILL_C_
#define _FALCON_LESION_FILL_C_

#include "falcon.h"
#include "falcon_morph.h"



/***************************************************************
 *
 * lesion filling
 *
 * from Mishkin's lab meeting presentation
 *
 ***************************************************************/

int niik_image_fill_lesion(nifti_image *img, nifti_image *les_mask, nifti_image *wm_mask) {
  if(!niik_image_fill_lesion_with_matrix(img,les_mask,NULL,wm_mask,NULL)) {
    fprintf(stderr,"ERROR: niik_image_fill_lesion_with_matrix\n");
    return 0;
  }
  return 1;
}

int niik_image_fill_lesion_with_matrix(nifti_image *img, nifti_image *les_mask, niikmat *lmat, nifti_image *wm_mask, niikmat *wmat)
/* -fills lesions with NAWM intensities
 * -img is the t1w image
 * -les_mask is the mask of lesions (t2 or t1)
 * -lmat is the matrix transforming the lesion mask to img
 * -if lmat is null then les_mask and img must be in the same space
 * -wm_mask is the NAWM mask
 * -wmat is the matrix that transforms wm_mask to img space
 */
{
  nifti_image
  *xfmimg=NULL,
   *tmpimg=NULL;
  double
  minvox,
  mean,stdv;
  int i;
  unsigned char
  *bimg;
  char fcname[64]="niik_image_fill_lesion_with_matrix";

  fprintf(stdout,"[%s] start\n",fcname);
  fprintf(stdout,"  image            %s\n",img->fname);
  fprintf(stdout,"  lesion           %s\n",les_mask->fname);
  fprintf(stdout,"  white matter     %s\n",wm_mask->fname);

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(wm_mask==NULL) {
    fprintf(stderr,"ERROR: white matter mask is null\n");
    return 0;
  }
  if(les_mask==NULL) {
    fprintf(stderr,"ERROR: lesion mask is null\n");
    return 0;
  }
  if(img->nvox!=wm_mask->nvox) {
    fprintf(stderr,"ERROR: WM #vox %i %i\n",img->nvox,wm_mask->nvox);
    return 0;
  }
  if(img->nvox!=les_mask->nvox) {
    fprintf(stderr,"ERROR: lesion #vox %i %i\n",img->nvox,les_mask->nvox);
    return 0;
  }

  minvox=1.2*NIIK_DMIN(img->dx,NIIK_DMIN(img->dy,img->dz));
  if(wmat!=NULL) {
    NIIK_RET0(((xfmimg = niik_image_affine_transform_3d(wm_mask,img,wmat,NIIK_INTERP_NN))==NULL),
              fcname,"niik_image_affine_transform_3d");
  } else {
    NIIK_RET0(((xfmimg = niik_image_copy(wm_mask))==NULL),
              fcname,"niik_image_copy");
  }
  NIIK_RET0((!niik_image_morph_3d_radius(xfmimg,NIIK_MORPH_ERODE,minvox)),fcname,"niik_image_morph_3d_radius");
  mean = niik_image_get_mean(img,xfmimg);
  stdv = niik_image_get_stdv(img,xfmimg);
  xfmimg = niik_image_free(xfmimg);
  fprintf(stdout,"[niik_image_fill_lesion_with_matrix] WM %8.3f %8.3f\n",mean,stdv);

  /* get lesion mask */
  if(lmat!=NULL) {
    if((xfmimg = niik_image_affine_transform_3d(les_mask,img,lmat,NIIK_INTERP_NN))==NULL) {
      fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
      return 0;
    }
    if((bimg = niik_image_get_voxels_as_uint8_vector(xfmimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
      return 0;
    }
    xfmimg = niik_image_free(xfmimg);
  } else {
    if((bimg = niik_image_get_voxels_as_uint8_vector(les_mask))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
      return 0;
    }
  }

  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }

  for(i=0; i<img->nvox; i++) {
    if(bimg[i]==0) continue;
    niik_image_set_voxel(tmpimg,i,mean + stdv*niik_get_rand_normal());
  }

  /* -Kunio's original implementation in cleveland did not have this
   *  gaussian filter
   * -it change the distribution so it may not be appropriate
   */
  fprintf(stdout,"[niik_image_fill_lesion_with_matrix] gaussian filter\n");
  if(!niik_image_filter_gaussian_update(tmpimg,1,1)) {
    fprintf(stderr,"ERROR: niik_image_filter_gaussian_update\n");
    return 0;
  }

  /* update the image */
  for(i=0; i<img->nvox; i++) {
    if(bimg[i]==0) continue;
    niik_image_set_voxel(img,i,niik_image_get_voxel(tmpimg,i));
  }

  free(bimg);
  tmpimg = niik_image_copy(tmpimg);
  return 1;
} /* niik_image_fill_lesion_with_matrix */


int niik_image_fill_lesion_with_feature(nifti_image *img, nifti_image *les_mask, niikmat *lmat, nifti_image *wm_mask, niikmat *wmat)
/*
 *  *** This function is not complete ***
 *
 * -fills the lesion in the image with similar feature as normal-appearing white matter
 *
 *
 */
{
  nifti_image
  *xfmlm=NULL,
   *xfmwm=NULL;
  niikvec *fv;
  int
  i,
  kernel=2,
  verbose=1;

  if(verbose) fprintf(stdout,"[niik_image_fill_lesion_with_feature] start\n");
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(les_mask==NULL) {
    fprintf(stderr,"ERROR: les_mask is null\n");
    return 0;
  }
  if(wm_mask==NULL) {
    fprintf(stderr,"ERROR: wm_mask is null\n");
    return 0;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_fill_lesion_with_feature]   img = %s\n",img->fname);
    fprintf(stdout,"[niik_image_fill_lesion_with_feature]   les = %s\n",les_mask->fname);
    fprintf(stdout,"[niik_image_fill_lesion_with_feature]   WM  = %s\n",wm_mask->fname);
  }

  if(verbose) fprintf(stdout,"[niik_image_fill_lesion_with_feature] transform lesion mask\n");
  if((xfmlm = niik_image_affine_transform_3d(les_mask,img,lmat,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
    return 0;
  }

  if(verbose) fprintf(stdout,"[niik_image_fill_lesion_with_feature] transform WM mask\n");
  if((xfmwm = niik_image_affine_transform_3d(wm_mask,img,lmat,NIIK_INTERP_NN))==NULL) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d\n");
    return 0;
  }

  fv = niikvec_init(12);

  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(xfmlm,i)>0) {
      if(!niik_image_feature_voxel_type1(img,i,kernel,fv)) {
        fprintf(stderr,"ERROR: niik_image_feature_voxel_type1\n");
        return 0;
      }
    }
  }

  if(verbose>1) fprintf(stdout,"[niik_image_fill_lesion_with_feature] free variables\n");
  xfmlm = niik_image_free(xfmlm);
  xfmwm = niik_image_free(xfmwm);

  if(verbose>1) fprintf(stdout,"[niik_image_fill_lesion_with_feature] success\n");
  return 1;
} /* niik_image_fill_lesion_with_feature */


#endif /* _FALCON_LESION_FILL_C_ */
