/* Filename:     nifti1_kunio_bseg.c
 * Description:  brain segmentation functions
 * Author:       Kunio Nakamura
 * Date:         February 29, 2012
 *
 * Main functions:
 *      niik_image_bseg_simplest
 *      niik_image_bseg_basic
 *
 * Helpful functions:
 *     niik_image_bseg_basic_thresh
 *     niik_image_brain_mask_from_mni
 *     niik_image_bseg_get_brain_surface
 *
 *
 *
 */



#ifndef _FALCON_BSEG_C_
#define _FALCON_BSEG_C_

#include "falcon.h"
#include "falcon_cortex.h"
#include "falcon_morph.h"
#include "falcon_seg.h"


/************************************************************
 *
 * very basic segmentation methods
 *
 * niik_image_bseg_simplest
 *
 * niik_image_bseg_from_mni
 *
 ************************************************************/

nifti_image *niik_image_bseg_simplest(nifti_image *img,niikmat *regmni,double thresh,int mni_dilate,double radius)
/* the simplest brain segmentation function
 * 1. transform and dilate mni brain mask (dilated mni_dilate-times)
 * 2. threshold
 * 3. mask threshold image with mni mask
 * 4. erosion by radius
 * 5. seed fill
 * 6. dilation by radius
 */
{
  nifti_image
  *outimg,
  *mni_seg;
  char
  *FSLDIR,
  fname[1024],
  fcname[32]="niik_image_bseg_simplest";
  int
  i,iter,
  verbose=1;
  niikmat *invmat=NULL;
  if(verbose>=1) niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }
  /*
   * 1. get / transform / dilate MNI brain mask
   */
  if(verbose>=1) fprintf(stdout,"[%s] 1. get / transform / dilate MNI brain mask\n",fcname);
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv FSLDIR\n",fcname);
    return NULL;
  }
  sprintf(fname,"%s/data/standard/MNI152_T1_2mm_brain_mask.nii.gz",FSLDIR);
  if(verbose) fprintf(stdout,"[%s]    reading mni seg   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  for(iter=0; iter<mni_dilate; iter++) {
    if(verbose) fprintf(stdout,"[%s]    MNI mask dilation %i\n",fcname,iter+1);
    if(!niik_image_morph_3d_mask(mni_seg,NULL,NIIK_MORPH_DILATE,NIIK_MORPH_3D_KERNEL_DIAMOND,3)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_mask\n");
      return NULL;
    }
  } /* dilate few times */
  if(verbose>=1) fprintf(stdout,"[%s]    transform mni mask\n",fcname);
  if(regmni==NULL) {
    invmat=niikmat_identity(4,4);
  } else {
    invmat=niikmat_inverse(regmni);
  }
  if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  invmat=niikmat_free(invmat);
  /*
   * 2. threshold image and mask with dilated MNI brain mask
   */
  if(verbose) fprintf(stdout,"[%s] 2. threshold image\n",fcname);
  if(thresh<0 || niik_check_double_problem(thresh)) {
    thresh=2.2;
    if(!niik_image_bseg_basic_thresh_with_mask(img,mni_seg,&thresh,0)) {
      fprintf(stdout,"niik_image_bseg_basic_thresh_with_mask\n");
      return NULL;
    }
  } /* calculate threshold if undefined */
  if(verbose>=1) fprintf(stdout,"[%s]    threshold %8.3f\n",fcname,thresh);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_threshold(outimg,thresh)) {
    fprintf(stderr,"ERROR: niik_image_thresh\n");
    return NULL;
  }
  if(!niik_image_mask(outimg,mni_seg)) {
    fprintf(stderr,"ERROR: niik_image_mask\n");
    return NULL;
  }
  mni_seg=niik_image_free(mni_seg);
  /*
   * 3. erosion
   */
  if(verbose>=1) fprintf(stdout,"[%s] 3. erode  %8.3f\n",fcname,radius);
  if(!niik_image_morph_3d_radius(outimg,NIIK_MORPH_ERODE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s] 4. seedfill\n",fcname);
  if(!niik_image_seed_fill_from_middle(outimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill_from_middle\n");
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s] 5. dilate %8.3f\n",fcname,radius);
  if(!niik_image_morph_3d_radius(outimg,NIIK_MORPH_DILATE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
    return NULL;
  }
  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(img,i)<thresh)
      niik_image_set_voxel(outimg,i,0);
  }
  if(verbose) fprintf(stdout,"[%s] finished\n",fcname);
  return outimg;
} /* niik_image_bseg_simplest */


nifti_image *niik_image_bseg_from_mni(nifti_image *img,niikmat *regmni,double thresh,double radius)
/* segmentation using (1) MNI registration, (2) threshold, (3) radius (dilate mni mask) and (4) combine masks from (2) and (3) */
{
  nifti_image
  *outimg,
  *mni_seg;
  char
  *FSLDIR,
  fname[1024],
  fcname[512]="niik_image_bseg_from_mni";
  int
  verbose=1;
  niikmat *invmat=NULL;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"[%s] 1. get/transform/dilate mni brain mask\n",fcname);
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"ERROR: please setenv FSLDIR\n");
    return NULL;
  }
  sprintf(fname,"%s/data/standard/MNI152_T1_2mm_brain_mask.nii.gz",FSLDIR);
  if(verbose>=1) fprintf(stdout,"[%s] reading mni seg   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  if(radius>=0) {
    if(verbose>=1) fprintf(stdout,"[%s] MNI mask dilation %-8.3f\n",fcname,radius);
    if(!niik_image_morph_3d_radius(mni_seg,NIIK_MORPH_DILATE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
      return NULL;
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] transform mni mask\n",fcname);
  if(regmni==NULL) {
    invmat=niikmat_identity(4,4);
  } else {
    invmat=niikmat_inverse(regmni);
  }
  if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  invmat=niikmat_free(invmat);
  if(verbose>=1) fprintf(stdout,"[%s] 2. threshold image\n",fcname);
  if(thresh<0 || niik_check_double_problem(thresh)) {
    thresh=1.5;
    if(!niik_image_bseg_basic_thresh_with_mask(img,mni_seg,&thresh,0)) {
      fprintf(stdout,"niik_image_bseg_basic_thresh_with_mask\n");
      return NULL;
    }
  } /* calculate threshold if undefined */
  if(verbose>=1) fprintf(stdout,"[%s]    threshold %8.3f\n",fcname,thresh);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_threshold(outimg,thresh)) {
    fprintf(stderr,"ERROR: niik_image_thresh\n");
    return NULL;
  }
  if(!niik_image_mask(outimg,mni_seg)) {
    fprintf(stderr,"ERROR: niik_image_mask\n");
    return NULL;
  }
  mni_seg=niik_image_free(mni_seg);
  if(verbose>=1) fprintf(stdout,"[%s] finished\n",fcname);
  return outimg;
} /* niik_image_bseg_from_mni */


nifti_image *niik_image_bseg_basic(nifti_image *img,int mni_dilate_iter,double lothresh,double mdthresh,double hithresh,double hithresh2,double radius,niikmat *regmni,int flag_pv)
/* bseg = brain segmentation
 * 1. read/transform MNI masks
 * 2. calculate threshold values:
 *    -lo  absolute threshold
 *    -md  initial threshold (slightly above lothresh to separate from non-brain)
 *    -hi  absolute threshold
 * 3. combine mni_seg and threshold mask
 * 4. conditional opening
 * 5. partial volume correction
 */

{
  nifti_image
  *tmpimg=NULL,
   *mni_fil=NULL,
    *mni_roi_dil=NULL,
     *mni_roi=NULL,
      *mni_seg=NULL;
  const char *NIIKDIR=NULL;
  char
  fname[4096],
        fcname[32]="niik_image_bseg_basic";
  niikmat
  *regmat=NULL,
   *invmat=NULL;
  int
  outidx=1,
  i,j,k,n,nn,
  xdim,area,
  verbose=1;
  unsigned char
  *bimg=NULL,*broi=NULL;
  double
  dval,dran,
       *affpar=NULL,
        *dimg=NULL;

  if(verbose>=1) {
    fprintf(stdout,"[%s] approximate brain segmentaiton\n",fcname);
    niik_fc_display(fcname,1);
  }
  /* check environment and inputs */
  if(!(NIIKDIR=get_NIIKDIR()))
    return NULL;

  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return NULL;
  }

  /* get affine MNI registration */
  if(regmni!=NULL) {} /* regmni exists */
  else if(img->sform_code!=NIFTI_XFORM_MNI_152) {
    affpar=(double *)calloc(25,sizeof(double));
    affpar[7]=affpar[8]=affpar[9]=affpar[10]=1;
    if(!niik_aregister_align_mni_default1(img,NULL,affpar)) {
      fprintf(stderr,"ERROR: niik_aregister_align_mni_default1\n");
      return NULL;
    }
    if((invmat = niik_aregister_matrix_from_affpar(affpar))==NULL) {
      fprintf(stderr,"ERROR: niik_aregister_matrix_from_affpar\n");
      return NULL;
    }
    if(!niikmat_inverse_update(invmat)) {
      fprintf(stderr,"ERROR: niikmat_inverse\n");
      return NULL;
    }
    free(affpar);
  } else {
    if((invmat=niik_image_bseg_get_mni_affine_matrix_from_sform(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_bseg_get_mni_affine_matrix_from_sform\n");
      return NULL;
    }
  }

  if((dimg=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return NULL;
  }

  /************************************************
   * 1. read / transform MNI brain mask
   ************************************************/

  fprintf(stdout,"[%s] 1.   read and transform MNI images\n",fcname);
  fprintf(stdout,"[%s] 1a.  read MNI images\n",fcname);
  if(mni_seg==NULL) {
    sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_brain_mask.nii.gz",NIIKDIR);
    fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
    if((mni_seg=nifti_image_read(fname,1))==NULL) {
      fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
      return NULL;
    }
  }
  sprintf(fname,"%s/data/CLADA/MNI152_T1_1mm_brain_edge_ROI_3-9mm_dilate3mm_fil.nii.gz",NIIKDIR);
  fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
  if((mni_roi=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  sprintf(fname,"%s/data/CLADA/MNI152_T1_1mm_brain_edge_ROI_nonbrain_6mm_blur.nii.gz",NIIKDIR);
  fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
  if((mni_fil=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }

  if(invmat==NULL) {
    if((invmat=niikmat_inverse(regmni))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse\n");
      return NULL;
    }
  }

  /* transform mni masks */
  fprintf(stdout,"[%s] 1b.  transform MNI images\n",fcname);
  if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  if(!niik_image_affine_transform_3d_update(mni_roi,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  if(!niik_image_affine_transform_3d_update(mni_fil,img,invmat,NIIK_INTERP_LINEAR)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  mni_seg->scl_slope=0;
  mni_seg->scl_inter=0;

  if(!niik_image_type_convert(mni_roi,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(mni_roi,NIFTI_TYPE_UINT8)\n");
    return NULL;
  }
  broi=mni_roi->data;

  if(verbose>=3) fprintf(stdout,"    mask count = %i \n",niik_image_count_mask(mni_seg));
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_mniseg.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
    sprintf(fname,"tmp_bseg_mask%i_mniroi.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_roi);
    sprintf(fname,"tmp_bseg_mask%i_mnifil.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_fil);
  }

  /*
   * 2. THRESHOLD USING RIDLER or INPUT
   */
  fprintf(stdout,"[%s] 2.   threshold calculation\n",fcname);
  fprintf(stdout,"[%s] 2a.  lower threshold\n",fcname);
  if(niik_check_double_problem(lothresh)) {
    lothresh=1.5;
    if(!niik_image_bseg_basic_thresh(img,invmat,&lothresh,2)) {
      fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
      return NULL;
    }
  }
  /*if(!niik_image_thresh_ridler(img,mni_seg,&lothresh)){
    fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
    return NULL; } }*/
  fprintf(stdout,"[%s] 2b.  upper threshold\n",fcname);
  if(niik_check_double_problem(hithresh)) {
    hithresh =  niik_image_get_upper_threshold(img,mni_seg,0.0001);
    if(fabs(hithresh)<0.1) {
      fprintf(stderr,"[%s] ERROR: hi threshold is too low\n",fcname);
      return NULL;
    }
  }
  if(niik_check_double_problem(hithresh)) {
    fprintf(stderr,"ERROR: hithresh is not normal, %f\n",hithresh);
    return NULL;
  }
  if(niik_check_double_problem(hithresh2)) {
    hithresh2 =  niik_image_get_upper_threshold(img,mni_seg,0.05);
  }
  if(niik_check_double_problem(hithresh2)) {
    fprintf(stderr,"ERROR: hithresh is not normal, %f\n",hithresh);
    return NULL;
  }

  fprintf(stdout,"[%s] 2c.  medium threshold\n",fcname);
  /* 2012-06-02, Kunio
   * -the objective of the medium threshold is to use the threshold
   *  that is slightly higher than 'threshold' and to make sure
   *  that brain and non-brain is disconnected.
   * -original method is Ridler's threshold with 0.8 factor
   * --this method sometimes resulted in high values...
   * -next method is the trimmed median or trimmed ridler
   *
   * 2012-10-29, Kunio
   * -currently this is disabled by 'mdthresh=lothresh;'
   */
  if(!niik_check_double_problem(mdthresh)) {}
  else {
    mdthresh=1.4;
    if(!niik_image_thresh_ridler(img,mni_seg,&mdthresh)) {
      fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
      return NULL;
    }
  } /* calculate mdthresh */
  if(niik_check_double_problem(mdthresh)) {
    fprintf(stderr,"ERROR: mdthresh is not normal, %f\n",mdthresh);
    return NULL;
  }
  mdthresh=lothresh;

  if(mni_dilate_iter>0) {
    fprintf(stdout,"[%s] Dilation: MNI mask   x %i\n",fcname,mni_dilate_iter);
    mni_seg=niik_image_free(mni_seg);
    sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_brain_mask.nii.gz",NIIKDIR);
    fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
    if((mni_seg=nifti_image_read(fname,1))==NULL) {
      fprintf(stderr,"[%s] ERROR: nifti_image_read %s\n",fcname,fname);
      return NULL;
    }
    for(i=0; i<mni_dilate_iter; i++) {
      if(!niik_image_morph_3d_mask(mni_seg,NULL,NIIK_MORPH_DILATE,NIIK_MORPH_3D_KERNEL_DIAMOND,3)) {
        fprintf(stderr,"[%s] ERROR: niik_image_morph_3d_mask\n",fcname);
        return NULL;
      }
    }
    if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
      fprintf(stderr,"[%s] ERROR: niik_image_affine_transform_3d_update\n",fcname);
      return NULL;
    }
  }

  if(verbose>=3) fprintf(stdout,"    mask count = %i \n",niik_image_count_mask(mni_seg));
  fprintf(stdout,"    threshold         %8.3f\n",lothresh);
  fprintf(stdout,"    initial threshold %8.3f\n",mdthresh);
  fprintf(stdout,"    upper threshold   %8.3f (global) \n",hithresh);
  fprintf(stdout,"    upper threshold   %8.3f (regional) \n",hithresh2);

  for(i=0; i<img->nvox; i++) {
    if(dimg[i]<mdthresh) { /* initial threshold */
      niik_image_set_voxel(mni_seg,i,0);
    } else if(broi[i]>0) {
      if(dimg[i]>hithresh2) /* regional */
        niik_image_set_voxel(mni_seg,i,0);
    } else {
      if(dimg[i]>hithresh)  /* global */
        niik_image_set_voxel(mni_seg,i,0);
    }
  }

  if(verbose>=3) fprintf(stdout,"    mask count = %i \n",niik_image_count_mask(mni_seg));
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_thresh.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
  }

  /*
   * 3. REMOVE BRIGHT AREA
   */
  fprintf(stdout,"[%s] 3.   dilate bright area\n",fcname);
  fprintf(stdout,"    using %8.3f regional threshold\n",hithresh2);
  mni_roi_dil = niik_image_copy(mni_roi);
  niik_image_clear(mni_roi_dil);
  for(i=0; i<img->nvox; i++) {
    if(broi[i]>0) {
      if(dimg[i]>hithresh2) /* regional */
        niik_image_set_voxel(mni_roi_dil,i,1);
    } else {
      if(dimg[i]>hithresh) /* global */
        niik_image_set_voxel(mni_roi_dil,i,1);
    }
  }
  niik_image_morph_3d_radius(mni_roi_dil,NIIK_MORPH_CLOSE,3.2);
  niik_image_morph_3d_radius(mni_roi_dil,NIIK_MORPH_DILATE,2.4);
  if(verbose>=3) fprintf(stdout,"    mask count = %i \n",niik_image_count_mask(mni_seg));
  niik_image_maskout(mni_seg,mni_roi_dil);
  if(verbose>=3) fprintf(stdout,"    mask count = %i \n",niik_image_count_mask(mni_seg));
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_remove_bright.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
    sprintf(fname,"tmp_bseg_mask%i_remove_roi.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_roi_dil);
  }

  /* seed fill */
  if(regmni==NULL) {
    if((regmat=niikmat_inverse(invmat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse\n");
      return NULL;
    }
  } else {
    if((regmat=niikmat_copy(regmni))==NULL) {
      fprintf(stderr,"ERROR: niikmat_copy\n");
      return NULL;
    }
  }
  if(!niik_image_bseg_seed_fill_from_eroded_standard_mask(mni_seg,regmat,-10)) {
    fprintf(stderr,"ERROR: niik_image_bseg_seed_fill_from_eroded_standard_mask\n");
    return NULL;
  }

  /*
   * 4. OPENING MNI MASK
   *    erosion
   *    seed fill
   *    dilation
   *    threshold
   */

  if(niik_check_double_problem(radius)) radius = 3.6; /* default value */
  fprintf(stdout,"[%s] 4.   conditional regional opening %-8.3f\n",fcname,radius);
  fprintf(stdout,"[%s] 4a.  regional radius map\n",fcname);
  niik_image_type_convert(mni_fil,NIFTI_TYPE_FLOAT32);
  for(i=0; i<img->nvox; i++) {
    niik_image_mul_voxel(mni_fil,i,radius);
  }
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_Rmap.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_fil);
  }

  fprintf(stdout,"[%s] 4b.  regional erosion %6.3f\n",fcname,radius);
  if(!niik_image_morph_3d_radius_map(mni_seg,mni_fil,NIIK_MORPH_ERODE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return NULL;
  }
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_erode.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
  }

  fprintf(stdout,"[%s] 4c.  seed-fill\n",fcname);
  if(!niik_image_bseg_seed_fill_from_eroded_standard_mask(mni_seg,regmat,-10)) {
    fprintf(stderr,"ERROR: niik_image_bseg_seed_fill_from_eroded_standard_mask\n");
    return NULL;
  }

  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_seedfill.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
  }

  fprintf(stdout,"[%s] 4d.  regional dilation %6.3f\n",fcname,radius);
  if(!niik_image_morph_3d_radius_map(mni_seg,mni_fil,NIIK_MORPH_DILATE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return NULL;
  }
  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_dilate.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
  }

  fprintf(stdout,"[%s] 4e.  apply thresholds (%6.2f %6.2f)\n",fcname,lothresh,hithresh2);
  for(i=0; i<img->nvox; i++) {
    if(dimg[i]<lothresh)   /* lower thresh */
      niik_image_set_voxel(mni_seg,i,0);
    if(broi[i]>0) {
      if(dimg[i]>hithresh2) /* regional upper thresh */
        niik_image_set_voxel(mni_seg,i,0);
    } else {
      if(dimg[i]>hithresh) /* global upper thresh */
        niik_image_set_voxel(mni_seg,i,0);
    }
    if(niik_image_get_voxel(mni_roi_dil,i))
      niik_image_set_voxel(mni_seg,i,0);
  }
  niik_image_maskout(mni_seg,mni_roi_dil);

  if(verbose>1) {
    sprintf(fname,"tmp_bseg_mask%i_conop.nii.gz",outidx++);
    fprintf(stdout,"      writing tmpseg: %s\n",fname);
    niik_image_write(fname,mni_seg);
  }

  /* free filter map */
  mni_fil = niik_image_free(mni_fil);

  if(flag_pv) {
    tmpimg=niik_image_copy(mni_seg);
    bimg=mni_seg->data;
    broi=tmpimg->data;
    xdim=mni_seg->nx;
    area=xdim*mni_seg->ny;
    dran=lothresh/2.0;
    for(k=n=0; k<mni_seg->nz; k++) {
      for(j=0; j<mni_seg->ny; j++) {
        for(i=0; i<mni_seg->nx; n++,i++) {
          if(bimg[n]>0) {
            bimg[n]=255;
            continue;
          }
          nn=0;
          if(i>0) if(broi[n-   1]>0) nn++;
          if(j>0) if(broi[n-xdim]>0) nn++;
          if(k>0) if(broi[n-area]>0) nn++;
          if(i<mni_seg->nx-1) if(broi[n+   1]>0) nn++;
          if(j<mni_seg->ny-1) if(broi[n+xdim]>0) nn++;
          if(k<mni_seg->nz-1) if(broi[n+area]>0) nn++;
          if(nn==0) continue;
          dval = niik_image_get_voxel(img,n);
          if(dval<dran)     bimg[n]=0;
          else if(dval>lothresh) bimg[n]=255;
          else {
            bimg[n] = NIIK_UCMINMAX(floor(255*(dval-lothresh)/dran+0.5),0,255);
          }
        }
      }
    }
    fprintf(stdout,"[%s] bseg volume %8.2f mm with    (parial volume)\n",fcname,niik_image_get_voxel_size(mni_seg)*niik_image_get_sum(mni_seg,mni_seg)*0.001/255.0);
  }

  else {
    fprintf(stdout,"[%s] bseg volume %8.2f mm\n",fcname,niik_image_get_voxel_size(mni_seg)*niik_image_count_mask(mni_seg)*0.001);
  }

  invmat=niikmat_free(invmat);
  regmat=niikmat_free(regmat);
  nifti_image_free(mni_roi_dil);
  nifti_image_free(mni_roi);
  free(dimg);
  return mni_seg;
} /* niik_image_bseg_basic */


/******************************************************************************
 *
 * helpful functions
 *
 *
 *   niikmat *niik_image_bseg_get_mni_affine_matrix_from_sform(nifti_image *img);
 *
 *
 ******************************************************************************/

int niik_image_bseg_basic_thresh(nifti_image *img,niikmat *regmni,double *thresh,int method)
/* -estimate a threshold for brain mask
 * -thresh should be the factor for ridler threshold function [default = 1.5 - 2.2 or something like that] */
{
  nifti_image
  *mni_seg;
  const char *NIIKDIR=NULL;
  char
  fname[1024],
        *FSLDIR=NULL,
         fcname[512]="niik_image_bseg_basic_thresh";
  int verbose=0;
  if(verbose) fprintf(stdout,"[%s] start %s\n",fcname,img->fname);
  if((FSLDIR=getenv("FSLDIR"))==NULL) {
    fprintf(stderr,"ERROR: please setenv FSLDIR\n");
    return 0;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    return 0;
  }
  /*sprintf(fname,"%s/data/standard/MNI152_T1_2mm_brain_mask.nii.gz",FSLDIR);*/
  sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.25mm_dilate.nii.gz",NIIKDIR);
  if(verbose) fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  if(!niik_image_affine_transform_3d_update(mni_seg,img,regmni,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  if(!niik_image_bseg_basic_thresh_with_mask(img,mni_seg,thresh,method)) {
    fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh_with_mask\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"[%s] finish\n",fcname);
  mni_seg=niik_image_free(mni_seg);
  return 1;
}


int niik_image_bseg_basic_thresh_with_mask(nifti_image *img,nifti_image *mni_seg,double *thresh,int method)
/* -estimate a threshold for brain mask
 * -mni_seg is the transformed mask (same space as img)
 * -input thresh should be ridler threshold factor [default = 2.2]
 * -method is (1) otsu, (2) ridler within 10-90 percentile,
 *  and (otherwise) riddler otherwise
 * */
{
  nifti_image *maskimg=NULL;
  char fcname[512]="niik_image_bseg_basic_thresh_with_mask";
  double p[2];
  int i;
  int verbose=1;
  if(verbose) fprintf(stdout,"[%s] start %s\n",fcname,img->fname);
  switch(method) {
  default:
    if(verbose>0) fprintf(stdout,"[%s] ridler threshold\n",fcname);
    if(!niik_image_thresh_ridler(img,mni_seg,thresh)) {
      fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
      return 0;
    }
    break;
  case 1:
    if(verbose>0) fprintf(stdout,"[%s] otsu threshold\n",fcname);
    if(!niik_image_thresh_otsu(img,mni_seg,thresh)) {
      fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
      return 0;
    }
    break;
  case 2:
    if(verbose>0) fprintf(stdout,"[%s] ridler threshold 0.1-0.9\n",fcname);
    p[0] = niik_image_get_percentile(img,mni_seg,0.1);
    p[1] = niik_image_get_percentile(img,mni_seg,0.9);
    if(verbose>0) fprintf(stdout,"[%s]   0.1 - 0.9 = %9.3f %9.3f\n",fcname,p[0],p[1]);
    if((maskimg=niik_image_copy(mni_seg))==NULL) {
      fprintf(stderr,"[%s] niik_image_copy\n",fcname);
      return 0;
    }
    for(i=0; i<img->nvox; i++) {
      if     (niik_image_get_voxel(img,i)<p[0]) niik_image_set_voxel(maskimg,i,0);
      else if(niik_image_get_voxel(img,i)>p[1]) niik_image_set_voxel(maskimg,i,0);
    }
    if(!niik_image_thresh_ridler(img,maskimg,thresh)) {
      fprintf(stderr,"ERROR: niik_image_thresh_otsu \n");
      return 0;
    }
    maskimg=niik_image_free(maskimg);
    break;
  }
  if(verbose) fprintf(stdout,"[%s]   thresh=%.4f\n",fcname,*thresh);
  return 1;
}


nifti_image *niik_image_brain_mask_from_mni(nifti_image *img,niikmat *regmni,int radius,int verbose)
/* read a pre-defined mask and affine transform */
{
  niikmat *invmat=NULL;
  const char *NIIKDIR=NULL;
  char fname[4096],fcname[512]="niik_image_brain_mask_from_mni";
  nifti_image *mni_seg=NULL;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(regmni==NULL) {
    fprintf(stderr,"ERROR: regmni is null\n");
    return NULL;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    return NULL;
  }
  if(radius==0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.nii.gz",NIIKDIR);
  } else if(radius>0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_dilate.nii.gz",NIIKDIR,radius);
  } else if(radius<0) {
    radius=abs(radius);
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_erode.nii.gz",NIIKDIR,radius);
  }
  if(verbose>=1) fprintf(stdout,"[%s]    reading brain mask   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s]    transform mni mask\n",fcname);
  if(regmni==NULL) {
    invmat=niikmat_identity(4,4);
  } else {
    invmat=niikmat_inverse(regmni);
  }
  if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  niikmat_free(invmat);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return mni_seg;
}


int niik_image_bseg_conditional_opening(nifti_image *maskimg,double radius)
/* conditional opening
 * -erode by radius
 * -seed fill
 * -dilate by radius */
{
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return 0;
  }
  if(!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_ERODE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return 0;
  }
  if(!niik_image_seed_fill_from_middle(maskimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill_3d_from_middle\n");
    return 0;
  }
  if(!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_DILATE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return 0;
  }
  return 1;
}


int niik_image_bseg_conditional_opening_imask(nifti_image *maskimg,nifti_image *imask,double radius)
/* conditional opening
 * -erode by radius
 * -seed fill using imask
 * -dilate by radius */
{
  nifti_image *tmpimg=NULL;
  char fcname[64]="niik_image_bseg_conditional_opening_imask";
  int verbose=0;
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return 0;
  }
  if(imask==NULL) {
    fprintf(stderr,"ERROR: imask is null\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] start %.3f\n",fcname,radius);
  if(!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_ERODE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return 0;
  }
  if((tmpimg=niik_image_copy(imask))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] %i %i\n",fcname,niik_image_count_mask(maskimg),niik_image_count_mask(tmpimg));
  if(!niik_image_seed_fill_slow(maskimg,tmpimg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] %i %i\n",fcname,niik_image_count_mask(maskimg),niik_image_count_mask(tmpimg));
  niik_image_copy_data(tmpimg,maskimg);
  if(!niik_image_morph_3d_radius(maskimg,NIIK_MORPH_DILATE,radius)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] %i\n",fcname,niik_image_count_mask(maskimg));
  tmpimg=niik_image_free(tmpimg);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
}


int niik_image_bseg_seed_fill_from_eroded_standard_mask(nifti_image *maskimg,niikmat *regmat,int erode_radius)
/* applies seed fill from eroded standard brain mask
 * -regmat is from maskimg to standard image
 */
{
  nifti_image *tmpimg=NULL,*mni_seg=NULL;
  const char *NIIKDIR=NULL;
  char fname[4096],fcname[512]="niik_image_bseg_seed_fill_from_inside";
  int verbose=1;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(regmat==NULL) {
    fprintf(stderr,"ERROR: regmat is null\n");
    return 0;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    return 0;
  }
  erode_radius=-fabs(erode_radius);
  if(erode_radius==0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.nii.gz",NIIKDIR);
  } else if(erode_radius>0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_dilate.nii.gz",NIIKDIR,erode_radius);
  } else if(erode_radius<0) {
    erode_radius=fabs(erode_radius);
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_erode.nii.gz",NIIKDIR,erode_radius);
  }
  if(verbose>=1) fprintf(stdout,"[%s]    reading brain mask   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s]    transform mni mask\n",fcname);
  if(!niik_image_inverse_affine_transform_3d_update(mni_seg,maskimg,regmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  if(verbose>=2) {
    fprintf(stdout,"[%s] writing tmp_bseg_seed0.nii.gz\n",fcname);
    niik_image_write("tmp_bseg_seed0.nii.gz",mni_seg);
  }
  if((tmpimg=niik_image_copy(maskimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] seed fill\n",fcname);
  if(!niik_image_seed_fill_slow(maskimg,mni_seg,0)) {
    fprintf(stderr,"ERROR: niik_image_seed_fill_slow\n");
    return 0;
  }
  if(verbose>=2) {
    fprintf(stdout,"[%s] writing tmp_bseg_seed1.nii.gz\n",fcname);
    niik_image_write("tmp_bseg_seed1.nii.gz",maskimg);
    fprintf(stdout,"[%s] writing tmp_bseg_seed2.nii.gz\n",fcname);
    niik_image_write("tmp_bseg_seed2.nii.gz",mni_seg);
  }
  if(verbose>=1) fprintf(stdout,"[%s] update mask image\n",fcname);
  niik_image_copy_data(mni_seg,maskimg);
  niik_image_mask(maskimg,tmpimg);
  tmpimg=niik_image_free(tmpimg);
  mni_seg=niik_image_free(mni_seg);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
} /* niik_image_bseg_seed_fill_from_eroded_standard_mask */


niikmat *niik_image_bseg_get_mni_affine_matrix_from_sform(nifti_image *img) {
  const char *NIIKDIR=NULL;
  char fname[4096],fcname[128]="niik_image_bseg_get_mni_affine_matrix_from_sform";
  int verbose=2;
  nifti_image *mni_seg=NULL;
  niikmat *invmat=NULL;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_brain_mask.nii.gz",NIIKDIR);
  if(verbose>=2) fprintf(stdout,"[%s]   reading %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  fprintf(stdout,"[%s] img %s sform\n",fcname,img->fname);
  niikmat_display(niikmat_mat44_matrix(img->sto_xyz));
  invmat=niikmat_mat44_matrix(img->sto_xyz);
  niikmat_multiply_mat1_free2(invmat,niikmat_scale_matrix(1.0/img->dx,1.0/img->dy,1.0/img->dz));
  niikmat_multiply_mat2_free1(niikmat_mat44_matrix(mni_seg->sto_ijk),invmat);
  niikmat_multiply_mat2_free1(niikmat_scale_matrix(mni_seg->dx,mni_seg->dy,mni_seg->dz),invmat);
  fprintf(stdout,"[%s] img %s mni registration matrix\n",fcname,img->fname);
  if(verbose>=1) niikmat_display(invmat);
  if(!niikmat_inverse_update(invmat)) {
    fprintf(stderr,"ERROR: niikmat_inverse\n");
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return invmat;
}


int niik_image_bseg_remove_eyes(nifti_image *img,nifti_image *maskimg,niikmat *regmni,double radius,double hithresh,int along_brain_edge)
/* img is the t1w image
 * maskimg is the mask for brain
 * regmni is the affine mni registration matrix
 * radius is the dilation radius for eye mask
 *   -if positive, dilation
 *   -if negative, erosion
 * hithresh is the threshold for intensity
 * along_brain_edge changes to another mask (not just the eyes but also brain edges, thus more areas)
 */
{
  nifti_image
  *mnieye=NULL;
  const char *NIIKDIR=NULL;
  char
  fname[4096],
        fcname[512]="niik_image_bseg_remove_eye";
  int verbose=1,num,i;
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"[%s] ERROR: please setenv NIIKDIR\n",fcname);
    return 0;
  }
  /* remove bright eyes */
  if(along_brain_edge) {
    sprintf(fname,"%s/data/CLADA/MNI152_T1_1mm_brain_edge_ROI_3-9mm_dilate3mm_fil.nii.gz",NIIKDIR);
  } else {
    sprintf(fname,"%s/data/bseg/MNI152_T1_3mm_eyes.nii.gz",NIIKDIR);
  }
  if(verbose>=1) fprintf(stdout,"[%s] reading %s\n",fcname,fname);
  if((mnieye=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s]   resample 1x1x1\n",fcname);
  if(!niik_image_resample_3d_update(mnieye,1,1,1,-1,-1,-1,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_resample_3d_update\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s]   radius %.2f\n",fcname,radius);
  if(radius>0) {
    if(!niik_image_morph_3d_radius(mnieye,NIIK_MORPH_DILATE,radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
      return 0;
    }
  } else if(radius<0) {
    if(!niik_image_morph_3d_radius(mnieye,NIIK_MORPH_DILATE,-radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
      return 0;
    }
  }
  if(!niik_image_inverse_affine_transform_3d_update(mnieye,maskimg,regmni,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] maskout eyes   (%i)\n",fcname,niik_image_count_mask(mnieye));
  if(verbose>=2) {
    sprintf(fname,"tmp_bseg_eyes.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,mnieye);
  }
  num=niik_image_count_mask(maskimg);
  if(hithresh<0) {
    for(i=0; i<maskimg->nvox; i++) {
      if(niik_image_get_voxel(mnieye,i)>0) {
        niik_image_set_voxel(maskimg,i,0);
      }
    }
  } else {
    for(i=0; i<maskimg->nvox; i++) {
      if(niik_image_get_voxel(mnieye,i)>0) {
        if(niik_image_get_voxel(img,i)>hithresh) {
          niik_image_set_voxel(maskimg,i,0);
        }
      }
    }
  }
  mnieye=niik_image_free(mnieye);
  if(verbose>=1)
    fprintf(stdout,"[%s] removed %i voxels %.4f ml\n",fcname,
            num-niik_image_count_mask(maskimg),
            (num-niik_image_count_mask(maskimg))*0.001*niik_image_get_voxel_size(maskimg));
  return 1;
}


int niik_image_bseg_get_hithresh(nifti_image *img,niikmat *regmni,double radius,double *hithresh) {
  nifti_image
  *mni_seg=NULL;
  niikmat
  *invmat=NULL;
  const char* NIIKDIR=NULL;
  char
  fname[4096],
        fcname[512]="niik_image_bseg_get_hithresh";
  int verbose=1;
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    return 0;
  }
  /* remove bright eyes */
  sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_erode.nii.gz",NIIKDIR,(int)radius);
  if(verbose>=1) fprintf(stdout,"[%s] reading %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  invmat=niikmat_inverse(regmni);
  if(!niik_image_affine_transform_3d_update(mni_seg,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  invmat=niikmat_free(invmat);
  *hithresh=niik_image_get_percentile(img,mni_seg,*hithresh);
  mni_seg=niik_image_free(mni_seg);
  return 1;
}

/***************************************************************************
 *
 * niik_image_bseg_surface
 *
 *
 ***************************************************************************/

kobj *niik_image_bseg_get_brain_surface(niikmat *regmni,double imove,int verbose) {
  kobj
  *obj=NULL;
  kvert
  *v;
  char fcname[512]="niik_image_bseg_get_brain_surface";
  const char* NIIKDIR=NULL;
  char
  fname[4096];
  niikmat *invmat=NULL;
  int amove=0;
  if(verbose>=1) {
    fprintf(stdout,"[%s] start\n",fcname);
  }
  if(regmni==NULL) {
    fprintf(stderr,"ERROR: regmni is null\n");
    return NULL;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    return NULL;
  }
  if(imove==0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask.2mm.off.gz",NIIKDIR);
  } else if(imove<0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_erode.off.gz",NIIKDIR,(int)-imove);
  } else if(imove>0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_dilate.off.gz",NIIKDIR,(int)imove);
  }
  if(verbose>=1) fprintf(stdout,"[%s]   read %s\n",fcname,fname);
  if((obj=off_kobj_read_offply(fname))==NULL) {
    fprintf(stderr,"ERROR: off_kobj_read_off %s\n",fname);
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask.2mm.off.gz",NIIKDIR);
    amove=1;
    if(verbose>=1) fprintf(stdout,"[%s]   read %s\n",fcname,fname);
    if((obj=off_kobj_read_offply(fname))==NULL) {
      fprintf(stderr,"ERROR: off_kobj_read_off %s\n",fname);
      return NULL;
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] transform\n",fcname);
  invmat=niikmat_inverse(regmni);
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v=niikpt_affine_transform(invmat,v->v);
  }
  invmat=niikmat_free(invmat);
  if(verbose>=1) fprintf(stdout,"[%s] normal calculations\n",fcname);
  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);
  off_smooth_kobj_vert_normal(obj);
  if(amove) {
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v=niikpt_move_normal(v->v,v->normal,imove);
    }
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_smooth_kobj_vert_normal(obj);
  }
  if(verbose>=1) fprintf(stdout,"[%s] finished\n",fcname);
  return obj;
}



/****************************************************************************************
 *
 * surface/volume version
 *
 ****************************************************************************************/

nifti_image *niik_image_bseg_test1(nifti_image *img,niikmat *regmni,kobj *obj,
                                   double thresh,
                                   double lolim,double hilim,
                                   double radius,
                                   double iradius,
                                   double mni_dilate_radius,
                                   int maxiter,
                                   int using_surface,
                                   double close_brain_radius,
                                   double radius_for_lim_erosion,
                                   double curvature_thresh,double curvature_range)
/* -img is the t1w image like t1g or t1p
 * -regmni is the affine matrix from img to mni space
 * -obj is the brain surface object in img space
 * -thresh is the threshold; if negative, used to calculate the threshold (using -thresh as the factor)
 *    for example -1.5
 * -lolim is the lower limit; if negative, used to calculate the lower limit (using -lolim as the factor)
 *    for example -1.2
 * -hilim is the higher limit; same as above for negative number
 *    for example -0.99
 * -radius is the conditional opening radius
 *    for example 3.2
 * -iradius is the initial dilation radius for MNI mask
 *    for example 1.5
 * -mni_dilate_radius is the dilation radius for MNI mask used to mask at each dilation
 * --the dilated mask must include all brain voxels
 *   for example 12
 * -maxiter is the number of iterations
 *   for example 10
 * -using_surface is the flag to use the object, otherwise skips the surface-based analysis
 *  and uses brain-closing to fill in the sulcal and ventricular CSF
 *   for example 10
 * -close_brain_radius is the radius for closing
 * --it is only used when using_surface is zero
 *   for example 6.5
 * -radius_for_lim_erosion is the erosion radius for calculating the lower/upper limits
 *   for example, 2.0 - 5.0
 * -curvature_thresh is the threshold for curvature-based surface deformation
 *   for example, 3.0 - 5.0
 * -curvature_range is the range for threshold in curvature-based surface deformation
 *   for example, 1.0 - 2.0
 */
{
  kvert *v;
  nifti_image
  *mnieye=NULL,
   *mniroi=NULL,
    *bounds=NULL,
     *tmpimg=NULL,
      *outimg=NULL;
  const char *NIIKDIR=NULL;
  char
  fname[512],
        fcname[512]="niik_image_bseg_test1";
  niikvec *cm=NULL;
  niikmat *invmat;
  double
  ym,
  rgb[9];
  int
  idx[9],
      vindex,
      verbose=1,
      num,
      i,j;

  if(verbose>=1) {
    fprintf(stdout,"[%s] start\n",fcname);
    fprintf(stdout,"  image              : %s\n",img->fname);
    if(obj!=NULL)
      fprintf(stdout,"  surface object     : %s\n",obj->fname);
    if(thresh>0)
      fprintf(stdout,"  threshold          : %.3f\n",thresh);
    else
      fprintf(stdout,"  thresh factor      : %.3f\n",-thresh);
    if(lolim>0)
      fprintf(stdout,"  lower limit        : %.3f\n",lolim);
    else
      fprintf(stdout,"  lower factor       : %.3f\n",-lolim);
    if(hilim>0)
      fprintf(stdout,"  upper limit        : %.3f\n",hilim);
    else
      fprintf(stdout,"  upper factor       : %.3f\n",-hilim);
    fprintf(stdout,"  opening radius     : %.3f\n",radius);
    fprintf(stdout,"  initial dilation   : %.3f\n",iradius);
    fprintf(stdout,"  mni_dilate_radius  : %.3f\n",mni_dilate_radius);
    fprintf(stdout,"  # iterations       : %i\n",maxiter);
    if(!using_surface)
      fprintf(stdout,"  close brain radius : %.3f\n",close_brain_radius);
    fprintf(stdout,"  lim erosion radius : %.3f\n",radius_for_lim_erosion);
    if(using_surface)
      fprintf(stdout,"  curvature thresh   : %.3f +/- %.2f\n",curvature_thresh,curvature_range);
  }

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(regmni==NULL) {
    fprintf(stderr,"ERROR: regmni is null\n");
    return NULL;
  }
  if(using_surface) {
    if(obj==NULL) {
      fprintf(stderr,"ERROR: obj is null\n");
      return NULL;
    }
    cm=niikvec_init(obj->nvert);
  }

  if(verbose>=3) {
    if(obj!=NULL) {
      off_kobj_add_one_color(obj,0.7,0.4,0.25);
      sprintf(fname,"tmp_bseg0.off.gz");
      fprintf(stdout,"\twriting %s\n",fname);
      off_kobj_write_offply(fname,obj,0);
    }
  }

  /* calculate threshold */
  if(thresh<0) {
    thresh=-thresh;
    if(!niik_image_bseg_basic_thresh(img,regmni,&thresh,0)) {
      fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh\n");
      return NULL;
    }
  }
  fprintf(stdout,"[%s]   thresh = %.3f\n",fcname,thresh);

  if(verbose>=3) {
    bounds=niik_image_copy(img);
  }

  /* initial estimate of brain mask */
  if(verbose>=1) fprintf(stdout,"[%s] get the brain mask R=%.2f\n",fcname,iradius);
  if((outimg=niik_image_bseg_from_mni(img,regmni,thresh,iradius))==NULL) {
    fprintf(stderr,"ERROR: niik_image_bseg_from_mni\n");
    return NULL;
  }
  if((tmpimg=niik_image_copy(outimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }

  /* calculate lolim and hilim
   * -intensity range for eroded areas
   */
  if(!niik_image_morph_3d_radius(tmpimg,NIIK_MORPH_ERODE,radius_for_lim_erosion)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius\n");
    return NULL;
  }
  /*niik_image_display_stats(img,tmpimg);*/
  if(lolim<0) {
    if(verbose>=2) fprintf(stdout,"[%s]   lo lim\n",fcname);
    lolim=-lolim;
    if(!niik_image_bseg_basic_thresh(img,regmni,&lolim,0)) {
      fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh\n");
      return NULL;
    }
  }
  if(hilim<0) {
    if(verbose>=2) fprintf(stdout,"[%s]   hi lim\n",fcname);
    hilim=-hilim;
    hilim=niik_image_get_percentile(img,tmpimg,hilim);
  }
  if(verbose>=1) fprintf(stdout,"[%s]   expand-threshold = %8.1f %8.1f\n",fcname,lolim,hilim);

  /* get eye mask */
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    return NULL;
  }

  /* remove bright eyes */
  sprintf(fname,"%s/data/bseg/MNI152_T1_3mm_eyes.nii.gz",NIIKDIR);
  if(verbose>=1) fprintf(stdout,"[%s] reading %s\n",fcname,fname);
  if((mnieye=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  invmat=niikmat_inverse(regmni);
  if(!niik_image_affine_transform_3d_update(mnieye,img,invmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  invmat=niikmat_free(invmat);
  if(verbose>=1) fprintf(stdout,"[%s] maskout eyes   (%i)\n",fcname,niik_image_count_mask(mnieye));
  if(verbose>=2) {
    sprintf(fname,"tmp_bseg_eyes.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,mnieye);
  }
  num=niik_image_count_mask(outimg);
  for(i=0; i<img->nvox; i++) {
    if(niik_image_get_voxel(img,i)>hilim) {
      if(niik_image_get_voxel(mnieye,i)>0) {
        niik_image_set_voxel(outimg,i,0);
      }
    }
  }
  mnieye=niik_image_free(mnieye);
  if(verbose>=1) fprintf(stdout,"[%s] removed %i voxels %.4f ml\n",fcname,
                           num-niik_image_count_mask(outimg),
                           (num-niik_image_count_mask(outimg))*0.001*niik_image_get_voxel_size(outimg));

  /* mni roi */
  if(verbose>=2) fprintf(stdout,"[%s]   get dilated mask  %i\n",fcname,(int)mni_dilate_radius);
  if((mniroi=niik_image_brain_mask_from_mni(img,regmni,(int)floor(0.5+mni_dilate_radius),1))==NULL) {
    fprintf(stderr,"ERROR: niik_image_brain_mask_from_mni\n");
    return NULL;
  }

  if(verbose>=1) fprintf(stdout,"[%s] conditional opening\n",fcname);
  if(!niik_image_bseg_conditional_opening(outimg,radius)) {
    fprintf(stderr,"ERROR: niik_image_bseg_conditional_opening\n");
    return NULL;
  }

  if(0) {
    sprintf(fname,"tmp_bseg1_imask1.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,outimg);
  }


  /* surface deformation */
  if(using_surface) {
    fprintf(stdout,"[%s] surface deformation with curvature-based smoothing\n",fcname);
    for(i=0; i<40; i++) {
      niik_off_curvature_map_update(obj,cm,2);
      for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
        niik_image_interpolate_3d_linear_xyz_update(outimg,v->v,&ym);
        if(ym>0.5) {
          cm->v[vindex]=0.2*(1.0-NIIK_Heaviside(fabs(cm->v[vindex])-curvature_thresh,curvature_range));
        } else if(ym<0.5) {
          cm->v[vindex]=-0.2*(1.0-NIIK_Heaviside(fabs(cm->v[vindex])-curvature_thresh,curvature_range));
        }
      }
      for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
        v->v=niikpt_move_normal(v->v,v->normal,cm->v[vindex]);
      }
      niik_off_apply_surface_smoothing(obj,2,0.9);
    }
    niik_off_curvature_map_update(obj,cm,2);
  }

  if(verbose>=3) {
    rgb[0]=hilim*2;
    rgb[1]=rgb[2]=0;
    niik_image_boundary(bounds,outimg,rgb,0,1,0);
    sprintf(fname,"tmp_bounds0.nii.gz");
    fprintf(stdout,"\twriting %s\n",fname);
    niik_image_write(fname,bounds);
    if(obj!=NULL && using_surface) {
      off_kobj_add_one_color(obj,0.77,0.8,0.25);
      off_kobj_apply_color_map(obj,cm->v,-5,5,NIIK_COLORMAP_SPECTRAL);
      sprintf(fname,"tmp_bseg1.off.gz");
      fprintf(stdout,"\twriting %s\n",fname);
      off_kobj_write_offply(fname,obj,0);
    }
  }

  /* MAIN LOOP */
  for(j=0; j<maxiter; j++) {
    if(verbose>=2) fprintf(stdout,"[%s] iteration %3i\n",fcname,j+1);
    if(!niik_image_copy_data(outimg,tmpimg)) {
      fprintf(stderr,"ERROR: niik_image_copy_data\n");
      return NULL;
    }
    if(!niik_image_morph_3d_mask(outimg,NULL,NIIK_MORPH_DILATE,NIIK_MORPH_3D_KERNEL_DIAMOND,3)) {
      fprintf(stderr,"ERRO:R niik_image_morph_3d_mask\n");
      return NULL;
    }
    for(i=0; i<img->nvox; i++) {
      if(niik_image_get_voxel(img,i)<thresh) niik_image_set_voxel(outimg,i,0);
      if(niik_image_get_voxel(tmpimg,i)==0) {
        if(niik_image_get_voxel(img,i)<lolim)
          niik_image_set_voxel(outimg,i,0);
        else if(niik_image_get_voxel(img,i)>hilim)
          niik_image_set_voxel(outimg,i,0);
      }
    }
    niik_image_mask(outimg,mniroi);
    if(verbose>=2) fprintf(stdout,"[%s] %4i conditional opening\n",fcname,j+1);
    if(!niik_image_bseg_conditional_opening(outimg,radius)) {
      fprintf(stderr,"ERROR: niik_image_bseg_conditional_opening\n");
      return NULL;
    }
    if(verbose>=1) fprintf(stdout,"[%s] %4i volume %9.4f\n",fcname,j+1,niik_image_get_mask_vol(outimg));
    if(verbose>=3) {
      sprintf(fname,"tmp_bseg%i.nii.gz",j+1);
      fprintf(stdout,"\twriting %s\n",fname);
      niik_image_write(fname,outimg);
    }
  }

  /* surface deformation */
  if(using_surface) {
    fprintf(stdout,"[%s] surface deformation\n",fcname);
    for(i=0; i<40; i++) {
      niik_off_curvature_map_update(obj,cm,2);
      for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
        niik_image_interpolate_3d_linear_xyz_update(outimg,v->v,&ym);
        if(ym>0.5) {
          cm->v[vindex]=0.2*(1.0-NIIK_Heaviside(fabs(cm->v[vindex])-curvature_thresh,curvature_range));
        } else if(ym<0.5) {
          cm->v[vindex]=-0.2*(1.0-NIIK_Heaviside(fabs(cm->v[vindex])-curvature_thresh,curvature_range));
        }
      }
      for(v=obj->vert,vindex=0; v!=NULL; v=v->next,vindex++) {
        v->v=niikpt_move_normal(v->v,v->normal,cm->v[vindex]);
      }
      if(i<35) niik_off_apply_surface_smoothing(obj,2,0.9);
    }
    niik_off_curvature_map_update(obj,cm,2);
  }

  if(verbose>=3) {
    sprintf(fname,"tmp_bseg%i.nii.gz",maxiter+1);
    fprintf(stdout,"\twriting %s\n",fname);
    niik_image_write(fname,outimg);
    if(obj!=NULL && using_surface) {
      off_kobj_apply_color_map(obj,cm->v,-5,5,NIIK_COLORMAP_SPECTRAL);
      sprintf(fname,"tmp_bseg2.off.gz");
      fprintf(stdout,"\twriting %s\n",fname);
      off_kobj_write_offply(fname,obj,0);
    }
  }

  /* final mask */
  if(using_surface) {
    fprintf(stdout,"[%s] create mask form surface object\n",fcname);
    niik_image_clear(outimg);
    off_obj2img(outimg,obj,1);
    for(i=0; i<20; i++) {
      for(v=obj->vert; v!=NULL; v=v->next) {
        v->v=niikpt_move_normal(v->v,v->normal,-0.1);
      }
      off_obj2img(outimg,obj,1);
    }
    niik_image_flip_bgfg(outimg);
    idx[1] = idx[2] = idx[3] = 0;
    if(!niik_image_seed_fill_xyz(outimg,idx,0)) {
      fprintf(stderr,"ERROR: niik_image_seed_fill_xyz\n");
      return NULL;
    }
    niik_image_flip_bgfg(outimg);
  } else {
    fprintf(stdout,"[%s] create mask by closing brain\n",fcname);
    if(!niik_image_morph_close_brain(outimg,close_brain_radius,close_brain_radius)) {
      fprintf(stderr,"ERROR: niik_image_morph_close_brain\n");
      return NULL;
    }
  }

  if(verbose>=3) {
    rgb[0]=rgb[1]=hilim*2.0;
    rgb[2]=0;
    niik_image_boundary(bounds,outimg,rgb,0,1,0);
    sprintf(fname,"tmp_bounds1.nii.gz");
    fprintf(stdout,"\twriting %s\n",fname);
    niik_image_write(fname,bounds);
    sprintf(fname,"tmp_bseg.nii.gz");
    fprintf(stdout,"\twriting %s\n",fname);
    niik_image_write(fname,outimg);
  }

  tmpimg=niik_image_free(tmpimg);
  bounds=niik_image_free(bounds);
  cm=niikvec_free(cm);

  fprintf(stdout,"[%s] finish\n",fcname);
  return outimg;
} /* niik_image_bseg_test1 */


int niik_image_bseg_test2_thresh(nifti_image *img,niikmat *regmat,double radius,double percentile,double *vout) {
  const char *NIIKDIR=NULL;
  char fname[4096],fcname[512]="niik_image_bseg_test2_thresh";
  nifti_image *mni_seg=NULL;
  int verbose=0;
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(regmat==NULL) {
    fprintf(stderr,"ERROR: regmat is null\n");
    return 0;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    return 0;
  }
  if(radius==0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.nii.gz",NIIKDIR);
  } else if(radius>0) {
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_dilate.nii.gz",NIIKDIR,(int)radius);
  } else if(radius<0) {
    radius=fabs(radius);
    sprintf(fname,"%s/data/bseg/MNI152_T1_1mm_brain_mask_edit.%imm_erode.nii.gz",NIIKDIR,(int)radius);
  }
  if(verbose>=1) fprintf(stdout,"[%s]    reading brain mask   %s\n",fcname,fname);
  if((mni_seg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s]    transform mni mask\n",fcname);
  if(!niik_image_inverse_affine_transform_3d_update(mni_seg,img,regmat,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  *vout = niik_image_get_percentile(img,mni_seg,percentile);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  mni_seg=niik_image_free(mni_seg);
  return 0;
}


nifti_image *niik_image_bseg_test2(nifti_image *img,niikmat *regmni,double thresh,double hithresh,double radius,double mni_radius) {
  nifti_image
  *mni_seg=NULL,
   *mni_fil=NULL,
    *maskimg=NULL;
  const char *NIIKDIR=NULL;
  char
  fname[512],
        fcname[64]="niik_image_bseg_test2";
  double d;
  kobj *obj=NULL;
  int
  thresh_method=1,
  i,
  verbose=1;

  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    return NULL;
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(regmni==NULL) {
    fprintf(stderr,"ERROR: regmat is null\n");
    return NULL;
  }

  /* THRESHOLD */
  if(thresh<0) {
    thresh=-thresh;
    if(!niik_image_bseg_basic_thresh(img,regmni,&thresh,thresh_method)) {
      fprintf(stderr,"ERROR: niik_image_bseg_basic_thresh\n");
      return NULL;
    }
  }
  if(verbose>=1) fprintf(stdout,"[%s] threshold %5.1f\n",fcname,thresh);
  if((maskimg=niik_image_threshold_new(img,thresh))==NULL) {
    fprintf(stderr,"ERROR: niik_image_threshold_new\n");
    return NULL;
  }

  /* GETTING DILATED BRAIN MASK */
  if(verbose>=1) fprintf(stdout,"[%s] standard brain masking %5.1f\n",fcname,mni_radius);
  if((mni_seg=niik_image_brain_mask_from_mni(img,regmni,mni_radius,0))==NULL) {
    fprintf(stderr,"ERROR: niik_image_brain_mask_from_mni\n");
    return NULL;
  }
  niik_image_mask(maskimg,mni_seg);
  mni_seg=niik_image_free(mni_seg);

  /* GETTING ERODED BRAIN MASK */
  if(verbose>=1) fprintf(stdout,"[%s] standard brain mask addition %5.1f\n",fcname,mni_radius);
  if((mni_seg=niik_image_brain_mask_from_mni(img,regmni,-mni_radius,0))==NULL) {
    fprintf(stderr,"ERROR: niik_image_brain_mask_from_mni\n");
    return NULL;
  }
  for(i=0; i<img->nvox; i++)
    if(niik_image_get_voxel(mni_seg,i)>0)
      niik_image_set_voxel(maskimg,i,1);

  /* UPPER THRESHOLD
     for(d=0.1;d<=0.899;d+=0.1){
     fprintf(stdout,"[%s]   %3.0f%%  %6.1f\n",fcname,d*100.0,niik_image_get_percentile(img,mni_seg,d)); }
     for(d=0.9;d<=1.01;d+=0.01){
     fprintf(stdout,"[%s]   %3.0f%%  %6.1f\n",fcname,d*100.0,niik_image_get_percentile(img,mni_seg,d)); }*/
  if(hithresh<0) {
    d=-hithresh;
    hithresh=niik_image_get_percentile(img,mni_seg,d);
    fprintf(stdout,"[%s]   %3.0f%%  %6.1f\n",fcname,d*100.0,hithresh);
  }

  /* REMOVE BRIGHT THINGS */
  if(verbose>=1) fprintf(stdout,"[%s] remove bright things around the brain edge %5.1f\n",fcname,hithresh);
  if(!niik_image_bseg_remove_eyes(img,maskimg,regmni,4.0,hithresh,1)) {
    fprintf(stderr,"ERROR: niik_image_bseg_remove_eyes\n");
    return NULL;
  }
  if(verbose>=2) {
    sprintf(fname,"tmp_bseg2_mask_rm.nii.gz");
    fprintf(stdout,"[%s]   writing output %s\n",fcname,fname);
    niik_image_write(fname,maskimg);
  }

  /* CONDITIONAL EROSION */
  sprintf(fname,"%s/data/CLADA/MNI152_T1_1mm_brain_edge_ROI_nonbrain_6mm_blur.nii.gz",NIIKDIR);
  if(verbose>=2) fprintf(stdout,"[%s] reading %s\n",fcname,fname);
  if((mni_fil=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  if(!niik_image_inverse_affine_transform_3d_update(mni_fil,img,regmni,NIIK_INTERP_LINEAR)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return NULL;
  }
  if(!niik_image_type_convert(mni_fil,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  for(i=0; i<img->nvox; i++) {
    niik_image_mul_voxel(mni_fil,i,radius);
  }
  if(!niik_image_morph_3d_radius_map(maskimg,mni_fil,NIIK_MORPH_ERODE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return NULL;
  }
  if(!niik_image_bseg_seed_fill_from_eroded_standard_mask(maskimg,regmni,-10)) {
    fprintf(stderr,"ERROR: niik_image_bseg_seed_fill_from_eroded_standard_mask\n");
    return NULL;
  }
  if(!niik_image_morph_3d_radius_map(maskimg,mni_fil,NIIK_MORPH_DILATE)) {
    fprintf(stderr,"ERROR: niik_image_morph_3d_radius_mask\n");
    return NULL;
  }

  /*if(verbose>=1) fprintf(stdout,"[%s] conditional opening %5.3f\n",fcname,radius);
    if(!niik_image_bseg_conditional_opening_imask(maskimg,mni_seg,radius)){
    fprintf(stderr,"ERROR: niik_image_bseg_conditional_opening_imask\n");
    return NULL; } */

  mni_seg=niik_image_free(mni_seg);

  if((obj=niik_image_bseg_get_brain_surface(regmni,0,0))==NULL) {
    fprintf(stderr,"ERROR: niik_image_bseg_get_brain_surface\n");
    return NULL;
  }

  if(1) {
    off_kobj_add_one_color(obj,0.7,0.4,0.25);
    sprintf(fname,"tmp_bseg2_obj0.off.gz");
    fprintf(stdout,"\twriting %s\n",fname);
    off_kobj_write_offply(fname,obj,0);
  }

  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return maskimg;
} /* niik_image_bseg_test2 */





#endif /* _FALCON_BSEG_C_ */
