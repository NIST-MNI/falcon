/* Filename:     nifti1_kunio_venseg.c
 * Description:  ventricle segmentation functions
 * Author:       Kunio Nakamura
 * Date:         March 13, 2012
 */


#ifndef _FALCON_VENSEG_C_
#define _FALCON_VENSEG_C_

#include "falcon.h"

int g_niik_venseg_debug = 0;

void niik_image_segven_set_debug(int debug) {
  g_niik_venseg_debug = debug;
}
int niik_image_segven_get_debug() {
  return g_niik_venseg_debug;
}

int niik_image_segment_ventricle(nifti_image *image,nifti_image *segimg,nifti_image *venimg,
                                 nifti_image *venroi_low_thresh,nifti_image *venroi_no_dilate,
                                 double thresh,int maxiter)
/***************************************
 * ventricle segmentation function
 * -image is the input t1w image
 * -segimg is the input brain mask
 * -venimg is the initial ventricle mask (probably obtained from
 *  nonlinear transform of ventricle mask)
 *   -the output is replaced here ***
 *   -get the most accurate parts so mainly the body or atria of the
 *    lateral ventricles and exclue the posterior horns (maybe even
 *    temporal horns as well
 * -venroi_low_thresh is the area of low threshold
 *   -this should contain postior horns and temporal horns as well as
 *    mid-sagittal section
 * -venroi_no_dilate is the areas of no dilation
 *   -mainly the midsagittal section
 * -thresh is the threshold between CSF and tissue
 *   -automatically estimated if thresh is very large
 * -maxiter is the number of dilation iteration
 *   -for large ventricles, one may need more iterations
 *   -potentail error (voxel outside the brain was in 'venimg') can
 *    result in really wild output
 ***************************************/
{
  nifti_image
  *tmpimg=NULL;
  double
  cthresh,
  ibrain,iCSF,
  imax;
  float *fimg;
  int
  verbose=0,
  iter,
  xdim,area,
  i,j,k,n;
  unsigned char
  *bven,*btmp,*broi,
  *bndil; /* venroi_no_dilate */


  /* checking inputs */
  if(image==NULL) {
    fprintf(stderr,"ERROR: image is null\n");
    return 0;
  }
  if(segimg==NULL) {
    fprintf(stderr,"ERROR: segimg is null\n");
    return 0;
  }
  if(venimg==NULL) {
    fprintf(stderr,"ERROR: venimg is null\n");
    return 0;
  }
  if(venroi_low_thresh==NULL) {
    fprintf(stderr,"ERROR: venroi_low_thresh is null\n");
    return 0;
  }
  if(venroi_no_dilate ==NULL) {
    fprintf(stderr,"ERROR: venroi_no_dilate is null\n");
    return 0;
  }

  if(niik_image_cmp_dim(image,segimg)!=0) {
    fprintf(stderr,"ERROR: wrong image dimension: %3i %3i %3i (segimg)\n",segimg->nx,segimg->ny,segimg->nz);
    return 0;
  }
  if(niik_image_cmp_dim(image,venimg)!=0) {
    fprintf(stderr,"ERROR: wrong image dimension: %3i %3i %3i (venimg)\n",venimg->nx,venimg->ny,venimg->nz);
    return 0;
  }
  if(niik_image_cmp_dim(image,venroi_low_thresh)!=0) {
    fprintf(stderr,"ERROR: wrong image dimension: %3i %3i %3i (venroi low thresh)\n",
            venroi_low_thresh->nx,venroi_low_thresh->ny,venroi_low_thresh->nz);
    return 0;
  }
  if(niik_image_cmp_dim(image,venroi_no_dilate)!=0) {
    fprintf(stderr,"ERROR: wrong image dimension: %3i %3i %3i (venroi no dilate)\n",
            venroi_no_dilate->nx,venroi_no_dilate->ny,venroi_no_dilate->nz);
    return 0;
  }

  /* convert data types */
  if(!niik_image_type_convert(segimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert %s \n",segimg->fname);
    return 0;
  }
  if(!niik_image_type_convert(venimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert %s \n",venimg->fname);
    return 0;
  }
  if(!niik_image_type_convert(venroi_low_thresh,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert %s \n",venroi_low_thresh->fname);
    return 0;
  }
  if(!niik_image_type_convert(venroi_no_dilate,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert %s \n",venroi_no_dilate->fname);
    return 0;
  }

  if(verbose) {
    fprintf(stdout,"[niik_image_segment_ventricle] start\n");
    fprintf(stdout,"[niik_image_segment_ventricle] iter   %i\n",maxiter);
  }

  if((fimg = niik_image_get_voxels_as_float_vector(image))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_float_vector(image)\n");
    return 0;
  }

  /* calculate the threshold */
  if(niik_check_double_problem(thresh)) {
    imax   = niik_image_get_max(image,NULL);
    if(verbose) fprintf(stdout,"   [niik_image_segment_ventricle] max intensity %7.3f\n",imax);
    ibrain = niik_image_get_mode(image,segimg,0,imax,128,3);
    for(i=0; i<image->nvox; i++) {
      if(fimg[i]>ibrain) niik_image_set_voxel(venimg,i,0);
    }
    iCSF   = niik_image_get_mode(image,venimg,0,imax,128,3);
    thresh = (ibrain + iCSF)/2.0;
    if(verbose>=0) {
      fprintf(stdout,"[niik_image_segment_ventricle] brain intensity %7.3f\n",ibrain);
      fprintf(stdout,"[niik_image_segment_ventricle] CSF intensity   %7.3f\n",iCSF);
    }
  } /* calculate thresh if undefined */
  if(verbose>=0) {
    fprintf(stdout,"[niik_image_segment_ventricle] threshold       %7.3f\n",thresh);
  }

  /* modify the input ventricle mask */
  bven = venimg->data;
  broi  = venroi_low_thresh->data;
  bndil = venroi_no_dilate->data;
  xdim = image->nx;
  area = image->ny * xdim;
  if(verbose) fprintf(stdout,"[niik_image_segment_ventricle] threshold image to improve original ventricle mask\n");
  for(i=0; i<image->nvox; i++) {
    if(fimg[i]>thresh) bven[i]=0;
  }
  if(verbose>1) {
    fprintf(stdout,"[niik_image_segment_ventricle] write tmp_segven01_predilation.nii.gz\n");
    niik_image_write("tmp_segven01_predilation.nii.gz",venimg);
  }

  /* create a temp space */
  if((tmpimg = niik_image_copy(venimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  btmp = tmpimg->data;

  for(iter=0; iter<maxiter; iter++) {
    /* conditional dilation (6 neighbors)
     * -check for intensity
     * -ROI-dependent threshold
     */
    if((iter%5)==0) {
      fprintf(stdout,"\t  iter %3i %9.4f\n",iter+1,
              niik_image_count_mask(venimg) *
              niik_image_get_voxel_size(venimg) / 1000.0);
    }

    for(k=n=0; k<image->nz; k++) {
      for(j=0; j<image->ny; j++) {
        for(i=0; i<image->nx; n++,i++) {
          if(!bven[n]) continue;
          cthresh=thresh;
          if(broi[n]) cthresh=thresh*0.5;

          if(i>0) {
            if(!bven[n-1]) {
              if(!bndil[n-1]) {
                if(fimg[n-1]<cthresh) {
                  btmp[n-1]=1;
                }
              }
            }
          }
          if(j>0) {
            if(!bndil[n-xdim]) {
              if(!bven[n-xdim]) {
                if(fimg[n-xdim]<cthresh) {
                  btmp[n-xdim]=1;
                }
              }
            }
          }
          if(k>0) {
            if(!bndil[n-area]) {
              if(!bven[n-area]) {
                if(fimg[n-area]<cthresh) {
                  btmp[n-area]=1;
                }
              }
            }
          }

          if(i<image->nx-1) {
            if(!bndil[n+1]) {
              if(!bven[n+1]) {
                if(fimg[n+1]<cthresh) {
                  btmp[n+1]=1;
                }
              }
            }
          }
          if(j<image->ny-1) {
            if(!bndil[n+xdim]) {
              if(!bven[n+xdim]) {
                if(fimg[n+xdim]<cthresh) {
                  btmp[n+xdim]=1;
                }
              }
            }
          }
          if(k<image->nz-1) {
            if(!bndil[n+area]) {
              if(!bven[n+area]) {
                if(fimg[n+area]<cthresh) {
                  btmp[n+area]=1;
                }
              }
            }
          }
        }
      }
    }

    for(i=0; i<image->nvox; i++) bven[i]=btmp[i];
  }
  fprintf(stdout,"\t  iter %3i %9.4f\n",iter+1,
          niik_image_count_mask(venimg) *
          niik_image_get_voxel_size(venimg) / 1000.0);

  if(verbose>1) {
    fprintf(stdout,"[niik_image_segment_ventricle] write tmp_segven02_postdilation.gz\n");
    niik_image_write("tmp_segven02_postdilation.nii.gz",venimg);
  }

  tmpimg=niik_image_free(tmpimg);
  return 1;
}


nifti_image *niik_image_segment_ventricle_prep_from_mni_warp(nifti_image *refimg,nifti_image *warpimg,int datatype) {
  nifti_image *outimg;
  const char *NIIKDIR=NULL;
  char  fname[4096];
  int verbose=0;

  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return NULL;
  }
  if((NIIKDIR=get_NIIKDIR())==NULL) {
    fprintf(stderr,"ERROR: please setenv NIIKDIR\n");
    return NULL;
  }

  switch(datatype) {
  case 1:
    sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_ven_mask.nii.gz",NIIKDIR);
    break;
  case 2:
    sprintf(fname,"%s/data/CLADA/MNI152_T1_2mm_ven_lothresh_roi.nii.gz",NIIKDIR);
    break;
  case 3:
    sprintf(fname,"%s/data/CLADA/MNI152_T1_1mm_ven_fil_no_dilation_area.nii.gz",NIIKDIR);
    break;
  default:
    fprintf(stderr,"ERROR: unknown datatype, %i\n",datatype);
    return NULL;
  }

  if(verbose) fprintf(stdout,"[niik_image_segment_ventricle_prep_from_mni_warp] reading image %s\n",fname);
  if((outimg=nifti_image_read(fname,1))==NULL) {
    fprintf(stderr,"ERROR: nifti_image_read %s\n",fname);
    return NULL;
  }
  if(warpimg==NULL) return outimg;

  if(verbose) fprintf(stdout,"[niik_image_segment_ventricle_prep_from_mni_warp] warp image\n");
  if(!niik_image_apply_3d_warp_update(outimg,refimg,warpimg,NIIK_WARP_MAP_LOC,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_3d_warp_update(outimg,refimg,warpimg,NIIK_WARP_MAP_LOC,NIIK_INTERP_NN)\n");
    return NULL;
  }

  if(verbose) fprintf(stdout,"[niik_image_segment_ventricle_prep_from_mni_warp] finish\n");
  return outimg;
}


#endif /* _FALCON_VENSEG_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/