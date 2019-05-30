#pragma once
#ifndef __FALCON_SEG_H__
#define __FALCON_SEG_H__

#include "falcon_surfaces.h"

/******************************************************************
 *
 * falcon_bseg.c
 *
 ******************************************************************/

/* helpful functions (falcon_bseg.c) */
int niik_image_bseg_basic_thresh(nifti_image *img,niikmat *regmni,double *thresh,int method);
int niik_image_bseg_basic_thresh_with_mask(nifti_image *img,nifti_image *mni_seg,double *thresh,int method);
kobj *niik_image_bseg_get_brain_surface(niikmat *regmni,double imove,int verbose);
nifti_image *niik_image_bseg_from_mni(nifti_image *img,niikmat *regmni,double thresh,double radius);
int niik_image_bseg_remove_eyes(nifti_image *img,nifti_image *maskimg,niikmat *regmni,double radius,double hithresh,int along_brain_edge);
int niik_image_bseg_seed_fill_from_eroded_standard_mask(nifti_image *maskimg,niikmat *regmat,int erode_radius);
niikmat *niik_image_bseg_get_mni_affine_matrix_from_sform(nifti_image *img);


/* actual segmentation functions (falcon_bseg.c) */
nifti_image *niik_image_bseg_simplest(nifti_image *img,niikmat *regmni,double thresh,int mni_dilate,double radius);
nifti_image *niik_image_bseg_basic(nifti_image *img,int dilate_iter,double lothresh,double mdthresh,double hithresh,double hithresh2,double radius,niikmat *regmni,int flag_pv);

nifti_image *niik_image_bseg_test1(nifti_image *img,niikmat *regmni,kobj *obj,double thresh,double lolim,double hilim,double radius,double iradius,double mni_dilate_radius,int maxiter,int using_surface,double close_brain_radius,double radius_for_lim_erosion,double curvature_thresh,double curvature_range);

nifti_image *niik_image_bseg_test2(nifti_image *img,niikmat *regmni,double thresh,double hithresh,double radius,double mni_radius); /* this is really in testing ... */
nifti_image *niik_image_bseg_test3(nifti_image *img,niikmat *regmni,double thresh,double hithresh,double radius,double mni_radius);


/******************************************************************
 *
 * falcon_venseg.c
 *
 ******************************************************************/

void niik_image_segven_set_debug(int debug);
int  niik_image_segven_get_debug();

int niik_image_segment_ventricle(nifti_image *image,nifti_image *segimg,nifti_image *venimg,
                                 nifti_image *venroi_low_thresh,nifti_image *venroi_no_dilate,
                                 double thresh,int maxiter);
nifti_image *niik_image_segment_ventricle_prep_from_mni_warp(nifti_image *refimg,nifti_image *warpimg,int datatype);


#endif /*__FALCON_SEG_H__*/



/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/