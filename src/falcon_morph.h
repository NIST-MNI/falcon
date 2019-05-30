#pragma once
#ifndef __FALCON_MORPH_H__
#define __FALCON_MORPH_H__

/******************************************************************
 *
 * falcon_morph.c
 *
 ******************************************************************/
int niik_image_morph_3d_radius(nifti_image *img,int morphtype,double radius);
int niik_image_morph_3d_radius_mask(nifti_image *img,nifti_image *maskimg,int morph_type,double radius);
int niik_image_seed_fill(nifti_image *img,nifti_image *initimg,int flag_grad);
int niik_image_seed_fill_slow(nifti_image *img,nifti_image *initimg,int flag_grad);
int niik_image_seed_fill_xyz(nifti_image *img,int *idx,int flag_grad);
int niik_image_seed_fill_edge (nifti_image *img,int flag_grad);
int niik_image_seed_fill_from_middle (nifti_image *img,int flag_grad);
int niik_image_morph_close_brain(nifti_image *segimg,double dilate_radius,double erode_radius);
int niik_image_close_holes(nifti_image *img);

int niik_image_morph_3d_radius_map(nifti_image *img,nifti_image *Rmap,int morph_type);

int niik_image_morph_3d_mask(nifti_image *img,nifti_image *maskimg,int morph_type,int shape, int kdim);

int niik_image_morph_gray_dilate(nifti_image *img,double radius,double val);



#endif /*__FALCON_MORPH_H__*/

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/