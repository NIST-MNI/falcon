/*
# Developed by Simon Fristed Eskildsen, eskild@gmail.com
# Part of FALCON
#
# Copyright notice:
# This code is copyright Simon Fristed Eskildsen.
# It may not be copied, altered in any way or transmitted
# to others (unless explicitly stated otherwise) without
# the written permission of the author/developer. 
*/

#ifndef VOLUME_TOOLS_H
#define VOLUME_TOOLS_H

#include <volume_io.h>
#include "surface_tools_basic.h"
#include "volume_io_wrap.h"

#define RCC 0
#define BCC 1

#define UNVISITED 0
#define VISITED 1
#define DELETED 2
#define LEAFNODE 3

#undef BLACK
#undef WHITE
#undef RED
#define BLACK 0
#define WHITE 1
#define RED 2


struct k_node{
  struct k_node *representative;
  short component;
  byte colour;
  int slice;
  byte type;
  int num_edges;
  struct k_edge **edges;
  int size;
};


struct k_edge{
  struct k_node *a,*b;
  int weight;
  byte colour;
};

struct k_graph{
  struct k_node *V;
  struct k_edge *E;
  int numV;
  int numE;
};


point3D volume_cog( VIO_Volume *original);
point3D volume_cog_iso( VIO_Volume *original, VIO_Real isovalue);
point3D volume_cog_thresh( VIO_Volume *original, float thresh);
point3D volume_cog_threshold( VIO_Volume *original, VIO_Real threshold);
point3D volume_cog_segment( VIO_Volume *original, VIO_Real lower, VIO_Real upper);
int rotate_volume(VIO_Volume *original, VIO_Volume *rotated, int direction);
int rotate_volume_angle(VIO_Volume *original, VIO_Volume *rotated, float xrot, float yrot, float zrot);
int create_binary_volume(VIO_Volume *vol, VIO_Real *starts, int *sizes);
//int expand_volume(VIO_Volume *vol, int *expansion);
int expand_volume(VIO_Volume *vol, VIO_Volume *expand, point3D *exp);
int decrease_volume(VIO_Volume *vol, VIO_Volume *decreased, point3D *dec);
float volume_density(VIO_Volume *vol, VIO_Real isovalue, VIO_Real threshold, int mode);
int add_volumes(VIO_Volume *vol1, VIO_Volume *vol2, VIO_Volume *result);
VIO_Volume or_volumes(VIO_Volume vol1, VIO_Volume vol2);
VIO_Volume not_volume(VIO_Volume vol);
int gradient_scale_volume(VIO_Volume *gradient, VIO_Volume *wm_membership, VIO_Volume *gm_membership, VIO_Volume *csf_membership, VIO_Volume *result);
int normalizeVolume_max_min(VIO_Volume *vol, VIO_Volume *normalized);
int normalize_volume_max_min(VIO_Volume *vol);
int set_intensity_by_mask(VIO_Volume original, VIO_Volume mask, VIO_Real intensity);
//int threshold_volume(VIO_Volume *orig, VIO_Volume *out, VIO_Real min, VIO_Real max);
VIO_Volume threshold_volume(VIO_Volume orig, VIO_Real min, VIO_Real max);
VIO_Status minc_blur3D_volume(VIO_Volume data, double fwhmx, double fwhmy, double fwhmz, int ndim, int kernel_type);
int binary_boundings(VIO_Volume binary, int *lower, int *upper, VIO_Real iso);
int get_bounding_box(VIO_Volume *vol, point3D *c1, point3D *c2, VIO_Real thresh);
float brain_size( VIO_Volume *original);
int brain_size_3D(VIO_Volume *original, float *x, float *y, float *z);
int createEllipsoid(VIO_Volume volume, float a, float b, float c, float r, float x0, float y0, float z0);
int erode3D(VIO_Volume *input, VIO_Volume *output, int iso, int times, int coeff, VIO_Volume *mask);
int sub_volumes(VIO_Volume *vol1,VIO_Volume *vol2,VIO_Volume *out, int iso1, int iso2);
VIO_Volume sub_volumes_real(VIO_Volume vol1,VIO_Volume vol2);
//int defuzzyfi(VIO_Volume *wm, VIO_Volume *gm, VIO_Volume *output, float level, int crisp_value);
int defuzzyfi(VIO_Volume *wm, VIO_Volume *gm, VIO_Volume *output, VIO_Real level, VIO_Real label);
void dilate_voxel(Volume_wrap *wdilated, int x, int y, int z, byte fill);
int dilate3D(VIO_Volume *input, VIO_Volume *output, VIO_Real iso, int times, VIO_Volume *mask);
int dilate3D_new(VIO_Volume *input, VIO_Volume *output, VIO_Real iso, int times, int six_connect, VIO_Volume *mask);
int dilate_measurements(VIO_Volume *input, VIO_Volume *output, VIO_Real background, int times, int six_connect, VIO_Volume *mask);
//int group(VIO_Volume *volume, VIO_Real label, VIO_Real fill, VIO_Real bg_value, VIO_Volume *output);
VIO_Volume group(VIO_Volume volume, VIO_Real label, VIO_Real fill, VIO_Real bg_value);
int mask_voxels(VIO_Volume *mask,VIO_Volume *vol, int iso);
float mean_intensity( VIO_Volume *original );
float mean_intensity2(VIO_Volume *original, VIO_Volume *segmented, short label);
//int median_filter_curvature(VIO_Volume *volume,VIO_Volume *output, int times, int median_val, int ignore_bg, int six_connect, int orig_six_connect, VIO_Volume *curvature, int use_curvature, float curv_val);
int median_filter_curvature(VIO_Volume *volume,VIO_Volume *output, int times, int median_val, int ignore_bg, int six_connect, int orig_six_connect, VIO_Volume *curvature, int use_curvature, float curv_val, float bg_threshold, VIO_Volume *mask, int use_mask);
int median_filter_new(VIO_Volume *volume,VIO_Volume *output, int times, int median_val, int ignore_bg, int six_connect, int orig_six_connect);
int removeGray(VIO_Volume *gm, VIO_Volume *wm, int iso);
int maskCuttingPlane(VIO_Volume *volume,VIO_Volume  *output, int value);
inline float interpolate_volume_nn(VIO_Volume *volume, point3D *currentPoint);

/* used by topology correction */
void zero_volume(Volume_wrap *wvolume);
byte isconnected(Volume_wrap *wbinary,int x,int y,int z, byte iso);
int topological_number(point3D point, Volume_wrap *wbody, short iso);
int topological_number_short(point3D point, Volume_wrap *wbody, short iso);
void SCC_labeling(Volume_wrap *wbody, Volume_wrap *wresidue, short body_label, short residue_label);
byte isconnected_short(Volume_wrap *wbinary,int x,int y,int z, short iso);
byte are_strongly_connected(Volume_wrap *wbody,point3D *point1,point3D *point2, short iso);
byte are_weakly_connected(Volume_wrap *wresidue, Volume_wrap *wbody, point3D *point1, point3D *point2, short body_iso, short residue_iso);
int gradient3D(VIO_Volume *volume, VIO_Volume *length);

float interpolate_volume( VIO_Volume *volume, point3D *currentPoint );
float interpolate_volume_cubic(VIO_Volume *volume, point3D *currentPoint, BOOLEAN world);
int gpf(VIO_Volume *volume, float stddev, VIO_Volume *dx, VIO_Volume *dy, VIO_Volume *dz );

int split_volume(VIO_Volume *vol, VIO_Real iso, point3D *result);
int optimize_slice_position(VIO_Volume *vol, VIO_Real iso, int dimension, int slice, int searchspace, int *no_of_voxels);
int mask_by_slices(VIO_Volume *vol, int dimension, int slicestart, int sliceend);

VIO_Real min_intensity( VIO_Volume original );
VIO_Real max_intensity( VIO_Volume original );

int volume_nearest_label(VIO_Volume mask, VIO_Volume atlas, VIO_Volume out, int margin);

/* this is not implemented yet */
int rodell_eskildsen_measure(VIO_Volume *inner_distmap, VIO_Volume *outer_distmap, VIO_Volume *inner_chamfer, VIO_Volume *outer_chamfer, VIO_Volume *output, VIO_Volume *mask);

#endif
 

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
*/