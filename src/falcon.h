/* Filename:     nifti1_kunio.h
 * Description:  general nifti1 function header
 * Author:       Kunio Nakamura
 * Date:         February 23, 2012
 *
 */
#pragma once

#ifndef __FALCON_H__
#define __FALCON_H__

#include "falcon_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif /*HAVE_MALLOC_H*/

#include <string.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>

/* nifti1 headers */
#include "nifti1.h"
#include "nifti1_io.h"

/* TIFF  library */
#include "tiffio.h"

/* ZLIB */
#include <zlib.h>

#define NIIK_MINC_SUPPORT TRUE

/* EDGE FILES ARE COMPRESSED */
#define COMPRESS_OFF_EDGE_FILE 1

/* CONSTANTS */
#define NIIK_VAL_UCMAX 255u
#define NIIK_VAL_USMAX 65535u
#define NIIK_VAL_UIMAX 16777215u
#define NIIK_VAL_ULMAX 4294967295u
#define NIIK_VAL_CMAX 127
#define NIIK_VAL_SMAX 32767
#define NIIK_VAL_IMAX 8388607
#define NIIK_VAL_LMAX 2147483647

#define NIIKMAX 1e301
#define NIIKPI  3.141592653589793238462643383280

/* CONSTANTS */
#define NIIK_PI       3.141592653589793238462643383280
#define NIIK_PIhalf   1.57079632679489661923132169164
#define NIIK_PI2      6.283185307179586476925286766559
/*#define NIIK_DEGREE2RAD  0.017453292519943295769236907685
  #define NIIK_RAD2DEGREE 57.295779513082320876798154814096*/

/* INTERPOLATION TYPE DEFINITIONS */
#define NIIK_INTERP_NN 0   /* nearest neighbor */
#define NIIK_INTERP_LINEAR  1  /* trilinear interpolation */
#define NIIK_INTERP_BSPLINE 2  /* cubic b-spline */
#define NIIK_INTERP_SPLINE  3  /* cubic spline */
#define NIIK_INTERP_SINC    4  /* sinc */

/* REGISTRATION TYPE DEFINITIONS */
#define NIIK_REGISTER_UNKNOWN  0   /* unknown */
#define NIIK_REGISTER_SAD      1   /* sum of absolute difference */
#define NIIK_REGISTER_SSQ      2   /* sum of squared difference */
#define NIIK_REGISTER_CC       3   /* correlation coefficietn */
#define NIIK_REGISTER_NMI      4   /* normalized mutual information */

/* MORPHOLOGIC TYPE DEFINITIONS */
#define NIIK_MORPH_DILATE 1   /* dilation */
#define NIIK_MORPH_ERODE  2   /* erosion */
#define NIIK_MORPH_CLOSE  3   /* closing */
#define NIIK_MORPH_OPEN   4   /* opening */
#define NIIK_MORPH_2D_KERNEL_DIAMOND   101   /* kernel is 2d diamond */
#define NIIK_MORPH_2D_KERNEL_SQUARE    102   /* kernel is 2d square */
#define NIIK_MORPH_3D_KERNEL_DIAMOND   103   /* kernel is 3d diamond */
#define NIIK_MORPH_3D_KERNEL_SQUARE    104   /* kernel is 3d square */

/* NUMERICAL METHOD DEFINITIONS */
#define NIIK_NELDER_MEAD_COST_UNKNOWN   0  /* unknown */
#define NIIK_NELDER_MEAD_COST_RATIO     1  /* ratio */
#define NIIK_NELDER_MEAD_COST_ABS       2  /* absolute */
#define NIIK_NELDER_MEAD_COST_DIFF      3  /* difference */

/* WARP MAP TYPE */
#define NIIK_WARP_MAP_UNKNOWN      0  /* unknown */
#define NIIK_WARP_MAP_DISP         1  /* displacement map */
#define NIIK_WARP_MAP_LOC          2  /* location map (coordinates) */
#define NIIK_WARP_MAP_DISP_BSPLINE 3  /* displacement map by b-spline coefficients */

/* COLOR MAP TYPE */
#define NIIK_COLORMAP_UNKNOWN      0  /* unknown */
#define NIIK_COLORMAP_ATROPHY      1  /* atrophy color map */
#define NIIK_COLORMAP_SUMMER       2  /* summer color map */
#define NIIK_COLORMAP_SPECTRAL     3  /* spectral color map */
#define NIIK_COLORMAP_JACOBIAN     4  /* atrophy color map 2 */
#define NIIK_COLORMAP_GRAYSCALE    5  /* grayscale map */

#define NIIK_COLORMAP_NEG_ATROPHY  6  /* negative part of atrophy lut (blues)*/
#define NIIK_COLORMAP_POS_ATROPHY  7  /* positive part of atrophy lut (reds)*/


/* CURVE FITTING METHOD DEFINITIONS */
#define NIIK_CURVE_FIT_UNKNOWN   0  /* unknown */
#define NIIK_CURVE_FIT_ABSOLUTE  1  /* absolute */
#define NIIK_CURVE_FIT_SQUARE    2  /* square */

/* CLASSIFICATION METHOD DEFINITIONS */
#define NIIK_CLASSIFY_UNKNOWN   0  /* unknown */
#define NIIK_CLASSIFY_KMEANS    1  /* k-means clustering */
#define NIIK_CLASSIFY_FCM       2  /* fuzzy c-means clustering */

/*Convenient errror-handling macros*/
#define NIIK_RET0(test,fcname,msg)      if((test)>0) { fprintf(stdout,"[%s:%i:%s] ERROR: %s\n",__FILE__,__LINE__,__func__,msg); return 0; }
#define NIIK_RET1(test,fcname,msg)      if((test)>0) { fprintf(stdout,"[%s:%i:%s] ERROR: %s\n",__FILE__,__LINE__,__func__,msg); return 1; }
#define NIIK_RET0n(test,fcname,msg)     if((test)>0) { fprintf(stdout,"[%s:%i:%s] ERROR: %s\n",__FILE__,__LINE__,__func__,msg); return NULL; }
#define NIIK_EXIT(test,fcname,msg,code) if((test)>0) { fprintf(stdout,"[%s:%i:%s] ERROR: %s\n",__FILE__,__LINE__,__func__,msg); exit(code); }
#define NIIK_RETURN(test,msg,code)      if((test)>0) { fprintf(stdout,"[%s:%i:%s] ERROR: %s\n",__FILE__,__LINE__,__func__,msg); return(code); }

typedef struct {                   /** 4x4 matrix struct **/
  double m[4][4] ;
} mat44d ;

/* 3d point with additional number */
typedef struct {
  double x,y,z,w;
} niikpt;  /* most functions are in falcon_point.c */

/* Spherical coordinates */
typedef struct {
  double the,psi;
} niiksph;  /* most functions are in falcon_point.c */


typedef struct {
  double **m;
  int row,col; /* row-col matrix */
} niikmat;  /* general matrix */

typedef struct {
  double *v;
  int num; /* num-element vector */
} niikvec;  /* general vector */


typedef struct {
  niikvec **col;
  char **name; /*column name*/
  int ncol;    /* number of columns */
} niiktable;   /* table */



typedef struct {                    /* object for normalized mutual information */
  double nmi;                       /* normalized mutual information */
  double hmin[2],hmax[2],hran[2];   /* histogram min/max/range for image 1/2 */
  double **h1,**h2;                 /* histograms for image 1/2 and 2d histo: h1[0,1][variable=(hnum[0] and hnum[1]?] */
  int hnum[2];                      /* size of histogram arrays */
} nmi_obj;                     /* functions are in falcon_register.c */



#include "falcon_text_io.h"



/******************************************************************
 *
 * inpainting
 *
 ******************************************************************/


/******************************************************************
 *
 * falcon_rand.c
 *
 ******************************************************************/
double niik_get_rand();
double niik_get_rand_normal();


/******************************************************************
 *
 * nifti1_kunio.c
 *
 ******************************************************************/

int niik_check_double_problem(double val);

int niik_access_file(char *f);

/* basic functions (nifti1_kunio.c) */
#include "falcon_inline.h"

double NIIK_Heaviside(double val,double eps);
double NIIK_Heaviside11(double val,double eps);
double NIIK_DiracDelta(double val,double eps);
double NIIK_BoxMuller(double v1,double v2);
double NIIK_RiceRnd(double v,double s);
double NIIK_RiceRnd2(double v,double s);

double niik_linear_step(double val,double eps);
double niik_pvc(double val,double m1, double m2);
double niik_pv(double val,double m1, double m2);
int niik_triangular_number(int num);
double niik_legendre_func(double x,int deg);
int niik_legendre_func_with_param(double x,double *y,double *param,int deg);

void niik_disp_exec_info(char *progname,int major,int minor,int micro);

char ** niik_csstring_to_list(char *CSstr,int *num);
double *niik_csstring_to_double_list(char *CSstr,int *num);
int *   niik_csstring_to_int_list(char *CSstr,int *num);
int     niik_csstring_count(char *CSstr);
char ** niik_txtfile_to_list(char *filename,int *num) ;

void    niik_free_list(char **list,int num);

int niik_check_fsldir_exists();


/* colormap functions (nifti1_kunio.c) */
niikmat *niik_colormap_get(int ctype,int num);
niikmat *niik_colormap_get_atrophy(int num);
niikmat *niik_colormap_get_pos_atrophy(int num);
niikmat *niik_colormap_get_neg_atrophy(int num);
niikmat *niik_colormap_get_summer(int num);
niikmat *niik_colormap_get_spectral(int num);
niikmat *niik_colormap_get_jacobian(int num);
niikmat *niik_colormap_get_grayscale(int num);
int niik_image_convert_to_color_image(nifti_image *img);
nifti_image *niik_image_map_color(nifti_image *img,double imin,double imax,niikmat *cmap);
nifti_image *niik_image_map_color_overlay_mask(nifti_image *img,double imin,double imax,
    nifti_image *map,double map_min,double map_max,
    double bg_factor,
    nifti_image *maskimg,niikmat *cmap);

/* basic (nifti1_kunio.c) */
int niik_version_display(const char *fcname,int major,int minor,int micro);
int niik_fc_display(const char *fcname,int flag_start);

int niik_underscore_to_space(char *str);

char *niik_nifti1_xform_string( int code );
char *niik_numerical_method_string( int code );
char *niik_image_dim_string( int code );
int niik_image_cmp_dim(nifti_image *img1,nifti_image *img2);
int niik_image_count_mask(nifti_image *maskimg);
int niik_image_clear(nifti_image *img);
int niik_image_flip_bgfg(nifti_image *maskimg);
int niik_image_mask(nifti_image *img,nifti_image *maskimg);
int niik_image_maskout(nifti_image *img,nifti_image *maskimg);
int niik_image_add_masks(nifti_image *modified,nifti_image *unchange);
nifti_image *niik_image_read(const char *fname);
nifti_image **niik_image_read_multiple(char *fname,int *num);
nifti_image **niik_image_read_multiple_from_file(char *fname,int *num);
int niik_image_write(const char *fname,nifti_image *img);
nifti_image *niik_image_free(nifti_image *img);
double niik_image_get_voxel_size(nifti_image *img);
double niik_image_get_mask_vol(nifti_image *maskimg);
int niik_image_get_index(nifti_image *img,int x,int y,int z);
int niik_image_get_index_niikpt(nifti_image *img,niikpt p);
int niik_image_check_index_niikpt(nifti_image *img,niikpt p);
niikpt niik_image_get_pt_from_index(nifti_image *img,int idx,int *x,int *y,int *z);
int niik_image_display_stats(nifti_image *img,nifti_image *maskimg);
int niik_image_label_display_stats(nifti_image *img,nifti_image *labelimg,int label);
int niik_image_dim_display(nifti_image *img,char *s);
int niik_image_get_mask_crop_roi(nifti_image *maskimg,int *xmin,int *xmax,int *ymin,int *ymax,int *zmin,int *zmax);
int niik_image_crop(nifti_image *img,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax);
int niik_image_pad_3d(nifti_image *img,int xdim,int ydim,int zdim);

int niik_image_iscale(nifti_image *img,double imin,double imax,double omin,double omax);
int niik_image_append_history(nifti_image *img,const char *history_entry);
char* niik_create_minc_timestamp(int argc,char *argv[]);



/* type converters (nifti1_kunio.c) */
int niik_image_type_convert_scl(nifti_image * img,int datatype,int scl_flag); /* type convert with/without intensity scaling */
int niik_image_type_convert(nifti_image * img,int datatype); /* generic type converter */
int niik_image_auto_type_convert(nifti_image * img);

/* copy image (nifti1_kunio.c) */
nifti_image *niik_image_copy(nifti_image * src); /* copies image */
nifti_image *niik_image_copy_as_type(nifti_image *img,int datatype);
nifti_image *niik_image_init(int x,int y,int z,int t,int u,int v,int w,double dx,double dy,double dz,double dt,double du,double dv,double dw,int datatype);
nifti_image **niik_image_array(int num);

/* voxel info (nifti1_kunio.c) */
static inline double  niik_image_get_voxel(nifti_image *img,int p); /* get voxel value */
static inline double  niik_image_voxel_get(nifti_image *img,int p); /* get voxel value */
int niik_image_get_voxels_as_double_vector_update(nifti_image *img,double *v);

double *niik_image_get_voxels_as_double_vector(nifti_image *img);  /* get the voxels in double vector */
float  *niik_image_get_voxels_as_float_vector(nifti_image *img);
unsigned char *niik_image_get_voxels_as_uint8_vector(nifti_image *img);
long *niik_image_get_voxels_as_long_vector(nifti_image *img);
niikvec *niik_image_get_voxels_as_double_vector_within_mask(nifti_image *img,nifti_image *maskimg);
int niik_image_set_voxels_from_double_vector(nifti_image *img,double *v); /* set the voxels from double vector */
int niik_image_set_voxels_from_uint8_vector(nifti_image *img,unsigned char *v); /* set the voxels from uint8 vector */
int niik_image_set_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val);
static inline int niik_image_set_voxel(nifti_image *img,int p,double d); /* set voxel value */
static inline int niik_image_add_voxel(nifti_image *img,int p,double d);
static inline int niik_image_mul_voxel(nifti_image *img,int p,double val);
int niik_image_add_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val);
int niik_image_mul_voxels_ROI(nifti_image *img,int xmin,int ymin,int zmin,int xmax,int ymax,int zmax,double val);

int niik_image_copy_data(nifti_image *src,nifti_image *dest);
int niik_image_copy_ref_info(nifti_image *src,nifti_image *dest);

#include "falcon_image_inline.h"

/* filter (nifti1_kunio.c) */
int niik_image_filter_gaussian_update(nifti_image *img,int kdim,double FWHM);
nifti_image *niik_image_filter_gaussian(nifti_image *img,int kdim,double FWHM);
nifti_image *niik_image_3d_average_filter(nifti_image *img,int size);

/* stats (nifti1_kunio.c) */
double niik_image_get_mean(nifti_image *img,nifti_image *maskimg);
double niik_image_get_stdv(nifti_image *img,nifti_image *maskimg);
double niik_image_get_var(nifti_image *img,nifti_image *maskimg);
double niik_image_get_skew(nifti_image *img,nifti_image *maskimg);
double niik_image_get_min(nifti_image *img,nifti_image *maskimg);
double niik_image_get_max(nifti_image *img,nifti_image *maskimg);
double niik_image_get_sum(nifti_image *img,nifti_image *maskimg);
double niik_image_get_percentile(nifti_image *img,nifti_image *maskimg,double percent);
double niik_image_get_mode(nifti_image *img,nifti_image *maskimg,double dmin,double dmax,int num,int avgnum);
double niik_image_get_median(nifti_image *img,nifti_image *maskimg);

nifti_image *niik_image_label_to_mask(nifti_image *labelimg,double label);
double niik_image_get_mode_label(nifti_image *img,nifti_image *labelimg,int label,double dmin,double dmax,int num,int avgnum);
double niik_image_get_mean_label(nifti_image *img,nifti_image *labelimg,int label);
double niik_image_get_stdv_label(nifti_image *img,nifti_image *labelimg,int label);

int niik_image_get_max_index(nifti_image *img,nifti_image *maskimg);
int niik_image_get_min_index(nifti_image *img,nifti_image *maskimg);


/* multiple dimensions (nifti1_kunio.c) */
nifti_image **niik_image_unmerge3d(nifti_image *img);
nifti_image *niik_image_unmerge3d_extract_one(nifti_image *img,int n);
nifti_image *niik_image_combine(nifti_image **imglist,int nimg,int dimnum,double scale_max);
nifti_image *niik_image_montage(nifti_image **imglist,int nimg,int tile,int xlo,int ylo,int zlo,int xhi,int yhi,int zhi);
int niik_image_combine_and_write_as_type(char *filename,nifti_image **imglist,int nimg,char dim,double scale_max,int datatype);
int niik_image_combine_and_write(char *filename,nifti_image **imglist,int nimg,char dim,double scale_max);

/* multiple images (nifti1_kunio.c) */
nifti_image *niik_image_average_multiple(nifti_image **imglist,int nimg,double *faclist);
nifti_image *niik_image_maximize_multiple(nifti_image **imglist,int nimg,double *faclist);
nifti_image *niik_image_add_multiple(nifti_image **imglist,int nimg,double *faclist);
nifti_image *niik_image_mul_multiple(nifti_image **imglist,int nimg);
int niik_image_multiply_2_images(nifti_image *changed,nifti_image *mult);
int niik_image_add_2_images(nifti_image *changed,nifti_image *add);


/* tiff (nifti1_kunio.c) */
int niik_tiff_write_xslice(const char *fname,nifti_image *img,double dmin,double dmax,int xslice);
int niik_tiff_write_yslice(const char *fname,nifti_image *img,double dmin,double dmax,int yslice);
int niik_tiff_write_zslice(const char *fname,nifti_image *img,double dmin,double dmax,int zslice);

int niik_image_rotate(nifti_image *img,char *rstr);
int niik_image_restructure(nifti_image *img,char *rstr);

int niik_image_subsample(nifti_image *img,int sub_ratio);


int niikmat_write_as_linear_xfm(char *fname,niikmat *mat);
int niikmat_convert_from_xfm(niikmat *m,nifti_image *srcimg,nifti_image *tgtimg);
int niikmat_convert_to_xfm(niikmat *m,nifti_image *srcimg,nifti_image *tgtimg);
niikmat *niikmat_read_xfm(char *fname);

/* principal axes */
int niikmat_jacobi(niikmat *m,niikvec *eigval,niikmat *eigvec);
int niik_image_principal_axes(nifti_image *img,nifti_image *maskimg,niikmat *eigvec,niikvec *eigval);


/******************************************************************
 *
 * falcon_lesion_fill.c
 *
 * -2012-06-04, Kunio
 * -used to be in nifti_kunio c-file
 *
 *
 ******************************************************************/

int niik_image_fill_lesion(nifti_image *img, nifti_image *les_mask, nifti_image *wm_mask);
int niik_image_fill_lesion_with_matrix(nifti_image *img, nifti_image *les_mask, niikmat *lmat, nifti_image *wm_mask, niikmat *wmat);
int niik_image_fill_lesion_with_feature(nifti_image *img, nifti_image *les_mask, niikmat *lmat, nifti_image *wm_mask, niikmat *wmat); /* Incomplete */


/******************************************************************
 *
 * falcon_histogram.c
 *
 ******************************************************************/

niikmat *niik_image_histogram_auto(nifti_image *img,nifti_image *maskimg,int num);
int niik_image_histogram(nifti_image *img,nifti_image *maskimg,double *hx,double *hy,int num);
int niik_image_histogram_limits(nifti_image *img,nifti_image *maskimg,double *hx,double *hy,int num);
int niik_image_fit_gaussian(nifti_image *img,nifti_image *maskimg,int avgnum,double *mean,double *stdv,double *err);
double niik_image_fit_gaussian_obj_func(double *v);
double niik_image_histogram_optim_bin_size(nifti_image *img,nifti_image *maskimg,int method);

int niik_image_histogram_matching_test1(nifti_image *refimg,nifti_image *refmask,
                                        nifti_image *img,nifti_image *imgmask,int num);

int niik_image_histogram_matching_same_space(nifti_image *refimg,nifti_image *img);

int niik_image_fit_gaussian_tissues(nifti_image *img,nifti_image *maskimg,double imin,double imax,int num,int avgnum,
                                    double *mean,double *stdv,double *peak,double *err);

int niik_image_histogram_2d(nifti_image *ximg,nifti_image *yimg,nifti_image *maskimg,double xmin,double dx,double ymin,double dy,niikmat *histo,int verbose);

/******************************************************************
 *
 * falcon_median.c
 *
 ******************************************************************/

int niik_image_filter_median_radius(nifti_image *img,nifti_image *maskimg,double radius);
int niik_image_filter_median_radius_dilated_area(nifti_image *maskimg,double radius);
double niik_image_filter_median_radius_ijk(nifti_image *img,nifti_image *maskimg,int i,int j,int k,int kx,int ky,int kz,double radius);
double niik_median_quicksort_double_untouch(double *v,int num);
double niik_median_quicksort_double(double *v,int num);
void niik_median_quicksort_double_func(double *v,int nlo,int nhi,int num);

int niik_image_get_median2(nifti_image *img,nifti_image *mask,double *out);



/******************************************************************
 *
 * falcon_bla.c
 *
 ******************************************************************/
/* niik vector functions (falcon_bla.c) */
niikvec *niikvec_init(int num);
niikvec *niikvec_init_range(double xmin, double xmax, double delta);
niikvec *niikvec_init_from_ascii(char *str);
niikvec *niikvec_free(niikvec *v);
int niikvec_display(niikvec *v);
niikvec *niikvec_copy(niikvec *src);
int niikvec_copy_update(niikvec *src,niikvec *dest);
int niikvec_set_all(niikvec *v,double val);


/* niik matrix functions (falcon_bla.c) */
niikmat *niikmat_init(int row,int col);
niikmat *niikmat_rand(int row,int col);
niikmat *niikmat_free(niikmat *mat);

int niikmat_add_to_a(niikmat *a,niikmat *b);
int niikmat_kmul(niikmat *mat,double k);
int niikmat_kadd(niikmat *mat,double k);
int niikmat_set_all(niikmat *mat,double k);
niikmat *niikmat_identity(int row,int col);
int niikmat_identity_update(niikmat *mat);
int niikmat_display(niikmat *mat);
int mat44_display(mat44 mat);
int niikmat_display_msg(const char *msg, niikmat *mat);
niikmat *niikmat_transpose(niikmat *mat);
niikmat *niikmat_transpose_free(niikmat *mat);
niikmat *niikmat_reshape_free(niikmat *mat,int row,int col);
niikmat *niikmat_reshape(niikmat *mat,int row,int col);
niikmat *niikmat_copy(niikmat *mat);
int niikmat_copy_update(niikmat *src,niikmat *dest);
int niikmat_clear(niikmat *mat);
int niik_image_one(nifti_image *img);
int niik_image_calc_nvox(nifti_image *img);
int niik_image_calc_nvox3d(nifti_image *img);
int niik_image_val(nifti_image *img,double val);
int niik_image_gaussian_noise(nifti_image *img,double FWHM, double mean);
int niikmat_inverse_update(niikmat *mat);
niikmat *niikmat_inverse(niikmat *mat);
int niik_singular_value_decomposition(niikmat *a,double *w, niikmat *v);
int niikmat_pinv(niikmat *mat);
niikmat *niikmat_pseudo_inverse(niikmat *mat);
niikmat *niikmat_pseudo_inverse_free(niikmat *mat);
int niikmat_lu_decompose(niikmat *mat,niikmat *Umat,niikmat *Lmat);

/* niik matrix functions for affine transformations (falcon_bla.c)
 * [ 0] [ 1] [ 2] [ 3]            [0][0]  [0][1]  [0][2]   [0][3]
 * [ 4] [ 5] [ 6] [ 7]   --->>>   [1][0]  [1][1]  [1][2]   [1][3]
 * [ 8] [ 9] [10] [11]            [2][0]  [2][1]  [2][2]   [2][3]
 * [12] [13] [14] [15]            [3][0]  [3][1]  [3][2]   [3][3]
 */
int niikmat_scale_matrix_update(niikmat *mat,double sx,double sy,double sz);
int niikmat_shear_matrix_update(niikmat *mat,double sx,double sy,double sz);
int niikmat_rotate_matrix_update(niikmat *mat,double rx,double ry,double rz);
int niikmat_translate_matrix_update(niikmat *mat,double x,double y,double z);
int niikmat_mat44_matrix_update(niikmat *mat,mat44 m);
int niikmat_flip_matrix_update(niikmat *mat,char dir);

niikmat *niikmat_scale_matrix(double sx,double sy,double sz);
niikmat *niikmat_shear_matrix(double sx,double sy,double sz);
niikmat *niikmat_rotate_matrix(double rx,double ry,double rz);
niikmat *niikmat_translate_matrix(double x,double y,double z);
niikmat *niikmat_flip_matrix(char dir);
niikmat *niikmat_mat44_matrix(mat44 m);
mat44 niikmat_make_mat44(niikmat *mat);
int niikmat_affine_matrix(niikmat *mat,double rx,double ry,double rz,double tx,double ty,double tz,double sx,double sy,double sz,double hx,double hy,double hz,double px,double py,double pz);
niikmat *niikmat_affine_matrix_new(double rx,double ry,double rz,
                                   double tx,double ty,double tz,
                                   double sx,double sy,double sz,
                                   double hx,double hy,double hz,
                                   double px,double py,double pz);
niikmat *niikmat_affine_matrix_val(double m00,double m01,double m02,double m03,
                                   double m10,double m11,double m12,double m13,
                                   double m20,double m21,double m22,double m23,
                                   double m30,double m31,double m32,double m33);

niikmat *niikmat_multiply(niikmat *mat1,niikmat *mat2);  /* output = mat1 * mat2 */
int niikmat_multiply_mat1(niikmat *mat1,niikmat *mat2);  /* here output is replaced in mat1 */
int niikmat_multiply_mat2(niikmat *mat1,niikmat *mat2);
int niikmat_multiply_mat1_free2(niikmat *mat1,niikmat *mat2);  /* here output is replaced in mat1 and mat2 is freed */
int niikmat_multiply_mat2_free1(niikmat *mat1,niikmat *mat2);
niikmat *niikmat_multiply_free12(niikmat *mat1,niikmat *mat2);
niikmat *niikmat_multiply_free1(niikmat *mat1,niikmat *mat2);
niikmat *niikmat_multiply_free2(niikmat *mat1,niikmat *mat2);



/* niik vector sort functions (falcon_bla.c) */
int niikvec_sort(niikvec *v);
int niik_sort_double_vector(double *v,int num);
void niik_vector_sort_func(double *v,int nlo,int nhi,int num);
double niik_get_median_from_sorted_double_vector(double *v,int num);
double niik_get_percentile_from_double_vector(double *v,int num,double percent);
int niik_sort_double_vector_index(double *v,int *idx,int num);
void niik_vector_sort_index_func(double *v,int *idx,int nlo,int nhi,int num);

double niikvec_get_max(niikvec *v);
double niikvec_get_min(niikvec *v);


/* bspline functions (falcon_bla.c) */
double niik_bspline_eval(double *c,int num,double x);
int niik_bspline_calc_coeff_update(niikmat *A,double *v,int num);
int niik_bspline_update_A(niikmat *A);
int niik_bspline_update_A_precomp_matrix(niikmat *A);
int niik_get_mode_double_vector(double *x,double *y,int num,double *xi,double *yi);
int niik_get_mode_bspline_vector(double *v,int num,double *xi,double *yi); /* get the mode using bspline for sub-bin accuracy */
int niik_get_mode_bspline_vector_with_A_matrix(double *v,int num,double *xi,double *yi,niikmat *A_mat);

/* watershed functions for double vector (falcon_bla.c) */
int niik_get_watershed_from_tallest_peaks_double_vector(double *y,int num,int pnum);
int niik_get_watershed_double_vector_untouch(double *y,int *f,int num,int p);
int niik_get_watershed_double_vector_update(double *y,int *f,int num,int p);
int niik_get_watershed_double_vector(double *y,int num,int p);


/*generic function for reading a vector from a csv file*/
void *niik_read_vector_ex(const char *fname,int is_int, int *num, int column_id, const char *column_name);


/* generic double vector algorithms (falcon_bla.c) */
double *niik_calloc_double_vector(int num);
double **niik_calloc_double_matrix(int m,int n);

int niik_write_double_vector(const char *fname,double *v,int num);
int niik_write_double_vector_ex(const char *fname,double *v,int num,const char *colname);

double *niik_read_double_vector(const char *fname,int *num);
int niik_write_double_vectors(const char *fname,double **v,int num,int ncol);
int niik_write_double_vectors_ex(const char *fname,double **v,int num,int ncol,const char **colnames);
int niik_write_double_vector_binary(const char *fname,double *v,int num);
int niik_display_float_vector(float *f,int num);
int niik_display_double_vector(double *v,int num);
int niik_runavg_double_vector(double *v,int num,int anum);
int niik_get_max_index_double_vector(double *v,int num);
int niik_get_min_index_double_vector(double *v,int num);
int niik_get_min_index_float_vector(float *v,int num);
int niik_get_max_index_float_vector(float *v,int num);
double niik_get_max_from_double_vector(double *v,int num);
double niik_get_min_from_double_vector(double *v,int num);
double niik_get_mean_from_double_vector(double *v,int num);
double niik_get_stdv_from_double_vector(double *v,int num);
double niik_get_var_from_double_vector(double *v,int num);
double niik_get_sum_from_double_vector(double *v,int num);
int niik_get_trimmed_average_from_double_vector(double *v,int num,double percent,double *out);
int niik_set_zero_for_double_vector(double *v,int num);

int niik_get_moments_from_double_vector(double *v,int num,double *mean,double *stdv,double *skew,double *kurt,double *adev);


double niik_interp1d_linear_in_double_vector(double *v,int num,double x); /* 1d interpolation */
int niik_central_difference_double_vector(double *v,int num); /* central difference */

int niik_display_stats_for_double_vector(double *v,int num);

/* generic int vector algorithms (falcon_bla.c) */
int *niik_read_int_vector(const char *fname,int *num);
int niik_write_int_vector(const char *fname,int *v,int num);
int niik_write_int_vector_ex(const char *fname,int *v,int num,const char *colname);
int niik_display_int_vector(int *v,int num);
double niik_sum_int_vector(int *v,int num);
int niik_count_nonzero_from_int_vector(int *v,int num);
int niik_count_positive_int_from_int_vector(int *v,int num);
int niik_count_zero_from_int_vector(int *v,int num);
int niik_get_max_from_int_vector(int *v,int num);
int niik_get_min_from_int_vector(int *v,int num);
int niik_get_max_index_from_int_vector(int *v,int num);
int niik_get_min_index_from_int_vector(int *v,int num);


/* affine matrix decomposition (falcon_bla.c) */
int niikmat_decompose_affine(niikmat *mat,double *par,int dof);
double niikmat_decompose_affine_obj_func(double *v);

/* halfway matrix calculation  (falcon_bla.c) */
int niikmat_halfway_matrix(niikmat *afmat, double *afpar,int dof);
double niikmat_halfway_matrix_obj_func(double *v);

niikmat *niikmat_halfway_matrix2(niikmat *A);

/* slope calculation (falcon_bla.c) */
int niik_slope_intercept_from_double_vector(double *x, double *y, int num, double *slope, double *intercept);
int niik_slope_from_double_vector(double *x, double *y, int num, double *slope);
int niik_slope_from_double_vector_least_absolute_difference(double *x,double *y,int num,double *par,int ndeg);

int niik_polynomial_fit(double *x, double *y, int num, int ndeg,double *par);

/* correlation coefficient (falcon_bla.c) */

int niik_double_vector_calc_corrcoef(double *v,double *w,int num,double *out);
int niikvec_calc_corrcoef(niikvec *v,niikvec *w,double *out);
int niikvec_cross_correlation(niikvec *v,niikvec *w,niikvec *out);

/* average matrices  (falcon_bla.c) */
double *niik_average_affine_param(niikmat **A,int num);
niikmat *niikmat_average_from_affine_param(niikmat **A,int num);





/******************************************************************
 *
 * falcon_nelder_mead.c
 *
 ******************************************************************/

int niik_nelder_mead(niikmat *p,int ndim,double *tol,int cost_method,double (* pfn)(),int *maxiter);
/* new multi-seed / multi-level optimization function */
int niik_nelder_mead_multi_level(niikmat *s,niikmat *p,int *nseed,int nlevel,int ndim,double *tol,int cost_method,double (* pfn)(),int *maxiter,int verbose);


/******************************************************************
 *
 * falcon_point.c
 *
 ******************************************************************/

niikpt niikpt_rand();

niikpt niikpt_val(double x, double y, double z, double w);
void niikpt_disp(niikpt p);
char *niikpt_display_xyz(niikpt p);
inline static niikpt niikpt_add(niikpt p,niikpt q);
inline static niikpt niikpt_avg(niikpt p,niikpt q);
inline static niikpt niikpt_wavg(niikpt p,niikpt q,double wp);
inline static niikpt niikpt_avg3(niikpt p,niikpt q,niikpt r);
inline static niikpt niikpt_wavg3(niikpt p,niikpt q,niikpt r, double wp, double wq, double wr);
inline static niikpt niikpt_sub(niikpt p,niikpt q);
inline static niikpt niikpt_mul(niikpt p,niikpt q);
inline static niikpt niikpt_div(niikpt p,niikpt q);
inline static double niikpt_dot(niikpt p,niikpt q);
inline static double niikpt_mag2(niikpt p);
inline static double niikpt_mag (niikpt p);
inline static double niikpt_distance(niikpt p,niikpt q);
inline static niikpt niikpt_cross(niikpt p,niikpt q);
inline static niikpt niikpt_kmul(niikpt p,double k);
inline static niikpt niikpt_move_normal(niikpt p,niikpt n,double k);
niikpt niikpt_interp_tri(niikpt p1,niikpt p2,niikpt p3,double t1,double t2);
double niikpt_plane_to_point_distance(niikpt pt,niikpt pt_on_plane,niikpt plane_normal);
double niikpt_plane_to_point_distance2(niikpt pt,niikpt plane_eq);
niikpt niikpt_plane_eq_to_plane_normal(niikpt plane);
niikpt niikpt_3pts_to_plane_eq(niikpt p1,niikpt p2,niikpt p3);
niikpt niikpt_closest_point_on_plane_to_point(niikpt pt,niikpt plane);
double niikpt_line_to_point_distance(niikpt pt,niikpt p1,niikpt p2);
niikpt niikpt_closest_point_on_line_to_point_distance(niikpt pt,niikpt p1,niikpt p2);
niikpt niikpt_line_plane_intersection(niikpt p1,niikpt p2,niikpt p3,niikpt L1,niikpt L2);
niikpt niikpt_closest_point_on_triangle_to_point(niikpt pt,niikpt p1,niikpt p2,niikpt p3);

niikpt *niikpt_alloc(int num);
niikpt **niikpt_matrix(int col,int row);
niikpt niikpt_image_get_pixdim(nifti_image *img);
int niikpt_point_is_on_triangle(niikpt pt,niikpt t1,niikpt t2,niikpt t3,niikpt v1,niikpt v2,niikpt v3,niikpt normal);


niikpt niikpt_image_get_centroid(nifti_image * img,nifti_image *maskimg);
int niikpt_image_get_centroid_update(nifti_image * img,nifti_image *maskimg,niikpt *pt);
int niikpt_image_get_max_ijk_position(nifti_image *maskimg,niikpt *pmax);
int niikpt_image_get_min_ijk_position(nifti_image *maskimg,niikpt *pmin);
double niikpt_image_get_bounding_sphere(nifti_image * img,double threshold,niikpt ctr);


/* spherical coordinates falcon_point.c */
niikpt niikpt_euc2sph(niikpt p);
niikpt niikpt_sph2euc(niikpt p);

#include "falcon_point_inline.h"
#include "falcon_sph_inline.h"



/******************************************************************
 *
 * falcon_interpolate.c
 *
 ******************************************************************/

const char * niik_interpolate_string(int interp);
const char * niik_interpolate_name(int interp);
double niik_image_interpolate_3d(nifti_image *img,niikpt p,int interp);
double niik_image_interpolate_3d_ijk(nifti_image *img,niikpt p,int interp);
double niik_image_interpolate_3d_nn(nifti_image *img,niikpt p);
double niik_image_interpolate_3d_nn_ijk(nifti_image *img,niikpt p);
unsigned char niik_image_interpolate_uint8_image_3d_nn(nifti_image *img,double x, double y,double z);
double niik_image_interpolate_3d_linear(nifti_image *img,niikpt p);
double niik_image_interpolate_3d_linear_ijk(nifti_image *img,niikpt p);
double niik_image_interpolate_double_image_3d_nn(nifti_image *img,double x, double y,double z);
int niik_image_interpolate_convert_3d_bspline_coeff(nifti_image *img);
nifti_image *niik_image_interpolate_convert_3d_bspline_coeff2(nifti_image *img);
int niik_image_interpolate_inverse_3d_bspline_coeff(nifti_image *coeffimg, nifti_image *img);
double niik_image_interpolate_3d_bspline(nifti_image *bsplimg,niikpt p);
double niik_image_interpolate_3d_bspline_ijk(nifti_image *bsplimg,niikpt p);
double niik_image_interpolate_float_image_3d_linear(nifti_image *img,niikpt p);
double niik_image_interpolate_float_image_3d_linear_ijk(nifti_image *img,niikpt p);

/* input image can be more than 3d but interpolation is 3d (falcon_interpolate.c */
int niik_image_interpolate_3d_ijk_update(nifti_image *img,niikpt p,int interp,double *v);
int niik_image_interpolate_3d_xyz_update(nifti_image *img,niikpt p,int interp,double *v);
int niik_image_interpolate_3d_nn_ijk_update(nifti_image *img,niikpt p,double *v);
int niik_image_interpolate_3d_linear_ijk_update(nifti_image *img,niikpt p,double *v);
int niik_image_interpolate_3d_linear_xyz_update(nifti_image *img,niikpt p,double *v);
int niik_image_interpolate_3d_bspline_ijk_update(nifti_image *bsplimg,niikpt p,double *v);
int niik_image_interpolate_3d_bspline_xyz_update(nifti_image *bsplimg,niikpt p,double *v);


int niik_image_interpolate_cost(nifti_image *img,nifti_image *refimg,niikmat *afmat,double *cost);

/* interpolate along the normal */
int niik_image_interp_along_normal(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v,niikvec *out);
int niik_image_interp_along_normal_update(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v);
int niik_image_interp_along_normal_double_vector(nifti_image *img,int interp,niikpt pt,niikpt normal,double *x,double *out,int num);




/******************************************************************
 *
 * falcon_affine.c
 *
 ******************************************************************/

int niik_image_inverse_affine_transform_3d_update(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp);

int niik_image_affine_transform_3d_update(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp);
int niik_image_resample_3d_update(nifti_image *img,double dx, double dy, double dz,int nx, int ny,int nz,int interp);
nifti_image *niik_image_affine_transform_3d(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp);
nifti_image *niik_image_resample_3d(nifti_image *img,double dx, double dy, double dz,int nx, int ny,int nz,int interp);

int niik_image_copy_pixdim_info(nifti_image *src,nifti_image *dest);
int niik_image_copy_dim_info(nifti_image *src,nifti_image *dest);

int niik_image_update_sform(nifti_image *img,niikmat *afmni);


/******************************************************************
 *
 * falcon_aregister.c
 *
 ******************************************************************/

int niik_image_aregister_set_g_img(nifti_image *img,int idx);
int niik_image_aregister_set_g_reg_method(int method);
int niik_image_aregister_set_g_img_to_null();
int niik_image_aregister_set_g_invmat(int imat);
int niik_image_aregister_get_g_invmat();
int niik_image_aregister_set_g_verbose_errfunc(int v);

char *niik_aregister_method_string( int reg_method );
int niik_aregister_display_affine( double *affpar );
niikmat *niik_aregister_matrix_from_affpar(double *affpar);
int      niik_aregister_matrix_from_affpar_update( niikmat *afmat, double *affpar );
double   niik_aregister_align_mni_maximize_symmetry(double *v);
int niik_aregister_align_mni(nifti_image *mni_img, nifti_image *mni_seg, nifti_image *img, nifti_image *seg, double *affpar, int cost_method, double filFWHM, double sample);
int niik_aregister_align_mni_default1(nifti_image *img, nifti_image *seg, double *affpar);
int niik_aregister_align_mni_predefined_imgs(nifti_image *img, nifti_image *seg, double *affpar,
    int cost_method, double filFWHM, double sample);

int niik_image_aregister(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,nifti_image *movseg,
                         double *affpar,double *perturb,int register_method,double sample,double filFWHM);
int niik_image_aregister_imat(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,nifti_image *movseg,niikmat *imat,
                              double *affpar,double *perturb,int register_method,double sample,double filFWHM);

int niik_aregister_nmi_obj_display_var();
int niik_aregister_nmi_obj_display(nmi_obj *nmiobj);
int niik_aregister_nmi_obj_update_var(double min1,double max1,double min2,double max2,int num1,int num2);
int niik_aregister_nmi_obj_update(nmi_obj *obj,double min1,double max1,double min2,double max2,int num1,int num2);
void niik_aregister_nmi_obj_free(nmi_obj *obj);
nmi_obj *niik_aregister_nmi_obj_init(double min1,double max1,double min2,double max2,int num1,int num2);
nmi_obj *niik_aregister_nmi_obj_alloc();
nmi_obj *niik_aregister_g_nmi_obj_get();
int niik_aregister_g_nmi_obj_set(nmi_obj *nmiobj);


/******************************************************************
 *
 * falcon_aregister2.c
 *
 * -2012-05-01, Kunio
 * -under-development, to replace falcon_aregister.c
 *
 *
 ******************************************************************/


void niik_aregister2_set_verbose( int verbose );
int  niik_aregister2_get_verbose();

nmi_obj *niik_aregister2_nmi_obj_alloc();
void niik_aregister2_nmi_obj_free(nmi_obj *obj);
nmi_obj *niik_aregister2_nmi_obj_init(double min1,double max1,double min2,double max2,int num1,int num2);
void niik_image_aregister_set_nmi_obj(nmi_obj *obj);
nmi_obj *niik_image_aregister_get_nmi_obj();

char *niik_aregister2_method_string( int reg_method );

int niik_image_aregister2_cost_cc(nifti_image *refimg,nifti_image *refmask,nifti_image *img,double *cc);
int niik_image_aregister2_cost_nmi(nifti_image *refimg,nifti_image *refmask,nifti_image *img,nmi_obj *nmi);
int niik_image_aregister2_hemispheric_symmetry_optiization(nifti_image *img,double FWHM,niikpt imgctr,double *ry, double *rz);

int niik_image_aregister2_test1(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int dof,int areg_cost,nmi_obj *nmiobj);
/*int niik_image_aregister2_test2(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,nmi_obj *nmiobj);*/
int niik_image_aregister2_test2(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,
                                int nlevel,double fFWHM,double fmaxiter,int fnseed,double fdelta,niikvec *idaffpar,nmi_obj *nmiobj);
int niik_image_aregister2_test_range(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,int areg_cost,
                                     int nlevel,double fFWHM,double fmaxiter,int fnseed,double fdelta,niikvec *range,niikvec *idaffpar,nmi_obj *nmiobj);


int niik_image_aregister2(nifti_image *refimg,nifti_image *refseg,nifti_image *movimg,niikmat *imat,
                          double *affpar,double *daffpar,int maxiter, double tol,int register_method);
int niik_image_aregister2_multilevel(nifti_image *refimg,nifti_image *refseg,nifti_image *img,niikmat *afmat,double *affpar,
                                     int areg_cost,nmi_obj *nmiobj,
                                     int nlevel,niikvec *FWHM, niikvec *delta, niikvec *tol,int *nseed,int *maxiter,niikmat *daffpar,niikmat *seed);
int niik_image_aregister2_multilevel_multimask(nifti_image *refimg,nifti_image **refseg,nifti_image *img,niikmat *afmat,double *affpar,
    int areg_cost,nmi_obj *nmiobj,
    int nlevel,niikvec *FWHM, niikvec *delta, niikvec *tol,int *nseed,int *maxiter,
    niikmat *daffpar,niikmat *seed);
int niik_image_aregister2_multiseed(nifti_image *refimg,nifti_image *refseg,
                                    nifti_image *img,niikmat *afmat,
                                    double *affpar,int areg_cost,nmi_obj *nmiobj,
                                    double FWHM, double delta, double tol,
                                    int maxiter, double *daffpar, niikmat *seed);


/******************************************************************
 *
 * falcon_nregister.c
 *
 ******************************************************************/

char *niik_warp_type_string( int code );
int niik_image_apply_3d_warp_update(nifti_image *img,nifti_image *refimg,nifti_image *warp_img,int warp_type,int interp);
int niik_image_nregister_invert_nonlinear_map(nifti_image *warpimg,int warp_type);
int niik_image_nregister_invert_displacement_map(nifti_image *warpimg);
int niik_image_nregister_invert_displacement_map_iterative(nifti_image *warpimg,int maxiter);

/* functions to register (falcon_nregister_bspline.c) */
int niik_image_nregister_bspline(nifti_image *refimg,nifti_image *refseg,
                                 nifti_image *movimg,nifti_image *movseg,
                                 nifti_image *wbsimg,
                                 int flag_sobel,
                                 double *idelta,double *wdelta,double *filFWHM,int *maxiter,
                                 double *maxdfm,int nlevel,int localsize,int ndir,
                                 double reg_dist,double grad_weight,double dfm_fc);
int niik_image_nregister_bspline_update_warp_coeff(nifti_image *img,nifti_image *refimg,double delta,double regularization_size);


/* functions to register (falcon_nregister_demons.c) */
int niik_image_nregister_demons_update(nifti_image *warpimg,nifti_image *gradimg,
                                       nifti_image *refimg,nifti_image *maskimg,
                                       nifti_image *tgtimg,niikmat *tgtmat,
                                       nifti_image *movimg,niikmat *movmat);
int niik_image_nregister_demons(nifti_image *warpimg,nifti_image *refimg,nifti_image *maskimg,
                                nifti_image *tgtimg,niikmat *tgtmat,nifti_image *tgtgrad,
                                nifti_image *movimg,niikmat *movmat,
                                int maxiter,   /* for intra-subject registration, small number is enough */
                                double wFWHM);
int niik_image_nregister_intrasubject_apply_warp(nifti_image *warpimg,nifti_image *refimg,nifti_image *movimg,niikmat *movmat,int interp);


/* non-essential (not well written) functions (falcon_nregister.c) */
int niik_image_3d_warp_with_ref_update(nifti_image *img,nifti_image *refimg,nifti_image *pos_img,int interp);
int niik_image_3d_warp_update(nifti_image *img,nifti_image *pos_img,int interp);
nifti_image *niik_image_combine_warp(nifti_image *movimg,nifti_image *refimg,nifti_image *warpimg,niikmat *premat, niikmat *postmat,int warp_map_type);


/* functions to register (falcon_nregister_intrasubject_template.c) */
int niik_image_nregister_average_with_fov(nifti_image *refimg,nifti_image **warpimg,nifti_image **imglist,niikmat **matlist,int num,int interp);
int niik_image_nregister_intrasubject_template(nifti_image *refimg,nifti_image *maskimg,
    nifti_image **warpimglist,nifti_image **imglist,niikmat **matlist,int num,
    niikvec *gFWHM,niikvec *wFWHM,
    niikvec *maxiter,
    int nlevel);
int niik_image_nregister_intrasubject_template_test1(nifti_image *refimg,nifti_image *maskimg,
    nifti_image **warpimglist,nifti_image **imglist,niikmat **matlist,int num);



/******************************************************************
 *
 * falcon_bounds.c
 *
 ******************************************************************/

int niik_image_boundary(nifti_image *img,nifti_image *maskimg,double *dval, int flag_dim_check, int flag_color, int flag_inclusive);

nifti_image *niik_image_mask_add_red_color_uint8(nifti_image *img,nifti_image *maskimg);
nifti_image *niik_image_mask_add_red_color(nifti_image *img,nifti_image *maskimg,double scale_max);

/******************************************************************
 *
 * falcon_threshold.c
 *
 ******************************************************************/

int niik_image_threshold(nifti_image *img,double thresh);
nifti_image *niik_image_threshold_new(nifti_image *img,double thresh);
int niik_image_thresh_otsu(nifti_image * img,nifti_image *maskimg,double *thresh);
int niik_image_thresh_ridler(nifti_image * img,nifti_image *maskimg,double *thresh);
int niik_image_thresh_kittler(nifti_image * img,nifti_image *maskimg,double *thresh);
double niik_image_get_upper_threshold(nifti_image *img,nifti_image *maskimg,double fac);


/******************************************************************
 *
 * falcon_sobel.c
 *
 ******************************************************************/

int          niik_image_sobel_filter_update(nifti_image *img,char dir);
niikpt       niik_image_sobel_filter_niikpt(nifti_image *img,niikpt p);
niikpt       niik_image_sobel_filter_voxel_par(nifti_image *img,int p,double a1,double a2,double a3);
niikpt       niik_image_sobel_filter_voxel(nifti_image *img,int p);
nifti_image *niik_image_sobel_filter(nifti_image *img,char dir);
nifti_image **niik_image_sobel_filters(nifti_image *img);
nifti_image **niik_image_sobel_filters_with_mag(nifti_image *img);
nifti_image *niik_image_sobel_filters_with_mag_single_output(nifti_image *img);

nifti_image *niik_image_gauss_sobel_filters_with_mag_single_output(nifti_image *img,double FWHM);
nifti_image *niik_image_gauss_sobel_filter(nifti_image *img,double FWHM);

nifti_image *niik_image_non_maximum_suppression(nifti_image *img,nifti_image *maskimg,double FWHM,double mean, double stdv,double thresh,double pthresh);
nifti_image *niik_image_divergence(nifti_image **grad_with_mag,int norm);


/******************************************************************
 *
 * falcon_laplace_map.c
 *
 ******************************************************************/

nifti_image *niik_image_laplace_map(nifti_image *Simg, nifti_image *Limg);



/******************************************************************
 *
 * falcon_distance_map.c
 *
 ******************************************************************/

nifti_image *niik_image_distance_map(nifti_image *img,double max_dist);


/******************************************************************
 *
 * falcon_spmat.c
 *
 * -functions for sparse matrix
 *
 ******************************************************************/

int niik_spmat_construct(double **A,int N,double *sa,int *ija);
int niik_spmat_reverse(double **A,int N,double *sa,int *ija);
int niik_spmat_op_mul_mat_vec(double *sa,int *ija,int N,double *x,double *y);
double niik_spmat_op_vec_dot_product(double *a,double *b,int N);
int niik_spmat_conj_grad(double *sa,int *ija,int N,double *x,double *b);
void niik_spmat_display_dig(double **A,int n1,int n2,int dig);


/******************************************************************
 *
 * falcon_nbcsr.c
 *
 ******************************************************************/

int niik_aregister_nbcr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg, double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag);
int niik_aregister_nbcsr(nifti_image *refimg, nifti_image *refseg, nifti_image *movimg, nifti_image *movseg, double *affpar, int cost_method, double filFWHM, double sample,int rigid_flag);

void nifti1_kunio_nbcsr_turn_on_debug();
void nifti1_kunio_nbcsr_turn_off_debug();






/* falcon_inorm.c */

int niik_image_linear_normalization(nifti_image *refimg, nifti_image *modimg, nifti_image *modmask, niikmat *mod2ref);
int niik_image_linear_normalization_get_factor(nifti_image *refimg, nifti_image *modimg, nifti_image *modmask, niikmat *mod2ref,double *slope);



/* falcon_dbc.c */

nifti_image *niik_image_dbc(nifti_image *refimg, nifti_image *img, nifti_image *maskimg, niikmat *img2ref, double radius, double lim);
nifti_image *niik_image_dbc_with_scaling(nifti_image *refimg, nifti_image *img, nifti_image *maskimg,
    niikmat *img2ref, double radius, double lim, int flag_scale);
nifti_image *niik_image_dbc_with_scaling_resample(nifti_image *refimg, nifti_image *img, nifti_image *maskimg,
    niikmat *img2ref, double radius, double lim, int flag_scale, double *xyz) ;

/*
 * falcon_feature_map.c
 *
 * -functions for feature maps
 */

double niik_image_voxel_gabor_filter(nifti_image *img,int x,int y,int z,double lambdax,double lambday,double lambdaz,double sigmax,double sigmay,double sigmaz,double rx,double ry,double rz);

int niik_image_feature_voxel_type1(nifti_image *img,int voxel,int kernel,niikvec *fv);
nifti_image *niik_image_feature_map_type1(nifti_image *img,int kernel);
nifti_image *niik_image_feature_map_type2(nifti_image *img);
nifti_image *niik_image_feature_map_type3(nifti_image *img,int kernel,int method);
nifti_image *niik_image_feature_map_type4(nifti_image *img);
nifti_image **niik_image_feature_map_mean_var(nifti_image *img,int kernel);

int niik_image_feature_voxel_type2(nifti_image *img,int voxel,niikvec *fv);
int niik_image_feature_voxel_type1(nifti_image *img,int voxel,int kernel,niikvec *fv);

nifti_image *niik_image_feature_map_type4_multi(nifti_image **imglist,int nimg);

/*
 * falcon_bimodal_fit.c
 */

int niik_bimodal_fit_double_vector(double *x,double *y,int num,double *mean,double *stdv,double *peak,double *err,int verbose);
int niik_image_bimodal_fit(nifti_image *img,nifti_image *mask,double imin,double delta,double imax,int num,double *mean,double *stdv,double *peak,double *err);



/*
 * falcon_PABIC.c
 */

double niik_image_pabic_valley(double d, double s);
nifti_image *niik_image_pabic_bias(nifti_image *img,nifti_image *maskimg,int nd,double dist);



/*
 * falcon_cortex_edge.c
 */

int niik_image_niikcortex_atrophy_edge(nifti_image *avgimg,nifti_image **imglist,niikmat **matrixlist,nifti_image *maskimg,niikmat *tstat,double xc_thresh,niikmat *atv,int *ijk,int verbose);


/******************************************************************
 *
 * falcon_jacobian_map.c
 *
 * -2012-09-28, Kunio
 * -not validated as of 2012-10-05
 * -2012-11-14, Kunio: appears to be working properly
 *
 *
 ******************************************************************/

nifti_image *niik_image_jacobian_map(nifti_image *warpimg,nifti_image *maskimg,int warptype);


/******************************************************************
 *
 * falcon_sienax.c
 *
 * -2012-10-05, Kunio
 * -not validated as of 2012-10-05
 *
 *
 ******************************************************************/

int niik_image_siena(nifti_image *templateimg,nifti_image *maskimg,nifti_image **imglist,niikmat **matrixlist,int numimglist,double ulim,double uran,int coloroutput,niikmat *out);



/******************************************************************
 *
 * falcon_pseudoT2.c
 *
 * -2012-11-07, Kunio
 * -not validated as of 2012-11-07
 *
 *
 ******************************************************************/

int niik_image_calculate_pseudoT2(nifti_image *PDimg,double PD_TE,double PD_TR,nifti_image *T2img,double T2_TE,double T2_TR,double omax,nifti_image *maskimg,nifti_image *outimg);


/******************************************************************
 *
 * falcon_interpacket.c
 *
 * -used to simulate inter-packet motion
 *
 ******************************************************************/

nifti_image *niik_image_interpacket_motion_correction(nifti_image *img,int num,int dof);
nifti_image *niik_image_merge_interpacket(nifti_image **imglist,int num);
nifti_image **niik_image_split_interpacket(nifti_image *img,int num);


/******************************************************************
 *
 * falcon_mnc2nii.c
 * falcon_nii2mnc.c
 *
 ******************************************************************/

nifti_image *niik_image_read_minc(const char *minc_name);
int niik_image_write_minc(char *outname,nifti_image *nii_ptr);

nifti_image *niik_image_read_minc2(const char *minc_name);
int niik_image_write_minc2(char *outname,nifti_image *nii_ptr);


int niik_mat44_to_cosines_start_step(mat44 mat,double dx,double dy,double dz,double *ocosx,double *ocosy,double *ocosz,double *ostart,double *ostep);



/******************************************************************
 *
 * falcon_segment.c
 *
 ******************************************************************/

int niik_image_classify_FGMM(nifti_image *img,nifti_image *mask,int nclass,niikmat *stats,int maxiter);



/******************************************************************
 *
 * falcon_nlseg.c
 *
 ******************************************************************/

int niik_image_nlseg_extract_3dpatch(void *img,void *patch,int datatype,int x,int y,int z,int size,int sx,int sy,int sz);
int niik_image_nlseg_extract_3dpatch_mean_var(void *img,void *patch,int datatype,int x,int y,int z,int size,int sx,int sy,int sz,float *M,float *V);



#endif /* __FALCON_H__ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/