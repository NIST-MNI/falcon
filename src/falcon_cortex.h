#pragma once
#ifndef __FALCON_CORTEX_H__
#define __FALCON_CORTEX_H__

#include "falcon.h"
#include "falcon_surfaces.h"


enum {
  CORTEX_ICS=0, /*inner cortical surface (interface between WM and GM) aka white surface*/
  CORTEX_OCS=1  /*outer cortical surface (interface between GM and CSF) aka pial surface*/
};

enum {
  WEIGHT_DT     =0,             /* differential equation solver step*/
  WEIGHT_SURFACE=1,             /* surface weights */
  WEIGHT_THICKNESS_SMOOTHNESS=2,/* thickness smoothness weights */
  WEIGHT_IMAGE=3,               /* image weights */
  WEIGHT_BRAIN_MASK=4,          /* brain mask weights */
  WEIGHT_VENTRICLE=5,           /* ventricle mask weights */
  WEIGHT_LESION=6,              /* lesion mask weights */
  WEIGHT_PROXIMITY=7,           /* proximity weights */
  WEIGHT_AVOID=8,               /* weight for avoid mask */
  WEIGHT_ABS_THICKINESS=9,      /* absolute thickness */
  WEIGHT_GRADIENT=10,           /* gradient term */
  WEIGHT_CURVATURE=11,          /* surface curvature term*/
  WEIGHT_DGRADIENT=12,          /* 2nd order gradient term */
  WEIGHT_MIN_THICKNESS=13,      /* minimal thickness */
  WEIGHT_MAX_THICKNESS=14,      /* minimal thickness */
  WEIGHT_THICKNESS_SIGMA=15,    /* thickness falloff function  */
  WEIGHT_UPDATE_SIGMA=16,       /* update field smoothing  */
  WEIGHT_SURFACE_SIMPLE=17,     /* Simple surface smoothing  */
  WEIGHT_PRIOR=18,              /* Weight of prior  */
  WEIGHT_MIX=19,                /* Weight of tissue intensities mix (i.e GM & WM or GM & CSF)*/
  WEIGHT_COUNT                  /* total number of weght terms*/
};


/* NUMERICAL METHOD DEFINITIONS */
enum {
  NIIK_NUM_METHOD_UNKNOWN=0,         /* unknown, problem */
  NIIK_NUM_METHOD_FORWARD_EULER=1,   /* forward euler */
  NIIK_NUM_METHOD_BACKWARD_EULER=2,  /* backward euler */
  NIIK_NUM_METHOD_MIDPOINT=3,        /* midpoint */
  NIIK_NUM_METHOD_RUNGE_KUTTA=4,     /* runge-kutta */
  NIIK_NUM_METHOD_MAX                /* */
};

typedef unsigned char indicator_t;


typedef struct {
  nifti_image *t1img;
  nifti_image *nonctx_mask;
  nifti_image *avoid_mask;
  nifti_image *gwi_mask;
  nifti_image *csf_mask;
  nifti_image *brain_mask;
  nifti_image *lesion_mask;
  nifti_image **prior;

  kobj **ctx;
  double delta;
  double apply_step;
  niikmat *weight;
  int *ctx_label;
  int iter;
  int iter2;
  int remesh;
  int remesh_fail_max;
  int use_mf;

  double gradient_FWHM[2];
  double divergence_FWHM[2];
  double prior_FWHM[3];

  double gradient_lambda[2];
  double *tissue_val;
  int bbox_depth;
  int cortex_id;
  int numerical_method;
  double proximity_min_distance;
  double tolerance;
  int debug_keep_tmp;
  niikmat *regmat;

  double *debug_pt;

  /*created during execution*/
  nifti_image ***grad_img;
  nifti_image ***grad_prior;
  nifti_image **div_prior;

  nifti_image **div_img;
  bbox *bb;
  double *thklist;
  double *mf_list;
  double *dd[2],*rd[2];
  double mflim;
  double mf_lim[2]; /* multiplication factor lo/hi limits */

  off_curvature_t crv[2];
  kvert ***vmat;
  int nvert;

  /*indicator helper function, using in paralell apply deformation only*/
  indicator_t *indicator; 

  niikpt **dfmlist; /* calculated deformation */
  
  /*CONVERGENCE*/
  FILE *convergence_log;

  /*maximum deformation distance*/
  double *dfm_limit;
  /*actual distance travelled*/
  double *travel;

  /*DEBUG*/
  FILE *debug_trace_log;
  int trace_pt_id;
} niikcortex_deform;


/*TODO: join into niikcortex_deform?*/
typedef struct {
  nifti_image *img;
  nifti_image ***grad_img;
  nifti_image ***grad_prior;
  nifti_image **div_prior;
  nifti_image **div_img;
  nifti_image *brain_mask;
  nifti_image *gwi_mask;
  nifti_image *ventricle_mask;
  nifti_image *avoid_mask;
  nifti_image *lesion_mask;
  nifti_image **prior;
  kvert ***vmat;
  bbox *bb;
  niikmat *deform_weights;

  double *thklist;
  double *mf_list;
  double mean_ics_value;
  double range_ics_value;
  double mean_ocs_value;
  double range_ocs_value;
  double intWM;

  double* gradient_lambda;
  double prox_min_dist;

  off_curvature_t *crv;
  niikcortex_deform *dfm;
} refine_info;



char *niik_numerical_method_string( int code ) ;

/* cortical deformation (falcon_cortex.c) */
int niikcortex_deform_cortex(niikcortex_deform * dfm);

niikcortex_deform *niikcortex_deform_init();
niikcortex_deform *niikcortex_deform_free(niikcortex_deform *dfm);

/* measurement of cortical thickness (falcon_cortex.c) */
double niikcortex_calc_area_weighted_thickness(kobj *ics, kobj *ocs, double *thk, int filter_type, int filter_size);
int niikcortex_calc_thickness(kobj *ics, kobj *ocs, double *thk, double *phi, double *the, int filter_type, int filter_size);
int niikcortex_add_thickness_color(kobj *ics, kobj *ocs, double *thk, double omin, double omax);
int niikcortex_add_color(kobj *obj, double *var, double omin, double omax,int color_map_type,int color_levels);
int niikcortex_add_color_discrete(kobj *obj, double *var, double omin, double omax, int color_map_type,int color_levels);

/* falcon_cortex_initics.c */
int    niikcortex_initics_expand(nifti_image *brain_img,double max_size,double step_size,kobj *obj);
int    niikcortex_initics_shrink(nifti_image *gwi_img,nifti_image *lap_map,nifti_image *dist_map,
                                 int maxiter,double  initlen, double finlen, double *dfm_step,int dfm_iter,
                                 kobj *obj,niikpt check_pt);
kobj  *niikcortex_initics_atimg(nifti_image *at_ics_img,double elen);
double niikcortex_initics_shrink_check_deform(nifti_image *img,bbox *bb,kvert *v,double step,double thresh);
int    niikcortex_initics_shrink_check_deform_interp(nifti_image *img,kvert *v,double thresh);


/* falcon_cortex_initocs.c */
int niikcortex_initocs_expand(nifti_image *img,nifti_image *brain_mask,double *step_size,
                              int stepiter,double intGM,double intWM,double intCSF,double intICS,double intOCS,
                              double gthresh, double max_cth,double init_cth,
                              int smooth_iter_maxcth,int smooth_iter,int maxiter2,
                              float smooth_percent,kobj *ics, kobj *ocs,nifti_image *border,int alt_mode);





/******************************************************************
 *
 * falcon_off_shrinkwrap.c
 *
 ******************************************************************/

int off_shrinkwrap_kobj_simple(nifti_image *img,kobj *obj,double len);
int off_shrinkwrap_kobj(nifti_image *img,kobj *obj,double len,int maxiter,double stepsize,int debug,double smooth);
int off_shrinkwrap_kobj_bbox(nifti_image *img,kobj *obj,double len,int maxiter,double stepsize,int flag_bbox,int debug,double smooth);
int off_shrinkwrap_kobj_bbox_remesh(nifti_image *img,kobj *obj,double len,int maxiter,double stepsize,int flag_bbox,int flag_remesh,int debug,double smooth);
int off_kobj_balloon(kobj *obj, int maxiter, double stepsize, int debug);


/******************************************************************
 *
 * falcon_cortex_*.c
 *
 ******************************************************************/

int niikcortex_get_debug();
void niikcortex_set_debug(int debug);

int niik_image_average_with_fov(nifti_image *refimg,nifti_image **imglist,niikmat **matlist,int num,int interp);



/* falcon_cortex.c */
int niikcortex_estimate_tissue_values(nifti_image *img,nifti_image *brain_mask,
                                      nifti_image *csf_mask,nifti_image *GWI_mask,
                                      double *iCSF, double *iWM, double *iGM,
                                      double *iBrain,double *intICS,double *intOCS,
                                      double *ranICS,double *ranOCS,double mixICS,double mixOCS);

int niikcortex_make_fuzzy(nifti_image *input,nifti_image *output,
                          double value,double range);


int niikcortex_off_count_intersection(bbox *bb,kobj *ics, kobj *ocs);

int niikcortex_off_correct_self_intersections(bbox *bb,kobj *obj[],int maxiter);

nifti_image *niikcortex_gwi_mask_get_CLADA_mask(char *filename,nifti_image *refimg,nifti_image *warpimg,int interp);

int *niikcortex_non_cortex_label(int *ctx_label, nifti_image *nctx_img,kobj *ics, kobj *ocs,double dran,double delta,int iter);
nifti_image *niikcortex_gwi_mask(nifti_image *img,
                                 nifti_image *brain_mask,
                                 nifti_image *ven_mask, double ven_dilate,
                                 nifti_image *wm_mask,
                                 double dgm_dilate,
                                 nifti_image *avoid_mask,
                                 nifti_image *lesion_mask,
                                 nifti_image *warpimg,double vessel_radius,
                                 double median_radius);

nifti_image *niikcortex_modify_image(nifti_image *img,
                                     nifti_image *brain_mask,
                                     nifti_image *ven_mask, double ven_dilate, double ven_blur,
                                     nifti_image *avoid_mask,
                                     nifti_image *dgm_mask, double dgm_dilate, double dgm_blur,
                                     nifti_image *lesion_mask,
                                     double fillval);

nifti_image *niikcortex_wm_mask(nifti_image *img,nifti_image *maskimg,nifti_image *warpimg,double mag,double thresh);
/* white matter segmentation: testing version
 * -includes the spatially varying threshold
 */



/*******************************************************************
 * 
 * falcon_slice.c
 * 
 * *****************************************************************/
typedef struct 
{
    int slice_dir; /*0 - x, 1 -y , 2- z*/
    int slice_num;
    int *slice;
    double dmin;
    double dmax;
    double scale;
    const char *fpattern;
    int compression;
} extract_slice_info;

typedef struct
{
  extract_slice_info slices_x;
  extract_slice_info slices_y;
  extract_slice_info slices_z;

  int dump_slices;
  int dump_surfaces;
  char *prefix;
} cortex_tracing_info;

int niik_image_write_slices_obj(extract_slice_info *info,
        nifti_image *img, bbox *off_bb);


/*******************************************************************
 * 
 * falcon_cortex_debug.c
 * 
 * *****************************************************************/
int  falcon_tracing_init(nifti_image * img, cortex_tracing_info *info);
void falcon_tracing_free(cortex_tracing_info *info);
void falcon_tracing_dump(cortex_tracing_info *info, int iter, const char *task, nifti_image * img,bbox *bb);
void falcon_tracing_dump_objects(cortex_tracing_info *info, int iter, const char *task,kobj *obj[],int nobj);



/**/


#endif /*__FALCON_CORTEX_H__*/


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/