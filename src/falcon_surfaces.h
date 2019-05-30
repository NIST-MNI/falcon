#pragma once
#ifndef __FALCON_SURFACES_H__
#define __FALCON_SURFACES_H__


/********************************************************
 *
 *  handling off files with vertex/triangle(face)/edge
 *
 ********************************************************/

typedef struct  vertex_struct kvert; /* vertex */
typedef struct  face_struct kface;   /* face = triangle */
typedef struct  edge_struct kedge;   /* edge */
typedef struct  obj_struct kobj;     /* polyhedron */

struct face_struct { /* always triangle */
  int index;         /* index = should be positive */
  kvert *vert[3];    /* triangle's vertices */
  kedge *edge[3];    /* triangle edges */
  niikpt normal;        /* triangle normal */
  niikpt pmin,pmax;     /* min/max values */
  double *color;     /* color if usd, or NULL */
  kface *prev,*next; /* double linked list */
};                   /* functions in falcon_off_*.c use this */

struct edge_struct {
  int   index;         /* index = should be positive */
  kvert *endpts[2];  /* end points of the edge */
  kface *adjface[2]; /* adjacent triangles */
  kedge *prev,*next; /* double linked list */
};                   /* functions in falcon_off_*.c use this */

struct vertex_struct {
  int index;          /* index = should be positive */
  niikpt v;           /* 3d coordinates */
  niikpt normal;      /* vertex normal */
  niiksph sph;        /* spherical coordinates */
  int nei;            /* number of neighbors */
  kvert **neivert;    /* list of pointers to the neighboring vertice */
  kface **neiface;    /* list of pointers to the neighboring triangles */
  kedge **neiedge;    /* list of pointers to the neighboring edges */
  double *color;      /* vertex color (need this?) */
  kvert *prev,*next;  /* double linked list */
  int  idata;         /* cortex 0  - ics 1 - ocs, when needed  */
};                    /* functions in falcon_off_*.c use this */

struct obj_struct {
  char *fname;             /* filename */
  int nface,nedge,nvert;   /* number of face,edge,vert */
  kface *face;             /* head pointer to the list of triangles */
  kedge *edge;             /* head pointer to the list of edges */
  kvert *vert;             /* head pointer to the list of vertices */
  int color;               /* flag for color usage */
  int spherecoo;           /* spherical coordinates system: 1 = use, 0 = don't use */
  int n_comments;          /* number of comment lines*/
  char ** comment;         /* comments extracted from the file (metadata)*/
};                         /* functions in falcon_off_*.c use this */

typedef struct {
  double bmax,delta;          /* bunding box is isotropic */
  niikpt origin;              /* origin location */
  int depth,                  /* octree depth */
      xdim,ydim,zdim,           /* bb is 3D */
      area,nvox;
  kface ***data;              /* data contains list of list of faces */
  int *ndata;                 /* size of data */
} bbox;                       /* functions in falcon_off_*.c use this */

typedef struct {
  int ndata;                 /* size of data */
  int *idata;                /* list of # data */
  kface ***data;             /* data contains list of faces */
  double dist;               /* min distance */
} bboxlut;                   /* functions in falcon_off_box_lut.c use this */


typedef struct {
  double bmax,delta;          /* bunding box is isotropic */
  niiksph origin;             /* origin location */
  int depth;                  /* tree depth */
  int psi_dim,the_dim,        /* bb is 3D */
      area,nvox;
  kface ***data;              /* data contains list of list of faces */
  int *ndata;                 /* size of data */
} bbox_sph;                   /*  */


typedef struct off_curvature_tmp {
  double *cornerareas[3];
  double *pointareas;
  niikpt *pdir1;
  niikpt *pdir2;
  double *curv1;
  double *curv2;
  double *curv12;
} off_curvature_t;


/******************************************************************
 *
 * falcon_off_bbox_lut.c
 *
 ******************************************************************/
bboxlut *off_bboxlut_calloc();
bboxlut *off_bboxlut_free(bboxlut *bb);
int off_bboxlut_update_objs(bboxlut *bb,kobj **objs,int nobjs);



/******************************************************************
 *
 * falcon_off.c
 *
 ******************************************************************/

kobj *off_obj_init();
kobj *off_kobj_copy(kobj *src);

kvert *off_vert_init();
kvert *off_vert_init_with_v(double x,double y,double z);
kface *off_face_init();
kedge *off_edge_init();
void off_kface_free(kface *f);
void off_kedge_free(kedge *e);
void off_kvert_free(kvert *v);
kobj *off_kobj_free(kobj *obj);

int off_kvert_test_link(kvert *vert);
int off_kface_test_link(kface *face);
int off_kedge_test_link(kedge *edge);

kvert **off_vert_list(int num);
kface **off_face_list(int num);
kedge **off_edge_list(int num);
kvert ***off_vert_matrix(int num1, int num2);

void off_kvert_remove(kvert *v,kobj *obj);
void off_kface_remove(kface *f,kobj *obj);
void off_kedge_remove(kedge *e,kobj *obj);

int off_kobj_write_off(const char * fname,kobj *obj,int write_vertex_normal);
kobj *off_kobj_read_off(const char * fname);
int off_kobj_write_offe(const char *fname,kobj *obj);
int off_kobj_read_offe(const char *fname,kobj *obj);
int off_kobj_show_edge_info(kobj *obj);
int off_kobj_read_off_make_edge(kobj *obj);

int off_kobj_update_vertex_nei(kvert *vert,kedge *edge);
void off_kobj_update_num(kobj *obj);
void off_kobj_update_nvert(kobj *obj);
void off_kobj_update_nface(kobj *obj);
void off_kobj_update_nedge(kobj *obj);

void off_kobj_update_all_index(kobj *obj);
void off_kobj_update_vert_index(kobj *obj);
void off_kobj_update_face_index(kobj *obj);
void off_kobj_update_edge_index(kobj *obj);

int off_kobj_gray_kface(kobj *obj);
int off_kobj_add_color(kobj *obj);
int off_kobj_add_white_color(kobj *obj);
int off_kobj_remove_color(kobj *obj);
int off_kobj_add_green_color(kobj *obj);
int off_kobj_add_blue_color(kobj *obj);
int off_kobj_add_red_color(kobj *obj);
int off_kobj_add_one_color(kobj *obj,double R,double G, double B);
int off_kobj_apply_color_map(kobj *obj,double *var, double dmin,double dmax,int colormap);


kobj *off_create_icosahedron();
/*kobj *off_create_plane(const niikpt *v,double border);HACK right now, look in falcon_midsag.c*/
kobj *off_create_rect(const niikpt *v);

int off_subdiv_kobj(kobj *obj);

void off_swap_kvert(kvert *e1,kvert *e2);
void off_swap_kface(kface *f1,kface *f2);
void off_swap_kedge(kedge *e1,kedge *e2);

int off_match_face_vert_edge(kface *f);

void off_display_vert_info(kvert *v,int type);
void off_display_edge_info(kedge *e,int type);
void off_display_face_info(kface *f,int type);
int off_display_kobj_edge_stats(kobj *obj);

niikpt off_avg_point_on_edge(kedge *e);
void off_update_kface_pminmax(kface *f);
void off_update_kvert_pminmax(kvert *v);
void off_update_kobj_kface_pminmax(kobj *obj);
niikpt off_calc_kobj_pmin(kobj *obj);
niikpt off_calc_kobj_pmax(kobj *obj);

void off_calc_facelist_bounds(kface **facelist,int num,niikpt *pmin,niikpt *pmax);
void off_calc_facelist_bounds_sph(kface **facelist,int num,niiksph *pmin,niiksph *pmax);


void off_update_kface_normal(kface *f);
void off_update_kvert_normal(kvert *v);
int off_update_kobj_face_normal(kobj *obj);
int off_update_kobj_vert_normal(kobj *obj);
int off_smooth_kobj_vert_normal(kobj *obj);

double off_get_kobj_mean_edge_length(kobj *obj);
int off_kobj_display_valency_stats(kobj *obj);

double off_get_kobj_area(kobj *obj);
double off_get_kobj_mean_tri_area(kobj *obj);
niikpt off_get_kobj_global_min_coord(kobj *obj);
niikpt off_get_kobj_global_max_coord(kobj *obj);

int off_get_local_kvert(kvert **vlist,int num);
int off_get_local_kvert_check(kvert **vlist,int vnum,int num);
kvert **off_get_local_kvert_list(kvert *vert,int *num,int local_num);
int off_get_local_kface(kface **flist,int num);
niikpt niikpt_kvert_local_average(kvert *v,int n);
niikpt niikpt_kvert_local_average2(kvert *v);
niikpt niikpt_kvert_local_average_with_tangential_relaxation(kvert *v,int use_tangential_relaxation);
niikpt niikpt_kvert_local_average_with_tangential_relaxation2(kvert *v,double lambda);

niikpt niikpt_kvert_simple_local_avg(kvert *v);


int off_surface_smooth_using_vert(kobj *obj,double *var,int vnei,double wself);
int off_surface_median_smooth_using_vert(kobj *obj,double *var,int vnei);
int off_surface_trimmed_average_smooth_using_vert(kobj *obj,double *var,int vnei,double trim);
int off_surface_field_smooth_using_vert(kobj *obj, niikpt *var, double sigma);
int off_surface_field_smooth_using_vert_thresholded(kobj *obj, niikpt *var, double sigma,double threshold);
int off_surface_gauss_smooth_using_vert(kobj *obj, double *var, double sigma,int iter);



kobj *off_make_sphere_from_icosahedron(double elen,double radius, niikpt ctr);

int niik_affine_transform_off(kobj *obj,niikmat *afmat);

double off_get_kobj_volume(kobj *obj);

/* curvature maps (falcon_off.c) */
int niik_off_curvature_vert(kvert *v,int nei,double *out);
int niik_off_curvature_vert_sph(kvert *v,int nei,double *out);
int niik_off_curvature_map_update(kobj *obj,niikvec *cm,int nei);

int niik_off_apply_surface_smoothing(kobj *obj,int nei,double delta);

kvert *niik_off_find_closest_vertex(kobj *obj,niikpt p);

/* spherical coordinates */
int off_kobj_read_off_sph (const char *fname,kobj *obj);
int off_kobj_write_off_sph(const char *fname,kobj *obj);
int off_avg_smooth_spherical_coordinates(kobj *obj);

/*metadata functions*/
int off_kobj_add_comment(kobj *obj,const char *comment);


/******************************************************************
 *
 * falcon_off_cuberille.c
 *
 ******************************************************************/

kobj *off_cuberille_kobj(nifti_image *img,int flag_all_obj);
kvert *off_cuberille_kobj_search_kvert(kobj *obj,kvert *vtail,int x,int y, int z);
int niik_image_correct_for_cuberille(nifti_image *img,int flag);


/******************************************************************
 *
 * falcon_off_remesh.c
 *
 ******************************************************************/

int off_remesh_kobj(kobj *obj,double elen,int maxiter,int wmark);
int off_remesh_kobj_ex(kobj *obj,double emin,double emax, int maxiter,int wmark,int smooth);
int off_relax_kobj(kobj *obj,double elen,int maxiter,int wmark);

int off_remesh_kedge_collapse_clean_up(kobj *obj,kvert **vrm,kface **frm,kedge **erm);
int off_remesh_kedge_collapse2(kobj *obj,kedge *edge);
int off_remesh_kobj_collapse_correction(kobj *obj);
int off_remesh_kedge_collapse_correction(kobj *obj,kedge *edge);
int off_remesh_kedge_collapse(kobj *obj,kedge *edge,kvert **remove_vert,kface **remove_face,kedge **remove_edge);
kvert *off_remesh_kedge_split(kobj *obj,kedge *edge);
int off_remesh_kedge_flip(kobj *obj,kedge *edge);
int off_remesh_kobj_count_global_valence6(kobj *obj);
int off_remesh_kobj_count_global_valence57(kobj *obj);
int off_remesh_kobj_count_global_valence_variance(kobj *obj);
int off_remesh_kobj_valence_minimization_kedge(kobj *obj,kedge *edge);
int off_remesh_kobj_valence_minimization(kobj *obj,int wmark);
int off_remesh_kvert_relocate(kvert *vert);
int off_remesh_kvert_relocate_lambda(kvert *vert,double lambda);
int off_kobj_test(kobj *obj);
int off_kobj_test_common_neighbor(kedge *edge);
int off_remesh_kedge_local_info(kedge *edge,kvert **v,kface **f,kedge **e);

int off_remesh_dual_kobj(kobj *objs[],off_curvature_t crv[], int maxiter);




/******************************************************************
 *
 * falcon_off_bbox.c
 *
 ******************************************************************/

bbox *off_bbox_calloc();
bbox *off_bbox_free(bbox *bb);
bbox *off_bbox_init(int depth,double bmax);
int off_create_bbox_from_kface_list(bbox *bb,kface **facelist,int num);
int off_create_bbox_from_multiple_kobj(bbox *bb,kobj **obj,int num_obj);
int off_create_bbox_from_kobj(bbox *bb,kobj *obj);
int off_create_bbox_octree(bbox *bb,kface **facelist,int num,int cdepth,niikpt pmin,niikpt pmax);
int off_check_tri_tri_intersect(kface *f1,kface *f2);
int off_check_tri_tri_intersect_with_isectline(niikpt v1,niikpt v2,niikpt v3,
    niikpt u1,niikpt u2,niikpt u3,
    int *coplanar,
    niikpt *isect1, niikpt *isect2);
int off_check_self_intersection_kvert(bbox *bb,kvert *v);
int off_check_self_intersection_kface(bbox *bb,kface *f);
int off_count_self_intersection(bbox *bb,kobj *obj);
int off_check_self_intersection_from_bbox(bbox *bb);
int off_count_self_intersection_add_color(bbox *bb,kobj *obj,int flag_add);
int off_correct_self_intersection(bbox *bb,kobj *obj);
int off_correct_self_intersection_with_elen(bbox *bb,kobj *obj,double elen);
void off_display_self_intersection_kface(kobj *obj);
int off_check_bbox_self_intersection(bbox *bb);
void off_display_bbox_info(bbox *bb);

/******************************************************************
 *
 * falcon_off_bbox_sph.c
 *
 ******************************************************************/
bbox_sph *off_bbox_sph_calloc();
bbox_sph *off_bbox_sph_free(bbox_sph *bb);
bbox_sph *off_bbox_sph_init(int depth,double bmax);
int off_create_bbox_sph_from_kface_list(bbox_sph *bb,
                                        kface **facelist,
                                        int num);
int off_create_bbox_sph_from_multiple_kobj(bbox_sph *bb,kobj **obj,int num_obj);
int off_create_bbox_sph_from_kobj(bbox_sph *bb,kobj *obj);
void off_display_bbox_sph_info(bbox_sph *bb);


/******************************************************************
 *
 * falcon_off_obj2img.c
 *
 ******************************************************************/

int off_obj2img(nifti_image *img,kobj *obj,double pval);
int off_obj2img_use_qform(nifti_image *img,kobj *obj,double pval);
int off_obj2img_color(nifti_image *img,kobj *obj,double pval);

/******************************************************************
 *
 * falcon_off_fsasc2off.c
 *
 ******************************************************************/

int niik_off_fsasc2off(char *iname,char *oname);


/******************************************************************
 *
 * falcon_off_ocv.c
 *
 ******************************************************************/

int niik_off_outer_contour_object(nifti_image *maskimg,kobj *obj,double elen);


/******************************************************************
 *
 * falcon_off_curvature.c
 *
 ******************************************************************/
int niik_off_pointareas(kobj *obj,double *pointareas,double *cornerareas[3]);

int niik_off_curvature_map_rusinkiewicz(kobj *obj,
                                        double *pointareas,
                                        double *cornerareas[3],
                                        niikpt *pdir1, niikpt *pdir2,
                                        double *curv1, double *curv2,double *curv12);


int init_curvature(off_curvature_t *crv,kobj *obj);

int re_init_curvature(off_curvature_t *crv, kobj *obj);

int free_curvature(off_curvature_t *crv);

int update_curvature(off_curvature_t *crv,kobj *obj);

double mean_curvature(off_curvature_t *crv,kobj *obj);

/******************************************************************
 *
 * falcon_ply.c
 *
 ******************************************************************/
int niik_off_write_ply(const char *fname, kobj *obj, const double *meas, int output_sph, int output_normal, int output_edge);
int off_kobj_write_ply(const char *fname, kobj *obj, int write_vertex_normal);
int off_kobj_write_ply_ex(const char *fname,kobj *obj,
       int output_normal,int output_sph,int output_edge, int output_color,
       int n_meas,const char **meas_name, const double **meas);

kobj *off_kobj_read_ply(const char * fname);
kobj *off_kobj_read_ply_ex(const char * fname, int *n_meas, char ***meas_name, double ***meas);

/*compatibility functions*/
kobj *off_kobj_read_offply(const char * fname);
kobj *off_kobj_read_offply_ex(const char * fname, int *n_meas, char ***meas_name, double ***meas);

int off_kobj_write_offply(const char * fname,kobj *obj, int write_vertex_normal);





#endif /*__FALCON_SURFACES_H__*/
/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/