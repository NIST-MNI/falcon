/* FILENAME:     niik_skull_registration.c
 * DESCRIPTION:  Kunio's nifti1 skull-based registration program
 * AUTHOR:       Kunio Nakamura
 * DATE:         January 29, 2013
 */

#include "falcon.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define MAJOR_VERSION (1)
#define MINOR_VERSION (0)
#define MICRO_VERSION (0)

static char *prog_version[] = {
  "  niik_skull_registration history\n"
  "  0.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -initial version\n"
  "  1.0.0 knakamura@mrs.bic.mcgill.ca\n"
  "  -different version\n"
};

static char *prog_describe[] = {
  "  [niik_skull_registration] description\n"
  "   1. reads input images, brain masks, and initial <xfm>\n"
  "     -xfm is converted to niik internal matrix format, aka niikmat (same as FSL)\n"
  "     -that is, niikmat's world coordinate is [Wx,Wy,Wz] = [Pixel_x,Pixel_y,Pixel_z] * [Vx,Vy,Vz]\n"
  "   2. create a sphere, centered at the centroid of <mov_brain> with a radius of\n"
  "      <icosahedron_radius> and edge length <icosahedron_elen>\n"
  "   3. shrink-wrap the sphere to fit <mov_brain>\n"
  "     -shrink-wrap has 2 passes: first with coarse <shrink_elen1> followed by\n"
  "      smaller <shrink_elen1> with smaller iteration <shrink_iter2>\n"
  "   4. remesh for smooth well-distributed surface model <remesh_iter>\n"
  "   5. find the minimal intensity value outside the brain surface\n"
  "   6. median-filter the distance to the location of minimal intensity on the surface space\n"
  "   7. transform the shrink-wrap/remeshed surface (from 2) using the initial <init.xfm>\n"
  "   8. repeat steps 3&4 with the transformed surface\n"
  "   9. calculates the updated transformation xfm using singular value decomposition\n"
  "  10. convert calculated matrix to xfm format and write it\n"
};

static char *prog_help[] = {
  "  niik_skull_registration:\n"
  "\n"
  "  optional usage:\n"
  "  -u -help --help                   : show this usage\n"
  "  --version                         : show version info\n"
  "  -icosahedron_radius=<radius>      : initial sphere radius (to shrink-wrap) based on icosahedron [default=100]\n"
  "  -icosahedron_elen=<elen>          : initial sphere edge length [default=3]\n"
  "  -search_distance_step=<step>      : search step distance [default=0.99]\n"
  "  -search_distance_neg=<neg>        : search distance in negative direction [default=0.5]\n"
  "  -search_distance_pos=<pos>        : search distance in positive direction [default=5.0]\n"
  "  -shrink_iter1=<iter>              : number of iteration for first shrink-wrap [default=20]\n"
  "  -shrink_iter2=<ite2>              : number of iteration for second shrink-wrap [default=10]\n"
  "  -shrink_step1=<step1>             : shrink-wrap step size for first pass  [default=10]\n"
  "  -shrink_step2=<step2>             : shrink-wrap step size for second pass [default=5]\n"
  "  -shrink_elen1=<elen1>             : shrink-wrap edge length step size for first pass [default=5]\n"
  "  -shrink_elen2=<elen2>             : shrink-wrap edge length step size for sercond pass [default=3.5]\n"
  "                                    : this is the final edge length\n"
  "  -remesh_iter=<iter>               : remesh iteration [default=15]\n"
  "  -median_filter_size=<size>        : surface median filter on skull distance map [default=3]\n"
  "                                    : this will be related to shrink_elen2 (such as edge length and neighborhood size)\n"
};

void usage() {
  fprintf(stdout,"niik_skull_registration\n");
  fprintf(stdout,"  usage: [options] <mov> <mov_brain> <target> <target_brain> <init.xfm> <out.xfm>\n\n");
  /*(  fprintf(stdout,"\n  use nifti format for images and brain masks\n");*/
}


int niikmat_linalg_SV_solve(niikmat *mat,double *v) {
  gsl_matrix
  *A, *V;
  gsl_vector *S,*work,*x,*b;
  int m,n,verbose=0;
  char fcname[32]="niikmat_linalg_SV_solve";
  NIIK_RET0(((A = gsl_matrix_alloc(mat->row,mat->col))==NULL),fcname,"gsl_matrix_alloc for A");
  for(m=0; m<mat->row; m++)
    for(n=0; n<mat->col; n++)
      gsl_matrix_set (A,m,n,mat->m[m][n]);
  NIIK_RET0(((V = gsl_matrix_alloc(mat->col,mat->col))==NULL),fcname,"gsl_matrix_alloc for V");
  NIIK_RET0(((S = gsl_vector_alloc(mat->col))==NULL),fcname,"gsl_matrix_alloc for S");
  NIIK_RET0(((work = gsl_vector_alloc(mat->col))==NULL),fcname,"gsl_vector_alloc for work");
  NIIK_RET0((gsl_linalg_SV_decomp(A,V,S,work)),fcname,"gsl_linalg_SV_decomp");
  NIIK_RET0(((x = gsl_vector_alloc(mat->col))==NULL),fcname,"gsl_vector_alloc for x");
  NIIK_RET0(((b = gsl_vector_alloc(mat->row))==NULL),fcname,"gsl_vector_alloc for b");
  for(n=0; n<mat->row; n++)
    gsl_vector_set (b,n,v[n]);
  NIIK_RET0((gsl_linalg_SV_solve(A,V,S,b,x)),fcname,"gsl_linalg_SV_solve");
  if(verbose>=1) gsl_vector_fprintf(stdout,x,"%12.6f");
  for(n=0; n<mat->col; n++)
    v[n] = gsl_vector_get (x,n);
  gsl_vector_free(work);
  gsl_vector_free(x);
  gsl_vector_free(b);
  gsl_vector_free(S);
  gsl_matrix_free(A);
  gsl_matrix_free(V);
  return 1;
}


niikmat *niikmat_affine_matrix_from_matching_landmarks(niikpt *x,niikpt *y,int num) {
  int verbose=0;
  char fcname[64]="niikmat_affine_matrix_from_matching_landmarks";
  niikmat *M=NULL,*out=NULL;
  int m,n;
  double *b;
  if(verbose>=1) niik_fc_display(fcname,1);
  NIIK_RET0(((M=niikmat_init(num,4))==NULL),
            fcname,"niikmat_init for M");
  NIIK_RET0(((out=niikmat_identity(4,4))==NULL),
            fcname,"niikmat_identity for out");
  b=(double *)calloc(num,sizeof(double));
  for(m=0; m<num; m++) {
    b[m]=y[m].x;
    M->m[m][0] = x[m].x;
    M->m[m][1] = x[m].y;
    M->m[m][2] = x[m].z;
    M->m[m][3] = 1;
  }
  /* X-direction */
  if(verbose>=1) fprintf(stdout,"[%s] x-direction\n",fcname);
  for(m=0; m<num; m++) {
    b[m]=y[m].x;
  }
  NIIK_RET0((!niikmat_linalg_SV_solve(M,b)),fcname,"niik_linalg_SV_solve X");
  for(n=0; n<4; n++) {
    out->m[0][n]=b[n];
  }
  if(verbose>=2) niikmat_display(out);
  /* Y-direction */
  if(verbose>=1) fprintf(stdout,"[%s] y-direction\n",fcname);
  for(m=0; m<num; m++) {
    b[m]=y[m].y;
  }
  NIIK_RET0((!niikmat_linalg_SV_solve(M,b)),fcname,"niik_linalg_SV_solve Y");
  for(n=0; n<4; n++) {
    out->m[1][n]=b[n];
  }
  if(verbose>=2) niikmat_display(out);
  /* Z-direction */
  if(verbose>=1) fprintf(stdout,"[%s] z-direction\n",fcname);
  for(m=0; m<num; m++) {
    b[m]=y[m].z;
  }
  NIIK_RET0((!niikmat_linalg_SV_solve(M,b)),fcname,"niik_linalg_SV_solve Z");
  for(n=0; n<4; n++) {
    out->m[2][n]=b[n];
  }
  if(verbose>=1) niikmat_display(out);
  if(verbose>=1) niik_fc_display(fcname,1);
  return out;
}


int niik_image_skull_registration_fit_skull(nifti_image *img,kobj *obj,double *warp_distance_map,
    double search_distance_neg,
    double search_distance_pos,
    double search_distance_step,
    int median_filter_size,
    int use_gradient) {
  char fcname[64]="niik_image_skull_registration_fit_skull";
  int n,m,num,verbose=1;
  kvert *v=NULL;
  double d,voxval,maxval,max_distance,minval,min_distance,*vlist=NULL;
  if(verbose>=1) niik_fc_display(fcname,1);
  if(verbose>=1) fprintf(stdout,"[%s] search for skull (point of min intensity)\n",fcname);
  NIIK_RET0((!off_update_kobj_face_normal(obj)),fcname,"off_update_kobj_face_normal");
  NIIK_RET0((!off_update_kobj_vert_normal(obj)),fcname,"off_update_kobj_vert_normal");
  NIIK_RET0((!off_smooth_kobj_vert_normal(obj)),fcname,"off_smooth_kobj_vert_normal");
  NIIK_RET0((!off_smooth_kobj_vert_normal(obj)),fcname,"off_smooth_kobj_vert_normal");
  NIIK_RET0((!off_smooth_kobj_vert_normal(obj)),fcname,"off_smooth_kobj_vert_normal");
  if(use_gradient) {
    fprintf(stdout,"[%s] using gradient\n",fcname);
    for(d=search_distance_neg,num=0; d<=search_distance_pos+search_distance_step*0.1; d+=search_distance_step) {
      num++;
    }
    vlist=(double *)calloc(num,sizeof(double));
  }
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    for(d=search_distance_neg,max_distance=min_distance=0,minval=1e11,maxval=-1e11,m=0; d<=search_distance_pos+search_distance_step*0.1; d+=search_distance_step,m++) {
      if(use_gradient) {
        NIIK_RET0((!niik_image_interpolate_3d_xyz_update(img,
                   niikpt_move_normal(v->v,v->normal,d),
                   NIIK_INTERP_LINEAR,
                   &vlist[m])),
                  fcname,"niik_image_interpolate_3d_xyz_update");
        if(m>0) {
          if((vlist[m]-vlist[m-1])>maxval) {
            maxval=(vlist[m]-vlist[m-1]);
            max_distance=d;
          }
        }
      } else {
        NIIK_RET0((!niik_image_interpolate_3d_xyz_update(img,
                   niikpt_move_normal(v->v,v->normal,d),
                   NIIK_INTERP_NN,
                   &voxval)),
                  fcname,"niik_image_interpolate_3d_xyz_update");
        if(minval>voxval) {
          minval = voxval;
          min_distance = d;
        }
      }
    } /* search region along the normal */
    if(use_gradient)
      warp_distance_map[n] = max_distance;
    else
      warp_distance_map[n] = min_distance;
  } /* each vertex */
  NIIK_RET0((!off_surface_median_smooth_using_vert(obj,warp_distance_map,median_filter_size)),
            fcname,"off_surface_median_smooth_using_vert");
  if(vlist!=NULL) {
    free(vlist);
    vlist=NULL;
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}


int niik_image_skull_registration(nifti_image *img1,nifti_image *bm1,
                                  nifti_image *img2,nifti_image *bm2,
                                  niikmat *mat,
                                  double icosahedron_elen,
                                  double icosahedron_radius,
                                  double shrink_elen1,double shrink_elen2,
                                  double shrink_step1,double shrink_step2,
                                  int shrink_iter1,int shrink_iter2,
                                  int remesh_iter,
                                  int median_filter_size,
                                  double search_distance_neg,
                                  double search_distance_pos,
                                  double search_distance_step)
/* img1 and bm1 are moving
 * img2 and bm2 are target
 * bm is brain mask
 * mat is the initial matrix (niikmat format)
 * img1 and img2 sto_xyz must be correct
 */
{
  char fcname[64]="niik_image_skull_registration";
  kobj *obj=NULL;
  niikpt
  **vlist=NULL,
    ctr;
  double
  sumerr=0,
  *warp_distance_map=NULL;
  int
  n,
  verbose=2;
  kvert
  *v=NULL;
  niikmat *tmpmat=NULL;

  if(verbose>=1) {
    niik_fc_display(fcname,1);
    fprintf(stdout,"[%s] parameters: \n",fcname);
    fprintf(stdout,"  icosahedron_radius                    %.1f:\n",icosahedron_radius);
    fprintf(stdout,"  icosahedron_elen                      %.2f:\n",icosahedron_elen);
    fprintf(stdout,"  shrink 1:\n");
    fprintf(stdout,"    shrink_elen1                        %.2f\n",shrink_elen1);
    fprintf(stdout,"    shrink_step1                        %.2f\n",shrink_step1);
    fprintf(stdout,"    shrink_iter1                        %i\n",shrink_iter1);
    fprintf(stdout,"  shrink 2:\n");
    fprintf(stdout,"    shrink_elen2                        %.2f\n",shrink_elen2);
    fprintf(stdout,"    shrink_step2                        %.2f\n",shrink_step2);
    fprintf(stdout,"    shrink_iter2                        %i\n",shrink_iter2);
    fprintf(stdout,"  remesh_iter                            %i\n",remesh_iter);
    fprintf(stdout,"  median_filter_size                     %i\n",median_filter_size);
    fprintf(stdout,"  search_distance_neg                    %.2f\n",search_distance_neg);
    fprintf(stdout,"  search_distance_pos                    %.2f\n",search_distance_pos);
    fprintf(stdout,"  search_distance_step                   %.2f\n",search_distance_step);
  }

  NIIK_RET0( (img1==NULL),fcname,"img1 is null");
  NIIK_RET0( (img2==NULL),fcname,"img2 is null");
  NIIK_RET0((!img1->sform_code),fcname,"img1 sform is missing");
  NIIK_RET0((!img2->sform_code),fcname,"img2 sform is missing");

  ctr=niikpt_image_get_centroid(bm1,NULL);
  if(verbose>=1) fprintf(stdout,"[%s] centroid %8.3f %8.3f %8.3f\n",fcname,ctr.x,ctr.y,ctr.z);
  NIIK_RET0(((obj=off_make_sphere_from_icosahedron(icosahedron_elen,icosahedron_radius,ctr))==NULL),fcname,"off_make_sphere_from_icosahedron");
  NIIK_RET0((!off_kobj_remove_color(obj)),fcname,"off_kobj_remove_color");

  if(verbose>=2)
    NIIK_RET0((!off_kobj_write_offply("test1.off",obj,0)),
              fcname,
              "off_kobj_write_off");

  if(verbose>=1) fprintf(stdout,"[%s] do shrink-wrap 1\n",fcname);
  NIIK_RET0((!off_shrinkwrap_kobj_bbox_remesh(bm1,obj,shrink_elen1,shrink_iter1,shrink_step1,1,1,0)),
            fcname,
            "off_shrinkwrap_kobj_bbox_remesh");
  if(verbose>=1) fprintf(stdout,"[%s] do shrink-wrap 2\n",fcname);
  NIIK_RET0((!off_shrinkwrap_kobj_bbox_remesh(bm1,obj,shrink_elen2,shrink_iter2,shrink_step2,1,1,0)),
            fcname,
            "off_shrinkwrap_kobj_bbox_remesh");
  NIIK_RET0((!off_remesh_kobj(obj,shrink_elen2,remesh_iter,0)),fcname,"off_remesh_kobj");
  if(verbose>=2) fprintf(stdout,"[%s] vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  if(verbose>=1) fprintf(stdout,"[%s] create vertex lists\n",fcname);
  vlist=(niikpt **)calloc(2,sizeof(niikpt *));
  for(n=0; n<2; n++) vlist[n]=(niikpt *)calloc(obj->nvert,sizeof(niikpt));
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++)
    vlist[0][n]=v->v;
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++)
    vlist[1][n]=niikpt_affine_transform(mat,v->v);

  NIIK_RET0((!off_kobj_add_one_color(obj,0.6,0.6,0.2)),fcname,"off_kobj_add_one_color");
  if(verbose>=2)  {
    fprintf(stdout,"[%s] write off test2_img1.off\n",fcname);
    NIIK_RET0((!off_kobj_write_offply("test2_img1.off",obj,0)),
              fcname,
              "off_kobj_write_off");
  }
  if(verbose>=2) fprintf(stdout,"[%s] vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  NIIK_RET0(((warp_distance_map = (double *)calloc(obj->nvert,sizeof(double)))==NULL),fcname,"calloc for warp_distance_map");

  if(verbose>=1) fprintf(stdout,"[%s] skull fitting for img1\n",fcname);
  NIIK_RET0((!niik_image_skull_registration_fit_skull(img1,obj,warp_distance_map,
             search_distance_neg,
             search_distance_pos,
             search_distance_step,
             median_filter_size,1)),
            fcname,"niik_image_skull_registration_fit_skull (img1)");
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    vlist[0][n] = v->v = niikpt_move_normal(v->v,v->normal,warp_distance_map[n]);
  }

  NIIK_RET0((!off_kobj_add_one_color(obj,0.4,0.9,0.9)),fcname,"off_kobj_add_one_color");
  if(verbose>=2)
    NIIK_RET0((!off_kobj_write_offply("test3_img1.off",obj,0)),
              fcname,
              "off_kobj_write_off");

  if(verbose>=1) fprintf(stdout,"[%s] skull fitting for img2\n",fcname);
  if(verbose>=1) fprintf(stdout,"[%s]   do shrink-wrap 2\n",fcname);
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    v->v = vlist[1][n];
  }
  NIIK_RET0((!off_shrinkwrap_kobj_bbox_remesh(bm2,obj,shrink_elen2,1,2.0,1,0,0)),
            fcname,
            "off_shrinkwrap_kobj_bbox_remesh");
  if(verbose>=2) fprintf(stdout,"[%s] vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  NIIK_RET0((!off_kobj_add_one_color(obj,0.6,0.6,0.2)),fcname,"off_kobj_add_one_color");
  if(verbose>=2)  {
    fprintf(stdout,"[%s] write off test2_img2.off\n",fcname);
    NIIK_RET0((!off_kobj_write_offply("test2_img2.off",obj,0)),
              fcname,
              "off_kobj_write_off");
  }

  for(n=0; n<obj->nvert; n++) warp_distance_map[n]=0;
  NIIK_RET0((!niik_image_skull_registration_fit_skull(img2,obj,warp_distance_map,
             search_distance_neg,
             search_distance_pos,
             search_distance_step,
             median_filter_size,1)),
            fcname,"niik_image_skull_registration_fit_skull (img2)");
  for(v=obj->vert,n=0; v!=NULL; v=v->next,n++) {
    vlist[1][n] = v->v = niikpt_move_normal(v->v,v->normal,warp_distance_map[n]);
  }

  NIIK_RET0((!off_kobj_add_one_color(obj,0.4,0.9,0.9)),fcname,"off_kobj_add_one_color");
  if(verbose>=2)
    NIIK_RET0((!off_kobj_write_offply("test3_img2.off",obj,0)),
              fcname,
              "off_kobj_write_off") ;

  /* estimate the new matrix */
  NIIK_RET0(((tmpmat=niikmat_affine_matrix_from_matching_landmarks(vlist[0],vlist[1],obj->nvert))==NULL),
            fcname,"niikmat_affine_matrix_from_matching_landmarks");
  niikmat_display(tmpmat);

  for(n=0,sumerr=0; n<obj->nvert; n++) {
    sumerr += niikpt_distance(niikpt_affine_transform(mat,vlist[0][n]),vlist[1][n]);
  }
  fprintf(stdout,"[%s] error  = %.10f old\n",fcname,sumerr/obj->nvert);
  for(n=0,sumerr=0; n<obj->nvert; n++) {
    sumerr += niikpt_distance(niikpt_affine_transform(tmpmat,vlist[0][n]),vlist[1][n]);
  }
  fprintf(stdout,"[%s] error  = %.10f new\n",fcname,sumerr/obj->nvert);

  NIIK_RET0((!niikmat_copy_update(tmpmat,mat)),fcname,"niikmat_copy_update");

  for(n=0; n<2; n++) free(vlist[n]);
  free(vlist);
  free(warp_distance_map);
  off_kobj_free(obj);
  tmpmat=niikmat_free(tmpmat);

  niik_fc_display(fcname,0);
  return 1;
}

int main(int argc,char *argv[],char *envp[]) {
  nifti_image
  *tmpimg[2],
  **maskimg=NULL,
    **img=NULL;
  int
  n,nc,sc;
  char
  *outcheck=NULL,
   fcname[32]="niik_skull_registration";
  niikmat *imat=NULL;
  FILE *fp=NULL;
  double
  icosahedron_elen=3,
  icosahedron_radius=100,
  shrink_step1=10,
  shrink_step2=5,
  shrink_elen1=10.5,
  shrink_elen2=3.5,
  search_distance_neg=-0.5,
  search_distance_pos=5.0,
  search_distance_step=0.99;
  int
  remesh_iter=15,
  median_filter_size=3,
  shrink_iter1=20,
  shrink_iter2=10;

  if(argc==1) {
    usage();
    exit(0);
  }

  nc=sc=1;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        fprintf(stdout,"%s\n",*prog_version);
        exit(0);
      } else if(!strncmp(argv[nc],"--describe",10)) {
        fprintf(stdout,"%s",*prog_describe);
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      } else if(!strncmp(argv[nc],"-help",5)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-search_distance_step=",22)) {
        search_distance_step = atof(argv[nc] + 21);
      } else if(!strncmp(argv[nc],"-search_distance_neg=",21)) {
        search_distance_neg = atof(argv[nc] + 21);
      } else if(!strncmp(argv[nc],"-search_distance_pos=",21)) {
        search_distance_pos = atof(argv[nc] + 21);
      } else if(!strncmp(argv[nc],"-median_filter_size=",20)) {
        median_filter_size = atoi(argv[nc] + 20);
      } else if(!strncmp(argv[nc],"-icosahedron_radius=",20)) {
        icosahedron_radius = atof(argv[nc] + 20);
      } else if(!strncmp(argv[nc],"-icosahedron_elen=",18)) {
        icosahedron_elen = atof(argv[nc] + 18);
      }

      else if(!strncmp(argv[nc],"-shrink_step1=",14)) {
        shrink_step1 = atof(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-shrink_step2=",14)) {
        shrink_step2 = atof(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-shrink_elen1=",14)) {
        shrink_elen1 = atof(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-shrink_elen2=",14)) {
        shrink_elen2 = atof(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-shrink_iter1=",14)) {
        shrink_iter1 = atoi(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-shrink_iter2=",14)) {
        shrink_iter2 = atoi(argv[nc] + 14);
      } else if(!strncmp(argv[nc],"-remesh_iter=",13)) {
        remesh_iter = atoi(argv[nc] + 13);
      }

      else if(!strncmp(argv[nc],"-outcheck=",10)) {
        outcheck = argv[nc]+10;
      }

      else if(!strncmp(argv[nc],"-u",2)) {
        fprintf(stdout,"%s\n",*prog_help);
        usage();
        exit(0);
      }

      else {
        fprintf(stderr,"[%s] ERROR: unknown option %s\n",fcname,argv[nc]);
        exit(0);
      }
      nc++;
    } else {
      argv[sc++]=argv[nc++];
    }
  } /* reading options (while) */
  argc=sc;

  if(argc>7) {
    fprintf(stderr,"[%s] ERROR: too many arguments\n",fcname);
    exit(0);
  } else if(argc<7) {
    fprintf(stderr,"[%s] ERROR: too few arguments\n",fcname);
    exit(0);
  }

  niik_fc_display(fcname,1);
  NIIK_EXIT(((    img=(nifti_image **)calloc(2,sizeof(nifti_image *)))==NULL),fcname,"calloc for img",9);
  NIIK_EXIT(((maskimg=(nifti_image **)calloc(2,sizeof(nifti_image *)))==NULL),fcname,"calloc for maskimg",9);

  fprintf(stdout,"[%s] reading input image:    %s\n",fcname,argv[1]);
  NIIK_EXIT(((img[0]=niik_image_read(argv[1]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img[0],NIFTI_TYPE_FLOAT32,1)),
            fcname,
            "niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading input mask:     %s\n",fcname,argv[2]);
  NIIK_EXIT(((maskimg[0]=niik_image_read(argv[2]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert(maskimg[0],NIFTI_TYPE_UINT8)),
            fcname,
            "niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading input image:    %s\n",fcname,argv[3]);
  NIIK_EXIT(((img[1]=niik_image_read(argv[3]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert_scl(img[1],NIFTI_TYPE_FLOAT32,1)),
            fcname,
            "niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading input mask:     %s\n",fcname,argv[4]);
  NIIK_EXIT(((maskimg[1]=niik_image_read(argv[4]))==NULL),
            fcname,
            "reading niik_image_read",9);
  NIIK_EXIT((!niik_image_type_convert(maskimg[1],NIFTI_TYPE_UINT8)),
            fcname,
            "niik_image_type_convert",9);

  fprintf(stdout,"[%s] reading init xfm:       %s\n",fcname,argv[5]);
  NIIK_EXIT(((imat = niikmat_read_xfm(argv[5]))==NULL),
            fcname,
            "niikmat_read_xfm",9);

  if(imat!=NULL) {
    fprintf(stdout,"[%s] using initial matrix\n",fcname);
    niikmat_display(imat);
    fprintf(stdout,"[%s] convert niikmat\n",fcname);
    NIIK_EXIT((!niikmat_convert_from_xfm(imat,img[0],img[1])),fcname,"niikmat_convert_from_xfm",9);
    /*niikmat_display(imat); */
    NIIK_RET0(((tmpimg[0]=niik_image_affine_transform_3d(img[0],img[1],imat,NIIK_INTERP_LINEAR))==NULL),fcname,"niik_image_affine_transform_3d");
    NIIK_RET0(((tmpimg[1]=niik_image_copy(img[1]))==NULL),fcname,"niik_image_copy");
    fprintf(stdout,"[%s] write output tmp.nii.gz\n",fcname);
    NIIK_RET0((!niik_image_combine_and_write("tmp.nii.gz",tmpimg,2,'t',140)),fcname,"niik_image_combine_and_write");
    for(n=0; n<2; n++) tmpimg[n]=niik_image_free(tmpimg[n]);
  }

  fprintf(stdout,"[%s] do skull-based registration\n",fcname);
  NIIK_EXIT((!niik_image_skull_registration(img[0],maskimg[0],
             img[1],maskimg[1],
             imat,
             icosahedron_elen,
             icosahedron_radius,
             shrink_elen1,shrink_elen2,
             shrink_step1,shrink_step2,
             shrink_iter1,shrink_iter2,
             remesh_iter,
             median_filter_size,
             search_distance_neg,
             search_distance_pos,
             search_distance_step)),
            fcname,
            "niik_image_skull_registration",9);

  if(outcheck!=NULL) {
    NIIK_RET0(((tmpimg[0]=niik_image_affine_transform_3d(img[0],img[1],imat,NIIK_INTERP_LINEAR))==NULL),fcname,"niik_image_affine_transform_3d");
    NIIK_RET0(((tmpimg[1]=niik_image_copy(img[1]))==NULL),fcname,"niik_image_copy");
    fprintf(stdout,"[%s] write output %s\n",fcname,outcheck);
    NIIK_RET0((!niik_image_combine_and_write(outcheck,tmpimg,2,'t',140)),fcname,"niik_image_combine_and_write");
    for(n=0; n<2; n++) tmpimg[n]=niik_image_free(tmpimg[n]);
  }

  NIIK_EXIT((!niikmat_convert_to_xfm(imat,img[0],img[1])),fcname,"niikmat_convert_to_xfm",9);
  niikmat_display(imat);

  fprintf(stdout,"[%s] writing output xfm:     %s\n",fcname,argv[6]);
  NIIK_RET0(((fp = fopen(argv[6],"w"))==NULL),fcname,"can not open file");
  fprintf(fp,"MNI Transform File\n");
  fprintf(fp,"Transform_Type = Linear;\n");
  fprintf(fp,"Linear_Transform = \n");
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f\n",
          imat->m[0][0],imat->m[0][1],imat->m[0][2],imat->m[0][3]);
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f\n",
          imat->m[1][0],imat->m[1][1],imat->m[1][2],imat->m[1][3]);
  fprintf(fp,"%15.10f %15.10f %15.10f %15.10f ;\n",
          imat->m[2][0],imat->m[2][1],imat->m[2][2],imat->m[2][3]);
  fclose(fp);

  imat = niikmat_free(imat);
  for(n=0; n<2; n++) {
    img[n]=niik_image_free(img[n]);
    maskimg[n]=niik_image_free(maskimg[n]);
  }
  free(img);
  free(maskimg);
  niik_fc_display(fcname,0);
  exit(0);
} /* niik_skull_registration */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/