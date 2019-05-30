/* Filename:     nifti1_kunio_off_ocv.c
 * Description:  outer contour object creation
 * Author:       Kunio Nakamura
 * Date:         October 28, 2012
 */

#include "falcon.h"
#include "falcon_cortex.h"

static kvert *g_niikcortex_initics_check_v=NULL;


int niik_off_outer_contour_object(nifti_image *maskimg,kobj *obj,double elen);
int niik_off_shrink_wrap(nifti_image *maskimg,kobj *obj,double len,double thresh,int maxiter,double stepsize1,double stepsize2,
                         int flag_bbox,int flag_remesh,int flag_xdetail);

/************************************************************
 *
 * outer contour object
 *
 ************************************************************/

int niik_off_outer_contour_object(nifti_image *maskimg,kobj *obj,double elen) {

  char
  fcname[512]="niik_off_outer_contour_object";
  double stepsize1=0.3, stepsize2=3.0;
  int maxiter=10;

  niik_fc_display(fcname,1);
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }

  fprintf(stdout,"[%s] initial remesh %9.4f\n",fcname,elen);
  if(!off_remesh_kobj(obj,elen,15,0)) {
    fprintf(stderr,"[%s] ERROR: off_remesh_kobj\n",fcname);
    return 0;
  }

  if(!niik_off_shrink_wrap(maskimg,obj,elen,0.5,maxiter,stepsize1,stepsize2,1,1,1)) {
    fprintf(stderr,"[%s] ERROR: niik_off_shrink_wrap\n",fcname);
    return 0;
  }

  niik_fc_display(fcname,0);
  return 1;
} /* niik_off_outer_contour_object */


int niik_off_shrink_wrap(nifti_image *maskimg,kobj *obj,double len,double thresh,int maxiter,double stepsize1,double stepsize2,
                         int flag_bbox,int flag_remesh,int flag_xdetail)
/*
 * shrink-wrap function
 * -based on the one from nifti1_kunio_off_shrinkwrap.c
 *
 */

{

  bbox *bb=NULL;
  kvert *v;
  niikpt origpt;
  double
  tol=0.01,
  stepsum,
  err;
  int
  iter1,
  iter2,
  midx,
  n,
  ct;
  char fcname[64]="niik_off_shrink_wrap";
  int verbose=0;

  /* CHECK INPUTS */
  if(verbose>=1) niik_fc_display(fcname,1);
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(obj==NULL) {
    fprintf(stderr,"[%s] ERROR: obj is null\n",fcname);
    return 0;
  }

  /* BOUNDING BOX */
  if(flag_bbox) {
    if(verbose>=1) fprintf(stdout,"[%s] bounding box\n",fcname);
    bb=off_bbox_init(6,320);
    if(verbose>=1) fprintf(stdout,"[%s] bbox = %f   %2i\n",fcname,bb->delta,bb->depth);
  } else {
    if(verbose>=1) fprintf(stdout,"[%s] not using bounding box\n",fcname);
  }

  /* check for self intersection */
  if(flag_bbox) {
    if((n=off_count_self_intersection(bb,obj))>0) {
      fprintf(stderr,"[%s] ERROR: self-intersection, %i\n",fcname,n);
      return 0;
    } else {
      if(verbose>=1) fprintf(stdout,"[%s] no self-intersection\n",fcname);
    }
  }

  /* this is the outer loop */
  for(iter2=1; iter2<=maxiter; iter2++) {

    if(verbose>=1) fprintf(stdout,"[%s] ShrinkWrapIter %3i\n",fcname,iter2);
    if(flag_bbox) {
      if(verbose>=2) fprintf(stdout,"[%s] bounding box add triangles\n",fcname);
      off_create_bbox_from_kobj(bb,obj);
    }

    /* update face and vert normal + pmin/pmax */
    if(verbose>=2) fprintf(stdout,"[%s] update normal/pmin/pmax\n",fcname);
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_update_kobj_kface_pminmax(obj);
    off_kobj_update_all_index(obj);

    /* w-mark for remeshing */
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v.w=1.0;
    }

    /* this is the inner loop */
    for(stepsum=0,iter1=1; stepsum<stepsize2; stepsum+=stepsize1,iter1++) {

      for(v=obj->vert,ct=0; v!=NULL; v=v->next) {

        origpt = v->v;

        /* vertex is in the object
         * -move out! */
        if(niik_image_interpolate_3d_linear(maskimg,v->v)>thresh) {
          v->v = niikpt_move_normal(origpt,v->normal,2*stepsize1);
          if(flag_bbox) {
            off_update_kvert_pminmax(v);
            if(off_check_self_intersection_kvert(bb,v)) {
              midx=1;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else {
              midx=2;
              v->v.w = 0;
              ct++;
            }
          } /* using bounding box */
        } /* vertex is in the object */

        /* vertex is outside the object & not using bbox */
        else if(!flag_bbox) {
          if(niik_image_interpolate_3d(maskimg,v->v,NIIK_INTERP_LINEAR)>thresh) {
            midx=3;
            v->v=origpt;
            off_update_kvert_pminmax(v);
          } else {
            midx=5;
            v->v.w = 0;
            ct++;
          }
        } /* not using bounding box */

        else {
          switch(flag_xdetail) {
          default:
          case 0: /* not detailed object intersection check */
            v->v = niikpt_move_normal(origpt,v->normal,-stepsize1);
            off_update_kvert_pminmax(v);
            if(off_check_self_intersection_kvert(bb,v)) {
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else {
              v->v.w=0;
              ct++;
            }
            break;
          case 1: /* detailed object intersection check */
            v->v = niikpt_move_normal(origpt,v->normal,-stepsize1);
            off_update_kvert_pminmax(v);
            if(!niikcortex_initics_shrink_check_deform_interp(maskimg,v,thresh)) {
              midx=1;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else {
              midx=2;
              v->v.w=0;
              ct++;
            }
            break;
          } /* switch */
        } /* using bounding box */

        if(verbose>=3) fprintf(stdout,"[%s] midx = %i\n",fcname,midx);
      } /* each vertex */

      err=(double)ct/obj->nvert;
      if(verbose>=1) fprintf(stdout,"[%s] ShrinkWrapIter %3i  %5.2f %4.2e   %8i / %8i\n",fcname,iter1,stepsum,err,ct,obj->nvert);
      if(err<=tol) {
        if(verbose>=1) fprintf(stdout,"[%s] ShrinkWrapIter %3i  %5.2f %4.2e   %8i / %8i  <tol %8.3g, break loop>\n",fcname,iter1,stepsum,err,ct,obj->nvert,err);
        break;
      }
      if(ct<2) {
        if(verbose>=1) fprintf(stdout,"[%s] ShrinkWrapIter %3i  %5.2f %4.2e   %8i / %8i  <count, break loop>\n",fcname,iter1,stepsum,err,ct,obj->nvert);
        break;
      }
    } /* inner loop */

    /* remesh */
    if(flag_remesh) {
      if(verbose>=1) fprintf(stdout,"[%s] remesh #v %i\n",fcname,obj->nvert);
      if(!off_remesh_kobj(obj,len,8,0)) {
        fprintf(stderr,"[%s] ERROR: off_remesh_kobj\n",fcname);
        return 0;
      }
    }

    if(flag_bbox) {
      if(!off_correct_self_intersection(bb,obj)) {
        fprintf(stderr,"[%s] ERROR: off_correct_self_intersection\n",fcname);
        return 0;
      }
      if((n=off_count_self_intersection(bb,obj))>0) {
        fprintf(stderr,"[%s] ERROR: self-intersection, %i\n",fcname,n);
        return 0;
      } else {
        if(verbose>=1) fprintf(stdout,"[%s] no self-intersection\n",fcname);
      }
    } /* using bounding box */

  } /* outer loop */

  bb=off_bbox_free(bb);
  return 1;
} /* niik_off_shrink_wrap */


int niikcortex_initics_shrink_check_deform_interp(nifti_image *img,kvert *v,double thresh)
/* -returns zero if the vertex is on the object
 * -otherwise return non-zero
 */
{
  char fcname[128]="niikcortex_initics_shrink_check_deform_interp";
  niikpt
  face_pt_list[64],
               edge_pt_list[64];
  double d;
  int
  n,
  xsc=0;
  if(niik_image_interpolate_3d_linear(img,v->v)>thresh) {
    return 0;
  }
  /*
   * check for mid-edge points
   */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"[%s] check mid-edge points\n",fcname);
  }
  #pragma omp parallel for reduction(+:xsc) private(n)
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    edge_pt_list[n] = niikpt_avg(v->neiedge[n]->endpts[0]->v,v->neiedge[n]->endpts[1]->v);
    if(niik_image_interpolate_3d_linear(img,edge_pt_list[n])>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for mid-face points
   */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"[%s] check mid-face points\n",fcname);
  }
  #pragma omp parallel for reduction(+:xsc) private(n)
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    face_pt_list[n] = niikpt_avg3(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v);
    if(niik_image_interpolate_3d_linear(img,face_pt_list[n])>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for mid (mid-face)-vertex points
   */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"[%s] check mid-face points\n",fcname);
  }
  #pragma omp parallel for reduction(+:xsc) private(n)
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(niik_image_interpolate_3d_linear(img,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[0]->v))>thresh) {
      xsc++;
      continue;
    }
    if(niik_image_interpolate_3d_linear(img,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[1]->v))>thresh) {
      xsc++;
      continue;
    }
    if(niik_image_interpolate_3d_linear(img,niikpt_avg(face_pt_list[n],v->neiface[n]->vert[2]->v))>thresh) {
      xsc++;
      continue;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for 1/4 and 3/4-edge points
   */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"[%s] check along 1/4 3/4 edge points\n",fcname);
  }
  #pragma omp parallel for reduction(+:xsc)  private(n)
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(niik_image_interpolate_3d_linear(img,niikpt_avg(edge_pt_list[n],v->neiedge[n]->endpts[0]->v))>thresh) {
      xsc++;
    }
    if(niik_image_interpolate_3d_linear(img,niikpt_avg(edge_pt_list[n],v->neiedge[n]->endpts[1]->v))>thresh) {
      xsc++;
    }
  }
  if(xsc) {
    return 0;
  }
  /*
   * check for 3/4-edge points
   */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"[%s] check along 3/4\n",fcname);
  }
  for(n=0; n<v->nei; n++) {
    if(xsc) continue;
    if(v->neiface[n]->vert[0]==v) {
      for(d=0.2; d<0.999; d+=0.2) {
        if(niik_image_interpolate_3d_linear(img,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[1]->v,v->neiface[n]->vert[2]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
    if(v->neiface[n]->vert[1]==v) {
      for(d=0.2; d<0.998; d+=0.2) {
        if(niik_image_interpolate_3d_linear(img,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[2]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
    if(v->neiface[n]->vert[2]==v) {
      for(d=0.2; d<0.998; d+=0.2) {
        if(niik_image_interpolate_3d_linear(img,niikpt_wavg(v->v,niikpt_wavg(v->neiface[n]->vert[0]->v,v->neiface[n]->vert[1]->v,d),0.25))>thresh) {
          xsc++;
        }
      }
    }
  }
  if(xsc) {
    return 0;
  }
  return 1;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/