/* Filename:     nifti1_kunio_off_shrinkwrap.c
* Description:  object functions for shrinkwrap
* Author:       Kunio Nakamura
* Date:         March 1, 2012
*/
#include "falcon.h"
#include "falcon_cortex.h"

int off_shrinkwrap_kobj_simple(nifti_image *img,kobj *obj,double len) {
  if(!off_shrinkwrap_kobj(img,obj,len,30 ,10,0,0.0 )) {
    fprintf(stderr,"ERROR: off_shrinkwrap_kobj \n");
    return 0;
  }
  return 1;
}

int off_shrinkwrap_kobj(nifti_image *img,kobj *obj,double len,int maxiter,double stepsize,int debug, double smooth) {
  if(!off_shrinkwrap_kobj_bbox(img,obj,len,maxiter,stepsize,1,debug,smooth)) {
    fprintf(stderr,"ERROR: off_shrinkwrap_kobj_bbox\n");
    return 0;
  }
  return 1;
}

int off_shrinkwrap_kobj_bbox(nifti_image *img,kobj *obj, double len, int maxiter, double stepsize, int flag_bbox, int debug, double smooth) {
  if(!off_shrinkwrap_kobj_bbox_remesh(img,obj,len,maxiter,stepsize,flag_bbox,1,debug,smooth)) {
    fprintf(stderr,"ERROR: off_shrinkwrap_kobj_bbox\n");
    return 0;
  }
  return 1;
}

int off_shrinkwrap_kobj_bbox_remesh(nifti_image *img,kobj *obj, double len, int maxiter, double stepsize, 
                                    int flag_bbox, int flag_remesh, int debug, double smooth)
/*
* shrink-wrap function
* -may/may not use bounding box
* -outer loop iterates maxiter times and updates bounding box,
* -inner loop iterates up to stepsize (with step = 1/10 of stepsize )
*  and moves inwards while checking the voxel intensity and self-
*  intersections
*
*/
{
  nifti_image *tmpimg = NULL;
  bbox *bb=NULL;
  kvert *v;
  niikpt origpt;

  double tol=0.01,
         stepsum,
         err,
         step;
  
  int     iter1,
          iter2,
          midx,
          n,
          ct;

  char fname[512];
  const char *fcname=__func__;

  int verbose = niik_verbose();
  cortex_tracing_info trace;
  int debug_tracing = falcon_tracing_init(img, &trace);

  if(img == NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",__func__);
    return 0;
  }
  if(obj == NULL) {
    fprintf(stderr,"[%s] ERROR: obj is null\n",__func__);
    return 0;
  }

  if(verbose>1) fprintf(stdout,"[%s] start off_shrinkwrap_kobj_simple\n",__func__);

  if(debug) {
    verbose=2;
    /* make a copy image for uint8 datatype */
    if((tmpimg = niik_image_copy(img))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_copy\n",__func__);
      return 0;
    }
    if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",__func__);
      return 0;
    }
  }

  if(flag_bbox) {
    if(verbose>1) fprintf(stdout,"[%s] bounding box\n", __func__ );
    bb=off_bbox_init(7,320);
    if(verbose>1) fprintf(stdout,"[%s] bbox = delta=%f depth=%2i\n", __func__ ,bb->delta,bb->depth);
  } else {
    if(verbose>1) fprintf(stdout,"[%s] not using bounding box\n", __func__);
  }


  /* determine the step size based on the image voxel size ?*/
  /*step = NIIK_DMIN(img->dx,NIIK_DMIN(img->dy,img->dz)) * 0.8;*/
  step = stepsize / 10.0;
  if(verbose>1) fprintf(stdout,"[%s] step %8.4f\n", __func__,step);

  /* check for self intersection */
  if(flag_bbox) {
    if((n=off_count_self_intersection(bb,obj))>0) {
      fprintf(stderr,"ERROR: [%s] self-intersection, %i\n", __func__ ,n);
      return 0;
    } else {
      if(verbose>=1) fprintf(stdout,"[%s]     no self-intersection\n", __func__);
    }
  }

  /* this is the outer loop */
  for(iter2=1; iter2<=maxiter; iter2++) {

    if(verbose>1) fprintf(stdout,"[%s] ShrinkWrapIter %3i\n", __func__,iter2);
    /* update face and vert normal + pmin/pmax */
    if(verbose>1) fprintf(stdout,"[%s] update normal/pmin/pmax\n", __func__);

    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_update_kobj_kface_pminmax(obj);
    off_kobj_update_all_index(obj);

    if(flag_bbox) {
      if(verbose>1) fprintf(stdout,"[%s] bounding box add triangles\n", __func__);
      off_create_bbox_from_kobj(bb,obj);

      if(debug_tracing)
        falcon_tracing_dump(&trace,iter2,"shrinkwrap",img,bb);
    }

    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v.w=1.0;
    }

    for(stepsum=0,iter1=1; stepsum<stepsize; stepsum+=step,iter1++) {

      for(v=obj->vert, ct=0; v!=NULL; v=v->next) {
        niikpt update_normal=v->normal;
        niikpt update_smooth=niikpt_zero();

        if(smooth>0.0) {
          int nei;
          for(nei=0; nei<v->nei; nei++)
            update_smooth = niikpt_add(update_smooth,v->neivert[nei]->v);

          update_smooth = niikpt_kmul(update_smooth, 1.0/(double)v->nei);
          update_smooth = niikpt_sub(update_smooth, v->v);
        }

        origpt = v->v;

        /* vertex is in the object
        * expand (direction of the normal) + smoothing
        * */
        if(niik_image_interpolate_3d_linear(img, v->v) >= 0.5) {

          v->v = niikpt_move_normal(origpt, update_smooth, 2.0*step*smooth       );
          v->v = niikpt_move_normal(v->v,   update_normal, 2.0*step*(1.0-smooth) );

          if(flag_bbox) {
            off_update_kvert_pminmax(v);

            /*check for self intersect*/
            if(off_check_self_intersection_kvert( bb, v)) {
              midx = 1;
              v->v = origpt;
              off_update_kvert_pminmax(v);
            } else {
              midx = 2;
              v->v.w = 0;
              ct++;
            }
          } /* using bounding box */
        }  /* vertex is in the object */
        else { /* vertex is outside the object, shrink + smooth */

          v->v = niikpt_move_normal(origpt, update_smooth,      step*smooth);
          v->v = niikpt_move_normal(v->v,   update_normal, -1.0*step*(1.0-smooth));

          if(flag_bbox) { /* using bounding box */
            off_update_kvert_pminmax(v);

            /*do not step into object*/
            if(niik_image_interpolate_3d_linear(img, v->v) >= 0.5) {
              midx=3;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            /*check for self intersect*/
            } else if(off_check_self_intersection_kvert(bb,v)) {
              midx=4;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else { /*actually moved successfully*/
              midx=5;
              v->v.w = 0;
              ct++;
            }
          } else { /* not using bounding box */
            if(niik_image_interpolate_3d_linear(img, v->v) >= 0.5) {
              midx=3;
              v->v=origpt;
              off_update_kvert_pminmax(v);
            } else {
              midx=5;
              v->v.w = 0;
              ct++;
            }
          }
          if(verbose>2) fprintf(stdout,"[%s] midx = %i\n", __func__ ,midx);
        }
      } /* each vertex */
      err = (double)ct/obj->nvert;

      if(verbose>0) fprintf(stdout,"[%s] ShrinkWrapIter %3i  %5.2f %4.2e   %8i / %8i\n", __func__, iter1 , stepsum, err , ct, obj->nvert);
      if(err <= tol ) break;
      if(ct<2) break;
    } /* inner loop */


    /* there should not be self-intersection here
    * ** unless the step size and bbox size are bad
    if((n=off_count_self_intersection_add_color(bb,obj,1))>0){
      off_display_self_intersection_kface(obj);
      fprintf(stderr,"ERROR: self-intersection, %i\n",n);
      sprintf(fname,"tmp_iter%i-prob.off",iter2);
      fprintf(stdout,"  writing %s\n",fname);
      if(!off_kobj_write_offply(fname,obj,0)){
        fprintf(stderr,"ERROR: off_kobj_write_off \n");
        exit(0); }
      return 0; }
      else { fprintf(stdout,"        no self-intersection\n"); }
    */

    /* remesh */
    if(flag_remesh) {
      if(!off_remesh_kobj_ex(obj,len/2,len*2,5,0,1)) {
        fprintf(stderr,"[%s] ERROR: off_remesh_kobj \n", __func__);
        return 0;
      }
    }

    if(flag_bbox) {
      if(!off_correct_self_intersection(bb,obj)) {
        fprintf(stderr,"[%s] ERROR: off_correct_self_intersection\n",__func__);
        return 0;
      }
      if((n=off_count_self_intersection(bb,obj))>0) {
        fprintf(stderr,"[%s] ERROR: self-intersection, %i\n",__func__,n);
        return 0;
      } else {
        if(verbose>=1) fprintf(stdout,"[%s] no self-intersection\n",__func__);
      }
    } /* using bounding box */

    if(debug) {
      sprintf(fname,"tmp_iter%03i.ply",iter2);
      fprintf(stdout,"  writing %s\n",fname);
      if(!off_kobj_write_offply(fname,obj,0)) {
        fprintf(stderr,"[%s] ERROR: off_kobj_write_off \n",__func__);
        exit(0);
      }
      off_obj2img(tmpimg,obj,iter2+10);
    }
  } /* larger loop */

  if(debug) {
    fprintf(stdout,"\twriting image:    tmp_shrink_wrap.mnc\n");
    niik_image_write("tmp_shrink_wrap.mnc",tmpimg);
    nifti_image_free( tmpimg );
    tmpimg=NULL;
  }

  if(flag_bbox)
    bb=off_bbox_free(bb);

  falcon_tracing_free(&trace);
  return 1;
} /* off_shrinkwrap_kobj */


/**************************************************************
*
* balloon
* expand surface (move along normals) avoiding self intersections
*
**************************************************************/
int off_kobj_balloon(kobj *obj, int maxiter, double stepsize, int debug) {
  const char *fcname=__func__;
  int iter=0;
  bbox *bb;
  niikpt origpt;
  kvert *v;
  if(obj==NULL) {
    fprintf(stderr,"[%s] ERORR: obj is null\n",fcname);
    return 0;
  }

  if((bb=off_bbox_init(7,320))==NULL) {
    fprintf(stderr,"[%s] ERROR: off_bbox_init\n",fcname);
    return 0;
  }
  off_create_bbox_from_kobj(bb,obj);
  
  for(iter=0; iter<maxiter; iter++) {
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_smooth_kobj_vert_normal(obj);
    off_update_kobj_kface_pminmax(obj);
    off_kobj_update_all_index(obj);

    for(v=obj->vert; v!=NULL; v=v->next) {
      origpt=v->v;
      v->v = niikpt_move_normal(origpt, v->normal, stepsize);
      off_update_kvert_pminmax(v);
      if(off_check_self_intersection_kvert(bb,v)) {
        v->v=origpt;
        off_update_kvert_pminmax(v);
      }
    }
  } /* iter */
  bb=off_bbox_free(bb);
  return 1;
}


/*
kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/