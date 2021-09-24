/* Filename:     nifti1_kunio_cortex_deformation.c
 * Description:  functions for deformation
 * Author:       Kunio Nakamura
 * Date:         July 4, 2013
 */

#include "falcon.h"
#include "falcon_cortex.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif

#define _USE_PRL_DEFORM

static const double g_niikcortex_min_deform_distance = 1e-4;
static const double g_niikcortex_min_deform_distance2 = 1e-3;
extern int g_niik_cortex_debug;
extern int  niikcortex_get_debug();
extern void niikcortex_set_debug(int debug);



/**********************************************************************
 *
 * niikcortex_deformation
 *
 **********************************************************************/


static int g_niikcortex_deform_calc_deformation_vertex_index=-1; /*VF: danger, can go out of array bounds*/

int niikcortex_deform_apply_deformation(bbox *bb, kobj *obj, double step, niikpt *dfmlist, double *actdfm)
/* -applies deformation
 * -bb is the updated bounding box system
 * -obj is the object to deform
 * -step is the deformation size
 * --if individual deformation magnitude is smaller than step, then
 *   applies magnitude of the deformation (whichever is smaller)
 * -dfmlist is the list of deformation
 *    -dfmlist[*].w must contain the magnitude of xyz deformation
 * -actdfm is the actual deformation distance and is added
 */
{
  kvert *v;
  double dist,cstep,dd;
  niikpt old_v;
  int vidx,num;
  const int verbose=0;
  const char *fcname=__func__;

  if(verbose>=1 && g_niikcortex_deform_calc_deformation_vertex_index>=0 ) {
    fprintf(stdout,"[%s] start  vertex %i actual %-8.4f residual %-8.4f\n",fcname,
            g_niikcortex_deform_calc_deformation_vertex_index,
            actdfm[g_niikcortex_deform_calc_deformation_vertex_index],
            dfmlist[g_niikcortex_deform_calc_deformation_vertex_index].w);
  }

  for(v=obj->vert,vidx=num=0; v!=NULL; v=v->next,vidx++) {
    dd=dfmlist[vidx].w;
    if(dd<g_niikcortex_min_deform_distance) {
      continue;
    } else if(dd<step) {
      cstep = dd;
    } else {
      cstep = step;
    }
    old_v = v->v;
    v->v = niikpt_move_normal(v->v, dfmlist[vidx], cstep);
    off_update_kvert_pminmax(v);

    if(bb!=NULL) {
      if(off_check_self_intersection_kvert(bb,v)) { /*VF: 85% of the program spends here*/
        v->v=old_v;
        off_update_kvert_pminmax(v);
        continue;
      }
    }
    num++;
    dist = dfmlist[vidx].w-cstep;
    /*dfmlist[vidx]=niikpt_kmul(dfmlist[vidx],dist/dfmlist[vidx].w);*/
    dfmlist[vidx].w = dist;
    actdfm[vidx] += cstep;
  } /* each vertex */

  if(verbose>=1) {
    if(g_niikcortex_deform_calc_deformation_vertex_index>=0) {
      fprintf(stdout,"[%s] finish vertex %i actual %-8.4f residual %-8.4f\n",fcname,
              g_niikcortex_deform_calc_deformation_vertex_index,
              actdfm[g_niikcortex_deform_calc_deformation_vertex_index],
              dfmlist[g_niikcortex_deform_calc_deformation_vertex_index].w);
    }
    if(verbose>=2)
      fprintf(stdout,"[%s] deformation %8i / %i\n",fcname,num,obj->nvert);
  }
  return 1;
} /* niikcortex_deform_apply_deformation */


int niikcortex_deform_apply_deformation_prl(bbox *bb,kobj *obj,double step,niikpt *dfmlist,double *actdfm,indicator_t *indicator,int vfilter)
/* -applies deformation using OpenMP parallelization
 * -bb is the updated bounding box system
 * -obj is the object to deform
 * -step is the deformation size
 * --if individual deformation magnitude is smaller than step, then
 *   applies magnitude of the deformation (whichever is smaller)
 * -dfmlist is the list of deformation
 *    -dfmlist[*].x,y,z - unit vector
 *    -dfmlist[*].w must contain the magnitude of xyz deformation
 * -actdfm is the actual deformation distance and is added
 */
{
  int num;
  const int verbose=0;
  const int stride=8;
  const char *fcname=__func__;
  int i,c,s;
  if(verbose>=1 && g_niikcortex_deform_calc_deformation_vertex_index>=0 ) {
    fprintf(stdout,"[%s] start  vertex %i actual %-8.4f residual %-8.4f\n",fcname,
            g_niikcortex_deform_calc_deformation_vertex_index,
            actdfm[g_niikcortex_deform_calc_deformation_vertex_index],
            dfmlist[g_niikcortex_deform_calc_deformation_vertex_index].w);
  }

  /*initialize indicator*/
  memset(indicator, 0, obj->nvert*sizeof(indicator_t));
  num=0;
  for(c=0; c<stride; c++) {
    int snum=0;
    
    /*distribute tasks in a predictable fashion*/
    #pragma omp parallel for reduction(+:snum)
    for(i=c; i<bb->zdim; i+=stride) {
      int j;
      for(j=0; j<bb->ydim; j++) {
        int k;
        int index = j*bb->xdim + i*bb->area;
        for(k=0; k<bb->xdim; k++,index++) {
          /*iterate over all voxels in cell*/

          int n;
          for(n=0; n<bb->ndata[index]; n++) {
            kface *cf=bb->data[index][n];
            int l;
            if(cf->vert[0]->idata != vfilter) continue; /*wrong cortex*/

            for(l=0; l<3; l++) {
              double dd,cstep,dist;
              niikpt old_v;
              int vidx;
              int skip;
              kvert *v;
              v=cf->vert[l];

              vidx=v->index-1;
              /*skip over processed vertex*/
              /*#pragma omp atomic read*/
              skip=indicator[vidx];
              if(skip) continue;

              dd=dfmlist[vidx].w;
              if (dd>=g_niikcortex_min_deform_distance2) {
                if(dd<step) {
                  cstep = dd;
                } else {
                  cstep = step;
                }

                old_v = v->v;
                v->v = niikpt_move_normal(v->v, dfmlist[vidx], cstep);
                off_update_kvert_pminmax(v);

                if(!isfinite(v->v.x) || !isfinite(v->v.y) || !isfinite(v->v.z) ||
                  off_check_self_intersection_kvert(bb,v)) { /*VF: 85% of the program spends here*/
                  v->v=old_v;
                  off_update_kvert_pminmax(v);
                } else {
                  snum++;
                  dist = dfmlist[vidx].w-cstep;
                  //dfmlist[vidx]=niikpt_kmul(dfmlist[vidx], dist/dfmlist[vidx].w );
                  dfmlist[vidx].w = dist;
                  actdfm[vidx] += cstep;
                }
              }
              /*#pragma omp atomic write*/
              indicator[vidx]=1;
            }
          }
        }
      }
    }
    num+=snum;
  }

  if(verbose>=1) {
    if(g_niikcortex_deform_calc_deformation_vertex_index>=0) {
      fprintf(stdout,"[%s] finish vertex %i actual %-8.4f residual %-8.4f\n",fcname,
              g_niikcortex_deform_calc_deformation_vertex_index,
              actdfm[g_niikcortex_deform_calc_deformation_vertex_index],
              dfmlist[g_niikcortex_deform_calc_deformation_vertex_index].w);
    }
    if(verbose>=2)
      fprintf(stdout,"[%s] deformation %8i / %i\n",fcname,num,obj->nvert);
  }
  return 1;
} /* niikcortex_deform_apply_deformation */


/********************************************************
 *
 * calculation of deformation for each vertex
 *
 ********************************************************/
niikpt niikcortex_deform_calc_deformation_vertex_surface_smoothness(kvert *v) {
  return niikpt_sub(niikpt_kvert_local_average_with_tangential_relaxation2(v,0.999),v->v);
}

niikpt niikcortex_deform_calc_deformation_vertex_surface_smoothness_simple(kvert *v) {
  return niikpt_sub(niikpt_kvert_simple_local_avg(v),v->v);
}

niikpt niikcortex_deform_calc_deformation_vertex_thickness_smoothness_term(kvert *v,kvert *vo,double *thickness) {
  /**
   * Move point to maintain constant thickness 
   */
  niikpt diff;
  double dist,proj;
  int i;
  double avg_thickness=0.0;

  diff=niikpt_unit(niikpt_sub(v->v, vo->v));
  
  /*local thickness average*/
  for(i=0;i<v->nei;i++)
    avg_thickness+=thickness[v->neivert[i]->index-1];
  avg_thickness/=v->nei;
  /*move in the opposite direction to minimize the difference*/
  return niikpt_kmul(diff, (avg_thickness-thickness[v->index-1])*0.5);
}

/**
*/
niikpt niikcortex_deform_calc_deformation_prior_term(kvert *v, nifti_image **prior, nifti_image ***grad_prior, int cortex_id)
{
    double prior_obj = niik_image_interpolate_3d_linear(     prior[cortex_id*2], v->v);
    double prior_bck = 0.0;
    double project;

    niikpt grad_prior_obj;
    niikpt grad_prior_bck;

    grad_prior_obj.x = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2][1], v->v);
    grad_prior_obj.y = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2][2], v->v);
    grad_prior_obj.z = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2][3], v->v);

    if(prior[cortex_id*2+1]){
        prior_bck=niik_image_interpolate_3d_linear(     prior[cortex_id*2+1], v->v);

        grad_prior_bck.x = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2+1][1], v->v);
        grad_prior_bck.y = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2+1][2], v->v);
        grad_prior_bck.z = niik_image_interpolate_3d_linear(grad_prior[cortex_id*2+1][3], v->v);
    } else {
        prior_bck = 1.0 - prior_obj;

        grad_prior_bck.x = -grad_prior_obj.x;
        grad_prior_bck.y = -grad_prior_obj.y;
        grad_prior_bck.z = -grad_prior_obj.z;
    }

    /*TODO: normalize gradients ?*/
    grad_prior_obj   = niikpt_unit(grad_prior_obj);  /*only interested in the direction*/
    grad_prior_bck   = niikpt_unit(grad_prior_bck);  /*only interested in the direction*/

    /*should balance at 0.5 p value between WM and GM   or GM and CSF*/
    project = ( niikpt_dot(v->normal, grad_prior_obj) * prior_obj +
                niikpt_dot(v->normal, grad_prior_bck) * prior_bck ) / 2.0;

    return niikpt_kmul(v->normal, -1 * project);
}


/* VF: Balooning to stay on the edge of the tissue "class", and in the area with high intensity gradient
*/
niikpt niikcortex_deform_calc_deformation_vertex_image_term(
  kvert *v, nifti_image *img, nifti_image ***grad_img,
  double *gradient_lambda,double mean_ics_value, double range_ics_value,
  double mean_ocs_value,double range_ocs_value, int cortex_id, double *gray, double *grad) {
  niikpt local_grad;
  double project;

  *gray = niik_image_interpolate_3d_linear(     img, v->v);

  /*move in the right direction*/
  local_grad.x = niik_image_interpolate_3d_linear(grad_img[cortex_id][1], v->v);
  local_grad.y = niik_image_interpolate_3d_linear(grad_img[cortex_id][2], v->v);
  local_grad.z = niik_image_interpolate_3d_linear(grad_img[cortex_id][3], v->v);
  local_grad   = niikpt_unit(local_grad); /*only interested in the direction*/

  project      = niikpt_dot(v->normal, local_grad);

  /*we want to move in the direction of the gradient that would bring the *gray closer to the mean_XXX_value*/

  switch(cortex_id) {
  case CORTEX_ICS: /* WHITE MATTER SURFACE */
    *grad = niik_image_interpolate_3d_linear(grad_img[0][0], v->v);

    return niikpt_kmul(v->normal, -1 * project *
                       gradient_lambda[0] / (*grad+gradient_lambda[0]) *
                       NIIK_Heaviside11(*gray-mean_ics_value, range_ics_value));

  case CORTEX_OCS: /* PIAL SURFACE */
    *grad = niik_image_interpolate_3d_linear(grad_img[1][0], v->v);
    return niikpt_kmul(v->normal, -1 * project *
                       gradient_lambda[1] / (*grad+gradient_lambda[1]) *
                       NIIK_Heaviside11(*gray-mean_ocs_value, range_ocs_value));
  }
  return niikpt_zero();
} /* niikcortex_deform_calc_deformation_vertex_image_term */


niikpt niikcortex_deform_calc_deformation_vertex_gradient_term(kvert *v, nifti_image *img, nifti_image **grad_img, niikpt *grad, int cortex_id) 
{
  grad->x = niik_image_interpolate_3d_linear(grad_img[1], v->v);
  grad->y = niik_image_interpolate_3d_linear(grad_img[2], v->v);
  grad->z = niik_image_interpolate_3d_linear(grad_img[3], v->v);

  *grad = niikpt_unit(*grad);

  /* ICS goes towards WM (positive sign), ocs goes outside*/
  return niikpt_kmul(v->normal, niikpt_dot(v->normal, *grad)*(cortex_id==CORTEX_ICS?1:-1)); /**/

} /* niikcortex_deform_calc_deformation_vertex_gradient_term */


niikpt niikcortex_deform_calc_deformation_vertex_flux_term(kvert *v, nifti_image *img, nifti_image *div_img, int cortex_id) 
{
  double d = niik_image_interpolate_3d_linear(div_img, v->v);

  /*clamp between -1 and 1 */
  if(d<-1.0)     d=-1.0;
  else if(d>1.0) d= 1.0;

  return niikpt_kmul(v->normal, -d); 
} /* niikcortex_deform_calc_deformation_vertex_flux_term */


niikpt niikcortex_deform_calc_abs_thickness_term(kvert *v,kvert *vo,
                      double max_thickness , double min_thickness,
                      double sigma,int cortex_id)
/*move away from vo if thickness is less then min_thickness and towards if thickness if more then max_thickness*/
{
  niikpt diff;
  double dist,proj;
  double attraction,repulsion;
  int dir;

  diff = niikpt_sub(vo->v, v->v);
  proj = niikpt_dot(v->normal, niikpt_unit(diff)); /*project to normal direction?*/
  dist = niikpt_mag(diff);

  /*dir = (cortex_id==CORTEX_ICS?1:-1);*/

  attraction = NIIK_Heaviside(dist-max_thickness,sigma)*proj;
  repulsion  = NIIK_Heaviside(min_thickness-dist,sigma)*proj;

  return niikpt_kmul(v->normal, attraction-repulsion);
}


int niikcortex_deform_calc_deformation_vertex_proximity_term_2sides_test(kvert *v, kface *f, double *min_dist)
/* testing for proximity term
 * min_dist[0] - below (negative vertical distance)
 * min_dist[1] - above (positive vertical distance)
 * -updates the min_dist and vdist (signed vertical distance)
 */
{
  const char *fcname="niikcortex_deform_calc_deformation_vertex_proximity_term_2sides_test";
  const int verbose=0;
  niikpt vec[3],pt;
  double t;
  if(verbose>0) niik_fc_display(fcname,1);

  /* calculate distance from the point to the triangle's plane */
  t = niikpt_det4(1,1,1,0, f->vert[0]->v.x, f->vert[1]->v.x, f->vert[2]->v.x,v->normal.x, f->vert[0]->v.y,f->vert[1]->v.y,f->vert[2]->v.y,v->normal.y,f->vert[0]->v.z,f->vert[1]->v.z,f->vert[2]->v.z,v->normal.z);
  if(fabs(t)<1e-5) return 0;

  t = -niikpt_det4(1,1,1,1,f->vert[0]->v.x,f->vert[1]->v.x, f->vert[2]->v.x,v->v.x, f->vert[0]->v.y,f->vert[1]->v.y,f->vert[2]->v.y,v->v.y,f->vert[0]->v.z,f->vert[1]->v.z,f->vert[2]->v.z,v->v.z) / t;

  if     (t< 0.0 && t<min_dist[0]) return 0;
  else if(t>=0.0 && t>min_dist[1]) return 0;

  pt=niikpt_move_normal(v->v, v->normal, t);

  vec[0]=niikpt_sub(f->vert[1]->v,f->vert[0]->v);
  vec[1]=niikpt_sub(f->vert[2]->v,f->vert[1]->v);
  vec[2]=niikpt_sub(f->vert[0]->v,f->vert[2]->v);

  if(!niikpt_point_is_on_triangle(pt,f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
                                  vec[0],vec[1],vec[2],f->normal))
    return 0;

  if( t<0.0 )
    min_dist[0] = t;
   else
    min_dist[1] = t;
  return 1;
} /* niikcortex_deform_calc_deformation_vertex_proximity_term_2sides_test */


int niikcortex_deform_calc_deformation_vertex_proximity_term_2sides(kvert *v, double *prox_min_dist, bbox *bb, niikpt *pprox)
/* calculates the deformation due to proximal surfaces
 * -but for both sides (front and back) 
 * -it will be the balance between 2 sides
 * it will try to maintain distance of at least *prox_min_dist from each side, and balance between two if unable 
 *
 *   v   is the vertex of interest
 *   prox_min_dist    is the maximum distance of interest
 *   bb  is the bounding box system
 *   pprox   is the output proximity-term displacement
 *
 */
{
  const char *fcname=__func__;
  double
    eps=0.001,
    dmin[2],vdis,min_dist,
    xmin,ymin,zmin,xmax,ymax,zmax,dval;

  int i,j,k,m,n,idx;
  kface *f;
  niikpt vmin,vmax;

  int verbose=1;

  if(verbose>0) {
    if(v->index-1!=g_niikcortex_deform_calc_deformation_vertex_index) {
      verbose=0;
    }
  }

  if(verbose>0) niik_fc_display(fcname,1);
  *pprox=niikpt_zero();
  min_dist = prox_min_dist[0];

  vmin.x = v->v.x - min_dist;
  vmin.y = v->v.y - min_dist;
  vmin.z = v->v.z - min_dist;

  vmax.x = v->v.x + min_dist;
  vmax.y = v->v.y + min_dist;
  vmax.z = v->v.z + min_dist;

  xmin = NIIK_IMAX(0, (v->v.x-min_dist-bb->origin.x)/bb->delta-1);
  ymin = NIIK_IMAX(0, (v->v.y-min_dist-bb->origin.y)/bb->delta-1);
  zmin = NIIK_IMAX(0, (v->v.z-min_dist-bb->origin.z)/bb->delta-1);

  xmax = NIIK_IMIN((v->v.x+min_dist-bb->origin.x)/bb->delta+1, bb->xdim-1);
  ymax = NIIK_IMIN((v->v.y+min_dist-bb->origin.y)/bb->delta+1, bb->ydim-1);
  zmax = NIIK_IMIN((v->v.z+min_dist-bb->origin.z)/bb->delta+1, bb->zdim-1);

  /*dval = dmin[0] = dmin[1] = min_dist + 0.001;*/
  dmin[0] = -min_dist;
  dmin[1] = min_dist;
  /* for each bbox */
  for(k=zmin; k<=zmax; k++) {
    for(j=ymin; j<=ymax; j++) {
      n=xmin + j*bb->xdim + k*bb->area;
      for(i=xmin; i<=xmax; n++,i++) {
        for(m=0; m<bb->ndata[n]; m++) {
          int nn;
          f = bb->data[n][m];
          /*face is too far*/

          if( vmax.x<f->pmin.x || vmax.y<f->pmin.y || vmax.z<f->pmin.z || vmin.x>f->pmax.x || vmin.y>f->pmax.y || vmin.z>f->pmax.z
               /*we found a corner of the current face*/
             || v==f->vert[0] || v==f->vert[1] || v==f->vert[2]) continue;

          /*check the neighours*/
          for(nn=0; nn<3; nn++) {
            int nei;
            for(nei=0; nei<f->vert[nn]->nei; nei++) if(v==f->vert[nn]->neivert[nei]) break;
            if(nei<f->vert[nn]->nei) break;
          }
          if(nn<3) continue;

          if(niikcortex_deform_calc_deformation_vertex_proximity_term_2sides_test(v,f,dmin)) {
            if(verbose>0) {
              fprintf(stdout,"[%s] dist = %9.5f %9.5f  \n",fcname,dmin[0],dmin[1]);
            }
          } /* new min */
        } /* each face */
      }
    }
  } /* each bbox */

  /*balance between two*/
  dval = min_dist/(-dmin[0] + eps) - min_dist/(dmin[1] + eps);

  if     (dval < -1.0) dval=-1.0;
  else if(dval >  1.0) dval=1.0;

  *pprox = niikpt_kmul(v->normal, dval);

  prox_min_dist[0] = -dmin[0];
  prox_min_dist[1] = dmin[1];

  if(verbose>0) {
    fprintf(stdout,"[%s] dist = %9.5f %9.5f    -> %9.5f\n",fcname,dmin[0],dmin[1],dval);
  }
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_calc_deformation_vertex_proximity_term_2sides */


int niikcortex_deform_calc_deformation_vertex_curvature_term(kvert *v, off_curvature_t *crv, niikpt *pcurv) {
  double curv1 = crv->curv1[v->index-1];
  double curv2 = crv->curv2[v->index-1];

  /*Using mean curvature right now*/
  *pcurv = niikpt_kmul(v->normal, -1.0*(curv1+curv2)/2.0 ); 
  return 1;
}


niikpt niikcortex_deform_calc_deformation_vertex(refine_info *dfm, kvert *vi, kvert *vo, int cortex_id, int vidx, double mfactor)
/* -computes deformation for a given vertex /cortex
 */
{
  niikpt
      normal,/*Vertex normal*/
      psum,  /*sum of all forces*/
      pavg,  /*smoothing force (local average)*/
      ptsurf,
      pssurf,
      pimag,
      pgrad,
      pmask,
      pvent,
      pavoid,
      plesm,
      pprox,
      pcurv,
      pdgrad,
      pthick,
      psmooth,
      pprior;
    
  double
        dmin2[2],
        gray,grad,
        bval=0, vval=0, lval=0, vcrm=0;
  kvert *v[2];

  v[0]=vi;
  v[1]=vo;

  //vidx = dfm->vi->index-1;
  bval = vval = lval = vcrm = 0;
  pprior= psmooth= psum = pavg = ptsurf = pssurf = pimag = pmask = pvent = plesm = pprox = pgrad = pavoid = pcurv = pdgrad = pthick = niikpt_zero();
  
  /* if(thklist[vidx]<0.01) {
     if(g_niikcortex_deform_calc_deformation_vertex_index==vidx) {
     fprintf(stdout,"[niikcortex_deform_calc_deformation_vertex] thklist is less than 0.01mm\n"); }
     return psum; } */

  /* Do not apply certain forces if this one is in effect ?*/
  /* brain mask term , surface stays inside of the brain */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_BRAIN_MASK]>0.0) {
    bval = niik_image_interpolate_3d_linear(dfm->brain_mask,v[cortex_id]->v);
    pmask = niikpt_kmul(v[cortex_id]->normal, bval-1.0);
  }

  /* ventricle term, surfaces traverse ventricles (?) */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_VENTRICLE]>0.0) {
    vval = niik_image_interpolate_3d_linear(dfm->ventricle_mask,v[cortex_id]->v);
    pvent = niikpt_kmul(v[cortex_id]->normal, vval);
  }

  /* avoid term, avoid mask */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_AVOID]>0.0) {
    if(dfm->avoid_mask!=NULL) {
      vcrm = niik_image_interpolate_3d_linear(dfm->avoid_mask,v[cortex_id]->v);
      pavoid = niikpt_kmul(v[cortex_id]->normal, -vcrm);
    }
  }

  /* lesion term */
  if(dfm->lesion_mask!=NULL) {
    if(dfm->deform_weights->m[cortex_id][WEIGHT_LESION]>0.0) {
      lval = niik_image_interpolate_3d_linear(dfm->lesion_mask,v[cortex_id]->v);
      plesm = niikpt_kmul(v[cortex_id]->normal, lval);
    }
  }

  /* tangential surface smoothness term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE]>0.0) {
    ptsurf = niikcortex_deform_calc_deformation_vertex_surface_smoothness(
              v[cortex_id]);
  }

  /* simple surface smoothness term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE_SIMPLE]>0.0) {
    pssurf = niikcortex_deform_calc_deformation_vertex_surface_smoothness_simple(
              v[cortex_id]);
  }

  /* another smoothness term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_THICKNESS_SMOOTHNESS]>0.0 ){
    psmooth = niikcortex_deform_calc_deformation_vertex_thickness_smoothness_term(
              v[cortex_id], v[1-cortex_id], dfm->thklist);
  }

  /* image term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_IMAGE]>0.0) {
    pimag = niikcortex_deform_calc_deformation_vertex_image_term(
              v[cortex_id], dfm->img, dfm->grad_img, dfm->gradient_lambda,
              dfm->mean_ics_value*mfactor, dfm->range_ics_value, dfm->mean_ocs_value*mfactor, dfm->range_ocs_value,
              cortex_id, &gray, &grad);
  }

  /* image gradient term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_GRADIENT] > 0.0) {
    pgrad = niikcortex_deform_calc_deformation_vertex_gradient_term(
              v[cortex_id], dfm->img, dfm->grad_img[cortex_id], &normal, cortex_id);
  }

  /* flux gradient term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_DGRADIENT] > 0.0) {
    pdgrad = niikcortex_deform_calc_deformation_vertex_flux_term(
               v[cortex_id], dfm->img, dfm->div_img[cortex_id], cortex_id);
  }

  /* proximity term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_PROXIMITY]>0.0) {
    dmin2[0] = dmin2[1] = dfm->prox_min_dist;
    niikcortex_deform_calc_deformation_vertex_proximity_term_2sides(v[cortex_id], dmin2, dfm->bb, &pprox);
  } /* proximity deformation */

  /* curvature term */
  if(dfm->deform_weights->m[cortex_id][WEIGHT_CURVATURE]>0.0) {
    niikcortex_deform_calc_deformation_vertex_curvature_term(v[cortex_id], &dfm->crv[cortex_id], &pcurv);
  } /* curvature deformation */

  /* abs thickness term*/
  if(dfm->deform_weights->m[cortex_id][WEIGHT_ABS_THICKINESS]>0.0){
    pthick=niikcortex_deform_calc_abs_thickness_term(
        v[cortex_id],v[1-cortex_id], 
        dfm->deform_weights->m[cortex_id][WEIGHT_MAX_THICKNESS], 
        dfm->deform_weights->m[cortex_id][WEIGHT_MIN_THICKNESS], 
        dfm->deform_weights->m[cortex_id][WEIGHT_THICKNESS_SIGMA],
        cortex_id ); 
  }

  /* prior term */
  if(dfm->prior[cortex_id*2]!=NULL && dfm->deform_weights->m[cortex_id][WEIGHT_PRIOR]>0.0) {
    pprior=niikcortex_deform_calc_deformation_prior_term(v[cortex_id], dfm->prior, dfm->grad_prior, cortex_id );
  } /* prior deformation */


  /* vidx is shifted by 1*/
  if(dfm->dfm->debug_trace_log && dfm->dfm->trace_pt_id==vidx)
  {
    fprintf(dfm->dfm->debug_trace_log,  
"%d,%f,%f,%f,%f,%f,%f,\
%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\
%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\
%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
cortex_id,
v[cortex_id]->v.x,v[cortex_id]->v.y,v[cortex_id]->v.z,v[cortex_id]->normal.x,v[cortex_id]->normal.y,v[cortex_id]->normal.z,\
ptsurf.x,ptsurf.y,ptsurf.z,pimag.x,pimag.y,pimag.z,pgrad.x,pgrad.y,pgrad.z,pmask.x,pmask.y,pmask.z,
pvent.x,pvent.y,pvent.z,plesm.x,plesm.y,plesm.z,pprox.x,pprox.y,pprox.z,pavoid.x,pavoid.y,pavoid.z,
pcurv.x,pcurv.y,pcurv.z,pdgrad.x,pdgrad.y,pdgrad.z,pthick.x,pthick.y,pthick.z);
  }/*TODO: add pssurf*/

  /* combine deformation with weighting */
  psum.x =
    ptsurf.x * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE] +
    pssurf.x * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE_SIMPLE] +
    pimag.x  * dfm->deform_weights->m[cortex_id][WEIGHT_IMAGE] +
    pgrad.x  * dfm->deform_weights->m[cortex_id][WEIGHT_GRADIENT] +
    pmask.x  * dfm->deform_weights->m[cortex_id][WEIGHT_BRAIN_MASK] +
    pvent.x  * dfm->deform_weights->m[cortex_id][WEIGHT_VENTRICLE] +
    plesm.x  * dfm->deform_weights->m[cortex_id][WEIGHT_LESION] +
    pprox.x  * dfm->deform_weights->m[cortex_id][WEIGHT_PROXIMITY] +
    pavoid.x   * dfm->deform_weights->m[cortex_id][WEIGHT_AVOID] +
    pcurv.x  * dfm->deform_weights->m[cortex_id][WEIGHT_CURVATURE] +
    pdgrad.x * dfm->deform_weights->m[cortex_id][WEIGHT_DGRADIENT] +
    pthick.x * dfm->deform_weights->m[cortex_id][WEIGHT_ABS_THICKINESS]+
    psmooth.x* dfm->deform_weights->m[cortex_id][WEIGHT_THICKNESS_SMOOTHNESS]+
    pprior.x * dfm->deform_weights->m[cortex_id][WEIGHT_PRIOR];

  psum.y =
    ptsurf.y * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE] +
    pssurf.y * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE_SIMPLE] +
    pimag.y  * dfm->deform_weights->m[cortex_id][WEIGHT_IMAGE] +
    pgrad.y  * dfm->deform_weights->m[cortex_id][WEIGHT_GRADIENT] +
    pmask.y  * dfm->deform_weights->m[cortex_id][WEIGHT_BRAIN_MASK] +
    pvent.y  * dfm->deform_weights->m[cortex_id][WEIGHT_VENTRICLE] +
    plesm.y  * dfm->deform_weights->m[cortex_id][WEIGHT_LESION] +
    pprox.y  * dfm->deform_weights->m[cortex_id][WEIGHT_PROXIMITY] +
    pavoid.y   * dfm->deform_weights->m[cortex_id][WEIGHT_AVOID]+
    pcurv.y  * dfm->deform_weights->m[cortex_id][WEIGHT_CURVATURE] +
    pdgrad.y * dfm->deform_weights->m[cortex_id][WEIGHT_DGRADIENT] +
    pthick.y * dfm->deform_weights->m[cortex_id][WEIGHT_ABS_THICKINESS]+
    psmooth.y* dfm->deform_weights->m[cortex_id][WEIGHT_THICKNESS_SMOOTHNESS] +
    pprior.y * dfm->deform_weights->m[cortex_id][WEIGHT_PRIOR];


  psum.z =
    ptsurf.z * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE] +
    pssurf.z * dfm->deform_weights->m[cortex_id][WEIGHT_SURFACE_SIMPLE] +
    pimag.z  * dfm->deform_weights->m[cortex_id][WEIGHT_IMAGE] +
    pgrad.z  * dfm->deform_weights->m[cortex_id][WEIGHT_GRADIENT] +
    pmask.z  * dfm->deform_weights->m[cortex_id][WEIGHT_BRAIN_MASK] +
    pvent.z  * dfm->deform_weights->m[cortex_id][WEIGHT_VENTRICLE] +
    plesm.z  * dfm->deform_weights->m[cortex_id][WEIGHT_LESION] +
    pprox.z  * dfm->deform_weights->m[cortex_id][WEIGHT_PROXIMITY] +
    pavoid.z   * dfm->deform_weights->m[cortex_id][WEIGHT_AVOID]+
    pcurv.z  * dfm->deform_weights->m[cortex_id][WEIGHT_CURVATURE]+
    pdgrad.z * dfm->deform_weights->m[cortex_id][WEIGHT_DGRADIENT]+
    pthick.z * dfm->deform_weights->m[cortex_id][WEIGHT_ABS_THICKINESS]+
    psmooth.z* dfm->deform_weights->m[cortex_id][WEIGHT_THICKNESS_SMOOTHNESS] +
    pprior.z * dfm->deform_weights->m[cortex_id][WEIGHT_PRIOR];


  return psum;
} /* niikcortex_deform_calc_deformation_vertex */


int niikcortex_deform_calc_deformation( niikcortex_deform * dfm  )
/* -computes the deformation for each vertex
 * -updates dfmlist
 * -uses niikcortex_deform_calc_deformation_vertex
 * -method of integration is determined by num_method
 * output is in dfmlist array, with unit vectors and magnitude stored in .w component
 */
{
  int n,cidx,vidx;
  double dt;
    niikpt
    **op, ***k;
  struct tm *stm;
  time_t ctm;
  const char *fcname=__func__;
  char  tmstr[256];
  const int verbose=0;

  /*TODO: move this up one level?*/
  refine_info  dfm_refine;
  memset(&dfm_refine,0,sizeof(refine_info));

  dfm_refine.img = dfm->t1img;
  dfm_refine.bb  = dfm->bb;
  dfm_refine.brain_mask = dfm->brain_mask;
  dfm_refine.avoid_mask = dfm->avoid_mask;
  dfm_refine.ventricle_mask = dfm->csf_mask;
  dfm_refine.crv = dfm->crv;
  dfm_refine.deform_weights = dfm->weight;
  dfm_refine.div_img = dfm->div_img;
  dfm_refine.grad_img = dfm->grad_img;
  dfm_refine.grad_prior = dfm->grad_prior;
  dfm_refine.prior = dfm->prior;
  dfm_refine.gradient_lambda = dfm->gradient_lambda;
  dfm_refine.lesion_mask = dfm->lesion_mask;
  dfm_refine.mean_ics_value = dfm->tissue_val[4];
  dfm_refine.mean_ocs_value = dfm->tissue_val[5];
  dfm_refine.range_ics_value = dfm->tissue_val[6];
  dfm_refine.range_ocs_value = dfm->tissue_val[7];
  dfm_refine.intWM = dfm->tissue_val[1];
  dfm_refine.vmat = dfm->vmat;
  dfm_refine.thklist = dfm->thklist;
  dfm_refine.dfm = dfm;
  dfm_refine.prox_min_dist = dfm->proximity_min_distance;
  

  dt = dfm->weight->m[CORTEX_ICS][WEIGHT_DT];
  ctm=time(NULL);
  if(verbose>=2) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    switch(dfm->cortex_id) {
    case 1:
      fprintf(stdout,"[%s] start white %s\n",fcname,tmstr);
      break;
    case 2:
      fprintf(stdout,"[%s] start pial %s\n",fcname,tmstr);
      break;
    case 3:
      fprintf(stdout,"[%s] start both %s\n",fcname,tmstr);
      break;
    default:
      fprintf(stderr,"[%s] ERROR: unknown cortex_id %i\n",fcname,dfm->cortex_id);
      return 0;
    }
  }

  switch(dfm->numerical_method) {
  case NIIK_NUM_METHOD_FORWARD_EULER:
    if(verbose>0) fprintf(stdout,"[%s] Forward Euler method\n",fcname);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;

        dfm->dfmlist[cidx][vidx] = niikcortex_deform_calc_deformation_vertex(&dfm_refine, dfm->vmat[CORTEX_ICS][vidx], dfm->vmat[CORTEX_OCS][vidx],cidx,vidx, dfm->use_mf?dfm->mf_list[vidx]:1.0 );
      } /* each cortex */
    } /* each vertex */
    break;

  case NIIK_NUM_METHOD_BACKWARD_EULER:  /* backward euler */
    if(verbose) fprintf(stdout,"[%s] Backward Euler method\n",fcname);
    op = niikpt_matrix(2, dfm->ctx[0]->nvert);
    if(verbose>2) fprintf(stdout,"[%s]   backup vertex positions\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        op[cidx][vidx]=dfm->vmat[cidx][vidx]->v;
      }
    }

    if(verbose>=2) fprintf(stdout,"[%s]   k1 calculation [Backward Euler]\n",fcname);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        dfm->dfmlist[cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex (Backward Eulder) */
    } /* each vertex (Backward Eulder) */
    /* (Backward Eulder)
     * dfm = dt Func(t,Ctx1)
     * Ctx1 = Ctx + k1
     */
    if(verbose>=2) fprintf(stdout,"[%s]   k2 preparation [Backward Euler]\n",fcname);

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = niikpt_move_normal(op[cidx][vidx],dfm->dfmlist[cidx][vidx],dt);
      }
      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      off_update_kobj_kface_pminmax(dfm->ctx[cidx]);

      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
    }
    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Backward Euler]\n",fcname);
    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        dfm->dfmlist[cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */
    /* update dfmlist
     * put the original vertex back */
    if(verbose>=2) fprintf(stdout,"[%s]   output and input update [Backward Euler]\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = op[cidx][vidx];
      }
    }
    /* free memory */
    if(verbose>2) fprintf(stdout,"[%s] free memory\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      free(op[cidx]);
    }
    free(op);
    break;

  case NIIK_NUM_METHOD_MIDPOINT:  /* midpoint */
    if(verbose>0) fprintf(stdout,"[%s] Mid-point method\n",fcname);
    op = niikpt_matrix(2,dfm->ctx[0]->nvert);
    if(verbose>=2) fprintf(stdout,"[%s]   backup vertex positions [Midpoint method]\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        op[cidx][vidx]=dfm->vmat[cidx][vidx]->v;
      }
    }

    /* k1 = dCtx = dt Func(t,Ctx) */
    if(verbose>1) fprintf(stdout,"[%s]   k1 calculation [Midpoint method]\n",fcname);
    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        dfm->dfmlist[cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /*
     * k2 = dCtx2 = dt Func(t,Ctx1)
     * Ctx1 = Ctx + 0.5 * k1
     */
    if(verbose>=3) fprintf(stdout,"[%s]   k2 preparation [Midpoint method]\n",fcname);

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = niikpt_move_normal(op[cidx][vidx],dfm->dfmlist[cidx][vidx],0.5*dt);
      }
      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      off_update_kobj_kface_pminmax(dfm->ctx[cidx]);

      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);

    }
    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Midpoint method]\n",fcname);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        dfm->dfmlist[cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /* update dfmlist
     * put the original vertex back */
    if(verbose>=2) fprintf(stdout,"[%s]   output and input update [Midpoint method]\n",fcname);

    for(cidx=0; cidx<2; cidx++) {
      if(cidx==0 && dfm->cortex_id==2) continue;
      if(cidx==1 && dfm->cortex_id==1) continue;
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = op[cidx][vidx];
      }
    }

    /* free memory */
    if(verbose>=2) fprintf(stdout,"[%s] free memory [Midpoint method]\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      free(op[cidx]);
    }
    free(op);
    break;


  case NIIK_NUM_METHOD_RUNGE_KUTTA:

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);

    if(verbose>=1) fprintf(stdout,"[%s] Runge-Kutta 4th order method %s\n",fcname,tmstr);
    op = niikpt_matrix(2,dfm->ctx[0]->nvert);
    k = (niikpt ***)calloc(4,sizeof(niikpt **));

    for(n=0; n<4; n++) {
      if((k[n] = niikpt_matrix(2,dfm->ctx[0]->nvert))==NULL) {
        fprintf(stderr,"[%s] ERROR: niikpt_matrix\n",fcname);
        return 0;
      }
    }

    if(verbose>=2) fprintf(stdout,"[%s]   backup vertex positions [Runge-Kutta]\n",fcname);
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        op[cidx][vidx]=dfm->vmat[cidx][vidx]->v;
      }
    }

    /* k1 = dCtx = dt Func(t,Ctx) */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k1 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[0][0][vidx] = k[0][1][vidx] = niikpt_zero();
        continue;
      }
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        k[0][cidx][vidx] = niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /*
     * k2 = dCtx2 = dt Func(t,Ctx1)
     * Ctx1 = Ctx + 0.5 * k1
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k2 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = niikpt_move_normal(op[cidx][vidx],k[0][cidx][vidx],0.5*dt);
      }

      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      off_update_kobj_kface_pminmax(dfm->ctx[cidx]);

      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
    }

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[1][0][vidx] = k[1][1][vidx] = niikpt_zero();
        continue;
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;

        k[1][cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /*
     * k3 = dCtx3 = dt Func(t,Ctx2)
     * Ctx2 = Ctx1 + 0.5 * k2
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k3 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = niikpt_move_normal(op[cidx][vidx],k[1][cidx][vidx],0.5*dt);
      }
      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      off_update_kobj_kface_pminmax(dfm->ctx[cidx]);

      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
    } /* each surface */

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k3 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private(cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[2][0][vidx] = k[2][1][vidx] = niikpt_zero();
        continue;
      }
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        k[2][cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /*
     * k4 = dCtx4 = dt Func(t,Ctx3)
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k4 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v = niikpt_move_normal(op[cidx][vidx],k[2][cidx][vidx],0.5*dt);
      }
      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      off_update_kobj_kface_pminmax(dfm->ctx[cidx]);
      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
    } /* each surface */

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k4 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[3][0][vidx] = k[3][1][vidx] = niikpt_zero();
        continue;
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;
        k[3][cidx][vidx]=niikcortex_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],dfm->vmat[CORTEX_OCS][vidx],cidx,vidx,dfm->use_mf?dfm->mf_list[vidx]:1.0);
      } /* each cortex */
    } /* each vertex */

    /* update dfmlist
     * put the original vertex back */
    if(verbose>=3) {
      ctm=time(NULL);
      stm=localtime(&ctm);
      strftime(tmstr,256,"%Y-%m-%d %T",stm);
      fprintf(stdout,"[%s]   output and input update [Runge-Kutta] %s\n",fcname,tmstr);
    }

    #pragma omp parallel for private (vidx) num_threads(2)
    for(cidx=0; cidx<2; cidx++) {
      if(cidx==0 && dfm->cortex_id==2) continue;
      if(cidx==1 && dfm->cortex_id==1) continue;

      #pragma omp parallel for
      for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        dfm->vmat[cidx][vidx]->v   = op[cidx][vidx];

        dfm->dfmlist[cidx][vidx].x = (k[0][cidx][vidx].x + 2*k[1][cidx][vidx].x + 2*k[2][cidx][vidx].x + k[3][cidx][vidx].x) / 6.0 * dt;
        dfm->dfmlist[cidx][vidx].y = (k[0][cidx][vidx].y + 2*k[1][cidx][vidx].y + 2*k[2][cidx][vidx].y + k[3][cidx][vidx].y) / 6.0 * dt;
        dfm->dfmlist[cidx][vidx].z = (k[0][cidx][vidx].z + 2*k[1][cidx][vidx].z + 2*k[2][cidx][vidx].z + k[3][cidx][vidx].z) / 6.0 * dt;

        // if( fabs(dfmlist[cidx][vidx].x)<g_niikcortex_min_deform_distance &&
        //     fabs(dfmlist[cidx][vidx].y)<g_niikcortex_min_deform_distance &&
        //     fabs(dfmlist[cidx][vidx].z)<g_niikcortex_min_deform_distance) {
        //   dfmlist[cidx][vidx].w = 0.0;
        // } else {
        //   dfmlist[cidx][vidx].w = niikpt_mag(dfmlist[cidx][vidx]);
        // }
      }
    }

    /* free memory */
    if(verbose>=2) {
      ctm=time(NULL);
      stm=localtime(&ctm);
      strftime(tmstr,256,"%Y-%m-%d %T",stm);
      fprintf(stdout,"[%s] free memory [Runge-Kutta] %s\n",fcname,tmstr);
    }

    for(n=0; n<4; n++) {
      for(cidx=0; cidx<2; cidx++) {
        free(k[n][cidx]);
      }
      free(k[n]);
    }
    free(k);
    for(cidx=0; cidx<2; cidx++) {
      free(op[cidx]);
    }
    free(op);
    break;

  default:
    fprintf(stderr,"[%s] ERROR: not implemented yet %i\n",fcname,dfm->numerical_method);
    return 0;
  }

  /*convert dfmlist into unit vectors with .w showing magnitude*/
  for(cidx=0; cidx<2; cidx++) {
    if(cidx==0 && dfm->cortex_id==2) continue;
    if(cidx==1 && dfm->cortex_id==1) continue;
    if(dfm->weight->m[cidx][WEIGHT_UPDATE_SIGMA]>0.0) {
      off_surface_field_smooth_using_vert(dfm->ctx[cidx], dfm->dfmlist[cidx], dfm->weight->m[cidx][WEIGHT_UPDATE_SIGMA]);
    }
  }

  #pragma omp parallel for private (vidx) num_threads(2)
  for(cidx=0; cidx<2; cidx++) {
    if(cidx==0 && dfm->cortex_id==2) continue;
    if(cidx==1 && dfm->cortex_id==1) continue;

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[cidx]->nvert; vidx++) {
        if( fabs(dfm->dfmlist[cidx][vidx].x)<g_niikcortex_min_deform_distance &&
            fabs(dfm->dfmlist[cidx][vidx].y)<g_niikcortex_min_deform_distance &&
            fabs(dfm->dfmlist[cidx][vidx].z)<g_niikcortex_min_deform_distance) {
          dfm->dfmlist[cidx][vidx].w = 0.0;
        } else {
          double d=niikpt_mag(dfm->dfmlist[cidx][vidx]);
          dfm->dfmlist[cidx][vidx] = niikpt_kmul(dfm->dfmlist[cidx][vidx], 1.0/d); /*make a unit vector with length in w*/
          dfm->dfmlist[cidx][vidx].w = d;
        }
    }
  }

  if(verbose>=1) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    fprintf(stdout,"[%s] exit function %s\n",fcname,tmstr);
  }
  return 1;
} /* niikcortex_deform_calc_deformation */

/*TODO: adapt this code to work without assumption that object intensity is always max (?)*/
int niikcortex_deform_cortex_calc_multiplication_factor(nifti_image *img, 
    nifti_image *brain_mask, kobj **ctx, double intWM, double mflim, double *mf_list) {
  char *pstr;
  const char *fcname=__func__;
  double  xmin=-3.0,
          xmax=-0.1,
          xdel=0.2,
          mf_lim[2];
  kvert *vi,*vo;
  niikpt pt,normal;
  int n,
    vidx,nvert,
    verbose=niik_verbose();
    niikvec *x,*y;

  if(verbose>=1) niik_fc_display(fcname,1);

  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((ctx==NULL),fcname,"cortical model is null");
  NIIK_RET0((ctx[0]==NULL),fcname,"WM surface is null");
  NIIK_RET0((ctx[1]==NULL),fcname,"Pial surface is null");
  NIIK_RET0((mf_list==NULL),fcname,"mf_list is null");

  nvert=ctx[0]->nvert;
  if(nvert!=ctx[1]->nvert) {
    fprintf(stderr,"[%s] ERROR: #v is different between WM[%i] and Pial[%i]\n",fcname,ctx[0]->nvert,ctx[1]->nvert);
    return 0;
  }
  /* CALCULATE MULTIPLICATION FACTOR */
  mf_lim[0] = (mflim>1.0)?1.0/mflim:mflim;
  mf_lim[1] = 1.0 / mf_lim[0];
  if(verbose>=1) fprintf(stdout,"[%s] scaling factor limits: %8.5f %8.5f\n",fcname,mf_lim[0],mf_lim[1]);
  if(verbose>=1) fprintf(stdout,"[%s] WM intensity           %8.2f\n",fcname,intWM);
  if((x=niikvec_init_range(xmin,xmax,xdel))==NULL) {
    fprintf(stderr,"[%s] ERROR: niikvec_init_range\n",fcname);
    return 0;
  }
  NIIK_RET0(((y=niikvec_init(x->num))==NULL),fcname,"niikvec_init");

  if(verbose>=1) fprintf(stdout,"[%s] interp vector          %5.2f %5.2f %5.2f %5i\n",fcname,x->v[0],x->v[1]-x->v[0],x->v[x->num-1],x->num);
  for(vi=ctx[0]->vert,vo=ctx[1]->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
    normal = niikpt_unit(niikpt_sub(vo->v,vi->v));
    if(fabs(normal.x)+fabs(normal.y)+fabs(normal.z)<1e-3) {
      mf_list[vidx]=1.0;
      continue;
    }
    if(!niik_image_interp_along_normal(img,NIIK_INTERP_LINEAR,vi->v,normal,x,y)) {
      fprintf(stderr,"[%s] ERROR: niik_image_interp_along_normal, vidx=%i\n",fcname,vidx);
      return 0;
    }
    /*debug*/
    if(g_niikcortex_deform_calc_deformation_vertex_index==vidx) {
      for(n=0; n<x->num; n++) {
        pstr=niikpt_display_xyz(niikpt_move_normal(vi->v,normal,x->v[n]));
        fprintf(stdout,"[%s] %5.1f %s   %8.2f %3.0f\n",
                fcname,x->v[n],pstr,
                niik_image_interpolate_3d_linear(img,pt),
                niik_image_interpolate_3d_nn(brain_mask,pt));
        free(pstr);
      }
    } /* display info */
    
    /*TODO: this is implicit rule that WM is the highest intensity ?*/
    mf_list[vidx] = NIIK_DMINMAX(niik_get_max_from_double_vector(y->v,x->num) / intWM, mf_lim[0], mf_lim[1]);
  }
  if(verbose>0) fprintf(stdout,"[niikcortex_deform_cortex]   smoothing multiplication factor\n");
  if(!off_surface_smooth_using_vert(ctx[0],mf_list,2,0.5)) {
    fprintf(stderr,"[%s] ERROR: off_surface_smooth_using_vert\n",fcname);
    return 0;
  }

  if(verbose>0)
    fprintf(stdout,"[%s]   %9.6f +/- %9.6f  (%9.6f %9.6f %9.6f %9.6f)\n",
            fcname,
            niik_get_mean_from_double_vector(mf_list,nvert),
            niik_get_stdv_from_double_vector(mf_list,nvert),
            niik_get_min_from_double_vector (mf_list,nvert),
            niik_get_percentile_from_double_vector(mf_list,nvert,0.25),
            niik_get_percentile_from_double_vector(mf_list,nvert,0.75),
            niik_get_max_from_double_vector(mf_list,nvert));

  if(g_niikcortex_deform_calc_deformation_vertex_index>=0) {
    fprintf(stdout,"[%s] scaling factor[%i] %9.4f\n",
            fcname,
            g_niikcortex_deform_calc_deformation_vertex_index,
            mf_list[g_niikcortex_deform_calc_deformation_vertex_index]);
  }
  x=niikvec_free(x);
  y=niikvec_free(y);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_cortex_calc_multiplication_factor */


/**************
 * allocate or reallocate transient structures
 * 
 * 
 * */
int niicortex_deform_prepare(niikcortex_deform * dfm)
{
  const char* fcname=__func__;
  kvert *vi,*vo;
  int vidx,cidx;
  dfm->nvert = dfm->ctx[0]->nvert;

  if(dfm->nvert!=dfm->ctx[1]->nvert) {
    fprintf(stderr,"[%s] ERROR: different #vertex %i %i\n",fcname,dfm->ctx[0]->nvert,dfm->ctx[1]->nvert);
    return 0;
  }

  /* nonctx labels
   * -ctx_label is a list of zeros and non-zeros for non-cortex and cortex, respectively
   */
  if(dfm->nonctx_mask!=NULL) {
    NIIK_RET0(((dfm->ctx_label = niikcortex_non_cortex_label(dfm->ctx_label,dfm->nonctx_mask,dfm->ctx[0],dfm->ctx[1],0.0,0.1,2))==NULL),
              fcname,"niikcortex_non_cortex_label");
  } else {
    int n;
    dfm->ctx_label = (int *)realloc(dfm->ctx_label,dfm->nvert*sizeof(int));
    for(n=0; n<dfm->nvert; n++) dfm->ctx_label[n] = 1;
  } 

  dfm->indicator = (indicator_t*)realloc(dfm->indicator,dfm->nvert*sizeof(indicator_t));
 
  for(cidx=0; cidx<2; cidx++) {
    dfm->dfmlist[cidx] = (niikpt *)realloc(dfm->dfmlist[cidx], dfm->nvert*sizeof(niikpt));
    dfm->vmat[cidx]    = (kvert **)realloc(dfm->vmat[cidx], dfm->nvert*sizeof(kvert *));
    dfm->rd[cidx]  =     (double*) realloc(dfm->rd[cidx], dfm->nvert*sizeof(double));
    dfm->dd[cidx]  =     (double*) realloc(dfm->dd[cidx], dfm->nvert*sizeof(double));

    re_init_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
  }

  for(vi=dfm->ctx[CORTEX_ICS]->vert,vo=dfm->ctx[CORTEX_OCS]->vert,vidx=0;
      vi!=NULL;
      vi=vi->next,vo=vo->next,vidx++) {
    /*assign labels*/
    vi->idata=CORTEX_ICS;
    vo->idata=CORTEX_OCS;
    dfm->vmat[0][vidx]=vi;
    dfm->vmat[1][vidx]=vo;
  }
 
  dfm->thklist = (double*) realloc(dfm->thklist, dfm->nvert*sizeof(double));

  /* CALCULATE MULTIPLICATION FACTOR */
  dfm->mf_lim[0] = (dfm->mflim>1.0)?1.0/dfm->mflim:dfm->mflim;
  dfm->mf_lim[1] = 1.0 / dfm->mf_lim[0];
  dfm->mf_list = (double *)realloc(dfm->mf_list, dfm->nvert*sizeof(double));

  if(dfm->use_mf) {
    if(!niikcortex_deform_cortex_calc_multiplication_factor(dfm->t1img, dfm->brain_mask, dfm->ctx, dfm->tissue_val[1], dfm->mf_lim[0], dfm->mf_list)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_deform_cortex_calc_multiplication_factor\n",fcname);
      return 0;
    }
  }

  // for(vi=dfm->ctx[0]->vert,vo=dfm->ctx[1]->vert,vidx=0; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
  //   niikpt normal = niikpt_unit(niikpt_sub(vo->v,vi->v));
  //   double dist,dmax;
  //   for(dist=0.1,dmax=0.0; dist<=2.01; dist+=0.1) {
  //     niikpt pt = niikpt_move_normal(vi->v,normal,-dist);
  //     if(niik_image_interpolate_3d_nn(dfm->brain_mask,pt)>0) {
  //       dmax = NIIK_DMAX(dmax, niik_image_interpolate_3d_linear(dfm->t1img, pt));
  //     }
  //   } /* step back into white matter */
  //   dfm->mf_list[vidx] = NIIK_DMINMAX(dmax / dfm->tissue_val[1], dfm->mf_lim[0], dfm->mf_lim[1]); /*intWM=dfm->tissue_val[1]*/
  //   /* mf_list[vidx] = 0.5 + 0.5*mf_list[vidx]; */
  // }

  // fprintf(stdout,"[%s] mflist:  %9.6f +/- %9.6f  (%9.6f %9.6f %9.6f %9.6f)\n",fcname,
  //           niik_get_mean_from_double_vector(dfm->mf_list,dfm->nvert),
  //           niik_get_stdv_from_double_vector(dfm->mf_list,dfm->nvert),
  //           niik_get_min_from_double_vector (dfm->mf_list,dfm->nvert),
  //           niik_get_percentile_from_double_vector(dfm->mf_list,dfm->nvert,0.25),
  //           niik_get_percentile_from_double_vector(dfm->mf_list,dfm->nvert,0.75),
  //           niik_get_max_from_double_vector(dfm->mf_list,dfm->nvert));


  // if(!off_surface_smooth_using_vert(dfm->ctx[0],dfm->mf_list,2,0.5)) {
  //   fprintf(stderr,"[%s] ERROR: off_surface_smooth_using_vert\n",fcname);
  //   return 0;
  // }

#if _OPENMP>=201307
  #pragma omp simd
#endif
  for(vidx=0; vidx<dfm->nvert; vidx++) {
    dfm->thklist[vidx] = niikpt_distance(dfm->vmat[0][vidx]->v,dfm->vmat[1][vidx]->v);
  }

  return 1;
}


/************************************************************
 *
 * int niikcortex_deform_cortex( ... )
 *
 * -deforms the cortical model
 *
 * see also:
 *   int    niikcortex_deform_calc_deformation  ( ... )
 *   niikpt niikcortex_deform_calc_deformation_vertex ( ... )
 *   int    niikcortex_deform_apply_deformation ( ... )
 *
 *
 ************************************************************/
int niikcortex_deform_cortex(niikcortex_deform * dfm)
/*                            )*/
/* -main deformation function
 * -deforms the cortical model (ics and ocs)
 * -returns 1 for success and zero for failure
 *
 * -img is T1W MRI
 * -lesion_mask is optional
 * -other masks are required and should be clear
 * -ctx_label is a list of zeros and non-zeros for non-cortex and cortex, respectively
 * -deform_weights are the deformation weighting
 * -tol is the tolerance (stop criterion) usually 1%
 * -maxiter is the maximum iteration (another stop criterion)
 * -num_method is the numerical method, 4th order Runge-Kutta is recommended (slower)
 * -tissue_values intensity values
 * -FWHM for calculating the gradient image
 * -gradient_lambda weighting for gradient
 * -dfm_ctx is a flag for cortex deformation
 *   1 = deform white surface only
 *   2 = deform pial surface only
 *   3 = deform both
 *   otherwise, error
 *
 */
{
  nifti_image
    **grad_img[4],
    *div_img[2],
    *blur_img[4],
    *tmpimg;
  nifti_image **grad_prior[4];
  bbox *bb;
  kface *f;
  kvert *vi,*vo;
    
  niikmat *invmat;

  niikpt pt,normal;
  struct tm *stm;
  time_t ctm,ref_time;
  char tmstr[256];

  double
    dist,
    dval,
    dmin,dmax,dmean,dstdv,
    dmax_dfm, /* cumulative maximum deformation for bbox update */
    dratio[2],
    dt,
    intWM,
    mean_ics_value,range_ics_value,
    mean_ocs_value,range_ocs_value;

  int
    n,nlo,nhi,
    ci[2],
    xsc=0,
    *this_label,
    cidx,vidx,
    iter,iter2,
    verbose=niik_verbose(),
    remesh_fail;
  
  int cortex_id=dfm->cortex_id;
  int debug_tracing=0;
  const char * falcon_trace_log=NULL;
  const char * falcon_trace_log_id=NULL;
  cortex_tracing_info trace;

  /*curvature storage*/
  /*niikmat *deform_weights=dfm->weight;*/

  char fname[2048];
  const char *fcname=__func__;

  debug_tracing = falcon_tracing_init(dfm->t1img, &trace);

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  if(verbose>=1) niik_fc_display(fcname,1);
  NIIK_RET0((     dfm->t1img==NULL),fcname,"t1img is null");
  NIIK_RET0((dfm->brain_mask==NULL),fcname,"brain mask is null");
  NIIK_RET0((  dfm->csf_mask==NULL),fcname,"csf mask is null");
  if(cortex_id>3 || cortex_id<1) {
    fprintf(stderr,"[%s] ERROR: invalid cortex_id, %i\n",fcname, cortex_id);
    return 0;
  }


  mean_ics_value  = dfm->tissue_val[4];
  mean_ocs_value  = dfm->tissue_val[5];
  range_ics_value = dfm->tissue_val[6];
  range_ocs_value = dfm->tissue_val[7];
  intWM           = dfm->tissue_val[1];

  dt              = dfm->weight->m[0][0];

  if(verbose>0) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"  t1w image            %s\n",dfm->t1img->fname);
    fprintf(stdout,"  brain mask           %s\n",dfm->brain_mask->fname);
    fprintf(stdout,"  white surface        %s\n",dfm->ctx[CORTEX_ICS]->fname);
    fprintf(stdout,"  pial surface         %s\n",dfm->ctx[CORTEX_OCS]->fname);

    if(dfm->prior[0])
      fprintf(stdout,"  WM prior             %s\n", dfm->prior[0]->fname);
    if(dfm->prior[1])
      fprintf(stdout,"  GM prior             %s\n", dfm->prior[1]->fname);
    if(dfm->prior[2])
      fprintf(stdout,"  CSF prior            %s\n", dfm->prior[2]->fname);

    fprintf(stdout,"  surface vfe          %i %i %i\n",dfm->ctx[CORTEX_ICS]->nvert,dfm->ctx[CORTEX_ICS]->nface,dfm->ctx[CORTEX_ICS]->nedge);
    fprintf(stdout,"  deform time step     %-7.4f\n",dt);
    fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n",dfm->apply_step);
    fprintf(stdout,"  surface sm weights   %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_SURFACE],dfm->weight->m[CORTEX_OCS][WEIGHT_SURFACE]);
    fprintf(stdout,"  image weights        %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_IMAGE],dfm->weight->m[CORTEX_OCS][WEIGHT_IMAGE]);
    fprintf(stdout,"  thickness sm weights %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_THICKNESS_SMOOTHNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_THICKNESS_SMOOTHNESS]);
    fprintf(stdout,"  brainmask weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_BRAIN_MASK],dfm->weight->m[CORTEX_OCS][WEIGHT_BRAIN_MASK]);
    fprintf(stdout,"  ventricle weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_VENTRICLE],dfm->weight->m[CORTEX_OCS][WEIGHT_VENTRICLE]);
    fprintf(stdout,"  avoid weights        %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_AVOID],dfm->weight->m[CORTEX_OCS][WEIGHT_AVOID]);
    fprintf(stdout,"  lesion mask weights  %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_LESION],dfm->weight->m[CORTEX_OCS][WEIGHT_LESION]);
    fprintf(stdout,"  proximity weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_PROXIMITY],dfm->weight->m[CORTEX_OCS][WEIGHT_PROXIMITY]);
    fprintf(stdout,"  abs thickness        %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_ABS_THICKINESS],dfm->weight->m[CORTEX_OCS][WEIGHT_ABS_THICKINESS]);
    fprintf(stdout,"  gradient weights     %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_GRADIENT],dfm->weight->m[CORTEX_OCS][WEIGHT_GRADIENT]);
    fprintf(stdout,"  curvature weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_CURVATURE],dfm->weight->m[CORTEX_OCS][WEIGHT_CURVATURE]);

    fprintf(stdout,"  min thicknes         %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_MIN_THICKNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_MIN_THICKNESS]);
    fprintf(stdout,"  max thicknes         %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_MAX_THICKNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_MAX_THICKNESS]);
    fprintf(stdout,"  thicknes  sigma      %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_THICKNESS_SIGMA],dfm->weight->m[CORTEX_OCS][WEIGHT_THICKNESS_SIGMA]);
    fprintf(stdout,"  prior weight         %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_PRIOR],dfm->weight->m[CORTEX_OCS][WEIGHT_PRIOR]);

    fprintf(stdout,"  use mf               %s\n", dfm->use_mf?"Yes":"No");
    fprintf(stdout,"  max iter             %i\n", dfm->iter);
    fprintf(stdout,"  max remesh fail ctr  %i\n", dfm->remesh_fail_max);
    fprintf(stdout,"  white matter         %9.2f\n",intWM);
    fprintf(stdout,"  white surface value  %9.2f %9.2f\n",mean_ics_value,range_ics_value);
    fprintf(stdout,"  pial surface value   %9.2f %9.2f\n",mean_ocs_value,range_ocs_value);
    fprintf(stdout,"  proximity min dist   %-7.2f\n",dfm->proximity_min_distance);
    if(dfm->gradient_lambda[0]>=0) {
      fprintf(stdout,"  ics gradient lambda      %7.3f\n",dfm->gradient_lambda[0]);
    } else {
      fprintf(stdout,"  ics gradient lambda      TBD\n");
    }
    if(dfm->gradient_lambda[1]>=0) {
      fprintf(stdout,"  ocs gradient lambda      %7.3f\n",dfm->gradient_lambda[1]);
    } else {
      fprintf(stdout,"  ocs gradient lambda      TBD\n");
    }
    fprintf(stdout,"  tolerance            %-7.4f\n",dfm->tolerance);
    if(dfm->ctx_label!=NULL)
      fprintf(stdout,"  non-cortex label     %-i / %i\n",(int)niik_count_zero_from_int_vector(dfm->ctx_label,dfm->nvert),dfm->nvert);
    else
      fprintf(stdout,"  non-cortex label     not used\n");

    switch(cortex_id) {
    case 1:
      fprintf(stdout,"  deform cortex        white surface\n");
      break;
    case 2:
      fprintf(stdout,"  deform cortex        pial surface\n");
      break;
    case 3:
      fprintf(stdout,"  deform cortex        white and pial surfaces\n");
      break;
    default:
      fprintf(stderr,"[%s] ERROR: unknown cortex_id, %i\n",fcname,cortex_id);
      return 0;
    }
    fprintf(stdout,"  numerical method     %s\n", niik_numerical_method_string(dfm->numerical_method));
    fflush(stdout);
  }

  if(dfm->debug_pt!=NULL)
    if(dfm->debug_pt[0]>0)
      fprintf(stdout,"  debug position       %9.2f %9.2f %9.2f\n",dfm->debug_pt[1]*dfm->t1img->dx, dfm->debug_pt[2]*dfm->t1img->dy, dfm->debug_pt[3]*dfm->t1img->dz);

  if((falcon_trace_log = getenv("FALCON_TRACE_LOG")) && (falcon_trace_log_id = getenv("FALCON_TRACE_LOG_ID")))
  {
    dfm->debug_trace_log=fopen(falcon_trace_log,"w");
    dfm->trace_pt_id=atoi(falcon_trace_log_id)-1; /*the code inside is 0-based*/
    fprintf(dfm->debug_trace_log,"cidx,\
x,y,z,nx,ny,nz,\
psurf.x,psurf.y,psurf.z,pimag.x,pimag.y,pimag.z,\
pgrad.x,pgrad.y,pgrad.z,pmask.x,pmask.y,pmask.z,\
pvent.x,pvent.y,pvent.z,plesm.x,plesm.y,plesm.z,\
pprox.x,pprox.y,pprox.z,pavoid.x,pavoid.y,pavoid.z,\
pcurv.x,pcurv.y,pcurv.z,pdgrad.x,pdgrad.y,pdgrad.z,\
pthick.x,pthick.y,pthick.z\n");
  }


  if(verbose>0) fprintf(stdout,"[%s] gradient image\n",fcname);

  for(n=0; n<2; n++) {
    NIIK_RET0(((blur_img[n  ] = niik_image_copy_as_type(dfm->t1img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
    NIIK_RET0(((blur_img[n+2] = niik_image_copy_as_type(dfm->t1img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");

    if( dfm->gradient_FWHM[n]>0.0 ) {
      if(verbose>1) fprintf(stdout,"[%s] gradient  gaussian filter\n",fcname);
      NIIK_RET0((!niik_image_filter_gaussian_update(blur_img[n],11,dfm->gradient_FWHM[n])),fcname,"niik_image_filter_gaussian_update");
    } else if(verbose>1) {
      fprintf(stdout,"[%s]   no gradient gaussian filter\n",fcname);
    }

    if( dfm->divergence_FWHM[n]>0.0 ) {
      if(verbose>1) fprintf(stdout,"[%s] divergence  gaussian filter\n",fcname);
      NIIK_RET0((!niik_image_filter_gaussian_update(blur_img[n+2],11,dfm->divergence_FWHM[n])),fcname,"niik_image_filter_gaussian_update");
    } else if(verbose>1) {
      fprintf(stdout,"[%s]   no gradient gaussian filter\n",fcname);
    }

    if(verbose>1) fprintf(stdout,"[%s]   sobel filter\n",fcname);
    NIIK_RET0(((grad_img[n  ] = niik_image_sobel_filters_with_mag(blur_img[n  ]))==NULL), fcname, "niik_image_sobel_filters_with_mag");
    NIIK_RET0(((grad_img[n+2] = niik_image_sobel_filters_with_mag(blur_img[n+2]))==NULL), fcname, "niik_image_sobel_filters_with_mag");

    /*divergence of the gradient field*/
    NIIK_RET0(((div_img[n] = niik_image_divergence(grad_img[n+2],0))==NULL),fcname,"niik_image_divergence");
  }

  for(n=0; n<4; n++)
    blur_img[n] = niik_image_free(blur_img[n]);

  for(n=0; n<4; n++) {
    if(dfm->prior[n]) {
        nifti_image *blur_prior;
        NIIK_RET0(((blur_prior = niik_image_copy_as_type(dfm->prior[n], NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");

        if(dfm->prior_FWHM[n/2]>0.0) { 
          if(verbose>1) fprintf(stdout,"[%s] gradient  gaussian filter\n",fcname);
          NIIK_RET0((!niik_image_filter_gaussian_update(blur_prior,11, dfm->prior_FWHM[n/2] )),fcname,"niik_image_filter_gaussian_update");
        } else if(verbose>1) {
          fprintf(stdout,"[%s]   no gradient gaussian filter\n",fcname);
        }

        if(verbose>1) fprintf(stdout,"[%s]   sobel filter\n",fcname);
        NIIK_RET0(((grad_prior[n] = niik_image_sobel_filters_with_mag(blur_prior))==NULL),fcname,"niik_image_sobel_filters_with_mag");

        niik_image_free(blur_prior);
    } else {
        grad_prior[n]=NULL;
    }
  }

  if(verbose>1) fprintf(stdout,"[%s]   sobel filtered\n",fcname);

  for(n=0; n<2; n++) {
    if(dfm->gradient_lambda[n]<0) {
      if(verbose>1) fprintf(stdout,"[%s] gradient lambda calculation\n",fcname);
      dfm->gradient_lambda[n] = niik_image_get_percentile(grad_img[n][0],dfm->brain_mask,0.7);
      fprintf(stdout,"[%s] gradient lambda  %d   %7.3f\n",fcname, n, dfm->gradient_lambda[n]);
    } else {
      if(verbose>1) fprintf(stdout,"[%s] skip gradient lambda calculation\n",fcname);
    }
  }

  if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] allocate memory\n");

  /*make sure we have normals*/
  for(cidx=0; cidx<2; cidx++) {
    off_update_kobj_face_normal(dfm->ctx[cidx]);
    off_update_kobj_vert_normal(dfm->ctx[cidx]);
    off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
  }

  NIIK_RET0(!niicortex_deform_prepare(dfm),fcname,"niicortex_deform_prepare");

  if(verbose>=3) fprintf(stdout,"[%s] bbox initialization\n",fcname);
  bb=off_bbox_init(dfm->bbox_depth,320);

  if(verbose>=1) fprintf(stdout,"[%s] bbox initialization %7.4f %i\n",fcname,bb->delta,bb->depth);

  /* VERTEX */
  if(dfm->debug_pt && dfm->debug_pt[0]>0) {
    pt.x = dfm->debug_pt[1] * dfm->t1img->dx; /*voxel to world here?*/
    pt.y = dfm->debug_pt[2] * dfm->t1img->dy;
    pt.z = dfm->debug_pt[3] * dfm->t1img->dz;
    for(vi=dfm->ctx[0]->vert,vo=dfm->ctx[1]->vert,dmin=1e9; vi!=NULL; vi=vi->next,vo=vo->next,vidx++) {
      dval = niikpt_distance(vi->v,pt);
      if(dval<dmin) {
        dmin = dval;
        g_niikcortex_deform_calc_deformation_vertex_index=vi->index;
      }
      dval = niikpt_distance(vo->v,pt);
      if(dval<dmin) {
        dmin = dval;
        g_niikcortex_deform_calc_deformation_vertex_index=vi->index;
      }
    }
  } /* debug vertex */

  /* calculate the inverse matrix */
  if(dfm->regmat==NULL) {
    if((invmat=niikmat_identity(4,4))==NULL) {
      fprintf(stderr,"[%s] ERROR: niikmat_identity\n",fcname);
      return 0;
    }
  } else {
    if((invmat=niikmat_inverse(dfm->regmat))==NULL) {
      fprintf(stderr,"[%s] ERROR: niikmat_inverse\n",fcname);
      return 0;
    }
  }
  dfm->grad_img=grad_img;
  dfm->div_img=div_img;
  dfm->grad_prior=grad_prior;

  if(verbose>0) fprintf(stdout,"[niikcortex_deform_cortex]   smoothing multiplication factor\n");

  if(verbose>0 && dfm->use_mf)
    fprintf(stdout,"[niikcortex_deform_cortex]   %9.6f +/- %9.6f  (%9.6f %9.6f %9.6f %9.6f)\n",
            niik_get_mean_from_double_vector(dfm->mf_list,dfm->nvert),
            niik_get_stdv_from_double_vector(dfm->mf_list,dfm->nvert),
            niik_get_min_from_double_vector (dfm->mf_list,dfm->nvert),
            niik_get_percentile_from_double_vector(dfm->mf_list,dfm->nvert,0.25),
            niik_get_percentile_from_double_vector(dfm->mf_list,dfm->nvert,0.75),
            niik_get_max_from_double_vector(dfm->mf_list,dfm->nvert));

  if(g_niikcortex_deform_calc_deformation_vertex_index>=0 && dfm->use_mf) {
    if(verbose>0) fprintf(stdout,"   scaling factor %9.4f\n",dfm->mf_list[g_niikcortex_deform_calc_deformation_vertex_index]);
  }

  if(verbose>3 && dfm->use_mf) {
    if(!niikcortex_add_color(dfm->ctx[0],dfm->mf_list,0.8,1.2,NIIK_COLORMAP_SPECTRAL,50)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_add_color\n",fcname);
      return 0;
    }

    if(!niikcortex_add_color(dfm->ctx[1],dfm->mf_list,0.8,1.2,NIIK_COLORMAP_SPECTRAL,50)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_add_color\n",fcname);
      return 0;
    }

    /* write output */
    if(verbose>3 && get_DEBUG_PREFIX()) {
      sprintf(fname,"%s_tmp_dfm_white.ply",get_DEBUG_PREFIX());
      fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
      off_kobj_write_offply(fname,dfm->ctx[0],0);
      sprintf(fname,"%s_tmp_dfm_pial.ply",get_DEBUG_PREFIX());
      fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
      off_kobj_write_offply(fname,dfm->ctx[1],0);
    }
    off_kobj_add_one_color(dfm->ctx[0],1,1,0); /* yellow for white matter-surface */
    off_kobj_add_one_color(dfm->ctx[1],1,0,0); /* red for pial-surface */
  }
  /* end of multiplication factor calculation */

  if(verbose>0) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i cth %8.6f +/- %8.6f  min/25%%/med/75%%/max %6.4f %6.4f %6.4f %6.4f %6.4f\n",
                          iter,
                          niik_get_mean_from_double_vector(dfm->thklist,dfm->nvert),
                          niik_get_stdv_from_double_vector(dfm->thklist,dfm->nvert),
                          niik_get_min_from_double_vector(dfm->thklist,dfm->nvert),
                          niik_get_percentile_from_double_vector(dfm->thklist,dfm->nvert,0.25),
                          niik_median_quicksort_double_untouch(dfm->thklist,dfm->nvert),
                          niik_get_percentile_from_double_vector(dfm->thklist,dfm->nvert,0.75),
                          niik_get_max_from_double_vector(dfm->thklist,dfm->nvert));


  /* check surface intersection before deformation */
  xsc = niikcortex_off_count_intersection(bb,dfm->ctx[0], dfm->ctx[1]);
  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);
  /*TODO: fix?*/

  if(dfm->convergence_log)
    fprintf(dfm->convergence_log,"iteration,wtime,mean_ics,std_ics,min_ics,max_ics,mean_ocs,std_ocs,min_ocs,max_ocs\n");

  /************************************************
   *
   * MAIN LOOP STARTS HERE
   *
   ************************************************/
  ref_time=time(NULL);
  if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] main loop\n");
  for(iter=0,dmax_dfm=0.0,remesh_fail=0; iter<dfm->iter; iter++,dmax_dfm+=(dfm->apply_step * dfm->iter2)) {
    /*
     * prepare for deformation
     * -normal, cortical thickness, and bounding box
     */

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    fprintf(stdout,"[%s] iter %-4i of %-4i %s\n",fcname,iter+1,dfm->iter,tmstr);
    if(verbose>2) fprintf(stdout,"[%s] iter %-4i update normal\n",fcname,iter+1);

    #pragma omp parallel for
    for(cidx=0; cidx<2; cidx++) {
      off_update_kobj_face_normal(dfm->ctx[cidx]);
      off_update_kobj_vert_normal(dfm->ctx[cidx]);
      off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
    }

    for(cidx=0; cidx<2; cidx++) {
      if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        fprintf(stdout,"iter %-4i cortex %i mean curvature=%8.6f\n",iter,cidx,mean_curvature(&dfm->crv[cidx],dfm->ctx[cidx]));
    }

    if(dmax_dfm>=bb->delta*0.45 ) { /* FOR TESTING, remove || 1 if debugged properly */
      if(verbose>1) fprintf(stdout,"[%s] iter %-4i bbox updating\n",fcname,iter+1);
      NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");
      dmax_dfm=0;
    }

    if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i update thickness\n",iter+1);

#if _OPENMP>=201307
    #pragma omp simd
#endif
    for(vidx=0; vidx<dfm->nvert; vidx++) {
      dfm->thklist[vidx] = niikpt_distance(dfm->vmat[0][vidx]->v,dfm->vmat[1][vidx]->v);
    }

    /*TODO: make this parameter (5 iterations)*/
    if(!(iter+1)%5 && dfm->use_mf) {
      if(!niikcortex_deform_cortex_calc_multiplication_factor(dfm->t1img,dfm->brain_mask,dfm->ctx,dfm->tissue_val[1],dfm->mf_lim[0],dfm->mf_list)) {
        fprintf(stderr,"[%s] ERROR: niikcortex_deform_cortex_calc_multiplication_factor\n",fcname);
        return 0;
      }
    }



    /******************************************************
     * calculate deformation
     ******************************************************/

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>1) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i calculate deformation %s\n",iter+1,tmstr);

    dfm->bb=bb;
    dfm->cortex_id=cortex_id;

    NIIK_RET0((!niikcortex_deform_calc_deformation(dfm)),
              fcname,"niikcortex_deform_calc_deformation");
    
    if(dfm->convergence_log)
    {
      double dmin_ics=dfm->dfmlist[0][0].w;
      double dmax_ics=dfm->dfmlist[0][0].w;
      double dmean_ics=0.0,dstdv_ics=0.0;

      double dmin_ocs=dfm->dfmlist[1][0].w;
      double dmax_ocs=dfm->dfmlist[1][0].w;
      double dmean_ocs=0.0, dstdv_ocs=0.0;

      for(vidx=0; vidx<dfm->nvert; vidx++) {
          double d0=dfm->dfmlist[0][vidx].w;
          double d1=dfm->dfmlist[1][vidx].w;

          if     (dmin_ics>d0) dmin_ics=d0;
          else if(dmax_ics<d0) dmax_ics=d0;
          dmean_ics += d0;
          dstdv_ics += NIIK_SQ(d0);

          if     (dmin_ocs>d1) dmin_ocs=d1;
          else if(dmax_ocs<d1) dmax_ocs=d1;
          dmean_ocs += d1;
          dstdv_ocs += NIIK_SQ(d1);
      }
      dmean_ics /= dfm->nvert;
      dstdv_ics = sqrt(dstdv_ics / dfm->nvert - dmean_ics * dmean_ics);
      
      dmean_ocs /= dfm->nvert;
      dstdv_ocs = sqrt(dstdv_ocs / dfm->nvert - dmean_ocs * dmean_ocs);

      fprintf(dfm->convergence_log,"%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",iter, (int)(ctm-ref_time),
              dmean_ics,dstdv_ics,dmin_ics,dmax_ics, dmean_ocs,dstdv_ocs,dmin_ocs,dmax_ocs);
    }

    if(verbose>1) fprintf(stdout,"[%s] iter %-4i calculate deformation done %s\n",fcname,iter+1,tmstr);

    if(verbose>1 && g_niikcortex_deform_calc_deformation_vertex_index>=0) {
      vidx = g_niikcortex_deform_calc_deformation_vertex_index % dfm->ctx[0]->nvert;
      fprintf(stdout,"  vertex   %9.4f %9.4f %9.4f |    %9.4f %9.4f %9.4f\n",
              dfm->vmat[0][vidx]->v.x,dfm->vmat[0][vidx]->v.y,dfm->vmat[0][vidx]->v.z,
              dfm->vmat[1][vidx]->v.x,dfm->vmat[1][vidx]->v.y,dfm->vmat[1][vidx]->v.z);
      fprintf(stdout,"  normal   %9.4f %9.4f %9.4f |    %9.4f %9.4f %9.4f\n",
              dfm->vmat[0][vidx]->normal.x,dfm->vmat[0][vidx]->normal.y,dfm->vmat[0][vidx]->normal.z,
              dfm->vmat[1][vidx]->normal.x,dfm->vmat[1][vidx]->normal.y,dfm->vmat[1][vidx]->normal.z);
      if(dfm->ctx_label[g_niikcortex_deform_calc_deformation_vertex_index]==0) {
        fprintf(stdout,"  deform   non-cortex\n");
      } else {
        fprintf(stdout,"  deform   %9.4f %9.4f %9.4f |    %9.4f %9.4f %9.4f\n",
                dfm->dfmlist[0][vidx].x,dfm->dfmlist[0][vidx].y,dfm->dfmlist[0][vidx].z,
                dfm->dfmlist[1][vidx].x,dfm->dfmlist[1][vidx].y,dfm->dfmlist[1][vidx].z); /*TODO: update to use .w for magnitude*/
      }
    } /* show vertex motion */

    /* show deformation stats */
    if(verbose>1) {
      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && cortex_id==2) continue;
        if(cidx==1 && cortex_id==1) continue;
        dmin=dmax=dfm->dfmlist[cidx][0].w;
        dmean=dstdv=0;
        for(vidx=nlo=nhi=0; vidx<dfm->nvert; vidx++) {
          if     (dmin>dfm->dfmlist[cidx][vidx].w) {
            dmin=dfm->dfmlist[cidx][vidx].w;
            nlo=vidx;
          } else if(dmax<dfm->dfmlist[cidx][vidx].w) {
            dmax=dfm->dfmlist[cidx][vidx].w;
            nhi=vidx;
          }
          dmean += dfm->dfmlist[cidx][vidx].w;
          dstdv += dfm->dfmlist[cidx][vidx].w * dfm->dfmlist[cidx][vidx].w;
        }
        dmean /= dfm->nvert;
        dstdv = sqrt(dstdv / dfm->nvert - dmean * dmean);
        if(!cidx)
          fprintf(stdout,"[%s]      estimated deformation  [white]: %9.5f +/- %9.5f (%9.5f, %6.3f) vidx: %6i %6i\n",fcname,dmean,dstdv,dmin,dmax,nlo,nhi);
        else
          fprintf(stdout,"[%s]      estimated deformation   [pial]: %9.5f +/- %9.5f (%9.5f, %6.3f) vidx: %6i %6i\n",fcname,dmean,dstdv,dmin,dmax,nlo,nhi);
      }
    }

    /******
     * Add random noise (HACK)
     */
    #if 0
    NIIK_RET0((!niikcortex_deform_randomize_deformation(dfm,dfmlist,iter)),
              fcname,"niikcortex_deform_randomize_deformation");
    #endif 

    /*************************************
     * apply deformation
     *************************************/
    strftime(tmstr,256,"%Y-%m-%d %T",localtime(&ctm));
    if(verbose>0) {
      fprintf(stdout,"[%s] iter %-4i apply deformation %s\n",fcname,iter+1,tmstr);
    }
    dratio[0]=dratio[1]=0;
    for(cidx=0; cidx<2; cidx++) {
      niik_set_zero_for_double_vector(dfm->dd[cidx],dfm->nvert);
    }

    NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");

    if(debug_tracing)
    {
      /*VF: maybe do it more effeciently?*/
      NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[0],1,1,0)),fcname,"off_kobj_add_one_color"); /* yellow for white matter-surface */
      NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[1],1,0,0)),fcname,"off_kobj_add_one_color"); /* red for pial-surface */

      falcon_tracing_dump(&trace,iter,"refine",dfm->t1img,bb);
      falcon_tracing_dump_objects(&trace, iter,"refine", dfm->ctx, 2 );
    }

    for(iter2=0,dist=dfm->apply_step; iter2<dfm->iter2; iter2++) {

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && cortex_id==2) continue;
        if(cidx==1 && cortex_id==1) continue;

        #ifndef _USE_PRL_DEFORM
        NIIK_RET0((!niikcortex_deform_apply_deformation(bb,dfm->ctx[cidx],dist,dfm->dfmlist[cidx],dfm->dd[cidx])),fcname,"niikcortex_deform_apply_deformation ");
        #else 
        NIIK_RET0((!niikcortex_deform_apply_deformation_prl(bb,dfm->ctx[cidx],dist,dfm->dfmlist[cidx],dfm->dd[cidx],dfm->indicator,cidx)),
                  fcname,"niikcortex_deform_apply_deformation ");
        #endif
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && cortex_id==2) continue;
        if(cidx==1 && cortex_id==1) continue;
        dmin=dmax=dfm->dd[cidx][0];
        dmean=dstdv=0;

        for(vidx=0; vidx<dfm->nvert; vidx++) {
          dfm->rd[cidx][vidx] = dfm->dfmlist[cidx][vidx].w; /* deformation distance calculated but not actually applied */
          if     (dmin>dfm->dd[cidx][vidx]) dmin=dfm->dd[cidx][vidx];
          else if(dmax<dfm->dd[cidx][vidx]) dmax=dfm->dd[cidx][vidx];
          dmean += dfm->dd[cidx][vidx];
          dstdv += NIIK_SQ(dfm->dd[cidx][vidx]);
        }
        dmean /= dfm->nvert;
        dstdv = sqrt(dstdv / dfm->nvert - dmean * dmean);

        if(verbose>1 && iter2==dfm->iter2-1) { /* show deformation stats */
          if(cidx==0) {
            fprintf(stdout,"[%s]   actual deformation [%i white]: %9.5f +/- %9.5f (%9.5f, %6.4f)\n",
                    fcname,cidx,dmean,dstdv,dmin,dmax);
          }
          if(cidx==1) {
            fprintf(stdout,"[%s]   actual deformation [%i  pial]: %9.5f +/- %9.5f (%9.5f, %6.4f)\n",
                    fcname,cidx,dmean,dstdv,dmin,dmax);
          }
          if(verbose>1) {
            fprintf(stdout,"[%s]          remaining deformation: %9.5f +/- %9.5f (%9.5f, %9.5f) Q=(%9.5f, %9.5f)\n",fcname,
                    niik_get_mean_from_double_vector(dfm->rd[cidx],dfm->nvert),
                    niik_get_stdv_from_double_vector(dfm->rd[cidx],dfm->nvert),
                    niik_get_min_from_double_vector(dfm->rd[cidx],dfm->nvert),
                    niik_get_max_from_double_vector(dfm->rd[cidx],dfm->nvert),
                    niik_get_percentile_from_double_vector(dfm->rd[cidx],dfm->nvert,0.25),
                    niik_get_percentile_from_double_vector(dfm->rd[cidx],dfm->nvert,0.75));
          }  /* verbose */
        } /* verbose */

      } /* each surface */
    } /* dist */

    for(cidx=ci[0]=ci[1]=0; cidx<2; cidx++) {
      if(cidx==0 && cortex_id==2) continue;
      if(cidx==1 && cortex_id==1) continue;
      for(vidx=0; vidx<dfm->nvert; vidx++) {
        if(dfm->dd[cidx][vidx]<g_niikcortex_min_deform_distance)
          ci[cidx]++;
      }
    }
    fprintf(stdout,"[%s] deform below threshold: %9.5f %9.5f %%\n",fcname, 100.0*ci[0]/dfm->nvert, 100.0*ci[1]/dfm->nvert);

    /*TODO: make this a parameter, also target elen*/
    /*Potentially this can make self-intersection */
    /*WARNING: this doesn't seem to work properly yet*/
    /*TODO: add curvature constraints ?*/
    if(dfm->remesh>0 && iter%dfm->remesh==(dfm->remesh-1) && remesh_fail<dfm->remesh_fail_max) {  
      int xsc;
      int r;
      kobj *backup[2];
      /*make a backup surfaces*/
      backup[0]=off_kobj_copy(dfm->ctx[0]);
      backup[1]=off_kobj_copy(dfm->ctx[1]);

      /*need curvature info*/
      for(cidx=0;cidx<2;cidx++)
      {
          /*if(dfm->weight->m[0][WEIGHT_CURVATURE]<=0.0 || dfm->weight->m[1][WEIGHT_CURVATURE]<=0.0)*/
          off_update_kobj_face_normal(dfm->ctx[cidx]);
          off_update_kobj_vert_normal(dfm->ctx[cidx]);
          off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
          update_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
      }

      fprintf(stdout,"[%s] Remeshing\n",fcname);

      for(r=0;r<5;r++)
      {
        off_remesh_dual_kobj(dfm->ctx, dfm->crv, 1 ); /*try not to do any smoothing of the surface*/
        /*split into independent stages: splitting, collapsing, valence optimization*/
        /*need to re-calculate everything...*/
        NIIK_RET0(!niicortex_deform_prepare(dfm),fcname,"niicortex_deform_prepare");
        for(cidx=0;cidx<2;cidx++) {
          off_update_kobj_face_normal(dfm->ctx[cidx]);
          off_update_kobj_vert_normal(dfm->ctx[cidx]);
          off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
          update_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
        }

        NIIK_RET0((!niikcortex_off_correct_self_intersections(bb, dfm->ctx, 25)),fcname,"niikcortex_off_correct_self_intersections");

        NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");
        xsc = niikcortex_off_count_intersection(bb, dfm->ctx[0], dfm->ctx[1]);
        fprintf(stdout,"[%s] after remesh surface intersection %i\n",fcname,xsc);
        if(!xsc) break;
      }
      if(xsc && get_POSTMORTEM_PREFIX())
      {
        sprintf(fname,"%s_failed_selfintersect_dfm%03i_ics.ply",get_POSTMORTEM_PREFIX(),iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[0],0);

        sprintf(fname,"%s_failed_selfintersect_dfm%03i_ocs.ply",get_POSTMORTEM_PREFIX(),iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[1],0);
      }
      if(xsc) {
        /*rollback!*/
        fprintf(stdout,"[%s] rolling back \n",__func__);
        off_kobj_free(dfm->ctx[0]);
        off_kobj_free(dfm->ctx[1]);
        dfm->ctx[0]=backup[0];
        dfm->ctx[1]=backup[1];
        remesh_fail++;
      } else {
        off_kobj_free(backup[0]);
        off_kobj_free(backup[1]);
      }

      NIIK_RET0(!niicortex_deform_prepare(dfm),fcname,"niicortex_deform_prepare");
      /*need curvature info*/
      for(cidx=0;cidx<2;cidx++)
      {
          /*if(dfm->weight->m[0][WEIGHT_CURVATURE]<=0.0 || dfm->weight->m[1][WEIGHT_CURVATURE]<=0.0)*/
          off_update_kobj_face_normal(dfm->ctx[cidx]);
          off_update_kobj_vert_normal(dfm->ctx[cidx]);
          off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
          update_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
      }
      NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");

      dmax_dfm=0.0;
    }

    /*TODO: use tol to stop iterations here maybe?*/

    /******************************************************
     * update cortical thickness and show
     ******************************************************/
    for(vidx=0; vidx<dfm->nvert; vidx++) {
      dfm->thklist[vidx] = niikpt_distance(dfm->vmat[0][vidx]->v,dfm->vmat[1][vidx]->v);
    }
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>1) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i cth %10.6f +/- %10.6f   min/25%%/med/75%%/max %6.4f %6.4f %6.4f %6.4f %6.4f %s\n",
                            iter+1,
                            niik_get_mean_from_double_vector(dfm->thklist,dfm->nvert),
                            niik_get_stdv_from_double_vector(dfm->thklist,dfm->nvert),
                            niik_get_min_from_double_vector(dfm->thklist,dfm->nvert),
                            niik_get_percentile_from_double_vector(dfm->thklist,dfm->nvert,0.25),
                            niik_median_quicksort_double_untouch(dfm->thklist,dfm->nvert),
                            niik_get_percentile_from_double_vector(dfm->thklist,dfm->nvert,0.75),
                            niik_get_max_from_double_vector(dfm->thklist,dfm->nvert),
                            tmstr);

    /*****************************************************
     * write temporary output
     *****************************************************/
    if(verbose>3 && (iter%5)==0 && get_DEBUG_PREFIX()) {
      if(cortex_id%2) {
        NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[0],1,1,0)),fcname,"off_kobj_add_one_color ics");
        sprintf(fname,"%s_tmp_dfm%03i_ics.ply",get_DEBUG_PREFIX(),iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[0],0);
      }
      if(cortex_id>=2) {
        NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[1],1,0,0)),fcname,"off_kobj_add_one_color ocs");
        sprintf(fname,"%s_tmp_dfm%03i_ocs.ply",get_DEBUG_PREFIX(),iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[1],0);
      }
    } /* writing temp output */

    /******************************************************
     * check for the stop condition
     * -to be implemented properly
     ******************************************************/

    /*xsc = niikcortex_off_count_intersection(bb, white_surface, pial_surface);
    fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);*/

  } /* iteration */

  if(dfm->debug_trace_log) fclose(dfm->debug_trace_log);

  xsc = niikcortex_off_count_intersection(bb, dfm->ctx[0], dfm->ctx[1]);
  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);

  off_bbox_free(bb);

  for(n=0; n<4; n++) {
    int i;
    for(i=0; i<4; i++)
      grad_img[n][i] = niik_image_free(grad_img[n][i]);
    free(grad_img[n]);
  }

  for(n=0; n<3; n++) {
    int i;
    if(grad_prior[n]) {
        for(i=0; i<4; i++)
          grad_prior[n][i] = niik_image_free(grad_prior[n][i]);
        free(grad_prior[n]);
      }
  }

  for(n=0; n<2; n++) {
    div_img[n]=niik_image_free(div_img[n]);
  }
  invmat = niikmat_free(invmat);
  falcon_tracing_free(&trace);

  if(verbose>1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_cortex */


niikcortex_deform *niikcortex_deform_init() {
  const char *fcname=__func__;
  niikcortex_deform *dfm=NULL;
  int n;
  if((dfm=(niikcortex_deform *)calloc(1,sizeof(niikcortex_deform)))==NULL) {
    fprintf(stderr,"[%s] ERROR: could not allocate memory for dfm\n",fcname);
    return NULL;
  }

  dfm->t1img =
    dfm->nonctx_mask =
      dfm->avoid_mask =
        dfm->gwi_mask =
          dfm->csf_mask =
            dfm->brain_mask =
              dfm->lesion_mask =
                NULL;

  NIIK_RET0(((dfm->ctx = (kobj **)calloc(2,sizeof(kobj *)))==NULL),fcname,"could not allocate memory for ctx");
  dfm->ctx[0]=dfm->ctx[1]=NULL;/*not needed because calloc*/
  NIIK_RET0(((dfm->prior = (nifti_image **)calloc(4,sizeof(nifti_image *)))==NULL),fcname,"could not allocate memory for ctx");
  dfm->prior[0]=dfm->prior[1]=dfm->prior[2]=dfm->prior[3]=NULL; /*not needed because calloc*/
  dfm->delta = 0.5;
  dfm->apply_step = 0.2;
  dfm->weight=niikmat_init(2,WEIGHT_COUNT);
  dfm->weight->m[CORTEX_ICS][ WEIGHT_SURFACE]     = 1.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_SURFACE]     = 1.0;  /* surface weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_THICKNESS_SMOOTHNESS] = 0.5;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_THICKNESS_SMOOTHNESS] = 0.5;  /* thickness smoothness weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_IMAGE]       = 1.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_IMAGE]       = 1.0;  /* image weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_BRAIN_MASK ] = 1.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_BRAIN_MASK ] = 1.0;  /* brain mask weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_VENTRICLE ]  = 1.5;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_VENTRICLE ]  = 1.5;  /* ventricle mask weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_LESION ]     = 1.5;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_LESION ]     = 1.5;  /* lesion mask weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_PROXIMITY ]  = 0.2;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_PROXIMITY ]  = 0.2;  /* proximity weights */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_AVOID ] = 1.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_AVOID ] = 1.0;  /* avoid mask */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_ABS_THICKINESS ] = 0.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_ABS_THICKINESS ] = 0.0;  /* absolute thickness */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_GRADIENT ]   = 0.1;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_GRADIENT ]   = 0.1;  /* gradient term */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_DGRADIENT ]  = 0.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_DGRADIENT ]  = 0.0;  /* divergence term */
  dfm->weight->m[CORTEX_ICS][ WEIGHT_CURVATURE ] = 0.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_CURVATURE ] = 0.0;     /* curvature term */

  /*Update smoothing terms, for now here*/
  dfm->weight->m[CORTEX_ICS][ WEIGHT_UPDATE_SIGMA ] = 0.0;
  dfm->weight->m[CORTEX_OCS][ WEIGHT_UPDATE_SIGMA ] = 0.0;  /* field update smoothing */

  /*Thinkness constraints, for now I am putting them here*/
  dfm->weight->m[CORTEX_ICS][ WEIGHT_MIN_THICKNESS ] = 0.0; 
  dfm->weight->m[CORTEX_OCS][ WEIGHT_MIN_THICKNESS ] = 0.0;  
  dfm->weight->m[CORTEX_ICS][ WEIGHT_MAX_THICKNESS ] = 3.0; 
  dfm->weight->m[CORTEX_OCS][ WEIGHT_MAX_THICKNESS ] = 3.0; 

  dfm->weight->m[CORTEX_ICS][ WEIGHT_THICKNESS_SIGMA ] = 1.0;  
  dfm->weight->m[CORTEX_OCS][ WEIGHT_THICKNESS_SIGMA ] = 1.0;  

  dfm->weight->m[CORTEX_ICS][ WEIGHT_PRIOR ] = 0.0; /* priors term */
  dfm->weight->m[CORTEX_OCS][ WEIGHT_PRIOR ] = 0.0;

  dfm->weight->m[CORTEX_ICS][ WEIGHT_MIX ] = 0.5; /* equal mix for ICS and OCS for tissue types */
  dfm->weight->m[CORTEX_OCS][ WEIGHT_MIX ] = 0.5;

  dfm->ctx_label=NULL;
  dfm->iter=100;
  dfm->iter2=5;
  dfm->remesh=0;
  dfm->remesh_fail_max=4; /*hackish way (?) maybe use different strategy?*/
  dfm->use_mf=1; 

  dfm->tolerance=0.1; /*TODO: make use of this*/
  dfm->gradient_FWHM[0]=dfm->gradient_FWHM[1]=1.0;
  dfm->divergence_FWHM[0]=dfm->divergence_FWHM[1]=1.0;
  dfm->prior_FWHM[0]=dfm->prior_FWHM[1]=1.0;
  dfm->gradient_lambda[0]=dfm->gradient_lambda[1]=-1;   /* negative -> auto */

  NIIK_RET0n(((dfm->tissue_val=(double *)calloc(12,sizeof(double)))==NULL),fcname,"calloc for tissue_val");
  for(n=0; n<12; n++)
    dfm->tissue_val[n]=-1;

  dfm->bbox_depth=7; /*VF: this seemm to need tweaking...*/
  dfm->cortex_id=3; /* 1=white, 2=pial, 3=both */
  dfm->numerical_method = NIIK_NUM_METHOD_FORWARD_EULER;
  dfm->proximity_min_distance = 0.6;
  dfm->debug_keep_tmp = 0;
  dfm->mflim=0.1; /*TODO: check ???*/
  NIIK_RET0n(((dfm->dfmlist = (niikpt **)calloc(2,sizeof(niikpt *)))==NULL),fcname,"calloc for dfmlist");
  NIIK_RET0n(((dfm->vmat = (kvert ***)calloc(2,sizeof(kvert **)))==NULL),fcname,"calloc for vmat");

  dfm->dfm_limit = NULL;
  dfm->travel = NULL;

  return dfm;
} /* niikcortex_deform_init */

niikcortex_deform *niikcortex_deform_free(niikcortex_deform *dfm) {
  int n,cidx;
  if(dfm==NULL) return NULL;

  dfm->t1img          = niik_image_free(dfm->t1img);
  dfm->nonctx_mask    = niik_image_free(dfm->nonctx_mask);
  dfm->avoid_mask     = niik_image_free(dfm->avoid_mask);
  dfm->gwi_mask       = niik_image_free(dfm->gwi_mask);
  dfm->csf_mask       = niik_image_free(dfm->csf_mask);
  dfm->brain_mask     = niik_image_free(dfm->brain_mask);
  dfm->lesion_mask    = niik_image_free(dfm->lesion_mask);
  dfm->weight         = niikmat_free( dfm->weight );

  for(cidx=0; cidx<2; cidx++) {
    free(dfm->vmat[cidx]);
    free(dfm->dfmlist[cidx]);
    free_curvature(&dfm->crv[cidx]);
    free(dfm->rd[cidx]);
    free(dfm->dd[cidx]);
  }

  free(dfm->indicator);
  free(dfm->dfmlist);
  free(dfm->vmat);
  free(dfm->thklist);
  free(dfm->mf_list);

  if(dfm->dfm_limit)
    free(dfm->dfm_limit);

  if(dfm->travel) /*do it for two surfaces*/
    free(dfm->travel);

  if(dfm->ctx!=NULL) {
    for(n=0; n<2; n++) 
      dfm->ctx[n]=off_kobj_free(dfm->ctx[n]);
  }
  if(dfm->ctx_label) {
    free(dfm->ctx_label);
    dfm->ctx_label=NULL;
  }
  free(dfm->tissue_val);
  free(dfm->ctx);

  for(n=0;n<4;n++) {
    if(dfm->prior[n]!=NULL) 
      dfm->prior[n] = niik_image_free(dfm->prior[n]);
  }
  free(dfm->prior);

  free(dfm);
  return NULL;
} /* niikcortex_deform_free */



int niikcortex_remesh_cortex(niikcortex_deform * dfm)
{
  int xsc,cidx,i;
  const char* fcname=__func__;
  bbox *bb;
  
  fprintf(stdout,"[%s] Remeshing\n",fcname);

  for(i=0;i<10;i++) {
    /*need curvature info*/
    #if 1
    for(cidx=0; cidx<2; cidx++)
    {
        re_init_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
        /*if(dfm->weight->m[0][WEIGHT_CURVATURE]<=0.0 || dfm->weight->m[1][WEIGHT_CURVATURE]<=0.0)*/
        off_update_kobj_face_normal(dfm->ctx[cidx]);
        off_update_kobj_vert_normal(dfm->ctx[cidx]);
        off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
        update_curvature(&dfm->crv[cidx], dfm->ctx[cidx]);
    }
    #endif
    off_remesh_dual_kobj(dfm->ctx,dfm->crv, 1 ); /*try not to do any smoothing of the surface*/

    bb=off_bbox_init(dfm->bbox_depth,320);

    NIIK_RET0((!niikcortex_off_correct_self_intersections(bb, dfm->ctx, 20)),fcname,"niikcortex_off_correct_self_intersections");

    NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");

    xsc = niikcortex_off_count_intersection(bb, dfm->ctx[0], dfm->ctx[1]);
    fprintf(stdout,"[%s] after remesh surface intersection %i\n",fcname,xsc);
  }
  off_bbox_free(bb);
  return 1;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
