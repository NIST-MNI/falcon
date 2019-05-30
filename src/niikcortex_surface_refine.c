/* FILENAME:     niikcortex_refine.c
 * DESCRIPTION:  Kunio's refine cortical surfaces to fit the data
 * AUTHOR:       Kunio Nakamura
 * DATE:         March 24, 2013
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

static const double g_niikmesh_min_deform_distance = 1e-4;
static const double g_niikmesh_min_deform_distance2 = 1e-3;
static int g_niikmesh_deform_calc_deformation_vertex_index=-1; /*DEBUG*/


/*defined in falcon_cortex_deformation.c*/
int niikcortex_deform_apply_deformation_prl(bbox *bb,kobj *obj,double step,niikpt *dfmlist,double *actdfm,indicator_t *indicator,int vfilter);
int niikcortex_deform_apply_deformation(bbox *bb, kobj *obj, double step, niikpt *dfmlist, double *actdfm);


int niikmesh_estimate_tissue_values(nifti_image *img,
               nifti_image * roi_mask,
               nifti_image * obj_mask,
               double *obj,
               double *back,
               double *range_obj,
               double *range_back)
{
  nifti_image *bck_mask=NULL;
  double dmin,dmax;

  dmax = niik_image_get_max(img, roi_mask);
  dmin = niik_image_get_min(img, roi_mask);


  NIIK_RET0(((bck_mask = niik_image_copy_as_type(roi_mask, NIFTI_TYPE_UINT8))==NULL),__func__,"niik_image_copy_as_type for roi_mask -> tmpimg");
  NIIK_RET0((!niik_image_maskout(bck_mask,obj_mask)),fcname,"niik_image_maskout bck_mask,obj_mask");
  /*subtract mask*/

  /*get mode of object*/
  *obj = niik_image_get_mode(img,obj_mask,dmin,dmax,128,3);

  /*get mode of background around object*/
  *back= niik_image_get_mode(img,bck_mask,dmin,dmax,128,3);

  /*TODO: find better idea?*/
  *range_obj = *range_back = fabs(*obj-*back)/2.0;

  bck_mask=niik_image_free(bck_mask);
  return 1;
}

/*******************************************************************/

/********************************************************
 *
 * calculation of deformation for each vertex
 *
 ********************************************************/
niikpt niikmesh_deform_calc_deformation_vertex_surface_smoothness(kvert *v) {
  return niikpt_sub(niikpt_kvert_local_average_with_tangential_relaxation2(v,0.999),v->v);
}

niikpt niikmesh_deform_calc_deformation_vertex_surface_smoothness_simple(kvert *v) {
  return niikpt_sub(niikpt_kvert_simple_local_avg(v),v->v);
}

/**
*/
niikpt niikmesh_deform_calc_deformation_prior_term(kvert *v, nifti_image **prior, nifti_image ***grad_prior)
{
    double prior_obj = niik_image_interpolate_3d_linear(     prior[0], v->v);
    double prior_bck = 0.0;
    double project;

    niikpt grad_prior_obj;
    niikpt grad_prior_bck;

    grad_prior_obj.x = niik_image_interpolate_3d_linear(grad_prior[0][1], v->v);
    grad_prior_obj.y = niik_image_interpolate_3d_linear(grad_prior[0][2], v->v);
    grad_prior_obj.z = niik_image_interpolate_3d_linear(grad_prior[0][3], v->v);

    if(prior[1]){
        prior_bck = niik_image_interpolate_3d_linear(     prior[1], v->v);

        grad_prior_bck.x = niik_image_interpolate_3d_linear(grad_prior[1][1], v->v);
        grad_prior_bck.y = niik_image_interpolate_3d_linear(grad_prior[1][2], v->v);
        grad_prior_bck.z = niik_image_interpolate_3d_linear(grad_prior[1][3], v->v);
    } else {
        prior_bck = 1.0 - prior_obj;

        grad_prior_bck.x = -grad_prior_obj.x;
        grad_prior_bck.y = -grad_prior_obj.y;
        grad_prior_bck.z = -grad_prior_obj.z;
    }
    grad_prior_obj   = niikpt_unit(grad_prior_obj);  /*only interested in the direction*/
    grad_prior_bck   = niikpt_unit(grad_prior_bck);  /*only interested in the direction*/

    /*TODO: normalize gradients ?*/

    /*should balance at 0.5 p value between background and object*/
    project = ( niikpt_dot(v->normal, grad_prior_obj) * prior_obj +
                niikpt_dot(v->normal, grad_prior_bck) * prior_bck ) / 2.0;

    return niikpt_kmul(v->normal, -1 * project);
}

/* VF: Balooning to stay on the edge of the tissue "class",
*/
niikpt niikmesh_deform_calc_deformation_vertex_image_term(
  kvert *v, nifti_image *img, nifti_image ***grad_img,
  double *gradient_lambda,
  double mean_ics_value,
  double range_ics_value,
  double *gray,
  double *grad)
{
  niikpt local_grad;
  double project;

  *gray = niik_image_interpolate_3d_linear(     img, v->v);

  /*move in the right direction*/
  local_grad.x = niik_image_interpolate_3d_linear(grad_img[0][1], v->v);
  local_grad.y = niik_image_interpolate_3d_linear(grad_img[0][2], v->v);
  local_grad.z = niik_image_interpolate_3d_linear(grad_img[0][3], v->v);

  local_grad   = niikpt_unit(local_grad); /*only interested in the direction*/

  project      = niikpt_dot(v->normal, local_grad);

  /*we want to move in the direction of the gradient that would bring the *gray closer to the mean_XXX_value*/
  *grad = niik_image_interpolate_3d_linear(grad_img[0][0], v->v);

  return niikpt_kmul(v->normal, -1.0 * project *
                      gradient_lambda[0] / (*grad+gradient_lambda[0]) *
                      NIIK_Heaviside11(*gray-mean_ics_value, range_ics_value));
} /* niikmesh_deform_calc_deformation_vertex_image_term */


niikpt niikmesh_deform_calc_deformation_vertex_gradient_term(kvert *v, nifti_image *img, nifti_image **grad_img, niikpt *normal) {
  normal->x = niik_image_interpolate_3d_linear(grad_img[1], v->v);
  normal->y = niik_image_interpolate_3d_linear(grad_img[2], v->v);
  normal->z = niik_image_interpolate_3d_linear(grad_img[3], v->v);

  *normal = niikpt_unit(*normal);

  /* ICS goes towards WM (positive sign)*/
  return niikpt_kmul(v->normal, niikpt_dot(v->normal, *normal)); /**/

} /* niikmesh_deform_calc_deformation_vertex_gradient_term */


niikpt niikmesh_deform_calc_deformation_vertex_flux_term(kvert *v, nifti_image *img, nifti_image *div_img,int sign) {
  double d = niik_image_interpolate_3d_linear(div_img, v->v);

  /*clamp between -1 and 1 */
  if(d<-1.0)     d=-1.0;
  else if(d>1.0) d= 1.0;

  return niikpt_kmul(v->normal, d*sign); 
} /* niikmesh_deform_calc_deformation_vertex_flux_term */



int niikmesh_deform_calc_deformation_vertex_proximity_term_2sides_test(kvert *v, kface *f, double *min_dist)
/* testing for proximity term
 * min_dist[0] - below (negative vertical distance)
 * min_dist[1] - above (positive vertical distance)
 * -updates the min_dist and vdist (signed vertical distance)
 */
{
  const int verbose=0;
  niikpt vec[3],pt;
  double t;

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
} /* niikmesh_deform_calc_deformation_vertex_proximity_term_2sides_test */


int niikmesh_deform_calc_deformation_vertex_proximity_term_2sides(kvert *v, double *prox_min_dist, bbox *bb, niikpt *pprox)
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
  double
    eps=0.001,
    dmin[2],vdis,min_dist,
    xmin,ymin,zmin,xmax,ymax,zmax,dval;

  int i,j,k,m,n,idx;
  kface *f;
  niikpt vmin,vmax;

  int verbose=1;

  if(verbose>0) {
    if(v->index-1!=g_niikmesh_deform_calc_deformation_vertex_index) {
      verbose=0;
    }
  }

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

          if(niikmesh_deform_calc_deformation_vertex_proximity_term_2sides_test(v,f,dmin)) {
            if(verbose>0) {
              fprintf(stdout,"[%s] dist = %9.5f %9.5f  \n",__func__,dmin[0],dmin[1]);
            }
          } /* new min */
        } /* each face */
      }
    }
  } /* each bbox */

  /*balance between two*/
  /*TODO: clamp between -1 and 1*/
  dval = min_dist/(-dmin[0] + eps) - min_dist/(dmin[1] + eps);

  if     (dval < -1.0) dval=-1.0;
  else if(dval >  1.0) dval=1.0;

  *pprox = niikpt_kmul(v->normal, dval);

  prox_min_dist[0] = -dmin[0];
  prox_min_dist[1] = dmin[1];

  return 1;
} /* niikmesh_deform_calc_deformation_vertex_proximity_term_2sides */


int niikmesh_deform_calc_deformation_vertex_curvature_term(kvert *v, off_curvature_t *crv, niikpt *pcurv) {
  double curv1=crv->curv1[v->index-1];
  double curv2=crv->curv2[v->index-1];

  /*Using mean curvature right now*/
  *pcurv=niikpt_kmul(v->normal, -1.0*(curv1+curv2)/2.0 );
  return 1;
}

niikpt mesh_deform_calc_deformation_vertex(refine_info  *dfm, kvert* v, int vidx)
{
  /*Calculates deformation force for one vertex*/
  niikpt
      normal,/*Vertex normal*/
      psum,  /*sum of all forces*/
      pavg,  /*smoothing force (local average)*/
      ptsurf,
      pssurf,
      pimag,
      pgrad,
      pmask,
      pcrm,
      plesm,
      pprox,
      pcurv,
      pdgrad,
      pprior;
    
  double
        dmin2[2],
        gray,grad,
        bval=0, vval=0, lval=0, vcrm=0;

  bval = vval = lval = vcrm = 0;
  psum = pavg = ptsurf = pssurf = pimag = pmask =  plesm = pprox = pgrad = pcrm = pcurv = pdgrad = pprior = niikpt_zero();
  
  /* Do not apply certain forces if this one is in effect ?*/
  /* brain mask term , surface stays inside of the brain */
  if(dfm->deform_weights->m[0][WEIGHT_BRAIN_MASK] != 0.0) {
    bval = niik_image_interpolate_3d_linear(dfm->brain_mask,v->v);
    pmask = niikpt_kmul(v->normal, bval-1.0);
  }

  /* lesion term */
  if(dfm->lesion_mask!=NULL) {
    if(dfm->deform_weights->m[0][WEIGHT_LESION] != 0.0) {
      lval = niik_image_interpolate_3d_linear(dfm->lesion_mask,v->v);
      plesm = niikpt_kmul(v->normal, lval);
    }
  }

  /* tangential surface smoothness term */
  if(dfm->deform_weights->m[0][WEIGHT_SURFACE] != 0.0) {
    ptsurf = niikmesh_deform_calc_deformation_vertex_surface_smoothness(v);
  }

  /* simple surface smoothness term */
  if(dfm->deform_weights->m[0][WEIGHT_SURFACE_SIMPLE] != 0.0) {
    pssurf = niikmesh_deform_calc_deformation_vertex_surface_smoothness_simple(v);
  }

  /* image term */
  if(dfm->deform_weights->m[0][WEIGHT_IMAGE] != 0.0) {
    pimag = niikmesh_deform_calc_deformation_vertex_image_term(
              v, dfm->img, dfm->grad_img, dfm->gradient_lambda,
              dfm->mean_ics_value * dfm->mf_list[v->index-1],
              dfm->range_ics_value,
              &gray, &grad);
  }

  /* image gradient term */
  if(dfm->deform_weights->m[0][WEIGHT_GRADIENT] != 0.0) {
    pgrad = niikmesh_deform_calc_deformation_vertex_gradient_term(
              v, dfm->img, dfm->grad_img[0], &normal);
  }

  /* flux gradient term */
  if(dfm->deform_weights->m[0][WEIGHT_DGRADIENT] != 0.0) {
    pdgrad = niikmesh_deform_calc_deformation_vertex_flux_term(
               v, dfm->img, dfm->div_img[0], (dfm->range_ocs_value>0.0?-1:1) ) ;
  }

  /* proximity term */
  if(dfm->deform_weights->m[0][WEIGHT_PROXIMITY] != 0.0) {
    dmin2[0] = dmin2[1] = dfm->prox_min_dist;
    niikmesh_deform_calc_deformation_vertex_proximity_term_2sides(v, dmin2, dfm->bb, &pprox);
  } /* proximity deformation */

  /* curvature term */
  if(dfm->deform_weights->m[0][WEIGHT_CURVATURE] != 0.0) {
    niikmesh_deform_calc_deformation_vertex_curvature_term(v, &dfm->crv[0], &pcurv);
  } /* curvature deformation */

  /* prior term */
  if(dfm->prior[0]!=NULL && dfm->deform_weights->m[0][WEIGHT_PRIOR] != 0.0) {
    pprior=niikmesh_deform_calc_deformation_prior_term(v, dfm->prior, dfm->grad_prior);
  } /* prior deformation */

  /* combine deformation with weighting */
  psum.x =
    ptsurf.x * dfm->deform_weights->m[0][WEIGHT_SURFACE] +
    pssurf.x * dfm->deform_weights->m[0][WEIGHT_SURFACE_SIMPLE] +
    pimag.x  * dfm->deform_weights->m[0][WEIGHT_IMAGE] +
    pgrad.x  * dfm->deform_weights->m[0][WEIGHT_GRADIENT] +
    pmask.x  * dfm->deform_weights->m[0][WEIGHT_BRAIN_MASK] +
    plesm.x  * dfm->deform_weights->m[0][WEIGHT_LESION] +
    pprox.x  * dfm->deform_weights->m[0][WEIGHT_PROXIMITY] +
    pcrm.x   * dfm->deform_weights->m[0][WEIGHT_CEREBELLUM] +
    pcurv.x  * dfm->deform_weights->m[0][WEIGHT_CURVATURE] +
    pdgrad.x * dfm->deform_weights->m[0][WEIGHT_DGRADIENT] +
    pprior.x * dfm->deform_weights->m[0][WEIGHT_PRIOR];

  psum.y =
    ptsurf.y * dfm->deform_weights->m[0][WEIGHT_SURFACE] +
    pssurf.y * dfm->deform_weights->m[0][WEIGHT_SURFACE_SIMPLE] +
    pimag.y  * dfm->deform_weights->m[0][WEIGHT_IMAGE] +
    pgrad.y  * dfm->deform_weights->m[0][WEIGHT_GRADIENT] +
    pmask.y  * dfm->deform_weights->m[0][WEIGHT_BRAIN_MASK] +
    plesm.y  * dfm->deform_weights->m[0][WEIGHT_LESION] +
    pprox.y  * dfm->deform_weights->m[0][WEIGHT_PROXIMITY] +
    pcrm.y   * dfm->deform_weights->m[0][WEIGHT_CEREBELLUM]+
    pcurv.y  * dfm->deform_weights->m[0][WEIGHT_CURVATURE] +
    pdgrad.y * dfm->deform_weights->m[0][WEIGHT_DGRADIENT] +
    pprior.y * dfm->deform_weights->m[0][WEIGHT_PRIOR];

  psum.z =
    ptsurf.z * dfm->deform_weights->m[0][WEIGHT_SURFACE] +
    pssurf.z * dfm->deform_weights->m[0][WEIGHT_SURFACE_SIMPLE] +
    pimag.z  * dfm->deform_weights->m[0][WEIGHT_IMAGE] +
    pgrad.z  * dfm->deform_weights->m[0][WEIGHT_GRADIENT] +
    pmask.z  * dfm->deform_weights->m[0][WEIGHT_BRAIN_MASK] +
    plesm.z  * dfm->deform_weights->m[0][WEIGHT_LESION] +
    pprox.z  * dfm->deform_weights->m[0][WEIGHT_PROXIMITY] +
    pcrm.z   * dfm->deform_weights->m[0][WEIGHT_CEREBELLUM]+
    pcurv.z  * dfm->deform_weights->m[0][WEIGHT_CURVATURE]+
    pdgrad.z * dfm->deform_weights->m[0][WEIGHT_DGRADIENT] +
    pprior.z * dfm->deform_weights->m[0][WEIGHT_PRIOR];

  return psum;
} /* niikmesh_deform_calc_deformation_vertex */




int mesh_deform_calc_deformation( niikcortex_deform * dfm  )
/* -computes the deformation for each vertex
 * -updates dfmlist
 * -uses mesh_deform_calc_deformation_vertex
 * -method of integration is determined by num_method
 * output is in dfmlist array, with unit vectors and magnitude stored in .w component
 */
{
  int n,cidx=0,vidx;
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
  dfm_refine.gwi_mask   = dfm->gwi_mask;
  dfm_refine.crv = dfm->crv;
  dfm_refine.deform_weights = dfm->weight;
  dfm_refine.mf_list = dfm->mf_list;
  dfm_refine.div_img = dfm->div_img;
  dfm_refine.grad_img = dfm->grad_img;
  dfm_refine.prior  = dfm->prior;
  dfm_refine.grad_prior = dfm->grad_prior;
  dfm_refine.gradient_lambda = dfm->gradient_lambda;
  dfm_refine.lesion_mask = dfm->lesion_mask;
  dfm_refine.prox_min_dist = dfm->proximity_min_distance;

  dfm_refine.mean_ics_value  =(dfm->tissue_val[0]+dfm->tissue_val[1])/2;
  dfm_refine.range_ics_value = dfm->tissue_val[2];
  dfm_refine.range_ocs_value = dfm->tissue_val[0] - dfm->tissue_val[1];

  /*HACK: using max intensity */
  dfm_refine.intWM = dfm->tissue_val[0]>dfm->tissue_val[1] ? dfm->tissue_val[0]:dfm->tissue_val[1];
  dfm_refine.vmat = dfm->vmat;
  dfm_refine.dfm = dfm;

  dt = dfm->weight->m[0][WEIGHT_DT];
  ctm=time(NULL);
  if(verbose>=2) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
  }

  switch(dfm->numerical_method) {
  case NIIK_NUM_METHOD_FORWARD_EULER:
    if(verbose>0) fprintf(stdout,"[%s] Forward Euler method\n",fcname);

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      dfm->dfmlist[0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine, dfm->vmat[CORTEX_ICS][vidx], vidx);
    } /* each vertex */
    break;

  case NIIK_NUM_METHOD_BACKWARD_EULER:  /* backward euler */
    if(verbose) fprintf(stdout,"[%s] Backward Euler method\n",fcname);
    op = niikpt_matrix(1, dfm->ctx[0]->nvert);
    if(verbose>2) fprintf(stdout,"[%s]   backup vertex positions\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      op[cidx][vidx] = dfm->vmat[cidx][vidx]->v;
    }

    if(verbose>=2) fprintf(stdout,"[%s]   k1 calculation [Backward Euler]\n",fcname);

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      dfm->dfmlist[0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex (Backward Eulder) */
    /* (Backward Eulder)
     * dfm = dt Func(t,Ctx1)
     * Ctx1 = Ctx + k1
     */
    if(verbose>=2) fprintf(stdout,"[%s]   k2 preparation [Backward Euler]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        dfm->vmat[0][vidx]->v = niikpt_move_normal(op[0][vidx],dfm->dfmlist[0][vidx],dt);
    }
    off_update_kobj_face_normal(dfm->ctx[0]);
    off_update_kobj_vert_normal(dfm->ctx[0]);
    off_smooth_kobj_vert_normal(dfm->ctx[0]);
    off_update_kobj_kface_pminmax(dfm->ctx[0]);

    if(dfm->weight->m[0][WEIGHT_CURVATURE]>0.0)
      update_curvature(&dfm->crv[0], dfm->ctx[0]);

    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Backward Euler]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = niikpt_zero();
        continue;
      }
      dfm->dfmlist[0][vidx]=mesh_deform_calc_deformation_vertex(&dfm_refine, dfm->vmat[CORTEX_ICS][vidx], vidx);
    } /* each vertex */
    /* update dfmlist
     * put the original vertex back */
    if(verbose>=2) fprintf(stdout,"[%s]   output and input update [Backward Euler]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        dfm->vmat[0][vidx]->v = op[0][vidx];
    }
    /* free memory */
    if(verbose>2) fprintf(stdout,"[%s] free memory\n",fcname);

    free(op[0]);
    free(op);
    break;

  case NIIK_NUM_METHOD_MIDPOINT:  /* midpoint */
    if(verbose>0) fprintf(stdout,"[%s] Mid-point method\n",fcname);
    op = niikpt_matrix(1, dfm->ctx[0]->nvert);
    if(verbose>=2) fprintf(stdout,"[%s]   backup vertex positions [Midpoint method]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        op[0][vidx]=dfm->vmat[0][vidx]->v;
    }

    /* k1 = dCtx = dt Func(t,Ctx) */
    if(verbose>1) fprintf(stdout,"[%s]   k1 calculation [Midpoint method]\n",fcname);
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      dfm->dfmlist[0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /*
     * k2 = dCtx2 = dt Func(t,Ctx1)
     * Ctx1 = Ctx + 0.5 * k1
     */
    if(verbose>=3) fprintf(stdout,"[%s]   k2 preparation [Midpoint method]\n",fcname);


    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      dfm->vmat[0][vidx]->v = niikpt_move_normal(op[0][vidx],dfm->dfmlist[0][vidx],0.5*dt);
    }
    off_update_kobj_face_normal(dfm->ctx[0]);
    off_update_kobj_vert_normal(dfm->ctx[0]);
    off_smooth_kobj_vert_normal(dfm->ctx[0]);
    off_update_kobj_kface_pminmax(dfm->ctx[0]);

    if(dfm->weight->m[0][WEIGHT_CURVATURE]>0.0)
      update_curvature(&dfm->crv[0], dfm->ctx[0]);

    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Midpoint method]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }
      dfm->dfmlist[0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /* update dfmlist
     * put the original vertex back */
    if(verbose>=2) fprintf(stdout,"[%s]   output and input update [Midpoint method]\n",fcname);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      dfm->vmat[0][vidx]->v = op[0][vidx];
    }

    /* free memory */
    if(verbose>=2) fprintf(stdout,"[%s] free memory [Midpoint method]\n",fcname);
    free(op[0]);
    free(op);
    break;

  case NIIK_NUM_METHOD_RUNGE_KUTTA:

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);

    if(verbose>=1) fprintf(stdout,"[%s] Runge-Kutta 4th order method %s\n",fcname,tmstr);
    op = niikpt_matrix(1,dfm->ctx[0]->nvert);
    k = (niikpt ***)calloc(4,sizeof(niikpt **));

    for(n=0; n<4; n++) {
      if((k[n] = niikpt_matrix(2,dfm->ctx[0]->nvert))==NULL) {
        fprintf(stderr,"[%s] ERROR: niikpt_matrix\n",fcname);
        return 0;
      }
    }

    if(verbose>=2) fprintf(stdout,"[%s]   backup vertex positions [Runge-Kutta]\n",fcname);
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      op[0][vidx]=dfm->vmat[0][vidx]->v;
    }

    /* k1 = dCtx = dt Func(t,Ctx) */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k1 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[0][0][vidx] = k[0][1][vidx] = niikpt_zero();
        continue;
      }
      k[0][0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /*
     * k2 = dCtx2 = dt Func(t,Ctx1)
     * Ctx1 = Ctx + 0.5 * k1
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k2 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        dfm->vmat[0][vidx]->v = niikpt_move_normal(op[0][vidx],k[0][0][vidx],0.5*dt);
    }

    off_update_kobj_face_normal(dfm->ctx[0]);
    off_update_kobj_vert_normal(dfm->ctx[0]);
    off_smooth_kobj_vert_normal(dfm->ctx[0]);
    off_update_kobj_kface_pminmax(dfm->ctx[0]);

    if(dfm->weight->m[0][WEIGHT_CURVATURE]>0.0)
      update_curvature(&dfm->crv[0], dfm->ctx[0]);

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k2 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[1][0][vidx] = k[1][1][vidx] = niikpt_zero();
        continue;
      }

      k[1][0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /*
     * k3 = dCtx3 = dt Func(t,Ctx2)
     * Ctx2 = Ctx1 + 0.5 * k2
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k3 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      dfm->vmat[0][vidx]->v = niikpt_move_normal(op[0][vidx],k[1][0][vidx],0.5*dt);
    }
    off_update_kobj_face_normal(dfm->ctx[0]);
    off_update_kobj_vert_normal(dfm->ctx[0]);
    off_smooth_kobj_vert_normal(dfm->ctx[0]);
    off_update_kobj_kface_pminmax(dfm->ctx[0]);

    if(dfm->weight->m[0][WEIGHT_CURVATURE]>0.0)
      update_curvature(&dfm->crv[0],dfm->ctx[0]);

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k3 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[2][0][vidx] = k[2][1][vidx] = niikpt_zero();
        continue;
      }
      k[2][0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /*
     * k4 = dCtx4 = dt Func(t,Ctx3)
     */
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=3) fprintf(stdout,"[%s]   k4 preparation [Runge-Kutta] %s\n",fcname,tmstr);

    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        dfm->vmat[0][vidx]->v = niikpt_move_normal(op[0][vidx],k[2][0][vidx],0.5*dt);
    }
    off_update_kobj_face_normal(dfm->ctx[0]);
    off_update_kobj_vert_normal(dfm->ctx[0]);
    off_smooth_kobj_vert_normal(dfm->ctx[0]);
    off_update_kobj_kface_pminmax(dfm->ctx[0]);
    if(dfm->weight->m[0][WEIGHT_CURVATURE]>0.0)
      update_curvature(&dfm->crv[0],dfm->ctx[0]);

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>=2) fprintf(stdout,"[%s]   k4 calculation [Runge-Kutta] %s\n",fcname,tmstr);

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        k[3][0][vidx] = k[3][1][vidx] = niikpt_zero();
        continue;
      }

      k[3][0][vidx] = mesh_deform_calc_deformation_vertex(&dfm_refine,dfm->vmat[CORTEX_ICS][vidx],vidx);
    } /* each vertex */

    /* update dfmlist
     * put the original vertex back */
    if(verbose>=3) {
      ctm=time(NULL);
      stm=localtime(&ctm);
      strftime(tmstr,256,"%Y-%m-%d %T",stm);
      fprintf(stdout,"[%s]   output and input update [Runge-Kutta] %s\n",fcname,tmstr);
    }

    #pragma omp parallel for
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
        dfm->vmat[0][vidx]->v   = op[0][vidx];

        dfm->dfmlist[0][vidx].x = (k[0][0][vidx].x + 2*k[1][0][vidx].x + 2*k[2][0][vidx].x + k[3][0][vidx].x) / 6.0 * dt;
        dfm->dfmlist[0][vidx].y = (k[0][0][vidx].y + 2*k[1][0][vidx].y + 2*k[2][0][vidx].y + k[3][0][vidx].y) / 6.0 * dt;
        dfm->dfmlist[0][vidx].z = (k[0][0][vidx].z + 2*k[1][0][vidx].z + 2*k[2][0][vidx].z + k[3][0][vidx].z) / 6.0 * dt;
    }

    /* free memory */
    if(verbose>=2) {
      ctm=time(NULL);
      stm=localtime(&ctm);
      strftime(tmstr,256,"%Y-%m-%d %T",stm);
      fprintf(stdout,"[%s] free memory [Runge-Kutta] %s\n",fcname,tmstr);
    }

    for(n=0; n<4; n++) {
      free(k[n][0]);
      free(k[n]);
    }
    free(k);

    free(op[0]);
    free(op);
    break;

  default:
    fprintf(stderr,"[%s] ERROR: not implemented yet %i\n",fcname,dfm->numerical_method);
    return 0;
  }

  /*convert dfmlist into unit vectors with .w showing magnitude*/
  if(dfm->weight->m[0][WEIGHT_UPDATE_SIGMA]>0.0) {
    off_surface_field_smooth_using_vert(dfm->ctx[0],dfm->dfmlist[0],dfm->weight->m[0][WEIGHT_UPDATE_SIGMA]);
  }
  
  #pragma omp parallel for
  for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
    if( fabs(dfm->dfmlist[0][vidx].x)<g_niikmesh_min_deform_distance &&
        fabs(dfm->dfmlist[0][vidx].y)<g_niikmesh_min_deform_distance &&
        fabs(dfm->dfmlist[0][vidx].z)<g_niikmesh_min_deform_distance) {
      dfm->dfmlist[0][vidx].w = 0.0;
    } else {
      double d=niikpt_mag(dfm->dfmlist[0][vidx]);
      dfm->dfmlist[0][vidx] = niikpt_kmul(dfm->dfmlist[0][vidx],1.0/d); /*make a unit vector with length in w*/
      dfm->dfmlist[0][vidx].w = d;
    }
  }

  if(verbose>=1) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    fprintf(stdout,"[%s] exit function %s\n",fcname,tmstr);
  }
  return 1;
} /* niikmesh_deform_calc_deformation */


/*TODO: adapt this code to work without assumption that object intensity is always max (?)*/

int niikmesh_deform_mesh_calc_multiplication_factor(
    nifti_image *img, nifti_image *brain_mask, kobj **ctx,
    double intWM, double mflim, double *mf_list) {
  char *pstr;
  const char *fcname=__func__;

  /*traverse symmetrically*/
  double  xmin =-2.0,
          xmax = 2.0,
          xdel = 0.2,
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
  NIIK_RET0((ctx[0]==NULL),fcname,"mesh is null");
  NIIK_RET0((mf_list==NULL),fcname,"mf_list is null");

  nvert=ctx[0]->nvert;

  /* CALCULATE MULTIPLICATION FACTOR */
  mf_lim[0] = (mflim>1.0)?1.0/mflim:mflim;
  mf_lim[1] = 1.0 / mf_lim[0];

  if(verbose>=1) fprintf(stdout,"[%s] scaling factor limits: %8.5f %8.5f\n",fcname,mf_lim[0],mf_lim[1]);
  if(verbose>=1) fprintf(stdout,"[%s] object intensity       %8.2f\n", fcname, intWM);
  if((x=niikvec_init_range(xmin, xmax, xdel))==NULL) {
    fprintf(stderr,"[%s] ERROR: niikvec_init_range\n",fcname);
    return 0;
  }
  NIIK_RET0(((y=niikvec_init(x->num))==NULL),fcname,"niikvec_init");

  if(verbose>=1) fprintf(stdout,"[%s] interp vector          %5.2f %5.2f %5.2f %5i\n",fcname, x->v[0], x->v[1]-x->v[0], x->v[x->num-1], x->num);

  for(vi=ctx[0]->vert,vidx=0; vi!=NULL; vi=vi->next,vidx++) {
    if(!niik_image_interp_along_normal(img, NIIK_INTERP_LINEAR, vi->v, vi->normal, x, y)) {
      fprintf(stderr,"[%s] ERROR: niik_image_interp_along_normal, vidx=%i\n",fcname,vidx);
      return 0;
    }
    /*TODO: this is implicit rule that WM is the highest intensity ?*/
    mf_list[vidx] = NIIK_DMINMAX( niik_get_max_from_double_vector(y->v, x->num) / intWM, mf_lim[0], mf_lim[1]);
  }
  if(verbose>0) fprintf(stdout,"[%s]   smoothing multiplication factor\n",__func__);
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

  if(g_niikmesh_deform_calc_deformation_vertex_index>=0) {
    fprintf(stdout,"[%s] scaling factor[%i] %9.4f\n",
            fcname,
            g_niikmesh_deform_calc_deformation_vertex_index,
            mf_list[g_niikmesh_deform_calc_deformation_vertex_index]);
  }
  x=niikvec_free(x);
  y=niikvec_free(y);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niikmesh_deform_mesh_calc_multiplication_factor */


/**************
 * allocate or reallocate transient structures
 * 
 * 
 * */
int niicortex_mesh_deform_prepare(niikcortex_deform * dfm)
{
  kvert *vi,*vo;
  int vidx;
  dfm->nvert = dfm->ctx[0]->nvert;

  /* nonctx labels
   * -ctx_label is a list of zeros and non-zeros for non-cortex and cortex, respectively
   */
  if(dfm->nonctx_mask!=NULL) {
    /*TODO: finish this!*/
  } else {
    int n;
    dfm->ctx_label = (int *)realloc(dfm->ctx_label,dfm->nvert*sizeof(int));
    for(n=0; n<dfm->nvert; n++) dfm->ctx_label[n] = 1;
  } 

  dfm->indicator = (indicator_t*)realloc(dfm->indicator,dfm->nvert*sizeof(indicator_t));
 
  dfm->dfmlist[0] = (niikpt *)realloc(dfm->dfmlist[0], dfm->nvert*sizeof(niikpt));
  dfm->vmat[0]    = (kvert **)realloc(dfm->vmat[0], dfm->nvert*sizeof(kvert *));
  dfm->rd[0]  =     (double*) realloc(dfm->rd[0], dfm->nvert*sizeof(double));
  dfm->dd[0]  =     (double*) realloc(dfm->dd[0], dfm->nvert*sizeof(double));

  re_init_curvature(&dfm->crv[0], dfm->ctx[0]);

  for(vi=dfm->ctx[CORTEX_ICS]->vert,vidx=0;
      vi!=NULL;
      vi=vi->next,vidx++) {
    /*assign labels*/
    vi->idata=CORTEX_ICS;
    dfm->vmat[0][vidx]=vi;
  }
 
  /* CALCULATE MULTIPLICATION FACTOR */
  dfm->mf_lim[0] = (dfm->mflim>1.0)?1.0/dfm->mflim:dfm->mflim;
  dfm->mf_lim[1] = 1.0 / dfm->mf_lim[0];
  dfm->mf_list = (double *)realloc(dfm->mf_list,dfm->nvert*sizeof(double));
  /*debug*/
  memset(dfm->mf_list,0,dfm->nvert*sizeof(double));

  /*HACK for tissue value normalization*/
  if(!niikmesh_deform_mesh_calc_multiplication_factor(dfm->t1img, dfm->brain_mask, dfm->ctx,
    dfm->tissue_val[0]>dfm->tissue_val[1]?dfm->tissue_val[0]:dfm->tissue_val[1],
    dfm->mf_lim[0],dfm->mf_list)) {
    fprintf(stderr,"[%s] ERROR: niikmesh_deform_mesh_calc_multiplication_factor\n", __func__);
    return 0;
  }

  return 1;
}


int niikmesh_refine( niikcortex_deform *dfm) {
  nifti_image
    **grad_img[2],
    *div_img[1],
    *blur_img[2],
    *tmpimg;

  nifti_image **grad_prior[4];

  bbox *bb;
  kface *f;
  kvert *vi,*vo;

  niikmat *invmat;

  niikpt pt,normal;
  struct tm *stm;
  time_t ref_time;
  char tmstr[256];

  double
    dist,
    dval,
    dmin,dmax,dmean,dstdv,
    dmax_dfm, /* cumulative maximum deformation for bbox update */
    dratio[2],
    dt,
    intWM,
    mean_ics_value,range_ics_value;

  int
    n,nlo,nhi,
    ci[2],
    xsc=0,
    *this_label,
    cidx,vidx,
    iter,iter2,
    verbose=1 /*niik_verbose()*/,
    remesh_fail;

  int debug_tracing=0;
  const char * falcon_trace_log=NULL;
  const char * falcon_trace_log_id=NULL;
  cortex_tracing_info trace;

  char fname[512];
  const char *fcname=__func__;

  debug_tracing = falcon_tracing_init(dfm->t1img, &trace);

  ref_time=time(NULL);
  stm=localtime(&ref_time);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  if(verbose>=1) niik_fc_display(fcname,1);
  NIIK_RET0((     dfm->t1img==NULL), fcname,"t1img is null");
  NIIK_RET0((dfm->brain_mask==NULL), fcname,"brain mask is null");
  NIIK_RET0((  dfm->gwi_mask==NULL), fcname,"gwi_mask mask is null");

  dt              = dfm->weight->m[0][0];

  if(verbose>0) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"  t1w image            %s\n",dfm->t1img->fname);
    fprintf(stdout,"  brain mask           %s\n",dfm->brain_mask->fname);
    if(dfm->prior[0])
        fprintf(stdout,"  prior object         %s\n",dfm->prior[0]->fname);
    if(dfm->prior[1])
        fprintf(stdout,"  prior background       %s\n",dfm->prior[1]->fname);
    fprintf(stdout,"  surface              %s\n",dfm->ctx[CORTEX_ICS]->fname);
    fprintf(stdout,"  surface vfe          %i %i %i\n",dfm->ctx[CORTEX_ICS]->nvert,dfm->ctx[CORTEX_ICS]->nface,dfm->ctx[CORTEX_ICS]->nedge);
    fprintf(stdout,"  deform time step     %-7.4f\n",dt);
    fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n",dfm->apply_step);
    fprintf(stdout,"  surface sm weight   %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_SURFACE]);
    fprintf(stdout,"  image weight        %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_IMAGE]);
    fprintf(stdout,"  prior weight        %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_PRIOR]);
    fprintf(stdout,"  mask weight         %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_BRAIN_MASK]);
    fprintf(stdout,"  proximity weight    %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_PROXIMITY]);
    fprintf(stdout,"  gradient weight     %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_GRADIENT]);
    fprintf(stdout,"  curvature weight    %-7.4f \n",dfm->weight->m[CORTEX_ICS][WEIGHT_CURVATURE]);

    fprintf(stdout,"  max iter             %i\n", dfm->iter);
    fprintf(stdout,"  proximity min dist   %-7.2f\n",dfm->proximity_min_distance);
    if(dfm->gradient_lambda[0]>=0) {
      fprintf(stdout,"  gradient lambda      %7.3f\n",dfm->gradient_lambda[0]);
    } else {
      fprintf(stdout,"  gradient lambda      TBD\n");
    }
    fprintf(stdout,"  tolerance            %-7.4f\n",dfm->tolerance);
    if(dfm->ctx_label!=NULL)
      fprintf(stdout,"  non-cortex label     %-i / %i\n",(int)niik_count_zero_from_int_vector(dfm->ctx_label, dfm->nvert),dfm->nvert);
    else
      fprintf(stdout,"  non-cortex label     not used\n");

    fprintf(stdout,"  numerical method     %s\n",niik_numerical_method_string(dfm->numerical_method));
    fflush(stdout);
  }

  if(dfm->debug_pt!=NULL)
    if(dfm->debug_pt[0]>0)
      fprintf(stdout,"  debug position       %9.2f %9.2f %9.2f\n",dfm->debug_pt[1]*dfm->t1img->dx,dfm->debug_pt[2]*dfm->t1img->dy,dfm->debug_pt[3]*dfm->t1img->dz);

    if(verbose>0) fprintf(stdout,"[%s] gradient image\n",fcname);

    NIIK_RET0(((blur_img[0]   = niik_image_copy_as_type(dfm->t1img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
    NIIK_RET0(((blur_img[1]   = niik_image_copy_as_type(dfm->t1img,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");

    if(dfm->gradient_FWHM[0]>0.0) {
      if(verbose>1) fprintf(stdout,"[%s] gradient  gaussian filter\n",fcname);
      NIIK_RET0((!niik_image_filter_gaussian_update(blur_img[0],11,dfm->gradient_FWHM[0])),fcname,"niik_image_filter_gaussian_update");
      /*debug*/
      /*
      printf("Save:debug_blur.mnc blurred with fwhm:%f\n",dfm->gradient_FWHM[0]);
      NIIK_EXIT((!niik_image_write("debug_blur.mnc",blur_img[0])),__func__,"niik_image_write",9);*/

    } else if(verbose>1) {
      fprintf(stdout,"[%s]   no gradient gaussian filter\n",fcname);
    }

    if(dfm->divergence_FWHM[0]>0.0) {
      if(verbose>1) fprintf(stdout,"[%s] divergence  gaussian filter\n",fcname);
      NIIK_RET0((!niik_image_filter_gaussian_update(blur_img[1],11,dfm->divergence_FWHM[0])),fcname,"niik_image_filter_gaussian_update");
    } else if(verbose>1) {
      fprintf(stdout,"[%s]   no gradient gaussian filter\n",fcname);
    }

    if(verbose>1) fprintf(stdout,"[%s]   sobel filter\n",fcname);
    NIIK_RET0(((grad_img[0]  = niik_image_sobel_filters_with_mag(blur_img[0]))==NULL),fcname,"niik_image_sobel_filters_with_mag");
    NIIK_RET0(((grad_img[1]  = niik_image_sobel_filters_with_mag(blur_img[1]))==NULL),fcname,"niik_image_sobel_filters_with_mag");

    /*divergence of the gradient field*/
    NIIK_RET0(((div_img[0] = niik_image_divergence(grad_img[1],0))==NULL),fcname, "niik_image_divergence");

    /*debug*/
    /*NIIK_EXIT((!niik_image_write("tmp_divergence.mnc",div_img[0])),__func__,"niik_image_write",9);*/

    /*debug*/

  for(n=0; n<2; n++)
    blur_img[n]=niik_image_free(blur_img[n]);


  for(n=0; n<3; n++) {
    if(dfm->prior[n]) {
        char tmp_fname[1024];
        int k;
        nifti_image *blur_prior;

        NIIK_RET0(((blur_prior = niik_image_copy_as_type(dfm->prior[n],NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");

        if(dfm->prior_FWHM[0]>0.0) { /*TODO: use different FWHM for ICS adn OCS*/
          if(verbose>1) fprintf(stdout,"[%s] gradient  gaussian filter\n",fcname);
          NIIK_RET0((!niik_image_filter_gaussian_update(blur_prior,11,dfm->prior_FWHM[0])),fcname,"niik_image_filter_gaussian_update");

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

  if(dfm->gradient_lambda[0]<0) {
      if(verbose>1) fprintf(stdout,"[%s] gradient lambda calculation\n",fcname);
      dfm->gradient_lambda[0] = niik_image_get_percentile(grad_img[0][0], dfm->brain_mask, 0.7); /*TODO: figure out if we need it*/
      if(dfm->gradient_lambda[0]<1e-7) dfm->gradient_lambda[0]=1e-7;
      /*debug*/
      /*NIIK_EXIT((!niik_image_write("debug_gradient_mag.mnc", grad_img[0][0])),__func__,"niik_image_write",9);*/
      fprintf(stdout,"[%s] gradient lambda      %7.5f\n",fcname, dfm->gradient_lambda[0]);
  } else {
      if(verbose>1) fprintf(stdout,"[%s] skip gradient lambda calculation\n",fcname);
  }

  if(verbose>2) fprintf(stdout,"[%s] allocate memory\n",__func__);

  if(dfm->dfm_limit)
  {
    /*assume that it contains one item*/
    /*TODO: make it work with preset map*/
    double _dfm_limit=dfm->dfm_limit[0];
    dfm->dfm_limit=(double*)realloc(dfm->dfm_limit, sizeof(double)*dfm->ctx[CORTEX_ICS]->nvert);

    for(vidx=0;vidx<dfm->ctx[CORTEX_ICS]->nvert;vidx++) dfm->dfm_limit[vidx]=_dfm_limit;

    dfm->travel = (double*)calloc(dfm->ctx[CORTEX_ICS]->nvert, sizeof(double));
  }
  /*make sure we have normals*/
  off_update_kobj_face_normal(dfm->ctx[0]);
  off_update_kobj_vert_normal(dfm->ctx[0]);
  off_smooth_kobj_vert_normal(dfm->ctx[0]);

  NIIK_RET0(!niicortex_mesh_deform_prepare(dfm), fcname, "niicortex_mesh_deform_prepare");

  if(verbose>=3) fprintf(stdout,"[%s] bbox initialization\n",fcname);
  bb=off_bbox_init(dfm->bbox_depth,320);

  if(verbose>=1) fprintf(stdout,"[%s] bbox initialization %7.4f %i\n",fcname,bb->delta,bb->depth);

  /* VERTEX */
  if(dfm->debug_pt && dfm->debug_pt[0]>0) {
    pt.x = dfm->debug_pt[1] * dfm->t1img->dx; /*voxel to world here?*/
    pt.y = dfm->debug_pt[2] * dfm->t1img->dy;
    pt.z = dfm->debug_pt[3] * dfm->t1img->dz;
    for(vi=dfm->ctx[0]->vert; vi!=NULL; vi=vi->next,vidx++) {
      dval = niikpt_distance(vi->v,pt);
      if(dval<dmin) {
        dmin = dval;
        g_niikmesh_deform_calc_deformation_vertex_index=vi->index;
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

  if(verbose>0)
    fprintf(stdout,"[%s]   %9.6f +/- %9.6f  (%9.6f %9.6f %9.6f %9.6f)\n",
            fcname,
            niik_get_mean_from_double_vector(dfm->mf_list, dfm->nvert),
            niik_get_stdv_from_double_vector(dfm->mf_list, dfm->nvert),
            niik_get_min_from_double_vector (dfm->mf_list, dfm->nvert),
            niik_get_percentile_from_double_vector(dfm->mf_list, dfm->nvert,0.25),
            niik_get_percentile_from_double_vector(dfm->mf_list, dfm->nvert,0.75),
            niik_get_max_from_double_vector(dfm->mf_list,dfm->nvert));

  if(g_niikmesh_deform_calc_deformation_vertex_index>=0) {
    if(verbose>0) fprintf(stdout,"   scaling factor %9.4f\n",dfm->mf_list[g_niikmesh_deform_calc_deformation_vertex_index]);
  }

  if(verbose>3) {
    if(!niikcortex_add_color(dfm->ctx[0],dfm->mf_list,0.8,1.2,NIIK_COLORMAP_SPECTRAL)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_add_color\n",fcname);
      return 0;
    }

    /* write output */
    if(verbose>3 && get_DEBUG_PREFIX()) {
      sprintf(fname,"%s_tmp_refine.ply",get_DEBUG_PREFIX());
      fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
      off_kobj_write_offply(fname,dfm->ctx[0],0);
    }
    off_kobj_add_one_color(dfm->ctx[0],1,0,0); /* red  */
  }
  /* end of multiplication factor calculation */

  /* check surface intersection before deformation */
  /*
  xsc = niikcortex_off_count_intersection(bb,dfm->ctx[0], dfm->ctx[1]);
  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);*/
  /*TODO: fix?*/
  if(dfm->convergence_log)
    fprintf(dfm->convergence_log,"iteration,wtime,mean_d,std_d,min_d,max_d\n");

  /************************************************
   *
   * MAIN LOOP STARTS HERE
   *
   ************************************************/
  if(verbose>2) fprintf(stdout,"[%s] main loop\n",__func__);
  ref_time=time(NULL);

  for(iter=0, dmax_dfm=0.0; iter<dfm->iter; iter++,dmax_dfm+=(dfm->apply_step * dfm->iter2)) {
    time_t ctm;
    /*
     * prepare for deformation
     * -normal, cortical thickness, and bounding box
     */
    cidx=0;
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    fprintf(stdout,"[%s] iter %-4i of %-4i %s\n",fcname,iter+1,dfm->iter,tmstr);
    if(verbose>2) fprintf(stdout,"[%s] iter %-4i update normal\n",fcname,iter+1);

    off_update_kobj_face_normal(dfm->ctx[cidx]);
    off_update_kobj_vert_normal(dfm->ctx[cidx]);
    off_smooth_kobj_vert_normal(dfm->ctx[cidx]);
    if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);

    if(dfm->weight->m[cidx][WEIGHT_CURVATURE]>0.0)
       fprintf(stdout,"iter %-4i mesh mean curvature=%8.6f\n",iter,mean_curvature(&dfm->crv[cidx],dfm->ctx[cidx]));

    if(dmax_dfm>=bb->delta*0.45 ) { /* FOR TESTING, remove || 1 if debugged properly */
      if(verbose>1) fprintf(stdout,"[%s] iter %-4i bbox updating\n",fcname,iter+1);
      NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,1)),fcname,"off_create_bbox_from_multiple_kobj");
      dmax_dfm=0;
    }

    if(verbose>2) fprintf(stdout,"[niikmesh_deform_cortex] iter %-4i update thickness\n",iter+1);

    /*TODO: make this parameter (5 iterations)*/
    if(!(iter+1)%5) {
      if(!niikmesh_deform_mesh_calc_multiplication_factor(dfm->t1img,dfm->brain_mask,dfm->ctx,dfm->tissue_val[0],dfm->mf_lim[0],dfm->mf_list)) {
        fprintf(stderr,"[%s] ERROR: niikmesh_deform_mesh_calc_multiplication_factor\n",fcname);
        return 0;
      }
    }

    /******************************************************
     * calculate deformation
     ******************************************************/

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>1) fprintf(stdout,"[%s] iter %-4i calculate deformation %s\n",fcname, iter+1,tmstr);

    dfm->bb=bb;
    dfm->cortex_id=0;

    NIIK_RET0((!mesh_deform_calc_deformation(dfm)), fcname,"mesh_deform_calc_deformation");

    if(dfm->dfm_limit)
    {
      for(vidx=0; vidx<dfm->nvert; vidx++) {
        /*check if we achieved threshold*/
        if( (dfm->travel[vidx]+dfm->dfmlist[0][vidx].w)>=dfm->dfm_limit[vidx])
        {
            dfm->dfmlist[0][vidx].w=dfm->dfm_limit[vidx]-dfm->travel[vidx];

            if(dfm->dfmlist[0][vidx].w<0.0)
                dfm->dfmlist[0][vidx].w=0.0;
            /*stop if we reached the limit*/
        }
        /*account for the movement*/
      }
    }

    if(dfm->convergence_log)
    {
      double dmin=dfm->dfmlist[0][0].w;
      double dmax=dfm->dfmlist[0][0].w;
      double dmean=0.0,dstdv=0.0;

      for(vidx=0; vidx<dfm->nvert; vidx++) {
          if     (dmin>dfm->dfmlist[0][vidx].w) dmin=dfm->dfmlist[0][vidx].w;
          else if(dmax<dfm->dfmlist[0][vidx].w) dmax=dfm->dfmlist[0][vidx].w;
          dmean += dfm->dfmlist[0][vidx].w;
          dstdv += NIIK_SQ(dfm->dfmlist[0][vidx].w);
      }
      dmean /= dfm->nvert;
      dstdv = sqrt(dstdv / dfm->nvert - dmean * dmean);

      fprintf(dfm->convergence_log,"%d,%d,%f,%f,%f,%f\n",iter, (int)(ctm-ref_time), dmean, dstdv, dmin, dmax);
    }
    if(verbose>1) fprintf(stdout,"[%s] iter %-4i calculate deformation done %s\n",fcname,iter+1,tmstr);

    /* show deformation stats */

    /*************************************
     * apply deformation
     *************************************/
    strftime(tmstr,256,"%Y-%m-%d %T",localtime(&ctm));
    if(verbose>0) {
      fprintf(stdout,"[%s] iter %-4i apply deformation %s\n",fcname,iter+1,tmstr);
    }
    dratio[0]=dratio[1]=0;
    niik_set_zero_for_double_vector(dfm->dd[0],dfm->nvert);

    NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,1)),fcname,"off_create_bbox_from_multiple_kobj");

    if(debug_tracing)
    {
      /*VF: maybe do it more effeciently?*/
      NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[0],1,0,0)),fcname,"off_kobj_add_one_color"); /* yellow for white matter-surface */
      falcon_tracing_dump(&trace, iter,"refine", dfm->t1img, bb);
      falcon_tracing_dump_objects(&trace, iter,"refine", dfm->ctx, 1 );
    }

    for(iter2=0,dist=dfm->apply_step; iter2<dfm->iter2; iter2++) {

      #ifndef _USE_PRL_DEFORM
      NIIK_RET0((!niikcortex_deform_apply_deformation(bb,dfm->ctx[0], dist, dfm->dfmlist[0], dfm->dd[0])),
                fcname,"niikcortex_deform_apply_deformation ");
      #else
      NIIK_RET0((!niikcortex_deform_apply_deformation_prl(bb,dfm->ctx[0], dist, dfm->dfmlist[0], dfm->dd[0], dfm->indicator, 0)),
                fcname,"niikcortex_deform_apply_deformation ");
      #endif

      dmin=dmax=dfm->dd[0][0];
      dmean=dstdv=0;

      /*This is just for tracing progress, can be removed*/
      for(vidx=0; vidx<dfm->nvert; vidx++) {
          dfm->rd[0][vidx] = dfm->dfmlist[0][vidx].w; /* deformation distance calculated but not actually applied */
          if     (dmin>dfm->dd[0][vidx]) dmin=dfm->dd[0][vidx];
          else if(dmax<dfm->dd[0][vidx]) dmax=dfm->dd[0][vidx];
          dmean += dfm->dd[0][vidx];
          dstdv += NIIK_SQ(dfm->dd[0][vidx]);
      }
      dmean /= dfm->nvert;
      dstdv = sqrt(dstdv / dfm->nvert - dmean * dmean);

      if(verbose>1 && iter2==(dfm->iter2-1)) { /* show deformation stats */
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
    } /* apply step */

    cidx=ci[0]=ci[1]=0;
    for(vidx=0; vidx<dfm->nvert; vidx++) {

        if(dfm->travel)
            dfm->travel[vidx]+=dfm->dd[cidx][vidx]; /*apply actual travel distance*/

        if(dfm->dd[cidx][vidx]<g_niikmesh_min_deform_distance)
          ci[cidx]++;
    }
    fprintf(stdout,"[%s] deform below threshold: %9.5f %%\n",fcname, 100.0*ci[0]/dfm->nvert);


     /******************************************************
     * check for the stop condition
     * -to be implemented properly
     ******************************************************/

    /*xsc = niikcortex_off_count_intersection(bb, white_surface, pial_surface);
    fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);*/
    if(ci[0]==dfm->nvert) {
        fprintf(stdout,"[%s] No movement, stopping \n",fcname);
        break;
    }
  } /* iteration */

//  xsc = niikcortex_off_count_intersection(bb, dfm->ctx[0], dfm->ctx[1]);
//  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);

//  if(verbose>3) {
//    fprintf(stdout,"[niikcortex_deform_cortex] write temp files\n");
//    off_kobj_write_offply("tmp_dfm_final_ics.off.gz",dfm->ctx[0],0);
//    off_kobj_write_offply("tmp_dfm_final_ocs.off.gz",dfm->ctx[1],0);
//  }


  off_bbox_free(bb);

  for(n=0; n<2; n++) {
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

  div_img[0]=niik_image_free(div_img[0]);

  invmat = niikmat_free(invmat);
  falcon_tracing_free(&trace);

  if(verbose>1) niik_fc_display(fcname,0);
  return 1;
}

void prog_history() {
  fprintf(stdout,"[falcon_mesh_refine] history\n");
  fprintf(stdout,"  version  0.0.0  February 26, 2019, Vladimir S. FONOV <vladimir.fonov@gmail.com>\n");
}

void usage() {
  fprintf(stdout,"falcon_mesh_refine\n");
  fprintf(stdout,"  usage: [options] <img> <global_mask> <object_mask> <in.off/ply> <out.off/ply>\n\n");
  fprintf(stdout,"\n");
  fprintf(stdout,"  optional usage:\n");
  fprintf(stdout,"  -u -help --help                   : show this usage\n");
  fprintf(stdout,"  --version                         : show version info\n");
  fprintf(stdout,"  -debug-keep-tmp                   : keep debug files (obsolete)\n");
  fprintf(stdout,"  -log <file>                       : write convergence log to file\n");
  fprintf(stdout,"  -travel <file>                    : write amount of deformation applied\n");
  fprintf(stdout,"  PDE Solver:                       \n");
  fprintf(stdout,"  -backward-euler                   : use backward euler method\n");
  fprintf(stdout,"  -forward-euler                    : use forward euler method (DEFAULT)\n");
  fprintf(stdout,"  -rungekutta                       : use runge-kutta method\n");
  fprintf(stdout,"  -midpoint                         : use midpoint method\n");
  fprintf(stdout,"  Processing options:               \n");
  fprintf(stdout,"  -lesion-mask   <img>              : lesion mask\n");
  fprintf(stdout,"  -nonctx-mask   <img>              : mask away non-cortex\n");
  fprintf(stdout,"  -priorobj      <img>              : priors for object\n");
  fprintf(stdout,"  -priorbkg      <img>              : priors for background\n");
  fprintf(stdout,"  -dfmlimit      <max>              : deformation maximum\n");
  fprintf(stdout,"  Optimization weights \n");
  fprintf(stdout,"  -wprox    <white>          : proximity weights, default 0.2 \n");
  fprintf(stdout,"  -wprior   <white>          : prior weights, default 0.0 \n");
  fprintf(stdout,"  -wimag    <white>          : image weights, default 1.0 \n");
  fprintf(stdout,"  -wgrad    <white>          : gradient weights, default 0.0 \n");
  fprintf(stdout,"  -wsmooth  <white>          : tangential surface smoothness weights, default 1.0 \n");
  fprintf(stdout,"  -wssmooth <white>          : simple surface smoothness weights, default 0.0 \n");
  fprintf(stdout,"  -wmask    <white>          : brain mask weights, default 1.0 \n");
  fprintf(stdout,"  -wles     <white>          : lesion mask weights, default 1.5 \n");
  fprintf(stdout,"  -wcurv    <white>          : curvature weights, default 0.0 \n");
  fprintf(stdout,"  -wflux    <white>          : flux term , default 0.0 \n");
  fprintf(stdout,"  Update smoothing \n");
  fprintf(stdout,"  -supdate  <white>          : update smoothing sigma, default 0.0 \n");
  fprintf(stdout,"  Proximity distance constraints \n");
  fprintf(stdout,"  -pmin <val>                : minimum proximity distance, default 0.6 \n");

  fprintf(stdout,"  Additional optimizer parameters \n");
  fprintf(stdout,"  -depth <n>                 : Quad-tree depth , default 7\n");
  fprintf(stdout,"  -gradient-FWHM <f>         : blurring kernels for image gradient calculation (voxels), default 1.0 \n");
  fprintf(stdout,"  -divergence-FWHM <f>       : blurring kernels for image divergence calculation (voxels), default 1.0 \n");
  fprintf(stdout,"  -prior-FWHM <f>            : blurring kernels for prior gradient calculation (voxels), default 1.0 \n");
  fprintf(stdout,"  -delta <f>                 : time-step (default 0.5)\n");
  fprintf(stdout,"  -apply <f>                 : apply-step (default 0.2)\n");
  fprintf(stdout,"  -iter  <n>                 : maximum number of iterations, default 100\n");
  fprintf(stdout,"  -iter2 <n>                 : maximum number of sub-iterations, default 5\n");
  fprintf(stdout,"  -remesh <n>                : remesh every nth iteration, default 0 (disabled)\n");
  fprintf(stdout,"\n\n");
  fprintf(stdout,"  Environment variables: \n");
  fprintf(stdout,"   FALCON_TRACE - debug tracing prefix\n");
  fprintf(stdout,"     output will be named  <prefix>_<stage>_<iteration>_<object>.ply - for surface \n");
  fprintf(stdout,"     or  <prefix>_<stage>_<iteration>_<x|y|z>_<slice>.tiff - for volume slice \n");
  fprintf(stdout,"   FALCON_TRACE_SCALE - debug tracing image scale, default 1.0 \n");
  fprintf(stdout,"   FALCON_TRACE_X,FALCON_TRACE_Y,FALCON_TRACE_Z - (n1)[,n2[,n3]] - slices to trace (voxel coordinates) \n");
  fprintf(stdout,"   FALCON_TRACE_SURF - set to any value to dump surfaces \n");
}

int niikmesh_deform_process(niikcortex_deform *dfm) {
  const char *fcname=__func__;
  niik_fc_display(fcname,1);

  NIIK_RET0((dfm==NULL),fcname,"missing self");

  fprintf(stdout,"[%s] parameters\n",fcname);
  fprintf(stdout,"  t1w image            %s\n", dfm->t1img->fname);
  fprintf(stdout,"  Global mask             %s\n", dfm->brain_mask->fname);
  fprintf(stdout,"  Object mask             %s\n", dfm->gwi_mask->fname);
  if(dfm->prior[0])
    fprintf(stdout,"  obj prior             %s\n", dfm->prior[0]->fname);
  if(dfm->prior[1])
    fprintf(stdout,"  bck prior             %s\n", dfm->prior[1]->fname);

  fprintf(stdout,"  surface              %s\n", dfm->ctx[0]->fname);
  fprintf(stdout,"  surface vfe          %i %i %i\n", dfm->ctx[0]->nvert, dfm->ctx[0]->nface, dfm->ctx[0]->nedge);
  fprintf(stdout,"  deform time step     %-7.4f\n", dfm->delta);
  fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n", dfm->apply_step);
  fprintf(stdout,"  surface sm weight    %-7.4f \n", dfm->weight->m[0][WEIGHT_SURFACE]);
  fprintf(stdout,"  image intensity weight %-7.4f \n", dfm->weight->m[0][WEIGHT_IMAGE]);
  fprintf(stdout,"  brainmask weights    %-7.4f \n", dfm->weight->m[0][WEIGHT_BRAIN_MASK]);
  fprintf(stdout,"  lesion mask weights  %-7.4f \n", dfm->weight->m[0][WEIGHT_LESION]);
  fprintf(stdout,"  proximity weights    %-7.4f \n", dfm->weight->m[0][WEIGHT_PROXIMITY]);
  fprintf(stdout,"  gradient weights     %-7.4f \n", dfm->weight->m[0][WEIGHT_GRADIENT]);
  fprintf(stdout,"  curvature weights    %-7.4f \n", dfm->weight->m[0][WEIGHT_CURVATURE]);
  fprintf(stdout,"  flux weights         %-7.4f \n", dfm->weight->m[0][WEIGHT_DGRADIENT]);
  fprintf(stdout,"  update sigma         %-7.4f \n", dfm->weight->m[0][WEIGHT_UPDATE_SIGMA]);
  fprintf(stdout,"  prior weight         %-7.4f \n", dfm->weight->m[0][WEIGHT_PRIOR]);

  fprintf(stdout,"  max iter             %i\n", dfm->iter);
  fprintf(stdout,"  max iter2            %i\n", dfm->iter2);
  fprintf(stdout,"  remesh               %i\n", dfm->remesh);
  fprintf(stdout,"  object mean intensity %9.2f\n", dfm->tissue_val[0]);
  fprintf(stdout,"  background intensity  %9.2f\n", dfm->tissue_val[1]);
  fprintf(stdout,"  intensity range       %9.2f\n", dfm->tissue_val[3]);

  fprintf(stdout,"  proximity min dist   %-7.2f\n", dfm->proximity_min_distance);
  fprintf(stdout,"  quadtree depth       %i\n", dfm->bbox_depth);
  fprintf(stdout,"  gradient FWHM    %-7.4f\n", dfm->gradient_FWHM[0]);
  fprintf(stdout,"  divergence FWHM  %-7.4f\n", dfm->divergence_FWHM[0]);
  fprintf(stdout,"  prior FWHM       %-7.4f\n", dfm->prior_FWHM[0]);

  if(dfm->gradient_lambda[0]>=0) {
    fprintf(stdout,"  surface gradient lambda      %7.3f\n",dfm->gradient_lambda[0]);
  } else {
    fprintf(stdout,"  surface gradient lambda      TBD\n");
  }

  fprintf(stdout,"  tolerance            %-7.4f\n", dfm->tolerance);

  fprintf(stdout,"  numerical method     %s\n", niik_numerical_method_string(dfm->numerical_method));

  niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_process */


int main(int argc,char *argv[],char *envp[]) {
  niikcortex_deform *dfm=NULL;
  const char *fcname="niikcortex_mesh_refine";
  char *travel_out=NULL;
  int n=0,i=0,nc=1,sc=1;
  double *ijk;

  nifti_image *debugimg;
  char* timestamp=niik_create_minc_timestamp(argc,argv);

  if(argc==1) {
    usage();
    exit(0);
  }

#ifdef HAVE_OPENMP
  fprintf(stderr,"[%s] Using OpenMP, max number of threads=%d\n",fcname,omp_get_max_threads());
#endif

  NIIK_EXIT(((dfm=niikcortex_deform_init())==NULL),fcname,"niikcortex_deform_init",1);
  niik_version_display(fcname, NIIK_MAJOR_VERSION, NIIK_MINOR_VERSION, NIIK_MICRO_VERSION);
  niik_fc_display(fcname,1);

  ijk=(double *)calloc(9,sizeof(double));
  ijk[0]=0;

  while(nc<argc) {
    if(argv[nc][0]=='-') {
      if(!strncmp(argv[nc],"--version",9)) {
        prog_history();
        exit(0);
      } else if(!strncmp(argv[nc],"--help",6)) {
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-backward-euler",15)) {
        dfm->numerical_method=NIIK_NUM_METHOD_BACKWARD_EULER;
      } /* backward-euler method */

      else if(!strncmp(argv[nc],"-forward-euler",14)) {
        dfm->numerical_method=NIIK_NUM_METHOD_FORWARD_EULER;
      } /* forward-euler method */

      else if(!strncmp(argv[nc],"-midpoint",14)) {
        dfm->numerical_method=NIIK_NUM_METHOD_MIDPOINT;
      } /* forward-euler method */

      else if(!strncmp(argv[nc],"-rungekutta",11)) {
        dfm->numerical_method=NIIK_NUM_METHOD_RUNGE_KUTTA;
      } /* runge-kutta (4th order) */

      else if(!strncmp(argv[nc],"-debug-keep-tmp",15)) {
        dfm->debug_keep_tmp=1;
      }

      else if(!strncmp(argv[nc],"-gradient-FWHM",13)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -gradient-FWHM\n",fcname);
          exit(1);
        }
        dfm->gradient_FWHM[0] = atof(argv[nc]);
      } /* gradient-FWHM */

      else if(!strncmp(argv[nc],"-divergence-FWHM",16)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -divergence-FWHM\n",fcname);
          exit(1);
        }
        dfm->divergence_FWHM[0] = atof(argv[nc]);
      } /* divergence-FWHM */

      else if(!strncmp(argv[nc],"-prior-FWHM",11)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argments for -prior-FWHM\n",fcname);
          exit(1);
        }
        dfm->prior_FWHM[0] = atof(argv[nc]);
      } /* prior-FWHM */

      else if(!strncmp(argv[nc],"-lesion-mask",12)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading lesion    %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->lesion_mask = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-priorobj",9)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading prior   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->prior[0] = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-priorbkg",9)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        fprintf(stdout,"[%s] reading prior   %s\n",fcname,argv[nc]);
        NIIK_EXIT(((dfm->prior[1] = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      }

      else if(!strncmp(argv[nc],"-delta",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->delta = atof(argv[nc]);
      } /* delta */

      else if(!strncmp(argv[nc],"-apply",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->apply_step = atof(argv[nc]);
      } /* delta */

      else if(!strncmp(argv[nc],"-iter2",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->iter2 = atoi(argv[nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-iter",5)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->iter = atoi(argv[nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-remesh",7)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->remesh = atoi(argv[nc]);
      } /* iter */

      else if(!strncmp(argv[nc],"-nonctx-mask",12)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        NIIK_EXIT(((dfm->nonctx_mask = niik_image_read(argv[nc]))==NULL),fcname,"niik_image_read",9);
      } /* nonctx-mask */

      else if(!strncmp(argv[nc],"-wprox",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprox <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_PROXIMITY] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wgrad",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wgrad <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_GRADIENT] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wflux",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wflux <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_DGRADIENT] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wimag",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wimag <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_IMAGE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wcurv",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wcurv <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_CURVATURE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wsmooth",8)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wsmooth <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_SURFACE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wssmooth",9)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wssmooth <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_SURFACE_SIMPLE] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wmask",6)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wmask <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_BRAIN_MASK] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wles",5)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wles <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_LESION] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-wprior",7)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -wprior <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_PRIOR] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-supdate",8)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -supdate <white> \n",fcname);
          exit(1);
        }
        dfm->weight->m[0][WEIGHT_UPDATE_SIGMA] = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-pmin",5)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -pmin <dist>\n",fcname);
          exit(1);
        }
        dfm->proximity_min_distance = atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-depth",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s)\n",fcname);
          exit(1);
        }
        dfm->bbox_depth = atoi(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-log",4)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -log <fname>\n",fcname);
          exit(1);
        }
        dfm->convergence_log=fopen(argv[nc],"w");
      }

      else if(!strncmp(argv[nc],"-dfmlimit",9)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -dfmlimit <f>\n",fcname);
          exit(1);
        }
        dfm->dfm_limit=(double*)calloc(1,sizeof(double));
        dfm->dfm_limit[0]=atof(argv[nc]);
      }

      else if(!strncmp(argv[nc],"-travel",7)) {
        if((nc++)>=argc) {
          fprintf(stderr,"[%s] ERROR: missing argment(s), -travel <fname>\n",fcname);
          exit(1);
        }
        travel_out=strdup(argv[nc]);
      }


      else if(!strncmp(argv[nc],"-help",5)) {
        usage();
        exit(0);
      }

      else if(!strncmp(argv[nc],"-u",2)) {
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

  if(argc<6) {
    fprintf(stderr,"[%s] ERROR: too few argments (%i != 6)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  } else if(argc>6) {
    fprintf(stderr,"[%s] ERROR: too many argments (%i != 6)\n",fcname,argc);
    for(i=0; i<argc; i++) fprintf(stderr," %s",argv[i]);
    fprintf(stderr,"\n");
    exit(1);
  }
  dfm->debug_pt=ijk;

  fprintf(stdout,"[%s] reading t1w image    %s\n",fcname,argv[1]);
  NIIK_EXIT(((dfm->t1img=niik_image_read(argv[1]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading mask image   %s\n",fcname,argv[2]);
  NIIK_EXIT(((dfm->brain_mask=niik_image_read(argv[2]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading object image  %s\n",fcname,argv[3]);
  NIIK_EXIT(((dfm->gwi_mask=niik_image_read(argv[3]))==NULL),fcname,"niik_image_read",9);
  fprintf(stdout,"[%s] reading init object  %s\n",fcname,argv[4]);
  NIIK_EXIT(((dfm->ctx[0]=off_kobj_read_offply(argv[4]))==NULL),fcname,"off_kobj_read_off",9);

  NIIK_EXIT((!niik_image_type_convert(dfm->brain_mask,  NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, brain_mask",1);
  NIIK_EXIT((!niik_image_type_convert(dfm->gwi_mask,    NIFTI_TYPE_UINT8   )),fcname,"niik_image_convert, gwi_mask",1);

  /* get the tissue values (niikcortex-deform) */
  if(dfm->tissue_val[0]>0  && dfm->tissue_val[1]>0 ) {
    fprintf(stdout,"[%s] using predefined tissue intensities\n",fcname);
  } else {
    fprintf(stdout,"[%s] niikcortex_estimate_tissue_values\n",fcname);
    NIIK_EXIT((!niikmesh_estimate_tissue_values(dfm->t1img,
               dfm->brain_mask,
               dfm->gwi_mask,
               &dfm->tissue_val[0],
               &dfm->tissue_val[1],
               &dfm->tissue_val[2],
               &dfm->tissue_val[3])),
              fcname,"niikmesh_estimate_tissue_values",9);
  }

  /* dt */
  dfm->weight->m[CORTEX_ICS][WEIGHT_DT] =
    dfm->weight->m[CORTEX_OCS][WEIGHT_DT] =
      dfm->delta;

  /*Make sure outputs are flushed before long process starts*/
  fflush(stdout);
  fflush(stderr);

  /*display all paramters*/
  niikmesh_deform_process(dfm);

  /* cortex deformation (niikcortex-deform) */
  NIIK_EXIT((!niikmesh_refine(dfm)), fcname, "niikmesh_refine",9);

  /*append metadata*/
  NIIK_RET1((!off_kobj_add_comment(dfm->ctx[0], timestamp)),fcname,"off_kobj_add_comment");

  fprintf(stdout,"[%s] yellow color for  surface\n",fcname);
  NIIK_RET1((!off_kobj_add_one_color(dfm->ctx[0], 0.8,0.8,0)),fcname,"off_kobj_add_one_color");

  /* writing output for pial surface */
  fprintf(stdout,"[%s] writing surface:   %s\n",fcname,argv[5]);
  NIIK_RET1((!off_kobj_write_offply(argv[5],dfm->ctx[0], 0)),fcname,"off_kobj_write_off");

  if(dfm->convergence_log) fclose(dfm->convergence_log);

  if(travel_out && dfm->travel ){
    /*write output file*/
    NIIK_RET1((!niik_write_double_vector_ex(travel_out, dfm->travel, dfm->ctx[0]->nvert, "travel")),fcname,"niik_write_double_vector_ex");


  }
  if(travel_out)
    free(travel_out);

  dfm=niikcortex_deform_free(dfm);

  free(ijk);
  niik_fc_display(fcname,0);
  free(timestamp);
  return 0;
} /* main */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
