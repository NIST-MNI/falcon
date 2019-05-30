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

/*from falcon_cortex_deformation.c*/
int niikcortex_deform_apply_deformation_prl(bbox *bb,kobj *obj,double step,niikpt *dfmlist,double *actdfm,indicator_t *indicator,int vfilter);
niikpt niikcortex_deform_calc_deformation_vertex_surface_smoothness_simple(kvert *v);


/**************
 * allocate or reallocate transient structures
 * 
 * 
 * */

int niicortex_deform_prepare_smoothing(niikcortex_deform * dfm)
{
  const char* fcname="niicortex_deform_prepare";
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

  return 1;
}


niikpt niikcortex_smooth_calc_deformation_vertex(refine_info *dfm, kvert *vi, kvert *vo, int cortex_id, int vidx, int iter)
/* -computes deformation for a given vertex /cortex
 */
{
  niikpt
      normal,/*Vertex normal*/
      psum,  /*sum of all forces*/
      psmooth;  /*smoothing force (local average)*/
  
  /*TODO: make it a parameter!*/
  /*double smoothing_terms[]={0.33,-0.34};*/

  kvert *v[2];

  v[0]=vi;
  v[1]=vo;

  psum = psmooth = niikpt_zero();
  
  psmooth = niikcortex_deform_calc_deformation_vertex_surface_smoothness_simple(
              v[cortex_id]);


  /* combine deformation with weighting */
  psum.x = psmooth.x * dfm->deform_weights->m[cortex_id][iter%2];
  psum.y = psmooth.y * dfm->deform_weights->m[cortex_id][iter%2];
  psum.z = psmooth.z * dfm->deform_weights->m[cortex_id][iter%2];

  return psum;
} /* niikcortex_deform_calc_deformation_vertex */


int niikcortex_deform_calc_deformation_smoothing( niikcortex_deform * dfm,int iter)
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
  const char *fcname="niikcortex_deform_calc_deformation";
  char  tmstr[256];
  const int verbose=3;

  /*TODO: move this up one level?*/
  refine_info  dfm_refine;
  memset(&dfm_refine,0,sizeof(refine_info));

  dfm_refine.img = dfm->t1img;
  dfm_refine.bb  = dfm->bb;
  dfm_refine.brain_mask = dfm->brain_mask;
  dfm_refine.cerebellum_mask = dfm->cerebellum_mask;
  dfm_refine.ventricle_mask = dfm->csf_mask;
  dfm_refine.crv = dfm->crv;
  dfm_refine.deform_weights = dfm->weight;
  dfm_refine.div_img = dfm->div_img;
  dfm_refine.grad_img = dfm->grad_img;
  dfm_refine.gradient_lambda = dfm->gradient_lambda;
  dfm_refine.lesion_mask = dfm->lesion_mask;
  dfm_refine.mean_ics_value = dfm->tissue_val[4];
  dfm_refine.mean_ocs_value = dfm->tissue_val[5];
  dfm_refine.range_ics_value = dfm->tissue_val[6];
  dfm_refine.range_ocs_value = dfm->tissue_val[7];
  dfm_refine.intWM = dfm->tissue_val[1];
  dfm_refine.vmat = dfm->vmat;
  dfm_refine.dfm = dfm;

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

    #pragma omp parallel for private (cidx)
    for(vidx=0; vidx<dfm->ctx[0]->nvert; vidx++) {
      if(!dfm->ctx_label[vidx]) {
        dfm->dfmlist[0][vidx] = dfm->dfmlist[1][vidx] = niikpt_zero();
        continue;
      }

      for(cidx=0; cidx<2; cidx++) {
        if(cidx==0 && dfm->cortex_id==2) continue;
        if(cidx==1 && dfm->cortex_id==1) continue;

        dfm->dfmlist[cidx][vidx] = niikcortex_smooth_calc_deformation_vertex(&dfm_refine, dfm->vmat[CORTEX_ICS][vidx], dfm->vmat[CORTEX_OCS][vidx], cidx, vidx, iter);
      } /* each cortex */
    } /* each vertex */

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
          dfm->dfmlist[cidx][vidx] = niikpt_kmul(dfm->dfmlist[cidx][vidx],1.0/d); /*make a unit vector with length in w*/
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
int niikcortex_smooth_cortex(niikcortex_deform * dfm)
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
  const char *fcname="niikcortex_smooth_cortex";
  bbox *bb;
  kface *f;
  kvert *vi,*vo;
    
  niikpt pt,normal;
  struct tm *stm;
  time_t ctm;
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
    verbose=niik_verbose();
  
  int cortex_id=dfm->cortex_id;
  int debug_tracing=0;
  const char * falcon_trace_log=NULL;
  const char * falcon_trace_log_id=NULL;
  cortex_tracing_info trace;

  char fname[512];

  debug_tracing=falcon_tracing_init(dfm->t1img, &trace);

  ctm=time(NULL);
  stm=localtime(&ctm);
  strftime(tmstr,256,"%Y-%m-%d %T",stm);

  if(verbose>=1) niik_fc_display(fcname,1);
  if(cortex_id>3 || cortex_id<1) {
    fprintf(stderr,"[%s] ERROR: invalid cortex_id, %i\n",fcname, cortex_id);
    return 0;
  }

  dt              = dfm->weight->m[0][0];

  if(verbose>0) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"  white surface        %s\n",dfm->ctx[CORTEX_ICS]->fname);
    fprintf(stdout,"  pial surface         %s\n",dfm->ctx[CORTEX_OCS]->fname);
    fprintf(stdout,"  surface vfe          %i %i %i\n",dfm->ctx[CORTEX_ICS]->nvert,dfm->ctx[CORTEX_ICS]->nface,dfm->ctx[CORTEX_ICS]->nedge);
    fprintf(stdout,"  deform time step     %-7.4f\n",dt);
    fprintf(stdout,"  deform apply step    %-7.4f    for each deform-apply\n",dfm->apply_step);
    fprintf(stdout,"  surface sm weights   %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_SURFACE],dfm->weight->m[CORTEX_OCS][WEIGHT_SURFACE]);
    fprintf(stdout,"  image weights        %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_IMAGE],dfm->weight->m[CORTEX_OCS][WEIGHT_IMAGE]);
    fprintf(stdout,"  thickness sm weights %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_THICKNESS_SMOOTHNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_THICKNESS_SMOOTHNESS]);
    fprintf(stdout,"  brainmask weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_BRAIN_MASK],dfm->weight->m[CORTEX_OCS][WEIGHT_BRAIN_MASK]);
    fprintf(stdout,"  ventricle weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_VENTRICLE],dfm->weight->m[CORTEX_OCS][WEIGHT_VENTRICLE]);
    fprintf(stdout,"  cerebellum weights   %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_CEREBELLUM],dfm->weight->m[CORTEX_OCS][WEIGHT_CEREBELLUM]);
    fprintf(stdout,"  lesion mask weights  %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_LESION],dfm->weight->m[CORTEX_OCS][WEIGHT_LESION]);
    fprintf(stdout,"  proximity weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_PROXIMITY],dfm->weight->m[CORTEX_OCS][WEIGHT_PROXIMITY]);
    fprintf(stdout,"  abs thickness        %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_ABS_THICKINESS],dfm->weight->m[CORTEX_OCS][WEIGHT_ABS_THICKINESS]);
    fprintf(stdout,"  gradient weights     %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_GRADIENT],dfm->weight->m[CORTEX_OCS][WEIGHT_GRADIENT]);
    fprintf(stdout,"  curvature weights    %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_CURVATURE],dfm->weight->m[CORTEX_OCS][WEIGHT_CURVATURE]);

    fprintf(stdout,"  min thicknes         %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_MIN_THICKNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_MIN_THICKNESS]);
    fprintf(stdout,"  max thicknes         %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_MAX_THICKNESS],dfm->weight->m[CORTEX_OCS][WEIGHT_MAX_THICKNESS]);
    fprintf(stdout,"  thicknes  sigma      %-7.4f %-7.4f\n",dfm->weight->m[CORTEX_ICS][WEIGHT_THICKNESS_SIGMA],dfm->weight->m[CORTEX_OCS][WEIGHT_THICKNESS_SIGMA]);

    fprintf(stdout,"  max iter             %i\n", dfm->iter);
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
    fflush(stdout);
  }

  if(dfm->debug_pt!=NULL)
    if(dfm->debug_pt[0]>0)
      fprintf(stdout,"  debug position       %9.2f %9.2f %9.2f\n",dfm->debug_pt[1]*dfm->t1img->dx,dfm->debug_pt[2]*dfm->t1img->dy,dfm->debug_pt[3]*dfm->t1img->dz);

  if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] allocate memory\n");

  NIIK_RET0(!niicortex_deform_prepare_smoothing(dfm),fcname,"niicortex_deform_prepare_smoothing");

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

  if(verbose>3) {
    if(!niikcortex_add_color(dfm->ctx[0],dfm->mf_list,0.8,1.2,NIIK_COLORMAP_SPECTRAL)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_add_color\n",fcname);
      return 0;
    }

    if(!niikcortex_add_color(dfm->ctx[1],dfm->mf_list,0.8,1.2,NIIK_COLORMAP_SPECTRAL)) {
      fprintf(stderr,"[%s] ERROR: niikcortex_add_color\n",fcname);
      return 0;
    }

    /* write output */
    if(verbose>3) {
      sprintf(fname,"tmp_dfm_white.off.gz");
      fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
      off_kobj_write_offply(fname,dfm->ctx[0],0);
      sprintf(fname,"tmp_dfm_pial.off.gz");
      fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
      off_kobj_write_offply(fname,dfm->ctx[1],0);
    }
    off_kobj_add_one_color(dfm->ctx[0],1,1,0); /* yellow for white matter-surface */
    off_kobj_add_one_color(dfm->ctx[1],1,0,0); /* red for pial-surface */
  }
  /* end of multiplication factor calculation */

  /* check surface intersection before deformation */
  xsc = niikcortex_off_count_intersection(bb,dfm->ctx[0], dfm->ctx[1]);
  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);
  /*TODO: fix?*/

  /************************************************
   *
   * MAIN LOOP STARTS HERE
   *
   ************************************************/
  if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] main loop\n");
  for(iter=0,dmax_dfm=0; iter<dfm->iter; iter++,dmax_dfm+=(dfm->apply_step * dfm->iter2)) {
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
      {
        update_curvature(&dfm->crv[cidx],dfm->ctx[cidx]);
        fprintf(stdout,"iter %-4i cortex %i mean curvature=%8.6f\n",iter,cidx,mean_curvature(&dfm->crv[cidx],dfm->ctx[cidx]));
      }
    }

    if(dmax_dfm>=bb->delta*0.45 ) { /* FOR TESTING, remove || 1 if debugged properly */
      if(verbose>1) fprintf(stdout,"[%s] iter %-4i bbox updating\n",fcname,iter+1);
      NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,dfm->ctx,2)),fcname,"off_create_bbox_from_multiple_kobj");
      dmax_dfm=0;
    }

    if(verbose>2) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i update thickness\n",iter+1);

    /******************************************************
     * calculate deformation
     ******************************************************/

    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);
    if(verbose>1) fprintf(stdout,"[niikcortex_deform_cortex] iter %-4i calculate deformation %s\n",iter+1,tmstr);

    dfm->bb=bb;
    dfm->cortex_id=cortex_id;

    NIIK_RET0((!niikcortex_deform_calc_deformation_smoothing(dfm,iter)),
              fcname,"niikcortex_deform_calc_deformation");

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
        dmin=dmax=dfm->dfmlist[cidx][cidx].w;
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

    if(debug_tracing && dfm->t1img)
    {
      /*VF: maybe do it more effeciently?*/
      NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[0],1,1,0)),fcname,"off_kobj_add_one_color"); /* yellow for white matter-surface */
      NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[1],1,0,0)),fcname,"off_kobj_add_one_color"); /* red for pial-surface */

      falcon_tracing_dump(&trace,iter,"smooth",dfm->t1img,bb);
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

    /*****************************************************
     * write temporary output
     *****************************************************/
    if(verbose>3 && (iter%5)==0) {
      if(cortex_id%2) {
        NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[0],1,1,0)),fcname,"off_kobj_add_one_color ics");
        sprintf(fname,"tmp_dfm%03i_ics.off.gz",iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[0],0);
      }
      if(cortex_id>=2) {
        NIIK_RET0((!off_kobj_add_one_color(dfm->ctx[1],1,0,0)),fcname,"off_kobj_add_one_color ocs");
        sprintf(fname,"tmp_dfm%03i_ocs.off.gz",iter+1);
        fprintf(stdout,"[%s] write temp files %s\n",fcname,fname);
        off_kobj_write_offply(fname,dfm->ctx[1],0);
      }
    } /* writing temp output */

    
  } /* iteration */

  if(dfm->debug_trace_log) fclose(dfm->debug_trace_log);

  xsc = niikcortex_off_count_intersection(bb, dfm->ctx[0], dfm->ctx[1]);
  fprintf(stdout,"[%s] surface intersection %i\n",fcname,xsc);

  if(verbose>3) {
    fprintf(stdout,"[%s] write temp files\n",fcname);
    off_kobj_write_offply("tmp_dfm_final_ics.off.gz",dfm->ctx[0],0);
    off_kobj_write_offply("tmp_dfm_final_ocs.off.gz",dfm->ctx[1],0);
  }

  off_bbox_free(bb);
  falcon_tracing_free(&trace);
  if(verbose>1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_deform_cortex */


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
