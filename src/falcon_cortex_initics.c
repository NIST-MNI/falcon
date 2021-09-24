/* Filename:     nifti1_kunio_cortex_initics.c
 * Description:  functions for crating initial white matter surface
 * Author:       Kunio Nakamura
 * Date:         March 25, 2012
 */

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1

#endif

#include "falcon.h"
#include "falcon_cortex.h"

int niikcortex_initics_shrink(nifti_image *gwi_img,nifti_image *lap_map,nifti_image *dist_map,int maxiter,double initlen,double finlen,double *dfm_step,int dfm_iter,kobj *obj,niikpt check_pt);
kobj *niikcortex_init_ics_atimg(nifti_image *at_ics_img,double elen);


static kvert *g_niikcortex_initics_check_v=NULL;


/**********************************************************************
 *
 * niikcortex_initics_expand
 *
 * -expands from atlas-based object to brain surface
 *
 **********************************************************************/

int niikcortex_initics_expand(nifti_image *brain_img,   /* atlas-based initial white matter surface image */
                              double max_size,
                              double step_size,
                              kobj *obj) {              /* input and output here */
  bbox *bb=NULL;
  kvert *v;
  niikpt origpt;
  double dist;
  int
  iter=0,
  verbose=0;
  const char *fcname=__func__;

  niik_fc_display(fcname,1);
  NIIK_RET0((brain_img==NULL),fcname,"brain_img is null");

  if(verbose>=0) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"    mask     %s\n",brain_img->fname);
    fprintf(stdout,"    step     %7.3f\n",step_size);
    fprintf(stdout,"    max      %7.3f\n",max_size);
  }
  bb=off_bbox_init(7,320);

  for(dist=0,iter=0; dist<=max_size; dist+=step_size,iter++) {
    if(verbose>=0) fprintf(stdout,"[%s] balloon %-3i %7.3f %7.3f\n",fcname,iter+1,dist,step_size);
    if(verbose>=2) fprintf(stdout,"[%s] update face/vert normals\n",fcname);

    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_smooth_kobj_vert_normal(obj);

    if(verbose>1) fprintf(stdout,"[%s] update bbox\n",fcname);
    off_update_kobj_kface_pminmax(obj);
    off_create_bbox_from_kobj(bb,obj);

    if(verbose>1) fprintf(stdout,"[%s] balloon\n",fcname);
    for(v=obj->vert; v!=NULL; v=v->next) {
      origpt=v->v;

      /* check for brain mask */
      if(niik_image_interpolate_3d_nn(brain_img,v->v)==0) {
        continue;
      }

      /* move */
      v->v = niikpt_move_normal(v->v,v->normal,step_size);

      /* check for intersection */
      off_update_kvert_pminmax(v);
      if(off_check_self_intersection_kvert(bb,v)) {
        v->v=origpt;
        off_update_kvert_pminmax(v);
        continue;
      }
    } /* each vertex */

    if(dist<max_size/2.0 && 0) {
      for(v=obj->vert; v!=NULL; v=v->next) {
        origpt = v->v;
        v->v=niikpt_kvert_local_average(v,1);
        off_update_kvert_pminmax(v);
        if(off_check_self_intersection_kvert(bb,v)) {
          v->v=origpt;
          off_update_kvert_pminmax(v);
          continue;
        }
      }
    }

    off_update_kobj_kface_pminmax(obj);
    off_create_bbox_from_kobj(bb,obj);
    if(verbose>1) fprintf(stdout,"[%s] # self-intersection %i\n",fcname,off_count_self_intersection(bb,obj));
  } /* iteratoins */

  bb=off_bbox_free(bb);

  niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_initics_expand */




/**********************************************************************
 *
 * niikcortex_initics_atimg
 *
 * -shrink-wrap based on atlas-based initial wihte matter surface image
 * -currently this function is unused because shrink-wrap function works
 *  sufficiently
 *
 **********************************************************************/

kobj *niikcortex_initics_atimg(nifti_image *at_ics_img,   /* atlas-based initial white matter surface image */
                               double elen) {             /* target edge length */
  kobj *obj;
  niikpt ctr;
  int verbose=1;
  const char *fcname=__func__;

  if(verbose) fprintf(stdout,"[%s] start function\n",fcname);
  if(at_ics_img==NULL) {
    fprintf(stderr,"ERROR: at_ics_img is null\n");
    return 0;
  }
  ctr=niikpt_image_get_centroid(at_ics_img,NULL);
  if(niik_check_double_problem(ctr.x)) {
    fprintf(stderr,"ERROR: niikpt_image_get_centroid(maskimg,NULL)\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"[%s] cetnroid = %7.3f %7.3f %7.3f\n",fcname,ctr.x,ctr.y,ctr.z);
  if((obj=off_make_sphere_from_icosahedron(3,300,ctr))==NULL) {
    fprintf(stderr,"ERROR: off_make_sphere_from_icosahedron(3,300,ctr)\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"[%s] initial spherical object\n",fcname);
  if(!off_shrinkwrap_kobj_simple(at_ics_img,obj,elen)) {
    fprintf(stderr,"ERROR: off_shrinkwrap_kobj_simple(img,obj,elen)\n");
    exit(0);
  }
  if(verbose) fprintf(stdout,"[%s] shrink-wrap\n",fcname);
  if(verbose) fprintf(stdout,"[%s] exit function\n",fcname);
  return obj;
}



/**********************************************
 *
 * remesh function with angular resampling
 *
 **********************************************/

int niikcortex_initics_shrink_remesh_angular_region(kobj *obj,double elen,int maxiter,int wmark,niikpt ctr, niikpt sph, double elen2,double athresh)
/* -ctr is the center of object
 * -elen2 is the target regional edge length
 * -par1 is the first angle
 * -par2 is the second angle
 * -athresh is the angular threshold
 */
{
  double emin,emax,dval,dang;
  kedge *e,*ne,*erm[9];
  kface *f,*frm[9];
  kvert *v,*vrm[2];
  int n,iter;
  int verbose=0;
  const char *fcname=__func__;  
  emax=4.0/3*elen;
  emin=4.0/5*elen;
  if(verbose) {
    fprintf(stdout,"[%s] %8.4f %8.4f %8.4f\n",fcname,elen,emin,emax);
    fprintf(stdout,"    regional elen = %8.4f\n",elen2);
    fprintf(stdout,"    vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }
  for(iter=1; iter<=maxiter; iter++) {
    /***************************
     * SPLIT EDGES
     ****************************/
    if(verbose) {
      fprintf(stdout,"    split edges  %8.4f \n",emax);
    }
    for(e=obj->edge; e!=NULL; e=e->next) {
      if(wmark)
        if(e->endpts[0]->v.w>0 || e->endpts[1]->v.w>0)
          continue;
      dval = niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
      dang = niikpt_angle_between_vectors(sph,niikpt_unit(niikpt_sub(e->endpts[0]->v,ctr)));
      if(fabs(dang)>athresh) {
        if(dval<emax) continue;
      } else {
        if(dval<4.0/3.0*elen2) continue;
      }
      if((v=off_remesh_kedge_split(obj,e))==NULL) {
        fprintf(stderr,"ERROR: off_remesh_kedge_split %i\n",e->index);
        return 0;
      }
    } /* edge split */
    /***************************
     * COLLAPSE EDGES
    ****************************/
    if(verbose) {
      fprintf(stdout,"    collapse edges  %8.4f \n",emin);
    }
    for(e=obj->edge; e!=NULL; e=ne) {
      ne=e->next;
      if(wmark)
        if(e->endpts[0]->v.w>0 || e->endpts[1]->v.w>0)
          continue;
      dval = niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
      dang = niikpt_angle_between_vectors(sph,niikpt_unit(niikpt_sub(e->endpts[0]->v,ctr)));
      if(dang>180) dang-=180;
      if(fabs(dang)>athresh) {
        if(dval>emin) continue;
      } else {
        if(dval>4.0/5.0*elen2) continue;
      }
      if(verbose>2) fprintf(stdout,"\t  collapse %i  %8.5f\n",e->index,dval);
      if((n=off_remesh_kedge_collapse(obj,e,vrm,frm,erm))==0) {
        fprintf(stderr,"ERROR: off_remesh_kedge_collapse \n");
        return 0;
      }
      /* the edge was not collapsed */
      if(n<0) continue;
      /* mark the vertex for next step */
      erm[0]->endpts[0]->v.w=erm[0]->endpts[1]->v.w=0;
      /* the edge was actually collapsed */
      v=erm[0]->endpts[0];
      if(verbose>3) {
        fprintf(stdout,"\t\tcollapse clean up\n");
        fprintf(stdout,"\t\t  v[%i] \n",vrm[0]->index);
        fprintf(stdout,"\t\t  f[%i %i] \n",frm[0]->index,frm[1]->index);
        fprintf(stdout,"\t\t  e[%i %i %i] \n",erm[0]->index,erm[1]->index,erm[2]->index);
      }
      /* make sure that the next edge is not removed */
      for(ne=e->next; ne!=NULL; ne=ne->next) {
        if(ne==erm[0]) continue;
        if(ne==erm[1]) continue;
        if(ne==erm[2]) continue;
        break;
      }
      if(verbose>3) fprintf(stdout,"\t\tcollapse clean up !\n");
      /* remove and free */
      off_remesh_kedge_collapse_clean_up(obj,vrm,frm,erm);
      vrm[0]=NULL;
      frm[0]=frm[1]=NULL;
      erm[0]=erm[1]=erm[2]=NULL;
      if(verbose>3) fprintf(stdout,"\t\tcollapse clean up done\n");
    } /* each edge for collapsing */
    /* correct 3 neighbors */
    if(!off_remesh_kobj_collapse_correction(obj)) {
      fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction \n");
      return 0;
    }
    /***************************
     * VALENCE MINIMIATION
    ****************************/
    if(verbose) fprintf(stdout,"    minimize the valence\n");
    for(n=0; n<3; n++) {
      if(!off_remesh_kobj_valence_minimization(obj,wmark)) {
        fprintf(stderr,"ERROR: off_remesh_kobj_valence_minimization\n");
        return 0;
      }
    }
    if(verbose) {
      n = off_remesh_kobj_count_global_valence57(obj);
      off_kobj_update_num(obj);
      fprintf(stdout,"        valence   = %3.0f%% = %i / %i\n",100.0*n/obj->nvert,n,obj->nvert);
      if(verbose>1) {
        fprintf(stdout,"    check status and update index \n");
        for(v=obj->vert,n=1; v!=NULL; v=v->next) {
          v->index=n++;
        }
        obj->nvert=n-1;
        for(f=obj->face,n=1; f!=NULL; f=f->next) {
          f->index=n++;
        }
        obj->nface=n-1;
        for(e=obj->edge,n=1; e!=NULL; e=e->next) {
          e->index=n++;
        }
        obj->nedge=n-1;
        fprintf(stdout,"      vfe = %i %i %i \n",obj->nvert,obj->nface,obj->nedge);
        fprintf(stdout,"      mean elen = %9.5f \n",off_get_kobj_mean_edge_length(obj));
        if(!off_kobj_test(obj)) {
          fprintf(stderr,"ERROR: test obj\n");
          exit(0);
        }
      }
    }
    if(!off_remesh_kobj_collapse_correction(obj)) {
      fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction\n");
      return 0;
    }
    /***************************
     * TANGENTIAL RELAXATION
    ****************************/
    if(verbose) fprintf(stdout,"    relocation: tangential relaxation\n");
    for(n=0; n<3; n++) {
      if(wmark) {
        for(v=obj->vert; v!=NULL; v=v->next)
          if(v->v.w==0) off_remesh_kvert_relocate(v);
      } else {
        for(v=obj->vert; v!=NULL; v=v->next)
          off_remesh_kvert_relocate(v);
      }
    }
    if(verbose)fprintf(stdout,"           remesh vfe = %i %i %i \n",obj->nvert,obj->nface,obj->nedge);
  } /* iteration */
  for(v=obj->vert,n=1; v!=NULL; v=v->next) {
    v->index=n++;
  }
  obj->nvert=n-1;
  for(f=obj->face,n=1; f!=NULL; f=f->next) {
    f->index=n++;
  }
  obj->nface=n-1;
  for(e=obj->edge,n=1; e!=NULL; e=e->next) {
    e->index=n++;
  }
  obj->nedge=n-1;
  if(verbose) fprintf(stdout,"          remeshed to vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  if(obj->color) {
    off_kobj_add_color(obj);
  }
  return 1;
} /* niikcortex_initics_shrink_remesh_angular_region */




/************************************************************
 *
 * MAIN FUNCTION FOR INITIAL WHITE SURFACE SHRINK-WRAP
 *
 *************************************************************/

int niikcortex_initics_shrink(nifti_image *gwi_img,      /* subject's gray-white interface mask image -> convert to distance */
                              nifti_image *lap_map,      /* laplace map */
                              nifti_image *dist_map,      /* laplace map */
                              int maxiter,               /* max iteration */
                              double initlen,            /* initial edge length */
                              double finlen,             /* final edge length */
                              double *dfm_step,          /* deformation step: [0]=initial, [1]=final, linearly interpolated inbetween (around 0.3?) */
                              int dfm_iter,              /* deformation iteration */
                              kobj *obj,                 /* white matter object -- updated here */
                              niikpt check_pt) {
  const char *fcname=__func__;
  kvert *v,*check_v=NULL;
  double
    phi,psi,
    dval,dsum,dmin,
    thresh = 0.5,
    angular_thresh = 30, /* regional edge-spliting's threshold in degrees */
    currstepsum=0,
    currstep, /* current deformation step function of (iter,maxiter,dfm_step[0-1]) */
    elen,elen2;
  int
    dphi,dpsi,
    vpos,
    iter,iter2,
    wsum,
    verbose=2;
  bbox *bb=NULL;
  niikpt
    objctr,
    origpt,
    vec,
    sf;
  niikmat
    *colormap;
  char
    fname[512];
  struct tm *stm;
  time_t ctm;
  char tmstr[256];

  if(verbose>=1) niik_fc_display(fcname,1);
  if(gwi_img==NULL) {
    fprintf(stderr,"ERROR: gwi_img is null\n");
    return 0;
  }
  if(lap_map==NULL) {
    fprintf(stderr,"ERROR: lap_map is null\n");
    return 0;
  }

  if(initlen<finlen) {
    fprintf(stderr,"ERROR: init len, %7.3f, is smaller than final len, %7.3f\n",initlen,finlen);
    return 0;
  }

  if(verbose>=1) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"  Initial object         %s\n",obj->fname);
    fprintf(stdout,"  GM-WM interface mask   %s\n",gwi_img->fname);
    fprintf(stdout,"  Laplace map            %s\n",lap_map->fname);
    /*if(dist_map!=NULL)
      fprintf(stdout,"  Distance map           %s\n",dist_map->fname);*/
    fprintf(stdout,"  Max iter               %i\n",maxiter);
    fprintf(stdout,"  Deform step            %-7.3f -> %-7.3f\n",dfm_step[0],dfm_step[1]);
    fprintf(stdout,"  Deform iter            %i\n",dfm_iter);
    fprintf(stdout,"  Edge length            %-7.3f -> %7.3f\n",initlen,finlen);
  }

  if(verbose>1) fprintf(stdout,"[%s] set up bounding box\n",fcname);
  bb=off_bbox_init(7,320); /* 5 or 6? Kunio 2012-09-12 */
  if(verbose>=1) fprintf(stdout,"[%s] bbox %8.4f\n",fcname,bb->delta);

  objctr = niikpt_image_get_centroid(gwi_img,NULL);
  if(verbose>=1) fprintf(stdout,"[%s] object's center: %7.3f %7.3f %7.3f\n",fcname,objctr.x,objctr.y,objctr.z);
  dpsi = 800 / maxiter;
  dphi = 400 / maxiter;
  if(verbose>=1) fprintf(stdout,"[%s] angular step: %3i %3i\n",fcname,dpsi,dphi);

  if((colormap=niik_colormap_get(NIIK_COLORMAP_SUMMER,maxiter))==NULL) {
    fprintf(stderr,"ERROR: niik_colormap_get(NIIK_COLORMAP_SUMMER,maxiter)\n");
    return 0;
  }


  /*************************
   * MAIN LOOP
   *************************/

  for(iter=1; iter<=maxiter; iter++) {
    ctm=time(NULL);
    stm=localtime(&ctm);
    strftime(tmstr,256,"%Y-%m-%d %T",stm);

    currstep = niik_pv((double)iter/maxiter,dfm_step[0],dfm_step[1]);
    elen = exp(-(iter-1.0)/maxiter*6.0)*(initlen-finlen)+finlen;
    if(iter==maxiter) elen=finlen;
    if(verbose>=1) fprintf(stdout,"[%s] ITERATION %4i %6.3f %s\n",fcname,iter,elen,tmstr);

    /* normal vector calculation */
    if(verbose>1) fprintf(stdout,"[%s] %4i update face/vert normals\n",fcname,iter);
    off_update_kobj_face_normal(obj);
    off_update_kobj_vert_normal(obj);
    off_smooth_kobj_vert_normal(obj);

    /* balloon sometimes */
    if(((iter%5)==4) && iter<maxiter*0.8 && iter>maxiter*0.2) {
      off_update_kobj_kface_pminmax(obj);
      off_create_bbox_from_kobj(bb,obj);
      if(verbose>1) fprintf(stdout,"[%s] %4i balloon\n",fcname,iter);
      for(iter2=0; iter2<4; iter2++) {
        for(v=obj->vert; v!=NULL; v=v->next) {
          origpt=v->v;
          /* go outward */
          v->v=niikpt_move_normal(v->v,v->normal,currstep*0.25);
          off_update_kvert_pminmax(v);
          if(off_check_self_intersection_kvert(bb,v))  {
            v->v=origpt;
            off_update_kvert_pminmax(v);
            continue;
          }
        } /* for each vertex */
        off_update_kobj_face_normal(obj);
        off_update_kobj_vert_normal(obj);
        off_smooth_kobj_vert_normal(obj);
        if(verbose>2) fprintf(stdout,"[%s] %4i balloon finished\n",fcname,iter);
      }
    } /* sometimes (iter%5)==0) */

    /* normal from laplace map */
    if(iter%2 && iter>5) {
      if(verbose>1) fprintf(stdout,"[%s] %4i laplace-normal\n",fcname,iter);
      for(v=obj->vert; v!=NULL; v=v->next) {
        vpos = niik_image_get_index_niikpt(lap_map,v->v);
        if(vpos<0) continue;
        sf = niik_image_sobel_filter_voxel(lap_map,vpos);
        if(fabs(sf.x) + fabs(sf.y) + fabs(sf.z) > 1e-3) {
          v->normal = niikpt_unit(sf);
        }
      }
    }

    /* set w-mark to NOT do remesh */
    if(verbose>1) fprintf(stdout,"[%s] %4i set w-marks\n",fcname,iter);
    for(v=obj->vert; v!=NULL; v=v->next) {
      v->v.w=1;
    }

    /* bounding box */
    if(verbose>1) fprintf(stdout,"[%s] %4i update bbox\n",fcname,iter);
    off_update_kobj_kface_pminmax(obj);
    off_create_bbox_from_kobj(bb,obj);

    /* deformation */
    if(verbose>=2) fprintf(stdout,"[niikcortex_initics_shrink] %4i start deformation (step=%0.4f)\n",iter,currstep);
    if(dist_map!=NULL) {
      thresh = cos(19.0*iter/maxiter) * exp(-3.0*iter/maxiter);
    }
    if(verbose>=1) fprintf(stdout,"[niikcortex_initics_shrink] %4i deform thresh = %7.3f\n",iter,thresh);
    if(check_pt.w>0) {
      for(v=obj->vert,dmin=1e6; v!=NULL; v=v->next) {
        dval = niikpt_distance(check_pt,v->v);
        if(dval<dmin) {
          dmin=dval;
          check_v=v;
        }
      }
      fprintf(stderr,"[niikcortex_initics_shrink] %8.3f %8.3f %8.3f\n",check_v->v.x,check_v->v.y,check_v->v.z);
    } else {
      check_v=NULL;
    }
    g_niikcortex_initics_check_v=check_v;
    for(iter2=0,currstepsum=0; iter2<dfm_iter; currstepsum+=currstep,iter2++) {
      if(currstepsum>=bb->delta*0.6) break;
      /* if(verbose>2) fprintf(stdout,"[niikcortex_initics_shrink] %4i deformation %i\n",iter,iter2); */
      if(dist_map!=NULL) { /* using distance map */
        for(v=obj->vert,dsum=0; v!=NULL; v=v->next) {
          dval = niikcortex_initics_shrink_check_deform(dist_map,bb,v,currstep,thresh);
          if(v==check_v) {
            fprintf(stderr,"[niikcortex_initics_shrink] %8.3f %8.3f %8.3f  |  %8.3f  |  %i\n",v->v.x,v->v.y,v->v.z,dval,v->index);
          }
          if(fabs(dval)>0.2) v->v.w=0;
          dsum += dval;
        }
        for(v=obj->vert,dval=0; v!=NULL; v=v->next) {
          dval+=(v->v.w==0);
        }
        if(verbose>=2) fprintf(stdout,"[niikcortex_initics_shrink] %4i deform %3i %8.2f | #v (no motion) %i / %i\n",iter,iter2,dsum,(int)dval,obj->nvert);
      } /* using distance map */
      else { /* no distance map */
        for(v=obj->vert,dsum=0; v!=NULL; v=v->next) {
          dval = niikcortex_initics_shrink_check_deform(gwi_img,bb,v,currstep,thresh);
          if(v==check_v) {
            fprintf(stderr,"[niikcortex_initics_shrink] %8.3f %8.3f %8.3f  |  %8.3f  |  %i\n",v->v.x,v->v.y,v->v.z,dval,v->index);
          }
          if(fabs(dval)>1e-5) v->v.w=0;
          dsum += dval;
        }
        for(v=obj->vert,dval=0; v!=NULL; v=v->next) {
          dval+=(v->v.w==0);
        }
        if(verbose>=2) fprintf(stdout,"[niikcortex_initics_shrink] %4i deform %3i %8.2f | #v (no motion) %i / %i\n",iter,iter2,dsum,(int)dval,obj->nvert);
      } /* not using distance map */
    } /* inner deformation loop */

    if(verbose>=1) {
      if(!off_kobj_add_one_color(obj,0,1,0)) {
        fprintf(stderr,"ERROR: off_kobj_add_one_color\n");
        return 0;
      }
      sprintf(fname,"tmp_initics_shrink%ip.ply",iter);
      fprintf(stdout,"[niikcortex_initics_shrink] %4i write output file: %s\n",iter,fname);
      off_kobj_write_offply(fname,obj,0);
    }


    /************************************
     * remesh
     ************************************/
    for(v=obj->vert,wsum=0; v!=NULL; v=v->next) {
      wsum+=(v->v.w==0);
    }
    if(iter==maxiter) {
      if(verbose>1) fprintf(stdout,"[%s] final remesh %7.3f\n",fcname,finlen);
      NIIK_RET0((!off_remesh_kobj(obj,finlen,20,0)),fcname,"off_remesh_kobj");
      NIIK_RET0((!off_display_kobj_edge_stats(obj)),fcname,"off_display_kobj_edge_stats");
    } else if(iter>NIIK_IMAX(10,maxiter*0.1) && iter<(maxiter*0.9)) {
      psi = NIIK_DEGREE2RAD((dpsi * iter +  90)%360);
      phi = NIIK_DEGREE2RAD((dphi * iter + 270)%360);
      vec.x = cos(psi) * sin(phi);
      vec.y = sin(psi) * sin(phi);
      vec.z =           -cos(phi);
      vec.w = 0;
      elen2 = NIIK_DMAX(elen*0.3,NIIK_DMAX(NIIK_DMIN(gwi_img->dx,NIIK_DMIN(gwi_img->dy,gwi_img->dz))*0.3,0.3));
      if(verbose>1) fprintf(stdout,"[niikcortex_initics_shrink] %4i remesh  %i / %i | %5.2f %5.2f | %5.2f %5.2f %5.2f from %5.2f %5.2f\n",
                              iter,wsum,obj->nvert,elen,elen2,
                              vec.x,vec.y,vec.z,
                              NIIK_RAD2DEGREE(psi),NIIK_RAD2DEGREE(phi));
      if(!niikcortex_initics_shrink_remesh_angular_region(obj,elen,5,1,objctr,vec,elen2,angular_thresh)) {
        fprintf(stderr,"ERROR: niikcortex_initics_shrink_remesh_angular_region\n");
        return 0;
      }
    } else {
      if(verbose>1) fprintf(stdout,"[niikcortex_initics_shrink] %4i remesh  %i / %i\n",iter,wsum,obj->nvert);
      if(!off_remesh_kobj(obj,elen,6,1)) {
        fprintf(stderr,"ERROR: off_remesh_kobj(obj,elen,10,1)\n");
        return 0;
      }
    } /* REMESH */

    /* SELF-INTERSECTION CORRECTION */
    if(!off_correct_self_intersection_with_elen(bb,obj,elen)) {
      fprintf(stderr,"ERROR: off_correct_self_intersection(bb,obj)\n");
      return 0;
    }

    if(verbose>2 && (iter%2)==0) {
      //if(verbose>=2 && (iter%10==0)) {
      if(!off_kobj_add_one_color(obj,colormap->m[iter-1][0],colormap->m[iter-1][1],colormap->m[iter-1][2])) {
        fprintf(stderr,"ERROR: off_kobj_add_one_color(obj,colormap->m[iter-1][0],colormap->m[iter-1][1],colormap->m[iter-1][2])\n");
        return 0;
      }
      sprintf(fname,"tmp_initics_shrink%i.ply",iter);
      fprintf(stdout,"[niikcortex_initics_shrink] %4i write output file: %s\n",iter,fname);
      off_kobj_write_offply(fname,obj,0);
    }

  } /* iteration */
  /* off_kobj_write_offply("tmp_initics_shrink.ply",obj,0); */

  colormap = niikmat_free(colormap);
  bb=off_bbox_free(bb);

  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_init_shrink */




double niikcortex_initics_shrink_check_deform(nifti_image *img,bbox *bb,kvert *v,double step,double thresh)
/*
 * -moves the vertex position by v->normal * step
 * -check intensity on the neighboring triangles
 * -check self-intersection
 * -if no problem, then keeps that position
 * -otherwise, vertex does not move
 * -retruns step for the actual distance moved
 *
 */
{
  niikpt
  origpt;
  origpt=v->v;
  /* go outward */
  if(niikcortex_initics_shrink_check_deform_interp(img,v,thresh)==0) {
    /* check there is surface intersection */
    if(v==g_niikcortex_initics_check_v) {
      fprintf(stdout,"\t\tmove outward\n");
    }
    v->v=niikpt_move_normal(v->v,v->normal,step);
    off_update_kvert_pminmax(v);
    if(off_check_self_intersection_kvert(bb,v))  {
      v->v=origpt;
      off_update_kvert_pminmax(v);
      return 0;
    }
    return step;
  }
  /* go inward */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"\t\tmove inward %.3f\n",thresh);
  }
  v->v=niikpt_move_normal(v->v,v->normal,-step);
  if(niikcortex_initics_shrink_check_deform_interp(img,v,thresh)==0) {
    v->v=origpt;
    return 0;
  }

  /* check there is surface intersection */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"\t\tno object\n");
  }
  off_update_kvert_pminmax(v);
  if(off_check_self_intersection_kvert(bb,v))  {
    v->v=origpt;
    off_update_kvert_pminmax(v);
    return 0;
  }
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"\t\tno triangle intersection\n");
  }

  /* check if there WILL be surface intersection */
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"\t\tno intersection\n");
  }
  v->v=niikpt_move_normal(origpt,v->normal,-NIIK_DMINMAX(5*step,1,2));
  off_update_kvert_pminmax(v);
  if(off_check_self_intersection_kvert(bb,v))  {
    v->v=origpt;
    off_update_kvert_pminmax(v);
    return 0;
  }
  if(v==g_niikcortex_initics_check_v) {
    fprintf(stdout,"\t\tno triangle intersection\n");
  }
  v->v=niikpt_move_normal(origpt,v->normal,-step);
  off_update_kvert_pminmax(v);
  return step;
} /* int niikcortex_initics_shrink_wrap_deform */


/*
 kate: space-indent on; indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
 */