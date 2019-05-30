/* Filename:     nifti1_kunio_cortex_initocs.c
 * Description:  functions for crateing initial pial surface
 * Author:       Kunio Nakamura
 * Date:         March 25, 2012
 */

#include "falcon.h"
#include "falcon_morph.h"
#include "falcon_cortex.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#endif


typedef unsigned char indicator_t;



/**********************************************************************
 *
 * niikcortex_initocs_expand
 *
 * -expands from atlas-based object to brain surface
 *
 **********************************************************************/

int niikcortex_initocs_expand_remove_intersection(bbox *bb,kobj *ics,kobj *ocs,kface *f,double dfm)
/* returns nonzero if successfully removed the intersection */
{
  const char *fcname="niikcortex_initocs_expand_remove_intersection";
  int
  iter,
  n0,nv,sn;
  niikpt
  normal,
  origpt;
  double
  maxdist,
  dist;
  int verbose=0;

  if(verbose>0) niik_fc_display(fcname,1);
  if(verbose>1) {
    if((n0 = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
      fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
      return 0;
    }
    if(!n0) return 1;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   intersections %i\n",n0);
  }
  maxdist = 5.0 * dfm;

  /* try small deformations along the normal */
  for(nv=0; nv<3; nv++) { /* each vertex */
    origpt = f->vert[nv]->v;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   vert %i [%i]\n",f->vert[nv]->index,nv);
    normal = f->vert[nv]->normal;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   normal %7.3f %7.3f %7.3f\n",normal.x,normal.y,normal.z);
    for(dist=dfm/2.0,sn=1; dist<=maxdist; dist+=dfm/2.0,sn=(sn>0)?-1:1) {
      f->vert[nv]->v = niikpt_move_normal(origpt,normal,dist*sn);
      if(off_check_self_intersection_kvert(bb,f->vert[nv])==0) {
        if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   vert %i [%i] rm dist = %7.4f\n",f->vert[nv]->index,nv,dist*sn);
        return 1;
      }
    }
    f->vert[nv]->v = origpt;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   dist did not work with v->normal\n");
  }

  /* change normal: try small deformations along the normal */
  normal = f->normal;
  for(nv=0; nv<3; nv++) { /* each vertex */
    origpt = f->vert[nv]->v;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   normal %7.3f %7.3f %7.3f\n",normal.x,normal.y,normal.z);
    for(dist=dfm/2.0,sn=1; dist<=maxdist; dist+=dfm/2.0,sn=(sn>0)?-1:1) {
      f->vert[nv]->v = niikpt_move_normal(origpt,normal,dist*sn);
      if(off_check_self_intersection_kvert(bb,f->vert[nv])==0) {
        if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   vert %i [%i] rm dist = %7.4f\n",f->vert[nv]->index,nv,dist*sn);
        return 1;
      }
    }
    f->vert[nv]->v = origpt;
    if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   dist did not work with f->normal\n");
  } /* each vertex */

  /* random normal: try small deformations along the normal */
  for(iter=0; iter<100; iter++) {
    normal = niikpt_unit(niikpt_rand());
    for(nv=0; nv<3; nv++) { /* each vertex */
      origpt = f->vert[nv]->v;
      if(verbose)fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   random normal %7.3f %7.3f %7.3f iter %i\n",normal.x,normal.y,normal.z,iter+1);
      for(dist=dfm/2.0,sn=1; dist<=maxdist; dist+=dfm/2.0,sn=(sn>0)?-1:1) {
        f->vert[nv]->v = niikpt_move_normal(origpt,normal,dist*sn);
        if(off_check_self_intersection_kvert(bb,f->vert[nv])==0) {
          if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   vert %i [%i] rm dist = %7.4f\n",f->vert[nv]->index,nv,dist*sn);
          if(verbose>1) {
            if((n0 = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
              fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
              return 0;
            }
            fprintf(stderr,"[niikcortex_initocs_expand_remove_intersection] intersections: %i\n",n0);
          }
          return 1;
        }
      }
      f->vert[nv]->v = origpt;
      if(verbose) fprintf(stdout,"[niikcortex_initocs_expand_remove_intersection]   dist did not work with random normal vtx %i %i\n",f->vert[nv]->index,nv);
    } /* each vertex */
  } /* iteration */

  if(verbose>1) {
    if((n0 = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
      fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
      return 0;
    }
    fprintf(stderr,"[niikcortex_initocs_expand_remove_intersection] intersections: %i\n",n0);
  }
  if(verbose>0) niik_fc_display(fcname,0);
  return 0;
}

int niikcortex_initocs_expand(nifti_image *img,          /* t1w image */
                              nifti_image *brain_mask,   /* brain mask */
                              double *step_size,         /* variable step size */
                              int stepiter,              /* size of step_size vector */
                              double intGM,              /* gray matter intensity */
                              double intWM,              /* white matter intensity */
                              double intCSF,             /* cerebrospinal fluid intensity */
                              double intICS,             /* white matter surface intensity */
                              double intOCS,             /* pial surface intensity */
                              double gthresh,            /* gradient threshold */
                              double max_cth,            /* maximum allowed cortical thickness */
                              double init_cth,           /* initial cortical thickness */
                              int smooth_iter_maxcth,    /* smooth number of iterations (for max_cth) */
                              int smooth_iter,           /* smooth number of iterations for each step */
                              int maxiter2,              /* maximum number of iteration */
                              float smooth_percent,      /* smooth percent for smooth_iter */
                              kobj *ics,                 /* white matter surface */
                              kobj *ocs,                 /* pial surface is updated */
                              nifti_image *border)       /* limit for GM expansion*/
{              

  const char *fcname=__func__;
  char fname[512];
  kobj **ctx;
  bbox *bb=NULL;
  niikpt origpt;
  kvert
    *vi,*vo;
  kface
    *fi,*fo;
  kvert **_vi,**_vo;
  double
    *thk_list,
    dist,sumdist,
    cth,
    xx,
    *x,*y,*dy,dx,xmin,xmax,
    *mcth_list;
  int *yb;
  int
    iter,
    iter2,
    vindex,
    nx,n0,
    m,n,nn,num,
    xsc,
    verbose=niik_verbose();
  int i,c;
  struct tm *stm;
  time_t ctm;

  const int stride=5;
  indicator_t *indicator; /*indicator helper function, using in paralell apply deformation only*/
  char tmstr[256];
  cortex_tracing_info trace;
  int debug_tracing=falcon_tracing_init(img,&trace);

  if(verbose>=1) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((ics==NULL),fcname,"ics is null");
  NIIK_RET0((ocs==NULL),fcname,"ocs is null");
  NIIK_RET0((brain_mask==NULL),fcname,"brain_mask is null");

  if(ics->nvert != ocs->nvert) {
    fprintf(stderr,"[%s] ERROR: #vertex is different %i %i\n",fcname,ics->nvert,ocs->nvert);
    return 0;
  }

  if(verbose) {
    fprintf(stdout,"[%s] parameters\n",fcname);
    fprintf(stdout,"    image         %s\n",img->fname);
    fprintf(stdout,"    brainmask     %s\n",brain_mask->fname);
    fprintf(stdout,"    step          ");
    for(n=0; n<stepiter; n++) {
      fprintf(stdout,"%7.3f ",step_size[n]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"    max iter      %-i / step\n",maxiter2);
    fprintf(stdout,"    g-thresh      %-7.3f\n",gthresh);
    fprintf(stdout,"    smooth iter maxcth  %-i\n",smooth_iter_maxcth);
    fprintf(stdout,"    smooth iter   %-i\n",smooth_iter);
    fprintf(stdout,"    init thick    %-7.3f\n",init_cth);
    fprintf(stdout,"    max thick     %-7.3f\n",max_cth);
    fprintf(stdout,"    obj vfe       %-i %-i %-i\n",ics->nvert,ics->nface,ics->nedge);
  }

  NIIK_RET0(((indicator = (indicator_t*)malloc(ocs->nvert*sizeof(indicator_t)))==NULL),fcname,"malloc");
  NIIK_RET0(((_vi = (kvert**)calloc(ocs->nvert,sizeof(kvert*)))==NULL),fcname,"calloc");
  NIIK_RET0(((_vo = (kvert**)calloc(ocs->nvert,sizeof(kvert*)))==NULL),fcname,"calloc");
  NIIK_RET0(((ctx=(kobj **)calloc(2,sizeof(kobj *)))==NULL),fcname,"calloc");

  ctx[0] = ics;
  ctx[1] = ocs;
  NIIK_RET0(((bb=off_bbox_init(7,320))==NULL),fcname,"off_bbox_init");
  if(verbose>=1) fprintf(stdout,"[%s] bbox delta = %.3f\n",fcname,bb->delta);

  NIIK_RET0(((thk_list=(double *)calloc(ics->nvert,sizeof(double)))==NULL),fcname,"calloc");

  if(verbose) fprintf(stdout,"[%s] check white surface intersection\n",fcname);
  if((xsc=off_count_self_intersection_add_color(bb,ics,1))<0) {
    fprintf(stderr,"[%s] ERROR: off_count_self_intersection_add_color(bb,ics,1)\n",fcname);
    return 0;
  }
  if(xsc) {
    fprintf(stderr,"[%s] ERROR: ics has surface intersections %i\n",fcname,xsc);
    return 0;
    NIIK_RET0((!off_correct_self_intersection(bb,ics)),fcname,"off_correct_self_intersection");
  } else {
    fprintf(stderr,"[%s]   no surface interaction %i\n",fcname,xsc);
  }/* check intersection */

  if(verbose>0) fprintf(stdout,"[%s] initialize cortex gradient vector\n",fcname);
  xmin = -3;
  xmax = 12;
  dx = 0.1;
  /*number of steps (?)*/
  for(nx=0,xx=xmin; xx<=(xmax+0.001); xx+=dx) {
    nx++;
  }
  if(verbose>0) fprintf(stdout,"[%s] vector %6.2f : %6.2f : %6.2f  | %i\n",fcname,xmin,dx,xmax,nx);
  x = niik_calloc_double_vector(nx);
  y = niik_calloc_double_vector(nx);
  dy= niik_calloc_double_vector(nx);
  yb= (int*)calloc(nx,sizeof(int));
  
  for(n=0; n<nx; n++) {
    x[n] = xmin + dx*n;
    yb[n] = 0;
  }
  n0=(int) floor(-xmin/dx+0.5);
  if(verbose>0) fprintf(stdout,"[%s] vector %6.2f : %6.2f : %6.2f  | %i n0 %i\n",fcname,xmin,dx,xmax,nx,n0);


  /* overall processing:
   * 1. initial outer cortical surface
   * 1.a. move surface outward by small amount
   * 1.b. if self-interesction is found,
   *    try moving in random directions (small steps)
   * 1.c. repeat until no self-intersection
   * 2. estimate the maximum cortical thickness
   *    -brain mask
   *    -max length 3mm
   *    -self-intersection
   *    -intensity
   *    -first derivative
   * 2.b. smooth the max cortical thickness
   * 3. apply balloon
   * 3.a. expand with small steps
   * 3.b. smooth (not too much)
   * 3.c. repeat
   */


  /******************************************
   *
   * INITIAL PIAL SURFACE
   *
   ******************************************/

  off_update_kobj_face_normal(ics);
  off_update_kobj_vert_normal(ics);

  /*off_smooth_kobj_vert_normal(ics);*/
  for(fi=ics->face,fo=ocs->face; fi!=NULL; fi=fi->next,fo=fo->next) {
    fo->normal=fi->normal;
  }
  for(vi=ics->vert,vo=ocs->vert; vi!=NULL; vi=vi->next,vo=vo->next) {
    vo->normal=vi->normal;
  }

  /*create index*/
  for(vi=ics->vert,vo=ocs->vert,vindex=0; vi!=NULL; vi=vi->next,vo=vo->next,vindex++) {
    _vi[vindex]=vi;
    _vo[vindex]=vo;
    vo->idata=1;
    vi->idata=0;
  }

  if(verbose>=1) fprintf(stdout,"[%s] initialize pial surface at white surface\n",fcname);
  for(vi=ics->vert,vo=ocs->vert; vi!=NULL; vi=vi->next,vo=vo->next) {
    vo->v=niikpt_move_normal(vi->v,vi->normal,init_cth);
  }

  if(verbose) fprintf(stdout,"[%s] check cortex intersections\n",fcname);
  if((xsc = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
    fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
    return 0;
  }

  if(verbose>1 && xsc>0) {
    if(verbose>2) {
      fprintf(stdout,"[niikcortex_initocs_expand] write tmp_niikcortex_initocs_ics.off\n");
      off_kobj_write_offply("tmp_niikcortex_initocs_ics.off",ics,0);
      fprintf(stdout,"[niikcortex_initocs_expand] write tmp_niikcortex_initocs_ocs.off\n");
      off_kobj_write_offply("tmp_niikcortex_initocs_ocs.off",ocs,0);
    }
    for(fo=ocs->face; fo!=NULL; fo=fo->next) {
      if(fo->color[0]==fo->color[1]) continue;
      fprintf(stdout,"\nface %8i   normal = %7.2f %7.2f %7.2f \n",fo->index,fo->normal.x,fo->normal.y,fo->normal.z);
      for(m=0; m<3; m++) {
        fprintf(stdout,  "vert %8i   vertex = %7.2f %7.2f %7.2f \n",fo->vert[m]->index,fo->vert[m]->v.x,fo->vert[m]->v.y,fo->vert[m]->v.z);
      }
    }
  } /* xsc write output */

  if(xsc) { /* correct intersections */
    nn = NIIK_IMAX(100,xsc);
    for(iter=0; iter<nn && xsc>0; iter++) {
      fprintf(stdout,"[%s] removing intersection iter %i\n",fcname,iter);
      /* correct each surface separately */
      for(fo=ocs->face; fo!=NULL; fo=fo->next) {
        if(fo->color[0]==fo->color[1]) continue;
        break;
      }
      if(fo==NULL) break;
      while(fo!=NULL) {
        if(niikcortex_initocs_expand_remove_intersection(bb,ics,ocs,fo,init_cth*2)) break;
        /* do the next one */
        for(fo=fo->next; fo!=NULL; fo=fo->next) {
          if(fo->color[0]==fo->color[1]) continue;
          break;
        }
        if(fo==NULL) break;
        niikcortex_initocs_expand_remove_intersection(bb,ics,ocs,fo,init_cth*2);
      }
      if(iter>=5) { /* after 5 iterations, try to change white surface in addition to pial surface */
        fprintf(stdout,"[niikcortex_initocs_expand] removing intersection iter %i [white]\n",iter);
        /* correct each surface separately */
        for(fi=ics->face; fi!=NULL; fi=fi->next) {
          if(fi->color[0]==fi->color[1]) continue;
          off_remesh_kvert_relocate(fi->vert[0]);
          off_remesh_kvert_relocate(fi->vert[1]);
          off_remesh_kvert_relocate(fi->vert[2]);
        }
      } /* white surface perturbation */
      if((xsc = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
        fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
        return 0;
      }
      fprintf(stdout,"[%s]   remaining intersections [iter %i] %i\n",fcname,iter,xsc);
    } /* each iteration */
  } /* xsc */
  if(verbose>0) fprintf(stdout,"[%s] after initial surface, check intersections = %i\n",fcname,xsc);

  if(verbose>0) fprintf(stdout,"[%s] regional max cortical thickness\n",fcname);


  /*
   * CALCULATE MAXIMUM CORTICAL THICKNESS
   * per vertex, by sampling intensities in the normal direction
   */
  NIIK_RET0(((mcth_list = niik_calloc_double_vector(ics->nvert))==NULL),fcname,"niik_calloc_double_vector");

  for(vi=ics->vert,vo=ocs->vert,vindex=0; vi!=NULL; vi=vi->next,vo=vo->next,vindex++) {
    for(n=0; n<nx; n++) {
      niikpt _pos=niikpt_move_normal(vo->v, vi->normal, x[n] );
      dy[n] = y[n] = niik_image_interpolate_3d_linear(img,_pos);
     
      if(border)
        yb[n] = niik_image_interpolate_3d_nn(border,_pos);
    }
    niik_runavg_double_vector(dy,nx,nx/20);
    niik_central_difference_double_vector(dy,nx);
    for(n=0; n<nx; n++) {
      dy[n] /= dx;
    }

    if(verbose>5) {
      for(n=0; n<nx; n++) {
        fprintf(stdout,"%3i %6.2f %6.2f %8.3f \n",n,x[n],y[n],dy[n]);
      }
    }

    /* find the very max */
    for(n=n0; n<nx; n++) {
      if(x[n]>=max_cth) {
        if(verbose>5) {
          fprintf(stdout,"\tmax cth %.5f\n",max_cth);
        }
        break;
      }
      if(y[n]<intCSF || yb[n]>0)   {
        if(verbose>5) {
          fprintf(stdout,"\tCSF %.2f < %.2f\n",y[n],intCSF);
        }
        break;
      }
      if(y[n]>intOCS && y[n]<intGM && dy[n]>gthresh) {
        if(verbose>5) {
          fprintf(stdout,"\ty=%.2f dy=%.2f\n",y[n],dy[n]);
        }
        break;
      }
      if(x[n]>=3 && y[n]>=intWM) {
        if(verbose>5) {
          fprintf(stdout,"\tx=%5.2f\n",x[n]);
        }
        break;
      }
    }
    if(verbose>5) {
      fprintf(stdout,"%3i %6.2f %6.2f %8.3f   max\n",n,x[n],y[n],dy[n]);
    }

    /* pull back if really low */
    if(y[n]<intCSF||yb[n]>0) {
      for(m=n; m>n0; m--) {
        if(y[m]>intCSF && y[n]<intGM) {
          if(m==(nx-1)) {
            continue;
          }
          if(dy[m]<dy[m-1] && dy[m]<dy[m+1])
            break;
        }
      }
      if(m>n0) n=m;
    }
    /* pull back if gradient is not the local max */
    else if(dy[n]>gthresh) {
      for(m=n; m>n0; m--) {
        if(y[m]>intGM) continue;
        if(dy[m]<=0) break;
      }
      if(m>n0) n=m;
    }

    mcth_list[vindex] = NIIK_DMIN(x[n],max_cth);
    if(verbose>5) fprintf(stdout,"%3i %6.2f %6.2f %8.3f \n",n,x[n],y[n],dy[n]);
  } /* each vertex */
  free(x);
  free(y);
  free(dy);
  free(yb);

  if(verbose>0) fprintf(stdout,"[%s] regional max cortical thickness done\n",fcname);


  /* SMOOTHING -MAX CORTICAL THICKNESS */
  if(verbose) {
    fprintf(stdout,"[niikcortex_initocs_expand] surface smoothing (regional max cortical thickness) %i\n",smooth_iter_maxcth);
    /*niik_display_stats_for_double_vector(mcth_list,ics->nvert);*/
  }
  for(n=0; n<smooth_iter_maxcth; n++) {
    NIIK_RET0((!off_surface_smooth_using_vert(ics,mcth_list,3,1)),fcname,"off_surface_smooth_using_vert");
  }
  if(verbose) {
    fprintf(stdout,"[niikcortex_initocs_expand] after surface smoothing\n");
    niik_display_stats_for_double_vector(mcth_list,ics->nvert);
  }


  /************************************************
   *
   * APPLY BALLOON
   *   main loop
   *
   ************************************************/

  if(verbose>1) fprintf(stdout,"[%s] apply balloon\n",fcname);

  for(iter=0; iter<stepiter; iter++) {
    if(verbose>2) {
      sprintf(fname,"tmp_initocs_%i.off",iter);
      fprintf(stdout,"[%s] iter %i  writing output for testing %s\n",fcname,iter,fname);
      NIIK_RET0((!off_kobj_add_one_color(ocs,1,0.2,0.2)),fcname,"off_kobj_add_one_color");
      NIIK_RET0((!off_kobj_write_offply(fname,ocs,0)),fcname,"off_kobj_write_off");
    }

    /* UPDATE THICKNESS */
    for(vo=ocs->vert,vi=ics->vert,vindex=0; vo!=NULL; vo=vo->next,vi=vi->next,vindex++) {
      thk_list[vindex] = niikpt_distance(vi->v,vo->v);
    }
    if(verbose>=1)
      fprintf(stdout,"[%s] iter %i step %6.3f  cth = %9.5f\n",fcname,iter,step_size[iter],
              niik_get_mean_from_double_vector(thk_list,ics->nvert));

    /* reset the w-flags to deform */
    for(vo=ocs->vert; vo!=NULL; vo=vo->next) {
      vo->v.w=0;
    }

    /**********************************************************
     *
     * BALLOON
     *
     **********************************************************/

    for(dist=0, sumdist=bb->delta, iter2=0;
        ( dist <= (max_cth+0.001) ) && iter2 < maxiter2 ;
        iter2++, dist+=step_size[iter], sumdist+=step_size[iter] ) {
      if(sumdist>bb->delta*0.5) {
        if(verbose>2) fprintf(stdout,"[%s] %i  bbox update\n",fcname,iter);

        #pragma omp parallel for private(n)
        for(n=0; n<2; n++) {
          off_update_kobj_face_normal(ctx[n]);
          off_update_kobj_vert_normal(ctx[n]);
          off_smooth_kobj_vert_normal(ctx[n]);
          off_update_kobj_kface_pminmax(ctx[n]);
        }
        NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,ctx,2)),
                  fcname,"off_create_bbox_from_multiple_kobj");
        sumdist=0;
      }  /* bounding box */
      else {
        if(verbose>2) fprintf(stdout,"[%s] %i  no bbox update  %.3f\n",fcname,iter,sumdist);
      }

      ctm=time(NULL);
      stm=localtime(&ctm);
      strftime(tmstr,256,"%Y-%m-%d %T",stm);
      if(verbose>2) fprintf(stdout,"[%s] deformation %i %s\n",fcname,iter2,tmstr);


      /**/
      memset(indicator,0,ocs->nvert);

      num=0;
      for(c=0; c<stride; c++) {
        /*distribute tasks in a predictable fashion*/
        #pragma omp parallel for reduction(+:num) private(i) shared(bb,indicator,_vi,_vo)
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
                if(cf->vert[0]->idata != 1) continue; /*wrong cortex*/
                for(l=0; l<3; l++) {
                  kvert *vi,*vo;
                  int vindex;
                  double cth;
                  niikpt normal, origpt;

                  vo=cf->vert[l];
                  vindex=vo->index-1;
                  vi=_vi[vindex];

                  if(indicator[vindex]) continue;
                  indicator[vindex]=1;

                  if(vo->v.w>0) continue;

                  origpt = vo->v;
                  vo->v = niikpt_move_normal(vo->v,vo->normal,step_size[iter]);
                  cth = niikpt_distance(vi->v,vo->v);

                  /* check for the max cortical thickness */
                  if(cth>mcth_list[vindex]) {
                    vo->v = origpt;
                    vo->v.w = 1;
                    continue;
                  }

                  if(cth>max_cth) {
                    vo->v = origpt;
                    vo->v.w = 1;
                    continue;
                  }

                  /* check for brain mask */
                  if(niik_image_interpolate_3d_nn(brain_mask,vo->v)==0) {
                    vo->v = origpt;
                    vo->v.w = 1;
                    continue;
                  }

                  /* check for border, if present*/
                  if(border)
                  {
                    if(niik_image_interpolate_3d_nn(border,vo->v)>0) {
                      vo->v = origpt;
                      vo->v.w = 1;
                      continue;
                    }
                  }

                  /* check for intersection 0.5mm */
                  vo->v = niikpt_move_normal(origpt,vo->normal,0.5);
                  off_update_kvert_pminmax(vo);
                  if(off_check_self_intersection_kvert(bb,vo)) {
                    vo->v=origpt;
                    off_update_kvert_pminmax(vo);
                    continue;
                  }

                  /* check for intersection defined step */
                  vo->v = niikpt_move_normal(origpt,vo->normal,step_size[iter]);
                  off_update_kvert_pminmax(vo);
                  if(off_check_self_intersection_kvert(bb,vo)) {
                    vo->v=origpt;
                    off_update_kvert_pminmax(vo);
                    continue;
                  }
                  num++;
                }
              }
            }
          } /* for each vertex */
        }
      }

      /***********************************************************
       *
       * SMOOTHING
       *
       ***********************************************************/

      /* if(verbose>=2) fprintf(stdout,"[%s]   bbox update\n",fcname);
         #pragma omp parallel for
         for(n=0;n<2;n++) {
         off_update_kobj_face_normal(ctx[n]);
         off_update_kobj_vert_normal(ctx[n]);
         off_smooth_kobj_vert_normal(ctx[n]);
         off_update_kobj_kface_pminmax(ctx[n]); }
         NIIK_RET0((!off_create_bbox_from_multiple_kobj(bb,ctx,2)),fcname,"off_create_bbox_from_multiple_kobj"); */

      if(smooth_percent>0.0) {
        for(n=0; n<smooth_iter; n++) {
          if(verbose>1) {
            for(vo=ocs->vert,vi=ics->vert,vindex=0; vo!=NULL; vo=vo->next,vi=vi->next,vindex++) {
              thk_list[vindex] = niikpt_distance(vi->v,vo->v);
            }
            fprintf(stdout,"[%s]   smooth %i.%i   cth = %6.4f (%6.4f %6.4f)\n",
                    fcname,iter,n,
                    niik_get_mean_from_double_vector(thk_list,ics->nvert),
                    niik_get_min_from_double_vector(thk_list,ics->nvert),
                    niik_get_max_from_double_vector(thk_list,ics->nvert));
          }

          for(vo=ocs->vert,vi=ics->vert,vindex=0; vo!=NULL; vo=vo->next,vi=vi->next,vindex++) {
            origpt = vo->v;
            vo->v = niikpt_wavg(vo->v, niikpt_kvert_local_average2(vo), smooth_percent);
            cth   = niikpt_distance(vi->v,vo->v);
            if(cth>max_cth) {
              vo->v=origpt;
              continue;
            }
            off_update_kvert_pminmax(vo);
            if(off_check_self_intersection_kvert(bb,vo)) {
              vo->v=origpt;
              off_update_kvert_pminmax(vo);
              continue;
            }
          } /* FOR EACH VERTEX */
        } /* SMOOTH ITERATION */
      }

      if(verbose>=2) {
        for(vo=ocs->vert,vi=ics->vert,vindex=0; vo!=NULL; vo=vo->next,vi=vi->next,vindex++) {
          thk_list[vindex] = niikpt_distance(vi->v,vo->v);
        }

        fprintf(stdout,"[%s] iter %i.%i step %6.3f cumstep %6.3f cth = %6.4f (%6.4f %6.4f)  %i\n",
                fcname,iter,iter2,
                step_size[iter],dist,
                niik_get_mean_from_double_vector(thk_list,ics->nvert),
                niik_get_min_from_double_vector (thk_list,ics->nvert),
                niik_get_max_from_double_vector (thk_list,ics->nvert),
                num);
      }

      if(debug_tracing)
      {
        NIIK_RET0((!off_kobj_add_one_color(ics,1,1,0)),fcname,"off_kobj_add_one_color"); /* yellow for white matter-surface */
        NIIK_RET0((!off_kobj_add_one_color(ocs,1,0,0)),fcname,"off_kobj_add_one_color"); /* red for pial-surface */

        falcon_tracing_dump(&trace,iter2+iter*maxiter2,"initocs",img,bb);
        falcon_tracing_dump_objects(&trace, iter2+iter*maxiter2 ,"initocs", &ocs, 1 );
      }
      /* can be removed after enough debugging */
      if(verbose>1) fprintf(stderr,"[%s]   test intersections\n",fcname);
      if((n = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
        fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
        return 0;
      }
      if(n) {
        fprintf(stderr,"[%s] ERROR: intersections %i\n",fcname,n);
      }
      /* up to here */

      if(!num) {
        if(verbose>1) fprintf(stdout,"[%s] iter %i  no motion\n",fcname,iter);
        break;
      } else if(100.0*num/ocs->nvert < 0.1) {
        if(verbose>1) fprintf(stdout,"[%s] iter %i  no motion (#moving vertex [%i] < 0.1%%)\n",fcname,iter,num);
        break;
      }

    } /* repeat the balloon until no change */

    /* can be removed after enough debugging */
    if((n = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
      fprintf(stderr,"ERROR: niikcortex_off_count_intersection(bb,ics,ocs)\n");
      return 0;
    }
    if(n) {
      fprintf(stderr,"[%s] ERROR: intersections %i\n",fcname,n);
    }
    /* up to here */
  } /* iter, different step sizes */


  if((n = niikcortex_off_count_intersection(bb,ics,ocs))<0) {
    fprintf(stderr,"[%s] ERROR: niikcortex_off_count_intersection\n",fcname);
    return 0;
  }
  if(n) {
    fprintf(stderr,"[%s] ERROR: found intersections %i\n",fcname,n);
  }

  if(verbose>=1) fprintf(stdout,"[%s] add red color\n",fcname);
  if(!off_kobj_add_one_color(ocs,1,0,0)) {
    fprintf(stderr,"[%s] ERROR: off_kobj_add_one_color(ocs,1,0,0)\n",fcname);
    return 0;
  }

  bb=off_bbox_free(bb);
  free(mcth_list);
  free(ctx);
  free(thk_list);
  free(indicator);
  free(_vi);
  free(_vo);
  falcon_tracing_free(&trace);
  
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niikcortex_initocs_expand */


/*
 kate: space-indent on; indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/