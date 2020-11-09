/* Filename:     nifti1_kunio_off_remesh.c
* Description:  object functions for remeshing
* Author:       Kunio Nakamura
* Date:         March 1, 2012
*
* geomview's off file
* -remeshing ( collapse_edge / split_edge / flip_edge / relocate )
*
* -this function is very difficult to validate
*
*/

#include "falcon.h"
#include "falcon_surfaces.h"

int off_kobj_test_common_neighbor(kedge *edge);

/****************************
 * reinitialize element ids
 * 
 * 
 * **************************/
void  off_obj_renumber(kobj *obj)
{
  kedge *e;
  kface *f;
  kvert *v;
  int n;

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
}


/**********************************************
* remesh function
*
* -remeshes the object so that the edge length is close to elen
*
* -elen is the target edge length
* -maxiter should be 5-10
* -wmark, if nonzero, will do selective remeshing
*  where v->v.w is used to decide where we would remsh
*  so that if v->v.w is nonzero, function skips remeshing
*  and     if v->v.w is zero, function remeshes
* -a point is usually initialized with p.w = 0
*  so by default, remeshing occurs
*
*
**********************************************/
int off_remesh_kobj(kobj *obj,double elen,int maxiter,int wmark) {
  return off_remesh_kobj_ex(obj,4.0/5*elen,4.0/3*elen,maxiter,wmark,1);
}

/**********************************************
* remesh function
*
* -remeshes the object so that the edge length is close to elen
*
* -emin is the minimum length after which edge will be collapsed
* -emaxn is the maximim length after which the edge will be split in tow
* -maxiter should be 5-10
* -wmark, if nonzero, will do selective remeshing
*  where v->v.w is used to decide where we would remsh
*  so that if v->v.w is nonzero, function skips remeshing
*  and     if v->v.w is zero, function remeshes
* -a point is usually initialized with p.w = 0
*  so by default, remeshing occurs
* -smooth - do smoothing after remeshing
*
*
**********************************************/
int off_remesh_kobj_ex(kobj *obj,double emin,double emax, int maxiter,int wmark,int smooth) 
{
  const char *fcname=__func__;

  double dval;
  kedge *e,*ne,*erm[9];
  kface *f,*frm[9];
  kvert *v,*vrm[2];
  int n,iter;
  const int verbose=0;

  /*emax=4.0/3*elen;
  emin=4.0/5*elen;*/
  if(fabs(emin)<0.001 || fabs(emax)<0.001||emax<emin)
    abort();

  if(verbose>=1) {
    fprintf(stdout,"\n[%s] start: %8.4f %8.4f \n",fcname,emin,emax);
    fprintf(stdout,"[%s]   vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);
  }

  for(iter=1; iter<=maxiter; iter++) {
    if(verbose>=1) {
      fprintf(stdout,"[%s] remesh iteration %3i\n",fcname,iter);
      fprintf(stdout,"[%s]   checking object\n",fcname);
      if(!off_kobj_test(obj)) {
        fprintf(stderr,"ERROR: off_kobj_test during split\n");
        exit(0);
      }
    }

    /***************************
    * SPLIT EDGES
    ****************************/
    if(verbose>=1) {
      fprintf(stdout,"[%s] split edges  %8.4f \n",fcname,emax);
    }

    for(e=obj->edge; e!=NULL; e=e->next) {
      if(wmark)
        if(e->endpts[0]->v.w>0 || e->endpts[1]->v.w>0)
          continue;

      dval = niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
      /* fprintf(stdout,"\t\t%i  %8.5f\n",e->index,dval); */
      if(dval<emax) continue;
      if((v=off_remesh_kedge_split(obj,e))==NULL) {
        fprintf(stderr,"ERROR: off_remesh_kedge_split %i\n",e->index);
        return 0;
      }

      if(verbose>=2) {
        fprintf(stdout,"[%s]   testing obj e = %i\n",fcname,e->index);
        if(!off_kobj_test(obj)) {
          fprintf(stderr,"ERROR: off_kobj_test during split\n");
          exit(0);
        }
      }
    } /* edge split */


    if(verbose>=1) {
      /* check the status */
      off_obj_renumber(obj);
      fprintf(stdout,"[%s]   vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);
      fprintf(stdout,"[%s]   mean elen = %9.5f \n",fcname,off_get_kobj_mean_edge_length(obj));
      if(verbose>=2) {
        fprintf(stdout,"[%s]    test obj after split\n",fcname);
        if(!off_kobj_test(obj)) {
          fprintf(stderr,"ERROR: off_kobj_test\n");
          exit(0);
        }
      }
    }


    /***************************
    * COLLAPSE EDGES
    ****************************/

    if(verbose>=1) {
      fprintf(stdout,"[%s] collapse edges  %8.4f \n",fcname,emin);
    }
    for(e=obj->edge; e!=NULL; e=ne) {
      ne=e->next;

      if(wmark) {
        /* 2012-06-26, Kunio
        * changed from && to || */
        if(e->endpts[0]->v.w>0 || e->endpts[1]->v.w>0)
          continue;
      }

      /* skip if long edge */
      dval=niikpt_distance(e->endpts[0]->v,e->endpts[1]->v);
      if(dval>emin) continue;

      if(verbose>=2) fprintf(stdout,"[%s]    collapse %i  %8.5f\n",fcname,e->index,dval);
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
      if(verbose>=4) {
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
      if(verbose>=4) fprintf(stdout,"\t\tcollapse clean up !\n");
      /* remove and free */
      off_remesh_kedge_collapse_clean_up(obj,vrm,frm,erm);
      vrm[0]=NULL;
      frm[0]=frm[1]=NULL;
      erm[0]=erm[1]=erm[2]=NULL;
      if(verbose>=4) fprintf(stdout,"\t\tcollapse clean up done\n");

      /* remove this after debugging
      * -check the object for everything at each step  */
      if(verbose>=3) {
        if(!off_kobj_test(obj)) {
          fprintf(stderr,"ERROR: off_kedge_test\n");
          fprintf(stderr,"       test off_remesh_kedge_collapse\n");
          exit(0);
        }
      }
    } /* each edge for collapsing */

    /* correct 3 neighbors */
    if(!off_remesh_kobj_collapse_correction(obj)) {
      fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction \n");
      return 0;
    }

    if(verbose>=2) {
      off_obj_renumber(obj);
      fprintf(stdout,"      vfe = %i %i %i \n",obj->nvert,obj->nface,obj->nedge);
      fprintf(stdout,"      mean elen = %9.5f \n",off_get_kobj_mean_edge_length(obj));
      if(!off_kobj_test(obj)) {
        fprintf(stderr,"ERROR: test obj\n");
        exit(0);
      }
    }

    /***************************
    * VALENCE MINIMIATION
    ****************************/
    if(verbose>=2) {
      fprintf(stdout,"\t  valence6  = %i\n",off_remesh_kobj_count_global_valence6(obj));
      fprintf(stdout,"\t  valence57 = %i\n",off_remesh_kobj_count_global_valence57(obj));
      fprintf(stdout,"\t  valence   = %i  variance\n",off_remesh_kobj_count_global_valence_variance(obj));
    }


    if(verbose>=1) fprintf(stdout,"[%s] minimize the valence\n",fcname);
    for(n=0; n<3; n++) {
      if(!off_remesh_kobj_valence_minimization(obj,wmark)) {
        fprintf(stderr,"ERROR: off_remesh_kobj_valence_minimization\n");
        return 0;
      }
    }

    if(verbose>=1) {
      n = off_remesh_kobj_count_global_valence57(obj);
      off_kobj_update_num(obj);
      fprintf(stdout,"[%s]   valence = %3.0f%% = %i / %i\n",fcname,100.0*n/obj->nvert,n,obj->nvert);

      if(verbose>=2) {
        fprintf(stdout,"    check status and update index \n");
        off_obj_renumber(obj);
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
    if(smooth) 
    {
      if(verbose>=1) fprintf(stdout,"[%s] relocation: tangential relaxation\n",fcname);
      for(n=0; n<3; n++) {
        if(wmark) {
          for(v=obj->vert; v!=NULL; v=v->next) {
            if(v->v.w==0) off_remesh_kvert_relocate(v);
            if(v->next==NULL) break;
          }
          for(v=v; v!=NULL; v=v->prev) {
            if(v->v.w==0) off_remesh_kvert_relocate(v);
          }
        } else {
          for(v=obj->vert; v!=NULL; v=v->next) off_remesh_kvert_relocate(v);
        }
      }

      if(verbose>=2) {
        fprintf(stdout,"    check status and update index \n");
        off_obj_renumber(obj);

        fprintf(stdout,"  vfe = %i %i %i \n",obj->nvert,obj->nface,obj->nedge);
        fprintf(stdout,"    mean elen = %9.5f \n",off_get_kobj_mean_edge_length(obj));

        if(!off_kobj_test(obj)) {
          fprintf(stderr,"ERROR: test obj\n");
          exit(0);
        }
      }
      if(obj->spherecoo==1) off_avg_smooth_spherical_coordinates(obj);
    }

    if(verbose>=1) fprintf(stdout,"[%s] remesh vfe = %i %i %i \n",fcname,obj->nvert,obj->nface,obj->nedge);
  } /* iteration */

  /* TESTING FOR FACE ORIENTATION
  * -I tested for the normal direction using a large sphere
  * -This way I can test the surface orientation by taking the surface normal
  * -I can copy/paste the below loop to right after split and collapse
  *  to test each routine. --but not after flip as it can cause opposite
  *  orientation in a valid manner
  for(f=obj->face,n=1;f!=NULL;f=f->next){
    off_update_kface_normal(f);
    if(niikpt_dot(niikpt_unit(f->vert[0]->v),f->normal)<-0.5) {
      fprintf(stdout,"%9i: %12.6f | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f\n",
              f->index,niikpt_dot(niikpt_unit(f->vert[0]->v),f->normal),
              f->normal.x,f->normal.y,f->normal.z,
              f->vert[0]->v.x,f->vert[0]->v.y,f->vert[0]->v.z);
      exit(0); } }*/

  off_obj_renumber(obj);
  if(verbose>=1) fprintf(stdout,"[%s] remeshed to vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);

  if(obj->color) {
    if(verbose>=1) fprintf(stdout,"[%s] adding color (memory)\n",fcname);
    if(!off_kobj_add_color(obj)) {
      fprintf(stderr,"ERROR: off_kobj_add_color\n");
      return 0;
    }
  }

  if(verbose) fprintf(stdout,"[%s] finished\n",fcname);
  return 1;
} /* int off_remesh_kobj_ex */


int off_relax_kobj(kobj *obj, double elen,int maxiter,int wmark) {
  const char *fcname=__func__;
  double emin,emax,dval;
  kedge *e,*ne,*erm[9];
  kface *f,*frm[9];
  kvert *v,*vrm[2];
  int n,iter;
  const int verbose=0;

  emax=4.0/3*elen;
  emin=4.0/5*elen;
  if(fabs(elen)<0.001)
    abort();

  if(verbose>=1) {
    fprintf(stdout,"\n[%s] start: %8.4f %8.4f %8.4f \n",fcname,elen,emin,emax);
    fprintf(stdout,"[%s]   vfe = %i %i %i\n",fcname,obj->nvert,obj->nface,obj->nedge);
  }

  for(iter=1; iter<=maxiter; iter++) {

    /***************************
    * TANGENTIAL RELAXATION
    ****************************/
    if(verbose>=1) fprintf(stdout,"[%s] relocation: tangential relaxation\n",fcname);
    for(n=0; n<3; n++) {
      if(wmark) {
        for(v=obj->vert; v!=NULL; v=v->next) {
          if(v->v.w==0) off_remesh_kvert_relocate(v);
          if(v->next==NULL) break;
        }
        for(v=v; v!=NULL; v=v->prev) {
          if(v->v.w==0) off_remesh_kvert_relocate(v);
        }
      } else {
        for(v=obj->vert; v!=NULL; v=v->next)
          off_remesh_kvert_relocate(v);
      }
    }

    if(obj->spherecoo==1) off_avg_smooth_spherical_coordinates(obj);

    if(verbose>=1) fprintf(stdout,"[%s] relax vfe = %i %i %i \n",fcname,obj->nvert,obj->nface,obj->nedge);
  } /* iteration */

  if(verbose) fprintf(stdout,"[%s] finished\n",fcname);
  return 1;
} /* int off_remesh_kobj(kobj *obj,double elen) */


/******************************************************
*
* functions related to off_remesh_kedge_collapse
*
*   off_remesh_kedge_collapse_check_kedge
*   off_remesh_kedge_collapse_clean_up
*   off_remesh_kedge_collapse2
*
*
*
******************************************************/
int off_remesh_kedge_collapse_clean_up(kobj *obj,kvert **vrm,kface **frm,kedge **erm) {
  const int verbose=0;
  if(verbose) fprintf(stdout,"\tremove v %i\n",vrm[0]->index);
  off_kvert_remove(vrm[0],obj);
  if(verbose) fprintf(stdout,"\tfree   v\n");
  off_kvert_free(vrm[0]);
  if(verbose) fprintf(stdout,"\tremove f %i %i\n",frm[0]->index,frm[1]->index);
  off_kface_remove(frm[0],obj);
  off_kface_remove(frm[1],obj);
  if(verbose) fprintf(stdout,"\tfree   f\n");
  off_kface_free(frm[0]);
  off_kface_free(frm[1]);
  if(verbose) fprintf(stdout,"\tremove e %i %i %i\n",erm[0]->index,erm[1]->index,erm[2]->index);
  off_kedge_remove(erm[0],obj);
  off_kedge_remove(erm[1],obj);
  off_kedge_remove(erm[2],obj);
  if(verbose) fprintf(stdout,"\tfree   e\n");
  off_kedge_free(erm[0]);
  off_kedge_free(erm[1]);
  off_kedge_free(erm[2]);
  return 1;
} /* int off_remesh_kedge_collapse_clean_up(kobj *obj,kvert **vrm,kface **frm,kedge **erm) */

int off_remesh_kedge_collapse2(kobj *obj,kedge *edge) {
  kvert *vrm[2];
  kface *frm[5];
  kedge *erm[5];

  if(!off_remesh_kedge_collapse(obj,edge,vrm,frm,erm)) {
    fprintf(stderr,"ERROR: off_remesh_kedge_collapse \n");
    return 0;
  }
  off_kvert_remove(vrm[0],obj);
  off_kvert_free(vrm[0]);
  off_kface_remove(frm[0],obj);
  off_kface_remove(frm[1],obj);
  off_kface_free(frm[0]);
  off_kface_free(frm[1]);
  off_kedge_remove(erm[0],obj);
  off_kedge_remove(erm[1],obj);
  off_kedge_remove(erm[2],obj);
  off_kedge_free(erm[0]);
  off_kedge_free(erm[1]);
  off_kedge_free(erm[2]);
  return 1;
}



/***************************************************************************
*
* off_remesh_kedge_collapse_correction
*
* -corrects the special case where there are 3 common neighbors by 2 points on the same edge
* -this happens in uneven distributed triangles
* -solution to this problem is to add another vertex on the edge that is not on
*  the adjface surface
*
* see also
*
*
*
***************************************************************************/
int off_remesh_kobj_collapse_correction(kobj *obj) {
  kedge *e;
  int nei,cnt=1;
  int verbose=0;
  while(cnt) {
    for(e=obj->edge,cnt=0; e!=NULL; e=e->next) {
      nei=off_kobj_test_common_neighbor(e);
      if(nei==2) continue;
      cnt++;
      if(!off_remesh_kedge_collapse_correction(obj,e)) {
        fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction \n");
        return 0;
      }
    }
    if(verbose) fprintf(stderr,"       off_remesh_kobj_collapse_correction # = %i\n",cnt);
  }
  return 1;
}


int off_remesh_kedge_collapse_correction(kobj *obj,kedge *edge) {
  kface *f[4];
  kvert *v[4];
  kedge *e=NULL;
  int n,m;
  int verbose=0;
  const char *fcname="off_remesh_kedge_collapse_correction";

  if(verbose) fprintf(stdout,"[%s]  off_remesh_kedge_collapse_correction edge = %i \n",fcname,edge->index);

  v[0]=v[1]=v[2]=v[3]=NULL;

  /* adjacent face */
  f[0]=edge->adjface[0];
  f[1]=edge->adjface[1];

  /* end points */
  v[0]=edge->endpts[0];
  v[1]=edge->endpts[1];

  /* make sure we didn't flip */
  for(m=0; m<3; m++)
    if(f[0]->vert[m]==v[1])
      if(f[0]->vert[(m+1)%3]!=v[0]) {
        v[0]=edge->endpts[1];
        v[1]=edge->endpts[0];
      }

  /* other points */
  for(m=0; m<2; m++) {
    for(n=0; n<3; n++) {
      if(f[m]->vert[n]==v[0]) continue;
      if(f[m]->vert[n]==v[1]) continue;
      v[m+2]=f[m]->vert[n];
      break;
    }
    if(v[m+2]==NULL) {
      fprintf(stderr,"[%s] ERROR: v[%i] could not be found \n",fcname,m+2);
      return 0;
    }
  }

  /* there should be additional point */
  for(n=0,e=NULL; n<edge->endpts[0]->nei&&e==NULL; n++) {
    for(m=0; m<edge->endpts[1]->nei&&e==NULL; m++) {
      if(edge->endpts[0]->neivert[n]==edge->endpts[1]->neivert[m]) {
        if(edge->endpts[1]->neivert[m]==v[2]) continue;
        if(edge->endpts[1]->neivert[m]==v[3]) continue;
        e = edge->endpts[1]->neiedge[m];
      }
    }
  }

  if(e!=NULL) {
    if(verbose) {
      fprintf(stdout,"[%s] %i \n",fcname,e->index);
    }
  } else if(e==NULL) {
    /* problem situation
    * -writes a smoothed output for me to visually inspect */
    fprintf(stderr,"[%s]  no edge!\n",fcname);
    fprintf(stderr,"  looking for v[3] %i \n",v[3]->index);
    for(n=0,e=NULL; n<v[2]->nei; n++) {
      fprintf(stderr,"     v[2]->neivert[%i] = %i\n",n,v[2]->neivert[n]->index);
    }
    off_kobj_add_color(obj);
    for(n=0; n<2; n++) {
      edge->adjface[n]->color[0]=1;
      edge->adjface[n]->color[2]=edge->adjface[n]->color[3]=0;
    }
    for(v[0]=obj->vert; v[0]!=NULL; v[0]=v[0]->next) off_remesh_kvert_relocate(v[0]);
    for(v[0]=obj->vert; v[0]!=NULL; v[0]=v[0]->next) off_remesh_kvert_relocate(v[0]);
    for(v[0]=obj->vert; v[0]!=NULL; v[0]=v[0]->next) off_remesh_kvert_relocate(v[0]);
    fprintf(stderr,"  writing tmp_prob.ply for error checking, look for red vertex \n");
    off_obj_renumber(obj);

    if(!off_kobj_write_offply("tmp_prob.ply",obj,0)) {
      fprintf(stderr,"[%s] ERROR: off_kobj_write_off \n",fcname);
    }
    return 0;
  }

  /* split this edge */
  if(off_remesh_kedge_split(obj,e)==NULL) {
    fprintf(stderr,"ERROR: off_remesh_kedge_split %i\n",e->index);
    return 0;
  }

  return 1;
}


/******************************************************
*
* off_remesh_kedge_collapse
*
* -actual function to do the edge collapsing
* -things that were removed are in the remove_* variables
* -they can be freed outside this function
* -function (off_remesh_kedge_collapse2) does freeing for you
* -function (off_remesh_kedge_collapse_clean_up) frees these variables after
*  this function
* -returns zero for error
* -return 1 for success
* -return -1 if no edge was actually collapsed
* -v[1] is to be removed
*
* see also
*
* int off_remesh_kedge_collapse2(kobj *obj,kedge *edge);
* int off_remesh_kedge_collapse_clean_up(kobj *obj,kvert **vrm,kface **frm,kedge **erm);
*
*
******************************************************/
int off_remesh_kedge_collapse(kobj *obj,kedge *edge,kvert **remove_vert,kface **remove_face,kedge **remove_edge) {
  kface *f[4];
  kvert *v[4];
  kedge *e[512];
  int n,m,nedge;
  niikpt midpt;
  niiksph sphpt;
  int verbose=0;

  if(edge==NULL) {
    fprintf(stderr,"ERROR: edge is a null pointer\n");
    return 0;
  }

  midpt=niikpt_avg(edge->endpts[0]->v,edge->endpts[1]->v);
  if(obj->spherecoo==1) {
    sphpt=niiksph_avg_lim(edge->endpts[0]->sph,edge->endpts[1]->sph,NIIK_PI);
  }

  if(verbose) {
    off_kobj_update_num(obj);
    fprintf(stdout,"\n  start off_remesh_kedge_collapse e=%i\n",edge->index);
    fprintf(stdout,"  vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }

  for(n=0; n<2; n++) f[n]=NULL;
  for(n=0; n<4; n++) v[n]=NULL;
  for(n=0; n<512; n++) e[n]=NULL;

  /* get the original configuration */
  if(off_remesh_kedge_local_info(edge,v,f,e)==0) {
    fprintf(stderr,"ERROR: off_remesh_kedge_local_info \n");
    return 0;
  }

  /* check if valence is not too small */
  if(v[0]->nei+v[1]->nei-4 <= 3) return -1;
  if(v[2]->nei-1 <= 3) return -1;
  if(v[3]->nei-1 <= 3) return -1;
  if(v[0]->nei+v[1]->nei>=13) return -1;

  /* check for common neighbros */
  if(off_kobj_test_common_neighbor(edge)>2) {
    if(!off_remesh_kedge_collapse_correction(obj,edge)) {
      fprintf(stderr,"ERROR: off_remesh_kedge_collapse_correction \n");
      return 0;
    }
    return -1;
  }

  /* -now we have f0-2, e0-4, v0-2
  * -we need to find f2-3 and v1's neighboring edges */
  if(verbose>2) fprintf(stdout,"\tfind f2-3\n");
  if     (e[1]->adjface[0]==f[0]) f[2]=e[1]->adjface[1];
  else if(e[1]->adjface[1]==f[0]) f[2]=e[1]->adjface[0];
  else {
    fprintf(stderr,"ERROR: could not find f[1]\n");
    return 0;
  }
  if     (e[2]->adjface[0]==f[1]) f[3]=e[2]->adjface[1];
  else if(e[2]->adjface[1]==f[1]) f[3]=e[2]->adjface[0];
  else {
    fprintf(stderr,"ERROR: could not find f[2]\n");
    return 0;
  }


  /* we then find the extra edges */
  for(n=0,m=4; n<v[1]->nei; n++) {
    if(v[1]->neiedge[n]==edge) continue;
    if(v[1]->neiedge[n]==e[1]) continue;
    if(v[1]->neiedge[n]==e[2]) continue;
    e[m++]=v[1]->neiedge[n];
    if(verbose>2) fprintf(stdout,"      edge [%i] = %i | %i\n",m-1,e[m-1]->index,n);
  }
  nedge=m;


  if(verbose) {
    fprintf(stdout,"    adjface %3i %3i\n",f[0]->index,f[1]->index);
    fprintf(stdout,"    face    %3i %3i %3i %3i\n",f[0]->index,f[1]->index,f[2]->index,f[3]->index);
    fprintf(stdout,"    vertex  %3i %3i %3i %3i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
    fprintf(stdout,"    edge [%3i] ",nedge);
    for(m=0; m<nedge; m++) fprintf(stdout,"%3i ",e[m]->index);
    fprintf(stdout,"\n");

    for(n=0; n<2; n++) {
      fprintf(stdout,"    vertex[%i]\n",n);
      off_display_vert_info(v[n],10);
    }
    for(n=0; n<4; n++) {
      fprintf(stdout,"    face[%i]\n",n);
      off_display_face_info(f[n],10);
    }
    fprintf(stdout,"    nedge = %i \n",nedge);
    for(n=0; n<nedge; n++) {
      fprintf(stdout,"    edge[%i]\n",n);
      off_display_edge_info(e[n],10);
    }

    for(n=0; n<2; n++) {
      fprintf(stdout,"  check v[%i]'s neivert's nei\n",n);
      for(m=0; m<v[n]->nei; m++) {
        fprintf(stdout,"    v[%i]->neivert[%i] (%i)->nei = %i\n",n,m,v[n]->neivert[m]->index,v[n]->neivert[m]->nei);
      }
    }
  } /* verbose */



  /* modify edges */
  if(verbose) {
    fprintf(stdout,"  modify edges e[0]=%i, e[3]=%i\n",e[0]->index,e[3]->index);
  }
  if     (e[0]->adjface[0]==f[0]) e[0]->adjface[0]=f[2];
  else if(e[0]->adjface[1]==f[0]) e[0]->adjface[1]=f[2];
  else {
    fprintf(stderr,"ERROR: could not find f[0] in e[0]->adjface %i %i\n",e[0]->adjface[0]->index,e[0]->adjface[1]->index);
    return 0;
  }
  if     (e[3]->adjface[0]==f[1]) e[3]->adjface[0]=f[3];
  else if(e[3]->adjface[1]==f[1]) e[3]->adjface[1]=f[3];
  else {
    fprintf(stderr,"ERROR: could not find f[1] in e[3]->adjface %i %i\n",e[3]->adjface[0]->index,e[3]->adjface[1]->index);
    return 0;
  }

  /* modify more edges */
  if(verbose) {
    fprintf(stdout,"  modify more edges\n");
  }
  for(n=4; n<nedge; n++) {
    if     (e[n]->endpts[0]==v[1]) e[n]->endpts[0]=v[0];
    else if(e[n]->endpts[1]==v[1]) e[n]->endpts[1]=v[0];
    else {
      fprintf(stderr,"ERROR: unknown edge in the list\n");
      fprintf(stderr,"  e[%i]\n",n);
      fprintf(stderr,"  e[%i] %i \n",n,e[n]->index);
      fprintf(stderr,"  e[%i] %i endpts %i\n",n,e[n]->index,e[n]->endpts[0]->index);
      fprintf(stderr,"  e[%i] %i endpts %i %i\n",n,e[n]->index,e[n]->endpts[0]->index,e[n]->endpts[1]->index);
      return 0;
    }
  }


  /* modify faces */
  if(verbose) {
    fprintf(stdout,"  modify faces\n");
  }
  for(n=0; n<v[1]->nei; n++) {
    for(m=0; m<3; m++) {
      if(v[1]->neiface[n]->vert[m]==v[1]) {
        if(verbose>2) {
          fprintf(stdout,"    v[1]->neiface[%i]->vert[%i] (%i) is changed to v=%i (neiface=%i)\n",n,m,
                  v[1]->neiface[n]->vert[m]->index,v[0]->index,v[1]->neiface[n]->index);
        }
        v[1]->neiface[n]->vert[m]=v[0];
        break;
      }
    }
  }
  for(m=0; m<3; m++) {
    if(f[2]->edge[m]==e[1]) {
      f[2]->edge[m]=e[0];
      break;
    }
  }
  for(m=0; m<3; m++) {
    if(f[3]->edge[m]==e[2]) {
      f[3]->edge[m]=e[3];
      break;
    }
  }


  /* update neigbors */
  if(verbose) {
    fprintf(stdout,"  update v0 neighbors\n");
  }
  if(!off_kobj_update_vertex_nei(v[0],e[0])) {
    fprintf(stderr,"ERROR at %d: off_kobj_update_vertex_nei(v=%i,e=%i)\n",__LINE__,v[0]->index,e[0]->index);
    return 0;
  }
  if(verbose) {
    fprintf(stdout,"  update v0->neivert neighbors\n");
  }
  for(n=0; n<v[0]->nei; n++) {
    if(!off_kobj_update_vertex_nei(v[0]->neivert[n],v[0]->neiedge[n])) {
      fprintf(stderr,"ERROR at %d: off_kobj_update_vertex_nei(v=%i,e=%i)\n",__LINE__,v[0]->neivert[n]->index,v[0]->neiedge[n]->index);
      return 0;
    }
  }

  /* remove things
  * -but don't free them yet */
  if(verbose) {
    fprintf(stdout,"  modify remove faces %i %i\n",f[0]->index,f[1]->index);
  }
  remove_face[0]=f[0];
  remove_face[1]=f[1];
  /*off_kface_remove(f[0],obj);
    off_kface_remove(f[1],obj);
    off_kface_free(f[0]);
    off_kface_free(f[1]);*/

  if(verbose) {
    fprintf(stdout,"  modify remove edges %i %i %i\n",edge->index,e[1]->index,e[2]->index);
  }
  remove_edge[0]=edge;
  remove_edge[1]=e[1];
  remove_edge[2]=e[2];
  /*off_kedge_remove(edge,obj);
    off_kedge_remove(e[1],obj);
    off_kedge_remove(e[2],obj);
    off_kedge_free(edge);
    off_kedge_free(e[1]);
    off_kedge_free(e[2]); */

  if(verbose) {
    fprintf(stdout,"  modify remove vertex %i\n",v[1]->index);
  }
  remove_vert[0]=v[1];
  /*off_kvert_remove(v[1],obj);
    off_kvert_free(v[1]);*/

  /* update the midpoint */
  v[0]->v=midpt;
  if(obj->spherecoo==1) {
    v[0]->sph=sphpt;
  }


  if(verbose>1) {
    fprintf(stdout,"  remaining edges\n");
    fprintf(stdout,"e[0] = ");
    off_display_edge_info(e[0],10);
    fprintf(stdout,"e[0] = adjface[0]: ");
    off_display_face_info(e[0]->adjface[0],10);
    fprintf(stdout,"e[0] = adjface[1]: ");
    off_display_face_info(e[0]->adjface[1],10);
    fprintf(stdout,"e[3] = ");
    off_display_edge_info(e[3],10);
    fprintf(stdout,"e[3] = adjface[0]: ");
    off_display_face_info(e[3]->adjface[0],10);
    fprintf(stdout,"e[3] = adjface[1]: ");
    off_display_face_info(e[3]->adjface[1],10);
    fprintf(stdout,"  remaining faces\n");
    fprintf(stdout,"f[2] = ");
    off_display_face_info(f[2],10);
    fprintf(stdout,"f[3] = ");
    off_display_face_info(f[3],10);

    fprintf(stdout,"  collapsed vertex\n");
    off_display_vert_info(v[0],10);

    fprintf(stdout,"  collapsed vertex's neiface\n");
    for(n=0; n<v[0]->nei; n++) {
      off_display_face_info(v[0]->neiface[n],10);
    }
    fprintf(stdout,"  collapsed vertex's neiedge\n");
    for(n=0; n<v[0]->nei; n++) {
      off_display_edge_info(v[0]->neiedge[n],10);
    }
    fprintf(stdout,"  collapsed vertex's neivert\n");
    for(n=0; n<v[0]->nei; n++) {
      off_display_vert_info(v[0]->neivert[n],10);
    }
  }


  if(verbose) {
    off_kobj_update_num(obj);
    fprintf(stdout,"  vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }

  /*  if(f[2]==f[3]) {fprintf(stdout,"special case\n");exit(0);}*/

  return 1;

} /* int off_remesh_kedge_collapse(kobj *obj,kedge *edge,kvert **remove_vert,kface **remove_face,kedge **remove_edge) */



/******************************************************
*
* kvert *off_remesh_kedge_split(kobj *obj,kedge *edge);
*
* splits edge (second version)
*
* -part of remeshing
* -returns the new vertex
*
* -adds 1 new vertex, 3 new edges, and 3 new faces
* -recycles edge and faces so there's no removel
*
*******************************************************/
kvert *off_remesh_kedge_split(kobj *obj,kedge *edge) {
  kface *f[2],*nf[2],*ff;
  kvert *v[4],*nv,*vv;
  kedge *e[4],*ne[3],*ee;
  int n;
  const int verbose=0;

  if(edge==NULL) {
    fprintf(stderr,"ERROR: edge is a null pointer\n");
    return NULL;
  }

  if(verbose) {
    off_kobj_update_num(obj);
    fprintf(stdout,"\nstart off_remesh_kedge_split e=%i\n",edge->index);
    fprintf(stdout,"  vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }

  for(n=0; n<2; n++) f[n]=NULL;
  for(n=0; n<4; n++) v[n]=NULL;
  for(n=0; n<4; n++) e[n]=NULL;

  /* get the original configuration */
  if(off_remesh_kedge_local_info(edge,v,f,e)==0) {
    fprintf(stderr,"ERROR: off_remesh_kedge_local_info \n");
    return 0;
  }

  if(verbose) {
    fprintf(stdout,"  adjface %7i %7i\n",f[0]->index,f[1]->index);
    fprintf(stdout,"  vertex  %7i %7i %7i %7i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
    fprintf(stdout,"  edge    %7i %7i %7i %7i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
  }

  /* we'll not remove anything for this version
  * and add 4 new faces, 4 new edges, and 1 vertex (returned) */

  /* new vertex */
  if(verbose) fprintf(stdout,"  make new vertex\n");
  nv = off_vert_init();
  nv->v = niikpt_avg(v[0]->v,v[1]->v); /* point is the mid point of the given edge  */
  nv->index=-1;
  if(obj->spherecoo==1) {
    nv->sph = niiksph_avg_lim(v[0]->sph, v[1]->sph, NIIK_PI);
  } /* spherecoordinates */

  /* new faces */
  if(verbose) fprintf(stdout,"  make new faces\n");
  for(n=0; n<2; n++) {
    nf[n]=off_face_init();
    nf[n]->index=-(n+1);
  }

  /* new edges */
  if(verbose) fprintf(stdout,"  make new edges\n");
  for(n=0; n<3; n++) {
    ne[n]=off_edge_init();
    ne[n]->index=-(n+1);
  }

  /* update edge info */
  n=0;
  ne[n]->adjface[0]=nf[0];
  ne[n]->adjface[1]=nf[1];
  n=1;
  ne[n]->adjface[0]= f[0];
  ne[n]->adjface[1]=nf[0];
  n=2;
  ne[n]->adjface[0]= f[1];
  ne[n]->adjface[1]=nf[1];
  n=0;
  ne[n]->endpts[0]=nv;
  ne[n]->endpts[1]=v[1];
  n=1;
  ne[n]->endpts[0]=nv;
  ne[n]->endpts[1]=v[2];
  n=2;
  ne[n]->endpts[0]=nv;
  ne[n]->endpts[1]=v[3];
  edge->endpts[0]=v[0];
  edge->endpts[1]=nv;

  if     (e[1]->adjface[0]==f[0]) e[1]->adjface[0]=nf[0];
  else if(e[1]->adjface[1]==f[0]) e[1]->adjface[1]=nf[0];
  if     (e[2]->adjface[0]==f[1]) e[2]->adjface[0]=nf[1];
  else if(e[2]->adjface[1]==f[1]) e[2]->adjface[1]=nf[1];

  /* update face info */
  if(verbose) fprintf(stdout,"  update face\n");
  for(n=0; n<3; n++)
    if(f[0]->vert[n]==v[1]) {
      f[0]->vert[n]=nv;
      break;
    }
  for(n=0; n<3; n++)
    if(f[0]->edge[n]==e[1]) {
      f[0]->edge[n]=ne[1];
      break;
    }
  if(verbose) {
    fprintf(stderr,"\tf[0]-> v=%i %i %i | edge %i(pt=%i,%i) %i(pt=%i,%i) %i(pt=%i,%i)\n",
            f[0]->vert[0]->index,f[0]->vert[1]->index,f[0]->vert[2]->index,
            f[0]->edge[0]->index,f[0]->edge[0]->endpts[0]->index,f[0]->edge[0]->endpts[1]->index,
            f[0]->edge[1]->index,f[0]->edge[1]->endpts[0]->index,f[0]->edge[1]->endpts[1]->index,
            f[0]->edge[2]->index,f[0]->edge[2]->endpts[0]->index,f[0]->edge[2]->endpts[1]->index);
  }

  for(n=0; n<3; n++)
    if(f[1]->vert[n]==v[1]) {
      f[1]->vert[n]=nv;
      break;
    }
  for(n=0; n<3; n++)
    if(f[1]->edge[n]==e[2]) {
      f[1]->edge[n]=ne[2];
      break;
    }
  if(verbose) {
    fprintf(stderr,"\tf[1]-> v=%i %i %i | edge %i(pt=%i,%i) %i(pt=%i,%i) %i(pt=%i,%i)\n",
            f[1]->vert[0]->index,f[1]->vert[1]->index,f[1]->vert[2]->index,
            f[1]->edge[0]->index,f[1]->edge[0]->endpts[0]->index,f[1]->edge[0]->endpts[1]->index,
            f[1]->edge[1]->index,f[1]->edge[1]->endpts[0]->index,f[1]->edge[1]->endpts[1]->index,
            f[1]->edge[2]->index,f[1]->edge[2]->endpts[0]->index,f[1]->edge[2]->endpts[1]->index);
  }

  nf[0]->vert[0]=nv;
  nf[0]->vert[1]=v[2];
  nf[0]->vert[2]=v[1];
  nf[0]->edge[0]=ne[1];
  nf[0]->edge[1]= e[1];
  nf[0]->edge[2]=ne[0];

  if(verbose) {
    fprintf(stderr,"\tnf[0]-> v=%i %i %i | edge %i(pt=%i,%i) %i(pt=%i,%i) %i(pt=%i,%i)\n",
            nf[0]->vert[0]->index,nf[0]->vert[1]->index,nf[0]->vert[2]->index,
            nf[0]->edge[0]->index,nf[0]->edge[0]->endpts[0]->index,nf[0]->edge[0]->endpts[1]->index,
            nf[0]->edge[1]->index,nf[0]->edge[1]->endpts[0]->index,nf[0]->edge[1]->endpts[1]->index,
            nf[0]->edge[2]->index,nf[0]->edge[2]->endpts[0]->index,nf[0]->edge[2]->endpts[1]->index);
  }

  nf[1]->vert[0]=nv;
  nf[1]->vert[1]=v[1];
  nf[1]->vert[2]=v[3];
  nf[1]->edge[0]=ne[0];
  nf[1]->edge[1]= e[2];
  nf[1]->edge[2]=ne[2];

  if(verbose) {
    fprintf(stderr,"\tnf[1]-> v=%i %i %i | edge %i(pt=%i,%i) %i(pt=%i,%i) %i(pt=%i,%i)\n",
            nf[1]->vert[0]->index,nf[1]->vert[1]->index,nf[1]->vert[2]->index,
            nf[1]->edge[0]->index,nf[1]->edge[0]->endpts[0]->index,nf[1]->edge[0]->endpts[1]->index,
            nf[1]->edge[1]->index,nf[1]->edge[1]->endpts[0]->index,nf[1]->edge[1]->endpts[1]->index,
            nf[1]->edge[2]->index,nf[1]->edge[2]->endpts[0]->index,nf[1]->edge[2]->endpts[1]->index);
    fprintf(stderr,"\tedge    (%i) -> adjface=%i %i | endpts %i %i \n",
            edge->index,
            edge->adjface[0]->index,edge->adjface[1]->index,
            edge->endpts[0]->index,edge->endpts[1]->index);
    for(n=0; n<4; n++)
      fprintf(stderr,"\te [%i] (%i) -> adjface=%i %i | endpts %i %i \n",n,
              e[0]->index,
              e[0]->adjface[0]->index,e[0]->adjface[1]->index,
              e[0]->endpts[0]->index,e[0]->endpts[1]->index);
    for(n=0; n<3; n++)
      fprintf(stderr,"\tne[%i] (%i) -> adjface=%i %i | endpts %i %i \n",n,
              ne[0]->index,
              ne[0]->adjface[0]->index,ne[0]->adjface[1]->index,
              ne[0]->endpts[0]->index,ne[0]->endpts[1]->index);
  }


  /* update neighbors */
  if(verbose) fprintf(stdout,"  update neighbors\n");
  if(!off_kobj_update_vertex_nei(nv,edge)) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(nv,edge)\n");
    return NULL;
  }

  if(verbose>1) fprintf(stdout,"    update neighbors\n");
  if(!off_kobj_update_vertex_nei(v[0],edge)) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[0],edge)\n");
    return NULL;
  }

  if(verbose>1) fprintf(stdout,"    update neighbors\n");
  if(!off_kobj_update_vertex_nei(v[1],ne[0])) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[1],ne[0])\n");
    return NULL;
  }

  if(verbose>1) fprintf(stdout,"    update neighbors\n");
  if(!off_kobj_update_vertex_nei(v[2],ne[1])) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[2],ne[1])\n");
    return NULL;
  }

  if(verbose>1) fprintf(stdout,"    update neighbors\n");
  if(!off_kobj_update_vertex_nei(v[3],ne[2])) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[3],ne[2])\n");
    return NULL;
  }

  if(verbose>1) fprintf(stdout,"    update neighbors done \n");


  /* insert new things */
  if(verbose) fprintf(stdout,"  insert vertex\n");
  nv->next = v[0]->next;
  if(nv->next!=NULL) nv->next->prev = nv;
  nv->prev = v[0];
  v[0]->next = nv;

  if(verbose) fprintf(stdout,"  connect new edges\n");
  ne[0]->next = ne[1];
  ne[1]->next = ne[2];
  ne[2]->prev = ne[1];
  ne[1]->prev = ne[0];

  if(verbose) fprintf(stdout,"  insert new edges\n");
  ne[2]->next = e[3]->next;
  if(e[3]->next!=NULL) e[3]->next->prev=ne[2];
  e[3]->next = ne[0];
  ne[0]->prev = e[3];

  if(verbose) fprintf(stdout,"  connect new faces\n");
  nf[0]->next = nf[1];
  nf[1]->prev = nf[0];

  if(verbose) fprintf(stdout,"  insert new faces after f[1]=%i\n",f[1]->index);
  nf[1]->next = f[1]->next;
  if(nf[1]->next!=NULL) nf[1]->next->prev = nf[1];
  f[1]->next = nf[0];
  nf[0]->prev = f[1];


  if(verbose) {
    if(verbose>1) {
      fprintf(stdout,"  forward \n");
      for(vv=obj->vert; vv!=NULL; vv=vv->next) {
        fprintf(stdout,"v = %7i\n",vv->index);
        if(vv->next==NULL) break;
      }
      for(ff=obj->face; ff!=NULL; ff=ff->next) {
        fprintf(stdout,"f = %7i\n",ff->index);
        if(ff->next==NULL) break;
      }
      for(ee=obj->edge; ee!=NULL; ee=ee->next) {
        fprintf(stdout,"e = %7i\n",ee->index);
        if(ee->next==NULL) break;
      }

      fprintf(stdout,"  backward \n");
      for(vv=vv; vv!=NULL; vv=vv->prev) {
        fprintf(stdout,"v = %7i\n",vv->index);
      }
      for(ff=ff; ff!=NULL; ff=ff->prev) {
        fprintf(stdout,"f = %7i\n",ff->index);
      }
      for(ee=ee; ee!=NULL; ee=ee->prev) {
        fprintf(stdout,"e = %7i\n",ee->index);
      }
    }

    off_kobj_update_num(obj);
    fprintf(stdout,"  vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);

    if(!off_kobj_test(obj)) {
      fprintf(stderr,"ERROR: test failed\n");
      exit(0);
    }
  } /* verbose */

  return nv;
} /* kvert *off_remesh_kedge_split(kobj *obj,kedge *edge) */



/******************************************************
*
* off_remesh_kedge_flip: flips and edge
*
* -no addition or removal
*
******************************************************/
int off_remesh_kedge_flip(kobj *obj,kedge *edge) {
  kvert *v[4];
  kedge *e[4];
  kface *f[2];
  int n;
  const int verbose=0;

  if(edge==NULL) {
    fprintf(stderr,"ERROR: edge is a null pointer\n");
    return 0;
  }

  if(verbose) {
    off_kobj_update_num(obj);
    fprintf(stdout,"\nstart off_remesh_kedge_flip  e=%i\n",edge->index);
    fprintf(stdout,"  vfe = %i %i %i\n",obj->nvert,obj->nface,obj->nedge);
  }

  for(n=0; n<2; n++) f[n]=NULL;
  for(n=0; n<4; n++) v[n]=NULL;
  for(n=0; n<4; n++) e[n]=NULL;

  /* get the original configuration */
  if(off_remesh_kedge_local_info(edge,v,f,e)==0) {
    fprintf(stderr,"ERROR: off_remesh_kedge_local_info \n");
    return 0;
  }

  /*      fprintf(stdout,"  local info\n");
      fprintf(stdout,"  edge = %i\n",edge->index);
      fprintf(stdout,"  vert = %i %i %i %i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
      fprintf(stdout,"  edge = %i %i %i %i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
      fprintf(stdout,"  face = %i %i\n",f[0]->index,f[1]->index);
      off_display_face_info(f[0],10);
      off_display_face_info(f[1],10);*/

  for(n=0; n<4; n++)
    if(e[n]==NULL) {
      fprintf(stderr,"ERROR: edge was not assigned\n");
      fprintf(stdout,"       off_remesh_kedge_flip e=%i\n",edge->index);
      return 0;
    }
  if(verbose) {
    fprintf(stdout,"  adjface %7i %7i\n",f[0]->index,f[1]->index);
    fprintf(stdout,"  vertex  %7i %7i %7i %7i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
    fprintf(stdout,"  edge    %7i %7i %7i %7i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
  }

  /*
  * edge to do
  * -edge's end points are different
  * -e[0] and e[2] have different adjface
  *
  * face to do
  * -vertex and edge are different
  *
  * vertex to do
  * -update neighbor list
  */

  edge->endpts[0]=v[2];
  edge->endpts[1]=v[3];
  if     (e[0]->adjface[0]==f[0]) e[0]->adjface[0]=f[1];
  else if(e[0]->adjface[1]==f[0]) e[0]->adjface[1]=f[1];
  else {
    fprintf(stderr,"ERROR: could not find f[0] in e[0]->adjface[?]\n");
    return 0;
  }
  if     (e[2]->adjface[0]==f[1]) e[2]->adjface[0]=f[0];
  else if(e[2]->adjface[1]==f[1]) e[2]->adjface[1]=f[0];
  else {
    fprintf(stderr,"ERROR: could not find f[1] in e[2]->adjface[?]\n");
    return 0;
  }


  f[0]->vert[0]=v[1];
  f[0]->vert[1]=v[3];
  f[0]->vert[2]=v[2];
  f[0]->edge[0]=e[2];
  f[0]->edge[1]=edge;
  f[0]->edge[2]=e[1];
  f[1]->vert[0]=v[2];
  f[1]->vert[1]=v[3];
  f[1]->vert[2]=v[0];
  f[1]->edge[0]=edge;
  f[1]->edge[1]=e[3];
  f[1]->edge[2]=e[0];

  if(verbose>1) {
    fprintf(stdout,"  update v=%i e=%i\n",v[0]->index,e[0]->index);
  }
  if(!off_kobj_update_vertex_nei(v[0],e[0])) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[0],e[0])\n");
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"  update v=%i e=%i\n",v[1]->index,e[1]->index);
  }
  if(!off_kobj_update_vertex_nei(v[1],e[1])) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[1],e[1])\n");
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"  update v=%i e=%i\n",v[2]->index,edge->index);
  }
  if(!off_kobj_update_vertex_nei(v[2],edge)) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[2],edge)\n");
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"  update v=%i e=%i\n",v[3]->index,edge->index);
  }
  if(!off_kobj_update_vertex_nei(v[3],edge)) {
    fprintf(stderr,"ERROR: off_kobj_update_vertex_nei(v[3],edge)\n");
    return 0;
  }

  if(verbose) {
    fprintf(stdout,"  adjface %7i %7i\n",f[0]->index,f[1]->index);
    fprintf(stdout,"  vertex  %7i %7i %7i %7i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
    fprintf(stdout,"  edge    %7i %7i %7i %7i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
    if(!off_kobj_test(obj)) {
      fprintf(stderr,"ERROR: test failed\n");
      exit(0);
    }
  }

  /*  for(n=0;n<2;n++){
    off_update_kface_normal(f[n]);
    if(niikpt_dot(niikpt_unit(f[n]->vert[0]->v),f[n]->normal)<-0.5) {
      fprintf(stdout,"%9i: %12.6f | %5.1f %5.1f %5.1f | %5.1f %5.1f %5.1f | %1i\n",
              f[n]->index,niikpt_dot(niikpt_unit(f[n]->vert[0]->v),f[n]->normal),
              f[n]->normal.x,f[n]->normal.y,f[n]->normal.z,
              f[n]->vert[0]->v.x,f[n]->vert[0]->v.y,f[n]->vert[0]->v.z,n);

      fprintf(stdout,"  local info\n");
      fprintf(stdout,"  edge = %i\n",edge->index);
      fprintf(stdout,"  vert = %i %i %i %i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
      fprintf(stdout,"  edge = %i %i %i %i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
      fprintf(stdout,"  face = %i %i\n",f[0]->index,f[1]->index);
      off_display_face_info(f[0],10);
      off_display_face_info(f[1],10);
      exit(0); }
    fprintf(stdout,"  normal direction is OK! \n");
    } */

  return 1;
} /* int off_remesh_kedge_flip(kobj *obj,kedge *edge) */





/******************************************************
*
* off_remesh_kobj_valence_minimization
*
* -no addition or removal
* -minimize the valence for each vertex
*
*
* functions
*   off_remesh_kobj_valence_minimization
*   off_remesh_kobj_valence_minimization_kedge
*   off_remesh_kobj_count_global_valence6
*   off_remesh_kobj_count_global_valence57
*   off_remesh_kobj_count_global_valence_variance
*
*
*
******************************************************/
int off_remesh_kobj_count_global_valence6(kobj *obj)
/* counts the number of vertex that have exactly 6 neighbors */
{
  kvert *v;
  int n;
  for(v=obj->vert,n=0; v!=NULL; v=v->next) {
    n+=(v->nei==6);
  }
  return n;
}

int off_remesh_kobj_count_global_valence57(kobj *obj)
/* counts the number of vertices that have 5,6, or 7 neighbors */
{
  kvert *v;
  int n;
  for(v=obj->vert,n=0; v!=NULL; v=v->next) {
    n+=(v->nei>=5 && v->nei<=7);
  }
  return n;
}


int off_remesh_kobj_count_global_valence_variance(kobj *obj) {
  kvert *v;
  int n;
  for(v=obj->vert,n=0; v!=NULL; v=v->next) {
    n+=(v->nei-6)*(v->nei-6);
  }
  return n;
}


int off_remesh_kobj_valence_minimization_kedge(kobj *obj,kedge *edge) {
  kface *f[2];
  kvert *v[4];
  int m,n,neisum,neisum2;
  const int verbose=0;

  for(m=0; m<2; m++) {
    f[m]=edge->adjface[m];
    v[m]=edge->endpts[m];
  }

  /* if valence is too small then don't decrease any more */
  if(v[0]->nei<=4) return 0;
  if(v[1]->nei<=4) return 0;

  /* find 4 vertices around these faces */
  for(m=0; m<2; m++) {
    for(n=0; n<3; n++) {
      if(f[m]->vert[n]==v[0]) continue;
      if(f[m]->vert[n]==v[1]) continue;
      v[m+2]=f[m]->vert[n];
      break;
    }
  }
  if(verbose>1) fprintf(stdout,"\t  %i %i %i %i \n",v[0]->nei,v[1]->nei,v[2]->nei,v[3]->nei);

  /* calculate the sum of valence */
  for(n=neisum=0; n<4; n++) {
    m = v[n]->nei - 6;
    neisum+=m*m*m*m;
  }
  if(m==0) return 0; /* can't get any better; m<=1 didn't work */

  /* calculate the sum of valence after flipping */
  for(n=neisum2=0; n<4; n++) {
    if(n<2) m = v[n]->nei - 7;
    else    m = v[n]->nei - 5;
    neisum2+=m*m*m*m;
  } /* neisum2+=m*m; }*/
  if(neisum<=neisum2) return 0;

  if(verbose) {
    fprintf(stdout,"\tedge %i flip (%i -> %i)! \n",edge->index,neisum,neisum2);
    fprintf(stdout,"\t  %i %i %i %i \n",v[0]->nei,v[1]->nei,v[2]->nei,v[3]->nei);
  }

  /* flip if valence is smaller */
  if(!off_remesh_kedge_flip(obj,edge)) {
    fprintf(stderr,"ERROR: off_remesh_kedge_flip %i\n",edge->index);
    return -1;
  }

  if(verbose) {
    fprintf(stdout,"\t  %i %i %i %i \n",v[0]->nei,v[1]->nei,v[2]->nei,v[3]->nei);
  }
  return 1;
}


int off_remesh_kobj_valence_minimization(kobj *obj,int wmark) {
  kedge *e;
  int n;
  const int verbose=0;
  if(verbose) fprintf(stdout,"  off_remesh_kobj_valence_minimization\n");
  /* visit all edges */
  for(e=obj->edge; e!=NULL; e=e->next) {
    if(wmark)
      if(e->endpts[0]->v.w>0 || e->endpts[1]->v.w>0) {
        if(e->next==NULL) break;
        continue;
      }
    if((n=off_remesh_kobj_valence_minimization_kedge(obj,e))<0) {
      fprintf(stderr,"ERROR: off_remesh_kobj_valence_minimization_kedge %i \n",e->index);
      return 0;
    }
    if(n) {
      e->endpts[0]->v.w=e->endpts[1]->v.w=0.0;
    }
    if(e->next==NULL) break;
  }

  for(e=e; e->prev!=NULL; e=e->prev) {
    if(wmark)
      if(e->endpts[0]->v.w>0 && e->endpts[1]->v.w>0)
        continue;

    if(off_remesh_kobj_valence_minimization_kedge(obj,e)<0) {
      fprintf(stderr,"ERROR: off_remesh_kobj_valence_minimization_kedge %i \n",e->index);
      return 0;
    }
  }
  return 1;
}




/******************************************************
*
* off_remesh_kvert_relocate
*
* -no addition or removal
* -tangential relaxation
* - will move the edge by 0.95 towards the weighted average of neigbours, 
* - in the direction of the norm
* 
*
******************************************************/
int off_remesh_kvert_relocate(kvert *vert) {
  return off_remesh_kvert_relocate_lambda(vert,0.95);
} /* off_remesh_kvert_relocate */


int off_remesh_kvert_relocate_lambda(kvert *vert,double lambda) {
  niikpt g,p;
  double a;
  int n;

  /* update normal just once */
  g=niikpt_zero();

  for(n=0; n<vert->nei; n++) {
    a = niikpt_area2(vert->neiface[n]->vert[0]->v,vert->neiface[n]->vert[1]->v,vert->neiface[n]->vert[2]->v);
    p = niikpt_avg3 (vert->neiface[n]->vert[0]->v,vert->neiface[n]->vert[1]->v,vert->neiface[n]->vert[2]->v);
    g.x += a * p.x;
    g.y += a * p.y;
    g.z += a * p.z;
    g.w += a;
  }
  g.x /= g.w;
  g.y /= g.w;
  g.z /= g.w;

  for(n=0; n<vert->nei; n++)
    off_update_kface_normal(vert->neiface[n]);

  off_update_kvert_normal(vert);

  g.w = vert->v.w; /* make sure that v.w flag stays the same */
  vert->v = niikpt_move_normal(g, vert->normal, lambda*niikpt_dot(vert->normal, niikpt_sub(vert->v,g)));

  return 1;
} /* off_remesh_kvert_relocate */





/*********************************************************************************
*
* int off_kobj_test(kobj *obj);
*
* -test the object's validity
* -it's not a complete test
* -check the following
* --Is the vertex its neighbor's neighbor?
* --Valence is too small? (<=3?)
* --End points should have two common neighbors
* --Is the backward link broken? for vert,face,edge?
*
*********************************************************************************/
int off_kobj_test(kobj *obj) {
  const char *fcname=__func__;
  kedge *e;
  kvert *v;
  kface *f;
  int n,m,hit;

  /* is is connected to other vertices? */
  fprintf(stdout,"[%s] Each vertex is connected to at least one other vertex\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) {
    if(v->nei==0) {
      fprintf(stderr,"\n\nERROR: lonely vertex\n");
      fprintf(stderr,"ERROR: v[%i]\n",v->index);
      off_display_vert_info(v,10);
      for(n=0; n<v->nei; n++) {
        off_display_vert_info(v->neivert[n],10);
      }
      for(f=obj->face; f!=NULL; f=f->next) {
        if(f->vert[0]==v) {
          off_display_face_info(f,10);
        }
        if(f->vert[1]==v) {
          off_display_face_info(f,10);
        }
        if(f->vert[2]==v) {
          off_display_face_info(f,10);
        }
      }
      for(e=obj->edge; e!=NULL; e=e->next) {
        if(e->endpts[0]==v) {
          off_display_edge_info(e,10);
        }
        if(e->endpts[1]==v) {
          off_display_edge_info(e,10);
        }
      }
      return 0;
    }
  }

  /* Am I my neighbor's neighbor? */
  fprintf(stdout,"[%s] Check neighbor's neighbors include the vertex itself\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) {
    hit=0;
    for(n=0; n<v->nei; n++) {
      for(m=0; m<v->neivert[n]->nei; m++) {
        if(v->neivert[n]->neivert[m]==v)
          hit=1;
      }
    }
    if(!hit) {
      fprintf(stderr,"\n\nERROR: neighbor's neighbor test\n");
      fprintf(stderr,"ERROR: v[%i]\n",v->index);
      off_display_vert_info(v,10);
      for(n=0; n<v->nei; n++) {
        off_display_vert_info(v->neivert[n],10);
      }
      return 0;
    }
  }

  /* Correspondence between neiedge and neivert */
  fprintf(stdout,"[%s] Check the vertex's conneting edge has itself as the endpoint\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) {
    hit=0;
    for(n=0; n<v->nei; n++) {
      if     (v->neiedge[n]->endpts[0]==v) break;
      else if(v->neiedge[n]->endpts[1]==v) break;
      else {
        fprintf(stderr,"ERROR: v[%i]->neivert[%i] edge and vert mismatch\n",
                v->index,n);
        return 0;
      }
    }
  }

  /* Valence is too small */
  fprintf(stdout,"[%s] Check valence is too small <=3\n",fcname);
  for(v=obj->vert; v!=NULL; v=v->next) {
    if(v->nei<=3) {
      fprintf(stderr,"ERROR: valence is too small\n");
      off_display_vert_info(v,10);
      return 0;
    }
  }

  /* Endpoint should have only 2 common neighbors */
  fprintf(stdout,"[%s] Checking the neighboring\n",fcname);


  /* Endpoint should have only 2 common neighbors */
  fprintf(stdout,"[%s] 2 common neighbors by 2 endpoints on an edge\n",fcname);
  for(e=obj->edge; e!=NULL; e=e->next) {
    if((hit=off_kobj_test_common_neighbor(e))!=2) {
      fprintf(stdout,"\nERROR: edge[%i]'s endpts(%i,%i) have %i common neighbors\n",
              e->index,e->endpts[0]->index,e->endpts[1]->index,
              hit);
      for(n=0; n<e->endpts[0]->nei; n++) {
        for(m=0; m<e->endpts[1]->nei; m++) {
          if(e->endpts[0]->neivert[n]==e->endpts[1]->neivert[m]) {
            off_display_vert_info(e->endpts[0]->neivert[n],10);
          }
        }
      }
      fprintf(stdout,"\n");
      for(n=0; n<2; n++) {
        fprintf(stdout,"e[%i]->endpts[%i]\n",e->index,n);
        off_display_vert_info(e-> endpts[n],20);
        off_display_face_info(e->adjface[n],20);
      }
      off_kobj_gray_kface(obj);
      e->adjface[0]->color[0]=e->adjface[1]->color[0]=1.0;
      for(n=1; n<=4; n++)
        e->adjface[0]->color[n]=e->adjface[1]->color[n]=0.0;
      fprintf(stdout,"\n");
      /*fprintf(stdout,"vertex around this edge:\n");
      for(v=obj->vert;v!=NULL;v=v->next){
        if     (niikpt_distance(v->v,e->endpts[0]->v)<1.2){
          off_display_vert_info(v,20); }
        else if(niikpt_distance(v->v,e->endpts[1]->v)<1.2){
          off_display_vert_info(v,20); }
          }*/
      fprintf(stderr,"  writing tmp_prob.off.gz for error checking \n");
      if(!off_kobj_write_offply("tmp_prob.off.gz",obj,0)) {
        fprintf(stderr,"ERROR: off_kobj_write_off \n");
        exit(0);
      }
      return 0;
    }
  }

  fprintf(stdout,"[%s] testing links (face)\n",fcname);
  if(!off_kface_test_link(obj->face)) {
    fprintf(stderr,"ERROR: face link broken\n");
    return 0;
  }
  fprintf(stdout,"[%s] testing links (vert)\n",fcname);
  if(!off_kvert_test_link(obj->vert)) {
    fprintf(stderr,"ERROR: vert link is broken\n");
    return 0;
  }
  fprintf(stdout,"[%s] testing links (edge)\n",fcname);
  if(!off_kedge_test_link(obj->edge)) {
    fprintf(stderr,"ERROR: edge link is broken\n");
    return 0;
  }
  return 1;
}



int off_kobj_test_common_neighbor(kedge *edge) {
  int hit=0,m,n;
  for(n=0; n<edge->endpts[0]->nei; n++) {
    for(m=0; m<edge->endpts[1]->nei; m++) {
      if(edge->endpts[0]->neivert[n]==edge->endpts[1]->neivert[m]) {
        hit++;
      }
    }
  }
  return hit;
}



/***************************************************************************
*
* int off_remesh_kedge_local_info(kedge *edge,kvert **v,kface **f,kedge **e);
*
* -get the local information
* -edge is the input
* -v[0,1] are the edge's endpts
* -f[0,1] are the edge's adjfaces
* -v[2] is the common neighboring vertex on f[0] side
* -v[3] is the common neighboring vertex on f[1] side
* -e[0] is between v[0] and v[2]
* -e[1] is between v[1] and v[2]
* -e[2] is between v[1] and v[3]
* -e[3] is between v[0] and v[3]
*
* -these aren't done b/c only collapse need them
* --f[2] is e[1]'s adjface (but not f[0])
* --f[3] is e[2]'s adjface (but not f[1])
*
* -checks for the number of common neighboring vertex
* --if not two, returns 0
* -checks other things, and returns 0 for unsuccessful cases
*
*
****************************************************************************/
int off_remesh_kedge_local_info(kedge *edge,kvert **v,kface **f,kedge **e) {
  int m,n;

  /* adjacent face */
  f[0]=edge->adjface[0];
  f[1]=edge->adjface[1];

  /* end points */
  v[0]=edge->endpts[0];
  v[1]=edge->endpts[1];

  /* make sure we didn't flip */
  for(m=0; m<3; m++)
    if(f[0]->vert[m]==v[1])
      if(f[0]->vert[(m+1)%3]!=v[0]) {
        v[0]=edge->endpts[1];
        v[1]=edge->endpts[0];
      }


  /* 2 common neighboring vertices
  * -check if there are actually two
  *
  * -this should be checked elsewhere
  *
  *
  for(n=0,nei=0;n<edge->endpts[0]->nei;n++){
    for(m=0;m<edge->endpts[1]->nei;m++){
      if(edge->endpts[0]->neivert[n]==edge->endpts[1]->neivert[m]){
        nei++; }
    } }
  if(nei!=2) return 0;
  */

  for(m=0; m<2; m++) {
    for(n=0; n<3; n++) {
      if(f[m]->vert[n]==v[0]) continue;
      if(f[m]->vert[n]==v[1]) continue;
      v[m+2]=f[m]->vert[n];
      break;
    }
  }

  for(n=0; n<3; n++)
    if(f[0]->edge[n]==edge) continue;
    else if(f[0]->edge[n]->endpts[0]==v[0]) {
      e[0]=f[0]->edge[n];
      break;
    } else if(f[0]->edge[n]->endpts[1]==v[0]) {
      e[0]=f[0]->edge[n];
      break;
    }
  if(n==3) {
    fprintf(stderr,"ERROR: could not find e[0] \n");
    return 0;
  }

  for(n=0; n<3; n++)
    if(f[0]->edge[n]==edge) continue;
    else if(f[0]->edge[n]->endpts[0]==v[1]) {
      e[1]=f[0]->edge[n];
      break;
    } else if(f[0]->edge[n]->endpts[1]==v[1]) {
      e[1]=f[0]->edge[n];
      break;
    }
  if(n==3) {
    fprintf(stderr,"ERROR: could not find e[1] \n");
    return 0;
  }

  for(n=0; n<3; n++)
    if(f[1]->edge[n]==edge) continue;
    else if(f[1]->edge[n]->endpts[0]==v[0]) {
      e[3]=f[1]->edge[n];
      break;
    } else if(f[1]->edge[n]->endpts[1]==v[0]) {
      e[3]=f[1]->edge[n];
      break;
    }
  if(n==3) {
    fprintf(stderr,"ERROR: could not find e[3] \n");
    return 0;
  }

  for(n=0; n<3; n++)
    if(f[1]->edge[n]==edge) continue;
    else if(f[1]->edge[n]->endpts[0]==v[1]) {
      e[2]=f[1]->edge[n];
      break;
    } else if(f[1]->edge[n]->endpts[1]==v[1]) {
      e[2]=f[1]->edge[n];
      break;
    }

  if(n==3) {
    fprintf(stderr,"ERROR: could not find e[2]\n");
    fprintf(stderr,"  f[1] v=%i %i %i   e=%i (%i %i) %i (%i %i) %i (%i %i)\n",
            f[1]->vert[0]->index,f[1]->vert[1]->index,f[1]->vert[2]->index,
            f[1]->edge[0]->index,f[1]->edge[0]->endpts[0]->index,f[1]->edge[0]->endpts[1]->index,
            f[1]->edge[1]->index,f[1]->edge[1]->endpts[0]->index,f[1]->edge[1]->endpts[1]->index,
            f[1]->edge[2]->index,f[1]->edge[2]->endpts[0]->index,f[1]->edge[2]->endpts[1]->index);
    return 0;
  }

  /*  fprintf(stdout,"  local info\n");
  fprintf(stdout,"  edge = %i\n",edge->index);
  fprintf(stdout,"  vert = %i %i %i %i\n",v[0]->index,v[1]->index,v[2]->index,v[3]->index);
  fprintf(stdout,"  edge = %i %i %i %i\n",e[0]->index,e[1]->index,e[2]->index,e[3]->index);
  fprintf(stdout,"  face = %i %i\n",f[0]->index,f[1]->index);
  off_display_face_info(f[0],10);
  off_display_face_info(f[1],10);
  exit(0);*/

  return 1;

} /* off_remesh_kedge_local_info */



/**********************************************
* remesh function
*
* -remeshes two objects so that the edge length is close to elen
* -both objects are remeshed in the same way, keeping 1-to-1 correspondence
*
* -elen is the target edge length
* -maxiter should be 5-10
* -wmark, if nonzero, will do selective remeshing
*  where v->v.w is used to decide where we would remsh
*  so that if v->v.w is nonzero, function skips remeshing
*  and     if v->v.w is zero, function remeshes
* -a point is usually initialized with p.w = 0
*  so by default, remeshing occurs
*
*
**********************************************/
int off_remesh_dual_kobj(kobj *objs[], off_curvature_t crv[], int maxiter) {
  const char *fcname=__func__;
  double emin=0.5,emax=3.0; /*TODO: add parameter?*/

  kedge *e[2], *ne[2];
  kedge *erm[2][9];
  kface *f[2], *frm[2][9];
  kvert *v[2], *vrm[2][2];
  int nobj;
  int n,iter;
  const int verbose=0;
  double elen = (off_get_kobj_mean_edge_length(objs[0])+off_get_kobj_mean_edge_length(objs[1]))/2.0;
  int before_remesh = objs[0]->nvert;
  
  emin = elen/10.0;
  emax = elen*2.0;

  if(verbose>=1) {
    for(nobj=0;nobj<2;nobj++) {
      fprintf(stdout,"[%s] obj:%d  vfe = %i %i %i\n",fcname,nobj,objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);
    }
  }

#if 1
  /*store max abs curvature as .w element in the vertex*/
  for(nobj=0;nobj<2;nobj++) {
    kvert *v;
    for(v=objs[nobj]->vert; v!=NULL; v=v->next) {
      /*v->index=n++;*/
      /*if(v->index==0) abort;*/
      v->v.w = NIIK_DMAX(fabs(crv[nobj].curv1[v->index-1]),
                         fabs(crv[nobj].curv2[v->index-1]));
    }
  }
#endif 

  for(iter=1; iter<=maxiter; iter++) {
    if(verbose>=1) {
      fprintf(stdout,"[%s] remesh iteration %3i\n",fcname,iter);
      fprintf(stdout,"[%s]   checking object\n",fcname);
      if( !off_kobj_test(objs[0]) ) {
        fprintf(stderr,"ERROR: off_kobj_test during split\n");
        exit(0);
      }
      if( !off_kobj_test(objs[1]) ) {
        fprintf(stderr,"ERROR: off_kobj_test during split\n");
        exit(0);
      }
    }

    /***************************
    * SPLIT EDGES
    * 
    ****************************/
    if(verbose>=1) {
      fprintf(stdout,"[%s] split edges  \n",fcname);
    }
    
    for(e[0]=objs[0]->edge,e[1]=objs[1]->edge; e[0]!=NULL; e[0]=e[0]->next,e[1]=e[1]->next) {
      int i;
      double min_elen,max_elen, max_curvature;

      #if 0
      if(e[0]->endpts[0]->index == 0 || e[0]->endpts[1]->index == 0 ||
         e[1]->endpts[0]->index == 0 || e[1]->endpts[1]->index == 0 ) /*new point , just created*/
        max_curvature = 0.0; /*don't split, just continue*/
      else
        max_curvature = NIIK_DMAX(NIIK_DMAX( e[0]->endpts[0]->v.w, e[0]->endpts[1]->v.w),
                                  NIIK_DMAX( e[1]->endpts[0]->v.w, e[1]->endpts[1]->v.w));

      #endif

      min_elen = NIIK_DMIN( niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ),
                            niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ) );

      max_elen = NIIK_DMAX( niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ),
                            niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ) );

      if(max_elen < emax ) continue;  /*don't make edges too short?*/
      if(min_elen < emin*2 ) continue; /*don't make edges too short?*/
      /*if(min_elen*max_curvature < 10.0 && min_elen < emax ) continue;*/ /*split edges on high curvature areas*/

      if(verbose>=2) fprintf(stdout,"[%s]  split %i  %8.5f \n", fcname, e[0]->index, min_elen);

      if( (v[0]=off_remesh_kedge_split(objs[0],e[0]))==NULL ||
          (v[1]=off_remesh_kedge_split(objs[1],e[1]))==NULL ) {
        fprintf(stderr,"ERROR: off_remesh_kedge_split %i\n",e[0]->index);
        return 0;
      }

      if(verbose>=2) {
        fprintf(stdout,"[%s]   testing obj e = %i\n",fcname,e[0]->index);
        if(!off_kobj_test(objs[0])||
           !off_kobj_test(objs[1])) {
          fprintf(stderr,"ERROR: off_kobj_test during split\n");
          exit(0);
        }
      }
    } /* edge split */

    if(verbose>=1) {
      /* check the status */
      for(nobj=0;nobj<2;nobj++) {
        kface *f;
        kvert *v;
        kedge *e;
        for(v=objs[nobj]->vert,n=1; v!=NULL; v=v->next) {
          v->index=n++;
        }
        objs[nobj]->nvert=n-1;
        for(f=objs[nobj]->face,n=1; f!=NULL; f=f->next) {
          f->index=n++;
        }
        objs[nobj]->nface=n-1;
        for(e= objs[nobj]->edge,n=1; e!=NULL; e=e->next) {
          e->index=n++;
        }
        objs[nobj]->nedge=n-1;
        fprintf(stdout,"[%s]   vfe = %i %i %i\n",fcname,objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);
        fprintf(stdout,"[%s] obj:%d  mean elen = %9.5f \n",fcname,nobj,off_get_kobj_mean_edge_length(objs[nobj]));
        if(verbose>=2) {
          fprintf(stdout,"[%s]    test obj after split\n",fcname);
          if(!off_kobj_test(objs[nobj])) {
            fprintf(stderr,"ERROR: off_kobj_test\n");
            exit(0);
          }
        }
      }
    }


    /***************************
    * COLLAPSE EDGES
    ****************************/
    if(verbose>=1) {
      fprintf(stdout,"[%s] collapse edges  %8.4f \n",fcname,emin);
    }
    ne[0]=ne[1]=NULL;
    for(e[0]=objs[0]->edge,e[1]=objs[1]->edge,n=0; e[0]!=NULL; e[0]=ne[0],e[1]=ne[1],n++) {
      int i;
      int n1,n2;
      double min_elen,max_elen,max_curvature;

      ne[0]=e[0]->next;
      ne[1]=e[1]->next;

      min_elen = NIIK_DMIN( niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ),
                            niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ) );

      max_elen = NIIK_DMAX( niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ),
                            niikpt_distance(e[0]->endpts[0]->v, e[0]->endpts[1]->v ) );

      if(min_elen > emin ) continue;  /*don't make edges too long?*/
      if(max_elen > emax/2 ) continue; /*don't make edges too long?*/

      if(e[0]->endpts[0]->index==0 || e[0]->endpts[1]->index==0 ||
         e[1]->endpts[0]->index==0 || e[1]->endpts[1]->index==0 ) /*new point , just created*/
        max_curvature=1e10; /*don't collapse, pretend the curvature is high*/
      else
        max_curvature = NIIK_DMAX(NIIK_DMAX( e[0]->endpts[0]->v.w, e[0]->endpts[1]->v.w),
                                  NIIK_DMAX( e[1]->endpts[0]->v.w, e[1]->endpts[1]->v.w));

      if( min_elen*max_curvature>2.0 ) continue; /*don't collapse edges in high curvature areas*/

      if(verbose>=2) fprintf(stdout,"[%s]  collapse %i  %8.5f \n", fcname, e[0]->index, min_elen);

      if((n1=off_remesh_kedge_collapse(objs[0], e[0], vrm[0], frm[0], erm[0]))==0 ||
         (n2=off_remesh_kedge_collapse(objs[1], e[1], vrm[1], frm[1], erm[1]))==0 ) {
        fprintf(stderr,"ERROR: off_remesh_kedge_collapse \n");
        return 0;
      }

      /* the edge was not collapsed */
      /* TODO: check if only one edge collapsed?*/
      /* shouldn't hapen since topology is the same */
      if( n1<0 && n2<0 ) continue;

      /* the edge was actually collapsed */
      v[0]=erm[0][0]->endpts[0];
      v[1]=erm[1][0]->endpts[0];

      if(verbose>=4) {
        fprintf(stdout,"\t\tcollapse clean up\n");
        fprintf(stdout,"\t\t  v[%i] \n",vrm[0][0]->index);
        fprintf(stdout,"\t\t  f[%i %i] \n",frm[0][0]->index,frm[0][1]->index);
        fprintf(stdout,"\t\t  e[%i %i %i] \n",erm[0][0]->index,erm[0][1]->index,erm[0][2]->index);

        fprintf(stdout,"\t\t  v[%i] \n",vrm[1][0]->index);
        fprintf(stdout,"\t\t  f[%i %i] \n",frm[1][0]->index,frm[1][1]->index);
        fprintf(stdout,"\t\t  e[%i %i %i] \n",erm[1][0]->index,erm[1][1]->index,erm[1][2]->index);
      }

      /* make sure that the next edge is not removed */
      for(ne[0]=e[0]->next, ne[1]=e[1]->next; ne[0]!=NULL; ne[0]=ne[0]->next, ne[1]=ne[1]->next) {
        if(ne[0]==erm[0][0] || ne[0]==erm[0][1]|| ne[0]==erm[0][2]||
           ne[1]==erm[1][0] || ne[1]==erm[0][1]|| ne[1]==erm[1][2] ) continue;
        break;
      }

      if(verbose>=4) fprintf(stdout,"\t\tcollapse clean up !\n");
      /* remove and free */
      off_remesh_kedge_collapse_clean_up(objs[0],vrm[0],frm[0],erm[0]);
      off_remesh_kedge_collapse_clean_up(objs[1],vrm[1],frm[1],erm[1]);
      vrm[0][0]=NULL;
      frm[0][0]=frm[0][1]=NULL;
      erm[0][0]=erm[0][1]=erm[0][2]=NULL;
      vrm[1][0]=NULL;
      frm[1][0]=frm[1][1]=NULL;
      erm[1][0]=erm[1][1]=erm[1][2]=NULL;
      if(verbose>=4) fprintf(stdout,"\t\tcollapse clean up done\n");

      /* remove this after debugging
      * -check the object for everything at each step  */
      if(verbose>=3) {
        if(!off_kobj_test(objs[0]) || !off_kobj_test(objs[1])) {
          fprintf(stderr,"ERROR: off_kedge_test\n");
          fprintf(stderr,"       test off_remesh_kedge_collapse\n");
          exit(0);
        }
      }
    } /* each edge for collapsing */

    /* correct 3 neighbors */
    if( !off_remesh_kobj_collapse_correction(objs[0])||
        !off_remesh_kobj_collapse_correction(objs[1])   ) {
      fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction \n");
      return 0;
    }

    if(verbose>=1) {
      for(nobj=0;nobj<2;nobj++) {
        kvert *v;
        kedge *e;
        kface *f;
        if(verbose>=2) {
          for(v=objs[nobj]->vert,n=1; v!=NULL; v=v->next) {
            v->index=n++;
          }
          objs[nobj]->nvert=n-1;
          for(f=objs[nobj]->face,n=1; f!=NULL; f=f->next) {
            f->index=n++;
          }
          objs[nobj]->nface=n-1;
          for(e=objs[nobj]->edge,n=1; e!=NULL; e=e->next) {
            e->index=n++;
          }
          objs[nobj]->nedge=n-1;
        }
        off_kobj_update_num(objs[nobj]);
        fprintf(stdout,"      vfe = %i %i %i \n",objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);
        fprintf(stdout,"  obj:%d    mean elen = %9.5f \n",nobj,off_get_kobj_mean_edge_length(objs[nobj]));

        if(verbose>=2) {
          if(!off_kobj_test(objs[nobj])) {
            fprintf(stderr,"ERROR: test obj\n");
            exit(0);
          }
        }
      }
    }

    /***************************
    * VALENCE MINIMIATION
    * depends only on topology, which is the same for both surfaces
    ****************************/
   for(nobj=0;nobj<2;nobj++) {
    if(verbose>=2) {
        fprintf(stdout,"\tobj:%d  valence6  = %i\n",nobj,off_remesh_kobj_count_global_valence6(objs[nobj]));
        fprintf(stdout,"\tobj:%d  valence57 = %i\n",nobj,off_remesh_kobj_count_global_valence57(objs[nobj]));
        fprintf(stdout,"\tobj:%d  valence   = %i  variance\n",nobj,off_remesh_kobj_count_global_valence_variance(objs[nobj]));
    }

    if(verbose>=1) fprintf(stdout,"[%s] minimize the valence\n",fcname);
    for(n=0; n<3; n++) {
      if(!off_remesh_kobj_valence_minimization(objs[nobj],0)) { /*TODO: this actually can create intersections!*/
        fprintf(stderr,"ERROR: off_remesh_kobj_valence_minimization\n");
        return 0;
      }
    }
   }

   if(verbose>=1) {
      for(nobj=0;nobj<2;nobj++) {
        n = off_remesh_kobj_count_global_valence57(objs[nobj]);
        off_kobj_update_num(objs[nobj]);
        fprintf(stdout,"[%s] obj:%d  valence = %3.0f%% = %i / %i\n",fcname,nobj,100.0*n/objs[nobj]->nvert,n,objs[nobj]->nvert);

        if(verbose>=2) {
          kvert *v;
          kedge *e;
          kface *f;
          fprintf(stdout,"    check status and update index \n");
          for(v=objs[nobj]->vert,n=1; v!=NULL; v=v->next) {
            v->index=n++;
          }
          objs[nobj]->nvert=n-1;
          for(f=objs[nobj]->face,n=1; f!=NULL; f=f->next) {
            f->index=n++;
          }
          objs[nobj]->nface=n-1;
          for(e=objs[nobj]->edge,n=1; e!=NULL; e=e->next) {
            e->index=n++;
          }
          objs[nobj]->nedge=n-1;
          fprintf(stdout," obj:%d     vfe = %i %i %i \n",nobj,objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);
          fprintf(stdout," obj:%d     mean elen = %9.5f \n",nobj, off_get_kobj_mean_edge_length(objs[nobj]));
          if(!off_kobj_test(objs[nobj])) {
            fprintf(stderr,"ERROR: test obj\n");
            exit(0);
          }
        }
      }
    }

    for(nobj=0;nobj<2;nobj++) {
      if(!off_remesh_kobj_collapse_correction(objs[nobj])) {
        fprintf(stderr,"ERROR: off_remesh_kobj_collapse_correction\n");
        return 0;
      }
    }

    if(verbose>=2) {
      for(nobj=0;nobj<2;nobj++) {
        kvert *v;
        kedge *e;
        kface *f;

        fprintf(stdout,"    check status and update index \n");
        for(v=objs[nobj]->vert,n=1; v!=NULL; v=v->next) {
          v->index=n++;
        }
        objs[nobj]->nvert=n-1;
        for(f=objs[nobj]->face,n=1; f!=NULL; f=f->next) {
          f->index=n++;
        }
        objs[nobj]->nface=n-1;
        for(e=objs[nobj]->edge,n=1; e!=NULL; e=e->next) {
          e->index=n++;
        }
        objs[nobj]->nedge=n-1;
        fprintf(stdout,"  vfe = %i %i %i \n",objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);
        fprintf(stdout,"    mean elen = %9.5f \n",off_get_kobj_mean_edge_length(objs[nobj]));

        if(!off_kobj_test(objs[nobj])) {
          fprintf(stderr,"ERROR: test obj\n");
          exit(0);
        }
      }
    }
  } /* iteration */

  for(nobj=0;nobj<2;nobj++) {
    kvert *v;
    kedge *e;
    kface *f;

    for(v=objs[nobj]->vert,n=1; v!=NULL; v=v->next) {
      v->index=n++;
    }
    objs[nobj]->nvert=n-1;
    for(f=objs[nobj]->face,n=1; f!=NULL; f=f->next) {
      f->index=n++;
    }
    objs[nobj]->nface=n-1;
    for(e=objs[nobj]->edge,n=1; e!=NULL; e=e->next) {
      e->index=n++;
    }
    objs[nobj]->nedge=n-1;
    if(verbose>=1) fprintf(stdout,"[%s] remeshed to vfe = %i %i %i\n",fcname,objs[nobj]->nvert,objs[nobj]->nface,objs[nobj]->nedge);

    if(objs[nobj]->color) {
      if(verbose>=1) fprintf(stdout,"[%s] adding color (memory)\n",fcname);
      if(!off_kobj_add_color(objs[nobj])) {
        fprintf(stderr,"ERROR: off_kobj_add_color\n");
        return 0;
      }
    }
  }

  fprintf(stdout,"remeshed %d to %d\n",before_remesh, objs[0]->nvert);

  if(verbose) fprintf(stdout,"[%s] finished\n",fcname);
  return 1;
} 




/*
kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
