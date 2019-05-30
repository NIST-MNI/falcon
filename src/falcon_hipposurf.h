#ifndef _FALCON_HIPPOSURF_H_
#define _FALCON_HIPPOSURF_H_

#include "falcon.h"

typedef struct {
  nifti_image *hippomask;
  kobj *lh; // left hippo
  kobj *rh; // right hippo
  char outname[4096];
  int verbose;
  float elen;
} niik_hipposurf; // class for hippocampal surface model

niik_hipposurf *niik_hipposurf_init()
/* constructor
 */
{
  niik_hipposurf *h;
  NIIK_RETURN(((h=(niik_hipposurf *)calloc(1,sizeof(niik_hipposurf)))==NULL),"calloc for niik_hipposurf",0);
  h->verbose=1;
  h->lh=h->rh=NULL;
  h->hippomask=NULL;
  h->outname[0]=0;
  h->elen=0.5;
  return h;
} /* niik_hipposurf_init */

niik_hipposurf *niik_hipposurf_free(niik_hipposurf *h) {
  if(h==NULL) return NULL;
  if(h->hippomask!=NULL) h->hippomask=niik_image_free(h->hippomask);
  if(h->lh!=NULL) h->lh=off_kobj_free(h->lh);
  if(h->rh!=NULL) h->rh=off_kobj_free(h->rh);
  free(h);
  return NULL;
} /* niik_hipposurf_free */


int niik_hipposurf_make_surf(niik_hipposurf *h) {
  nifti_image
  *himgs[2],
  *lhimg=NULL,
   *rhimg=NULL;
  int i,n;
  double dv;
  kobj *hippo[2];
  niikpt ctr[2];
  double radius=10;
  niikmat *mat=NULL;
  niikvec *vec=NULL;
  kvert *v;

  if(h->verbose>0) niik_fc_display((char *)__func__,1);
  NIIK_RETURN((h==NULL),"hipposurf obj is null",0);
  NIIK_RETURN((h->hippomask==NULL),"hippo mask is null",0);


  /**********************************************************/
  if(h->verbose>1) {
    fprintf(stdout,"[%s:%i:%s] creating left/right hippo masks\n",__FILE__,__LINE__,__func__);
  }
  NIIK_RETURN(((lhimg=niik_image_copy_as_type(h->hippomask,NIFTI_TYPE_UINT8))==NULL),"niik_image_copy_as_type",0);
  NIIK_RETURN(((rhimg=niik_image_copy_as_type(h->hippomask,NIFTI_TYPE_UINT8))==NULL),"niik_image_copy_as_type",0);
  niik_image_clear(lhimg);
  niik_image_clear(rhimg);
  for(i=0; i<h->hippomask->nvox; i++) {
    dv=niik_image_get_voxel(h->hippomask,i);
    if(dv>1.5 && dv<2.5) {
      niik_image_set_voxel(rhimg,i,1);
    } else if(dv>3.5 && dv<4.5) {
      niik_image_set_voxel(lhimg,i,1);
    }
  }
  if(niik_image_count_mask(lhimg)==0) {
    for(i=0; i<h->hippomask->nvox; i++) {
      dv=niik_image_get_voxel(h->hippomask,i);
      if(dv>2.5 && dv<3.5) {
        niik_image_set_voxel(rhimg,i,1);
      } else if(dv>15.5 && dv<16.5) {
        niik_image_set_voxel(lhimg,i,1);
      }
    }
  }
  NIIK_RETURN((niik_image_count_mask(lhimg)==0),"no left hc",0);
  NIIK_RETURN((niik_image_count_mask(rhimg)==0),"no right hc",0);

  himgs[0]=lhimg;
  himgs[1]=rhimg;

  if(h->verbose>2) {
    fprintf(stdout,"[%s:%i:%s] writing images\n",__FILE__,__LINE__,__func__);
    niik_image_write("lhimg.nii.gz",lhimg);
    niik_image_write("rhimg.nii.gz",rhimg);
  }


  /**********************************************************/
  if(h->verbose>1) {
    fprintf(stdout,"[%s:%i:%s] creating left/right hippo objects\n",__FILE__,__LINE__,__func__);
  }
  for(n=0; n<2; n++) {
    ctr[n]=niikpt_image_get_centroid(himgs[n],NULL);
    fprintf(stdout,"[%s] centroid = %7.3f %7.3f %7.3f\n",__func__,ctr[n].x,ctr[n].y,ctr[n].z);
    mat=niikmat_affine_matrix_new(60,0,0,  0,0,0,  2,2,5,  0,0,0,   ctr[n].x,ctr[n].y,ctr[n].z);
    NIIK_RETURN(((hippo[n]=off_make_sphere_from_icosahedron(2,radius,ctr[n]))==NULL),"off_make_sphere_from_icosahedron",0);
    for(v=hippo[n]->vert; v!=NULL; v=v->next) {
      v->v=niikpt_affine_transform(mat,v->v);
    }
  }
  h->lh=hippo[0];
  h->rh=hippo[1];

  if(h->verbose>2) { // write the initial elliptical model
    for(n=0; n<2; n++) {
      vec=niikvec_init(hippo[n]->nvert);
      for(v=hippo[n]->vert,i=0; v!=NULL; v=v->next,i++) {
        vec->v[i]=v->sph.the;
      }
      fprintf(stdout,"%i  %f %f\n",n,niikvec_get_min(vec),niikvec_get_max(vec));
      NIIK_RETURN((!off_kobj_apply_color_map(hippo[n],vec->v,
                                             niikvec_get_min(vec),niikvec_get_max(vec),
                                             NIIK_COLORMAP_SPECTRAL)),
                  "off_kobj_apply_color_map",1);
      vec=niikvec_free(vec);
    }
    off_kobj_write_off("lh0.off",h->lh,0);
    off_kobj_write_off("rh0.off",h->rh,0);
  }

  #pragma omp parallel for
  for(n=0; n<2; n++) {
    /* shrink-wrap roughly -> finely, then remesh to target */
    off_shrinkwrap_kobj_bbox_remesh(himgs[n],hippo[n],
                                    1.5,20,2,1,1,0);
    off_shrinkwrap_kobj_bbox_remesh(himgs[n],hippo[n],
                                    h->elen,20,2,1,1,0);
    off_shrinkwrap_kobj_bbox_remesh(himgs[n],hippo[n],h->elen/1.5,10,2,1,1,0);
    off_remesh_kobj(hippo[n],h->elen,10,0);
  }


  /* put colors */
  for(n=0; n<2; n++) {
    vec=niikvec_init(hippo[n]->nvert);
    for(v=hippo[n]->vert,i=0; v!=NULL; v=v->next,i++) {
      vec->v[i]=v->sph.the;
    }
    //fprintf(stdout,"%i  %f %f\n",n,niikvec_get_min(vec),niikvec_get_max(vec));
    NIIK_RETURN((!off_kobj_apply_color_map(hippo[n],vec->v,
                                           //0,NIIK_PI,
                                           niikvec_get_min(vec),niikvec_get_max(vec),
                                           NIIK_COLORMAP_SPECTRAL)),
                "off_kobj_apply_color_map",1);
    vec=niikvec_free(vec);
  }


  /**********************************************************/
  // WRITING OBJECTS
  if(h->verbose>2) {
    fprintf(stdout,"[%s:%i:%s] writing objects\n",__FILE__,__LINE__,__func__);
    off_kobj_write_off("lh.off",h->lh,0);
    off_kobj_write_off("rh.off",h->rh,0);
  }

  lhimg=niik_image_free(lhimg);
  rhimg=niik_image_free(rhimg);
  if(h->verbose>0) niik_fc_display((char *)__func__,0);
  return 1;
}

double angle_diff(double v,double w) {
  double d;
  d=v-w;
  if(fabs(d)<NIIK_PI2) return d;
  if(d>0) return d-NIIK_PI2;
  return d+NIIK_PI2;
}

double angle_diff2(double the1,double psi1,double the2,double psi2) {
  /*return NIIK_RAD2DEGREE(niikpt_distance(niikpt_val(sin(the1)*cos(psi1),sin(the1)*sin(psi1),cos(the1),1),
    niikpt_val(sin(the2)*cos(psi2),sin(the2)*sin(psi2),cos(the2),1)));*/
  return niikpt_angle_between_vectors(niikpt_val(sin(the1)*cos(psi1),
                                      sin(the1)*sin(psi1),
                                      cos(the1),
                                      1),
                                      niikpt_val(sin(the2)*cos(psi2),
                                          sin(the2)*sin(psi2),
                                          cos(the2),
                                          1));
}

kvert *niik_off_resample_parametric_surface_find_vert1(kobj *obj,double the,double psi) {
  kvert *v,*v1=NULL;
  double d,d1;
  for(v=v1=obj->vert,d1=10000; v!=NULL; v=v->next) {
    d=angle_diff2(the,psi,v->sph.the,v->sph.psi);
    if(d<d1) {
      d1=d;
      v1=v;
    }
    if(0)fprintf(stdout,"    %9.1f %9.1f %9.1f %9.1f | %4.0f\n",
                   NIIK_RAD2DEGREE(the),NIIK_RAD2DEGREE(psi),
                   NIIK_RAD2DEGREE(v->sph.the),NIIK_RAD2DEGREE(v->sph.psi),d);
  }
  return v1;
}

kvert *niik_off_resample_parametric_surface_find_vert(kobj *obj,kvert *vi,double the,double psi) {
  kvert *v=NULL,*vlist[1024];
  int nei,neimin,num;
  double es[1024],emin;
  int const verbose=0;
  if(verbose>0)NIIK_RETURN((obj==NULL),"obj is null",0);
  if(vi==NULL) v=obj->vert;
  else v=vi;
  nei=0;
  while(nei>=0) {
    neimin=-1;
    //emin=fabs(angle_diff(v->sph.the,the)) + fabs(angle_diff(v->sph.psi,psi));
    emin=angle_diff2(the,psi,v->sph.the,v->sph.psi);
    for(nei=0; nei<v->nei; nei++) {
      //es[nei] = fabs(angle_diff(v->neivert[nei]->sph.the,the)) + fabs(angle_diff(v->neivert[nei]->sph.psi,psi));
      //es[nei] = acos(cos(v->neivert[nei]->sph.the)*cos(v->sph.the) + sin(v->neivert[nei]->sph.the)*sin(v->sph.the)*cos(v->neivert[nei]->sph.psi-v->sph.psi));
      es[nei] = angle_diff2(the,psi,v->neivert[nei]->sph.the,v->neivert[nei]->sph.psi);
      if(verbose>0)fprintf(stdout,"  nei %i %6.3f |  %6.2f %6.2f   %6.2f %6.2f\n",nei,es[nei],the,psi,v->neivert[nei]->sph.the,v->neivert[nei]->sph.psi);
      if(es[nei]<emin) {
        emin=es[nei];
        neimin=nei;
      }
    }
    if(neimin<0) {
      num=1;
      vlist[0]=v;
      NIIK_RETURN(((num=off_get_local_kvert_check(vlist,1024,num))<0),"off_get_local_kvert_check",NULL);
      NIIK_RETURN(((num=off_get_local_kvert_check(vlist,1024,num))<0),"off_get_local_kvert_check",NULL);
      NIIK_RETURN(((num=off_get_local_kvert_check(vlist,1024,num))<0),"off_get_local_kvert_check",NULL);
      emin=angle_diff2(the,psi,v->sph.the,v->sph.psi);
      neimin=-1;
      for(nei=0; nei<num; nei++) {
        es[nei] = angle_diff2(the,psi,vlist[nei]->sph.the,vlist[nei]->sph.psi);
        if(verbose>0)fprintf(stdout,"  nei %i %6.3f |  %6.2f %6.2f   %6.2f %6.2f\n",nei,es[nei],the,psi,vlist[nei]->sph.the,vlist[nei]->sph.psi);
        if(es[nei]<emin) {
          emin=es[nei];
          neimin=nei;
        }
      }
      if(neimin<0) {
        return v;
      } else {
        //fprintf(stdout,"[%s] running second niik_off_resample_parametric_surface_find_vert\n",__func__);
        return niik_off_resample_parametric_surface_find_vert(obj,vlist[neimin],the,psi);
      }
    } else  {
      if(verbose>0)fprintf(stdout,"  nei %i %6.3f |  %6.2f %6.2f   %6.2f %6.2f *\n\n",neimin,es[neimin],the,psi,v->neivert[neimin]->sph.the,v->neivert[neimin]->sph.psi);
      v=v->neivert[neimin];
    }
  }
  return NULL;
}


kface *niik_off_resample_parametric_surface_find_face(kobj *obj,kvert *v,double the,double psi) {
  int nei,neimin,n;
  double es[48];
  for(nei=neimin=0; nei<v->nei; nei++) {
    es[nei]=0;
    for(n=0; n<3; n++) {
      es[nei] +=
        fabs(angle_diff(v->neiface[nei]->vert[n]->sph.the,the)) +
        fabs(angle_diff(v->neiface[nei]->vert[n]->sph.psi,psi));
    }
    if(es[nei]<es[neimin]) {
      neimin=nei;
    }
  }
  return v->neiface[neimin];
}

int niik_off_resample_parametric_surface(kobj *obj,kobj *ref,niikmat *vin,niikmat *vout) {
  kvert *vo,*vr;
  int i,n;
  int const verbose=0;
  niik_fc_display((char *)__func__,1);
  NIIK_RETURN((obj==NULL),"obj is null",0);
  NIIK_RETURN((ref==NULL),"ref is null",0);
  NIIK_RETURN((vin ==NULL),"vin is null",0);
  NIIK_RETURN((vout==NULL),"vout is null",0);
  NIIK_RETURN((vin->row!=vout->row),"wrong #row",0);
  NIIK_RETURN((vout->col!=ref->nvert),"wrong #row",0);
  NIIK_RETURN((vin ->col!=obj->nvert),"wrong #row",0);
  if(verbose>0)fprintf(stdout,"[%s:%i:%s] re-index\n",__FILE__,__LINE__,__func__);
  for(vo=obj->vert,i=0; vo!=NULL; vo=vo->next,i++) {
    vo->index=i;
  }
  for(vr=ref->vert,i=0; vo!=NULL; vr=vr->next,i++) {
    vr->index=i;
  }

  vo=NULL;
  if(verbose>0)fprintf(stdout,"[%s:%i:%s] analysis\n",__FILE__,__LINE__,__func__);
  for(vr=ref->vert,i=0; vr!=NULL; vr=vr->next,i++) {
    if(verbose>0)fprintf(stdout,"[%s:%i:%s] analysis %i\n",__FILE__,__LINE__,__func__,i);
    vo=niik_off_resample_parametric_surface_find_vert1(obj,vr->sph.the,vr->sph.psi);
    if(verbose>0) {
      fprintf(stdout,"[%s:%i:%s] analysis %i %i\n",__FILE__,__LINE__,__func__,i,vo->index);
      if(verbose>1) {
        fprintf(stdout,"%6.2f %6.2f %6.2f    %6.2f %6.2f %6.2f  |   %8.1f %8.1f      %8.1f %8.1f  |  %6.2f %6.2f %6.2f %8.1f %8.1f | %8.1f *\n",
                vr->v.x,vr->v.y,vr->v.z,
                vo->v.x,vo->v.y,vo->v.z,
                NIIK_RAD2DEGREE(vr->sph.the),NIIK_RAD2DEGREE(vr->sph.psi),NIIK_RAD2DEGREE(vo->sph.the),NIIK_RAD2DEGREE(vo->sph.psi),
                vo->v.x-vr->v.x,
                vo->v.y-vr->v.y,
                vo->v.z-vr->v.z,
                NIIK_RAD2DEGREE(vo->sph.the)-NIIK_RAD2DEGREE(vr->sph.the),
                NIIK_RAD2DEGREE(vo->sph.psi)-NIIK_RAD2DEGREE(vr->sph.psi),
                NIIK_RAD2DEGREE(niikpt_distance(niikpt_val(sin(vo->sph.the)*cos(vo->sph.psi),sin(vo->sph.the)*sin(vo->sph.psi),cos(vo->sph.the),1),
                                                niikpt_val(sin(vr->sph.the)*cos(vr->sph.psi),sin(vr->sph.the)*sin(vr->sph.psi),cos(vr->sph.the),1))));
      }
    }

    for(n=0; n<vin->row; n++) {
      vout->m[n][i]=vin->m[n][vo->index];
    }
  }
  niik_fc_display((char *)__func__,0);
  return 1;
} /* niik_off_resample_parametric_surface */

int niik_hipposurf_principal_axes(kobj *obj,niikpt *pa,niikpt ctr) {
  niikmat *cov=NULL;
  kvert *v;
  niikvec *eigval;
  niikmat *eigvec;
  int j;
  cov=niikmat_init(3,3);
  eigval=niikvec_init(3);
  eigvec=niikmat_init(3,3);
  for(v=obj->vert; v!=NULL; v=v->next) {
    cov->m[0][0] += (v->v.x-ctr.x)*(v->v.x-ctr.x);
    cov->m[0][1] += (v->v.x-ctr.x)*(v->v.y-ctr.y);
    cov->m[0][2] += (v->v.x-ctr.x)*(v->v.z-ctr.z);
    cov->m[1][0] += (v->v.y-ctr.y)*(v->v.x-ctr.x);
    cov->m[1][1] += (v->v.y-ctr.y)*(v->v.y-ctr.y);
    cov->m[1][2] += (v->v.y-ctr.y)*(v->v.z-ctr.z);
    cov->m[2][0] += (v->v.z-ctr.z)*(v->v.x-ctr.x);
    cov->m[2][1] += (v->v.z-ctr.z)*(v->v.y-ctr.y);
    cov->m[2][2] += (v->v.z-ctr.z)*(v->v.z-ctr.z);
  }
  niikmat_kmul(cov,1.0/obj->nvert);
  NIIK_RETURN((!niikmat_jacobi(cov,eigval,eigvec)),"niikmat_jacobi",0);
  /*niikmat_display(eigvec);
    niikvec_display(eigval);*/
  if(fabs(eigval->v[0])>fabs(eigval->v[1]) && fabs(eigval->v[0])>fabs(eigval->v[2])) {
    j=0;
    pa[0]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    if(fabs(eigval->v[1])>fabs(eigval->v[2])) {
      j=1;
      pa[1]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
      j=2;
      pa[2]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    }
  } else if(fabs(eigval->v[1])>fabs(eigval->v[0]) && fabs(eigval->v[1])>fabs(eigval->v[2])) {
    j=1;
    pa[0]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    if(fabs(eigval->v[0])>fabs(eigval->v[2])) {
      j=0;
      pa[1]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
      j=2;
      pa[2]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    }
  } else {
    j=2;
    pa[0]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    if(fabs(eigval->v[0])>fabs(eigval->v[1])) {
      j=0;
      pa[1]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
      j=1;
      pa[2]=niikpt_val(eigvec->m[0][j],eigvec->m[1][j],eigvec->m[2][j],1);
    }
  }
  pa[3]=niikpt_val(eigval->v[0],eigval->v[1],eigval->v[2],0);
  eigval=niikvec_free(eigval);
  eigvec=niikmat_free(eigvec);
  cov=niikmat_free(cov);
  return 1;
}

int niik_hipposurf_extract_feature(kobj *obj,niikmat *mat) {
  int i,j,nei;
  kvert *v;
  niikpt ctr,*pa;
  double dsum,dvol,dvolabs,dval,dsph,dvec;
  const char *fcname="niik_hipposurf_extract_feature";

  off_update_kobj_face_normal(obj);
  off_update_kobj_vert_normal(obj);
  off_smooth_kobj_vert_normal(obj);
  off_update_kobj_kface_pminmax(obj);

  for(v=obj->vert,ctr=niikpt_zero(); v!=NULL; v=v->next) {
    ctr=niikpt_add(ctr,v->v);
  }
  ctr=niikpt_kmul(ctr,1.0/obj->nvert);
  pa=(niikpt *)calloc(4,sizeof(niikpt));
  niik_hipposurf_principal_axes(obj,pa,ctr);
  fprintf(stdout,"[%s] principal axis = %9.4f %9.4f %9.4f\n",__func__,pa[0].x,pa[0].y,pa[0].z);
  for(i=0; i<3; i++) pa[i]=niikpt_unit(pa[i]);

  fprintf(stdout,"[%s] resampling for features\n",fcname);
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    mat->m[0][i]=v->v.x; // coordinate x
    mat->m[1][i]=v->v.y; // coordinate y
    mat->m[2][i]=v->v.z; // coordinate z
    mat->m[3][i]=niikpt_distance(v->v,ctr); // distance from centroid
    for(nei=0,dsum=dvec=0,dvol=dvolabs=0; nei<v->nei; nei++) { // local surface area
      // areas
      mat->m[4][i]+=niikpt_area2(v->neiface[nei]->vert[0]->v,
                                 v->neiface[nei]->vert[1]->v,
                                 v->neiface[nei]->vert[2]->v);
      dsph = acos(cos(v->neivert[nei]->sph.the)*cos(v->sph.the) + sin(v->neivert[nei]->sph.the)*sin(v->sph.the)*cos(v->neivert[nei]->sph.psi-v->sph.psi));
      dsum += dsph;
      dvec += niikpt_distance(v->v,v->neivert[nei]->v) / dsph;
      // volumes
      dval = niikpt_det4(v->neiface[nei]->vert[0]->v.x,v->neiface[nei]->vert[1]->v.x,v->neiface[nei]->vert[2]->v.x,ctr.x,
                         v->neiface[nei]->vert[0]->v.y,v->neiface[nei]->vert[1]->v.y,v->neiface[nei]->vert[2]->v.y,ctr.y,
                         v->neiface[nei]->vert[0]->v.z,v->neiface[nei]->vert[1]->v.z,v->neiface[nei]->vert[2]->v.z,ctr.z,
                         1,1,1,1);
      dvol+=dval;
      dvolabs+=fabs(dval);
    }
    mat->m[4][i]/=nei;
    mat->m[5][i]=mat->m[4][i]/dsum;
    niik_off_curvature_vert(v, 5,&mat->m[6][i]);
    niik_off_curvature_vert(v,11,&mat->m[7][i]);
    // niik_off_curvature_vert_meek2000(v,&mat->m[6][i]);
    mat->m[8][i]=dvol;
    mat->m[9][i]=dvolabs;
    mat->m[10][i]=dvec/nei;
    mat->m[11][i]=niikpt_dot(niikpt_sub(v->v,ctr),pa[0]);
    mat->m[12][i]=niikpt_dot(niikpt_sub(v->v,ctr),pa[1]);
    mat->m[13][i]=niikpt_dot(niikpt_sub(v->v,ctr),pa[2]);
    mat->m[14][i]=NIIK_SSQ(mat->m[12][i],mat->m[13][i]);
  }

  // normalize the curvatures
  j=6;
  for(v=obj->vert,dsum,i=0; v!=NULL; v=v->next,i++) {
    dsum+=fabs(mat->m[j][i]);
  }
  dsum/=i;
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    mat->m[j][i]/=dsum;
  }
  j=7;
  for(v=obj->vert,dsum,i=0; v!=NULL; v=v->next,i++) {
    dsum+=mat->m[j][i];
  }
  dsum/=i;
  for(v=obj->vert,i=0; v!=NULL; v=v->next,i++) {
    mat->m[j][i]/=dsum;
  }

  i=6;
  for(j=0; j<5; j++) off_surface_smooth_using_vert(obj,mat->m[i],1,0);
  i=7;
  for(j=0; j<5; j++) off_surface_smooth_using_vert(obj,mat->m[i],1,0);

  return 1;
}

#endif
