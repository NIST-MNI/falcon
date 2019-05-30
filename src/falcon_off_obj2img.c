/* Filename:     nifti1_kunio_off_obj2img.c
 * Description:  object functions for creating bounds image
 * Author:       Kunio Nakamura
 * Date:         March 1, 2012
 */


#ifndef _FALCON_OFF_OBJ2IMG_C_
#define _FALCON_OFF_OBJ2IMG_C_

#include "falcon.h"
#include "falcon_surfaces.h"

/*****************************************************
 *
 * off_obj2img
 *
 * -object to image (boundary)
 *
 *****************************************************/

int off_obj2img_use_qform(nifti_image *img,kobj *obj,double pval) {
  niikmat *afmat;
  kvert *v;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(!img->qform_code) {
    fprintf(stderr,"ERROR: img->qform_code is zero\n");
    return 0;
  }
  if((afmat = niikmat_mat44_matrix(img->qto_ijk))==NULL) {
    fprintf(stderr,"ERROR: niikmat_mat44_matrix(img->qto_ijk)\n");
    return 0;
  }
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v = niikpt_affine_transform(afmat,v->v);
  }
  if(!off_obj2img(img,obj,pval)) {
    fprintf(stderr,"ERROR: off_obj2img(img,obj,pval)\n");
    return 0;
  }
  if(!niikmat_inverse_update(afmat)) {
    fprintf(stderr,"ERROR: niikmat_inverse_update(afmat)\n");
    return 0;
  }
  for(v=obj->vert; v!=NULL; v=v->next) {
    v->v = niikpt_affine_transform(afmat,v->v);
  }
  afmat=niikmat_free(afmat);
  return 1;
}

int off_obj2img(nifti_image *img,kobj *obj,double pval)
/* create a boundary image
 * off_obj2img_color works better */
{
  kface *f;
  kedge *e;
  kvert *v;
  niikpt
  p,vec;
  int
  n;
  double
  delta,
  dist,dval,dval2,
  dp[3],
  *idata;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(img->nu>1) {
    fprintf(stderr,"ERROR: img->nu is %i > 1\n",img->nu);
    return 0;
  }
  if(img->nt>1) {
    fprintf(stderr,"ERROR: img->nt is %i > 1\n",img->nt);
    return 0;
  }
  if((idata=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  delta = NIIK_DMIN(img->dx,NIIK_DMIN(img->dy,img->dz)) * 0.2;
  /* put value on the vertex */
  for(v=obj->vert; v!=NULL; v=v->next) {
    if(v->v.x<0) continue;
    if(v->v.y<0) continue;
    if(v->v.z<0) continue;
    n=niik_image_check_index_niikpt(img,v->v);
    if(n>=0) idata[n]=pval;
  }
  /* put value along the edge */
  for(e=obj->edge; e!=NULL; e=e->next) {
    vec  = niikpt_sub(e->endpts[0]->v,e->endpts[1]->v);
    dist = niikpt_mag(vec);
    vec  = niikpt_unit(vec);
    for(dval=delta; dval<=dist; dval+=delta) {
      p = niikpt_move_normal(e->endpts[1]->v,vec,dval);
      if(p.x<0) continue;
      n=niik_image_check_index_niikpt(img,p);
      if(n>=0) idata[n]=pval;
    }
  }
  for(f=obj->face; f!=NULL; f=f->next) {
    p=niikpt_avg3(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v);
    n=niik_image_check_index_niikpt(img,p);
    if(n>=0) idata[n]=pval;
    /* -we want to normalize so that 1.0 = average length
     * -estimate the step size by 1.0 / dval * delta
     * -approximately 5 steps par voxel */
    dval = (niikpt_distance(f->vert[0]->v,f->vert[1]->v) + niikpt_distance(f->vert[0]->v,f->vert[2]->v))/2.0;
    dval = 1.0 / dval * delta;
    dval2 = niikpt_distance(f->vert[1]->v,f->vert[2]->v);
    dval2 = 1.0 / dval2 * delta;
    for(dp[0]=delta; dp[0]<1.0; dp[0]+=dval2) {
      for(dp[1]=delta; dp[1]<1.0; dp[1]+=dval) {
        /* fprintf(stderr,"\n  %5.1f %5.1f %5.1f   %5.1f %5.1f %5.1f     %5.1f %5.1f %5.1f\n",
        	f->vert[0]->v.x,f->vert[0]->v.y,f->vert[0]->v.z,
        	f->vert[1]->v.x,f->vert[1]->v.y,f->vert[1]->v.z,
        	f->vert[2]->v.x,f->vert[2]->v.y,f->vert[2]->v.z);*/
        p = niikpt_interp_tri(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,dp[0],dp[1]);
        /* fprintf(stderr,"  %5.1f %5.1f %5.1f  from %5.2f %5.2f\n",p.x,p.y,p.z,dp[0],dp[1]); */
        n=niik_image_check_index_niikpt(img,p);
        if(n>=0) idata[n]=pval;
      }
    }
  }
  /* put the values in the image */
  if(!niik_image_set_voxels_from_double_vector(img,idata)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(idata);
  return 1;
}


int off_obj2img_color(nifti_image *img,kobj *obj,double pval)
/* 2012-03-13 Kunio
 * -create a boundary image with colors
 * -color info comes from the object, pval is used to scale
 * -if the object does not have a color, then grayscale image is written
 */
{
  char fcname[32]="off_obj2img_color";
  bbox *bb;
  kface *f;
  /*kedge *e;*/
  niikpt
  tri[3],
      is1,is2,
      p;
  int
  verbose=0,
  cop,
  i,j,k,nn,n,mm=-1,m,size;
  double
  dval,dd,delta,
       *idata,
       *dx,*dy,*dz;
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0((obj==NULL),fcname,"obj is null");
  NIIK_RET0((img->nu>1),fcname,"img->nu greater than 1");
  NIIK_RET0((img->nt>1),fcname,"img->nt greater than 1");
  delta=NIIK_DMIN(img->dx,NIIK_DMIN(img->dy,img->dz))*0.1;
  /*fprintf(stdout,"\tdelta = %8.4f\n",delta);
  for(e=obj->edge,dval=0;e!=NULL;e=e->next){
    dval=NIIK_DMAX(dval,niikpt_distance(e->endpts[0]->v,e->endpts[1]->v)); }
  fprintf(stdout,"\telenM = %8.4f\n",dval);
  dd=delta/dval;
  fprintf(stdout,"\tddist = %8.4f\n",dd);*/
  /*
   * OBJECT WITHOUT COLOR
   */
  if(!obj->color) {
    size=img->nx*img->ny*img->nz;
    if((idata=niik_image_get_voxels_as_double_vector(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
      return 0;
    }
    off_update_kobj_kface_pminmax(obj);
    /* prepare bounding box */
    if((bb=off_bbox_init(7,360))==NULL) {
      fprintf(stderr,"ERROR: off_bbox_init\n");
      return 0;
    }
    if(!off_create_bbox_from_kobj(bb,obj)) {
      fprintf(stderr,"ERROR: off_create_bbox_from_kobj(bb,obj)\n");
      return 0;
    }
    /* find triangle-triangle intersections */
    if(verbose>0) fprintf(stdout,"[%s] find tri-tri intersect for colors  %8.3f\n",fcname,pval);
    /* z-direction processing */
    for(k=0; k<img->nz; k++) {
      tri[0]=tri[1]=tri[2]=niikpt_zero();
      tri[0].x=tri[0].y=0;
      tri[1].x=tri[2].y=NIIK_DMAX(img->nx*img->dx,img->ny*img->dy)*3;
      tri[0].z=tri[1].z=tri[2].z=k*img->dz;
      for(nn=0; nn<bb->nvox; nn++) {
        for(n=0; n<bb->ndata[nn]; n++) {
          f=bb->data[nn][n];
          if(f->color==NULL) continue;
          if(f->pmax.z<tri[0].z) continue;
          if(f->pmin.z>tri[0].z) continue;
          if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
              tri[0],tri[1],tri[2],
              &cop,&is1,&is2)) {
            continue;
          }
          if(cop) continue;
          if((m=niik_image_check_index_niikpt(img,is2))>=0) {
            idata[m]=pval;
          }
          dd=delta/niikpt_distance(is1,is2);
          for(dval=0; dval<=1.0; dval+=dd) {
            p=niikpt_wavg(is1,is2,dval);
            m=niik_image_check_index_niikpt(img,p);
            if(mm==m) continue;
            if(m<0) continue;
            idata[m]=pval;
            mm=m;
          } /* add color along intersect line */
        }
      } /* each face in bounding box */
    }
    /* y-direction processing */
    for(j=0; j<img->ny; j++) {
      tri[0]=tri[1]=tri[2]=niikpt_zero();
      tri[0].x=tri[0].z=0;
      tri[1].x=tri[2].z=NIIK_DMAX(img->nx*img->dx,img->nz*img->dz)*3;
      tri[0].y=tri[1].y=tri[2].y=j*img->dy;
      for(nn=0; nn<bb->nvox; nn++) {
        for(n=0; n<bb->ndata[nn]; n++) {
          f=bb->data[nn][n];
          if(f->pmax.y<tri[0].y) continue;
          if(f->pmin.y>tri[0].y) continue;
          if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
              tri[0],tri[1],tri[2],
              &cop,&is1,&is2)) {
            continue;
          }
          if(cop) continue;
          if((m=niik_image_check_index_niikpt(img,is2))>=0) {
            idata[m]=pval;
          }
          dd=delta/niikpt_distance(is1,is2);
          for(dval=0; dval<=1.0; dval+=dd) {
            p=niikpt_wavg(is1,is2,dval);
            m=niik_image_check_index_niikpt(img,p);
            if(mm==m) continue;
            if(m<0) continue;
            idata[m]=pval;
            mm=m;
          } /* add color along intersect line */
        }
      } /* each face in bounding box */
    }
    /* x-direction processing */
    for(i=0; i<img->nx; i++) {
      tri[0]=tri[1]=tri[2]=niikpt_zero();
      tri[0].y=tri[0].z=0;
      tri[1].y=tri[2].z=NIIK_DMAX(img->ny*img->dy,img->nz*img->dz)*3;
      tri[0].x=tri[1].x=tri[2].x=i*img->dx;
      for(nn=0; nn<bb->nvox; nn++) {
        for(n=0; n<bb->ndata[nn]; n++) {
          f=bb->data[nn][n];
          if(f->pmax.x<tri[0].x) continue;
          if(f->pmin.x>tri[0].x) continue;
          if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
              tri[0],tri[1],tri[2],
              &cop,&is1,&is2)) {
            continue;
          }
          if(cop) continue;
          if((m=niik_image_check_index_niikpt(img,is2))>=0) {
            idata[m]=pval;
          }
          dd=delta/niikpt_distance(is1,is2);
          for(dval=0; dval<=1.0; dval+=dd) {
            p=niikpt_wavg(is1,is2,dval);
            m=niik_image_check_index_niikpt(img,p);
            if(mm==m) continue;
            if(m<0) continue;
            idata[m]=pval;
            mm=m;
          } /* add color along intersect line */
        }
      } /* each face in bounding box */
    }
    /* put the values in the image */
    if(!niik_image_set_voxels_from_double_vector(img,idata)) {
      fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
      return 0;
    }
    free(idata);
    bb=off_bbox_free(bb);
    return 1;
  } /* OBJECT WITHOUT COLOR */

  /*
   * OBJECT WITH COLOR
   */
  if(!niik_image_convert_to_color_image(img)) {
    fprintf(stderr,"ERROR: niik_image_convert_to_color_image(img)\n");
    return 0;
  }
  if((idata=niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  size=img->nx*img->ny*img->nz;
  dx=idata;
  dy=dx+size;
  dz=dy+size;
  off_update_kobj_kface_pminmax(obj);
  /* prepare bounding box */
  if((bb=off_bbox_init(7,360))==NULL) {
    fprintf(stderr,"ERROR: off_bbox_init\n");
    return 0;
  }
  if(!off_create_bbox_from_kobj(bb,obj)) {
    fprintf(stderr,"ERROR: off_create_bbox_from_kobj(bb,obj)\n");
    return 0;
  }
  /* find triangle-triangle intersections */
  fprintf(stdout,"  find tri-tri intersect for colors  %8.3f\n",pval);
  /* z-direction processing */
  for(k=0; k<img->nz; k++) {
    tri[0]=tri[1]=tri[2]=niikpt_zero();
    tri[0].x=tri[0].y=0;
    tri[1].x=tri[2].y=NIIK_DMAX(img->nx*img->dx,img->ny*img->dy)*3;
    tri[0].z=tri[1].z=tri[2].z=k*img->dz;
    for(nn=0; nn<bb->nvox; nn++) {
      for(n=0; n<bb->ndata[nn]; n++) {
        /*fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i\n",k,ii,jj,kk,nn,n);*/
        f=bb->data[nn][n];
        if(f->color==NULL) continue;
        if(f->pmax.z<tri[0].z) continue;
        if(f->pmin.z>tri[0].z) continue;
        if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
            tri[0],tri[1],tri[2],
            &cop,&is1,&is2)) {
          continue;
        }
        if(cop) continue;
        /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop\n",k,ii,jj,kk,nn,n); */
        if((m=niik_image_check_index_niikpt(img,is2))>=0) {
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
        }
        dd=delta/niikpt_distance(is1,is2);
        for(dval=0; dval<=1.0; dval+=dd) {
          p=niikpt_wavg(is1,is2,dval);
          m=niik_image_check_index_niikpt(img,p);
          if(mm==m) continue;
          if(m<0) continue;
          /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop; voxel %8i\n",k,ii,jj,kk,nn,n,m); */
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
          mm=m;
        } /* add color along intersect line */
      }
    } /* each face in bounding box */
  }
  /* y-direction processing */
  for(j=0; j<img->ny; j++) {
    tri[0]=tri[1]=tri[2]=niikpt_zero();
    tri[0].x=tri[0].z=0;
    tri[1].x=tri[2].z=NIIK_DMAX(img->nx*img->dx,img->nz*img->dz)*3;
    tri[0].y=tri[1].y=tri[2].y=j*img->dy;
    for(nn=0; nn<bb->nvox; nn++) {
      for(n=0; n<bb->ndata[nn]; n++) {
        /*fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i\n",k,ii,jj,kk,nn,n);*/
        f=bb->data[nn][n];
        if(f->pmax.y<tri[0].y) continue;
        if(f->pmin.y>tri[0].y) continue;
        if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
            tri[0],tri[1],tri[2],
            &cop,&is1,&is2)) {
          continue;
        }
        if(cop) continue;
        /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop\n",k,ii,jj,kk,nn,n); */
        if((m=niik_image_check_index_niikpt(img,is2))>=0) {
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
        }
        dd=delta/niikpt_distance(is1,is2);
        for(dval=0; dval<=1.0; dval+=dd) {
          p=niikpt_wavg(is1,is2,dval);
          m=niik_image_check_index_niikpt(img,p);
          if(mm==m) continue;
          if(m<0) continue;
          /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop; voxel %8i\n",k,ii,jj,kk,nn,n,m); */
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
          mm=m;
        } /* add color along intersect line */
      }
    } /* each face in bounding box */
  }
  /* x-direction processing */
  for(i=0; i<img->nx; i++) {
    tri[0]=tri[1]=tri[2]=niikpt_zero();
    tri[0].y=tri[0].z=0;
    tri[1].y=tri[2].z=NIIK_DMAX(img->ny*img->dy,img->nz*img->dz)*3;
    tri[0].x=tri[1].x=tri[2].x=i*img->dx;
    for(nn=0; nn<bb->nvox; nn++) {
      for(n=0; n<bb->ndata[nn]; n++) {
        /*fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i\n",k,ii,jj,kk,nn,n);*/
        f=bb->data[nn][n];
        if(f->pmax.x<tri[0].x) continue;
        if(f->pmin.x>tri[0].x) continue;
        if(!off_check_tri_tri_intersect_with_isectline(f->vert[0]->v,f->vert[1]->v,f->vert[2]->v,
            tri[0],tri[1],tri[2],
            &cop,&is1,&is2)) {
          continue;
        }
        if(cop) continue;
        /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop\n",k,ii,jj,kk,nn,n); */
        if((m=niik_image_check_index_niikpt(img,is2))>=0) {
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
        }
        dd=delta/niikpt_distance(is1,is2);
        for(dval=0; dval<=1.0; dval+=dd) {
          p=niikpt_wavg(is1,is2,dval);
          m=niik_image_check_index_niikpt(img,p);
          if(mm==m) continue;
          if(m<0) continue;
          /* fprintf(stdout,"depth = %3i; bb %3i %3i %3i; idx %6i; nf %5i; cop; voxel %8i\n",k,ii,jj,kk,nn,n,m); */
          dx[m]=pval*f->color[0];
          dy[m]=pval*f->color[1];
          dz[m]=pval*f->color[2];
          mm=m;
        } /* add color along intersect line */
      }
    } /* each face in bounding box */
  }
  /* put the values in the image */
  if(!niik_image_set_voxels_from_double_vector(img,idata)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    return 0;
  }
  free(idata);
  bb=off_bbox_free(bb);
  return 1;
}


#endif /* _FALCON_OFF_OBJ2IMG_C_ */
