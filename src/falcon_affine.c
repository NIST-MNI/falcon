/* Filename:     nifti1_kunio_affine.c
 * Description:  affine and resample functions
 * Author:       Kunio Nakamura
 * Date:         February 27, 2012
 */

#ifndef _FALCON_AFFINE_C_
#define _FALCON_AFFINE_C_

#include "falcon.h"


int niik_image_affine_transform_3d_update(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp)
/* niik_image_affine_transform_3d_update
 * -transform img to the space of refimg
 *  voxel size and image dimension are matched
 * -img is replaced with the output
 * -if interp is NIIK_INTERP_BSPLINE, then b-spline coefficients are calculated from img
 * --but img is unchanged
 *
 * see also niik_image_affine_transform_3d_update */
{
  nifti_image
  *tmpimg=NULL,
   *bsplimg=NULL,
    *outimg=NULL;
  niikmat
  *v=NULL,
   *invmat=NULL;
  niikpt pt,qt;
  int
  b,bb,
  size,
  i,j,k,n,
  verbose=0;
  char fcname[64]="niik_image_affine_transform_3d_update";
  if(   img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  /*if(   img->ndim>3) { fprintf(stderr,"ERROR: img has more than 3d\n"); return 0; }
    if(refimg->ndim>3) { fprintf(stderr,"ERROR: refimg has more than 3d\n"); return 0; }*/
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) created outimg\n");
  /* copy important data from reference
  if(!niik_image_copy_ref_info(refimg,img)){
    fprintf(stderr,"ERROR: niik_image_copy_ref_info\n");
    reutrn 0; }*/
  /*img->ndim = img->dim[0] = refimg->ndim;*/
  img->nx = img->dim[1] = refimg->nx;
  img->ny = img->dim[2] = refimg->ny;
  img->nz = img->dim[3] = refimg->nz;
  img->dx = img->pixdim[1] = refimg->dx;
  img->dy = img->pixdim[2] = refimg->dy;
  img->dz = img->pixdim[3] = refimg->dz;
  img->nvox = niik_image_calc_nvox(img);
  img->sform_code = refimg->sform_code;
  img->qform_code = refimg->qform_code;
  /*img-> freq_dim = refimg-> freq_dim;
    img->phase_dim = refimg->phase_dim;
    img->slice_dim = refimg->slice_dim;
    img->slice_start = refimg->slice_start;
    img->slice_end = refimg->slice_end;
    img->slice_duration = refimg->slice_duration;*/
  if(img->qform_code) {
    img->quatern_b = refimg->quatern_b;
    img->quatern_c = refimg->quatern_c;
    img->quatern_d = refimg->quatern_d;
    img->qoffset_x = refimg->qoffset_x;
    img->qoffset_x = refimg->qoffset_y;
    img->qoffset_z = refimg->qoffset_z;
    img->qfac      = refimg->qfac;
    img->qto_xyz   = refimg->qto_xyz;
    img->qto_ijk   = refimg->qto_ijk;
  }
  if(img->sform_code) {
    img->sto_xyz   = refimg->sto_xyz;
    img->sto_ijk   = refimg->sto_ijk;
  }
  img->toffset    = refimg->toffset;
  img->xyz_units  = refimg->xyz_units;
  img->time_units = refimg->time_units;
  free(img->data);
  if(verbose>=1) fprintf(stdout,"[%s] calloc %i %i\n",fcname,img->nvox,img->nbyper);
  if((img->data = (void *)calloc(img->nvox*img->nbyper,1))==NULL) {
    fprintf(stderr,"ERROR: malloc for image data\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) prepared output image\n");
  if(regmat==NULL) {
    if((invmat=niikmat_identity(4,4))==NULL) {
      fprintf(stderr,"ERROR: niikmat_identity\n");
      return 0;
    }
  } else {
    if((invmat=niikmat_inverse(regmat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse_matrix\n");
      return 0;
    }
  }
  if(verbose) {
    fprintf(stdout,"-d (niik_image_affine_transform_3d) inverse matrix\n");
    niikmat_display(invmat);
  }
  if(verbose)   fprintf(stdout,"-d (niik_image_affine_transform_3d) modify matrix for voxel spacing\n");
  if(interp==NIIK_INTERP_BSPLINE) {
    if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) bspline coefficients\n");
    if((bsplimg=niik_image_copy(outimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) bspline coefficient\n");
    if(!niik_image_interpolate_convert_3d_bspline_coeff(bsplimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n");
      return 0;
    }
    tmpimg=bsplimg;
  } else
    tmpimg=outimg;
  niikmat_multiply_mat2_free1(niikmat_scale_matrix(1.0/tmpimg->dx,1.0/tmpimg->dy,1.0/tmpimg->dz),invmat);
  if(verbose) {
    fprintf(stdout,"-d (niik_image_affine_transform_3d) include voxel spacing\n");
    niikmat_display(invmat);
  }
  if(!niik_image_clear(img)) {
    fprintf(stderr,"ERROR: niik_image_clear\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) go thru voxels %s\n",niik_interpolate_string(interp));
  size=img->nx*img->ny*img->nz;
  bb=img->nvox/size;
  bb=(bb<1)?1:bb;
  v=niikmat_init(img->nz,12);

  #pragma omp parallel for private (i,j,n,pt,qt,b)
  for(k=0; k<img->nz; k++) {
    pt.z = k*img->dz;
    n=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {
      pt.y = j*img->dy;
      for(i=0; i<img->nx; n++,i++) {
        pt.x = i*img->dx;
        qt = niikpt_affine_transform(invmat,pt);
        if(qt.x<0) continue;
        if(qt.y<0) continue;
        if(qt.z<0) continue;
        if(!niik_image_interpolate_3d_ijk_update(tmpimg,qt,interp,v->m[k])) {
          fprintf(stderr,"ERROR: niik_image_interpolate_3d_ijk_update\n");
          continue;
        }
        for(b=0; b<bb; b++) {
          niik_image_set_voxel(img,n+b*size,v->m[k][b]);
        }
        /*if(!niik_image_set_voxel(img,n,niik_image_interpolate_3d_ijk(tmpimg,qt,interp))){
          fprintf(stderr,"ERROR: niik_image_set_voxel\n");
          return 0; }*/
      }
    }
  } /* 3d part */
  bsplimg=niik_image_free(bsplimg);
  outimg=niik_image_free(outimg);
  invmat=niikmat_free(invmat);
  v=niikmat_free(v);
  return 1;
}


int niik_image_inverse_affine_transform_3d_update(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp) {
  niikmat *invmat=NULL;
  if(regmat==NULL) {
    fprintf(stderr,"ERROR: regmat is null\n");
    return 0;
  }
  invmat=niikmat_inverse(regmat);
  if(!niik_image_affine_transform_3d_update(img,refimg,invmat,interp)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  invmat=niikmat_free(invmat);
  return 1;
}


nifti_image *niik_image_affine_transform_3d(nifti_image *img,nifti_image *refimg,niikmat *regmat,int interp)
/* see also niik_image_affine_transform_3d_update */
{
  nifti_image
  *tmpimg=NULL,
   *bsplimg=NULL,
    *outimg=NULL;
  niikmat *invmat=NULL;
  niikpt pt,qt;
  double d;
  int
  i,j,k,n,
  verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return NULL;
  }
  if(   img->ndim>3) {
    fprintf(stderr,"ERROR: img has more than 3d\n");
    return NULL;
  }
  if(refimg->ndim>3) {
    fprintf(stderr,"ERROR: refimg has more than 3d\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) start\n");
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) created outimg\n");
  /* copy important data from reference */
  outimg->ndim = outimg->dim[0] = refimg->ndim;
  outimg->nx = outimg->dim[1] = refimg->nx;
  outimg->ny = outimg->dim[2] = refimg->ny;
  outimg->nz = outimg->dim[3] = refimg->nz;
  outimg->dx = outimg->pixdim[1] = refimg->dx;
  outimg->dy = outimg->pixdim[2] = refimg->dy;
  outimg->dz = outimg->pixdim[3] = refimg->dz;
  outimg->nvox = outimg->nx*outimg->ny*outimg->nz;
  outimg->sform_code = refimg->sform_code;
  outimg->qform_code = refimg->qform_code;
  /*outimg-> freq_dim = refimg-> freq_dim;
    outimg->phase_dim = refimg->phase_dim;
    outimg->slice_dim = refimg->slice_dim;
    outimg->slice_start = refimg->slice_start;
    outimg->slice_end = refimg->slice_end;
    outimg->slice_duration = refimg->slice_duration;*/
  if(outimg->qform_code) {
    outimg->quatern_b = refimg->quatern_b;
    outimg->quatern_c = refimg->quatern_c;
    outimg->quatern_d = refimg->quatern_d;
    outimg->qoffset_x = refimg->qoffset_x;
    outimg->qoffset_x = refimg->qoffset_y;
    outimg->qoffset_z = refimg->qoffset_z;
    outimg->qfac      = refimg->qfac;
    outimg->qto_xyz   = refimg->qto_xyz;
    outimg->qto_ijk   = refimg->qto_ijk;
  }
  if(outimg->sform_code) {
    outimg->sto_xyz   = refimg->sto_xyz;
    outimg->sto_ijk   = refimg->sto_ijk;
  }
  outimg->toffset    = refimg->toffset;
  outimg->xyz_units  = refimg->xyz_units;
  outimg->time_units = refimg->time_units;
  free(outimg->data);
  if((outimg->data = (void *)calloc(outimg->nvox*outimg->nbyper,1))==NULL) {
    fprintf(stderr,"ERROR: malloc for image data\n");
    return NULL;
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) prepared output image\n");
  if(regmat==NULL) {
    if((invmat=niikmat_identity(4,4))==NULL) {
      fprintf(stderr,"ERROR: niikmat_identity\n");
      return NULL;
    }
  } else {
    if((invmat=niikmat_inverse(regmat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse_matrix\n");
      return NULL;
    }
  }
  if(interp==NIIK_INTERP_BSPLINE) {
    if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) bspline coefficients\n");
    if((bsplimg=niik_image_copy(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return NULL;
    }
    if(!niik_image_interpolate_convert_3d_bspline_coeff(bsplimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n");
      return NULL;
    }
    tmpimg=bsplimg;
  } else
    tmpimg=img;
  if(verbose) {
    fprintf(stdout,"-d (niik_image_affine_transform_3d) modify matrix for voxel spacing\n");
  }
  niikmat_multiply_mat2_free1(niikmat_scale_matrix(1.0/tmpimg->dx,1.0/tmpimg->dy,1.0/tmpimg->dz),invmat);
  if(verbose) {
    fprintf(stdout,"-d (niik_image_affine_transform_3d) include voxel spacing\n");
    niikmat_display(invmat);
  }
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) ref max val %8.3f \n",niik_image_get_max(refimg,NULL));
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) img max val %8.3f \n",niik_image_get_max(img,NULL));
  if(verbose) fprintf(stdout,"-d (niik_image_affine_transform_3d) go thru voxels %s\n",niik_interpolate_string(interp));
  #pragma omp parallel for private (i,j,n,pt,qt,d)
  for(k=0; k<outimg->nz; k++) {
    pt.z = k*outimg->dz;
    n=k*outimg->nx*outimg->ny;
    for(j=0; j<outimg->ny; j++) {
      pt.y = j*outimg->dy;
      for(i=0; i<outimg->nx; n++,i++) {
        pt.x = i*outimg->dx;
        if(!niik_image_set_voxel(outimg,n,0)) {
          fprintf(stderr,"ERROR: niik_image_set_voxel\n");
          continue;
        }
        qt = niikpt_affine_transform(invmat,pt);
        if(qt.x<0) continue;
        if(qt.y<0) continue;
        if(qt.z<0) continue;
        d=niik_image_interpolate_3d_ijk(tmpimg,qt,interp);
        if(!niik_image_set_voxel(outimg,n,d)) {
          fprintf(stderr,"ERROR: niik_image_set_voxel\n");
          continue;
        }
      }
    }
  }
  nifti_image_free(bsplimg);
  return outimg;
}

nifti_image *niik_image_resample_3d(nifti_image *img,double dx, double dy, double dz,int nx, int ny,int nz,int interp) {
  nifti_image *outimg=NULL;
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_resample_3d_update(outimg,dx, dy, dz,nx, ny,nz,interp)) {
    fprintf(stderr,"ERROR: niik_image_resample_3d_update(img,dx, dy, dz,nx, ny,nz,interp)\n");
    outimg=niik_image_free(outimg);
    return NULL;
  }
  return outimg;
}

int niik_image_resample_3d_update(nifti_image *img,double dx, double dy, double dz,int nx, int ny,int nz,int interp)
/* -img is updated with resample with dx,dy,dz and nx,ny,nz
 * -dx,dy,dz will be the new voxel delta
 * -nx,ny,nz will be the new image dimension (if negative, automatic)
 */
{
  nifti_image
  *tmpimg=NULL,
   *bsplimg=NULL,
    *outimg=NULL;
  niikpt pt;
  double
  **vlist=NULL;
  niikmat *afmat=NULL;
  int
  b=0,nb=0,
  area=0,size=0,
  i=0,j=0,k=0,n=0,
  verbose=2;

  const char *fcname="niik_image_resample_3d_update";
  pt=niikpt_zero();
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  /*if(   img->ndim>3) { fprintf(stderr,"ERROR: img has more than 3d\n"); return 0; }*/
  if(img->ndim<3) {
    fprintf(stderr,"[%s] ERROR: img needs at least 3 dimensions\n",fcname);
    return 0;
  }
  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;

  }
  if(verbose>=1) fprintf(stdout,"[%s]   created outimg\n",fcname);
  /* copy important data from reference */
  /* img->ndim = img->dim[0] = 3; */
  img->dx = img->pixdim[1] = dx;
  img->dy = img->pixdim[2] = dy;
  img->dz = img->pixdim[3] = dz;

  if(nx<0) img->nx = img->dim[1] = (int)floor(outimg->nx*outimg->dx/dx+0.5);
  else     img->nx = img->dim[1] = nx;
  if(ny<0) img->ny = img->dim[2] = (int)floor(outimg->ny*outimg->dy/dy+0.5);
  else     img->ny = img->dim[2] = ny;
  if(nz<0) img->nz = img->dim[3] = (int)floor(outimg->nz*outimg->dz/dz+0.5);
  else     img->nz = img->dim[3] = nz;

  img->nvox = niik_image_calc_nvox(img);
  if(verbose) {
    fprintf(stdout,"\t  image dim = %7i %7i %7i   <-   %7i %7i %7i\n",img->nx,img->ny,img->nz,outimg->nx,outimg->ny,outimg->nz);
    fprintf(stdout,"\t  image vox = %7.4f %7.4f %7.4f   <-   %7.4f %7.4f %7.4f\n",img->dx,img->dy,img->dz,outimg->dx,outimg->dy,outimg->dz);
  }
  free(img->data);
  if(verbose>=1) fprintf(stdout,"[%s]   calloc %i %i\n",fcname,img->nvox,img->nbyper);
  if((img->data = (void *)calloc(img->nvox*img->nbyper,1))==NULL) {
    fprintf(stderr,"ERROR: malloc for image data\n");
    return 0;
  }
  if(verbose>1) fprintf(stdout,"-d (niik_image_resample_3d) prepared output image\n");
  if(verbose) {
    fprintf(stdout,"\t  resampling with %s\n",niik_interpolate_string(interp));
  }
  if(interp==NIIK_INTERP_BSPLINE) {
    if(verbose>1) fprintf(stdout,"-d (niik_image_resample_3d) bspline coefficients\n");
    if((bsplimg=niik_image_copy(outimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(!niik_image_interpolate_convert_3d_bspline_coeff(bsplimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n");
      return 0;
    }
    tmpimg=bsplimg;
  } else {
    tmpimg=outimg;
  } /* BSPLINE OR ELSE */
  if(verbose>1) fprintf(stdout,"[%s] img max val %8.3f\n",fcname,niik_image_get_max(img,NULL));
  if(verbose>1) fprintf(stdout,"[%s] go thru voxels %s\n",fcname,niik_interpolate_string(interp));
  area=img->nx*img->ny;
  for(n=4,nb=1; n<=img->ndim; n++) {
    nb *= img->dim[n];
  }

  size= img->nx * img->ny * img->nz;
  NIIK_RET0((nb<1),fcname,"wrong nb");
  if(verbose>=1) fprintf(stdout,"[%s] area = %i, size = %i, nb = %i\n",fcname,area,size,nb);
  if(verbose>=1) fprintf(stdout,"[%s]   calloc %i %i\n",fcname,img->nz,(int)sizeof(double *));
  if((vlist = (double **)calloc(img->nz,sizeof(double *)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc for vlist, (%i,%i)\n",fcname,img->nz,(int)sizeof(double *));
    return 0;
  }
  for(k=0; k<img->nz; k++) {
    if((vlist[k] = (double *)calloc(nb,sizeof(double)))==NULL) {
      fprintf(stderr,"[%s] ERROR: calloc (%i,%i)  for k=%i\n",fcname,nb,(int)sizeof(double),k);
      return 0;
    }
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] clear img\n",fcname);
  }
  NIIK_RET0((!niik_image_clear(img)),fcname,"niik_image_clear failed");
  /* 2012-04-07 Kunio
   * -fixed memory issues with parallel processing here */
  img->qform_code = NIFTI_XFORM_UNKNOWN;
  afmat = niikmat_mat44_matrix(img->sto_xyz);
  niikmat_multiply_mat1_free2(afmat,niikmat_scale_matrix(img->dx/outimg->dx,
                              img->dy/outimg->dy,
                              img->dz/outimg->dz));
  img->sto_xyz=niikmat_make_mat44(afmat);
  niikmat_inverse_update(afmat);
  img->sto_ijk=niikmat_make_mat44(afmat);

  #pragma omp parallel for private (i,j,pt,n,b)
  for(k=0; k<img->nz; k++) {
    pt.w = 0;
    /*pt.z = k*img->dz;*/
    n=k*area;
    for(j=0; j<img->ny; j++) {
      /*pt.y = j*img->dy;*/
      for(i=0; i<img->nx; n++,i++) {
        /*pt.x = i*img->dx;*/

        /*if(pt.x<0) continue;
        if(pt.y<0) continue;
        if(pt.z<0) continue;*/

        niik_index_to_world(img,i,j,k,&pt);

        if(!niik_image_interpolate_3d_xyz_update(tmpimg,pt,interp,vlist[k])) {
          fprintf(stderr,"ERROR: niik_image_interpolate_3d_xyz_update(tmpimg,pt,interp,vlist)\n");
          continue;
        }
        /*d=niik_image_interpolate_3d_ijk(tmpimg,pt,interp);*/
        for(b=0; b<nb; b++) {
          niik_image_set_voxel(img,n+b*size,vlist[k][b]);
        }
      }
    }  /* i,j */
  } /* k */
  free(vlist);
  /* img->sform_code = NIFTI_XFORM_UNKNOWN; */
  bsplimg=niik_image_free(bsplimg);
  outimg=niik_image_free(outimg);
  return 1;
} /* niik_image_resample_3d_update */


int niik_image_copy_dim_info(nifti_image *src,nifti_image *dest)
/* copy dim info from src to dest including all nx,ny,nz... dim[1-7], and nvox */
{
  if(src ==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  if(dest==NULL) {
    fprintf(stderr,"ERROR: dest is null\n");
    return 0;
  }
  if(src->dim[0]>0) dest->nx = dest->dim[1] = src->nx;
  if(src->dim[0]>1) dest->ny = dest->dim[2] = src->ny;
  if(src->dim[0]>2) dest->nz = dest->dim[3] = src->nz;
  if(src->dim[0]>3) dest->nt = dest->dim[4] = src->nt;
  if(src->dim[0]>4) dest->nu = dest->dim[5] = src->nu;
  if(src->dim[0]>5) dest->nv = dest->dim[6] = src->nv;
  if(src->dim[0]>6) dest->nw = dest->dim[7] = src->nw;
  dest->nvox = src->nvox;
  return 1;
}


int niik_image_copy_pixdim_info(nifti_image *src,nifti_image *dest)
/* copy dim info from src to dest including all dx,dy,dz... pixdim[1-7] */
{
  if(src ==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  if(dest==NULL) {
    fprintf(stderr,"ERROR: dest is null\n");
    return 0;
  }
  if(src->dim[0]>0) dest->dx = dest->pixdim[1] = src->dx;
  if(src->dim[0]>1) dest->dy = dest->pixdim[2] = src->dy;
  if(src->dim[0]>2) dest->dz = dest->pixdim[3] = src->dz;
  if(src->dim[0]>3) dest->dt = dest->pixdim[4] = src->dt;
  if(src->dim[0]>4) dest->du = dest->pixdim[5] = src->du;
  if(src->dim[0]>5) dest->dv = dest->pixdim[6] = src->dv;
  if(src->dim[0]>6) dest->dw = dest->pixdim[7] = src->dw;
  return 1;
}


#endif /* _FALCON_AFFINE_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/