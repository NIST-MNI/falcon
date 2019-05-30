/* Filename:     nifti1_kunio_nregister.c
 * Description:  nonlinear registration function
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 *
 * Important functions
 *
 *   char *niik_warp_type_string( int code );
 *   int niik_image_apply_3d_warp_update(nifti_image *img,nifti_image *out_img,nifti_image *warp_img,int warp_type,int interp);
 *
 *
 */

#ifndef _FALCON_NREGISTER_C_
#define _FALCON_NREGISTER_C_

#include "falcon.h"



char *niik_warp_type_string( int code ) {
  switch(code) {
  case NIIK_WARP_MAP_DISP:
    return "displacement map";
  case NIIK_WARP_MAP_LOC:
    return "location map";
  case NIIK_WARP_MAP_DISP_BSPLINE:
    return "b-spline displacement map";
  default:
    break;
  }
  return "unknown";
}


/***************************************************
 *
 * APPLICATION OF WARP
 *
 ***************************************************/


int niik_image_apply_3d_warp_update(nifti_image *img,nifti_image *refimg,nifti_image *warp_img,int warp_type,int interp)
/* applies the warp
 * -img is the pre-warped image, and warped image is REPLACED HERE
 *   -if img has more than 3d, then output has the same number of dimensions
 *    but the x,y,z dimensions/voxel size come from refimg
 * -refimg serves as the template/reference image
 * -warp_img is the warp field (displacement or location map)
 * -warp_type is the type warp: NIIK_WARP_MAP_*
 * -interp is the image interpolation type (not the warp)
 */
{
  char fcname[64]="niik_image_apply_3d_warp_update";
  nifti_image
  *tmpimg=NULL; /* for b-spline coefficient calculation */
  niikpt pt,qt;
  niikmat *bt,*ct;
  int
  warp_interp=NIIK_INTERP_LINEAR,
  verbose=1,
  nimg,size,
  ni,mi,n,i,j,k;
  if(verbose) {
    fprintf(stdout,"[niik_image_apply_3d_warp] start\n");
    fprintf(stdout,"[niik_image_apply_3d_warp]   warp type         %s\n",niik_warp_type_string(warp_type));
    fprintf(stdout,"[niik_image_apply_3d_warp]   interp type       %s\n",niik_interpolate_string(interp));
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  if(warp_img==NULL) {
    fprintf(stderr,"ERROR: warp_img is null\n");
    return 0;
  }
  /* assign tmpimg with bspline or regular image */
  if(interp==NIIK_INTERP_BSPLINE) {
    if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp] bspline coeff calculation\n");
    if((tmpimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL) {
      fprintf(stderr,"ERROR: niik_image_type_convert(tmpimg)\n");
      return 0;
    }
    if(!niik_image_interpolate_convert_3d_bspline_coeff(tmpimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_convert_3d_bspline_coeff(tmpimg)\n");
      return 0;
    }
  } else {
    if((tmpimg=niik_image_copy(img))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
  }
  for(n=4,ni=1; n<img->ndim; n++) {
    ni *= img->dim[n];
  }
  nimg=ni;
  if(!niik_image_affine_transform_3d_update(img,refimg,NULL,NIIK_INTERP_NN)) {
    fprintf(stderr,"ERROR: niik_image_affine_transform_3d_update\n");
    return 0;
  }
  size=img->nx * img->ny * img->nz;
  for(n=4,mi=1; n<img->ndim; n++) {
    mi *= img->dim[n];
  }
  if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp]   #img %i\n",nimg);
  if(nimg!=mi) {
    /* update the image size beyond dimension 3 */
    if(verbose) fprintf(stdout,"[niik_image_apply_3d_warp] changing img's size\n");
    if(!niik_image_copy_ref_info(tmpimg,img)) {
      fprintf(stderr,"ERROR: niik_image_copy_ref_info\n");
      return 0;
    }
    free(img->data);
    if((img->data=(void *)calloc(img->nvox,img->nbyper))==NULL) {
      fprintf(stderr,"ERROR: caloc \n");
      return 0;
    }
  }
  /* for parallel processing */
  bt = niikmat_init(img->nz,3);
  ct = niikmat_init(img->nz,ni);
  if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp]   clear out image\n");
  niik_image_clear(img);
  /* actual image interpolation */
  switch(warp_type) {
  case NIIK_WARP_MAP_DISP_BSPLINE:
    if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp] NIIK_WARP_MAP_DISP_BSPLINE\n");
    warp_interp = NIIK_INTERP_BSPLINE;
  case NIIK_WARP_MAP_DISP:
    if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp] NIIK_WARP_MAP_DISP\n");
    #pragma omp parallel for private(i,j,n,pt,qt,ni)
    for(k=0; k<img->nz; k++) {
      pt.w = 0;
      pt.z = k*img->dz;
      n = k * img->nx * img->ny;
      for(j=0; j<img->ny; j++) {
        pt.y = j*img->dy;
        for(i=0; i<img->nx; n++,i++) {
          pt.x = i*img->dx;
          if(!niik_image_interpolate_3d_xyz_update(warp_img,pt,warp_interp,bt->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
            continue;
          }
          /* warp displacement */
          qt = niikpt_val(pt.x+bt->m[k][0],pt.y+bt->m[k][1],pt.z+bt->m[k][2],0);
          if(!niik_image_interpolate_3d_xyz_update(tmpimg,qt,interp,ct->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
            continue;
          }
          for(ni=0; ni<nimg; ni++) {
            if(!niik_image_set_voxel(img,n+ni*size,ct->m[k][ni])) {
              fprintf(stderr,"ERROR: niik_image_set_voxel\n");
              continue;
            }
          }
        }
      }
    }
    break;
  case NIIK_WARP_MAP_LOC:
    if(verbose>1) fprintf(stdout,"[niik_image_apply_3d_warp] NIIK_WARP_MAP_LOC\n");
    #pragma omp parallel for private(i,j,n,pt,qt,ni)
    for(k=0; k<img->nz; k++) {
      pt.z = k*img->dz;
      n = k * img->nx * img->ny;
      for(j=0; j<img->ny; j++) {
        pt.y = j*img->dy;
        for(i=0; i<img->nx; n++,i++) {
          pt.x = i*img->dx;
          if(!niik_image_interpolate_3d_xyz_update(warp_img,pt,warp_interp,bt->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
            continue;
          }
          /* warp displacement */
          qt = niikpt_val(bt->m[k][0],bt->m[k][1],bt->m[k][2],0);
          if(!niik_image_interpolate_3d_xyz_update(tmpimg,qt,interp,ct->m[k])) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
            continue;
          }
          for(ni=0; ni<nimg; ni++) {
            if(!niik_image_set_voxel(img,n+ni*size,ct->m[k][ni])) {
              fprintf(stderr,"ERROR: niik_image_set_voxel\n");
              continue;
            }
          }
        }
      }
    }
    break;
  default:
    fprintf(stderr,"ERROR: Unknown warp_type %i\n",warp_type);
    return 0;
  }
  tmpimg = niik_image_free(tmpimg);
  bt=niikmat_free(bt);
  ct=niikmat_free(ct);
  if(verbose) fprintf(stdout,"[niik_image_apply_3d_warp] successful\n");
  return 1;
} /* niik_image_apply_3d_warp_update */




nifti_image *niik_image_combine_warp(nifti_image *movimg,nifti_image *refimg,nifti_image *warpimg,
                                     niikmat *premat, niikmat *postmat,int warp_map_type)
/* -the order is
 *    movimg -> premat -> (fwd)warp -> postmat -> refimg
 * -therefore, for each voxel in refimg, we need coordinates from movimg
 * -The images should be registered by
 *  1. movimg -> premat affine registration
 *  2. Premat(movimg) -> warp with warpimg
 *  3. Warp(Premat(movimg)) -> registered image
 *
 * 2012-06-01, Kunio
 * -I'm not sure if it is working "accurately"
 * -It is working for now... so don't touch unless I know what I'm doing...
 */
{
  nifti_image
  *outimg=NULL,
   *tmpimg=NULL;
  niikmat *postmatinv=NULL,*prematinv=NULL;
  niikpt
  p,q,r,w;
  float *fx,*fy,*fz;
  double v[4];
  int
  warp_interp=NIIK_INTERP_LINEAR,
  dimadd=0,
  size3d,
  i,j,k,n,
  verbose=1;

  if(verbose>1) fprintf(stdout,"[niik_image_combine_warp] start\n");
  if(movimg==NULL)  {
    fprintf(stderr,"ERROR: movimg is null\n");
    return NULL;
  }
  if(refimg==NULL)  {
    fprintf(stderr,"ERROR: refimg is null\n");
    return NULL;
  }
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg is null\n");
    return NULL;
  }
  /*if(premat==NULL) { fprintf(stderr,"ERROR: premat is null\n"); return NULL; }
    if(postmat==NULL) { fprintf(stderr,"ERROR: postmat is null\n"); return NULL; }*/
  if(verbose>1) fprintf(stdout,"[niik_image_combine_warp] prepare direction image\n");

  p.w=q.w=r.w=w.w=0;

  /* output positional image */
  if((outimg=niik_image_copy(movimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  free(outimg->data);
  outimg->ndim=outimg->dim[0]=5;
  outimg->nx=outimg->dim[1]=outimg->nx+dimadd*2;
  outimg->ny=outimg->dim[2]=outimg->ny+dimadd*2;
  outimg->nz=outimg->dim[3]=outimg->nz+dimadd*2;
  outimg->nt=outimg->dim[4]=1;
  outimg->nu=outimg->dim[5]=3;
  size3d=outimg->nx*outimg->ny*outimg->nz;
  outimg->nvox=outimg->nx*outimg->ny*outimg->nz*outimg->nu;
  if((outimg->data=(void *)calloc(outimg->nvox,outimg->nbyper))==NULL) {
    fprintf(stderr,"ERROR: calloc for outimg->data\n");
    return NULL;
  }
  fx=(float *)outimg->data;
  fy=fx+size3d;
  fz=fy+size3d;
  if(premat!=NULL) {
    if((prematinv = niikmat_inverse(premat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse\n");
      return NULL;
    }
  }
  if(postmat!=NULL) {
    if((postmatinv = niikmat_inverse(postmat))==NULL) {
      fprintf(stderr,"ERROR: niikmat_inverse\n");
      return NULL;
    }
  }
  v[0]=v[1]=v[2]=0;

  /* warp */
  if(verbose) fprintf(stdout,"[niik_image_combine_warp] positional image\n");
  for(k=n=0; k<outimg->nz; k++) {
    p.z = outimg->dz*(k-dimadd);
    for(j=0; j<outimg->ny; j++) {
      p.y = outimg->dy*(j-dimadd);
      for(i=0; i<outimg->nx; n++,i++) {
        p.x = outimg->dx*(i-dimadd);
        fx[n]=p.x;
        fy[n]=p.y;
        fz[n]=p.z;
      }
    }
  }
  if(0) {
    fprintf(stdout,"writing output image\n");
    niik_image_write("tmp_dirimg-1.nii.gz",outimg);
  }

  /* final warp image */
  if((tmpimg=niik_image_copy(refimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert\n");
    return NULL;
  }
  free(tmpimg->data);
  tmpimg->ndim=tmpimg->dim[0]=5;
  tmpimg->nx=tmpimg->dim[1]=tmpimg->nx;
  tmpimg->ny=tmpimg->dim[2]=tmpimg->ny;
  tmpimg->nz=tmpimg->dim[3]=tmpimg->nz;
  tmpimg->nt=tmpimg->dim[4]=1;
  tmpimg->nu=tmpimg->dim[5]=3;
  size3d=tmpimg->nx*tmpimg->ny*tmpimg->nz;
  tmpimg->nvox=tmpimg->nx*tmpimg->ny*tmpimg->nz*tmpimg->nu;
  if((tmpimg->data=(void *)calloc(tmpimg->nvox,tmpimg->nbyper))==NULL) {
    fprintf(stderr,"ERROR: calloc for tmpimg->data\n");
    return NULL;
  }

  /* warp */
  if(verbose) fprintf(stdout,"[niik_image_combine_warp] nonlinear warp\n");

  switch(warp_map_type) {
  case NIIK_WARP_MAP_LOC: /* position map */
    if(verbose>1) fprintf(stdout,"[niik_image_combine_warp]   warp using location map\n");
    fprintf(stderr,"ERROR: not implemented yet\n");
    return NULL;
    break;

  case NIIK_WARP_MAP_DISP_BSPLINE: /* bspline displacement map */
    warp_interp=NIIK_INTERP_BSPLINE;
    if(verbose>1) fprintf(stdout,"[niik_image_combine_warp]   warp using b-spline displacement map\n");

  case NIIK_WARP_MAP_DISP: /* displacement map */
    if(verbose>1) fprintf(stdout,"[niik_image_combine_warp]   warp using displacement map\n");
    n=tmpimg->nx*tmpimg->ny*tmpimg->nz;
    fx=(float *)tmpimg->data;
    fy=fx+n;
    fz=fy+n;
    for(k=n=0; k<tmpimg->nz; k++) {
      p.z = tmpimg->dz*(k-dimadd);
      for(j=0; j<tmpimg->ny; j++) {
        p.y = tmpimg->dy*(j-dimadd);
        for(i=0; i<tmpimg->nx; n++,i++) {
          p.x = tmpimg->dx*(i-dimadd);
          if(postmatinv!=NULL) {
            q = niikpt_affine_transform(postmatinv,p);
          } else
            q = p;
          if(verbose>2)fprintf(stdout,"[%3i %3i %3i] %9.3f %9.3f %9.3f   %9.3f %9.3f %9.3f\n",i,j,k,p.x,p.y,p.z,q.x,q.y,q.z);
          if(!niik_image_interpolate_3d_xyz_update(warpimg,q,warp_interp,v)) {
            fprintf(stderr,"ERROR: niik_image_interpolate_3d_xyz_update\n");
            return NULL;
          }
          r.x=q.x+v[0];
          r.y=q.y+v[1];
          r.z=q.z+v[2];
          if(prematinv!=NULL)
            w = niikpt_affine_transform(prematinv,r);
          else
            w = r;
          fx[n]=w.x;
          fy[n]=w.y;
          fz[n]=w.z;
        }
      }
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown warp map type %i\n",warp_map_type);
    return NULL;
  }

  tmpimg->cal_min = niik_image_get_min(tmpimg,NULL);
  tmpimg->cal_max = niik_image_get_max(tmpimg,NULL);
  prematinv=niikmat_free(prematinv);
  postmatinv=niikmat_free(postmatinv);

  /* fprintf(stdout,"writing output image\n");
     niik_image_write("tmp_dirimg-2.nii.gz",tmpimg);
     niik_image_apply_3d_warp_update(movimg,refimg,tmpimg,NIIK_WARP_MAP_LOC,NIIK_INTERP_LINEAR);
     niik_image_write("tmp_dirimg-3.nii.gz",movimg); */

  return tmpimg;
}




/************************************************************************
 *
 * INVERSION OF NONLINEAR WARP
 *
 * -not exact, approximateion
 *
 ************************************************************************/


int niik_image_nregister_invert_nonlinear_map(nifti_image *warpimg,int warp_type)
/* invert the warp map (generic version)
 * -warpimg is a warp map of type warp_type
 * -warpimg is replaced with the inverted warp
 */
{
  nifti_image *tmpimg;
  int verbose=1;
  if(verbose) fprintf(stdout,"[niik_image_nregister_invert_nonlinear_map] start\n");
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg is null\n");
    return 0;
  }
  switch(warp_type) {
  case NIIK_WARP_MAP_DISP_BSPLINE:
    if(verbose) fprintf(stdout,"[niik_image_nregister_invert_nonlinear_map]   bspline displacement map\n");
    if((tmpimg=niik_image_copy(warpimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_copy\n");
      return 0;
    }
    if(!niik_image_interpolate_inverse_3d_bspline_coeff(tmpimg,warpimg)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_inverse_3d_bspline_coeff\n");
      return 0;
    }
    tmpimg=niik_image_free(tmpimg);
  case NIIK_WARP_MAP_DISP:
    if(verbose) fprintf(stdout,"[niik_image_nregister_invert_nonlinear_map]   displacement map\n");
    if(!niik_image_nregister_invert_displacement_map(warpimg)) {
      fprintf(stderr,"ERROR: niik_image_nregister_invert_displacement_map\n");
      return 0;
    }
    break;
  case NIIK_WARP_MAP_LOC:
    if(verbose) fprintf(stdout,"[niik_image_nregister_invert_nonlinear_map]   location map\n");
    fprintf(stderr,"ERROR: not implemented yet\n");
    return 0;
  default:
    fprintf(stderr,"ERROR: unknown warp type %i\n",warp_type);
    return 0;
  }
  return 1;
} /* niik_image_nregister_invert_nonlinear_map */


int niik_image_nregister_invert_displacement_map_iterative(nifti_image *warpimg,int maxiter)
/* Simple estimation of inverse
 * Reference :
 * Mingli Chen, et al. 2008. A simple fixed-point approach to invert a deformation field. Medical Physics. 35(1), 81-88.
 * -10 iterations is enough.
 */
{
  nifti_image
  *tmpimg=NULL;
  char fcname[64]="niik_image_nregister_invert_displacement_map_iterative";
  int i,dim3,iter;
  niikpt p;
  double v[9];
  float *fi[3];
  niik_fc_display(fcname,1);
  NIIK_RET0((!niik_image_type_convert(warpimg,NIFTI_TYPE_FLOAT32)),fcname,"niik_image_type_convert");
  NIIK_RET0(((tmpimg = niik_image_copy_as_type(warpimg,NIFTI_TYPE_FLOAT32))==NULL),fcname,"niik_image_copy_as_type");
  NIIK_RET0((!niik_image_clear(warpimg)),fcname,"niik_image_clear");
  dim3=warpimg->nvox/3;
  fi[0]=(float *)warpimg->data;
  fi[1]=fi[0]+dim3;
  fi[2]=fi[1]+dim3;
  p.w=0;
  for(iter=0; iter<maxiter; iter++) {
    for(i=0; i<dim3; i++) {
      p.x=warpimg->dx*(i%warpimg->nx);
      p.y=warpimg->dy*((i/warpimg->nx)%warpimg->ny);
      p.z=warpimg->dz*((i/warpimg->nx)/warpimg->ny);
      niik_image_interpolate_3d_xyz_update(warpimg,p,NIIK_INTERP_LINEAR,v);
      p.x+=v[0];
      p.y+=v[1];
      p.z+=v[2];
      niik_image_interpolate_3d_xyz_update(tmpimg,p,NIIK_INTERP_LINEAR,v);
      fi[0][i]=-v[0];
      fi[1][i]=-v[1];
      fi[2][i]=-v[2];
    }
  }
  tmpimg=niik_image_free(tmpimg);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_nregister_invert_displacement_map_iterative */


int niik_image_nregister_invert_displacement_map(nifti_image *warpimg)
/* invert the warp displacement map
 * -warpimg is a displacement map
 * -warpimg is replaced with the inverted warp
 * -warp should be close to zero near the image edge
 *
 * -for each voxel,
 *    calculate the inverse warp position
 *    put the inverse warp at the inverse warp position
 *    dilate the warp and its weighting (euclidian distance)
 * -divide by the weighting
 */
{
  nifti_image
  *maskimg=NULL,
   *tmpimg=NULL;
  int
  verbose=2,
  /*warp_interp=NIIK_INTERP_LINEAR,*/
  size3,
  ii,jj,kk,nn,
  ilo,ihi,jlo,jhi,klo,khi,
  i,j,k,m,n;
  double
  plim=1e-2,
  *fx,*fy,*fz,*f[3],
  *wx,*wy,*wz,*w[3];
  niikpt p,q,r;
  niikmat *bt;
  float *fimg;
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] start\n");
  if(warpimg==NULL) {
    fprintf(stderr,"ERROR: warpimg si null\n");
    return 0;
  }
  if(warpimg->nu!=3) {
    fprintf(stderr,"ERROR: warpimg needs nu=3 (%i) \n",warpimg->nu);
    return 0;
  }
  if(warpimg->datatype!=NIFTI_TYPE_FLOAT64) {
    fprintf(stderr,"ERROR: warpimg is not double float (%s)\n",nifti_datatype_string(warpimg->datatype));
    return 0;
  }
  if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] copy image\n");
  if((tmpimg=niik_image_copy_as_type(warpimg,NIFTI_TYPE_FLOAT64))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return 0;
  }
  niik_image_clear(warpimg);
  if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] maskimg\n");
  if((maskimg=niik_image_copy_as_type(warpimg,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
    return 0;
  }
  niik_image_clear(maskimg);
  fimg=maskimg->data;
  bt = niikmat_init(warpimg->nz,3);
  size3 = warpimg->nx * warpimg->ny * warpimg->nz;
  w[0] = wx = warpimg->data;
  w[1] = wy = wx + size3;
  w[2] = wz = wy + size3;
  f[0] = fx = tmpimg->data;
  f[1] = fy = fx + size3;
  f[2] = fz = fy + size3;
  /* initial guess */
  if(verbose) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] main loop\n");
  #pragma omp parallel for private(i,j,n,p,q,ii,jj,kk,ilo,ihi,jlo,jhi,klo,khi,r,nn,m)
  for(k=0; k<warpimg->nz; k++) {
    p.z=k*warpimg->dz;
    n = k*warpimg->nx*warpimg->ny;
    for(j=0; j<warpimg->ny; j++) {
      p.y=j*warpimg->dy;
      for(i=0; i<warpimg->nx; n++,i++) {
        p.x=i*warpimg->dx;
        /*if(!niik_image_interpolate_3d_xyz_update(tmpimg,p,warp_interp,bt->m[k])){
          fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
          continue; } */
        for(m=0; m<3; m++) {
          bt->m[k][m]=f[m][n];
        }
        if(fabs(bt->m[k][0])<plim) {
          if(fabs(bt->m[k][1])<plim) {
            if(fabs(bt->m[k][2])<plim) {
              continue;
            }
          }
        }
        q=niikpt_val(p.x+bt->m[k][0],p.y+bt->m[k][1],p.z+bt->m[k][2],0);
        /*fprintf(stdout,"[%3i %3i %3i] %7.2f %7.2f %7.2f -> %7.2f %7.2f %7.2f\n",i,j,k,
          p.x,p.y,p.z,q.x,q.y,q.z); */
        ilo=NIIK_IMAX(q.x-2,0);
        ihi=NIIK_IMIN(q.x+2,warpimg->nx-1);
        jlo=NIIK_IMAX(q.y-2,0);
        jhi=NIIK_IMIN(q.y+2,warpimg->ny-1);
        klo=NIIK_IMAX(q.z-2,0);
        khi=NIIK_IMIN(q.z+2,warpimg->nz-1);
        /*fprintf(stdout,"[%3i %3i %3i] %7.2f %7.2f %7.2f -> %7.2f %7.2f %7.2f  %3i %3i %3i\n",i,j,k,
          p.x,p.y,p.z,q.x,q.y,q.z,ii,jj,kk); */
        for(kk=klo; kk<=khi; kk++) {
          r.z = NIIK_SQ(kk*warpimg->dz-q.z);
          for(jj=jlo; jj<=jhi; jj++) {
            r.y = NIIK_SQ(jj*warpimg->dy-q.y);
            for(ii=ilo; ii<=ihi; ii++) {
              nn = ii + jj * warpimg->nx + kk * warpimg->nx * warpimg->ny;
              r.x = NIIK_SQ(ii*warpimg->dx-q.x);
              r.w = NIIK_GaussPDF(sqrt(r.x+r.y+r.z),1.5);
              /*fprintf(stdout,"[%3i %3i %3i] %7.2f %7.2f %7.2f -> %7.2f %7.2f %7.2f  [%3i %3i %3i]  %7.4f %7.4f %7.4f | %7.4f\n",i,j,k,
              p.x,p.y,p.z,q.x,q.y,q.z,ii,jj,kk,r.x,r.y,r.z,r.w); */
              fimg[nn] += r.w;
              for(m=0; m<3; m++) {
                w[m][nn] += -bt->m[k][m] * r.w;
              }
            }
          }
        }
      }
    }
  } /* each voxel */
  for(i=0; i<size3; i++) {
    if(fabs(fimg[i])<plim) continue;
    for(m=0; m<3; m++) {
      w[m][i] /= fimg[i];
    }
  }
  if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] free memory\n");
  tmpimg=niik_image_free(tmpimg);
  maskimg=niik_image_free(maskimg);
  bt=niikmat_free(bt);
  if(verbose>1) fprintf(stdout,"[niik_image_nregister_bspline_invert_displacement_map] successful finish\n");
  return 1;
} /* int niik_image_nregister_invert_displacement_map */

#endif /* _FALCON_NREGISTER_C_ */
