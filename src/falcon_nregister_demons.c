/* Filename:     nifti1_kunio_nregister_demons.c
 * Description:  demons-based nonlinear registration functions
 * Author:       Kunio Nakamura
 * Date:         April 9, 2012
 */

#ifndef _FALCON_NREGISTER_DEMONS_C_
#define _FALCON_NREGISTER_DEMONS_C_

#include "falcon.h"

int niik_image_nregister_demons_update(nifti_image *warpimg,nifti_image *gradimg,
                                       nifti_image *refimg,nifti_image *maskimg,
                                       nifti_image *tgtimg,niikmat *tgtmat,
                                       nifti_image *movimg,niikmat *movmat)
/* optical flow registration function
 * -updates warpfield (warpimg)
 * -refimg and maskimg are in the same space as warpimg
 * -tgtimg is transformed to refimg via tgmat 4x4 affine matrix
 * -movimg is transformed to refimg via movmat 4x4 affine matrix
 * -movimg is the moving image and tgtimg is the target image
 */

{
  char
  fcname[64]="niik_image_nregister_demons_update";
  int
  i,j,k,n,
  xdim,ydim,zdim,size,
  verbose=0;
  niikpt
  pt,qt,p1,p2;
  unsigned char
  *bimg=NULL;
  niikmat
  *ww=NULL,*ws=NULL,*gv=NULL,
   *tgtinv=NULL,
    *movinv=NULL;
  nifti_image
  *imgs[2];
  float
  *wx=NULL,*wy=NULL,*wz=NULL;
  imgs[0]=imgs[1]=NULL;
  if(verbose>=1) niik_fc_display(fcname,1);
  /* CHECK INPUTS */
  if(warpimg==NULL) {
    fprintf(stderr,"[%s] ERROR: warpimg is null\n",fcname);
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(tgtimg==NULL) {
    fprintf(stderr,"[%s] ERROR: tgtimg is null\n",fcname);
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"[%s] ERROR: movimg is null\n",fcname);
    return 0;
  }
  if(tgtmat==NULL) {
    fprintf(stderr,"[%s] ERROR: tgtmat is null\n",fcname);
    return 0;
  }
  if(movmat==NULL) {
    fprintf(stderr,"[%s] ERROR: movmat is null\n",fcname);
    return 0;
  }
  if(warpimg->nx!=refimg->nx) {
    fprintf(stderr,"[%s] ERROR: x-dimension for warpimg and refimg, %i %i\n",
            fcname,warpimg->nx,refimg->nx);
    return 0;
  }
  if(warpimg->ny!=refimg->ny) {
    fprintf(stderr,"[%s] ERROR: x-dimension for warpimg and refimg, %i %i\n",
            fcname,warpimg->ny,refimg->ny);
    return 0;
  }
  if(warpimg->nz!=refimg->nz) {
    fprintf(stderr,"[%s] ERROR: x-dimension for warpimg and refimg, %i %i\n",
            fcname,warpimg->nz,refimg->nz);
    return 0;
  }
  if(niik_image_cmp_dim(refimg,maskimg)!=0) {
    fprintf(stderr,"[%s] ERROR: image dimension for refimg and maskimg\n",fcname);
    return 0;
  }
  if(warpimg->nu!=3) {
    fprintf(stderr,"[%s] ERROR: warp image nu is not 3, %i\n",fcname,warpimg->nu);
    return 0;
  }
  if(warpimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: warp image datatype is not float32, %s\n",
            fcname,nifti_datatype_string(warpimg->datatype));
    return 0;
  }
  if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"[%s] ERROR: maskimg datatype is not uint8, %s\n",
            fcname,nifti_datatype_string(maskimg->datatype));
    return 0;
  }
  xdim=warpimg->nx;
  ydim=warpimg->ny;
  zdim=warpimg->nz;
  size=xdim*ydim*zdim;
  wx=warpimg->data;
  wy=wx+size;
  wz=wy+size;
  ww=niikmat_init(zdim,3);
  ws=niikmat_init(zdim,4);
  gv=niikmat_init(zdim,4);
  bimg=maskimg->data;
  /*niik_image_clear(refimg);*/
  tgtinv=niikmat_inverse(tgtmat);
  movinv=niikmat_inverse(movmat);
  imgs[0]=tgtimg;
  imgs[1]=movimg;
  if(verbose>=2) fprintf(stdout,"[%s] start main loop\n",fcname);
  #pragma omp parallel for private(i,j,n,pt,qt,p1,p2)
  for(k=0; k<zdim; k++) {
    pt.z=k*refimg->dz;
    n=k*xdim*ydim;
    for(j=0; j<ydim; j++) {
      pt.y=j*refimg->dy;
      for(i=0; i<xdim; n++,i++) {
        if(!bimg[n]) continue;
        pt.x=i*refimg->dx;
        niik_image_interpolate_3d_xyz_update(warpimg,pt,NIIK_INTERP_LINEAR,ww->m[k]);
        qt.x=pt.x+ww->m[k][0];
        qt.y=pt.y+ww->m[k][1];
        qt.z=pt.z+ww->m[k][2];
        p1=niikpt_affine_transform(tgtinv,pt);
        p2=niikpt_affine_transform(movinv,qt);
        niik_image_interpolate_3d_xyz_update(gradimg,p1,NIIK_INTERP_LINEAR,ws->m[k]);
        gv->m[k][0]=niik_image_interpolate_3d_linear(imgs[0],p1);
        gv->m[k][1]=niik_image_interpolate_3d_linear(imgs[1],p2);
        /*fprintf(stdout,"\t[%3i %3i %3i] %5.0f %5.0f | %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f : %7.3f\n",i,j,k,
          gv->m[k][0],gv->m[k][1],ww->m[k][0],ww->m[k][1],ww->m[k][2],ws->m[k][1],ws->m[k][2],ws->m[k][3],ws->m[k][0]);*/
        gv->m[k][2]=gv->m[k][1]-gv->m[k][0];
        gv->m[k][3]=gv->m[k][2] / ( ws->m[k][0]*ws->m[k][0] + NIIK_SQ(gv->m[k][2]) );
        /*niik_image_set_voxel(refimg,n,gv->m[k[2]);*/
        wx[n]-=gv->m[k][3]*ws->m[k][1];
        wy[n]-=gv->m[k][3]*ws->m[k][2];
        wz[n]-=gv->m[k][3]*ws->m[k][3];
        if(i==108)
          if(j==50)
            if(k==150)
              fprintf(stdout,"\t[%3i %3i %3i] %5.0f %5.0f | %7.3f %7.3f %7.3f | %5.0f %5.0f %5.0f : %5.0f | %7.3f %7.3f %7.3f\n",
                      i,j,k,gv->m[k][0],gv->m[k][1],ww->m[k][0],ww->m[k][1],ww->m[k][2],ws->m[k][1],ws->m[k][2],ws->m[k][3],ws->m[k][0],wx[n],wy[n],wz[n]);
      }
    }
  }
  ww=niikmat_free(ww);
  ws=niikmat_free(ws);
  gv=niikmat_free(gv);
  tgtinv=niikmat_free(tgtinv);
  movinv=niikmat_free(movinv);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_nregister_demons_update */


int niik_image_nregister_intrasubject_apply_warp(nifti_image *warpimg,nifti_image *refimg,
    nifti_image *movimg,niikmat *movmat,int interp)
/* -applies warp to movimg with movmat before warping
 * -refimg is modified/updated
 * -works with the functions in nifti1_kunio_nregister_intrasubject_template.c
 */
{
  nifti_image *tmpimg=NULL; /* for bspline */
  char fcname[128]="niik_image_nregister_intrasubject_apply_warp";
  int i,j,k,n;
  niikpt pt,qt;
  double ww[5];
  niikmat *movinv=NULL;
  pt.w=qt.w=0;
  if(warpimg==NULL) {
    fprintf(stderr,"[%s] ERROR: warpimg is null\n",fcname);
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"[%s] ERROR: movimg is null\n",fcname);
    return 0;
  }
  if(movmat==NULL) {
    fprintf(stderr,"[%s] ERROR: movmat is null\n",fcname);
    return 0;
  }
  if((movinv=niikmat_inverse(movmat))==NULL) {
    fprintf(stderr,"[%s] ERROR: niikmat_inverse\n",fcname);
    return 0;
  }
  if(interp==NIIK_INTERP_BSPLINE) {
    if((tmpimg=niik_image_interpolate_convert_3d_bspline_coeff2(movimg))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_interpolate_convert_3d_bspline_coeff2\n",fcname);
      return 0;
    }
  } else {
    tmpimg=movimg;
  }
  for(k=n=0; k<refimg->nz; k++) {
    pt.z=k*refimg->dz;
    for(j=0; j<refimg->ny; j++) {
      pt.y=j*refimg->dy;
      for(i=0; i<refimg->nx; n++,i++) {
        pt.x=i*refimg->dx;
        niik_image_interpolate_3d_xyz_update(warpimg,pt,NIIK_INTERP_LINEAR,ww);
        qt.x=pt.x+ww[0];
        qt.y=pt.y+ww[1];
        qt.z=pt.z+ww[2];
        niik_image_set_voxel(refimg,n,niik_image_interpolate_3d(tmpimg,niikpt_affine_transform(movinv,qt),interp));
      }
    }
  }
  movinv=niikmat_free(movinv);
  if(tmpimg!=movimg) {
    tmpimg=niik_image_free(tmpimg);
  }
  return 1;
} /* niik_image_nregister_intrasubject_apply_warp */


int niik_image_nregister_demons(nifti_image *warpimg,nifti_image *refimg,nifti_image *maskimg,
                                nifti_image *tgtimg,niikmat *tgtmat,nifti_image *tgtgrad,
                                nifti_image *movimg,niikmat *movmat,
                                int maxiter,   /* for intra-subject registration, small number is enough */
                                double wFWHM)  /* 2?  */

{

  char
  fname[4096],
        fcname[64]="niik_image_nregister_demons";
  int
  iter,
  verbose=1;
  nifti_image
  *gradimg=NULL;

  niik_fc_display(fcname,1);

  /* CHECK INPUTS */
  if(verbose>=2) fprintf(stdout,"[%s] check inputs\n",fcname);
  if(warpimg==NULL) {
    fprintf(stderr,"[%s] ERROR: warpimg is null\n",fcname);
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"[%s] ERROR: refimg is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(tgtimg==NULL) {
    fprintf(stderr,"[%s] ERROR: tgtimg is null\n",fcname);
    return 0;
  }
  if(movimg==NULL) {
    fprintf(stderr,"[%s] ERROR: movimg is null\n",fcname);
    return 0;
  }
  if(tgtmat==NULL) {
    fprintf(stderr,"[%s] ERROR: tgtmat is null\n",fcname);
    return 0;
  }
  if(movmat==NULL) {
    fprintf(stderr,"[%s] ERROR: movmat is null\n",fcname);
    return 0;
  }

  if(warpimg->nx!=refimg->nx) {
    fprintf(stderr,"[%s] ERROR: x-dimension for warpimg and refimg, %i %i\n",fcname,warpimg->nx,refimg->nx);
    return 0;
  }
  if(warpimg->ny!=refimg->ny) {
    fprintf(stderr,"[%s] ERROR: y-dimension for warpimg and refimg, %i %i\n",fcname,warpimg->ny,refimg->ny);
    return 0;
  }
  if(warpimg->nz!=refimg->nz) {
    fprintf(stderr,"[%s] ERROR: z-dimension for warpimg and refimg, %i %i\n",fcname,warpimg->nz,refimg->nz);
    return 0;
  }
  if(niik_image_cmp_dim(refimg,maskimg)!=0) {
    fprintf(stderr,"[%s] ERROR: image dimension for refimg and maskimg\n",fcname);
    return 0;
  }
  if(warpimg->nu!=3) {
    fprintf(stderr,"[%s] ERROR: warp image nu is not 3, %i\n",fcname,warpimg->nu);
    return 0;
  }
  if(warpimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: warp image datatype is not float32, %s\n",fcname,nifti_datatype_string(warpimg->datatype));
    return 0;
  }

  if(verbose>=1) fprintf(stdout,"[%s] gradient image from target image\n",fcname);
  if(tgtgrad!=NULL) {
    if(verbose>=1) fprintf(stdout,"[%s] no sobel fliter\n",fcname);
    gradimg=tgtgrad;
  } else {
    if(verbose>=1) fprintf(stdout,"[%s] sobel fliter %5.1f\n",fcname,1.5);
    if((gradimg=niik_image_gauss_sobel_filters_with_mag_single_output(tgtimg,1.5))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_sobel_filters_with_mag_single_output\n",fcname);
      return 0;
    }
  }
  /*
   * main loop
   */
  if(verbose>=2) fprintf(stdout,"[%s] start main loop\n",fcname);
  for(iter=0; iter<maxiter; iter++) {
    if(!niik_image_nregister_demons_update(warpimg,gradimg,refimg,maskimg,tgtimg,tgtmat,movimg,movmat)) {
      fprintf(stderr,"[%s] ERROR: niik_image_nregister_demons_update\n",fcname);
      return 0;
    }
    if(wFWHM>1e-4) {
      if(!niik_image_filter_gaussian_update(warpimg,wFWHM*2,wFWHM)) {
        fprintf(stderr,"[%s] ERROR: niik_image_filter_gaussian_update\n",fcname);
        return 0;
      }
    }
    if(verbose>=10) {
      fprintf(stdout,"[%s] writing output image\n",fcname);
      sprintf(fname,"tmp_nregimg_warp%i.nii.gz",iter+1);
      niik_image_write(fname,warpimg);
      sprintf(fname,"tmp_nregimg_diff%i.nii.gz",iter+1);
      niik_image_write(fname,refimg);
    }
  } /* iteration */
  /* clean up and finish */
  if(gradimg!=tgtgrad) {
    gradimg=niik_image_free(gradimg);
  }
  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_nregister_demons */

#endif /* _FALCON_NREGISTER_DEMONS_C_ */
