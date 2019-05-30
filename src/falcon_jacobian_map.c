/* Filename:     nifti1_kunio_jacobian_map.c
 * Description:  Jacobian Map
 * Author:       Kunio Nakamura
 * Date:         September 29, 2012
 */

#ifndef _FALCON_JACOBIAN_MAP_C_
#define _FALCON_JACOBIAN_MAP_C_

#include "falcon.h"

nifti_image *niik_image_jacobian_map(nifti_image *warpimg,nifti_image *maskimg,int warptype)
/* calculate Jacobian map
 *
 */
{

  nifti_image *jmap=NULL;
  char fcname[64]="niik_image_jacobian_map";
  int
  xdim,area,size,
       i,j,k,n;
  niikmat *M;
  float *fx,*fy,*fz,*fj;

  if(warpimg==NULL) {
    fprintf(stderr,"[%s] ERROR: warpimg is null\n",fcname);
    return NULL;
  }

  if((jmap=niik_image_copy_as_type(warpimg,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }

  M=niikmat_init(3,3);
  jmap->dim[0]=jmap->ndim=3;
  jmap->nvox=jmap->nx*jmap->ny*jmap->nz;
  fj=jmap->data;

  if(!niik_image_type_convert(warpimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return NULL;
  }

  xdim = warpimg->nx;
  area = warpimg->nx*warpimg->ny;
  size = warpimg->nx*warpimg->ny*warpimg->nz;
  fx=warpimg->data;
  fy=fx+size;
  fz=fy+size;

  switch(warptype) {
  case NIIK_WARP_MAP_DISP:

    fprintf(stdout,"[%s] using displacement map\n",fcname);
    for(k=n=0; k<warpimg->nz; k++) {
      for(j=0; j<warpimg->ny; j++) {
        for(i=0; i<warpimg->nx; n++,i++) {

          if(i==0) {
            M->m[0][0]=(fx[n+1]-fx[n  ])/warpimg->dx;
            M->m[1][0]=(fy[n+1]-fy[n  ])/warpimg->dx;
            M->m[2][0]=(fz[n+1]-fz[n  ])/warpimg->dx;
          } else if(i==warpimg->nx-1) {
            M->m[0][0]=(fx[n  ]-fx[n-1])/warpimg->dx;
            M->m[1][0]=(fy[n  ]-fy[n-1])/warpimg->dx;
            M->m[2][0]=(fz[n  ]-fz[n-1])/warpimg->dx;
          } else {
            M->m[0][0]=(fx[n+1]-fx[n-1])/warpimg->dx/2.0;
            M->m[1][0]=(fy[n+1]-fy[n-1])/warpimg->dx/2.0;
            M->m[2][0]=(fz[n+1]-fz[n-1])/warpimg->dx/2.0;
          }

          if(j==0) {
            M->m[0][1]=(fx[n+xdim]-fx[n  ])/warpimg->dy;
            M->m[1][1]=(fy[n+xdim]-fy[n  ])/warpimg->dy;
            M->m[2][1]=(fz[n+xdim]-fz[n  ])/warpimg->dy;
          } else if(j==warpimg->ny-1) {
            M->m[0][1]=(fx[n  ]-fx[n-xdim])/warpimg->dy;
            M->m[1][1]=(fy[n  ]-fy[n-xdim])/warpimg->dy;
            M->m[2][1]=(fz[n  ]-fz[n-xdim])/warpimg->dy;
          } else {
            M->m[0][1]=(fx[n+xdim]-fx[n-xdim])/warpimg->dy/2.0;
            M->m[1][1]=(fy[n+xdim]-fy[n-xdim])/warpimg->dy/2.0;
            M->m[2][1]=(fz[n+xdim]-fz[n-xdim])/warpimg->dy/2.0;
          }

          if(k==0) {
            M->m[0][2]=(fx[n+area]-fx[n  ])/warpimg->dz;
            M->m[1][2]=(fy[n+area]-fy[n  ])/warpimg->dz;
            M->m[2][2]=(fz[n+area]-fz[n  ])/warpimg->dz;
          } else if(k==warpimg->nz-1) {
            M->m[0][2]=(fx[n  ]-fx[n-area])/warpimg->dz;
            M->m[1][2]=(fy[n  ]-fy[n-area])/warpimg->dz;
            M->m[2][2]=(fz[n  ]-fz[n-area])/warpimg->dz;
          } else {
            M->m[0][2]=(fx[n+area]-fx[n-area])/warpimg->dz/2.0;
            M->m[1][2]=(fy[n+area]-fy[n-area])/warpimg->dz/2.0;
            M->m[2][2]=(fz[n+area]-fz[n-area])/warpimg->dz/2.0;
          }

          M->m[0][0] += 1.0;
          M->m[1][1] += 1.0;
          M->m[2][2] += 1.0;

          fj[n] = niikpt_det3(M->m[0][0],M->m[0][1],M->m[0][2],
                              M->m[1][0],M->m[1][1],M->m[1][2],
                              M->m[2][0],M->m[2][1],M->m[2][2]);
          /*fj[n] =  niikpt_det3(M->m[0][0],M->m[1][0],M->m[2][0],
            M->m[0][1],M->m[1][1],M->m[2][1],
            M->m[0][2],M->m[1][2],M->m[2][2]);*/
          /*fj[n] = nifti_mat33_determ(M);*/
        }
      }
    }
    break;
  case NIIK_WARP_MAP_LOC:
    break;
  case NIIK_WARP_MAP_DISP_BSPLINE:
    break;
  }

  M=niikmat_free(M);

  if(maskimg!=NULL) {
    fprintf(stdout,"[%s] mask size : %i\n",fcname,niik_image_count_mask(maskimg));
    fprintf(stdout,"[%s] JACOBIAN INTEGRATION %19.9f using %s\n",fcname,
            niik_image_get_mean(jmap,maskimg),maskimg->fname);
  }

  return jmap;
} /* niik_image_jacobian_map */

#endif /* _FALCON_JACOBIAN_MAP_C_ */
