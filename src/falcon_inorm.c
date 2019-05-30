/* Filename:     nifti1_kunio_inorm.c
 * Description:  functions for intensity normalization
 * Author:       Kunio Nakamura
 * Date:         April 11, 2012
 *
 */

#ifndef _FALCON_INORM_C_
#define _FALCON_INORM_C_

#include "falcon.h"

int niik_image_linear_normalization(nifti_image *refimg, nifti_image *modimg, nifti_image *modmask, niikmat *mod2ref) {
  const char *fcname="niik_image_linear_normalization";
  double
  x,y,
  S,Y0,E,
  sx,sy,sxx,sxy,syy;
  niikpt
  pfov,
  p,q;
  int
  verbose=0,
  m,n,i,j,k;
  unsigned char
  *bimg=NULL;
  niikmat *afmat=NULL;

  NIIK_RET0((refimg==NULL),fcname,"refimg is null");
  NIIK_RET0((modimg==NULL),fcname,"refimg is null");

  sx = sy = sxx = sxy = syy = 0;
  if(modmask!=NULL) {
    if(modimg->nvox!=modmask->nvox) {
      fprintf(stderr,"ERROR: different nvox %i %i\n",modimg->nvox,modmask->nvox);
      return 0;
    }
    if((bimg = niik_image_get_voxels_as_uint8_vector(modmask))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
      return 0;
    }
  } else {
    if((bimg = (unsigned char *)calloc(modimg->nvox,sizeof(char)))==NULL) {
      fprintf(stderr,"ERROR: calloc\n");
      return 0;
    }
    for(i=0; i<modimg->nvox; i++) {
      bimg[i]=1;
    }
  }

  if(mod2ref==NULL) {
    if((afmat = niikmat_identity(4,4))==NULL) {
      fprintf(stderr,"ERROR: niikmat_identity\n");
      return 0;
    }
  } else {
    if((afmat = niikmat_copy(mod2ref))==NULL) {
      fprintf(stderr,"ERROR: niikmat_copy\n");
      return 0;
    }
  }

  pfov = niikpt_val(refimg->dx*(refimg->nx-1),
                    refimg->dy*(refimg->ny-1),
                    refimg->dz*(refimg->nz-1),
                    1);
  p.w=0.0;
  for(k=m=n=0; k<modimg->nz; k++) {
    p.z = k*modimg->dz;
    for(j=0; j<modimg->ny; j++) {
      p.y = j*modimg->dy;
      for(i=0; i<modimg->nx; n++,i++) {
        p.x = i*modimg->dx;
        if(!bimg[n]) continue;
        q = niikpt_affine_transform(afmat,p);
        if(q.x<0) continue;
        else if(q.x>pfov.x) continue;
        if(q.y<0) continue;
        else if(q.y>pfov.y) continue;
        if(q.z<0) continue;
        else if(q.z>pfov.z) continue;
        x = niik_image_interpolate_3d(modimg,p,NIIK_INTERP_NN);
        y = niik_image_interpolate_3d(refimg,q,NIIK_INTERP_LINEAR);
        sx += x;
        sy += y;
        sxx += x*x;
        syy += y*y;
        sxy += x*y;
        m++;
      }
    }
  }

  /*S = (m*sxy - sx*sy) / (m*sxx - sx*sx);*/
  S = sxy / sxx;
  Y0 = sy/m - S*sx/m;
  if(verbose>=0) {
    fprintf(stdout,"[%s] mean mod = %12.6f\n",fcname,sx/m);
    fprintf(stdout,"[%s] mean ref = %12.6f\n",fcname,sy/m);
    fprintf(stdout,"[%s] SLOPE %12.8f\n",fcname,S);
    fprintf(stdout,"[%s] INTERCEPT %12.8f\n",fcname,Y0);
    for(k=m=n=0,E=0; k<modimg->nz; k++) {
      p.z = k*modimg->dz;
      for(j=0; j<modimg->ny; j++) {
        p.y = j*modimg->dy;
        for(i=0; i<modimg->nx; n++,i++) {
          p.x = i*modimg->dx;
          if(!bimg[n]) continue;
          q = niikpt_affine_transform(afmat,p);
          if(q.x<0) continue;
          else if(q.x>pfov.x) continue;
          if(q.y<0) continue;
          else if(q.y>pfov.y) continue;
          if(q.z<0) continue;
          else if(q.z>pfov.z) continue;
          x = niik_image_interpolate_3d(modimg,p,NIIK_INTERP_NN);
          y = niik_image_interpolate_3d(refimg,q,NIIK_INTERP_LINEAR);
          E += fabs(y-x*S-Y0);
        }
      }
    }
    fprintf(stdout,"[%s] ERROR %12.8f\n",fcname,E);
  }

  /* apply scaling */
  for(i=0; i<modimg->nvox; i++) {
    y = S * niik_image_get_voxel(modimg,i) + Y0;
    niik_image_set_voxel(modimg,i,y);
  }

  if(verbose>0) fprintf(stdout,"[%s] free\n",fcname);
  afmat = niikmat_free(afmat);
  free(bimg);
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
}



#endif /* _FALCON_INORM_C_ */
