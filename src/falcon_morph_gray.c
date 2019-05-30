/* Filename:     nifti1_kunio_morph_gray.c
 * Description:  grayscale-morphology
 * Author:       Kunio Nakamura
 * Date:         August 2, 2013
 */

#ifndef _FALCON_MORPH_GRAY_C_
#define _FALCON_MORPH_GRAY_C_

#include "falcon.h"

int niik_image_morph_gray_dilate(nifti_image *img,double radius,double val)

{
  char fcname[32]="niik_image_morph_gray_dilate";
  int verbose=1;
  int
  im,jm,km,
  ii,jj,kk,xydim,
  a,b,c,
  i,j,k,m,n;
  nifti_image
  *tmpimg=NULL;
  double
  cval,
  r2,
  dx,dy,dz,
  *dimg=NULL;

  if(verbose>1) niik_fc_display(fcname,1);

  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0(((tmpimg=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL),fcname,"niik_image_copy_as_type");
  r2=radius*radius;
  im=(int)floor(0.5+radius/img->dx);
  jm=(int)floor(0.5+radius/img->dy);
  km=(int)floor(0.5+radius/img->dz);
  dimg=tmpimg->data;
  xydim=img->nx*img->ny;

  for(k=m=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; m++,i++) {
        cval = niik_image_get_voxel(img,m);

        for(kk=-km; kk<=km; kk++) {
          c=kk+k;
          if(c<0) continue;
          if(c>=img->nz) break;
          dz=kk*img->dz;
          dz*=dz;
          if(dz>r2) continue;

          for(jj=-jm; jj<=jm; jj++) {
            b=jj+j;
            if(b<0) continue;
            if(b>=img->ny) break;;
            dy=jj*img->ny;
            dy*=dy;
            if(dz+dy>r2) continue;

            for(ii=-im,n=jj*img->nx+kk*xydim-im; ii<=im; n++,ii++) {
              a=ii+i;
              if(a<0) continue;
              if(a>=img->nx) break;
              dx=ii*img->dx;
              dx*=dx;
              if(dx+dy+dz>r2) continue;

              if(cval<dimg[m+n]+val)
                cval = dimg[m+n]+val;
            }
          }
        } /* neighbourhood */
        niik_image_set_voxel(img,m,cval);
      }
    }
  } /* for each voxel */

  tmpimg=niik_image_free(tmpimg);
  if(verbose>1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_morph_gray_dilate */


#endif  /* _FALCON_MORPH_GRAY_C_ */
