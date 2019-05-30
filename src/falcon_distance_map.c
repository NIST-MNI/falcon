/* Filename:     nifti1_kunio_distance_map.c
 * Description:  distance map function
 * Author:       Kunio Nakamura
 * Date:         April 2, 2012
 */

#ifndef _FALCON_DISTANCE_MAP_C_
#define _FALCON_DISTANCE_MAP_C_

#include "falcon.h"



/******************************************************************
 * distance map
 *
 * -calculates the distance map from mask image (Positive = object and Zero = background)
 * -max_dist sets the limit of distance (because the distance near the edge
 *  is frequently needed and not far away from the edge (as in level set method)
 * -compared to my previous versions, the edge is assumed to be between voxel
 * -so the voxel at the edge has distance of half pixel size
 *
 ******************************************************************/

nifti_image *niik_image_distance_map(nifti_image *img,double max_dist) {
  nifti_image
  *tmpimg=NULL,
   *bakimg=NULL,
    *outimg=NULL;
  int
  kd,xdim,dim12,
  ii,jj,kk,
  m,n,i,j,k;
  double
  *dimg,*kimg,
  dval,max2;
  unsigned char *bimg;
  int verbose=1;

  if(verbose>=1) {
    niik_fc_display("niik_image_distance_map",1);
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }

  dval = NIIK_DMIN(img->pixdim[1],NIIK_DMIN(img->pixdim[2],img->pixdim[3]));
  if(max_dist<dval) {
    fprintf(stderr,"ERROR: max_dist %f is too small\n",max_dist);
    return 0;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] copy image \n");
  }
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] type convert image\n");
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT64)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(outimg,%s) \n",nifti_datatype_string(NIFTI_TYPE_FLOAT64));
    return NULL;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] type convert image\n");
  }
  if((tmpimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return NULL;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] type convert image\n");
  }
  if(!niik_image_type_convert(tmpimg,NIFTI_TYPE_UINT8)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(tmpimg,%s) \n",nifti_datatype_string(NIFTI_TYPE_UINT8));
    return NULL;
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] type convert image\n");
  }
  if((bakimg=niik_image_copy(outimg))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy \n");
    return NULL;
  }

  bimg = tmpimg -> data;
  kimg = bakimg -> data;
  dimg = outimg -> data;
  xdim = outimg -> dim[1];
  dim12 = xdim * outimg -> dim[2];

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] clear image\n");
  }
  max2 = max_dist*max_dist;
  for(i=0; i<outimg->nvox; i++) dimg[i]=max2;
  for(i=0; i<outimg->nvox; i++) bimg[i]=(bimg[i]>0);

  /* x-dir */
  kd=floor(max_dist/img->dx+1.5);
  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] x-dir %i\n",kd);
  }
  for(k=n=0; k<img->dim[3]; k++) {
    for(j=0; j<img->dim[2]; j++) {
      for(i=0; i<img->dim[1]; n++,i++) {
        for(ii=1,m=kd; ii<=kd; ii++) {
          if(ii+i<0) continue;
          if(ii+i>=img->dim[1]) continue;
          /* if(i==114 && j==153 && k==93) fprintf(stdout,"%3i %3i %3i %5.2f %4i \n",i,j,k,dimg[n],ii); */
          if(bimg[n+ii]!=bimg[n]) {
            dimg[n]=NIIK_DMIN(NIIK_SQ(img->pixdim[1]*(ii-0.5)),dimg[n]);
            m=ii;
            break;
          }
        }
        for(ii=-1; ii>=-m; ii--) {
          if(ii+i<0) continue;
          if(ii+i>=img->dim[1]) continue;
          /* if(i==114 && j==153 && k==93) fprintf(stdout,"%3i %3i %3i %5.2f %4i \n",i,j,k,dimg[n],ii);*/
          if(bimg[n+ii]!=bimg[n]) {
            dimg[n]=NIIK_DMIN(NIIK_SQ(img->pixdim[1]*(ii-0.5)),dimg[n]);
          }
        }
      }
    }
  }  /* x-dir */
  for(i=0; i<img->nvox; i++) {
    kimg[i]=dimg[i];
  }


  /* y-dir */
  kd=floor(max_dist/img->dy+1.5);
  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] y-dir %i\n",kd);
  }
  for(k=0; k<img->dim[3]; k++) {
    for(i=0; i<img->dim[1]; n++,i++) {
      n=k*dim12+i;
      for(j=0; j<img->dim[2]; n+=xdim,j++) {
        for(jj=-kd,m=n-kd*xdim; jj<=kd; m+=xdim,jj++) {
          if(jj+j<0) continue;
          if(jj==0)  continue;
          if(jj+j>=img->dim[2]) break;
          dval = NIIK_SQ(img->pixdim[2]*(jj-0.5));
          /*if(i==114 && j==153 && k==93) fprintf(stdout,"%3i %3i %3i %5.2f %4i \n",i,j,k,dimg[n],jj);*/
          if(dval>=dimg[n]) continue;
          if(bimg[m]!=bimg[n]) {
            dimg[n]=NIIK_DMIN(dval,dimg[n]);
          } else if(kimg[m]>kimg[n]) continue;
          else
            dimg[n]=NIIK_DMIN(dimg[n],kimg[m]+dval);
        }
      }
    }
  } /* y-dir */
  for(i=0; i<img->nvox; i++) {
    kimg[i]=dimg[i];
  }


  /* z-dir */
  kd=floor(max_dist/img->dz+1.5);
  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] z-dir %i\n",kd);
  }
  for(i=0; i<img->dim[1]; n++,i++) {
    for(j=0; j<img->dim[2]; n+=xdim,j++) {
      n=j*xdim+i;
      for(k=0; k<img->dim[3]; n+=dim12,k++) {
        for(kk=-kd,m=n-kd*dim12; kk<=kd; m+=dim12,kk++) {
          if(kk+k<0) continue;
          if(kk==0)  continue;
          if(kk+k>=img->dim[3]) break;
          dval = NIIK_SQ(img->pixdim[3]*(kk-0.5));
          /*if(i==114 && j==153 && k==93) fprintf(stdout,"%3i %3i %3i %5.2f %4i \n",i,j,k,dimg[n],kk);*/
          if(dval>=dimg[n]) continue;
          if(bimg[m]!=bimg[n]) {
            dimg[n]=NIIK_DMIN(dval,dimg[n]);
          } else if(kimg[m]>kimg[n]) continue;
          else
            dimg[n]=NIIK_DMIN(dimg[n],kimg[m]+dval);
        }
      }
    }
  } /* z-dir */

  for(i=0; i<img->nvox; i++) {
    dimg[i]=sqrt(dimg[i])*((bimg[i]>0)?1:-1);
  }

  if(verbose>1) {
    fprintf(stdout,"[niik_image_distance_map] convert outimg\n");
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"ERROR: niik_image_type_convert(outimg,%s) \n",nifti_datatype_string(NIFTI_TYPE_FLOAT32));
    return NULL;
  }

  nifti_image_free(bakimg);
  nifti_image_free(tmpimg);
  if(verbose>=1) {
    niik_fc_display("niik_image_distance_map",0);
  }
  return outimg;
}


#endif /* NIFTI1_KUNIO_DISTANCE_MAP_C_ */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/