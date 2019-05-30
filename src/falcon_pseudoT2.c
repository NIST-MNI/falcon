/* Filename:     nifti1_kunio_psuedoT2.c
 * Description:  pseudoT2 calculation functions
 * Author:       Kunio Nakamura
 * Date:         November 7, 2012
 */

#ifndef _FALCON_PSEUDOT2_C_
#define _FALCON_PSEUDOT2_C_

#include "falcon.h"

#define EPS 0.001

int niik_image_calculate_pseudoT2(nifti_image *PDimg,double PD_TE,double PD_TR,nifti_image *T2img,double T2_TE,double T2_TR,double omax,nifti_image *maskimg,nifti_image *outimg) {
  char fcname[64]="niik_image_calculate_pseudoT2";
  int
  verbose=0,
  i;
  float *fPD=NULL,*fT2=NULL,*fout=NULL;
  unsigned char *bimg=NULL;

  niik_fc_display(fcname,1);

  if(verbose>=1) fprintf(stdout,"[%s] check inputs\n",fcname);
  if(PDimg==NULL) {
    fprintf(stderr,"[%s] ERROR: PDimg is null\n",fcname);
    return 0;
  }
  if(T2img==NULL) {
    fprintf(stderr,"[%s] ERROR: T2img is null\n",fcname);
    return 0;
  }
  if(niik_image_cmp_dim(PDimg,T2img)!=0) {
    fprintf(stderr,"[%s] ERROR: wrong dimensions\n",fcname);
    return 0;
  }
  if(maskimg!=NULL) {
    if(niik_image_cmp_dim(PDimg,maskimg)!=0) {
      fprintf(stderr,"[%s] ERROR: wrong dimensions (mask)\n",fcname);
      return 0;
    }
  }
  if(niik_image_cmp_dim(PDimg,outimg)!=0) {
    fprintf(stderr,"[%s] ERROR: wrong dimensions (output)\n",fcname);
    return 0;
  }
  if(!niik_image_type_convert(outimg,NIFTI_TYPE_FLOAT32)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }

  if(verbose>=1) fprintf(stdout,"[%s] get data\n",fcname);
  fout=outimg->data;
  if((fPD=niik_image_get_voxels_as_float_vector(PDimg))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_float_vector\n",fcname);
    return 0;
  }
  if((fT2=niik_image_get_voxels_as_float_vector(T2img))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_float_vector\n",fcname);
    return 0;
  }


  if(maskimg!=NULL) {
    if(verbose>=1) fprintf(stdout,"[%s] convert mask to uin8\n",fcname);
    if(!niik_image_type_convert(maskimg,NIFTI_TYPE_UINT8)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
      return 0;
    }
    bimg=maskimg->data;
    for(i=0; i<PDimg->nvox; i++) {
      fout[i]=0;
      if(bimg[i]==0) {
        continue;
      }
      fout[i]=log((EPS+fPD[i])/(EPS+fT2[i]));
      if(fout[i]<EPS) continue;
      fout[i]=(T2_TE - PD_TE) / (EPS + fout[i]);
      if(fout[i]>omax) fout[i]=omax;
    }
  } /* using mask */

  else {
    fprintf(stdout,"[%s] WARNING: no mask !\n",fcname);
    for(i=0; i<PDimg->nvox; i++) {
      fout[i]=log((EPS+fPD[i])/(EPS+fT2[i]));
      if(fout[i]<EPS) continue;
      fout[i]=(T2_TE - PD_TE) / (EPS + fout[i]);
    }
  } /* using mask */

  if(verbose>=1) fprintf(stdout,"[%s] free memory\n",fcname);
  free(fPD);
  free(fT2);

  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_calculate_pseudoT2 */

#undef EPS

#endif /* _FALCON_PSEUDOT2_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/