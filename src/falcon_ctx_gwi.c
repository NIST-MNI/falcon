/* Filename:     nifti1_kunio_ctx_gwi.c
 * Description:  Gray matter - white matter interface mask
 * Author:       Kunio Nakamura
 * Date:         February 29, 2012
 */


#ifndef _FALCON_CTX_GWI_C_
#define _FALCON_CTX_GWI_C_

#include "falcon.h"


nifti_image *niik_image_ctx_gwi(nifti_image **imglist) {
  nifti_image
  *img,*segimg,*venimg,*wmmask,*dgmmask,*cerebellummask,*brainstemmask,
  *outimg=NULL;

  fprintf(stdout,"[niik_image_ctx_gwi] start\n");

  if(imglist==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }

  img=imglist[0];
  segimg=imglist[1];
  wmmask=imglist[2];
  venmask=imglist[3];
  dgmmask=imglist[4];
  cerebellummask=imglist[5];
  brainstemmask=imglist[6];


  fprintf(stdout,"[niik_image_ctx_gwi] finish\n");
  return outimg;
} /* niik_image_ctx_gwi */


#endif
