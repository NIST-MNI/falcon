/* FILENAME:     nifti1_kunio_segment.c
 * DESCRIPTION:  general segmentation/classification function
 * REFERENCE:    Cuadra MB, Cammoun L, Butz T, Cuisenaire O, Thiran JP, 2005. Comparison and Validation of Tissue Modelization and Statistical Classification Methods in T1-Weighted MR Brain Images. IEEE Transactions on Medical Imaging. 24(12); 1548-65.
 * Author:       Kunio Nakamura
 * Date:         August 24, 2012
 */

#ifndef _FALCON_SEGMENT_C_
#define _FALCON_SEGMENT_C_

#include "falcon.h"

int niik_image_classify_FGMM(nifti_image *img,nifti_image *mask,int nclass,niikmat *stats,int maxiter)
/* -img is the grayscale image to be classified
 * -mask is the mask (ROI) to be classified; output is replaced here
 * -nclass is the number of classes
 * -stats is the initial states
 *  stats->m[nclass][0] = probability
 *  stats->m[nclass][1] = means
 *  stats->m[nclass][2] = stdv */
{
  char
  fcname[64]="niik_image_classify_FGMM";
  int
  iter,
  i,j,n,
  nvox,
  verbose=1;
  niikmat *pr;
  niikvec *vec;
  double dval,dsum;

  if(verbose>=1) fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(mask==NULL) {
    fprintf(stderr,"ERROR: mask is null\n");
    return 0;
  }
  if(stats==NULL) {
    fprintf(stderr,"ERROR: stats is null\n");
    return 0;
  }

  if(verbose>=2) {
    fprintf(stdout,"[%s] initial stats\n",fcname);
    niikmat_display(stats);
  }

  nvox = niik_image_count_mask(mask);
  if(verbose>=2) fprintf(stdout,"[%s] #voxel %i\n",fcname,nvox);
  pr = niikmat_init(nclass,nvox);

  for(iter=0; iter<maxiter; iter++) {

    /* E-step */
    for(i=j=0; i<img->nvox; i++) {
      if(niik_image_get_voxel(mask,i)==0) continue;
      for(n=0; n<nclass; n++) {
        pr->m[n][j] =
          /*stats->m[n][0] **/
          NIIK_GaussPDF(niik_image_get_voxel(img,i)-stats->m[n][1],
                        stats->m[n][2]);
      }
      if(verbose>=2) {
        if(j==300) {
          if(nclass==3)
            fprintf(stdout,"pr->m[300] = %8.3f -> %8.3f %8.3f %8.3f\n",niik_image_get_voxel(img,i),pr->m[0][300],pr->m[1][300],pr->m[2][300]);
          else if(nclass==2)
            fprintf(stdout,"pr->m[300] = %8.3f -> %8.3f %8.3f\n",niik_image_get_voxel(img,i),pr->m[0][300],pr->m[1][300]);
        }
      }
      for(n=0,dsum=1e-3; n<nclass; n++) {
        dsum+=pr->m[n][j];
      }
      for(n=0; n<nclass; n++) {
        pr->m[n][j]/=dsum;
      }
      j++;
    }  /* each voxel */

    if(verbose>=2) {
      if(nclass==3)
        fprintf(stdout,"pr->m[300] =             %8.3f %8.3f %8.3f\n",pr->m[0][300],pr->m[1][300],pr->m[2][300]);
      else if(nclass==2)
        fprintf(stdout,"pr->m[300] =             %8.3f %8.3f\n",pr->m[0][300],pr->m[1][300]);
    }

    /* M-step */
    for(n=0; n<nclass; n++) {
      stats->m[n][0]=stats->m[n][1]=stats->m[n][2]=0;
      for(i=j=0; i<img->nvox; i++) {
        if(niik_image_get_voxel(mask,i)==0) continue;
        stats->m[n][0] += pr->m[n][j];
        stats->m[n][1] += pr->m[n][j] * niik_image_get_voxel(img,i);
        stats->m[n][2] += pr->m[n][j] * NIIK_SQ(niik_image_get_voxel(img,i));
        j++;
      }
      stats->m[n][1] /= stats->m[n][0];
      stats->m[n][2] = sqrt(stats->m[n][2] / stats->m[n][0] - NIIK_SQ(stats->m[n][1]));
      stats->m[n][0] /= j;
      if(verbose>=1) fprintf(stdout,"[%s] class %i   %8.3f %8.3f %8.3f\n",fcname,n,stats->m[n][0],stats->m[n][1],stats->m[n][2]);
    } /* each voxel, M-step */
  } /* iteration */


  if(verbose>=1) fprintf(stdout,"[%s] create a class image  (includes all voxels)\n",fcname);
  vec=niikvec_init(nclass);
  for(i=0; i<img->nvox; i++) {
    for(n=0; n<nclass; n++) {
      pr->m[n][0] =
        NIIK_GaussPDF(niik_image_get_voxel(img,i)-stats->m[n][1],
                      stats->m[n][2]);
    }
    for(n=1,j=0; n<nclass; n++) {
      if(pr->m[j][0]<pr->m[n][0]) {
        j=n;
      }
    }
    if(niik_image_get_voxel(mask,i)>0)
      vec->v[j]+=1;
    niik_image_set_voxel(mask,i,j+1);
  }
  dval = niik_image_get_voxel_size(mask) / 1000.0;
  for(n=0; n<nclass; n++) {
    vec->v[n] *= dval;
  }
  if(verbose>=1) fprintf(stdout,"[%s] volume within mask: ",fcname);
  niikvec_display(vec);

  pr=niikmat_free(pr);
  vec=niikvec_free(vec);
  if(verbose>=1) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
} /* niik_image_classify_FGMM */

#endif /* _FALCON_SEGMENT_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/