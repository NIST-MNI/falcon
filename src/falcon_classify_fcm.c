/* Filename:     nifti1_kunio_classify_fcm.c
 * Description:  fuzzy c-means classification function
 * Author:       Kunio Nakamura
 * Date:         January 2, 2013
 *
 *
 */

#ifndef _FALCON_CLASSIFY_FCM_C_
#define _FALCON_CLASSIFY_FCM_C_

#include "falcon.h"


niik_classify_fcm *niik_classify_fcm_init(nifti_image *img,nifti_image *mask,int num_class,int iter,int verbose) {
  niik_classify_fcm *C=NULL;
  char fcname[64]="niik_classify_fcm_init";
  int i=0;
  if(verbose>=0) niik_fc_display(fcname,1);
  C = (niik_classify_fcm *)calloc(1,sizeof(niik_classify_fcm));
  C->num_class = num_class;
  if(verbose>=1) fprintf(stdout,"[%s] copy mask\n",fcname);
  if((C->mask = niik_image_copy_as_type(mask,NIFTI_TYPE_UINT8))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s] copy img\n",fcname);
  if((C->img = niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  if((C->primg=niik_image_copy_as_type(mask,NIFTI_TYPE_FLOAT64))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type, %i\n",fcname,i);
    return NULL;
  }
  free(C->primg->data);
  C->primg->ndim=C->primg->dim[0]=5;
  C->primg->nt=C->primg->dim[4]=1;
  C->primg->nu=C->primg->dim[5]=C->num_class;
  C->primg->nvox=C->primg->nx*C->primg->ny*C->primg->nz*C->primg->nt*C->primg->nu;
  C->primg->data=(void *)calloc(C->primg->nvox,sizeof(double));
  C->mean = niikmat_init(img->nu,num_class);
  C->stdv = niikmat_init(img->nu,num_class);
  C->iter = iter;
  if(verbose>=0) niik_fc_display(fcname,0);
  return C;
} /* niik_classify_fcm_init */

niik_classify_fcm *niik_classify_fcm_free(niik_classify_fcm *C) {
  if(C==NULL) return NULL;
  if(C->mask!=NULL) C->mask=niik_image_free(C->mask);
  if(C->img !=NULL) C->img =niik_image_free(C->img );
  if(C->mean!=NULL) C->mean=niikmat_free(C->mean);
  if(C->stdv!=NULL) C->stdv=niikmat_free(C->stdv);
  free(C);
  return NULL;
} /* niik_classify_fcm_free */

int niik_classify_fcm_display(niik_classify_fcm *C) {
  char fcname[64]="niik_classify_fcm_display";
  int n,m;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: class is null\n",fcname);
    return 0;
  }
  niik_fc_display(fcname,1);
  fprintf(stdout,"  mask          :   %s\n",C->mask->fname);
  fprintf(stdout,"  img           :   %s\n",C->img->fname);
  fprintf(stdout,"  iter          :   %i\n",C->iter);
  fprintf(stdout,"  # img         :   %i\n",C->img->nu);
  fprintf(stdout,"  # class       :   %i\n",C->num_class);
  for(n=0; n<C->img->nu; n++) {
    fprintf(stdout,"  mean %3i      :   ",n);
    for(m=0; m<C->num_class; m++) {
      fprintf(stdout,"%12.9f %8.5f | ",C->mean->m[n][m],C->stdv->m[n][m]);
    }
    fprintf(stdout,"\n");
  }
  return 1;
} /* niik_classify_fcm_display */

int niik_classify_fcm_initialize_means(niik_classify_fcm *C,int verbose) {
  char fcname[64]="niik_classify_fcm_initialize_means";
  int ni,nc,i,j,nvox,nm;
  unsigned char *bimg=NULL;
  double
  mean,stdv,
       *dimg=NULL;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  bimg=C->mask->data;
  dimg=C->img->data;
  nvox=C->img->nx*C->img->ny*C->img->nz;
  nm=niik_image_count_mask(C->mask);
  for(ni=0; ni<C->img->nu; ni++) {
    C->mean->m[ni][0]=0;
    stdv = 0;
    for(i=0,j=ni*nvox; i<nvox; i++,j++) {
      if(bimg[i]==0) continue;
      C->mean->m[ni][0] += dimg[j];
      stdv += dimg[j]*dimg[j];
    }
    C->mean->m[ni][0] /= nm;
    mean = C->mean->m[ni][0];
    stdv = sqrt(stdv/nm - C->mean->m[ni][0]);
    for(nc=0; nc<C->num_class; nc++) {
      C->mean->m[ni][nc] = (nc-C->num_class/2.0)*stdv + mean;
    }
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_fcm_initialize_means */

int niik_classify_fcm_assign_probability(niik_classify_fcm *C,int verbose) {
  char fcname[64]="niik_classify_fcm_assign_classes";
  int
  ni,nc,nk,
  i,nvox;
  unsigned char *bimg=NULL;
  double
  *dlist,dsum,
  *dpr=NULL,
   *dimg=NULL;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  nvox=C->img->nx*C->img->ny*C->img->nz;
  bimg=C->mask->data;
  dimg=C->img->data;
  dpr =C->primg->data;
  dlist=(double *)calloc(C->num_class,sizeof(double));
  for(i=0; i<nvox; i++) {
    if(bimg[i]==0) continue;
    for(nc=0; nc<C->num_class; nc++) {
      dlist[nc]=0;
      for(ni=0; ni<C->img->nu; ni++) {
        dlist[nc]+=NIIK_SQ(dimg[i+ni*nvox]-C->mean->m[ni][nc]);
      }
      dlist[nc]=sqrt(dlist[nc]);
    }
    /*fprintf(stdout,"[%3i %3i %3i] %8.3f %8.3f %8.3f\n",i%C->img->nx,(i/C->img->nx)%C->img->ny,i/C->img->nx/C->img->ny,dlist[0],dlist[1],dlist[2]); */
    for(nc=0; nc<C->num_class; nc++) {
      for(nk=0,dsum=1e-6; nk<C->num_class; nk++) {
        dsum+=NIIK_SQ(dlist[nc]/(1e-6+dlist[nk]));
      }
      dpr[i+nc*nvox]=1.0/dsum;
    }
    /*fprintf(stdout,"[%3i %3i %3i] %8.3f %8.3f %8.3f | %8.3f %8.3f %8.3f\n",
            i%C->img->nx,(i/C->img->nx)%C->img->ny,i/C->img->nx/C->img->ny,dlist[0],dlist[1],dlist[2],
            dpr[i],dpr[i+nvox],dpr[i+nvox*2]); */
  }
  free(dlist);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_fcm_assign_classes */

int niik_classify_fcm_update_means(niik_classify_fcm *C,int verbose) {
  char fcname[64]="niik_classify_fcm_update_means";
  int
  ni,nc,
  i,nvox;
  unsigned char *bimg=NULL;
  double
  *dlist,dsum,
  *dpr=NULL,
   *dimg=NULL;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  nvox=C->img->nx*C->img->ny*C->img->nz;
  bimg=C->mask->data;
  dimg=C->img->data;
  dpr =C->primg->data;
  dlist=(double *)calloc(C->num_class,sizeof(double));
  for(nc=0; nc<C->num_class; nc++) {
    for(ni=0; ni<C->img->nu; ni++) {
      C->mean->m[ni][nc]=1e-6;
      C->stdv->m[ni][nc]=1e-6;
      dsum=1e-6;
      for(i=0; i<nvox; i++) {
        if(bimg[i]==0) continue;
        C->mean->m[ni][nc]+=dimg[i+ni*nvox]*dpr[i+nc*nvox];
        dsum+=dpr[i+nc*nvox];
      }
      C->mean->m[ni][nc]/=dsum;
      C->stdv->m[ni][nc]=1e-6;
      dsum=1e-6;
      for(i=0; i<nvox; i++) {
        if(bimg[i]==0) continue;
        C->stdv->m[ni][nc]+=NIIK_SQ(dimg[i+ni*nvox]-C->mean->m[ni][nc])*dpr[i+nc*nvox];
        dsum+=dpr[i+nc*nvox];
      }
      C->stdv->m[ni][nc]=sqrt(C->stdv->m[ni][nc]/dsum);
    } /* each dimension */
  } /* each class */
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_update_means */

int niik_classify_fcm_iterate(niik_classify_fcm *C,int verbose) {
  char fcname[64]="niik_classify_fcm_iterate";
  int
  iter=0;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  for(iter=0; iter<C->iter; iter++) {
    niik_classify_fcm_assign_probability(C,1);
    niik_classify_fcm_update_means(C,1);
    niik_classify_fcm_display(C);
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_update_means */


#endif /*  _FALCON_CLASSIFY_FCM_C_ */
