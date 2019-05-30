/* Filename:     nifti1_kunio_classify_kmeans.c
 * Description:  kmeans classification function
 * Author:       Kunio Nakamura
 * Date:         January 1, 2013
 *
 *
 */

#ifndef _FALCON_CLASSIFY_KMEANS_C_
#define _FALCON_CLASSIFY_KMEANS_C_

#include "falcon.h"


niik_classify_kmeans *niik_classify_kmeans_init(nifti_image *img,nifti_image *mask,int num_class,int iter,int verbose) {
  niik_classify_kmeans *K=NULL;
  char fcname[64]="niik_classify_kmeans_init";
  if(verbose>=0) niik_fc_display(fcname,1);
  K = (niik_classify_kmeans *)calloc(1,sizeof(niik_classify_kmeans));
  if(verbose>=1) fprintf(stdout,"[%s] copy mask\n",fcname);
  if((K->mask = niik_image_copy_as_type(mask,NIFTI_TYPE_UINT8))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  if(verbose>=1) fprintf(stdout,"[%s] copy img\n",fcname);
  if((K->img = niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT64))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }
  K->num_class = num_class;
  K->mean = niikmat_init(img->nu,num_class);
  K->wts  = niikmat_init(img->nu,num_class);
  NIIK_RET0((!niikmat_set_all(K->wts,1)),fcname,"niikmat_set_all(K->wts,1)");
  K->iter = iter;
  if(verbose>=0) niik_fc_display(fcname,0);
  return K;
} /* niik_classify_kmeans_init */

niik_classify_kmeans *niik_classify_kmeans_free(niik_classify_kmeans *K) {
  if(K==NULL) return NULL;
  if(K->mask!=NULL) {
    K->mask=niik_image_free(K->mask);
  }
  if(K->img !=NULL) {
    K->img =niik_image_free(K->img);
  }
  if(K->mean!=NULL) {
    K->mean=niikmat_free(K->mean);
  }
  free(K);
  return NULL;
} /* niik_classify_kmeans_free */

int niik_classify_kmeans_display(niik_classify_kmeans *K) {
  char fcname[64]="niik_classify_kmeans_display";
  int n,m;
  if(K==NULL) {
    fprintf(stderr,"[%s] ERROR: K is null\n",fcname);
    return 0;
  }
  niik_fc_display(fcname,1);
  fprintf(stdout,"  mask          :   %s\n",K->mask->fname);
  fprintf(stdout,"  img           :   %s\n",K->img->fname);
  fprintf(stdout,"  iter          :   %i\n",K->iter);
  fprintf(stdout,"  # img         :   %i\n",K->img->nu);
  fprintf(stdout,"  # class       :   %i\n",K->num_class);
  for(n=0; n<K->img->nu; n++) {
    fprintf(stdout,"  weights %3i    :   ",n);
    for(m=0; m<K->num_class; m++) {
      fprintf(stdout,"%12.9f ",K->wts->m[n][m]);
    }
    fprintf(stdout,"\n");
  }
  for(n=0; n<K->img->nu; n++) {
    fprintf(stdout,"  mean %3i      :   ",n);
    for(m=0; m<K->num_class; m++) {
      fprintf(stdout,"%12.9f ",K->mean->m[n][m]);
    }
    fprintf(stdout,"\n");
  }
  return 1;
} /* niik_classify_kmeans_display */


int niik_classify_kmeans_initialize_means(niik_classify_kmeans *K,int verbose) {
  char fcname[64]="niik_classify_kmeans_initialize_means";
  int ni,nc,i,j,nvox,nm;
  unsigned char *bimg=NULL;
  double
  mean,stdv,
       *dimg=NULL;
  if(K==NULL) {
    fprintf(stderr,"[%s] ERROR: K is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  bimg=K->mask->data;
  dimg=K->img->data;
  nvox=K->img->nx*K->img->ny*K->img->nz;
  nm=niik_image_count_mask(K->mask);
  for(ni=0; ni<K->img->nu; ni++) {
    K->mean->m[ni][0]=0;
    stdv = 0;
    for(i=0,j=ni*nvox; i<nvox; i++,j++) {
      if(bimg[i]==0) continue;
      K->mean->m[ni][0] += dimg[j];
      stdv += dimg[j]*dimg[j];
    }
    K->mean->m[ni][0] /= nm;
    mean = K->mean->m[ni][0];
    stdv = sqrt(stdv/nm - K->mean->m[ni][0]);
    for(nc=0; nc<K->num_class; nc++) {
      K->mean->m[ni][nc] = (nc-K->num_class/2.0)*stdv*0.1 + mean;
    }
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_initialize_means */


int niik_classify_kmeans_assign_classes(niik_classify_kmeans *K,int verbose) {
  char fcname[64]="niik_classify_kmeans_assign_classes";
  int
  ni,nc,
  i,nvox;
  unsigned char *bimg=NULL;
  double
  *dlist,
  *dimg=NULL;
  if(K==NULL) {
    fprintf(stderr,"[%s] ERROR: K is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  nvox=K->img->nx*K->img->ny*K->img->nz;
  bimg=K->mask->data;
  dimg=K->img->data;
  dlist=(double *)calloc(K->num_class,sizeof(double));
  for(i=0; i<nvox; i++) {
    if(bimg[i]==0) continue;
    for(nc=0; nc<K->num_class; nc++) {
      dlist[nc]=0;
      for(ni=0; ni<K->img->nu; ni++) {
        dlist[nc]+=K->wts->m[ni][nc]*NIIK_SQ(dimg[i+ni*nvox]-K->mean->m[ni][nc]);
      }
    }
    nc=niik_get_min_index_double_vector(dlist,K->num_class);
    bimg[i]=nc+1;
  }
  free(dlist);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_assign_classes */


int niik_classify_kmeans_update_means(niik_classify_kmeans *K,int verbose) {
  char fcname[64]="niik_classify_kmeans_update_means";
  int
  ni,nc,nn,
  i,j,nvox;
  unsigned char *bimg=NULL;
  double
  *dimg=NULL;
  if(K==NULL) {
    fprintf(stderr,"[%s] ERROR: K is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  nvox=K->img->nx*K->img->ny*K->img->nz;
  bimg=K->mask->data;
  dimg=K->img->data;
  for(nc=0; nc<K->num_class; nc++) {
    nn=nc+1;
    for(ni=0; ni<K->img->nu; ni++) {
      K->mean->m[ni][nc]=1e-6;
      for(i=j=0; i<nvox; i++) {
        if(bimg[i]!=nn) continue;
        K->mean->m[ni][nc]+=dimg[i+ni*nvox];
        j++;
      }
      K->mean->m[ni][nc]/=j;
    } /* each dimension */
  } /* each class */
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_update_means */


int niik_classify_kmeans_iterate(niik_classify_kmeans *K,int verbose) {
  char fcname[64]="niik_classify_kmeans_iterate";
  int iter;
  if(K==NULL) {
    fprintf(stderr,"[%s] ERROR: K is null\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  for(iter=0; iter<K->iter; iter++) {
    fprintf(stdout,"[%s] iteration %i\n",fcname,iter+1);
    niik_classify_kmeans_assign_classes(K,verbose);
    niik_classify_kmeans_update_means(K,verbose);
    /*if(verbose>=1) niik_classify_kmeans_display(K);*/
  }
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_classify_kmeans_iterate */

#endif /*  _FALCON_CLASSIFY_KMEANS_C_ */
