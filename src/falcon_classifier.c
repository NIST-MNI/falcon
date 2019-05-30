/* Filename:     nifti1_kunio_classifier.c
 * Description:  classification functions
 * Author:       Kunio Nakamura
 * Date:         December 28, 2012
 */

#ifndef _FALCON_CLASSIFIER_C_
#define _FALCON_CLASSIFIER_C_

#include "falcon.h"

niik_naive_bayesian_classifier *niik_naive_bayesian_classifier_init(int num_img,int num_class) {
  char fcname[64]="niik_naive_bayesian_classifier_init";
  niik_naive_bayesian_classifier *C;
  if((C=(niik_naive_bayesian_classifier *)calloc(1,sizeof(niik_naive_bayesian_classifier)))==NULL) {
    fprintf(stderr,"[%s] ERROR: calloc\n",fcname);
    return NULL;
  }
  C->num_img = num_img;
  C->num_class = num_class;
  C->mask=NULL;
  C->img  =(nifti_image **)calloc(num_img,sizeof(nifti_image *));
  C->primg=(nifti_image **)calloc(num_class,sizeof(nifti_image *));
  C->mean = niikmat_init(num_img,num_class);
  C->stdv = niikmat_init(num_img,num_class);
  return C;
} /* niik_naive_bayesian_classifier_init */


int niik_naive_bayesian_classifier_display(niik_naive_bayesian_classifier *C) {
  char fcname[64]="niik_naive_bayesian_classifier_display";
  int i,j;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  fprintf(stdout,"[%s] naive bayesian classification\n",fcname);
  fprintf(stdout,"  # img              : %i\n",C->num_img);
  fprintf(stdout,"  # class            : %i\n",C->num_class);
  for(i=0; i<C->num_img; i++) {
    if(C->img[i]==NULL) continue;
    fprintf(stdout,"  image %3i          : %s\n",i,C->img[i]->fname);
  }
  for(i=0; i<C->num_img; i++) {
    fprintf(stdout,"  image %3i          : ",i);
    for(j=0; j<C->num_class; j++) {
      fprintf(stdout,"%9.3f %9.3f  | ",C->mean->m[i][j],C->stdv->m[i][j]);
    }
    fprintf(stdout,"\n");
  }
  return 1;
} /* niik_naive_bayesian_classifier_display */


int niik_naive_bayesian_classifier_alloc_primg(niik_naive_bayesian_classifier *C,nifti_image *refimg) {
  char fcname[64]="niik_naive_bayesian_classifier_alloc_primg";
  int i,n;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  for(i=0; i<C->num_class; i++) {
    if(C->primg[i]!=NULL) {
      if(C->mask!=NULL) {
        if((n=niik_image_cmp_dim(C->mask,C->primg[i]))!=0) {
          fprintf(stderr,"[%s] ERROR: different image dimension, %s\n",fcname,niik_image_dim_string(n));
          return 0;
        }
      }
    } else {
      if((C->primg[i]=niik_image_copy_as_type(refimg,NIFTI_TYPE_FLOAT32))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
        return 1;
      }
    }
  } /* each class */
  for(i=0; i<C->num_class; i++) {
    if(!niik_image_type_convert(C->primg[i],NIFTI_TYPE_FLOAT32)) {
      fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
      return 0;
    }
  }
  return 1;
} /* niik_naive_bayesian_classifier_alloc_primg */



/*
 * niik_naive_bayesian_classifeir
 *   for brain tissue classification
 *
 * prefix = niik_naive_bayesian_classfier_brain_classify
 *
 */

int niik_naive_bayesian_classifier_brain_classify_add_t1w(niik_naive_bayesian_classifier *C,
    nifti_image *t1w) {
  char fcname[64]="niik_naive_bayesian_classifier_brain_classify_add_t1w";
  int n;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: classifier->mask is null\n",fcname);
    return 0;
  }
  for(n=0; n<C->num_img; n++) {
    if(C->img[n]==NULL) break;
  }
  if(n>=C->num_img) {
    fprintf(stderr,"[%s] ERROR: img vector is full\n",fcname);
    return 0;
  }
  C->img[n]=t1w;
  if(C->num_class!=5) {
    fprintf(stderr,"[%s] ERROR: not brain classification\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: brain mask is missing\n",fcname);
    return 0;
  }
  /* UPDATING MEAN:   BG, CSF, GM, WM, LES */
  C->mean->m[n][0]=niik_image_get_percentile(t1w,C->mask,0.02);
  C->mean->m[n][1]=niik_image_get_percentile(t1w,C->mask,0.10);
  C->mean->m[n][2]=niik_image_get_percentile(t1w,C->mask,0.35);
  C->mean->m[n][3]=niik_image_get_percentile(t1w,C->mask,0.85);
  C->mean->m[n][4]=niik_image_get_percentile(t1w,C->mask,0.50);
  C->stdv->m[n][0]=C->mean->m[n][0]/1.5;
  C->stdv->m[n][1]=C->mean->m[n][1]/1.5;
  C->stdv->m[n][2]=(C->mean->m[n][3]-C->mean->m[n][2])/1.2;
  C->stdv->m[n][3]=(C->mean->m[n][3]-C->mean->m[n][2])/1.8;
  C->stdv->m[n][4]=(C->mean->m[n][3]-C->mean->m[n][2])/0.8;
  return 1;
}

int niik_naive_bayesian_classifier_brain_classify_add_t2w(niik_naive_bayesian_classifier *C,
    nifti_image *img) {
  char fcname[64]="niik_naive_bayesian_classifier_brain_classify_add_t2w";
  int n;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: classifier->mask is null\n",fcname);
    return 0;
  }
  for(n=0; n<C->num_img; n++) {
    if(C->img[n]==NULL) break;
  }
  if(n>=C->num_img) {
    fprintf(stderr,"[%s] ERROR: img vector is full\n",fcname);
    return 0;
  }
  C->img[n]=img;
  if(C->num_class!=5) {
    fprintf(stderr,"[%s] ERROR: not brain classification\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: brain mask is missing\n",fcname);
    return 0;
  }
  /* UPDATING MEAN:   BG, CSF, GM, WM, LES */
  C->mean->m[n][0]=niik_image_get_percentile(img,C->mask,0.02);
  C->mean->m[n][1]=niik_image_get_percentile(img,C->mask,0.90);
  C->mean->m[n][2]=niik_image_get_percentile(img,C->mask,0.40);
  C->mean->m[n][3]=niik_image_get_percentile(img,C->mask,0.15);
  C->mean->m[n][4]=niik_image_get_percentile(img,C->mask,0.98);
  C->stdv->m[n][0]=C->mean->m[n][0]/1.5;
  C->stdv->m[n][1]=(C->mean->m[n][1]-C->mean->m[n][2])/1.2;
  C->stdv->m[n][2]=(C->mean->m[n][1]-C->mean->m[n][2])/2.2;
  C->stdv->m[n][3]=(C->mean->m[n][1]-C->mean->m[n][2])/1.5;
  C->stdv->m[n][4]=(C->mean->m[n][4]-C->mean->m[n][2])/2.0;
  return 1;
}

int niik_naive_bayesian_classifier_brain_classify_add_pdw(niik_naive_bayesian_classifier *C,
    nifti_image *img) {
  char fcname[64]="niik_naive_bayesian_classifier_brain_classify_add_pdw";
  int n;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: classifier->mask is null\n",fcname);
    return 0;
  }
  for(n=0; n<C->num_img; n++) {
    if(C->img[n]==NULL) break;
  }
  if(n>=C->num_img) {
    fprintf(stderr,"[%s] ERROR: img vector is full\n",fcname);
    return 0;
  }
  C->img[n]=img;
  if(C->num_class!=5) {
    fprintf(stderr,"[%s] ERROR: not brain classification\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: brain mask is missing\n",fcname);
    return 0;
  }
  /* UPDATING MEAN:   BG, CSF, GM, WM, LES */
  C->mean->m[n][0]=niik_image_get_percentile(img,C->mask,0.02);
  C->mean->m[n][1]=niik_image_get_percentile(img,C->mask,0.70);
  C->mean->m[n][2]=niik_image_get_percentile(img,C->mask,0.60);
  C->mean->m[n][3]=niik_image_get_percentile(img,C->mask,0.15);
  C->mean->m[n][4]=niik_image_get_percentile(img,C->mask,0.98);
  C->stdv->m[n][0]=C->mean->m[n][0]/1.5;
  C->stdv->m[n][1]=(C->mean->m[n][2]-C->mean->m[n][3])/0.6;
  C->stdv->m[n][2]=(C->mean->m[n][2]-C->mean->m[n][3])/1.5;
  C->stdv->m[n][3]=(C->mean->m[n][2]-C->mean->m[n][3])/1.5;
  C->stdv->m[n][4]=(C->mean->m[n][4]-C->mean->m[n][2])/2.0;
  return 1;
}

int niik_naive_bayesian_classifier_brain_classify_add_flr(niik_naive_bayesian_classifier *C,
    nifti_image *img) {
  char fcname[64]="niik_naive_bayesian_classifier_brain_classify_add_flr";
  int n;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: classifier->mask is null\n",fcname);
    return 0;
  }
  for(n=0; n<C->num_img; n++) {
    if(C->img[n]==NULL) break;
  }
  if(n>=C->num_img) {
    fprintf(stderr,"[%s] ERROR: img vector is full\n",fcname);
    return 0;
  }
  C->img[n]=img;
  if(C->num_class!=5) {
    fprintf(stderr,"[%s] ERROR: not brain classification\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: brain mask is missing\n",fcname);
    return 0;
  }
  /* UPDATING MEAN:   BG, CSF, GM, WM, LES */
  C->mean->m[n][0]=niik_image_get_percentile(img,C->mask,0.02);
  C->mean->m[n][1]=niik_image_get_percentile(img,C->mask,0.05);
  C->mean->m[n][2]=niik_image_get_percentile(img,C->mask,0.52);
  C->mean->m[n][3]=niik_image_get_percentile(img,C->mask,0.48);
  C->mean->m[n][4]=niik_image_get_percentile(img,C->mask,0.99);
  C->stdv->m[n][0]=C->mean->m[n][0]/1.5;
  C->stdv->m[n][1]=C->mean->m[n][1]/1.5;
  C->stdv->m[n][2]=(C->mean->m[n][2]-C->mean->m[n][1])/1.2;
  C->stdv->m[n][3]=(C->mean->m[n][2]-C->mean->m[n][1])/1.2;
  C->stdv->m[n][4]=(C->mean->m[n][4]-C->mean->m[n][2])/0.8;
  return 1;
}

int niik_naive_bayesian_classifier_brain_classify_add_images(niik_naive_bayesian_classifier *C,
    nifti_image *t1w,
    nifti_image *t1c,
    nifti_image *t2w,
    nifti_image *pdw,
    nifti_image *flr) {
  char fcname[64]="niik_naive_bayesian_classifier_brain_classify_add_images";
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(t1w!=NULL) {
    if(!niik_naive_bayesian_classifier_brain_classify_add_t1w(C,t1w)) {
      fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_t1w\n",fcname);
      return 0;
    }
  } /* add t1w */
  if(t2w!=NULL) {
    if(!niik_naive_bayesian_classifier_brain_classify_add_t2w(C,t2w)) {
      fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_t2w\n",fcname);
      return 0;
    }
  } /* add t2w */
  if(pdw!=NULL) {
    if(!niik_naive_bayesian_classifier_brain_classify_add_pdw(C,pdw)) {
      fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_pdw\n",fcname);
      return 0;
    }
  } /* add pdw */
  if(flr!=NULL) {
    if(!niik_naive_bayesian_classifier_brain_classify_add_flr(C,flr)) {
      fprintf(stderr,"[%s] ERROR: niik_naive_bayesian_classifier_brain_classify_add_flr\n",fcname);
      return 0;
    }
  } /* add pdw */
  return 1;
}

int niik_naive_bayesian_classifier_brain_classify_calc_probability(niik_naive_bayesian_classifier *C) {
  int m,n,i,verbose=0;
  double pv[88],psum;
  char fcname[128]="niik_naive_bayesian_classifier_brain_classify_calc_probability";
  unsigned char *bimg=NULL;
  niikmat *tmpmat=NULL;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: C->mask is null\n",fcname);
    return 0;
  }
  bimg=niik_image_get_voxels_as_uint8_vector(C->mask);
  tmpmat=niikmat_init(C->num_img,C->num_class);
  if(0) {
    for(m=0; m<C->num_class; m++)
      for(n=0; n<C->num_img; n++)
        tmpmat->m[n][m]=sqrt(2.0*NIIK_PI*C->stdv->m[n][m]*C->stdv->m[n][m]);
  } else {
    niikmat_set_all(tmpmat,1);
  }
  for(i=0; i<C->mask->nvox; i++) {
    if(bimg[i]==0) continue;
    if(verbose>=1) {
      fprintf(stdout,"[%3i %3i %3i]\n",
              i%C->mask->nx,(i/C->mask->nx)%C->mask->ny,i/C->mask->nx/C->mask->ny);
    }
    for(m=0; m<C->num_class; m++) {
      pv[m]=1;
      for(n=0; n<C->num_img; n++) {
        pv[m] *= NIIK_GaussPDF(niik_image_get_voxel(C->img[n],i)-C->mean->m[n][m],
                               C->stdv->m[n][m]) / tmpmat->m[n][m];
        if(verbose>=1) {
          fprintf(stdout,"[%3i %3i %3i | %2i %2i] %9.3f | %9.3f %9.3f  ->  %9.5f   %9.5f\n",
                  i%C->mask->nx,(i/C->mask->nx)%C->mask->ny,i/C->mask->nx/C->mask->ny,n,m,
                  niik_image_get_voxel(C->img[n],i), C->mean->m[n][m], C->stdv->m[n][m],
                  NIIK_GaussPDF(niik_image_get_voxel(C->img[n],i)-C->mean->m[n][m],C->stdv->m[n][m]) / tmpmat->m[n][m],
                  pv[m]);
        }
      }
    }
    for(m=0,psum=1e-6; m<C->num_class; m++) {
      psum+=pv[m];
    }
    for(m=0; m<C->num_class; m++) {
      niik_image_set_voxel(C->primg[m],i,pv[m]);
    }
    if(verbose>=1) {
      fprintf(stdout,"[%3i %3i %3i] %9.3f %9.3f %9.3f %9.3f %9.3f\n",
              i%C->mask->nx,(i/C->mask->nx)%C->mask->ny,i/C->mask->nx/C->mask->ny,
              niik_image_get_voxel(C->primg[0],i),
              niik_image_get_voxel(C->primg[1],i),
              niik_image_get_voxel(C->primg[2],i),
              niik_image_get_voxel(C->primg[3],i),
              niik_image_get_voxel(C->primg[4],i));
      fprintf(stdout,"[%3i %3i %3i] %9.3f %9.3f %9.3f %9.3f %9.3f\n",
              i%C->mask->nx,(i/C->mask->nx)%C->mask->ny,i/C->mask->nx/C->mask->ny,
              pv[0],pv[1],pv[2],pv[3],pv[4]);
    }
    for(m=0; m<C->num_class; m++) {
      niik_image_set_voxel(C->primg[m],i,pv[m]/psum);
    }
    if(verbose>=1) {
      if(i>=414844) {
        fprintf(stdout,"[%3i %3i %3i] \n",
                i%C->mask->nx,(i/C->mask->nx)%C->mask->ny,i/C->mask->nx/C->mask->ny);
        exit(0);
      }
    }
  }
  free(bimg);
  tmpmat=niikmat_free(tmpmat);
  return 1;
} /* niik_naive_bayesian_classifier_brain_classify_calc_probability */


int niik_naive_bayesian_classifier_brain_classify_calc_ML(niik_naive_bayesian_classifier *C) {
  int m,n,nn,i,verbose=0;
  double pv[88],dmax,dval,dsum;
  char fcname[128]="niik_naive_bayesian_classifier_brain_classify_calc_probability";
  unsigned char *bimg=NULL;
  niikmat *tmpmat=NULL;
  if(C==NULL) {
    fprintf(stderr,"[%s] ERROR: C is null\n",fcname);
    return 0;
  }
  if(C->mask==NULL) {
    fprintf(stderr,"[%s] ERROR: C->mask is null\n",fcname);
    return 0;
  }
  bimg=niik_image_get_voxels_as_uint8_vector(C->mask);
  for(i=0; i<C->mask->nvox; i++) {
    if(bimg[i]==0) continue;
    m=0;
    dmax = niik_image_get_voxel(C->primg[0],i);
    for(n=1; n<C->num_class; n++) {
      dval = niik_image_get_voxel(C->primg[n],i);
      if(dmax<dval) {
        m=n;
        dmax=dval;
      }
    }
    bimg[i] = m+1;
  }

  tmpmat=niikmat_init(C->num_img,C->num_class);
  for(n=0; n<C->num_class; n++) {
    nn=n+1;
    for(m=0; m<C->num_img; m++) {
      C->mean->m[m][n]=1e-6;
      C->stdv->m[m][n]=1e-6;
      for(i=0; i<C->mask->nvox; i++) {
        if(bimg[i]==nn) {
          C->mean->m[m][n]+= niik_image_get_voxel(C->img[m],i) * niik_image_get_voxel(C->primg[n],i);
          tmpmat->m[m][n] += niik_image_get_voxel(C->primg[n],i);
        }
      }  /* each voxel */
      C->mean->m[m][n] = C->mean->m[m][n] / tmpmat->m[m][n];
      for(i=0; i<C->mask->nvox; i++) {
        if(bimg[i]==nn) {
          C->stdv->m[m][n]+= NIIK_SQ(niik_image_get_voxel(C->img[m],i)-C->mean->m[m][n]) * niik_image_get_voxel(C->primg[n],i);
        }
      }  /* each voxel */
      C->stdv->m[m][n] = sqrt(C->stdv->m[m][n] / tmpmat->m[m][n]);
    } /* each image */
  } /* each class */
  tmpmat=niikmat_free(tmpmat);

  niik_naive_bayesian_classifier_display(C);
  free(bimg);
  return 1;
} /* niik_naive_bayesian_classifier_brain_classify_calc_probability */



int niik_classifier_run_FCM(niik_naive_bayesian_classifier *C) {

}

#endif

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/