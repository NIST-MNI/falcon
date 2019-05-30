/* Filename:     nifti1_kunio_threshold.c
 * Description:  threshold functions
 * Author:       Kunio Nakamura
 * Date:         February 29, 2012
 *
 * niik_image_threshold       to create a mask
 * niik_image_thresh_otsu     otsu threshold method
 * niik_image_thresh_ridler   ridler threshold method
 * niik_image_thresh_kittler  kittler threshold method
 * niik_image_get_upper_threshold   upper threshold calculation from histogram
 */

#ifndef _FALCON_THRESHOLD_C_
#define _FALCON_THRESHOLD_C_

#include "falcon.h"

int niik_image_threshold(nifti_image *img,double thresh) {
  double *dimg;
  int i;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if((dimg = niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_get_voxels_as_double_vector\n");
    return 0;
  }
  if(img->scl_slope!=0.0l) {
    for(i=0; i<img->nvox; i++) {
      dimg[i] = dimg[i] * img->scl_slope + img->scl_inter;
    }
    img->scl_slope = img->scl_inter = 0.0;
  }
  for(i=0; i<img->nvox; i++) {
    dimg[i] = (dimg[i]>=thresh);
  }
  if(!niik_image_set_voxels_from_double_vector(img,dimg)) {
    fprintf(stderr,"ERROR: niik_image_set_voxels_from_double_vector\n");
    free(dimg);
    return 0;
  }
  if(!niik_image_type_convert_scl(img,NIFTI_TYPE_UINT8,0)) {
    fprintf(stderr,"ERROR: niik_image_type_convert_scl\n");
    free(dimg);
    return 0;
  }
  free(dimg);
  return 1;
}

nifti_image *niik_image_threshold_new(nifti_image *img,double thresh)
/* count the number of voxels and returns it */
{
  nifti_image
  *outimg=NULL;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy\n");
    return NULL;
  }
  if(!niik_image_threshold(outimg,thresh)) {
    fprintf(stderr,"ERROR: niik_image_thresh\n");
    return NULL;
  }
  return outimg;
} /* niik_image_thresh_new */

#define NIFTI_K_OTSU_NUM 11

int niik_image_thresh_otsu(nifti_image * img,nifti_image *maskimg,double *thresh)
/* otsu threshold */
{
  int n,m,iter,i,num;
  double
  m1,m2,dm,w1,w2,svec[NIFTI_K_OTSU_NUM],
  fgray,thresh2,dt,tvec[NIFTI_K_OTSU_NUM];
  unsigned char *bimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(maskimg==NULL) {
    num=NIFTI_K_OTSU_NUM-1;
    tvec[  0] = niik_image_get_min(img,NULL)+1;
    tvec[num] = niik_image_get_max(img,NULL)-1;
    for(iter=0; iter<10; iter++) {
      dt = ( tvec[num] - tvec[0] ) / (NIFTI_K_OTSU_NUM-1.0);
      for(n=m=0; n<NIFTI_K_OTSU_NUM; n++) {
        tvec[n] = thresh2 = tvec[0] + dt * n;
        for(i=0,m1=m2=w1=w2=0; i<img->nvox; i++) {
          fgray=niik_image_get_voxel(img,i);
          if( fgray < thresh2 ) {
            m1 += fgray;
            w1 += 1;
          } else {
            m2 += fgray;
            w2 += 1;
          }
        }
        m1 /= w1;
        m2 /= w2;
        w1 = w1 / (w1 + w2);
        w2 = 1.0 - w1;
        dm = m1 - m2;
        svec[n] = w1 * w2 * dm * dm;
        if(svec[m]<svec[n]) {
          m=n;
        }
      }
      if(m==0)        {
        tvec[  0]=tvec[0];
        tvec[num]=tvec[  1];
      } else if(m==num) {
        tvec[num]=tvec[9];
        tvec[num]=tvec[num];
      } else            {
        tvec[0]=tvec[m-1];
        tvec[num]=tvec[m+1];
      }
      if(fabs(tvec[0]-tvec[num])<1e-3)  {
        *thresh = tvec[m];
        return 1;
      }
    }
    *thresh = tvec[m];
    return 1;
  }
  if(niik_image_cmp_dim(img,maskimg)) {
    fprintf(stderr,"ERORR: nifti_k_cmp_nifti_image_dim \n");
    return 0;
  }
  bimg = niik_image_get_voxels_as_uint8_vector(maskimg);
  num=NIFTI_K_OTSU_NUM-1;
  tvec[  0] = niik_image_get_min(img,maskimg)+1;
  tvec[num] = niik_image_get_max(img,maskimg)-1;
  for(iter=0; iter<10; iter++) {
    dt = ( tvec[num] - tvec[0] ) / (NIFTI_K_OTSU_NUM-1.0);
    for(n=m=0; n<NIFTI_K_OTSU_NUM; n++) {
      tvec[n] = thresh2 = tvec[0] + dt * n;
      for(i=0,m1=m2=w1=w2=0; i<img->nvox; i++) {
        if(!bimg[i]) {
          continue;
        }
        fgray=niik_image_get_voxel(img,i);
        if( fgray < thresh2 ) {
          m1 += fgray;
          w1 += 1;
        } else {
          m2 += fgray;
          w2 += 1;
        }
      }
      m1 /= w1;
      m2 /= w2;
      w1 = w1 / (w1 + w2);
      w2 = 1.0 - w1;
      dm = m1 - m2;
      svec[n] = w1 * w2 * dm * dm;
      if(svec[m]<svec[n]) {
        m=n;
      }
    }
    if(m==0)        {
      tvec[  0]=tvec[0];
      tvec[num]=tvec[  1];
    } else if(m==num) {
      tvec[num]=tvec[9];
      tvec[num]=tvec[num];
    } else            {
      tvec[0]=tvec[m-1];
      tvec[num]=tvec[m+1];
    }
    if(fabs(tvec[0]-tvec[num])<1e-3)  {
      *thresh = tvec[m];
      free(bimg);
      return 1;
    }
  }
  *thresh = tvec[m];
  free(bimg);
  return 1;
}

#undef NIFTI_K_OTSU_NUM


int niik_image_thresh_ridler(nifti_image * img,nifti_image *maskimg,double *thresh)
/*
 * niik_image_thresh_ridler
 *
 * -calculates threshold using Ridler's method (similar to it)
 * -at input thresh is the weighting factor ~[0.5,4.0] around 2.0
 * --increasing
 * -at output, thresh becomes the threshold value
 */
{
  double *dimg = NULL;
  unsigned char *bimg = NULL;
  int i,j,no,nb;
  double omean,bmean,fac,cthresh,pthresh;
  int verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if((dimg = niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: nifti_k_get_voxel_values_vector\n");
    return 0;
  }
  fac = *thresh;
  if(verbose) fprintf(stdout,"  fac %7.1f\n",fac);
  cthresh = pthresh = niik_get_mean_from_double_vector(dimg,img->nvox);
  if(maskimg==NULL) {
    for(j=0; j<10; j++) {
      bmean = omean = 0;
      for(i=no=nb=0; i<img->nvox; i++) {
        if(dimg[i]>cthresh) {
          omean+=dimg[i];
          no++;
        } else {
          bmean+=dimg[i];
          nb++;
        }
      }
      bmean/=nb;
      omean/=no;
      cthresh = (omean+fac*bmean) / (1.0+fac);
      if(pthresh==cthresh) {
        break;
      }
      pthresh = cthresh;
    }
  } else {
    bimg = niik_image_get_voxels_as_uint8_vector(maskimg);
    for(j=0; j<10; j++) {
      bmean = omean = 0;
      for(i=no=nb=0; i<img->nvox; i++) {
        if(!bimg[i]) continue;
        else if(dimg[i]>cthresh) {
          omean+=dimg[i];
          no++;
        } else                     {
          bmean+=dimg[i];
          nb++;
        }
      }
      bmean/=nb;
      omean/=no;
      cthresh = (omean+fac*bmean) / (1.0+fac);
      if(pthresh==cthresh) {
        break;
      }
      if(verbose) fprintf(stdout,"  thresh %7.1f   (%7.1f %7.1f)  %i %i\n",cthresh,omean,bmean,nb,no);
      pthresh = cthresh;
    }
  }
  *thresh = pthresh;
  if(dimg!=NULL) free(dimg);
  if(bimg!=NULL) free(bimg);
  return 1;
} /* niik_image_thresh_ridler */

int niik_image_thresh_kittler(nifti_image * img,nifti_image *maskimg,double *thresh) {
  const int N=6;
  double *dimg = NULL;
  unsigned char *bimg = NULL;
  int
  iter,maxiter=10,
       i,j,n,nlo,no,nb;
  double
  tvec[N],svec[N],dt,
       osum,bsum,ossq,bssq,omean,bmean,ostdv,bstdv;
  int verbose=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if((dimg = niik_image_get_voxels_as_double_vector(img))==NULL) {
    fprintf(stderr,"ERROR: nifti_k_get_voxel_values_vector\n");
    return 0;
  }
  if(maskimg!=NULL) {
    bimg = niik_image_get_voxels_as_uint8_vector(maskimg);
  }
  tvec[  0] = dimg[niik_get_min_index_double_vector(dimg,img->nvox)] + 1;
  tvec[N-1] = dimg[niik_get_max_index_double_vector(dimg,img->nvox)] - 1;
  for(j=0; j<N; j++) {
    svec[j]=NIIKMAX;
  }
  for(iter=0; iter<maxiter; iter++) {
    if(verbose) fprintf(stdout,"  iteration %i\n",iter);
    dt = (tvec[N-1]-tvec[0]) / (N-1.0);
    for(n=1; n<N-1; n++) {
      tvec[n]=dt*n+tvec[0];
    }
    if(verbose) fprintf(stdout,"    %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n",tvec[0],tvec[1],tvec[2],tvec[3],tvec[4],tvec[5]);
    for(j=0; j<N; j++) {
      if(niik_check_double_problem(svec[j])) continue;
      osum = bsum = ossq = bssq = omean = bmean = 0;
      for(i=no=nb=0; i<img->nvox; i++) {
        if(maskimg!=NULL) {
          if(!bimg[i]) continue;
        } else if(dimg[i]>tvec[j]) {
          osum+=dimg[i];
          ossq+=dimg[i]*dimg[i];
          no++;
        } else {
          bsum+=dimg[i];
          bssq+=dimg[i]*dimg[i];
          nb++;
        }
      }
      bmean = bsum / nb;
      omean = osum / no;
      if(verbose) fprintf(stdout,"      mean  %6.2f %6.2f\n",bmean,omean);
      bstdv = sqrt(bssq/nb - bmean * bmean);
      ostdv = sqrt(ossq/no - omean * omean);
      if(verbose) fprintf(stdout,"      stdv  %6.2f %6.2f\n",bstdv,ostdv);
      omean = (double)no / (no+nb);
      bmean = (double)nb / (no+nb);
      if(verbose) fprintf(stdout,"      prob  %6.2f %6.2f\n",bmean,omean);
      svec[j] = omean * log(ostdv+1e-6) + bmean * log(bstdv+1e-6) - omean*log(omean+1e-6) - bmean * log(bmean+1e-6);
      if(verbose) fprintf(stdout,"  thresh %6.2f    %9.6f \n",tvec[j],svec[j]);
    } /* j=0-N */
    /* choose the best one */
    nlo = niik_get_min_index_double_vector(svec,N);
    if(iter==maxiter-1) {
      break;
    } else if(nlo==0) {
      tvec[N-1] = tvec[1];
      svec[N-1] = svec[1];
      for(n=1; n<N-1; n++) svec[n]=NIIKMAX;
    } else if(nlo==N-1) {
      tvec[0] = tvec[N-2];
      svec[0] = svec[N-2];
      for(n=1; n<N-1; n++) svec[n]=NIIKMAX;
    } else {
      tvec[0]   = tvec[nlo-1];
      svec[0]   = svec[nlo-1];
      tvec[N-1] = tvec[nlo+1];
      svec[N-1] = svec[nlo+1];
      for(n=1; n<N-1; n++) svec[n]=NIIKMAX;
    }
  } /* iterations */
  *thresh = tvec[nlo];
  if(dimg!=NULL) free(dimg);
  if(bimg!=NULL) free(bimg);
  if(verbose) fprintf(stdout,"thresh = %7.2f \n",*thresh);
  return 1;
} /* niik_image_thresh_kittler */





/* -i wanted to make a function that detected upper threshold
 *  but this has never been successful
 * -this version is incomplete */
double niik_image_get_upper_threshold(nifti_image *img,nifti_image *maskimg,double fac) {
  double *hx,*hy,dmin,dmax,dran,thresh;
  int n,num,nhi;
  FILE *fp;
  int verbose=0;
  char fcname[64]="niik_image_get_upper_threshold";
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return NIIKMAX;
  }
  dmin = niik_image_get_min(img,NULL);
  dmax = niik_image_get_max(img,NULL);
  num  = 180;
  hx=(double *)calloc(num,sizeof(double));
  hy=(double *)calloc(num,sizeof(double));
  dran = (dmax-dmin)/(num-1.0);
  for(n=0; n<num; n++) {
    hx[n]=n*dran+dmin;
  }
  if(!niik_image_histogram(img,maskimg,hx,hy,num)) {
    fprintf(stderr,"ERROR: niik_image_histogram\n");
    return NIIKMAX;
  }
  niik_runavg_double_vector(hy,num,4);
  nhi = niik_get_max_index_double_vector(hy,num);
  nhi = nhi+0;
  if(verbose) fprintf(stdout,"    nhi %i  %6.0f @ %6.1f\n",nhi,hy[nhi],hx[nhi]);
  thresh = hy[nhi] / 2.0;
  for(n=nhi; n<num; n++) {
    if(hy[n]<thresh)
      break;
  }
  if(verbose) fprintf(stdout,"    nmi %i  %6.0f @ %6.1f\n",n,hy[n],hx[n]);
  thresh = hy[nhi] * fac;
  if(verbose) fprintf(stdout,"    histogram thresh %6.0f \n",thresh);
  for(n=nhi; n<num; n++) {
    if(hy[n] < thresh) break;
    if(hy[n] <= 2) break;
  }
  if(n>=num) {
    n=num-1;
  }
  thresh = hx[n];
  if(verbose) fprintf(stdout,"    thresh %6.0f   (%i,%6.1f) \n",thresh,n,hy[n]);
  if(verbose) {
    fprintf(stdout,"writing histo.txt\n");
    fp = fopen("histo.txt","w");
    for(n=0; n<num; n++) fprintf(fp,"%-9.3f %-9.3f\n",hx[n],hy[n]);
    fclose(fp);
  }
  free(hx);
  free(hy);
  return thresh;
}


#endif /* _FALCON_THRESHOLD_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/