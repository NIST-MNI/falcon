/* FILENAME:     nifti1_kunio_PABIC.c
 * DESCRIPTION:  PABIC: PArametric BIas Correction
 * REFERENCE:    Styner M, Brechbuhler C, Szekely G, Gerig G. 2000. Parametric estimate of intensity inhomogeneities applied to MRI. IEEE Transactions on Medical Imaging. 19(3); 153-65.
 * Author:       Kunio Nakamura
 * Date:         August 24, 2012
 *
 * main function
 *   nifti_image *niik_image_pabic_bias(nifti_image *img,nifti_image *maskimg,int nd,double dist);
 */


#ifndef _FALCON_PABIC_C_
#define _FALCON_PABIC_C_

#include "falcon.h"

int g_niik_image_pabic_bias_nd=0; /* degree */
nifti_image *g_niik_image_pabic_bias=NULL;
nifti_image *g_niik_image_pabic_img=NULL;
nifti_image *g_niik_image_pabic_mask=NULL;
niikmat *g_niik_image_pabic_stats=NULL;

int niik_image_pabic_bias_ndim(int nd);
int niik_image_pabic_bias_img_update(nifti_image *biasimg,nifti_image *maskimg,double *par,int nd);
double niik_image_pabic_valley(double d, double s) {
  d*=d;  /* basic cost function */
  s*=s;
  return d/(d+3*s);
}

double niik_image_pabic_bias_valley(double gray,niikmat *stats)
/* cost basis function for PABIC algorithm
 * -includes both Gray and White matter components
 */
{
  int n,verbose=0;
  double dv=1;
  if(verbose>=1) fprintf(stdout,"[niik_image_pabic_bias_valley] gray %f\n",gray);
  for(n=0; n<stats->row; n++) {
    if(verbose>=1) fprintf(stdout,"[niik_image_pabic_bias_valley] %i %f %f\n",n,dv,niik_image_pabic_valley(gray-stats->m[n][0],stats->m[n][1]));
    dv*=niik_image_pabic_valley(gray-stats->m[n][0],stats->m[n][1]);
  }
  if(verbose>=1)fprintf(stdout,"[niik_image_pabic_bias_valley] fin %f\n",dv);
  return dv;
}

double niik_image_pabic_bias_cost_func(double *p)
/* cost function for PABIC algorithm */
{
  static int iter=0;
  static double gerr=1e10;
  int i;
  double dval,dout=0;
  unsigned char *bimg;
  if(p==NULL) {
    iter=0;
    gerr=1e10;
    return 0;
  }
  if(!niik_image_pabic_bias_img_update(g_niik_image_pabic_bias,g_niik_image_pabic_mask,p,g_niik_image_pabic_bias_nd)) {
    fprintf(stderr,"[niik_image_pabic_bias_cost_func] ERROR: niik_image_pabic_bias_img_updata\n");
    exit(0);
  }
  bimg=g_niik_image_pabic_mask->data;
  for(i=0; i<g_niik_image_pabic_bias->nvox; i++) {
    if(bimg[i]==0) continue;
    dval = niik_image_get_voxel(g_niik_image_pabic_img,i) * niik_image_get_voxel(g_niik_image_pabic_bias,i);
    dout += niik_image_pabic_bias_valley(dval,g_niik_image_pabic_stats);
  } /* each voxel */
  dout/=niik_image_count_mask(g_niik_image_pabic_mask);
  iter++;
  if(gerr>dout) {
    gerr=dout;
    fprintf(stdout,"%9i %12.7f ",iter,dout);
    niik_display_double_vector(p,niik_image_pabic_bias_ndim(g_niik_image_pabic_bias_nd));
  }
  return dout;
} /* niik_image_pabic_bias_cost_func */

int niik_image_pabic_bias_ndim(int nd)
/* calculate the number of parameters from the number of degrees (nd) */
{
  if(nd==0) return 1;
  return 3*nd + niik_triangular_number(nd-2) + niik_image_pabic_bias_ndim(nd-1);
}

int niik_image_pabic_bias_img_update(nifti_image *biasimg,nifti_image *maskimg,double *par,int nd) {
  int
  i,j,k,
  u,v,w,m,n;
  float *fimg=NULL;
  char fcname[64]="niik_image_pabic_bias_img_update";
  double lyz,lz;
  unsigned char *bimg;
  niikvec *lx;
  int verbose=0;
  if(verbose>=1)niik_fc_display(fcname,1);
  if(!niik_image_clear(biasimg)) {
    fprintf(stderr,"[niik_image_pabic_bias_img_update] ERROR: niik_image_clear\n");
    return 0;
  }
  fimg=biasimg->data;
  lx=niikvec_init(biasimg->nx);
  if(maskimg!=NULL) {
    bimg=maskimg->data;
    for(u=m=0; u<=nd; u++) {
      for(i=0; i<biasimg->nx; i++) {
        lx->v[i]=niik_legendre_func(2.0*i/(biasimg->nx-1.0),u);
      }
      for(v=0; v<=nd; v++) {
        for(w=0; w<=nd; w++) {
          if(u+v+w>nd) continue;
          if(verbose>=1) fprintf(stdout,"[niik_image_pabic_bias_img_update] %2i,%2i,%2i | %2i | %8.3f\n",u,v,w,m,par[m]);
          #pragma omp parallel for private(n,i,j,lz,lyz)
          for(k=0; k<biasimg->nz; k++) {
            n=k*biasimg->nx*biasimg->ny;
            lz=niik_legendre_func(2.0*k/(biasimg->nz-1.0),w);
            for(j=0; j<biasimg->ny; j++) {
              lyz=lz*niik_legendre_func(2.0*j/(biasimg->ny-1.0),v);
              for(i=0; i<biasimg->nx; n++,i++) {
                if(bimg[n]==0) continue;
                fimg[n] += par[m] * lx->v[i] * lyz;
              }
            }
          }
          m++;
        }
      }
    }
    lx=niikvec_free(lx);
  } /* with maskimg */
  else {
    for(u=m=0; u<=nd; u++) {
      for(i=0; i<biasimg->nx; i++) {
        lx->v[i]=niik_legendre_func(2.0*i/(biasimg->nx-1.0),u);
      }
      for(v=0; v<=nd; v++) {
        for(w=0; w<=nd; w++) {
          if(u+v+w>nd) continue;
          if(verbose>=1) fprintf(stdout,"[niik_image_pabic_bias_img_update] %2i,%2i,%2i | %2i | %8.3f\n",u,v,w,m,par[m]);
          #pragma omp parallel for private(n,i,j,lz,lyz)
          for(k=0; k<biasimg->nz; k++) {
            n=k*biasimg->nx*biasimg->ny;
            lz=niik_legendre_func(2.0*k/(biasimg->nz-1.0),w);
            for(j=0; j<biasimg->ny; j++) {
              lyz=lz*niik_legendre_func(2.0*j/(biasimg->ny-1.0),v);
              for(i=0; i<biasimg->nx; n++,i++) {
                fimg[n] += par[m] * lx->v[i] * lyz;
              }
            }
          }
          m++;
        }
      }
    }
    lx=niikvec_free(lx);
  } /* without maskimg */
  if(verbose>=1)niik_fc_display(fcname,0);
  return 1;
} /* niik_image_pabic_bias_img_update */


/************************************************************************************
 *
 * main function
 *
 ************************************************************************************/

nifti_image *niik_image_pabic_bias(nifti_image *img,nifti_image *maskimg,int nd,double dist)
/* -img is the image to be corrected, but not corrected here!
 * -maskimg is the mask for img
 * -nd is the # of degrees (2-3)
 * -dist is the sampling distance
 * -returns the bias field image
 */
{
  nifti_image
  *tmpimg=NULL,
   *tmpmask=NULL,
    *gradimg=NULL,
     *biasimg=NULL;
  double
  mean[9],stdv[9],peak[9],err,
       grad_thresh,
       grad_percent = 0.5,
       atol=1e-5,(* pfn)();
  char
  fname[512],
        fcname[64]="niik_image_pabic_bias";
  niikmat
  *FGMM_stats,
  *p;
  int
  i,
  maxiter=99999,
  ndim,
  verbose=1;

  if(verbose>=1) niik_fc_display(fcname,1);
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return NULL;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is null\n");
    return NULL;
  }

  /* find ndegree
   * 0 ->  1 = 1
   * 1 ->  4 = 1 + 3 (x,y,z)
   * 2 -> 10 = 1 + 3 (x,y,z) + 6 (x2,xy,xz,y2,yz,z2)
   * 3 -> 20 = 1 + 3 (x,y,z) + 6 (x2,xy,xz,y2,yz,z2) + 10 (x3,x2y,x2z,xy2,xyz,xz2,y3,y2z,yz2,z3)
   * 4 -> 35 = 1 + 3 (x,y,z) + 6 (x2,xy,xz,y2,yz,z2) + 10 (x3,x2y,x2z,xy2,xyz,xz2,y3,y2z,yz2,z3) +
   *           15 (x4,x3y,x3z,x2y2,x2yz,x2z2,xy3,xy2z,xyz2,xz3,y4,y3z,y2z2,yz3,z4)
   * 5 -> 56 = 1 + 3 + 6 + 10 + 15 + 21
   * 6 -> 84 = 1 + 3 + 6 + 10 + 15 + 21 + 28
   */
  g_niik_image_pabic_bias_nd=nd;
  if(verbose>=1) fprintf(stdout,"[%s] #degree: %i\n",fcname,nd);
  ndim=niik_image_pabic_bias_ndim(nd) * 3;
  if(verbose>=1) fprintf(stdout,"[%s] #parameters: %i\n",fcname,ndim);

  /* ESTIMATE TISSUE INTENSITIES */
  if(!niik_image_bimodal_fit(img,maskimg,0,0,0,-1,mean,stdv,peak,&err)) {
    fprintf(stderr,"ERROR: niik_image_bimodal_fit\n");
    return NULL;
  }
  fprintf(stdout,"[%s] bimodal fit  %12.6f\n",fcname,err);
  fprintf(stdout,"[%s] distr1 %12.6f %12.6f %12.6f\n",fcname,mean[0],stdv[0],peak[0]);
  fprintf(stdout,"[%s] distr2 %12.6f %12.6f %12.6f\n",fcname,mean[1],stdv[1],peak[1]);

  /* ESTIMATE TISSUE INTENSITIES */
  FGMM_stats=niikmat_init(2,3);
  FGMM_stats->m[0][0]=0.65;
  FGMM_stats->m[1][0]=0.35;
  FGMM_stats->m[0][1]=niik_image_get_percentile(img,maskimg,0.35);
  FGMM_stats->m[1][1]=niik_image_get_percentile(img,maskimg,0.85);
  FGMM_stats->m[0][2]=fabs(FGMM_stats->m[1][1] - FGMM_stats->m[0][1]) / 2.0;
  FGMM_stats->m[1][2]=fabs(FGMM_stats->m[1][1] - FGMM_stats->m[0][1]) / 2.0;
  tmpimg=niik_image_copy(maskimg);
  if(!niik_image_classify_FGMM(img,tmpimg,2,FGMM_stats,15)) {
    fprintf(stderr,"[%s] ERROR: niik_image_bimodal_fit\n",fcname);
    return NULL;
  }
  tmpimg=niik_image_free(tmpimg);
  fprintf(stdout,"[%s] FGMM classification\n",fcname);
  fprintf(stdout,"[%s] GM     %12.6f %12.6f %12.6f\n",fcname,FGMM_stats->m[0][0],FGMM_stats->m[0][1],FGMM_stats->m[0][2]);
  fprintf(stdout,"[%s] WM     %12.6f %12.6f %12.6f\n",fcname,FGMM_stats->m[1][0],FGMM_stats->m[1][1],FGMM_stats->m[1][2]);


  /* SOBEL FILTER */
  if(verbose>=1) fprintf(stdout,"[%s] gradient filter\n",fcname);
  if((gradimg = niik_image_gauss_sobel_filter(img,1.0))==NULL) {
    fprintf(stderr,"ERROR: niik_image_gauss_sobel_filter\n");
    return NULL;
  }

  /* MODIFY AND CREATE MASK */
  if(verbose>=1) fprintf(stdout,"[%s] using gradient percentile  %5.3f\n",fcname,grad_percent);
  grad_thresh = niik_image_get_percentile(gradimg,maskimg,grad_percent);
  if(verbose>=1) fprintf(stdout,"[%s] gradient threshold  %8.3f\n",fcname,grad_thresh);

  if((tmpmask = niik_image_copy_as_type(maskimg,NIFTI_TYPE_UINT8))==NULL) {
    fprintf(stderr,"ERROR: niik_image_copy_as_type\n");
    return NULL;
  }
  for(i=0; i<maskimg->nvox; i++) {
    if(niik_image_get_voxel(gradimg,i)>grad_thresh)
      niik_image_set_voxel(tmpmask,i,0);
  }
  gradimg=niik_image_free(gradimg);
  if(verbose>=1) fprintf(stdout,"[%s] mask size %i\n",fcname,niik_image_count_mask(tmpmask));

  if(verbose>=2) {
    sprintf(fname,"tmp_pabic_mask.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,tmpmask);
  }

  /* PREPARE FOR PABIC OPTIMIZATION */
  if(verbose>=2) fprintf(stdout,"[%s] prepare for optimization\n",fcname);
  atol = 1e-6;
  pfn = niik_image_pabic_bias_cost_func;
  p=niikmat_rand(ndim+1,ndim);
  niikmat_kadd(p,-0.5);
  niikmat_kmul(p,0.3);
  for(i=0; i<ndim; i++) p->m[0][i]=0;
  for(i=0; i<ndim+1; i++) p->m[i][0]+=1.0;
  if(verbose>=2) niikmat_display(p);

  g_niik_image_pabic_img=img;
  g_niik_image_pabic_mask=tmpmask;
  g_niik_image_pabic_stats=niikmat_init(2,2);
  g_niik_image_pabic_stats->m[0][0]=mean[0];
  g_niik_image_pabic_stats->m[1][0]=mean[1];
  g_niik_image_pabic_stats->m[0][1]=stdv[0];
  g_niik_image_pabic_stats->m[1][1]=stdv[1];

  g_niik_image_pabic_stats->m[0][0]=FGMM_stats->m[0][1];
  g_niik_image_pabic_stats->m[1][0]=FGMM_stats->m[1][1];
  g_niik_image_pabic_stats->m[0][1]=FGMM_stats->m[0][2];
  g_niik_image_pabic_stats->m[1][1]=FGMM_stats->m[1][2];
  niikmat_display(g_niik_image_pabic_stats);

  if((biasimg=g_niik_image_pabic_bias=niik_image_copy_as_type(img,NIFTI_TYPE_FLOAT32))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy_as_type\n",fcname);
    return NULL;
  }

  /* run optimization */
  if(verbose>=1) fprintf(stdout,"[%s] nelder mead\n",fcname);
  if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"ERROR: niik_nelder_mead\n");
    return 0;
  }

  /* post-processing
   * -display stats
   * -update the bias field for the entire image coverage
   */
  fprintf(stdout,"[%s] optimized %i %12.8f\n",fcname,maxiter,atol);

  niik_image_pabic_bias_img_update(biasimg,NULL,p->m[0],nd);

  if(verbose>=2) {
    sprintf(fname,"tmp_pabic_test.nii.gz");
    fprintf(stdout,"[%s] writing %s\n",fcname,fname);
    niik_image_write(fname,biasimg);
  }

  if(verbose>=2) fprintf(stdout,"[%s] free memory\n",fcname);
  tmpmask=niik_image_free(tmpmask);
  p=niikmat_free(p);
  FGMM_stats=niikmat_free(FGMM_stats);
  if(verbose>=1) niik_fc_display(fcname,0);
  return biasimg;
} /* niik_image_pabic_bias */


#endif /* _FALCON_PABIC_C_ */
