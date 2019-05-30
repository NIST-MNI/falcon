/* Filename:     nifti1_kunio_cortex_edge.c
 * Description:  functions for edge-based cortical atrophy measurement
 * Author:       Kunio Nakamura
 * Date:         August 28, 2012
 */

#ifndef _FALCON_CORTEX_EDGE_C_
#define _FALCON_CORTEX_EDGE_C_

#include "falcon.h"

const int nimg=2;

int niik_image_niikcortex_atrophy_edge_measure_shift(niikvec *x,niikmat *y,double *pial_shift,double *white_shift,double csf_mean,double csf_stdv,double gm_mean,double gm_stdv,double wm_mean,double wm_stdv,double pial_mean,double pial_range,double white_mean,double white_range,double xc_thresh,nifti_image *avgimg,niikmat *bsp_matrix,niikvec *max_shift_vec,int voxel,int verbose) {
  niikmat
  **pti=NULL,   /* intensity-based tissue probability */
    *dy=NULL;   /* first derivative of intensity */
  int
  i,j,k,
  m,n,nn,n0,num,num2;
  double
  xi,yi,xp,yp,
  dx,
  dmax,dval,
  /*ymax=0,*/
  inorm=1,
  ydiv[nimg],
  dymean[nimg],dystdv[nimg],dypeak[nimg],dyerr;
  niikvec **v,*t;
  char
  fname[512],
        fcname[258]="niik_image_niikcortex_atrophy_edge_measure_shift";

  if(verbose>=1) niik_fc_display(fcname,1);
  if(x==NULL) {
    fprintf(stderr,"[%s] ERROR: x is null\n",fcname);
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",fcname);
    return 0;
  }
  num=x->num;
  if(num!=y->col) {
    fprintf(stderr,"[%s] ERROR: size is different between x and y, %i %i\n",fcname,num,y->row);
    return 0;
  }

  if(verbose>=4) {
    fprintf(stdout,"[%s] args\n",fcname);
    fprintf(stdout,"[%s] CSF_mean      %12.6f\n",fcname,csf_mean);
    fprintf(stdout,"[%s] CSF_stdv      %12.6f\n",fcname,csf_stdv);
    fprintf(stdout,"[%s] GM_mean       %12.6f\n",fcname,gm_mean);
    fprintf(stdout,"[%s] GM_stdv       %12.6f\n",fcname,gm_stdv);
    fprintf(stdout,"[%s] WM_mean       %12.6f\n",fcname,wm_mean);
    fprintf(stdout,"[%s] WM_stdv       %12.6f\n",fcname,wm_stdv);
    fprintf(stdout,"[%s] Pial_mean     %12.6f\n",fcname,pial_mean);
    fprintf(stdout,"[%s] Pial_range    %12.6f\n",fcname,pial_range);
    fprintf(stdout,"[%s] White_mean    %12.6f\n",fcname,white_mean);
    fprintf(stdout,"[%s] xc_thresh     %12.6f\n",fcname,xc_thresh);
  }


  /* find the center */
  for(n=1,n0=0; n<num; n++) {
    if(fabs(x->v[n])<fabs(x->v[n0])) n0=n;
  }
  if(verbose>=2) fprintf(stdout,"[%s] centre   %7.4f %i\n",fcname,x->v[n0],n0);

  /* intensity-based tissue probabilities */
  pti=(niikmat **)calloc(nimg,sizeof(niikmat *));
  pti[0]=niikmat_init(4,num);
  pti[1]=niikmat_init(4,num);

  t=niikvec_init(num);
  for(n=0; n<num; n++) {
    t->v[n] = exp(-pow(x->v[n]/10.0,14)*2.0);
  }

  if(verbose>=2) fprintf(stdout,"[%s] surface edge stat\n",fcname);
  if(verbose>=1) fprintf(stdout,"[%s] pial surface edge:  %8.4f +/- %8.4f\n",fcname,pial_mean,pial_range);
  if(verbose>=1) fprintf(stdout,"[%s] white surface edge: %8.4f +/- %8.4f\n",fcname,white_mean,white_range);

  inorm=1;
  dx=x->v[1]-x->v[0];
  dy=niikmat_init(nimg,num);

  if(verbose>=2) fprintf(stdout,"[%s] processing image\n",fcname);
  for(m=0; m<nimg; m++) {
    for(n=0; n<num; n++) {
      dy->m[m][n]=y->m[m][n];
    }
  }

  if(verbose>=2) fprintf(stdout,"[%s] central difference\n",fcname);
  for(m=0; m<nimg; m++) {
    if(1) { /* non-decreasing function */
      for(n=n0+1; n<num; n++) {
        if(dy->m[m][n]>gm_mean) break;
      }
      /* fprintf(stdout,"  white %i %f\n",n,x->v[n]); */
      for(; n<num; n++) {
        dy->m[m][n]=(dy->m[m][n-1]>dy->m[m][n])?dy->m[m][n-1]:dy->m[m][n];
      }
      for(n=n0-1; n>=0; n--) {
        if(dy->m[m][n]<gm_mean) break;
      }
      /* fprintf(stdout,"  pial  %i %f\n",n,x->v[n]); */
      for(; n>=0; n--) {
        dy->m[m][n]=(dy->m[m][n+1]<dy->m[m][n])?dy->m[m][n+1]:dy->m[m][n];
      }
    }

    if(!niik_central_difference_double_vector(dy->m[m],num)) {
      fprintf(stderr,"[%s] ERROR: niik_central_difference_double_vector\n",fcname);
      continue;
    }
    for(n=0; n<num; n++) {
      dy->m[m][n]/=dx;
    }
    if(verbose>=10) {
      niik_bimodal_fit_double_vector(x->v,dy->m[m],num,dymean,dystdv,dypeak,&dyerr,0);
      fprintf(stdout,"[%s]   distr1 %8.4f %8.4f %8.4f\n",fcname,dypeak[0],dymean[0],dystdv[0]);
      fprintf(stdout,"[%s]   distr2 %8.4f %8.4f %8.4f\n",fcname,dypeak[1],dymean[1],dystdv[1]);
    }
  }

  if(verbose>=2) fprintf(stdout,"[%s] splitting pial / white areas\n",fcname);
  for(m=0; m<nimg; m++) {
    ydiv[m]=-100;
    /* modulating gradients */

    /* FIND RELEVANT PARTITIONS */
    if(verbose>=2) fprintf(stdout,"[%s] finding 2 relevant partitions\n",fcname);
    for(n=nn=0,dmax=-1e10; n<n0; n++) {
      if(y->m[m][n]<gm_mean) {
        dval=NIIK_GaussPDF(y->m[m][n]-pial_mean,pial_range);
        if(dy->m[m][n]*dval>dmax) {
          dmax=dy->m[m][n]*dval;
          nn=n;
        }
      }
    }
    n=nn;
    if(verbose>=2) fprintf(stdout,"[%s]   search from: %i %6.2f %6.1f   %6.4f\n",fcname,n,x->v[n],dy->m[m][n],NIIK_GaussPDF(y->m[m][n]-pial_mean,pial_range));

    if(n<=0) continue;
    for(n=n; n<num; n++) {
      if(x->v[n]>3.5) break;
      if(dy->m[m][n]>150) continue;
      if(dy->m[m][n]<=0) break;
      if(dy->m[m][n-1]>dy->m[m][n])
        if(dy->m[m][n+1]>dy->m[m][n])
          break;
    }
    if(x->v[n]>3.5) {
      if(verbose>=2) fprintf(stdout,"[%s]   partition position: %i %6.2f %6.1f   NOT FOUND!\n",fcname,n,x->v[n],dy->m[m][n]);
      ydiv[m]=-100;
    } else {
      if(verbose>=2) fprintf(stdout,"[%s]   partition position: %i %6.2f %6.1f\n",fcname,n,x->v[n],dy->m[m][n]);
      ydiv[m]=x->v[n];
    }
  }

  /*ydiv[0]=ydiv[1]=-100;*/

  /* fprintf(stdout,"ydiv %12.5f %12.5f \n",ydiv[0],ydiv[1]);*/
  if     (ydiv[0]<-99 && ydiv[1]>-99)
    ydiv[0]=ydiv[1];
  else if(ydiv[1]<-99 && ydiv[0]>-99)
    ydiv[1]=ydiv[0];
  else if(fabs(ydiv[0]-ydiv[1])>2.0)
    ydiv[0]=ydiv[1]=NIIK_ABS_DMIN(ydiv[0],ydiv[1]);

  if(ydiv[0]<-99 && ydiv[1]<-99) {
    if(verbose>=2) fprintf(stdout,"[%s] intensity-based tissue probability\n",fcname);
    for(m=0; m<nimg; m++) {
      for(n=0; n<num; n++) {
        pti[m]->m[0][n] = NIIK_Heaviside(inorm*y->m[m][n]-csf_mean,csf_stdv) * exp(-pow((inorm*y->m[m][n] - pial_mean)/ pial_range,8.0)/2.0);   /*pial*/
        /*pti[m]->m[0][n] = NIIK_Heaviside(inorm*y->m[m][n]-csf_mean,csf_stdv) * (1.0-NIIK_Heaviside(inorm*y->m[m][n] - pial_mean, pial_range)); *//* pial */
        pti[m]->m[1][n] = NIIK_Heaviside(inorm*y->m[m][n]-white_mean,fabs(wm_mean-white_mean));  /*white*/
      }
    } /* use intensity-based tissue probabilities */
  } /* each image */

  else {
    if(verbose>=2) fprintf(stdout,"[%s] partition-based tissue probability\n",fcname);
    for(m=0; m<nimg; m++) {
      /* calculate the probability based on this partition */
      for(n=0; n<num; n++) {
        if(x->v[n]<ydiv[m]) {
          pti[m]->m[0][n]=1.0;
          pti[m]->m[1][n]=0.0;
        } else if(x->v[n]>ydiv[m]) {
          pti[m]->m[1][n]=1.0;
          pti[m]->m[0][n]=0.0;
        }
      }
    }
  }

  if(verbose>=2) fprintf(stdout,"[%s] modulated gradient\n",fcname);
  for(m=0; m<nimg; m++) {
    for(n=0; n<num; n++) {
      pti[m]->m[2][n] = dy->m[m][n] * pti[m]->m[0][n] * t->v[n];
      pti[m]->m[3][n] = dy->m[m][n] * pti[m]->m[1][n] * t->v[n];
    }
  } /* each time point */

  /* show info */
  if(verbose>=2) {
    for(m=0; m<nimg; m++) {
      /* fprintf(stdout,"[%s] display intermediate variables\n",fcname);*/
      fprintf(stdout,"  t     = [");
      niikvec_display(t);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  x     = [");
      niikvec_display(x);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  y%i   = [",m);
      niik_display_double_vector( y->m[m],num);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  dy%i  = [",m);
      niik_display_double_vector(dy->m[m],num);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  pti%i = [",m);
      niik_display_double_vector(pti[m]->m[0],num);
      fprintf(stdout,";\n");
      fprintf(stdout,"  ");
      niik_display_double_vector(pti[m]->m[1],num);
      fprintf(stdout,"];\n");
      fprintf(stdout,"  dym%i = [",m);
      niik_display_double_vector(pti[m]->m[2],num);
      fprintf(stdout,";\n");
      fprintf(stdout,"  ");
      niik_display_double_vector(pti[m]->m[3],num);
      fprintf(stdout,"];\n");
      dmax=niik_get_max_from_double_vector(y->m[m],num);
      if     (dmax>1000) {
        dmax = floor(dmax/1000.0)*1000.0;
      } else if(dmax>100 ) {
        dmax = floor(dmax/100.0 )*100.0;
      } else if(dmax>10  ) {
        dmax = floor(dmax/10.0  )*10.0;
      }
    }
  } /* show info */


  /* edge detection */
  v=(niikvec **)calloc((1+nimg),sizeof(niikvec *));
  for(m=0; m<nimg; m++) {
    v[m]=niikvec_init(num);
  }
  num2=2*num+1;
  v[nimg]=niikvec_init(num2);

  /* pial surface shift */
  for(n=0; n<num; n++) {
    v[0]->v[n]=(pti[0]->m[2][n]<50)?0:pti[0]->m[2][n];
    v[1]->v[n]=(pti[1]->m[2][n]<50)?0:pti[1]->m[2][n];
  }

  niikvec_cross_correlation(v[0],v[1],v[nimg]);
  if(max_shift_vec->num!=num2) {
    fprintf(stderr,"[%s] ERROR: num2 != max_shift_vec size, %i %i\n",fcname,num2,max_shift_vec->num);
    return 0;
  }
  for(n=0; n<num2; n++) {
    v[nimg]->v[n] *= max_shift_vec->v[n];
  }
  n=niik_get_max_index_double_vector(v[nimg]->v,num2);
  if(v[nimg]->v[n]>xc_thresh) {
    niik_get_mode_bspline_vector_with_A_matrix(v[nimg]->v,num2,&xi,&yi,bsp_matrix);
    xi = x->v[0]*2.0 + xi*dx - dx;
    xp = xi;
    yp = yi;
  } else {
    xp=0;
    yp=0;
  }

  if(verbose>=2) {
    fprintf(stdout,"  xx  = [");
    for(n=0; n<num2; n++) {
      fprintf(stdout,"%5.2f ",x->v[0]*2.0 + n*dx - dx);
    }
    fprintf(stdout,"];\n");
    fprintf(stdout,"  xc1 = [");
    for(n=0; n<num2; n++) {
      fprintf(stdout,"%5.2f ",v[nimg]->v[n]);
    }
    fprintf(stdout,"];\n");
  }

  /* white surface shift */
  for(n=0; n<num; n++) {
    v[0]->v[n]=(pti[0]->m[3][n]<50)?0:pti[0]->m[3][n];
    v[1]->v[n]=(pti[1]->m[3][n]<50)?0:pti[1]->m[3][n];
  }
  niikvec_cross_correlation(v[0],v[1],v[nimg]);
  for(n=0; n<num2; n++) {
    v[nimg]->v[n] *= max_shift_vec->v[n];
  }
  n=niik_get_max_index_double_vector(v[nimg]->v,num2);
  if(v[nimg]->v[n]>xc_thresh) {
    niik_get_mode_bspline_vector(v[nimg]->v,num2,&xi,&yi);
    /* fprintf(stdout,"[%s] peak %12.2f  at  %9.5f   white\n" ,fcname,yi,xi); */
    xi = x->v[0]*2.0 + xi*dx - dx;
  } else {
    xi = 0;
    yi = 0;
  }

  if(verbose>=2) {
    fprintf(stdout,"  xc2 = [");
    for(n=0; n<num2; n++) {
      fprintf(stdout,"%5.2f ",v[nimg]->v[n]);
    }
    fprintf(stdout,"];\n");
    fprintf(stdout,"  px = [%9.5f %9.5f];\n",xp,xi);
    fprintf(stdout,"  py = [%9.5f %9.5f];\n",yp,yi);
    if(avgimg==NULL) {
      sprintf(fname,"unknown image");
      i=j=k=0;
    } else {
      niik_image_get_pt_from_index(avgimg,voxel,&i,&j,&k);
      sprintf(fname,"%s",avgimg->fname);
      niik_underscore_to_space(fname);
    }
    fprintf(stdout,"\nfigure(1); \nsubplot(4,1,1); \nplot(x,y0,x,y1); \nlegend('m0','m12'); \ntitle('%s voxel %i %i %i (shift=%5.4f)','fontsize',14); \nylabel('Intensity Profile','fontsize',14);\nsubplot(4,1,2);\nplot(x,t,x,pti0(1,:),x,pti0(2,:),x,pti1(1,:),x,pti1(2,:));\naxis([%.1f %.1f 0 1.4]);\nlegend('Distance (mm)',  'pial (base)',  'white (base)',  'pial (FU)',  'white (FU)');\nxlabel('Distance along the edge normal (mm)','fontsize',14);\nylabel('Probabilities','fontsize',14);\ntitle('%s voxel %i %i %i','fontsize',14);\nsubplot(4,1,3);\nplot(x,dy0,x,dy1);\nlegend('m0','m12');\nylabel('First Derivative','fontsize',14);\nsubplot(4,1,4);\nplot(x,dym0(1,:),'bo-',x,dym0(2,:),'ro-',x,dym1(1,:),'bx-',x,dym1(2,:),'rx-');\nlegend('m0 -pial','m0 -white','m12 -pial','m12 -white');\nylabel('Modified Derivative','fontsize',14);\nxlabel('Distance along the edge normal (mm)','fontsize',14);\nnx=find(abs(xx)<=8);\nfigure(2);\nplot(xx(nx),xc1(nx),xx(nx),xc2(nx),px,py,'ro');\nlegend('pial','white','peaks');\nxlabel('Distance along the edge normal (mm)','fontsize',14);\nylabel('Cross-correlation','fontsize',14);\ntitle('%s voxel %i %i %i','fontsize',14);\n",fname,i,j,k,xp-xi,x->v[0],x->v[x->num-1],fname,i,j,k,fname,i,j,k);
  }

  if(verbose>=1) {
    fprintf(stdout,"[%s] peak %12.2f  at  %9.5f   pial\n",fcname,yp,xp);
    fprintf(stdout,"[%s] peak %12.2f  at  %9.5f   white\n",fcname,yi,xi);
  }
  if(verbose>=1) fprintf(stdout,"[%s] shift   %9.5f\n",fcname,xp-xi);
  *pial_shift = xp;
  *white_shift = xi;

  for(m=0; m<nimg+1; m++) {
    v[m]=niikvec_free(v[m]);
  }
  free(v);
  dy=niikmat_free(dy);
  t=niikvec_free(t);
  for(m=0; m<nimg; m++)
    pti[m]=niikmat_free(pti[m]);
  free(pti);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_niikcortex_atrophy_edge_measure_shift */



int niik_image_niikcortex_atrophy_edge_voxel(nifti_image *avgimg,nifti_image **bsimg,niikmat **invmatlist,int interp,double csf_mean,double csf_stdv,double gm_mean,double gm_stdv,double wm_mean,double wm_stdv,double pial_mean,double pial_range,double white_mean,double white_range,double xc_thresh,double xmin,double xmax, double xdelta,niikmat *bsp_matrix,niikvec *max_shift_vec,int voxel,double *pial_shift,double *white_shift,int verbose) {
  char fcname[64]="niik_image_niikcortex_atrophy_edge_voxel";
  niikpt pt,ps,pp[nimg],norm[nimg];
  int i,j,k,m;
  niikvec *x=NULL;
  niikmat *y=NULL;
  pt = niik_image_get_pt_from_index(avgimg,voxel,&i,&j,&k);
  if(verbose>=1) fprintf(stdout,"[%s] voxel = [%3i %3i %3i]\n",fcname,i,j,k);
  ps=niikpt_add(pt,niikpt_unit(niik_image_sobel_filter_voxel(avgimg,voxel)));
  if(verbose>=2) {
    ps=niikpt_unit(niik_image_sobel_filter_voxel(avgimg,voxel));
    fprintf(stdout,"[%s] normal = [%.3f %.3f %.3f]\n",fcname,ps.x,ps.y,ps.z);
    ps=niikpt_add(pt,niikpt_unit(niik_image_sobel_filter_voxel(avgimg,voxel)));
  }
  x=niikvec_init_range(xmin,xmax,xdelta);
  y=niikmat_init(nimg,x->num);
  for(m=0; m<nimg; m++) {
    if(verbose>=2) fprintf(stdout,"[%s] interpolation\n",fcname);
    pp[m]   = niikpt_affine_transform(invmatlist[m],pt);
    norm[m] = niikpt_unit(niikpt_sub(niikpt_affine_transform(invmatlist[m],ps),pp[m]));
    if(verbose>=3) fprintf(stdout,"\t%4i    %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f\n",m,
                             pp[m].x / bsimg[m]->dx,  pp[m].y  / bsimg[m]->dy,   pp[m].z / bsimg[m]->dz,
                             norm[m].x,norm[m].y,norm[m].z);
    if(!niik_image_interp_along_normal_double_vector(bsimg[m],interp,pp[m],norm[m],x->v,y->m[m],x->num)) {
      fprintf(stderr,"[%s] ERROR: niik_image_interp_along_normal_double_vector\n",fcname);
      return 0;
    }
    niik_runavg_double_vector(y->m[m],y->col,5);
  } /* get intensities */
  /***************************
   * do the edge shift calculation
   ***************************/
  if(verbose>=2) fprintf(stdout,"[%s] shift calculation\n",fcname);
  if(!niik_image_niikcortex_atrophy_edge_measure_shift(x,y,pial_shift,white_shift,csf_mean,csf_stdv,gm_mean,gm_stdv,wm_mean,wm_stdv,pial_mean,pial_range,white_mean,white_range,
      xc_thresh,avgimg,bsp_matrix,max_shift_vec,voxel,(verbose>=3)*2)) {
    fprintf(stderr,"[%s] ERROR: niik_image_niikcortex_atrophy_edge_measure_shift\n",fcname);
    return 0;
  }
  y=niikmat_free(y);
  x=niikvec_free(x);
  return 1;
} /* niik_image_niikcortex_atrophy_edge_voxel */



int niik_image_niikcortex_atrophy_edge(nifti_image *avgimg,nifti_image **imglist,niikmat **matrixlist,nifti_image *maskimg,niikmat *tstat,double xc_thresh,niikmat *atv,int *ijk,int verbose)
/* MAIN FUNCTION
 * -img is the averaged image to calculate the gradient
 * -imglist is the list of images (baseline, follow-up)
 * -matrixlist is the list of matrices that would transform images (imglist) to average image (img)
 * -maskimg is a binary image describing the edge voxels (nonzero for edge)
 * -tstat is a matrix of tissue statistics tstat->m[CSF|GM|WM][MEAN|STDV]
 * -xc_thresh is the threshold for cross-correlation (not sure what the right value should be... 2012-09-05, Kunio)
 * -verbose the level of commenting (0-2)
 */
{
  char fcname[64]="niik_image_niikcortex_atrophy_edge";
  nifti_image
  **bsimg=NULL;
  int
  interp=NIIK_INTERP_BSPLINE,
  nvox,
  num,i,j,k,nn,m,n;
  unsigned char *bimg=NULL;
  niikvec
  *max_shift_vec;
  niikmat
  **invmatlist=NULL,
    *bsp_matrix=NULL,
     *cm=NULL;
  double
  d[3],
  xmax=8,xdelta=0.2,
  max_shift=3,
  max_shift_range=1,
  pial_mean,pial_range,
  white_mean,white_range;

  if(verbose>=1) niik_fc_display(fcname,1);
  if(avgimg==NULL) {
    fprintf(stderr,"[%s] ERROR: avgimg is null\n",fcname);
    return 0;
  }
  if(imglist==NULL) {
    fprintf(stderr,"[%s] ERROR: imglist is null\n",fcname);
    return 0;
  }
  if(matrixlist==NULL) {
    fprintf(stderr,"[%s] ERROR: matrixlist is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  if(niik_image_cmp_dim(avgimg,maskimg)!=0) {
    fprintf(stderr,"[%s] ERROR: niik_image_cmp_dim %s %s\n",fcname,avgimg->fname,maskimg->fname);
    return 0;
  }
  if(maskimg->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"[%s] ERROR: maskimg is not NIFTI_TYPE_UINT8\n",fcname);
    return 0;
  }

  /* bspline matrix for histogram interpolation */
  num=(xmax*2)/xdelta+1;
  num=num*2+1;
  fprintf(stdout,"[%s] bspline %i\n",fcname,num);
  bsp_matrix=niikmat_init(num,num);
  niik_bspline_update_A(bsp_matrix);

  invmatlist=(niikmat **)calloc(nimg,sizeof(niikmat *));
  for(m=0; m<nimg; m++) {
    if((invmatlist[m] = niikmat_inverse(matrixlist[m]))==NULL) {
      fprintf(stderr,"[%s] ERROR: niikmat_inverse\n",fcname);
      return 0;
    }
  }

  pial_mean   = 0.3*tstat->m[0][0] + 0.7*tstat->m[1][0];
  white_mean  = 0.5*tstat->m[2][0] + 0.5*tstat->m[1][0];
  pial_range  = fabs(tstat->m[0][0] - tstat->m[1][0]) / 2.0;
  white_range = fabs(tstat->m[2][0] - tstat->m[1][0]) / 2.0;

  max_shift_vec=niikvec_init(num);
  for(n=0; n<num; n++) {
    max_shift_vec->v[n]=1.0-NIIK_Heaviside(fabs(n*xdelta - xmax*2 - xdelta)-max_shift,max_shift_range);
  }
  /*niikvec_display(max_shift_vec);*/

  nvox = niik_image_count_mask(maskimg);

  if(ijk[0]>0) {
    fprintf(stdout,"[%s] checking a voxel %i %i %i\n",fcname,ijk[1],ijk[2],ijk[3]);
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] xcor thresh: %12.1f\n",fcname,xc_thresh);
    fprintf(stdout,"[%s] tissue stats: CSF   %8.3f +/- %8.3f\n",fcname,tstat->m[0][0],tstat->m[0][1]);
    fprintf(stdout,"[%s] tissue stats: GM    %8.3f +/- %8.3f\n",fcname,tstat->m[1][0],tstat->m[1][1]);
    fprintf(stdout,"[%s] tissue stats: WM    %8.3f +/- %8.3f\n",fcname,tstat->m[2][0],tstat->m[2][1]);
    fprintf(stdout,"[%s] edge stats:   CSF   %8.3f +/- %8.3f\n",fcname,pial_mean,pial_range);
    fprintf(stdout,"[%s] edge stats:   WM    %8.3f +/- %8.3f\n",fcname,white_mean,white_range);
    if(verbose>=2) fprintf(stdout,"[%s] mask #voxel: %i\n",fcname,nvox);
  }

  if(verbose>=1) fprintf(stdout,"[%s] bspline coefficient\n",fcname);
  bsimg=(nifti_image **)calloc(nimg,sizeof(nifti_image *));
  if(interp==NIIK_INTERP_BSPLINE) {
    #pragma omp parallel for
    for(m=0; m<nimg; m++) {
      if((bsimg[m]=niik_image_copy(imglist[m]))==NULL) {
        fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
        continue;
      }
      if(!niik_image_interpolate_convert_3d_bspline_coeff(bsimg[m])) {
        fprintf(stderr,"[%s] ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n",fcname);
        continue;
      }
    } /* b-spline coefficients */
    for(m=0; m<nimg; m++) {
      if(bsimg[m]==NULL) {
        fprintf(stderr,"[%s] ERROR: bsimg is null %i\n",fcname,m);
        return 0;
      }
    }
  } else {
    for(m=0; m<nimg; m++) {
      bsimg[m]=imglist[m];
    }
  }

  if(verbose>=1) fprintf(stdout,"[%s] main loop\n",fcname);
  bimg=maskimg->data;

  if(ijk[0]>0) {
    n=ijk[1]+ijk[2]*avgimg->nx+ijk[3]*avgimg->nx*avgimg->ny;
    if(!niik_image_niikcortex_atrophy_edge_voxel(avgimg,bsimg,invmatlist,interp,tstat->m[0][0],tstat->m[0][1],tstat->m[1][0],tstat->m[1][1],tstat->m[2][0],tstat->m[2][1],pial_mean,pial_range,white_mean,white_range,xc_thresh,-xmax,xmax,xdelta,bsp_matrix,max_shift_vec,n,&d[0],&d[1],3)) {
      fprintf(stderr,"[%s] ERROR: niik_image_niikcortex_atrophy_edge_voxel\n",fcname);
      return 0;
    }
    d[2]=d[0]-d[1];
    /*fprintf(stdout,"%i [%3i %3i %3i] %8.3f (=%8.3f + %8.3f)\n",n,ijk[1],ijk[2],ijk[3],d[2],d[0],d[1]);*/
    return 0;
  } /* only this voxel */

  for(n=nn=0; n<avgimg->nvox; n++) {
    if(bimg[n]==0) continue;
    niik_image_get_pt_from_index(avgimg,n,&i,&j,&k);
    if(!niik_image_niikcortex_atrophy_edge_voxel(avgimg,bsimg,invmatlist,interp,tstat->m[0][0],tstat->m[0][1],tstat->m[1][0],tstat->m[1][1],tstat->m[2][0],tstat->m[2][1],pial_mean,pial_range,white_mean,white_range,xc_thresh,-xmax,xmax,xdelta,bsp_matrix,max_shift_vec,n,&d[0],&d[1],0)) {
      fprintf(stderr,"[%s] ERROR: niik_image_niikcortex_atrophy_edge_voxel\n",fcname);
      return 0;
    }
    d[2]=d[0]-d[1];
    /* fprintf(stdout,"%i [%3i %3i %3i] %8.3f (=%8.3f + %8.3f)\n",n,i,j,k,d[2],d[0],d[1]); */
    atv->m[1][nn]=d[0];
    atv->m[2][nn]=d[1];
    atv->m[0][nn]=d[2];
    /*if(n==140+137*avgimg->nx+73*avgimg->nx*avgimg->ny){
      fprintf(stdout,"%i [%3i %3i %3i] %8.3f (=%8.3f + %8.3f)\n",n,i,j,k,atv->m[0][nn],atv->m[1][nn],atv->m[2][nn]);
      fprintf(stdout,"%i [%3i %3i %3i] %8.3f (=%8.3f + %8.3f)\n",n,i,j,k,d[2],d[0],d[1]);
      exit(0); }*/
    nn++;
  } /* each voxel */

  fprintf(stdout,"--- pial motion---\n");
  niik_display_stats_for_double_vector(atv->m[1],nn);
  fprintf(stdout,"--- white motion---\n");
  niik_display_stats_for_double_vector(atv->m[2],nn);
  fprintf(stdout,"--- cortical atrophy---\n");
  niik_display_stats_for_double_vector(atv->m[0],nn);

  /*niik_write_double_vector("tmp_cat.txt",at->m[0],nn);
    fprintf(stdout,"[%s] writing tmp.nii.gz\n",fcname);
    niik_image_write("tmp.nii.gz",cimg);*/

  cm=niikmat_free(cm);
  max_shift_vec=niikvec_free(max_shift_vec);
  for(m=0; m<nimg; m++) {
    invmatlist[m]=niikmat_free(invmatlist[m]);
    if(interp==NIIK_INTERP_BSPLINE)
      bsimg[m]=niik_image_free(bsimg[m]);
  }
  free(bsimg);
  free(invmatlist);
  bsp_matrix=niikmat_free(bsp_matrix);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}

#endif /* _FALCON_CORTEX_EDGE_C_ */
