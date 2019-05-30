/* Filename:     nifti1_kunio_bimodal_fit.c
 * Description:  fitting bimodal distribution
 * Author:       Kunio Nakamura
 * Date:         August 20, 2012
 *
 *
 *
 *
 * int niik_image_bimodal_fit(nifti_image *img,nifti_image *mask,double imin,double delta,double imax,int num,double *mean,double *stdv,double *peak,double *err);
 *   -fitting the image histogram
 *
 *
 * int niik_bimodal_fit_double_vector(double *x,double *y,int num,double *mean,double *stdv,double *peak,double *err,int verbose);
 *   -fitting from a vector (x,y)
 *
 *
 * int niik_bimodal_fit(niikvec *x,niikvec *y,double *mean,double *stdv,double *peak,double *err,int verbose);
 *   -fitting function
 *
 *
 * double niik_bimodal_fit_cost(double *par);
 *   -fitting cost function
 *
 *
 */


#ifndef _FALCON_BIMODAL_FIT_C_
#define _FALCON_BIMODAL_FIT_C_

#include "falcon.h"

int niik_bimodal_fit(niikvec *x,niikvec *y,double *mean,double *stdv,double *peak,double *err,int verbose);
double niik_bimodal_fit_cost(double *par);

/* OPTIMIAZATION VARIABLES */

niikvec *g_niik_bimodal_fit_vec[2];
double g_niik_bimodal_fit_mean_ratio_limits[2]= {3,1}; /* 0=limit value, 1=limit range */
double g_niik_bimodal_fit_peak_ratio_limits[2]= {3,1}; /* 0=limit value, 1=limit range */
int g_niik_bimodal_fit_positive_peak=1;  /* 1=enforce peak is a positive number */

int niik_image_bimodal_fit_set_mean_ratio_limits(double value,double range) {
  g_niik_bimodal_fit_mean_ratio_limits[0]=value;
  g_niik_bimodal_fit_mean_ratio_limits[1]=range;
  return 1;
}
double niik_image_bimodal_fit_get_mean_ratio_limits_value() {
  return g_niik_bimodal_fit_mean_ratio_limits[0];
}
double niik_image_bimodal_fit_get_mean_ratio_limits_range() {
  return g_niik_bimodal_fit_mean_ratio_limits[1];
}
int niik_image_bimodal_fit_set_peak_ratio_limits(double value,double range) {
  g_niik_bimodal_fit_peak_ratio_limits[0]=value;
  g_niik_bimodal_fit_peak_ratio_limits[1]=range;
  return 1;
}
double niik_image_bimodal_fit_get_peak_ratio_limits_value() {
  return g_niik_bimodal_fit_peak_ratio_limits[0];
}
double niik_image_bimodal_fit_get_peak_ratio_limits_range() {
  return g_niik_bimodal_fit_peak_ratio_limits[1];
}


int niik_image_bimodal_fit(nifti_image *img,nifti_image *mask,double imin,double delta,double imax,int num,double *mean,double *stdv,double *peak,double *err)
/* bimodal fit using image with mask
 * -if imin>=imax, then imin/imax are determined automatically
 */
{
  char fcname[64]="niik_image_bimodal_fit";
  niikvec *x,*y;
  int verbose=0;
  if(verbose>=1) niik_fc_display(fcname,1);
  /* input-checking */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  if(mask==NULL) {
    fprintf(stderr,"ERROR: mask is null\n");
    return 0;
  }
  /* set initial values if not defined */
  if(num<0) num=200;
  if(niik_check_double_problem(imin) && niik_check_double_problem(imax)) imin=imax=0;
  if(imin>=imax) {
    imin=niik_image_get_min(img,mask);
    imax=niik_image_get_max(img,mask);
    delta=(imax-imin)/(num-1);
  }
  if(niik_check_double_problem(delta)) {
    delta = (imax-imin)/(num-1);
  }
  if(verbose>=1) fprintf(stdout,"[%s]   x range: %5.2f %5.2f %5.2f   %i\n",fcname,imin,delta,imax,num);
  /*x=niikvec_init_range(imin-delta*20,imax+delta*20,delta);*/
  x=niikvec_init_range(imin,imax,delta);
  y=niikvec_init(x->num);
  num=x->num;
  if(verbose>=2) fprintf(stdout,"[%s] make histogram\n",fcname);
  if(!niik_image_histogram(img,mask,x->v,y->v,x->num)) {
    fprintf(stderr,"ERROR: niik_image_histogram\n");
    x=niikvec_free(x);
    y=niikvec_free(y);
    return 0;
  }
  if(verbose>=2) {
    niikvec_write("tmp_bimodal_x.txt",x);
    niikvec_write("tmp_bimodal_y.txt",y);
  }
  if(!niik_bimodal_fit(x,y,mean,stdv,peak,err,verbose)) {
    fprintf(stderr,"ERROR: niik_bimodal_fit\n");
    x=niikvec_free(x);
    y=niikvec_free(y);
    return 0;
  }
  x=niikvec_free(x);
  y=niikvec_free(y);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_image_bimodal_fit */


int niik_bimodal_fit_double_vector(double *x,double *y,int num,double *mean,double *stdv,double *peak,double *err,int verbose) {
  niikvec *vx,*vy;
  vx=(niikvec *)calloc(1,sizeof(niikvec));
  vy=(niikvec *)calloc(1,sizeof(niikvec));
  vx->num=vy->num=num;
  vx->v=x;
  vy->v=y;
  if(!niik_bimodal_fit(vx,vy,mean,stdv,peak,err,verbose)) {
    fprintf(stderr,"[niik_bimodal_fit_double_vector] ERROR: niik_bimodal_fit\n");
    return 0;
  }
  free(vx); /* not niikvec_free */
  free(vy);
  return 1;
} /* niik_bimodal_fit_double_vector */


int niik_bimodal_fit(niikvec *x,niikvec *y,double *mean,double *stdv,double *peak,double *err,int verbose)
/* fits bimodal function by 2 Gaussian distributions
 * -given x and y, calculate parameters mean, stdv, and peak to fit y
 */
{
  char
  fcname[65]="niik_bimodal_fit";
  int
  nlevel=3,
  ndim=6,
  *nseed,
  m,n,nn,
  n10[2],num;
  niikmat *p=NULL,*s=NULL;
  double
  opeak,
  dval,
  delta,xrange,
  *atol,(* pfn)();
  int *maxiter;
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
  if(y->num!=num) {
    fprintf(stderr,"[%s] ERROR: num is different between x %i and y %i\n",fcname,num,y->num);
    return 0;
  }
  /* initialization
   * -find 2 peaks */
  nn=niik_get_max_index_double_vector(y->v,num);
  opeak=y->v[nn];
  if(verbose>=2) fprintf(stdout,"[%s] overal peak     %8.4f at %8.4f\n",fcname,y->v[nn],x->v[nn]);
  /* Full width at 10% of max */
  for(n=nn; n>0; n--) {
    if(y->v[n]<y->v[nn]*0.1) break;
  }
  n10[0]=n+1;
  if(verbose>=2) fprintf(stdout,"[%s] left 10%% max   %8.4f at %8.4f\n",fcname,y->v[n10[0]],x->v[n10[0]]);
  for(n=nn; n<num; n++) {
    if(y->v[n]<y->v[nn]*0.1) break;
  }
  n10[1]=n-1;
  if(verbose>=2) fprintf(stdout,"[%s] right 10%% max  %8.4f at %8.4f\n",fcname,y->v[n10[1]],x->v[n10[1]]);
  /* find the closest full width at 50% max */
  for(n=n10[0]; n<nn; n++) {
    if(y->v[n]>y->v[nn]*0.5) break;
  }
  peak[0]=y->v[n];
  mean[0]=x->v[n];
  stdv[0]=fabs(x->v[n]-x->v[n10[0]]);
  for(n=n10[1]; n>nn; n--) {
    if(y->v[n]>y->v[nn]*0.5) break;
  }
  peak[1]=y->v[n];
  mean[1]=x->v[n];
  stdv[1]=fabs(x->v[n]-x->v[n10[1]]);
  if(fabs(stdv[0])==0) stdv[0]=stdv[1];
  if(fabs(stdv[1])==0) stdv[1]=stdv[0];
  if(verbose>=2) {
    fprintf(stdout,"[%s] initialization\n",fcname);
    fprintf(stdout,"[%s]   distr1 %8.4f %8.4f %8.4f\n",fcname,peak[0],mean[0],stdv[0]);
    fprintf(stdout,"[%s]   distr2 %8.4f %8.4f %8.4f\n",fcname,peak[1],mean[1],stdv[1]);
  }
  /* set variables for optimization */
  g_niik_bimodal_fit_vec[0]=x;
  g_niik_bimodal_fit_vec[1]=y;
  delta=x->v[1]-x->v[0];
  xrange=x->v[x->num-1]-x->v[0];
  pfn = niik_bimodal_fit_cost;
  maxiter = (int *)calloc(nlevel,sizeof(int));
  nseed   = (int *)calloc(nlevel,sizeof(int));
  atol = (double *)calloc(nlevel,sizeof(double));
  maxiter[0]=100;
  maxiter[1]=1000;
  maxiter[2]=10000;
  nseed[0]=850;
  nseed[1]=20;
  nseed[2]=5;
  atol[0]=1e-3;
  atol[1]=1e-4;
  atol[2]=1e-5;
  s=niikmat_rand(nseed[0],ndim);
  /*niikmat_display(s);*/
  /* get realistic initial estiamtes
   * -means and stdvs are randomly calculated (it's really pseudo-random)
   * -peaks are calculated according to the means
   * -if peaks are too small (10% of overall maximum), then both mean and peak are rejected and
   *  a new random number is generated for mean
   */
  m=0;
  s->m[m][0] = peak[0];
  s->m[m][1] = mean[0];
  s->m[m][2] = stdv[0];
  s->m[m][3] = peak[1];
  s->m[m][4] = mean[1];
  s->m[m][5] = stdv[1];
  for(m=1; m<50; m++) {
    s->m[m][0] = niik_pv(s->m[m][0],peak[0]*0.5,peak[0]*2.0);
    s->m[m][1] = niik_pv(s->m[m][1],mean[0],mean[1]);
    s->m[m][2] = niik_pv(s->m[m][2],stdv[0]*0.5,stdv[0]*2.0);
    s->m[m][3] = niik_pv(s->m[m][3],peak[1]*0.5,peak[1]*2.0);
    s->m[m][4] = niik_pv(s->m[m][4],mean[0],mean[1]);
    s->m[m][5] = niik_pv(s->m[m][5],stdv[1]*0.5,stdv[1]*2.0);
  }
  for(m=50; m<nseed[0]; m++) {
    s->m[m][1] = niik_pv(s->m[m][1],x->v[0],x->v[x->num-1]);
    dval = y->v[(int)NIIK_DMINMAX((s->m[m][1]-x->v[0])/delta+0.5,0,x->num-1)];
    s->m[m][0] = niik_pv(s->m[m][0],dval*0.8,dval*1.2);
    while(s->m[m][0]<0.3*opeak) {
      s->m[m][1] = niik_pv(niik_get_rand(),x->v[0],x->v[x->num-1]);
      dval = y->v[(int)NIIK_DMINMAX((s->m[m][1]-x->v[0])/delta+0.5,0,x->num-1)];
      s->m[m][0] = niik_pv(niik_get_rand(),dval*0.8,dval*1.2);
    }
    s->m[m][4] = niik_pv(s->m[m][4],x->v[0],x->v[x->num-1]);
    dval = y->v[(int)NIIK_DMINMAX((s->m[m][4]-x->v[0])/delta+0.5,0,x->num-1)];
    s->m[m][3] = niik_pv(s->m[m][3],dval*0.8,dval*1.2);
    while(s->m[m][3]<0.5*opeak) {
      s->m[m][4] = niik_pv(niik_get_rand(),x->v[0],x->v[x->num-1]);
      dval = y->v[(int)NIIK_DMINMAX((s->m[m][4]-x->v[0])/delta+0.5,0,x->num-1)];
      s->m[m][3] = niik_pv(niik_get_rand(),dval*0.8,dval*1.2);
    }
    s->m[m][2] = niik_pv(s->m[m][2],delta*num*0.5,delta*num*0.1);
    s->m[m][5] = niik_pv(s->m[m][5],delta*num*0.5,delta*num*0.1);
  }
  /*niikmat_display(s);*/
  p=niikmat_init(nlevel,ndim);
  p->m[0][0] = p->m[0][3] = peak[0]*0.5;
  p->m[0][1] = p->m[0][4] = xrange*0.0;
  p->m[0][2] = p->m[0][5] = delta*num*0.5;
  p->m[1][0] = p->m[1][3] = peak[0]*0.2;
  p->m[1][1] = p->m[1][4] = xrange*0.2;
  p->m[1][2] = p->m[1][5] = delta*num*0.2;
  p->m[2][0] = p->m[2][3] = peak[0]*0.1;
  p->m[2][1] = p->m[2][4] = xrange*0.1;
  p->m[2][2] = p->m[2][5] = delta*num*0.1;
  /* multi-level / multi-seed optimization */
  if(!niik_nelder_mead_multi_level(s,p,nseed,nlevel,ndim,atol,NIIK_NELDER_MEAD_COST_RATIO,pfn,maxiter,0)) {
    fprintf(stderr,"[%s] ERROR: niik_nelder_mead_multi_level\n",fcname);
    return 0;
  }
  if(s->m[0][1]>s->m[0][4]) {
    peak[1] = fabs(s->m[0][0]);
    peak[0] = fabs(s->m[0][3]);
    mean[1] = s->m[0][1];
    mean[0] = s->m[0][4];
    stdv[1] = fabs(s->m[0][2]);
    stdv[0] = fabs(s->m[0][5]);
  } else {
    peak[0] = fabs(s->m[0][0]);
    peak[1] = fabs(s->m[0][3]);
    mean[0] = s->m[0][1];
    mean[1] = s->m[0][4];
    stdv[0] = fabs(s->m[0][2]);
    stdv[1] = fabs(s->m[0][5]);
  }
  err[0] = atol[0];
  p=niikmat_free(p);
  s=niikmat_free(s);
  free(nseed);
  free(maxiter);
  free(atol);
  if(verbose>=2) fprintf(stdout," distr   %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                           peak[0],mean[0],stdv[0],peak[1],mean[1],stdv[1]);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
} /* niik_bimodal_fit */


double niik_bimodal_fit_cost(double *par) {
  int n;
  niikvec *x,*y;
  double delta,out=0,sum=0;
  static int iter=0;
  static double E=1e9;
  /* initialization */
  if(par==NULL) {
    E=1e9;
    iter=0;
  }
  /* get the vectors */
  x=g_niik_bimodal_fit_vec[0];
  y=g_niik_bimodal_fit_vec[1];
  delta=x->v[1]-x->v[0];
  if(g_niik_bimodal_fit_positive_peak) {
    par[0]=fabs(par[0]);
    par[3]=fabs(par[3]);
  }
  out+=NIIK_Heaviside(fabs(par[2]/par[5])-g_niik_bimodal_fit_mean_ratio_limits[0],g_niik_bimodal_fit_mean_ratio_limits[1])*1e5;
  out+=NIIK_Heaviside(fabs(par[5]/par[2])-g_niik_bimodal_fit_mean_ratio_limits[0],g_niik_bimodal_fit_mean_ratio_limits[1])*1e5;
  out+=NIIK_Heaviside(fabs(par[0]/par[3])-g_niik_bimodal_fit_peak_ratio_limits[0],g_niik_bimodal_fit_peak_ratio_limits[1])*1e5;
  out+=NIIK_Heaviside(fabs(par[3]/par[0])-g_niik_bimodal_fit_peak_ratio_limits[0],g_niik_bimodal_fit_peak_ratio_limits[1])*1e5;
  out+=(1.0-NIIK_Heaviside(par[1]-x->v[0]-delta*10,delta*5))*1e9;
  out+=(1.0-NIIK_Heaviside(par[4]-x->v[0]-delta*10,delta*5))*1e9;
  out+=NIIK_Heaviside(par[1]-x->v[x->num-1]+delta*10,delta*5)*1e9;
  out+=NIIK_Heaviside(par[4]-x->v[x->num-1]+delta*10,delta*5)*1e9;
  if(0) {
    fprintf(stdout,"mean ratio :   %9.5f -> %9.5f or %9.5f -> %9.5f\n",
            fabs(par[2]/par[5]),NIIK_Heaviside(fabs(par[2]/par[5])-g_niik_bimodal_fit_mean_ratio_limits[0],g_niik_bimodal_fit_mean_ratio_limits[1]),
            fabs(par[5]/par[2]),NIIK_Heaviside(fabs(par[5]/par[2])-g_niik_bimodal_fit_mean_ratio_limits[0],g_niik_bimodal_fit_mean_ratio_limits[1]));
    fprintf(stdout,"peak ratio :   %9.5f -> %9.5f or %9.5f -> %9.5f\n",
            fabs(par[0]/par[3]),NIIK_Heaviside(fabs(par[0]/par[3])-g_niik_bimodal_fit_peak_ratio_limits[0],g_niik_bimodal_fit_peak_ratio_limits[1]),
            fabs(par[3]/par[0]),NIIK_Heaviside(fabs(par[3]/par[0])-g_niik_bimodal_fit_peak_ratio_limits[0],g_niik_bimodal_fit_peak_ratio_limits[1]));
    fprintf(stdout,"position   :   %9.5f -> %9.5f or %9.5f -> %9.5f    |   %9.5f %9.5f %9.5f\n",
            par[1],(1.0-NIIK_Heaviside(par[1]-x->v[0]-delta*10,delta*5)),
            par[4],(1.0-NIIK_Heaviside(par[4]-x->v[0]-delta*10,delta*5)),x->v[0],delta,x->v[0]+delta*10.0);
    fprintf(stdout,"position   :   %9.5f -> %9.5f or %9.5f -> %9.5f\n",
            par[1],NIIK_Heaviside(par[1]-x->v[x->num-1]+delta*10,delta*5),
            par[4],NIIK_Heaviside(par[4]-x->v[x->num-1]+delta*10,delta*5));
    fprintf(stdout,"out = %.3e\n",out);
  }
  /* calculate errors */
  for(n=0; n<x->num; n++) {
    sum += NIIK_SQ( par[0] * NIIK_GaussPDF(x->v[n]-par[1],par[2]) +
                    par[3] * NIIK_GaussPDF(x->v[n]-par[4],par[5]) -
                    y->v[n] );
  }
  sum/=x->num;
  out+=sum;
  /* display if needed */
  if(out<E) {
    E=out;
    if(0)fprintf(stdout,"%6i  distr %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f | %12.5f\n",iter,par[0],par[1],par[2],par[3],par[4],par[5],out);
  }
  iter++;
  return out;
} /* niik_bimodal_fit_cost */


#endif /* _FALCON_BIMODAL_FIT_C_ */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/