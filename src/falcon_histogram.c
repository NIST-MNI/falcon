/* Filename:     nifti1_kunio_histogram.c
 * Description:  histogram functions
 * Author:       Kunio Nakamura
 * Date:         March 18, 2012
 */

#ifndef _FALCON_HISTOGRAM_C_
#define _FALCON_HISTOGRAM_C_

#include "falcon.h"
#include "falcon_morph.h"

/* global variables for this file */
double *g_niik_fit_gaussian_vec[3];
int     g_niik_fit_gaussian_num=0;
int     g_niik_fit_gaussian_mix_num=0;


niikmat *niik_image_histogram_auto(nifti_image *img,nifti_image *maskimg,int num) {
  char fcname[64]="niik_image_histogram_auto";
  niikmat *out=NULL;
  double *x,*y;
  int verbose=0;
  if(verbose>0) niik_fc_display(fcname,1);
  NIIK_RET0((img==NULL),fcname,"img is null");
  NIIK_RET0(((out=niikmat_init(2,num))==NULL),fcname,"niikmat_init");
  x=out->m[0];
  y=out->m[1];
  x[0]=niik_image_get_min(img,maskimg);
  x[num-1]=niik_image_get_max(img,maskimg);
  NIIK_RET0((!niik_image_histogram(img,maskimg,x,y,num)),fcname,"niik_image_histogram");
  if(verbose>0) niik_fc_display(fcname,0);
  return out;
} /* niik_image_histogram_auto */


int niik_image_histogram(nifti_image *img,nifti_image *maskimg,double *hx,double *hy,int num) {
  double
  dval,dmin,dmax,dran;
  int
  i,n;
  unsigned char
  *bimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if( hx==NULL) {
    fprintf(stderr,"ERROR:  hx is a null pointer\n");
    return 0;
  }
  if( hy==NULL) {
    fprintf(stderr,"ERROR:  hy is a null pointer\n");
    return 0;
  }
  dmin = hx[0];
  dmax = hx[num-1];
  dran = (dmax-dmin)/(num-1.0);
  for(n=0; n<num; n++) hy[n]=0;
  for(n=0; n<num; n++) hx[n]=dmin+dran*n;
  if(maskimg!=NULL) {
    if(img->nvox!=maskimg->nvox) {
      fprintf(stderr,"ERROR: different image nvox %i %i\n",img->nvox,maskimg->nvox);
      fprintf(stderr,"       img  = %i %s \n",img->nvox,img->fname);
      fprintf(stderr,"       mask = %i %s \n",maskimg->nvox,maskimg->fname);
      return 0;
    }
    if((bimg = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
      return 0;
    }
  } else {
    bimg = (unsigned char *)calloc(img->nvox,sizeof(char));
    for(i=0; i<img->nvox; i++) bimg[i]=1;
  }
  for(i=0; i<img->nvox; i++) {
    if(bimg[i]) {
      dval=niik_image_get_voxel(img,i);
      if(niik_check_double_problem(dval)) {
        fprintf(stderr,"ERROR: niik_image_get_voxel\n");
        return 0;
      }
      n=(int)floor((dval-dmin)/dran+0.5);
      /* 2012-03-30, Kunio, corrected */
      if(n<0) n=0;
      else if(n>=num) {
        n=num-1;
      }
      hy[n]+=1.0;
    }
  }
  free(bimg);
  return 1;
}

int niik_image_histogram_limits(nifti_image *img,nifti_image *maskimg,double *hx,double *hy,int num)
/* below hx[0] and above hx[num-1] are not included
 */
{
  double
  dval,dmin,dmax,dran;
  int
  i,n;
  unsigned char
  *bimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if( hx==NULL) {
    fprintf(stderr,"ERROR:  hx is a null pointer\n");
    return 0;
  }
  if( hy==NULL) {
    fprintf(stderr,"ERROR:  hy is a null pointer\n");
    return 0;
  }
  dmin = hx[0];
  dmax = hx[num-1];
  dran = (dmax-dmin)/(num-1.0);
  for(n=0; n<num; n++) hy[n]=0;
  for(n=0; n<num; n++) hx[n]=dmin+dran*n;
  if(maskimg!=NULL) {
    if(img->nvox!=maskimg->nvox) {
      fprintf(stderr,"ERROR: different image nvox %i %i\n",img->nvox,maskimg->nvox);
      return 0;
    }
    if((bimg = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
      fprintf(stderr,"ERROR: niik_image_get_voxels_as_uint8_vector\n");
      return 0;
    }
  } else {
    bimg = (unsigned char *)calloc(img->nvox,sizeof(char));
    for(i=0; i<img->nvox; i++) bimg[i]=1;
  }
  for(i=0; i<img->nvox; i++) {
    if(bimg[i]==0) continue;
    dval=niik_image_get_voxel(img,i);
    if(niik_check_double_problem(dval)) {
      fprintf(stderr,"ERROR: niik_image_get_voxel\n");
      return 0;
    }
    n=(int)floor((dval-dmin)/dran+0.5);
    /*fprintf(stdout,"%9i %9i %9.5f\n",i,n,dval);*/
    if(n<0) continue;
    else if(n>=num) continue;
    hy[n]+=1.0;
  }
  free(bimg);
  return 1;
}

int niik_image_histogram_2d(nifti_image *ximg,nifti_image *yimg,nifti_image *maskimg,double xmin,double dx,double ymin,double dy,niikmat *histo,int verbose) {
  char fcname[64]="niik_image_histogram_2d";
  unsigned char *bimg=NULL;
  int i,n,nn[2],mm[2],nm[2];
  double delta[2],vmin[2],val[2],xx[2],dd[2];
  if(verbose>=1) niik_fc_display(fcname,1);
  if(ximg==NULL) {
    fprintf(stderr,"[%s] ERROR: ximg is null\n",fcname);
    return 0;
  }
  if(yimg==NULL) {
    fprintf(stderr,"[%s] ERROR: yimg is null\n",fcname);
    return 0;
  }
  if(niik_image_cmp_dim(ximg,yimg)!=0) {
    fprintf(stderr,"[%s] ERROR: ximg [%i: %ix%ix%i] and yimg [%i: %ix%ix%i]\n",fcname,
            ximg->ndim,ximg->nx,ximg->ny,ximg->nz,
            yimg->ndim,yimg->nx,yimg->ny,yimg->nz);
    return 0;
  }
  if(maskimg!=NULL) {
    bimg=niik_image_get_voxels_as_uint8_vector(maskimg);
  } else {
    bimg=(unsigned char *)calloc(ximg->nvox,sizeof(char));
    for(i=0; i<ximg->nvox; i++) bimg[i]=1;
  }
  niikmat_clear(histo);
  vmin[0] = xmin;
  delta[0] = dx;
  nm[0]=histo->col;
  vmin[1] = ymin;
  delta[1] = dy;
  nm[1]=histo->row;
  for(i=0; i<ximg->nvox; i++) {
    if(bimg[i]==0) continue;
    val[0] = niik_image_get_voxel(ximg,i);
    val[1] = niik_image_get_voxel(yimg,i);
    for(n=0; n<2; n++) {
      xx[n] = (val[n]-vmin[n])/delta[n];
      nn[n] = (int)floor(xx[n]);
      dd[n] = xx[n]-nn[n];
      mm[n] = nn[n]+1;
      if(nn[n]>=nm[n]) nn[0]=-1;
      if(mm[n]>=nm[n]) mm[0]=-1;
    }
    if(nn[0]>=0) {
      if(nn[1]>=0) {
        histo->m[nn[0]][nn[1]] += (1.0-dd[0]) * (1.0-dd[1]);
      }
      if(mm[1]>=0) {
        histo->m[nn[0]][mm[1]] += (1.0-dd[0]) * dd[1];
      }
    }
    if(mm[0]>=0) {
      if(nn[1]>=0) {
        histo->m[mm[0]][nn[1]] += dd[0] * (1.0-dd[1]);
      }
      if(mm[1]>=0) {
        histo->m[mm[0]][mm[1]] += dd[0] * dd[1];
      }
    }
  } /* each voxel */
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}

int niik_image_fit_gaussian(nifti_image *img,nifti_image *maskimg,int avgnum,double *mean,double *stdv,double *err)
/* -estimates the mean and stdev by fitting a normal distribution to
 *  the image histogram
 * -returns 1 for success or 0 for error
 * -outputs are placed in mean and stdv
 * -error is put in err
 */
{
  int
  maxiter=200,
  i,n,num;
  niikmat *p;
  double
  *hy,*hx,
  xi,yi,
  dval,
  dmin,dmax,dran,
  (* pfn)(),tol;
  int verbose=0;
  /* FILE *fp; */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"ERROR: maskimg is a null pointer\n");
    return 0;
  }
  dmin = niik_image_get_min(img,maskimg);
  dmax = niik_image_get_max(img,maskimg);
  dmin-=20; /* artificial amounts */
  dmax+=20;
  num=(dmax-dmin)/3.0;
  dran=(dmax-dmin)/(num-1.0);
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) image intensity values:  min %f  max  %f  ran %f   hist num %i\n",dmin,dmax,dran,num);
  hx = (double *)calloc(num,sizeof(double));
  hy = (double *)calloc(num,sizeof(double));
  for(i=0; i<num; i++) {
    hx[i]=dmin+dran*i;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) histogram x %8.2f %8.2f %8.2f \n",hx[0],dran,hx[num-1]);
  if(!niik_image_histogram(img,maskimg,hx,hy,num)) {
    fprintf(stderr,"ERROR: niik_image_histogram\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) created histogram\n");
  /* global variable */
  g_niik_fit_gaussian_vec[0] = hx;
  g_niik_fit_gaussian_vec[1] = hy;
  g_niik_fit_gaussian_num = num;
  if(avgnum<0) {
    if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) created histogram\n");
    if(!niik_runavg_double_vector(hy,num,num/20)) {
      fprintf(stderr,"ERROR: niik_runavg_double_vector\n");
      return 0;
    }
  } else {
    if(!niik_runavg_double_vector(hy,num,avgnum)) {
      fprintf(stderr,"ERROR: niik_runavg_double_vector\n");
      return 0;
    }
  }
  if(verbose) {
    fprintf(stdout,"-v (niik_image_fit_gaussian) histogram x/y\n");
    niik_write_double_vector("tmp_hx.txt",hx,num);
    niik_write_double_vector("tmp_hy.txt",hy,num);
  }
  if(!niik_get_mode_bspline_vector(hy,num,&xi,&yi)) {
    fprintf(stderr,"ERROR: niik_get_mode_bspline_vector\n");
    return 0;
  }
  xi=dran*xi+hx[0];
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) histogram mode: %f @ %f\n",yi,xi);
  n=niik_get_max_index_double_vector(hy,num);
  for(i=n; i<num; i++) {
    if(hy[i]<hy[n]/2.0) break;
  }
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) stop position %7.4f (@ %2.0f)\n",hy[i],hx[i]);
  if(verbose) fprintf(stdout,"-v (niik_image_fit_gaussian) starting params: mean=%f stdv=%f peak=%f (@ %2.0f)\n",xi,hx[i]-hx[n],yi,hx[n]);
  p=niikmat_init(4,3);
  p->m[0][0]=xi;
  p->m[0][1]=hx[i]-hx[n];
  p->m[0][2]=yi;
  for(n=1; n<4; n++) {
    p->m[n][0]=p->m[0][0]*(1+0.5*(1==n));
    p->m[n][1]=p->m[0][1]*(1+2.5*(2==n));
    p->m[n][2]=p->m[0][2]+(dran*n*(num/10.0));
  }
  pfn = niik_image_fit_gaussian_obj_func;
  tol = 1e-4;
  if(!niik_nelder_mead(p,3,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
    return 0;
  }
  *mean = p->m[0][0];
  *stdv = p->m[0][1];
  dval  = p->m[0][2];
  *err = tol;
  if(verbose) {
    fprintf(stdout,"   mean = %15.8f \n",*mean);
    fprintf(stdout,"   stdv = %15.8f \n",*stdv);
    fprintf(stdout,"   peak = %15.8f \n",dval);
    fprintf(stdout,"   err  = %15.8f \n",*err);
  }
  /* fp=fopen("test.txt","w");
  for(i=0;i<num;i++){
    fprintf(fp,"%f %f %f\n",
      hx[i],hy[i],
      NIIK_GaussPDF(hx[i]-p->m[0][0],p->m[0][1])*p->m[0][2]);
  }
  fclose(fp); */
  return 1;
}

int niik_fit_gaussian_double_vector(double *hy,double *hx,int num,double *mean,double *stdv,double *peak,double *err) {
  int
  maxiter=200,
  i,n;
  niikmat *p;
  double
  xi,yi,
  dran,
  (* pfn)(),tol;
  int verbose=0;
  char fcname[32]="niik_fit_gaussian_double_vector";
  if(hx==NULL) {
    fprintf(stderr,"[%s] ERROR: hx is a null pointer\n",fcname);
    return 0;
  }
  if(hy==NULL) {
    fprintf(stderr,"[%s] ERROR: hy is a null pointer\n",fcname);
    return 0;
  }
  if(verbose>=1) niik_fc_display(fcname,1);
  /* global variable */
  g_niik_fit_gaussian_vec[0] = hx;
  g_niik_fit_gaussian_vec[1] = hy;
  g_niik_fit_gaussian_num = num;
  /* estimate initial condition */
  if(!niik_get_mode_bspline_vector(hy,num,&xi,&yi)) {
    fprintf(stderr,"ERROR: niik_get_mode_bspline_vector\n");
    return 0;
  }
  dran=(hx[num-1]-hx[0])/(num-1.0);
  xi=dran*xi+hx[0];
  if(verbose>=1) fprintf(stdout,"[%s] histogram mode: %f @ %f\n",fcname,yi,xi);
  n=niik_get_max_index_double_vector(hy,num);
  for(i=n; i<num; i++) {
    if(hy[i]<hy[n]/2.0) break;
  }
  if(verbose>=1) fprintf(stdout,"[%s] stop position %7.4f (@ %2.0f)\n",fcname,hy[i],hx[i]);
  if(verbose>=1) fprintf(stdout,"[%s] starting params: mean=%f stdv=%f peak=%f (@ %2.0f)\n",fcname,xi,hx[i]-hx[n],yi,hx[n]);
  p=niikmat_init(4,3);
  p->m[0][0]=xi;
  p->m[0][1]=hx[i]-hx[n];
  p->m[0][2]=yi;
  for(n=1; n<4; n++) {
    p->m[n][0]=p->m[0][0]*(1+0.5*(1==n));
    p->m[n][1]=p->m[0][1]*(1+2.5*(2==n));
    p->m[n][2]=p->m[0][2]+(dran*n*(num/10.0));
  }
  pfn = niik_image_fit_gaussian_obj_func;
  tol = 1e-4;
  if(!niik_nelder_mead(p,3,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
    return 0;
  }
  *mean = p->m[0][0];
  *stdv = p->m[0][1];
  *peak = p->m[0][2];
  *err = tol;
  if(verbose>=1) {
    fprintf(stdout,"   mean = %15.8f \n",*mean);
    fprintf(stdout,"   stdv = %15.8f \n",*stdv);
    fprintf(stdout,"   peak = %15.8f \n",*peak);
    fprintf(stdout,"   err  = %15.8f \n",*err);
    fprintf(stdout,"[%s] y3=%1.5f*exp(-(x-%1.5f).^2/(2*%1.5f^2));\n",fcname,*peak,*mean,*stdv);
  }
  p=niikmat_free(p);
  if(verbose>=1) niik_fc_display(fcname,0);
  return 1;
}  /* niik_fit_gaussian_double_vector */

double niik_image_fit_gaussian_obj_func(double *v)
/* gaussian fit's objective function */
{
  double mean,stdv,peak,dout;
  int n,num;
  mean = v[0];
  stdv = v[1];
  peak = v[2];
  num  = g_niik_fit_gaussian_num;
  for(n=0,dout=0; n<num; n++) {
    dout += fabs(NIIK_GaussPDF(g_niik_fit_gaussian_vec[0][n]-mean,stdv)*peak-g_niik_fit_gaussian_vec[1][n]);
  }
  return dout / num;
}

double niik_fit_gaussian_mixture_obj_func(double *v)
/* gaussian fit's objective function */
{
  double mean,stdv,peak,d,dout,*p;
  int i,n,num;
  num  = g_niik_fit_gaussian_num;
  for(n=0,dout=0; n<num; n++) {
    for(i=0,p=v,d=0; i<g_niik_fit_gaussian_mix_num; i++,p=p+3) {
      mean = p[0];
      stdv = p[1];
      peak = p[2];
      d += peak*NIIK_GaussPDF(g_niik_fit_gaussian_vec[0][n]-mean,stdv);
    }
    dout += fabs(d-g_niik_fit_gaussian_vec[1][n]);
  }
  /*if(1){
    fprintf(stdout,"\t%9.5f %9.5f %9.5f | %9.5f %9.5f %9.5f | %9.5f %9.5f %9.5f | %2.12f\n",
    v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],dout/num); }*/
  return dout / num;
} /* niik_fit_gaussian_mixture_obj_func */


/**** end of stat functions ****/


double niik_image_histogram_optim_bin_size(nifti_image *img,nifti_image *maskimg,int method) {
  double
  dmin,dmax,dvar,
       q1,q3;
  int num,hnum;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(maskimg==NULL) {
    num = img->nvox;
  } else {
    num=niik_image_count_mask(maskimg);
  }
  switch(method) {
  case 1:
    q1 = niik_image_get_percentile(img,maskimg,0.25);
    q3 = niik_image_get_percentile(img,maskimg,0.75);
    dmin=niik_image_get_min(img,maskimg);
    dmax=niik_image_get_max(img,maskimg);
    fprintf(stdout,"  min,q1,q3,max = %9.4f,%9.4f,%9.4f,%9.4f\n",dmin,q1,q3,dmax);
    hnum = pow(2.0*(q3-q1)*num,1.0/3.0);
    return (dmax-dmin) / hnum;
  case 2:
    dvar=niik_image_get_var(img,maskimg);
    return 3.5 * sqrt(dvar) / pow(num,1.0/3.0);
  default:
    fprintf(stderr,"ERROR: niik_image_histogram_optim_bin_size\n");
    return NIIKMAX;
  }
  return NIIKMAX;
}



int niik_image_histogram_matching_test1(nifti_image *refimg,nifti_image *refmask,
                                        nifti_image *img,nifti_image *imgmask,int num)

{
  nifti_image
  *imgs[2],*masks[2];
  char fcname[64]="niik_image_histogram_matching_test1";
  double
  dsum,
  dx,
  *slope,*inter,*thresh,
  dval,
  *hx[2],*hy[2],
  dmax[3],dmin[3];
  int
  i,j,ni[2],n,hnum,
  verbose=2;
  FILE *fp=NULL;

  if(verbose>0) {
    fprintf(stdout,"[%s] start\n",fcname);
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is null\n");
    return 0;
  }
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }

  imgs[0]=refimg;
  imgs[1]=img;
  masks[0]=refmask;
  masks[1]=imgmask;

  /* extreme values */
  for(n=0; n<2; n++) {
    dmin[n] = niik_image_get_min(imgs[n],masks[n]);
    dmax[n] = niik_image_get_max(imgs[n],masks[n]);
    if(verbose>1) fprintf(stdout,"[%s] %i  %12.3f %12.3f\n",fcname,n,dmin[n],dmax[n]);
  }

  /* HISTOGARMS */
  dmin[2]=NIIK_DMIN(dmin[0],dmin[1]);
  dmax[2]=NIIK_DMAX(dmax[0],dmax[1]);
  hnum=200; /*(dmax[2]-dmin[2])/5.0;*/
  if(verbose>1) fprintf(stdout,"[%s] hnum = %i\n",fcname,hnum);
  for(n=0; n<2; n++) {
    if(verbose>2) fprintf(stdout,"[%s] histogram %i\n",fcname,n);
    hx[n]=(double *)calloc(hnum,sizeof(double));
    hy[n]=(double *)calloc(hnum,sizeof(double));
    hx[n][0]     =NIIK_DMIN(dmin[0],dmin[1]);
    hx[n][hnum-1]=NIIK_DMAX(dmax[0],dmax[1]);
    if(!niik_image_histogram(imgs[n],masks[n],hx[n],hy[n],hnum)) {
      fprintf(stderr,"ERROR: niik_image_histogram\n");
      return 0;
    }
  } /* histogram */

  if(verbose>1) fprintf(stdout,"[%s] cumulative probability\n",fcname);
  for(n=0; n<2; n++) {
    for(i=0,dsum=0; i<hnum; i++) {
      dsum+=hy[n][i];
    }
    for(i=0; i<hnum; i++) {
      hy[n][i]/=dsum;
    }
    for(i=1; i<hnum; i++) {
      hy[n][i]+=hy[n][i-1];
    }
  }

  /* check if histogram is more or less ok */
  if(hy[0][0]>=0.5) {
    fprintf(stderr,"[%s] ERROR: background (<%5.2f) is too big (ref)\n",fcname,hx[0][1]);
    return 0;
  }
  if(hy[1][0]>=0.5) {
    fprintf(stderr,"[%s] ERROR: background (<%5.2f) is too big (img)\n",fcname,hx[1][1]);
    return 0;
  }

  /*fp=fopen("tmp_histo_matching_img1.txt","w");*/
  if(fp!=NULL) {
    for(i=0; i<hnum; i++) {
      fprintf(fp,"%8.3f %12.3f\n",hx[0][i],hy[0][i]);
    }
    fclose(fp);
    fp=fopen("tmp_histo_matching_img2.txt","w");
    for(i=0; i<hnum; i++) {
      fprintf(fp,"%8.3f %12.3f\n",hx[1][i],hy[1][i]);
    }
    fclose(fp);
  }

  for(i=0; (i<hnum)&&(verbose>2); i++) {
    fprintf(stdout,"%8.3f %12.3f %12.3f\n",hx[0][i],hy[0][i],hy[1][i]);
  }

  if(verbose>0) {
    fprintf(stdout,"[%s]   num=%i\n",fcname,num);
  }
  slope=(double *)calloc(num,sizeof(double));
  inter=(double *)calloc(num,sizeof(double));
  thresh=(double *)calloc(num+1,sizeof(double));
  thresh[0]=hx[0][0];

  if(verbose>1) fprintf(stdout,"idx  percent   X_ref    Y_ref :    X_img    Y_img  ->    slope    inter\n");
  dx=1.0/num;
  ni[0]=ni[1]=0;
  for(j=1; j<num; j++) {
    for(n=0; n<2; n++) {
      dmin[n]=hx[n][ni[n]];
      for(i=0; i<hnum; i++) {
        if(hy[n][i]>dx*j) break;
      }
      ni[n]=i;
    }
    if(j==1) {
      dmin[0]=dmin[1]=hx[0][0];
    }
    slope[j] = (hx[0][ni[0]]-dmin[0]) / (hx[1][ni[1]]-dmin[1]);
    inter[j] = dmin[0] - slope[j]*dmin[1];
    if(verbose>1)
      fprintf(stdout,"%-2i %8.3f %8.3f %8.3f : %8.3f %8.3f  -> %8.3f %8.3f\n",
              j,j*dx,hx[0][ni[0]],hy[0][ni[0]],hx[1][ni[1]],hy[1][ni[1]],
              slope[j],inter[j]);
    thresh[j]=hx[1][ni[1]];
  }

  if(verbose>0) {
    fprintf(stdout,"[%s]   apply\n",fcname);
  }
  for(i=0; i<img->nvox; i++) {
    dval=niik_image_voxel_get(img,i);
    if(dval>thresh[num-1]) dval=dval*slope[num-1]+inter[num-1];
    else if(dval<thresh[0]) dval=dval*slope[0]+inter[0];
    else {
      for(j=1; j<num; j++) {
        if(thresh[j-1]<dval && dval<thresh[j]) break;
      }
      dval=dval*slope[j]+inter[j];
    }
    niik_image_set_voxel(img,i,dval);
  }

  free(thresh);
  free(slope);
  free(inter);
  for(n=0; n<2; n++) {
    free(hx[n]);
    free(hy[n]);
  }
  if(verbose>0) fprintf(stdout,"[%s] finish\n",fcname);
  return 1;
}



/****************************************************************
 *
 * histogram-based tissue characterization
 *
 ****************************************************************/
int niik_image_fit_gaussian_tissues(nifti_image *img,nifti_image *maskimg,double imin,double imax,int num,int avgnum,
                                    double *mean,double *stdv,double *peak,double *err)
/* -img is the input t1-weighted image
 * -maskimg is the brain parenchymal mask
 * -imin imax are the min and max intensities for histogram
 * -num is the histogram #bins, bin size is determined by imin,imax and num
 * -avgnum is the histogram averaging kernel size
 * -mean stdv are the 3-element vectors for GM,WM,&CSF
 * -err is the fit error
 */
{

  char fcname[32]="niik_image_fit_gaussian_tissues";
  nifti_image *maskroi=NULL;
  niikmat *histomat=NULL,*p=NULL;
  double dran,dval,tol=1e-6,(* pfn)();
  int i,j,ndim,maxiter=1e3,
               verbose=0;

  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is null\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }

  niik_fc_display(fcname,1);

  if((maskroi=niik_image_copy(maskimg))==NULL) {
    fprintf(stderr,"[%s] niik_image_copy\n",fcname);
    return 0;
  }
  if(!niik_image_morph_close_brain(maskroi,5.5,5.5)) {
    fprintf(stderr,"[%s] ERROR: niik_image_morph_close_brain\n",fcname);
    return 0;
  }

  if(imin>imax) {
    imin=niik_image_get_min(img,maskroi);
    imax=niik_image_get_max(img,maskroi);
    dran=imax-imin;
    imin-=dran*0.2;
    imax+=dran*0.2;
  }
  dran=(imax-imin)/(num-1);
  histomat=niikmat_init(4,num);
  for(i=0; i<num; i++) {
    histomat->m[0][i]=imin+dran*i;
  }
  fprintf(stdout,"  %9.5f : %9.5f : %9.5f   %i\n",imin,dran,imax,num);

  /*
   * bimodal gaussian fitting
   */
  fprintf(stdout,"[%s] bimodal fitting for brain parenchyma\n",fcname);
  niik_image_histogram(img,  maskimg,histomat->m[0],histomat->m[1],num);
  if(avgnum>0) niik_runavg_double_vector(histomat->m[1],num,avgnum);
  if(!niik_bimodal_fit_double_vector(histomat->m[0],histomat->m[1],num,mean,stdv,peak,&dval,0)) {
    fprintf(stderr,"[%s] ERROR: niik_bimodal_fit_double_vector\n",fcname);
    return 0;
  }
  if(verbose>=1)
    fprintf(stdout,"[%s] bimodal fit: y1=%1.5f*exp(-(x-%1.5f).^2/(2*%1.5f^2)); y2=%1.5f*exp(-(x-%1.5f).^2/(2*%1.5f^2));\n",fcname,
            peak[0],mean[0],stdv[0],
            peak[1],mean[1],stdv[1]);

  /*
   * single gaussian fitting
   */
  fprintf(stdout,"[%s] single fitting for csf\n",fcname);
  if(!niik_image_maskout(maskroi,maskimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_maskout\n",fcname);
    return 0;
  }
  niik_image_histogram(img,maskroi,histomat->m[0],histomat->m[2],num);
  if(avgnum>0) niik_runavg_double_vector(histomat->m[2],num,avgnum);

  if(!niik_fit_gaussian_double_vector(histomat->m[2],histomat->m[0],num,mean+2,stdv+2,peak+2,&dval)) {
    fprintf(stderr,"[%s] ERROR: niik_image_fit_gaussian\n",fcname);
    return 0;
  }

  if(verbose>=1) {
    for(i=0; i<3*0; i++) {
      fprintf(stdout,"  distribution %i:  %9.3f %9.3f %9.3f\n",i,peak[i],mean[i],stdv[i]);
    }
  }

  /*
   * combined fitting
   */
  fprintf(stdout,"[%s] combined fitting\n",fcname);
  if(!niik_image_add_masks(maskroi,maskimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_add_masks\n",fcname);
    return 0;
  }
  niik_image_histogram(img,maskroi,histomat->m[0],histomat->m[3],num);
  if(avgnum>0) niik_runavg_double_vector(histomat->m[3],num,avgnum);
  ndim=9;
  g_niik_fit_gaussian_mix_num=3;
  g_niik_fit_gaussian_vec[0] = histomat->m[0];
  g_niik_fit_gaussian_vec[1] = histomat->m[3];
  p = niikmat_init(ndim+1,ndim);
  for(i=0; i<ndim+1; i++) {
    j=0;
    p->m[i][j++] = mean[0];
    p->m[i][j++] = stdv[0];
    p->m[i][j++] = peak[0];
    p->m[i][j++] = mean[1];
    p->m[i][j++] = stdv[1];
    p->m[i][j++] = peak[1];
    p->m[i][j++] = mean[2];
    p->m[i][j++] = stdv[2];
    p->m[i][j++] = peak[2];
  }
  for(i=1; i<ndim+1; i++) {
    for(j=0; j<ndim; j++) {
      p->m[i][j] *= 0.9+0.2*niik_get_rand();
    }
  }
  /*niikmat_display(p);*/
  pfn = niik_fit_gaussian_mixture_obj_func;
  tol = 1e-4;
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&maxiter)) {
    fprintf(stderr,"[%s] ERROR: nifti_k_nelder_mead\n",fcname);
    return 0;
  }
  j=0;
  mean[0]=p->m[0][j++];
  stdv[0]=p->m[0][j++];
  peak[0]=p->m[0][j++];
  mean[1]=p->m[0][j++];
  stdv[1]=p->m[0][j++];
  peak[1]=p->m[0][j++];
  mean[2]=p->m[0][j++];
  stdv[2]=p->m[0][j++];
  peak[2]=p->m[0][j++];
  p=niikmat_free(p);
  *err=tol;

  for(i=0; i<num; i++) {
    fprintf(stdout,"%9.4f %9.0f %9.0f %9.0f    %9.3f %9.3f %9.3f %9.3f\n",
            histomat->m[0][i],histomat->m[1][i],histomat->m[2][i],histomat->m[1][i]+histomat->m[2][i],
            peak[0]*NIIK_GaussPDF(histomat->m[0][i]-mean[0],stdv[0]),
            peak[1]*NIIK_GaussPDF(histomat->m[0][i]-mean[1],stdv[1]),
            peak[2]*NIIK_GaussPDF(histomat->m[0][i]-mean[2],stdv[2]),
            peak[0]*NIIK_GaussPDF(histomat->m[0][i]-mean[0],stdv[0])+
            peak[1]*NIIK_GaussPDF(histomat->m[0][i]-mean[1],stdv[1])+
            peak[2]*NIIK_GaussPDF(histomat->m[0][i]-mean[2],stdv[2]));
  }

  i=2;
  fprintf(stdout,"  CSF distribution :  %9.3f %9.3f %9.3f\n",peak[i],mean[i],stdv[i]);
  i=0;
  fprintf(stdout,"  GM distribution  :  %9.3f %9.3f %9.3f\n",peak[i],mean[i],stdv[i]);
  i=1;
  fprintf(stdout,"  WM distribution  :  %9.3f %9.3f %9.3f\n",peak[i],mean[i],stdv[i]);

  histomat=niikmat_free(histomat);
  maskroi=niik_image_free(maskroi);
  niik_fc_display(fcname,0);
  return 1;
} /* niik_image_fit_gaussian_tissues */


#endif /* _FALCON_HISTOGRAM_C_ */

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/