/* Filename:     falcon_bla.c
 * Description:  basic linear algebra functions by Kunio
 * Author:       Kunio Nakamura
 * Date:         February 23, 2012
 */

#include "falcon.h"


int svdcmp(double **a, int m, int n, double *w, double **v);
int niik_singular_value_decomposition(niikmat *a,double *w, niikmat *v);

/***************************************************
 *
 * VECTOR FUNCTIONS
 *
 ***************************************************/

niikvec *niikvec_init(int num) {
  niikvec *v;
  if((v=(niikvec *)calloc(1,sizeof(niikvec)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc 1\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((v->v=(double *)calloc(num,sizeof(double)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc 2\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  v->num=num;
  return v;
}

niikvec *niikvec_init_range(double xmin, double xmax, double delta) {
  niikvec *v;
  int n,num;
  double x;
  for(x=xmin,num=0; x<=xmax+delta/2.0; x+=delta) num++;
  if((v=niikvec_init(num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikvec_init\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  for(n=0; n<num; n++) {
    v->v[n]=n*delta+xmin;
  }
  return v;
}

niikvec *niikvec_init_from_ascii(char *str) {
  int num=0,n;
  char **CSlist=NULL;
  niikvec *v=NULL;
  if((CSlist=niik_csstring_to_list(str,&num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_csstring_to_list %s\n",__FILE__,__LINE__,__func__,str);
    return NULL;
  }
  v=niikvec_init(num);
  for(n=0; n<num; n++) {
    v->v[n]=atof(CSlist[n]);
    free(CSlist[n]);
  }
  free(CSlist);
  return v;
}  /* niikvec_init_from_ascii */

niikvec *niikvec_free(niikvec *v) {
  if(v==NULL) return NULL;
  if(v->v!=NULL) free(v->v);
  free(v);
  return NULL;
}

int niikvec_display(niikvec *v) {
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<v->num; n++) {
    fprintf(stdout,"%12.4lf ",v->v[n]);
  }
  fprintf(stdout,"\n");
  return 1;
}

niikvec *niikvec_read(const char *fname) {
  FILE *fp;
  niikvec *v;
  int num;
  double d;
  if(fname==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fname is null\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((fp=fopen(fname,"r"))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: fopen %s\n",__FILE__,__LINE__,__func__,fname);
    return NULL;
  }
  num=0;
  while(!feof(fp)) {
    if((fscanf(fp,"%lf ",&d))!=1) {
      fprintf(stderr,"[%s:%i:%s] ERROR: fscanf \n",__FILE__,__LINE__,__func__);
      return NULL;
    }
    num++;
  }
  rewind(fp);
  v=niikvec_init(num);
  num=0;
  while(!feof(fp)) {
    if((fscanf(fp,"%lf ",&v->v[num++]))!=1) {
      fprintf(stderr,"[%s:%i:%s] ERROR: fscanf \n",__FILE__,__LINE__,__func__);
      niikvec_free(v);
      return NULL;
    }
  }
  return v;
}


niikvec *niikvec_copy(niikvec *src) {
  niikvec *out;
  if(src==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: src is null\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  out=niikvec_init(src->num);
  if(!niikvec_copy_update(src,out)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikvec_copy_update\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  return out;
}

int niikvec_copy_update(niikvec *src,niikvec *dest) {
  int i;
  if(src==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: src is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(dest==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: dest is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(src->num!=dest->num) {
    fprintf(stderr,"[%s:%i:%s] ERROR: wrong size %i %i\n",__FILE__,__LINE__,__func__,src->num,dest->num);
    return 0;
  }
  for(i=0; i<src->num; i++) {
    dest->v[i]=src->v[i];
  }
  return 1;
}

int niikvec_set_all(niikvec *v,double val) {
  int n;
  for(n=0; n<v->num; n++) {
    v->v[n]=val;
  }
  return 0;
}




/*
 *  sort functions
 *
 *  -similar to niik_median_quicksort_double in nifti1_kunio_median.c
 *  -this function sorts the entire list while niik_median_quicksort_double sorts
 *   halfway for the median calculaation
 *  -niikvec_sort is the main function
 *  -niik_vector_sort_func is the recursive internal function that is accessed from
 *   niikvec_sort
 */

int niikvec_sort(niikvec *v) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  niik_vector_sort_func(v->v,0,v->num-1,v->num);
  return 1;
}

double niik_get_percentile_from_double_vector(double *v,int num,double percent) {
  int n;
  double *w,x;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  if((w = niik_calloc_double_vector(num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_calloc_double_vector(num)\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  for(n=0; n<num; n++) {
    w[n]=v[n];
  }
  if(!niik_sort_double_vector(w,num)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_sort_double_vector\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  x=w[(int)floor((num-1)*percent+0.5)];
  free(w);
  return x;
}

int niik_sort_double_vector(double *v,int num) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  niik_vector_sort_func(v,0,num-1,num);
  return 1;
}

void niik_vector_sort_func(double *v,int nlo,int nhi,int num) {
  double x;
  int i,j;
  x = v[nhi];
  i=nlo;
  j=nhi;
  do {
    while(v[i]<x) {
      i++;
    }
    while(v[j]>x) {
      j--;
    }
    if(i<=j) {
      if(v[i]!=v[j]) {
        NIIK_DSWAP(&v[j],&v[i]);
      }
      i++;
      j--;
    }
  } while(i<=j);
  if(nlo<=j) niik_vector_sort_func(v,nlo,j,num);
  if(i<=nhi) niik_vector_sort_func(v,i,nhi,num);
}

double niik_get_median_from_sorted_double_vector(double *v,int num) {
  if(num%2) {
    return v[num/2];
  }
  return (v[num/2] + v[num/2-1])/2.0;
}


/* return index for sorting
 * Kunio Nakamura, January 12, 2014
 */

int niik_sort_double_vector_index(double *v,int *idx,int num)
/* return the index of sorted array
 * v is un-changed
 * -v and idx are num-long arrays
 */
{
  int i,verbose=0;
  double *w=NULL;
  if(verbose>0) niik_fc_display(__func__,1);
  NIIK_RET0((v==NULL),__func__,"v is null");
  NIIK_RET0((idx==NULL),__func__,"idx is null");
  NIIK_RET0(((w=(double *)calloc(num,sizeof(double)))==NULL),__func__,"calloc");
  for(i=0; i<num; i++) {
    idx[i]=i;
    w[i]=v[i];
  }
  niik_vector_sort_index_func(w,idx,0,num-1,num);
  free(w);
  if(verbose>0) niik_fc_display(__func__,0);
  return 1;
}

void niik_vector_sort_index_func(double *v,int *idx,int nlo,int nhi,int num) {
  double x;
  int i,j;
  x = v[nhi];
  i=nlo;
  j=nhi;
  do {
    while(v[i]<x) {
      i++;
    }
    while(v[j]>x) {
      j--;
    }
    if(i<=j) {
      if(v[i]!=v[j]) {
        NIIK_DSWAP(&v[j],&v[i]);
        NIIK_ISWAP(&idx[j],&idx[i]);
      }
      i++;
      j--;
    }
  } while(i<=j);
  if(nlo<=j) niik_vector_sort_index_func(v,idx,nlo,j,num);
  if(i<=nhi) niik_vector_sort_index_func(v,idx,i,nhi,num);
  return;
}



/***************************************************
 *
 * GENERAL DOUBLE VECTOR FUNCTIONS
 *
 ***************************************************/

double *niik_calloc_double_vector(int num) {
  double *v;
  if(num<=0) {
    fprintf(stderr,"[%s:%i:%s] ERROR: num is invalid\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((v=(double *)calloc(num,sizeof(double)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc failure\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  return v;
}

double **niik_calloc_double_matrix(int m,int n) {
  double **v;
  int i;
  if(m<=0) {
    fprintf(stderr,"[%s:%i:%s] ERROR: m is invalid\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if(n<=0) {
    fprintf(stderr,"[%s:%i:%s] ERROR: n is invalid\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((v=(double **)calloc(m,sizeof(double *)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: calloc failure\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  for(i=0; i<m; i++) {
    v[i]=niik_calloc_double_vector(n);
  }
  return v;
}


double niik_get_max_from_double_vector(double *v,int num) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  return v[niik_get_max_index_double_vector(v,num)];
}

double niikvec_get_max(niikvec *v) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  return niik_get_max_from_double_vector(v->v,v->num);
}

int niik_get_max_index_double_vector(double *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]<v[n]) m=n;
  }
  return m;
}

double niik_get_min_from_double_vector(double *v,int num) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  return v[niik_get_min_index_double_vector(v,num)];
}

double niikvec_get_min(niikvec *v) {
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  return niik_get_min_from_double_vector(v->v,v->num);
}

int niik_get_min_index_double_vector(double *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]>v[n]) m=n;
  }
  return m;
}

double niik_get_mean_from_double_vector(double *v,int num) {
  double d;
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  for(n=0,d=0; n<num; n++) {
    d+=v[n];
  }
  return d/num;
}

double niik_get_stdv_from_double_vector(double *v,int num) {
  return sqrt(niik_get_var_from_double_vector(v,num));
}

double niik_get_var_from_double_vector(double *v,int num) {
  double d,d2;
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  for(n=0,d=d2=0; n<num; n++) {
    d+=v[n];
    d2+=v[n]*v[n];
  }
  return sqrt(d2/num - NIIK_SQ(d/num));
}

double niik_get_sum_from_double_vector(double *v,int num) {
  double d;
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  for(n=0,d=0; n<num; n++) {
    d+=v[n];
  }
  return d;
}

int niik_get_moments_from_double_vector(double *v,int num,double *mean,double *stdv,double *skew,double *kurt,double *adev)
/* based on numerical recipe book */
{
  int
  n,
  verbose=0;
  double s,p,ep,var;
  if(verbose>1) fprintf(stdout,"[niik_get_moments_from_double_vector] start\n");
  *mean = *stdv = *skew = *kurt = *adev = var = 0;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(num==0)  {
    fprintf(stderr,"[%s:%i:%s] ERROR: num is zero\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0,*mean=0; n<num; n++) {
    *mean+=v[n];
  }
  *mean = (*mean) / num;
  if(verbose>0) fprintf(stdout,"mean = %f\n",*mean);
  *stdv = *skew = *kurt = *adev = 0;
  ep = s = p = 0;
  for(n=0; n<num; n++) {
    *adev += fabs(s=v[n]-(*mean));
    ep += s;
    var += (p=s*s);
    *skew += (p *= s);
    *kurt += (p *= s);
  }
  *adev = (*adev) / num;
  if(verbose>0) fprintf(stdout,"adev = %f\n",*adev);
  if(num<=1) {
    return 1;
  }
  var = (var-ep*ep/num)/(num-1);
  if(verbose>0) fprintf(stdout,"var = %f\n",var);
  if(var<1e-5) {
    return 1;
  }
  *stdv=sqrt(var);
  if(verbose>0) fprintf(stdout,"stdv = %f\n",*stdv);
  *skew = *skew / (num*(var)*(*stdv));
  if(verbose>0) fprintf(stdout,"skew = %f\n",*skew);
  *kurt=(*kurt)/(num*(var)*(var))-3.0;
  if(verbose>0) fprintf(stdout,"kurt = %f\n",*kurt);
  if(verbose>1) fprintf(stdout,"[niik_get_moments_from_double_vector] finish\n");
  return 1;
}

int niik_get_trimmed_average_from_double_vector(double *v,int num,double percent,double *out)
/* -get trimmed average from double vector
 * -percent is the upper and lower trimming percent (0-0.5)
 */
{
  niikvec *w=NULL;
  double d;
  int m,n,nn,verbose=1;
  if((w=niikvec_init(num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikvec_init\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) w->v[n]=v[n];
  if(!niikvec_sort(w)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikvec_sort\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(verbose>1) {
    fprintf(stdout,"\tall data\n");
    for(n=0; n<num; n++) fprintf(stdout,"%5i %9.4f\n",n,w->v[n]);
  }
  nn=(int)floor(num*percent);
  for(n=nn,d=0,m=0; n<num-nn; n++) {
    d+=w->v[n];
    m++;
  }
  d/=(num-nn-nn);
  if(verbose>0) {
    fprintf(stdout,"\tm = %i, %i\n",m,num-nn-nn);
    if(verbose>1) {
      fprintf(stdout,"\ttrimmed data\n");
      for(n=nn; n<num-nn; n++) fprintf(stdout,"%5i %9.4f\n",n,w->v[n]);
      fprintf(stdout,"\td = %.9f\n",d);
    }
  }
  *out=d;
  return 1;
} /* niik_get_trimmed_average_from_double_vector */

int niik_get_min_index_float_vector(float *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]>v[n]) m=n;
  }
  return m;
}

int niik_get_max_index_float_vector(float *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]<v[n]) m=n;
  }
  return m;
}

int niik_runavg_double_vector(double *v,int num,int anum) {
  double *w;
  int n,m,nlo,nhi,nn;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if((w=niik_calloc_double_vector(num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_calloc_double_vector\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) {
    w[n]=v[n];
  }
  /*fprintf(stdout,"  num %i    anum %i\n",num,anum);*/
  for(n=0; n<num; n++) {
    nlo=NIIK_IMAX(n-anum,0);
    nhi=NIIK_IMIN(n+anum,num-1);
    v[n]=0;
    for(m=nlo,nn=0; m<=nhi; m++,nn++) {
      v[n]+=w[m];
    }
    v[n]/=nn;
  }
  free(w);
  return 1;
}

int niik_set_zero_for_double_vector(double *v,int num) {
  int n;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) v[n]=0;
  return 1;
}

double niik_interp1d_linear_in_double_vector(double *v,int num,double x)
/* -linearly interpolate 1d vector
 * -returns the value in vector v at x where x is between 0 and num-1 */
{
  int n;
  double d;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return NIIKMAX;
  }
  if(x<0) {
    fprintf(stderr,"[%s:%i:%s] ERROR: x(%7.2f) < 0\n",__FILE__,__LINE__,__func__,x);
    return NIIKMAX;
  }
  if(x>num-1) {
    fprintf(stderr,"[%s:%i:%s] ERROR: x(%7.2f) > num-1\n",__FILE__,__LINE__,__func__,x);
    return NIIKMAX;
  }
  n=(int)floor(x);
  if(n==num-1) return v[n];
  d=x-n;
  /*fprintf(stdout,"  %7.5f -> %7.3f x %7.3f + %7.3f x %7.3f   near %i\n",x,(1.0-d),v[n],d,v[n+1],n);*/
  return v[n+1]*d + v[n]*(1.0-d);
} /* niik_interp1d_linear_in_double_vector */

int niik_central_difference_double_vector(double *v,int num) {
  int n,verbose=0;
  double *w=NULL;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(verbose>=2) fprintf(stdout,"[%s:%i:%s]  num %i\n",__FILE__,__LINE__,__func__,num);
  if((w=(double *)calloc(num,sizeof(double)))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR calloc\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(n=0; n<num; n++) {
    w[n]=v[n];
  }
  for(n=1; n<num-1; n++) {
    v[n]=(w[n+1]-w[n-1])/2.0;
  }
  v[0] = w[1] - w[0];
  v[num-1] = w[num-1] - w[num-2];
  free(w);
  if(verbose>=2) fprintf(stdout,"[niik_central_difference_double_vector] finish\n");
  return 1;
}

int niik_display_stats_for_double_vector(double *v,int num) {
  fprintf(stdout,"[niik_display_stats_for_double_vector]\n");
  fprintf(stdout,"  length %-i\n",num);
  fprintf(stdout,"    mean %-15.9f\n",niik_get_mean_from_double_vector(v,num));
  fprintf(stdout,"    stdv %-15.9f\n",niik_get_stdv_from_double_vector(v,num));
  fprintf(stdout,"     min %-15.9f\n",niik_get_min_from_double_vector(v,num));
  fprintf(stdout,"  median %-15.9f\n",niik_median_quicksort_double_untouch(v,num));
  fprintf(stdout,"     max %-15.9f\n",niik_get_max_from_double_vector(v,num));
  return 1;
}

int niikvec_legendre(niikvec *x, niikvec *y,double *param,int ndeg)
/* if ndeg = 0 then y = constant
 * y is updated
 */
{
  int n;
  if(x==NULL) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: x is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",__func__);
    return 0;
  }
  if(param==NULL) {
    fprintf(stderr,"[%s] ERROR: param is null\n",__func__);
    return 0;
  }
  if(x->num!=y->num) {
    fprintf(stderr,"[%s] ERROR: x and y num is different %i %i\n",__func__,x->num,y->num);
    return 0;
  }
  for(n=0; n<x->num; n++) {
    if(!niik_legendre_func_with_param(x->v[n],y->v+n,param,ndeg)) {
      fprintf(stderr,"[%s] ERROR: niik_legendre_func_with_param\n",__func__);
      return 0;
    }
  }
  return 1;
} /* niikvec_legendre */

/***************************************************
 *
 * GENERAL INTEGER VECTOR FUNCTIONS
 *
 ***************************************************/


double niik_sum_int_vector(int *v,int num) {
  int n;
  double dout=0;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return 0;
  }
  for(n=0; n<num; n++) {
    dout+=v[n];
  }
  return dout;
}

int niik_count_nonzero_from_int_vector(int *v,int num) {
  int n,m;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return 0;
  }
  for(n=m=0; n<num; n++) {
    m+=(v[n]!=0);
  }
  return m;
}

int niik_count_positive_int_from_int_vector(int *v,int num) {
  int n,m;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return 0;
  }
  for(n=m=0; n<num; n++) {
    m+=(v[n]>0);
  }
  return m;
}

int niik_count_zero_from_int_vector(int *v,int num) {
  int n,m;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return 0;
  }
  for(n=m=0; n<num; n++) {
    m+=(v[n]==0);
  }
  return m;
}

int niik_get_max_index_from_int_vector(int *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]<v[n]) m=n;
  }
  return m;
}

int niik_get_min_index_from_int_vector(int *v,int num) {
  int m,n;
  if(v==NULL) {
    fprintf(stderr,"ERROR: v is null\n");
    return -1;
  }
  for(n=m=0; n<num; n++) {
    if(v[m]>v[n]) m=n;
  }
  return m;
}

int niik_get_max_from_int_vector(int *v,int num) {
  return v[niik_get_max_index_from_int_vector(v,num)];
}

int niik_get_min_from_int_vector(int *v,int num) {
  return v[niik_get_min_index_from_int_vector(v,num)];
}








/***************************************************
 *
 * B-SPLINE FUNCTIONS
 *
 * -first calculate A matrix  (niik_bspline_update_A)
 * -second calculate the b-spline coefficients (niik_bspline_calc_coeff)
 * -finally evaluate the b-splines for interpolations (niik_bspline_eval)
 *
 ***************************************************/

double niik_bspline_eval(double *c,int num,double x)
/* given the list of b-spline coefficients, interpolate at x and returns it  */
{
  double
    w[4],d;
  int
    degree=3,md=1,
    ii,ei,m,v;
  v=floor(x);
  ii=v-md;
  ei=ii+degree;
  d=x-v;
  /*fprintf(stdout,"  d  %8.3f\n",d);
    fprintf(stdout,"  v  %i\n",v);
    fprintf(stdout,"  ii %i\n",ii);
    fprintf(stdout,"  ei %i\n",ei);*/
  w[3]=d*d*d/6.0l;
  w[0]=1.0l/6 + (d*d-d) / 2.0l - w[3];
  w[2]=d+w[0]-2.0l*w[3];
  w[1]=1.0l-w[0]-w[2]-w[3];
  ei=NIIK_IMIN(ei,num-1);
  for(m=ii,d=0,v=0; m<=ei; v++,m++) {
    d+=w[v]*c[m];
  }
  return d;
}

int niik_bspline_calc_coeff_update(niikmat *A,double *v,int num)
/* -A must be pre-computed */
{
  double *ctmp;
  int n,ki,ko,k;
  ctmp=niik_calloc_double_vector(num);
  for(n=0; n<num; n++) {
    ctmp[n]=0;
    ki=NIIK_IMINMAX(n-20,0,num-1);
    ko=NIIK_IMINMAX(n+20,0,num-1);
    for(k=ki; k<=ko; k++) {
      ctmp[n]+=A->m[n][k]*v[k];
    }
  }
  for(n=0; n<num; n++) {
    v[n]=ctmp[n];
  }
  free(ctmp);
  return 1;
}

int niik_bspline_update_A(niikmat *A)
/* calculate matrix
 * A must be pre-defined  */
{
  int n;
  if(A->row<150) {
    for(n=1; n<A->row-1; n++) {
      A->m[n][n-1]=1.0l/6;
      A->m[n][n+1]=1.0l/6;
      A->m[n][n  ]=2.0l/3;
    }
    A->m[0][0]=A->m[A->row-1][A->row-1]=2.0l/3;
    A->m[0][1]=A->m[A->row-1][A->row-2]=1.0l/6;
    niikmat_inverse_update(A);
  } else {
    if(!niik_bspline_update_A_precomp_matrix(A)) {
      fprintf(stderr,"ERROR: niik_bspline_update_A_precomp_matrix\n");
      return 0;
    }
  }
  return 1;
}

int niik_bspline_update_A_precomp_matrix(niikmat *A)
/* pre-computed matrix for computing bspline coefficient
 * -must be sufficiently large
 * -otherwise, matrix inversion should be possible */
{
  int i,j,k,num;
  double m[512]= { 1.607695155,-0.430780618,0.115427319,-0.030928657,0.008287309,-0.002220578,0.000595002,-0.00015943,
                   0.000042719,-0.000011447,0.000003067,-0.000000822,0.00000022,-0.000000059,0.000000016,-0.000000004,
                   0.000000001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   -0.430780618,1.723122473,-0.461709275,0.123714627,-0.033149235,
                   0.008882311,-0.002380008,0.000637721,-0.000170877,0.000045786,
                   -0.000012268,0.000003287,-0.000000881,0.000000236,-0.000000063,
                   0.000000017,-0.000000005,0.000000001,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0.115427319,-0.461709275,
                   1.731409782,-0.463929853,0.124309629,-0.033308665,0.00892503,
                   -0.002391455,0.000640788,-0.000171699,0.000046007,-0.000012327,
                   0.000003303,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   -0.030928657,0.123714627,-0.463929853,1.732004784,-0.464089283,
                   0.124352349,-0.033320111,0.008928097,-0.002392276,0.000641009,
                   -0.000171758,0.000046022,-0.000012332,0.000003304,-0.000000885,
                   0.000000237,-0.000000064,0.000000017,-0.000000005,0.000000001,
                   0,0,0,0,0,0,0,0,0,0,0,0,0.008287309,-0.033149235,0.124309629,
                   -0.464089283,1.732047503,-0.46410073,0.124355416,-0.033320933,
                   0.008928317,-0.002392335,0.000641024,-0.000171762,0.000046023,
                   -0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,
                   0.000000017,-0.000000005,0.000000001,0,0,0,0,0,0,0,0,0,0,0,
                   -0.002220578,0.008882311,-0.033308665,0.124352349,-0.46410073,
                   1.73205057,-0.464101552,0.124355636,-0.033320992,0.008928333,
                   -0.00239234,0.000641025,-0.000171762,0.000046024,-0.000012332,
                   0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,0,0,0,0,0,0,0,0.000595002,
                   -0.002380008,0.00892503,-0.033320111,0.124355416,-0.464101552,
                   1.732050791,-0.464101611,0.124355652,-0.033320996,0.008928334,
                   -0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,
                   0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,0,0,0,0,0,0,-0.00015943,0.000637721,
                   -0.002391455,0.008928097,-0.033320933,0.124355636,-0.464101611,
                   1.732050806,-0.464101615,0.124355653,-0.033320997,0.008928334,
                   -0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,
                   0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,0,0,0,0,0,0.000042719,-0.000170877,
                   0.000640788,-0.002392276,0.008928317,-0.033320992,0.124355652,
                   -0.464101615,1.732050807,-0.464101615,0.124355653,-0.033320997,
                   0.008928334,-0.00239234,0.000641026,-0.000171762,0.000046024,
                   -0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,
                   0.000000017,-0.000000005,0.000000001,0,0,0,0,0,0,0,-0.000011447,
                   0.000045786,-0.000171699,0.000641009,-0.002392335,0.008928333,
                   -0.033320996,0.124355653,-0.464101615,1.732050808,-0.464101615,
                   0.124355653,-0.033320997,0.008928334,-0.00239234,0.000641026,
                   -0.000171762,0.000046024,-0.000012332,0.000003304,-0.000000885,
                   0.000000237,-0.000000064,0.000000017,-0.000000005,0.000000001,
                   0,0,0,0,0,0,0.000003067,-0.000012268,0.000046007,-0.000171758,
                   0.000641024,-0.00239234,0.008928334,-0.033320997,0.124355653,
                   -0.464101615,1.732050808,-0.464101615,0.124355653,-0.033320997,
                   0.008928334,-0.00239234,0.000641026,-0.000171762,0.000046024,
                   -0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,
                   0.000000017,-0.000000005,0.000000001,0,0,0,0,0,-0.000000822,
                   0.000003287,-0.000012327,0.000046022,-0.000171762,0.000641025,
                   -0.00239234,0.008928334,-0.033320997,0.124355653,-0.464101615,
                   1.732050808,-0.464101615,0.124355653,-0.033320997,0.008928334,
                   -0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,
                   0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,0,0.00000022,-0.000000881,0.000003303,
                   -0.000012332,0.000046023,-0.000171762,0.000641026,-0.00239234,
                   0.008928334,-0.033320997,0.124355653,-0.464101615,1.732050808,
                   -0.464101615,0.124355653,-0.033320997,0.008928334,-0.00239234,
                   0.000641026,-0.000171762,0.000046024,-0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,
                   -0.000000005,0.000000001,0,0,0,-0.000000059,0.000000236,-0.000000885,0.000003304,-0.000012332,0.000046024,
                   -0.000171762,0.000641026,-0.00239234,0.008928334,-0.033320997,0.124355653,-0.464101615,1.732050808,-0.464101615,
                   0.124355653,-0.033320997,0.008928334,-0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,0.000003304,
                   -0.000000885,0.000000237,-0.000000064,0.000000017,-0.000000005,0.000000001,0,0,0.000000016,-0.000000063,0.000000237,
                   -0.000000885,0.000003304,-0.000012332,0.000046024,-0.000171762,0.000641026,-0.00239234,0.008928334,-0.033320997,
                   0.124355653,-0.464101615,1.732050808,-0.464101615,0.124355653,-0.033320997,0.008928334,-0.00239234,0.000641026,
                   -0.000171762,0.000046024,-0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,-0.000000005,
                   0.000000001,0,-0.000000004,0.000000017,-0.000000064,0.000000237,-0.000000885,0.000003304,-0.000012332,0.000046024,-0.000171762,
                   0.000641026,-0.00239234,0.008928334,-0.033320997,0.124355653,-0.464101615,1.732050808,-0.464101615,0.124355653,-0.033320997,0.008928334,
                   -0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,-0.000000005,0.000000001
                 };

  double v[33]= { 0.000000001,-0.000000005,0.000000017,-0.000000064,0.000000237,-0.000000885,0.000003304,-0.000012332,0.000046024,-0.000171762,
                  0.000641026,-0.00239234,0.008928334,-0.033320997,0.124355653,-0.464101615,1.732050808,-0.464101615,0.124355653,-0.033320997,0.008928334,
                  -0.00239234,0.000641026,-0.000171762,0.000046024,-0.000012332,0.000003304,-0.000000885,0.000000237,-0.000000064,0.000000017,-0.000000005,0.000000001
                };

  num=A->row;
  if(num<128) {
    fprintf(stderr,"ERROR: matrix is too small\n");
    return 0;
  }
  /* beginning */
  for(i=k=0; i<16; i++) {
    for(j=0; j<32; j++) {
      A->m[i][j]=m[k++];
    }
  }
  /* fprintf(stdout,"final k = %i\n",k); */
  /* middle */
  for(i=16; i<num-16; i++) {
    for(j=i-16,k=0; j<=i+16; j++) {
      A->m[i][j]=v[k++];
    }
  }
  /* fprintf(stdout,"final k = %i\n",k); */
  /* ending */
  for(i=num-16,k=511; i<num; i++) {
    for(j=num-32; j<num; j++) {
      A->m[i][j]=m[k--];
    }
  }
  /*fprintf(stdout,"final k = %i\n",k);
    fprintf(stdout,"final i,j = %i %i\n",i,j);*/
  return 1;
}


/**********************************************************
 *
 * watershed functions
 **********************************************************/

int niik_get_watershed_from_tallest_peaks_double_vector(double *y,int num,int pnum) {
  int m,n,nn;
  int *v=NULL;
  v=(int *)calloc(num,sizeof(int));
  for(n=0; n<num; n++) v[n]=0;
  for(m=0; m<pnum; m++) {
    for(n=0; n<num; n++) {
      if(v[n]>0) v[n]+=1;
    }
    for(n=nn=0; n<num-1; n++) {
      if(y[n]<y[n+1]) continue;
      if(y[n]<y[n-1]) continue;
      if(v[n]>0) continue;
      if(!nn) nn=n;
      else if(y[n]>y[nn])
        nn=n;
    }
    if(nn) {
      if(!niik_get_watershed_double_vector_update(y,v,num,nn)) {
        fprintf(stderr,"[%s] niik_get_watershed_double_vector_update\n",__func__);
        return 0;
      }
    }
  } /* for pnum */
  /*for(n=0;n<num;n++){ fprintf(stdout,"%i %i %f\n",n,v[n],y[n]); }*/
  for(n=0; n<num; n++) {
    y[n] *= (v[n]>0);
  }
  free(v);
  v=NULL;
  return 1;
} /* niik_get_watershed_from_tallest_peaks_double_vector */

int niik_get_watershed_double_vector_untouch(double *y,int *f,int num,int p)
/* -update f so that 0 is outside the watershed region
 * -it does not update the watershed region
 * -y is the vector of interest; to do watershed processing
 * -num is the vector size for y and f
 * -p is the position of peak (0,num-1)
 */
{
  int n,nn;
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",__func__);
    return 0;
  }
  nn=num-1;
  for(n=p; n<nn; n++) {
    if(y[n]<=y[n+1])
      break;
  }
  for( ; n<num; n++) {
    f[n]=0;
  }
  for(n=p; n>0; n--) {
    if(y[n]<=y[n-1])
      break;
  }
  for( ; n>=0; n--) {
    f[n]=0;
  }
  return 1;
} /* niik_get_watershed_double_vector_untouch */

int niik_get_watershed_double_vector_update(double *y,int *f,int num,int p)
/* -update f so that 1 is in the watershed region
 * -it does not update outside the watershed region
 * -opposite of niik_get_watershed_double_vector_untouch
 * -y is the vector of interest; for watershed processing
 * -num is the vector size for y and f
 * -p is the position of peak (0,num-1)
 */
{
  int n;
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",__func__);
    return 0;
  }
  f[p]=1;
  for(n=p+1; n<num; n++) {
    if(y[n-1]<=y[n]) break;
    f[n]=1;
  }
  for(n=p-1; n>=0; n--) {
    if(y[n+1]<=y[n]) break;
    f[n]=1;
  }
  return 1;
} /* niik_get_watershed_double_vector_untouch */

int niik_get_watershed_double_vector(double *y,int num,int p)
/* get the N-large distributions in y */
{
  int n,nn;
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",__func__);
    return 0;
  }
  nn=num-1;
  for(n=p; n<nn; n++) {
    if(y[n]<=y[n+1])
      break;
  }
  for( ; n<num; n++) {
    y[n]=0;
  }
  for(n=p; n>0; n--) {
    if(y[n]<=y[n-1])
      break;
  }
  for( ; n>=0; n--) {
    y[n]=0;
  }
  return 1;
} /* niik_get_watershed_double_vector */


int niik_get_mode_double_vector(double *x,double *y,int num,double *xi,double *yi)
/* get the peak position (no sub-sampling) */
{
  int n;
  if(x==NULL) {
    fprintf(stderr,"[%s] ERROR: x is null\n",__func__);
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",__func__);
    return 0;
  }
  if((n=niik_get_max_index_double_vector(y,num))<0) {
    fprintf(stderr,"[%s] ERROR: niik_get_max_index_double_vector\n",__func__);
    return 0;
  }
  *xi=x[n];
  *yi=y[n];
  return 1;
} /* niik_get_mode_double_vector */

int niik_get_mode_bspline_vector(double *v,int num,double *xi,double *yi)
/* -A matrix is not pre-defined
 * -can be slow
 */
{
  if(!niik_get_mode_bspline_vector_with_A_matrix(v,num,xi,yi,NULL)) {
    fprintf(stderr,"[%s] ERROR: niik_get_mode_bspline_vector_with_A_matrix\n",__func__);
    return 0;
  }
  return 1;
}

int niik_get_mode_bspline_vector_with_A_matrix(double *v,int num,double *xi,double *yi,niikmat *A_mat)
/* -calculates the mode using b-spline for sub-bin size
 * -A matrix is pre-defined
 * -can be faster is A matrix size did not change */
{
  niikmat *A=NULL;
  niikvec *coeff=NULL;
  int i,n,m;
  double dx,xvec[6],yvec[6];
  int verbose=0;
  if(verbose>=1) niik_fc_display(__func__,1);
  if(v==NULL) {
    fprintf(stderr,"[%s] ERROR: v is null\n",__func__);
    return 0;
  }
  /* find the max position */
  n=niik_get_max_index_double_vector(v,num);
  /* if the max is at the end, no spline to estimate the max */
  if(n==0 || n==num-1) {
    *xi=(double)n;
    *yi=v[n];
    return 1;
  }
  /* set min/max */
  xvec[0]=n-1.0;
  xvec[5]=n+1.0;
  if(verbose>=1) fprintf(stdout,"[%s] %2i %6.2f %6.2f\n",__func__,n,xvec[0],xvec[5]);
  /* spline prep */
  if(A_mat==NULL) {
    A=niikmat_init(num,num);
    niik_bspline_update_A(A);
  } else {
    A=A_mat;
  }
  coeff=niikvec_init(num);
  for(i=0; i<num; i++) coeff->v[i]=v[i];
  niik_bspline_calc_coeff_update(A,coeff->v,num);
  /* find max in finer resolution */
  for(i=0; i<9; i++) {
    dx=xvec[5]-xvec[0];
    for(m=0; m<6; m++) {
      xvec[m] = xvec[0]+dx*m/5.0;
      yvec[m] = niik_bspline_eval(coeff->v,num,xvec[m]);
      if(verbose>=1) fprintf(stdout,"[%s] %2i %9.5f  %12.8f\n",__func__,m,xvec[m],yvec[m]);
    }
    n=niik_get_max_index_double_vector(yvec,6);
    if(verbose>=1) fprintf(stdout,"[%s] %8.5f %12.8f\n\n",__func__,xvec[n],yvec[n]);
    if(n<0) {
      fprintf(stderr,"[%s] ERROR: n is invalid %i\n",__func__,n);
      A=niikmat_free(A);
      coeff=niikvec_free(coeff);
      return 0;
    }
    if(n>5) {
      fprintf(stderr,"[%s] ERROR: n is invalid %i\n",__func__,n);
      A=niikmat_free(A);
      coeff=niikvec_free(coeff);
      return 0;
    }
    xvec[0]=xvec[n]-dx/5.0;
    xvec[5]=xvec[n]+dx/5.0;
  }
  /* update the value again */
  n=niik_get_max_index_double_vector(yvec,6);
  yvec[n] = niik_bspline_eval(coeff->v,num,xvec[n]);
  if(verbose>=1) fprintf(stdout,"[%s] %8.5f %12.8f\n",__func__,xvec[n],yvec[n]);
  /* update the outputs */
  *xi=xvec[n];
  *yi=yvec[n];
  /* free variables */
  if(A_mat==NULL) A=niikmat_free(A);
  niikvec_free(coeff);
  if(verbose>=1) niik_fc_display(__func__,0);
  return 1;
} /* int niik_get_mode_bspline_vector(double *v,int num,double *xi,double *yi); */


/****** end of bspline functions ******/






/******************************************************
 *
 * GENERAL MATRIX
 *
 ******************************************************/

niikmat *niikmat_init(int row,int col) {
  niikmat *mat;
  int i;
  if(row<=0) {
    fprintf(stderr,"[%s] ERROR: row is invalid\n",__func__);
    return NULL;
  }
  if(col<=0) {
    fprintf(stderr,"[%s] ERROR: col is invalid\n",__func__);
    return NULL;
  }
  mat=(niikmat *)calloc(1,sizeof(niikmat));
  mat->m=(double **)calloc(row,sizeof(double *));
  for(i=0; i<row; i++) {
    mat->m[i]=(double *)calloc(col,sizeof(double));
  }
  mat->row=row;
  mat->col=col;
  return mat;
}

niikmat *niikmat_rand(int row,int col) {
  niikmat *mat;
  int i,j;
  if(row<=0) {
    fprintf(stderr,"[niikmat_rand] ERROR: row is invalid\n");
    return NULL;
  }
  if(col<=0) {
    fprintf(stderr,"[niikmat_rand] ERROR: col is invalid\n");
    return NULL;
  }
  mat=(niikmat *)calloc(1,sizeof(niikmat));
  mat->m=(double **)calloc(row,sizeof(double *));
  for(i=0; i<row; i++) {
    mat->m[i]=(double *)calloc(col,sizeof(double));
    for(j=0; j<col; j++) {
      mat->m[i][j]=niik_get_rand();
    }
  }
  mat->row=row;
  mat->col=col;
  return mat;
}

niikmat *niikmat_free(niikmat *mat) {
  int i,verbose=0;
  if(mat==NULL) return NULL;
  if(verbose) fprintf(stdout,"-v (niikmat_free) start A = %i-by-%i\n",mat->row,mat->col);
  for(i=0; i<mat->row; i++) {
    if(verbose) fprintf(stdout,"-v (niikmat_free) check row %i\n",i);
    if(mat->m[i]!=NULL)  {
      free(mat->m[i]);
    }
    if(verbose) fprintf(stdout,"-v (niikmat_free) finish row %i\n",i);
  }
  if(verbose) fprintf(stdout,"-v (niikmat_free) free matrix\n");
  free(mat->m);
  if(verbose) fprintf(stdout,"-v (niikmat_free) free object\n");
  free(mat);
  if(verbose) fprintf(stdout,"-v (niikmat_free) return null\n");
  return NULL;
}




int niikmat_kmul(niikmat *mat,double k) {
  int i,j;
  char fcname[24]="niikmat_kmul";
  if(mat==NULL) {
    fprintf(stderr,"[%s] ERROR: mat is null\n",fcname);
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      mat->m[i][j]*=k;
    }
  }
  return 1;
}

int niikmat_kadd(niikmat *mat,double k) {
  int i,j;
  char fcname[24]="niikmat_kadd";
  if(mat==NULL) {
    fprintf(stderr,"[%s] ERROR: mat is null\n",fcname);
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      mat->m[i][j]+=k;
    }
  }
  return 1;
}

int niikmat_set_all(niikmat *mat,double k) {
  int i,j;
  char fcname[24]="niikmat_set_all";
  if(mat==NULL) {
    fprintf(stderr,"[%s] ERROR: mat is null\n",fcname);
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      mat->m[i][j]=k;
    }
  }
  return 1;
}

niikmat *niikmat_identity(int row,int col) {
  niikmat *mat;
  int i,j;
  if((mat=niikmat_init(row,col))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  for(i=0; i<row; i++) {
    for(j=0; j<col; j++) {
      mat->m[i][j]=(i==j?1.0:0.0);
    }
  }
  return mat;
}

int niikmat_identity_update(niikmat *mat) {
  int i,j;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      mat->m[i][j]=(i==j?1.0:0.0);
    }
  }
  return 1;
}

int niikmat_display_msg(const char *msg, niikmat *mat) {
  if(mat==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(msg!=NULL)
    fprintf(stdout,"%s",msg);
  if(!niikmat_display(mat)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_display(mat)\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  return 1;
}

int niikmat_display(niikmat *mat) {
  int i,j;
  if(mat==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(mat->row<=0) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: row is invalid\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(mat->col<=0) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: col is invalid\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    fprintf(stdout,"\t    ");
    for(j=0; j<mat->col; j++) {
      fprintf(stdout,"%15.9f ",mat->m[i][j]);
    }
    fprintf(stdout,"\n");
  }
  return 1;
}

int mat44_display(mat44 mat) {
  niikmat *M;
  M=niikmat_mat44_matrix(mat);
  niikmat_display(M);
  M=niikmat_free(M);
  return 1;
}

niikmat *niikmat_copy(niikmat *mat) {
  niikmat *out;
  int i,j;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return NULL;
  }
  if((out=niikmat_init(mat->row,mat->col))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      out->m[i][j]=mat->m[i][j];
    }
  }
  return out;
}

int niikmat_copy_update(niikmat *src,niikmat *dest) {
  int i,j;
  if(src==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  if(dest==NULL) {
    fprintf(stderr,"ERROR: dest is null\n");
    return 0;
  }
  if(src->row!=dest->row) {
    fprintf(stderr,"ERROR: rows are different %i %i\n",src->row,dest->row);
    return 0;
  }
  if(src->col!=dest->col) {
    fprintf(stderr,"ERROR: cols are different %i %i\n",src->col,dest->col);
    return 0;
  }
  for(i=0; i<src->row; i++) {
    for(j=0; j<src->col; j++) {
      dest->m[i][j]=src->m[i][j];
    }
  }
  return 1;
}

int niikmat_clear(niikmat *mat) {
  int i,j;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: src is null\n");
    return 0;
  }
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      mat->m[i][j]=0;
    }
  }
  return 1;
}

niikmat *niikmat_transpose_free(niikmat *mat) {
  niikmat *out=NULL;
  NIIK_RET0(((out=niikmat_transpose(mat))==NULL),__func__,"niikmat_transpose");
  mat=niikmat_free(mat);
  return out;
}

niikmat *niikmat_transpose(niikmat *mat) {
  niikmat *out;
  int i,j;
  if(mat==NULL) {
    fprintf(stderr,"[niikmat_transpose] ERROR: mat is null\n");
    return NULL;
  }
  out=niikmat_init(mat->col,mat->row);
  for(i=0; i<mat->row; i++) {
    for(j=0; j<mat->col; j++) {
      out->m[j][i]=mat->m[i][j];
    }
  }
  return out;
}

niikmat *niikmat_reshape_free(niikmat *mat,int row,int col) {
  niikmat *out=NULL;
  int i,j,n=0,m=0;
  NIIK_RETURN((row*col!=mat->row*mat->col),"wrong num of elements",NULL);
  out=niikmat_init(row,col);
  for(i=0; i<row; i++) {
    for(j=0; j<col; j++) {
      out->m[i][j] = mat->m[m][n++];
      if(n==mat->col) {
        m++;
        n=0;
      }
    }
  }
  mat=niikmat_free(mat);
  return out;
}

niikmat *niikmat_reshape(niikmat *mat,int row,int col) {
  niikmat *out=NULL;
  int i,j,n=0,m=0;
  NIIK_RETURN((row*col!=mat->row*mat->col),"wrong num of elements",NULL);
  mat=niikmat_init(row,col);
  for(i=0; i<row; i++) {
    for(j=0; j<col; j++) {
      out->m[i][j] = mat->m[m][n++];
      if(n==mat->col) {
        m++;
        n=0;
      }
    }
  }
  return out;
}

int niikmat_inverse_update(niikmat *mat) {
  niikmat *tmp=NULL;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  if((tmp=niikmat_inverse(mat))==NULL) {
    fprintf(stderr,"ERROR: niikmat_inverse_matrix\n");
    return 0;
  }
  if(!niikmat_copy_update(tmp,mat)) {
    fprintf(stderr,"ERROR: niikmat_copy_update\n");
    return 0;
  }
  tmp=niikmat_free(tmp);
  return 1;
}

niikmat *niikmat_inverse(niikmat *mat)
/* matrix inversion
 *
 * Revision
 * -2012-08-29, Kunio
 * -corrected memory leak problem*/
{
  niikmat *Lmat=NULL,*Umat=NULL,*out=NULL;
  int *pvec=NULL;
  double *v=NULL;
  int
  verbose=0,
  i,j,k,N;
  if(verbose) fprintf(stdout,"-d (niikmat_inverse) start\n");
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return NULL;
  }
  if(mat->row!=mat->col) {
    fprintf(stderr,"ERROR: not a square matrix %i %i\n",mat->row,mat->col);
    return NULL;
  }
  if(mat->row>2100) { /* arbitrary limit */
    fprintf(stderr,"ERROR: matix is too large for this function, %i\n",mat->row);
    return NULL;
  }
  N = mat->row;
  out=niikmat_copy(mat);
  Umat=niikmat_init(mat->row,mat->col);
  Lmat=niikmat_init(mat->row,mat->col);
  pvec=(int *)calloc(mat->row,sizeof(int *));
  /* permutation matrix (vector format) */
  for(i=0; i<N; i++) {
    pvec[i]=i;
  }
  /*niik_display_int_vector(pvec,N);*/
  for(i=0; i<N; i++) {
    for(j=k=i; j<N; j++) {
      if(fabs(out->m[j][i])>fabs(out->m[k][i])) {
        k=j;
      }
    }
    if(k>=N) continue;
    if(fabs(out->m[i][i])>fabs(out->m[k][i])) continue;
    for(j=0; j<N; j++) {
      NIIK_DSWAP(&out->m[k][j],&out->m[i][j]);
    }
    NIIK_ISWAP(&pvec[k],&pvec[i]);
    /*niik_display_int_vector(pvec,N);
      fprintf(stdout,"\n%i",i);
      niikmat_display(mat);*/
  }
  /* niik_display_int_vector(pvec,N);
     niikmat_display(mat);*/
  if(verbose) {
    fprintf(stdout,"-d (niikmat_inverse) lu_decompose\n\norig");
    niikmat_display(mat);
  }
  if(!niikmat_lu_decompose(out,Umat,Lmat)) {
    fprintf(stderr,"ERROR: niikmat_lu_decompose\n");
    return NULL;
  }
  if(verbose) {
    fprintf(stdout,"-d (niikmat_inverse) results of lu_decompose\n");
    fprintf(stdout,"up");
    niikmat_display(Umat);
    fprintf(stdout,"lo");
    niikmat_display(Lmat);
  }
  v=niik_calloc_double_vector(N);
  for(k=0; k<N; k++) { /* by column */
    if(verbose) fprintf(stdout,"  k = %i\n",k);
    for(i=0; i<N; i++) {
      v[i]=(i==k);
    }
    if(verbose) {
      fprintf(stdout,"start:    ");
      niik_display_double_vector(v,N);
    }
    /* forward substitution method */
    for(i=0; i<N; i++) {
      for(j=i-1; j>=0; j--) {
        v[i]-=Lmat->m[i][j]*v[j];
      }
    }
    if(verbose) {
      fprintf(stdout,"fwd fin:  ");
      niik_display_double_vector(v,N);
    }
    /* backward substitution method */
    v[N-1]=v[N-1]/Umat->m[N-1][N-1];
    for(i=N-2; i>=0; i--) {
      for(j=i+1; j<N; j++) {
        v[i]-=Umat->m[i][j]*v[j];
      }
      v[i]/=Umat->m[i][i];
    }
    if(verbose) {
      fprintf(stdout,"bac fin: ");
      niik_display_double_vector(v,N);
    }
    /* put values into output matrix */
    for(i=0; i<N; i++) {
      out->m[i][pvec[k]]=v[i];
    }
    if(verbose) {
      fprintf(stdout,"out");
      niikmat_display(out);
    }
  }
  free(pvec);
  Lmat=niikmat_free(Lmat);
  Umat=niikmat_free(Umat);
  free(v);
  return out;
}

int niikmat_pinv(niikmat *mat)
/*
 * pseudo-inverse
 * using singular value decomposition
 * updates mat
 */
{
  char fcname[64]="niikmat_pinv";
  niikmat *U,*V,*W;
  double *w;
  long i,j;
  int const verbose=0;
  U=niikmat_copy(mat);
  V=niikmat_init(mat->col,mat->col);
  w=niik_calloc_double_vector(mat->col);
  if(verbose>0) fprintf(stdout,"[%s:%i:%s] niik_singular_value_decomposition\n",__FILE__,__LINE__,fcname);
  if(!niik_singular_value_decomposition(U,w,V)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_singular_value_decomposition\n",__FILE__,__LINE__,fcname);
    return 0;
  }
  if(verbose>=1) {
    NIIK_EXIT((!niikmat_display_msg("U = \n",U)),fcname,"display U",1);
    fprintf(stdout,"[%s] w = \n",fcname);
    niik_display_double_vector(w,mat->col);
    NIIK_EXIT((!niikmat_display_msg("V = \n",V)),fcname,"display V",1);
  }
  W=niikmat_init(mat->col,mat->col); /* W+ */
  for(i=0; i<mat->col; i++)
    for(j=0; j<mat->col; j++)
      if(i==j) W->m[i][j]=1.0/w[j];
      else W->m[i][j]=0;
  free(w);
  NIIK_RET0(((U=niikmat_transpose_free(U))==NULL),fcname,"niikmat_transpose_free");
  NIIK_EXIT(((U=niikmat_multiply_free12(W,U))==NULL),fcname,"niikmat_multiply_free12 W*U",1);
  if(verbose>=1) NIIK_EXIT((!niikmat_display_msg("WU = \n",U)),fcname,"display WU",1);
  NIIK_EXIT(((U = niikmat_multiply_free12(V,U))==NULL),fcname,"niikmat_multiply_free12 V*W*U",1);
  if(verbose>=1) NIIK_EXIT((!niikmat_display_msg("pinv = \n",U)),fcname,"display pinv",1);
  for(i=0; i<mat->row; i++)
    if(mat->m[i]!=NULL)
      free(mat->m[i]);
  free(mat->m);
  mat->m=U->m;
  mat->row=U->row;
  mat->col=U->col;
  free(U);
  return 1;
} /* niikmat_pinv */

niikmat *niikmat_pseudo_inverse(niikmat *mat) {
  niikmat *tmp=NULL,*out=NULL;
  if((tmp=niikmat_copy(mat))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_copy\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((out=niikmat_pseudo_inverse_free(tmp))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_pseudo_inverse\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  return out;
} /* niikmat_pseudo_inverse */

niikmat *niikmat_pseudo_inverse_free(niikmat *mat) {
  niikmat *tmp=NULL,*out=NULL;
  if((tmp=niikmat_copy(mat))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_copy\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((out=niikmat_pseudo_inverse_free(tmp))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_pseudo_inverse\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  mat=niikmat_free(mat);
  return out;
} /* niikmat_pseudo_inverse */

double dpythag(double a, double b)
/* returns sqrt(a^2+b^2) */
{
  double aa,bb;
  aa=fabs(a);
  bb=fabs(b);
  if(aa>bb) return aa*sqrt(1.0+NIIK_SQ(bb/aa));
  else if(bb==0.0) return 0.0;
  return bb*sqrt(1.0+NIIK_SQ(aa/bb));
}

int niik_singular_value_decomposition(niikmat *a,double *w, niikmat *v)
/* -singular value decomposition where a = U w trans(v)
 * -a is the input matrix MxN, compact form
 * -w is the N-length vector, allocate memory before this function, updated on output
 * -v is the NxN matrix, allocate memory before this function, updated on output
 * -a is updated on output as U
 */
{
  char fcname[64]="niik_singular_value_decomposition";
  double **A,**V;
  int i,j;
  int const verbose=0;
  if(verbose>=1) {
    fprintf(stdout,"row = %i\n",a->row);
    fprintf(stdout,"col = %i\n",a->col);
  }
  A=(double **)calloc(a->row+1,sizeof(double *));
  for(i=1; i<=a->row; i++) {
    A[i]=(double *)calloc(a->col+1,sizeof(double));
    for(j=1; j<=a->col; j++) {
      A[i][j] = a->m[i-1][j-1];
    }
  }
  V=(double **)calloc(a->col+1,sizeof(double *));
  for(i=1; i<=a->col; i++) {
    V[i]=(double *)calloc(a->col+1,sizeof(double));
  }
  if(verbose>0) fprintf(stdout,"[%s] svdcmp\n",fcname);
  svdcmp(A,a->row,a->col,w-1,V);
  for(i=1; i<=a->col; i++) {
    for(j=1; j<=a->col; j++) {
      v->m[i-1][j-1]=V[i][j];
    }
  }
  for(i=1; i<=a->col; i++) {
    free(V[i]);
  }
  free(V);
  for(i=1; i<=a->row; i++) {
    for(j=1; j<=a->col; j++) {
      a->m[i-1][j-1]=A[i][j];
    }
    free(A[i]);
  }
  free(A);
  return 1;
}

int svdcmp(double **a, int m, int n, double *w, double **v)
/* from numerical recipe book */
{
  int flag,i,its,j,jj,k,l,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  int verbose=0;
  if(verbose>=1) {
    fprintf(stdout,"  a (input) = \n");
    for(i=1; i<=m; i++) {
      for(j=1; j<=n; j++) {
        fprintf(stdout,"%9.5f ",a[i][j]);
      }
      fprintf(stdout,"\n");
    }
  }
  rv1=(double *)calloc(1+n,sizeof(double));
  g=scale=anorm=0.0;
  for (i=1; i<=n; i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i; k<=m; k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i; k<=m; k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        /* g = -SIGN(sqrt(s),f);*/
        if(f<0) g=sqrt(s);
        else    g=-sqrt(s);
        h=f*g-s;
        a[i][i]=f-g;
        for (j=l; j<=n; j++) {
          for (s=0.0,k=i; k<=m; k++) s += a[k][i]*a[k][j];
          f=s/h;
          for (k=i; k<=m; k++) a[k][j] += f*a[k][i];
        }
        for (k=i; k<=m; k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l; k<=n; k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l; k<=n; k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        /*g = -SIGN(sqrt(s),f);*/
        if(f<0) g= sqrt(s);
        else    g=-sqrt(s);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l; k<=n; k++) rv1[k]=a[i][k]/h;
        for (j=l; j<=m; j++) {
          for (s=0.0,k=l; k<=n; k++) s += a[j][k]*a[i][k];
          for (k=l; k<=n; k++) a[j][k] += s*rv1[k];
        }
        for (k=l; k<=n; k++) a[i][k] *= scale;
      }
    }
    anorm=NIIK_DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n; i>=1; i--) {
    if (i < n) {
      if (g) {
        for (j=l; j<=n; j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l; j<=n; j++) {
          for (s=0.0,k=l; k<=n; k++) s += a[i][k]*v[k][j];
          for (k=l; k<=n; k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l; j<=n; j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=NIIK_IMIN(m,n); i>=1; i--) {
    l=i+1;
    g=w[i];
    for (j=l; j<=n; j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l; j<=n; j++) {
        for (s=0.0,k=l; k<=m; k++) s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
        for (k=i; k<=m; k++) a[k][j] += f*a[k][i];
      }
      for (j=i; j<=m; j++) a[j][i] *= g;
    } else for (j=i; j<=m; j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n; k>=1; k--) {
    for (its=1; its<=30; its++) {
      flag=1;
      for (l=k; l>=1; l--) {
        nm=l-1;
        if ((double)(fabs(rv1[l])+anorm) == anorm) {
          flag=0;
          break;
        }
        if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l; i<=k; i++) {
          f=s*rv1[i];
          rv1[i]=c*rv1[i];
          if ((double)(fabs(f)+anorm) == anorm) break;
          g=w[i];
          h=dpythag(f,g);
          w[i]=h;
          h=1.0/h;
          c=g*h;
          s = -f*h;
          for (j=1; j<=m; j++) {
            y=a[j][nm];
            z=a[j][i];
            a[j][nm]=y*c+z*s;
            a[j][i]=z*c-y*s;
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=1; j<=n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 30) {
        fprintf(stderr,"ERROR: no convergence in 30 svdcmp iterations");
        return 0;
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=dpythag(f,1.0);
      /*f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;*/
      if(f<0)
        f=((x-z)*(x+z)+h*((y/(f-g))-h))/x;
      else
        f=((x-z)*(x+z)+h*((y/(f+g))-h))/x;
      c=s=1.0;
      for (j=l; j<=nm; j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=dpythag(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g = g*c-x*s;
        h=y*s;
        y *= c;
        for (jj=1; jj<=n; jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=dpythag(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=c*g+s*y;
        x=c*y-s*g;
        for (jj=1; jj<=m; jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free(rv1);
  if(0) {
    fprintf(stdout,"  w = \n");
    for(i=1; i<=n; i++) {
      fprintf(stdout,"%9.5f ",w[i]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"  u = \n");
    for(i=1; i<=m; i++) {
      for(j=1; j<=n; j++) {
        fprintf(stdout,"%9.5f ",a[i][j]);
      }
      fprintf(stdout,"\n");
    }
    fprintf(stdout,"  v = \n");
    for(i=1; i<=n; i++) {
      for(j=1; j<=n; j++) {
        fprintf(stdout,"%9.5f ",v[i][j]);
      }
      fprintf(stdout,"\n");
    }
    exit(0);
  }
  return 1;
}





int niikmat_lu_decompose(niikmat *mat,niikmat *Umat,niikmat *Lmat)
/* LU decomposition
 * -mat is the input
 * -outputs are replaced in Umat and Lmat
 * -Umat and Lmat must have the same dimension as mat
 *  and their memory should be allocated before this function
 *
 * -currently there's a problem when the diagonal terms are zero
 */
{
  int i,j,k,N;
  double d;
  int verbose=0;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  if(mat->row!=mat->col) {
    fprintf(stderr,"ERROR: not a square matrix %i %i\n",mat->row,mat->col);
    return 0;
  }
  N = mat->row;
  if(!niikmat_copy_update(mat,Umat)) {
    fprintf(stderr,"ERROR: niikmat_copy_update\n");
    return 0;
  }
  if(!niikmat_identity_update(Lmat)) {
    fprintf(stderr,"ERROR: niikmat_identity_update\n");
    return 0;
  }
  if(verbose) {
    fprintf(stdout,"-d (niikmat_lu_decompose) start\n");
    niikmat_display(mat);
  }
  for(i=0; i<N-1; i++) { /* col */
    for(j=i+1; j<N; j++) {
      Lmat->m[j][i] = d = Umat->m[j][i] / Umat->m[i][i];
      if(Umat->m[i][i]==0.0l) Lmat->m[j][i] = d = 0;
      /*fprintf(stdout,"%i %i  %9.3f\n",i,j,d);*/
      for(k=i; k<N; k++) {
        Umat->m[j][k] = Umat->m[j][k] - Umat->m[i][k] * d;
      }
      if(verbose) {
        fprintf(stdout,"\n\n%i %i up",i,j);
        niikmat_display(Umat);
        fprintf(stdout,"%i %i lo",i,j);
        niikmat_display(Lmat);
      }
    }
  }
  /*Lmat->m[N-1][N-1]=Umat->m[N-1][N-1]/Lmat->m[N-1][N-1];*/
  /* Umat->m[N-1][N-1]=1; */
  if(verbose) {
    fprintf(stdout,"-d (niikmat_lu_decompose) results\nup");
    niikmat_display(Umat);
    fprintf(stdout,"lo");
    niikmat_display(Lmat);
    fprintf(stdout,"-d (niikmat_lu_decompose) end\n");
  }
  return 1;
} /* niikmat_lu_decomposition */





/************************************************************
 *
 * AFFINE MATRIX
 *
 ************************************************************/

int niikmat_scale_matrix_update(niikmat *mat,double sx,double sy,double sz) {
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  if(!niikmat_clear(mat)) {
    fprintf(stderr,"ERROR: niikmat_clear\n");
    return 0;
  }
  mat->m[0][0] = sx;
  mat->m[1][1] = sy;
  mat->m[2][2] = sz;
  mat->m[3][3] = 1.;
  return 1;
}

niikmat *niikmat_scale_matrix(double sx,double sy,double sz) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  mat->m[0][0] = sx;
  mat->m[1][1] = sy;
  mat->m[2][2] = sz;
  mat->m[3][3] = 1.;
  return mat;
}

int niikmat_shear_matrix_update(niikmat *mat,double sx,double sy,double sz) {
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  if(!niikmat_identity_update(mat)) {
    fprintf(stderr,"ERROR: niikmat_identity_update\n");
    return 0;
  }
  mat->m[0][1] = sx;
  mat->m[0][2] = sy;
  mat->m[1][2] = sz;
  return 1;
}

niikmat *niikmat_shear_matrix(double sx,double sy,double sz) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_shear_matrix_update(mat,sx,sy,sz)) {
    fprintf(stderr,"ERROR: niikmat_shear_matrix_update\n");
    return NULL;
  }
  return mat;
}

int niikmat_rotate_matrix_update(niikmat *mat,double rx,double ry,double rz) {
  double sx,sy,sz,cx,cy,cz;
  rx=NIIK_DEGREE2RAD(rx);
  ry=NIIK_DEGREE2RAD(ry);
  rz=NIIK_DEGREE2RAD(rz);
  sx = sin(rx);
  sy = sin(ry);
  sz = sin(rz);
  cx = cos(rx);
  cy = cos(ry);
  cz = cos(rz);
  mat->m[0][0]=cy*cz;
  mat->m[0][1]=-cy*sz;
  mat->m[0][2]=sy;
  mat->m[1][0]=cx*sz+sx*sy*cz;
  mat->m[1][1]=cx*cz-sx*sy*sz;
  mat->m[1][2]=-sx*cy;
  mat->m[2][0]=sx*sz-cx*sy*cz;
  mat->m[2][1]=sx*cz+cx*sy*sz;
  mat->m[2][2]=cx*cy;
  mat->m[3][0]=mat->m[3][1]=mat->m[3][2]=0;
  mat->m[3][3]=1;
  return 1;
}

niikmat *niikmat_rotate_matrix(double rx,double ry,double rz) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_rotate_matrix_update(mat,rx,ry,rz)) {
    fprintf(stderr,"ERROR: niikmat_rotate_matrix_update\n");
    return NULL;
  }
  return mat;
}

int niikmat_translate_matrix_update(niikmat *mat,double x,double y,double z) {
  if(!niikmat_identity_update(mat)) {
    fprintf(stderr,"ERROR: niikmat_identity_update\n");
    return 0;
  }
  mat->m[0][3]=x;
  mat->m[1][3]=y;
  mat->m[2][3]=z;
  return 1;
}

niikmat *niikmat_translate_matrix(double x,double y,double z) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_translate_matrix_update(mat,x,y,z)) {
    fprintf(stderr,"ERROR: niikmat_translate_matrix_update\n");
    return NULL;
  }
  return mat;
}

int niikmat_flip_matrix_update(niikmat *mat,char dir) {
  if(!niikmat_identity_update(mat)) {
    fprintf(stderr,"ERROR: niikmat_identity_update\n");
    return 0;
  }
  switch(dir) {
  case 'x':
    mat->m[0][0]=-1;
    break;
  case 'y':
    mat->m[1][1]=-1;
    break;
  case 'z':
    mat->m[2][2]=-1;
    break;
  }
  return 1;
}

niikmat *niikmat_flip_matrix(char dir) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_flip_matrix_update(mat,dir)) {
    fprintf(stderr,"ERROR: niikmat_flip_matrix_update\n");
    return NULL;
  }
  return mat;
}

int niikmat_mat44_matrix_update(niikmat *mat,mat44 m) {
  int i,j;
  if(!niikmat_identity_update(mat)) {
    fprintf(stderr,"ERROR: niikmat_identity_update\n");
    return 0;
  }
  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      mat->m[i][j]=m.m[i][j];
  return 1;
}

niikmat *niikmat_mat44_matrix(mat44 m) {
  niikmat *mat;
  if((mat=niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_mat44_matrix_update(mat,m)) {
    fprintf(stderr,"ERROR: niikmat_mat44_matrix_update\n");
    return NULL;
  }
  return mat;
}

mat44 niikmat_make_mat44(niikmat *mat) {
  mat44 m;
  int i,j;
  m.m[3][3]=0;
  if(mat==NULL)   {
    fprintf(stderr,"[niikmat_make_mat44] ERROR: mat is null\n");
    return m;
  }
  if(mat->row!=4) {
    fprintf(stderr,"[niikmat_make_mat44] ERROR: #row is not 4\n");
    return m;
  }
  if(mat->col!=4) {
    fprintf(stderr,"[niikmat_make_mat44] ERROR: #col is not 4\n");
    return m;
  }
  for(i=0; i<4; i++)
    for(j=0; j<4; j++)
      m.m[i][j]=mat->m[i][j];
  return m;
}

niikmat *niikmat_affine_matrix_val(double m00,double m01,double m02,double m03,
                                   double m10,double m11,double m12,double m13,
                                   double m20,double m21,double m22,double m23,
                                   double m30,double m31,double m32,double m33) {
  niikmat *afmat=NULL;
  if((afmat = niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  afmat->m[0][0]=m00;
  afmat->m[0][1]=m01;
  afmat->m[0][2]=m02;
  afmat->m[0][3]=m03;
  afmat->m[1][0]=m10;
  afmat->m[1][1]=m11;
  afmat->m[1][2]=m12;
  afmat->m[1][3]=m13;
  afmat->m[2][0]=m20;
  afmat->m[2][1]=m21;
  afmat->m[2][2]=m22;
  afmat->m[2][3]=m23;
  afmat->m[3][0]=m30;
  afmat->m[3][1]=m31;
  afmat->m[3][2]=m32;
  afmat->m[3][3]=m33;
  return afmat;
}

niikmat *niikmat_affine_matrix_new(double rx,double ry,double rz,
                                   double tx,double ty,double tz,
                                   double sx,double sy,double sz,
                                   double hx,double hy,double hz,
                                   double px,double py,double pz) {
  niikmat *afmat=NULL;
  if((afmat = niikmat_init(4,4))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  if(!niikmat_affine_matrix(afmat,rx,ry,rz,tx,ty,tz,sx,sy,sz,hx,hy,hz,px,py,pz)) {
    fprintf(stderr,"ERROR: niikmat_affine_matrix\n");
    return NULL;
  }
  return afmat;
}

int niikmat_affine_matrix(niikmat *mat,
                          double rx,double ry,double rz,
                          double tx,double ty,double tz,
                          double sx,double sy,double sz,
                          double hx,double hy,double hz,
                          double px,double py,double pz) {
  niikmat *tmp;
  int verbose=0;
  if(mat==NULL) {
    fprintf(stderr,"ERROR: mat is null\n");
    return 0;
  }
  if(mat->row!=4) {
    fprintf(stderr,"ERROR: mat is not 4-by-4 (%i,%i)\n",mat->row,mat->col);
    return 0;
  }
  if(mat->col!=4) {
    fprintf(stderr,"ERROR: mat is not 4-by-4 (%i,%i)\n",mat->row,mat->col);
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) translate matrix update\n");
  if(!niikmat_translate_matrix_update(mat,-px,-py,-pz)) {
    fprintf(stderr,"ERROR: niikmat_translate_matrix_update\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) scale matrix\n");
  tmp=niikmat_scale_matrix       (     sx, sy, sz);
  if(tmp==NULL) {
    fprintf(stderr,"ERROR: niikmat_scale_matrix\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) multiply matrix\n");
  niikmat_multiply_mat2(tmp,mat);
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) shear matrix update\n");
  niikmat_shear_matrix_update    (tmp, hx, hy, hz);
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) multiply matrix\n");
  niikmat_multiply_mat2(tmp,mat);
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) rotate matrix\n");
  niikmat_rotate_matrix_update   (tmp, rx, ry, rz);
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) multiply matrix\n");
  niikmat_multiply_mat2_free1(tmp,mat);
  if(verbose) fprintf(stderr,"-d (niikmat_affine_matrix) translate + pivot (adjust)\n");
  mat->m[0][3] += tx + px;
  mat->m[1][3] += ty + py;
  mat->m[2][3] += tz + pz;
  return 1;
}

double niikmat_mag(niikmat *A) {
  int i,j;
  double d=1e-7;
  if(A==NULL) {
    fprintf(stderr,"ERROR: input matrix is null\n");
    return -1;
  }
  for(i=0; i<A->row; i++) {
    for(j=0; j<A->col; j++) {
      d+=A->m[i][j]*A->m[i][j];
    }
  }
  if(d<1e-8) return 0;
  return sqrt(d);
}

int niikmat_add_to_a(niikmat *a,niikmat *b) {
  int i,j;
  NIIK_RETURN((a==NULL),"a is null",0);
  NIIK_RETURN((b==NULL),"a is null",0);
  NIIK_RETURN((a->row!=b->row),"#row is different",0);
  NIIK_RETURN((a->col!=b->col),"#col is different",0);
  for(i=0; i<a->row; i++) {
    for(j=0; j<a->col; j++) {
      a->m[i][j]+=b->m[i][j];
    }
  }
  return 1;
}


/****************************************
 * matrix multiplicatoin
 *
 * -general format
 * -niikmat4 may be no longer needed
 *
 ****************************************/

niikmat *niikmat_multiply(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2 */
{
  niikmat *out;
  int i,j,k;
  char fcname[32]="niikmat_multiply";
  if(mat1==NULL) {
    fprintf(stderr,"[%s] ERROR: mat1 is null\n",fcname);
    return NULL;
  }
  if(mat2==NULL) {
    fprintf(stderr,"[%s] ERROR: mat2 is null\n",fcname);
    return NULL;
  }
  if(mat1->col!=mat2->row) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat1.col(%i) and mat2.row(%i) are different\n",__FILE__,__LINE__,__func__,mat1->col,mat2->row);
    return NULL;
  }
  if((out=niikmat_init(mat1->row,mat2->col))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return NULL;
  }
  for(i=0; i<mat1->row; i++) {
    for(j=0; j<mat2->col; j++) {
      out->m[i][j]=0;
      for(k=0; k<mat1->col; k++) {
        out->m[i][j]+=mat1->m[i][k]*mat2->m[k][j];
      }
    }
  }
  return out;
}

niikmat *niikmat_multiply_free12(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2
 * free mat1 and mat2 */
{
  niikmat *out=NULL;
  char fcname[32]="niikmat_multiply_free12";
  if((out=niikmat_multiply(mat1,mat2))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_multiply\n",__FILE__,__LINE__,fcname);
    return NULL;
  }
  mat1=niikmat_free(mat1);
  mat2=niikmat_free(mat2);
  return out;
}

niikmat *niikmat_multiply_free1(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2
 * free mat1 */
{
  niikmat *out=NULL;
  char fcname[32]="niikmat_multiply_free1";
  if((out=niikmat_multiply(mat1,mat2))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_multiply\n",__FILE__,__LINE__,fcname);
    return NULL;
  }
  niikmat_free(mat1);
  return out;
}

niikmat *niikmat_multiply_free2(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2
 * free mat2 */
{
  niikmat *out=NULL;
  char fcname[32]="niikmat_multiply_free2";
  if((out=niikmat_multiply(mat1,mat2))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_multiply\n",__FILE__,__LINE__,fcname);
    return NULL;
  }
  mat2=niikmat_free(mat2);
  return out;
}

int niikmat_multiply_mat1(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2
 * mat1 is updated
 * mat2 is not freed; to free, use niikmat_multiply_mat1_free2  */
{
  niikmat *out;
  int i,j,k;
  const char *fcname="niikmat_multiply_mat1";
  if(mat1==NULL) {
    fprintf(stderr,"[%s] ERROR: mat1 is null\n",fcname);
    return 0;
  }
  if(mat2==NULL) {
    fprintf(stderr,"[%s] ERROR: mat2 is null\n",fcname);
    return 0;
  }
  if(mat1->col!=mat2->row) {
    fprintf(stderr,"ERROR: mat1.col(%i) and mat2.row(%i) are different\n",mat1->col,mat2->row);
    return 0;
  }
  if((out=niikmat_init(mat1->row,mat2->col))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_init\n",__FILE__,__LINE__,fcname);
    return 0;
  }
  for(i=0; i<mat1->row; i++) {
    for(j=0; j<mat2->col; j++) {
      out->m[i][j]=0;
      for(k=0; k<mat1->col; k++) {
        out->m[i][j]+=mat1->m[i][k]*mat2->m[k][j];
      }
    }
  }
  if(out->row==mat1->row && out->col==mat1->col) {
    for(i=0; i<out->row; i++) {
      for(j=0; j<out->col; j++) {
        mat1->m[i][j]=out->m[i][j];
      }
    }
    niikmat_free(out);
    return 1;
  } else {
    for(i=0; i<mat1->row; i++) if(mat1->m[i]!=NULL) free(mat1->m[i]);
    free(mat1->m);
    mat1->m=out->m;
    mat1->row=out->row;
    mat1->col=out->col;
    free(out);
  }
  return 1;
}

int niikmat_multiply_mat2(niikmat *mat1,niikmat *mat2)
/* output = mat2 * mat1
 * mat2 is updated
 * mat1 is not freed; to free, use niikmat_multiply_mat2_free1  */
{
  niikmat *out;
  int i,j,k,verbose=0;
  if(mat1==NULL) {
    fprintf(stderr,"ERROR: mat1 is null\n");
    return 0;
  }
  if(mat2==NULL) {
    fprintf(stderr,"ERROR: mat2 is null\n");
    return 0;
  }
  if(mat1->col!=mat2->row) {
    fprintf(stderr,"ERROR: mat1.col(%i) and mat2.row(%i) are different\n",mat1->col,mat2->row);
    return 0;
  }
  if((out=niikmat_init(mat1->row,mat2->col))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return 0;
  }
  for(i=0; i<mat1->row; i++) {
    for(j=0; j<mat2->col; j++) {
      out->m[i][j]=0;
      for(k=0; k<mat1->col; k++) {
        out->m[i][j]+=mat1->m[i][k]*mat2->m[k][j];
      }
    }
  }
  if(out->row==mat2->row && out->col==mat2->col) {
    if(verbose>0)fprintf(stdout,"\trow/col match\n");
    for(i=0; i<out->row; i++) {
      for(j=0; j<out->col; j++) {
        mat2->m[i][j]=out->m[i][j];
      }
    }
    niikmat_free(out);
    return 1;
  } else {
    if(verbose>0) {
      fprintf(stdout,"\trow/col different\n");
      niikmat_display(out);
    }
    for(i=0; i<mat2->row; i++) free(mat2->m[i]);
    free(mat2->m);
    mat2->m=out->m;
    mat2->row=out->row;
    mat2->col=out->col;
    free(out);
  }
  return 1;
}

int niikmat_multiply_mat1_free2(niikmat *mat1,niikmat *mat2)
/* output = mat1 * mat2
 * mat1 is updated
 * mat2 is freed  */
{
  niikmat *out;
  int i,j,k;
  if(mat1==NULL) {
    fprintf(stderr,"ERROR: mat1 is null\n");
    return 0;
  }
  if(mat2==NULL) {
    fprintf(stderr,"ERROR: mat2 is null\n");
    return 0;
  }
  if(mat1->col!=mat2->row) {
    fprintf(stderr,"ERROR: mat1.col(%i) and mat2.row(%i) are different\n",mat1->col,mat2->row);
    return 0;
  }
  if((out=niikmat_init(mat1->row,mat2->col))==NULL) {
    fprintf(stderr,"ERROR: niikmat_init\n");
    return 0;
  }
  for(i=0; i<mat1->row; i++) {
    for(j=0; j<mat2->col; j++) {
      out->m[i][j]=0;
      for(k=0; k<mat1->col; k++) {
        out->m[i][j]+=mat1->m[i][k]*mat2->m[k][j];
      }
    }
  }
  if(out->row==mat1->row && out->col==mat1->col) {
    for(i=0; i<out->row; i++) {
      for(j=0; j<out->col; j++) {
        mat1->m[i][j]=out->m[i][j];
      }
    }
    niikmat_free(mat2);
    niikmat_free(out);
    return 1;
  } else {
    if(mat1->m!=NULL) {
      for(i=0; i<mat1->row; i++) if(mat1->m[i]!=NULL) free(mat1->m[i]);
      free(mat1->m);
    }
    mat1->m=out->m;
    mat1->row=out->row;
    mat1->col=out->col;
    niikmat_free(mat2);
    free(out);
  }
  return 1;
}

int niikmat_multiply_mat2_free1(niikmat *mat1,niikmat *mat2)
/* output = mat2 * mat1
 * mat2 is updated
 * mat1 is freed  */
{
  niikmat *out;
  int i,j,k;
  if(mat1==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat1 is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(mat2==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat2 is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(mat1->col!=mat2->row) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat1.col(%i) and mat2.row(%i) are different\n",__FILE__,__LINE__,__func__,mat1->col,mat2->row);
    return 0;
  }
  if((out=niikmat_init(mat1->row,mat2->col))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_init\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(i=0; i<mat1->row; i++) {
    for(j=0; j<mat2->col; j++) {
      out->m[i][j]=0.0;
      for(k=0; k<mat1->col; k++) {
        out->m[i][j]+=mat1->m[i][k]*mat2->m[k][j];
      }
    }
  }
  if(out->row==mat2->row && out->col==mat2->col) {
    for(i=0; i<out->row; i++) {
      for(j=0; j<out->col; j++) {
        mat2->m[i][j]=out->m[i][j];
      }
    }
    niikmat_free(mat1);
    niikmat_free(out);
    return 1;
  } else {
    if(mat2!=NULL) {
      for(i=0; i<mat2->row; i++) if(mat2->m[i]!=NULL) free(mat2->m[i]);
      free(mat2->m);
    }
    mat2->m=out->m;
    mat2->row=out->row;
    mat2->col=out->col;
    niikmat_free(mat1);
    free(out);
  }
  return 1;
}

/***
 *** end of sub-section of niikmat's affine things
 *** (within section for niikmat function)
 ***
 ***/


/**** end of niikmat functions *****/




/*******************************************************
 *
 * DECOMPOSE MATRIX
 *
 *******************************************************/

niikmat *g_niikmat_decompose_affine_obj_func_mat;
int g_niikmat_decompose_affine_obj_func_dof;
double *g_niikmat_decompose_affine_obj_func_ipar;

int niikmat_decompose_affine(niikmat *mat,double *par,int dof) {
  int
  verbose=0,
  m,n,
  ndim,iter=165536;
  niikmat *p=NULL,*matout=NULL;
  double
  dout=0,
  tol=1e-5,(* pfn)()=NULL;
  if(mat==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: mat is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(par==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: par is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) start\n");
  g_niikmat_decompose_affine_obj_func_mat  = mat;
  g_niikmat_decompose_affine_obj_func_dof  = dof;
  if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) amoeba matrix\n");
  ndim=225;
  p=niikmat_init(ndim+1,ndim);
  for(n=7; n<=10; n++) {
    p->m[0][n] = 1;
  }
  for(n=14; n<=16; n++) {
    p->m[0][n] = par[n];
  }
  for(n=1; n<=16; n++) {
    p->m[1][n] = 1;
  }
  for(m=2; m<=ndim; m++) {
    p->m[m][0]=0;
    for(n=1; n<=3; n++) {
      p->m[m][n] = niik_get_rand() * 60.0 - 30.0;
    }
    for(n=4; n<=6; n++) {
      p->m[m][n] = niik_get_rand() * 40 - 20;
    }
    for(n=7; n<=10; n++) {
      p->m[m][n] = niik_get_rand() * 0.4 + 0.8;
    }
    for(n=11; n<=13; n++) {
      p->m[m][n] = niik_get_rand() * 0.4 - 0.2;
    }
    for(n=14; n<=16; n++) {
      p->m[m][n] = par[n];
    }
  }
  p->m[52][7]=-1;
  p->m[50][8]=-1;
  p->m[51][9]=-1;
  p->m[52][10]=-1;
  p->m[54][8]=p->m[54][9]=-1;
  p->m[55][8]=p->m[55][10]=-1;
  p->m[56][9]=p->m[56][10]=-1;
  g_niikmat_decompose_affine_obj_func_ipar=par;
  /* niikmat_display(p); */
  if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) obj func\n");
  pfn = niikmat_decompose_affine_obj_func;
  if(dof<12) {
    if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) start nelder-mead\n");
    if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&iter)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
      return 0;
    }
  } else {
    if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) start nelder-mead\n");
    if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
      return 0;
    }
  }
  if(verbose>1) fprintf(stdout,"-v (niikmat_decompose_affine) results\n");
  for(n=0; n<17; n++) {
    par[n] = p->m[0][n];
  }
  switch(dof) {
  case 6:
    par[7]=par[8]=par[9]=par[10]=1;
    par[11]=par[12]=par[13]=0;
    break;
  case 7:
    par[8]=par[9]=par[10]=1;
    par[11]=par[12]=par[13]=0;
    break;
  case 9:
    par[7]=1;
    par[11]=par[12]=par[13]=0;
    break;
  case 12:
    par[7]=1;
    break;
  }
  if(verbose>1) {
    if(!niik_aregister_display_affine(par)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_display_affine\n",__FILE__,__LINE__,__func__);
      return 0;
    }
  }
  if((matout = niik_aregister_matrix_from_affpar(par))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  /* niikmat_display(matout); */
  if(verbose) fprintf(stdout,"   error = %15.9f | iter = %i \n",tol/10000.0,iter);
  /* subroutine for translations */
  par[4]+=mat->m[0][3]-matout->m[0][3];
  par[5]+=mat->m[1][3]-matout->m[1][3];
  par[6]+=mat->m[2][3]-matout->m[2][3];
  if(!niik_aregister_matrix_from_affpar_update(matout,par)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(verbose) {
    niikmat_display(matout);
    fprintf(stdout,"      r  = %15.9f %15.9f %15.9f \n",par[1],par[2],par[3]); /* in degrees */
    fprintf(stdout,"      t  = %15.9f %15.9f %15.9f \n",par[ 4],par[ 5],par[ 6]);
    fprintf(stdout,"      s  = %15.9f %15.9f %15.9f   includes global scaling of %-15.9f\n",par[ 7]*par[8],par[7]*par[ 9],par[7]*par[10],par[7]);
    fprintf(stdout,"      sh = %15.9f %15.9f %15.9f \n",par[11],par[12],par[13]);
    fprintf(stdout,"      c  = %15.9f %15.9f %15.9f \n",par[14],par[15],par[16]);
  }
  if(!niik_aregister_matrix_from_affpar_update(matout,par)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(m=0; m<4; m++)
    for(n=0; n<4; n++)
      dout+=fabs(mat->m[m][n]-matout->m[m][n]);
  /* niikmat_display(matout); */
  if(verbose) {
    fprintf(stdout,"   error = %15.9f    iter %i\n",dout,iter);
    niikmat_display(matout);
  }
  p=niikmat_free(p);
  matout=niikmat_free(matout);
  return 1;
} /* int niikmat_decompose_affine(niikmat *mat,double *par,int dof) */

double niikmat_decompose_affine_obj_func(double *v) {
  niikmat *m;
  int
  i,j,
  verbose=0;
  double out=0,tmp[25];
  static double optout=1e9;
  static int iter=0;
  if(v==NULL) {
    optout=1e9;
    iter=0;
    return 0;
  }
  for(i=0; i<17; i++) tmp[i]=v[i];
  for(i=14; i<=16; i++) {
    v[i] = tmp[i] = g_niikmat_decompose_affine_obj_func_ipar[i];
  }
  switch(g_niikmat_decompose_affine_obj_func_dof) {
  case 6:
    tmp[7]=tmp[8]=tmp[9]=tmp[10]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 7:
    tmp[8]=tmp[9]=tmp[10]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 9:
    tmp[7]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 12:
    tmp[7]=1;
    break;
  default:
    fprintf(stderr,"[%s:%i:%s] ERROR: unknown dof %i\n",__FILE__,__LINE__,__func__,g_niikmat_decompose_affine_obj_func_dof);
    fprintf(stderr,"ERROR: during niikmat_decompose_affine_obj_func\n");
    exit(0);
  }
  if((m = niik_aregister_matrix_from_affpar(tmp))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    fprintf(stderr,"ERROR: during niikmat_decompose_affine_obj_func\n");
    exit(0);
  }
  iter++;
  if(fabs(tmp[1])>360) out+=fabs(tmp[1]-360);
  if(fabs(tmp[2])>360) out+=fabs(tmp[2]-360);
  if(fabs(tmp[3])>360) out+=fabs(tmp[3]-360);
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      out += 10000.0*fabs ( g_niikmat_decompose_affine_obj_func_mat->m[i][j] -
                            m->m[i][j]);
  if(verbose>1) fprintf(stdout,"  %15.9f %9i\n",out,iter);
  if(optout>out) {
    optout=out;
    if(verbose) {
      fprintf(stdout,"  %15.9f %9i *\n",out,iter);
      /*      niik_aregister_display_affine(tmp);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[0][0],m->m[0][1],m->m[0][2],m->m[0][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[1][0],m->m[1][1],m->m[1][2],m->m[1][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[2][0],m->m[2][1],m->m[2][2],m->m[2][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[3][0],m->m[3][1],m->m[3][2],m->m[3][3]); */
    }
  }
  m=niikmat_free(m);
  return out;
} /* niikmat_decompose_affine_obj_func */



/******* end of decompose affine matrix *******/







/*******************************************************
 *
 * AVERAGING MATRIX
 *
 * October 13, 2012
 * -removed gsl again, because gsl appears to be wrong and
 *  I can't correct it
 *
 *
 * June 27, 2012
 * -working and added gsl
 * -their eigenvector appear to be wrong...
 * -this section is incomplete... use matlab to check
 *  the results...
 *
 *******************************************************/

niikmat *niikmat_average_from_affine_param(niikmat **A,int num)
/* calculate average matrix from multiple matrices
 * using affine decomposition */
{
  double *p;
  niikmat *m;
  if(A==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: input list of matrices is null\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  fprintf(stdout,"[niikmat_average_from_affine_param] start\n");
  if((p=niik_average_affine_param(A,num))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_average_affine_param\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if((m=niik_aregister_matrix_from_affpar(p))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  fprintf(stdout,"[niikmat_average_from_affine_param] final matrix\n");
  niikmat_display(m);
  free(p);
  return m;
}

double *niik_average_affine_param(niikmat **A,int num)
/* average matrices using affine parameters */
{
  int m,n,dof=12;
  double **p,*out;
  if(A==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: input list of matrices is null\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  if(num==0) {
    fprintf(stderr,"ERROR: input list of matrices is null\n");
    return NULL;
  }
  for(n=0; n<num; n++) {
    if(A[n]->row!=4) {
      fprintf(stderr,"[%s:%i:%s] ERROR: matrix [%i] is not 4x4 matrix, %i x %i\n",__FILE__,__LINE__,__func__,n,A[n]->row,A[n]->col);
      return NULL;
    }
    if(A[n]->col!=4) {
      fprintf(stderr,"[%s:%i:%s] ERROR: matrix [%i] is not 4x4 matrix, %i x %i\n",__FILE__,__LINE__,__func__,n,A[n]->row,A[n]->col);
      return NULL;
    }
  }
  p=(double **)calloc(num,sizeof(double *));
  /*  #pragma omp parallel for private(m)*/
  for(n=0; n<num; n++) {
    p[n]=(double *)calloc(25,sizeof(double));
    for(m=7; m<=10; m++) p[n][m]=1;
    fprintf(stdout,"[niik_average_affine_param] matrix %i\n",n+1);
    if(!niikmat_decompose_affine(A[n],p[n],dof)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_decompose_affine\n",__FILE__,__LINE__,__func__);
      continue;
    }
    niik_aregister_display_affine( p[n] );
  }
  out=(double *)calloc(25,sizeof(double));
  for(m=0; m<18; m++) {
    for(n=0; n<num; n++) {
      out[m]+=p[n][m];
    }
    out[m]/=num;
  }
  for(n=0; n<num; n++) {
    free(p[n]);
  }
  free(p);
  fprintf(stdout,"[niik_average_affine_param] average\n");
  niik_aregister_display_affine( out );
  return out;
}

/******* end of averaging matrices *******/


/*******************************************************
 *
 * HALFWAY MATRIX
 *
 *******************************************************/

niikmat *niikmat_halfway_matrix2(niikmat *A) {
  niikmat *X=NULL,*Y=NULL,*iX=NULL,*iY=NULL,*XX;
  const double eps=1e-2;
  double d=100;
  int i,j;
  if(A==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: input matrix is null\n",__FILE__,__LINE__,__func__);
    return NULL;
  }
  X=niikmat_copy(A);
  Y=niikmat_identity(A->row,A->col);
  iX=niikmat_copy(A);
  iY=niikmat_copy(A);
  while(d>eps) {
    niikmat_copy_update(X,iX);
    niikmat_copy_update(Y,iY);
    niikmat_inverse_update(iX);
    niikmat_inverse_update(iY);
    for(i=0; i<A->row; i++) {
      for(j=0; j<A->col; j++) {
        X->m[i][j]=0.5*(X->m[i][j]+iY->m[i][j]);
        Y->m[i][j]=0.5*(Y->m[i][j]+iX->m[i][j]);
      }
    }
    XX=niikmat_multiply(X,X);
    for(i=0; i<A->row; i++) {
      for(j=0; j<A->col; j++) {
        XX->m[i][j]=XX->m[i][j]-A->m[i][j];
      }
    }
    d=niikmat_mag(XX);
    XX=niikmat_free(XX);
    fprintf(stdout,"[niikmat_halfway_matrix2] e %.9f\n",d);
  }
  Y=niikmat_free(Y);
  iX=niikmat_free(iX);
  iY=niikmat_free(iY);
  return X;
}


niikmat *g_niikmat_halfway_matrix_obj_func_mat;
int g_niikmat_halfway_matrix_obj_func_dof;
int g_niikmat_halfway_matrix_obj_func_3x3;

int niikmat_halfway_matrix(niikmat *afmat, double *afpar,int dof)
/* -calculates the halfway matrix
 * -returns the parameters to calculate the forward matrix (1/2 of afmat)
 * -inverse can be calculated by inverting the forward matrix
 * -afmat = inv(fwd) * fwd */
{
  double niikmat_halfway_matrix_obj_func(double *v);
  niikmat
  *tmpmat,
  *p;
  int
  verbose=1,
  m,n,
  iter=1e5,
  ndim=225;
  double tol=1e-5,(* pfn)()=NULL;

  if(verbose) fprintf(stdout,"[niikmat_halfway_matrix] start\n");
  if(afmat==NULL)   {
    fprintf(stderr,"[%s:%i:%s]  ERROR: afmat is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(afmat->col!=4) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: afmat->col is not 4, %i\n",__FILE__,__LINE__,__func__,afmat->col);
    return 0;
  }
  if(afmat->row!=4) {
    fprintf(stderr,"[%s:%i:%s]  ERROR: afmat->row is not 4, %i\n",__FILE__,__LINE__,__func__,afmat->row);
    return 0;
  }

  g_niikmat_halfway_matrix_obj_func_mat = afmat;
  p=niikmat_init(ndim+1,ndim);

  /* initialize the first pass
   * -optimize within 3x3 matrix
   */
  for(n=7; n<=10; n++) {
    p->m[0][n] = 1;
  }
  for(n=14; n<=16; n++) {
    p->m[0][n] = afpar[n];
  }
  for(n=1; n<=16; n++) {
    p->m[1][n] = 1;
  }
  for(m=2; m<=ndim; m++) {
    p->m[m][0]=0;
    for(n=1; n<=3; n++) {
      p->m[m][n] = niik_get_rand() * 60.0 - 30.0;
    }
    for(n=4; n<=6; n++) {
      p->m[m][n] = niik_get_rand() * 40 - 20;
    }
    for(n=7; n<=10; n++) {
      p->m[m][n] = niik_get_rand() * 0.4 + 0.8;
    }
    for(n=11; n<=13; n++) {
      p->m[m][n] = niik_get_rand() * 0.4 - 0.2;
    }
    for(n=14; n<=16; n++) {
      p->m[m][n] = afpar[n];
    }
  }
  pfn=niikmat_halfway_matrix_obj_func;
  g_niikmat_halfway_matrix_obj_func_dof = dof;
  g_niikmat_halfway_matrix_obj_func_3x3 = 1;

  if(verbose) fprintf(stdout,"[niikmat_halfway_matrix] first pass (3x3)\n");
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_RATIO,pfn,&iter)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
    return 0;
  }

  /* initialize for the second pass
   */
  if(verbose) fprintf(stdout,"[niikmat_halfway_matrix] second pass (translation)\n");
  g_niikmat_halfway_matrix_obj_func_3x3 = 0;
  tol = 1e-7;
  iter=1e4;
  for(m=2; m<=ndim; m++) {
    for(n=4; n<=6; n++) {
      p->m[m][n] = niik_get_rand() * 100.0 - 50.0;
    }
  }
  if(!niik_nelder_mead(p,ndim,&tol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
    return 0;
  }

  /* get the parameters */
  for(n=0; n<18; n++) {
    afpar[n]=p->m[0][n];
  }
  switch(dof) {
  case 6:
    afpar[7]=afpar[8]=afpar[9]=afpar[10]=1;
    afpar[11]=afpar[12]=afpar[13]=0;
    break;
  case 7:
    afpar[8]=afpar[9]=afpar[10]=1;
    afpar[11]=afpar[12]=afpar[13]=0;
    break;
  case 9:
    afpar[7]=1;
    afpar[11]=afpar[12]=afpar[13]=0;
    break;
  case 12:
    afpar[7]=1;
    break;
  default:
    return 0;
  }

  niikmat_display(afmat);
  niik_aregister_display_affine(afpar);

  tmpmat=niikmat_init(4,4);
  niik_aregister_matrix_from_affpar_update(tmpmat,afpar);

  if(verbose) {
    fprintf(stdout,"[niikmat_halfway_matrix] fwd halfway matrix\n");
    niikmat_display(tmpmat);
  }

  niikmat_multiply_mat1_free2(tmpmat,niikmat_copy(tmpmat));
  if(verbose) {
    fprintf(stdout,"[niikmat_halfway_matrix] reconstructed full matrix\n");
    niikmat_display(tmpmat);
  }
  tmpmat = niikmat_free(tmpmat);

  g_niikmat_halfway_matrix_obj_func_mat=NULL;
  return 1;
}

double niikmat_halfway_matrix_obj_func(double *v) {
  niikmat *m;
  int
  i,j,
  verbose=0;
  double out=0,tmp[25];
  static double optout=1e9;
  static int iter=0;
  if(v==NULL) {
    optout=1e9;
    iter=0;
    return 0;
  }
  for(i=0; i<17; i++) tmp[i]=v[i];
  switch(g_niikmat_halfway_matrix_obj_func_dof) {
  case 6:
    tmp[7]=tmp[8]=tmp[9]=tmp[10]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 7:
    tmp[8]=tmp[9]=tmp[10]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 9:
    tmp[7]=1;
    tmp[11]=tmp[12]=tmp[13]=0;
    break;
  case 12:
    tmp[7]=1;
    break;
  default:
    fprintf(stderr,"[%s:%i:%s] ERROR: unknown dof %i\n",__FILE__,__LINE__,__func__,g_niikmat_halfway_matrix_obj_func_dof);
    fprintf(stderr,"ERROR: during niikmat_decompose_affine_obj_func\n");
    exit(0);
  }
  if((m=niik_aregister_matrix_from_affpar(tmp))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_aregister_matrix_from_affpar\n",__FILE__,__LINE__,__func__);
    fprintf(stderr,"ERROR: during niikmat_decompose_affine_obj_func\n");
    exit(0);
  }
  iter++;
  if(fabs(tmp[1])>360) out+=fabs(tmp[1]-360);
  if(fabs(tmp[2])>360) out+=fabs(tmp[2]-360);
  if(fabs(tmp[3])>360) out+=fabs(tmp[3]-360);
  niikmat_multiply_mat1_free2(m,niikmat_copy(m));
  if(g_niikmat_halfway_matrix_obj_func_3x3) {
    for(i=0; i<3; i++)
      for(j=0; j<3; j++)
        out += 10000.0*fabs ( g_niikmat_halfway_matrix_obj_func_mat->m[i][j] -
                              m->m[i][j]);
  } else if(!g_niikmat_halfway_matrix_obj_func_3x3) {
    for(i=0; i<4; i++)
      out += fabs ( g_niikmat_halfway_matrix_obj_func_mat->m[i][3] -
                    m->m[i][3]);
  }
  if(verbose>1) fprintf(stdout,"  %15.9f %9i\n",out,iter);
  if(optout>out) {
    optout=out;
    if(verbose) {
      fprintf(stdout,"  %15.9f %9i *\n",out,iter);
      /*      niik_aregister_display_affine(tmp);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[0][0],m->m[0][1],m->m[0][2],m->m[0][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[1][0],m->m[1][1],m->m[1][2],m->m[1][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[2][0],m->m[2][1],m->m[2][2],m->m[2][3]);
      fprintf(stdout,"  %15.9f %15.9f %15.9f %15.9f\n",m->m[3][0],m->m[3][1],m->m[3][2],m->m[3][3]); */
    }
  }
  m=niikmat_free(m);
  return out;
} /* niikmat_halfway_matrix_obj_func */


/******* end of halfway matrix *******/




/******************************************
 *
 * slope calculation
 *
 ******************************************/

int niik_slope_intercept_from_double_vector(double *x, double *y, int num, double *slope, double *intercept) {
  double
  xy=0,
  ux=0,
  uy=0,
  x2=0;
  int n;
  for(n=0; n<num; n++) {
    xy += x[n]*y[n];
    ux += x[n];
    uy += y[n];
    x2 += x[n]*x[n];
  }
  xy/=num;
  ux/=num;
  uy/=num;
  x2/=num;
  *slope = (xy - ux*uy) / (x2 - ux*ux);
  *intercept = uy - *slope * ux;
  return 1;
}

int niik_slope_from_double_vector(double *x, double *y, int num, double *slope)
/* -calculates the slope of y - x
 * -to calculate the intercept, Intercept = Mean(y) - slope * Mean(x)
 */
{
  double
  xy=0,
  ux=0,
  uy=0,
  x2=0;
  int n;
  for(n=0; n<num; n++) {
    xy += x[n]*y[n];
    ux += x[n];
    uy += y[n];
    x2 += x[n]*x[n];
  }
  xy/=num;
  ux/=num;
  uy/=num;
  x2/=num;
  *slope = (xy - ux*uy) / (x2 - ux*ux);
  return 1;
} /* simple linear regression */


int niik_slope_from_double_vector_least_absolute_difference_cost_function_num=0;
int niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim=2;
double *niik_slope_from_double_vector_least_absolute_difference_cost_function_var[2];

double niik_slope_from_double_vector_least_absolute_difference_cost_function(double *p) {
  static double dd=1e30;
  int m,n;
  double d=0,y,z;
  if(p==NULL) {
    dd=1e30;
    return 0;
  }
  for(n=0; n<niik_slope_from_double_vector_least_absolute_difference_cost_function_num; n++) {
    for(m=0,y=0,z=1; m<niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim; m++) {
      y += z * p[m];
      z *= niik_slope_from_double_vector_least_absolute_difference_cost_function_var[0][n];
    }
    d += fabs(y-niik_slope_from_double_vector_least_absolute_difference_cost_function_var[1][n]);
  }
  d=d/niik_slope_from_double_vector_least_absolute_difference_cost_function_num;
  if(d<dd) {
    dd=d;
    fprintf(stdout,"%3i",niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim);
    for(m=0; m<niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim; m++) {
      fprintf(stdout,"%18.6f ",p[m]);
    }
    fprintf(stdout," -> %24.9f\n",d);
  }
  return d;
} /* niik_slope_from_double_vector_least_absolute_difference_cost_function */

int niik_slope_from_double_vector_least_absolute_difference(double *x,double *y,int num,double *par,int ndeg) {
  int m,n,ndim,iter;
  niikmat *p;
  double
  atol,
  (* pfn)();
  ndim=50+ndeg;
  p = niikmat_init(ndim+1,ndim);
  pfn = niik_slope_from_double_vector_least_absolute_difference_cost_function;
  atol=1e-4;
  niik_slope_from_double_vector_least_absolute_difference_cost_function_var[0]=x;
  niik_slope_from_double_vector_least_absolute_difference_cost_function_var[1]=y;
  niik_slope_from_double_vector_least_absolute_difference_cost_function_num=num;

  switch(ndeg) {
  default:
    for(m=0; m<ndim+1; m++) {
      for(n=0; n<ndeg; n++) {
        p->m[m][n]=10*niik_get_rand()-5;
      }
    }
    iter=1e4;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n");
    break;
  case 1:
    for(m=0; m<ndim+1; m++) {
      for(n=0; n<ndeg; n++) {
        p->m[m][n]=10*niik_get_rand()-5;
      }
    }
    iter=1e2;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);
    break;

  case 2:
    for(m=0; m<ndim+1; m++) {
      for(n=0; n<ndeg; n++) {
        p->m[m][n]=10*niik_get_rand()-5;
      }
    }
    iter=1e3;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"[%s:%i:%s] ERROR: nifti_k_nelder_mead\n",__FILE__,__LINE__,__func__);
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<ndeg; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);
    break;

  case 3:
    /* initial estimate with 2 parameters */
    for(m=0; m<ndim+1; m++) {
      for(n=0; n<2; n++) {
        p->m[m][n]=10*niik_get_rand()-5;
      }
    }
    iter=5e2;
    niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim=2;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);

    niik_slope_from_double_vector_least_absolute_difference_cost_function(NULL);
    for(m=0; m<ndim+1; m++) {
      p->m[m][n]=10*niik_get_rand()-5;
    }
    iter=5e3;
    niik_slope_from_double_vector_least_absolute_difference_cost_function_ndim=ndeg;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<ndeg; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);

    niik_slope_from_double_vector_least_absolute_difference_cost_function(NULL);
    for(m=0; m<ndim+1; m++) {
      p->m[m][n]*=1.05;
    }
    iter=5e3;
    atol*=0.1;
    if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_ABS,pfn,&iter)) {
      fprintf(stderr,"ERROR: nifti_k_nelder_mead\n");
      return 0;
    }
    fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
    for(n=0; n<ndeg; n++) {
      fprintf(stdout,"%12.6f ",par[n]);
    }
    fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);

    break;

  case 4:
  case 5:
    break;
  }

  fprintf(stdout,"[niik_slope_from_double_vector_least_abs_diff] parameters: ");
  for(n=0; n<ndeg; n++) {
    fprintf(stdout,"%12.6f ",par[n]);
  }
  fprintf(stdout,"\n[niik_slope_from_double_vector_least_abs_diff] iter %i\n",iter);

  p=niikmat_free(p);
  return 1;
}


/******** end of slope *******/






/******************************************
 *
 * polynomial fitting function
 *
 * -not validated, 2012-06-01, Kunio
 *
 ******************************************/

int niik_polynomial_fit(double *x, double *y, int num, int ndeg,double *par) {
  niikmat *mat,*imat;
  int m,n;
  if(x==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: x is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: y is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if((mat = niikmat_init(num,ndeg+1))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_init\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(m=0; m<mat->row; m++) {
    mat->m[m][0]=1;
    for(n=1; n<mat->col; n++) {
      mat->m[m][n]=mat->m[m][n-1]*x[n];
    }
  }
  if((imat = niikmat_pseudo_inverse(mat))==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikmat_pseudo_inverse\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  for(m=0; m<mat->row; m++) {
    par[m]=0;
    for(n=0; n<mat->col; n++) {
      par[m]+=imat->m[n][m]*y[m];
    }
  }
  mat=niikmat_free(mat);
  imat=niikmat_free(imat);
  return 1;
} /* niik_polynomial_fit */




/******************************************
 *
 * correlation coefficient function
 *
 *
 ******************************************/

int niik_double_vector_calc_corrcoef(double *v,double *w,int num,double *out) {
  niikvec *V,*W;
  int n;
  V=niikvec_init(num);
  W=niikvec_init(num);
  for(n=0; n<num; n++) {
    V->v[n]=v[n];
    W->v[n]=w[n];
  }
  if(!niikvec_calc_corrcoef(V,W,out)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niikvec_calc_corrcoef\n",__FILE__,__LINE__,__func__);
    V=niikvec_free(V);
    W=niikvec_free(W);
    return 0;
  }
  V=niikvec_free(V);
  W=niikvec_free(W);
  return 1;
} /* niik_double_vector_calc_corrcoef */

int niikvec_calc_corrcoef(niikvec *v,niikvec *w,double *out) {
  int n,num;
  double
  xsum,ysum,xssq,yssq,xysum,wsum;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(w==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: w is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(v->num!=w->num) {
    fprintf(stderr,"[%s:%i:%s] ERROR: num is different %i %i\n",__FILE__,__LINE__,__func__,v->num,w->num);
    return 0;
  }
  num=v->num;
  xsum=ysum=xssq=yssq=xysum=wsum=0;
  for(n=0; n<num; n++) {
    xsum  += v->v[n];
    ysum  += w->v[n];
    xssq  += v->v[n] * v->v[n];
    yssq  += w->v[n] * w->v[n];
    xysum += v->v[n] * w->v[n];
    wsum  += 1.0;
  }
  *out = (wsum * xysum - xsum * ysum) / sqrt(wsum*xssq-xsum*xsum) / sqrt(wsum*yssq-ysum*ysum);
  return 1;
} /* niikvec_calc_corrcoef */


int niikvec_cross_correlation(niikvec *v,niikvec *w,niikvec *out) {
  int m,n,nm;
  if(v==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: v is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(w==NULL) {
    fprintf(stderr,"[%s:%i:%s] ERROR: w is null\n",__FILE__,__LINE__,__func__);
    return 0;
  }
  if(v->num!=w->num) {
    fprintf(stderr,"[%s:%i:%s] ERROR: num is different %i %i\n",__FILE__,__LINE__,__func__,v->num,w->num);
    return 0;
  }
  if(v->num*2+1!=out->num) {
    fprintf(stderr,"[%s:%i:%s] ERROR: num is incorrect %i %i\n",__FILE__,__LINE__,__func__,v->num,out->num);
    return 0;
  }
  for(n=0; n<out->num; n++) {
    out->v[out->num-n-1]=0;
    for(m=0; m<v->num; m++) {
      nm=n+m-(out->num-1)/2;
      if(nm<0) continue;
      if(nm>=w->num) continue;
      /*fprintf(stdout,"%3i f(%3i) %12.4f  g(%3i) %12.4f\n",n,m,v->v[m],nm,w->v[nm]);*/
      out->v[out->num-n-1]+=v->v[m]*w->v[nm];
    }
  }
  return 1;
}


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/