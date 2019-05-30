/* Filename:     nifti1_kunio_curve_fitting.c
 * Description:  basic linear algebra functions by Kunio
 *               for curve fitting
 * Author:       Kunio Nakamura
 * Date:         October 13, 2012
 * DatParent:    nifti1_kunio_bla.c
 */

#ifndef _FALCON_CURVE_FITTING_C_
#define _FALCON_CURVE_FITTING_C_

#include "falcon.h"

int     niik_curve_fitting_method=NIIK_CURVE_FIT_SQUARE;
int     niik_curve_fitting_least_absolute_difference_cost_function_num=0;
int     niik_curve_fitting_least_absolute_difference_cost_function_ndim=2;
double *niik_curve_fitting_least_absolute_difference_cost_function_var[2];

double niik_curve_fitting_least_absolute_difference_cost_function(double *p) {
  static double dd=1e30;
  static int iter=0;
  int m,n,num,ndim,verbose=0;
  double d=0,y,z,*x,*v;
  if(p==NULL) {
    dd=1e30;
    iter=0;
    return 0;
  }
  num =niik_curve_fitting_least_absolute_difference_cost_function_num;
  ndim=niik_curve_fitting_least_absolute_difference_cost_function_ndim;
  x   =niik_curve_fitting_least_absolute_difference_cost_function_var[0];
  v   =niik_curve_fitting_least_absolute_difference_cost_function_var[1];
  if(verbose>=1) fprintf(stdout,"\npar %12.9f %12.9f %12.9f\n",p[0],p[1],p[2]);
  for(n=0; n<num; n++) {
    for(m=0,y=0,z=1; m<ndim; m++) {
      y += z*p[m];
      z *= x[n];
    }
    switch(niik_curve_fitting_method) {
    default:
    case NIIK_CURVE_FIT_UNKNOWN:
    case NIIK_CURVE_FIT_SQUARE:
      d += NIIK_SQ(y-v[n]);
      break;
    case NIIK_CURVE_FIT_ABSOLUTE:
      d += fabs(y-v[n]);
      break;
    }
    if(verbose>=1) fprintf(stdout,"\t\t%i %5.1f %12.9f %12.9f  %12.9f\n",n,x[n],v[n],y,fabs(y-v[n]));
  }
  d=d/num;
  iter++;
  if(d<dd) {
    dd=d;
    if(verbose>=1) {
      fprintf(stdout,"%3i",niik_curve_fitting_least_absolute_difference_cost_function_ndim);
      for(m=0; m<ndim; m++) {
        fprintf(stdout,"%18.6f ",p[m]);
      }
      fprintf(stdout," -> %24.9f  %12i\n",d,iter);
    }
  }
  return d;
} /* niik_curve_fitting_least_absolute_difference_cost_function */

int niik_curve_fitting_polynomial(double *x,double *y,int num,double *par,int ndeg,int method) {
  int m,n,ndim,iter;
  niikmat *p;
  double
  ymean,
  atol,
  (* pfn)();
  char fcname[64]="niik_curve_fittin_least_absolute_difference";

  if(x==NULL) {
    fprintf(stderr,"[%s] ERROR: x is null\n",fcname);
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[%s] ERROR: y is null\n",fcname);
    return 0;
  }
  if(par==NULL) {
    fprintf(stderr,"[%s] ERROR: par is null\n",fcname);
    return 0;
  }

  niik_fc_display(fcname,1);

  /* initialize optimization parameters */
  ndim=5+ndeg;
  p = niikmat_init(ndim+1,ndim);
  pfn = niik_curve_fitting_least_absolute_difference_cost_function;
  niik_curve_fitting_least_absolute_difference_cost_function_var[0]=x;
  niik_curve_fitting_least_absolute_difference_cost_function_var[1]=y;
  niik_curve_fitting_least_absolute_difference_cost_function_num=num;
  niik_curve_fitting_least_absolute_difference_cost_function_ndim=ndeg;
  niik_curve_fitting_method=method;
  ymean=niik_get_mean_from_double_vector(y,num);
  fprintf(stdout,"\tymean = %12.9f\n",ymean);

  for(m=0; m<ndim+1; m++) {
    p->m[m][0]=ymean + ymean * (0.5*niik_get_rand()-0.25);
    for(n=1; n<ndeg; n++) {
      p->m[m][n]=10*niik_get_rand()-5;
      p->m[m][n]/=pow(10.0,(double)n);
    }
  }
  p->col=ndeg;
  niikmat_display(p);
  p->col=ndim;
  iter=1e6;
  atol=1e-20;
  if(!niik_nelder_mead(p,ndim,&atol,NIIK_NELDER_MEAD_COST_DIFF,pfn,&iter)) {
    fprintf(stderr,"[%s] ERROR: nifti_k_nelder_mead\n",fcname);
    return 0;
  }
  fprintf(stdout,"[%s] error: %12.9f\n",fcname,atol);
  for(n=0; n<ndeg; n++) {
    par[n]=p->m[0][n];
  }

  fprintf(stdout,"[%s] parameters: ",fcname);
  for(n=0; n<ndeg; n++) {
    fprintf(stdout,"%12.6f ",par[n]);
  }
  fprintf(stdout,"\n[%s] iter %i\n",fcname,iter);

  p=niikmat_free(p);
  niik_fc_display(fcname,0);
  return 1;
}

#endif /* _FALCON_CURVE_FITTING_C_ */
