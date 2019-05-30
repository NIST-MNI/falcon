/* Filename:     nifti1_kunio_nelder_mead.c
 * Description:  nelder-mead simplex downhill optimization function
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 */

#ifndef _FALCON_NELDER_MEAD_C_
#define _FALCON_NELDER_MEAD_C_

#include "falcon.h"

#define NIIK_NM_SWAP(a,b) { swap=(a); (a)=(b); (b)=swap; }

/* local functions */
int niik_nelder_mead_order(niikmat *p,double *y);


/* global variables */
double alpha = 1.0;
double phi   = 2.0;
double sigma = 0.5;

/**************************************************************************
 *
 * int niik_nelder_mead(niikmat *p,int ndim,double tol,double (* pfn)(),int maxiter);
 *
 * -main nelder-mead function
 * -p is a matrix of parameters (ndim+1-by-ndim matrix)
 * -ndim is the number of variables
 * -tol is the tolerance, replaced with the minimized cost
 * -pfn is the pointer to the function
 * -maxiter is the maximum iteration, replaced with actual iteration
 * -cost_method is the cost method
 *   1 = 2 * abs(min_cost-max_cost)/(min_cost+max_cost)
 *   2 = min_cost
 *
 **************************************************************************/

int niik_nelder_mead(niikmat *p,int ndim,double *tol,int cost_method,double (* pfn)(),int *maxiter) {
  double *y,*g,*r,*e,*c,yr,ye,yc,se;
  int m,n,ndim1,iter;
  int verbose=0;

  if(verbose>=1) {
    fprintf(stdout,"-v%i (niik_nelder_mead) start\n",verbose);
    fprintf(stdout,"  ndim = %-i\n",ndim);
    fprintf(stdout,"  tol  = %0.1e\n",*tol);
    fprintf(stdout,"  iter = %-i\n",*maxiter);
  }

  if(verbose>=1) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation\n");
  ndim1 = ndim + 1;
  if(verbose>=3) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation y %i\n",ndim1);
  y = (double *)calloc(ndim1,sizeof(double));
  if(verbose>=3) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation g\n");
  g = (double *)calloc(ndim,  sizeof(double));
  if(verbose>=3) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation r\n");
  r = (double *)calloc(ndim,  sizeof(double));
  if(verbose>=3) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation r\n");
  e = (double *)calloc(ndim,  sizeof(double));
  if(verbose>=3) fprintf(stdout,"-v  (niik_nelder_mead) memory allocation c\n");
  c = (double *)calloc(ndim,  sizeof(double));

  if(verbose>=1) fprintf(stdout,"-v  (niik_nelder_mead) initial errors\n");
  for(n=0; n<ndim1; n++) {
    if(verbose>=2) fprintf(stdout,"-v  (niik_nelder_mead) initial error %i\n",n);
    y[n] = pfn(p->m[n]); /* initial estimates */
    if(verbose>=1) fprintf(stdout,"-v  (niik_nelder_mead) initial error %i %.7f\n",n,y[n]);
  }
  niik_nelder_mead_order(p,y);

  if(verbose>=2) {
    fprintf(stdout,"-v%i (niik_nelder_mead) initialization\n",verbose);
    for(n=0; n<ndim1; n++) {
      fprintf(stdout,"  %12.6f | ",y[n]);
      for(m=0; m<ndim; m++) {
        fprintf(stdout,"%7.3f ",p->m[n][m]);
      }
      fprintf(stdout,"\n");
    }
  }

  for(iter=1; iter<=*maxiter; iter++) {

    switch(cost_method) {
    default:
    case NIIK_NELDER_MEAD_COST_RATIO:
      se = 2.0*fabs(y[0]-y[ndim])/(fabs(y[0])+fabs(y[ndim])); /* similar to NR */
      if(verbose>=2)
        fprintf(stdout,"  se = %19.9f   %19.9f %19.9f\n",se,y[0],y[ndim]);
      break;
    case NIIK_NELDER_MEAD_COST_ABS:
      se = y[0];
      if(verbose>=2)
        fprintf(stdout,"  se = %19.9f   %19.9f\n",se,y[0]);
      break;
    case NIIK_NELDER_MEAD_COST_DIFF:
      /* Kunio, 2012-10-20
       * -simple difference in case error is close to zero */
      se = fabs(y[0]-y[ndim]);
      break;
    }

    /* finish iterations */
    if(se < *tol) {
      break;
    }

    niik_nelder_mead_order(p,y);
    if(verbose>2) {
      fprintf(stdout,"w %12.6f | ",y[ndim]);
      for(m=0; m<ndim; m++) {
        fprintf(stdout,"%7.3f ",p->m[ndim][m]);
      }
      fprintf(stdout,"\n");
    }

    /* center of gravity */
    for(n=0; n<ndim; n++) {
      g[n]=0;
      for(m=0; m<ndim; m++) {
        g[n] += p->m[m][n];
      }
      g[n] /= ndim;
    }

    if(verbose>2) {
      fprintf(stdout,"g %12.6f | ",0*y[ndim]);
      for(m=0; m<ndim; m++) {
        fprintf(stdout,"%7.3f ",g[m]);
      }
      fprintf(stdout,"\n");
    }

    /* try reflection */
    if(verbose>=2) fprintf(stdout,"-v%i (niik_nelder_mead) reflection\n",verbose);
    for(n=0; n<ndim; n++) {
      r[n] = g[n] + alpha * (g[n] - p->m[ndim][n]);
    }
    yr = pfn(r);
    if(verbose>=3) {
      fprintf(stdout,"R %12.6f | ",yr);
      for(m=0; m<ndim; m++) {
        fprintf(stdout,"%7.3f ",r[m]);
      }
      fprintf(stdout,"\n");
    }

    /* reflection is better than the best */
    if(yr<y[0]) {
      /* try expansion */
      if(verbose>=2) fprintf(stdout,"-v%i (niik_nelder_mead) expansion\n",verbose);
      for(n=0; n<ndim; n++) {
        e[n] = g[n] + phi * (g[n] - p->m[ndim][n]);
      }
      ye = pfn(e);
      if(verbose>=3) {
        fprintf(stdout,"E %12.6f | ",ye);
        for(m=0; m<ndim; m++) {
          fprintf(stdout,"%7.3f ",e[m]);
        }
        fprintf(stdout,"\n");
      }
      if(ye<y[0]) { /* expansion was even better! */
        y[ndim]=ye;
        for(n=0; n<ndim; n++) p->m[ndim][n]=e[n];
      } else { /* expansion was not better */
        y[ndim] = yr;
        for(n=0; n<ndim; n++) p->m[ndim][n]=r[n];
      }
      continue;
    }

    /* reflection was better than the second worse */
    else if(yr<y[ndim-1]) {
      y[ndim] = yr;
      for(n=0; n<ndim; n++) p->m[ndim][n] = r[n];
      continue;
    }

    /* contraction --reflection was not better than the second worst */
    else {
      if(verbose>=2) fprintf(stdout,"-v%i (niik_nelder_mead) contraction\n",verbose);
      for(n=0; n<ndim; n++) {
        c[n] = p->m[ndim][n] + sigma * (g[n] - p->m[ndim][n]);
      }
      yc = pfn(c);
      if(verbose>=3) {
        fprintf(stdout,"C %12.6f | ",yc);
        for(m=0; m<ndim; m++) {
          fprintf(stdout,"%7.3f ",c[m]);
        }
        fprintf(stdout,"\n");
      }
      if(yc<y[ndim]) {
        y[ndim]=yc;
        for(n=0; n<ndim; n++) p->m[ndim][n]=c[n];
        continue;
      }
    }

    /* contraction was not better -> reduction */
    if(verbose>=2) fprintf(stdout,"-v%i (niik_nelder_mead) reduction\n",verbose);
    for(m=1; m<ndim1; m++) {
      for(n=0; n<ndim; n++) {
        p->m[m][n] = p->m[0][n] + sigma * (p->m[m][n] - p->m[0][n]);
      }
      y[m] = pfn(p->m[m]);
      if(verbose>=3) {
        fprintf(stdout,"c %12.6f | ",y[m]);
        for(n=0; n<ndim; n++) {
          fprintf(stdout,"%7.3f ",p->m[m][n]);
        }
        fprintf(stdout,"\n");
      }
    }

  } /* iteration */

  niik_nelder_mead_order(p,y);

  *maxiter = iter;
  *tol = y[0];

  free(c);
  free(e);
  free(r);
  free(g);
  free(y);

  return 1;
}


/**************************************************************************
 *
 * int niik_nelder_mead_order(niikmat *p,double *y,int ndim);
 *
 * -re-orders the parameters (p) and cost values (y)
 *  from min [0] to max [ndim]
 *
 **************************************************************************/

int niik_nelder_mead_order(niikmat *p,double *y) {
  int m,n,k,ndim,ndim1;
  double swap;
  if(p==NULL) {
    fprintf(stderr,"[niik_nelder_mead_order] ERROR: p is null\n");
    return 0;
  }
  if(y==NULL) {
    fprintf(stderr,"[niik_nelder_mead_order] ERROR: y is null\n");
    return 0;
  }
  ndim  = p->col;
  ndim1 = p->row;
  for(m=0; m<ndim1; m++) {
    for(n=m+1; n<ndim1; n++) {
      if(y[m]>y[n]) {
        NIIK_NM_SWAP(y[m],y[n]);
        for(k=0; k<ndim; k++) NIIK_NM_SWAP(p->m[m][k],p->m[n][k]);
      }
    }
  }
  return 1;
} /* niik_nelder_mead_order */


int niik_nelder_mead_multi_level(niikmat *s,niikmat *p,int *nseed,int nlevel,int ndim,double *tol,int cost_method,double (* pfn)(),int *maxiter,int verbose)
/* -multi-level / multi-seed optimization function
 * -s   : list of seed points for the first level
 *        nseed-by-ndim matrix
nn * -p   : perturbation at each level
 *        nlevel-by-ndim matrix
 * -nseed  : list of number of seed points
 *           nlevel-long vector
 * -nlevel : number of levels
 * -ndim   : number of parameters
 * -tol    : list of tolerance
 *           nlevel-long vector
 * -cost_method  : method for cost calculation
 * -pfn          : pointer to the function
 * -maxiter      : list of maximum iterations
 *                 nlevel-long vector
 */
{
  niikmat *t=NULL;
  niikvec *v=NULL;
  int
  i,j,
  nl,ns,
  srow,
  cmaxiter,amaxiter;
  double ctol,atol;
  char fcname[64]="niik_nelder_mead_multi_level";
  if(verbose>=1) niik_fc_display(fcname,1);
  if(s==NULL) {
    fprintf(stderr,"[%s] ERROR: s is null\n",fcname);
    return 0;
  }
  if(p==NULL) {
    fprintf(stderr,"[%s] ERROR: p is null\n",fcname);
    return 0;
  }
  if(nseed==NULL) {
    fprintf(stderr,"[%s] ERROR: nseed is null\n",fcname);
    return 0;
  }
  if(tol==NULL) {
    fprintf(stderr,"[%s] ERROR: tol is null\n",fcname);
    return 0;
  }
  if(pfn==NULL) {
    fprintf(stderr,"[%s] ERROR: pfn is null\n",fcname);
    return 0;
  }
  if(maxiter==NULL) {
    fprintf(stderr,"[%s] ERROR: maxiter is null\n",fcname);
    return 0;
  }
  if(s->row<nseed[0]) {
    fprintf(stderr,"[%s] ERROR: nseed is too large %i (s) < %i\n",fcname,s->row,nseed[0]);
    return 0;
  }
  if(p->row!=nlevel) {
    fprintf(stderr,"[%s] ERROR: nlevel did not match %i (p), %i\n",fcname,p->row,nlevel);
    return 0;
  }
  if(p->col!=ndim) {
    fprintf(stderr,"[%s] ERROR: ndim did not match %i (p), %i\n",fcname,p->row,ndim);
    return 0;
  }
  if(verbose>=1) {
    fprintf(stdout,"[%s] # level        : %i\n",fcname,nlevel);
    fprintf(stdout,"[%s] # param        : %i\n",fcname,ndim);
    fprintf(stdout,"[%s] # seed         : ",fcname);
    for(nl=0; nl<nlevel; nl++) {
      fprintf(stdout,"%i ",nseed[nl]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"[%s] tol            : ",fcname);
    for(nl=0; nl<nlevel; nl++) {
      fprintf(stdout,"%.3g ",tol[nl]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"[%s] max iter       : ",fcname);
    for(nl=0; nl<nlevel; nl++) {
      fprintf(stdout,"%i ",maxiter[nl]);
    }
    fprintf(stdout,"\n");
  }
  t=niikmat_init(ndim+1,ndim);
  v=niikvec_init(nseed[0]);
  srow=s->row;
  for(nl=0; nl<nlevel; nl++) {
    /*if(verbose>=1) fprintf(stdout,"[%s] level %i\n",fcname,nl);*/
    for(ns=0; ns<nseed[nl]; ns++) {
      for(i=0; i<ndim+1; i++) {
        for(j=0; j<ndim; j++) {
          t->m[i][j] = s->m[ns][j] + p->m[nl][j] * niik_get_rand();
        }
      }
      if(verbose>=3) {
        fprintf(stdout,"[%s] initial parameters\n",fcname);
        niikmat_display(t);
        fprintf(stdout,"\n");
      }
      /* First pass */
      cmaxiter=maxiter[nl];
      ctol=tol[nl];
      if(verbose>=1) fprintf(stdout,"[%s] level %-3i seed %-5i\n",fcname,nl,ns);
      if(!niik_nelder_mead(t,ndim,&ctol,cost_method,pfn,&cmaxiter)) {
        fprintf(stderr,"[%s] ERROR: niik_nelder_mead\n",fcname);
        return 0;
      }
      if(verbose>=2) {
        fprintf(stdout,"[%s] improved  %20.8f  %9i    ",fcname,ctol,cmaxiter);
        niik_display_double_vector(t->m[0],ndim);
      }
      amaxiter=cmaxiter;
      atol=ctol;
      /* Second pass */
      cmaxiter=maxiter[nl];
      ctol=tol[nl];
      for(i=1; i<ndim; i++) {
        for(j=0; j<ndim; j++) {
          t->m[i][j] = s->m[ns][j] + p->m[nl][j] * niik_get_rand();
        }
      }
      for(j=0; j<ndim; j++) {
        t->m[ndim][j] = s->m[ns][j];
      }
      if(!niik_nelder_mead(t,ndim,&ctol,cost_method,pfn,&cmaxiter)) {
        fprintf(stderr,"[%s] ERROR: niik_nelder_mead\n",fcname);
        return 0;
      }
      if(verbose>=2) {
        fprintf(stdout,"[%s] optimized %20.8f  %9i    ",fcname,ctol,cmaxiter);
        niik_display_double_vector(t->m[0],ndim);
      }
      cmaxiter+=amaxiter;
      ctol=atol;
      /* sort */
      v->v[ns]=ctol;
      for(j=0; j<ndim; j++)
        s->m[ns][j]=t->m[0][j];
      s->row=ns+1;
      niik_nelder_mead_order(s,v->v);
      s->row=srow;
      /*for(n=0;n<=ns;n++){
        fprintf(stdout,"%5i %20.5f : ",n,v->v[n]);
        niik_display_double_vector(s->m[n],s->col); }
        fprintf(stdout,"\n\n");*/
    } /* each seed */
    if(verbose>=2) {
      for(ns=0; ns<nseed[nl]; ns++) {
        fprintf(stdout,"%5i %20.5f : ",ns,v->v[ns]);
        niik_display_double_vector(s->m[ns],s->col);
      }
      fprintf(stdout,"\n\n");
    }
  } /* level */
  if(verbose>=1) niik_fc_display(fcname,0);
  tol[0]=v->v[0];
  t=niikmat_free(t);
  v=niikvec_free(v);
  return 1;
} /* niik_nelder_mead_multi_level */


#endif /* _FALCON_NELDER_MEAD_C_ */
