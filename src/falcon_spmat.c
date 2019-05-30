/* Filename:     nifti1_kunio_spmat.c
 * Description:  functions for sparse matrix
 * Author:       Kunio Nakamura
 * Date:         March 4, 2012
 *
 * -Functions for sparse matrix using the new Yale Sparse Matrix Format.
 *  It is the one used in Numerical Recipe.
 * -Conjugate gradient method is here.
 * -So it can solve system of linear equations like the Laplace equations.
 */
#include "falcon.h"


int niik_spmat_construct(double **A,int N,double *sa,int *ija)
/* A is a N-by-N matrix
 * sa is replaced with new values
 * ija is alo replaced
 * This is the new Yale Sparse Matrix System Format
 * -It is the same as in Numerical Recipe */
{
  int i,j,m;
  /* diagonal terms */
  for(i=0; i<N; i++) {
    sa[i]=A[i][i];
  }
  /* off-diagonal terms */
  for(i=0,m=N+1; i<N; i++) {
    ija[i]=m;
    for(j=0; j<N; j++) {
      if(i==j) continue; /* skip diagonal terms */
      else if(fabs(A[i][j])<=1e-10) continue; /* skip zero terms */
      sa[m]=A[i][j];
      ija[m]=j;
      m++;
    }
  }
  ija[i]=m;
  return 1;
}



int niik_spmat_reverse(double **A,int N,double *sa,int *ija)
/* -instead of creating sa and ija from A
 *  this function creates A from sa and ija.
 * -it's a reverse
 * -it's primary purpose is to check sa and ija
 *  created directly from data (not through matrix A */
{
  int i,j,m;
  /* diagonal terms */
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      A[i][j]=0;
  for(i=0; i<N; i++) {
    A[i][i]=sa[i];
  }
  /* off-diagonal terms */
  for(i=0,m=N+1; i<N; i++) {
    while(m<ija[i+1]) {
      A[i][ija[m]]=sa[m];
      m++;
    }
  }
  return 1;
}



int niik_spmat_op_mul_mat_vec(double *sa,int *ija,int N,double *x,double *y)
/* -multiple sparse matrix and vector
 *  y = Ax
 *  where A is represnted by sa and ija
 * -the output is replaced in y
 * -all of these variables need to have memory allocated */
{
  int i,j;
  for(i=0,j=N+1; i<N; i++)  {
    y[i]=sa[i]*x[i];
    /* fprintf(stdout,"%i %10.2f  = %10.2f x %10.2f   j=%i  d\n",i,y[i],sa[i],x[i],j); */
    while(j<ija[i+1]) {
      y[i]+=sa[j]*x[ija[j]];
      /* fprintf(stdout,"%i %10.2f += %10.2f x %10.2f   j=%i\n",i,y[i],sa[j],x[ija[j]],j); */
      j++;
    }
  }
  return 1;
}



double niik_spmat_op_vec_dot_product(double *a,double *b,int N)
/* -computes the square of vectors
 *  i.e., output_scaler = transpose(a)b */
{
  double dval=0.0;
  int n;
  for(n=0; n<N; n++) {
    dval+=a[n]*b[n];
  }
  return dval;
}


int niik_spmat_conj_grad(double *sa,int *ija,int N,double *x,double *b)
/* -main conjugate gradient method using sparse matrix A (sa and ija)
 * -given b, initial x, and matrix A in sparse format (sa and ija)
 *  updates x iteratively
 */
{
  int i,iter;
  double
  rs1,rs2,alpha,beta,
      *a,*p,*r;
  int verbose=niik_verbose();
  char fcname[32]="niik_spmat_conj_grad";

  if(verbose>0) niik_fc_display(fcname,1);

  NIIK_RET0(((r=(double *)calloc(N,sizeof(double)))==NULL),fcname,"memory allocation for r");
  NIIK_RET0(((p=(double *)calloc(N,sizeof(double)))==NULL),fcname,"memory allocation for p");
  NIIK_RET0(((a=(double *)calloc(N,sizeof(double)))==NULL),fcname,"memory allocation for a");
  niik_spmat_op_mul_mat_vec(sa,ija,N,x,r);

  for(i=0; i<N; i++) r[i]=p[i]=b[i]-r[i];

  if(verbose>1) {
    fprintf(stdout,"r0  = ");
    niik_display_double_vector(r,N);
    fprintf(stdout,"x0  = ");
    niik_display_double_vector(x,N);
  }

  rs1=niik_spmat_op_vec_dot_product(r,r,N);
  if(verbose>1) {
    fprintf(stdout,"rs1 = %12.6f\n",rs1);
  }

  for(iter=1; iter<N; iter++) {

    niik_spmat_op_mul_mat_vec(sa,ija,N,p,a);
    alpha = rs1 / niik_spmat_op_vec_dot_product(p,a,N);
    if(verbose>1) {
      fprintf(stdout,"rs1 = %12.6f \n",rs1);
      fprintf(stdout,"pap = %12.6f \n",niik_spmat_op_vec_dot_product(p,a,N));
      fprintf(stdout,"al  = %12.6f \n",alpha);
    }

    /* update x */
    for(i=0; i<N; i++) {
      x[i] = x[i] + alpha*p[i];
    }
    if(verbose>1) {
      fprintf(stdout,"x%i  = ",iter);
      niik_display_double_vector(x,N);
    }

    /* update r */
    for(i=0; i<N; i++) {
      r[i] = r[i] - alpha*a[i];
    }
    if(verbose>1) {
      fprintf(stdout,"r%i  = ",iter);
      niik_display_double_vector(r,N);
    }

    rs2 = niik_spmat_op_vec_dot_product(r,r,N);
    beta=rs2/rs1;
    if(verbose>1) {
      fprintf(stdout,"rs2  = %12.6f \n",rs2);
      fprintf(stdout,"bet  = %12.6f \n",beta);
    }

    /* update p */
    for(i=0; i<N; i++) {
      p[i] = r[i] + beta*p[i];
    }
    if( verbose>1 ) {
      fprintf(stdout,"p%i  = ",iter);
      niik_display_double_vector(p,N);
    }

    if(sqrt(rs2)<1e-10 ) {
      if(verbose>1)        fprintf(stdout,"[niik_spmat_conj_grad] ConjGrad Converge:  %6.3e\n",sqrt(rs2));
      break;
    }

    if( verbose>1 )        fprintf(stdout,"[niik_spmat_conj_grad] ConjGrad Iter %4i: %6.3e\n",iter,sqrt(rs2));
    rs1=rs2;

    if( verbose>1 )        fprintf(stdout,"\n");
  }

  free(r);
  free(p);
  free(a);
  if(verbose>0) niik_fc_display(fcname,0);
  return 1;
} /* niik_spmat_conj_grad */


void niik_spmat_display_dig(double **A,int n1,int n2,int dig) {
  int i,j;
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {
      switch(dig) {
      case 0:
        fprintf(stdout,"%3.0f ",A[i][j]);
        break;
      case 1:
        fprintf(stdout,"%4.1f ",A[i][j]);
        break;
      case 2:
        fprintf(stdout,"%6.2f ",A[i][j]);
        break;
      default:
      case 3:
        fprintf(stdout,"%7.3f ",A[i][j]);
        break;
      case 4:
        fprintf(stdout,"%8.4f ",A[i][j]);
        break;
      case 5:
        fprintf(stdout,"%9.5f ",A[i][j]);
        break;
      case 6:
        fprintf(stdout,"%10.6f ",A[i][j]);
        break;
      case 7:
        fprintf(stdout,"%11.7f ",A[i][j]);
        break;
      case 8:
        fprintf(stdout,"%12.8f ",A[i][j]);
        break;
      case 9:
        fprintf(stdout,"%13.9f ",A[i][j]);
        break;
      }
    }
    fprintf(stdout,"\n");
  }
  return;
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/