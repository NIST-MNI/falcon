#ifndef _FALCON_ANN_H_
#define _FALCON_ANN_H_

typedef struct {
  niikmat *alpha;     // ni x nh matrix
  niikmat *beta;      // nh x nk matrix
  int num;            // # data
  int ni;             // # inputs
  int nh;             // # hidden node
  int nk;             // # outputs
  niikmat *in;        // input data, num * ni
  niikmat *out;       // output data, num * nk
  int verbose;
} niik_simple_ann; // class

niik_simple_ann *niik_simple_ann_init() {
  niik_simple_ann *a=NULL;
  NIIK_RET0(((a=(niik_simple_ann *)calloc(1,sizeof(niik_simple_ann)))==NULL),__func__,"calloc");
  a->alpha=NULL;
  a->beta=NULL;
  a->in=NULL;
  a->out=NULL;
  a->num=a->ni=a->nh=a->nk=0;
  a->verbose=1;
  return a;
}

niik_simple_ann *niik_simple_ann_init_num(int num,int ni,int nh,int nk) {
  niik_simple_ann *a=NULL;
  NIIK_RET0(((a=(niik_simple_ann *)calloc(1,sizeof(niik_simple_ann)))==NULL),__func__,"calloc");
  a->alpha=niikmat_init(ni,nh);
  a->beta=niikmat_init(nh,nk);
  a->in=niikmat_init(num,ni);
  a->out=niikmat_init(num,nk);
  a->verbose=1;
  return a;
}

niik_simple_ann *niik_simple_ann_data_from_fann_file(char *fname) {
  niik_simple_ann *a=NULL;
  FILE *fp;
  int n,i,k;
  NIIK_RET0(((a=(niik_simple_ann *)calloc(1,sizeof(niik_simple_ann)))==NULL),__func__,"calloc");
  NIIK_RET0(((fp=fopen(fname,"r"))==NULL),__func__,"fopen");
  fscanf(fp,"%i %i %i",&a->num,&a->ni,&a->nk);
  a->nh=5;
  a->in=niikmat_init(a->num,a->ni);
  a->out=niikmat_init(a->num,a->nk);
  a->alpha=niikmat_init(a->ni,a->nh);
  a->beta=niikmat_init(a->nh,a->nk);
  for(n=0; n<a->num; n++) {
    for(i=0; i<a->ni; i++) {
      fscanf(fp,"%lf ",&a->in->m[n][i]);
    }
    for(k=0; k<a->nk; k++) {
      fscanf(fp,"%lf ",&a->out->m[n][k]);
    }
  }
  fclose(fp);
  return a;
}

niik_simple_ann *niik_simple_ann_free(niik_simple_ann *a) {
  if(a==NULL) return NULL;
  free(a);
  return NULL;
}


int nmdim=0;
niik_simple_ann *nma=NULL;

double niik_simple_ann_train_dumb_nm_func(double *v) {
  int i,j,k,n;
  static double **b,**h;
  static double **alpha=NULL;
  static double **beta;
  double e=0;
  static int iter=0;
  static double emin=1e31;
  static double *eli=NULL;

  //niik_fc_display((char *)__func__,1);
  if(alpha==NULL) {
    b=(double **)calloc(nma->num,sizeof(double *));
    h=(double **)calloc(nma->num,sizeof(double *));
    for(n=0; n<nma->num; n++) {
      b[n]=(double *)calloc(nma->nk,sizeof(double));
      h[n]=(double *)calloc(nma->nh,sizeof(double));
    }
    eli=(double *)calloc(nma->num,sizeof(double));
    alpha=(double **)calloc(nma->ni,sizeof(double *));
    beta =(double **)calloc(nma->nh,sizeof(double *));
  }
  for(i=0; i<nma->ni; i++) {
    alpha[i]=v+i*nma->nh;
    //niik_display_double_vector(alpha[i],nma->nh);
  }
  for(i=0; i<nma->nh; i++) {
    beta[i]=v+nma->ni*nma->nh+i*nma->nk;
    //niik_display_double_vector(beta[i],nma->nk);
  }


  #pragma omp parallel for private(i,j,k) ordered schedule(static)
  for(n=0; n<nma->num; n++) {
    //    fprintf(stdout,"%i\n",n);
    eli[n]=0;
    for(k=0; k<nma->nk; k++) {
      b[n][k]=0;
      for(j=0; j<nma->nh; j++) {
        h[n][j]=0;
        for(i=0; i<nma->ni; i++) {
          h[n][j]+=alpha[i][j]*nma->in->m[n][i];
        }
        // b[n][k]+=beta[j][k]*h[n][j];
        b[n][k]+=beta[j][k]*NIIK_Heaviside(h[n][j],2);
      }
      //eli[n]+=(nma->out->m[n][k]-b[n][k])*(nma->out->m[n][k]-b[n][k]);
      //eli[n]+=NIIK_SQ(nma->out->m[n][k]-NIIK_Heaviside(b[n][k],2));
      eli[n]+=fabs(nma->out->m[n][k]-NIIK_Heaviside(b[n][k],2));
    }
  }

  for(n=0,e=0; n<nma->num; n++) {
    e+=eli[n];
  }

  /*  free(b);
      free(h);
      free(alpha);
      free(beta);
      free(eli);*/
  e/=nma->num;
  iter++;
  if(emin>e) {
    fprintf(stdout,"e = %9.6f   %9i\n",e,iter);
    emin=e;
    //niik_display_double_vector(v,nmdim);
    niik_write_double_vector("tmp_v.txt",v,nmdim);
  }
  return e;
}

int niik_simple_ann_train_dumb(niik_simple_ann *a) {
  int i,j,k,n;
  niikvec *h,*b;
  double e=0;
  niikmat *nmp=NULL;
  double nmtol;
  int nmcost=NIIK_NELDER_MEAD_COST_RATIO;
  double (*nmpfn)();
  int nmmax=5e6;

  niik_fc_display((char *)__func__,1);
  fprintf(stdout,"[%s:%i:%s] ni,nh,nk %i %i %i %i\n",__FILE__,__LINE__,__func__,a->num,a->ni,a->nh,a->nk);
  h=niikvec_init(a->nh);
  b=niikvec_init(a->nk);

  nmdim=a->ni*a->nh+a->nh*a->nk;
  nmp=niikmat_rand(nmdim+1,nmdim);
  niikmat_kadd(nmp,-0.5);
  niikmat_kmul(nmp,0.1);
  for(i=0; i<nmdim; i++)
    nmp->m[i][i]+=1;
  //nmp=niikmat_identity(nmdim+1,nmdim);
  //niikmat_display(nmp); exit(0);

  nmpfn=niik_simple_ann_train_dumb_nm_func;
  nma=a;

  //fprintf(stdout,"%i %f\n",n,niik_simple_ann_train_dumb_nm_func(nmp->m[n]));
  niik_nelder_mead(nmp,nmdim,&nmtol,nmcost,nmpfn,&nmmax);
  n=0;
  fprintf(stdout,"%i %i %f\n",n,nmmax,niik_simple_ann_train_dumb_nm_func(nmp->m[n]));
  niik_display_double_vector(nmp->m[n],nmdim);
  exit(0);

  for(n=0; n<a->num; n++) {
    for(k=0; k<a->nk; k++) {
      b->v[k]=0;
      for(j=0; j<a->nh; j++) {
        h->v[j]=0;
        for(i=0; i<a->ni; i++) {
          h->v[j]+=a->alpha->m[i][j]*a->in->m[n][i];
        }
        b->v[k]+=a->beta->m[j][k]*h->v[j];
      }
      e+=(a->out->m[n][k]-b->v[k])*(a->out->m[n][k]-b->v[k]);
      //fprintf(stdout,"e = %0.1f %0.2f %0.2f\n",e,a->out->m[n][k],b->v[k]);
    }
  }

  fprintf(stdout,"e = %0.1f\n",e);

  h=niikvec_free(h);
  b=niikvec_free(b);
  return 1;
}

#endif
