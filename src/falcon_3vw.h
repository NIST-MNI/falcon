#ifndef _FALCON_3VW_H_
#define _FALCON_3VW_H_

#include "falcon.h"

typedef struct {
  nifti_image *img;
  niikmat *stx;
  niikpt ctr;    // point of interest in stx space
  double width;  // output
  int nmethod;   // 0 = single measurement, 1 = weighted average 5x5 neighborhood
  double xran;   // sample range
  double dx;     // sampling step size
  int verbose;
} niik_image_3vw; /* class for measuring 3ventricular width */

niik_image_3vw *niik_image_3vw_init() {
  niik_image_3vw *v=NULL;
  NIIK_RET0(((v=(niik_image_3vw *)calloc(1,sizeof(niik_image_3vw)))==NULL),__func__,"calloc in niik_image_3vw");
  v->img=NULL;
  v->stx=NULL;
  v->ctr=niikpt_val(95,116,75,1);
  v->width=-1;
  v->xran=20;
  v->dx=1;
  v->nmethod=0;
  v->verbose=1;
  return v;
}
niik_image_3vw *niik_image_3vw_free(niik_image_3vw *v) {
  NIIK_RET0n((v==NULL),__func__,"v is null");
  if(v->img!=NULL) v->img=niik_image_free(v->img);
  if(v->stx!=NULL) v->stx=niikmat_free(v->stx);
  free(v);
  return NULL;
}

int niik_image_3vw_update_measure_width(niik_image_3vw *v,niikpt p,niikpt normal,double *width) {
  niikvec
  *yy=NULL,
   *y=NULL;
  double
  yp,yn,
  xmin,dx,xp,xn;
  int n,num,nmid;

  if(v->verbose>1)fprintf(stdout,"[%s:%i:%s] %8.3f %8.3f %8.3f\n",__FILE__,__LINE__,__func__,p.x,p.y,p.z);
  xmin=-fabs(v->xran);
  dx=v->dx;
  num=v->xran*2/dx+1;
  if(v->verbose>1) fprintf(stdout,"[%s:%i:%s] %8.3f %8.3f %8.3f  %5i\n",__FILE__,__LINE__,__func__,xmin,dx,v->xran,num);
  y=niikvec_init(num);
  yy=niikvec_init(num);
  // interpolation
  for(n=0; n<num; n++) {
    yy->v[n]=y->v[n]=niik_image_interpolate_3d(v->img,niikpt_move_normal(p,normal,n*dx+xmin),NIIK_INTERP_LINEAR);
  }

  // zero-crossing with smoothing then without smoothing
  if(v->verbose>1) niikvec_display(y);
  niik_runavg_double_vector(y->v,y->num,num/10);
  if(v->verbose>1) niikvec_display(y);
  nmid=(y->num-1)/2;
  niik_central_difference_double_vector(y->v,y->num);
  if(v->verbose>1)niikvec_display(y);
  for(n=0; n<num/5; n++) {
    if     (y->v[nmid-1]<0 && y->v[nmid+1]>0) break;
    else if(y->v[nmid-1]>0 && y->v[nmid+1]<0) break;
    else if(y->v[nmid]<0) nmid++;
    else if(y->v[nmid]>0) nmid--;
  }
  for(n=0; n<num; n++) {
    y->v[n]=yy->v[n];
  }
  niik_central_difference_double_vector(y->v,y->num);
  if(v->verbose>1)niikvec_display(y);
  for(n=0; n<num/15; n++) {
    if     (y->v[nmid-1]<0 && y->v[nmid+1]>0) break;
    else if(y->v[nmid-1]>0 && y->v[nmid+1]<0) break;
    else if(y->v[nmid]<0) nmid++;
    else if(y->v[nmid]>0) nmid--;
  }
  for(n=0; n<num; n++) {
    y->v[n]=fabs(y->v[n]);
  }
  if(v->verbose>1) fprintf(stdout,"  min = %f at %.1f\n",yy->v[nmid],nmid*dx+xmin);

  // mode right and left
  niik_get_mode_bspline_vector(y->v+nmid,num-nmid,&xp,&yp);
  xp=(xp+nmid)*dx+xmin;
  if(v->verbose>1)fprintf(stdout,"[%s] peak+ %9.4f at %9.4f\n",__func__,yp,xp);
  niik_get_mode_bspline_vector(y->v,nmid+1,&xn,&yn);
  xn=xn*dx+xmin;
  if(v->verbose>1)fprintf(stdout,"[%s] peak- %9.4f at %9.4f\n",__func__,yn,xn);

  // width (difference)
  *width=(xp-xn);
  if(v->verbose>1)fprintf(stdout,"[%s] width %9.4f\n",__func__,*width);
  y=niikvec_free(y);
  yy=niikvec_free(yy);
  return 1;
}

int niik_image_3vw_update(niik_image_3vw *v)
/* main function */
{
  char fcname[32]="niik_image_3vw_update";
  niikpt p,pp,normal,ps[25];
  double w[25],dm;
  niikmat *inv;
  int i,j,n;

  if(v->verbose>0)niik_fc_display(fcname,1);
  NIIK_RET0((v==NULL),fcname,"v is null");
  NIIK_RET0((v->img==NULL),fcname,"v->img is null");
  NIIK_RET0((v->width>0),fcname,"width is already calculated");

  for(i=n=0; i<5; i++) {
    for(j=0; j<5; j++,n++) {
      ps[n]=niikpt_zero();
      ps[n].y=i-2.0;
      ps[n].z=j-2.0;
    }
  }

  inv=niikmat_inverse(v->stx);
  for(n=0; n<25; n++)
    ps[n]=niikpt_affine_transform(inv,niikpt_add(ps[n],v->ctr));
  p=niikpt_affine_transform(inv,v->ctr);
  pp=niikpt_affine_transform(inv,niikpt_add(v->ctr,niikpt_val(1,0,0,0)));
  if(v->verbose>0) fprintf(stdout,"[%s] center %8.3f %8.3f %8.3f\n",fcname,p.x/v->img->dx,p.y/v->img->dy,p.z/v->img->dz);
  normal=niikpt_unit(niikpt_sub(pp,p));
  if(v->verbose>0) fprintf(stdout,"[%s] normal %8.3f %8.3f %8.3f\n",fcname,normal.x,normal.y,normal.z);
  inv=niikmat_free(inv);

  switch(v->nmethod) {
  case 0:
    niik_image_3vw_update_measure_width(v,p,normal,&v->width);
    break;
  case 1:
    niik_image_3vw_update_measure_width(v,p,normal,&dm);
    for(n=0,dm*=5; n<25; n++) {
      niik_image_3vw_update_measure_width(v,ps[n],normal,w+n);
      dm+=w[n];
    }
    v->width=dm/30.0;
    break;
  }

  if(v->verbose>0)niik_fc_display(fcname,0);
  return 1;
}


#endif

/*
 kate: space-indent on; hl c;indent-width 4; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/