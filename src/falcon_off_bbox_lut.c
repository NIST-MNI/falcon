/* Filename:     nifti1_kunio_off_bbox_lut.c
 * Description:  object functions lut for faces around vertices
 * Author:       Kunio Nakamura
 * Date:         May 7, 2014
 *
 */

#include "falcon.h"

/*
 * initialize the bounding box system
 */

bboxlut *off_bboxlut_calloc() {
  bboxlut *B;
  B = (bboxlut *)calloc(1,sizeof(bbox));
  B->ndata=0;
  B->dist=0;
  B->data=NULL;
  B->idata=NULL;
  return B;
}

bboxlut *off_bboxlut_free(bboxlut *bb) {
  int n;
  if(bb==NULL) return NULL;
  if(bb->idata!=NULL) free(bb->idata);
  if(bb->data!=NULL) {
    for(n=0; n<bb->ndata; n++) {
      if(bb->data[n]!=NULL) free(bb->data[n]);
    }
    free(bb->data);
  }
  free(bb);
  return NULL;
}


int off_bboxlut_update_objs(bboxlut *bb,kobj **objs,int nobjs) {
  const char *fcname="off_bboxlut_update_objs";
  int m,n,max_nf,nn,i,j,k,minidx[4],maxidx[4],index;
  kvert *v;
  kface **flist,*f;
  niikpt pmin,pmax;
  bbox *b;
  niik_fc_display(fcname,1);
  b=off_bbox_init(7,320);
  NIIK_RET0((!off_create_bbox_from_multiple_kobj(b,objs,nobjs)),fcname,"off_create_bbox_from_multiple_kobj");
  NIIK_RET0((bb==NULL),fcname,"bb is null");
  NIIK_RET0((objs==NULL),fcname,"objs is null");
  NIIK_RET0((nobjs<=0),fcname,"#objs is not accepted");
  for(n=0; n<bb->ndata; n++) {
    if(bb->data[n]!=NULL) {
      free(bb->data[n]);
      bb->data[n]=NULL;
    }
  }
  for(n=max_nf=bb->ndata=0; n<nobjs; n++) {
    off_update_kobj_kface_pminmax(objs[n]);
    max_nf=(max_nf<objs[n]->nvert)?objs[n]->nvert:max_nf;
    bb->ndata+=objs[n]->nvert;
  }
  flist=(kface **)calloc(max_nf,sizeof(kface *));
  bb->data=(kface ***)calloc(bb->ndata,sizeof(kface **));
  bb->idata=(int *)calloc(bb->ndata,sizeof(int));
  for(n=m=0; n<nobjs; n++) {
    for(v=objs[n]->vert; v!=NULL; v=v->next,m++) {
      bb->idata[m]=0;
      pmin.x=v->v.x-bb->dist;
      pmin.y=v->v.y-bb->dist;
      pmin.z=v->v.z-bb->dist;
      pmax.x=v->v.x+bb->dist;
      pmax.y=v->v.y+bb->dist;
      pmax.z=v->v.z+bb->dist;

      minidx[1]=NIIK_IMINMAX(floor(pmin.x/b->delta-1.5),0,b->xdim-1);
      minidx[2]=NIIK_IMINMAX(floor(pmin.y/b->delta-1.5),0,b->ydim-1);
      minidx[3]=NIIK_IMINMAX(floor(pmin.z/b->delta-1.5),0,b->zdim-1);
      maxidx[1]=NIIK_IMINMAX(floor(pmax.x/b->delta+1.5),0,b->xdim-1);
      maxidx[2]=NIIK_IMINMAX(floor(pmax.y/b->delta+1.5),0,b->ydim-1);
      maxidx[3]=NIIK_IMINMAX(floor(pmax.z/b->delta+1.5),0,b->zdim-1);

      for(k=minidx[3]; k<=maxidx[3]; k++) {
        for(j=minidx[2]; j<=maxidx[2]; j++) {
          index = minidx[1] + j*b->xdim + k*b->area;
          for(i=minidx[1]; i<=maxidx[1]; index++,i++) {
            for(nn=0; nn<b->ndata[index]; nn++) {
              f=b->data[index][nn];
              if(f->pmax.x < pmin.x) continue;
              if(f->pmax.y < pmin.y) continue;
              if(f->pmax.z < pmin.z) continue;
              if(f->pmin.x > pmax.x) continue;
              if(f->pmin.y > pmax.y) continue;
              if(f->pmin.z > pmax.z) continue;
              if(f->vert[0]==v) continue;
              if(f->vert[1]==v) continue;
              if(f->vert[2]==v) continue;
              if(niikpt_distance(niikpt_closest_point_on_triangle_to_point(v->v,f->vert[0]->v,f->vert[1]->v,f->vert[2]->v),v->v)<bb->dist) {
                flist[bb->idata[m]++]=f;
              } /* added to temp list */
            }
          }
        }
      }
      bb->data[m]=(kface **)calloc(bb->idata[m],sizeof(kface *));
      for(nn=0; nn<bb->idata[m]; nn++) {
        bb->data[m][nn]=flist[nn];
      }
      //fprintf(stdout,"[%s] %9i   %-i \n",fcname,m,bb->idata[m]);
    } /* each vertex */
  } /* nobjs */
  off_bbox_free(b);
  free(flist);
  niik_fc_display(fcname,0);
  return 1;
}

