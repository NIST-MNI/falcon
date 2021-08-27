/* Filename:     nifti1_kunio_off_bbox.c
 * Description:  object functions for bounding box
 * Author:       Kunio Nakamura
 * Date:         March 1, 2012
 *
 *
 * geomview's off file
 * -bounding box system
 * -useful for fast triangle-triangle intersection test
 * -useful for fast point-triangle distance calculation
 *
 *
 * -bounding box uses octree system
 * -without octree, memory requirement would be too high
 *  or it would take a very long time
 *
 *
 *
 */
#include <math.h>

#include "falcon.h"
#include "falcon_surfaces.h"
#include "falcon_tritri_inline.h"

#ifndef   MAX
#define MAX(a,b) (a>b?a:b)
#endif

#ifndef   MIN
#define MIN(a,b) (a<b?a:b)
#endif

#ifndef   MAX3
#define MAX3(a,b,c) MAX(MAX(a,b),c)
#endif

#ifndef   MIN3
#define MIN3(a,b,c) MIN(MIN(a,b),c)
#endif

#ifndef   MAX4
#define MAX4(a,b,c,d) MAX(MAX(MAX(a,b),c),d)
#endif

#ifndef   MIN4
#define MIN4(a,b,c,d) MIN(MIN(MIN(a,b),c),d)
#endif


/*
 * initialize the bounding box system
 */

bbox *off_bbox_calloc() {
  bbox *B;
  B = (bbox *)calloc(1,sizeof(bbox));
  B->delta = B->bmax = 0;
  B->xdim=B->ydim=B->zdim=B->area=B->nvox=0;
  B->depth=0;
  B->data = NULL;
  B->ndata = NULL;
  B->origin.x=B->origin.y=B->origin.z=0.0;
  B->origin.w=1.0;
  return B;
}


bbox *off_bbox_free(bbox *bb) {
  int n;
  if(bb==NULL) return NULL;
  if(bb->data!=NULL) {
    for(n=0; n<bb->nvox; n++) {
      if(bb->data[n]!=NULL) free(bb->data[n]);
    }
    free(bb->data);
    free(bb->ndata);
  }
  free(bb);
  return NULL;
}


/*
 * off_bbox_init_depth
 *
 * -initializes the bouding box with depth and bmax
 * -depth 5 with bmax = 320 --> delta = 10mm
 * -depth 6 with bmax = 320 --> delta = 5mm
 * -depth 7 with bmax = 320 --> delta = 2.5mm
 */

bbox *off_bbox_init(int depth,double bmax) {
  bbox *bb;
  int n;
  bb = off_bbox_calloc();
  bb->depth=depth;
  n=1;

  bb->zdim=bb->ydim=bb->xdim=bb->xdim=(n<<depth);

  bb->area=bb->xdim*bb->ydim;
  bb->nvox=bb->area*bb->zdim;

  /*VF: hack: make it symmetric around the origin*/
  bb->origin.x=bb->origin.y=bb->origin.z= - bmax;
  bb->bmax  = bmax * 2;

  bb->delta = bb->bmax / bb->xdim;

  /*fprintf(stdout,"  off_bbox_init %3i %3i %3i  %8.4f\n",bb->xdim,bb->ydim,bb->zdim,bb->delta);*/

  bb->ndata = (int *)calloc(bb->nvox,sizeof(int));
  bb->data  = (kface ***)calloc(bb->nvox,sizeof(kface **));
  return bb;
}


/*
 * functions to update bbox system
 *
 * -sets up the bounding box system from list of triangles
 *
 *
 *   int off_create_bbox_from_kface_list(bbox *bb,kface **facelist,int num);
 *   int off_create_bbox_from_kobj(bbox *bb,kobj *obj);
 *
 *
 */

int off_create_bbox_from_kface_list(bbox *bb,
                                    kface **facelist,
                                    int num) {
  niikpt pmin,pmax;
  int n;
  double width;
  int vcount;

  if(bb==NULL)      {
    fprintf(stderr,"ERROR: bb is a null pointer\n");
    return 0;
  }
  if(facelist==NULL) {
    fprintf(stderr,"ERROR: facelist is a null pointer\n");
    return 0;
  }

  for( n=0; n<bb->nvox; n++ ) {
    if( bb->ndata[n]>0 ) {
      free(bb->data[n]);
      bb->data[n] = NULL;
      bb->ndata[n]=0;
    }
  }

  off_calc_facelist_bounds(facelist,num,&pmin,&pmax);

  width=MAX4(pmax.x-pmin.x, pmax.y-pmin.y, pmax.z-pmin.z, 1.0); /*VF: avoid unreasonably small width*/
  /*update bbox*/
  bb->origin=pmin;
  bb->bmax=width;
  
  bb->delta=width / (double)(bb->xdim);

  /*making cubic bounding box, TODO: allow for non-uniform delta*/
  pmax.x=pmin.x+width;
  pmax.y=pmin.y+width;
  pmax.z=pmin.z+width;

  off_create_bbox_octree(bb, facelist, num, bb->depth, pmin, pmax);

  return 1;
}


/* for multiple objects (like the pial/WM surfaces of cortex) */
int off_create_bbox_from_multiple_kobj(bbox *bb, kobj **obj, int num_obj) {
  const int verbose=0;
  int       n,nface;
  kface **facelist,*f;

  if(bb==NULL) {
    fprintf(stderr,"ERROR: bb is null\n");
    return 0;
  }
  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is null\n");
    return 0;
  }
  if(!num_obj) {
    fprintf(stderr,"ERROR: num_obj is zero\n");
    return 0;
  }
  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) start\n");
  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) pminmax\n");

  /*#pragma omp parallel for*/
  for(n=0; n<num_obj; n++) {
    off_update_kobj_kface_pminmax(obj[n]);
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) count faces\n");

  nface=0;
  for(n=0; n<num_obj; n++) {
    nface += obj[n]->nface;
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) calloc facelist %i\n",nface);

  facelist = (kface **)calloc(nface,sizeof(kface *));

  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) copy faces\n");

  for(n=nface=0; n<num_obj; n++) {
    for(f=obj[n]->face; f!=NULL; f=f->next) {
      facelist[nface++] = f;
    }
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) create from the list\n");

  if(!off_create_bbox_from_kface_list(bb,facelist,nface)) {
    fprintf(stderr,"ERROR: off_create_bbox_from_kface_list\n");
    return 0;
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_from_multiple_kobj) successful finish\n");

  free(facelist);

  return 1;
}


int off_create_bbox_from_kobj(bbox *bb,kobj *obj) {
  int n;
  kface **tmplist,*f;
  const int verbose=0;
  niikpt pmin,pmax;

  if(bb==NULL) return 0;
  if(obj==NULL) return 0;
  if(bb->nvox==0) return 0;

  off_kobj_update_num(obj);
  off_update_kobj_kface_pminmax(obj);

  if(verbose) fprintf(stdout,"  off_create_bbox_from_kobj\n");

  /* make enough temp space */
  if(verbose>1) fprintf(stdout,"    copy obj->face to facelist\n");
  if((tmplist = (kface **)calloc(obj->nface,sizeof(kface *)))==NULL) {
    fprintf(stderr,"ERROR: calloc tmplist\n");
    return 0;
  }

  for(f=obj->face,n=0; f!=NULL; f=f->next) {
    tmplist[n++]=f;
  }

  /* use off_create_bbox_from_kface_list */
  off_create_bbox_from_kface_list(bb,tmplist,obj->nface);

  free(tmplist);
  return 1;
}  /* int off_create_bbox_from_kobj(bbox *bb,kobj *obj) */


/*
 * main bounding box function
 *
 * -recursive function
 * -updates bbox so that bbox.data contains list of faces
 * -thus, data[index] is a list of triangles that falls between
 *  delta x and delta (x+1) where x = index % bbox.xdim
 *  delta y and delta (y+1) where y = floor((index % bbox.area) / bbox.xdim)
 *  delta z and delta (z+1) where z = floor(index/bbox.area)
 *
 */

int off_create_bbox_octree(bbox *bb, kface **facelist,
                           int num,
                           int cdepth,
                           niikpt pmin,
                           niikpt pmax) {
  kface ***newlist,*f;
  int
      index,
      n,*nlist;
  double
      xmid,ymid,zmid;
      niikpt fmid;
  const int verbose=0;

  if(verbose>1 && cdepth==0)
    fprintf(stdout,"  off_create_bbox_octree %i (%5.1f %5.1f %5.1f) (%5.1f %5.1f %5.1f)\n",cdepth,
            pmin.x, pmin.y, pmin.z, pmax.x, pmax.y, pmax.z );

  /* put triangles into bounding box system */
  if(cdepth==0) {
    int i,j,k;

    if(verbose>3) fprintf(stdout,"    put triangles into bbox\n");

    i=NIIK_IMINMAX((int)floor( (pmin.x - bb->origin.x) / bb->delta + 0.1 ), 0, bb->xdim-1);
    j=NIIK_IMINMAX((int)floor( (pmin.y - bb->origin.y) / bb->delta + 0.1 ), 0, bb->ydim-1);
    k=NIIK_IMINMAX((int)floor( (pmin.z - bb->origin.z) / bb->delta + 0.1 ), 0, bb->zdim-1);

    index = i+ j* bb->xdim + k* bb->area;

    if(verbose>3)
      fprintf(stdout,"      %3f %3f %3f  %3f %3f %3f  %3i %3i %3i| %8i | %8i| %4i %4i %4i\n",
              ((pmin.x - bb->origin.x) / bb->delta ),
              ((pmin.y - bb->origin.y) / bb->delta ),
              ((pmin.z - bb->origin.z) / bb->delta ),
              pmin.x,pmin.y,pmin.z,
              i,j,k,
              index,num,
              bb->xdim,bb->ydim,bb->zdim
             );

    if(bb->ndata[index] != 0) {
      fprintf(stderr,"Found cell already populated: %f,%f,%f\n",pmin.x,pmin.y,pmin.z);
      abort();
    }

    bb->ndata[index] = num;
    bb->data[index]=(kface **)calloc(num,sizeof(kface *));

    for(n=0; n<num; n++) {
      bb->data[index][n] = facelist[n];
    }
    return 1;
  }

  /* make eight lists */
  newlist = (kface ***)calloc(8,sizeof(kface **));
  nlist = (int *)calloc(8,sizeof(int));

  for(n=0; n<8; n++)  {
    if(verbose>2) fprintf(stdout,"    calloc %i\n",n);
    newlist[n] = (kface **)calloc(num,sizeof(kface *));
    nlist[n] = 0;
  }

  xmid = (pmin.x + pmax.x)/2.0;
  ymid = (pmin.y + pmax.y)/2.0;
  zmid = (pmin.z + pmax.z)/2.0;

  /* put each face into appropriate octree component */
  if(verbose>3) fprintf(stdout,"    go thru each triangle\n");

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 0: <xmid, <ymid, <zmid */
    if(f->pmin.x <= xmid &&
        f->pmin.y <= ymid &&
        f->pmin.z <= zmid )
      newlist[0][nlist[0]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];
    /* ROI 1: >xmid, <=ymid, <=zmid */
    if(f->pmax.x > xmid  &&
        f->pmin.y <= ymid &&
        f->pmin.z <= zmid )
      newlist[1][nlist[1]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 2: <=xmid, >ymid, <=zmid */
    if(f->pmin.x <= xmid &&
        f->pmax.y > ymid  &&
        f->pmin.z <= zmid)
      newlist[2][nlist[2]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 3: >xmid, >ymid, <=zmid */
    if(f->pmax.x > xmid &&
        f->pmax.y > ymid &&
        f->pmin.z <= zmid)
      newlist[3][nlist[3]++] = f;
  }
  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 4: <=xmid, <=ymid, >zmid */
    if(f->pmin.x <= xmid &&
        f->pmin.y <= ymid &&
        f->pmax.z > zmid)
      newlist[4][nlist[4]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 5: >xmid, <=ymid, >zmid */
    if(f->pmax.x > xmid &&
        f->pmin.y <= ymid &&
        f->pmax.z > zmid)
      newlist[5][nlist[5]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 6: Mxmid, >ymid, >zmid */
    if(f->pmin.x <= xmid &&
        f->pmax.y > ymid &&
        f->pmax.z > zmid)
      newlist[6][nlist[6]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];

    /* ROI 7: Mxmid, >ymid, >zmid */
    if(f->pmax.x > xmid &&
        f->pmax.y > ymid &&
        f->pmax.z > zmid)
      newlist[7][nlist[7]++] = f;
  }

  if(verbose>3) {
    niikpt p,q;
    for(n=0; n<8; n++) {
      if(n%2)     {
        p.x=xmid;
        q.x=pmax.x;
      } else        {
        p.x=pmin.x;
        q.x=xmid;
      }

      if((n%4)<2) {
        p.y=pmin.y;
        q.y=ymid;
      } else        {
        p.y=ymid;
        q.y=pmax.y;
      }

      if(n<4)     {
        p.z=pmin.z;
        q.z=zmid;
      } else        {
        p.z=zmid;
        q.z=pmax.z;
      }

      fprintf(stdout,"    %2i %9i (%5.1f %5.1f %5.1f) (%5.1f %5.1f %5.1f)\n",n,nlist[n],
              p.x,p.y,p.z,q.x,q.y,q.z);
    }
  }

  for(n=0; n<8; n++) {
    niikpt p,q;
    if(nlist[n]!=0) {

      if(n&1)     {
        p.x=xmid;
        q.x=pmax.x;
      } else        {
        p.x=pmin.x;
        q.x=xmid;
      }

      if(n&2)     {
        p.y=ymid;
        q.y=pmax.y;
      } else        {
        p.y=pmin.y;
        q.y=ymid;
      }

      if(n&4)     {
        p.z=zmid;
        q.z=pmax.z;
      } else        {
        p.z=pmin.z;
        q.z=zmid;
      }

      p.w=q.w=0;

      off_create_bbox_octree(bb, newlist[n], nlist[n], cdepth-1, p, q);
    }
  }

  for(n=0; n<8; n++)
    free(newlist[n]);

  free(newlist);
  free(nlist);

  return 1;
} /* int off_create_bbox_octree(bbox *bb,kface **facelist,int num,int cdepth,niikpt pmin,niikpt pmax) */




/*
 * -checks for intersection between triangle (f) and other
 *  triangles in the area (defined by bounding box, bb)
 *
 */

int off_check_tri_tri_intersect(kface * f1,kface * f2)
/* checks both ways */
{
  double
      v0[3],v1[3],v2[3],
      u0[3],u1[3],u2[3];

  v0[0]=f1->vert[0]->v.x;
  v0[1]=f1->vert[0]->v.y;
  v0[2]=f1->vert[0]->v.z;
  
  v1[0]=f1->vert[1]->v.x;
  v1[1]=f1->vert[1]->v.y;
  v1[2]=f1->vert[1]->v.z;

  v2[0]=f1->vert[2]->v.x;
  v2[1]=f1->vert[2]->v.y;
  v2[2]=f1->vert[2]->v.z;

  u0[0]=f2->vert[0]->v.x;
  u0[1]=f2->vert[0]->v.y;
  u0[2]=f2->vert[0]->v.z;

  u1[0]=f2->vert[1]->v.x;
  u1[1]=f2->vert[1]->v.y;
  u1[2]=f2->vert[1]->v.z;

  u2[0]=f2->vert[2]->v.x;
  u2[1]=f2->vert[2]->v.y;
  u2[2]=f2->vert[2]->v.z;

  return niik_tri_tri_intersect(v0,v1,v2, u0,u1,u2) ||
         niik_tri_tri_intersect(u0,u1,u2, v0,v1,v2);

} /* off_check_tri_tri_intersect */


int off_check_tri_tri_intersect_with_isectline(niikpt v1,niikpt v2,niikpt v3,
    niikpt u1,niikpt u2,niikpt u3,
    int *coplanar,
    niikpt *isect1, niikpt *isect2) {

  int outflag;
  double
    V0[3],V1[3],V2[3],
    U0[3],U1[3],U2[3],
    I1[3],I2[3];
    
  V0[0] = v1.x;
  V0[1] = v1.y;
  V0[2] = v1.z;

  V1[0] = v2.x;
  V1[1] = v2.y;
  V1[2] = v2.z;

  V2[0] = v3.x;
  V2[1] = v3.y;
  V2[2] = v3.z;

  U0[0] = u1.x;
  U0[1] = u1.y;
  U0[2] = u1.z;

  U1[0] = u2.x;
  U1[1] = u2.y;
  U1[2] = u2.z;

  U2[0] = u3.x;
  U2[1] = u3.y;
  U2[2] = u3.z;

  outflag=niik_tri_tri_intersect_with_isectline(V0,V1,V2,U0,U1,U2,coplanar,I1,I2);
  
  isect1->x = I1[0];
  isect1->y = I1[1];
  isect1->z = I1[2];
  isect2->x = I2[0];
  isect2->y = I2[1];
  isect2->z = I2[2];
  return outflag;
} /* off_check_tri_tri_intersect_with_isectline */


int off_check_self_intersection_kvert(bbox *bb,kvert *v) {
  int n,m=0;

  /*#pragma omp parallel for reduction(+:m) private(n)*/ /*TODO: potentially remove this one*/
  for(n=0; n<v->nei; n++) {
    if( m>0 ) break;/*continue  break - is breaking omp?*/

    if( off_check_self_intersection_kface(bb, v->neiface[n]) )
      m++;
  }
  return m;
} /* off_check_self_intersection_kvert */



static inline int off_check_faces_intersect(kface* f1, kface *f2)
{
    if(
       f1->pmax.x < f2->pmin.x || /*spatially separate*/
       f1->pmax.y < f2->pmin.y ||
       f1->pmax.z < f2->pmin.z ||
       f2->pmax.x < f1->pmin.x ||
       f2->pmax.y < f1->pmin.y ||
       f2->pmax.z < f1->pmin.z ||

       f1->vert[0]==f2->vert[0]|| /*neigbours*/
       f1->vert[0]==f2->vert[1]||
       f1->vert[0]==f2->vert[2]||
       f1->vert[1]==f2->vert[0]||
       f1->vert[1]==f2->vert[1]||
       f1->vert[1]==f2->vert[2]||
       f1->vert[2]==f2->vert[0]||
       f1->vert[2]==f2->vert[1]||
       f1->vert[2]==f2->vert[2] )
       return 0;

    return off_check_tri_tri_intersect(f1,f2);
}


int off_check_self_intersection_kface(bbox *bb, kface *f)
/* returns nonzero if there IS intersection
 * returns zero if there is NO intersection */
{

  int i,j,k,index,n;
  int minidx[4],maxidx[4];
  const int verbose=0;
  int testa=0;
  int intersect=0;

  minidx[1]=NIIK_IMINMAX(floor( (f->pmin.x - bb->origin.x)/bb->delta - 1),0,bb->xdim-1);
  minidx[2]=NIIK_IMINMAX(floor( (f->pmin.y - bb->origin.y)/bb->delta - 1),0,bb->ydim-1);
  minidx[3]=NIIK_IMINMAX(floor( (f->pmin.z - bb->origin.z)/bb->delta - 1),0,bb->zdim-1);

  maxidx[1]=NIIK_IMINMAX(ceil( (f->pmax.x - bb->origin.x)/bb->delta + 1),0,bb->xdim-1);
  maxidx[2]=NIIK_IMINMAX(ceil( (f->pmax.y - bb->origin.y)/bb->delta + 1),0,bb->ydim-1);
  maxidx[3]=NIIK_IMINMAX(ceil( (f->pmax.z - bb->origin.z)/bb->delta + 1),0,bb->zdim-1);

  if(verbose>0) {
    for(k=minidx[3],n=0; k<=maxidx[3]; k++) {
      for(j=minidx[2]; j<=maxidx[2]; j++) {
        index = minidx[1] + j*bb->xdim + k*bb->area;
        for(i=minidx[1]; i<=maxidx[1]; index++,i++) {
          n+=bb->ndata[index];
        }
      }
    }
    fprintf(stdout,"[off_check_self_intersection_kface]    num = %i\n",n);
  }

  /*#pragma omp parallel for private(i,j,n,index) reduction(+:intersect) */ /*reduction(+:testa) */
#if _OPENMP>=201307
  #pragma omp simd
#endif
  for(k=minidx[3]; k<=maxidx[3]; k++) {
    for(j=minidx[2]; j<=maxidx[2]; j++) {
      if( intersect>0 ) continue;
      index = minidx[1] + j*bb->xdim + k*bb->area;
      for(i=minidx[1]; i<=maxidx[1]; index++,i++) {

        for(n=0; n<bb->ndata[index]; n++) {
          kface *cf;
          cf=bb->data[index][n];

          if( off_check_faces_intersect(f,cf)) {
            intersect++;
            /*VF: TODO: break here?*/
          }
        }
      }
    }
  }
  return intersect; 
} /* off_check_self_intersection_kface */



/*
 * function to check self-intersection
 *
 *
 * off_count_self_intersection
 * -basic function for counting # of self-intersecting triangles
 *
 *
 * off_count_self_intersection_add_color
 * -similar to off_count_self_intersection
 * -the self-intersecting triangle becomes red (others are gray 0.66)
 * -color coding occurs when flag_add is nonzero
 *
 */


int off_count_self_intersection(bbox *bb,kobj *obj) {
  return off_count_self_intersection_add_color(bb,obj,0);
}


int off_check_self_intersection_from_bbox(bbox *bb)
/* makes the face red for self-intersections */
{
  int   n;
  if(bb==NULL) {
    fprintf(stderr,"ERROR: bbox is null\n");
    return -1;
  }

  /*#pragma omp parallel for private(i,j)*/
  for(n=0; n<bb->nvox; n++) {
    int j;
    for(j=0; j<bb->ndata[n]; j++) {
      if(bb->data[n][j]->color==NULL) {
        bb->data[n][j]->color=(double *)calloc(4,sizeof(double));
      } else {
        /*set to black */
        /*TODO: maybe use something else?*/
        bb->data[n][j]->color[0]=
        bb->data[n][j]->color[1]=
        bb->data[n][j]->color[2]=
        bb->data[n][j]->color[3]=0;
      }
      off_update_kface_pminmax(bb->data[n][j]);
    }
  }

  /*#pragma omp parallel for private(n)*/
  for(n=0; n<bb->nvox; n++) {
    int i;
    for(i=0; i<bb->ndata[n]; i++) {
      int j;
      kface *f1;
      f1=bb->data[n][i];

      for(j=i+1; j<bb->ndata[n]; j++) {
        kface *f2;
        f2=bb->data[n][j];

        if( off_check_faces_intersect(f1,f2)) {
          /*marking both faces red*/
          f1->color[0] = 1.0;
          f1->color[1] = f1->color[2] = f1->color[3] = 0.0;

          f2->color[0] = 1.0;
          f2->color[1] = f2->color[2]  = f2->color[3] = 0.0;
        }
      }
    }
  }

  return 1;
} /* int off_check_self_intersection_from_bbox(bbox *bb) */


int off_count_self_intersection_add_color(bbox *bb, kobj *obj, int flag_add)
/* -bb can be null if not created before this function
 * -obj is the object
 * -flag_add if non-zero, puts red surface for self intersection
 * -returns the number of surface intersection triangles */
{
  bbox *thisbb;
  unsigned char *sivec; /* self-intersection results here */
  int  n;
  int  cnt=0,
       nobb=0;

  if(obj==NULL) {
    fprintf(stderr,"ERROR: obj is a null pointer\n");
    return -1;
  }
  if(bb==NULL) {
    nobb=1;
    thisbb=off_bbox_init(7,320);
  } else {
    nobb=0;
    thisbb=bb;
  }

  if(!off_create_bbox_from_kobj(thisbb,obj)) {
    fprintf(stderr,"ERROR: off_create_bbox_from_kobj\n");
    return -1;
  }

  if(flag_add) {
    off_kobj_gray_kface(obj);
  }

  off_kobj_update_face_index(obj);

  sivec = (unsigned char *)calloc(obj->nface+1, sizeof(unsigned char));

  /*#pragma omp parallel for private(n,i,j)*/
  for(n=0; n<thisbb->nvox; n++) {
    int i;
    for(i=0; i<thisbb->ndata[n]; i++) {
      int j;
      kface *f1;

      f1 = thisbb->data[n][i];
      for(j=i+1; j<thisbb->ndata[n]; j++) {
        kface *f2;
        f2 = thisbb->data[n][j];

        if(sivec[f1->index]>0 && sivec[f2->index]>0) continue;

        if(off_check_faces_intersect(f1,f2)) {
          sivec[f1->index] = 1;
          sivec[f2->index] = 1;
          if(flag_add) {
            f1->color[0]=f2->color[0]=1.0;
            f1->color[1]=f2->color[1]=f1->color[2]=f2->color[2]=f1->color[3]=f2->color[3]=0.0;
          }
        }
      }
    }
  }

#if _OPENMP>=201307
//  #pragma omp simd
#endif
  for(n=1,cnt=0; n<=obj->nface; n++) {
    cnt+=sivec[n];
  }

  free(sivec);
  if( nobb ) {
     thisbb=off_bbox_free(thisbb);
  }
  return cnt;
} /* off_count_self_intersection_add_color */



int off_correct_self_intersection(bbox *bb,kobj *obj)
/* -corrects (removes) surface intersections
 * -smoothes (relocates vertices)
 * -remeshes if necessary
 */
{
  int n;
  if((n=off_correct_self_intersection_with_elen(bb,obj,-1))<0) {
    fprintf(stderr,"ERROR: off_correct_self_intersection_with_elen(bb,obj,-1)\n");
    return -1;
  }
  return n;
}


int off_correct_self_intersection_with_elen(bbox *bb,kobj *obj,double elen)
/* -corrects (removes) surface intersections
 * -smoothes (relocates vertices)
 * -remeshes if necessary
 * -elen is for remeshing
 * --if negative, then use the mean from the object
 */
{
  kface **flist,*f;
  kvert *v;
  int iter = 0,
      cnt,n,m;
  int verbose = niik_verbose();

  cnt = off_count_self_intersection_add_color(bb,obj,1);

  while(cnt>0) {
    if(verbose>1) fprintf(stdout,"[off_correct_self_intersection] iter %3i   intersections %-i\n",iter,cnt);
    iter++;
    flist = off_face_list(cnt);

    /* find the surface with self-intersections */
    for(f=obj->face,n=0; f!=NULL; f=f->next) {
      if(f->color[0]>f->color[1]) {
        /* fprintf(stdout,"  face %i   %i\n",f->index,n); */
        flist[n++]=f;
      }
    }
    /* smoothes with relocate function */
    for(n=0; n<cnt; n++) {
      for(m=0; m<3; m++) {
        off_remesh_kvert_relocate(flist[n]->vert[m]);
      }
    }

    /* do remeshing (iter>10) */
    if((iter%2)==0 || iter>10) {
      for(v=obj->vert; v!=NULL; v=v->next) {
        v->v.w=1;
      }
      for(n=0; n<cnt; n++) {
        for(m=0; m<3; m++) {
          flist[n]->vert[m]->v.w=0;
        }
      }
      for(m=15; m<iter; m++) {
        for(v=obj->vert; v!=NULL; v=v->next) {
          if(v->v.w>0) continue;
          for(n=0; n<v->nei; n++) {
            v->neivert[n]->v.w=0;
          }
        }
      }
      if(elen<0) {
        elen = off_get_kobj_mean_edge_length(obj);
      }
      elen+=0.2*elen*(niik_get_rand()-0.5);
      if(!off_remesh_kobj(obj,elen,5,1)) {
        fprintf(stderr,"ERROR: off_remesh_kobj\n");
        exit(0);
      }
    }
    free(flist);

    cnt = off_count_self_intersection_add_color(bb,obj,1);
    /* fprintf(stdout,"      off_count_self_intersection_add_color %i\n",cnt); */
    if(iter>50) {
      fprintf(stderr,"ERROR: can't correct self-intersection\n");
      for(f=obj->face,n=0; f!=NULL; f=f->next) {
        if(f->color[0]>f->color[1]) {
          fprintf(stderr,"\n");
          fprintf(stderr,"%3i face[%i] %7.3f %7.3f %7.3f edge\n",n,f->index,
                  niikpt_distance(f->edge[0]->endpts[0]->v,f->edge[0]->endpts[1]->v),
                  niikpt_distance(f->edge[1]->endpts[0]->v,f->edge[1]->endpts[1]->v),
                  niikpt_distance(f->edge[2]->endpts[0]->v,f->edge[2]->endpts[1]->v));
          fprintf(stderr,"%3i face[%i] %7.3f %7.3f %7.3f\n",n,f->index,f->vert[0]->v.x,f->vert[0]->v.y,f->vert[0]->v.z);
          fprintf(stderr,"%3i face[%i] %7.3f %7.3f %7.3f\n",n,f->index,f->vert[1]->v.x,f->vert[1]->v.y,f->vert[1]->v.z);
          fprintf(stderr,"%3i face[%i] %7.3f %7.3f %7.3f\n",n,f->index,f->vert[2]->v.x,f->vert[2]->v.y,f->vert[2]->v.z);
        }
      }
      for(n=0; n<30; n++) {
        for(v=obj->vert; v!=NULL; v=v->next) {
          off_remesh_kvert_relocate(v);
        }
      }
      fprintf(stderr,"  writing tmp_prob.off for error checking \n");
      if(!off_kobj_write_offply("tmp_prob.off",obj,0)) {
        fprintf(stderr,"ERROR: off_kobj_write_off \n");
        exit(0);
      }
      return 0;
    }
  } /* while (cnt>0) */
  return 1;
}


/* assumes that self-intersection test has been done previously with _add_color */
void off_display_self_intersection_kface(kobj *obj) {
  kface *f;
  for(f=obj->face; f!=NULL; f=f->next) {
    if(f->color[0]>f->color[1]) {
      off_display_face_info(f,20);
    }
  }
  return;
}


/* assuming that we already have everything in bb
 * checks every data in bb */
int off_check_bbox_self_intersection(bbox *bb) {

  int n;
  int verbose=niik_verbose();

  for(n=0; n<bb->nvox; n++) {
    int i;
    for(i=0; i<bb->ndata[n]; i++) {
      int j;
      kface *f1 = bb->data[n][i];
      for(j=i+1; j<bb->ndata[n]; j++) {
        kface *f2 = bb->data[n][j];

        if(off_check_faces_intersect(f1,f2)) {
          if(verbose>1) {
            fprintf(stdout,"  selfintersection found in %i [%3i %3i %3i] (delta=%5.1f) between f(%i) and f(%i)\n",
                    n,n%bb->xdim,(int)floor((n%bb->area)/bb->xdim),n/bb->area,bb->delta,f1->index,f2->index);
            off_display_face_info(f1,20);
            off_display_face_info(f2,20);
          }
          return 1;
        }
      }
    }
  }
  return 0;
}


void off_display_bbox_info(bbox *bb) {
  int n,num;
  fprintf(stdout,"bbox info\n");
  fprintf(stdout,"  max = %0.4f   delta = %0.6f\n",bb->bmax,bb->delta);
  fprintf(stdout,"  dim = %i %i %i  |  xy=%i | xyz=%i\n",bb->xdim,bb->ydim,bb->zdim,bb->area,bb->nvox);
  for(n=num=0; n<bb->nvox; n++) {
    if(bb->ndata[n]>0) num++;
  }
  fprintf(stdout,"  ndata  %i has at least one triangle\n",num);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/