/* Filename:     nifti1_kunio_off_bbox_sph.c
 * Description:  object functions for spherical bounding box
 * Author:       Kunio Nakamura , Vladimir FONOV
 * Date:         Aug 22, 2016
 *
 *
 *
 * -bounding box uses tree system
 *
 */
#include "falcon.h"
#include "falcon_surfaces.h"


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


int off_create_bbox_sph_octree(bbox_sph *bb, kface **facelist,
                               int num,
                               int cdepth,
                               niiksph pmin,
                               niiksph pmax);


static inline void face_sph_bounds(kface *face, niiksph *pmin, niiksph *pmax) {
  *pmin=niiksph_min(face->vert[0]->sph, niiksph_min(face->vert[1]->sph, face->vert[2]->sph));
  *pmax=niiksph_max(face->vert[0]->sph, niiksph_max(face->vert[1]->sph, face->vert[2]->sph));
  if(pmin->the<0) pmin->the+=NIIK_PI2;
  if(pmin->psi<0) pmin->psi+=NIIK_PI2;

}


/*
 * initialize the bounding box system
 */

bbox_sph *off_bbox_sph_calloc() {
  bbox_sph *B;
  B          = (bbox_sph *)calloc(1,sizeof(bbox_sph));
  B->delta   = B->bmax = 0;
  B->psi_dim = B->the_dim=B->area=B->nvox=0;
  B->depth   = 0;
  B->data    = NULL;
  B->ndata   = NULL;
  B->origin  = niiksph_zero();
  return B;
}


bbox_sph *off_bbox_sph_free(bbox_sph *bb) {
  int n;
  if(bb == NULL) return NULL;
  if(bb->data != NULL) {
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
 * off_bbox_sph_init_depth
 *
 * -initializes the bouding box with depth and bmax
 * -depth 5 with bmax = 320 --> delta = 10mm
 * -depth 6 with bmax = 320 --> delta = 5mm
 * -depth 7 with bmax = 320 --> delta = 2.5mm
 */

bbox_sph *off_bbox_sph_init(int depth,double bmax) {
  bbox_sph *bb;
  int         n=1;
  bb          = off_bbox_sph_calloc();
  bb->depth   = depth;

  bb->bmax    = bmax;

  bb->origin  = niiksph_zero();

  bb->the_dim = bb->psi_dim = bb->psi_dim = (n<<depth);

  bb->nvox    = bb->area=bb->psi_dim*bb->the_dim;

  bb->delta   = bb->bmax / bb->psi_dim;

  bb->ndata   = (int *)calloc(bb->nvox,sizeof(int));

  bb->data    = (kface ***)calloc(bb->nvox,sizeof(kface **));
  return bb;
}


/*
 * functions to update bbox_sph system
 *
 * -sets up the bounding box system from list of triangles
 *
 *
 *   int off_create_bbox_sph_from_kface_list(bbox_sph *bb,kface **facelist,int num);
 *   int off_create_bbox_sph_from_kobj(bbox_sph *bb,kobj *obj);
 *
 *
 */
int off_create_bbox_sph_from_kface_list(bbox_sph *bb,
                                        kface **facelist,
                                        int num) {
  niiksph pmin,pmax;
  int n;
  double width;
  int vcount;

  if(bb==NULL)       {
    fprintf(stderr,"ERROR: bb is a null pointer\n");
    return 0;
  }
  if(facelist==NULL) {
    fprintf(stderr,"ERROR: facelist is a null pointer\n");
    return 0;
  }

  for(n=0; n<bb->nvox; n++) {
    if(bb->ndata[n]>0) {
      free(bb->data[n]);
      bb->data[n] = NULL;
      bb->ndata[n]=0;
    }
  }
  /*MAYBE always use 0 - PI*2 ? */
  /*off_calc_facelist_bounds_sph(facelist,num,&pmin,&pmax);*/

  /*width=MAX3(pmax.psi-pmin.psi, pmax.the-pmin.the, 1.0); /*VF: avoid unreasonably small width*/
  width=NIIK_PI2;
  /*update bbox_sph*/
  bb->origin=niiksph_zero();
  pmin=niiksph_zero();
  /*making square bounding box*/
  pmax.the=pmin.the+width;
  pmax.psi=pmin.psi+width;

  bb->bmax =width;
  bb->delta=width / (double)(bb->psi_dim);

  off_create_bbox_sph_octree(bb, facelist, num, bb->depth, pmin, pmax);

  return 1;
}


/* for multiple objects (like the pial/WM surfaces of cortex) */
int off_create_bbox_sph_from_multiple_kobj(bbox_sph *bb,kobj **obj,int num_obj) {
  const int  verbose=0;
  int   n,nface;
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
  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) start\n");
  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) pminmax\n");

  /*for(n=0;n<num_obj;n++) { off_update_kobj_kface_pminmax(obj[n]); }*/

  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) count faces\n");

  for(n=nface=0; n<num_obj; n++) {
    nface += obj[n]->nface;
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) calloc facelist %i\n",nface);

  facelist = (kface **)calloc(nface,sizeof(kface *));

  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) copy faces\n");

  for(n=0,nface=0; n<num_obj; n++) {
    for(f=obj[n]->face; f!=NULL; f=f->next) {
      facelist[nface++] = f;
    }
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) create from the list\n");

  if(!off_create_bbox_sph_from_kface_list(bb,facelist,nface)) {
    fprintf(stderr,"ERROR: off_create_bbox_sph_from_kface_list\n");
    return 0;
  }

  if(verbose) fprintf(stdout,"-v (off_create_bbox_sph_from_multiple_kobj) successful finish\n");

  free(facelist);

  return 1;
}


int off_create_bbox_sph_from_kobj(bbox_sph *bb,kobj *obj) {
  int n;
  kface **tmplist,*f;
  const int verbose=0;
  niiksph pmin,pmax;

  if(bb==NULL) return 0;
  if(obj==NULL) return 0;
  if(bb->nvox==0) return 0;

  off_kobj_update_num(obj);
  /*off_update_kobj_kface_pminmax(obj);*/

  if(verbose) fprintf(stdout,"  off_create_bbox_sph_from_kobj\n");

  /* make enough temp space */
  if(verbose>1) fprintf(stdout,"    copy obj->face to facelist\n");
  if((tmplist = (kface **)calloc(obj->nface,sizeof(kface *)))==NULL) {
    fprintf(stderr,"ERROR: calloc tmplist\n");
    return 0;
  }

  for(f=obj->face,n=0; f!=NULL; f=f->next) {
    tmplist[n++]=f;
  }

  /* use off_create_bbox_sph_from_kface_list */
  off_create_bbox_sph_from_kface_list(bb,tmplist,obj->nface);

  free(tmplist);
  return 1;
}  /* int off_create_bbox_sph_from_kobj(bbox_sph *bb,kobj *obj) */


/*
 * main bounding box function
 *
 * -recursive function
 * -updates bbox_sph so that bbox_sph.data contains list of faces
 * -thus, data[index] is a list of triangles that falls between
 *  delta x and delta (x+1) where x = index % bbox_sph.psi_dim
 *  delta y and delta (y+1) where y = floor((index % bbox_sph.area) / bbox_sph.psi_dim)
 *  delta z and delta (z+1) where z = floor(index/bbox_sph.area)
 *
 */

int off_create_bbox_sph_octree(bbox_sph *bb, kface **facelist,
                               int num,
                               int cdepth,
                               niiksph pmin,
                               niiksph pmax) {
  kface ***newlist,*f;
  int
  index,
  n,*nlist;
  double
  psi_mid,the_mid;
  niiksph fmid;
  const int verbose=0;
  int k=0;

  if(verbose>1 && cdepth==0)
    fprintf(stdout,"  off_create_bbox_sph_octree %i (%5.1f %5.1f ) (%5.1f %5.1f )\n",cdepth,
            pmin.psi, pmin.the, pmax.psi, pmax.the );

  /* put triangles into bounding box system */
  if(cdepth==0) {
    int i,j,k;

    if(verbose>3) fprintf(stdout,"    put triangles into bbox_sph\n");

    i=NIIK_IMINMAX((int)floor( (pmin.psi - bb->origin.psi) / bb->delta + 0.5 ), 0, bb->psi_dim-1);
    j=NIIK_IMINMAX((int)floor( (pmin.the - bb->origin.the) / bb->delta + 0.5 ), 0, bb->the_dim-1);
    
    index = i + j * bb->psi_dim ;

    if(verbose>3)
      fprintf(stdout,"      %3f %3f %3f %3f %3i %3i | %8i | %8i| %4i %4i\n",
              ((pmin.psi - bb->origin.psi) / bb->delta ),
              ((pmin.the - bb->origin.the) / bb->delta ),
              pmin.psi, pmin.the,
              i, j,
              index, num,
              bb->psi_dim, bb->the_dim
             );

    if(bb->ndata[index] != 0) {
      fprintf(stderr,"Found cell already populated: %f,%f %i,%i %f  o(%f,%f)\n",pmin.psi,pmin.the, i,j, bb->delta, bb->origin.psi,bb->origin.the);
      abort();
    }

    bb->ndata[index] = num;
    bb->data[index] = (kface **)calloc(num,sizeof(kface *));

    for(n=0; n<num; n++) {
      bb->data[index][n] = facelist[n];
    }
    return 1;
  }

  /* make four lists */
  newlist = (kface ***)calloc(4,sizeof(kface **));
  nlist   = (int *)    calloc(4,sizeof(int));

  for(n=0; n<4; n++)  {
    if(verbose>2) fprintf(stdout,"    calloc %i\n",n);
    newlist[n] = (kface **)calloc(num,sizeof(kface *));
    nlist[n] = 0;
  }

  psi_mid = (pmin.psi + pmax.psi)/2.0;
  the_mid = (pmin.the + pmax.the)/2.0;

  /* put each face into appropriate octree component */
  if(verbose>3) fprintf(stdout,"    go thru each triangle\n");

  for(n=0; n<num; n++) {
    kface *f = facelist[n];
    niiksph fmin,fmax;
    face_sph_bounds(f, &fmin, &fmax);
    /* ROI 0: <psi_mid, <the_mid */
    if( fmin.psi <= psi_mid &&
        fmin.the <= the_mid )
      newlist[0][nlist[0]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];
    niiksph fmin,fmax;
    face_sph_bounds(f, &fmin, &fmax);
    /* ROI 1: >psi_mid, <=the_mid */
    if( fmax.psi > psi_mid  &&
        fmin.the <= the_mid )
      newlist[1][nlist[1]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];
    niiksph fmin,fmax;
    face_sph_bounds(f, &fmin, &fmax);

    /* ROI 2: <=psi_mid, >the_mid */
    if( fmin.psi <= psi_mid &&
        fmax.the > the_mid )
      newlist[2][nlist[2]++] = f;
  }

  for(n=0; n<num; n++) {
    kface *f = facelist[n];
    niiksph fmin,fmax;
    face_sph_bounds(f, &fmin, &fmax);

    /* ROI 3: >psi_mid, >the_mid */
    if( fmax.psi > psi_mid &&
        fmax.the > the_mid )
      newlist[3][nlist[3]++] = f;
  }

  if(verbose>3) {
    niiksph p,q;
    for(n=0; n<4; n++) {
      if(n%2)     {
        p.psi=psi_mid;
        q.psi=pmax.psi;
      } else        {
        p.psi=pmin.psi;
        q.psi=psi_mid;
      }

      fprintf(stdout,"    %2i %9i (%5.1f %5.1f) (%5.1f %5.1f)\n",
              n,nlist[n],
              p.psi,p.the,
              q.psi,q.the);
    }
  }

  for(n=0; n<4; n++) {
    niiksph p,q;
    if(nlist[n]!=0) {

      if(n&1)     {
        p.psi=psi_mid;
        q.psi=pmax.psi;
      } else        {
        p.psi=pmin.psi;
        q.psi=psi_mid;
      }

      if(n&2)     {
        p.the=the_mid;
        q.the=pmax.the;
      } else        {
        p.the=pmin.the;
        q.the=the_mid;
      }

      off_create_bbox_sph_octree(bb, newlist[n], nlist[n], cdepth-1, p, q);
    }
  }

  for(n=0; n<4; n++)
    free(newlist[n]);

  free(newlist);
  free(nlist);

  return 1;
} /* int off_create_bbox_sph_octree(bbox_sph *bb,kface **facelist,int num,int cdepth,niikpt pmin,niikpt pmax) */

void off_display_bbox_sph_info(bbox_sph *bb) {
  int n,num;
  fprintf(stdout,"bbox_sph info\n");
  fprintf(stdout,"  max = %0.4f   delta = %0.6f\n",bb->bmax,bb->delta);
  fprintf(stdout,"  dim = %i %i |  nvox=%i\n",bb->psi_dim,bb->the_dim,bb->nvox);
  for(n=num=0; n<bb->nvox; n++) {
    if(bb->ndata[n]>0) num++;
  }
  fprintf(stdout,"  ndata  %i has at least one triangle\n",num);
}


void off_calc_facelist_bounds_sph(kface **facelist,int num,niiksph *pmin,niiksph *pmax) {
  int i;
  /*TODO: paralelelize?*/
  if(num<1) return;
  face_sph_bounds(facelist[0],pmin,pmax);
  for(i=1; i<num; i++) {
    niiksph fmin,fmax;
    face_sph_bounds(facelist[i],&fmin,&fmax);
    *pmin=niiksph_min(*pmin,fmin);
    *pmax=niiksph_max(*pmax,fmax);
  }
}


/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
