/* Filename:     nifti1_kunio_point.c
 * Description:  point functions
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 */

#ifndef _FALCON_POINT_C_
#define _FALCON_POINT_C_

#include "falcon.h"

niikpt niikpt_rand() {
  niikpt p;
  p.x=(niik_get_rand()-0.5)*2.0;
  p.y=(niik_get_rand()-0.5)*2.0;
  p.z=(niik_get_rand()-0.5)*2.0;
  p.w=0;
  return p;
}

niikpt niikpt_val(double x, double y, double z, double w) {
  niikpt p;
  p.x=x;
  p.y=y;
  p.z=z;
  p.w=w;
  return p;
}
void niikpt_disp(niikpt p) {
  fprintf(stdout,"%f %f %f %f\n",p.x,p.y,p.z,p.w);
}
char *niikpt_display_xyz(niikpt p) {
  char *str;
  str=(char *)calloc(100,sizeof(char));
  sprintf(str,"%8.5f %8.5f %8.5f",p.x,p.y,p.z);
  return str;
}


/* interpolating on a triangle */
niikpt niikpt_interp_tri(niikpt p1,niikpt p2,niikpt p3,double t1,double t2) {
  niikpt e1,e2;
  e1=niikpt_wavg(p1,p2,t2);
  e2=niikpt_wavg(p1,p3,t2);
  return niikpt_wavg(e1,e2,t1);
}

/* signed distance from the plane to the point */
double niikpt_plane_to_point_distance(niikpt pt,niikpt pt_on_plane,niikpt plane_normal) {
  return ( niikpt_dot(pt,plane_normal) - niikpt_dot(pt_on_plane,plane_normal) ) / sqrt(1e-12+niikpt_dot(plane_normal,plane_normal));
}
double niikpt_plane_to_point_distance2(niikpt pt,niikpt plane_eq) {
  return (plane_eq.x*pt.x + plane_eq.y*pt.y + plane_eq.z*pt.z + plane_eq.w) / niikpt_mag(plane_eq);
}


niikpt niikpt_plane_eq_to_plane_normal(niikpt plane) {
  return niikpt_unit(plane);
}
niikpt niikpt_3pts_to_plane_eq(niikpt p1,niikpt p2,niikpt p3) {
  niikpt n;
  n=niikpt_cross(niikpt_sub(p2,p1),niikpt_sub(p3,p1));
  n.w=-niikpt_dot(p1,n);
  return n;
}

niikpt niikpt_closest_point_on_plane_to_point(niikpt pt,niikpt plane) {
  double d;
  d = niikpt_plane_to_point_distance2(pt,plane);
  return niikpt_move_normal(pt,niikpt_plane_eq_to_plane_normal(plane),-d);
}

/* line = between p1 and p2
 * point = pt */
double niikpt_line_to_point_distance(niikpt pt,niikpt p1,niikpt p2) {
  return niikpt_mag(niikpt_cross(niikpt_sub(pt,p1),niikpt_sub(pt,p2))) / niikpt_mag(niikpt_sub(p2,p1));
}

niikpt niikpt_closest_point_on_line_to_point_distance(niikpt pt,niikpt p1,niikpt p2) {
  double t;
  niikpt v;
  v = niikpt_sub(p2,p1);
  t=-niikpt_dot(niikpt_sub(p1,pt),v) / niikpt_mag2(v);
  return niikpt_move_normal(p1,v,t);
}

/* -finds the intersection of plane (defined by 3 points = p1,p2,p3)
 *  and a line (defined by L1 and L2)
 */
niikpt niikpt_line_plane_intersection(niikpt p1,niikpt p2,niikpt p3,niikpt L1,niikpt L2) {
  niikpt p,v;
  double t;
  v = niikpt_sub(L2,L1);
  t = niikpt_det4(1,1,1,0,p1.x,p2.x,p3.x,v.x,p1.y,p2.y,p3.y,v.y,p1.z,p2.z,p3.z,v.z);
  if(fabs(t)<1e-12) {
    return niikpt_problem();
  }
  t = niikpt_det4(1,1,1,1,p1.x,p2.x,p3.x,L1.x,p1.y,p2.y,p3.y,L1.y,p1.z,p2.z,p3.z,L1.z) / t;
  p = niikpt_move_normal(L1,v,-t);
  return p;
}


niikpt niikpt_closest_point_on_triangle_to_point(niikpt pt,niikpt p1,niikpt p2,niikpt p3) {
  /* triangle = p1,p2,p3
   * point = pt */
  char fcname[64]="niikpt_closest_point_on_triangle_to_point";
  double d;
  niikpt v[3],plane,norm;
  int idx=0;
  /* get the plane normal and plane eq
     plane = niikpt_3pts_to_plane_eq(p1,p2,p3); */
  v[0]=niikpt_sub(p2,p1);
  v[1]=niikpt_sub(p3,p2);
  v[2]=niikpt_sub(p1,p3);
  norm = plane = niikpt_cross(v[0],niikpt_sub(p3,p1));
  plane.w = -niikpt_dot(p1,norm);
  norm = niikpt_unit(norm);
  /* distance from triangle's plane to pt */
  d = niikpt_plane_to_point_distance2(pt,plane);
  /* projected point on plane */
  pt = niikpt_move_normal(pt,norm,-d);
  /* direction relative to the triangle edges */
  if(niikpt_dot(norm,niikpt_cross(v[0],niikpt_sub(pt,p1)))>0) idx+=1;
  if(niikpt_dot(norm,niikpt_cross(v[1],niikpt_sub(pt,p2)))>0) idx+=2;
  if(niikpt_dot(norm,niikpt_cross(v[2],niikpt_sub(pt,p3)))>0) idx+=4;
  /* fprintf(stdout,"idx = %i\n",idx); */
  switch(idx) {
  case 1: /* close to p3 */
    return p3;
  case 2: /* close to p1 */
    return p1;
  case 3: /* close to line p1-p3 */
    return niikpt_closest_point_on_line_to_point_distance(pt,p1,p3);
  case 4: /* close to p2 */
    return p2;
  case 5: /* close to line p2-p3 */
    return niikpt_closest_point_on_line_to_point_distance(pt,p2,p3);
  case 6: /* close to line p1-p2 */
    return niikpt_closest_point_on_line_to_point_distance(pt,p1,p2);
  case 7: /* in the triangle */
    return pt;
  default:
    fprintf(stderr,"[%s] ERROR: unknown error \n",fcname);
    fprintf(stderr,"[%s] ERROR: niikpt_point_to_triangle_distance\n",fcname);
    return niikpt_zero();
  }
  return niikpt_zero();
}

niikpt *niikpt_alloc(int num) {
  niikpt *p;
  if((p=(niikpt *)calloc(num,sizeof(niikpt)))==NULL) {
    fprintf(stderr,"ERROR: calloc for niikpt_alloc \n");
    return NULL;
  }
  return p;
}

niikpt **niikpt_matrix(int col,int row) {
  niikpt **p;
  int n;
  if((p=(niikpt **)calloc(col,sizeof(niikpt *)))==NULL) {
    fprintf(stderr,"ERROR: calloc for niikpt_alloc \n");
    return NULL;
  }
  for(n=0; n<col; n++) {
    if((p[n] = niikpt_alloc(row))==NULL) {
      fprintf(stderr,"ERROR: niikpt_alloc, %i\n",n);
      return NULL;
    }
  }
  return p;
}

niikpt niikpt_image_get_pixdim(nifti_image *img) {
  niikpt delta;
  delta.x=img->pixdim[1];
  delta.y=img->pixdim[2];
  delta.z=img->pixdim[3];
  delta.w=0;
  return delta;
}


/* -returns 1 is a point is on triangle
 * -assumes that pt is on the same plane
 * -triangle is defined by t1, t2, t3
 * -edges are defined by v1=t2-t1, v2=t3-t2, and v3=t1-t3
 * -triangle normal is normal
 */
int niikpt_point_is_on_triangle(niikpt pt,niikpt t1,niikpt t2,niikpt t3,niikpt v1,niikpt v2,niikpt v3,niikpt normal) {
  niikpt v;
  /* vertex 1 */
  v=niikpt_sub(pt,t1);
  if(niikpt_mag(v)<1e-8) return 1;
  if(niikpt_dot(niikpt_cross(v1,niikpt_unit(v)),normal) < -0.3 ) return 0;
  /* vertex 2 */
  v=niikpt_sub(pt,t2);
  if(niikpt_mag(v)<1e-8) return 1;
  if(niikpt_dot(niikpt_cross(v2,niikpt_unit(v)),normal) < -0.3 ) return 0;
  /* vertex 3 */
  v=niikpt_sub(pt,t3);
  if(niikpt_mag(v)<1e-8) return 1;
  if(niikpt_dot(niikpt_cross(v3,niikpt_unit(v)),normal) < -0.3 ) return 0;
  return 1;
}

/* pixel spacing is multiplied */
int niikpt_image_get_centroid_update(nifti_image * img,nifti_image *maskimg,niikpt *pt) {
  const char *fcname="niikpt_image_get_centroid_update";
  niikpt p;
  int i,j,k,n;
  double dval,ds[9];
  unsigned char *bimg;
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }
  if(maskimg==NULL) {
    bimg = (unsigned char *)calloc(img->nvox,sizeof(char));
    for(i=0; i<img->nvox; i++) bimg[i]=1;
  } else {
    if(niik_image_cmp_dim(img,maskimg)!=0) {
      fprintf(stderr,"[%s] ERROR: image [%ix%ix%i] and mask [%ix%ix%i] don't have the same dimension\n",fcname,
              img->nx,img->ny,img->nz,maskimg->nx,maskimg->ny,maskimg->nz);
      return 0;
    }
    if((bimg = niik_image_get_voxels_as_uint8_vector(maskimg))==NULL) {
      fprintf(stderr,"[%s] ERROR: niik_image_get_voxels_as_uint8_vector\n",fcname);
      return 0;
    }
  }
  ds[0]=ds[1]=ds[2]=ds[3]=0;
  for(k=n=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; n++,i++) {
        if(!bimg[n]) continue;
        dval = niik_image_get_voxel(img,n);
        if(niik_check_double_problem(dval)) {
          fprintf(stderr,"ERROR: niik_image_get_voxel, %i,%i,%i \n",i,j,k);
          return 0;
        }
        ds[0]+=dval;
        ds[1]+=dval*i;
        ds[2]+=dval*j;
        ds[3]+=dval*k;
      }
    }
  }
  for(i=1; i<=3; i++) {
    ds[i]/=ds[0];
    /*ds[i]*=img->pixdim[i]; */

  }
  p.x = ds[1];
  p.y = ds[2];
  p.z = ds[3];
  p.w = ds[0];
  free(bimg);
  bimg=NULL;
  niik_ijk_to_world(img,&p,pt);
  /**pt=p;*/
  return 1;
} /* niikpt_image_get_centroid_update */

niikpt niikpt_image_get_centroid(nifti_image * img,nifti_image *maskimg) {
  niikpt p;
  if(!niikpt_image_get_centroid_update(img,maskimg,&p)) {
    fprintf(stderr,"[%s:%i:%s] ERROR: niik_image_get_centroid_update\n",__FILE__,__LINE__,__func__);
    return niikpt_problem();
  }
  return p;
} /* niikpt_image_get_centroid */

double niikpt_image_get_bounding_sphere(nifti_image * img,double threshold,niikpt ctr) {
  const char *fcname="niikpt_image_get_bounding_sphere";
  niikpt p,w;
  int i,j,k,n;
  double dval;
  unsigned char *bimg;
  double max_r2=0.0;
  if(img==NULL)     {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }

  for(k=n=0; k<img->nz; k++) {
    for(j=0; j<img->ny; j++) {
      for(i=0; i<img->nx; n++,i++) {
        double r2;
        dval=niik_image_get_voxel(img,n);
        if(niik_check_double_problem(dval)) {
          fprintf(stderr,"ERROR: niik_image_get_voxel, %i,%i,%i \n",i,j,k);
          return 0;
        }

        if(dval<threshold) continue;

        p.x = i;
        p.y = j;
        p.z = k;
        niik_ijk_to_world(img,&p,&w);
        r2=niikpt_distance2(ctr,w);
        if(r2>max_r2) max_r2=r2;
      }
    }
  }

  return sqrt(max_r2);
}


int niikpt_image_get_min_ijk_position(nifti_image *maskimg,niikpt *pmin) {
  char fcname[64]="niik_image_get_min_ijk_position";
  int i,j,k,n;
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  pmin->x=pmin->y=pmin->z=1e10;
  for(k=n=0; k<maskimg->nz; k++) {
    for(j=0; j<maskimg->ny; j++) {
      for(i=0; i<maskimg->nx; n++,i++) {
        if(niik_image_get_voxel(maskimg,n)>0) {
          if(pmin->x>i) pmin->x=(double)i;
          if(pmin->y>j) pmin->y=(double)j;
          if(pmin->z>k) pmin->z=(double)k;
        }
      }
    }
  }
  return 1;
} /* niikpt_image_get_min_ijk_position */

int niikpt_image_get_max_ijk_position(nifti_image *maskimg,niikpt *pmax) {
  char fcname[64]="niik_image_get_max_ijk_position";
  int i,j,k,n;
  if(maskimg==NULL) {
    fprintf(stderr,"[%s] ERROR: maskimg is null\n",fcname);
    return 0;
  }
  pmax->x=pmax->y=pmax->z=-1e10;
  for(k=n=0; k<maskimg->nz; k++) {
    for(j=0; j<maskimg->ny; j++) {
      for(i=0; i<maskimg->nx; n++,i++) {
        if(niik_image_get_voxel(maskimg,n)>0) {
          if(pmax->x<i) pmax->x=(double)i;
          if(pmax->y<j) pmax->y=(double)j;
          if(pmax->z<k) pmax->z=(double)k;
        }
      }
    }
  }
  return 1;
} /* niikpt_image_get_max_ijk_position */


/******************************************
 *
 * spherical coordinates
 *
 *******************************************/
niikpt niikpt_euc2sph(niikpt p) {
  niikpt q;
  q.x = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  q.y = atan2(p.y,p.x);
  q.z = acos(p.z/q.x); //theta
  q.w = 0;
  return q;
}

niikpt niikpt_sph2euc(niikpt p) {
  niikpt q;
  q.x = p.x * sin(p.z) * cos(p.y);
  q.y = p.x * sin(p.z) * sin(p.y);
  q.z = p.x * cos(p.z);
  q.w = 0;
  return q;
}





#endif
/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/