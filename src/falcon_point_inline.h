#pragma once

#ifndef __FALCON_POINT_INLINE_H__
#define __FALCON_POINT_INLINE_H__


/*point operations*/
inline static niikpt niikpt_zero()   {
  niikpt p;
  p.x=p.y=p.z=p.w=0;
  return p;
}
inline static niikpt niikpt_nan()    {
  niikpt p;
  p.x=p.y=p.z=p.w=FP_NAN;
  return p;
}
inline static niikpt niikpt_problem() {
  niikpt p;
  p.x=p.y=p.z=p.w=NIIKMAX;
  return p;
}
inline static niikpt niikpt_array(double* a) {
  niikpt p;
  p.x=a[0];
  p.y=a[1];
  p.z=a[2];
  p.w=0.0;
  return p;
}


inline static niikpt niikpt_unit(niikpt p) {
  double d;
  if(fabs(p.x)+fabs(p.y)+fabs(p.z)<1e-10) return niikpt_zero(); /*maybe specify 1e-10 as epsilon ?*/
  d = niikpt_mag(p);
  p.x/=d;
  p.y/=d;
  p.z/=d;
  return p;
}

inline static niikpt niikpt_unit_normal(niikpt p,niikpt q,niikpt r) {
  niikpt v1,v2,n;
  v1=niikpt_sub(q,p);
  if(fabs(v1.x)+fabs(v1.y)+fabs(v1.z)<1e-12) return niikpt_zero();
  v2=niikpt_sub(r,p);
  if(fabs(v2.x)+fabs(v2.y)+fabs(v2.z)<1e-12) return niikpt_zero();
  n=niikpt_cross(niikpt_unit(v1),niikpt_unit(v2));
  if(fabs(n.x)+fabs(n.y)+fabs(n.z)<1e-12) return niikpt_zero();
  return niikpt_unit(n);
}


inline static niikpt niikpt_add(niikpt p,niikpt q) {
  p.x+=q.x;
  p.y+=q.y;
  p.z+=q.z;
  return p;
}
inline static niikpt niikpt_avg(niikpt p,niikpt q) {
  p.x=(p.x+q.x)/2.0;
  p.y=(p.y+q.y)/2.0;
  p.z=(p.z+q.z)/2.0;
  return p;
}

inline static niikpt niikpt_wavg(niikpt p,niikpt q,double wp) {
  double wq;
  wq = 1.0 - wp;
  p.x = p.x*wp + q.x*wq;
  p.y = p.y*wp + q.y*wq;
  p.z = p.z*wp + q.z*wq;
  return p;
}
inline static niikpt niikpt_avg3(niikpt p,niikpt q,niikpt r) {
  p.x=(p.x+q.x+r.x)/3.0;
  p.y=(p.y+q.y+r.y)/3.0;
  p.z=(p.z+q.z+r.z)/3.0;
  return p;
}
inline static niikpt niikpt_wavg3(niikpt p,niikpt q,niikpt r, double wp, double wq, double wr) {
  double s;
  s=wp+wq+wr;
  wp/=s;
  wq/=s;
  wr/=s;
  p.x=wp*p.x+wq*q.x+wr*r.x;
  p.y=wp*p.y+wq*q.y+wr*r.y;
  p.z=wp*p.z+wq*q.z+wr*r.z;
  return p;
}
inline static niikpt niikpt_sub(niikpt p,niikpt q) {
  p.x-=q.x;
  p.y-=q.y;
  p.z-=q.z;
  return p;
}
inline static niikpt niikpt_mul(niikpt p,niikpt q) {
  p.x*=q.x;
  p.y*=q.y;
  p.z*=q.z;
  return p;
}

inline static niikpt niikpt_div(niikpt p,niikpt q) {
#ifdef _DEBUG
  if(q.x==0) {
    fprintf(stderr,"WARNING: division by 0 x\n");
  }
  if(q.y==0) {
    fprintf(stderr,"WARNING: division by 0 y\n");
  }
  if(q.z==0) {
    fprintf(stderr,"WARNING: division by 0 z\n");
  }
#endif
  p.x/=q.x;
  p.y/=q.y;
  p.z/=q.z;
  return p;
}

inline static double niikpt_dot(niikpt p,niikpt q) {
  return p.x*q.x + p.y*q.y + p.z*q.z;
}
inline static double niikpt_mag2(niikpt p) {
  return p.x*p.x + p.y*p.y + p.z*p.z;
}
inline static double niikpt_mag (niikpt p) {
  return sqrt(niikpt_mag2(p));
}
inline static double niikpt_distance2(niikpt p,niikpt q) {
  return NIIK_SQ(p.x-q.x)+NIIK_SQ(p.y-q.y)+NIIK_SQ(p.z-q.z);
}
inline static double niikpt_distance(niikpt p,niikpt q) {
  return sqrt(niikpt_distance2(p,q));
}
inline static niikpt niikpt_cross(niikpt p,niikpt q)   {
  niikpt r;
  r.x = p.y*q.z - p.z*q.y;
  r.y = p.z*q.x - p.x*q.z;
  r.z = p.x*q.y - p.y*q.x;
  r.w=0;
  return r;
}

inline static niikpt niikpt_kmul(niikpt p,double k)                {
  p.x*=k;
  p.y*=k;
  p.z*=k;
  return p;
}
inline static niikpt niikpt_move_normal(niikpt p,niikpt n,double k) {
  p.x+=k*n.x;
  p.y+=k*n.y;
  p.z+=k*n.z;
  return p;
}


inline static double niikpt_angle_between_vectors(niikpt v1,niikpt v2) {
  return NIIK_RAD2DEGREE(acos(niikpt_dot(v1,v2)/(1e-12+niikpt_mag(v1)*niikpt_mag(v2))));
}
inline static double niikpt_area2(niikpt p1,niikpt p2,niikpt p3) {
  return niikpt_mag(niikpt_cross(niikpt_sub(p3,p1),niikpt_sub(p3,p2)));
}

inline static niikpt niikpt_min(niikpt p,niikpt q) {
  p.x=(p.x<q.x)?p.x:q.x;
  p.y=(p.y<q.y)?p.y:q.y;
  p.z=(p.z<q.z)?p.z:q.z;
  return p;
}
inline static niikpt niikpt_max(niikpt p,niikpt q) {
  p.x=(p.x>q.x)?p.x:q.x;
  p.y=(p.y>q.y)?p.y:q.y;
  p.z=(p.z>q.z)?p.z:q.z;
  return p;
}

/* for determinants (det2, det3, det4) matrix is like this
 * | a b |
 * | c d |
 *
 * | a b c |
 * | d e f |
 * | g h i |
 *
 * | a b c d |
 * | e f g h |
 * | i j k l |
 * | m n o p |
 *
 */
inline static double niikpt_det2(double a, double b, double c, double d) {
  return a*d - b*c;
}

inline static double niikpt_det3(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
  return a*niikpt_det2(e,f,h,i) - b*niikpt_det2(d,f,g,i) + c*niikpt_det2(d,e,g,h);
  /*double out=0;
    if(fabs(a)>1e-20) out =a*niikpt_det2(e,f,h,i);
    if(fabs(b)>1e-20) out-=b*niikpt_det2(d,f,g,i);
    if(fabs(c)>1e-20) out+=c*niikpt_det2(d,e,g,h);
    return out; */
}

inline static double niikpt_det4(double a, double b, double c, double d, double e, double f, double g, double h,
                                 double i, double j, double k, double l,  double m, double n, double o, double p ) {
  return a*niikpt_det3(f,g,h,j,k,l,n,o,p) - b*niikpt_det3(e,g,h,i,k,l,m,o,p) + c*niikpt_det3(e,f,h,i,j,l,m,n,p) - d*niikpt_det3(e,f,g,i,j,k,m,n,o);
  /*  double out=0;
      if(fabs(a)>1e-20) out =a*niikpt_det3(f,g,h,j,k,l,n,o,p);
  if(fabs(b)>1e-20) out-=b*niikpt_det3(e,g,h,i,k,l,m,o,p);
  if(fabs(c)>1e-20) out+=c*niikpt_det3(e,f,h,i,j,l,m,n,p);
  if(fabs(d)>1e-20) out-=d*niikpt_det3(e,f,g,i,j,k,m,n,o);
  return out; */

}

inline static niikpt niikpt_affine_transform( niikmat *m, niikpt p ) {
  niikpt q;
  /* assumption: m is good
     if(m==NULL) return niikpt_problem();
     else if(m->row!=4) return niikpt_problem();
     else if(m->col!=4) return niikpt_problem(); */
  q.x=m->m[0][0]*p.x+m->m[0][1]*p.y+m->m[0][2]*p.z+m->m[0][3];
  q.y=m->m[1][0]*p.x+m->m[1][1]*p.y+m->m[1][2]*p.z+m->m[1][3];
  q.z=m->m[2][0]*p.x+m->m[2][1]*p.y+m->m[2][2]*p.z+m->m[2][3];
  q.w=0;
  return q;
}

inline static niikpt niikpt_affine_transform_m44( mat44 m, niikpt p ) {
  niikpt q;
  q.x=m.m[0][0]*p.x+m.m[0][1]*p.y+m.m[0][2]*p.z+m.m[0][3];
  q.y=m.m[1][0]*p.x+m.m[1][1]*p.y+m.m[1][2]*p.z+m.m[1][3];
  q.z=m.m[2][0]*p.x+m.m[2][1]*p.y+m.m[2][2]*p.z+m.m[2][3];
  q.w=0;
  return q;
}

/*coordinate conversion function*/
inline static void niik_world_to_ijk(nifti_image *img, const niikpt *pw, niikpt *p) {
  /*use img->sto_xyz*/
  *p=niikpt_affine_transform_m44(img->sto_ijk,*pw);
}

inline static void niik_ijk_to_world(nifti_image *img, const niikpt *p, niikpt *pw) {
  *pw=niikpt_affine_transform_m44(img->sto_xyz,*p);
}

inline static void niik_index_to_world(nifti_image *img, int i, int j, int k, niikpt *p) {
  niikpt pw;
  /*use img->sto_ijk*/
  pw.x = i ;
  pw.y = j ;
  pw.z = k ;

  *p=niikpt_affine_transform_m44(img->sto_xyz,pw);
}

inline static int niik_ijk_to_offset(nifti_image *img, int i, int j, int k) {
  return i + j*img->nx + k*img->nx*img->ny;
}


inline static void niik_ijk_round(const niikpt *p,int* coord) {
  coord[1] = (int)floor(p->x+0.5);
  coord[2] = (int)floor(p->y+0.5);
  coord[3] = (int)floor(p->z+0.5);
}


#endif /*__FALCON_POINT_INLINE_H__*/
/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
