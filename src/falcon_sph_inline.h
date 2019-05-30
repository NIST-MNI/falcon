#pragma once

#ifndef __FALCON_SPH_INLINE_H__
#define __FALCON_SPH_INLINE_H__


inline static niiksph niiksph_zero()   {
  niiksph p;
  p.psi=p.the=0;
  return p;
}
inline static niiksph niiksph_nan()    {
  niiksph p;
  p.psi=p.the=FP_NAN;
  return p;
}
inline static niiksph niiksph_problem() {
  niiksph p;
  p.psi=p.the=NIIKMAX;
  return p;
}

inline static niiksph niiksph_norm(niiksph p) {
  if(p.the<0.0) p.the+=NIIK_PI2;
  if(p.psi<0.0) p.psi+=NIIK_PI2;

  if(p.the>NIIK_PI2) p.the-=NIIK_PI2;
  if(p.psi>NIIK_PI2) p.psi-=NIIK_PI2;

  /*TODO: figure out how to normalize ambiguity between the and psi */
  return p;
}



inline static niiksph niiksph_norm_sym(niiksph p) {
  if(p.the< (-NIIK_PI)) p.the+=NIIK_PI2;
  if(p.psi< (-NIIK_PI)) p.psi+=NIIK_PI2;

  if(p.the> NIIK_PI) p.the-=NIIK_PI2;
  if(p.psi> NIIK_PI) p.psi-=NIIK_PI2;

  /*TODO: figure out how to normalize ambiguity between the and psi */
  return p;
}

inline static niiksph niiksph_add(niiksph p,niiksph q) {
  p.psi+=q.psi;
  p.the+=q.the;
  return p;
}
inline static niiksph niiksph_sub(niiksph p,niiksph q) {
  p.psi-=q.psi;
  p.the-=q.the;
  return p;
}
inline static double  niiksph_dot(niiksph p,niiksph q) {
  return p.psi*q.psi+p.the*q.the;
}
inline static niiksph niiksph_avg(niiksph p,niiksph q) {
  p.psi=(p.psi+q.psi)/2.0;
  p.the=(p.the+q.the)/2.0;
  return p;
}

inline static niiksph niiksph_min(niiksph p,niiksph q) {
  p.psi=(p.psi<q.psi)? p.psi:q.psi;
  p.the=(p.the<q.the)? p.the:q.the;
  return p;
}
inline static niiksph niiksph_max(niiksph p,niiksph q) {
  p.psi=(p.psi>q.psi)?p.psi:q.psi;
  p.the=(p.the>q.the)?p.the:q.the;
  return p;
}

inline static niikpt  niiksph_xyz(niiksph p) {
  niikpt pt;
  pt.x = sin(p.the) * cos(p.psi);
  pt.y = sin(p.the) * sin(p.psi);
  pt.z = cos(p.the);

  return pt;
}

inline static niiksph  niikxyz_sph(niikpt pt) {
  niiksph p;
  pt=niikpt_unit(pt);

  p.the=acos(pt.z);
  p.psi=atan2(pt.y,pt.x);

  return niiksph_norm(p);
}

inline static double off_avg_radians(double r1, double r2, double lim)
/* add r1 and r2 but considers the limits like PI,
 * lim is half of the limit */
{
  double d;
  if(fabs(r2-r1)>=lim) {
    d=(r1+r2)/2.0+lim;
  } else {
    d=(r1+r2)/2.0;
  }
  if(d<0) d+=NIIK_PI2;
  if(d>NIIK_PI2) d-=NIIK_PI2;
  return d;
}

inline static niiksph niiksph_avg_lim(niiksph p,niiksph q, double lim) {
  niiksph r;
  r.psi=off_avg_radians(p.psi,q.psi,lim);
  r.the=off_avg_radians(p.the,q.the,lim);
  return r;
}


#endif /*__FALCON_SPH_INLINE_H__*/

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/
