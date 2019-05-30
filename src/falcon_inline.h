#pragma once

#ifndef __FALCON_INLINE_H__
#define __FALCON_INLINE_H__


static inline char NIIK_CMIN(char a,char b)                             {
  if(a<b) return a;
  return b;
}
static inline int  NIIK_IMIN(int a,int b)                               {
  if(a<b) return a;
  return b;
}
static inline long NIIK_LMIN(long a,long b)                             {
  if(a<b) return a;
  return b;
}
static inline float  NIIK_FMIN(float a,float b)                         {
  if(a<b) return a;
  return b;
}
static inline double NIIK_DMIN(double a,double b)                       {
  if(a<b) return a;
  return b;
}
static inline unsigned char NIIK_UCMIN(unsigned char a,unsigned char b) {
  if(a<b) return a;
  return b;
}
static inline unsigned int  NIIK_UIMIN(unsigned int a,unsigned int b)   {
  if(a<b) return a;
  return b;
}
static inline unsigned long NIIK_ULMIN(unsigned long a,unsigned long b) {
  if(a<b) return a;
  return b;
}

static inline char NIIK_CMAX(char a,char b)                             {
  if(a<b) return b;
  return a;
}
static inline int  NIIK_IMAX(int a,int b)                               {
  if(a<b) return b;
  return a;
}
static inline long NIIK_LMAX(long a,long b)                             {
  if(a<b) return b;
  return a;
}
static inline float  NIIK_FMAX(float a,float b)                         {
  if(a<b) return b;
  return a;
}
static inline double NIIK_DMAX(double a,double b)                       {
  if(a<b) return b;
  return a;
}
static inline unsigned char NIIK_UCMAX(unsigned char a,unsigned char b) {
  if(a<b) return b;
  return a;
}
static inline unsigned int  NIIK_UIMAX(unsigned int a,unsigned int b)   {
  if(a<b) return b;
  return a;
}
static inline unsigned long NIIK_ULMAX(unsigned long a,unsigned long b) {
  if(a<b) return b;
  return a;
}

static inline char   NIIK_CMINMAX(char v,char min,char max)    {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline int    NIIK_IMINMAX(int v,int min,int max)       {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline long   NIIK_LMINMAX(long v,long min,long max)    {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline float  NIIK_FMINMAX(float v,float min,float max) {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline double NIIK_DMINMAX(double v,double min,double max) {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline unsigned char NIIK_UCMINMAX(unsigned char v,unsigned char min,unsigned char max) {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline unsigned int  NIIK_UIMINMAX(unsigned int v, unsigned int min, unsigned int max)  {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}
static inline unsigned long NIIK_ULMINMAX(unsigned long v,unsigned long min,unsigned long max) {
  if(v<min) return min;
  else if(v>max) return max;
  return v;
}

static inline void NIIK_CSWAP(char *a,char *b) {
  char c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_ISWAP(int  *a,int  *b) {
  int  c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_LSWAP(long *a,long *b) {
  long c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_UCSWAP(unsigned char *a,unsigned char *b) {
  unsigned char c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_UISWAP(unsigned int  *a,unsigned int  *b) {
  unsigned int  c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_ULSWAP(unsigned long *a,unsigned long *b) {
  unsigned long c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_FSWAP(float  *a,float  *b) {
  float  c = *a;
  *a=*b;
  *b=c;
}
static inline void NIIK_DSWAP(double *a,double *b) {
  double c = *a;
  *a=*b;
  *b=c;
}

static inline double NIIK_ABS_DMIN(double a,double b) {
  if(fabs(a)<=fabs(b)) return a;
  return b;
}
static inline double NIIK_ABS_DMAX(double a,double b) {
  if(fabs(a)>=fabs(b)) return a;
  return b;
}

static inline double NIIK_SQ(double v)                   {
  return v*v;
}
static inline double NIIK_SSQ(double v,double w)         {
  return sqrt(v*v+w*w);
}
static inline double NIIK_GaussPDF(double val,double std) {
  return exp(-NIIK_SQ(val/std)/2.0);
}


/* conversion functions */
static inline double NIIK_FWHM2SIGMA(double FWHM)   {
  return FWHM/2.35482004503095;
}
static inline double NIIK_SIGMA2FWHM(double sigma)  {
  return sigma*2.35482004503095;
}
static inline double NIIK_RAD2DEGREE(double radian) {
  return radian*57.295779513082320876798154814096;
}
static inline double NIIK_DEGREE2RAD(double degree) {
  return degree/57.295779513082320876798154814096;
}




#endif //__FALCON_INLINE_H__

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/