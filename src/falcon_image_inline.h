#pragma once

#ifndef __FALCON_IMAGE_INLINE_H__
#define __FALCON_IMAGE_INLINE_H__

/* image accees functions*/
static inline double niik_image_get_voxel(nifti_image *img,int p);

static inline double niik_image_voxel_get(nifti_image *img,int p) {
  return niik_image_get_voxel(img,p);
}

static inline double niik_image_get_voxel(nifti_image *img,int p) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  const int  verbose=0;
#ifdef _DEBUG
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_get_voxel] ERROR: img is null\n");
    return NIIKMAX;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_get_voxel] ERROR: img->data is null\n");
    return NIIKMAX;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxel) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return NIIKMAX;
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxel) switch\n");
#endif
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    return (double)cptr[p];
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    return (double)sptr[p];
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    return (double)iptr[p];
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    return (double)lptr[p];
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    return (double)ucptr[p];
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    return (double)usptr[p];
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    return (double)uiptr[p];
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    return (double)ulptr[p];
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    return (double)fptr[p];
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    return (double)dptr[p];
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    return (double)ldptr[p];
  case NIFTI_TYPE_COMPLEX64:
    return NIIKMAX;
  case NIFTI_TYPE_COMPLEX128:
    return NIIKMAX;
  case NIFTI_TYPE_COMPLEX256:
    return NIIKMAX;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
#ifdef _DEBUG
    fprintf(stderr,"ERROR: unknown datatype\n");
#endif
    return NIIKMAX;
  }
#ifdef _DEBUG
  if(verbose) fprintf(stderr,"-d (niik_image_get_voxel) return error\n");
#endif
  return NIIKMAX;
}

static inline int niik_image_set_voxel(nifti_image *img,int p,double d) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  const int verbose=0;
#ifdef _DEBUG
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_set_voxel] ERROR: img is null\n");
    return 0;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_set_voxel] ERROR: img->data is null\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_set_voxel) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_set_voxel) switch\n");
#endif
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    cptr[p]=(char )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_CMAX));
    return 1;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    sptr[p]=(short)floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_SMAX));
    return 1;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    iptr[p]=(int  )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_IMAX));
    return 1;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    lptr[p]=(long )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_LMAX));
    return 1;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    ucptr[p]=(unsigned char )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_UCMAX));
    return 1;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    usptr[p]=(unsigned short)floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_USMAX));
    return 1;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    uiptr[p]=(unsigned int  )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_UIMAX));
    return 1;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    ulptr[p]=(unsigned long )floor(NIIK_DMINMAX(d+0.5,0,NIIK_VAL_ULMAX));
    return 1;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    fptr[p] =d;
    return 1;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    dptr[p] =d;
    return 1;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    ldptr[p]=d;
    return 1;
  case NIFTI_TYPE_COMPLEX64:
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    return 0;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
#ifdef _DEBUG
    fprintf(stderr,"ERROR: unknown datatype %i %s\n",img->datatype,nifti_datatype_string(img->datatype));
#endif
    return 0;
  }
#ifdef _DEBUG
  if(verbose) fprintf(stderr,"-d (niik_image_set_voxel) return error\n");
#endif
  return 0;
}

static inline int niik_image_add_voxel(nifti_image *img,int p,double d) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  const int verbose=0;
#ifdef _DEBUG
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_add_voxel] ERROR: img is null\n");
    return 0;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_add_voxel] ERROR: img->data is null\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) switch\n");
#endif
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    cptr[p]=(char )floor(NIIK_DMINMAX(cptr[p]+d+0.5,0,NIIK_VAL_CMAX));
    return 1;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    sptr[p]=(short)floor(NIIK_DMINMAX(sptr[p]+d+0.5,0,NIIK_VAL_SMAX));
    return 1;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    iptr[p]=(int  )floor(NIIK_DMINMAX(iptr[p]+d+0.5,0,NIIK_VAL_IMAX));
    return 1;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    lptr[p]=(long )floor(NIIK_DMINMAX(lptr[p]+d+0.5,0,NIIK_VAL_LMAX));
    return 1;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    ucptr[p]=(unsigned char )floor(NIIK_DMINMAX(ucptr[p]+d+0.5,0,NIIK_VAL_UCMAX));
    return 1;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    usptr[p]=(unsigned short)floor(NIIK_DMINMAX(usptr[p]+d+0.5,0,NIIK_VAL_USMAX));
    return 1;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    uiptr[p]=(unsigned int  )floor(NIIK_DMINMAX(uiptr[p]+d+0.5,0,NIIK_VAL_UIMAX));
    return 1;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    ulptr[p]=(unsigned long )floor(NIIK_DMINMAX(ulptr[p]+d+0.5,0,NIIK_VAL_ULMAX));
    return 1;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    fptr[p]=fptr[p]+d;
    return 1;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    dptr[p]=dptr[p]+d;
    return 1;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    ldptr[p]=ldptr[p]+d;
    return 1;
  case NIFTI_TYPE_COMPLEX64:
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    return 0;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
#ifdef _DEBUG
    fprintf(stderr,"ERROR: unknown datatype\n");
#endif
    return 0;
  }
#ifdef _DEBUG
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) return error\n");
#endif
  return 0;
} /* niik_image_add_voxel */

static inline int niik_image_mul_voxel(nifti_image *img,int p,double val) {
  unsigned char *ucptr;
  char *cptr;
  unsigned short *usptr;
  short *sptr;
  unsigned int *uiptr;
  int *iptr;
  float *fptr;
  double *dptr;
  long *lptr ;
  unsigned long *ulptr;
  long double *ldptr;
  const int verbose=0;
#ifdef _DEBUG
  if(img      ==NULL) {
    fprintf(stderr,"[niik_image_mul_voxel] ERROR: img is null\n");
    return 0;
  }
  if(img->data==NULL) {
    fprintf(stderr,"[niik_image_mul_voxel] ERROR: img->data is null\n");
    return 0;
  }
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) img %s\n",nifti_datatype_string(img->datatype));
  if(img->datatype == NIFTI_TYPE_COMPLEX64 ||
      img->datatype == NIFTI_TYPE_COMPLEX128 ||
      img->datatype == NIFTI_TYPE_COMPLEX256 ||
      img->datatype == NIFTI_TYPE_RGBA32 ) return 0;
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) switch\n");
#endif
  switch(img->datatype) {
  case NIFTI_TYPE_INT8:
    cptr=(char  *)img->data;
    cptr[p]=(char )floor(NIIK_DMINMAX(cptr[p]*val+0.5,0,NIIK_VAL_CMAX));
    return 1;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    sptr[p]=(short)floor(NIIK_DMINMAX(sptr[p]*val+0.5,0,NIIK_VAL_SMAX));
    return 1;
  case NIFTI_TYPE_INT32:
    iptr=(int   *)img->data;
    iptr[p]=(int  )floor(NIIK_DMINMAX(iptr[p]*val+0.5,0,NIIK_VAL_IMAX));
    return 1;
  case NIFTI_TYPE_INT64:
    lptr=(long  *)img->data;
    lptr[p]=(long )floor(NIIK_DMINMAX(lptr[p]*val+0.5,0,NIIK_VAL_LMAX));
    return 1;
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char  *)img->data;
    ucptr[p]=(unsigned char )floor(NIIK_DMINMAX(ucptr[p]*val+0.5,0,NIIK_VAL_UCMAX));
    return 1;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    usptr[p]=(unsigned short)floor(NIIK_DMINMAX(usptr[p]*val+0.5,0,NIIK_VAL_USMAX));
    return 1;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int   *)img->data;
    uiptr[p]=(unsigned int  )floor(NIIK_DMINMAX(uiptr[p]*val+0.5,0,NIIK_VAL_UIMAX));
    return 1;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long  *)img->data;
    ulptr[p]=(unsigned long )floor(NIIK_DMINMAX(ulptr[p]*val+0.5,0,NIIK_VAL_ULMAX));
    return 1;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float  *)img->data;
    fptr[p]=fptr[p]*val;
    return 1;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    dptr[p]=dptr[p]*val;
    return 1;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    ldptr[p]=ldptr[p]*val;
    return 1;
  case NIFTI_TYPE_COMPLEX64:
    return 0;
  case NIFTI_TYPE_COMPLEX128:
    return 0;
  case NIFTI_TYPE_COMPLEX256:
    return 0;
  case NIFTI_TYPE_RGB24:
  case NIFTI_TYPE_RGBA32:
  default:
#ifdef _DEBUG
    fprintf(stderr,"ERROR: unknown datatype\n");
#endif
    return 0;
  }
#ifdef _DEBUG
  if(verbose) fprintf(stderr,"-d (niik_image_add_voxel) return error\n");
#endif
  return 0;
} /* niik_image_mul_voxel */


#endif /* __FALCON_IMAGE_INLINE_H__*/

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/