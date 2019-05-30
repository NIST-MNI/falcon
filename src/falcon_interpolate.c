/* Filename:     nifti1_kunio_interpolate.c
 * Description:  interpolation functions
 * Author:       Kunio Nakamura
 * Date:         February 24, 2012
 *
 * important functions
 *
 *   double niik_image_interpolate_3d(nifti_image *img,niikpt p,int interp);
 *   double niik_image_interpolate_3d_ijk(nifti_image *img,niikpt p,int interp);
 *   int niik_image_interpolate_3d_ijk_update(nifti_image *img,niikpt p,int interp,double *v);
 *   int niik_image_interpolate_3d_xyz_update(nifti_image *img,niikpt p,int interp,double *v);
 *
 *
 * 2012-06-17, Kunio
 * -added functions
 *
 *   int niik_image_interp_along_normal(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v,niikvec *out);
 *   int niik_image_interp_along_normal_update(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v);
 */

#ifndef _FALCON_INTERP_C_

#include "falcon.h"

/* internal function */
int niik_image_interpolate_3d_bspline_ijk_update_double_type(nifti_image *bsplimg,niikpt p,double *v);
int niik_image_interpolate_3d_bspline_ijk_update_float_type(nifti_image *bsplimg,niikpt p,double *v);


/******************************************************************
 *
 * general interpolation function for 3d image
 *
 ******************************************************************/

double niik_image_interpolate_3d(nifti_image *img,niikpt p,int interp) {
  double d;
  switch(interp) {
  case NIIK_INTERP_NN:
    d=niik_image_interpolate_3d_nn(img,p);
    return d;
  case NIIK_INTERP_LINEAR:
    d=niik_image_interpolate_3d_linear(img,p);
    return d;
  case NIIK_INTERP_BSPLINE:
    d=niik_image_interpolate_3d_bspline(img,p);
    /*if(niik_check_double_problem(d)){
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline\n");
      return NIIKMAX; }*/
    return d;
  default:
    fprintf(stderr,"ERROR: unknown inteprolation type %i \n",interp);
    break;
  } /* switch (interp) */
  return 0;
} /* int niik_image_interpolate_3d */

double niik_image_interpolate_3d_ijk(nifti_image *img,niikpt p,int interp) {
  double d;
  switch(interp) {
  case NIIK_INTERP_NN:
    d=niik_image_interpolate_3d_nn_ijk(img,p);
    return d;
  case NIIK_INTERP_LINEAR:
    d=niik_image_interpolate_3d_linear_ijk(img,p);
    /*if(niik_check_double_problem(d)){
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear\n");
      return NIIKMAX; }*/
    return d;
  case NIIK_INTERP_BSPLINE:
    d=niik_image_interpolate_3d_bspline_ijk(img,p);
    return d;
  default:
    fprintf(stderr,"ERROR: unknown inteprolation type %i \n",interp);
    break;
  } /* switch (interp) */
  return 0;
} /* int niik_image_interpolate_3d_ijk */

int niik_image_interpolate_3d_ijk_update(nifti_image *img,niikpt p,int interp,double *v) {
  switch(interp) {
  case NIIK_INTERP_NN:
    if(!niik_image_interpolate_3d_nn_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_nn_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  case NIIK_INTERP_LINEAR:
    if(!niik_image_interpolate_3d_linear_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  case NIIK_INTERP_BSPLINE:
    if(!niik_image_interpolate_3d_bspline_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown inteprolation type %i \n",interp);
    break;
  } /* switch (interp) */
  return 1;
} /* int niik_image_interpolate_3d_ijk_update */

int niik_image_interpolate_3d_xyz_update(nifti_image *img,niikpt pw,int interp,double *v) {
  niikpt p;
  niik_world_to_ijk(img,&pw,&p);

  switch(interp) {
  case NIIK_INTERP_NN:
    if(!niik_image_interpolate_3d_nn_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_nn_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  case NIIK_INTERP_LINEAR:
    if(!niik_image_interpolate_3d_linear_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  case NIIK_INTERP_BSPLINE:
    if(!niik_image_interpolate_3d_bspline_ijk_update(img,p,v)) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear_ijk_update(img,pt,v)\n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown inteprolation type %i \n",interp);
    break;
  } /* switch (interp) */
  return 1;
} /* int niik_image_interpolate_3d_ijk_update */


/******************************************************************
 * string for name of interpolation
 *
 ******************************************************************/
const char * niik_interpolate_string(int interp) {
  return niik_interpolate_name(interp);
}

const char * niik_interpolate_name(int interp) {
  switch(interp) {
  case NIIK_INTERP_NN:
    return "nearest neighbor";
  case NIIK_INTERP_LINEAR:
    return "trilinear";
  case NIIK_INTERP_SPLINE:
    return "cubic spline";
  case NIIK_INTERP_BSPLINE:
    return "b-spline";
  case NIIK_INTERP_SINC:
    return "sinc";
  default:
    break;
  }
  return "unknown";
}



/******************************************************************
 *
 * nearest neighbor interpolation
 *
 * 2012-06-18, Kunio
 * -fixed a bug for ijk functions (rounding was not there)
 *
 ******************************************************************/

double niik_image_interpolate_3d_nn(nifti_image *img,niikpt pw) {
  niikpt p;
  niik_world_to_ijk(img,&pw,&p);
  return niik_image_interpolate_3d_nn_ijk(img,p);
}

double niik_image_interpolate_3d_nn_ijk(nifti_image *img,niikpt p) {
  int
  n,voffset,
  coord[7];
  double d=0;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return NIIKMAX;
  }
  if(img->ndim>3) {
    fprintf(stderr,"ERROR: img has more than 3 dimensions, %i\n",img->ndim);
    return NIIKMAX;
  }
  /* get the voxel index while checking for img.dim */
  niik_ijk_round(&p,coord);

  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    if(coord[n]<0) return 0; /* background is zero */
    if(coord[n]>=img->dim[n]) return 0; /* background is zero */
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* get the voxel intensity */
  if(img->dim[0]<=3) {
    d=niik_image_get_voxel(img,coord[0]);
    /*if(niik_check_double_problem(d)) {
      fprintf(stderr,"ERROR: niik_image_get_voxel\n");
      return NIIKMAX; } */
  }
  return d;
} /* niik_image_interpolate_3d_nn_ijk */

int niik_image_interpolate_3d_nn_ijk_update(nifti_image *img,niikpt p,double *v)
/* img is the input image
 * p is the point of interest
 * v is replaced with the output (in a vector format for more than 3d image)
 */
{
  int
  mm,m,n,voffset,
  coord[7];
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  /* get the voxel index while checking for img.dim */
  niik_ijk_round(&p,coord);

  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    if(coord[n]<0)  {
      mm=(img->nvox/img->nx/img->ny/img->nz);
      for(m=0; m<mm; m++) v[m]=0;
      return 1; /* background is zero */
    }

    if(coord[n]>=img->dim[n]) {
      mm=(img->nvox/img->nx/img->ny/img->nz);
      for(m=0; m<mm; m++) v[m]=0;
      return 1; /* background is zero */
    }
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* get the voxel intensity */
  for(m=0,n=coord[0]; n<img->nvox; n+=voffset,m++) {
    v[m]=niik_image_get_voxel(img,n);
  }
  return 1;
} /* niik_image_interpolate_3d_nn_ijk_update */



/*
 * faster version of nearest neighbor interpolation
 *
 * -return 0 with problem
 * -out of image is 0
 * -assumes that img is uint8 type
 * -assumes image is 3D
 */

unsigned char niik_image_interpolate_uint8_image_3d_nn(nifti_image *img,double x, double y,double z) {
  int
  n,voffset,
  coord[7];
  niikpt pw,p;
  unsigned char *bimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  } else if(img->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"ERROR: datatype is not NIFTI_TYPE_UINT8 \n");
    return 0;
  } else if(img->ndim>3) {
    fprintf(stderr,"ERROR: img has more than 3 dimensions, %i\n",img->ndim);
    return 0;
  }
  /* get the voxel index while checking for img.dim */
  pw.x=x;
  pw.y=y;
  pw.z=z;
  niik_world_to_ijk(img,&pw,&p);
  niik_ijk_round(&p,coord);

  /* calculate the index */
  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    if(coord[n]<0) return 0;
    if(coord[n]>=img->dim[n]) return 0;
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* get the voxel intensity */
  bimg = img->data;
  return bimg[coord[0]];
} /* niik_image_interpolate_uint8_image_3d_nn */

double niik_image_interpolate_double_image_3d_nn(nifti_image *img,double x, double y,double z) {
  int
  n,voffset,
  coord[7];
  niikpt pw,p;
  double *bimg;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  } else if(img->datatype!=NIFTI_TYPE_UINT8) {
    fprintf(stderr,"ERROR: datatype is not NIFTI_TYPE_UINT8 \n");
    return 0;
  } else if(img->ndim>3) {
    fprintf(stderr,"ERROR: img has more than 3 dimensions, %i\n",img->ndim);
    return 0;
  }

  /* get the voxel index while checking for img.dim */
  pw.x=x;
  pw.y=y;
  pw.z=z;
  niik_world_to_ijk(img,&pw,&p);
  niik_ijk_round(&p,coord);

  /* calculate the index */
  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    if(coord[n]<0) return 0;
    if(coord[n]>=img->dim[n]) return 0;
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* get the voxel intensity */
  bimg = img->data;
  return bimg[coord[0]];
} /* niik_image_interpolate_double_image_3d_nn */






/******************************************************************
 *
 * trilinear interpolation
 *
 * -assume intensity=0 outside the image
 *
 ******************************************************************/

double niik_image_interpolate_3d_linear(nifti_image *img,niikpt pw) {
  niikpt p;
  niik_world_to_ijk(img,&pw,&p);
  return niik_image_interpolate_3d_linear_ijk(img,p);
}

double niik_image_interpolate_3d_linear_ijk(nifti_image *img,niikpt p) {
  unsigned char  *ucptr;
  unsigned short *usptr;
  unsigned int   *uiptr;
  unsigned long  *ulptr;
  char  *cptr;
  short *sptr;
  int   *iptr;
  long  *lptr;
  float *fptr;
  double *dptr;
  long double *ldptr;
  int
  n,voffset,
  idxlist[9],
  coord[7];
  double
  dval,
  pos[9],
  dd[7],id[7],
  vallist[9];
  /* check inputs */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return NIIKMAX;
  }
  /* get the voxel index while checking for img.dim */
  pos[1]=p.x;
  pos[2]=p.y;
  pos[3]=p.z;
  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    coord[n] = floor(pos[n]); /* floor */
    if(coord[n]<-1) return 0; /* -1 is OK */
    if(coord[n]>=img->dim[n]) return 0;
    dd[n] = pos[n] - coord[n];
    id[n] = 1.0 - dd[n];
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* fprintf(stdout,"%3i %3i %3i   %i \n",coord[1],coord[2],coord[3],coord[0]); */
  idxlist[0]=coord[0];
  idxlist[1]=idxlist[0]+1;
  idxlist[2]=idxlist[1]+img->nx;
  idxlist[3]=idxlist[2]-1;
  idxlist[4]=idxlist[0]+img->nx*img->ny;
  idxlist[5]=idxlist[4]+1;
  idxlist[6]=idxlist[5]+img->nx;
  idxlist[7]=idxlist[6]-1;
  if(coord[1]<0)          idxlist[0]=idxlist[3]=idxlist[4]=idxlist[7]=-1;
  if(coord[1]>=img->nx-1) idxlist[1]=idxlist[2]=idxlist[5]=idxlist[6]=-1;
  if(coord[2]<0)          idxlist[0]=idxlist[1]=idxlist[4]=idxlist[5]=-1;
  if(coord[2]>=img->ny-1) idxlist[2]=idxlist[3]=idxlist[6]=idxlist[7]=-1;
  if(coord[3]<0)          idxlist[0]=idxlist[1]=idxlist[2]=idxlist[3]=-1;
  if(coord[3]>=img->nz-1) idxlist[4]=idxlist[5]=idxlist[6]=idxlist[7]=-1;
  if(img->dim[0]<=3) {
    img->nu=img->nt=1;
  }
  /* 2012-04-09 Kunio
   * -do not use niik_image_get_voxel
   * * * * * * * * * * * * * * * * * * * *
   for(n=0;n<8;n++){
   if(idxlist[n]<0) vallist[n]=0;
   else {
   vallist[n]=niik_image_get_voxel(img,idxlist[n]);
   if(niik_check_double_problem(vallist[n])){
   fprintf(stderr,"ERROR: niik_image_get_voxel\n");
   return NIIKMAX; } } }*/
  switch(img->datatype) {
  case NIFTI_TYPE_UINT8:
    ucptr=(unsigned char *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=ucptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_UINT16:
    usptr=(unsigned short *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=usptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_UINT32:
    uiptr=(unsigned int *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=uiptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_UINT64:
    ulptr=(unsigned long *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=ulptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_INT8:
    cptr=(char *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=cptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_INT16:
    sptr=(short *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=sptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_INT32:
    iptr=(int *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=iptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_INT64:
    lptr=(long *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=lptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    fptr=(float *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=fptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_FLOAT64:
    dptr=(double *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=dptr[idxlist[n]];
      }
    }
    break;
  case NIFTI_TYPE_FLOAT128:
    ldptr=(long double *)img->data;
    for(n=0; n<8; n++) {
      if(idxlist[n]<0) vallist[n]=0;
      else {
        vallist[n]=ldptr[idxlist[n]];
      }
    }
    break;
  default:
    fprintf(stderr,"ERROR: unknown datatype %i\n",img->datatype);
    return NIIKMAX;
  }
  dval =
    vallist[0] * id[1] * id[2] * id[3] +
    vallist[1] * dd[1] * id[2] * id[3] +
    vallist[2] * dd[1] * dd[2] * id[3] +
    vallist[3] * id[1] * dd[2] * id[3] +
    vallist[4] * id[1] * id[2] * dd[3] +
    vallist[5] * dd[1] * id[2] * dd[3] +
    vallist[6] * dd[1] * dd[2] * dd[3] +
    vallist[7] * id[1] * dd[2] * dd[3];
  /* fprintf(stdout,"out = %f \n",dval); */
  return dval;
} /* double niik_image_interpolate_3d_linear_ijk(nifti_image *img,niikpt p) */

int niik_image_interpolate_3d_linear_xyz_update(nifti_image *img,niikpt pw,double *v) {
  niikpt p;
  niik_world_to_ijk(img,&pw,&p);

  if(!niik_image_interpolate_3d_linear_ijk_update(img,p,v)) {
    fprintf(stderr,"ERROR: niik_image_interpolate_3d_linear_ijk_update\n");
    return 0;
  }
  return 1;
}

int niik_image_interpolate_3d_linear_ijk_update(nifti_image *img,niikpt p,double *v)
/* img is the input image
 * p is the point of interest
 * v is replaced with the output (in a vector format for more than 3d image)
 */
{
  int
  m,mm,n,nn,k,voffset,
  idxlist[9],
  coord[7];
  double
  pos[9],
      dd[7],id[7],
      vallist[9];
  /* check inputs */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer \n");
    return 0;
  }
  /* get the voxel index while checking for img.dim */
  pos[1]=p.x;
  pos[2]=p.y;
  pos[3]=p.z;
  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    coord[n] = floor(pos[n]); /* floor */
    if(coord[n]<-1) {
      mm=(img->nvox/img->nx/img->ny/img->nz);
      for(m=0; m<mm; m++) v[m]=0;
      return 1;
    } /* -1 is OK */
    if(coord[n]>=img->dim[n]) {
      mm=(img->nvox/img->nx/img->ny/img->nz);
      for(m=0; m<mm; m++) v[m]=0;
      return 1;
    }
    dd[n] = pos[n] - coord[n];
    id[n] = 1.0 - dd[n];
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /*fprintf(stdout,"%3i %3i %3i   %i \n",coord[1],coord[2],coord[3],coord[0]);
    fprintf(stdout,"voffset %3i\n",voffset);*/
  mm=(img->nvox/img->nx/img->ny/img->nz);
  for(m=0,nn=coord[0]; m<mm; nn+=voffset,m++) {
    idxlist[0]=nn;
    idxlist[1]=idxlist[0]+1;
    idxlist[2]=idxlist[1]+img->nx;
    idxlist[3]=idxlist[2]-1;
    idxlist[4]=idxlist[0]+img->nx*img->ny;
    idxlist[5]=idxlist[4]+1;
    idxlist[6]=idxlist[5]+img->nx;
    idxlist[7]=idxlist[6]-1;
    if(coord[1]<0)          idxlist[0]=idxlist[3]=idxlist[4]=idxlist[7]=-1;
    if(coord[1]>=img->nx-1) idxlist[1]=idxlist[2]=idxlist[5]=idxlist[6]=-1;
    if(coord[2]<0)          idxlist[0]=idxlist[1]=idxlist[4]=idxlist[5]=-1;
    if(coord[2]>=img->ny-1) idxlist[2]=idxlist[3]=idxlist[6]=idxlist[7]=-1;
    if(coord[3]<0)          idxlist[0]=idxlist[1]=idxlist[2]=idxlist[3]=-1;
    if(coord[3]>=img->nz-1) idxlist[4]=idxlist[5]=idxlist[6]=idxlist[7]=-1;
    for(n=k=0; n<8; n++) {
      if(idxlist[n]<0) {
        vallist[n]=0;
        k++;
      } else {
        vallist[n]=niik_image_get_voxel(img,idxlist[n]);
        /*if(niik_check_double_problem(vallist[n])){
          fprintf(stderr,"ERROR: niik_image_get_voxel\n");
          return 0; } */
        /*fprintf(stdout,"m=%i n=%i val=%f\n",m,n,vallist[n]);*/
      }
    }
    if(k==8) {
      v[m]=0;
      break;
    }
    v[m] =
      vallist[0] * id[1] * id[2] * id[3] +
      vallist[1] * dd[1] * id[2] * id[3] +
      vallist[2] * dd[1] * dd[2] * id[3] +
      vallist[3] * id[1] * dd[2] * id[3] +
      vallist[4] * id[1] * id[2] * dd[3] +
      vallist[5] * dd[1] * id[2] * dd[3] +
      vallist[6] * dd[1] * dd[2] * dd[3] +
      vallist[7] * id[1] * dd[2] * dd[3];
  } /* each 3d slab */
  return 1;
} /* int niik_image_interpolate_3d_linear_ijk_update(nifti_image *img,niikpt p,double *v) */

double niik_image_interpolate_float_image_3d_linear(nifti_image *img,niikpt pw)
/* p is xyz coordinates and not ijk */
{
  int
  n,voffset,
  idxlist[9],
  coord[7];
  niikpt p;
  float *fimg;
  double
  dval,
  pos[9],
  dd[7],id[7],
  vallist[9];
  /* check inputs */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(img->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"ERROR: img is not float32\n");
    return NIIKMAX;
  }
  /* get the voxel index while checking for img.dim */
  niik_world_to_ijk(img,&pw,&p);
  pos[1]=p.x;
  pos[2]=p.y;
  pos[3]=p.z;

  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    dd[n] = pos[n];

    coord[n] = floor(dd[n]); /* floor */
    if(coord[n]<-1) return 0; /* -1 is OK */
    if(coord[n]>=img->dim[n]) return 0;
    dd[n] = dd[n] - coord[n];
    id[n] = 1.0 - dd[n];
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* fprintf(stdout,"%3i %3i %3i   %i \n",coord[1],coord[2],coord[3],coord[0]); */
  idxlist[0]=coord[0];
  idxlist[1]=idxlist[0]+1;
  idxlist[2]=idxlist[1]+img->nx;
  idxlist[3]=idxlist[2]-1;
  idxlist[4]=idxlist[0]+img->nx*img->ny;
  idxlist[5]=idxlist[4]+1;
  idxlist[6]=idxlist[5]+img->nx;
  idxlist[7]=idxlist[6]-1;
  if(coord[1]<0)          idxlist[0]=idxlist[3]=idxlist[4]=idxlist[7]=-1;
  if(coord[1]>=img->nx-1) idxlist[1]=idxlist[2]=idxlist[5]=idxlist[6]=-1;
  if(coord[2]<0)          idxlist[0]=idxlist[1]=idxlist[4]=idxlist[5]=-1;
  if(coord[2]>=img->ny-1) idxlist[2]=idxlist[3]=idxlist[6]=idxlist[7]=-1;
  if(coord[3]<0)          idxlist[0]=idxlist[1]=idxlist[2]=idxlist[3]=-1;
  if(coord[3]>=img->nz-1) idxlist[4]=idxlist[5]=idxlist[6]=idxlist[7]=-1;
  if(img->dim[0]<=3) {
    img->nu=img->nt=1;
  }
  fimg=img->data;
  for(n=0; n<8; n++) {
    if(idxlist[n]<0) vallist[n]=0;
    else {
      vallist[n]=fimg[idxlist[n]];
      /*if(niik_check_double_problem(vallist[n])){
      fprintf(stderr,"ERROR: niik_image_get_voxel\n");
      return NIIKMAX; } */
    }
  }
  dval =
    vallist[0] * id[1] * id[2] * id[3] +
    vallist[1] * dd[1] * id[2] * id[3] +
    vallist[2] * dd[1] * dd[2] * id[3] +
    vallist[3] * id[1] * dd[2] * id[3] +
    vallist[4] * id[1] * id[2] * dd[3] +
    vallist[5] * dd[1] * id[2] * dd[3] +
    vallist[6] * dd[1] * dd[2] * dd[3] +
    vallist[7] * id[1] * dd[2] * dd[3];
  /* fprintf(stdout,"out = %f \n",dval); */
  return dval;
} /* double niik_image_interpolate_3d_linear(nifti_image *img,niikpt p) */

double niik_image_interpolate_float_image_3d_linear_ijk(nifti_image *img,niikpt p) {
  int
  n,voffset,
  idxlist[9],
  coord[7];
  float *fimg;
  double
  dval,
  pos[9],
  dd[7],id[7],
  vallist[9];
  /* check inputs */
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return NIIKMAX;
  }
  if(img->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"ERROR: img is not float32\n");
    return NIIKMAX;
  }
  /* get the voxel index while checking for img.dim */
  pos[1]=p.x;
  pos[2]=p.y;
  pos[3]=p.z;
  for(n=1,coord[0]=0,voffset=1; n<=3; n++) {
    coord[n] = floor(pos[n]); /* floor */
    if(coord[n]<-1) return 0; /* -1 is OK */
    if(coord[n]>=img->dim[n]) return 0;
    dd[n] = pos[n] - coord[n];
    id[n] = 1.0 - dd[n];
    coord[0] += coord[n] * voffset;
    voffset *= img->dim[n];
  }
  /* fprintf(stdout,"%3i %3i %3i   %i \n",coord[1],coord[2],coord[3],coord[0]); */
  idxlist[0]=coord[0];
  idxlist[1]=idxlist[0]+1;
  idxlist[2]=idxlist[1]+img->nx;
  idxlist[3]=idxlist[2]-1;
  idxlist[4]=idxlist[0]+img->nx*img->ny;
  idxlist[5]=idxlist[4]+1;
  idxlist[6]=idxlist[5]+img->nx;
  idxlist[7]=idxlist[6]-1;
  if(coord[1]<0)          idxlist[0]=idxlist[3]=idxlist[4]=idxlist[7]=-1;
  if(coord[1]>=img->nx-1) idxlist[1]=idxlist[2]=idxlist[5]=idxlist[6]=-1;
  if(coord[2]<0)          idxlist[0]=idxlist[1]=idxlist[4]=idxlist[5]=-1;
  if(coord[2]>=img->ny-1) idxlist[2]=idxlist[3]=idxlist[6]=idxlist[7]=-1;
  if(coord[3]<0)          idxlist[0]=idxlist[1]=idxlist[2]=idxlist[3]=-1;
  if(coord[3]>=img->nz-1) idxlist[4]=idxlist[5]=idxlist[6]=idxlist[7]=-1;
  if(img->dim[0]<=3) {
    img->nu=img->nt=1;
  }
  fimg=img->data;
  for(n=0; n<8; n++) {
    if(idxlist[n]<0) vallist[n]=0;
    else {
      vallist[n]=fimg[idxlist[n]];
      /*if(niik_check_double_problem(vallist[n])){
      fprintf(stderr,"ERROR: niik_image_get_voxel\n");
      return NIIKMAX; } */
    }
  }
  dval =
    vallist[0] * id[1] * id[2] * id[3] +
    vallist[1] * dd[1] * id[2] * id[3] +
    vallist[2] * dd[1] * dd[2] * id[3] +
    vallist[3] * id[1] * dd[2] * id[3] +
    vallist[4] * id[1] * id[2] * dd[3] +
    vallist[5] * dd[1] * id[2] * dd[3] +
    vallist[6] * dd[1] * dd[2] * dd[3] +
    vallist[7] * id[1] * dd[2] * dd[3];
  /* fprintf(stdout,"out = %f \n",dval); */
  return dval;
} /* double niik_image_interpolate_3d_linear(nifti_image *img,niikpt p) */




/******************************************************************
 * b-spline interpolation
 *
 * -assume intensity=0 outside the image
 * Reference:     Thevenaz, Interpolation Revisited, IEEE TMI, 2000 (19)7: 739-758.
 ******************************************************************/

nifti_image *niik_image_interpolate_convert_3d_bspline_coeff2(nifti_image *img) {
  nifti_image *outimg=NULL;
  const char * fcname="niik_image_interpolate_convert_3d_bspline_coeff2";
  if((outimg=niik_image_copy(img))==NULL) {
    fprintf(stderr,"[%s] ERROR: niik_image_copy\n",fcname);
    return NULL;
  }
  if(!niik_image_interpolate_convert_3d_bspline_coeff(outimg)) {
    fprintf(stderr,"[%s] ERROR: niik_image_interpolate_convert_3d_bspline_coeff\n",fcname);
    return NULL;
  }
  return outimg;
} /* niik_image_interpolate_convert_3d_bspline_coeff2 */


int niik_image_interpolate_convert_3d_bspline_coeff(nifti_image *img)
/* creates b-spline's coefficient image
 * -more than 3d is OK but coefficients are calculated by 3D basis
 */
{
  char fcname[129]="niik_image_interpolate_convert_3d_bspline_coeff";
  int num,i,j,k,u,t,v,w,ut,n,dim12,size;
  double *vv,*imgdata;
  niikmat *A;
  int verbose=0;

  if(verbose>=1)  fprintf(stdout,"[%s] start\n",fcname);
  if(img==NULL) {
    fprintf(stderr,"[%s] ERROR: img is a null pointer\n",fcname);
    return 0;
  }

  if(verbose>=1) fprintf(stdout,"[%s] convert to double\n",fcname);
  if(!niik_image_type_convert(img,NIFTI_TYPE_FLOAT64)) {
    fprintf(stderr,"[%s] ERROR: niik_image_type_convert\n",fcname);
    return 0;
  }

  imgdata=img->data;
  dim12 = img->nx*img->ny;
  size  = dim12 * img->nz;
  if(verbose>=1) fprintf(stdout,"[%s] img = %i %i %i %i %i\n",fcname,
                           img->nx,img->ny,img->nz,img->nt,img->nu);

  num=NIIK_IMAX(img->nx,NIIK_IMAX(img->ny,img->nz));
  vv=(double *)malloc(num*sizeof(double));

  if(verbose>=1) fprintf(stdout,"[%s] v size = %i\n",fcname,num);
  /*fprintf(stdout,"  %i %i %i %i %i %i %i \n",
    img->dim[1],img->dim[2],img->dim[3],img->dim[4],img->dim[5],img->dim[6],img->dim[7]);*/

  if(verbose>=1) fprintf(stdout,"[%s] main loop\n",fcname);
  for(w=ut=0; w<img->dim[7]; w++) {
    for(v=0; v<img->dim[6]; v++) {
      for(u=0; u<img->dim[5]; u++) {
        for(t=0; t<img->dim[4]; ut++,t++) {
          if(verbose>=2) fprintf(stdout,"[%s] ut loop %i\n",fcname,ut);

          /* along x */
          /* prepare A for x-dir */
          if(verbose>=1) fprintf(stderr,"[%s]  x-dir\n",fcname);
          A=niikmat_init(img->nx,img->nx);
          if(!niik_bspline_update_A(A)) {
            fprintf(stderr,"[%s] ERROR: niik_bspline_update for x-dir\n",fcname);
            return 0;
          }
          if(verbose>=1) fprintf(stderr,"[%s]  x-dir %i %i *\n",fcname,(int)A->row,(int)A->col);
          n=0;
          for(k=0; k<img->nz; k++) {
            for(j=0; j<img->ny; j++) {
              for(i=0; i<img->nx; n++,i++) {
                vv[i]=imgdata[n + ut*size];
              }
              if(verbose>=3) fprintf(stderr,"[%s]     yz=%i %i | %i %i\n",fcname,j,k,img->ny,img->nz);

              if(!niik_bspline_calc_coeff_update(A,vv,img->nx)) {
                fprintf(stderr,"[%s] ERROR: niik_bspline_calc_coeff_update\n",fcname);
                return 0;
              }
              if(verbose>=2) fprintf(stderr,"[%s]     yz=%i %i | %i %i\n",fcname,j,k,img->ny,img->nz);

              for(i=0,n-=img->nx; i<img->nx; n++,i++) {
                imgdata[n+ut*size]=vv[i];
              }
            }
          }
          if(verbose>1) fprintf(stderr,"\t\tx-dir free A %i %i\n",(int)A->row,(int)A->col);
          niikmat_free(A);

          /* along y */
          /* prepare A for y-dir */
          if(verbose) fprintf(stderr,"\t\ty-dir\n");
          A=niikmat_init(img->ny,img->ny);
          if(verbose) fprintf(stderr,"\t\ty-dir bspline update %i\n",img->ny);
          if(!niik_bspline_update_A(A)) {
            fprintf(stderr,"ERROR: niik_bspline_update for y-dir\n");
            return 0;
          }
          if(verbose) fprintf(stderr,"\t\ty-dir *\n");
          for(k=0; k<img->nz; k++) {
            for(i=0; i<img->nx; i++) {
              for(j=0,n=k*dim12+i; j<img->ny; n+=img->nx,j++) {
                vv[j]=imgdata[n+ut*size];
              }
              if(!niik_bspline_calc_coeff_update(A,vv,img->ny))\
              {
                fprintf(stderr,"ERROR: niik_bspline_calc_coeff_update\n");
                return 0;
              }

              for(j=0,n=k*dim12+i; j<img->ny; n+=img->nx,j++) {
                imgdata[n+ut*size]=vv[j];
              }
            }
          }
          niikmat_free(A);

          /* along z */
          /* prepare A for z-dir */
          if(verbose) fprintf(stderr,"\t\tz-dir\n");
          A=niikmat_init(img->nz,img->nz);
          if(!niik_bspline_update_A(A)) {
            fprintf(stderr,"ERROR: niik_bspline_update for x-dir\n");
            return 0;
          }
          if(verbose) fprintf(stderr,"\t\tz-dir *\n");
          for(i=0; i<img->nx; i++) {
            for(j=0; j<img->ny; j++) {
              for(k=0,n=i+j*img->nx; k<img->nz; n+=dim12,k++) {
                vv[k]=imgdata[n+ut*size];
              }
              if(!niik_bspline_calc_coeff_update(A,vv,img->nz)) {
                fprintf(stderr,"ERROR: niik_bspline_calc_coeff_update\n");
                return 0;
              }
              for(k=0,n=i+j*img->nx; k<img->nz; n+=dim12,k++) {
                imgdata[n+ut*size]=vv[k];
              }
            }
          }
          niikmat_free(A);
        }
      }
    }
  } /* u and t */
  free(vv);

  if(verbose) fprintf(stdout,"-v (niik_image_interpolate_convert_3d_bspline_coeff) finished\n");
  return 1;
} /* int niik_image_interpolate_convert_3d_bspline_coeff */

int niik_image_interpolate_inverse_3d_bspline_coeff(nifti_image *coeffimg, nifti_image *img)
/* -img is replaced */
{
  int i,j,k,m,nn,n,s,no=1;
  niikpt p;
  niikmat *bt;
  if(coeffimg==NULL) {
    fprintf(stderr,"ERROR: coeffimg is null\n");
    return 0;
  }
  if(img==NULL)     {
    fprintf(stderr,"ERROR: img is null\n");
    return 0;
  }
  /* img must have same nt, nu, nv, nw */
  s=img->nx*img->ny*img->nz;

  if(coeffimg->ndim > img->ndim) {
    img->ndim = img->dim[0] = coeffimg->ndim;
    img->nt   = img->dim[4] = coeffimg->nt;
    img->nu   = img->dim[5] = coeffimg->nu;
    img->nv   = img->dim[6] = coeffimg->nv;
    img->nw   = img->dim[7] = coeffimg->nw;
    img->dt   = img->pixdim[4] = coeffimg->dt;
    img->du   = img->pixdim[5] = coeffimg->du;
    img->dv   = img->pixdim[6] = coeffimg->dv;
    img->dw   = img->pixdim[7] = coeffimg->dw;

    for(n=img->nvox=1; n<=img->ndim; n++) {
      img->nvox*=img->dim[n];
    }
    /*fprintf(stdout,"  %i   %3i %3i %3i   %3i %3i %3i \n",img->ndim,img->nx,img->ny,img->nz,img->nt,img->nu,img->nv);
      fprintf(stdout,"  nvox = %i -> %i\n",img->nvox,img->nvox*img->nbyper);*/
    free(img->data);
    img->data=(void *)calloc(img->nvox,img->nbyper);
  }

  bt=niikmat_init(img->nz,250);
  p.w=0;
  /*#pragma omp parallel for private(i,j,m,n,p,nn) reduction(+:no)*/
  for(k=0; k<img->nz; k++) {

    n=k*img->nx*img->ny;
    for(j=0; j<img->ny; j++) {

      for(i=0; i<img->nx; n++,i++) {

        /*TODO: check with Kunio if this is correct*/
        niik_index_to_world(img,i,j,k,&p);

        if(!niik_image_interpolate_3d_bspline_xyz_update(coeffimg, p, bt->m[k])) {
          fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_xyz_update\n");
          no=0;
        }
        for(m=n,nn=0; m<img->nvox; m+=s,nn++) {
          if(!niik_image_set_voxel(img,m,bt->m[k][nn])) {
            fprintf(stderr,"ERROR: niik_image_set_voxel\n");
            no=0;

          }
        }
      }
    }
  }
  bt=niikmat_free(bt);
  return no;
}

double niik_image_interpolate_3d_bspline(nifti_image *bsplimg, niikpt pw) {
  niikpt p;
  niik_world_to_ijk(bsplimg,&pw,&p);
  return niik_image_interpolate_3d_bspline_ijk(bsplimg,p);
}

double niik_image_interpolate_3d_bspline_ijk(nifti_image *bsplimg,niikpt p)
/* interpolating a 3D image with b-spline */
{
  double outval;
  const char * fcname="niik_image_interpolate_3d_bspline_ijk";
  if(bsplimg==NULL) {
    fprintf(stderr,"[%s] ERROR: bsplimg is a null pointer\n",fcname);
    return 0;
  }
  if(bsplimg->ndim>3) {
    fprintf(stderr,"[%s] ERROR: bsplimg is not 3d\n",fcname);
    fprintf(stderr,"  please use 'niik_image_interpolate_3d_bspline_xyz_update'\n");
    return 0;
  }
  if(!niik_image_interpolate_3d_bspline_ijk_update(bsplimg,p,&outval)) {
    fprintf(stderr,"[%s] ERROR: niik_image_interpolate_3d_bspline_xyz_update\n",fcname);
    return 0;
  }
  return outval;
} /* double niik_image_interpolate_3d_bspline_ijk(nifti_image *bsplimg,niikpt p);  */

int niik_image_interpolate_3d_bspline_xyz_update(nifti_image *bsplimg,niikpt pw,double *v) {
  niikpt p;
  niik_world_to_ijk(bsplimg,&pw,&p);

  if(!niik_image_interpolate_3d_bspline_ijk_update(bsplimg,p,v)) {
    fprintf(stderr,"ERROR: niik_image_interpolate_3d_bspline_ijk_update(bsplimg,p,v)\n");
    return 0;
  }
  return 1;
}

int niik_image_interpolate_3d_bspline_ijk_update(nifti_image *bsplimg,niikpt p,double *v)
/* -wrapper to decide float32 or float64
 * -2012-10-27, Kunio */
{
  switch(bsplimg->datatype) {
  case NIFTI_TYPE_FLOAT64:
    if(!niik_image_interpolate_3d_bspline_ijk_update_double_type(bsplimg,p,v)) {
      fprintf(stderr,"[niik_image_interpolate_3d_bspline_ijk_update] ERROR: niik_image_interpolate_3d_bspline_ijk_update_double_type\n");
      return 0;
    }
    break;
  case NIFTI_TYPE_FLOAT32:
    if(!niik_image_interpolate_3d_bspline_ijk_update_float_type(bsplimg,p,v)) {
      fprintf(stderr,"[niik_image_interpolate_3d_bspline_ijk_update] ERROR: niik_image_interpolate_3d_bspline_ijk_update_float_type\n");
      return 0;
    }
    break;
  default:
    fprintf(stderr,"[niik_image_interpolate_3d_bspline_ijk_update] ERROR: image type is not float32 or float64, %s\n",
            nifti_datatype_string(bsplimg->datatype));
    return 0;
  }
  return 1;
} /* int niik_image_interpolate_3d_bspline_ijk_update(nifti_image *bsplimg,niikpt p,double *v) */

int niik_image_interpolate_3d_bspline_ijk_update_double_type(nifti_image *bsplimg,niikpt p,double *v)
/* -typically called from niik_image_interpolate_3d_bspline_ijk_update
 * -this can be internal function
 */
{
  const char * fcname="niik_image_interpolate_3d_bspline_ijk_upate_double_type";
  int
  area,size,i,j,k,ii[4],ei[4],m[4],n,ni,nj,nk,bb,b;
  double
  *imgdata,
  pos[4],
  Weight[4][4],
  w,w2,t0,t1;
  const int degree=3, md=1;
  if(bsplimg==NULL) {
    fprintf(stderr,"ERROR: bsplimg is a null pointer\n");
    return 0;
  }
  if(bsplimg->datatype!=NIFTI_TYPE_FLOAT64) {
    fprintf(stderr,"[%s] ERROR: bspline image must be double\n",fcname);
    return 0;
  }
  pos[1]=p.x;
  pos[2]=p.y;
  pos[3]=p.z;
  /* check for voxel position and calculate voxel space position */
  for(n=1; n<=3; n++) {
    if(pos[n]<-2.5)  {
      bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz; /*bsplimg->nt*bsplimg->nu*bsplimg->nv*bsplimg->nw;*/
      for(b=0; b<bb; b++) v[b]=0;
      return 1;
    } else if(pos[n]>(1.5+bsplimg->dim[n]) )  {
      bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz; /*bsplimg->nt*bsplimg->nu*bsplimg->nv*bsplimg->nw;*/
      for(b=0; b<bb; b++) v[b]=0;
      return 1;
    }
  }
  for(n=1; n<=3; n++) {
    m[n]=(int)floor(pos[n]);
    ii[n]=m[n]-md;
    ei[n]=ii[n]+degree;
    w=pos[n]-m[n];
    w2=w*w;
    Weight[n][3] = w2 * w / 6.0;
    Weight[n][0] = (1.0 / 6.0) + (w2-w) / 2.0 - Weight[n][3];
    Weight[n][2] = w + Weight[n][0] - 2.0 * Weight[n][3];
    Weight[n][1] = 1.0 - Weight[n][0] - Weight[n][2] - Weight[n][3];
    ei[n]=NIIK_IMIN(ei[n],bsplimg->dim[n]-1);
  }
  area = bsplimg->nx*bsplimg->ny;
  size = area*bsplimg->nz;
  /* perform interpolation */
  /*bb=bsplimg->nt*bsplimg->nu*bsplimg->nv*bsplimg->nw;*/
  bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz;
  bb=(bb<=0)?1:bb;
  for(b=0,imgdata=bsplimg->data; b<bb; b++,imgdata+=size) {
    v[b]=0;
    for (nk=0,k = ii[3]; k <= ei[3]; nk++,k++) {
      if(k<0) continue;
      t0=0;
      for (nj=0,j = ii[2]; j <= ei[2]; nj++,j++) {
        if(j<0) continue;
        t1=0;
        n = ii[1] + j*bsplimg->nx + k*area;
        for (ni=0,i = ii[1]; i <= ei[1]; n++,ni++,i++) {
          if(i<0) continue;
          t1 += Weight[1][ni] * imgdata[n];
        }
        t0 += Weight[2][nj] * t1;
      }
      v[b]+=Weight[3][nk] * t0;
    }
  } /* each 3d slab */
  return 1;
} /* int niik_image_interpolate_3d_bspline_ijk_update_double_type(nifti_image *bsplimg,niikpt p,double *v) */

int niik_image_interpolate_3d_bspline_ijk_update_float_type(nifti_image *bsplimg,niikpt p,double *v)
/* -typically called from niik_image_interpolate_3d_bspline_ijk_update
 * -this can be internal function
 * -the only differnece from niik_image_interpolate_3d_bspline_ijk_update_double_type
 *  is that imgdata is <float *>
 */
{
  const char * fcname="niik_image_interpolate_3d_bspline_ijk_upate_float_type";
  int
  area,size,i,j,k,ii[3],ei[3],m[3],n,ni,nj,nk,bb,b;
  float
  *imgdata,
  pos[3],
  Weight[3][4],
  w,w2,t0,t1;
  const int degree=3, md=1;
  if(bsplimg==NULL) {
    fprintf(stderr,"[%s] ERROR: bsplimg is a null pointer\n",fcname);
    return 0;
  }
  if(bsplimg->datatype!=NIFTI_TYPE_FLOAT32) {
    fprintf(stderr,"[%s] ERROR: bspline image must be float\n",fcname);
    return 0;
  }
  pos[0]=p.x;
  pos[1]=p.y;
  pos[2]=p.z;
  /* check for voxel position and calculate voxel space position */
  for(n=0; n<3; n++) {
    if(pos[n]<-2.5)  {
      bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz;
      for(b=0; b<bb; b++) v[b]=0;
      return 1;
    } else if(pos[n]>(1.5+bsplimg->dim[n+1]) )  {
      bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz;
      for(b=0; b<bb; b++) v[b]=0;
      return 1;
    }
  }
  for(n=0; n<3; n++) {
    m[n]=(int)floor(pos[n]);
    ii[n]=m[n]-md;
    w=pos[n]-m[n];
    w2=w*w;
    Weight[n][3] = w2 * w / 6.0;
    Weight[n][0] = (1.0 / 6.0) + (w2-w) / 2.0 - Weight[n][3];
    Weight[n][2] = w + Weight[n][0] - 2.0 * Weight[n][3];
    Weight[n][1] = 1.0 - Weight[n][0] - Weight[n][2] - Weight[n][3];
    ei[n]=NIIK_IMIN(ii[n]+degree,bsplimg->dim[n+1]-1);
  }
  area = bsplimg->nx*bsplimg->ny;
  size = area*bsplimg->nz;
  /* perform interpolation */
  bb=bsplimg->nvox/bsplimg->nx/bsplimg->ny/bsplimg->nz;
  bb=(bb<=0)?1:bb;
  for(b=0,imgdata=bsplimg->data; b<bb; b++,imgdata+=size) {
    v[b]=0;
    for (nk=0,k = ii[2]; k <= ei[2]; nk++,k++) {
      if(k<0) continue;
      t0=0;
      for (nj=0,j = ii[1]; j <= ei[1]; nj++,j++) {
        if(j<0) continue;
        t1=0;
        n = k*area + j*bsplimg->nx + ii[0];
        for (ni=0,i = ii[0]; i <= ei[0]; ni++,i++) {
          if(i<0) continue;
          t1 += Weight[0][ni] * imgdata[n+ni];
        }
        t0 += Weight[1][nj] * t1;
      }
      v[b]+=Weight[2][nk] * t0;
    }
  } /* each 3d slab */
  return 1;
} /* int niik_image_interpolate_3d_bspline_ijk_update_float_type */


int niik_image_interpolate_cost(nifti_image *img,nifti_image *refimg,niikmat *afmat,double *cost)
/* interpolation cost (sum of distance from the reference grid) */
{
  niikpt pt,qt;
  int i,j,k;
  const int verbose=0;
  double tmpout;
  if(verbose) {
    fprintf(stdout,"[niik_image_interpolate_cost] start\n");
  }
  if(img==NULL)   {
    fprintf(stderr,"ERROR: img is a null pointer\n");
    return 0;
  }
  if(refimg==NULL) {
    fprintf(stderr,"ERROR: refimg is a null pointer\n");
    return 0;
  }
  if(afmat==NULL) {
    fprintf(stderr,"ERROR: afmat is a null pointer\n");
    return 0;
  }

  for(k=0,tmpout=0; k<refimg->nz; k++) {

    /*pt.z = k * refimg->dz;*/

    for(j=0; j<refimg->ny; j++) {

      /*pt.y = j * refimg->dy;*/

      for(i=0; i<refimg->nx; i++) {

        /*pt.x = i * refimg->dx;*/
        niik_index_to_world(refimg,i,j,k,&pt);

        qt = niikpt_affine_transform(afmat,pt);

        /*qt.x /= img->dx;
        qt.y /= img->dy;
        qt.z /= img->dz;*/
        niik_world_to_ijk(img,&pt,&qt);

        qt.x = qt.x - floor(qt.x);
        qt.y = qt.y - floor(qt.y);
        qt.z = qt.z - floor(qt.z);

        tmpout += niikpt_mag(qt);
      }
    }
  }
  *cost=tmpout;
  if(verbose) {
    fprintf(stdout,"[niik_image_interpolate_cost]   cost = %12.4f\n",tmpout);
  }
  if(verbose) {
    fprintf(stdout,"[niik_image_interpolate_cost] successful exit\n");
  }
  return 1;
} /* niik_image_interpolate_cost */



/************************************************************
 *
 * niik_image_interp_along_normal
 *
 *
 ************************************************************/

int niik_image_interp_along_normal_double_vector(nifti_image *img,int interp,niikpt pt,niikpt normal,double *x,double *out,int num)
/* img is the image to be interpolated
 * interp is the type of interpolation
 * pt is the point of origin
 * normal is the normal from the origin
 * x is the vector of distances along the normal
 * out is the output vector of image intensities
 * num is the size of x- and out-vectors */
{
  int n;
  niikpt qt;

  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n" );
    return 0;
  }
  if(x==NULL)   {
    fprintf(stderr,"ERROR: x is null\n" );
    return 0;
  }
  if(out==NULL) {
    fprintf(stderr,"ERROR: out is null\n" );
    return 0;
  }
  /*TODO:double check coord conversion*/
  for(n=0; n<num; n++) {
    qt=niikpt_move_normal(pt,normal,x[n]);
    if(!niik_image_interpolate_3d_xyz_update(img,qt,interp,&out[n])) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_xyz_update\n");
      return 0;
    }
  }
  return 1;
} /* niik_image_interp_along_normal */

int niik_image_interp_along_normal(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v,niikvec *out)
/* img is the image to be interpolated
 * interp is the type of interpolation
 * pt is the point of origin
 * normal is the normal from the origin
 * v is the vector of distances along the normal
 * out is the output vector of image intensities */
{
  int n;
  niikpt qt;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n" );
    return 0;
  }
  if(v==NULL)   {
    fprintf(stderr,"ERROR: v is null\n" );
    return 0;
  }
  if(out==NULL) {
    fprintf(stderr,"ERROR: out is null\n" );
    return 0;
  }
  if(out->num!=v->num) {
    fprintf(stderr,"ERROR: vector lengths are different\n" );
    return 0;
  }
  /*#pragma omp parallel for*/
  for(n=0; n<v->num; n++) {
    qt=niikpt_move_normal(pt, normal, v->v[n]);

    if(!niik_image_interpolate_3d_xyz_update(img,qt,interp,&out->v[n])) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_xyz_update\n");
      continue; /*return 0;*/
    }
  }
  return 1;
} /* niik_image_interp_along_normal */

int niik_image_interp_along_normal_update(nifti_image *img,int interp,niikpt pt,niikpt normal,niikvec *v) {
  int n;
  niikpt qt;
  if(img==NULL) {
    fprintf(stderr,"ERROR: img is null\n" );
    return 0;
  }
  if(v==NULL)   {
    fprintf(stderr,"ERROR: v is null\n" );
    return 0;
  }
  /*#pragma omp parallel for*/
  for(n=0; n<v->num; n++) {
    qt=niikpt_move_normal(pt,normal,v->v[n]);
    if(!niik_image_interpolate_3d_xyz_update(img,qt,interp,&v->v[n])) {
      fprintf(stderr,"ERROR: niik_image_interpolate_3d_xyz_update\n");
      continue; /*return 0;*/
    }
  }
  return 1;
} /* niik_image_interp_along_normal_update */


#endif /* _FALCON_INTERP_C_ */

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8
*/