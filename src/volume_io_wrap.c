/*
# Developed by Simon Fristed Eskildsen, eskild@gmail.com
# Part of FALCON
#
# Copyright notice:
# This code is copyright Simon Fristed Eskildsen.
# It may not be copied, altered in any way or transmitted
# to others (unless explicitly stated otherwise) without
# the written permission of the author/developer. 
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "basic.h"
#include <volume_io.h>
#include <limits.h>
#include "volume_io_wrap.h"
#include <sys/types.h>
#include <signal.h>
#include "array_alloc.h"
#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

//#define DEBUG 1

int volumeToWrap(VIO_Volume *volume, Volume_wrap *wrap) {

  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;
  int i,j,k;

  get_volume_sizes(*volume,wrap->sizes);

  wrap->type = get_volume_nc_data_type(*volume,&wrap->sign);
  
  switch(wrap->type){
  case NC_BYTE:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_BYTE\n");
#endif
    wrap->type_size=NC_BYTE_SIZE;
    break;
  case NC_CHAR:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_CHAR\n");
#endif
    wrap->type_size=NC_CHAR_SIZE;
    break;
  case NC_NAT:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_NAT\n");
#endif
    wrap->type_size=NC_NAT_SIZE;
    break;
  case NC_SHORT:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_SHORT\n");
#endif
    wrap->type_size=NC_SHORT_SIZE;
    break;
  case NC_INT:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_INT\n");
#endif
    wrap->type_size=NC_INT_SIZE;
    break;
  case NC_FLOAT:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_FLOAT\n");
#endif
    wrap->type_size=NC_FLOAT_SIZE;
    break;
  case NC_DOUBLE:
#ifdef DEBUG
    fprintf(stderr,"Data type: NC_DOUBLE\n");
#endif
    wrap->type_size=NC_DOUBLE_SIZE;
    break;
  } 

  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = malloc(wrap->sizes[0]*sizeof(*bdata));
    bdata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**bdata));
    for(i=1;i<wrap->sizes[0];i++)
      bdata[i] = bdata[0] + i*wrap->sizes[1];
    
    bdata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***bdata));

    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	bdata[i][j] = bdata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];

    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  bdata[i][j][k] = get_volume_voxel_value(*volume,i,j,k,0,0);

    wrap->data = (void *)bdata;
    break;
  case NC_SHORT:
    sdata = malloc(wrap->sizes[0]*sizeof(*sdata));
    sdata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**sdata));
    for(i=1;i<wrap->sizes[0];i++)
      sdata[i] = sdata[0] + i*wrap->sizes[1];
    
    sdata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***sdata));
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	sdata[i][j] = sdata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  sdata[i][j][k] = get_volume_voxel_value(*volume,i,j,k,0,0);
    wrap->data = (void *)sdata;
    break;
  case NC_INT:
    idata = malloc(wrap->sizes[0]*sizeof(*idata));
    idata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**idata));
    for(i=1;i<wrap->sizes[0];i++)
      idata[i] = idata[0] + i*wrap->sizes[1];
    
    idata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***idata));
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
      idata[i][j] = idata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  idata[i][j][k] = get_volume_voxel_value(*volume,i,j,k,0,0);
    wrap->data = (void *)idata;
    break;
  case NC_FLOAT:
    fdata = malloc(wrap->sizes[0]*sizeof(*fdata));
    fdata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**fdata));
    for(i=1;i<wrap->sizes[0];i++)
      fdata[i] = fdata[0] + i*wrap->sizes[1];
    
    fdata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***fdata));
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	fdata[i][j] = fdata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  fdata[i][j][k] = get_volume_voxel_value(*volume,i,j,k,0,0);
    wrap->data = (void *)fdata;
    break;
  case NC_DOUBLE:
    ddata = malloc(wrap->sizes[0]*sizeof(*ddata));
    ddata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**ddata));
    for(i=1;i<wrap->sizes[0];i++)
      ddata[i] = ddata[0] + i*wrap->sizes[1];
    
    ddata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***ddata));
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
      ddata[i][j] = ddata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  ddata[i][j][k] = get_volume_voxel_value(*volume,i,j,k,0,0);
    wrap->data = (void *)ddata;
    break;
  }


  return STATUS_OK;

}

int volumeToWrap_real(VIO_Volume *volume, Volume_wrap *wrap) {
  double ***ddata;
  int i,j,k;

  get_volume_sizes(*volume,wrap->sizes);

  wrap->type = NC_DOUBLE;
  wrap->type_size=NC_DOUBLE_SIZE;

#ifdef DEBUG
    fprintf(stderr,"Wrapping data...\n");
    fprintf(stderr,"\tvolume data type: %d\n",get_volume_nc_data_type(*volume,&wrap->sign));
    fprintf(stderr,"\twrap data type: %d\n",wrap->type);
    fprintf(stderr,"\tsizes: %d %d %d\n",wrap->sizes[0],wrap->sizes[1],wrap->sizes[2]);
#endif
  
  ddata = malloc(wrap->sizes[0]*sizeof(*ddata));
  ddata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**ddata));
  for(i=1;i<wrap->sizes[0];i++)
    ddata[i] = ddata[0] + i*wrap->sizes[1];
  
  ddata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***ddata));
  for(i=0;i<wrap->sizes[0];i++)
    for(j=0;j<wrap->sizes[1];j++)
      ddata[i][j] = ddata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];
  for(i=0;i<wrap->sizes[0];i++)
    for(j=0;j<wrap->sizes[1];j++)
      for(k=0;k<wrap->sizes[2];k++)
	ddata[i][j][k] = get_volume_real_value(*volume,i,j,k,0,0);
  wrap->data = (void *)ddata;

  return STATUS_OK;

}

int wrapDataToVolume(VIO_Volume *volume, Volume_wrap *wrap){

  int i,j,k,sizes[3];
  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;

#ifdef DEBUG
    fprintf(stderr,"Unwrapping data...\n");
    fprintf(stderr,"\tvolume data type: %d\n",get_volume_nc_data_type(*volume,&wrap->sign));
    fprintf(stderr,"\twrap data type: %d\n",wrap->type);
    get_volume_sizes(*volume,sizes);
    fprintf(stderr,"\tvolume sizes: %d %d %d\n",sizes[0],sizes[1],sizes[2]);
    fprintf(stderr,"\twrap sizes: %d %d %d\n",wrap->sizes[0],wrap->sizes[1],wrap->sizes[2]);
#endif

  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
#ifdef DEBUG
    fprintf(stderr,"Data type: BYTE\n");
#endif
    bdata = (byte ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,bdata[i][j][k]);
    break;
  case NC_SHORT:
#ifdef DEBUG
    fprintf(stderr,"Data type: SHORT\n");
#endif
    sdata = (short ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++) {
	  set_volume_real_value(*volume,i,j,k,0,0,sdata[i][j][k]);
	}
#ifdef DEBUG
    fprintf(stderr,"Done assigning\n");
#endif    
    break;
  case NC_INT:
#ifdef DEBUG
    fprintf(stderr,"Data type: INT\n");
#endif
    idata = (int ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,idata[i][j][k]);
    break;
  case NC_FLOAT:
#ifdef DEBUG
    fprintf(stderr,"Data type: FLOAT\n");
#endif
    fdata = (float ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,fdata[i][j][k]);
    break;
  case NC_DOUBLE:
#ifdef DEBUG
    fprintf(stderr,"Data type: DOUBLE\n");
#endif
    ddata = (double ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,ddata[i][j][k]);
    break;
  }

  return STATUS_OK;

}

/* copy contents of a Volume_wrap to a VIO_Volume */
int wrapVoxelDataToVolume(VIO_Volume *volume, Volume_wrap *wrap){

  int i,j,k;
  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;

  if ((wrap->sizes[0] == 0) || (wrap->sizes[1] == 0) || (wrap->sizes[2] == 0)){
    fprintf(stderr,"ERROR! Trying to copy data from empty wrapped volume.\n");
    return STATUS_ERR;
  }

  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = (byte ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_voxel_value(*volume,i,j,k,0,0,bdata[i][j][k]);
    break;
  case NC_SHORT:
    sdata = (short ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++) {
	  set_volume_voxel_value(*volume,i,j,k,0,0,sdata[i][j][k]);
	}
    break;
  case NC_INT:
    idata = (int ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_voxel_value(*volume,i,j,k,0,0,idata[i][j][k]);
    break;
  case NC_FLOAT:
    fdata = (float ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_voxel_value(*volume,i,j,k,0,0,fdata[i][j][k]);
    break;
  case NC_DOUBLE:
    ddata = (double ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_voxel_value(*volume,i,j,k,0,0,ddata[i][j][k]);
    break;
  }

  return STATUS_OK;

}

/* copy contents of a Volume_wrap to a VIO_Volume */
int wrapVoxelDataToVolume_real(VIO_Volume *volume, Volume_wrap *wrap){

  int i,j,k;
  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;

  if ((wrap->sizes[0] == 0) || (wrap->sizes[1] == 0) || (wrap->sizes[2] == 0)){
    fprintf(stderr,"ERROR! Trying to copy data from empty wrapped volume.\n");
    return STATUS_ERR;
  }

  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = (byte ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,bdata[i][j][k]);
    break;
  case NC_SHORT:
    sdata = (short ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++) {
	  set_volume_real_value(*volume,i,j,k,0,0,sdata[i][j][k]);
	}
    break;
  case NC_INT:
    idata = (int ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,idata[i][j][k]);
    break;
  case NC_FLOAT:
    fdata = (float ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,fdata[i][j][k]);
    break;
  case NC_DOUBLE:
    ddata = (double ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  set_volume_real_value(*volume,i,j,k,0,0,ddata[i][j][k]);
    break;
  }

  return STATUS_OK;

}

void free_wrap(Volume_wrap *wrap) {

  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;

  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = (byte ***)wrap->data;
    free(bdata[0][0]);
    free(bdata[0]);
    free(bdata);
    break;
  case NC_SHORT:
    sdata = (short ***)wrap->data;
    free(sdata[0][0]);
    free(sdata[0]);
    free(sdata);
    break;
  case NC_INT:
    idata = (int ***)wrap->data;
    free(idata[0][0]);
    free(idata[0]);
    free(idata);
    break;
  case NC_FLOAT:
    fdata = (float ***)wrap->data;
    free(fdata[0][0]);
    free(fdata[0]);
    free(fdata);
    break;
  case NC_DOUBLE:
    ddata = (double ***)wrap->data;
    free(ddata[0][0]);
    free(ddata[0]);
    free(ddata);
    break;
  default:
    fprintf(stderr,"ERROR! free_wrap() is freeing nothing!");
  }
  
  wrap->sizes[0] = 0;
  wrap->sizes[1] = 0;
  wrap->sizes[2] = 0;
}

/* void *** alloc_data3D_old(int sizes[3],byte size_element) { */

/*   int i,j; */
/*   void ***data; */

/*   data = (void ***)malloc(sizes[0]*sizeof(void **)); */
/*   data[0] = (void **)malloc(sizes[0]*sizes[1]*sizeof(void *)); */
/*   for(i=1;i<sizes[0];i++) */
/*     data[i] = data[0] + i*sizes[1]; */
    
/*   data[0][0] = (void *)malloc(sizes[0]*sizes[1]*sizes[2]*size_element); */
/*   for(i=0;i<sizes[0];i++) */
/*     for(j=0;j<sizes[1];j++) */
/*       data[i][j] = data[0][0] +i*sizes[1]*sizes[2]*size_element + j*sizes[2]*size_element; */
  
/*   return data; */

/* } */

void ***alloc_data3D(int sizes[3], byte size_element)
{
  void ***iii, **ii, *i;
  int j,limit;
    
  iii = (void ***) malloc(sizes[0]*sizeof(void **));
  alloc_error_check(iii);
  ii = (void **) malloc(sizes[0] * sizes[1] * sizeof(void *));
  alloc_error_check(ii);
  i = (void *) malloc(sizes[0] * sizes[1] * sizes[2]*size_element);
  alloc_error_check(i);

  iii[0] = ii;
  for (j = 1; j < sizes[0]; j++) {
    iii[j] = iii[0] + j*sizes[1];
  }

  limit = sizes[0] * sizes[1];

  ii[0] = i;
  for (j = 1; j < limit; j++) {
    ii[j] = ii[0] + j*sizes[2]*size_element;
  }

  return iii;
}


void copy_wrap(Volume_wrap *original, Volume_wrap *copy) {

  void ***data, ***data_o;

  //data_o = (byte ***)original->data;
  data_o = original->data;

  copy->type = original->type;
  copy->type_size = original->type_size;
  copy->sizes[0] = original->sizes[0];  copy->sizes[1] = original->sizes[1];  copy->sizes[2] = original->sizes[2];
  copy->sign = original->sign;
  
  data = alloc_data3D(copy->sizes,copy->type_size);

  memcpy(data[0][0], data_o[0][0], copy->sizes[0]*copy->sizes[1]*copy->sizes[2]*copy->type_size);

  //copy->data = (void *)data;
  copy->data = data;

}

VIO_Real get_wrap_min(Volume_wrap *wrap){
  int i,j,k;
  VIO_Real min=FLT_MAX;
  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;


  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = (byte ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (bdata[i][j][k]<min)
	    min = bdata[i][j][k];
    break;
  case NC_SHORT:
    sdata = (short ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++) {
	  if (sdata[i][j][k]<min)
	    min = sdata[i][j][k];
	}
    break;
  case NC_INT:
    idata = (int ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (idata[i][j][k]<min)
	    min = idata[i][j][k];
    break;
  case NC_FLOAT:
    fdata = (float ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (fdata[i][j][k]<min)
	    min = fdata[i][j][k];
    break;
  case NC_DOUBLE:
    ddata = (double ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (ddata[i][j][k]<min)
	    min = ddata[i][j][k];
    break;
  }

  return min;
}

VIO_Real get_wrap_max(Volume_wrap *wrap){
  int i,j,k;
  VIO_Real max=0;
  byte ***bdata;
  short ***sdata;
  int ***idata;
  float ***fdata;
  double ***ddata;


  switch(wrap->type) {
  case NC_CHAR:
  case NC_BYTE:
  case NC_NAT:
    bdata = (byte ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (bdata[i][j][k]>max)
	    max = bdata[i][j][k];
    break;
  case NC_SHORT:
    sdata = (short ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++) {
	  if (sdata[i][j][k]>max)
	    max = sdata[i][j][k];
	}
    break;
  case NC_INT:
    idata = (int ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (idata[i][j][k]>max)
	    max = idata[i][j][k];
    break;
  case NC_FLOAT:
    fdata = (float ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (fdata[i][j][k]>max)
	    max = fdata[i][j][k];
    break;
  case NC_DOUBLE:
    ddata = (double ***)wrap->data;
    for(i=0;i<wrap->sizes[0];i++)
      for(j=0;j<wrap->sizes[1];j++)
	for(k=0;k<wrap->sizes[2];k++)
	  if (ddata[i][j][k]>max)
	    max = ddata[i][j][k];
    break;
  }

  return max;
}

void make_binary_wrap(Volume_wrap *wrap){
  VIO_Real max;
  byte ***bdata;
  double ***ddata;
  int i,j,k;

  max = get_wrap_max(wrap);

  ddata = (double ***)wrap->data;

  bdata = malloc(wrap->sizes[0]*sizeof(*bdata));
  bdata[0] = malloc(wrap->sizes[0]*wrap->sizes[1]*sizeof(**bdata));
  for(i=1;i<wrap->sizes[0];i++)
    bdata[i] = bdata[0] + i*wrap->sizes[1];
  bdata[0][0] = malloc(wrap->sizes[0]*wrap->sizes[1]*wrap->sizes[2]*sizeof(***bdata));
  for(i=0;i<wrap->sizes[0];i++)
    for(j=0;j<wrap->sizes[1];j++)
      bdata[i][j] = bdata[0][0] +i*wrap->sizes[1]*wrap->sizes[2] + j*wrap->sizes[2];

  for(i=0;i<wrap->sizes[0];i++)
    for(j=0;j<wrap->sizes[1];j++)
      for(k=0;k<wrap->sizes[2];k++)
	if (ddata[i][j][k]==max){
	  bdata[i][j][k] = 1;
	}else{
	  bdata[i][j][k] = 0;
	}
  
  free(ddata[0][0]);
  free(ddata[0]);
  free(ddata);

  wrap->type=NC_BYTE;
  wrap->type_size=NC_BYTE_SIZE;
  wrap->data=(void *)bdata;

}

void output_volume_wrap(Volume_wrap *wrap, VIO_STR filename){
  VIO_Volume vol;
  VIO_STR dim_names[3]={"zspace","yspace","xspace"};

  vol = create_volume(3,dim_names,wrap->type,FALSE,get_wrap_min(wrap),get_wrap_max(wrap));
  set_volume_sizes(vol,wrap->sizes);
  alloc_volume_data(vol);

  wrapVoxelDataToVolume(&vol,wrap);
  
  output_volume( filename, NC_UNSPECIFIED,0,0,0,vol,NULL,
		 (minc_output_options *)NULL);

  delete_volume(vol);
}

/*
 kate: space-indent on; hl c;indent-width 2; indent-mode cstyle;replace-tabs on;word-wrap-column 80;tab-width 8 
*/