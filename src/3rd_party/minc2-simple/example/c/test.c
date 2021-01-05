#include "minc2-simple.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc,char **argv)
{
  minc2_file_handle h;
  minc2_file_handle o;
  int ndim;
  int nelement;
  double *buffer;
  struct minc2_dimension *store_dims;
  
  if(argc<3)
  {
      fprintf(stderr,"Usage:%s <input.mnc> <output.mnc>\n",argv[0]);
      return 1;
  }
  minc2_allocate(&h);
  minc2_allocate(&o);
  
  if(minc2_open(h,argv[1])!=MINC2_SUCCESS)
  {
    fprintf(stderr,"Can't open %s for reading\n",argv[1]);
    return 1;
  }
    
  minc2_ndim(h,&ndim);
  minc2_nelement(h,&nelement);
  
  fprintf(stdout, "File:%s dimensions:%d \n",argv[1],ndim);
  
  minc2_get_store_dimensions(h,&store_dims);
  
  /*setup writing*/
  minc2_define(o,store_dims,MINC2_USHORT,MINC2_DOUBLE); /*writing to ushort volume, using double*/
  
  if(minc2_create(o,argv[2])!=MINC2_SUCCESS)
  {
    fprintf(stderr,"Can't open %s for writinf\n",argv[2]);
    return 1;
  }
  buffer=(double*)calloc(nelement,sizeof(double));
  /*reading full volume into memory buffer, using double data type, in file order*/
  minc2_setup_standard_order(h);
  minc2_setup_standard_order(o);
  
  if(minc2_load_complete_volume(h,buffer,MINC2_DOUBLE)!=MINC2_SUCCESS)
  {
    fprintf(stderr,"Error reading data from %s\n",argv[1]);
    return 1;
  }
  
  if(minc2_save_complete_volume(o,buffer,MINC2_DOUBLE)!=MINC2_SUCCESS)
  {
    fprintf(stderr,"Error writing data to %s\n",argv[2]);
    return 1;
  }
  free(buffer);
  minc2_close(h);
  minc2_close(o);
  
  /*deallocate*/
  minc2_free(h);
  minc2_free(o);
  return 0;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
