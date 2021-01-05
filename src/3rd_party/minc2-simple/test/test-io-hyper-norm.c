#include "minc2-simple.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

void print_dimension_info(struct minc2_dimension *dims)
{
  while(dims->id!=MINC2_DIM_END)
  {
    fprintf(stdout,"Dimension:%s length:%d start:%f step:%f",minc2_dim_type_name(dims->id),dims->length,dims->start,dims->step);
    if(dims->irregular)
      fprintf(stdout," irregular");
    if(dims->have_dir_cos)
      fprintf(stdout," Cosines: %f %f %f",dims->dir_cos[0],dims->dir_cos[1],dims->dir_cos[2]);
    
    fprintf(stdout,"\n");
    dims++;
  }
}

int hyperslab_iterate(int *_start,const struct minc2_dimension *dims,int skip)
{
  int d=0;
  while(skip>0) {
    _start++;
    dims++;
    skip--;
  }
  
  do
  {
    _start[d]++;
    if(_start[d]>=dims[d].length) {
      memset(_start,0,sizeof(int)*(d+1));/*reset counters*/
    } else {
      return 1;
    }
    d++;
  } while(dims[d].id!=MINC2_DIM_END);
  return 0;
}

int main(int argc,char **argv)
{
  minc2_file_handle h;
  minc2_file_handle o;
  int err=0;
  /*test basic read functionality*/
  if(argc<3)
  {
      fprintf(stderr,"Usage:%s <input.mnc> <output.mnc>\n",argv[0]);
      return 1;
  }
  minc2_allocate(&h);
  minc2_allocate(&o);
  
  if(minc2_open(h,argv[1])==MINC2_SUCCESS)
  {
    int ndim;
    int data_type=-1;
    int storage_type=-1;
    int nelement=-1;
    int i_slice_dim=0;
    int o_slice_dim=0;
    
    struct minc2_dimension *store_dims;
    struct minc2_dimension *repr_dims;
    clock_t start,diff;
    
    /*setup reading*/
    minc2_ndim(h,&ndim);
    minc2_data_type(h,&data_type);
    minc2_storage_data_type(h,&storage_type);
    minc2_slice_ndim(h,&i_slice_dim);
    
    fprintf(stdout,
            "File:%s dimensions:%d data type:%s storage type:%s\n",argv[1],ndim,
                                   minc2_data_type_name(data_type), minc2_data_type_name(storage_type));
    fprintf(stdout,"nelement:%d\n",nelement);
    minc2_get_store_dimensions(h,&store_dims);
    fprintf(stdout,"File order:\n");
    print_dimension_info(store_dims);
    fprintf(stdout,"Input slice dimensions:%d\n",i_slice_dim);
    
    if(ndim<3)
    {
      fprintf(stderr,"run test on 3D volume!\n");
      return 1;
    }
    
    /*setup writing*/
    minc2_define(o,store_dims,MINC2_USHORT,MINC2_FLOAT); /*writing to flaot volume, using float*/
    /*going to setup slice scaling*/
    minc2_set_scaling(o,0,1);

    if(minc2_create(o,argv[2])==MINC2_SUCCESS)
    {
      struct minc2_dimension *repr_dims_o;
      int nelement;
      int i;
      float *buffer;
      double s_min,g_min;
      double s_max,g_max;
      double g_avg;
      int    cnt;
      
      int *_start=(int*)calloc(ndim,sizeof(int));
      int *_count=(int*)calloc(ndim,sizeof(int));
      
      minc2_slice_ndim(o,&o_slice_dim);
      fprintf(stdout,"Output slice dimensions:%d\n",o_slice_dim);
    
      minc2_get_store_dimensions(o,&repr_dims_o);
      
      /*figure out the size of the two fastest varying (in standardized order) dimensions and allocate buffer*/
      /*will read one hyperslab at a time*/
      for(i=0,nelement=1;i<o_slice_dim;i++) {
        nelement*=repr_dims_o[i].length;
        _count[i]=repr_dims_o[i].length;
        _start[i]=0;
      }
      /**/
      for(i=o_slice_dim;i<ndim;i++) {
        _count[i]=1;
        _start[i]=0;
      }
      
      buffer=(float*)calloc(nelement,sizeof(float));
      
      /*reading  volume into memory buffer, using file data type, in standardized order, by hyperslab*/
      start=clock();
      g_avg=0.0;
      cnt=0;
      
      do 
      {
        if(minc2_read_hyperslab(h,_start,_count,buffer,MINC2_FLOAT)!=MINC2_SUCCESS)
        {
          fprintf(stderr,"Error reading hyperslab\n");
          err++;
        }
        s_max=buffer[0];
        s_min=buffer[0];
        if(!cnt)
        {
          g_min=s_min;
          g_max=s_max;
        }
        for(i=0;i<nelement;i++)
        {
          g_avg+=buffer[i];
          if(buffer[i]>s_max) s_max=buffer[i];
          if(buffer[i]<s_min) s_min=buffer[i];
          cnt++;
        }
        if(s_max>g_max) g_max=s_max;
        if(s_min<g_min) g_min=s_min;
        
        if(!err && minc2_set_slice_range(o,_start,s_min,s_max)!=MINC2_SUCCESS)
        {
          fprintf(stderr,"Error setting hyperslab range\n");
          err++;
        }
        
        if(!err && minc2_write_hyperslab(o,_start,_count,buffer,MINC2_FLOAT)!=MINC2_SUCCESS )
        {
          fprintf(stderr,"Error writing hyperslab\n");
          err++;
        }
      } while( hyperslab_iterate(_start,repr_dims_o,o_slice_dim) && !err);
      
      g_avg/=cnt;
      diff=clock();
      fprintf(stdout,"Avg:%lf min:%lf max:%lf time:%ld msec\n",g_avg,g_min,g_max, (diff * 1000 / CLOCKS_PER_SEC));
      
      free(buffer);
      /*copy metadata*/
      minc2_copy_metadata(h,o);
    } else {
      fprintf(stderr,"Can't open file %s for writing",argv[2]);
      err++;
    }
    minc2_close(h);
    minc2_close(o);
  } else {
    fprintf(stderr,"Can't open file %s for reading",argv[1]);
    err++;
  }
  minc2_free(h);
  minc2_free(o);
  return err;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
