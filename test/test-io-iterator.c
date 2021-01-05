#include "minc2-simple.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

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
    minc2_file_iterator_handle input_it;
    minc2_file_iterator_handle output_it;
    
    double g_min;
    double g_max;
    double g_avg;
    int    cnt;
    
    struct minc2_dimension *store_dims;
    clock_t start,diff;
    
    /*setup reading*/
    minc2_ndim(h,&ndim);
    minc2_data_type(h,&data_type);
    minc2_storage_data_type(h,&storage_type);
    minc2_slice_ndim(h,&i_slice_dim);
    minc2_get_store_dimensions(h,&store_dims);
    
    fprintf(stdout,"Input slice dimensions:%d \n",i_slice_dim);
    /*iterate over input volume to determine data range*/
    input_it=minc2_iterator_allocate0();
    start=clock();
    cnt=0;
    g_avg=0.0;
    g_min=1e10;
    g_max=-1e10;
    
    minc2_iterator_input_start(input_it,h,MINC2_DOUBLE);
    
    do {
      double v;
      minc2_iterator_get_values(input_it,&v);
      if(v<g_min) g_min=v;
      if(v>g_max) g_max=v;
      g_avg+=v;
      cnt++;
    } while(minc2_iterator_next(input_it)==MINC2_SUCCESS);
    
    diff=clock()-start;
    g_avg/=cnt;
    
    fprintf(stdout,"Reading: Avg:%lf min:%lf max:%lf time:%lf msec\n",g_avg,g_min,g_max, (diff * 1000.0 / CLOCKS_PER_SEC));
    
        
    /*setup writing*/
    minc2_define(o,store_dims,MINC2_USHORT,MINC2_DOUBLE); /*writing to ushort volume, using double*/
    
    minc2_set_scaling(o,0,1);
    if(minc2_create(o,argv[2])==MINC2_SUCCESS)
    {
      int cnt_o=0;
      /*going to setup global scaling*/
      /*minc2_set_volume_range(o,g_min,g_max);*/
      
      minc2_slice_ndim(o,&o_slice_dim);
      fprintf(stdout,"Output slice dimensions:%d\n",o_slice_dim);
      output_it=minc2_iterator_allocate0();  
      start=clock();
      minc2_iterator_input_start(input_it,h,MINC2_DOUBLE);
      minc2_iterator_output_start(output_it,o,MINC2_DOUBLE);
      do 
      {
        double v;
        minc2_iterator_get_values(input_it,&v);
        minc2_iterator_put_values(output_it,&v);
        cnt_o++;
        minc2_iterator_next(output_it); /*have to advance to make sure we flush data to disk*/
      } while(minc2_iterator_next(input_it)==MINC2_SUCCESS );
      
      g_avg/=cnt;
      diff=clock();
      fprintf(stdout,"Writing: time:%ld msec input_cnt=%d output_cnt=%d\n",(diff * 1000 / CLOCKS_PER_SEC),cnt,cnt_o);
      minc2_iterator_free(output_it);
      /*copy metadata*/
      minc2_copy_metadata(h,o);
    } else {
      fprintf(stderr,"Can't open file %s for writing",argv[2]);
      err++;
    }
    minc2_close(h);
    minc2_close(o);
    minc2_iterator_free(input_it);
  } else {
    fprintf(stderr,"Can't open file %s for reading",argv[1]);
    err++;
  }
  
  minc2_free(h);
  minc2_free(o);
  return err;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
