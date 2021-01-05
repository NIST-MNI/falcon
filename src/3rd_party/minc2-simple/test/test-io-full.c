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


void print_all_attributes(minc2_file_handle h)
{
  minc2_info_iterator_handle group_iterator=minc2_allocate_info_iterator();
  minc2_info_iterator_handle attr_iterator=minc2_allocate_info_iterator();
  
  if(minc2_start_group_iterator(h,group_iterator)==MINC2_SUCCESS)
  {
    while(minc2_iterator_group_next(group_iterator)==MINC2_SUCCESS)
    {
      if(minc2_start_attribute_iterator(h,minc2_iterator_group_name(group_iterator),attr_iterator)==MINC2_SUCCESS)
      {
        fprintf(stdout,"%s >>> \n",minc2_iterator_group_name(group_iterator));
        while(minc2_iterator_attribute_next(attr_iterator)==MINC2_SUCCESS)
        {
          int i;
          int attr_length=0;
          int attr_type=MINC2_UNKNOWN;
          void *attr;
          minc2_get_attribute_type(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),&attr_type);
          minc2_get_attribute_length(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),&attr_length);
          switch(attr_type)
          {
            case MINC2_STRING:
              attr=malloc(attr_length);
              minc2_read_attribute(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),attr,attr_length);
              fprintf(stdout,"\t\"%s\" = \"%s\"\n",minc2_iterator_attribute_name(attr_iterator),(char*)attr);
              free(attr);
              break;
            case MINC2_DOUBLE:
              attr=malloc(attr_length*sizeof(double));
              minc2_read_attribute(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),attr,attr_length);
              fprintf(stdout,"\t\"%s\" = %f",minc2_iterator_attribute_name(attr_iterator),*(double*)attr);
              for(i=1;i<attr_length;i++)
                fprintf(stdout," ,%f",((double*)attr)[i]);
              fprintf(stdout,"\n");
              free(attr);
              break;
            case MINC2_FLOAT:
              attr=malloc(attr_length*sizeof(float));
              minc2_read_attribute(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),attr,attr_length);
              fprintf(stdout,"\t\"%s\" = %f",minc2_iterator_attribute_name(attr_iterator),(double)*(float*)attr);
              for(i=1;i<attr_length;i++)
                fprintf(stdout," ,%f",(double)((float*)attr)[i]);
              fprintf(stdout,"\n");
              free(attr);
              break;
            case MINC2_UINT:
            case MINC2_INT:
              attr=malloc(attr_length*sizeof(int));
              minc2_read_attribute(h,minc2_iterator_group_name(group_iterator),minc2_iterator_attribute_name(attr_iterator),attr,attr_length);
              fprintf(stdout,"\t\"%s\" = %d",minc2_iterator_attribute_name(attr_iterator),(int)*(int*)attr);
              for(i=1;i<attr_length;i++)
                fprintf(stdout," ,%d",(int)((int*)attr)[i]);
              fprintf(stdout,"\n");
              free(attr);
              break;
            default:
              fprintf(stdout,"\t\"%s\" unsupported type %s\n",minc2_iterator_attribute_name(attr_iterator),minc2_data_type_name(attr_type));
          }
        }
        minc2_stop_info_iterator(attr_iterator);
        fprintf(stdout,"\n");
      }
    }
    minc2_stop_info_iterator(group_iterator);
  } else {
    fprintf(stderr,"Can't iterate over groups\n");
  }
  minc2_free_info_iterator(group_iterator);
  minc2_free_info_iterator(attr_iterator);
}

int main(int argc,char **argv)
{
  minc2_file_handle h;
  minc2_file_handle o;
  int err=0;
  char *timestamp;
  /*test basic read functionality*/
  if(argc<3)
  {
      fprintf(stderr,"Usage:%s <input.mnc> <output.mnc>\n",argv[0]);
      return 1;
  }
  minc2_allocate(&h);
  minc2_allocate(&o);
  
  timestamp=minc2_timestamp(argc,argv);
  
  if(minc2_open(h,argv[1])==MINC2_SUCCESS)
  {
    int ndim;
    int data_type=-1;
    int storage_type=-1;
    int nelement=-1;
    double *buffer;
    double f_avg,s_avg;
    double f_max,s_max;
    double f_min,s_min;
    struct minc2_dimension *store_dims;
    struct minc2_dimension *repr_dims;
    clock_t start,diff;
    int history_type;
    int history_len;
    char *history=NULL;
    char *new_history=NULL;
    
    /*setup reading*/
    minc2_ndim(h,&ndim);
    minc2_data_type(h,&data_type);
    minc2_storage_data_type(h,&storage_type);
    minc2_nelement(h,&nelement);
    
    fprintf(stdout,
            "File:%s dimensions:%d data type:%s storage type:%s\n",argv[1],ndim,
                                   minc2_data_type_name(data_type), minc2_data_type_name(storage_type));
    fprintf(stdout,"nelement:%d\n",nelement);
    minc2_get_store_dimensions(h,&store_dims);
    fprintf(stdout,"File order:\n");
    print_dimension_info(store_dims);

    if(minc2_get_attribute_type(h,"","history",&history_type)!=MINC2_SUCCESS)
    {
      fprintf(stderr,"Can't read :history type\n");
      return 1;
    }
    if(minc2_get_attribute_length(h,"","history",&history_len)!=MINC2_SUCCESS)
    {
      fprintf(stderr,"Can't read :history length\n");
      return 1;
    }
    
    history=malloc(history_len);
    if(minc2_read_attribute(h,"","history",history,history_len)!=MINC2_SUCCESS)
    {
      fprintf(stderr,"Can't read :history\n");
      return 1;
    }
    fprintf(stdout,":history type:%s history length:%d history:%s\n",minc2_data_type_name(history_type),history_len,history);
    new_history=malloc(history_len+strlen(timestamp));
    strcpy(new_history,history);
    strcat(new_history,timestamp);
    
    print_all_attributes(h);
    
    /*setup writing*/
    minc2_define(o,store_dims,MINC2_USHORT,data_type); /*writing to ushort volume, using double*/
    
    if(minc2_create(o,argv[2])==MINC2_SUCCESS)
    {
      buffer=(double*)calloc(nelement,sizeof(double));
      /*reading full volume into memory buffer, using double data type, in file order*/
      start=clock();
      if(minc2_load_complete_volume(h,buffer,MINC2_DOUBLE)==MINC2_SUCCESS)
      {
        int i;
        diff = clock() - start;
        f_avg=0.0;
        f_min=f_max=buffer[0];
        for(i=0;i<nelement;i++)
        {
          f_avg+=buffer[i];
          if(buffer[i]>f_max) f_max=buffer[i];
          if(buffer[i]<f_min) f_min=buffer[i];
        }
        f_avg/=nelement;
        fprintf(stdout,"Avg:%lf min:%lf max:%lf time:%ld msec\n",f_avg,f_min,f_max, (diff * 1000 / CLOCKS_PER_SEC));
      } else {
        fprintf(stderr,"Error reading data from %s\n",argv[1]);
        err++;
      }

      minc2_setup_standard_order(h);
      fprintf(stdout,"\nStandard order:\n");
      minc2_get_representation_dimensions(h,&repr_dims);
      print_dimension_info(repr_dims);
      
      if(minc2_write_attribute(o,"","history",new_history,history_len+strlen(timestamp),MINC2_STRING)!=MINC2_SUCCESS)
      {
        fprintf(stderr,"Can't write history entry\n");
        return 1;
      }
      
      start=clock();
      if(minc2_load_complete_volume(h,buffer,MINC2_DOUBLE)==MINC2_SUCCESS)
      {
        int i;
        diff = clock() - start;
        s_avg=0.0;
        s_min=s_max=buffer[0];
        
        for(i=0;i<nelement;i++)
        {
          s_avg+=buffer[i];
          if(buffer[i]>s_max) s_max=buffer[i];
          if(buffer[i]<s_min) s_min=buffer[i];
        }
        s_avg/=nelement;
        fprintf(stdout,"Avg:%lf min:%lf max:%lf time:%ld msec\n",s_avg,s_min,s_max, (diff * 1000 / CLOCKS_PER_SEC));
        
        minc2_setup_standard_order(o);
        start=clock();
        if(minc2_save_complete_volume(o,buffer,MINC2_DOUBLE)==MINC2_SUCCESS)
        {
          fprintf(stdout,"Saved values to %s using standard order time:%ld msec\n",argv[2],(diff * 1000 / CLOCKS_PER_SEC));
        } else {
          fprintf(stderr,"Error writing data to %s\n",argv[2]);
          err++;
        }

      } else {
        fprintf(stderr,"Error reading data from %s\n",argv[1]);
        err++;
      }
      
      free(buffer);
    } else {
      fprintf(stderr,"Can't open file %s for writing",argv[2]);
      err++;
    }
    free(new_history);
    free(history);
    minc2_close(h);
    minc2_close(o);
  } else {
    fprintf(stderr,"Can't open file %s for reading",argv[1]);
    err++;
  }
  free(timestamp);
  minc2_free(h);
  minc2_free(o);
  return err;
}


/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; */
