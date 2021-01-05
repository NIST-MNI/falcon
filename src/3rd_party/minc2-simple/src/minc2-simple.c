#include "minc2.h"
#include "minc_config.h"
#include <volume_io.h>

#include "minc2-simple.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <float.h>
/**
 * internal functions
 */
static int      _minc2_allocate_dimensions(minc2_file_handle h,int nDims);
static int      _minc2_cleanup_dimensions(minc2_file_handle h);
static mitype_t _minc2_type_to_mitype(int minc2_type);
static int      _mitype_to_minc2_type(mitype_t t);


/*these functions are defined in minc2-matrix-ops*/
/* ----------------------------- MNI Header -----------------------------------
@NAME       : build_transformation_matrix
@INPUT      : center, translations, scales, rotations
@OUTPUT     : *lt->mat - a linear transformation matrix
@RETURNS    : nothing
@DESCRIPTION: mat = (T)(C)(S)(SH)(R)(-C)
               the matrix is to be  PREmultiplied with a column vector (mat*colvec)
               when used in the application
---------------------------------------------------------------------------- */
void build_transformation_matrix(VIO_Transform *trans,
                                  double *center,
                                  double *translations,
                                  double *scales,
                                  double *shears,
                                  double *rotations);

/* extract parameters from linear transform
   trans = [scale][shear][rot]
         = [scale][shear][rz][ry][rx];
                                  */
VIO_BOOL extract2_parameters_from_matrix(VIO_Transform *trans,
                                         double *center,
                                         double *translations,
                                         double *scales,
                                         double *shears,
                                         double *rotations);


/**
 * internal representation of the minc file
 */
struct minc2_file {
  mihandle_t vol;

  int            ndims;
  int            store_type;  /*how data is stored in minc2 file*/
  int            data_type;   /*how data should be interpreted*/

  char         **dimension_name;
  misize_t      *dimension_size;
  double        *dimension_start;
  double        *dimension_step;

  struct         minc2_dimension *store_dims;
  struct         minc2_dimension *representation_dims;
 
  midimhandle_t *file_dims;
  midimhandle_t *apparent_dims;
  
  miboolean_t    slice_scaling_flag;
  miboolean_t    global_scaling_flag;
  
  miboolean_t    using_apparent_order;
  
  /*internal temporary data*/
  misize_t      *tmp_start;
  misize_t      *tmp_count;
};

/**
 * internale representation of an xfm file
 */
struct minc2_xfm_file
{
  VIO_General_transform xfm;
};


enum { ITERATOR_INFO_SIZE=256 };

struct minc2_info_iterator
{
  milisthandle_t _it;
  minc2_file_handle _minc2;

  char group_name[ITERATOR_INFO_SIZE];
  char attr_name[ITERATOR_INFO_SIZE];
};




/**
 * public functions
 */
int minc2_allocate(minc2_file_handle * h)
{
  *h=(minc2_file_handle)calloc(1,sizeof(struct minc2_file));
  return *h==NULL?MINC2_ERROR:MINC2_SUCCESS;
}

minc2_file_handle minc2_allocate0(void)
{
  minc2_file_handle h;
  if(minc2_allocate(&h)!=MINC2_SUCCESS)
    return NULL;
  return h;
}


int minc2_destroy(minc2_file_handle h)
{
  if(h->vol)
    minc2_close(h);
  return minc2_free(h);
}


int minc2_init(minc2_file_handle h)
{
  memset(h,0,sizeof(struct minc2_file));
  return MINC2_SUCCESS;
}


int minc2_free(minc2_file_handle h)
{
  if(!h)
    return MINC2_SUCCESS;
  _minc2_cleanup_dimensions(h);
  free(h);
  return MINC2_SUCCESS;
}


int minc2_open(minc2_file_handle h, const char * path)
{
  /*voxel valid range*/
  double valid_min,valid_max;
  
  /*real volume range, only awailable when slice scaling is off*/
  double volume_min=0.0,volume_max=1.0;
  miclass_t volume_data_class;
  mitype_t  store_type;
  int n_dims;
  int i;

  if ( miopen_volume(path, MI2_OPEN_READ, &h->vol) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't open minc file");
    return MINC2_ERROR;
  }
  
  if ( miget_volume_dimension_count(h->vol, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &n_dims)<0) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension count");
    return MINC2_ERROR;
  }
  
  if( _minc2_allocate_dimensions(h,n_dims)<0) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't allocate memory for dimensions");
    return MINC2_ERROR;
  }
  
  if ( miget_volume_dimensions(h->vol, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, h->ndims,
                               h->file_dims) < 0 ){
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension info");
    return MINC2_ERROR;
  }
  
  if ( miget_dimension_sizes( h->file_dims, h->ndims, h->dimension_size ) < 0 ){
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension sizes");
    return MINC2_ERROR;
  }

  if ( miget_dimension_separations(h->file_dims, MI_ORDER_FILE, h->ndims, h->dimension_step) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension step");
    return MINC2_ERROR;
  }

  if ( miget_dimension_starts(h->file_dims, MI_ORDER_FILE, h->ndims, h->dimension_start) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension step");
    return MINC2_ERROR;
  }
  
  if ( miget_data_type(h->vol, &store_type) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get minc data type");
    return MINC2_ERROR;
  }
  
  h->store_type=_mitype_to_minc2_type(store_type);
  
  if ( miget_slice_scaling_flag(h->vol, &h->slice_scaling_flag) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get slice scaling ");
    return MINC2_ERROR;
  }
  if(miget_volume_valid_range(h->vol,&valid_max,&valid_min) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get valid range");
    return MINC2_ERROR;
  }

  if( !h->slice_scaling_flag )
  {
    if( miget_volume_range(h->vol,&volume_max,&volume_min) < 0 ) {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get volume range");
      return MINC2_ERROR;
    }

    h->global_scaling_flag= !(volume_min == valid_min && volume_max == valid_max);
  }

  /*get dimension information*/
  for (i = 0; i < h->ndims; i++ )
  {
    char *      name;
    miboolean_t _sampling;
    
    /*const char *_sign="+";*/

    if ( miget_dimension_name(h->file_dims[i], &name) < 0 ) {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension name");
      return MINC2_ERROR;
    }

    h->dimension_name[i] = name;
    h->store_dims[h->ndims-i-1].length=h->dimension_size[i];
    
    /*TODO:Do not read this information for vector_dimension!*/
    if(miget_dimension_separation(h->file_dims[i],MI_FILE_ORDER,&h->store_dims[h->ndims-i-1].step)<0) {
      /*MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension step");
      return MINC2_ERROR;*/
      h->store_dims[h->ndims-i-1].step=0.0; /*set default value of 0, if there is no step size*/
    }
    
    /*TODO:Do not read this information for vector_dimension!*/
    if(miget_dimension_cosines(h->file_dims[i],&h->store_dims[h->ndims-i-1].dir_cos[0])==MI_NOERROR)
      h->store_dims[h->ndims-i-1].have_dir_cos=1;
    else
      h->store_dims[h->ndims-i-1].have_dir_cos=0;
    
    /*TODO:Do not read this information for vector_dimension!*/
    if(miget_dimension_start(h->file_dims[i],MI_FILE_ORDER,&h->store_dims[h->ndims-i-1].start)<0)
    {
      /*MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension start");
      return MINC2_ERROR;*/
      h->store_dims[h->ndims-i-1].start=0.0; /*set default value of 0, if there is no step size*/
    }

    if(miget_dimension_sampling_flag(h->file_dims[i],&_sampling)<0)
    {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get dimension sampling");
      return MINC2_ERROR;
    }
    
    h->store_dims[h->ndims-i-1].irregular=_sampling; /*documentation is wrong*/
    
    if(!strcmp(name,MIxspace) || !strcmp(name,MIxfrequency) ) /*this is X space*/
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_X;
    }
    else if(!strcmp(name,MIyspace) || !strcmp(name,MIyfrequency) ) /*this is Y
                                                                      space */
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_Y;
    }
    else if(!strcmp(name,MIzspace) || !strcmp(name,MIzfrequency) ) /*this is Z
                                                                    space*/
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_Z;
    }
    else if(!strcmp(name,MIvector_dimension) ) /*this is vector space*/
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_VEC;
    }
    else if(!strcmp(name,MItime) || !strcmp(name,MItfrequency) ) /*this is time
                                                                   space */
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_TIME;
    }
    else
    {
      h->store_dims[h->ndims-i-1].id=MINC2_DIM_TIME;
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported dimension type:%s",name);
      return MINC2_ERROR;
    }
  }
  /*mark the end of dimensions*/
  h->store_dims[h->ndims].id=MINC2_DIM_END;
  
  /*copy store to reprenetation dimension*/
  memcpy(h->representation_dims,h->store_dims,sizeof(struct minc2_dimension)*(h->ndims+1));

  if ( miget_data_class(h->vol, &volume_data_class) < 0 ) {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't get volume data class");
    return MINC2_ERROR;  
  }

  /* set the file data type*/
  if(h->slice_scaling_flag || h->global_scaling_flag)
  {
    switch ( h->store_type )
    {
      case MI_TYPE_FLOAT:
        h->data_type=MINC2_FLOAT;
        break;
      case MI_TYPE_DOUBLE:
        h->data_type=MINC2_DOUBLE;
        break;
      case MI_TYPE_FCOMPLEX:
        h->data_type=MINC2_FCOMPLEX;
        break;
      case MI_TYPE_DCOMPLEX:
        h->data_type=MINC2_DCOMPLEX;
        break;
      default:
        h->data_type=MINC2_FLOAT;
        break;
    } 
  }
  else /*not using normalization*/
  {
    switch ( h->store_type )
    {
      case MI_TYPE_BYTE:
        h->data_type=MINC2_BYTE;
        break;
      case MI_TYPE_UBYTE:
        h->data_type=MINC2_UBYTE;
        break;
      case MI_TYPE_SHORT:
        h->data_type=MINC2_SHORT;
        break;
      case MI_TYPE_USHORT:
        h->data_type=MINC2_USHORT;
        break;
      case MI_TYPE_INT:
        h->data_type=MINC2_INT;
        break;
      case MI_TYPE_UINT:
        h->data_type=MINC2_UINT;
        break;
      case MI_TYPE_FLOAT:
        h->data_type=MINC2_FLOAT;
        break;
      case MI_TYPE_DOUBLE:
        h->data_type=MINC2_DOUBLE;
        break;
      case MI_TYPE_SCOMPLEX:
        h->data_type=MINC2_SCOMPLEX;
        break;
      case MI_TYPE_ICOMPLEX:
        h->data_type=MINC2_ICOMPLEX;
        break;
      case MI_TYPE_FCOMPLEX:
        h->data_type=MINC2_FCOMPLEX;
        break;
      case MI_TYPE_DCOMPLEX:
        h->data_type=MINC2_DCOMPLEX;
        break;
      default:
        MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported file data type");
        return MINC2_ERROR;
    } 
  }

  switch ( volume_data_class )
  {
    case MI_CLASS_REAL:
    case MI_CLASS_INT:
    case MI_CLASS_LABEL: /* create an array of label names and values ?*/
/*      if(numberOfComponents == 1)
        {
        h->SetPixelType(SCALAR);
        }
      else
        {
        h->SetPixelType(VECTOR); 
        }*/
      break;
    case MI_CLASS_COMPLEX:
      /*h->SetPixelType(COMPLEX);*/
      /*numberOfComponents *= 2;*/
      break;
    default:
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported data class");
      return MINC2_ERROR;
  } //end of switch

  return MINC2_SUCCESS;
}




int minc2_slice_ndim(minc2_file_handle h,int *slice_ndim)
{
  if(h->slice_scaling_flag)
  {
    if( miget_slice_dimension_count(h->vol,MI_DIMCLASS_ANY, MI_DIMATTR_ALL, slice_ndim)<0)
      return MINC2_ERROR;
  } else {
    /*we don't have slice scaling?*/
    *slice_ndim= (h->ndims>2?2:h->ndims);
  }
  return MINC2_SUCCESS;
}


int minc2_setup_standard_order(minc2_file_handle h)
{
  /*int spatial_dimension=0;*/
  int usable_dimensions=0;
  int i;
  int dimension_indeces[5]={-1, -1, -1, -1, -1};
  
  if(!h->store_dims)
  {
    /*minc file is not opened or created yet*/
    return MINC2_ERROR;
  }
  
  /*create mapping*/
  for(i=0; i< h->ndims; i++)
  {
    switch(h->store_dims[i].id)
    {
      case MINC2_DIM_X:
        dimension_indeces[1]=i;
        break;
      case MINC2_DIM_Y:
        dimension_indeces[2]=i;
        break;
      case MINC2_DIM_Z:
        dimension_indeces[3]=i;
        break;
      case MINC2_DIM_VEC:
        dimension_indeces[0]=i;
        break;
      case MINC2_DIM_TIME:
        dimension_indeces[4]=i;
        break;
      default:
        /*error?*/
        MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported dimension");
        break;
    }
  }
  
  /*remap dimensions*/
  for(i=0; i<5; i++)
  {
    if( dimension_indeces[i]!=-1 )
    {
      h->apparent_dims[h->ndims-1-usable_dimensions]=h->file_dims[h->ndims-1-dimension_indeces[i]];
      
      /*always use positive, unless it is a vector dimension?*/
      if(i>0)
        miset_dimension_apparent_voxel_order(h->apparent_dims[h->ndims-1-usable_dimensions],MI_POSITIVE);
      
      h->representation_dims[usable_dimensions] = h->store_dims[dimension_indeces[i]];
      
      miget_dimension_separation(h->apparent_dims[h->ndims-1-usable_dimensions],MI_POSITIVE,&h->representation_dims[usable_dimensions].step);
      miget_dimension_start(     h->apparent_dims[h->ndims-1-usable_dimensions],MI_POSITIVE,&h->representation_dims[usable_dimensions].start);
      
      /*
      if(i>0 && i<4)
        spatial_dimension++;
      */
      usable_dimensions++;
    }
  }
  h->representation_dims[usable_dimensions].id=MINC2_DIM_END; /*mark the end*/

  /*Set apparent dimension order to the MINC2 api*/
  if(miset_apparent_dimension_order(h->vol, usable_dimensions, h->apparent_dims)<0)
    return MINC2_ERROR;  
  h->using_apparent_order=1;
  return MINC2_SUCCESS;
}

int minc2_close(minc2_file_handle h)
{
  if(h->vol)
  {
    if ( miclose_volume(h->vol) < 0 )
      return MINC2_ERROR;
    
    h->vol=0;
    
    return _minc2_cleanup_dimensions(h);
  } else {
    /*File was not open?*/
    return _minc2_cleanup_dimensions(h);
  }
}

int minc2_ndim(minc2_file_handle h,int *ndim)
{
  if(h->ndims)
  {
    *ndim=h->ndims;
  } else {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Dimensions undefined");
    return MINC2_ERROR;
  }
  return MINC2_SUCCESS;
}

int minc2_nelement(minc2_file_handle h,int *nelement)
{
  if(h->ndims && h->dimension_size)
  {
    int i;
    *nelement=1;
    for (i = 0; i < h->ndims; i++ )
    {
      *nelement *= h->dimension_size[i];
    }
    return MINC2_SUCCESS; 
  } else 
    return MINC2_ERROR;
}

int minc2_load_complete_volume( minc2_file_handle h,void *buffer,int representation_type)
{
  mitype_t buffer_type=MI_TYPE_UBYTE;
  int i;
  int err=MINC2_SUCCESS;

  if(h->using_apparent_order)
  {
    /*need to specify dimensions in apparent order, with minc2 convention that fasted dimensions are last*/
    for ( i = 0; i < h->ndims ; i++ )
    {
      h->tmp_start[i]=0;
      h->tmp_count[i]=h->representation_dims[h->ndims-i-1].length;
    }
  } else {
    /*will read information on file order*/
    for ( i = 0; i < h->ndims ; i++ )
    {
      h->tmp_start[i]=0;
      h->tmp_count[i]=h->store_dims[h->ndims-i-1].length;
    }
  }
  buffer_type=_minc2_type_to_mitype(representation_type);

  if ( miget_real_value_hyperslab(h->vol, buffer_type, h->tmp_start, h->tmp_count, buffer) < 0 )
    err=MINC2_ERROR;
  else
    err=MINC2_SUCCESS;
  
  return err;
}

#define \
_GET_BUFFER_MIN_MAX(type_out,buffer,buffer_length,buffer_min,buffer_max) \
  { \
    size_t _i;\
    const type_out *_buffer = (const type_out *)buffer; \
    buffer_min=buffer_max=_buffer[0]; \
    for(_i=0;_i<buffer_length;_i++,_buffer++)\
    {\
      if( *_buffer > buffer_max ) buffer_max=*_buffer; \
      if( *_buffer < buffer_min ) buffer_min=*_buffer; \
    }\
  }

#define \
_GET_BUFFER_MIN_MAX_PROT(type_out,buffer,buffer_length,buffer_min,buffer_max) \
  { \
    size_t _i;\
    const type_out *_buffer = (const type_out *)buffer; \
    buffer_min=buffer_max=FP_NAN; \
    for(_i=0;_i<buffer_length;_i++,_buffer++)\
    {\
      if( isfinite(*_buffer) ) {\
      if( *_buffer > buffer_max || !isfinite(buffer_max) ) buffer_max=*_buffer; \
      if( *_buffer < buffer_min || !isfinite(buffer_min) ) buffer_min=*_buffer; \
      } \
    }\
  }


int minc2_save_complete_volume( minc2_file_handle h,const void *buffer,int representation_type)
{
  mitype_t buffer_type=MI_TYPE_UBYTE;
  /*mitype_t file_store_type=MI_TYPE_UBYTE;*/
  
  int i;
  int err=MINC2_SUCCESS;
  size_t   buffer_length=1;
  double   buffer_min,buffer_max;

  if(h->slice_scaling_flag)
  {
    /*currently not supported*/
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Slice scaling in minc2_save_complete_volume is not supported yet");
    return MINC2_ERROR;
  }

  if(h->using_apparent_order)
  {
    /*need to specify dimensions in apparent order, with minc2 convention that fasted dimensions are last*/
    for ( i = 0; i < h->ndims ; i++ )
    {
      h->tmp_start[i]=0;
      h->tmp_count[i]=h->representation_dims[h->ndims-i-1].length;
      buffer_length*=h->tmp_count[i];
    }
  } else {
    /*will write information in file order*/
    for ( i = 0; i < h->ndims ; i++ )
    {
      h->tmp_start[i]=0;
      h->tmp_count[i]=h->store_dims[h->ndims-i-1].length;;
      buffer_length*=h->tmp_count[i];
    }
  }
  buffer_type=_minc2_type_to_mitype(representation_type);
  switch(representation_type )
  {
    case MINC2_UBYTE:
      _GET_BUFFER_MIN_MAX(unsigned char,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_BYTE:
      _GET_BUFFER_MIN_MAX(char,buffer, buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_USHORT:
      _GET_BUFFER_MIN_MAX(unsigned short,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_SHORT:
      _GET_BUFFER_MIN_MAX(short,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_UINT:
      _GET_BUFFER_MIN_MAX(unsigned int,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_INT:
      _GET_BUFFER_MIN_MAX(int,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_FLOAT:
      _GET_BUFFER_MIN_MAX_PROT(float,buffer,buffer_length,buffer_min,buffer_max);
      break;
    case MINC2_DOUBLE:
      _GET_BUFFER_MIN_MAX_PROT(double,buffer,buffer_length,buffer_min,buffer_max);
      break;
    default:
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported volume data type");
      return MINC2_ERROR;
  }
  if(minc2_set_volume_range(h,buffer_min,buffer_max)!=MINC2_SUCCESS)
    return MINC2_ERROR;

  if ( miset_real_value_hyperslab(h->vol, buffer_type, h->tmp_start, h->tmp_count, (void*)buffer ) < 0 )
    err=MINC2_ERROR;
  else
    err=MINC2_SUCCESS;
  
  return err;  
}

int minc2_set_scaling(minc2_file_handle h,int use_global_scaling,int use_slice_scaling)
{
  int err=MINC2_SUCCESS;

  if(use_global_scaling&&use_slice_scaling)
    /*can't have it both ways*/
    return MINC2_ERROR;

  h->global_scaling_flag=use_global_scaling;
  h->slice_scaling_flag=use_slice_scaling;

  return err;
}


int minc2_set_volume_range(minc2_file_handle h,
                           double value_min,
                           double value_max)
{
  int err=MINC2_SUCCESS;
  
  if( !h->global_scaling_flag )
  {
    
    if(miset_volume_valid_range( h->vol, value_max, value_min)<0) err=MINC2_ERROR;
    if(miset_volume_range(       h->vol, value_max, value_min)<0) err=MINC2_ERROR;
  }
  else // we are using scaling
  {
    if(miset_volume_range(h->vol,value_max,value_min)<0) err=MINC2_ERROR;
  }
  return err;
}

int minc2_set_slice_range(minc2_file_handle h,int *start,double value_min,double value_max)
{
  int i;
  for ( i = 0; i < h->ndims ; i++ )
  {
    h->tmp_start[i]=start[h->ndims-i-1];
  }
  if( miset_slice_range(h->vol,h->tmp_start, (size_t)h->ndims, value_max, value_min) < 0)
    return MINC2_ERROR;

  return MINC2_SUCCESS;
}


int minc2_write_hyperslab(minc2_file_handle h,int *start,int *count,const void* buffer,int representation_type)
{
  mitype_t buffer_type=_minc2_type_to_mitype(representation_type);
  int i;
  int err=MINC2_SUCCESS;
  
  /*need to specify dimensions with minc2 convention that fasted dimensions are last*/
  for ( i = 0; i < h->ndims ; i++ )
  {
    h->tmp_start[i]=start[h->ndims-i-1];
    h->tmp_count[i]=count[h->ndims-i-1];
  }

  if ( miset_real_value_hyperslab(h->vol, buffer_type, h->tmp_start, h->tmp_count, (void*)buffer ) < 0 )
    err=MINC2_ERROR;
  else
    err=MINC2_SUCCESS;

  return err;
}

int minc2_read_hyperslab(minc2_file_handle h, int *start, int *count, void* buffer, int representation_type)
{
  mitype_t buffer_type=_minc2_type_to_mitype(representation_type);
  int i;
  int err=MINC2_SUCCESS;

  /*need to specify dimensions with minc2 convention that fasted dimensions are last*/
  for ( i = 0; i < h->ndims ; i++ )
  {
    h->tmp_start[i]=start[h->ndims-i-1];
    h->tmp_count[i]=count[h->ndims-i-1];
  }

  if ( miget_real_value_hyperslab(h->vol, buffer_type, h->tmp_start, h->tmp_count, buffer) < 0 )
    err=MINC2_ERROR;
  else
    err=MINC2_SUCCESS;

  return err;
}

int minc2_data_type(minc2_file_handle h,int *_type)
{
  if(h->data_type>0)
  {
    *_type=h->data_type;
    return MINC2_SUCCESS;
  } else {
    /*not initialized!*/
    return MINC2_ERROR;
  }
}

int minc2_storage_data_type(minc2_file_handle h,int *_type)
{
  if(h->data_type>0)
  {
    *_type=(int)h->store_type;
    return MINC2_SUCCESS;
  } else {
    /*not initialized!*/
    return MINC2_ERROR;
  }
}

int minc2_get_representation_dimensions(minc2_file_handle h,struct minc2_dimension **dims)
{
  if(!h->representation_dims)
    return MINC2_ERROR;
  *dims=h->representation_dims;
  return MINC2_SUCCESS;
}

int minc2_get_store_dimensions(minc2_file_handle h,struct minc2_dimension **dims)
{
  if(!h->store_dims)
    return MINC2_ERROR;
  *dims=h->store_dims;
  return MINC2_SUCCESS;
}

int minc2_compare_voxel_dimensions(const struct minc2_dimension *one,const struct minc2_dimension *two)
{
  while(one->id!=MINC2_DIM_END && two->id!=MINC2_DIM_END)
  {
    if(one->id!=two->id)
      return MINC2_ERROR;
    if(one->length!=two->length)
      return MINC2_ERROR;

    one++;
    two++;
  };

  return one->id==MINC2_DIM_END && two->id==MINC2_DIM_END?MINC2_SUCCESS:MINC2_ERROR;
}


int minc2_compare_dimensions(const struct minc2_dimension *one,const struct minc2_dimension *two)
{
  int i;
  while(one->id!=MINC2_DIM_END && two->id!=MINC2_DIM_END)
  {
    if(one->id!=two->id)
      return MINC2_ERROR;

    if(one->length!=two->length)
      return MINC2_ERROR;

    if(fabs(one->step-two->step)>1e-6)
      return MINC2_ERROR;

    if(one->irregular!=two->irregular)
      return MINC2_ERROR;

    if(fabs(one->start-two->start)>1e-6)
      return MINC2_ERROR;

    if(one->have_dir_cos && two->have_dir_cos)
    {
      for(i=0;i<3;i++)
        if(fabs(one->dir_cos[i]-two->dir_cos[i])>1e-6)
          return MINC2_ERROR;
    }

    one++;
    two++;
  };

  return one->id==MINC2_DIM_END && two->id==MINC2_DIM_END?MINC2_SUCCESS:MINC2_ERROR;
}


int minc2_define(minc2_file_handle h, struct minc2_dimension *store_dims, int store_data_type,int data_type)
{
  int i;
  int ndims=0;
  struct minc2_dimension * dim;
  /*figure out number of dimension*/
  for(dim=store_dims;dim->id!=MINC2_DIM_END;dim++)
  {
    ndims++;
  }
  
  h->store_type=store_data_type; /*TODO: add mapping?*/
  h->data_type=data_type;
  
  /*TODO: add more cases*/
  if( ( h->store_type==MI_TYPE_FLOAT ||
        h->store_type==MI_TYPE_DOUBLE ) || 
        
      ( 
        ( h->store_type==MI_TYPE_BYTE  || h->store_type==MI_TYPE_INT  || h->store_type==MI_TYPE_SHORT ||
          h->store_type==MI_TYPE_UBYTE || h->store_type==MI_TYPE_UINT || h->store_type==MI_TYPE_USHORT ) && 
        ( h->data_type==MI_TYPE_BYTE   || h->data_type==MI_TYPE_INT   || h->data_type==MI_TYPE_SHORT ||
          h->data_type==MI_TYPE_UBYTE  || h->data_type==MI_TYPE_UINT  || h->data_type==MI_TYPE_USHORT )
      )
    )
  {
    h->slice_scaling_flag=0;
    h->global_scaling_flag=0;
  } else {
    h->slice_scaling_flag=0; /*TODO: use slice scaling sometimes?*/
    h->global_scaling_flag=1;
  }
  
  _minc2_allocate_dimensions(h,ndims);
  memcpy(h->store_dims         ,store_dims,sizeof(struct minc2_dimension)*(h->ndims+1));
  memcpy(h->representation_dims,store_dims,sizeof(struct minc2_dimension)*(h->ndims+1));
  
  for(dim=store_dims,i=0;dim->id!=MINC2_DIM_END;dim++,i++)
  {
    switch(dim->id)
    {
    case MINC2_DIM_X:
      micreate_dimension(MIxspace,MI_DIMCLASS_SPATIAL, 
                         dim->irregular?MI_DIMATTR_NOT_REGULARLY_SAMPLED:MI_DIMATTR_REGULARLY_SAMPLED, 
                         dim->length,
                         &h->file_dims[ndims-i-1] );
      break;
    case MINC2_DIM_Y:
      micreate_dimension(MIyspace,MI_DIMCLASS_SPATIAL, 
                         dim->irregular?MI_DIMATTR_NOT_REGULARLY_SAMPLED:MI_DIMATTR_REGULARLY_SAMPLED, 
                         dim->length,
                         &h->file_dims[ndims-i-1] );
      break;
    case MINC2_DIM_Z:
      micreate_dimension(MIzspace,MI_DIMCLASS_SPATIAL, 
                         dim->irregular?MI_DIMATTR_NOT_REGULARLY_SAMPLED:MI_DIMATTR_REGULARLY_SAMPLED, 
                         dim->length,
                         &h->file_dims[ndims-i-1] );
      break;
    case MINC2_DIM_TIME:
      micreate_dimension(MItime, MI_DIMCLASS_TIME, 
                         dim->irregular?MI_DIMATTR_NOT_REGULARLY_SAMPLED:MI_DIMATTR_REGULARLY_SAMPLED, 
                         dim->length,
                         &h->file_dims[ndims-i-1] );
      break;
    case MINC2_DIM_VEC:
      micreate_dimension(MIvector_dimension,MI_DIMCLASS_RECORD, MI_DIMATTR_REGULARLY_SAMPLED, 
                         dim->length,
                         &h->file_dims[ndims-i-1] );
      break;
    default:
      /*don't know this dimension type*/
      /*TODO: report error*/
      break;
    }
    miset_dimension_start(     h->file_dims[ndims-i-1],dim->start);
    miset_dimension_separation(h->file_dims[ndims-i-1],dim->step );
    if(dim->have_dir_cos)
      miset_dimension_cosines( h->file_dims[ndims-i-1],dim->dir_cos);
  }
  return MINC2_SUCCESS;
}


int minc2_create(minc2_file_handle h,const char * path)
{
  int err=MINC2_SUCCESS;
  /**/
  mivolumeprops_t hprops;
  
  if( minew_volume_props(&hprops) < 0)
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't set volume properties");
    return MINC2_ERROR;
  }
  
  /*TODO: move it to volume definition*/
  if(miget_cfg_present(MICFG_COMPRESS) && miget_cfg_int(MICFG_COMPRESS)>0  )
  {
    if(miset_props_compression_type(hprops, MI_COMPRESS_ZLIB)<0)
    {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't set compression");
      return MINC2_ERROR;
    }

    if(miset_props_zlib_compression(hprops,miget_cfg_int(MICFG_COMPRESS))<0)
    {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't set compression");
      return MINC2_ERROR;
    }
  }
  else
  {
    if(miset_props_compression_type(hprops, MI_COMPRESS_NONE)<0)
    {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Can't set compression");
      return MINC2_ERROR;
    }
  }

  if ( micreate_volume ( path, h->ndims, h->file_dims, h->store_type,
                         MI_CLASS_REAL, hprops, &h->vol )<0 ) /*change MI_CLASS_REAL to something else?*/
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Couldn't open file %s",path);
    return MINC2_ERROR;
  }

  /*have to set slice scaling flag before image is allocated*/
  if ( miset_slice_scaling_flag(h->vol, h->slice_scaling_flag )<0 )
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Couldn't set slice scaling");
    return MINC2_ERROR;
  }

  if (  micreate_volume_image ( h->vol ) <0 )
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Couldn't create image in file %s",path);
    return MINC2_ERROR;
  }

  if(h->global_scaling_flag || h->slice_scaling_flag)
  {
    switch(h->store_type)
    {
      case MI_TYPE_BYTE:
        if(miset_volume_valid_range(h->vol,SCHAR_MAX-1,SCHAR_MIN)<0) err=MINC2_ERROR;
        break;
      case MI_TYPE_UBYTE:
        if(miset_volume_valid_range(h->vol,UCHAR_MAX-1,0)<0) err=MINC2_ERROR;
        break;
      case MI_TYPE_SHORT:
        if(miset_volume_valid_range(h->vol,SHRT_MAX-1,SHRT_MIN)<0) err=MINC2_ERROR;
        break;
      case MI_TYPE_USHORT:
        if(miset_volume_valid_range(h->vol,USHRT_MAX-1,0)<0) err=MINC2_ERROR;
        break;
      case MI_TYPE_INT:
        if(miset_volume_valid_range(h->vol,INT_MAX-1,INT_MIN)<0) err=MINC2_ERROR;
        break;
      case MI_TYPE_UINT:
        if(miset_volume_valid_range(h->vol,UINT_MAX-1,0)<0) err=MINC2_ERROR;
        break;
      default:
        /*error*/
        MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported store data type");
        return MINC2_ERROR;
    }
  }

  return err;
}


int minc2_world_to_voxel(minc2_file_handle h,const double *world,double *voxel)
{
  if(!h->vol)
    return MINC2_ERROR;

  if(miconvert_world_to_voxel(h->vol, world, voxel)<0)
    return MINC2_ERROR;

  return MINC2_SUCCESS;
}


int minc2_voxel_to_world(minc2_file_handle h, const double *voxel, double *world)
{
  if(!h->vol)
    return MINC2_ERROR;

  if(miconvert_voxel_to_world(h->vol, voxel, world)<0)
    return MINC2_ERROR;

  return MINC2_SUCCESS;
}


int minc2_world_to_voxel_vec(minc2_file_handle h,int n,int stride, const double *world,double *voxel)
{
  int i;
  int ret=MINC2_SUCCESS;
  if(!h->vol)
    return MINC2_ERROR;

  for(i=0;i<n;i++ ) {
    const double* in_world=&world[i*stride];
    double *out_voxel=&voxel[i*stride];

    ret = ret || minc2_world_to_voxel(h, in_world, out_voxel);
  }

  return ret;
}


int minc2_voxel_to_world_vec(minc2_file_handle h, int n,int stride, const double *voxel, double *world)
{
  int i;
  int ret=MINC2_SUCCESS;
  if(!h->vol)
    return MINC2_ERROR;

  for(i=0;i<n;i++ ) {
    const double* in_voxel=&voxel[i*stride];
    double *out_world=&voxel[i*stride];

    ret =  ret || minc2_voxel_to_world(h, in_voxel, out_world);
  }

  return ret;
}



static int _minc2_cleanup_dimensions(minc2_file_handle h)
{
  int i;
  if( h->dimension_name )
  {
    for ( i = 0; i < h->ndims; i++ )
    {
      if(h->dimension_name[i]) 
        mifree_name( h->dimension_name[i] );
      h->dimension_name[i]=NULL;
    }
    free(h->dimension_name);
  }

  
  if(h->dimension_size)  free(h->dimension_size);
  if(h->dimension_start) free(h->dimension_start);
  if(h->dimension_step)  free(h->dimension_step);
  if(h->file_dims)       free(h->file_dims);
  if(h->apparent_dims)   free(h->apparent_dims);
  if(h->store_dims)      free(h->store_dims);
  if(h->representation_dims)free(h->representation_dims);
  if(h->tmp_start)       free(h->tmp_start);
  if(h->tmp_count)       free(h->tmp_count);
  
  h->dimension_name    = NULL;
  h->dimension_size    = NULL;
  h->dimension_start   = NULL;
  h->dimension_step    = NULL;
  h->file_dims         = NULL;
  h->apparent_dims     = NULL;
  h->store_dims        = NULL;
  h->representation_dims=NULL;
  h->using_apparent_order=0;
  h->tmp_count         =NULL;
  h->tmp_start         =NULL;
  
  return MINC2_SUCCESS;
}

static int _minc2_allocate_dimensions(minc2_file_handle h,int nDims)
{
  _minc2_cleanup_dimensions(h);

  h->ndims=nDims;

  h->dimension_name  = (char**)         calloc(h->ndims,sizeof(char*));
  h->dimension_size  = (misize_t*)      calloc(h->ndims,sizeof(misize_t));
  h->dimension_start = (double*)        calloc(h->ndims,sizeof(double));
  h->dimension_step  = (double*)        calloc(h->ndims,sizeof(double));
  h->file_dims       = (midimhandle_t*) calloc(h->ndims,sizeof(midimhandle_t));
  h->apparent_dims   = (midimhandle_t*) calloc(h->ndims,sizeof(midimhandle_t));
  h->store_dims      = (struct minc2_dimension*)calloc(h->ndims+1,sizeof(struct minc2_dimension));
  h->representation_dims= (struct minc2_dimension*)calloc(h->ndims+1,sizeof(struct minc2_dimension));

  h->tmp_start       = (misize_t *)calloc(h->ndims,sizeof(misize_t));
  h->tmp_count       = (misize_t *)calloc(h->ndims,sizeof(misize_t));

  /*TODO: check if memory was allocated?*/
  return MINC2_SUCCESS;
}

int minc2_copy_metadata(minc2_file_handle src,minc2_file_handle dst)
{
  milisthandle_t grplist;
  int err=0;
  if( !src->vol || !dst->vol)
    return MINC2_ERROR;

  if ( (milist_start(src->vol, "", 0, &grplist) ) == MI_NOERROR )
  {
      char           group_name[256];
      /*milisthandle_t attlist;*/
      while( milist_grp_next(grplist, group_name, sizeof(group_name) ) == MI_NOERROR )
      {
        if(micopy_attr(src->vol,group_name,dst->vol)<0)
          err++;
      }
    milist_finish(grplist);
  } else {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Error iterating through metadata");
    return MINC2_ERROR;
  }
  /*TODO: copy history attribute, because micopy_attr doesn't copy it*/


  return err>0?MINC2_ERROR:MINC2_SUCCESS;
}


const char * minc2_data_type_name(int minc2_type_id)
{
  switch(minc2_type_id )
    {
    case MINC2_UBYTE:
      return "unsigned char";
    case MINC2_BYTE:
      return "char";
    case MINC2_USHORT:
      return "unsigned short";
    case MINC2_SHORT:
      return "short";
    case MINC2_UINT:
      return "unsigned int";
    case MINC2_INT:
      return "int";
    case MINC2_FLOAT:
      return "float";
    case MINC2_DOUBLE:
      return "double";
    case MINC2_STRING:
      return "string";
    default:
      return "Unknown";
  }
}

const char * minc2_dim_type_name(int minc2_dim_id)
{
  switch(minc2_dim_id)
  {
    case MINC2_DIM_X:
      return "X";
    case MINC2_DIM_Y:
      return "Y";
    case MINC2_DIM_Z:
      return "Z";
    case MINC2_DIM_TIME:
      return "Time";
    case MINC2_DIM_VEC:
      return "Vector";
    case MINC2_DIM_UNKNOWN:
    case MINC2_DIM_END:
    default:
      return "Unknown";
  }
}

static mitype_t _minc2_type_to_mitype(int minc2_type)
{
  /*this is identity transform at the moment*/
  return (mitype_t)minc2_type;
  /*
    switch(representation_type )
  {
    case MINC2_UBYTE:
      buffer_type=MI_TYPE_UBYTE;
      break;
    case MINC2_BYTE:
      buffer_type=MI_TYPE_BYTE;
      break;
    case MINC2_USHORT:
      buffer_type=MI_TYPE_USHORT;
      break;
    case MINC2_SHORT:
      buffer_type=MI_TYPE_SHORT;
      break;
    case MINC2_UINT:
      buffer_type=MI_TYPE_UINT;
      break;
    case MINC2_INT:
      buffer_type=MI_TYPE_INT;
      break;
    case MINC2_FLOAT:
      buffer_type=MI_TYPE_FLOAT;
      break;
    case MINC2_DOUBLE:
      buffer_type=MI_TYPE_DOUBLE;
      break;
    default:
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported volume data type");
      free( start ); 
      free( count );
      return MINC2_ERROR;
  }
   */
}

static int _mitype_to_minc2_type(mitype_t t)
{
  /*this is identity transform at the moment*/
  return (int)t;
}

static int _minc2_type_size(int minc2_type_id)
{
  switch(minc2_type_id )
    {
    case MINC2_UBYTE:
      return sizeof(unsigned char);
    case MINC2_BYTE:
      return sizeof(char);
    case MINC2_USHORT:
      return sizeof(unsigned short);
    case MINC2_SHORT:
      return sizeof(short);
    case MINC2_UINT:
      return sizeof(unsigned int);
    case MINC2_INT:
      return sizeof(int);
    case MINC2_FLOAT:
      return sizeof(float);
    case MINC2_DOUBLE:
      return sizeof(double);
    case MINC2_STRING:
      return sizeof(char);
    default:
      return 1; /*ERROR?*/
  }
}


minc2_info_iterator_handle minc2_allocate_info_iterator(void)
{
  return calloc(1,sizeof(struct minc2_info_iterator));
}


int minc2_stop_info_iterator(minc2_info_iterator_handle it)
{
  if(!it)
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"NULL pointer");
    return MINC2_ERROR;
  }
  if(it->_it)
  {
    milist_finish(it->_it);
    it->_it=NULL;
  }
  return MINC2_SUCCESS;
}

int minc2_start_group_iterator(minc2_file_handle h,minc2_info_iterator_handle it)
{
  it->_minc2=h;
  *it->group_name=0;
  if(milist_start(it->_minc2->vol, "", 0, &it->_it) == MI_NOERROR)
    return MINC2_SUCCESS;
  it->_it=0;

  return MINC2_ERROR;
}

int minc2_start_attribute_iterator(minc2_file_handle h,const char* group,minc2_info_iterator_handle it)
{
  it->_minc2=h;
  strncpy(it->group_name,group,ITERATOR_INFO_SIZE);
  *it->attr_name=0;
  if(milist_start(it->_minc2->vol, it->group_name, 1, &it->_it) == MI_NOERROR)
    return MINC2_SUCCESS;
  it->_it=0;
  return MINC2_ERROR;
}

int minc2_iterator_group_next(minc2_info_iterator_handle it)
{
  if(!it)
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"NULL pointer");
    return MINC2_ERROR;
  }
  if(!it->_it)
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"Iterator not started");
    return MINC2_ERROR;
  }
  return milist_grp_next(it->_it, it->group_name, ITERATOR_INFO_SIZE ) == MI_NOERROR?MINC2_SUCCESS:MINC2_ERROR;
}


int minc2_iterator_attribute_next(minc2_info_iterator_handle it)
{
  if(!it)
  {
    MI_LOG_ERROR(MI2_MSG_GENERIC,"NULL pointer");
    return MINC2_ERROR;
  }
  if(!it->_it)
  {
            MI_LOG_ERROR(MI2_MSG_GENERIC,"Iterator not started");
    return MINC2_ERROR;
  }

  return milist_attr_next(it->_minc2->vol,it->_it,it->group_name,ITERATOR_INFO_SIZE,it->attr_name,ITERATOR_INFO_SIZE ) == MI_NOERROR?MINC2_SUCCESS:MINC2_ERROR;
}

const char* minc2_iterator_group_name(minc2_info_iterator_handle it)
{
  if(!it || !it->_it)
    return "";
  return it->group_name;
}

const char* minc2_iterator_attribute_name(minc2_info_iterator_handle it)
{
  if(!it || !it->_it)
    return "";
  return it->attr_name;
}


int minc2_free_info_iterator(minc2_info_iterator_handle it)
{
  int err=MINC2_SUCCESS;
  if(!it) return MINC2_ERROR;
  err=minc2_stop_info_iterator(it);
  free(it);
  return err;
}



int minc2_get_attribute_type(minc2_file_handle h,const char* group,const char* attr,int *minc2_type)
{
  mitype_t    att_data_type;
  *minc2_type=MINC2_UNKNOWN;
  if(miget_attr_type(h->vol,group,attr,&att_data_type)== MI_NOERROR)
  {
    *minc2_type=_mitype_to_minc2_type(att_data_type);
    return MINC2_SUCCESS;
  }

  return MINC2_ERROR;
}

int minc2_get_attribute_length(minc2_file_handle h,const char* group,const char* attr,int *attr_length)
{
  size_t      _att_length;
  *attr_length=0;
  if(miget_attr_length(h->vol,group,attr,&_att_length) == MI_NOERROR)
  {
    *attr_length=(int)_att_length;
    return MINC2_SUCCESS;
  }
  return MINC2_ERROR;

}

int minc2_read_attribute(minc2_file_handle h,const char* group,const char* attr,void *buf,int buf_size)
{
  mitype_t    att_data_type;
  if( miget_attr_type(h->vol,group,attr,&att_data_type)== MI_NOERROR &&
      miget_attr_values(h->vol,att_data_type,group,attr,buf_size,buf) == MI_NOERROR )
  {
    return MINC2_SUCCESS;
  }

  return MINC2_ERROR;
}

int minc2_write_attribute(minc2_file_handle h,const char* group,const char* attr,const void *buf,int buf_size,int minc2_type)
{
  mitype_t    att_data_type=_minc2_type_to_mitype(minc2_type);
  if(att_data_type==MI_TYPE_UNKNOWN) return MINC2_ERROR;

  if(miset_attr_values(h->vol,att_data_type,group,attr,buf_size,buf ) == MI_NOERROR)
  {
    return MINC2_SUCCESS;
  }
  return MINC2_ERROR;
}

int minc2_delete_attribute(minc2_file_handle h,const char* group,const char* attr)
{
  return midelete_attr ( h->vol, group, attr )==MI_NOERROR?MINC2_SUCCESS:MINC2_ERROR;
}


int minc2_delete_group(minc2_file_handle h,const char* group)
{
  return midelete_group ( h->vol, "",group )==MI_NOERROR?MINC2_SUCCESS:MINC2_ERROR;
}


char* minc2_timestamp(int argc,char **argv)
{
  char cur_time[200];
  time_t t;
  struct tm *tmp;
  char *out=NULL;
  int total_len=0;
  int i;

  t = time(NULL);
  tmp = localtime(&t);

  strftime(cur_time, sizeof(cur_time), "%a %b %d %T %Y>>>", tmp);
  total_len=strlen(cur_time);

  for (i=0; i<argc; i++) {
    total_len+=strlen(argv[i])+2;
  }

  out=malloc(total_len+1);
  strcpy(out,cur_time);
  /* Copy the program name and arguments */
  for (i=0; i<argc; i++) {
    strcat(out,argv[i]);
    strcat(out," ");
  }
  strcat(out,"\n");

  return out;
}


int minc2_xfm_allocate(minc2_xfm_file_handle * h)
{
  *h=(minc2_xfm_file_handle)calloc(1,sizeof(struct minc2_xfm_file));
  return *h==NULL?MINC2_ERROR:MINC2_SUCCESS;
}

minc2_xfm_file_handle minc2_xfm_allocate0(void)
{
  minc2_xfm_file_handle h;
  if(minc2_xfm_allocate(&h) != MINC2_SUCCESS)
    return NULL;
  return h;
}

int minc2_xfm_init(minc2_xfm_file_handle h)
{
  memset(h,0,sizeof(struct minc2_xfm_file));
  return MINC2_SUCCESS;
}


int minc2_xfm_free(minc2_xfm_file_handle h)
{
  if(!h)
    return MINC2_SUCCESS;
  free(h);
  return MINC2_SUCCESS;
}


int minc2_xfm_destroy(minc2_xfm_file_handle h)
{
  delete_general_transform(&h->xfm);
  return minc2_xfm_free(h);
}


int minc2_xfm_open(minc2_xfm_file_handle h,const char * path)
{
  if(input_transform_file((char*)path, &h->xfm)!=VIO_OK)
    return MINC2_ERROR;
  return MINC2_SUCCESS;
}

int minc2_xfm_save(minc2_xfm_file_handle h,const char * path)
{
  if(output_transform_file(path,(char*)"minc2-simple",&h->xfm)!=VIO_OK)
    return MINC2_ERROR;
  return MINC2_SUCCESS;
}

int minc2_xfm_transform_point(minc2_xfm_file_handle h,const double* in,double* out)
{
  return general_transform_point(&h->xfm,in[0],in[1],in[2],&out[0],&out[1],&out[2])==VIO_OK?MINC2_SUCCESS:MINC2_ERROR;
}

int minc2_xfm_transform_point_vec(minc2_xfm_file_handle h,int n,int stride,const double* in,double* out)
{
    int i;
    int ret=MINC2_SUCCESS;
    for(i=0;i<n;i++)
    {
        double *pnt_in=&in[i*stride];
        double *pnt_out=&out[i*stride];
        ret = ret ||  minc2_xfm_transform_point(h,pnt_in,pnt_out);
    }

    return ret;
}


int minc2_xfm_inverse_transform_point(minc2_xfm_file_handle h,const double* in,double* out)
{
  return general_inverse_transform_point(&h->xfm,in[0],in[1],in[2],&out[0],&out[1],&out[2])==VIO_OK?MINC2_SUCCESS:MINC2_ERROR;
}

int minc2_xfm_inverse_transform_point_vec(minc2_xfm_file_handle h,int n,int stride,const double* in,double* out)
{
    int i;
    int ret=MINC2_SUCCESS;
    for(i=0;i<n;i++)
    {
        double *pnt_in=&in[i*stride];
        double *pnt_out=&out[i*stride];
        ret = ret || minc2_xfm_inverse_transform_point(h,pnt_in,pnt_out);
    }

    return ret;
}


int minc2_xfm_invert(minc2_xfm_file_handle h)
{
  invert_general_transform(&h->xfm);
  return MINC2_SUCCESS;
}

int minc2_xfm_get_n_concat(minc2_xfm_file_handle h,int *n)
{
  /*heuristic check if transform is empty*/
  if(h->xfm.type==0 && h->xfm.linear_transform==0)
  {
    *n=0;
    return MINC2_SUCCESS;
  }

  switch(get_transform_type(&h->xfm))
  {
    default:
      /*unsupported*/
      *n=0;
      return MINC2_SUCCESS;
    case GRID_TRANSFORM:
    case LINEAR:
      *n=1;
      return MINC2_SUCCESS;
    case  CONCATENATED_TRANSFORM:
      *n=get_n_concated_transforms(&h->xfm);
      return MINC2_SUCCESS;
  }
}

int _minc2_xfm_type_convert(VIO_General_transform *_xfm)
{
  switch(get_transform_type(_xfm))
  {
    case THIN_PLATE_SPLINE:
      return MINC2_XFM_THIN_PLATE_SPLINE;
    case GRID_TRANSFORM:
      return MINC2_XFM_GRID_TRANSFORM;
    case LINEAR:
      return MINC2_XFM_LINEAR;
    case USER_TRANSFORM:
      return MINC2_XFM_USER_TRANSFORM;
    default :
      return MINC2_XFM_END; /* Unsupported ? */
  }
}

VIO_General_transform *_get_nth_transform(VIO_General_transform *_xfm,int n)
{
  if( get_transform_type(_xfm)==CONCATENATED_TRANSFORM &&
      n<get_n_concated_transforms(_xfm) )
  {
    return get_nth_general_transform(_xfm, n);
  } else if(n==0) {
    return _xfm;
  } else
    return NULL;

}

int minc2_xfm_get_n_type(minc2_xfm_file_handle h,int n,int *xfm_type)
{
  VIO_General_transform *_xfm=_get_nth_transform(&h->xfm, n);
  if(_xfm)
  {
    *xfm_type=_minc2_xfm_type_convert(_xfm);
    return MINC2_SUCCESS;
  } else {
    return MINC2_ERROR;
  }
}

int minc2_xfm_get_linear_transform(minc2_xfm_file_handle h,int n,double *matrix)
{
  int i,j;
  VIO_Transform *lin;
  VIO_General_transform *_xfm=_get_nth_transform(&h->xfm, n);
  if(_xfm)
  {
    lin=get_linear_transform_ptr(_xfm);

    for(j = 0; j < 4; ++j)
    {
      for(i = 0; i < 4; ++i)
      {
        matrix[j*4+i]=Transform_elem(*lin,j,i);
      }
    }
    return MINC2_SUCCESS;
  } else {
    return MINC2_ERROR;
  }
}

int minc2_xfm_get_grid_transform(minc2_xfm_file_handle h,int n,int *inverted,char **grid_file)
{

  VIO_General_transform *_xfm=_get_nth_transform(&h->xfm, n);
  if(_xfm)
  {
    *grid_file=strdup(_xfm->displacement_volume_file);
    *inverted=_xfm->inverse_flag;
    return MINC2_SUCCESS;
  } else {
    return MINC2_ERROR;
  }
}

int minc2_xfm_append_linear_transform(minc2_xfm_file_handle h,double *matrix)
{
  int n;
  int i,j;
  VIO_Transform lin;
  memset(&lin, 0, sizeof(VIO_Transform));

  for(j = 0; j < 4; ++j)
  {
    for(i = 0; i < 4; ++i)
    {
      Transform_elem(lin,j,i)=matrix[j*4+i];
    }
  }
  minc2_xfm_get_n_concat(h,&n);

  if(n==0) /*first transform*/
  {
    create_linear_transform(&h->xfm, &lin);
    return MINC2_SUCCESS;
  } else {
    VIO_General_transform lin_xfm;
    VIO_General_transform concated;

    memset(&lin_xfm, 0, sizeof(VIO_General_transform));
    create_linear_transform(&lin_xfm, &lin);

    concat_general_transforms( &h->xfm, &lin_xfm, &concated );
    delete_general_transform( &h->xfm );
    delete_general_transform( &lin_xfm );
    h->xfm = concated;
    return MINC2_SUCCESS;
  }
}

int minc2_xfm_append_grid_transform(minc2_xfm_file_handle h,const char * grid_path,int inv)
{
  int n;
  VIO_General_transform *_xfm;
  minc2_xfm_get_n_concat(h,&n);

  /*TODO: create another function that will work with volume*/
  if(n==0)
  {
    create_grid_transform_no_copy( &h->xfm, 0, 0 );
    if(inv) h->xfm.inverse_flag=TRUE;
    return MINC2_SUCCESS;
  }
  else
  {
    VIO_General_transform concated;
    VIO_General_transform nl_xfm;

    memset(&nl_xfm, 0, sizeof(VIO_General_transform));
    create_grid_transform_no_copy( &nl_xfm, 0, 0 );
    if(inv) nl_xfm.inverse_flag=TRUE;

    concat_general_transforms( &h->xfm, &nl_xfm, &concated );
    delete_general_transform( &h->xfm );
    delete_general_transform( &nl_xfm );

    h->xfm = concated;
    return MINC2_SUCCESS;
  }
}

int minc2_xfm_concat_xfm(minc2_xfm_file_handle h,minc2_xfm_file_handle o)
{
  int n;
  VIO_General_transform *_xfm;
  minc2_xfm_get_n_concat(h,&n);

  if(n==0)
  {
    copy_general_transform(&o->xfm,&h->xfm);
  } else {
    VIO_General_transform concated;
    concat_general_transforms( &h->xfm, &o->xfm, &concated );
    delete_general_transform( &h->xfm );
    h->xfm = concated;
  }
  return MINC2_SUCCESS;
}



int minc2_xfm_append_linear_param(minc2_xfm_file_handle h,
                              double *center,
                              double *translations,
                              double *scales,
                              double *shears,
                              double *rotations)
{
  int n;
  int i,j;
  VIO_Transform lin;
  memset(&lin, 0, sizeof(VIO_Transform));

  build_transformation_matrix(&lin,
                              center, translations,
                              scales, shears, rotations);
                              
  minc2_xfm_get_n_concat(h,&n);

  if(n==0) /*first transform*/
  {
    create_linear_transform(&h->xfm, &lin);
    return MINC2_SUCCESS;
  } else {
    VIO_General_transform lin_xfm;
    VIO_General_transform concated;

    memset(&lin_xfm, 0, sizeof(VIO_General_transform));
    create_linear_transform(&lin_xfm, &lin);

    concat_general_transforms( &h->xfm, &lin_xfm, &concated );
    delete_general_transform( &h->xfm );
    delete_general_transform( &lin_xfm );
    h->xfm = concated;
    return MINC2_SUCCESS;
  }
  
}

int minc2_xfm_extract_linear_param(minc2_xfm_file_handle h,
                             int n,
                             double *center,
                             double *translations,
                             double *scales,
                             double *shears,
                             double *rotations)
{
  int i,j;
  VIO_Transform *lin;
  VIO_General_transform *_xfm=_get_nth_transform(&h->xfm, n);
  if(_xfm)
  {
    lin=get_linear_transform_ptr(_xfm);
    
    return extract2_parameters_from_matrix(lin, center,
                                    translations,
                                    scales,
                                    shears,rotations)?MINC2_SUCCESS:MINC2_ERROR;
    
    
  } else {
    return MINC2_ERROR;
  }
}



/**
* minc2 file iterator functions
*/

struct minc2_file_iterator
{
  minc2_file_handle *_minc_file; /*< associated minc file*/
  int _ndim;                    /*< number of dimesnions (copy of _minc_file->ndims*/
  int _fnum;
 
  int *_index;                  /*< current voxel index */
  int *_start;                  /*< ROI start index*/
  int *_end;                    /*< ROI end index*/
  int *_count;                  /*< slice size index*/

  int _output_mode;             /*< are we dealing with output or input iterator */
  int _data_type;               /*< iterator data type */
                 
  int _slice_dimensions;        /*< number of slice dimensions */

  void *_buffer;                /*< slice buffer*/

  int _buffer_size;             /*< slice buffer size in elements */
  int _buffer_index;            /*< current buffer index */
  int _element_size;            /*< size of one element */

  double *_min;                 /*buffer to keep volume range*/
  double *_max;                 /*buffer to keep volume range*/
};


minc2_file_iterator_handle minc2_iterator_allocate0(void)
{
  minc2_file_iterator_handle h=calloc(1,sizeof(struct minc2_file_iterator));
  return h;
}

int minc2_iterator_free(minc2_file_iterator_handle h)
{
  if(!h) return MINC2_SUCCESS;

  if(h->_minc_file) free(h->_minc_file);
  if(h->_index) free(h->_index);
  if(h->_start) free(h->_start);
  if(h->_end)   free(h->_end);
  if(h->_buffer)free(h->_buffer);

  if(h->_min)   free(h->_min);
  if(h->_max)   free(h->_max);

  free(h);

  return MINC2_SUCCESS;
}

static int minc2_iterator_start(minc2_file_iterator_handle h,minc2_file_handle *m,int data_type,int fnumber)
{
  int i;
  h->_data_type=data_type;
  h->_fnum=fnumber;

  h->_minc_file=(minc2_file_handle*)realloc(h->_minc_file,sizeof(minc2_file_handle)*fnumber);

  memcpy(h->_minc_file,m,sizeof(minc2_file_handle)*fnumber);

  for(i=1;i<fnumber;i++)
  {
    if(minc2_compare_voxel_dimensions(h->_minc_file[0]->representation_dims,h->_minc_file[i]->representation_dims)!=MINC2_SUCCESS)
    {
      MI_LOG_ERROR(MI2_MSG_GENERIC,"Iterator: files have incompatible dimensions: %d and %d",0,i);
      return MINC2_ERROR;
    }
  }

  h->_ndim=h->_minc_file[0]->ndims;

  h->_index=(int*)realloc(h->_index,sizeof(int)*h->_ndim);
  h->_start=(int*)realloc(h->_start,sizeof(int)*h->_ndim);
  h->_end=  (int*)realloc(h->_end,  sizeof(int)*h->_ndim);
  h->_count=(int*)realloc(h->_count,sizeof(int)*h->_ndim);

  /**/
  h->_min=(double*)realloc(h->_min,sizeof(double)*fnumber);
  h->_max=(double*)realloc(h->_max,sizeof(double)*fnumber);

  /*for now use the whole volume */
  for ( i = 0; i < h->_ndim ; i++ )
  {
    h->_start[i]=0;
    h->_count[i]=1;
    h->_index[i]=0;

    if(h->_minc_file[0]->using_apparent_order)
      h->_end[i]=h->_minc_file[0]->representation_dims[i].length;
    else
      h->_end[i]=h->_minc_file[0]->store_dims[i].length;
   
  }

  for(i=0;i<fnumber;i++)
  {
    h->_min[i]=DBL_MAX;
    h->_max[i]=-DBL_MAX;
  }

  minc2_slice_ndim(h->_minc_file[0],&h->_slice_dimensions);
  if(h->_slice_dimensions<1)
    h->_slice_dimensions=1;

  h->_buffer_size=1;
  for(i=0;i<h->_slice_dimensions;i++)
  {
    h->_count[i]=h->_end[i]-h->_start[i];
    h->_buffer_size*=h->_count[i];
  }
  
  h->_element_size=_minc2_type_size(data_type);
  h->_buffer=(void*)realloc(h->_buffer,h->_element_size*h->_buffer_size*h->_fnum);
  h->_buffer_index=0;
  
  return MINC2_SUCCESS;
}

static int minc2_iterator_flush(minc2_file_iterator_handle h)
{
  int i,f;
  int r=1;

  for(f=0;f< h->_fnum; f++)
  {
    void *f_buffer=h->_buffer + h->_buffer_size*h->_element_size*f;

    if(h->_output_mode)
    {
      double buffer_min,buffer_max;
      switch(h->_data_type )
      {
        case MINC2_UBYTE:
          _GET_BUFFER_MIN_MAX(unsigned char,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_BYTE:
          _GET_BUFFER_MIN_MAX(char,f_buffer, h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_USHORT:
          _GET_BUFFER_MIN_MAX(unsigned short,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_SHORT:
          _GET_BUFFER_MIN_MAX(short,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_UINT:
          _GET_BUFFER_MIN_MAX(unsigned int,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_INT:
          _GET_BUFFER_MIN_MAX(int,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_FLOAT:
          _GET_BUFFER_MIN_MAX_PROT(float,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        case MINC2_DOUBLE:
          _GET_BUFFER_MIN_MAX_PROT(double,f_buffer,h->_buffer_size,buffer_min,buffer_max);
          break;
        default:
          MI_LOG_ERROR(MI2_MSG_GENERIC,"Unsupported volume data type");
          return MINC2_ERROR;
      }

      if( h->_minc_file[f]->slice_scaling_flag )
        if ( minc2_set_slice_range(h->_minc_file[f],h->_index,buffer_min,buffer_max)!=MINC2_SUCCESS )
          return MINC2_ERROR;

      r=minc2_write_hyperslab(h->_minc_file[f], h->_index,h->_count, f_buffer, h->_data_type)==MINC2_SUCCESS&&r;

      if(isfinite(buffer_min) && buffer_min<h->_min[f]) h->_min[f]=buffer_min;
      if(isfinite(buffer_max) && buffer_max>h->_max[f]) h->_max[f]=buffer_max;

    } else {
      r=minc2_read_hyperslab( h->_minc_file[f], h->_index,h->_count, f_buffer, h->_data_type)==MINC2_SUCCESS&&r;
    }
  }
  return r?MINC2_SUCCESS:MINC2_ERROR;
}


int minc2_iterator_input_start(minc2_file_iterator_handle h,minc2_file_handle m,int data_type)
{
  if(minc2_iterator_start(h,&m,data_type,1)!=MINC2_SUCCESS)
    return MINC2_ERROR;
  
  h->_output_mode=0;
  
  return minc2_iterator_flush(h); /*read first slice*/
  
}

int minc2_iterator_output_start(minc2_file_iterator_handle h,minc2_file_handle m,int data_type)
{
  if(minc2_iterator_start(h,&m,data_type,1)!=MINC2_SUCCESS)
    return MINC2_ERROR;

  h->_output_mode=1;
  return MINC2_SUCCESS;
}

int minc2_multi_iterator_input_start (minc2_file_iterator_handle h,minc2_file_handle *m,int data_type,int fnum)
{
  if(minc2_iterator_start(h,m,data_type,fnum)!=MINC2_SUCCESS)
    return MINC2_ERROR;

  h->_output_mode=0;

  return minc2_iterator_flush(h); /*read first slice*/
}

int minc2_multi_iterator_output_start(minc2_file_iterator_handle h,minc2_file_handle *m,int data_type,int fnum)
{
  if(minc2_iterator_start(h,m,data_type,fnum)!=MINC2_SUCCESS)
    return MINC2_ERROR;

  h->_output_mode=1;
  return MINC2_SUCCESS;
}

int minc2_iterator_next(minc2_file_iterator_handle h)
{
  int i=0;
  h->_buffer_index++;
  do {
    if( i == h->_slice_dimensions && h->_output_mode ) /*write last slice*/
    {
      int f;
      int r=1;
      h->_buffer_index=0;
      
      if(minc2_iterator_flush(h)!=MINC2_SUCCESS) 
        return MINC2_ERROR; /*I/O error, report somehow?*/
      
      for(f=0;f<h->_fnum;f++)
        r = r && ( minc2_set_volume_range(h->_minc_file[f], h->_min[f], h->_max[f] )==MINC2_SUCCESS);

      if(!r)
        return MINC2_ERROR;
    }
    
    h->_index[i]++;
    
    if(h->_index[i] < h->_end[i]) 
    {
      if( i == h->_slice_dimensions && !h->_output_mode )  /*read next slice*/
      {
        h->_buffer_index=0;
      
        if(minc2_iterator_flush(h)!=MINC2_SUCCESS) 
          return MINC2_ERROR; /*I/O error, report somehow?*/
      }
      return MINC2_SUCCESS;
      
    } else {
      if( i < (h-> _ndim-1) )
      {
        h->_index[i]=h->_start[i];
        i++;
      } else {
        return MINC2_ERROR; /*EOF*/
      }
    }
  } while( 1 );
  return MINC2_ERROR; /*should never get here*/
}

int minc2_iterator_get_values(minc2_file_iterator_handle h,void *val)
{
  int f;
  const void *_buffer=h->_buffer+h->_buffer_index*h->_element_size;
  void *_val=val;

  for(f=0;f<h->_fnum;f++)
  {
    memcpy(_val,_buffer,h->_element_size);
    _buffer+=h->_buffer_size*h->_element_size;
    _val+=h->_element_size;
  }
  return MINC2_SUCCESS;
}

int minc2_iterator_put_values(minc2_file_iterator_handle h,const void *val)
{
  int f;
  void *_buffer=h->_buffer+h->_buffer_index*h->_element_size;
  const void *_val=val;

  for(f=0;f<h->_fnum;f++)
  {
    memcpy(_buffer,_val,h->_element_size);
    _buffer+=h->_buffer_size*h->_element_size;
    _val+=h->_element_size;
  }
  return MINC2_SUCCESS;
}


/**
 * tag operations
 */
minc2_tags_handle  minc2_tags_allocate0(void)
{
  minc2_tags_handle m=calloc(1,sizeof(struct minc2_tags));

  return m;
}

int minc2_tags_free(minc2_tags_handle tags)
{
  int i,n;
  if(!tags) return MINC2_ERROR;

  if(tags->tags_volume1)
    free(tags->tags_volume1);

  if(tags->tags_volume2)
    free(tags->tags_volume2);

  if(tags->weights)
    free(tags->weights);

  if(tags->structure_ids)
    free(tags->structure_ids);

  if(tags->patient_ids)
    free(tags->patient_ids);

  if(tags->labels)
  {
    for(i=0;i<tags->n_tag_points;i++)
    {
      if(tags->labels[i])
      {
         free(tags->labels[i]);
      }
    }
    free(tags->labels);
  }
  free(tags);
  return MINC2_SUCCESS;
}


static void _convert_array_from_VIO(double *dst, VIO_Real  **src, int n_tag_points)
{
  int i=0;
  for(i=0;i<n_tag_points;i++)
  {
    *dst++=src[i][0];
    *dst++=src[i][1];
    *dst++=src[i][2];
  }
}

static void _convert_array_to_VIO(double *src, VIO_Real  **dst, int n_tag_points)
{
  int i=0;
  for(i=0;i<n_tag_points;i++)
  {
    dst[i][0]=*src++;
    dst[i][1]=*src++;
    dst[i][2]=*src++;
  }
}

static VIO_Real ** _allocate_vio_vectors(int n,int m)
{
  VIO_Real ** vol=NULL;
  int i;
  SET_ARRAY_SIZE( vol, 0, n, 10 );
  for(i=0;i<n;i++)
  {
    ALLOC( vol[i], m );
  }
  return vol;
}

static int minc2_tags_convert_from_VIO(minc2_tags_handle tags,
                    int n_volumes,
                    int n_tag_points,
                    VIO_Real  **tags_volume1,
                    VIO_Real  **tags_volume2,
                    VIO_Real  *weights,
                    int       *structure_ids,
                    int       *patient_ids,
                    VIO_STR    labels[]
                   )
{
  tags->n_volumes=n_volumes;
  tags->n_tag_points=n_tag_points;

  tags->tags_volume1=malloc(n_tag_points*3*sizeof(double));
  _convert_array_from_VIO(tags->tags_volume1, tags_volume1, n_tag_points);

  if(n_volumes>1 && tags_volume2){
    tags->tags_volume2=malloc(n_tag_points*3*sizeof(double));
    _convert_array_from_VIO(tags->tags_volume2, tags_volume2, n_tag_points);
  } else {
    tags->tags_volume2=NULL;
  }

  if(weights) {
    int i;
    tags->weights=malloc(n_tag_points*sizeof(double));
    for(i=0;i<n_tag_points;i++) tags->weights[i]=weights[i];
  } else {
    tags->weights=NULL;
  }

  if(structure_ids) {
    tags->structure_ids=malloc(n_tag_points*sizeof(int));
    memmove(tags->structure_ids, structure_ids, sizeof(int)*n_tag_points);
  } else {
    tags->structure_ids=NULL;
  }

  if(patient_ids) {
    tags->patient_ids=calloc(n_tag_points,sizeof(int));
    memmove(tags->patient_ids,patient_ids,sizeof(int)*n_tag_points);
  } else {
    tags->patient_ids=NULL;
  }

  if(labels) {
    int i;
    tags->labels=malloc(n_tag_points*sizeof(const char *));
    for(i=0;i<n_tag_points;i++) tags->labels[i]=strdup(labels[i]);
  } else {
    tags->labels=NULL;
  }

  /*TODO: check if all memoru properly allocated*/
  return MINC2_SUCCESS;
}

static int minc2_tags_convert_to_VIO(minc2_tags_handle tags,
                    int *n_volumes,
                    int *n_tag_points,
                    VIO_Real  ***tags_volume1,
                    VIO_Real  ***tags_volume2,
                    VIO_Real  **weights,
                    int       **structure_ids,
                    int       **patient_ids,
                    VIO_STR   *labels[]
                   )
{

  *n_volumes=tags->n_volumes;
  *n_tag_points=tags->n_tag_points;

  //tags->tags_volume1=malloc( tags->n_tag_points*3*sizeof(VIO_Real) );
  *tags_volume1=_allocate_vio_vectors(tags->n_tag_points,3);
  _convert_array_to_VIO(tags->tags_volume1, *tags_volume1, tags->n_tag_points);

  if(tags->n_volumes>1 && tags->tags_volume2){
    *tags_volume2=_allocate_vio_vectors(tags->n_tag_points,3);
    _convert_array_to_VIO(tags->tags_volume2, *tags_volume2, tags->n_tag_points);
  } else {
    *tags_volume2=NULL;
  }

  if(tags->weights) {
    int i;
    *weights=malloc(tags->n_tag_points*sizeof(VIO_Real));
    for(i=0;i<tags->n_tag_points;i++) (*weights)[i]=tags->weights[i];
  } else {
    *weights=NULL;
  }

  if(tags->structure_ids) {
    *structure_ids=malloc(tags->n_tag_points*sizeof(int));
    memmove(*structure_ids,tags->structure_ids,sizeof(int)*tags->n_tag_points);
  } else {
    *structure_ids=NULL;
  }

  if(tags->patient_ids) {
    *patient_ids=malloc(tags->n_tag_points*sizeof(int));
    memmove(*patient_ids,tags->patient_ids,sizeof(int)*tags->n_tag_points);
  } else {
    *patient_ids=NULL;
  }

  if(tags->labels) {
    int i;
    *labels=malloc(tags->n_tag_points*sizeof(const char *));
    for(i=0;i< tags->n_tag_points ;i++) {
      if(tags->labels[i]) (*labels)[i]=strdup(tags->labels[i]);
      else  (*labels)[i]=strdup("");
    }
  } else {
    *labels=NULL;
  }

  /*TODO: check if all memoru properly allocated*/
  return MINC2_SUCCESS;
}

int minc2_tags_load(minc2_tags_handle tags,const char *file)
{
  int       n_volumes;
  int       n_tag_points;
  VIO_Real      **tags_volume1;
  VIO_Real      **tags_volume2;
  VIO_Real      *weights;
  int       *structure_ids;
  int       *patient_ids;
  VIO_STR   *labels;
  int i;
  int ret=MINC2_ERROR;

  if ( input_tag_file((char *)file, &n_volumes, &n_tag_points,
       &tags_volume1, &tags_volume2, &weights, &structure_ids, &patient_ids, &labels ) != VIO_OK ) {
    return MINC2_ERROR;
  }
  ret=minc2_tags_convert_from_VIO(tags,n_volumes,n_tag_points,tags_volume1,tags_volume2,weights,structure_ids,patient_ids,labels);

  free_tag_points(
    n_volumes, n_tag_points,
    tags_volume1, tags_volume2,
    weights,structure_ids,
    patient_ids,labels );

  return ret;
}

int minc2_tags_save(minc2_tags_handle tags,const char *file)
{
  int         n_volumes;
  int         n_tag_points;
  VIO_Real  **tags_volume1;
  VIO_Real  **tags_volume2;
  VIO_Real   *weights;
  int         *structure_ids;
  int         *patient_ids;
  VIO_STR     *labels;

  int ret=MINC2_ERROR;

  ret=minc2_tags_convert_to_VIO(tags,&n_volumes, &n_tag_points, &tags_volume1, &tags_volume2,&weights,&structure_ids,&patient_ids,&labels);

  if ( output_tag_file((char *)file, "minc2-simple tag output", n_volumes, n_tag_points,
       tags_volume1, tags_volume2, weights, structure_ids, patient_ids, labels ) != VIO_OK ) {
    ret=MINC2_ERROR;
  }

  free_tag_points(
    n_volumes, n_tag_points,
    tags_volume1, tags_volume2,
    weights,structure_ids,
    patient_ids,labels );

  return ret;
}

int minc2_tags_init(minc2_tags_handle tags, int n_tag_points, int n_volumes, int have_weights, int have_strucure_ids, int have_patient_ids, int have_labels)
{
  tags->n_volumes=n_volumes;
  tags->n_tag_points=n_tag_points;

  tags->tags_volume1=malloc(n_tag_points*3*sizeof(double));

  /*TODO: check if the memory was already allocated!*/
  if(n_volumes>1){
    tags->tags_volume2=malloc(n_tag_points*3*sizeof(double));
  }

  if(have_weights) {
    tags->weights=malloc(n_tag_points*sizeof(double));
  }

  if(have_strucure_ids) {
    tags->structure_ids=calloc(n_tag_points,sizeof(int));
  }

  if(have_patient_ids) {
    tags->patient_ids=calloc(n_tag_points,sizeof(int));
  }

  if(have_labels) {
    tags->labels=calloc(n_tag_points,sizeof(const char *));
  }

  /*TODO: check memory allocation!*/
  return MINC2_SUCCESS;
}

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; remove-trailing-spaces modified; hl c*/
