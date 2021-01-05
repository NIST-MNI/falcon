/**
  * minc2 dimension types
  */ 
enum  minc2_dimensions {
  MINC2_DIM_UNKNOWN=0,
  MINC2_DIM_X,
  MINC2_DIM_Y,
  MINC2_DIM_Z,
  MINC2_DIM_TIME,
  MINC2_DIM_VEC,
  MINC2_DIM_END=255
};


/**
 * minc2 data storage type
 * compatible with minc2 API
 */
enum  minc2_type {
  MINC2_ORIGINAL = 0,     /**< MI_ORIGINAL_TYPE */
  MINC2_BYTE = 1,         /**< 8-bit signed integer */
  MINC2_SHORT = 3,        /**< 16-bit signed integer */
  MINC2_INT = 4,          /**< 32-bit signed integer */
  MINC2_FLOAT = 5,        /**< 32-bit floating point */
  MINC2_DOUBLE = 6,       /**< 64-bit floating point */
  MINC2_STRING = 7,       /**< ASCII string (?) */
  MINC2_UBYTE = 100,      /**< 8-bit unsigned integer */
  MINC2_USHORT = 101,     /**< 16-bit unsigned integer */
  MINC2_UINT = 102,       /**< 32-bit unsigned integer */
  MINC2_SCOMPLEX = 1000,  /**< 16-bit signed integer complex */
  MINC2_ICOMPLEX = 1001,  /**< 32-bit signed integer complex */
  MINC2_FCOMPLEX = 1002,  /**< 32-bit floating point complex */
  MINC2_DCOMPLEX = 1003,  /**< 64-bit floating point complex */
  MINC2_MAX_TYPE_ID,
  MINC2_UNKNOWN  = -1     /**< when the type is a record */
};


/**
  * XFM type
  */
enum  minc2_xfm {
  MINC2_XFM_LINEAR=1,
  MINC2_XFM_THIN_PLATE_SPLINE,
  MINC2_XFM_USER_TRANSFORM,
  MINC2_XFM_CONCATENATED_TRANSFORM,
  MINC2_XFM_GRID_TRANSFORM,
  MINC2_XFM_END
};


/**
 * minc2 dimension information
 */
struct minc2_dimension
{
  int    id;             /**< dimension id, see enum minc2_dimensions*/
  int    length;         /**< dimension length */
  int    irregular;      /**< flag to show irregular sampling */
  double step;           /**< dimension step   */
  double start;          /**< dimension start  */
  int    have_dir_cos;   /**< flag that dimension cosines is valid*/
  double dir_cos[3];     /**< direction cosines*/
};

/**
 *
 */
struct minc2_info_iterator;
typedef struct minc2_info_iterator* minc2_info_iterator_handle;


/**
 *
 */
struct minc2_file_iterator;
typedef struct minc2_file_iterator* minc2_file_iterator_handle;


/**
 * minc2 error codes, compatible with minc2 API
 */
enum { MINC2_SUCCESS=0,MINC2_ERROR=-1};

/** Opaque structure representing minc2 file
 * 
 */
struct minc2_file;
typedef struct minc2_file* minc2_file_handle; 

/** Opaque structure representing minc2 XFM file
 *
 */
struct minc2_xfm_file;
typedef struct minc2_xfm_file *minc2_xfm_file_handle;


/**
 * minc tag file info
 */
struct minc2_tags
{
  int         n_volumes;
  int         n_tag_points;

  double     *tags_volume1;
  double     *tags_volume2;

  double     *weights;
  int        *structure_ids;
  int        *patient_ids;

  char      **labels;
};

typedef struct minc2_tags *minc2_tags_handle;


/** 
 * allocate empty minc2 file structure, no need to call minc2_init after
 */
int minc2_allocate(minc2_file_handle * h);

/**
 * alternative version
 */
minc2_file_handle minc2_allocate0(void);


/** 
 * initialize minc2 file structure
 */
int minc2_init(minc2_file_handle h);


/**
 * deallocate minc2 file structure
 * will call standard free on it
 */
int minc2_free(minc2_file_handle h);


/**
 * close minc2 file if it's open,
 * then deallocate minc2 file structure
 */
int minc2_destroy(minc2_file_handle h);

/**
 * open existing file
 */
int minc2_open(minc2_file_handle h,const char * path);


/**
 * define a new minc2 volume, using provided storage dimension information and storage data type
 */
int minc2_define(minc2_file_handle h, struct minc2_dimension *store_dims, int store_data_type,int data_type);

/**
 * create a new file, owerwriting an existing one if needed
 */
int minc2_create(minc2_file_handle h,const char * path);

/**
 * close file, flushing data to disk
 * sets minccomplete flag too
 */
int minc2_close(minc2_file_handle h);

/**
 * query number of dimensions
 */
int minc2_ndim(minc2_file_handle h,int *ndim);

/**
 * query total number of voxels
 */
int minc2_nelement(minc2_file_handle h,int *nelement);

/**
 * query data type, used to represent data
 */
int minc2_data_type(minc2_file_handle h,int *_type);

/**
 * query data type, used to store data on disk
 */
int minc2_storage_data_type(minc2_file_handle h,int *_type);

/**
 * query number of slice dimensions 
 */
int minc2_slice_ndim(minc2_file_handle h,int *slice_ndim);

/**
 * Setup minc file for reading or writing information
 * in standardized order ( Vector dimension - X - Y -Z -TIME )
 * with positive steps
 */
int minc2_setup_standard_order(minc2_file_handle h);

/**
 * get dimension information in current representation format
 */
int minc2_get_representation_dimensions(minc2_file_handle h,struct minc2_dimension **dims);

/**
 * get dimension information in file format
 */
int minc2_get_store_dimensions(minc2_file_handle h,struct minc2_dimension **dims);

/**
 * Compare if volumes have compatible storage space (ignore spacing, direction cosines, etc)
 */
int minc2_compare_voxel_dimensions(const struct minc2_dimension *one,const struct minc2_dimension *two);

/**
 * Compare if volumes have compatible everything
 */
int minc2_compare_dimensions(const struct minc2_dimension *one,const struct minc2_dimension *two);


/**
 * Load complete volume into memory
 */
int minc2_load_complete_volume(minc2_file_handle h,void *buffer,int representation_type);

/**
 * Save complete volume into memory
 */
int minc2_save_complete_volume(minc2_file_handle h,const void *buffer,int representation_type);

/**
 * Specify flags to use scaling
 * this have to be set before minc2_create
 */
int minc2_set_scaling(minc2_file_handle h,int use_global_scaling,int use_slice_scaling);

/**
 * Specify volume range, only when using hyperslab writing
 * Implies no slice scaling 
 */
int minc2_set_volume_range(minc2_file_handle h,double value_min,double value_max);

/**
 * Specify slice range, only when using hyperslab writing with slice scaling
 * Implies no slice scaling
 */
int minc2_set_slice_range(minc2_file_handle h,int *start,double value_min,double value_max);

/**
 * convert world X,Y,Z coordinates to voxel indexes (also in X,Y,Z order)
 */
int minc2_world_to_voxel(minc2_file_handle h,const double *world,double *voxel);


/**
 * convert world X,Y,Z coordinates to voxel indexes (also in X,Y,Z order)
 * vectorized version
 */
int minc2_world_to_voxel_vec(minc2_file_handle h, int n, int stride, const double *world,double *voxel);


/**
 * convert voxel X,Y,Z indexes to world coordinates (also in X,Y,Z order)
 */
int minc2_voxel_to_world(minc2_file_handle h,const double *voxel,double *world);


/**
 * convert voxel X,Y,Z indexes to world coordinates (also in X,Y,Z order)
 * vectorized version
 */
int minc2_voxel_to_world_vec(minc2_file_handle h,int n, int stride, const double *voxel,double *world);


/**
 * transfer attributes from one volume to another
 */
int minc2_copy_metadata(minc2_file_handle src,minc2_file_handle dst);

/**
 * write hyperslab, using current dimension order
 * WARNING: no range checks performed!
 */
int minc2_write_hyperslab(minc2_file_handle h,int *start,int *count,const void* buffer,int representation_type);

/**
 * read hyperslab, using current dimension order
 * WARNING: no range checks performed!
 */
int minc2_read_hyperslab(minc2_file_handle h,int *start,int *count,void* buffer,int representation_type);

/**
 * return human-readable type name
 */
const char * minc2_data_type_name(int minc2_type_id);

/**
 * return human-readable dimension name
 */
const char * minc2_dim_type_name(int minc2_dim_id);


/**
 * get attribute type
 */
int minc2_get_attribute_type(minc2_file_handle h,const char* group,const char* attr,int *minc2_type);

/**
 * get attribute length
 */
int minc2_get_attribute_length(minc2_file_handle h,const char* group,const char* attr,int *attr_length);

/**
 * read attribute
 */
int minc2_read_attribute(minc2_file_handle h,const char* group,const char* attr,void *buf,int buf_size);

/**
 * write attribute
 */
int minc2_write_attribute(minc2_file_handle h,const char* group,const char* attr,const void *buf,int buf_size,int minc2_type);

/**
 * delete attribute
 */
int minc2_delete_attribute(minc2_file_handle h,const char* group,const char* attr);

/**
 * delete the whole group
 */
int minc2_delete_group(minc2_file_handle h,const char* group);

/**
 * Ititialize info iterator
 */
minc2_info_iterator_handle minc2_allocate_info_iterator(void);

/**
 * Free info iterator: stop iterator if needed and deallocate memory
 */
int minc2_free_info_iterator(minc2_info_iterator_handle it);

/**
 * Stop iterator: stop itarating , the iterator handle can be re-used for another time
 */
int minc2_stop_info_iterator(minc2_info_iterator_handle it);


/**
 * Start iterating over groups
 */
int minc2_start_group_iterator(minc2_file_handle h,minc2_info_iterator_handle group_it);

/**
 * Start iterating over attrobutes
 */
int minc2_start_attribute_iterator(minc2_file_handle h,const char* group,minc2_info_iterator_handle it);

/**
 * advance to next available group item
 */
int minc2_iterator_group_next(minc2_info_iterator_handle it);

/**
 * advance to next available attrobute item
 */
int minc2_iterator_attribute_next(minc2_info_iterator_handle it);

/**
 * Get current iterator contents
 */
const char* minc2_iterator_group_name(minc2_info_iterator_handle it);

/**
 * Get current iterator contents
 */
const char* minc2_iterator_attribute_name(minc2_info_iterator_handle it);

/**
 * generate timestamp
 */
char* minc2_timestamp(int argc,char **argv);

/**
 * allocate empty minc2 xfm file structure
 */
int minc2_xfm_allocate(minc2_xfm_file_handle * h);

/**
 * alternative version
 */
minc2_xfm_file_handle minc2_xfm_allocate0(void);

/**
 * initialize minc2 xfm file structure
 */
int minc2_xfm_init(minc2_xfm_file_handle h);

/**
 * deallocate minc2 xfm file structure
 * will call standard free on it
 */
int minc2_xfm_free(minc2_xfm_file_handle h);

/**
 * close minc2 xfm file if it's open,
 * then deallocate minc2 file structure
 */
int minc2_xfm_destroy(minc2_xfm_file_handle h);

/**
 * open existing XFM file
 */
int minc2_xfm_open(minc2_xfm_file_handle h,const char * path);

/**
 * save XFM  to file
 */
int minc2_xfm_save(minc2_xfm_file_handle h,const char * path);

/**
 * transform x,y,z coordinates
 */
int minc2_xfm_transform_point(minc2_xfm_file_handle h,const double* in,double* out);

/**
 * invert transform x,y,z coordinates
 */
int minc2_xfm_inverse_transform_point(minc2_xfm_file_handle h,const double* in,double* out);


/**
 * transform x,y,z coordinates for n points separated by stride
 */
int minc2_xfm_transform_point_vec(minc2_xfm_file_handle h,int n,int stride,const double* in,double* out);

/**
 * invert transform x,y,z coordinates for n points separated by stride
 */
int minc2_xfm_inverse_transform_point_vec(minc2_xfm_file_handle h,int n,int stride,const double* in,double* out);



/**
 * set flag to invert transform
 */
int minc2_xfm_invert(minc2_xfm_file_handle h);

/**
 * get number of concatenated transforms, return at least 1
 */
int minc2_xfm_get_n_concat(minc2_xfm_file_handle h,int *n);

/**
 * get type of nth transform
 */
int minc2_xfm_get_n_type(minc2_xfm_file_handle h,int n,int *xfm_type);

/**
 * extract n'th transform, if it is linear, as a 4x4 matrix , in a row-major fashion
 */
int minc2_xfm_get_linear_transform(minc2_xfm_file_handle h,int n,double *matrix);

/**
 * extract n'th transform, if it is nonlinear, as a reference to a grid file
 */
int minc2_xfm_get_grid_transform(minc2_xfm_file_handle h,int n,int *inverted,char **grid_file);

/**
 * Adds another  transform, a 4x4 matrix , in a row-major fashion
 */
int minc2_xfm_append_linear_transform(minc2_xfm_file_handle h,double *matrix);

/**
 * Adds another nonlinear grid transform
 */
int minc2_xfm_append_grid_transform(minc2_xfm_file_handle h,const char * grid_path,int inv);

/**
 * Concatenate another general xfm transform
 */
int minc2_xfm_concat_xfm(minc2_xfm_file_handle h,minc2_xfm_file_handle o);


/**
 * Generate linear transform based on parameters and append it
 */
int minc2_xfm_append_linear_param(minc2_xfm_file_handle h,
                              double *center,
                              double *translations,
                              double *scales,
                              double *shears,
                              double *rotations);

/**
* Extract linear parameters from the transform
* center have to be specied
*/
int minc2_xfm_extract_linear_param(minc2_xfm_file_handle h,
                             int n,
                             double *center,
                             double *translations,
                             double *scales,
                             double *shears,
                             double *rotations);


/**
 * Iterators
 */
minc2_file_iterator_handle minc2_iterator_allocate0(void);
int minc2_iterator_free(minc2_file_iterator_handle h);

int minc2_iterator_input_start (minc2_file_iterator_handle h,minc2_file_handle m,int data_type);
int minc2_iterator_output_start(minc2_file_iterator_handle h,minc2_file_handle m,int data_type);


int minc2_multi_iterator_input_start (minc2_file_iterator_handle h,minc2_file_handle *m,int data_type,int fnum);
int minc2_multi_iterator_output_start(minc2_file_iterator_handle h,minc2_file_handle *m,int data_type,int fnum);


int minc2_iterator_next(minc2_file_iterator_handle h);
/*int minc2_iterator_flush(minc2_file_iterator_handle h);*/

int minc2_iterator_get_values(minc2_file_iterator_handle h,void *val);
int minc2_iterator_put_values(minc2_file_iterator_handle h,const void *val);


/**
 * Tags io
 */

/**
 * allocate tag structure (empty)
 */
minc2_tags_handle minc2_tags_allocate0(void);

/**
 * free all memory resources
 * used in tags
 */
int minc2_tags_free(minc2_tags_handle tags);

/**
 * load tags from file
 */
int minc2_tags_load(minc2_tags_handle tags,const char *file);

/**
 * save tags to file
 */
int minc2_tags_save(minc2_tags_handle tags,const char *file);

/**
 * initialize minc tags with fixed number of tags
 */
int minc2_tags_init(minc2_tags_handle tags,int n_tag_points,int n_volumes,int have_weights,int have_strucure_ids,int have_patient_ids,int have_labels);

/**
 *
 */

/* kate: indent-mode cstyle; indent-width 2; replace-tabs on; remove-trailing-spaces modified; hl c*/
