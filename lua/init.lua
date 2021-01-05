-- lua module to read and write minc2 files
-- using minc2-simple c glue
-- using FFI 
local ffi = require("ffi")
require('torch')

-- contents of ../src/minc2-simple.h :
ffi.cdef[[
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
 * convert voxel X,Y,Z indexes to world coordinates (also in X,Y,Z order)
 */
int minc2_voxel_to_world(minc2_file_handle h,const double *voxel,double *world);


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

]]

local lib = ffi.load("minc2-simple") -- for now fixed path

minc2_file = {
    -- minc2 constants
    
    -- minc2 dimensions
    MINC2_DIM_UNKNOWN=ffi.C.MINC2_DIM_UNKNOWN,
    MINC2_DIM_X    = ffi.C.MINC2_DIM_X,
    MINC2_DIM_Y    = ffi.C.MINC2_DIM_Y,
    MINC2_DIM_Z    = ffi.C.MINC2_DIM_Z,
    MINC2_DIM_TIME = ffi.C.MINC2_DIM_TIME,
    MINC2_DIM_VEC  = ffi.C.MINC2_DIM_VEC,
    MINC2_DIM_END  = ffi.C.MINC2_DIM_END,
    
    -- minc2 data types
    MINC2_BYTE     = ffi.C.MINC2_BYTE ,
    MINC2_SHORT    = ffi.C.MINC2_SHORT ,
    MINC2_INT      = ffi.C.MINC2_INT ,
    MINC2_FLOAT    = ffi.C.MINC2_FLOAT ,
    MINC2_DOUBLE   = ffi.C.MINC2_DOUBLE ,
    MINC2_STRING   = ffi.C.MINC2_STRING ,
    MINC2_UBYTE    = ffi.C.MINC2_UBYTE ,
    MINC2_USHORT   = ffi.C.MINC2_USHORT ,
    MINC2_UINT     = ffi.C.MINC2_UINT ,
    MINC2_SCOMPLEX = ffi.C.MINC2_SCOMPLEX ,
    MINC2_ICOMPLEX = ffi.C.MINC2_ICOMPLEX ,
    MINC2_FCOMPLEX = ffi.C.MINC2_FCOMPLEX ,
    MINC2_DCOMPLEX = ffi.C.MINC2_DCOMPLEX ,
    MINC2_MAX_TYPE_ID=ffi.C.MINC2_MAX_TYPE_ID,
    MINC2_UNKNOWN  = ffi.C.MINC2_UNKNOWN  ,

    -- minc2 status
    MINC2_SUCCESS  = ffi.C.MINC2_SUCCESS,
    MINC2_ERROR    = ffi.C.MINC2_ERROR,

    -- XFM types
    MINC2_XFM_LINEAR=ffi.C.MINC2_XFM_LINEAR,
    MINC2_XFM_THIN_PLATE_SPLINE=ffi.MINC2_XFM_THIN_PLATE_SPLINE,
    MINC2_XFM_USER_TRANSFORM=ffi.MINC2_XFM_USER_TRANSFORM,
    MINC2_XFM_CONCATENATED_TRANSFORM=ffi.MINC2_XFM_CONCATENATED_TRANSFORM,
    MINC2_XFM_GRID_TRANSFORM=ffi.MINC2_XFM_GRID_TRANSFORM
}
minc2_file.__index = minc2_file

function minc2_file.new(path)
  local self = setmetatable({}, minc2_file)
  self._v=ffi.gc(lib.minc2_allocate0(),lib.minc2_destroy)
  if path~=nil then
      self:open(path)
  end
  return self
end

-- open existing minc2 file
function minc2_file:open(path)
    --print("Going to open:"..path)
    assert(path~=nil,"Provide minc2 file")
    assert( lib.minc2_open(self._v,path) == ffi.C.MINC2_SUCCESS )
end

-- close a minc2 file
function minc2_file:close()
    assert(lib.minc2_close(self._v)==ffi.C.MINC2_SUCCESS)
end

-- query number of dimensions
function minc2_file:ndim()
    local dd=ffi.new("int[1]")
    assert(lib.minc2_ndim(self._v,dd)==ffi.C.MINC2_SUCCESS)
    return dd[0]
end

-- provide descriptor of dimensions
function minc2_file:store_dims()
    local dims=ffi.new("struct minc2_dimension*[1]")
    assert(lib.minc2_get_store_dimensions(self._v,dims)==ffi.C.MINC2_SUCCESS)
    return dims[0]
end

-- provide descriptor of dimensions
function minc2_file:representation_dims()
    local dims=ffi.new("struct minc2_dimension*[1]")
    assert(lib.minc2_get_representation_dimensions(self._v,dims)==ffi.C.MINC2_SUCCESS)
    return dims[0]
end

-- provide volume size, using current representation
function minc2_file:volume_size()

    local _dims=self:representation_dims()
    local ndims=self:ndim()

    local sz=torch.LongStorage(ndims)
    local i

    for i=1,ndims do
        sz[i]=_dims[ndims-i].length
    end
    return sz
end

-- define a new volume
function minc2_file:define(dims, store_type, representation_type)
    assert(dims~=nil,"dims need to be defined")
    assert(store_type~=nil,"Store data type need to be set")
    assert(representation_type~=nil,"Data type need to be set")
    
    if type(dims)== "table" then
        --assume user didn't provide m2.MINC2_DIM_END
        local mydims={}
        for k, v in pairs(dims) do mydims[k] = v end
        mydims[#mydims+1]={id=minc2_file.MINC2_DIM_END }
        
        dims=ffi.new("struct minc2_dimension[?]",#mydims,mydims)
    end
    
    assert(lib.minc2_define(self._v,dims,store_type,representation_type)==ffi.C.MINC2_SUCCESS)
end

-- create a  new minc2 file
function minc2_file:create(path)
    assert( path~=nil )
    assert( lib.minc2_create(self._v, path ) == ffi.C.MINC2_SUCCESS)
end

function minc2_file:copy_metadata(another)
    assert(another)
    assert(lib.minc2_copy_metadata(another._v,self._v)==ffi.C.MINC2_SUCCESS)
end


function minc2_file:load_complete_volume(data_type)
    -- will be torch tensors
    -- require('torch')
    local data_type=data_type or ffi.C.MINC2_FLOAT
    -- local buf_len=ffi.new("int[1]")
    -- lib.minc2_nelement(self._v,buf_len)
    -- buf_len=buf_len[0]
    local buf=nil
    local _dims=self:representation_dims()
    local dims=torch.LongStorage(self:ndim())
    -- local nelements=1
    
    -- Torch tensor defines dimensions in a slowest first fashion
    for i=0,(self:ndim()-1) do 
        dims[self:ndim()-i]=_dims[i].length
        --nelements=nelements*_dims[i].length
    end
    
    if data_type==ffi.C.MINC2_BYTE then 
        buf=torch.CharTensor(dims)
    elseif data_type==ffi.C.MINC2_UBYTE then 
        buf=torch.ByteTensor(dims)
    elseif data_type==ffi.C.MINC2_SHORT then 
        buf=torch.ShortTensor(dims)
    elseif data_type==ffi.C.MINC2_USHORT then 
        buf=torch.ShortTensor(dims)
    elseif data_type==ffi.C.MINC2_INT then 
        buf=torch.IntTensor(dims)
    elseif data_type==ffi.C.MINC2_UINT then 
        buf=torch.IntTensor(dims)
    elseif data_type==ffi.C.MINC2_FLOAT then 
        buf=torch.FloatTensor(dims)
    elseif data_type==ffi.C.MINC2_DOUBLE then 
        buf=torch.DoubleTensor(dims)
    else
        error("Unsupported  yet")
    end
    assert( 
        lib.minc2_load_complete_volume(self._v, buf:storage():data(), data_type)==ffi.C.MINC2_SUCCESS 
    )

    return buf
end


function minc2_file:load_hyperslab(data_type, slab)
    if not slab then
        return self:load_complete_volume(data_type)
    else
        -- will be torch tensors
        local data_type=data_type or ffi.C.MINC2_FLOAT
        -- local buf_len=ffi.new("int[1]")
        -- lib.minc2_nelement(self._v,buf_len)
        -- buf_len=buf_len[0]
        local buf=nil
        local _dims=self:representation_dims()
        local ndims=self:ndim()
        local dims=torch.LongStorage(ndims)

        local slab_start=ffi.new("int[?]",ndims)
        local slab_count=ffi.new("int[?]",ndims)
        -- local nelements=1

        -- Torch tensor defines dimensions in a slowest first fashion, so as slab
        for i=0,(ndims-1) do
            if slab[i+1] ~= nil then
                if type(slab[i+1])=='table' then
                    if #slab[i+1]==2 then
                        slab_count[ndims-1-i]=slab[i+1][2]-slab[i+1][1]+1
                        slab_start[ndims-1-i]=slab[i+1][1]-1
                    else -- assume it's the whole dimension
                        slab_count[ndims-1-i]=_dims[ndims-i-1].length
                        slab_start[ndims-1-i]=0
                    end
                else -- assume it's a number
                    slab_count[ndims-1-i]=1
                    slab_start[ndims-1-i]=slab[i+1]-1
                end
            else
                slab_count[ndims-1-i]=_dims[ndims-i-1].length
                slab_start[ndims-1-i]=0
            end
            dims[i+1]=slab_count[ndims-1-i]
        end

        if data_type==ffi.C.MINC2_BYTE then
            buf=torch.CharTensor(dims)
        elseif data_type==ffi.C.MINC2_UBYTE then
            buf=torch.ByteTensor(dims)
        elseif data_type==ffi.C.MINC2_SHORT then
            buf=torch.ShortTensor(dims)
        elseif data_type==ffi.C.MINC2_USHORT then
            buf=torch.ShortTensor(dims)
        elseif data_type==ffi.C.MINC2_INT then
            buf=torch.IntTensor(dims)
        elseif data_type==ffi.C.MINC2_UINT then
            buf=torch.IntTensor(dims)
        elseif data_type==ffi.C.MINC2_FLOAT then
            buf=torch.FloatTensor(dims)
        elseif data_type==ffi.C.MINC2_DOUBLE then
            buf=torch.DoubleTensor(dims)
        else
            error("Unsupported  yet")
        end
        local i

        assert(
            lib.minc2_read_hyperslab(self._v, slab_start, slab_count, buf:storage():data(), data_type)==ffi.C.MINC2_SUCCESS
        )

        return buf
    end
end


function minc2_file:setup_standard_order()
    assert( lib.minc2_setup_standard_order(self._v) == ffi.C.MINC2_SUCCESS)
end


function minc2_file:save_complete_volume(buf)
    assert(buf~=nil)
    --local t=require('torch')
    local data_type=ffi.C.MINC2_FLOAT
    -- local s=buf:storage()
    local store_type=torch.type(buf)
    

    -- TODO: implement dimension checking!
    -- TODO: check if tensor is contigious
    -- TODO: figure out how to save non-contigious tensor
    
    if store_type == 'torch.CharTensor' then
        data_type=ffi.C.MINC2_BYTE
    elseif store_type=='torch.ByteTensor' then 
        data_type=ffi.C.MINC2_UBYTE
    elseif store_type=='torch.ShortTensor' then 
        data_type=ffi.C.MINC2_SHORT
    elseif store_type=='torch.ShortTensor' then 
        data_type=ffi.C.MINC2_USHORT
    elseif store_type=='torch.IntTensor' then 
        data_type=ffi.C.MINC2_INT
    elseif store_type == 'torch.IntTensor' then 
        data_type=ffi.C.MINC2_UINT
    elseif store_type == 'torch.FloatTensor' then 
        data_type=ffi.C.MINC2_FLOAT
    elseif store_type == 'torch.DoubleTensor' then 
        data_type=ffi.C.MINC2_DOUBLE
    else
        print(string.format("store_type=%s",store_type))
        error("Unsupported  yet")
    end
    
    assert(
        lib.minc2_save_complete_volume(self._v,buf:storage():data(),data_type)==ffi.C.MINC2_SUCCESS
        )
    return buf
end


function minc2_file:set_volume_range(rmin,rmax)
    assert( rmin ~= nil)
    assert( rmax ~= nil)
    assert( lib.minc2_set_volume_range(self._v,rmin,rmax) == ffi.C.MINC2_SUCCESS)
end

function minc2_file:save_hyperslab(buf, start)
    if not start then
        return self:save_complete_volume(buf)
    else
        assert(buf~=nil)

        -- will be torch tensors
        local data_type=ffi.C.MINC2_FLOAT

        local ndims =self:ndim()
        local dims = buf:size()

        local slab_start=ffi.new("int[?]",ndims)
        local slab_count=ffi.new("int[?]",ndims)
        -- local nelements=1

        -- Torch tensor defines dimensions in a slowest first fashion, so as slab
        for i=0,(ndims-1) do
            if start[i+1] ~= nil then
                if type(start[i+1])=='table' then
                    if #start[i+1]==2 then
                        slab_count[ndims-1-i]=start[i+1][2]-start[i+1][1]+1
                        slab_start[ndims-1-i]=start[i+1][1]-1
                    else -- assume it's the whole dimension
                        slab_count[ndims-1-i]=dims[i+1]
                        slab_start[ndims-1-i]=0
                    end
                else -- assume it's a number
                    slab_count[ndims-1-i]=dims[i+1]
                    slab_start[ndims-1-i]=start[i+1]-1
                end
            else
                slab_count[ndims-1-i]=dims[i+1]
                slab_start[ndims-1-i]=0
            end
        end


        print(string.format("slab_start=[%d,%d,%d]",slab_start[0],slab_start[1],slab_start[2]))
        print(string.format("slab_count=[%d,%d,%d]",slab_count[0],slab_count[1],slab_count[2]))

        local data_type=ffi.C.MINC2_FLOAT
        local store_type=torch.type(buf)

        -- TODO: implement dimension checking!
        -- TODO: check if tensor is contigious
        -- TODO: figure out how to save non-contigious tensor

        if store_type == 'torch.CharTensor' then
            data_type=ffi.C.MINC2_BYTE
        elseif store_type=='torch.ByteTensor' then
            data_type=ffi.C.MINC2_UBYTE
        elseif store_type=='torch.ShortTensor' then
            data_type=ffi.C.MINC2_SHORT
        elseif store_type=='torch.ShortTensor' then
            data_type=ffi.C.MINC2_USHORT
        elseif store_type=='torch.IntTensor' then
            data_type=ffi.C.MINC2_INT
        elseif store_type == 'torch.IntTensor' then
            data_type=ffi.C.MINC2_UINT
        elseif store_type == 'torch.FloatTensor' then
            data_type=ffi.C.MINC2_FLOAT
        elseif store_type == 'torch.DoubleTensor' then
            data_type=ffi.C.MINC2_DOUBLE
        else
            print(string.format("store_type=%s",store_type))
            error("Unsupported  yet")
        end
        local i

        assert(
            lib.minc2_write_hyperslab(self._v, slab_start, slab_count, buf:storage():data(), data_type)==ffi.C.MINC2_SUCCESS
        )

        return buf
    end
end


function minc2_file:read_attribute(group, attribute)

    local attr_type=ffi.new("int[1]")
    local attr_length=ffi.new("int[1]")
    
    -- assume that if we can't get attribute type, it's missing, return nil

    if lib.minc2_get_attribute_type(self._v,group, attribute, attr_type)~=ffi.C.MINC2_SUCCESS then
        return nil
    end

    assert(lib.minc2_get_attribute_length(self._v,group,attribute,attr_length)==ffi.C.MINC2_SUCCESS)

    if attr_type[0] == ffi.C.MINC2_STRING then
        local buf = ffi.new("uint8_t[?]", attr_length[0])
        assert(lib.minc2_read_attribute(self._v,group,attribute,buf,attr_length[0])==ffi.C.MINC2_SUCCESS);
        return ffi.string(buf, attr_length[0])
    else
        local buf
        local dims=attr_length[0]
        local data_type=attr_type[0]
        
        if data_type==ffi.C.MINC2_BYTE then 
            buf=torch.CharTensor(dims)
        elseif data_type==ffi.C.MINC2_UBYTE then 
            buf=torch.ByteTensor(dims)
        elseif data_type==ffi.C.MINC2_SHORT then 
            buf=torch.ShortTensor(dims)
        elseif data_type==ffi.C.MINC2_USHORT then 
            buf=torch.ShortTensor(dims)
        elseif data_type==ffi.C.MINC2_INT then 
            buf=torch.IntTensor(dims)
        elseif data_type==ffi.C.MINC2_UINT then 
            buf=torch.IntTensor(dims)
        elseif data_type==ffi.C.MINC2_FLOAT then 
            buf=torch.FloatTensor(dims)
        elseif data_type==ffi.C.MINC2_DOUBLE then 
            buf=torch.DoubleTensor(dims)
        else
            error("Unsupported  yet:"..data_type)
        end

        assert(lib.minc2_read_attribute(self._v,group,attribute,buf:storage():data(),attr_length[0])==ffi.C.MINC2_SUCCESS);
        
        return buf
    end
end

function minc2_file:write_attribute(group,attribute,value)
    -- local attr_type=ffi.new("int[1]")
    -- local attr_length=ffi.new("int[1]")
    local dtype=type(value)
    
    if dtype=="string" then
        attr_type=ffi.C.MINC2_STRING
        attr_length=#value+1
        assert(
            lib.minc2_write_attribute(self._v,group,attribute,ffi.cast("const char[]",value),#value+1,ffi.C.MINC2_STRING)==ffi.C.MINC2_SUCCESS
            )
    else
        local _value=value
        if dtype=="table" then 
            _value=torch.Tensor(value)
        elseif dtype=="number" then
            _value=torch.Tensor(1)
            _value[1]=value
        end
        
        local store_type=torch.type(_value)
        local data_type
        
        if store_type == 'torch.CharTensor' then
            data_type=ffi.C.MINC2_BYTE
        elseif store_type=='torch.ByteTensor' then 
            data_type=ffi.C.MINC2_UBYTE
        elseif store_type=='torch.ShortTensor' then 
            data_type=ffi.C.MINC2_SHORT
        elseif store_type=='torch.ShortTensor' then 
            data_type=ffi.C.MINC2_USHORT
        elseif store_type=='torch.IntTensor' then 
            data_type=ffi.C.MINC2_INT
        elseif store_type == 'torch.IntTensor' then 
            data_type=ffi.C.MINC2_UINT
        elseif store_type == 'torch.FloatTensor' then 
            data_type=ffi.C.MINC2_FLOAT
        elseif store_type == 'torch.DoubleTensor' then 
            data_type=ffi.C.MINC2_DOUBLE
        else
            error("Unsupported  yet:"..store_type)
        end
        
        assert(
            lib.minc2_write_attribute(self._v,group,attribute,_value:storage():data(),_value:storage():size(),data_type)==ffi.C.MINC2_SUCCESS
            )
    end
end

function minc2_file:metadata()
    local ret={}
    
    local group_iterator=ffi.gc(lib.minc2_allocate_info_iterator(), lib.minc2_free_info_iterator)
    local attr_iterator=ffi.gc(lib.minc2_allocate_info_iterator(),lib.minc2_free_info_iterator)

    assert(lib.minc2_start_group_iterator(self._v,group_iterator)==ffi.C.MINC2_SUCCESS)
    
    while lib.minc2_iterator_group_next(group_iterator)==ffi.C.MINC2_SUCCESS do
        local gname=lib.minc2_iterator_group_name(group_iterator)
        assert(lib.minc2_start_attribute_iterator(self._v, gname, attr_iterator)==ffi.C.MINC2_SUCCESS)
        local g={}
        
        while lib.minc2_iterator_attribute_next(attr_iterator)==ffi.C.MINC2_SUCCESS do
            local aname=lib.minc2_iterator_attribute_name(attr_iterator)
            g[ ffi.string(aname) ] = self:read_attribute(gname, aname)
        end
        
        ret[ ffi.string(lib.minc2_iterator_group_name(group_iterator)) ] = g
        lib.minc2_stop_info_iterator(attr_iterator)
    end
    lib.minc2_stop_info_iterator(group_iterator)
    return ret
end

function minc2_file:write_metadata(m)
    local group,g
    for group,g in pairs(m) do
        local attr,a
        for attr,a in pairs(g) do
            self:write_attribute(group,attr,a)
        end
    end
end

function minc2_file:world_to_voxel(xyz)
    -- convert world (X,Y,Z) coordinates to index (0-based) (i,j,k) coordinates
    local dtype=type(xyz)
    if dtype=="table" then -- return table too
        local _xyz=ffi.new("double[3]",xyz)
        local _ijk=ffi.new("double[3]")
        assert(lib.minc2_world_to_voxel(self._v,_xyz,_ijk)==ffi.C.MINC2_SUCCESS)
        return {_ijk[0],_ijk[1],_ijk[2]}
    else -- assume it is a torch tensor, return tensor
        local ijk=torch.Tensor(3)
        assert(lib.minc2_world_to_voxel(self._v,xyz:storage():data(),ijk:storage():data())==ffi.C.MINC2_SUCCESS)
        return ijk
    end
end

function minc2_file:voxel_to_world(ijk)
    -- convert index (0-based) (i,j,k) to world (X,Y,Z) coordinates
    local dtype=type(ijk)
    if dtype=="table" then -- return table too
        local _xyz=ffi.new("double[3]")
        local _ijk=ffi.new("double[3]",ijk)
        assert(lib.minc2_voxel_to_world(self._v,_ijk,_xyz)==ffi.C.MINC2_SUCCESS)
        return {_xyz[0],_xyz[1],_xyz[2]}
    else -- assume it is a torch tensor, return tensor
        local xyz=torch.Tensor(3)
        assert(lib.minc2_voxel_to_world(self._v,ijk:storage():data(),xyz:storage():data())==ffi.C.MINC2_SUCCESS)
        return xyz
    end
end


minc2_xfm = {
    -- XFM types
    MINC2_XFM_LINEAR                 = ffi.C.MINC2_XFM_LINEAR,
    MINC2_XFM_THIN_PLATE_SPLINE      = ffi.C.MINC2_XFM_THIN_PLATE_SPLINE,
    MINC2_XFM_USER_TRANSFORM         = ffi.C.MINC2_XFM_USER_TRANSFORM,
    MINC2_XFM_CONCATENATED_TRANSFORM = ffi.C.MINC2_XFM_CONCATENATED_TRANSFORM,
    MINC2_XFM_GRID_TRANSFORM         = ffi.C.MINC2_XFM_GRID_TRANSFORM
}
minc2_xfm.__index = minc2_xfm


function minc2_xfm.new(path)
  local self = setmetatable({}, minc2_xfm)
  self._v=ffi.gc(lib.minc2_xfm_allocate0(),lib.minc2_xfm_destroy)
  if path~=nil then
      self:open(path)
  end
  return self
end


function minc2_xfm:open(path)
    --print("Going to open:"..path)
    assert(path~=nil,"Provide minc2 file")
    assert( lib.minc2_xfm_open(self._v,path) == ffi.C.MINC2_SUCCESS )
end


function minc2_xfm:save(path)
    --print("Going to open:"..path)
    assert(path~=nil,"Provide minc2 file")
    assert( lib.minc2_xfm_save(self._v,path) == ffi.C.MINC2_SUCCESS )
end


function minc2_xfm:transform_point(xyz_in)
    local dtype=type(xyz_in)
    if dtype=="table" then -- return table too
        local _xyz_in=ffi.new("double[3]",xyz_in)
        local _xyz_out=ffi.new("double[3]")
        assert(lib.minc2_xfm_transform_point(self._v,_xyz_in,_xyz_out)==ffi.C.MINC2_SUCCESS)
        return {_xyz_out[0],_xyz_out[1],_xyz_out[2]}
    else -- assume it is a torch tensor, return tensor
        local xyz_out=torch.Tensor(3)
        assert(lib.minc2_xfm_transform_point(self._v,xyz_in:storage():data(),xyz_out:storage():data())==ffi.C.MINC2_SUCCESS)
        return xyz_out
    end
end


function minc2_xfm:inverse_transform_point(xyz_in)
    local dtype=type(xyz_in)
    if dtype=="table" then -- return table too
        local _xyz_in=ffi.new("double[3]",xyz_in)
        local _xyz_out=ffi.new("double[3]")
        assert(lib.minc2_xfm_inverse_transform_point(self._v,_xyz_in,_xyz_out)==ffi.C.MINC2_SUCCESS)
        return {_xyz_out[0],_xyz_out[1],_xyz_out[2]}
    else -- assume it is a torch tensor, return tensor
        local xyz_out=torch.Tensor(3)
        assert(lib.minc2_xfm_inverse_transform_point(self._v,xyz_in:storage():data(),xyz_out:storage():data())==ffi.C.MINC2_SUCCESS)
        return xyz_out
    end
end

function minc2_xfm:invert()
    assert(lib.minc2_xfm_invert(self._v)==ffi.C.MINC2_SUCCESS)
end

function minc2_xfm:get_n_concat()
    local n=ffi.new("int[1]")
    assert(lib.minc2_xfm_get_n_concat(self._v,n)==ffi.C.MINC2_SUCCESS)
    return n[0]
end

function minc2_xfm:get_n_type(n)
    local n=n or 0
    local t=ffi.new("int[1]")
    assert(lib.minc2_xfm_get_n_type(self._v,n,t)==ffi.C.MINC2_SUCCESS)
    return t[0]
end

function minc2_xfm:get_grid_transform(n)
    local n=n or 0
    local c_file=ffi.new("char*[1]")
    local inv=ffi.new("int[1]")
    assert(lib.minc2_xfm_get_grid_transform(self._v,n,inv,c_file)==ffi.C.MINC2_SUCCESS)
    local _file=ffi.string(c_file[0])
    ffi.C.free(c_file[0])
    return _file,inv[0]
end

function minc2_xfm:get_linear_transform(n)
    local n=n or 0
    local mat=torch.eye(4)
    assert(lib.minc2_xfm_get_linear_transform(self._v,n,mat:storage():data())==ffi.C.MINC2_SUCCESS)
    return mat
end

function minc2_xfm:get_linear_transform_param(n,center)
    local n=n or 0
    
    local _center=torch.Tensor(3):zero()
    if center ~= nil then
      _center=torch.Tensor(center)
    end
    
    local transform={
      center=_center,
      translations=torch.Tensor(3):zero(),
      scales=torch.Tensor(3):zero(),
      shears=torch.Tensor(3):zero(),
      rotations=torch.Tensor(3):zero(),
      invalid=false
    }
    
    if lib.minc2_xfm_extract_linear_param(self._v,n,
        transform.center:storage():data(),
        transform.translations:storage():data(),
        transform.scales:storage():data(),
        transform.shears:storage():data(),
        transform.rotations:storage():data()
       )~=ffi.C.MINC2_SUCCESS then
           transform.invalid=true
    else
        transform.invalid=false
    end
    return transform
end

function minc2_xfm:append_linear_transform(par)
    if torch.type(par)=="torch.DoubleTensor" then
        assert(lib.minc2_xfm_append_linear_transform(self._v,mat:storage():data())==ffi.C.MINC2_SUCCESS)
    else --assume it is parameters table
        assert(par and
            par.center and
            par.translations and
            par.scales and
            par.shears and
            par.rotations)

        assert(lib.minc2_xfm_append_linear_param(self._v,
            par.center:storage():data(),
            par.translations:storage():data(),
            par.scales:storage():data(),
            par.shears:storage():data(),
            par.rotations:storage():data()
            )==ffi.C.MINC2_SUCCESS)
    end
    return self
end

function minc2_xfm:append_grid_transform(grid_file,inv)
    local inv=inv or 0
    assert(lib.minc2_xfm_append_grid_transform(self._v,grid_file,inv)==ffi.C.MINC2_SUCCESS)
    return self
end

function minc2_xfm:concat_xfm(another)
    assert(lib.minc2_xfm_concat_xfm(self._v,another._v)==ffi.C.MINC2_SUCCESS)
end

-- this is all we have in the module
return { 
    -- minc2 file reader/writer
    minc2_file=minc2_file 
}

-- kate: indent-width 4; replace-tabs on; remove-trailing-space on; hl lua
