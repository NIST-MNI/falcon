from __future__ import print_function

from ._simple import ffi,lib
from .utils   import to_bytes,to_unicode
from .utils   import text_type
import six
import sys
import collections

class minc2_error(Exception):
    """
    minc2-generated error
    """
    def __init__(self, message=""):
        super(minc2_error, self).__init__(message)
    

class minc2_transform_parameters(object):
    """
    parameters describing affine transformation
    """
    def __init__(self):
        import numpy as np
        self.center=np.zeros(3)
        self.translations=np.zeros(3)
        self.scales=np.zeros(3)
        self.shears=np.zeros(3)
        self.rotations=np.zeros(3)
        self.invalid=False

    def __str__(self):
        return "center: {} {} {}\n".format(self.center[0],self.center[1],self.center[2])+ \
               "translations: {} {} {}\n".format(self.translations[0],self.translations[1],self.translations[2])+ \
               "scales: {} {} {}\n".format(self.scales[0],self.scales[1],self.scales[2])+ \
               "rotations: {} {} {}\n".format(self.rotations[0],self.rotations[1],self.rotations[2])+ \
               "shears: {} {} {}\n".format(self.shears[0],self.shears[1],self.shears[2])

    def __repr__(self):
        return self.__str__()


minc2_dim=collections.namedtuple('minc2_dim',['id','length', 'start', 'step', 'have_dir_cos', 'dir_cos'])


class minc2_file:
    """
    MINC2 file object (Volume on the disk, stored in .mnc file)

    """
    #: constants
    MINC2_DIM_UNKNOWN=lib.MINC2_DIM_UNKNOWN
    MINC2_DIM_X    = lib.MINC2_DIM_X
    MINC2_DIM_Y    = lib.MINC2_DIM_Y
    MINC2_DIM_Z    = lib.MINC2_DIM_Z
    MINC2_DIM_TIME = lib.MINC2_DIM_TIME
    MINC2_DIM_VEC  = lib.MINC2_DIM_VEC
    MINC2_DIM_END  = lib.MINC2_DIM_END
    
    #: minc2 data types
    MINC2_BYTE     = lib.MINC2_BYTE 
    MINC2_SHORT    = lib.MINC2_SHORT
    MINC2_INT      = lib.MINC2_INT 
    MINC2_FLOAT    = lib.MINC2_FLOAT 
    MINC2_DOUBLE   = lib.MINC2_DOUBLE 
    MINC2_STRING   = lib.MINC2_STRING 
    MINC2_UBYTE    = lib.MINC2_UBYTE 
    MINC2_USHORT   = lib.MINC2_USHORT 
    MINC2_UINT     = lib.MINC2_UINT 
    MINC2_SCOMPLEX = lib.MINC2_SCOMPLEX 
    MINC2_ICOMPLEX = lib.MINC2_ICOMPLEX 
    MINC2_FCOMPLEX = lib.MINC2_FCOMPLEX 
    MINC2_DCOMPLEX = lib.MINC2_DCOMPLEX 
    MINC2_MAX_TYPE_ID = lib.MINC2_MAX_TYPE_ID
    MINC2_UNKNOWN  = lib.MINC2_UNKNOWN  

    # minc2 status
    MINC2_SUCCESS  = lib.MINC2_SUCCESS,
    MINC2_ERROR    = lib.MINC2_ERROR

    # data types
    __minc2_to_numpy = {
            lib.MINC2_BYTE:   'int8',
            lib.MINC2_UBYTE:  'uint8',
            lib.MINC2_SHORT:  'int16',
            lib.MINC2_USHORT: 'uint16',
            lib.MINC2_INT:    'int32',
            lib.MINC2_UINT:   'uint32',
            lib.MINC2_FLOAT:  'float32',
            lib.MINC2_DOUBLE: 'float64',
        }
    __numpy_to_minc2 = {y: x for x, y in six.iteritems(__minc2_to_numpy)}

    __minc2_to_torch = {
            lib.MINC2_BYTE:   'torch.CharTensor',
            lib.MINC2_UBYTE:  'torch.ByteTensor',
            lib.MINC2_SHORT:  'torch.ShortTensor',
            lib.MINC2_USHORT: 'torch.ShortTensor', # WARNING: no support for unsigned short
            lib.MINC2_INT:    'torch.IntTensor',
            lib.MINC2_UINT:   'torch.IntTensor', # WARNING: no support for unsigned int
            lib.MINC2_FLOAT:  'torch.FloatTensor',
            lib.MINC2_DOUBLE: 'torch.DoubleTensor'
        }
    __torch_to_minc2 = {y:x for x,y in six.iteritems(__minc2_to_torch)}
    
    minc2_to_numpy=__minc2_to_numpy
    numpy_to_minc2=__numpy_to_minc2

    def __init__(self, path=None, standard=False, handle=None):
        """

        :param path: file path to open
        :param standard: after opening, immedeately switch to standard orientation ( see @setup_standard_order)
        :param handle: use library handle to already opened minc2 file
        """

        if handle is None:
            self._v = ffi.gc(lib.minc2_allocate0(), lib.minc2_destroy)
        else:
            self._v = handle

        if path is not None:
            self.open(path)
            if standard:
                self.setup_standard_order()

    def open(self, path):
        """
        Open existing minc2 file
        :param path: file path
        :return: None
        """
        if lib.minc2_open(self._v, to_bytes(path))!=lib.MINC2_SUCCESS:
            raise minc2_error("Can't open file:"+path)

    def close(self):
        """
        Close opened file, flush data on disk, free resources
        :return:
        """
        if lib.minc2_close(self._v)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error closing file")

    def ndim(self):
        """
        Query number of dimensions
        :return: integer number of dimensions
        """
        dd=ffi.new("int*", 0)
        lib.minc2_ndim(self._v, dd)
        return dd[0]

    def store_dims_(self):
        """
        internal function to query dimension structure
        :return: struct minc2_dimension
        """
        dims=ffi.new("struct minc2_dimension*[1]")
        if lib.minc2_get_store_dimensions(self._v,dims)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error defining dimensions")
        return dims[0]

    def minc2_dim_to_python_(self,d):
        """
        internal function to convert dimension descriptions
        :param d:  dimensions description in minc2 C library
        :return:  dimensions description in python format
        """
        import numpy as np
        dd = minc2_dim(id=d.id,length=d.length, start=d.start, step=d.step, have_dir_cos=d.have_dir_cos, dir_cos=np.zeros(3,np.float64))
        if d.have_dir_cos:
            ffi.memmove(ffi.cast("double [3]", dd.dir_cos.ctypes.data), d.dir_cos, 3*ffi.sizeof('double'))
        return dd

    def store_dims(self):
        """
        Description of file dimensions (as stored on disk)
        :return:  list of dimensions , in the order as stored on disk
        """
        d_=self.store_dims_()
        return [ self.minc2_dim_to_python_(d_[j]) for j in range(self.ndim())]

    def representation_dims_(self):
        """
        internal function
        :return:
        """
        dims=ffi.new("struct minc2_dimension*[1]")
        if lib.minc2_get_representation_dimensions(self._v,dims)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error getting dimension information")
        return dims[0]

    def representation_dims(self):
        """
        Description of file dimensions (current memory representation)
        :return:  list of dimensions , as visible to python
        """
        d_ = self.representation_dims_()
        return [ self.minc2_dim_to_python_(d_[j]) for j in range(self.ndim() ) ]

    def imitate(self, another, store_type=None, representation_type=None, path=None):
        """
        Generate new minc volume with the same parameters as another one

        :param another: another minc2_file object or path to a file
        :param store_type: data type for storage
        :param representation_type: preferred data type for representation
        :param path file path to generate
        :return: None
        """

        if not isinstance(another, minc2_file):
            another_ = minc2_file(another)
        else:
            another_ = another

        if store_type is None:
            store_type = another_.store_dtype()
        if representation_type is None:
            representation_type = another_.representation_dtype()

        # copy data
        self.define(another_.store_dims(), store_type=store_type, representation_type=representation_type)
        if path is not None:
            self.create(path)

    def define(self, dims, store_type=None, representation_type=None, slice_scaling=None, global_scaling=None, path=None):
        """
        Define new minc2 volume
        :param dims:  dimensions description (as will be stored on disk)
        :param store_type:  data format as stored on disk
        :param representation_type:  data format as visible to python
        :param slice_scaling:  use slice scaling
        :param global_scaling: use global scaling
        :param path: output file path
        :return:
        """
        _dims = dims
        if store_type is None:
            store_type=lib.MINC2_SHORT
        if representation_type is None:
            representation_type=store_type

        _store_type = store_type
        _representation_type = representation_type

        if isinstance(_store_type, six.string_types):
            _store_type = minc2_file.__numpy_to_minc2[_store_type]

        if isinstance(_representation_type, six.string_types):
            _representation_type = minc2_file.__numpy_to_minc2[_representation_type]

        if isinstance(dims, list ) or isinstance(dims,tuple):
            _dims = ffi.new("struct minc2_dimension[]", len(dims)+1)
            for i,j in enumerate(dims):
                if isinstance(j, minc2_dim):
                    _dims[i].id=j.id
                    _dims[i].length=j.length
                    _dims[i].start=j.start
                    _dims[i].step=j.step
                    _dims[i].have_dir_cos=j.have_dir_cos
                    if j.have_dir_cos: 
                        ffi.memmove(_dims[i].dir_cos, ffi.cast("double [3]", j.dir_cos.ctypes.data ), 3*ffi.sizeof('double'))
                else:
                    _dims[i]=j
            _dims[len(dims)]={'id':lib.MINC2_DIM_END}

        if lib.minc2_define(self._v, _dims, _store_type, _representation_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error defining new minc file")

        if slice_scaling is not None or global_scaling is not None:
            _slice_scaling=0
            _global_scaling=0

            if slice_scaling: _slice_scaling=1
            if global_scaling: _global_scaling=1

            if lib.minc2_set_scaling(self._v,_global_scaling,_slice_scaling )!=lib.MINC2_SUCCESS:
                raise minc2_error()
        if path is not None:
            self.create(path)

    def create(self, path):
        """
        Create new minc2 file, using predefined parameters
        :param path: file path
        :return:
        """
        if lib.minc2_create(self._v, to_bytes(path) )!=lib.MINC2_SUCCESS:
            raise minc2_error("Error creating file:"+path)
    
    def copy_metadata(self, another):
        """
        Copy minc2 metadata (headers)
        :param another: another volume
        :return:
        """
        if lib.minc2_copy_metadata(another._v,self._v)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error copying metadata")
    
    def load_complete_volume(self, data_type=None):
        """
        Load the whole volume as numpy ndarray
        :param data_type: python data required
        :return: numpy.ndarray
        """
        import numpy as np
        if data_type is None:
            data_type = np.dtype( self.representation_dtype() )
        buf=None
        _dims=self.representation_dims()
        # dims=torch.LongStorage(self:ndim())
        shape=list(range(self.ndim()))
        # numpy array  defines dimensions in a slowest first fashion
        for i in range(self.ndim()):
            shape[self.ndim()-i-1] = _dims[i].length

        dtype=None
        
        if data_type in minc2_file.__minc2_to_numpy:
            dtype = minc2_file.__minc2_to_numpy[data_type]
        elif data_type in minc2_file.__numpy_to_minc2:
            dtype = data_type
            data_type = minc2_file.__numpy_to_minc2[dtype]
        elif isinstance(data_type, np.dtype):
            dtype = data_type
            data_type = minc2_file.__numpy_to_minc2[dtype.name]
        else:
            raise minc2_error("Unsupported buffer data type:"+repr(data_type))

        buf = np.empty(shape, dtype, 'C')
        if lib.minc2_load_complete_volume(self._v, ffi.cast("void *", buf.ctypes.data) , data_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error loading volume")
        return buf

    def load_complete_volume_tensor(self, data_type=None):
        """
        Load the whole volume as pytorch tensor
        :param data_type: python data required
        :return: torch.Tensor
        """
        import torch
        if data_type is None:
            data_type=self.representation_dtype_tensor()
        buf=None
        _dims=self.representation_dims()
        # dims=torch.LongStorage(self:ndim())
        shape=list(range(self.ndim()))
        # numpy array  defines dimensions in a slowest first fashion
        for i in range(self.ndim()):
            shape[self.ndim()-i-1]=_dims[i].length

        dtype=None

        if data_type in minc2_file.__minc2_to_torch:
            dtype=eval(minc2_file.__minc2_to_torch[data_type])
        elif data_type in minc2_file.__torch_to_minc2:
            dtype=eval(data_type)
            data_type=minc2_file.__torch_to_minc2[data_type]
        else:
            raise minc2_error("Unsupported data type:"+repr(data_type))

        buf=dtype(*shape)

        if lib.minc2_load_complete_volume(self._v, ffi.cast("void *", buf.storage().data_ptr()) , data_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error loading volume")
        return buf

    def setup_standard_order(self):
        """
        Request library to use stamdard order: positive step sizes
        dimension order (slowest to fastest): time, Z, Y, X, Vector
        :return:
        """
        if lib.minc2_setup_standard_order(self._v)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error setting up standard volume order")

    def save_complete_volume(self, buf):
        """
        Dump whole numpy.ndarray into minc2 volume
        volume have to be open for writing and initialized (i.e .create or .imitate have been called)
        :param buf: numpy.ndarray
        :return: numpy.ndarray
        """
        import numpy as np
        data_type=lib.MINC2_FLOAT
        store_type=buf.dtype.name
        # TODO: make sure array is in "C" order
        
        assert(store_type in minc2_file.__numpy_to_minc2)
        # TODO: verify dimensions of the array

        data_type=minc2_file.__numpy_to_minc2[store_type]
        
        if lib.minc2_save_complete_volume(self._v,ffi.cast("void *", buf.ctypes.data),data_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error saving volume")
        return buf

    def save_complete_volume_tensor(self, buf):
        """
        Dump whole torch.Tensor into minc2 volume
        volume have to be open for writing and initialized (i.e .create or .imitate have been called)
        :param buf: torch.Tensor
        :return: torch.Tensor
        """
        #import torch
        data_type=lib.MINC2_FLOAT
        store_type=buf.type()

        assert(store_type in minc2_file.__torch_to_minc2)
        # TODO: verify dimensions of the array

        data_type=minc2_file.__torch_to_minc2[store_type]

        if lib.minc2_save_complete_volume(self._v,ffi.cast("void *", buf.storage().data_ptr()),data_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error saving volume")
        return buf

    def world_to_voxel(self, xyz):
        """
        Convert world coordinates (X,Y,Z) to voxel coordinates (i,j,k)

        :param xyz: - numpy array either length of 3 or 2D array with each row being X,Y,Z
        :return:    - either 1D array of (i,j,k) or 2D array with each row (i,j,k)
        """
        import numpy as np
        in_xyz=np.asarray(xyz,dtype=np.float64,order="C")
        if len(in_xyz.shape)==1: # single 3D array
            out_ijk = np.empty(3, np.float64, 'C')
            if lib.minc2_world_to_voxel(self._v, ffi.cast("double *", in_xyz.ctypes.data), ffi.cast("double *", out_ijk.ctypes.data))!=lib.MINC2_SUCCESS:
                raise minc2_error("Error world_to_voxel")
            return out_ijk
        else:
            out_ijk = np.empty(in_xyz.shape, np.float64, 'C')
            if lib.minc2_world_to_voxel_vec(self._v, in_xyz.shape[0], 3, ffi.cast("double *", in_xyz.ctypes.data), ffi.cast("double *", out_ijk.ctypes.data))!=lib.MINC2_SUCCESS:
                raise minc2_error("Error world_to_voxel_vec")
            return out_ijk

    def voxel_to_world(self, ijk):
        """
        Convert voxel coordinates (i,j,k) to world coordinates (X,Y,Z)

        :param ijk: - numpy array either length of 3 or 2D array with each row being i,j,k
        :return:    - either 1D array of (X,Y,Z) or 2D array with each row (X,Y,Z)
        """
        import numpy as np
        in_ijk=np.asarray(ijk, dtype=np.float64, order="C")
        if len(in_ijk.shape)==1: # single 3D array
            out_xyz = np.empty(3, np.float64, 'C')
            if lib.minc2_voxel_to_world(self._v, ffi.cast("double *", in_ijk.ctypes.data), ffi.cast("double *", out_xyz.ctypes.data)) != lib.MINC2_SUCCESS:
                raise minc2_error("Error in voxel_to_world")
            return out_xyz
        else:
            out_xyz = np.empty(in_ijk.shape, np.float64, 'C')
            if lib.minc2_world_to_voxel_vec(self._v, in_ijk.shape[0], 3, ffi.cast("double *", in_ijk.ctypes.data), ffi.cast("double *", out_xyz.ctypes.data)) != lib.MINC2_SUCCESS:
                raise minc2_error("Error in voxel_to_world_vec")
            return out_xyz


    def read_attribute(self, group, attribute):
        """
        Read minc2 header attribute
        :param group: attribute group name
        :param attribute: attribute name
        :return:
        """
        import numpy as np

        attr_type=ffi.new("int*",0)
        attr_length=ffi.new("int*",0)

        if isinstance(group, six.string_types ):
            group = to_bytes(group)
        if isinstance(attribute, six.string_types ):
            attribute = to_bytes(attribute)

        # assume that if we can't get attribute type, it's missing, return nil

        if lib.minc2_get_attribute_type(self._v, group, attribute, attr_type)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error getting attribute type {}:{}".format(group,attribute))

        if lib.minc2_get_attribute_length(self._v, group, attribute, attr_length)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error getting attribute length {}:{}".format(group,attribute))
            
        if attr_type[0] == lib.MINC2_STRING:
            buf = ffi.new("char[]", attr_length[0])
            if lib.minc2_read_attribute(self._v,group,attribute,buf,attr_length[0])!=lib.MINC2_SUCCESS:
                raise minc2_error("Error reading string attribute {}:{}".format(group,attribute))
            return to_unicode( ffi.string(buf, attr_length[0]) )
        else:
            data_type = attr_type[0]
            buf=None
            if data_type in minc2_file.__minc2_to_numpy:
                dtype=minc2_file.__minc2_to_numpy[data_type]
                shape=[attr_length[0]]
                buf=np.empty(shape,dtype,'C')
            else:
                raise minc2_error("Error determining attribute type {}:{}".format(group,attribute))

            if lib.minc2_read_attribute(self._v,group,attribute,ffi.cast("void *", buf.ctypes.data),attr_length[0])!=lib.MINC2_SUCCESS:
                raise minc2_error("Error reading attribute {}:{}".format(group,attribute))

            return buf

    def write_attribute(self, group, attribute, value):
        """
        Store attribute into minc2 file
        :param group:  group name
        :param attribute:  attribute name
        :param value:  attribute value
        :return:
        """
        if isinstance(group, six.string_types ):
            group=to_bytes(group)
        if isinstance(attribute, six.string_types ):
            attribute=to_bytes(attribute)

        if isinstance(value, six.string_types ):
            value=to_bytes(value)
            attr_type=lib.MINC2_STRING
            attr_length=len(value)

            if lib.minc2_write_attribute(self._v, group, attribute, value, attr_length+1,lib.MINC2_STRING)!=lib.MINC2_SUCCESS:
                raise minc2_error("Error writing attribute {}:{}".format(group,attribute))

        elif isinstance(value, six.binary_type): # assume it's already binary encoded
            attr_type = lib.MINC2_STRING
            attr_length = len(value)

            if lib.minc2_write_attribute(self._v, group, attribute, value, attr_length + 1,
                                         lib.MINC2_STRING) != lib.MINC2_SUCCESS:
                raise minc2_error("Error writing attribute {}:{}".format(group, attribute))
        else:
            import numpy as np
            data_type=lib.MINC2_FLOAT
            store_type=value.dtype.name

            assert(store_type in minc2_file.__numpy_to_minc2)
            data_type=minc2_file.__numpy_to_minc2[store_type]

            if lib.minc2_write_attribute(self._v, group, attribute, ffi.cast("void *", value.ctypes.data), value.size, data_type)!=lib.MINC2_SUCCESS:
                raise minc2_error("Error writing attribute {}:{}".format(group,attribute))

    def metadata(self):
        """
        Read complete metadata from minc2 volume into a dictionary
        :return: dict
        """
        ret={}

        group_iterator=ffi.gc(lib.minc2_allocate_info_iterator(), lib.minc2_free_info_iterator)
        attr_iterator=ffi.gc(lib.minc2_allocate_info_iterator(), lib.minc2_free_info_iterator)

        if lib.minc2_start_group_iterator(self._v,group_iterator)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error navigating groups")

        while lib.minc2_iterator_group_next(group_iterator)==lib.MINC2_SUCCESS:
            gname=lib.minc2_iterator_group_name(group_iterator)
            if lib.minc2_start_attribute_iterator(self._v, gname, attr_iterator)!=lib.MINC2_SUCCESS:
                raise minc2_error("Error iterating group:"+gname)
            g={}

            while lib.minc2_iterator_attribute_next(attr_iterator)==lib.MINC2_SUCCESS:
                aname=lib.minc2_iterator_attribute_name(attr_iterator)
                g[ to_unicode(ffi.string(aname)) ] = self.read_attribute(ffi.string(gname), ffi.string(aname))

            ret[ to_unicode( ffi.string(lib.minc2_iterator_group_name(group_iterator)) ) ] = g
            lib.minc2_stop_info_iterator(attr_iterator)

        lib.minc2_stop_info_iterator(group_iterator)

        # special attributes
        ret['']={}
        ret['']['history']=self.read_attribute('', 'history')
        ret['']['ident']=self.read_attribute('', 'ident')
        ret['']['minc_version']=self.read_attribute('', 'minc_version')

        return ret

    def write_metadata(self, m):
        """
        Dump complete dictionary as metadata info for minc2 volume (see @metadata)
        :param m:
        :return:
        """
        for group,g in six.iteritems(m):
            for attr,a in six.iteritems(g):
                self.write_attribute(group,attr,a)

    def store_dtype(self):
        """
        query storage datatype (disk representation)
        :return: numpy.dtype describing disk storage
        """
        _dtype = ffi.new("int*",0)
        if lib.minc2_storage_data_type(self._v,_dtype)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error setting store type")
        return minc2_file.__minc2_to_numpy[_dtype[0]]

    def representation_dtype(self):
        """
        query representation datatype (as python sees the volume)
        :return: numpy.dtype describing memory representation
        """
        _dtype = ffi.new("int*",0)
        if lib.minc2_data_type(self._v,_dtype)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error setting representation type")
        return minc2_file.__minc2_to_numpy[_dtype[0]]

    def representation_dtype_tensor(self):
        """
        query representation datatype (as python sees the volume)
        :return: torch.type describing memory representation
        """

        _dtype=ffi.new("int*",0)
        if lib.minc2_data_type(self._v,_dtype)!=lib.MINC2_SUCCESS:
            raise minc2_error("Error setting representation type")
        return minc2_file.__minc2_to_torch[_dtype[0]]


    def get_data(self):
        """
        shortcut to setup_standard_order;load_complete_volume
        :return: numpy.ndarray
        """
        self.setup_standard_order()
        return self.load_complete_volume()


    def get_tensor(self):
        """
        shortcut to setup_standard_order;load_complete_volume_tensor
        :return: torch.Tensor
        """
        self.setup_standard_order()
        return self.load_complete_volume_tensor()


    def set_data(self,new_data):
        """
        shortcut to setup_standard_order,save_complete_volume
        :param new_data: numpy.ndarray
        :return:
        """
        self.setup_standard_order()
        self.save_complete_volume(new_data)


    def set_tensor(self,new_data):
        """
        shortcut to setup_standard_order,save_complete_volume_tensor
        :param new_data: torch.Tensor
        :return:
        """
        self.setup_standard_order()
        self.save_complete_volume_tensor(new_data)

    def get_shape(self):
        """
        get volume shape (numpy style) as memory representation
        :return tuple of dimension sizes
        """
        _dims=self.representation_dims_()
        return tuple( _dims[i].length for i in range(self.ndim()) )

    def get_start(self):
        """
        get volume starts (numpy style)
        :return tuple of starts
        """
        _dims=self.representation_dims_()
        return tuple(_dims[i].start for i in range(self.ndim()))

    def get_step(self):
        """
        get volume steps (numpy style)
        :return tuple of step sizes
        """
        _dims=self.representation_dims_()
        return tuple(_dims[i].step for i in range(self.ndim()))

    # get/set whole volume
    data = property(get_data, set_data,None,"Complete volume")

    # get/set whole volume
    tensor = property(get_tensor, set_tensor,None,"Complete volume")

    # get volume shape
    shape = property(get_shape,None,None,"Volume shape (numpy style)")

    # get volume start
    start = property(get_start,None,None,"Volume start (numpy style)")

    # get volume step
    step = property(get_step,None,None,"Volume step (numpy style)")

    # hyperslab-based functions

    def set_volume_range(self, rmin, rmax):
        """
        Specify expected volume range (for volume-based intensity normalization)
        :param rmin: expected minimal value
        :param rmax: expected maximal value
        :return:
        """
        if lib.minc2_set_volume_range(self._v,rmin,rmax) != lib.MINC2_SUCCESS:
            raise minc2_error()

    def save_hyperslab(self, buf, start=None):
        """
        save a hyperslab into minc2 volume
        :param buf: numpy.ndarray with information for the hyperslab
        :param start: array-like list describing the offset
        :return buf
        """
        if start is None:
            return self.save_complete_volume(buf)
        else:
            #import numpy as np
            data_type  = lib.MINC2_FLOAT
            store_type = buf.dtype.name

            ndims = self.ndim()
            _dims = self.representation_dims_()
            #dims = buf.shape # can't be used if we are using array with missing dimensions

            # TODO: unsqueeze array

            slab_start=ffi.new("int[]",ndims)
            slab_count=ffi.new("int[]",ndims)

            for i in range(ndims):
                if start[i] is not None:
                    if isinstance(start[i], list) or isinstance(start[i], tuple):
                        if len(start[i])==2:
                            slab_count[ndims-1-i]=start[i][1]-start[i][0]
                            slab_start[ndims-1-i]=start[i][0]
                        else: # assume it's the whole dimension
                            slab_count[ndims-1-i]=_dims[ndims-i-1].length
                            slab_start[ndims-1-i]=0
                    else: # assume it's a number
                        slab_count[ndims-1-i]=1 # _dims[ndims-i-1].length
                        slab_start[ndims-1-i]=start[i]
                else:
                    slab_count[ndims-1-i]=_dims[ndims-i-1].length
                    slab_start[ndims-1-i]=0
            
            data_type=minc2_file.__numpy_to_minc2[store_type]

            if lib.minc2_write_hyperslab( self._v, slab_start, slab_count,
                                          ffi.cast("void *", buf.ctypes.data),
                                          data_type ) != lib.MINC2_SUCCESS:
                raise minc2_error("Error writing hyperslab")
            return buf

    def save_hyperslab_t(self, buf, start=None):
        """
        save a hyperslab into minc2 volume (for torch.Tensor)
        :param buf: torch.Tensor with information for the hyperslab
        :param start: array-like list describing the offset
        :return buf
        """

        if start is None:
            return self.save_complete_volume(buf)
        else:
            import torch
            store_type=buf.type()
            assert(store_type in minc2_file.__torch_to_minc2)
            # TODO: verify dimensions of the array

            #data_type = minc2_file.store_type[store_type]
            data_type=minc2_file.__torch_to_minc2[store_type]

            ndims =self.ndim()
            #dims = buf.size()
            _dims = self.representation_dims_()

            slab_start=ffi.new("int[]",ndims)
            slab_count=ffi.new("int[]",ndims)

            for i in range(ndims):
                if start[i] is not None:
                    if isinstance(start[i], list) or isinstance(start[i], tuple):
                        if len(start[i])==2:
                            slab_count[ndims-1-i]=start[i][1]-start[i][0]
                            slab_start[ndims-1-i]=start[i][0]
                        else: # assume it's the whole dimension (?)
                            slab_count[ndims-1-i]=_dims[ndims-1-i]
                            slab_start[ndims-1-i]=0
                    else: #assume it's a number
                        slab_count[ndims-1-i]=1 # dims[i]
                        slab_start[ndims-1-i]=start[i]
                else: # use the whole dimension (?)
                    slab_count[ndims-1-i]=_dims[ndims-1-i]
                    slab_start[ndims-1-i]=0

            if lib.minc2_write_hyperslab( self._v, slab_start, slab_count,
                                          ffi.cast("void *", buf.storage().data_ptr()),
                                          data_type ) != lib.MINC2_SUCCESS:
                raise minc2_error("Error writing hyperslab")
            return buf

    def load_hyperslab(self, slab=None, data_type=None):
        """
        Load hyperslab into memory
        :param slab: array of format ((dim1_start[,dim1_stop]),(dim2_start[,dim2_stop]),....) describing the hyperslab to read
        :param data_type: requested numpy datatype
        :return: numpy.ndarray
        """
        if data_type is None:
            data_type=self.representation_dtype()
        if slab is None:
            return self.load_complete_volume(data_type)
        else:
            import numpy as np
            buf = None
            _dims = self.representation_dims_()
            ndims = self.ndim()
            dims = [None]*ndims

            slab_start = ffi.new("int[]",ndims)
            slab_count = ffi.new("int[]",ndims)

            for i in range(ndims):
                if len(slab)>i and slab[i] is not None:
                    if isinstance(slab[i], list) or isinstance(slab[i], tuple):
                        if len(slab[i]) == 2:
                            slab_count[ndims-1-i] = slab[i][1]-slab[i][0]
                            slab_start[ndims-1-i] = slab[i][0]
                        else:# -- assume it's the whole dimension
                            slab_count[ndims-1-i] = _dims[ndims-i-1].length
                            slab_start[ndims-1-i] = 0
                    else: # assume it's a number
                        slab_count[ndims-1-i] = 1
                        slab_start[ndims-1-i] = slab[i]
                else:
                    slab_count[ndims-1-i] = _dims[ndims-i-1].length
                    slab_start[ndims-1-i] = 0
                dims[i] = slab_count[ndims-1-i]

            dtype = None

            if data_type in minc2_file.__minc2_to_numpy:
                dtype=minc2_file.__minc2_to_numpy[data_type]
            elif data_type in minc2_file.__numpy_to_minc2:
                dtype=data_type
                data_type=minc2_file.__numpy_to_minc2[dtype]
            elif isinstance(data_type, np.dtype):
                dtype=data_type
                data_type=minc2_file.__numpy_to_minc2[dtype.name]
            else:
                raise minc2_error("Unsupported data type:"+repr(data_type))

            buf = np.empty(dims, dtype, 'C')

            if lib.minc2_read_hyperslab( self._v, slab_start, slab_count,
                                         ffi.cast("void *", buf.ctypes.data),
                                         data_type ) != lib.MINC2_SUCCESS:
                raise minc2_error("Error reading hyperslab")
            return buf

    def load_hyperslab_t(self, slab=None, data_type=None):
        """
        Load hyperslab into memory, for torch
        :param slab: array of format ((dim1_start[,dim1_stop]),(dim2_start[,dim2_stop]),....) describing the hyperslab to read
        :param data_type: requested numpy datatype
        :return: torch.Tensor
        """
        if data_type is None:
            data_type=self.representation_dtype_tensor()
        if slab is None:
            return self.load_complete_volume_tensor(data_type)
        else:
            import torch
            buf = None
            _dims = self.representation_dims()
            ndims = self.ndim()
            dims = [None]*ndims

            slab_start = ffi.new("int[]",ndims)
            slab_count = ffi.new("int[]",ndims)

            for i in range(ndims):
                if len(slab)>i and slab[i] is not None:
                    if isinstance(slab[i], list) or isinstance(slab[i], tuple):
                        if len(slab[i]) == 2:
                            slab_count[ndims-1-i] = slab[i][1]-slab[i][0]
                            slab_start[ndims-1-i] = slab[i][0]
                        else:# -- assume it's the whole dimension
                            slab_count[ndims-1-i] = _dims[ndims-i-1].length
                            slab_start[ndims-1-i] = 0
                    else: # assume it's a number
                        slab_count[ndims-1-i] = 1
                        slab_start[ndims-1-i] = slab[i]
                else:
                    slab_count[ndims-1-i] = _dims[ndims-i-1].length
                    slab_start[ndims-1-i] = 0
                dims[i] = slab_count[ndims-1-i]

            dtype=None

            if data_type in minc2_file.__minc2_to_torch:
                dtype=eval(minc2_file.__minc2_to_torch[data_type])
            elif data_type in minc2_file.__torch_to_minc2:
                #print("evaluate:{}".format(data_type))
                dtype=eval(data_type)
                data_type=minc2_file.__torch_to_minc2[data_type]
            else:
                raise minc2_error("Unsupported data type:"+repr(data_type))

            buf=dtype(*dims)

            if lib.minc2_read_hyperslab( self._v, slab_start, slab_count,
                                         ffi.cast("void *", buf.storage().data_ptr()),
                                         data_type ) != lib.MINC2_SUCCESS:
                raise minc2_error("Error reading hyperslab")
            return buf

    def _slices_to_slab(self, args):
        # assume it's either a number, or slice
        args = args if isinstance(args, tuple) else (args,)
        slab=[]
        _dims = self.representation_dims()
        ndims = self.ndim()

        if len(args)!=self.ndim():
            raise minc2_error("Unsupported number of dimensions")

        for i,a in enumerate(args):
            if isinstance(a, slice):
                if a.step is not None and a.step!=1:
                    raise minc2_error()
                slab+=[(a.start if a.start is not None  else 0,
                        a.stop  if a.stop  is not None  else _dims[ndims-i-1].length )]
            elif isinstance(a, int):
                slab+=[a]
            else:
                raise minc2_error("Unsupported dimension type:"+repr(a))
        return slab

    def __getitem__(self,s):
        """
        numpy-style interface for reading hyperslabs from volume
        :param s: slice information
        :return: numpy.ndarray
        """
        idx=self._slices_to_slab(s)
        return self.load_hyperslab(idx).squeeze()

    def __setitem__(self, s, val):
        """
        numpy-style interface for writing hyperslabs into volume
        :param s: slice information
        :param val:  numpy.ndarray
        :return:
        """
        idx=self._slices_to_slab(s)
        return self.save_hyperslab(val,idx)

class minc2_xfm:
    """
    MINC2 .xfm file object
    
    """
    # constants
    MINC2_XFM_LINEAR                 = lib.MINC2_XFM_LINEAR
    MINC2_XFM_THIN_PLATE_SPLINE      = lib.MINC2_XFM_THIN_PLATE_SPLINE
    MINC2_XFM_USER_TRANSFORM         = lib.MINC2_XFM_USER_TRANSFORM
    MINC2_XFM_CONCATENATED_TRANSFORM = lib.MINC2_XFM_CONCATENATED_TRANSFORM
    MINC2_XFM_GRID_TRANSFORM         = lib.MINC2_XFM_GRID_TRANSFORM


    def __init__(self, path=None):
        """
        create new xfm object
        :param path: file path
        """
        self._v=ffi.gc(lib.minc2_xfm_allocate0(),lib.minc2_xfm_destroy)
        if path is not None:
            self.open(path)

    def open(self, path):
        """
        Open existing .xfm file
        :param path: file path
        :return:
        """
        assert path is not None,"Provide xfm file"
        assert lib.minc2_xfm_open(self._v,to_bytes(path)) == lib.MINC2_SUCCESS

    def save(self, path):
        """
        Save information into file
        :param path:
        :return:
        """
        assert path is not None,"Provide xfm file"
        assert(lib.minc2_xfm_save(self._v,to_bytes(path)) == lib.MINC2_SUCCESS)

    def transform_point(self, xyz_in):
        """
        Apply transformation to coordinates
        :param xyz_in:  either 1d array or 2d array
        :return:  1d or 3d array
        """
        import numpy as np
        _xyz_in=np.asarray(xyz_in,'float64','C')
        if len(_xyz_in.shape)==1:
            xyz_out=np.empty([3],'float64','C')
            assert(lib.minc2_xfm_transform_point(self._v,ffi.cast("double *", _xyz_in.ctypes.data),ffi.cast("double *", xyz_out.ctypes.data))==lib.MINC2_SUCCESS)
            return xyz_out
        else:
            xyz_out=np.empty(_xyz_in.shape,'float64','C')
            assert(lib.minc2_xfm_transform_point_vec(self._v,_xyz_in.shape[0], 3, ffi.cast("double *", _xyz_in.ctypes.data), ffi.cast("double *", xyz_out.ctypes.data))==lib.MINC2_SUCCESS)
            return xyz_out

    def inverse_transform_point(self,xyz_in):
        """
        Apply inverse transformation to coordinates
        :param xyz_in:  either 1d array or 2d array
        :return:  1d or 3d array
        """
        import numpy as np
        _xyz_in=np.asarray(xyz_in,'float64','C')
        if len(_xyz_in.shape)==1:
            xyz_out=np.empty([3],'float64','C')
            assert(lib.minc2_xfm_inverse_transform_point(self._v,ffi.cast("double *", _xyz_in.ctypes.data),ffi.cast("double *", xyz_out.ctypes.data))==lib.MINC2_SUCCESS)
            return xyz_out
        else:
            xyz_out=np.empty(_xyz_in.shape,'float64','C')
            assert(lib.minc2_xfm_inverse_transform_point_vec(self._v,_xyz_in.shape[0], 3, ffi.cast("double *", _xyz_in.ctypes.data), ffi.cast("double *", xyz_out.ctypes.data))==lib.MINC2_SUCCESS)
            return xyz_out

    def invert(self):
        """
        invert transform
        :return:
        """
        assert(lib.minc2_xfm_invert(self._v)==lib.MINC2_SUCCESS)

    def get_n_concat(self):
        """
        get number of sub-transforms
        :return: integer
        """
        n=ffi.new("int*",0)
        assert(lib.minc2_xfm_get_n_concat(self._v,n)==lib.MINC2_SUCCESS)
        return n[0]

    def get_n_type(self,n=0):
        """
        Get n'th transform type
        :param n: transform number
        :return: transform type
        """
        t=ffi.new("int*",0)
        assert(lib.minc2_xfm_get_n_type(self._v,n,t)==lib.MINC2_SUCCESS)
        return t[0]

    def get_grid_transform(self,n=0):
        """
        Extract non-linear grid transform
        :param n: subtransform number
        :return: (grid file, inversion flag)
        """
        c_file=ffi.new("char**")
        inv=ffi.new("int*",0)
        assert(lib.minc2_xfm_get_grid_transform(self._v,n,inv,c_file)==lib.MINC2_SUCCESS)
        _file=ffi.string(c_file[0])
        lib.free(c_file[0])
        return (_file,inv[0]!=0)

    def get_linear_transform(self,n=0):
        """
        Extract affine transform
        :param n: subtransform number
        :return: numpy.ndarray([4,4]) - matrix describing affine transform
        """
        import numpy as np
        _mat=np.empty([4,4],'float64','C')
        assert(lib.minc2_xfm_get_linear_transform(self._v,n,ffi.cast("double *", _mat.ctypes.data))==lib.MINC2_SUCCESS)
        return _mat

    def get_linear_transform_param(self,n=0,center=None):
        """
        Extract affine transform parameters
        :param n: subtransform number
        :param center: rotation center
        :return: minc2_transform_parameters linear transformation parameters
        """
        import numpy as np
        par=minc2_transform_parameters()

        if center is not None:
            par.center=np.asarray(center,'float64','C')
        # TODO: detect if extraction of parameters failed

        if lib.minc2_xfm_extract_linear_param(self._v,n,
                ffi.cast("double *", par.center.ctypes.data),
                ffi.cast("double *", par.translations.ctypes.data),
                ffi.cast("double *", par.scales.ctypes.data),
                ffi.cast("double *", par.shears.ctypes.data),
                ffi.cast("double *", par.rotations.ctypes.data)
            )!=lib.MINC2_SUCCESS:
            par.invalid=True
        return par


    def append_linear_transform(self, par):
        """
        Concatenate linear transformation
        :param par: if numpy.ndarray - should be 4x4 matrix describing affine transform, or it can be minc2_transform_parameters
        :return:
        """
        import numpy as np
        if isinstance(par,np.ndarray): # assume it's a matrix
            _mat=np.asarray(par,'float64','C')
            assert(lib.minc2_xfm_append_linear_transform(self._v,ffi.cast("double *", _mat.ctypes.data))==lib.MINC2_SUCCESS)
        else: # must be an object with parameters
            _par=minc2_transform_parameters()

            _par.center=np.asarray(par.center,'float64','C')
            _par.translations=np.asarray(par.translations,'float64','C')
            _par.scales=np.asarray(par.scales,'float64','C')
            _par.shears=np.asarray(par.shears,'float64','C')
            _par.rotations=np.asarray(par.rotations,'float64','C')

            assert(lib.minc2_xfm_append_linear_param(self._v,
                ffi.cast("double *", _par.center.ctypes.data),
                ffi.cast("double *", _par.translations.ctypes.data),
                ffi.cast("double *", _par.scales.ctypes.data),
                ffi.cast("double *", _par.shears.ctypes.data),
                ffi.cast("double *", _par.rotations.ctypes.data)
                )==lib.MINC2_SUCCESS)
        return self

    def append_grid_transform(self,grid_file,inv=False):
        """
        Concatenate nonlinear grid transform
        :param grid_file: grid file path
        :param inv:  inversion flag
        :return:
        """
        assert(lib.minc2_xfm_append_grid_transform(self._v,to_bytes(grid_file),inv)==lib.MINC2_SUCCESS)
        return self

    def concat_xfm(self,another):
        """
        Cobcatenate another minc2_xfm transform
        :param another: another minc2_xfm object
        :return:
        """
        assert(lib.minc2_xfm_concat_xfm(self._v,another._v)==lib.MINC2_SUCCESS)


class minc2_tags:
    """
    MINC2 tag object
    """
    def __init__(self, path=None,n_volumes=1):
        """
        Create object
        :param path:  existing .tag file
        :param n_volumes: number of expected volumes
        """
        import numpy as np
        self.n_volumes=n_volumes
        self.tag=[]
        self.weights=None
        self.structure_ids=None
        self.patient_ids=None
        self.labels=None

        if path is not None:
            self.load(path)

    def __len__(self):
        """
        number of tags or tag pairs
        :return: integer
        """
        if self.tag[0] is not None:
            return self.tag[0].shape[0]
        else:
            return 0

    def load(self, path):
        """
        load tags from a file
        :param path: file path
        """
        assert path is not None,"Provide tag file"
        _t=ffi.gc(lib.minc2_tags_allocate0(),lib.minc2_tags_free)
        assert(lib.minc2_tags_load(_t,to_bytes(path))==lib.MINC2_SUCCESS)
        import numpy as np
        # convert internal data structure to numpy arrays
        self.n_volumes=_t.n_volumes
        self.tag=[]

        self.tag+=[np.empty((_t.n_tag_points,3),'float64','C')]
        ffi.memmove(ffi.cast("double *",self.tag[0].ctypes.data), _t.tags_volume1, _t.n_tag_points*3*ffi.sizeof('double'))
        if self.n_volumes>1:
            self.tag+=[np.empty((_t.n_tag_points,3),'float64','C')]
            ffi.memmove( ffi.cast("double *",self.tag[1].ctypes.data), _t.tags_volume2, _t.n_tag_points*3*ffi.sizeof('double'))
        #
        if _t.weights is not None:
            self.weights=np.empty(_t.n_tag_points,'float64','C')
            ffi.memmove( ffi.cast("double *",self.weights.ctypes.data), _t.weights, _t.n_tag_points*ffi.sizeof('double'))
        else:
            self.weights=None
            
        if _t.structure_ids is not None:
            self.structure_ids=np.zeros(_t.n_tag_points, np.int_,'C')
            self.structure_ids[:]=np.frombuffer(ffi.buffer(_t.structure_ids, _t.n_tag_points*ffi.sizeof('int')), dtype='int32')
        else:
            self.structure_ids=None

        if _t.patient_ids is not None:
            self.patient_ids=np.empty(_t.n_tag_points,np.int_,'C')
            self.patient_ids[:]=np.frombuffer(ffi.buffer(_t.patient_ids, _t.n_tag_points*ffi.sizeof('int')), dtype='int32')
        else:
            self.patient_ids=None

        if _t.labels is not None:
            self.labels=[]
            for i in range(_t.n_tag_points):
                self.labels.append(ffi.string(_t.labels[i]))
        else:
            self.labels=None

    def save(self,path):
        """
        Save object into .tag file
        :param path: file path
        :return:
        """
        assert path is not None,"Provide tag file"
        _t=ffi.gc(lib.minc2_tags_allocate0(),lib.minc2_tags_free)

        import numpy as np
        import weakref
        # convert numpy arrays to internal data structure 
        assert self.n_volumes==1 or self.n_volumes==2 
        assert self.n_volumes == len(self.tag)

        if self.tag[0] is not None:
            weakdict = weakref.WeakKeyDictionary()

            shape=self.tag[0].shape
            assert shape[1]==3
            if self.n_volumes==2 :
                assert self.tag[1].shape == self.tag[0].shape

            assert lib.minc2_tags_init(_t,shape[0],self.n_volumes,1 if self.weights is not None else 0,1 if self.structure_ids is not None else 0,1 if self.patient_ids is not None else 0,1 if self.labels is not None else 0)==lib.MINC2_SUCCESS

            ffi.memmove(_t.tags_volume1, ffi.cast("double *",self.tag[0].ctypes.data),  _t.n_tag_points*3*ffi.sizeof('double'))
            if self.n_volumes>1:
                ffi.memmove( _t.tags_volume2, ffi.cast("double *",self.tag[1].ctypes.data),  _t.n_tag_points*3*ffi.sizeof('double'))
            #
            if self.weights is not None:
                ffi.memmove( _t.weights, ffi.cast("double *",self.weights.ctypes.data), _t.n_tag_points*ffi.sizeof('double'))
                
            if self.structure_ids is not None:
                np.frombuffer(ffi.buffer(_t.structure_ids, _t.n_tag_points*ffi.sizeof('int')), dtype='int32')[:] = \
                    self.structure_ids[:]

            if self.patient_ids is not None:
                np.frombuffer(ffi.buffer(_t.patient_ids, _t.n_tag_points*ffi.sizeof('int')), dtype='int32')[:] = \
                    self.patient_ids[:]

            weakdict[_t] = []
            if self.labels is not None:
                for i in range(shape[0]):
                    ss=ffi.new("char[]",to_bytes(self.labels[i]))
                    _t.labels[i] = ss # WARNING: have to make sure C code doesn't try to free these strings
                    weakdict[_t].append(ss)
                    
            # save to disk
            assert lib.minc2_tags_save(_t,to_bytes(path))==lib.MINC2_SUCCESS

            if self.labels is not None:
                  for i in range(shape[0]):
                      _t.labels[i] = ffi.cast("void *", 0)
            
            print("tag labels should be de-allocated now")
        else: # it's empty
            pass
        assert lib.minc2_tags_save(_t,to_bytes(path))==lib.MINC2_SUCCESS

class minc2_input_iterator:
    """
    minc2 input file iterator, can work with multiple input files
    and return a vector of values
    """
    def __init__(self, files=None, data_type=None):
        """
        initialize iterator
        :param files: list of minc2 files (either minc2 objects or paths)
        :param data_type: expected data type
        """
        self._i = ffi.gc(lib.minc2_iterator_allocate0(),lib.minc2_iterator_free)
        self._handles=[]
        self._dtype=None
        self._last=False
        if files is not None:
            self.open(files,data_type=data_type)
    
    def open(self,files,data_type=None):
        """
        open files
        :param files: list of minc2 files (either minc2 objects or paths)
        :param data_type: expected data type
        """
        if isinstance(files, six.string_types):
            files=(files,)

        for f in files:
            _h=ffi.gc(lib.minc2_allocate0(), lib.minc2_destroy)
            if lib.minc2_open(_h, to_bytes(f))!=lib.MINC2_SUCCESS:
                raise minc2_error("Can't open file:"+f)
            self._handles+=[_h]

        if data_type is None:
            _dtype = ffi.new("int*",0)
            if lib.minc2_data_type(self._handles[0],_dtype)!=lib.MINC2_SUCCESS:
                raise minc2_error("Error getting representation type")
            data_type=_dtype[0]
        
        import numpy as np
        self._dtype = minc2_file.minc2_to_numpy[data_type]
        self._val =  np.empty(len(self._handles), self._dtype, 'C')
        lib.minc2_multi_iterator_input_start(self._i,self._handles,data_type,len(self._handles))
        self._last=False

    def close(self):
        """
        Close all input file, close iterator
        :return:
        """
        for i in self._handles:
            lib.minc2_close(i)
        self._handles=[]

    def dim(self):
        """
        number of files being used
        :return:
        """
        return len(self._handles)

    def __iter__(self):
        """
        iterate
        :return:
        """
        return self

    def next(self):
        """
        read and advance
        """
        # minc2_iterator_get_values(input_minc_it,&voxels[0]);
        #_val= np.empty(shape, dtype, 'C')
        if self._last:
            raise StopIteration()

        lib.minc2_iterator_get_values(self._i,ffi.cast("void *", self._val.ctypes.data))

        self._last=lib.minc2_iterator_next(self._i)!=lib.MINC2_SUCCESS

        return self._val
            

    def __next__(self):
        return self.next()

    def val(self):
        """
        read the current value without advancing
        :return nunpy.ndarray
        """
        # TODO: make sure we can read
        lib.minc2_iterator_get_values(self._i,ffi.cast("void *", self._val.ctypes.data))
        return self._val

    def __del__(self):
        self.close()


class minc2_output_iterator:
    """
    minc2 output file iterator, can work with multiple output files
    and write a vector of values
    doesn't quite follow python semantics
    """
    def __init__(self,files=None, reference=None, data_type=None, store_type=None,slice_scaling=None, global_scaling=None,):
        self._i = ffi.gc(lib.minc2_iterator_allocate0(), lib.minc2_iterator_free)
        self._handles=[]
        self._last=False
        
        if files is not None:
            self.open(files, reference=reference, 
                      data_type=data_type, store_type=store_type,
                      slice_scaling=slice_scaling, global_scaling=global_scaling)

    def open(self,files,reference=None, data_type=None, store_type=None,slice_scaling=None, global_scaling=None,):
        self._handles=[]

        if isinstance(files, six.string_types):
            files=(files,)

        for f in files:
            _h=ffi.gc(lib.minc2_allocate0(), lib.minc2_destroy)
            self._handles+=[_h]

        _ref=None
        _dims = None

        if isinstance(reference, minc2_input_iterator ):
            reference = minc2_file(handle=reference._handles[0])

        elif isinstance(reference, six.string_types):
            reference = minc2_file(reference)

        if isinstance(reference, minc2_file ):
            if store_type is None:
                store_type = minc2_file.numpy_to_minc2[reference.store_dtype()]
            if data_type is None:
                data_type = minc2_file.numpy_to_minc2[reference.representation_dtype()]
            _dims=reference.store_dims()

        elif isinstance(reference,list) or isinstance(reference,tuple):
            # assume it's dimensions definition
            _dims=reference
            reference=None
        
        # default
        if data_type is None and store_type is None:
            data_type=store_type=lib.MINC2_FLOAT
        elif store_type is None:
            data_type=store_type
        elif data_type is None:
            store_type=data_type

        # create dimensions for all volumes
        __dims = ffi.new("struct minc2_dimension[]", len(_dims)+1)
        for i,j in enumerate(_dims):
            if isinstance(j, minc2_dim):
                __dims[i].id=j.id
                __dims[i].length=j.length
                __dims[i].start=j.start
                __dims[i].step=j.step
                __dims[i].have_dir_cos=j.have_dir_cos
                if j.have_dir_cos: 
                    ffi.memmove(__dims[i].dir_cos, ffi.cast("double [3]", j.dir_cos.ctypes.data ), 3*ffi.sizeof('double'))
            else:
                __dims[i]=j
        __dims[len(_dims)]={ 'id': lib.MINC2_DIM_END }

        for f,h in zip(files,self._handles):
            if lib.minc2_define(h, __dims, store_type, data_type)!=lib.MINC2_SUCCESS:
                raise minc2_error("Error defining new minc file ")

            if slice_scaling is not None or global_scaling is not None:
                _slice_scaling=0
                _global_scaling=0

                if slice_scaling: _slice_scaling=1
                if global_scaling: _global_scaling=1

                if lib.minc2_set_scaling(self._v,_global_scaling,_slice_scaling )!=lib.MINC2_SUCCESS:
                    raise minc2_error()

            if lib.minc2_create(h, to_bytes(f) )!=lib.MINC2_SUCCESS:
                raise minc2_error("Error creating file:"+f)

        self._dtype = minc2_file.minc2_to_numpy[data_type]
        import numpy as np
        self._val =  np.empty( len(self._handles), self._dtype, 'C')

        lib.minc2_multi_iterator_output_start(self._i, self._handles, data_type, len(self._handles))
        self._last=False

    def close(self):
        # close all input files
        for i in self._handles:
            lib.minc2_close(i)
        self._handles=[]

    def __del__(self):
        self.close()

    def dim(self):
        return len(self._handles)

    def __iter__(self):
        return self

    def next(self, value=None):
        """
        write value and adance
        """
        if self._last:
            raise StopIteration()

        if value is not None:
            self.set_value(value)

        self._last=lib.minc2_iterator_next(self._i)!=lib.MINC2_SUCCESS
            
    def __next__(self):
        return self.next()

    def set_value(self,v):
        self._val[:]=v
        lib.minc2_iterator_put_values(self._i,ffi.cast("void *", self._val.ctypes.data))
        return v

# kate: indent-width 4; replace-tabs on; remove-trailing-space on; hl python; show-tabs on
