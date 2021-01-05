from minc2_simple import minc2_file
from minc2_simple import minc2_xfm

import sys
import numpy as np


if __name__ == "__main__":
    #infile='/home/vfonov/data01/deep-seg/hc/jens_ec_lib_identity/local_model.mnc'
    infile='/opt/minc/share/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc'
    #m=minc2_file("/opt/minc/share/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc")
    m=minc2_file(infile)
    print("dims={}".format(m.ndim()))
    dims=m.store_dims()

    for i in range(m.ndim()): 
        print('Dimension {} length:{} id:{} start:{} step:{}'.format(i,dims[i].length, dims[i].id, dims[i].start, dims[i].step))

    o=minc2_file()
    # will create file with same dimensions

    my_dims=[
        { 'id':minc2_file.MINC2_DIM_X, 'length':193, 'start':96.0,   'step':-1.0},
        { 'id':minc2_file.MINC2_DIM_Y, 'length':229, 'start':-132.0, 'step':1.0},
        { 'id':minc2_file.MINC2_DIM_Z, 'length':193, 'start':-78.0,  'step':1.0}
    ]
    print("Will define new volume...")
    o.define(dims, minc2_file.MINC2_BYTE, minc2_file.MINC2_FLOAT)
    print("Will create new volume...")
    o.create('test_out.mnc')
    
    meta=m.metadata()
    
    print("Metadata:")
    print(repr(meta))
    
    print("History:")
    print(m.read_attribute("","history"))
    
    print("Copying metadata")
    o.copy_metadata(m)

    # going to read and write in standard order (XYZ)
    m.setup_standard_order()
    o.setup_standard_order()
    print("Loading from file...")

    # load into a c buffer
    data=m.load_complete_volume(minc2_file.MINC2_FLOAT)

    print("loaded array {} of size :{}".format(data.dtype,data.shape))

    # save from buffer to volume
    o.save_complete_volume(data)

    # not strictly needed , but will flush the data to disk immedeately
    o.close()
    # not strictly needed  either, the file will be close by garbage collection
    m.close()


    # working with slices, numpy style
    m=minc2_file(infile)
    o=minc2_file()
    o.define(m.store_dims(), minc2_file.MINC2_BYTE, minc2_file.MINC2_FLOAT)
    o.create('test_out2.mnc')

    m.setup_standard_order()
    o.setup_standard_order()
    d=m.representation_dims()

    _slice=m[d[0].length//2,:,:]

    o[0,:,:]=_slice
    
    o.close()
    m.close()

    

    
    
