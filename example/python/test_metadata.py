from minc2_simple import minc2_file
from minc2_simple import minc2_xfm

import sys
import numpy as np


if __name__ == "__main__":
    
    if len(sys.argv)<3:
        print("Usage: {} input.mnc output.mnc".format(sys.argv[0]))
        sys.exit(1)
    
    infile=sys.argv[1]
    outfile=sys.argv[2]
    
    m=minc2_file(infile)
    o=minc2_file()
    # will create file with same dimensions
    o.define(m.store_dims(), minc2_file.MINC2_BYTE, minc2_file.MINC2_FLOAT)
    
    print("Will create new volume...")
    o.create(outfile)
    
    meta=m.metadata()
    
    print("Metadata:")
    print(repr(meta))
    
    print("History:")
    print(m.read_attribute("","history"))
    
    print("Writing metadata")
    o.write_metadata(meta)

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
