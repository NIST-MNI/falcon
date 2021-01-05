# based on https://github.com/Mouse-Imaging-Centre/pyminc/blob/master/test/generatorTests.py

import unittest
import numpy as N
import os
import subprocess
import tempfile

from minc2_simple import minc2_file,minc2_xfm,minc2_error

def setUpModule():
    global outputFilename,emptyFilename,newFilename,inputFile_byte,inputFile_short,inputFile_int
    global inputFile_float,inputFile_double,inputFile_ubyte,inputFile_ushort,inputFile_uint
    global inputVector,input3DdirectionCosines
    
    
    outputFilename = tempfile.NamedTemporaryFile(prefix="test-out-", suffix=".mnc").name
    emptyFilename = tempfile.NamedTemporaryFile(prefix="test-empty-", suffix=".mnc").name
    newFilename = tempfile.NamedTemporaryFile(prefix="test-new-volume-", suffix=".mnc").name

    inputFile_byte = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_byte, '-osigned', '-obyte', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_short = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_short, '-osigned', '-oshort', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_int = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_int, '-oint', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_float = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_float, '-ofloat', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_double = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_double, '-odouble', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_ubyte = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_ubyte, '-ounsigned', '-obyte', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_ushort = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_ushort, '-ounsigned', '-oshort', '-input', '/dev/urandom', '100', '150', '125'])

    inputFile_uint = tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputFile_uint, '-ounsigned', '-oint', '-input', '/dev/urandom', '100', '150', '125'])


    inputVector = tempfile.NamedTemporaryFile(prefix="test-vector-", suffix=".mnc").name
    subprocess.check_call(['rawtominc', inputVector, '-input', '/dev/urandom', '3', '100', '150', '125',
                        '-dimorder', 'vector_dimension,xspace,yspace,zspace'])

    input3DdirectionCosines = tempfile.NamedTemporaryFile(prefix="test-3d-direction-cosines", suffix=".mnc").name
    subprocess.check_call(['rawtominc', input3DdirectionCosines, '-input', '/dev/urandom',
                           '100', '150', '125',
                        '-xstep', '1.0', '-ystep', '1.2', '-zstep', '0.9',
                        '-xstart', '-50.0', '-ystart', '-100.0', '-zstart', '-20.0',
                        '-xdircos',  "0.930537594119596",  "0.130822050579164",  "0.342031251514211" ,
                        '-ydircos', "-0.195843899594348",  "0.966963987997126",  "0.163174179662017",
                        '-zdircos', "-0.309385127921852", "-0.218824442010364",  "0.925417044472184"])


def tearDownModule():
    os.remove(inputFile_byte)
    os.remove(inputFile_short)
    os.remove(inputFile_int)
    os.remove(inputFile_float)
    os.remove(inputFile_double)
    os.remove(inputFile_ubyte)
    os.remove(inputFile_ushort)
    os.remove(inputFile_uint)
    os.remove(inputVector)
    os.remove(input3DdirectionCosines)
    
    if os.path.exists(outputFilename):
        os.remove(outputFilename)
    if os.path.exists(newFilename):
        os.remove(newFilename)
    

class minc2_file_io_numpy(unittest.TestCase):
    """test the minc2_file reading using numpy"""
    def testFromFileError(self):
        """attempting to load a garbage file should raise exception"""
        self.assertRaises(minc2_error, minc2_file, "garbage.mnc")
    def testFromFileDataTypeByte(self):
        """ensure byte data is read as float by default"""
        v = minc2_file(inputFile_byte)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeShort(self):
        """ensure short data is read as float by default"""
        v = minc2_file(inputFile_short)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeInt(self):
        """ensure int data is read as float by default"""
        v = minc2_file(inputFile_int)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeFloat(self):
        """ensure float data is read as float by default"""
        v = minc2_file(inputFile_float)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeDouble(self):
        """ensure double data is read as float"""
        v = minc2_file(inputFile_double)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float64')
    def testFromFileDataTypeUByte(self):
        """ensure unsigned byte data is read as float by default"""
        v = minc2_file(inputFile_ubyte)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeUShort(self):
        """ensure unsigned short data is read as float by default"""
        v = minc2_file(inputFile_ushort)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataTypeUInt(self):
        """ensure unsigned int data is read as float by default"""
        v = minc2_file(inputFile_uint)
        dt = v.representation_dtype()
        v.close()
        self.assertEqual(dt, 'float32')
    def testFromFileDataByte(self):
        """ensure that byte data is read correct with a precision of 8 decimals on a call to average()"""
        v = minc2_file(inputFile_byte)
        a = N.average(v.load_complete_volume('float64'))
        v.close()
        pipe = os.popen("mincstats -mean -quiet %s" % inputFile_byte, "r")
        output = float(pipe.read())
        pipe.close()
        self.assertAlmostEqual(a, output, 8)
    def testFromFileDataShort(self):
        """ensure that short data is read correct with a precision of 8 decimals on a call to average()"""
        v = minc2_file(inputFile_short)
        a = N.average(v.load_complete_volume('float64'))
        v.close()
        pipe = os.popen("mincstats -mean -quiet %s" % inputFile_short, "r")
        output = float(pipe.read())
        pipe.close()
        self.assertAlmostEqual(a, output, 8)
    def testFromFileDataInt(self):
        """ensure that int data is read correct with a precision of 8 decimals on a call to aveage()"""
        v = minc2_file(inputFile_int)
        a = N.average(v.load_complete_volume('float64'))
        v.close()
        pipe = os.popen("mincstats -mean -quiet %s" % inputFile_int, "r")
        output = float(pipe.read())
        pipe.close()
        self.assertAlmostEqual(a, output, 8)
    def testFromFileDataFloat(self):
        """ensure that float data is read correct with a precision of 8 decimals on a call to aveage()"""
        v = minc2_file(inputFile_float)
        a = N.average(v.load_complete_volume('float64'))
        v.close()
        pipe = os.popen("mincstats -mean -quiet %s" % inputFile_float, "r")
        output = float(pipe.read())
        pipe.close()
        self.assertAlmostEqual(a, output, 8)
    def testFromFileDataDouble(self):
        """ensure that double data is read correct with a precision of 8 decimals on a call to aveage()"""
        v = minc2_file(inputFile_double) 
        a = N.average(v.data)
        v.close()
        pipe = os.popen("mincstats -mean -quiet %s" % inputFile_double, "r")
        output = float(pipe.read())
        pipe.close()
        self.assertAlmostEqual(a, output, 8)
    def testDims(self):
        """Check data dimensions are correct"""
        v = minc2_file(inputFile_double) 
        v.setup_standard_order()
        dims=v.store_dims()

        self.assertEqual(len(dims), 3)
        # '100', '150', '125'    
        self.assertEqual(dims[0].id, minc2_file.MINC2_DIM_X) ## X
        self.assertEqual(dims[0].length, 125 )
        self.assertEqual(dims[1].id, minc2_file.MINC2_DIM_Y) ## Y
        self.assertEqual(dims[1].length, 150 )
        self.assertEqual(dims[2].id, minc2_file.MINC2_DIM_Z) ## X
        self.assertEqual(dims[2].length, 100 )

        
try:
    import torch # this is going to work only if torch is present
    
    class minc2_file_io_torch(unittest.TestCase):
        """test the minc2_file reading using pytorch tensor"""
        def testFromFileDataTypeByte(self):
            """ensure byte data is read as float by default"""
            v = minc2_file(inputFile_byte)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeShort(self):
            """ensure short data is read as float by default"""
            v = minc2_file(inputFile_short)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeInt(self):
            """ensure int data is read as float by default"""
            v = minc2_file(inputFile_int)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeFloat(self):
            """ensure float data is read as float by default"""
            v = minc2_file(inputFile_float)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeDouble(self):
            """ensure double data is read as float"""
            v = minc2_file(inputFile_double)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.DoubleTensor')
        def testFromFileDataTypeUByte(self):
            """ensure unsigned byte data is read as float by default"""
            v = minc2_file(inputFile_ubyte)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeUShort(self):
            """ensure unsigned short data is read as float by default"""
            v = minc2_file(inputFile_ushort)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataTypeUInt(self):
            """ensure unsigned int data is read as float by default"""
            v = minc2_file(inputFile_uint)
            dt = v.representation_dtype_tensor()
            v.close()
            self.assertEqual(dt, 'torch.FloatTensor')
        def testFromFileDataByte(self):
            """ensure that byte data is read correct with a precision of 8 decimals on a call to average()"""
            v = minc2_file(inputFile_byte)
            a = v.load_complete_volume_tensor('torch.DoubleTensor').mean().item()
            v.close()
            pipe = os.popen("mincstats -mean -quiet %s" % inputFile_byte, "r")
            output = float(pipe.read())
            pipe.close()
            self.assertAlmostEqual(a, output, 8)
        def testFromFileDataShort(self):
            """ensure that short data is read correct with a precision of 8 decimals on a call to average()"""
            v = minc2_file(inputFile_short)
            a = v.load_complete_volume_tensor('torch.DoubleTensor').mean().item()
            v.close()
            pipe = os.popen("mincstats -mean -quiet %s" % inputFile_short, "r")
            output = float(pipe.read())
            pipe.close()
            self.assertAlmostEqual(a, output, 8)
        def testFromFileDataInt(self):
            """ensure that int data is read correct with a precision of 8 decimals on a call to average()"""
            v = minc2_file(inputFile_int)
            a = v.load_complete_volume_tensor('torch.DoubleTensor').mean().item()
            v.close()
            pipe = os.popen("mincstats -mean -quiet %s" % inputFile_int, "r")
            output = float(pipe.read())
            pipe.close()
            self.assertAlmostEqual(a, output, 8)
        def testFromFileDataFloat(self):
            """ensure that float data is read correct with a precision of 8 decimals on a call to average()"""
            v = minc2_file(inputFile_float)
            a = v.load_complete_volume_tensor('torch.DoubleTensor').mean().item()
            v.close()
            pipe = os.popen("mincstats -mean -quiet %s" % inputFile_float, "r")
            output = float(pipe.read())
            pipe.close()
            self.assertAlmostEqual(a, output, 8)
        def testFromFileDataDouble(self):
            """ensure that double data is read correct with a precision of 8 decimals on a call to aveage()"""
            v = minc2_file(inputFile_double) 
            a = v.tensor.mean().item()
            v.close()
            pipe = os.popen("mincstats -mean -quiet %s" % inputFile_double, "r")
            output = float(pipe.read())
            pipe.close()
            self.assertAlmostEqual(a, output, 8)
except ImportError:
    pass

class TestWriteFileDataTypes(unittest.TestCase):
    ############################################################################
    # volumeFromDescription
    ############################################################################
    def testWriteDataAsByte(self):
        """ensure that a volume created by volumeFromDescription as byte is written out as such"""
        # TODO
        pass
    def testWriteDataAsShort(self):
        """ensure that a volume created by volumeFromDescription as short is written out as such"""
        # TODO
        pass
    def testWriteDataAsInt(self):
        """ensure that a volume created by volumeFromDescription as int is written out as such"""
        # TODO
        pass
    def testWriteDataAsFloat(self):
        """ensure that a volume created by volumeFromDescription as float is written out as such"""
        # TODO
        pass
    def testWriteDataAsDouble(self):
        """ensure that a volume created by volumeFromDescription as double is written out as such"""
        # TODO
        pass
    def testWriteDataAsUByte(self):
        """ensure that a volume created by volumeFromDescription as unsigned byte is written out as such"""
        # TODO
        pass
    def testWriteDataAsUShort(self):
        """ensure that a volume created by volumeFromDescription as unsigned short is written out as such"""
        # TODO
        pass
    def testWriteDataAsUInt(self):
        """ensure that a volume created by volumeFromDescription as unsigned int is written out as such"""
        # TODO
        pass

class minc2_file_hyperslabs_numpy(unittest.TestCase):
    """test getting and setting of hyperslabs"""
    def testGetHyperslab(self):
        """hyperslab should be same as slice from data array"""
        
        inputFile=inputFile_ushort
        #inputFile='/export01/data/vfonov/src1/minc2-simple/python/test_icbm.mnc'
        
        v = minc2_file(inputFile)
        v.setup_standard_order()
        sliceFromData_x = v.data[10,:,:]
        sliceFromData_y = v.data[:,10,:]
        sliceFromData_z = v.data[:,:,10]
        v.close()
        
        b = minc2_file(inputFile)
        b.setup_standard_order()
        hyperslab_x = b.load_hyperslab( [10, None, None] ).squeeze()
        hyperslab_y = b.load_hyperslab( [None, 10, None] ).squeeze()
        hyperslab_z = b.load_hyperslab( [None, None, 10] ).squeeze()
        b.close()

        self.assertEqual(N.average((sliceFromData_x-hyperslab_x)**2),0.0)
        self.assertEqual(N.average((sliceFromData_y-hyperslab_y)**2),0.0)
        self.assertEqual(N.average((sliceFromData_z-hyperslab_z)**2),0.0)

    def testSlicingGet(self):
        """volume slice should be same as slice from data array"""
        
        inputFile=inputFile_ushort
        
        v = minc2_file(inputFile)
        v.setup_standard_order()
        sliceFromData_x = v.data[10,:,:]
        sliceFromData_y = v.data[:,10,:]
        sliceFromData_z = v.data[:,:,10]
        v.close()
        
        b = minc2_file(inputFile)
        b.setup_standard_order()
        hyperslab_x = b[10, :, :]
        hyperslab_y = b[:, 10, :]
        hyperslab_z = b[:, :, 10]
        b.close()

        self.assertEqual(N.average((sliceFromData_x-hyperslab_x)**2),0.0)
        self.assertEqual(N.average((sliceFromData_y-hyperslab_y)**2),0.0)
        self.assertEqual(N.average((sliceFromData_z-hyperslab_z)**2),0.0)
        
    def testSetHyperslabFloat(self):
        """setting hyperslab should change underlying volume (float)"""
        
        # read some data from somwhere
        v  = minc2_file(inputFile_ushort)
        dims=v.store_dims()
        v.setup_standard_order()
        # use the whole volume
        hyperslab_a = v.data[10,:,:]
        hyperslab_a_ = v[10,:,:]
        v.close()

        print("Hyperslab:", hyperslab_a.shape )
        print("Hyperslab2:", hyperslab_a_.shape )
        print("dims:",dims)
        
        v2 = minc2_file()
        v2.define(dims,'float32','float32')
        v2.create(outputFilename)
        v2.setup_standard_order()
        
        # because we are saving float32 , we don't need slice normalization
        v2.save_hyperslab(hyperslab_a,   [10, None, None] )
        v2.close()

        v3  = minc2_file(outputFilename)
        hyperslab_b = v3.load_hyperslab( [10, None, None] )

        print(N.average((hyperslab_a-hyperslab_b)**2))

        self.assertEqual(N.average((hyperslab_a-hyperslab_b)**2),0.0)
        v3.close()

    def testSetSliceFloat(self):
        """volume slice setting should change underlying volume (float)"""
        
        # read some data from somwhere
        v  = minc2_file(inputFile_ushort)
        dims = v.store_dims()
        v.setup_standard_order()
        hyperslab_a = v.data[10,:,:]
        v.close()
        
        v2 = minc2_file()
        v2.define(dims,'float32','float32')
        v2.create(outputFilename)
        v2.setup_standard_order()

        # because we are saving float32 , we don't need slice normalization
        v2[10,:,:] = hyperslab_a
        v2.close()
        
        v3  = minc2_file(inputFile_ushort)
        v3.setup_standard_order()
        hyperslab_b = v3[10,:,:]
        v3.close()

        self.assertEqual(N.average((hyperslab_a-hyperslab_b)**2),0.0)

    def testSetHyperslabShort(self):
        """setting hyperslab should change underlying volume (short)"""
        
        # read some data from somwhere
        v  = minc2_file(inputFile_ushort)
        dims=v.store_dims()
        v.setup_standard_order()
        hyperslab_a = v.load_hyperslab( [10, None, None] )
        
        # try with normalization
        v2 = minc2_file()
        v2.define(dims,'uint16','float32') # , global_scaling=True
        v2.create(outputFilename)
        v2.set_volume_range(N.min(hyperslab_a),N.max(hyperslab_a))
        v2.setup_standard_order()
        
        # have to set slice normalization
        v2.save_hyperslab(hyperslab_a,   [10,None,None] )
        hyperslab_b = v2.load_hyperslab( [10, None, None] )
        self.assertAlmostEqual(N.average((hyperslab_a-hyperslab_b)**2),0.0,8)
        v2.close()
        v.close()
        
        
    def testHyperslabArray(self):
        """hyperslab should be reinsertable into volume"""
        if False:
            v = minc2_file(inputFile_ushort)
            v2 = minc2_file()
            v2.create(outputFilename)
            v2.close()
            v.close()

    try: # run tests if torch is present
        import torch
        
        class minc2_file_hyperslabs_torch(unittest.TestCase):
            """test getting and setting of hyperslabs"""
            def testGetHyperslab(self):
                """hyperslab should be same as slice from data array"""
                
                inputFile=inputFile_ushort
                #inputFile='/export01/data/vfonov/src1/minc2-simple/python/test_icbm.mnc'
                
                v = minc2_file(inputFile)
                v.setup_standard_order()
                sliceFromData_x = v.data[10,:,:]
                sliceFromData_y = v.data[:,10,:]
                sliceFromData_z = v.data[:,:,10]
                v.close()
                
                b = minc2_file(inputFile)
                b.setup_standard_order()
                hyperslab_x = b.load_hyperslab_t( [10, None, None] ).squeeze()
                hyperslab_y = b.load_hyperslab_t( [None, 10, None] ).squeeze()
                hyperslab_z = b.load_hyperslab_t( [None, None, 10] ).squeeze()
                b.close()

                self.assertEqual(torch.mean((sliceFromData_x-hyperslab_x)**2),0.0)
                self.assertEqual(torch.mean((sliceFromData_y-hyperslab_y)**2),0.0)
                self.assertEqual(torch.mean((sliceFromData_z-hyperslab_z)**2),0.0)
                
            def testSetHyperslabFloat(self):
                """setting hyperslab should change underlying volume (float)"""
                
                # read some data from somwhere
                v  = minc2_file(inputFile_ushort)
                dims=v.store_dims()
                v.setup_standard_order()
                hyperslab_a = v.load_hyperslab_t( [10, None, None] )
                
                v2 = minc2_file()
                v2.define(dims,'float32','float32')
                v2.create(outputFilename)
                v2.setup_standard_order()
                
                # because we are saving float32 , we don't need slice normalization
                v2.save_hyperslab_t(hyperslab_a,   [10,None,None] )
                hyperslab_b = v2.load_hyperslab_t( [10, None, None] )
                self.assertEqual(N.average((hyperslab_a-hyperslab_b)**2),0.0)
                v2.close()
                v.close()

            def testSetHyperslabShort(self):
                """setting hyperslab should change underlying volume (short)"""
                
                # read some data from somwhere
                v  = minc2_file(inputFile_ushort)
                dims=v.store_dims()
                v.setup_standard_order()
                hyperslab_a = v.load_hyperslab_t( [10, None, None] )
                
                # try with normalization
                v2 = minc2_file()
                v2.define(dims,'uint16','float32') # , global_scaling=True
                v2.create(outputFilename)
                v2.set_volume_range(torch.min(hyperslab_a),torch.max(hyperslab_a))
                v2.setup_standard_order()
                
                # have to set slice normalization
                v2.save_hyperslab_t(hyperslab_a,   [10,None,None] )
                hyperslab_b = v2.load_hyperslab_t( [10, None, None] )

                # compare results
                self.assertAlmostEqual(torch.mean((hyperslab_a-hyperslab_b)**2).item(),0.0,8)
                v2.close()
                v.close()
                
                
            def testHyperslabArray(self):
                """hyperslab should be reinsertable into volume"""
                if False:
                    v = minc2_file(inputFile_ushort)
                    v2 = minc2_file()
                    v2.create(outputFilename)
                    v2.close()
                    v.close()
    except ImportError:
        # torch is not available
        pass


class testVectorFiles(unittest.TestCase):
    """test reading and writing of vector files"""
    def testVectorRead(self):
        """make sure that a vector file can be read correctly"""
        v = minc2_file(inputVector)
        v.setup_standard_order()
        dims = v.representation_dims()
        self.assertEqual(dims[0].id, minc2_file.MINC2_DIM_VEC)
        v.close()
    def testVectorRead2(self):
        """make sure that volume has four dimensions"""
        v = minc2_file(inputVector)
        ndims = v.ndim()
        self.assertEqual(ndims, 4)
        data=v.data
        self.assertEqual(len(data.shape),4)
        v.close()
        
class testDirectionCosines(unittest.TestCase):
    """test that minc2_simple deals correctly with direction cosines"""
    def testDefaultDirCos3DVFF(self):
        """testing reading the direction cosines of a file with standard values (volumeFromFile)"""
        v = minc2_file(inputFile_ushort)
        #
        # This file was created without explicitly setting the direction cosines.
        # in that case, the attribute is not set altogether, so we should test
        # for it using the known defaults, because libminc does extract the correct
        # default values
        #
        v.setup_standard_order()
        dims = v.representation_dims()
        self.assertAlmostEqual(dims[0].dir_cos[0], 1.0, 8)
        self.assertAlmostEqual(dims[0].dir_cos[1], 0.0, 8)
        self.assertAlmostEqual(dims[0].dir_cos[2], 0.0, 8)
        
        self.assertAlmostEqual(dims[1].dir_cos[0], 0.0, 8)
        self.assertAlmostEqual(dims[1].dir_cos[1], 1.0, 8)
        self.assertAlmostEqual(dims[1].dir_cos[2], 0.0, 8)
        
        self.assertAlmostEqual(dims[2].dir_cos[0], 0.0, 8)
        self.assertAlmostEqual(dims[2].dir_cos[1], 0.0, 8)
        self.assertAlmostEqual(dims[2].dir_cos[2], 1.0, 8)
        
    def testNonDefaultDirCos3DVFF(self):
        """testing reading the direction cosines of a file with non-standard values (volumeFromFile)"""
        v = minc2_file(input3DdirectionCosines)
        v.setup_standard_order()
        dims = v.representation_dims()
        
        pipe = os.popen("mincinfo -attvalue xspace:direction_cosines %s" % input3DdirectionCosines, "r")
        from_file = pipe.read().rstrip().split(" ")
        pipe.close()
        
        self.assertAlmostEqual(dims[0].dir_cos[0], float(from_file[0]), 8)
        self.assertAlmostEqual(dims[0].dir_cos[1], float(from_file[1]), 8)
        self.assertAlmostEqual(dims[0].dir_cos[2], float(from_file[2]), 8)
        
        pipe = os.popen("mincinfo -attvalue yspace:direction_cosines %s" % input3DdirectionCosines, "r")
        from_file = pipe.read().rstrip().split(" ")
        pipe.close()
        self.assertAlmostEqual(dims[1].dir_cos[0], float(from_file[0]), 8)
        self.assertAlmostEqual(dims[1].dir_cos[1], float(from_file[1]), 8)
        self.assertAlmostEqual(dims[1].dir_cos[2], float(from_file[2]), 8)
        
        pipe = os.popen("mincinfo -attvalue zspace:direction_cosines %s" % input3DdirectionCosines, "r")
        from_file = pipe.read().rstrip().split(" ")
        pipe.close()
        self.assertAlmostEqual(dims[2].dir_cos[0], float(from_file[0]), 8)
        self.assertAlmostEqual(dims[2].dir_cos[1], float(from_file[1]), 8)
        self.assertAlmostEqual(dims[2].dir_cos[2], float(from_file[2]), 8)


class testWorlToVoxel(unittest.TestCase):
    """test that minc2_simple deals correctly with converting between coordinate systems"""


    def testWorldToVoxel(self):
        """testing world_to_voxel conversion of a file with non-standard values"""
        v = minc2_file(input3DdirectionCosines)
        v.setup_standard_order()
        xyz=N.array([50.0,-80.0,-30.0])
        ijk=v.world_to_voxel(xyz)
        self.assertAlmostEqual(ijk[0], -6.362013409627703453 , 8)
        self.assertAlmostEqual(ijk[1], 6.6280285942264356436 , 8)
        self.assertAlmostEqual(ijk[2], 75.806692060998855709 , 8)


    def testVoxelToWorld(self):
        """testing voxel_to_world conversion of a file with non-standard values"""
        v = minc2_file(input3DdirectionCosines)
        v.setup_standard_order()


        ijk=N.array([10,20,30])
        xyz=v.voxel_to_world(ijk)
        self.assertAlmostEqual(xyz[0],  -0.32337910608109865507, 8)
        self.assertAlmostEqual(xyz[1],  -73.698635237250869068,  8)
        self.assertAlmostEqual(xyz[2],  -29.421450173791534155,  8)

    def testWorldToVoxelVec(self):
        """Compare against binary world to voxel"""
        v = minc2_file(input3DdirectionCosines)
        v.setup_standard_order()

        x,y,z=N.meshgrid( N.linspace(-10,10,3),N.linspace(0,20,3),N.linspace(-5,15,3) )
        xyz=N.column_stack( ( N.ravel(x), N.ravel(y), N.ravel(z)))

        ijk=v.world_to_voxel(xyz)

        for i,x in enumerate(xyz):
            pipe = os.popen("worldtovoxel {} {} {} {}".format(input3DdirectionCosines,x[0], x[1], x[2]), "r")
            from_file = [float(i) for i in pipe.read().rstrip().split(" ")]
            pipe.close()
            for k in range(3):
                self.assertAlmostEqual(from_file[k], ijk[i,k] , 8)

if __name__ == "__main__":
    unittest.main()
