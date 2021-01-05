# based on https://github.com/Mouse-Imaging-Centre/pyminc/blob/master/test/generatorTests.py

import unittest
import numpy as np
import os
import subprocess
import tempfile

from minc2_simple import minc2_input_iterator,minc2_output_iterator,minc2_file

def setUpModule():
    global outputFilename
    global inputFile_float
    global outputAVG
    global dim_x,dim_y,dim_z

    dim_x=7
    dim_y=13
    dim_z=23

    outputFilename=[]
    inputFile_float=[]
    outputAVG = tempfile.NamedTemporaryFile(prefix="test-avg-", suffix=".mnc").name

    for i in range(4):
        f=tempfile.NamedTemporaryFile(prefix="test-", suffix=".mnc").name
        inputFile_float += [f]
        subprocess.check_call(['rawtominc', f, '-ofloat', '-input', '/dev/urandom', str(dim_x), str(dim_y), str(dim_z)])
        outputFilename += [tempfile.NamedTemporaryFile(prefix="test-out-", suffix=".mnc").name]
    
    # make average
    subprocess.check_call(['mincaverage','-float','-q']+inputFile_float+[outputAVG])
    

def tearDownModule():
    for f in inputFile_float:
        os.remove(f)
    os.remove(outputAVG)
    
    for f in outputFilename:
        if os.path.exists(f): os.remove(f)

class minc2_iterators(unittest.TestCase):
    """test the minc2_iterators """
    def testInput(self):
        
        """attempting to load something"""
        inp1=minc2_input_iterator(files=inputFile_float)
        inp2=minc2_input_iterator(files=outputAVG)
        cnt=0

        for i,j in zip(inp1,inp2):
            self.assertAlmostEqual(np.mean(i),j[0],6)
            cnt+=1
        # make sure we passed all voxels
        self.assertEqual(cnt,dim_x*dim_y*dim_z)

    def testOutput(self):
        """ensure byte data is read as float by default"""
        inp=minc2_input_iterator(files=outputAVG)
        out=minc2_output_iterator(files=outputFilename,reference=outputAVG,store_type=minc2_file.MINC2_FLOAT)
        cnt=0
        # write same value to all output files
        for i in inp:
            out.set_value([i[0]]*len(outputFilename))
            cnt+=1
            out.next()

        self.assertEqual(cnt,dim_x*dim_y*dim_z)
        # make sure all files are closed
        inp.close()
        out.close()

        # TODO: make sure we wrote correct info
        for f in outputFilename:
            pipe = os.popen("minccmp -q -rmse %s %s" % (f,outputAVG), "r")
            output = float(pipe.read())
            pipe.close()
            self.assertEqual(output, 0)

if __name__ == "__main__":
    unittest.main()
