# based on https://github.com/Mouse-Imaging-Centre/pyminc/blob/master/test/generatorTests.py

import unittest
import numpy as N
import os
import subprocess
import tempfile

from minc2_simple import minc2_file,minc2_xfm,minc2_error

def setUpModule():
    global outputXfmFilename1,outputXfmFilename2,outputXfmFilename3
    
     # testing for applying transformations to coordinates:
    outputXfmFilename1 = tempfile.NamedTemporaryFile(prefix="test-xfm-1", suffix=".xfm").name
    outputXfmFilename2 = tempfile.NamedTemporaryFile(prefix="test-xfm-2", suffix=".xfm").name
    outputXfmFilename3 = tempfile.NamedTemporaryFile(prefix="test-xfm-3", suffix=".xfm").name

    subprocess.check_call(["param2xfm", "-center", '2.21', '-3.765', '4.09', "-translation", '1.23', '6.4', '-7.8', "-scales", '0.2', '4.3', '-3', outputXfmFilename1])
    subprocess.check_call(["param2xfm", "-center", '-23.98', '0.46', '9.5', "-translation", '0.0', '-46', '89.3', "-scales", '10', '7.33', '84', outputXfmFilename2])
    subprocess.check_call(["xfmconcat", outputXfmFilename1, outputXfmFilename2, outputXfmFilename3])


def tearDownModule():
   
    os.remove(outputXfmFilename1)
    os.remove(outputXfmFilename2)
    os.remove(outputXfmFilename3)
   
class testXfmsAppliedToCoordinates(unittest.TestCase):
    """test that xfm files can be used to transform x,y,z coordinates"""
    def testForwardTransformSingleXfm(self):
        """testing coordinates transformed using the forward transform and a single transformation"""
        _xfm=minc2_xfm(outputXfmFilename1)
        out=_xfm.transform_point([6.68, 3.14, 7.00])
        self.assertAlmostEqual(out[0], 4.33400016486645, 8)
        self.assertAlmostEqual(out[1], 32.3265016365052, 8)
        self.assertAlmostEqual(out[2], -12.4399995803833, 8)
    
    def testInverseTransformSingleXfm(self):
        """testing coordinates transformed using the inverse transform and a single transformation"""
        
        _xfm=minc2_xfm(outputXfmFilename1)
        out=_xfm.inverse_transform_point([6.68, 3.14, 7.00])
        self.assertAlmostEqual(out[0], 18.4099990008772, 8)
        self.assertAlmostEqual(out[1], -3.64755821904214, 8)
        self.assertAlmostEqual(out[2], 0.520000139872233, 8)
    
    def testForwardTransformConcatenatedXfm(self):
        """testing coordinates transformed using the forward transform and a concatenated transformation"""
        
        _xfm=minc2_xfm(outputXfmFilename3)
        out=_xfm.transform_point([6.68, 3.14, 7.00])
        self.assertAlmostEqual(out[0],  259.159993714094, 8)
        self.assertAlmostEqual(out[1],  188.041454144745, 8)
        self.assertAlmostEqual(out[2], -1744.15997695923, 8)
    
    def testInverseTransformConcatenatedXfm(self):
        """testing coordinates transformed using the inverse transform and a concatenated transformation"""
        
        _xfm=minc2_xfm(outputXfmFilename3)
        out=_xfm.inverse_transform_point([6.68, 3.14, 7.00])
        self.assertAlmostEqual(out[0], -119.559994975925, 8)
        self.assertAlmostEqual(out[1], -2.72634880128239, 8)
        self.assertAlmostEqual(out[2], 0.0509524723840147, 8)
    
        
        
if __name__ == "__main__":
    unittest.main()
