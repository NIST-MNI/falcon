# based on https://github.com/Mouse-Imaging-Centre/pyminc/blob/master/test/generatorTests.py

import unittest
import numpy as N
import os
import subprocess
import tempfile

from minc2_simple import minc2_file, minc2_xfm, minc2_error, minc2_tags

def setUpModule():
    global tagFilename1,tagFilename2
    
     # testing for applying transformations to coordinates:
    tagFilename1 = tempfile.NamedTemporaryFile(prefix="test-tags1", suffix=".tag").name
    tagFilename2 = tempfile.NamedTemporaryFile(prefix="test-tags2", suffix=".tag").name
    
    with open(tagFilename1,'w') as f:
        f.write("""MNI Tag Point File
Volumes = 1;
% sample

Points =
 -39.1479525683508 -45.36582296518 -3.23530867558283 0 -1 -1 ""
 -39.1697845674545 -31.4945908372409 -12.6454579127804 0 -1 -1 ""
 -39.1528893087858 -12.736547262279 -19.5313607447827 0 -1 -1 ""
 44.4980182876921 -26.4856018308304 -18.2762593240494 0 -1 -1 ""
 44.5089190874764 -9.65617868177457 -26.413921582182 0 -1 -1 ""
 44.5352641099345 -45.9362159138925 -8.23509222318752 0 -1 -1 ""
 -0.750402788468333 -50.1286868142753 -73.5322832423614 0 -1 -1 ""
 -17.571498190414 0.402400648410242 33.4572970664624 0 -1 -1 ""
 -16.4704834167603 -12.6566404729527 33.5072391032967 0 -1 -1 ""
 -16.4715157056962 -26.2173239426784 33.5077521896347 0 -1 -1 ""
 19.4105560369935 2.06130787542515 33.5406287658606 0 -1 -1 ""
 18.3179199854826 -9.39715820391679 33.4929234623072 0 -1 -1 ""
 18.2770539201511 -22.9501031378546 33.4956185255666 0 -1 -1 ""
 18.2470549959536 35.4844009781165 -0.569014837963635 0 -1 -1 ""
 -16.1447884644039 35.392440668468 -0.708769405942937 0 -1 -1 ""
 -19.2135261139424 -88.530497477629 1.36036612200456 0 -1 -1 ""
 19.9738228151403 -82.2136875065283 0.761234963281745 0 -1 -1 "";
""")
    #
    
def tearDownModule():
    os.remove(tagFilename1)
    os.remove(tagFilename2)
   
class testTags(unittest.TestCase):
    """test that xfm files can be used to transform x,y,z coordinates"""
    def testTagsLoad(self):
        """Testing that we can open and read tags"""
        _tags=minc2_tags(tagFilename1)
        print("\nTag points:")
        print(_tags.tag[0])
        print("weights:")
        print(_tags.weights)
        print("structure_ids:")
        print(_tags.structure_ids)
        print("patient_ids")
        print(_tags.patient_ids)
        print("labels")
        print(_tags.labels)

    def testTagsSave(self):
        _tags=minc2_tags(tagFilename1)

        # save
        _tags.save(tagFilename2)
        import re

        with open(tagFilename1,'r') as f:
            tags1=f.read()

        with open(tagFilename2,'r') as f:
            tags2=f.read()

        # strip comment
        tags1=re.sub(r"^\%.*$","",tags1,flags=re.MULTILINE)
        tags2=re.sub(r"^\%.*$","",tags2,flags=re.MULTILINE)

        self.maxDiff=None
        self.assertEqual(tags1,tags2)
        


        
if __name__ == "__main__":
    unittest.main()
