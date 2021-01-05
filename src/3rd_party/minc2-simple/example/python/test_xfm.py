from minc2_simple import minc2_file
from minc2_simple import minc2_xfm
from minc2_simple import minc2_transform_parameters
import sys
import numpy as np

if __name__ == "__main__":
    # experiment with XFM files
    x=minc2_xfm()
    identity=np.eye(4)
    x.append_linear_transform(identity)
    x.save('identity.xfm')
    
    par=minc2_transform_parameters()
    par.translations=np.array([1,2,3])
    par.scales=np.array([1.1,1.0,0.9])
    
    y=minc2_xfm()
    y.append_linear_transform(par)
    y.save('test_par.xfm')
    z=minc2_xfm('test_par.xfm')
    par2=z.get_linear_transform_param(0)
    
    print(par)
    print(par.translations)
    print(par2)
    print(par2.translations)
    
    
    
