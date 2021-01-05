#! /usr/bin/env python

import argparse
import subprocess
import traceback
import os
import numpy as np

# needed for matrix log and exp
import scipy.linalg

# needed to read and write XFM files
from minc2_simple import minc2_xfm


def do_cmd(cmds,verbose=False):
    try:
        if not verbose:
            with open(os.devnull, "w") as fnull:
                p=subprocess.Popen(cmds, stdout=fnull, stderr=subprocess.PIPE)
        else:
            p=subprocess.Popen(cmds, stderr=subprocess.PIPE)
        
        (output,output_stderr)=p.communicate()
        outvalue=p.wait()
    except OSError:
        raise Exception("ERROR: command {} Error:{}!\nMessage: {}\n{}".format(str(cmds),str(outvalue),output_stderr,traceback.format_exc()))
    if not outvalue == 0:
        raise Exception("ERROR: command {} failed {}!\nMessage: {}\n{}".format(str(cmds),str(outvalue),output_stderr,traceback.format_exc()))
    return outvalue

def xfmavg(inputs, output, verbose=False):
    # TODO: handle inversion flag correctly
    all_linear=True
    all_nonlinear=True
    input_xfms=[]
    input_grids=[]
    
    for j in inputs:
        x=minc2_xfm(j)
        if x.get_n_concat()==1 and x.get_n_type(0)==minc2_xfm.MINC2_XFM_LINEAR:
            # this is a linear matrix
            input_xfms.append(np.asmatrix(x.get_linear_transform()))
        else:
            all_linear&=False
            # strip identity matrixes
            nl=[]
            _identity=np.asmatrix(np.identity(4))
            _eps=1e-6
            if x.get_n_type(0)==minc2_xfm.MINC2_XFM_LINEAR and x.get_n_type(1)==minc2_xfm.MINC2_XFM_GRID_TRANSFORM:
                if scipy.linalg.norm(_identity-np.asmatrix(x.get_linear_transform(0)) )>_eps: # this is non-identity matrix
                    all_nonlinear&=False
                else:
                    input_grids.append(x.get_grid_transform(1)[0])
            elif x.get_n_type(1)==minc2_xfm.MINC2_XFM_GRID_TRANSFORM:
                input_grids.append(x.get_grid_transform(0)[0])
                
    if all_linear:
        acc=np.asmatrix(np.zeros([4,4],dtype=np.complex))
        for i in input_xfms:
            print(i)
            acc+=scipy.linalg.logm(i)
            
        acc/=len(input_xfms)
        acc=np.asarray(scipy.linalg.expm(acc).real,'float64','C')
        
        print("result:",acc)
        
        x=minc2_xfm()
        x.append_linear_transform(acc)
        x.save(output)
        
    elif all_nonlinear:
        
        output_grid=output.rsplit('.xfm',1)[0]+'_grid_0.mnc'
        
        cmds=['mincaverage','-clob']
        cmds.extend(input_grids)
        cmds.append(output_grid)
        do_cmd(cmds,verbose=verbose)
        x=minc2_xfm()
        x.append_grid_transform(output_grid, False)
        x.save(output)
    else:
        raise Exception("Mixed XFM files provided as input")

def parse_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Average MNI .xfm files')

    parser.add_argument("--verbose",
                    action="store_true",
                    dest="verbose",
                    default=False,
                    help="Print verbose information" )
    
    parser.add_argument("--clobber",
                    action="store_true",
                    dest="clobber",
                    default=False,
                    help="Overwrite output" )

    parser.add_argument("inputs",
                        nargs='+',
                        help="Input xfm files")

    parser.add_argument("output",
                        help="Output xfm file")
    
    options = parser.parse_args()

    return options    


if __name__ == '__main__':
    options = parse_options()
    if not options.clobber and os.path.exists(options.output):
        raise Exception("File {} exists! Run with --clobber to overwrite".format(options.output))
        
    xfmavg(options.inputs,options.output,verbose=options.verbose)


# kate: space-indent on; indent-width 4; indent-mode python;replace-tabs on;word-wrap-column 80
