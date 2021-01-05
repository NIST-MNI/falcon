#! /usr/bin/env python

# standard library
import string
import os
import argparse
import pickle
import sys
from minc2_simple import minc2_file

# minc
#import minc

# numpy
import numpy as np

# scikit-learn
# from sklearn import svm
# from sklearn import neighbors
# from sklearn import ensemble
# from sklearn import linear_model
# from sklearn import tree
# from sklearn import naive_bayes
from sklearn import mixture

from sklearn.pipeline import Pipeline
from sklearn import preprocessing


def load_labels( infile ):
    #with minc2_file(infile) as m:
    m=minc2_file(infile)
    m.setup_standard_order()
    data=m.load_complete_volume(minc2_file.MINC2_INT)
    return data


def load_image( infile ):
    #with minc2_file(infile) as m:
    m=minc2_file(infile)
    m.setup_standard_order()
    data=m.load_complete_volume(minc2_file.MINC2_FLOAT)
    return data


def save_labels( outfile, reference, data, history=None ):
    # TODO: add history
    ref=minc2_file(reference)
    out=minc2_file()
    out.define(ref.store_dims(), minc2_file.MINC2_BYTE, minc2_file.MINC2_INT)
    out.create(outfile)
    out.copy_metadata(ref)
    out.setup_standard_order()
    out.save_complete_volume(data)

def parse_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='Run tissue classifier ')
    
    parser.add_argument('image',    
                        help="Run classifier on a set of given images",nargs='+')
    
    parser.add_argument('--output', 
                        help="Output image")
    
    parser.add_argument('--mask', 
                        help="Mask output results, set to 0 outside" )

    parser.add_argument('--trainmask', 
                        help="Mask output results, set to 0 outside" )

 
    parser.add_argument('-n',
                        type=int,
                        help="Number of classes",default=3)
    
    parser.add_argument('--debug', action="store_true",
                        dest="debug",
                        default=False,
                        help='Print debugging information' )
    
    parser.add_argument('--coord', action="store_true",
                        dest="coord",
                        default=False,
                        help='Use image coordinates as additional features' )
    
    parser.add_argument('--random', type=int,
                        dest="random",
                        help='Provide random state if needed' )
    
    options = parser.parse_args()
    
    return options

if __name__ == "__main__":
    #history=minc.format_history(sys.argv)
    history='TODO:history' #TODO
    options = parse_options()
    
   
    # load prior and input image
    if options.image is not None:
        if options.debug: print("Loading images...")
        
        images= [ load_image(i) for i in options.image ]

        if options.coord:
            # add features dependant on coordinates
            c=np.mgrid[0:images[0].shape[0] , 0:images[0].shape[1] , 0:images[0].shape[2]]
        
            # use with center at 0 and 1.0 at the edge, could have used preprocessing 
            images.append( ( c[0]-images[0].shape[0]/2.0)/ (images[0].shape[0]/2.0) )
            images.append( ( c[1]-images[0].shape[1]/2.0)/ (images[0].shape[1]/2.0) )
            images.append( ( c[2]-images[0].shape[2]/2.0)/ (images[0].shape[1]/2.0) )

        mask=None
        weights=None

        if options.mask is not None:
            mask=load_labels(options.mask)

        # if options.weights is not None:
        #     weights=load_image(options.weights)

        if options.debug: print("Done")
        
        clf=None
        
        if options.trainmask is not None:
            trainmask = load_labels(options.trainmask)
            training_X = np.column_stack( tuple( np.ravel( j[ trainmask>0] ) for j in images  ) )
        else:
            training_X = np.column_stack( tuple( np.ravel( j ) for j in images  ) )
    
        if options.debug: print("Fitting...")

        gmm = mixture.GaussianMixture(  n_components=options.n, covariance_type='full')
        gmm.fit( training_X )
        
        # sorting labels
        ren = np.argsort(gmm.means_[:,0])
        if options.debug: 
            print(gmm)
            print(gmm.means_[:,0])
            print(ren)

        ren += 1
        print(ren)
        
        if mask is not None:
            if options.debug: print("Using mask")
            out_cls=np.zeros_like(images[0], dtype=np.int32 )
            out_cls[mask>0] = np.take(ren, gmm.predict( np.column_stack( tuple( np.ravel( j[ mask>0 ] ) for j in images  ) ) ))
        else:
            out_cls = np.take(ren, gmm.predict( np.column_stack( tuple( np.ravel( j ) for j in images  ) ) ))
        
        save_labels(options.output, options.image[0], out_cls, history=history)
    else:
        print("Error in arguments")


# kate: indent-width 4; replace-tabs on; remove-trailing-space on; hl python; show-tabs on
