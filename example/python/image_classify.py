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
from sklearn import svm
from sklearn import neighbors
from sklearn import ensemble
from sklearn import linear_model
from sklearn import tree
from sklearn import naive_bayes
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
    
    parser.add_argument('prior',
                        help="classification prior")
    
    parser.add_argument('image',    
                        help="Run classifier on a set of given images",nargs='+')
    
    parser.add_argument('--output', 
                        help="Output image")
    
    parser.add_argument('--mask', 
                        help="Mask output results, set to 0 outside" )

    parser.add_argument('--weights',
                        help="Training weights" )

    parser.add_argument('--trainmask', 
                        help="Training mask" )

    parser.add_argument('--method',
                        choices=['SVM','lSVM','nuSVM','NN','RanForest','AdaBoost','tree','Logistic','Bayes','GM'],
                        default='lSVM',
                        help='Classification algorithm')
    
    parser.add_argument('--preprocess',
                        choices=['StandardScaler','Normalizer','MinMax'],
                        help='Apply pre-processing stage')
    
    parser.add_argument('-n',
                        type=int,
                        help="integer parameter for classifier",default=15)
    
    parser.add_argument('-C',
                        type=float,
                        help="Regularizing parameter for classifier")
    
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
    
    parser.add_argument('--save',
                        help='Save training results in a file')
    parser.add_argument('--load',
                        help='Load training results from a file')
    
    options = parser.parse_args()
    
    return options

if __name__ == "__main__":
    #history=minc.format_history(sys.argv)
    history='TODO:history' #TODO
    options = parse_options()
    
    
    #print(repr(options))
    
    # load prior and input image
    if (options.prior is not None or options.load is not None) and options.image is not None:
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

        if options.weights is not None:
            weights=load_image(options.weights)

        if options.debug: print("Done")
        
        clf=None
        
        if options.load is not None:
            with open(options.load, 'rb') as f:
                clf = pickle.load(f)
        else:
            prior=load_labels(options.prior)

            labels=list(np.unique(prior))
            counts=list(np.bincount(np.ravel(prior)))

            if 0 in labels:
                if options.debug: print("Label 0 will be discarded...")
                labels.remove(0)
                counts.pop(0) # assume it's first
        
            if options.debug: print("Available labels:{} counts: {} available images:{}".format(repr(labels),repr(counts),len(images)))
        
            if options.debug: print("Creating training dataset for classifier")
        
            training_weights=None

            if options.trainmask is not None:
                trainmask = load_labels(options.trainmask)
            
                training_X = np.column_stack( tuple( np.ravel( j[ np.logical_and(prior>0 , trainmask>0 ) ] ) for j in images  ) )
                training_Y = np.ravel( prior[ np.logical_and(prior>0 , trainmask>0 ) ] )
                if weights is not None:
                    training_weights = np.ravel( weights[ np.logical_and(prior>0 , trainmask>0 ) ] )
            else:
                training_X = np.column_stack( tuple( np.ravel( j[ prior>0 ] ) for j in images  ) )
                training_Y = np.ravel( prior[ prior>0 ] )
                if weights is not None:
                    training_weights = np.ravel( weights[ prior>0 ] )

        
            if options.debug: print("Fitting...")
        
            if options.method=="SVM":
                clf = svm.SVC()
            elif options.method=="nuSVM":
                clf = svm.NuSVC()
            elif options.method=='NN':
                clf = neighbors.KNeighborsClassifier(options.n)
            elif options.method=='RanForest':
                clf = ensemble.RandomForestClassifier(n_estimators=options.n,random_state=options.random)
            elif options.method=='AdaBoost':
                clf = ensemble.AdaBoostClassifier(n_estimators=options.n,random_state=options.random)
            elif options.method=='tree':
                clf = tree.DecisionTreeClassifier(random_state=options.random)
            elif options.method=='Logistic':
                clf = linear_model.LogisticRegression(C=options.C)
            elif options.method=='Bayes':
                clf = naive_bayes.GaussianNB()
            elif options.method=='GM':
                clf = mixture.BayesianGaussianMixture(n_components=len(labels))
            else:
                clf = svm.LinearSVC(C=options.C)
            
            if options.preprocess is not None:
                pp=None
                if options.preprocess=='StandardScaler':
                    pp = preprocessing.StandardScaler()
                elif options.preprocess=='Normalizer':
                    pp = preprocessing.Normalizer()
                else : # MinMax
                    pp = preprocessing.MinMaxScaler()

                clf=Pipeline([('preprocess',pp),
                              ('clf',clf)])

                if training_weights is not None:
                    clf.fit(training_X, training_Y, clf__sample_weight=training_weights)
                else:
                    clf.fit(training_X, training_Y)
            else:
                if training_weights is not None:
                    clf.fit(training_X, training_Y, sample_weight=training_weights)
                else:
                    clf.fit(training_X, training_Y)
        
        if options.debug: print(clf)
        
        if options.save is not None:
            with open(options.save,'wb') as f:
                pickle.dump(clf, f)
        
        if options.output is not None:
            if options.debug: print("Classifying...")
        
            out_cls=None
        
            if mask is not None:
                if options.debug: print("Using mask")
                out_cls=np.zeros_like(images[0], dtype=np.int32 )
                out_cls[mask>0]=clf.predict( np.column_stack( tuple( np.ravel( j[ mask>0 ] ) for j in images  ) ) )
            else:
                out_cls=clf.predict( np.column_stack( tuple( np.ravel( j ) for j in images  ) ) )
        
            if options.debug: print("Saving output...")
            
            #out=minc.Label(data=out_cls)
            #out.save(name=options.output, imitate=options.image[0],history=history)
            save_labels(options.output, options.image[0], out_cls, history=history)
    else:
        print "Error in arguments"


# kate: indent-width 4; replace-tabs on; remove-trailing-space on; hl python; show-tabs on
