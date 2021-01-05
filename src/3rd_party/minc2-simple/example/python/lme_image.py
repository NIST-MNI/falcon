#! /usr/bin/env python

import sys
import os
import csv

from minc2_simple import minc2_input_iterator,minc2_output_iterator,minc2_file

import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr
from rpy2.rinterface import RRuntimeError
import numpy as np

def load_csv(csv_file):
    '''Load csv file into a dictionary'''
    data={}
    # load CSV file 
    with open(input_csv,'r') as f:
        for r in csv.DictReader(f):
            for k in r.keys():
                try:
                    data[k].append(r[k])
                except KeyError:
                    data[k]=[r[k]]
    return data


# setup automatic conversion for numpy to Rpy
#numpy2ri.activate()

# import R objects
stats = importr('stats')
base  = importr('base')
#nlme  = importr('nlme')
lme   = importr('lmerTest')

# read the input data
input_csv='lng_t1nm.csv'
global_mask='mask.mnc'

# load CSV file
data=load_csv(input_csv)
# columns:
# signal,mask,subject,group,visit
# 

# setup R objects for performing linear modelling
subject = ro.FactorVector(data['subject'])
visit   = ro.FactorVector(data['visit'])
group   = ro.FactorVector(data['group'])

def run_lme(signal, mask):
    effects = ro.Formula('signal ~ group + visit +  (1|subject)')

    good_voxels=np.sum(mask>0.5)
    effects.environment["mask"]   = rm = ro.BoolVector(mask>0.5)
    effects.environment["signal"] = ro.FloatVector(signal).rx(rm)
    # assign variables 
    effects.environment["subject"]  = subject.rx(rm)
    effects.environment["visit"]    = visit.rx(rm)
    effects.environment["group"]    = group.rx(rm)
    
    # allocate space for output
    result=np.zeros(8)
    result[0]=good_voxels

    if good_voxels>4:
        try:
            # run linear mixed-effect model
            m = base.summary(lme.lmer(effects))
             # extract DF (for the visit)
            result[1]  = m.rx2('coefficients').rx(True,'df')[2]
            # extract coeffecients
            result[2:5] = m.rx2('coefficients').rx(True,'Estimate')[:]
            # extract t-values
            result[5:8] = m.rx2('coefficients').rx(True,'t value')[:]

        except RRuntimeError:
            # probably model didn't converge
            pass
    else:
        # not enough information
        pass

    return result




inp_signal=minc2_input_iterator(files=data['signal'],data_type=minc2_file.MINC2_DOUBLE)
inp_mask=minc2_input_iterator(  files=data['mask'],data_type=minc2_file.MINC2_DOUBLE)
inp_global_mask=minc2_input_iterator(  files=[global_mask],data_type=minc2_file.MINC2_DOUBLE)

# setup output iterator
out=minc2_output_iterator(files=["output_Count.mnc","output_df.mnc",
                                 "output_Intercept.mnc","output_group.mnc","output_visit.mnc",
                                 "output_Intercept_t.mnc","output_group_t.mnc","output_visit_t.mnc" ],
                            reference=data['signal'][0],data_type=minc2_file.MINC2_DOUBLE)
try:
    for signal,mask,gmask in zip(inp_signal,inp_mask,inp_global_mask):

        if gmask[0]>0.5:
            res=run_lme(signal,mask)
        else:
            res=np.zeros(8)

        out.next(res)

except StopIteration:
    print("Stopped early?")
    pass
# delete input iterator, free memory, close files, usually done automatically
inp_signal.close()
inp_mask.close()
# free up memory, close file not really needed here, usually done automatically
out.close()

# kate: space-indent on; indent-width 4; indent-mode python;replace-tabs on;word-wrap-column 80;show-tabs on;hl python
