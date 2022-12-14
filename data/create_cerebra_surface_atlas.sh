#! /bin/bash

# <extract surfaces on icbm template> 

# create mid surface
../build/src/igl/falcon_igl_mesh_avg mni_icbm152_t1_tal_nlin_sym_09c_out6.1_ocs-0.9.9.ply \
                    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_ics-0.9.9.ply \
                    -o mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9.ply

# split left-right
../build/src/igl/falcon_igl_mesh_split mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9.ply \
    mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-%d.ply 


for s in 0 1;do 
# sample
../build/src/igl/falcon_igl_field_sampler mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-${s}.ply \
    /data/data01/models/icbm152_model_09c/mni_icbm152_CerebrA_tal_nlin_sym_09c.mnc  cerebra_atlas_${s}.csv.gz \
    --header cerebra --labels --clobber

done
# resample to atlas space

../build/src/igl/falcon_igl_field_resample \
        -i cerebra_atlas_0.csv.gz \
        mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-0.ply \
        icbm152_model_09c/mni_icbm152_ics_sm_lt.ply \
        -o icbm152_model_09c/mni_icbm152_ics_sm_lt_atlas_cerebra.csv.gz  \
        --majority_invexp --SO3 --knn 3 --clobber

../build/src/igl/falcon_igl_field_resample \
        -i cerebra_atlas_1.csv.gz \
        mni_icbm152_t1_tal_nlin_sym_09c_out6.1_mid-0.9.9-1.ply \
        icbm152_model_09c/mni_icbm152_ics_sm_rt.ply \
        -o icbm152_model_09c/mni_icbm152_ics_sm_rt_atlas_cerebra.csv.gz  \
        --majority_invexp --SO3 --knn 3 --clobber

