#! /bin/bash

tempdir=$(mktemp -d --tmpdir)
trap "rm -rf $tempdir" 0 1 2 15

param2xfm $tempdir/identity.xfm

itk_resample --order 1 --like mni_icbm152_t1_tal_nlin_sym_09c.mnc \
    miccai2012_challenge/prior_ventricles.mnc \
    $tempdir/prior_ventricles.mnc 

minccalc -byte -labels -express 'A[0]>0.5?1:0' \
    $tempdir/prior_ventricles.mnc $tempdir/ventricles.mnc

./falcon_run_v2_d.sh \
    mni_icbm152_t1_tal_nlin_sym_09c.mnc \
    -omp 4 -debug -trace \
    -nl  $tempdir/identity.xfm -use_icbm \
    -brain mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc \
    -vent $tempdir/ventricles.mnc \
    -priors miccai2012_challenge/tissue_prior_3.mnc \
            miccai2012_challenge/tissue_prior_2.mnc \
            miccai2012_challenge/tissue_prior_1.mnc \
    mni_icbm152_t1_tal_nlin_sym_09c_out6.1
