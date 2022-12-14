#! /bin/bash


if [[ ! -e adni_model_3d_v2/to_icbm_0_NL.xfm ]];then

antsRegistration --float 1 --minc 1 \
    --collapse-output-transforms 1 \
    --dimensionality 3 \
    --initialize-transforms-per-stage 0 --write-composite-transform 1 \
    --interpolation Linear \
    --output "[adni_model_3d_v2/to_icbm_]" \
    --transform SyN[0.1,3.0,0.5] \
    --metric CC[adni_model_3d_v2/model_t1w.mnc,icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c.mnc,1,4] \
    --convergence [ 200x150x100x50, 1e-07, 20 ] \
    --smoothing-sigmas 4.0x2.0x1.0x0.0vox \
    --shrink-factors 8x4x2x1 --use-histogram-matching 1 \
    --winsorize-image-intensities [ 0.05, 0.95 ]  \
    --write-composite-transform 0 

fi

if [[ ! -e adni_model_3d_v2/atlas_CerebrA.mnc ]];then
itk_resample --transform adni_model_3d_v2/to_icbm_0_NL.xfm \
    --invert_transform \
    --like adni_model_3d_v2/model_t1w.mnc \
     icbm152_model_09c/mni_icbm152_CerebrA_tal_nlin_sym_09c_b.mnc \
     adni_model_3d_v2/atlas_CerebrA.mnc --baa --order 1 --clob \
     --labels
fi

# WARP the prior surfaces
for hemi in lt rt; do
    falcon_transform_surface \
            icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.ply \
            adni_model_3d_v2/to_icbm_0_inverse_NL.xfm   \
            adni_model_3d_v2/adni_model_3d_v2_init_ics_${hemi}.ply
done
