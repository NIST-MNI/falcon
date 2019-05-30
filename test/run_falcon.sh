#!/bin/bash

BINDIR=$1
MODELDIR=$2
INPUT=$3
OUTPUT=$4
WORKDIR=$5

if [[ -z $OUTPUT ]];then
 echo "Usage $0 BINDIR MODELDIR INPUT_PREFIX OUTPUT_PREFIX [WORKDIR]"
 exit 1
fi

VERBOSE=yes

####################################
# SCRIPT PARAMETERS
####################################
export OMP_NUM_THREADS=1 # multi-threading for parallel for niik
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS} # for itk

icbm_dir=$MODELDIR
# icbm_dir=${NIIKDIR}/data/icbm152_model_09c
stximg=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c.mnc
stxseg=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc

falcon_math=${BINDIR}/falcon_math

#setup file names and check that they are present

# STD linear
stx_t1w_mnc=${INPUT}.mnc
# BRAIN
stx_brain=${INPUT}_mask.mnc

# CLS
stx_cls=${INPUT}_cls.mnc
stx_lob=${INPUT}_lob.mnc
#stx_deep=${INPUT}_deep.mnc

# STD non-linear
nl_xfm=${INPUT}_nl.xfm

# VENTRICLES
stx_vent=${INPUT}_vent.mnc


# check input
function verify_files {
  while [ $# -gt 0 ]; do
    if [ ! -e $1 ];then 
      echo Missing $1
      exit 1
    fi
    shift
  done
}


function run_cmd {
  if [[ $VERBOSE == yes  ]];then
     echo Running: $@
     echo 
  fi
  $@
  if [[ $? != 0 ]]; then
    echo Failed: $@
    exit 1
  else
  echo 
  echo done
  echo 
  fi

  return 0
}

verify_files $stx_t1w_mnc $stx_brain $stx_cls $stx_lob $nl_xfm $stx_vent




img=$stx_t1w_mnc
ventmask=$stx_vent

if [ -z $WORKDIR ];then
    tempdir=$(mktemp -d -t FALCON.XXXXXX)
    trap "rm -rf $tempdir" 0 1 2 15
else
    tempdir=$WORKDIR
    mkdir -p $tempdir
fi

id=$(basename ${OUTPUT})
ffn=${tempdir}/${id}


# test
#run_cmd minccalc -labels -byte -express 'abs(A[0]-3)<0.5||abs(A[0]-9)<0.5?1:0' $stx_lob ${tempdir}/${id}_vent.mnc -clob
# experimental option
#ventmask=${tempdir}/${id}_vent.mnc

###############################################################
# 2. SEGMENTATION
###############################################################

wmmask=${ffn}_WM_mask.mnc
gwimask=${ffn}_GWI_mask.mnc


#if false;then

run_cmd itk_resample --labels --byte ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_deep2_4_cerebellum_dilate.mnc \
  ${ffn}_cerebellum_mask_tmp1.mnc --like ${img} --transform ${nl_xfm} --invert_transform  --clobber
# resample 
run_cmd itk_resample --like ${img} ${ffn}_cerebellum_mask_tmp1.mnc ${ffn}_cerebellum_mask_R.mnc --labels --byte --clobber
run_cmd mv -fv ${ffn}_cerebellum_mask_R.mnc ${ffn}_cerebellum_mask.mnc
run_cmd rm -fv ${ffn}_cerebellum_mask_tmp?.mnc  cerebellummask=${ffn}_cerebellum_mask.mnc

run_cmd $falcon_math maskout -in=$stx_brain -mask=${ffn}_cerebellum_mask.mnc -out=${ffn}_cerebral_brain_mask.mnc -uint8

echo "white matter segmentation"
WMmask=${ffn}_WM_mask
run_cmd minccalc  $stx_cls -express 'abs(A[0]-3)<0.5?1:0'  ${WMmask}_tmp1.mnc -clobber
run_cmd $falcon_math maskimg -in=${WMmask}_tmp1.mnc -mask=${stx_brain} -out=${WMmask}_tmp2.mnc 
run_cmd itk_resample --byte --labels --like ${img} ${ffn}_WM_mask_tmp2.mnc ${ffn}_WM_mask_R.mnc --clobber
run_cmd mv -fv ${ffn}_WM_mask_R.mnc ${ffn}_WM_mask.mnc
run_cmd rm -fv ${WMmask}_tmp?.mnc 


echo "  FALCON-filled mask"
run_cmd itk_resample --byte --labels ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_fill_mask.mnc ${ffn}_FALCON_fill_mask_tmp1.mnc --like ${img} --transform ${nl_xfm} --invert --clobber
run_cmd itk_resample --byte --labels --like ${img} ${ffn}_FALCON_fill_mask_tmp1.mnc ${ffn}_FALCON_fill_mask_R.mnc --clobber 
run_cmd mv -fv ${ffn}_FALCON_fill_mask_R.mnc ${ffn}_FALCON_fill_mask.mnc
run_cmd rm -fv ${ffn}_FALCON_fill_mask_tmp?.mnc

echo "  FALCON-removal mask"
run_cmd itk_resample --byte --labels ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_remove_mask.mnc ${ffn}_FALCON_remove_mask_tmp1.mnc --like ${img} --transform ${nl_xfm} --invert --clobber
run_cmd itk_resample --byte --labels --like ${img} ${ffn}_FALCON_remove_mask_tmp1.mnc ${ffn}_FALCON_remove_mask_R.mnc --clobber
run_cmd mv -fv ${ffn}_FALCON_remove_mask_R.mnc ${ffn}_FALCON_remove_mask.mnc
run_cmd rm -fv ${ffn}_FALCON_remove_mask_tmp?.mnc

echo "  brainstem-'cut' mask"
run_cmd itk_resample --byte --labels ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_brainstem.mnc ${ffn}_brainstem_mask_tmp1.mnc --like ${img} --transform ${nl_xfm} --invert --clobber
run_cmd itk_resample --byte --labels --like ${img} ${ffn}_brainstem_mask_tmp1.mnc ${ffn}_brainstem_mask_R.mnc  --clobber
run_cmd mv -fv ${ffn}_brainstem_mask_R.mnc ${ffn}_brainstem_mask.mnc
run_cmd rm -fv ${ffn}_brainstem_mask_tmp?.mnc 


echo "white matter - gray matter interface segmentation"
run_cmd minccalc -express '(A[0]+A[1]-A[2]-A[3]-A[4]+A[5])>0.5' ${ffn}_FALCON_fill_mask.mnc ${wmmask} ${ffn}_cerebellum_mask.mnc ${ffn}_brainstem_mask.mnc ${ffn}_FALCON_remove_mask.mnc ${ventmask} ${ffn}_GWI_mask_tmp1.mnc  -clobber

echo "  seedfill"
C=`$falcon_math world2voxel-int -in=${img} -xyz=0,14,9 | tail -n1 | awk '{print $2 "," $3 "," $4}'`
run_cmd $falcon_math seedfill -in=${ffn}_GWI_mask_tmp1.mnc -out=${ffn}_GWI_mask_tmp2.mnc -ijk=${C} 
echo "  close holes"
run_cmd $falcon_math closeholes -in=${ffn}_GWI_mask_tmp2.mnc -out=${ffn}_GWI_mask_tmp3.mnc 

echo "  close holes"
run_cmd $falcon_math median -in=${ffn}_GWI_mask_tmp3.mnc -out=${ffn}_GWI_mask_tmp4.mnc -radius=1.02 

run_cmd itk_resample --byte --labels  --like ${img} ${ffn}_GWI_mask_tmp4.mnc ${ffn}_GWI_mask_R.mnc 
mv -fv ${ffn}_GWI_mask_R.mnc ${ffn}_GWI_mask.mnc
rm -fv ${ffn}_GWI_mask_tmp?.mnc 

echo "non-ctx mask"
run_cmd itk_resample --byte --labels ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_nonctx_mask_2mm.mnc ${ffn}_FALCON_nonctx_mask.mnc --like ${img} --invert --clobber

echo "splitting into left/right"
run_cmd mincresample -byte ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_midsagittal_plane_CC4.mnc ${ffn}_FALCON_split_tmp1_.mnc -like ${img} -transform $nl_xfm -invert -clobber
run_cmd minccalc -express 'A[0]>0.1?1:0' ${ffn}_FALCON_split_tmp1_.mnc ${ffn}_FALCON_split_tmp1.mnc -clob
  
  #n=0
for n in `seq 0 1`; do 
  if [ ! -e ${ffn}_GWI_mask_rt_tmp1_${n}.mnc ];then
    run_cmd minccalc -express '(A[0]-A[1])>0.5' ${gwimask} ${ffn}_FALCON_split_tmp1.mnc ${ffn}_GWI_mask_rt_tmp1_${n}.mnc  -clobber
  fi
  
  C=`$falcon_math world2voxel-int -in=${img}  -xyz=24,-17,24 | tail -n1 | awk '{print $2 "," $3 "," $4}'`
  run_cmd $falcon_math seedfill -in=${ffn}_GWI_mask_rt_tmp1_${n}.mnc -out=${ffn}_GWI_mask_rt_tmp2_${n}.mnc -ijk=${C} 
  run_cmd $falcon_math thresh   -in=${ffn}_GWI_mask_rt_tmp2_${n}.mnc -out=${ffn}_GWI_mask_rt_tmp3_${n}.mnc -thresh=0.5 
  run_cmd itk_resample --byte --labels  --like ${img} ${ffn}_GWI_mask_rt_tmp3_${n}.mnc ${ffn}_GWI_mask_rt_R_${n}.mnc --clobber
  run_cmd mv -fv ${ffn}_GWI_mask_rt_R_${n}.mnc ${ffn}_GWI_mask_rt_${n}.mnc

  C=`$falcon_math world2voxel-int -in=${img} -xyz=-24,-17,24 | tail -n1 | awk '{print $2 "," $3 "," $4}'`
  run_cmd $falcon_math seedfill -in=${ffn}_GWI_mask_rt_tmp1_${n}.mnc -out=${ffn}_GWI_mask_lt_tmp2_${n}.mnc -ijk=${C} 
  run_cmd $falcon_math thresh -in=${ffn}_GWI_mask_lt_tmp2_${n}.mnc -out=${ffn}_GWI_mask_lt_tmp3_${n}.mnc -thresh=0.5 
  run_cmd itk_resample --byte --labels  --like ${img} ${ffn}_GWI_mask_lt_tmp3_${n}.mnc ${ffn}_GWI_mask_lt_R_${n}.mnc --clobber
  run_cmd mv -fv ${ffn}_GWI_mask_lt_R_${n}.mnc ${ffn}_GWI_mask_lt_${n}.mnc
  
  #run_kmd rm -fv ${ffn}_GWI_mask_?t_tmp?.mnc 
  LRcount=`mincstats -q -sum ${gwimask}`
  Lcount=`mincstats -q -sum ${ffn}_GWI_mask_lt_${n}.mnc`
  # Lpct=`fcalc -dig 0 $Lcount / $LRcount x 100 | awk '{print $1}'`
  Lpct=`echo $Lcount $LRcount | awk '{printf "%d\n", $1/$2*100}'`
  echo "  check left size = $Lpct = $Lcount / $LRcount"
  
  if [ ${Lpct} -gt 25 -a ${Lpct} -lt 75 ]; then echo "  left size is OK"; break; fi
  #run_kmd rm -fv ${ffn}_GWI_mask_?t_tmp?.mnc ${ffn}_GWI_mask_?t.mnc ${ffn}_GWI_mask_?t_R.mnc
  run_cmd mincmorph -success D ${ffn}_FALCON_split_tmp1_${n}.mnc ${ffn}_FALCON_split_tmp1_dilated_${n}.mnc -clobber
  run_cmd mv -fv ${ffn}_FALCON_split_tmp1_dilated_${n}.mnc ${ffn}_FALCON_split_tmp1_${n}.mnc
  if [ $n -eq 5 ]; then 
    #run_kmd rm -f ${ffn}_FALCON_split_tmp1.mnc
    echo ${ffn}_FALCON_split_tmp1.mnc
    echo "can't cut in 2 hemispheres... Please manually disconnect 2 hemispheres..." 
    echo "  Display ${img} -label ${gwimask}"
    exit 1
  fi
done
run_cmd mv ${ffn}_GWI_mask_lt_${n}.mnc ${ffn}_GWI_mask_lt.mnc
run_cmd mv ${ffn}_GWI_mask_rt_${n}.mnc ${ffn}_GWI_mask_rt.mnc

leftmask=${ffn}_GWI_mask_lt.mnc
rightmask=${ffn}_GWI_mask_rt.mnc

run_cmd itk_resample --byte --labels --like ${img} ${gwimask} $tempdir/gwimask.mnc --clobber
run_cmd mv -fv $tempdir/gwimask.mnc ${gwimask}
run_cmd $falcon_math bounds-color -in=${img} -imglist=${ffn}_brainstem_mask.mnc,${ffn}_FALCON_remove_mask.mnc,${ffn}_FALCON_fill_mask.mnc,${wmmask},${ventmask},${ffn}_cerebellum_mask.mnc,${stx_brain},${leftmask},${rightmask} -out=${ffn}_check_seg.nii.gz -vlist=210,30,210,90,60,255,60,90,255,255,120,80,220,150,150,120,120,255,60,255,60,180,255,20,255,180,20 

###############################################################
# 5. INITIAL INNER SURFACE
###############################################################
#    preparation
###############################################################

echo "initial inner surface"

for hemi in lt rt; do 
  echo "initial inner surface  $hemi"
  if [ -e ${ffn}_GWI_mask_init_ics_${hemi}.ply ]; then continue; fi
  run_cmd $falcon_math closebrain -in=${ffn}_GWI_mask_${hemi}.mnc -out=${ffn}_GWI_mask_${hemi}_close.mnc -radius=10.5 
  run_cmd $falcon_math dilate -in=${ffn}_GWI_mask_${hemi}_close.mnc -out=${ffn}_GWI_mask_${hemi}_dilate.mnc -radius=2 
  
  run_cmd $falcon_math laplacemap -imglist=${ffn}_GWI_mask_${hemi}.mnc,${ffn}_GWI_mask_${hemi}_dilate.mnc -out=${ffn}_GWI_${hemi}_Laplace_map.mnc -xyz=0.5,0.5,0.5 
  
  # not actually used!
  # run_cmd $falcon_math shrinkwrap -mask=${ffn}_GWI_mask_${hemi}_dilate.mnc -out=${ffn}_GWI_mask_${hemi}_dilate_obj.off -iter=10 -radius=35 -val=10 

  #run_kmd xfm2niikmat ${ffn}_stx.xfm $img $stximg ${ffn}_stx.mat 
  #run_kmd falcon_math applyaffine-obj -obj=${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.off -out=${ffn}_init_ics_${hemi}_from_icbm.off
  # copy all components
#   run_cmd cp ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.off ${ffn}_init_ics_${hemi}_from_icbm.off
#   run_cmd cp ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.offe ${ffn}_init_ics_${hemi}_from_icbm.offe
#   run_cmd cp ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.offsp ${ffn}_init_ics_${hemi}_from_icbm.offsp

### USING ICBM template
#  run_cmd ${BINDIR}/falcon_transform_off ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.off $nl_xfm  \
#                                  ${ffn}_init_ics_${hemi}_from_icbm.off --invert_transform --clob
#  run_cmd ${BINDIR}/falcon_cortex_initics ${img} ${ffn}_GWI_mask_${hemi}.mnc ${ffn}_GWI_${hemi}_Laplace_map.mnc \
#                                          ${ffn}_init_ics_${hemi}_from_icbm.off ${ffn}_GWI_mask_init_ics_${hemi}.off
###

   run_cmd ${BINDIR}/falcon_cortex_shrinkwrap ${ffn}_GWI_mask_${hemi}_dilate.mnc ${ffn}_GWI_mask_${hemi}_dilate2.ply --iter 50  --val 10 --elen 2 --clob #-radius=35
   run_cmd ${BINDIR}/falcon_cortex_shrinkwrap ${ffn}_GWI_mask_${hemi}_dilate.mnc ${ffn}_GWI_mask_${hemi}_dilate.ply  --iter 10  --val 2 --elen 1 --clob --obj=${ffn}_GWI_mask_${hemi}_dilate2.ply


   run_cmd ${BINDIR}/falcon_cortex_initics ${img} ${ffn}_GWI_mask_${hemi}.mnc ${ffn}_GWI_${hemi}_Laplace_map.mnc \
                                           ${ffn}_GWI_mask_${hemi}_dilate.ply ${ffn}_GWI_mask_init_ics_${hemi}.ply
done

run_cmd $falcon_math combine-obj -objlist=${ffn}_GWI_mask_init_ics_lt.ply,${ffn}_GWI_mask_init_ics_rt.ply -out=${ffn}_GWI_mask_init_ics.ply

###############################################################
# pial surface
###############################################################
echo "initial pial surface"
echo "  initial pial surface"

run_cmd ${BINDIR}/falcon_cortex_initocs \
      ${img} \
      ${ffn}_cerebral_brain_mask.mnc \
      ${ventmask} \
      ${gwimask} \
      ${ffn}_GWI_mask_init_ics.ply ${ffn}_GWI_mask_init_ocs.ply \
      --max-thick 3


###############################################################
# deform both surfaces
###############################################################

echo "refine both surfaces"
run_cmd ${BINDIR}/falcon_cortex_refine \
                  ${img} ${ffn}_cerebral_brain_mask.mnc \
                  ${ventmask} ${gwimask} \
                  ${ffn}_cerebellum_mask.mnc \
                  ${ffn}_brainstem_mask.mnc \
                  ${ffn}_GWI_mask_init_ics.ply \
                  ${ffn}_GWI_mask_init_ocs.ply \
                  ${OUTPUT}_ics.ply ${OUTPUT}_ocs.ply \
                  -nonctx-mask ${ffn}_FALCON_nonctx_mask.mnc \
                  -rungekutta
#fi

# split up to left and right surfaces
run_cmd ${BINDIR}/off_split ${OUTPUT}_ics.ply ${OUTPUT}_ics
run_cmd ${BINDIR}/off_split ${OUTPUT}_ocs.ply ${OUTPUT}_ocs

# calculate thickness
for s in 0 1;do
# remap to ICBM template
if [ $s == 0 ];then hemi=lt;else hemi=rt;fi

run_cmd ${BINDIR}/falcon_cortex_calc_thickness ${OUTPUT}_ics_${s}.ply ${OUTPUT}_ocs_${s}.ply ${OUTPUT}_thickness_${hemi}.txt

run_cmd ${BINDIR}/falcon_resample_field_sph \
          ${OUTPUT}_thickness_${hemi}.txt \
          ${OUTPUT}_ocs_${s}.ply \
          ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.ply \
          ${OUTPUT}_thickness_icbm_${hemi}.txt 
        
# convert to ply for easier visualisation
run_cmd ${BINDIR}/falcon_off2ply -sph ${icbm_dir}/data/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.ply \
        ${OUTPUT}_thickness_icbm_${hemi}.txt  \
        ${OUTPUT}_thickness_icbm_${hemi}.ply

        
done

echo -e "to check results: ${txtrst}"
echo "  NVK ${img} -obj=${ffn}_ics.ply,${ffn}_ocs.ply"
