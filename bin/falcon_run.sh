#! /bin/bash
# FILENAME:  falcon_run.sh
# AUTHOR:    Kunio Nakamura, Vladimir S. FONOV
# DATE:      May 01, 2019
#
# example of variables for Falcon
# export FALCON_HOME=<falcon home directory>
# export FALCON_DATA=$FALCON_HOME/share/falcon (optional)
# export FALCON_BIN=$FALCON_HOME/bin (optional)
# export FALCON_SCRIPTS=$FALCON_HOME/bin (optional)

####################################
# SCRIPT INFO
####################################

MAJOR_VERSION=0
MINOR_VERSION=9
MICRO_VERSION=1
ver=${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}

progname=$(basename $0)

pp=$$

set -o pipefail -E -e

####################################
# SCRIPT PARAMETERS
####################################

export OMP_NUM_THREADS=2 # multi-threading for parallel for niik
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS} # for itk

if [[ -z "$FALCON_HOME"  ]]; then
  FALCON_HOME="/opt/minc/1.9.17"
fi

if [[ -z "$FALCON_DATA"  ]]; then
  FALCON_DATA="${FALCON_HOME}/share/falcon"
fi

if [[ -z "$FALCON_BIN"  ]]; then
  FALCON_BIN="${FALCON_HOME}/bin"
fi

if [[ -z "$FALCON_SCRIPTS"  ]]; then
  FALCON_SCRIPTS="${FALCON_HOME}/bin"
fi


icbm_dir=$FALCON_DATA/data

####################################
# FUNCTIONS
####################################

function Usage {
  cat <<EOF

  ${progname} version ${ver} 

  ${progname} <input.mnc> <output_base> -brain <brain_mask.mnc>

  ${progname} parameters:
  -help                          :  show this usage

  ---- Required parameters ---
  -brain <brain mask.mnc>        :  brain mask

  --- Recomended parameters ---
  -nl <nl.xfm>                   :  nonlinear registration to icbm [default = None, will run ANTs]
  -omp <omp>                     :  change number of processors [default=1]
  -use_icbm                      :  use mesh from ICBM for initialization , default - use shrink-wrap
  -anlm                          :  apply anlm filter to input t1w scan (if not done before)
  -postprocess                   :  apply post-procesiing: resample mesh to atlas, calculate thickness

  --- Optional parameters  ---
  -vent <vent.mnc>               :  ventricle mask [default = None]
  -cerebellum <cerebellum.mnc>   :  cerebellum mask [default = None]
  -brainstem <brainstem.mnc>     :  brainstem mask [default = None]
  -cls <cls.mnc>                 :  tissue classification map
  -priors <WM> <GM> <CSF>        :  tissue priors, ( default: none )
  -nopriors                      :  don't use tissue priors, computed internally

  --- Optional parameters (don't touch, if you don't know what you are doing)   ---
  -gwimask <mask.mnc>            :  gray matter - white matter interface mask [default = None]
  -cerebellummask <mask.mnc>     :  cerebellum mask [default = None]
  -wmmask <mask.mnc>             :  white matter mask [default = None]
  -csfmask <csf.mnc>             :  CSF mask [default = None]
  -left <left.mnc>               :  left hemisphere mask [default = None]
  -right <right.mnc>             :  right hemisphere mask [default = None]
  -sides <left.mnc> <right.mnc>  :  left and right hemisphere masks [default = None]
  -variant  <var>                :  for debugging
  -trace                         :  produce traces of intermediate surfaces
  -hr                            :  produce high-resolution thickness measurement in ICBM space
  -smooth <fwhm>                 :  smooth thickness maps (in addition)
  -noremesh                      :  don't remesh initial mesh (for DEBUGGING only)
  -debug                         :  run in debug mode (keep temp, echo commands)
  -verbose                       :  echo commands
EOF
}

function emsg {
  echo $@ 
}


function Version {
  cat <<EOF

  [${progname}] version ${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}

EOF
}


####################################
# SCRIPT STARTS
####################################

if [[ $# -eq 0 ]]; then Usage; exit 1; fi

while  [[ $# -gt 0 ]]; do
  if   [[ $1 = -help ]]; then Usage; exit 1
  elif [[ $1 = -u ]]; then Usage; exit 1
  elif [[ $1 = -omp ]]; then export OMP_NUM_THREADS=$2; export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS}; shift 2
  elif [[ $1 = -vent    ]]; then vent=$2;shift 2;
  elif [[ $1 = -cerebellum ]]; then cerebellum=$2;shift 2;
  elif [[ $1 = -brainstem ]]; then brainstem=$2;shift 2;
  elif [[ $1 = -cls     ]]; then cls=$2;     shift 2;
  elif [[ $1 = -nl      ]]; then nlxfm=$2;   shift 2;
  elif [[ $1 = -gwimask ]]; then gwimask=$2; shift 2;
  elif [[ $1 = -wmmask  ]]; then wmmask=$2;  shift 2;
  elif [[ $1 = -csfmask ]]; then csfmask=$2;  shift 2;
  elif [[ $1 = -left    ]]; then leftmask=$2;shift 2;
  elif [[ $1 = -right   ]]; then rightmask=$2;shift 2;
  elif [[ $1 = -sides   ]]; then leftmask=$2; rightmask=$3;shift 3;
  elif [[ $1 = -priors  ]]; then prior_WM=$2; prior_GM=$3;prior_CSF=$4;shift 4;
  elif [[ $1 = -create_prior ]];then create_prior=yes; shift;
  elif [[ $1 = -cerebellummask ]]; then cerebellummask=$2;  shift 2;
  elif [[ $1 = -brain ]]; then brainmask=$2;  shift 2;
  elif [[ $1 = -variant ]]; then variant=$2;  shift 2;
  elif [[ $1 = -smooth ]]; then smooth=$2;  shift 2;
  elif [[ $1 = -use_icbm ]];then use_icbm=yes; shift;
  elif [[ $1 = -nopriors ]];then nopriors=yes; shift;
  elif [[ $1 = -trace ]];then trace=yes; shift;
  elif [[ $1 = -anlm ]];then use_anlm=yes; shift;
  elif [[ $1 = -postprocess ]];then postprocess=yes; shift;
  elif [[ $1 = -hr ]];then use_hr=yes; shift;
  elif [[ $1 = -noremesh ]];then noremesh="--noremesh"; shift;
  elif [[ $1 = -debug ]];then DEBUG=yes ; VERBOSE=yes; shift;
  elif [[ $1 = -verbose ]];then VERBOSE=yes; shift;
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

nargs=2
if   [[ ${#args[@]} -lt ${nargs} ]]; then echo -e "[${progname}] ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [[ ${#args[@]} -gt ${nargs} ]]; then echo -e "[${progname}] ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi

if [[ -z $brainmask ]];then
 echo "[${progname}] need -brain <brain> "
 exit 9
fi

alias date2='date +%Y-%m-%d\ %T'
shopt -s expand_aliases

img=${args[0]}
fn=${args[1]}
OUTPUT=$fn
nm=$(basename $img .mnc)


if [[ -n "$VERBOSE" ]];then
set -x
fi

if [[ -n "$DEBUG" ]];then
tempdir=${OUTPUT}/temp
mkdir -p $tempdir
#TODO: cleanup temp?
else
tempdir=$(mktemp -d --tmpdir)
trap "rm -rf $tempdir" 0 1 2 15
fi

icbm_template=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c.mnc
icbm_template_mask=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc

###############################################################
# 1B. NONLINEAR STANDARD SPACE REGISTRATION
###############################################################
if [[ -n "$nlxfm" ]]; then
  echo "  using nl xfm: $nlxfm"
else
  if [[ ! -e ${fn}_nlstx.xfm ]];then
    antsRegistration \
   --use-histogram-matching 1 \
   --minc -a --dimensionality 3 \
   --metric "CC[${img},${icbm_template},1,3,Regular,0.5]" \
   --masks "[${brainmask},${icbm_template_mask}]" \
   --convergence [150x50x50,1.e-9,20] \
   --shrink-factors   8x4x2 \
   --smoothing-sigmas 8x4x2 \
   --transform "SyN[.25,2,0.5]" \
   --output ${fn}_nlstx
  fi
  nlxfm=${fn}_nlstx.xfm
fi

if [[ ! -z "${use_anlm}" ]];then
  echo "Denoising..."
  if [[ ! -e ${fn}_anlm.mnc ]];then 
    itk_minc_nonlocal_filter --anlm --regularize 0.5 $img ${fn}_anlm.mnc
  fi
  scan=${fn}_anlm.mnc
else
  scan=$img
fi

if [[ -z "$cls" ]];then
  # use atropos to do tissue classification
  cls=${fn}_cls.mnc
  if [[ ! -e $cls ]];then
    for l in $(seq 1 3);do
      itk_resample  ${icbm_dir}/miccai2012_challenge/tissue_prior_$l.mnc ${tempdir}/prior_${l}.mnc --transform $nlxfm --invert_transform --order 1  --like $scan --short
    done

    if [[ -z  "$ANTSPATH" ]];then
      export ANTSPATH=${MINC_TOOLKIT}/bin
    fi

    # run version patched for better mnc support
    falcon_antsAtroposN4.sh -d 3 -a $scan -x $brainmask -c 3 -o ${tempdir}/atropos_ -s mnc  -p ${tempdir}/prior_%d.mnc -y 3 -m 1
    cp ${tempdir}/atropos_Segmentation.mnc ${fn}_cls.mnc

    for l in $(seq 1 3);do
      mincreshape -q -clob -short ${tempdir}/atropos_SegmentationPosteriors${l}.mnc ${fn}_p_${l}.mnc
    done
  fi

  if [[ -z "$nopriors" ]]; then
    if [[ -z "$prior_WM" ]];  then prior_WM=${fn}_p_3.mnc; fi
    if [[ -z "$prior_GM" ]];  then prior_GM=${fn}_p_2.mnc; fi
    if [[ -z "$prior_CSF" ]]; then prior_CSF=${fn}_p_1.mnc;fi
  fi
else
  itk_split_labels $cls "${tempdir}/${nm}_p_%d.mnc" --expit 1.0 --aa --normalize
fi

###############################################################
# 1B.  Tissue split
###############################################################
# split into GM WM CSF
itk_split_labels $cls "${tempdir}/${nm}_cls_%d.mnc" --byte
wmmask=${tempdir}/${nm}_cls_3.mnc
gmmask=${tempdir}/${nm}_cls_2.mnc
csfmask=${tempdir}/${nm}_cls_1.mnc

###############################################################
# 2. SEGMENTATION
###############################################################

for l in brainstem cerebellum deep_gm hc ventricles;do
  mincresample  ${icbm_dir}/miccai2012_challenge/prior_${l}.mnc \
        ${tempdir}/${nm}_prior_${l}.mnc \
        -like ${scan} \
        -transformation ${nlxfm} \
        -invert -q 
done

if [[ -n "${cerebellummask}" ]]; then 
  echo "  using cerebellum mask: $cerebellummask"
else
  cerebellummask=${fn}_cerebellum_mask.mnc
  if [[ ! -e $cerebellummask ]];then
  minccalc -q -byte -labels -express 'A[0]>=0.5' ${tempdir}/${nm}_prior_cerebellum.mnc $cerebellummask
  fi
fi

if [[ -n "${ventmask}" ]]; then 
  echo "  using vent mask: $ventmask"
else
  ventmask=${fn}_vent_mask.mnc
  if [ ! -e $ventmask ];then
  minccalc -q -byte -labels -express 'A[0]*A[1]>=0.5' ${tempdir}/${nm}_prior_ventricles.mnc ${fn}_p_1.mnc ${ventmask}
  fi
fi

echo "cerebral segmentation" # reove cerebellum and brainstem from the brain mask
if [[ ! -e ${fn}_cerebral_brain_mask.mnc ]]; then
  minccalc -q -byte -labels -express 'A[1]<0.1&&A[2]<0.1?A[0]:0' ${brainmask} ${tempdir}/${nm}_prior_cerebellum.mnc ${tempdir}/${nm}_prior_brainstem.mnc ${fn}_cerebral_brain_mask.mnc
fi

# merge deep GM,HC and ventricles
if [[ ! -e ${fn}_fill_mask.mnc ]]; then 
  #echo "  filled mask"
  #mincresample ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_fill_mask.mnc  ${fn}_fill_mask.mnc -like ${img} -nearest_neighbour -transformation ${nlxfm} -invert  -labels
  minccalc -q -byte -labels  -express '(A[0]+A[1]+A[2])>0.1?1:0' \
    ${tempdir}/${nm}_prior_ventricles.mnc ${tempdir}/${nm}_prior_deep_gm.mnc ${tempdir}/${nm}_prior_hc.mnc ${fn}_fill_mask.mnc 
fi

# HC  - removal mask
if [[ ! -e ${fn}_remove_mask.mnc ]]; then 
  #echo "  removal mask"
  #mincresample ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_remove_mask.mnc ${fn}_remove_mask.mnc -like ${img} -nearest_neighbour -transformation ${nlxfm} -invert -labels
  minccalc -q -byte -labels  -express 'A[0]>=0.5?1:0' ${tempdir}/${nm}_prior_hc.mnc ${fn}_remove_mask.mnc
fi

# brainstem-'cut' maskbrainstem
if [[ ! -e ${fn}_brainstem_mask.mnc ]]; then
  #mincresample ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_symbrainstem_09c_CLADA_brainstem.mnc  ${fn}_brainstem_mask.mnc -like ${img} -nearest_neighbour -transformation ${nlxfm} -invert -labels
  minccalc -q -byte -labels  -express 'A[0]>0.5?1:0' ${tempdir}/${nm}_prior_brainstem.mnc ${fn}_brainstem_mask.mnc
fi

echo "white matter - gray matter interface segmentation"

if [[ -n "${gwimask}" ]]; then
  echo "  using white-gray interface mask:  ${gwimask}"
elif [[ ! -e ${fn}_GWI_mask.mnc ]]; then
  itk_morph --exp 'D[2]' ${ventmask} ${tempdir}/${nm}_vent_d.mnc
  itk_morph --exp 'D[2]' ${fn}_brainstem_mask.mnc ${tempdir}/${nm}_brainstem_d.mnc

  minccalc -byte -labels -express '(A[0]+A[1]-A[2]-A[3]-A[4]+A[5])>0.5?1:0' \
    ${wmmask} \
    ${fn}_fill_mask.mnc  \
    ${fn}_cerebellum_mask.mnc \
    ${tempdir}/${nm}_brainstem_d.mnc \
    ${fn}_remove_mask.mnc \
    ${tempdir}/${nm}_vent_d.mnc \
    ${tempdir}/${nm}_GWI_mask_tmp1.mnc

  #  seedfill
  ${FALCON_BIN}/falcon_math seedfill -in=${tempdir}/${nm}_GWI_mask_tmp1.mnc -out=${tempdir}/${nm}_GWI_mask_tmp2.mnc -xyz=0,14,9
  #  close holes
  ${FALCON_BIN}/falcon_math closeholes -in=${tempdir}/${nm}_GWI_mask_tmp2.mnc -out=${tempdir}/${nm}_GWI_mask_tmp3.mnc 
  # median
  ${FALCON_BIN}/falcon_math median -in=${tempdir}/${nm}_GWI_mask_tmp3.mnc -out=${tempdir}/${nm}_GWI_mask_tmp4.mnc -radius=1.02

  # take largest connected component
  #mincmorph -successive 'GB[1:1:1:0]' ${fn}_GWI_mask_tmp3.mnc ${fn}_GWI_mask_tmp4.mnc

  # resample
  mincresample -like ${scan}  ${tempdir}/${nm}_GWI_mask_tmp4.mnc ${fn}_GWI_mask.mnc -nearest_neighbour -labels -byte
  rm -fv ${tempdir}/${nm}_GWI_mask_tmp?.mnc
  gwimask=${fn}_GWI_mask.mnc
else
  gwimask=${fn}_GWI_mask.mnc
fi

# "non-ctx mask" # TODO: replace with something else ?
if [[ ! -e ${fn}_nonctx_mask.mnc ]]; then
  mincresample ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_nonctx_mask_2mm.mnc \
            ${fn}_nonctx_mask.mnc \
            -like ${img} \
            -nearest_neighbour \
            -transformation ${nlxfm} \
            -invert \
            -labels -q
fi

#setup tracing
# TODO: make this a parameter?
if [[ -n "${trace}" ]];then
  mkdir -p ${OUTPUT}/${ver}
  export FALCON_TRACE=${OUTPUT}/${ver}/trace
  export FALCON_TRACE_X=43
  export FALCON_TRACE_Y=90
  export FALCON_TRACE_Z=128
  export FALCON_TRACE_SCALE=2.0
fi

# "splitting into left/right"
if [[  -n "${leftmask}"  ]] || [[ -n "${rightmask}" ]]; then
  if [[  -n "${leftmask}" ]]; then
    minccalc -byte -labels -express 'A[1]>0.5?A[0]:0' ${gwimask} ${leftmask} ${fn}_GWI_mask_lt.mnc
  fi
  if [[  -n "${rightmask}" ]]; then
    mincmask -byte -labels -express 'A[1]>0.5?A[0]:0' ${gwimask} ${rightmask} ${fn}_GWI_mask_rt.mnc
  fi
  if [[ ! -e ${fn}_GWI_mask_lt.mnc ]]; then
    minccalc -express 'A[0]>A[1]?1:0' ${gwimask} ${rightmask} ${fn}_GWI_mask_lt.mnc -labels -byte
  elif [[ ! -e ${fn}_GWI_mask_rt.mnc ]]; then
    minccalc -express 'A[0]>A[1]?1:0' ${gwimask} ${leftmask} ${fn}_GWI_mask_rt.mnc -labels -byte
  else
    echo "missing masks?"
    eixt 1
  fi
elif [[ ! -e ${fn}_GWI_mask_lt.mnc ]] || [[ ! -e ${fn}_GWI_mask_rt.mnc ]]; then
    ${FALCON_BIN}/falcon_midsag $gwimask --left ${tempdir}/${nm}_GWI_mask_lt.mnc --right ${tempdir}/${nm}_GWI_mask_rt.mnc
    minccalc -labels -byte 'A[0]>0.5?1:0' ${tempdir}/${nm}_GWI_mask_lt.mnc ${fn}_GWI_mask_lt.mnc
    minccalc -labels -byte 'A[0]>0.5?1:0' ${tempdir}/${nm}_GWI_mask_rt.mnc ${fn}_GWI_mask_rt.mnc
fi

leftmask=${fn}_GWI_mask_lt.mnc
rightmask=${fn}_GWI_mask_rt.mnc

# generate QC volume
#if [ ! -e ${fn}_check_seg.mnc ];then
#  ${FALCON_BIN}/falcon_math bounds-color \
#  -in=${img} \
#  -imglist=${fn}_brainstem_mask.mnc,${fn}_remove_mask.mnc,${fn}_fill_mask.mnc,${wmmask},${ventmask},${fn}_cerebellum_mask.mnc,${brainmask},${leftmask},${rightmask} \
#  -out=${fn}_check_seg.mnc \
#  -vlist=210,30,210,90,60,255,60,90,255,255,120,80,220,150,150,120,120,255,60,255,60,180,255,20,255,180,20
#fi


if [[ ! -e ${fn}_GWI_mask_init_ics.ply ]];then
  # "initial inner surface"
  for hemi in lt rt; do 
    # extract largest connected component first
    mincmorph -successive 'GK[0.5:1.5:0]' ${fn}_GWI_mask_${hemi}.mnc ${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc
    
    ${FALCON_BIN}/falcon_math closebrain -in=${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc    -out=${tempdir}/${nm}_GWI_mask_${hemi}_close.mnc  -radius=10.5
    ${FALCON_BIN}/falcon_math dilate     -in=${tempdir}/${nm}_GWI_mask_${hemi}_close.mnc -out=${tempdir}/${nm}_GWI_mask_${hemi}_dilate.mnc -radius=2

    ##### DEBUG ##########
    # Use a mask to restrict further processing
    # erode cc mask, to move surface in
    #falcon_math erode      -in=${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc -out=${tempdir}/${nm}_GWI_mask_${hemi}_cce.mnc -radius=1

    # create laplace field map at 0.5^3mm resolution
    ${FALCON_BIN}/falcon_math laplacemap -imglist=${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc,${tempdir}/${nm}_GWI_mask_${hemi}_dilate.mnc -out=${tempdir}/${nm}_GWI_${hemi}_Laplace_map.mnc -xyz=0.5,0.5,0.5 

    if [[ -n "${trace}" ]];then
      export FALCON_TRACE=${OUTPUT}/${ver}/trace_${hemi}
    fi

    if [[ "$hemi" == "lt" ]];then
      export FALCON_TRACE_X=43
    else
      export FALCON_TRACE_X=153
    fi

    if [[ -z "$use_icbm" ]];then
      ${FALCON_BIN}/falcon_cortex_shrinkwrap ${tempdir}/${nm}_GWI_mask_${hemi}_dilate.mnc  \
        ${tempdir}/${nm}_GWI_mask_${hemi}_dilate_obj.ply --iter 10 --radius 35 --val 10

      ${FALCON_BIN}/falcon_cortex_initics ${scan} \
        ${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc ${tempdir}/${nm}_GWI_${hemi}_Laplace_map.mnc \
        ${tempdir}/${nm}_GWI_mask_${hemi}_dilate_obj.ply ${tempdir}/${nm}_GWI_mask_init_ics_${hemi}.ply ${noremesh}
    else
      ${FALCON_BIN}/falcon_transform_surface \
        ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.ply ${nlxfm}  \
        ${tempdir}/${nm}_init_ics_${hemi}_from_icbm.ply --invert_transform --clob

      ${FALCON_BIN}/falcon_cortex_initics ${scan} \
        ${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc ${tempdir}/${nm}_GWI_${hemi}_Laplace_map.mnc \
        ${tempdir}/${nm}_init_ics_${hemi}_from_icbm.ply ${tempdir}/${nm}_GWI_mask_init_ics_${hemi}.ply ${noremesh}
    fi
  done

  ${FALCON_BIN}/falcon_math combine-obj -objlist=${tempdir}/${nm}_GWI_mask_init_ics_lt.ply,${tempdir}/${nm}_GWI_mask_init_ics_rt.ply -out=${tempdir}/${nm}_GWI_mask_init_ics_1.ply
  # fix potential errors
  ${FALCON_BIN}/falcon_surface_check ${tempdir}/${nm}_GWI_mask_init_ics_1.ply ${fn}_GWI_mask_init_ics.ply --fix  # ${tempdir}/${nm}_GWI_mask_init_ics_2.ply 
fi

if [[ -n "${trace}" ]];then
  export FALCON_TRACE=${OUTPUT}/${ver}/trace
  export FALCON_TRACE_X=43
fi

###############################################################
# pial surface
###############################################################
# "initial pial surface"

if [[ ! -e ${fn}_GWI_mask_init_ocs.ply ]];then
  ${FALCON_BIN}/falcon_cortex_initocs ${scan} \
        ${fn}_cerebral_brain_mask.mnc \
        ${ventmask} ${gwimask} \
        ${fn}_GWI_mask_init_ics.ply \
        ${fn}_GWI_mask_init_ocs.ply \
        --max-thick 2 \
        --border $csfmask
fi

# combine cerebellum and brainstem in one mask
if [[ ! -e ${fn}_cerebellum_and_brainstem.mnc ]]; then
      minccalc -labels -byte -express 'A[0]>0.5||A[1]>0.5?1:0' \
        ${cerebellummask} \
        ${fn}_brainstem_mask.mnc \
        ${fn}_cerebellum_and_brainstem.mnc
fi


###############################################################
# refine both surfaces
###############################################################
if [[ ! -e ${OUTPUT}_ics-${ver}.ply ]];then

  export FALCON_TRACE_X=43
  export FALCON_TRACE_Y=90
  export FALCON_TRACE_Z=128
  export FALCON_TRACE_SCALE=2.0

  # find slices with errors
  if [[ -n "${trace}" ]];then
    export FALCON_TRACE=${OUTPUT}/${ver}/trace
  fi

  if [[ -n "$prior_WM" ]];then
  priors="-priorwm $prior_WM -priorgm $prior_GM -priorcsf $prior_CSF"

  ${FALCON_BIN}/falcon_cortex_refine \
          ${scan} ${fn}_cerebral_brain_mask.mnc \
          ${ventmask} ${wmmask}  \
          ${fn}_cerebellum_and_brainstem.mnc \
          ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply \
          ${OUTPUT}_ics-${ver}.ply    ${OUTPUT}_ocs-${ver}.ply \
          -nonctx-mask ${fn}_nonctx_mask.mnc  \
          -gradient-FWHM    1.0 1.0 \
          -divergence-FWHM  1.0 1.0 \
          -prior-FWHM       1.0 1.0 \
          -wimag     0.8     0.8 \
          -wprior    0.2     0.2 \
          -wcurv     0.1     0.1 \
          -wsmooth   0.2     0.2 \
          -wssmooth  0.1     0.1 \
          -wtsmooth  1.0     1.0 \
          -wgrad     0.0     0.0 \
          -wflux     0.1     0.1 \
          -wprox     0.5     0.5 \
          -wbrain    1.0     1.0 \
          -wvent     1.0     1.0 \
          -wabs      0.2     0.2 \
          -tmin      0.0     0.0 \
          -tmax      5.0     5.0 \
          -tsigma    1.0     1.0 \
          -supdate   0.0     0.0 \
          -wmix      0.5     0.5 \
          -pmin      0.4    \
          -iter      100    \
          -delta     0.5    \
          -iter2     5      \
          -apply     0.2    \
          -remesh    5      \
          -log ${OUTPUT}_convergence-${ver}.csv \
          ${priors} \
          -rungekutta
  else
  ${FALCON_BIN}/falcon_cortex_refine \
          ${scan} \
          ${fn}_cerebral_brain_mask.mnc \
          ${ventmask} \
          ${wmmask}  \
          ${fn}_cerebellum_and_brainstem.mnc \
          ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply \
          ${OUTPUT}_ics-${ver}.ply    ${OUTPUT}_ocs-${ver}.ply \
          -nonctx-mask ${fn}_nonctx_mask.mnc  \
          -gradient-FWHM    1.0 1.0 \
          -divergence-FWHM  1.0 1.0 \
          -prior-FWHM       1.0 1.0 \
          -wimag     0.8     0.8 \
          -wprior    0.0     0.0 \
          -wcurv     0.1     0.1 \
          -wsmooth   0.2     0.2 \
          -wssmooth  0.1     0.1 \
          -wtsmooth  1.0     1.0 \
          -wgrad     0.0     0.0 \
          -wflux     0.1     0.1 \
          -wprox     0.5     0.5 \
          -wbrain    1.0     1.0 \
          -wvent     1.0     1.0 \
          -wabs      0.5     0.5 \
          -tmin      0.0     0.0 \
          -tmax      5.0     5.0 \
          -tsigma    1.0     1.0 \
          -supdate   0.0     0.0 \
          -wmix      0.5     0.5 \
          -pmin      0.5    \
          -iter      100    \
          -delta     0.5    \
          -iter2     5      \
          -apply     0.2    \
          -remesh    5      \
          -log ${OUTPUT}_convergence-${ver}.csv \
          -rungekutta 
  fi
fi

if [[ -z "$use_icbm" ]];then
  echo "WARNING: not using ICBM source mesh, surface registration will probably produce garbage"
fi

# split up to left and right surfaces
${FALCON_BIN}/falcon_surface_split ${OUTPUT}_ics-${ver}.ply ${tempdir}/${nm}_ics-${ver}
${FALCON_BIN}/falcon_surface_split ${OUTPUT}_ocs-${ver}.ply ${tempdir}/${nm}_ocs-${ver}

# process left and right separately
for s in 0 1;do
  if [[ $s == 0 ]];then hemi=lt;else hemi=rt;fi
  
  atlas=${icbm_dir}/icbm152_model_09c/mni_icbm152_ics_sm_${hemi}_atlas_dirty.csv.gz
  model=${icbm_dir}/icbm152_model_09c/mni_icbm152_ics_sm_${hemi}.ply
  atlas_hr=${icbm_dir}/icbm152_model_09c/mni_icbm152_ocs_${hemi}_atlas_dirty.txt

  if [[ ! -e ${tempdir}/${nm}_ics-${ver}_${s}.ply ]];then
    echo "Something went wrong, failed to split surfaces" 1>&2
    exit 1
  fi

  #thickness in native space + atlas
  if [[ ! -e ${OUTPUT}_thickness-${ver}_${hemi}.csv.gz ]];then
    # measure thickness
    ${FALCON_BIN}/falcon_cortex_calc_thickness \
      ${tempdir}/${nm}_ics-${ver}_${s}.ply \
      ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
      ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv

    # resample atlas back into subject's space
    ${FALCON_BIN}/falcon_resample_field_sph \
            $atlas \
            $model \
            ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
            ${tempdir}/${nm}_thickness-${ver}_${hemi}_atlas.csv  --clob --nearest

    paste -d ','  ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
                  ${tempdir}/${nm}_thickness-${ver}_${hemi}_atlas.csv \
                  | gzip -9 -c  > ${OUTPUT}_thickness-${ver}_${hemi}.csv.gz
  fi

  #thickness in ICBM space, updated
  if [[ ! -e ${OUTPUT}_thickness_icbm-${ver}_${hemi}.csv.gz ]];then
    ${FALCON_BIN}/falcon_resample_field_sph \
              ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
              ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
              ${model} \
              ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}.csv --clob

    # unpack atlas
    gunzip -c $atlas > ${tempdir}/atlas_icbm_${hemi}.csv

    paste -d ','  ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}.csv \
                  ${tempdir}/atlas_icbm_${hemi}.csv \
                  | gzip -9 -c  > ${OUTPUT}_thickness_icbm-${ver}_${hemi}.csv.gz
  fi

  if [[ -n "${smooth}" ]];then
    # smoothed thickness in native space
    if [[ ! -e ${OUTPUT}_thickness-${ver}_${hemi}_sm_${smooth}.csv.gz ]];then
      # measure thickness
      ${FALCON_BIN}/falcon_cortex_calc_thickness \
        ${tempdir}/${nm}_ics-${ver}_${s}.ply \
        ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
        ${OUTPUT}_thickness-${ver}_${hemi}_sm_${smooth}.csv.gz \
          --smooth ${smooth}
    fi

    #thickness in ICBM space
    if [[ ! -e ${OUTPUT}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv.gz ]];then
      ${FALCON_BIN}/falcon_resample_field_sph \
                ${OUTPUT}_thickness-${ver}_${hemi}_sm_${smooth}.csv.gz \
                ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
                ${model} \
                ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv --clob

      gunzip -c $atlas > ${tempdir}/atlas_icbm_${hemi}.csv
      paste -d ','  ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv \
                    ${tempdir}/atlas_icbm_${hemi}.csv \
                    | gzip -9 -c  > ${OUTPUT}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv.gz

    fi
  fi

  if [[ ! -z $use_hr ]] && [[ ! -e ${OUTPUT}_thickness_icbm_hr-${ver}_${hemi}.csv.gz ]];then
    # resample thickness into common (MNI-ICBM152) hires space
    ${FALCON_BIN}/falcon_resample_field_sph \
              ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
              ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
              ${icbm_dir}/icbm152_model_09c/mni_icbm152_ocs_${hemi}.ply \
              ${tempdir}/${nm}_thickness_icbm_hr-${ver}_${hemi}.csv --clob

    paste -d ','  ${tempdir}/${nm}_thickness_icbm_hr-${ver}_${hemi}.csv \
                  $atlas_hr | gzip -9  -c > ${OUTPUT}_thickness-${ver}_${hemi}.csv
  fi
done


# make a qc image
if [[ ! -e ${OUTPUT}_qc-${ver}.png ]];then
  ${FALCON_SCRIPTS}/falcon_slice_qc.sh \
    ${scan} ${OUTPUT}_ics-${ver}.ply ${OUTPUT}_ocs-${ver}.ply \
        ${OUTPUT}_qc-${ver}.png
fi


if [[ ! -e ${OUTPUT}_qc_ocs-${ver}.png ]];then
  ${FALCON_SCRIPTS}/falcon_off_qc_2.sh \
        ${tempdir}/${nm}_ocs-${ver}_0.ply  ${OUTPUT}_thickness-${ver}_lt.csv.gz \
        ${tempdir}/${nm}_ocs-${ver}_1.ply  ${OUTPUT}_thickness-${ver}_rt.csv.gz \
        ${OUTPUT}_qc_ocs-${ver}.png -min 0 -max 7
fi

if [[ -n "${smooth}" ]] && [[ ! -e ${OUTPUT}_qc_ocs-${ver}_sm_${smooth}.png ]];then
  ${FALCON_SCRIPTS}/falcon_off_qc_2.sh  \
        ${tempdir}/${nm}_ocs-${ver}_0.ply  ${OUTPUT}_thickness-${ver}_lt_sm_${smooth}.csv.gz \
        ${tempdir}/${nm}_ocs-${ver}_1.ply  ${OUTPUT}_thickness-${ver}_rt_sm_${smooth}.csv.gz \
        ${OUTPUT}_qc_ocs-${ver}_sm_${smooth}.png -min 0 -max 7
fi
