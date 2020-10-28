#! /bin/bash
# FILENAME:  falcon_run_v2.sh
# AUTHOR:    Kunio Nakamura, Vladimir S. FONOV
# DATE:      Dec 30, 2019
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
MICRO_VERSION=3
ver=${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}

progname=$(basename $0)

pp=$$

set -o pipefail -E -e 

####################################
# SCRIPT PARAMETERS
####################################

export OMP_NUM_THREADS=2 # multi-threading for parallel for niik
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS} # for itk

if [[ -z "$FALCON_HOME" ]];then
  FALCON_HOME="${MINC_TOOLKIT}"
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

if [[ -z $ANTSPATH ]];then 
  export ANTSPATH=$MINC_TOOLKIT/bin;
fi


icbm_dir=$FALCON_DATA/data

####################################
# FUNCTIONS
####################################

function Usage {
  cat <<EOF

  ${progname} version ${ver} 

  ${progname} <input.mnc> <output_base> 

  ${progname} parameters:
  -help                          :  show this usage

  --- Recomended parameters ---
  -nl <nl.xfm>                   :  nonlinear registration to icbm [default = None, will run ANTs]
  -omp <omp>                     :  change number of processors [default=1]
  -use_icbm                      :  use mesh from ICBM for initialization , default - use shrink-wrap
  -anlm                          :  apply anlm filter to input t1w scan (if not done before)
  -postprocess                   :  apply post-procesiing: resample mesh to atlas, calculate thickness
  -brain <brain mask.mnc>        :  brain mask
  

  --- Optional parameters  ---
  -vent <vent.mnc>               :  ventricle mask [default = None]
  -cerebellum <cerebellum.mnc>   :  cerebellum mask [default = None]
  -brainstem <brainstem.mnc>     :  brainstem mask [default = None]
  -nopriors                      :  don't use tissue priors, computed internally

  --- Optional parameters (don't touch, if you don't know what you are doing)   ---
  -cls <cls.mnc>                 :  tissue classification map <don't use>
  -priors <WM> <GM> <CSF>        :  tissue priors, ( default: none )
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
  -hc                            :  Don't exclude Hippocampus from cortical fitting procedure
  -cb                            :  fit surfaces on cerebellum (experimental)
  -atrophy                       :  for subjects with lots of atrophy (using AD template and priors)
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
  elif [[ $1 = -vent    ]]; then ventmask=$2;shift 2;
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
  elif [[ $1 = -cb ]];then PROCESS_CB=yes; shift;
  elif [[ $1 = -hc ]];then KEEP_HC=yes; shift;
  elif [[ $1 = -atrophy ]];then ATROPHY=yes; shift;
  elif [[ $1 = -ad ]];then ATROPHY=yes; shift;
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

nargs=2
if   [[ ${#args[@]} -lt ${nargs} ]]; then echo -e "[${progname}] ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [[ ${#args[@]} -gt ${nargs} ]]; then echo -e "[${progname}] ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
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



if [[ -n "$ATROPHY" ]];then
  prior_base=${icbm_dir}/miccai2012_challenge_ad
  prior=miccai2012_challenge_ad
  template=${icbm_dir}/adni_model_3d_v2/model_t1w.mnc
  template_mask=${icbm_dir}/adni_model_3d_v2/model_t1w_mask.mnc
else
  prior_base=${icbm_dir}/miccai2012_challenge
  prior=miccai2012_challenge
  template=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c.mnc
  template_mask=${icbm_dir}/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc
fi



function atropos {
  # 1. Ventricle CSF
  # 2. Cortical GM
  # 3. WM
  # 4. Subcortical GM
  # 5. Brainstem
  # 6. Cerebelum MW
  # 7. Cerebelum GM
  # 8. Additional CSF: extracerebral + 4th ventricle

  # output: supratentorial
  # 1. supratentorial CSF
  # 2. GM cortical
  # 3. WM
  # 4. GM subcortical

  # output: rest:
  # 5. Brainstem
  # 6. Cerebelum WM
  # 7. Cerebelum GM
  # 8. Additional CSF: extracerebral + 4th ventricle

  scan=$1
  mask=$2
  xfm=$3
  out=$4

  # optional, 0.25 chosen using cross-validation
  w1=${6:-0.25}
  w2=${7:-0.25}
  atropos_prior="Socrates[1]" # no real difference between Aristotle and Socrates
  atropos_model="Gaussian"

  tissue_prior=${prior_base}/tissue_prior8

  for l in $(seq 1 8);do
    itk_resample  ${tissue_prior}_$l.mnc ${tempdir}/prior_${l}.mnc --transform $xfm --invert_transform --order 1  --like $scan --float
  done

  minccalc -express '(A[0]+A[1]+A[2]+A[3])>0.05&&A[4]>0.5?1:0' -label -byte ${tempdir}/prior_{1,2,3,4}.mnc $mask ${tempdir}/supra_mask.mnc
  minccalc -express '(A[0]+A[1]+A[2]+A[3])<0.05&&A[4]>0.5?1:0' -label -byte ${tempdir}/prior_{1,2,3,4}.mnc $mask ${tempdir}/infra_mask.mnc

  # supratentorial tissue priors + all of CSF
  minccalc -express 'A[0]+A[1]' ${tempdir}/prior_1.mnc ${tempdir}/prior_8.mnc  ${tempdir}/supra_p_1.mnc # CSF
  mv ${tempdir}/prior_2.mnc     ${tempdir}/supra_p_2.mnc # cortical GM
  mv ${tempdir}/prior_3.mnc     ${tempdir}/supra_p_3.mnc # WM
  mv ${tempdir}/prior_4.mnc     ${tempdir}/supra_p_4.mnc # Subcortical GM

  # infratentorial tissue priors + CSF
  mv ${tempdir}/prior_5.mnc ${tempdir}/infra_p_1.mnc # brainstem
  mv ${tempdir}/prior_6.mnc ${tempdir}/infra_p_2.mnc # Cerebellum WM
  mv ${tempdir}/prior_7.mnc ${tempdir}/infra_p_3.mnc # Cerebellum GM
  mv ${tempdir}/prior_8.mnc ${tempdir}/infra_p_4.mnc # CSF

  # GM+WM+CSF+Subcortical Gray
  ${FALCON_BIN}/falcon_antsAtroposN4.sh -d 3 -a $scan -x ${tempdir}/supra_mask.mnc \
      -c 4 -o ${tempdir}/supra_ -s mnc  -p ${tempdir}/supra_p_%d.mnc \
      -y 3 \
      -m 2 -n 20  \
      -f 4 \
      -e "[200]" \
      -i "[50x50x50,0.0]" \
      -l "[0.1,1.0]" \
      -b "${tissue_prior}" \
      -q "${atropos_model}" \
      -w ${w1}

  # Brainstem + Cerebellum + Extracortical CSF
  ${FALCON_BIN}/falcon_antsAtroposN4.sh -d 3 -a $scan -x ${tempdir}/infra_mask.mnc \
      -c 4 -o ${tempdir}/infra_ -s mnc  -p ${tempdir}/infra_p_%d.mnc \
      -y 1 -y 2 -y 3 \
      -m 2 -n 20  \
      -f 4 \
      -e "[200]" \
      -i "[50x50x50,0.0]" \
      -b "${tissue_prior}" \
      -q "${atropos_model}" \
      -l "[0.1,1.0]"  \
      -w ${w2} # need to distinguish brainstem from cerebellum WM TODO: join ?

  cp ${tempdir}/supra_Segmentation.mnc ${out}_supra.mnc
  cp ${tempdir}/infra_Segmentation.mnc ${out}_infra.mnc

  # fuzzy classes
  for l in $(seq 1 4);do
    mincreshape -byte ${tempdir}/supra_SegmentationPosteriors${l}.mnc ${out}_supra_${l}.mnc
  done

  for l in $(seq 1 4);do
    mincreshape -byte ${tempdir}/infra_SegmentationPosteriors${l}.mnc ${out}_infra_${l}.mnc
  done

  # merging supra and infra -tentorial regions together 
  # all CSF is label 1 now
  minccalc -express 'A[0]>0.5?A[0]:(abs(A[1]-4)<0.5?1:A[1]>0?A[1]+4:0)' \
    -labels -byte \
      ${tempdir}/supra_Segmentation.mnc \
      ${tempdir}/infra_Segmentation.mnc \
      ${out}_seg.mnc
}


if [[ -z $brainmask ]];then
  ${FALCON_BIN}/falcon_antsBrainExtraction.sh -d 3 \
    -s mnc \
    -a ${img} \
    -e ${template} \
    -m ${template_mask}  \
    -o ${tempdir}/brain_
  cp ${tempdir}/brain_BrainExtractionMask.mnc  ${fn}_mask.mnc
  brainmask=${fn}_mask.mnc
fi


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
      --metric "CC[${img},${template},1,3,Regular,0.5]" \
      --masks    "[${brainmask},${template_mask}]" \
      --convergence [150x100x50,1.e-9,20] \
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
  cls=${fn}_atropos_seg.mnc
  if [[ ! -e $cls ]];then
    atropos  $scan $brainmask $nlxfm ${fn}_atropos
  fi

  prior_SUPRA_CSF=${fn}_atropos_supra_1.mnc
  prior_SUPRA_GM=${fn}_atropos_supra_2.mnc
  prior_SUPRA_WM=${fn}_atropos_supra_3.mnc
  prior_DEEP=${fn}_atropos_supra_4.mnc

  prior_BS=${fn}_atropos_infra_1.mnc
  prior_CB_WM=${fn}_atropos_infra_2.mnc
  prior_CB_GM=${fn}_atropos_infra_3.mnc
  prior_INFRA_CSF=${fn}_atropos_infra_4.mnc

  #CSF
  if [[ -z $prior_CSF ]];then
    prior_CSF=${fn}_atropos_CSF.mnc
    if [[ ! -e $prior_CSF ]];then
      minccalc -express 'clamp(A[0]+A[1],0,1)' $prior_SUPRA_CSF $prior_INFRA_CSF $prior_CSF
    fi
  fi
  #GM
  if [[ -z $prior_GM ]];then
    prior_GM=${fn}_atropos_GM.mnc
    if [[ ! -e $prior_GM ]];then
      minccalc -express 'clamp(A[0]+A[1],0,1)' $prior_SUPRA_GM $prior_CB_GM $prior_GM
    fi
  fi
  #WM
  if [[ -z $prior_WM ]];then
    prior_WM=${fn}_atropos_WM.mnc
    if [[ ! -e $prior_WM ]];then
      minccalc -express 'clamp(A[0]+A[1],0,1)' $prior_SUPRA_WM $prior_CB_WM $prior_WM
    fi
  fi
else
  echo "This version doesn't work with external cls"
  exit 1
fi

###############################################################
# 1B.  Tissue split
###############################################################
if [[ -z $csfmask ]];then
  csfmask=${tempdir}/${nm}_csfmask.mnc
  if [[ ! -e $csfmask ]];then
  minccalc -q -byte -labels -express 'A[0]>=0.5?1:0' $prior_CSF $csfmask
  fi
fi

if [[ -z $wmmask ]];then
  wmmask=${tempdir}/${nm}_wmmask.mnc
  if [[ ! -e $wmmask ]];then
  minccalc -q -byte -labels -express 'A[0]>=0.5?1:0' $prior_SUPRA_WM $wmmask
  fi
fi


###############################################################
# 2. SEGMENTATION
###############################################################
for l in hc ventricles;do # brainstem cerebellum deep_gm
  if [[ ! -e ${tempdir}/${nm}_prior_${l}.mnc ]];then
    mincresample  ${icbm_dir}/miccai2012_challenge/prior_${l}.mnc \
          ${tempdir}/${nm}_prior_${l}.mnc \
          -like ${scan} \
          -transformation ${nlxfm} \
          -invert -q 
  fi
done

if [[ -n "${cerebellummask}" ]]; then 
   echo "  using cerebellum mask: $cerebellummask"
else
   cerebellummask=${fn}_cerebellum_mask.mnc
   if [[ ! -e $cerebellummask ]];then
    minccalc -q -byte -labels -express 'A[0]+A[1]>=0.5?1:0' $prior_CB_WM $prior_CB_GM $cerebellummask
   fi
fi

 if [[ -n "${ventmask}" ]]; then 
   echo "  using vent mask: $ventmask"
 else
   ventmask=${fn}_vent_mask.mnc
   if [[ ! -e $ventmask ]]; then
    # have to be conservative with the threshold - better to exclude some tissue which is not CSF
    minccalc -q -byte -labels -express 'A[0]*A[1]>=0.1?1:0' ${tempdir}/${nm}_prior_ventricles.mnc $prior_CSF ${ventmask}
   fi
 fi

# echo "cerebral segmentation" # remove cerebellum and brainstem from the brain mask
if [[ ! -e ${fn}_cerebral_brain_mask.mnc ]]; then
   minccalc -q -byte -labels -express 'A[1]<0.1&&A[2]<0.1&&A[3]<0.1?A[0]:0' ${brainmask} $prior_BS $prior_CB_WM $prior_CB_GM ${fn}_cerebral_brain_mask.mnc
fi

# only needed for CB surface
if [[ ! -e ${fn}_brain_mask2.mnc ]]; then
   minccalc -q -byte -labels -express 'A[1]<0.1?A[0]:0' ${brainmask} $prior_BS ${fn}_brain_mask2.mnc
fi

if [[ ! -e ${fn}_brainstem_mask.mnc ]]; then
  minccalc -q -byte -labels  -express 'A[0]>=0.5?1:0' $prior_BS ${fn}_brainstem_mask.mnc
fi

ventmaskd=${ventmask}
# HACK for large atrophy, exclude area around ventricles
# TODO: figure out if maybe always use D[3]
if [[ -n "${ATROPHY}" ]];then
  ventmaskd=$tempdir/${nm}_vent_mask_d.mnc
  itk_morph --exp 'D[3]' ${ventmask} ${ventmaskd}
else
  ventmaskd=$tempdir/${nm}_vent_mask_d.mnc
  itk_morph --exp 'D[1]' ${ventmask} ${ventmaskd}
fi


# "non-ctx mask" # TODO: replace with something else ?
if [[ ! -e ${fn}_nonctx_mask-${ver}.mnc ]]; then
  mincresample ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_CLADA_nonctx_mask_2mm.mnc \
            $tempdir/${nm}_nonctx_mask.mnc \
            -like ${img} \
            -nearest_neighbour \
            -transformation ${nlxfm} \
            -invert \
            -labels -q

  # combine with ventricle mask, hippocampus masks and brainstem masks
  if [[ -z ${KEEP_HC} ]]; then
    minccalc -q -labels -byte -express 'A[0]>=0.5||A[1]>=0.5||A[2]>=0.5?1:0' \
      $tempdir/${nm}_nonctx_mask.mnc ${ventmaskd} $prior_BS \
      ${fn}_nonctx_mask-${ver}.mnc
  else
    minccalc -q -labels -byte -express 'A[0]>=0.5||A[1]>=0.5||A[2]>=0.5||A[3]>=0.5?1:0' \
      $tempdir/${nm}_nonctx_mask.mnc ${ventmaskd} ${tempdir}/${nm}_prior_hc.mnc $prior_BS \
      ${fn}_nonctx_mask-${ver}.mnc
  fi
fi

# echo "white matter - gray matter interface segmentation" , fills deep GM, Ventricles and optionally HC too 
if [[ -n "${gwimask}" ]]; then
   echo "  using white-gray interface mask:  ${gwimask}"
else
  gwimask=${fn}_GWI_mask.mnc
  if [[ ! -e ${gwimask} ]]; then
    # ${tempdir}/${nm}_prior_ventricles.mnc - 
    if [[ -z ${KEEP_HC} ]]; then
      minccalc -q -byte -labels  -express '(A[0]+A[1]+A[2]+A[3])>=0.5?1:0' $prior_SUPRA_WM ${ventmaskd} $prior_DEEP ${tempdir}/${nm}_prior_hc.mnc ${gwimask}
    else
      minccalc -q -byte -labels  -express '(A[0]+A[1]+A[2])>=0.5?1:0'      $prior_SUPRA_WM ${ventmaskd} $prior_DEEP ${gwimask}
    fi
  fi
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
    minccalc -byte -labels -express 'A[1]>=0.5?A[0]:0' ${gwimask} ${leftmask} ${fn}_GWI_mask_lt.mnc
  fi
  if [[  -n "${rightmask}" ]]; then
    mincmask -byte -labels -express 'A[1]>=0.5?A[0]:0' ${gwimask} ${rightmask} ${fn}_GWI_mask_rt.mnc
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
    minccalc -labels -byte -express 'A[0]>=0.5?1:0' ${tempdir}/${nm}_GWI_mask_lt.mnc ${fn}_GWI_mask_lt.mnc
    minccalc -labels -byte -express 'A[0]>=0.5?1:0' ${tempdir}/${nm}_GWI_mask_rt.mnc ${fn}_GWI_mask_rt.mnc
fi

leftmask=${fn}_GWI_mask_lt.mnc
rightmask=${fn}_GWI_mask_rt.mnc

# generate QC volume

if [[ ! -e ${fn}_GWI_mask_init_ics.ply ]];then
  # "initial inner surface"
  for hemi in lt rt; do
    if [[ ! -e ${tempdir}/${nm}_GWI_mask_init_ics_${hemi}.ply ]];then
      # extract largest connected component first
      mincmorph -successive 'GK[0.5:1.5:0]' ${fn}_GWI_mask_${hemi}.mnc ${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc
      
      ${FALCON_BIN}/falcon_math closebrain -in=${tempdir}/${nm}_GWI_mask_${hemi}_cc.mnc    -out=${tempdir}/${nm}_GWI_mask_${hemi}_close.mnc  -radius=10.5
      ${FALCON_BIN}/falcon_math dilate     -in=${tempdir}/${nm}_GWI_mask_${hemi}_close.mnc -out=${tempdir}/${nm}_GWI_mask_${hemi}_dilate.mnc -radius=2

      ##### DEBUG ##########
      # Use a mask to restrict further processing
      # erode cc mask, to move surface in

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
    fi
  done

  # cerebellum processing
  if [[ -n "$PROCESS_CB" ]]; then
    if [[ ! -e ${tempdir}/${nm}_GWI_mask_init_ics_CB.ply ]];then
      minccalc -labels -byte -express 'A[0]>=0.5?1:0' $prior_CB_WM ${fn}_GWI_mask_CB.mnc
      mincmorph -successive 'GK[0.5:1.5:0]' ${fn}_GWI_mask_CB.mnc ${tempdir}/${nm}_GWI_mask_CB_cc.mnc

      ${FALCON_BIN}/falcon_math closebrain -in=${tempdir}/${nm}_GWI_mask_CB_cc.mnc    \
        -out=${tempdir}/${nm}_GWI_mask_CB_close.mnc  -radius=10.5

      ${FALCON_BIN}/falcon_math dilate     -in=${tempdir}/${nm}_GWI_mask_CB_close.mnc \
        -out=${tempdir}/${nm}_GWI_mask_CB_dilate.mnc -radius=2
      # create laplace field map at 0.5^3mm resolution
      ${FALCON_BIN}/falcon_math laplacemap -imglist=${tempdir}/${nm}_GWI_mask_CB_cc.mnc,${tempdir}/${nm}_GWI_mask_CB_dilate.mnc -out=${tempdir}/${nm}_GWI_CB_Laplace_map.mnc -xyz=0.5,0.5,0.5 
      # now we don't have ICBM atlas yet
      # TODO: add ICBM template
      # if [[ -z "$use_icbm" ]];then

      ${FALCON_BIN}/falcon_cortex_shrinkwrap ${tempdir}/${nm}_GWI_mask_CB_dilate.mnc  \
        ${tempdir}/${nm}_GWI_mask_CB_dilate_obj.ply --iter 10 --radius 20 --val 10

      if [[ -n "${trace}" ]];then
        export FALCON_TRACE=${OUTPUT}/${ver}/trace_CB
        # TODO: choose  better values
        export FALCON_TRACE_X=43
        export FALCON_TRACE_Y=90
        export FALCON_TRACE_Z=128
      fi

      ${FALCON_BIN}/falcon_cortex_initics ${scan} \
        ${tempdir}/${nm}_GWI_mask_CB_cc.mnc ${tempdir}/${nm}_GWI_CB_Laplace_map.mnc \
        ${tempdir}/${nm}_GWI_mask_CB_dilate_obj.ply ${tempdir}/${nm}_GWI_mask_init_ics_CB.ply ${noremesh}
    fi
    # adding cerebellum whitematter to GWI mask 
    if [[ ! -e ${fn}_GWI_mask2.mnc ]]; then
      minccalc -labels -byte -express 'A[0]>0.5||A[1]>0.5?1:0' ${gwimask} ${tempdir}/${nm}_GWI_mask_CB_cc.mnc ${fn}_GWI_mask2.mnc
    fi
    gwimask=${fn}_GWI_mask2.mnc
  fi

  if [[ -n "$PROCESS_CB" ]]; then
    ${FALCON_BIN}/falcon_math combine-obj -objlist=${tempdir}/${nm}_GWI_mask_init_ics_lt.ply,${tempdir}/${nm}_GWI_mask_init_ics_rt.ply,${tempdir}/${nm}_GWI_mask_init_ics_CB.ply -out=${tempdir}/${nm}_GWI_mask_init_ics_1.ply
  else
    ${FALCON_BIN}/falcon_math combine-obj -objlist=${tempdir}/${nm}_GWI_mask_init_ics_lt.ply,${tempdir}/${nm}_GWI_mask_init_ics_rt.ply -out=${tempdir}/${nm}_GWI_mask_init_ics_1.ply
  fi

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
  if [[ -n "$PROCESS_CB" ]];then
    ${FALCON_BIN}/falcon_cortex_initocs ${scan} \
          ${fn}_brain_mask2.mnc \
          ${ventmask} ${gwimask} \
          ${fn}_GWI_mask_init_ics.ply \
          ${fn}_GWI_mask_init_ocs.ply \
          --max-thick 4.0 --delta 0.02  \
          --vlist 0.5,0.2,0.1 \
          --border $csfmask
  else
    ${FALCON_BIN}/falcon_cortex_initocs ${scan} \
          ${fn}_cerebral_brain_mask.mnc \
          ${ventmask} ${gwimask} \
          ${fn}_GWI_mask_init_ics.ply \
          ${fn}_GWI_mask_init_ocs.ply \
          --max-thick 4.0 --delta 0.02  \
          --vlist 0.5,0.2,0.1 \
          --border $csfmask
  fi
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

  # TODO: update cerebellum mask ?
  if [[ -n "$PROCESS_CB" ]];then
    # TODO: make it work without CB mask?
    # HACK:
    minccalc -q -labels -byte -express '0' ${fn}_brain_mask2.mnc ${tempdir}/${nm}_zero.mnc -clob
    # HACK2: combine priors 
    mincmath -max $prior_WM $prior_CB_WM ${fn}_prior_WM.mnc 
    mincmath -max $prior_GM $prior_CB_GM ${fn}_prior_GM.mnc 
    minccalc -express 'A[0]>0.5||A[1]>=0.5' $wmmask $prior_CB_WM ${fn}_WM_AND_CB.mnc

    priors="-priorwm ${fn}_prior_WM.mnc -priorgm ${fn}_prior_GM.mnc -priorcsf $prior_CSF"

    ${FALCON_BIN}/falcon_cortex_refine \
            ${scan} ${fn}_brain_mask2.mnc \
            ${ventmask} ${fn}_WM_AND_CB.mnc  \
            ${fn}_brainstem_mask.mnc \
            ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply \
            ${OUTPUT}_ics-${ver}.ply    ${OUTPUT}_ocs-${ver}.ply \
            -nonctx-mask      ${fn}_nonctx_mask-${ver}.mnc  \
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
            -iter      200    \
            -delta     0.5    \
            -iter2     5      \
            -apply     0.2    \
            -remesh    5      \
            -log ${OUTPUT}_convergence-${ver}.csv \
            ${priors} \
            -rungekutta
    else
      ${FALCON_BIN}/falcon_cortex_refine \
              ${scan} ${fn}_cerebral_brain_mask.mnc \
              ${ventmask} ${wmmask}  \
              ${fn}_cerebellum_and_brainstem.mnc \
              ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply \
              ${OUTPUT}_ics-${ver}.ply    ${OUTPUT}_ocs-${ver}.ply \
              -nonctx-mask      ${fn}_nonctx_mask-${ver}.mnc  \
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
              -iter      200    \
              -delta     0.5    \
              -iter2     5      \
              -apply     0.2    \
              -remesh    5      \
              -log ${OUTPUT}_convergence-${ver}.csv \
              ${priors} \
              -rungekutta
    fi
  else
    ${FALCON_BIN}/falcon_cortex_refine \
            ${scan} \
            ${fn}_cerebral_brain_mask.mnc \
            ${ventmask} \
            ${wmmask}  \
            ${fn}_cerebellum_and_brainstem.mnc \
            ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply \
            ${OUTPUT}_ics-${ver}.ply    ${OUTPUT}_ocs-${ver}.ply \
            -nonctx-mask ${fn}_nonctx_mask-${ver}.mnc  \
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
            -iter      200    \
            -delta     0.5    \
            -iter2     5      \
            -apply     0.2    \
            -remesh    5      \
            -log ${OUTPUT}_convergence-${ver}.csv \
            -rungekutta 
  fi
fi

# debugging 
if [[ `type -p RScript` > /dev/null ]]; then
    cat - >${OUTPUT}_convergence-${ver}.R <<END
library( ggplot2 )
conv <- read.csv( "${OUTPUT}_convergence-${ver}.csv" )
myPlot <- ggplot( conv, aes( x = iteration, y = mean_ics ) ) +
 geom_line( data = conv, aes( x = iteration, y = mean_ics ) , col="blue" ) +
 geom_line( data = conv, aes( x = iteration, y = mean_ocs ) , col="red" ) +
 theme( legend.position = "none" )
ggsave( filename = "${OUTPUT}_convergence-${ver}.png", plot = myPlot, width = 4, height = 3, units = 'in' )
END
`RScript  ${OUTPUT}_convergence-${ver}.R`
rm -f  ${OUTPUT}_convergence-${ver}.R
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
    ${FALCON_BIN}/falcon_igl_field_resample \
            -i $atlas \
            $model \
            ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
            -o ${tempdir}/${nm}_thickness-${ver}_${hemi}_atlas.csv  \
            --majority_invexp --SO3 --knn 3

    paste -d ','  ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
                  ${tempdir}/${nm}_thickness-${ver}_${hemi}_atlas.csv \
                  | gzip -9 -c  > ${OUTPUT}_thickness-${ver}_${hemi}.csv.gz
  fi

  #thickness in ICBM space, updated
  if [[ ! -e ${OUTPUT}_thickness_icbm-${ver}_${hemi}.csv.gz ]];then
    ${FALCON_BIN}/falcon_igl_field_resample \
              -i ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
              ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
              ${model} \
              -o ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}.csv \
              --knn 3 --SO3  --invexp

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
      ${FALCON_BIN}/falcon_igl_field_resample \
                -i ${OUTPUT}_thickness-${ver}_${hemi}_sm_${smooth}.csv.gz \
                ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
                ${model} \
                -o ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv \
                --knn 3 --SO3  --invexp

      gunzip -c $atlas > ${tempdir}/atlas_icbm_${hemi}.csv
      paste -d ','  ${tempdir}/${nm}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv \
                    ${tempdir}/atlas_icbm_${hemi}.csv \
                    | gzip -9 -c  > ${OUTPUT}_thickness_icbm-${ver}_${hemi}_sm_${smooth}.csv.gz
    fi
  fi

  if [[ ! -z $use_hr ]] && [[ ! -e ${OUTPUT}_thickness_icbm_hr-${ver}_${hemi}.csv.gz ]];then

    if [[ ! -e ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv ]];then
      ${FALCON_BIN}/falcon_cortex_calc_thickness \
        ${tempdir}/${nm}_ics-${ver}_${s}.ply \
        ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
        ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv
    fi


    # resample thickness into common (MNI-ICBM152) hires space
    ${FALCON_BIN}/falcon_igl_field_resample \
          -i ${tempdir}/${nm}_thickness-${ver}_${hemi}.csv \
          ${tempdir}/${nm}_ocs-${ver}_${s}.ply \
          ${icbm_dir}/icbm152_model_09c/mni_icbm152_ocs_${hemi}.ply \
          -o ${tempdir}/${nm}_thickness_icbm_hr-${ver}_${hemi}.csv \
          --knn 3 --SO3  --invexp

    paste -d ','  ${tempdir}/${nm}_thickness_icbm_hr-${ver}_${hemi}.csv \
                  $atlas_hr | gzip -9  -c > ${OUTPUT}_thickness_icbm_hr-${ver}_${hemi}.csv.gz
  fi
done


# make a qc image
if [[ ! -e ${OUTPUT}_qc-${ver}.png ]];then
  ${FALCON_SCRIPTS}/falcon_slice_qc.sh \
    ${scan} ${OUTPUT}_ics-${ver}.ply ${OUTPUT}_ocs-${ver}.ply \
        ${OUTPUT}_qc-${ver}.png
fi

if [[ ! -e ${OUTPUT}_qc_ocs_thickness-${ver}.png ]];then
  ${FALCON_SCRIPTS}/falcon_off_qc_2.sh \
        ${tempdir}/${nm}_ocs-${ver}_0.ply  ${OUTPUT}_thickness-${ver}_lt.csv.gz \
        ${tempdir}/${nm}_ocs-${ver}_1.ply  ${OUTPUT}_thickness-${ver}_rt.csv.gz \
        ${OUTPUT}_qc_ocs_thickness-${ver}.png \
        -min 0 -max 7
fi

if [[ -n "${smooth}" ]] && [[ ! -e ${OUTPUT}_qc_ocs_thickness-${ver}_sm_${smooth}.png ]];then
  ${FALCON_SCRIPTS}/falcon_off_qc_2.sh  \
        ${tempdir}/${nm}_ocs-${ver}_0.ply  ${OUTPUT}_thickness-${ver}_lt_sm_${smooth}.csv.gz \
        ${tempdir}/${nm}_ocs-${ver}_1.ply  ${OUTPUT}_thickness-${ver}_rt_sm_${smooth}.csv.gz \
        ${OUTPUT}_qc_ocs_thickness-${ver}_sm_${smooth}.png \
        -min 0 -max 7
fi

if [[ ! -e ${OUTPUT}_qc_ocs_lobe-${ver}.png ]];then
  ${FALCON_SCRIPTS}/falcon_off_qc_2.sh \
        ${tempdir}/${nm}_ocs-${ver}_0.ply  ${OUTPUT}_thickness-${ver}_lt.csv.gz \
        ${tempdir}/${nm}_ocs-${ver}_1.ply  ${OUTPUT}_thickness-${ver}_rt.csv.gz \
        ${OUTPUT}_qc_ocs_lobe-${ver}.png \
        -min 0 -max 255 \
        -column lobe -discrete -levels 256
fi
