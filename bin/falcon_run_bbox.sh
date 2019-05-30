#! /bin/bash
# SCRIPT:    Script to quickly test falcon parameters using a cutout of a scan (bbox)
# AUTHOR:    Vladimir S. FONOV
# DATE:      May 22, 2019
#
# export FALCON_BIN=<build directory>
# export FALCON_DATA=<>


####################################
# SCRIPT INFO
####################################

MAJOR_VERSION=0
MINOR_VERSION=1
MICRO_VERSION=8

progname=FALCON
pp=$$


####################################
# SCRIPT PARAMETERS
####################################

export OMP_NUM_THREADS=4 # multi-threading for parallel for niik
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS} # for itk

if [[ -z "$FALCON_BIN"  ]]; then
    echo "FALCON_BIN variable not set!"
    exit
fi
if [[ -z "$FALCON_DATA"  ]]; then
    echo "FALCON_DATA variable not set!"
    exit
fi
# if [[ -z "$BEASTDIR"  ]]; then
#     echo "BEASTDIR variable not set!"
#     exit
# fi
# if [[ -z "$VENTRICLEDIR"  ]]; then
#     echo "VENTRICLEDIR variable not set!"
#     exit
# fi

icbm_dir=$FALCON_DATA/data
# icbm_dir=${NIIKDIR}/data/icbm152_model_09c

variant=m81       # last variant

####################################
# FUNCTIONS
####################################

function Usage {
  cat <<EOF

  [${progname}] version ${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}

  [${progname}] usage: <input.mnc> <output_base>

  [${progname}] optional usage:
  -help                          :  show this usage
  -omp <omp>                     :  change number of processors [default=1]
  -vent <vent.mnc>               :  ventricle mask [default = None]
  -nl <nl.xfm>                   :  nonlinear registration to icbm [default = None]
  -gwimask <mask.mnc>            :  gray matter - white matter interface mask [default = None]
  -cerebellummask <mask.mnc>     :  cerebellum mask [default = None]
  -wmmask <mask.mnc>             :  white matter mask [default = None]
  -left <left.mnc>               :  left hemisphere mask [default = None]
  -right <right.mnc>             :  right hemisphere mask [default = None]
  -sides <left.mnc> <right.mnc>  :  left and right hemisphere masks [default = None]
  -variant  <var>                :  for debugging
  -priors <WM> <GM> <CSF>        :  tissue priors, ( default: none )
EOF
}

function emsg {
  echo $@
}


function Version {
  cat <<EOF

  [${progname}] version ${MAJOR_VERSION}.${MINOR_VERSION}.${MICRO_VERSION}

  0.0.0  November 17, 2013, Kunio Nakamura <knakamura@mrs.mni.mcgill.ca>
  -another test version

EOF
}

function run_cmd {
  echo Running: $@
  $@
  echo Done
  echo
  return 0
}


####################################
# SCRIPT STARTS
####################################

if [ $# -eq 0 ]; then Usage; exit 1; fi

while [ $# -gt 0 ]; do
  if   [ $1 = -help ]; then Usage; exit 1
  elif [ $1 = -u ]; then Usage; exit 1
  elif [ $1 = -omp ]; then export OMP_NUM_THREADS=$2; export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${OMP_NUM_THREADS}; shift 2
  elif [ $1 = -vent    ]; then ventmask=$2;shift 2;
  elif [ $1 = -cls     ]; then cls=$2;     shift 2;
  elif [ $1 = -nl      ]; then nlxfm=$2;   shift 2;
  elif [ $1 = -gwimask ]; then gwimask=$2; shift 2;
  elif [ $1 = -wmmask  ]; then wmmask=$2;  shift 2;
  elif [ $1 = -left    ]; then leftmask=$2;shift 2;
  elif [ $1 = -right   ]; then rightmask=$2;shift 2;
  elif [ $1 = -sides   ]; then leftmask=$2; rightmask=$3;shift 3;  
  elif [ $1 = -cerebellummask ]; then cerebellummask=$2;  shift 2; 
  elif [ $1 = -brain ]; then brainmask=$2;  shift 2; 
  elif [ $1 = -variant ]; then variant=$2;  shift 2; 
  elif [ $1 = -priors  ]; then prior_WM=$2; prior_GM=$3;prior_CSF=$4;shift 4;
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

nargs=2
if   [ ${#args[@]} -lt ${nargs} ]; then echo -e "${txtred}[${progname}] ERROR: too few arguments${txtrst}";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt ${nargs} ]; then echo -e "${txtred}[${progname}] ERROR: too many arguments${txtrst}"; echo "  args = ${args[@]}"; exit 9
fi


#niikmath=${FALCON_BIN}/bin/niikmath
niikmath=${FALCON_BIN}/falcon_math
midsag=${FALCON_BIN}/falcon_midsag
BINDIR=${FALCON_BIN}

alias date2='date +%Y-%m-%d\ %T'
shopt -s expand_aliases

img=${args[0]}
fn=${args[1]}
OUTPUT=$fn

gwimask=${fn}_GWI_mask.mnc
wmmask=${fn}_cls_3.mnc

set -e -x

echo "Denoising..."
if [ ! -e ${fn}_anlm.mnc ];then 
  run_cmd itk_minc_nonlocal_filter --anlm --regularize 0.5 $img ${fn}_anlm.mnc
fi


for hemi in lt rt; do 
  echo "initial inner surface  $hemi"
  
  if [ -e ${fn}_GWI_mask_init_ics_${hemi}.ply ]; then continue; fi
  
  # extract largest connected component first
  
  run_cmd mincmorph -successive 'GK[0.5:1.5:0]' ${fn}_GWI_mask_${hemi}.mnc ${fn}_GWI_mask_${hemi}_cc.mnc -clob
  
  run_cmd $niikmath closebrain -in=${fn}_GWI_mask_${hemi}_cc.mnc    -out=${fn}_GWI_mask_${hemi}_close.mnc  -radius=10.5
  run_cmd $niikmath dilate     -in=${fn}_GWI_mask_${hemi}_close.mnc -out=${fn}_GWI_mask_${hemi}_dilate.mnc -radius=2

  ##### DEBUG ##########
  # Use a mask to restrict further processing

  # create laplace field map at 0.5^3mm resolution
  run_cmd $niikmath laplacemap -imglist=${fn}_GWI_mask_${hemi}_cc.mnc,${fn}_GWI_mask_${hemi}_dilate.mnc -out=${fn}_GWI_${hemi}_Laplace_map.mnc -xyz=0.5,0.5,0.5 
  
  # not actually used!
  if [ -z $use_icbm ];then
    run_cmd ${BINDIR}/falcon_cortex_shrinkwrap ${fn}_GWI_mask_${hemi}_dilate.mnc  ${fn}_GWI_mask_${hemi}_dilate_obj.ply --iter 10 --radius 35 --val 10
    run_cmd ${BINDIR}/falcon_cortex_initics ${fn}_anlm.mnc \
                                        ${fn}_GWI_mask_${hemi}_cc.mnc ${fn}_GWI_${hemi}_Laplace_map.mnc \
                                        ${fn}_GWI_mask_${hemi}_dilate_obj.ply   ${fn}_GWI_mask_init_ics_${hemi}.ply
  else
    run_cmd ${BINDIR}/falcon_transform_off ${icbm_dir}/icbm152_model_09c/mni_icbm152_t1_tal_nlin_sym_09c_init_ics_${hemi}.ply ${nlxfm}  \
            ${fn}_init_ics_${hemi}_from_icbm.ply --invert_transform --clob
            
    run_cmd ${BINDIR}/falcon_cortex_initics ${fn}_anlm.mnc \
                                        ${fn}_GWI_mask_${hemi}_cc.mnc ${fn}_GWI_${hemi}_Laplace_map.mnc \
                                        ${fn}_init_ics_${hemi}_from_icbm.ply ${fn}_GWI_mask_init_ics_${hemi}.ply
  fi
done


if [ ! -e ${fn}_GWI_mask_init_ics.ply ];then
run_cmd $niikmath combine-obj -objlist=${fn}_GWI_mask_init_ics_lt.ply,${fn}_GWI_mask_init_ics_rt.ply -out=${fn}_GWI_mask_init_ics.ply
fi

###############################################################
# pial surface
###############################################################
echo "initial pial surface"
echo "  initial pial surface"

if [ ! -e ${fn}_GWI_mask_init_ocs.ply ];then
run_cmd ${BINDIR}/falcon_cortex_initocs ${fn}_anlm.mnc \
        ${fn}_cerebral_brain_mask.mnc ${ventmask} ${gwimask} \
        ${fn}_GWI_mask_init_ics.ply ${fn}_GWI_mask_init_ocs.ply
fi

###############################################################
# deform both surfaces
###############################################################

echo "Refine both surfaces"

# ${gwimask}

if [ ! -e ${OUTPUT}_ics_${variant}.ply ];then

mkdir -p ${OUTPUT}_trace_${variant}

export FALCON_TRACE=${OUTPUT}_trace_${variant}/trace
export FALCON_TRACE_X=40
export FALCON_TRACE_Y=35
export FALCON_TRACE_Z=70
export FALCON_TRACE_SCALE=2

# create csf and wm masks
# TODO: can repeat many times with different parameters

#export FALCON_TRACE_LOG=${OUTPUT}_trace_${variant}/trace.csv
#export FALCON_TRACE_LOG_ID=3920
if [ -n "$prior_WM"  ];then

  priors="-priorwm $prior_WM -priorgm $prior_GM -priorcsf $prior_CSF"

  run_cmd ${BINDIR}/falcon_cortex_refine \
        ${fn}_anlm.mnc ${fn}_cerebral_brain_mask.mnc \
        ${ventmask} ${wmmask}  \
        ${fn}_cerebellum_mask.mnc ${fn}_brainstem_mask.mnc \
        ${fn}_GWI_mask_init_ics.ply  ${fn}_GWI_mask_init_ocs.ply \
        ${OUTPUT}_ics_${variant}.ply ${OUTPUT}_ocs_${variant}.ply \
        -nonctx-mask ${fn}_nonctx_mask.mnc  \
          -gradient-FWHM    1.0 1.0 \
          -divergence-FWHM  1.0 1.0 \
          -prior-FWHM       1.0 1.0 \
          -wimag     0.7     0.7 \
          -wprior    0.2     0.3 \
          -wcurv     0.1     0.1 \
          -wsmooth   0.2     0.2 \
          -wssmooth  0.1     0.1 \
          -wtsmooth  0.5     0.5 \
          -wgrad     0.0     0.0 \
          -wflux     0.3     0.1 \
          -wprox     0.5     0.5 \
          -wbrain    1.0     1.0 \
          -wvent     1.0     1.0 \
          -wabs      0.2     0.2 \
          -tmin      0.0     0.0 \
          -tmax      5.0     5.0 \
          -tsigma    1.0     1.0 \
          -supdate   0.0     0.0 \
          -wmix      0.5     0.5 \
          -nomf             \
          -pmin      0.4    \
          -iter      100    \
          -delta     0.5    \
          -iter2     5      \
          -apply     0.2    \
          -remesh    0      \
          -rungekutta \
          -log ${OUTPUT}_convergence.csv \
          ${priors}

else

  run_cmd ${BINDIR}/falcon_cortex_refine \
        ${fn}_anlm.mnc ${fn}_cerebral_brain_mask.mnc \
        ${ventmask} ${wmmask}  \
        ${fn}_cerebellum_mask.mnc ${fn}_brainstem_mask.mnc \
        ${fn}_GWI_mask_init_ics.ply  ${fn}_GWI_mask_init_ocs.ply \
        ${OUTPUT}_ics_${variant}.ply ${OUTPUT}_ocs_${variant}.ply \
        -nonctx-mask ${fn}_nonctx_mask.mnc  \
          -gradient-FWHM     1.0 1.0 \
          -divergence-FWHM   1.0 1.0 \
          -wimag     0.8     0.7 \
          -wprior    0.2     0.3 \
          -wcurv     0.1     0.1 \
          -wsmooth   0.2     0.2 \
          -wssmooth  0.1     0.1 \
          -wtsmooth  0.5     0.5 \
          -wgrad     0.0     0.0 \
          -wflux     0.2     0.1 \
          -wprox     0.5     0.5 \
          -wbrain    1.0     1.0 \
          -wvent     1.0     1.0 \
          -wabs      0.2     0.2 \
          -tmin      0.0     0.0 \
          -tmax      5.0     5.0 \
          -tsigma    1.0     1.0 \
          -supdate   0.0     0.0 \
          -wmix      0.6     0.6 \
          -pmin      0.4    \
          -iter      100    \
          -delta     0.5    \
          -iter2     5      \
          -apply     0.2    \
          -remesh    0      \
          -log ${OUTPUT}_convergence.csv \
          -rungekutta 

fi


fi
