#! /bin/bash
set -e
#set -v
#set -x

BINDIR=$1
SAMPLE_BRAIN=$2
OUT=$3
GM=$4
BCK=$5

# optional
THICK=$6
TMP=$7
THRESHOLD=$8
CPLX=$9

VERBOSE=yes

if [ -z $BCK ];then
  echo "Usage $0 <bindir> <input> <output_prefix> <GM threshold> <BRAIN Threshold> [expected average thickness] [work directory] [threshold] [trace_flag] [cplx]"
  exit 1
fi

if [ -z $TMP ];then
  tempdir=$(mktemp -d -t FALCON.XXXXXX)
  trap "rm -rf $tempdir" 0 1 2 15
else
  tempdir=$TMP 
  if [ ! -e $tempdir ];then
    mkdir -p $tempdir
  fi
  # add tracing

  export FALCON_TRACE=$TMP/trace

  if [ ! -z ${CPLX} ];then
  export FALCON_TRACE_X=140
  export FALCON_TRACE_Y=100
  export FALCON_TRACE_Z=100
  export FALCON_TRACE_SCALE=2.0
  else # slices for the complex phantom
  export FALCON_TRACE_X=33
  export FALCON_TRACE_Y=47
  export FALCON_TRACE_Z=34
  export FALCON_TRACE_SCALE=2.0
  fi
fi

if [ -z $THRESHOLD ];then
 THRESHOLD=0.1
fi


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

verify_files  ${SAMPLE_BRAIN}

if [ -z $OMP_NUM_THREADS ];then OMP_NUM_THREADS=1;fi

echo OMP_NUM_THREADS=${OMP_NUM_THREADS}

#export NIIKVERBOSE=3

# create zero
run_cmd minccalc -express 0 -byte -clobber ${SAMPLE_BRAIN} ${tempdir}/zero.mnc
# create one
run_cmd minccalc -express 1 -byte -clobber ${SAMPLE_BRAIN} ${tempdir}/one.mnc

# extract GW interface
run_cmd ${BINDIR}/falcon_math thresh  -in=${SAMPLE_BRAIN} -thresh=${GM} -out=${tempdir}/GW.mnc

# extract 'brain'
run_cmd ${BINDIR}/falcon_math thresh  -in=${SAMPLE_BRAIN} -thresh=${BCK} -out=${tempdir}/BRAIN.mnc

# extract 'non zero'
run_cmd ${BINDIR}/falcon_math thresh  -in=${SAMPLE_BRAIN} -thresh=0.01 -out=${tempdir}/NONZERO.mnc


# close
run_cmd ${BINDIR}/falcon_math closebrain -in=${tempdir}/GW.mnc -out=${tempdir}/GW_close.mnc    -radius=1

# dilate WM
run_cmd ${BINDIR}/falcon_math dilate -in=${tempdir}/GW_close.mnc -out=${tempdir}/GW_dilate.mnc -radius=1

# dilate brain
# dilate
run_cmd ${BINDIR}/falcon_math dilate -in=${tempdir}/BRAIN.mnc -out=${tempdir}/BRAIN_dilate.mnc -radius=1


# extract CSF (everything above 0, outside of BRAIN
run_cmd ${BINDIR}/falcon_math maskout -img=${tempdir}/NONZERO.mnc -mask=${tempdir}/BRAIN_dilate.mnc -out=${tempdir}/CSF.mnc

# laplace
run_cmd ${BINDIR}/falcon_math laplacemap \
                 -imglist=${tempdir}/GW.mnc,${tempdir}/GW_dilate.mnc -out=${tempdir}/GW_laplace.mnc -xyz=0.5,0.5,0.5

# shrink-wrap
echo Shrink-Wrap
run_cmd ${BINDIR}/falcon_cortex_shrinkwrap ${tempdir}/GW_dilate.mnc \
  ${tempdir}/GW_dilate_obj.ply \
  --iter 50  --step 2.0 --elen 1.0 --clob

# initics
echo Init ICS
run_cmd ${BINDIR}/falcon_cortex_initics ${SAMPLE_BRAIN} ${tempdir}/GW.mnc ${tempdir}/GW_laplace.mnc \
          ${tempdir}/GW_dilate_obj.ply ${tempdir}/GW_init_ics.ply --clob

echo
echo extract pial surface
# extract pial surface
run_cmd ${BINDIR}/falcon_cortex_initocs ${SAMPLE_BRAIN} \
  ${tempdir}/BRAIN_dilate.mnc ${tempdir}/CSF.mnc ${tempdir}/GW.mnc \
  ${tempdir}/GW_init_ics.ply ${tempdir}/GW_init_ocs.ply \
  --iter 50 --clob --max-thick 2.0

echo
echo "refine (deform) both WM and pial surfaces"

# refine (deform) both WM and pial surfaces
run_cmd ${BINDIR}/falcon_cortex_refine ${SAMPLE_BRAIN} \
                       ${tempdir}/BRAIN_dilate.mnc  ${tempdir}/CSF.mnc  ${tempdir}/GW.mnc \
                       ${tempdir}/zero.mnc \
                       ${tempdir}/GW_init_ics.ply ${tempdir}/GW_init_ocs.ply \
                       ${OUT}_ics.ply ${OUT}_ocs.ply \
                      -gradient-FWHM    1.0 1.0 \
                      -divergence-FWHM  1.0 1.0 \
                      -wimag     0.5     0.5 \
                      -wcurv     0.1     0.1 \
                      -wsmooth   1.0     1.0 \
                      -wssmooth  0.1     0.1 \
                      -wtsmooth  1.0     1.0 \
                      -wgrad     0.0     0.0 \
                      -wflux     0.1     0.1 \
                      -wprox     0.5     0.5 \
                      -wbrain    1.0     1.0 \
                      -wvent     1.0     1.0 \
                      -wabs      0.0     0.0 \
                      -tmin      0.0     0.0 \
                      -tmax      5.0     5.0 \
                      -tsigma    1.0     1.0 \
                      -supdate   0.0     0.0 \
                      -pmin      0.5    \
                      -iter      100    \
                      -delta     0.5    \
                      -iter2     5      \
                      -apply     0.2    \
                      -remesh    0     \
                      -rungekutta


echo "Calculating cortical thickness"
# calculate thickness
run_cmd ${BINDIR}/falcon_cortex_calc_thickness ${OUT}_ics.ply ${OUT}_ocs.ply ${OUT}_thickness.txt.gz  --clob

# calculate average thickness using awk and compare against expected
if [ ! -z $THICK ];then
  zcat ${OUT}_thickness.txt.gz| awk "{s+=\$1;c+=1} END {th=s/c; if( (th-$THICK)>$THRESHOLD || (th-$THICK)<-$THRESHOLD ) {print \"Thickness: \" th \" expected $THICK\" ;exit 10 } else print \"Thickness OK\" }"
fi
