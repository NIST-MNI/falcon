#! /bin/bash

TOOLKIT=$1
bindir=$2
FG=$3
BG=$4
BIN=$5
DIST=$6
OUT=$7

if [ -z $OUT ];then
  echo Usage : $0 MINCTOOLKIT bindir FG BG BIN DIST OUT
  exit 1    
fi

source $TOOLKIT/minc-toolkit-config.sh

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
export OMP_NUM_THREADS=1

set -ex

#tempdir=$(mktemp -d -t FALCON.XXXXXX)
#trap "rm -rf $tempdir" 0 1 2 15
tempdir=$OUT/tmp
mkdir -p $tempdir

#TODO: add noise?
R_INIT=55
R=49

#TRACING
export FALCON_TRACE_X=80
export FALCON_TRACE_Y=100
export FALCON_TRACE_Z=100
export FALCON_TRACE_SCALE=2.0
export FALCON_TRACE=$tempdir/trace

# make a "global" mask
minccalc -labels -byte -express '1' $BIN $tempdir/all.mnc -clob

if [ ! -e $tempdir/signed_dist.mnc ];then
itk_distance --signed $BIN $tempdir/signed_dist.mnc
fi

if [ ! -e $tempdir/fg_blur.mnc ];then
fast_blur --fwhm 4.0 $FG $tempdir/fg_blur.mnc --float
fast_blur --fwhm 4.0 $BG $tempdir/bg_blur.mnc --float
fi

# make initial sphere
${bindir}/falcon_cortex_shrinkwrap --radius $R_INIT $tempdir/init_sphere.ply --clob
${bindir}/falcon_test_sphere $tempdir/init_sphere.ply $R_INIT --verbose --tolerance 0.5 

# follow inverse distance gradient
if true ;then
${bindir}/falcon_surface_refine \
  $DIST $tempdir/all.mnc $BIN $tempdir/init_sphere.ply $tempdir/dist_grad.ply \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.0 \
  -wsmooth 1.0 \
  -wssmooth 0.1 \
  -wgrad -1.0 \
  -iter 300

${bindir}/falcon_test_sphere $tempdir/dist_grad.ply $R --verbose --output $tempdir/dist_grad_err.ply --tolerance 0.5 
fi

if true ;then
#shrink-wrap
${bindir}/falcon_cortex_shrinkwrap $BIN --obj $tempdir/init_sphere.ply $tempdir/shrinkwrap.ply --clob
${bindir}/falcon_test_sphere $tempdir/shrinkwrap.ply $R --verbose
fi

if true;then
export FALCON_TRACE=$tempdir/trace_in

# shrink initial surface by 1mm
${bindir}/falcon_surface_refine \
  $tempdir/signed_dist.mnc $tempdir/all.mnc $BIN $tempdir/dist_grad.ply $tempdir/dist_grad_in_1mm.ply \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.0 \
  -wsmooth 1.0 \
  -wssmooth 0.1 \
  -wgrad -1.0 \
  -iter 100 \
  -dfmlimit 1.0 

${bindir}/falcon_test_sphere $tempdir/dist_grad_in_1mm.ply $(($R-1)) --verbose


# TODO: check that the surface is of radius R_INIT-1

export FALCON_TRACE=$tempdir/trace_out

${bindir}/falcon_surface_refine \
  $tempdir/signed_dist.mnc $tempdir/all.mnc $BIN $tempdir/dist_grad.ply $tempdir/dist_grad_out_1mm.ply \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.0 \
  -wsmooth 1.0 \
  -wssmooth 0.1 \
  -wgrad 1.0 \
  -iter 100 \
  -dfmlimit 1.0 

${bindir}/falcon_test_sphere $tempdir/dist_grad_out_1mm.ply $(($R+1)) --verbose

# TODO: check that the surface is of radius R_INIT+1
fi



export FALCON_TRACE=$tempdir/trace_out_to_fg
# check that we can go back to sphere
${bindir}/falcon_surface_refine \
  $FG $tempdir/all.mnc $BIN $tempdir/dist_grad_out_1mm.ply $tempdir/dist_grad_out_1mm_to_fg.ply  \
  -gradient-FWHM 1.0 -divergence-FWHM 2.0 \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.2 \
  -wflux 0.5 \
  -wsmooth 1.0 \
  -wssmooth 0.0 \
  -iter 100 
${bindir}/falcon_test_sphere $tempdir/dist_grad_out_1mm_to_fg.ply $R --verbose

export FALCON_TRACE=$tempdir/trace_out_to_bg
# check that we can go back to sphere
${bindir}/falcon_surface_refine \
  $BG $tempdir/all.mnc $BIN $tempdir/dist_grad_out_1mm.ply $tempdir/dist_grad_out_1mm_to_bg.ply  \
  -gradient-FWHM 1.0 -divergence-FWHM 2.0 \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.2 \
  -wflux 0.5 \
  -wsmooth 1.0 \
  -wssmooth 0.0 \
  -iter 100 
${bindir}/falcon_test_sphere $tempdir/dist_grad_out_1mm_to_bg.ply $R --verbose

export FALCON_TRACE=$tempdir/trace_in_to_fg
# check that we can go back to sphere
${bindir}/falcon_surface_refine \
  $FG $tempdir/all.mnc $BIN $tempdir/dist_grad_in_1mm.ply $tempdir/dist_grad_in_1mm_to_fg.ply  \
  -gradient-FWHM 1.0 -divergence-FWHM 2.0 \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.2 \
  -wflux 0.5 \
  -wsmooth 1.0 \
  -wssmooth 0.0 \
  -iter 100 
${bindir}/falcon_test_sphere $tempdir/dist_grad_in_1mm_to_fg.ply $R --verbose

export FALCON_TRACE=$tempdir/trace_in_to_bg
# check that we can go back to sphere
${bindir}/falcon_surface_refine \
  $BG $tempdir/all.mnc $BIN $tempdir/dist_grad_in_1mm.ply $tempdir/dist_grad_in_1mm_to_bg.ply  \
  -gradient-FWHM 1.0 -divergence-FWHM 2.0 \
  -rungekutta \
  -wprox 1.0 \
  -wimag 0.2 \
  -wflux 0.5 \
  -wsmooth 1.0 \
  -wssmooth 0.0 \
  -iter 100 
${bindir}/falcon_test_sphere $tempdir/dist_grad_in_1mm_to_bg.ply $R --verbose
