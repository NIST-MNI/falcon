#! /bin/bash

TOOLKIT=$1
bindir=$2
OUT=$3

if [ -z $OUT ];then
  echo Usage : $0 MINCTOOLKIT bindir OUT
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
R1=20
R2=5

#TRACING
export FALCON_TRACE_X=80
export FALCON_TRACE_Y=100
export FALCON_TRACE_Z=100
export FALCON_TRACE_SCALE=2.0
export FALCON_TRACE=$tempdir/trace


# make sphere 1
${bindir}/falcon_cortex_shrinkwrap --radius $R1 $tempdir/sphere1.ply --clobber --start 0.1 
${bindir}/falcon_test_sphere $tempdir/sphere1.ply $R1 --verbose --tolerance 0.5 

# make sphere 2
${bindir}/falcon_cortex_shrinkwrap --radius $R2 $tempdir/sphere2.ply --clobber --start 0.1 
${bindir}/falcon_test_sphere $tempdir/sphere2.ply $R2 --verbose --tolerance 0.5 

# move sphere 2
param2xfm -translation $R1 0 0 $tempdir/move_x.xfm -clob
${bindir}/falcon_transform_surface $tempdir/sphere2.ply $tempdir/move_x.xfm $tempdir/sphere2_shift.ply --clobber
${bindir}/falcon_test_sphere $tempdir/sphere2_shift.ply $R2 --verbose --tolerance 0.5 --center $R1,0,0

# join two spheres together
${bindir}/falcon_math combine-obj -objlist=$tempdir/sphere1.ply,$tempdir/sphere2_shift.ply \
  -out=$tempdir/sphere1_and_2.ply


# this should pass:
${bindir}/falcon_surface_check $tempdir/sphere1.ply $tempdir/sphere1_.ply --clobber
# this should also pass
${bindir}/falcon_surface_check $tempdir/sphere2_shift.ply $tempdir/sphere2_shift_.ply --clobber

# this should fail
time ${bindir}/falcon_surface_check $tempdir/sphere1_and_2.ply $tempdir/sphere1_and_2_.ply --clobber
echo $?
