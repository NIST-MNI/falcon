#! /bin/bash
set -e 
BINDIR=$1
TMP=$2

if [ -z $BINDIR ];then
  echo "Usage $0 <bindir> [work directory]"
  exit 1
fi


if [ -z $TMP ];then
  tempdir=$(mktemp --tmpdir -d FALCON_XXXXXX)
  trap "rm -rf $tempdir" 0 1 2 15
else
  tempdir=$TMP 
  if [ ! -e $tempdir ];then
    mkdir -p $tempdir
  fi
fi

#export NIIKVERBOSE=3



function verify_volume {
vol=$(mincstats -q -sum  $1 )

# make sure it's 21*21*21 voxels exactly
if [ x"$vol"x != x"$2"x ];then
 echo $(basename $1 .mnc) volume test fail!
 echo Expected $2, got $vol
 exit 1
fi

}

# let's create couple of phantoms first

COMMON="-background 0 -fill_value 0 -edge_value 0 -fill_value 1 -no_partial  -start -25 -25 -25 -nelements 50 50 50 -step 1 1 1 -byte -labels -clobber"

# prepare test dataset
make_phantom $COMMON -center 0 0 0 -width 20 20 20 -rectangle $tempdir/cube_20.mnc
verify_volume $tempdir/cube_20.mnc 9261

make_phantom $COMMON -center 0 0 0 -width 20 20 20 -ellipse $tempdir/ellipse_20.mnc
verify_volume $tempdir/ellipse_20.mnc 3191

make_phantom $COMMON -center 1 1 0 -width 0.5 0.5 0.5 -rectangle $tempdir/hole.mnc
verify_volume $tempdir/hole.mnc 1

minccalc -labels -express 'A[1]>0.5?0:A[0]' $tempdir/ellipse_20.mnc  $tempdir/hole.mnc $tempdir/ellipse_20_with_hole.mnc -clob
verify_volume $tempdir/ellipse_20_with_hole.mnc 3190

# test morphological operators
${BINDIR}/falcon_math dilate -in=${tempdir}/cube_20.mnc -out=${tempdir}/cube_20_dilate.mnc -radius=1
verify_volume $tempdir/cube_20_dilate.mnc 11907

${BINDIR}/falcon_math dilate -in=${tempdir}/ellipse_20.mnc -out=${tempdir}/ellipse_20_dilate.mnc -radius=1
verify_volume $tempdir/ellipse_20_dilate.mnc 4145

${BINDIR}/falcon_math erode -in=${tempdir}/cube_20_dilate.mnc -out=${tempdir}/cube_20_dilate_erode.mnc -radius=1
verify_volume $tempdir/cube_20_dilate_erode.mnc 9261

${BINDIR}/falcon_math erode -in=${tempdir}/ellipse_20_dilate.mnc -out=${tempdir}/ellipse_20_dilate_erode.mnc -radius=1
verify_volume $tempdir/ellipse_20_dilate_erode.mnc 3191

${BINDIR}/falcon_math close -in=${tempdir}/cube_20.mnc -out=${tempdir}/cube_20_close.mnc -radius=1
verify_volume $tempdir/cube_20_close.mnc 9261

${BINDIR}/falcon_math close -in=${tempdir}/ellipse_20.mnc -out=${tempdir}/ellipse_20_close.mnc -radius=1
verify_volume $tempdir/ellipse_20_close.mnc 3191

${BINDIR}/falcon_math open -in=${tempdir}/cube_20.mnc -out=${tempdir}/cube_20_open.mnc -radius=1
verify_volume $tempdir/cube_20_open.mnc 9025

${BINDIR}/falcon_math open -in=${tempdir}/ellipse_20.mnc -out=${tempdir}/ellipse_20_open.mnc -radius=1
verify_volume $tempdir/ellipse_20_open.mnc 3191

${BINDIR}/falcon_math closeholes -in=${tempdir}/cube_20.mnc -out=${tempdir}/cube_20_closeholes.mnc -radius=1
verify_volume $tempdir/cube_20_closeholes.mnc 9261

${BINDIR}/falcon_math closeholes -in=${tempdir}/ellipse_20.mnc -out=${tempdir}/ellipse_20_closeholes.mnc -radius=1
verify_volume $tempdir/ellipse_20_closeholes.mnc 3191

${BINDIR}/falcon_math closeholes -in=$tempdir/ellipse_20_with_hole.mnc -out=$tempdir/ellipse_20_with_hole_closeholes.mnc -radius=1
verify_volume $tempdir/ellipse_20_with_hole_closeholes.mnc 3191

${BINDIR}/falcon_math closebrain -in=$tempdir/ellipse_20_with_hole.mnc -out=$tempdir/ellipse_20_with_hole_closebrain.mnc -radius=1
verify_volume $tempdir/ellipse_20_with_hole_closebrain.mnc 3191

