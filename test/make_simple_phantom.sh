#! /bin/bash

TOOLKIT=$1
OUT_1=$2
OUT_2=$3
OUT_3=$4
OUT_4=$5


if [ -z $OUT_4 ];then
  echo Usage : $0 MINCTOOLKIT OUT_1 OUT_2 OUT_3 OUT_4
  exit 1
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

source $TOOLKIT/minc-toolkit-config.sh

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
unset MINC_FORCE_V2


tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

FG="-fill_value 1 -edge_value 1 -nelements 200 200 200 -background 0 -step 1 1 1 -continuous -start -100 -100 -100 -clobber -float"
#BG="-fill_value 10 -edge_value 10 -nelements 200 200 200 -background 110 -step 1 1 1 -continuous -start -100 -100 -100 -clobber -float -voxel_range 10 110 -real_range 10 110 "
BIN="-fill_value 1 -edge_value 1 -nelements 200 200 200 -background 0 -step 1 1 1  -discrete -start -100 -100 -100 -clobber -byte -no_partial"
set -x
make_phantom $FG  -ellipse -center 0 0 0     -width 100 100 100  $tempdir/fg.mnc -clobber
minccalc -express '10+100*A[0]'       $tempdir/fg.mnc       $OUT_1 -q -clob
minccalc -express '10+100*(1.0-A[0])' $tempdir/fg.mnc $OUT_2 -q -clob

#make_phantom $BG  -ellipse -center 0 0 0     -width 100 100 100  $OUT_2 -clobber
make_phantom $BIN -ellipse -center 0 0 0     -width 100 100 100  $OUT_3 -clobber -labels

itk_distance --signed $OUT_3 $tempdir/dist.mnc
minccalc -express 'abs(A[0])' $tempdir/dist.mnc $OUT_4 -clob -q



