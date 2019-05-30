#! /bin/bash


TOOLKIT=$1
OUT=$2



if [ -z $OUT ];then
  echo Usage : $0 MINCTOOLKIT OUT
  exit 1
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

source $TOOLKIT/minc-toolkit-config.sh

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1

WM="-fill_value 200  -edge_value 200 -nelements 200 200 200 -background 0 -step 1 1 1 -continuous -start -100 -100 -100 -clobber -float "
BCK="-fill_value 10  -edge_value 10  -nelements 200 200 200 -background 0 -step 1 1 1 -continuous -start -100 -100 -100 -clobber -float "

make_phantom $BCK -ellipse -center 0 0 0    -width 160 160 160  $tempdir/background.mnc    
make_phantom $WM -ellipse -center 0 0 0     -width 100 100 100  $tempdir/sphere_WM.mnc     
make_phantom $WM -rectangle -center 50 0 0  -width 100 100 100  $tempdir/rec_WM.mnc        
make_phantom $WM -ellipse -center 0 0 0     -width 80    5 100  $tempdir/ellipse_WM_1.mnc  
make_phantom $WM -ellipse -center 20 0 0    -width 100   3 100  $tempdir/ellipse_WM_2.mnc  
make_phantom $WM -ellipse -center 50 0 0    -width 5   20  20   $tempdir/ellipse_WM_3.mnc  
make_phantom $WM -ellipse -center 40 0 0    -width 3   15  20   $tempdir/ellipse_WM_4.mnc  
#wait

param2xfm -rotations 0 0 20 $tempdir/rot_z20.xfm -clob
param2xfm -rotations 0 0 -20 $tempdir/rot_z-20.xfm -clob

mincresample -transform $tempdir/rot_z20.xfm  -use_input_sampling $tempdir/ellipse_WM_1.mnc $tempdir/ellipse_WM_1_r1.mnc -clob -q -float  
mincresample -transform $tempdir/rot_z-20.xfm -use_input_sampling $tempdir/ellipse_WM_1.mnc $tempdir/ellipse_WM_1_r2.mnc -clob -q -float 
mincresample -transform $tempdir/rot_z20.xfm  -use_input_sampling $tempdir/ellipse_WM_4.mnc $tempdir/ellipse_WM_4_r1.mnc -clob -q -float 
mincresample -transform $tempdir/rot_z-20.xfm -use_input_sampling $tempdir/ellipse_WM_4.mnc $tempdir/ellipse_WM_4_r2.mnc -clob -q -float 
#wait


mincmath -max $tempdir/ellipse_WM_2.mnc $tempdir/ellipse_WM_1_r1.mnc $tempdir/ellipse_WM_1_r2.mnc $tempdir/ellipse_WM_3.mnc $tempdir/ellipse_WM_4_r1.mnc $tempdir/ellipse_WM_4_r2.mnc $tempdir/ellipses.mnc -clob -q

# minccalc -express 'clamp(A[0]-A[1],0,200)' $tempdir/sphere_WM.mnc $tempdir/rec_WM.mnc $tempdir/combination_WM_.mnc -clob -q 
# mincmath -max -q $tempdir/ellipses.mnc $tempdir/combination_WM_.mnc $tempdir/combination_WM.mnc

itk_g_morph --exp 'D[2]' $tempdir/ellipses.mnc $tempdir/combination_WM_d.mnc --clob

minccalc -express 'A[0]/2.0' $tempdir/combination_WM_d.mnc $tempdir/combination_WM_d2.mnc -q

mincmath -max $tempdir/ellipses.mnc $tempdir/combination_WM_d2.mnc $tempdir/background.mnc $OUT -q -clob
