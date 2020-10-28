#! /bin/bash
set -e -o pipefail

# default
cmap=-spectral
title=""
column=""

if [[ -z "$FALCON_HOME" ]];then
  FALCON_HOME="${MINC_TOOLKIT}"
fi

if [[ -z "$FALCON_DATA"  ]]; then
  FALCON_DATA="${FALCON_HOME}/share/falcon"
fi

if [[ -z "$FALCON_BIN"  ]]; then
  FALCON_BIN="${FALCON_HOME}/bin"
fi

function Usage {
  cat <<EOF
 Usage: 
  `basename $0` <input_lt.ply> <input_measurement_lt.txt/csv> <input_rt.ply> <input_measurement_rt.txt/csv> <output.png> \
    [foreground] [background] 
  --- Optional parameters  ---
  -min <m> 
  -max <m> 
  -title <title>
  -spectral - use spectral colour map (default)
  -atrophy  - use atrophy colour map
  -summer   - use summer colour map
  -jacobian - use jacobian colour map
  -gray     - use gray colour map
  -sphere   - output on the spherical map, using spherical coordinates instead of x,y,z
  -column  <n> - specify column name from txt file to use (default will use the first column)
  -discrete  - don't interpolate values 
  -levels <n> - number of different color levels 
EOF
}

####################################
# SCRIPT STARTS
####################################

if [ $# -eq 0 ]; then Usage; exit 1; fi

while [ $# -gt 0 ]; do
  if   [ $1 = -help ]; then Usage; exit 1
  elif [ $1 = -u ]; then Usage; exit 1
  elif [ $1 = -min  ]; then th_min=$2;  shift 2
  elif [ $1 = -max  ]; then th_max=$2;  shift 2
  elif [ $1 = -spectral  ]; then cmap=$1;  shift 
  elif [ $1 = -atrophy  ]; then cmap=$1;  shift 
  elif [ $1 = -summer  ]; then cmap=$1;  shift 
  elif [ $1 = -jacobian  ]; then cmap=$1;  shift 
  elif [ $1 = -gray  ]; then cmap=$1;  shift
  elif [ $1 = -sphere ];then sphere=--sphere;shift
  elif [ $1 = -title  ]; then title="$2";  shift 2
  elif [ $1 = -column  ]; then column="--column $2"; shift 2
  elif [ $1 = -discrete  ]; then discrete="--discrete"; shift
  elif [ $1 = -levels  ]; then levels="--levels $2"; shift 2
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

if   [ ${#args[@]} -lt 5 ]; then echo -e "ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt 7 ]; then echo -e "ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi

in_off_lt=${args[0]}
in_txt_lt=${args[1]}
in_off_rt=${args[2]}
in_txt_rt=${args[3]}
output=${args[4]}

if [ ${#args[@]} -gt 5 ]; then
foreground=${args[5]}
background=${args[6]}
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

if [ -z $foreground ];then
 foreground="white"
fi


if [ -z $background ];then
 background="black"
fi


function deal_with_csv {
  in=$1
  ext="${in##*.}"
  
  if [ "$ext" == "gz" ];then
    b=$(basename $1 .gz)
    gzip -d -c $1 > $tempdir/$b
    in=$tempdir/$b
    ext="${in##*.}"
  fi
  if [ "$ext" == "csv" ];then
    b=$(basename $1 .csv)
    # strip first line and take the first column (thickness)
    cat $in|tail -n +2 |cut -d , -f 1 > $tempdir/$b.txt
    echo $tempdir/$b.txt
  else
    echo $in
  fi
}

template=$FALCON_DATA/data/template.pov
colourbar=$FALCON_DATA/data/colourbar${cmap}.png

if [ ! -e $template ];then
 echo "Missing ${template}" 2>&1 
 exit 1
fi



# convert off to povray mesh, get the range
declare -a range

if [ -z $th_min ] || [ -z $th_max ];then
  declare -a range1
  declare -a range2
  range1=( $(falcon_csv_stats $in_txt_lt --range $column)  )
  range2=( $(falcon_csv_stats $in_txt_rt --range $column)  )
  
  range=($(echo "if(${range1[0]}<${range2[0]}) ${range1[0]} else ${range2[0]}"|bc -l) $(echo "if(${range1[1]}>${range2[1]}) ${range1[1]} else ${range2[1]}"|bc -l) )
 else
   range=($th_min $th_max)
fi

$FALCON_BIN/falcon_off2pov --colorize $in_txt_lt $column $in_off_lt $tempdir/brain_lt.pov --object brain_lt --clobber --metallic $sphere -$cmap --min ${range[0]} --max ${range[1]} $discrete $levels
$FALCON_BIN/falcon_off2pov --colorize $in_txt_rt $column $in_off_rt $tempdir/brain_rt.pov --object brain_rt --clobber --metallic $sphere -$cmap --min ${range[0]} --max ${range[1]} $discrete $levels


h=$(identify -format %h $colourbar)

pos1=$(($h+60))
pos2=30

convert -border 30 \
  -fill ${foreground} -background ${background} -bordercolor ${background} -font Times-Bold -pointsize 20 \
  -draw "text 2,$pos2 '${range[1]}'" \
  -draw "text 2,$pos1 '${range[0]}'"  \
  -alpha remove -alpha off  \
  $colourbar $tempdir/bar.miff


Width=320
Height=240

# modify template for different orientations
for angle in Front Back Top Bottom Left Right; do

    echo "#version 3.6;"> $tempdir/render_${angle}.pov

    case "$angle" in
        Top|Bottom)
            height=$Width
            width=$Height
            min=$Height

            echo "#include \"$tempdir/brain_lt.pov\"" >> $tempdir/render_${angle}.pov
            echo "#include \"$tempdir/brain_rt.pov\"" >> $tempdir/render_${angle}.pov
            sed -e "s/ANGLE/$angle/" -e "s/WIDTH/$width/" -e "s/HEIGHT/$height/" -e "s/MIN/$min/" $template >> $tempdir/render_${angle}.pov

            echo -e "object {\n  brain_lt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov
            echo -e "object {\n  brain_rt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov

            povray  $tempdir/render_${angle}.pov +o$tempdir/render_${angle}.png -D -GA +A +H${height} +W${width} -V  +UA -Q12 # >/dev/null  2>&1 
            convert $tempdir/render_${angle}.png -trim +repage $tempdir/render_${angle}.png
            convert $tempdir/render_${angle}.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}.png
            ;;
        Front|Back)
            height=$Height
            width=$Height
            min=$Height
            echo "#include \"$tempdir/brain_lt.pov\"" >> $tempdir/render_${angle}.pov
            echo "#include \"$tempdir/brain_rt.pov\"" >> $tempdir/render_${angle}.pov
            sed -e "s/ANGLE/$angle/" -e "s/WIDTH/$width/" -e "s/HEIGHT/$height/" -e "s/MIN/$min/" $template >> $tempdir/render_${angle}.pov
            echo -e "object {\n  brain_lt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov
            echo -e "object {\n  brain_rt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov

            povray $tempdir/render_${angle}.pov +o$tempdir/render_${angle}.png -D -GA +A +H${height} +W${width} -V  +UA -Q12 # >/dev/null  2>&1 
            convert $tempdir/render_${angle}.png -trim +repage $tempdir/render_${angle}.png
            convert $tempdir/render_${angle}.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}.png
            ;;

        Left|Right)
            height=$Height
            width=$Width
            min=$Height
            echo "#version 3.6;"> $tempdir/render_${angle}_lt.pov
            echo "#include \"$tempdir/brain_lt.pov\"" >> $tempdir/render_${angle}_lt.pov

            echo "#version 3.6;"> $tempdir/render_${angle}_rt.pov
            echo "#include \"$tempdir/brain_rt.pov\"" >> $tempdir/render_${angle}_rt.pov
            
            sed -e "s/ANGLE/$angle/" -e "s/WIDTH/$width/" -e "s/HEIGHT/$height/" -e "s/MIN/$min/" $template >> $tempdir/render_${angle}_lt.pov
            sed -e "s/ANGLE/$angle/" -e "s/WIDTH/$width/" -e "s/HEIGHT/$height/" -e "s/MIN/$min/" $template >> $tempdir/render_${angle}_rt.pov
            echo -e "object {\n  brain_lt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}_lt.pov
            echo -e "object {\n  brain_rt\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}_rt.pov

            povray $tempdir/render_${angle}_lt.pov +o$tempdir/render_${angle}_lt.png -D -GA +A +H${height} +W${width} -V  +UA -Q12 # >/dev/null  2>&1 
            povray $tempdir/render_${angle}_rt.pov +o$tempdir/render_${angle}_rt.png -D -GA +A +H${height} +W${width} -V  +UA -Q12 # >/dev/null  2>&1 
            convert $tempdir/render_${angle}_lt.png -trim +repage $tempdir/render_${angle}_lt.png
            convert $tempdir/render_${angle}_rt.png -trim +repage $tempdir/render_${angle}_rt.png
            
            convert $tempdir/render_${angle}_lt.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}_lt.png
            convert $tempdir/render_${angle}_rt.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}_rt.png
            
#             montage -geometry +0+0 -tile 2x -background ${background} -bordercolor ${background} -fill ${background} \
#                     $tempdir/render_${angle}_lt.png $tempdir/render_${angle}_rt.png $tempdir/render_${angle}.png
            ;;
    esac
done

montage \
 -geometry ${Width}x${Height}+5+5 \
 -tile 2x4 \
 -background ${background} \
 -fill   ${foreground} \
 -pointsize 20 \
 -font Times-Bold \
 -label Top     $tempdir/render_Top.png \
 -label Bottom  $tempdir/render_Bottom.png \
 -label "Left_LH"    $tempdir/render_Left_lt.png \
 -label "Right_LH"   $tempdir/render_Right_lt.png \
 -label "Left_RH"   $tempdir/render_Left_rt.png \
 -label "Right_RH"   $tempdir/render_Right_rt.png \
 -label Front   $tempdir/render_Front.png \
 -label Back    $tempdir/render_Back.png \
 $tempdir/pik.miff

#if [ ! -z "$title" ];then
#title="-title '${title}'"
#fi

montage \
 -geometry +0+0 \
 -tile 2x1 \
 -background ${background} \
 -fill   ${foreground} \
 -stroke ${foreground} \
 -title  "$title"  \
 $tempdir/pik.miff \
 $tempdir/bar.miff \
 $output
