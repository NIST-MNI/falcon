#! /bin/bash
set -e -o pipefail

cmap=-spectral
title=""

if [[ -z "$FALCON_DATA"  ]]; then
  FALCON_DATA="${FALCON_HOME}/share/falcon"
fi

if [[ -z "$FALCON_BIN"  ]]; then
  FALCON_BIN="${FALCON_HOME}/bin"
fi

function Usage {
  cat <<EOF
Usage: `basename $0` <input.off> <input_measurement.txt> <output.png> [foreground] [background]
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
EOF
}


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
  elif [ $1 = -column  ]; then column="--column $2";  shift 2
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

if   [ ${#args[@]} -lt 2 ]; then echo -e "ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt 4 ]; then echo -e "ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi

in_off=${args[0]}
in_txt=${args[1]}
output=${args[2]}

if [ ${#args[@]} -gt 2 ]; then
foreground=${args[3]}
background=${args[4]}
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

if [ -z $foreground ];then
 foreground="white"
fi

if [ -z $background ];then
 background="black"
fi

template=$FALCON_DATA/data/template.pov
colourbar=$FALCON_DATA/data/colourbar${cmap}.png

if [ ! -e $template ];then
 echo "Missing ${template}" 2>&1 
 exit 1
fi


# figure out measurement range for colouring
declare -a range
if [ -z $th_min ] || [ -z $th_max ];then
  range=( $(falcon_csv_stats $in_txt --range $column)  )
else
  range=($th_min $th_max)
fi

h=$(identify -format %h $colourbar)

pos1=$(($h+60))
pos2=30

convert -border 30 -fill ${foreground} -background ${background} -bordercolor ${background} -font Times-Bold -pointsize 20 -draw "text 2,$pos2 '${range[1]}'" -draw "text 2,$pos1 '${range[0]}'"  -alpha remove -alpha off  $colourbar $tempdir/bar.miff

# convert off to povray mesh
$FALCON_BIN/falcon_off2pov $in_off $tempdir/brain.pov --object brain --clobber $sphere \
  --colorize $in_txt --min ${range[0]} --max ${range[1]} -$cmap


Width=320
Height=240

# modify template for different orientations
for angle in Front Back Top Bottom Left Right; do

    case "$angle" in 
        Top|Bottom)
            height=$Width
            width=$Height
            min=$Height
            ;;
        Left|Right)
            height=$Height
            width=$Width
            min=$Height
            ;;
        Front|Back)
            height=$Height
            width=$Height
            min=$Height
            ;;
    esac
    
    echo "#include \"$tempdir/brain.pov\"" > $tempdir/render_${angle}.pov
    sed -e "s/ANGLE/$angle/" -e "s/WIDTH/$width/" -e "s/HEIGHT/$height/" -e "s/MIN/$min/" $template >> $tempdir/render_${angle}.pov
    echo -e "object {\n  brain\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov

    povray $tempdir/render_${angle}.pov +o$tempdir/render_${angle}.png -D -GA +A +H${height} +W${width} -V  +UA 
    convert $tempdir/render_${angle}.png -trim +repage $tempdir/render_${angle}.png
    convert $tempdir/render_${angle}.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}.png
done


montage \
 -geometry ${Width}x${Height}+5+5 \
 -tile 2x3 \
 -background ${background} \
 -fill   ${foreground} \
 -pointsize 20 \
 -font Times-Bold \
 -label Top     $tempdir/render_Top.png \
 -label Bottom  $tempdir/render_Bottom.png \
 -label Left    $tempdir/render_Left.png \
 -label Right   $tempdir/render_Right.png \
 -label Front   $tempdir/render_Front.png \
 -label Back    $tempdir/render_Back.png \
 $tempdir/pik.miff

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
