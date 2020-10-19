#! /bin/bash
set -e -o pipefail


if [[ -z "$FALCON_HOME" ]];then
  FALCON_HOME="${MINC_TOOLKIT}"
fi

if [[ -z "$FALCON_DATA"  ]]; then
  FALCON_DATA="${FALCON_HOME}/share/falcon"
fi

if [[ -z "$FALCON_BIN"  ]]; then
  FALCON_BIN="${FALCON_HOME}/bin"
fi

title=""
Width=320
Height=600
Distance=600


function Usage {
  cat <<EOF
Usage: `basename $0` <input.off/ply> <output.png/jpg/...> [foreground] [background]
  --- Optional parameters  ---
  -title <title>
  -width  <n>, default $Width
  -height <n>, default $Height
  -distance <n>, default $Distance
EOF
}


if [ $# -eq 0 ]; then Usage; exit 1; fi

while [ $# -gt 0 ]; do
  if   [ $1 = -help ]; then Usage; exit 1
  elif [ $1 = -u ]; then Usage; exit 1
  elif [ $1 = -title  ]; then title="$2";  shift 2
  elif [ $1 = -width  ]; then Width="$2";  shift 2
  elif [ $1 = -height  ]; then Height="$2";  shift 2
  elif [ $1 = -distance  ]; then Distance="$2";  shift 2
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

if   [ ${#args[@]} -lt 1 ]; then echo -e "ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt 4 ]; then echo -e "ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi

in_off=${args[0]}
output=${args[1]}

if [ ${#args[@]} -gt 2 ]; then
foreground=${args[2]}
background=${args[3]}
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

if [ -z $foreground ];then
 foreground="white"
fi

if [ -z $background ];then
 background="black"
fi

template=$FALCON_DATA/data/template_face.pov
#colourbar=$FALCON_DATA/data/colourbar${cmap}.png

if [ ! -e $template ];then
 echo "Missing ${template}" 2>&1 
 exit 1
fi

set -xe
# convert off to povray mesh
$FALCON_BIN/falcon_off2pov $in_off $tempdir/brain.pov \
  --object brain --clobber 
  

# modify template for different orientations
for angle in Left Front Right; do

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
    
    echo "#version 3.6;"> $tempdir/render_${angle}.pov
    echo "#include \"$tempdir/brain.pov\"" >> $tempdir/render_${angle}.pov
    sed -e "s/__ANGLE/$angle/" \
        -e "s/__WIDTH/$width/" \
        -e "s/__HEIGHT/$height/" \
        -e "s/__MIN/$min/" \
        -e "s/__DISTANCE/$Distance/" \
        $template >> $tempdir/render_${angle}.pov
    echo -e "object {\n  brain\n  rotate ObjectRotation${angle}\n}\n" >>$tempdir/render_${angle}.pov

    povray $tempdir/render_${angle}.pov +o$tempdir/render_${angle}.png -D -GA +J +R5 +A +H${height} +W${width} -V  +UA
    convert $tempdir/render_${angle}.png -trim +repage $tempdir/render_${angle}.png
    convert $tempdir/render_${angle}.png -background ${background} -bordercolor ${background} -fill ${background} -border 10  -alpha remove -alpha off  $tempdir/render_${angle}.png
done


montage \
 -geometry ${Width}x${Height}+5+5 \
 -tile 3x1 \
 -background ${background} \
 -fill   ${foreground} \
 -pointsize 20 \
 -font Times-Bold \
 $tempdir/render_Left.png \
 $tempdir/render_Front.png \
 $tempdir/render_Right.png \
 $tempdir/pik.miff

montage \
 -geometry +0+0 \
 -tile 1x1 \
 -background ${background} \
 -fill   ${foreground} \
 -stroke ${foreground} \
 -title  "$title"  \
 $tempdir/pik.miff \
 $output
