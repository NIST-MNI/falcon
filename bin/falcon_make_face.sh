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
foreground="white"
background="black"
threshold="--bimodal"

function Usage {
  cat <<EOF
Usage: `basename $0` <input.mnc> <output.png/jpg/...> [foreground] [background]
  --- Optional parameters  ---
  -title <title>
  -width  <n>, default $Width
  -height <n>, default $Height
  -distance <n>, default $Distance
  -threshold <f>, default using bimodalT
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
  elif [ $1 = -threshold  ]; then threshold="--threshold $2";  shift 2
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

if   [ ${#args[@]} -lt 1 ]; then echo -e "ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt 3 ]; then echo -e "ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi

in_mnc=${args[0]}
output=${args[1]}

if [ ${#args[@]} -gt 2 ]; then
foreground=${args[2]}
background=${args[3]}
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

set -xe

# extract surface using bimodalT
itk_morph $threshold --exp 'D[2] E[2]' $in_mnc $tempdir/mask.mnc

# extract surface using shrinkwrap
falcon_cortex_shrinkwrap  \
    $tempdir/mask.mnc \
    $tempdir/surface.ply \
    --start 5.0 --elen 3.0 --iter 40 --step 4.0 

# render face
$FALCON_BIN/falcon_face_render.sh \
  $tempdir/surface.ply \
  $output \
  $foreground \
  $background \
  -title "$title" \
  -width "$Width" \
  -height "$Height" \
  -distance "$Distance"
