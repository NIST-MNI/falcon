#! /bin/bash

progname=$(basename $0)

function Usage {
  cat <<EOF
   Surface QC tool
   Works best in icbm152 09c space
  Usage: ${progname} <input.mnc> <ics.ply> <ocs.ply> <output.jpg/png/tiff>

  --- Optional parameters  ---
  --min <f> - set image min
  --max <f> - set image max
  --pct <f> - set image max using perectile
EOF
}

if [ $# -eq 0 ]; then Usage; exit 1; fi

while [ $# -gt 0 ]; do
  if   [ $1 = --help ]; then Usage; exit 1
  elif [ $1 = -h ]; then Usage; exit 1
  elif [ $1 = --min  ]; then MIN="--min $2";  shift 2
  elif [ $1 = --max  ]; then MAX="--max $2";  shift 2
  elif [ $1 = --pct  ]; then PCT="--pct $2";  shift 2
  else
    args=( ${args[@]} $1 )
    shift
  fi
done

nargs=4
if   [ ${#args[@]} -lt ${nargs} ]; then echo -e "${progname} ERROR: too few arguments";  echo "  args = ${args[@]}"; exit 9
elif [ ${#args[@]} -gt ${nargs} ]; then echo -e "${progname} ERROR: too many arguments"; echo "  args = ${args[@]}"; exit 9
fi


t1w=${args[0]}
ics=${args[1]}
ocs=${args[2]}
out=${args[3]}

TMPDIR=$(mktemp -d --tmpdir)
trap "rm -rf $TMPDIR" 0 1 2 15
set -e

s=($(mincinfo -dimlength xspace -dimlength yspace -dimlength zspace $t1w) )
set -x

falcon_slice_extract $t1w -o $ics,$ocs $MIN $MAX $PCT \
 -y $(seq -s , 30 $(((${s[1]}-60+15)/18)) $((${s[1]}-29)) ) \
 -x $(seq -s , 20 $(((${s[0]}-40+15)/16)) $((${s[0]}-20)) ) \
 -z $(seq -s , 20 $(((${s[2]}-40+17)/18)) $((${s[2]}-20)) ) \
 $TMPDIR/pic_%s_%03d.png > /dev/null

 montage -geometry +0+0 -tile 8x2 \
  $TMPDIR/pic_x_???.png \
  $TMPDIR/pic_x.miff

 montage -geometry +0+0 -tile 9x2 \
  $TMPDIR/pic_y_???.png \
  $TMPDIR/pic_y.miff

  montage -geometry +0+0 -tile 9x2 \
  $TMPDIR/pic_z_???.png \
  $TMPDIR/pic_z.miff

montage \
-background black -fill white \
-geometry +0+0 -tile 1x3 \
 $TMPDIR/pic_y.miff \
 $TMPDIR/pic_x.miff \
 $TMPDIR/pic_z.miff \
  $out

