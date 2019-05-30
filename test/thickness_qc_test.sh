#! /bin/bash
set -x
#set -v

BINDIR=$1
SCRIPTDIR=$2
INPUT=$3
OUTPUT=$4

VERBOSE=yes

if [ -z $OUTPUT ];then
  echo "Usage $0 <bindir> <scriptdir> <input_prefix> <output_file>"
  exit 1
fi

tempdir=$(mktemp -d -t FALCON.XXXXXX)
trap "rm -rf $tempdir" 0 1 2 15

function verify_files {
  while [ $# -gt 0 ]; do
    if [ ! -e $1 ];then 
      echo Missing $1
      exit 1
    fi
    shift
  done
}

function run_cmd {
  if [[ $VERBOSE == yes  ]];then
     echo Running: $@
     echo 
  fi
  $@
  if [[ $? != 0 ]]; then
    echo Failed: $@
    exit 1
  else
  echo done
  fi

  return 0
}

in_surface=${INPUT}_ics.ply
in_thickness=${INPUT}_thickness.txt.gz

verify_files  $in_surface $in_thickness

export FALCON_DATA=$SCRIPTDIR/..
export FALCON_BIN=$BINDIR
export FALCON_SCRIPTS=$SCRIPTDIR

# run_cmd strips quote marks
run_cmd $SCRIPTDIR/falcon_off_qc.sh $in_surface $in_thickness $OUTPUT -min 1 -max 3 -title Complex_Phantom
