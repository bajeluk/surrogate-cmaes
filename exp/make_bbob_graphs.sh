#!/bin/bash

# make_bbob_graphs [SOURCE_DIR] [DEST_DIR]
# 

EXPID="exp_geneEC_08"
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
DATADIR="$CWD/../log/bbob_data/$EXPID"

if [ -n "$1" ]; then
  BBOB_RESULTS_DIR=$1
else
  BBOB_RESULTS_DIR="."
fi
if [ -n "$2" ]; then
  OUTPUT_DIR=$2
else
  PWD=`pwd`
  OUTPUT_DIR="$PWD/ppdata"
fi

if [ -d "$DATADIR" ]; then
  echo Directory with the BBOB data already exists:
  echo "$DATADIR"
  echo ""
  echo "!!! Omitting BBOB data transformation !!!"
  echo ""
else
  mkdir -p $DATADIR
  $CWD/split_bbob_results.sh $BBOB_RESULTS_DIR $DATADIR
fi

# computer dependent -> Lukas will probably rewrite it
#python $CWD/vendor/bbob_pproc/rungeneric.py --expensive $DATADIR/* -o $OUTPUT_DIR
mkdir -p $OUTPUT_DIR
cd $DATADIR
python $CWD/vendor/bbob_pproc/rungeneric.py --expensive --in-a-hurry 0 -o $OUTPUT_DIR *
