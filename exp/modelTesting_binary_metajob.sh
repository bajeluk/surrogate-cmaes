#!/bin/sh
#PBS -l select=1:ncpus=1:mem=1500mb:scratch_local=1gb

# it suppose the following variables set:
#
#   FUNC           -- list of integers of functions
#   DIM            -- list of integers of dimensions
#   INST           -- list of instances to process
#   OPTS           -- string with options to be eval()-ed
#   EXPID          -- string with the experiment name
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments

# MATLAB Runtime environment
# setting LD_LIBRARY_PATH moved into bash_settings.sh

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

MATLAB_BINARY_CALL="exp/metacentrum_testmodels"
DATASET="$EXPPATH_SHORT/$EXPID/dataset/DTS_005.mat"

export SCRATCHDIR
export LOGNAME

if [ -z "$FUNC" ] ; then
  echo "Error: Task FUNC numbers are not known"; exit 1
fi
if [ -z "$DIM" ] ; then
  echo "Error: Task DIM numbers are not known"; exit 1
fi
if [ -z "$INST" ]; then
  echo "Warning: INST is empty, default instances will be used."
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

# replace critical characters in $OPTS: '|' with ',' and "%" with "'"
OPTS=`echo $OPTS | tr '%|' "',"`

cd "$EXPPATH_SHORT/../.."
cp "$DATASET" "$SCRATCHDIR"
DATASET=$SCRATCHDIR/`basename "$DATASET"`

echo "====================="
echo -n "Current dir:    "; pwd
echo -n "Current node:   "; cat "$PBS_NODEFILE"
echo    '$HOME =         ' $HOME
echo    '$MCR_CACHE_ROOT = ' $MCR_CACHE_ROOT
echo    '$DATASET =      ' $DATASET
echo "====================="

######### CALL #########
#
echo '##############'
echo Will be called: $MATLAB_BINARY_CALL \"$EXPID\" \"$EXPPATH_SHORT\" \"$FUNC\" \"$DIM\" \"$INST\" \"$OPTS\" \"$DATASET\"
echo '##############'

$MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" "$FUNC" "$DIM" "$INST" "$OPTS" "$DATASET"
#
########################

if [ $? -eq 0 ]; then
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== FINISHED ===="
  rm -rf $SCRATCHDIR/*
  exit 0
else
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== ENDED WITH ERROR! ===="
  exit 1
fi
