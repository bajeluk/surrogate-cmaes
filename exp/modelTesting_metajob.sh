#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l scratch=4gb
#PBS -l matlab=1,matlab_Statistics_Toolbox=1,matlab_Optimization_Toolbox=1

# it suppose the following variables set:
#
#   FUNC           -- list of integers of functions
#   DIM            -- list of integers of dimensions
#   INST           -- list of instances to process
#   MATLAB_FCN     -- matlab script to be called
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments

module add matlab
MATLAB=matlab

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

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
if [ -z "$MATLAB_FCN" ] ; then
  echo "Error: Task MATLAB_FCN is not known"; exit 1
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

cd "$EXPPATH_SHORT/../.."

echo "====================="
echo -n "Current dir:    "; pwd
echo    '$HOME =        ' $HOME
echo "====================="

######### CALL #########
#
echo '##############'
echo Will be called: $MATLAB -nodisplay -singleCompThread -r "func=[$FUNC]; dims=[$DIM]; instances=[$INST]; ${MATLAB_FCN}; exit(0);";
echo '##############'

$MATLAB -nodisplay -singleCompThread -r "func=[$FUNC]; dims=[$DIM]; instances=[$INST]; ${MATLAB_FCN}; exit(0);";
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

