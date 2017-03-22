#!/bin/sh
#PBS -l select=1:ncpus=1:mem=1500mb:scratch_local=1gb

# it suppose the following variables set:
#
#   ID             -- integer idetificating this task's ID
#   EXPID          -- string with unique identifier of the whole experiment
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments
# 

# MATLAB Runtime environment
# setting LD_LIBRARY_PATH moved into bash_settings.sh

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

# this should be set in 'exp/bash_settings.sh'
# DEPLOY_DIR=deploy
# DEPLOY_ARCHIVE=${EXPID}_src.tar

export SCRATCHDIR
export LOGNAME
EXPPATH="$EXPPATH_SHORT/$EXPID"

if [ -z "$EXPID" ] ; then
  echo "Error: EXPID (experiment ID) is not known"; exit 1
fi
if [ -z "$ID" ] ; then
  echo "Error: Task ID number is not known"; exit 1
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

# clean up the lock-file(s) on exit
trap "rm -f $EXPPATH/calculating_$ID $EXPPATH/queued_$ID" TERM EXIT

#
# Un-pack the archive with current sources and change dir. to therein
#
DEPLOY_FILE="$EXPPATH_SHORT/../../$DEPLOY_DIR/$DEPLOY_ARCHIVE"
if [ ! -f "$DEPLOY_FILE" ]; then
  echo "Error: archive with current sources does not exist. It should be here:"
  echo "  $DEPLOY_FILE"
  exit 1
else
  # defined in exp/bash_settings.sh: RUNDIR="$SCRATCHDIR/surrogate-cmaes"
  mkdir -p "$RUNDIR"
  cd "$RUNDIR"
  echo "====================="
  echo "Unpacking the sources and the previously compiled binary..."
  tar -xf "$DEPLOY_FILE"
fi

# cd "$EXPPATH_SHORT/.."
# ulimit -t unlimited

echo "====================="
echo -n "Current dir:    "; pwd
echo -n "Current node:   "; cat $PBS_NODEFILE
echo    '$HOME =        ' $HOME
echo    '$MCR_CACHE_ROOT = ' $MCR_CACHE_ROOT
echo    "Will be called:" $MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" $ID
echo "====================="

######### CALL #########
#
$MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" $ID
#
# # this is for debug purposes: disable exit on error and direct call matlab
#
# module add matlab
# sed -i 's/^  try/  % try/;s/^  catch err/  return;\n  % catch err/;/ catch err/,/^end/s/^  end/  % end/' exp/bbob_test_01.m 
# matlab -nodisplay -r "dbstop if error; metacentrum_task_matlab('$EXPID','"$EXPPATH_SHORT"', $ID)";
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
