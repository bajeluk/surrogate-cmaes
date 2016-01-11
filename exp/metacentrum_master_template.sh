#!/bin/bash
# Metacentrum task manager for deployd Matlab-compiled 'metacentrum_task_matlab' binary

# usage:
#   ./metacentrum_master_template.sh EXPID META_QUEUE [ID1] [ID2]...
#
# where
#   EXPID       -- string with experiment's unique ID
#   META_QUEUE  -- string with walltime for Metacentrum (2h/4h/1d/2d/1w)
#   ID1,ID2,... -- integers defining numeric IDs of concrete experiments to run
#               all IDs are taken from file 'allids.txt' if no IDs supplied on
#               command-line (allids.txt is expected in experiment's directory
#               and should be produced by the experiment generator script)
#
# settings within this file:
#   EXPPATH_SHORT  -- $CWD/experiments

# EXPID = ExperimentID (string)
EXPID=$1

# QUEUE = Metacentrum walltime (2h/4h/1d/2d/1w) -- queue will be decided accordingly
QUEUE=$2

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Load global settings and variables
. $CWD/bash_settings.sh

# IDs of the tasks to be submitted (CWD == path where the current file is)
EXPPATH_SHORT="$CWD/experiments"
# SCRIPT=`basename ${BASH_SOURCE[0]}`
if [ $# -gt 2 ]; then
  shift; shift;
  IDS=$*
else
  IDS=`cat $EXPPATH_SHORT/$EXPID/allids.txt`
fi


#
# Packing of current sources
#
# check that all the deployed files and directories exist; $FILES_TO_DEPLOY is
#   set in exp/bash_settings.sh
for file in $FILES_TO_DEPLOY; do
  if [ ! -f $file -a ! -d $file ]; then
    echo "Error: the specified file to deploy '$file' does not exist. Exiting."
    exit 1;
  fi
done
#
lastdir=`pwd`
cd "$CWD/.."
mkdir -p "$DEPLOY_DIR"
if [ -f "$DEPLOY_DIR/$DEPLOY_ARCHIVE" ]; then
  echo "Warning: tar archive with the current EXPID already exists, skipping packaging"
  echo "=======  and this old archive will be used."
  echo " :: Delete the tar archive if you want to update the source package."
  echo " :: You can interrupt the Task scheduling within the next 5 seconds."
  echo ""
  ls -l "$DEPLOY_DIR/$DEPLOY_ARCHIVE"
  sleep 5
else
  # Packaging itself
  tar -cvf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" $FILES_TO_DEPLOY
  echo ""
fi
cd "$lastdir"
#
# End of packing of sources
#


export EXPID
export EXPPATH_SHORT
export ID

for ID in $IDS; do
  qsub -N "${EXPID}__${ID}" -l "walltime=$QUEUE" -v EXPID,ID,EXPPATH_SHORT $EXPPATH_SHORT/$EXPID/binary_task.sh
  if [ ! $? -eq 0 ] ; then
    echo "Nepodarilo se zadat ulohu segment ${ID}! Koncim."; exit 1
  else
    echo "Job ${EXPID} / ${ID} submitted to the '$QUEUE' queue."
    touch "$EXPPATH_SHORT/$EXPID/queued_$ID"
  fi
done
