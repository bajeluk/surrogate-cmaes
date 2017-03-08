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
  if [ $3 == "-k" ]; then
    MAXID=`cat $EXPPATH_SHORT/$EXPID/allids.txt | tr ' ' '\n' | tail -2 | head -1`
    echo "We will try to submit not-finished and not-running jobs up to nubmer ${MAXID}..."
    IDS=`$CWD/metacentrum_run_killed.sh $EXPID $MAXID`
    echo "We will submit the following IDs:"
    echo $IDS
  else
    shift; shift;
    IDS=$*
  fi
else
  IDS=`cat $EXPPATH_SHORT/$EXPID/allids.txt`
fi

# Ensure that 'scmaes_params.mat' exists
#
if [ ! -f "$EXPPATH_SHORT/$EXPID/scmaes_params.mat" ]; then
  echo "Warning: 'scmaes_params.mat' does not exist. I will create it by calling"
  echo ""
  echo "matlab -nodisplay -nojvm -r \"expInit('$EXPID'); exit(0);\""
  echo ""
  matlab -nodisplay -nojvm -r "expInit('$EXPID'); exit(0);"
  if [ $? != 0 ]; then
    echo "Matlab ended with error. I'm ending, too."
    exit 1
  fi
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
  # There is already a tar achive with the current EXPID
  echo "#####################"
  echo "Warning: tar archive with the current EXPID already exists, skipping packaging"
  echo "#######  and compiling and this old archive will be used."
  echo "         :: Delete the tar archive if you want to update the source package."
  echo "         :: You can interrupt the Task scheduling within the next 5 seconds."
  echo ""
  ls -l "$DEPLOY_DIR/$DEPLOY_ARCHIVE"
  sleep 5

else
  # The tar archive has to be made.
  # Do we need MCR re-compilation?
  if make -q; then
    # seems there is nothing new to compile!
    echo "#####################"
    echo "Warning: compiled binary '`pwd`/$MATLAB_BINARY_CALL' is up to date,"
    echo "######## so no compilation will happen! It is supposed to be OK."
    echo "         :: If not, delete that file and run this script again."
    echo "         :: You can interrupt the Task scheduling within the next 5 seconds."
    echo ""
    ls -l $MATLAB_BINARY_CALL
    sleep 5
    echo "====================="
    echo "OK, we are using the old binary."
    COMPILED=1
    FILES_TO_DEPLOY="$FILES_TO_DEPLOY $MATLAB_BINARY_CALL"
  else
    # yes, there is something new to compile
    COMPILED=0
  fi

  # Packaging itself
  echo "====================="
  echo "Packaging sources..."
  echo "  tar -cvf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" $FILES_TO_DEPLOY"
  tar -cvf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" $FILES_TO_DEPLOY
  echo "====================="

  if [ $COMPILED -eq 0 ]; then
    # We have to compile the packed sources into the binary
    echo "Preparing sources for compilation"
    COMPILE_DIR=$RUNDIR/compile_$$
    mkdir -p $COMPILE_DIR
    echo "extracting $DEPLOY_DIR/$DEPLOY_ARCHIVE into $COMPILE_DIR ..."
    tar -xf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" -C $COMPILE_DIR
    echo "Compiling..."
    cd $COMPILE_DIR
    module add matlab
    lasthome="$HOME"
    HOME=$COMPILE_DIR

    make

    if [ $? -gt 0 ]; then
      echo "Make failed. Removing '$COMPILE_DIR'." >&2
      rm -rf $COMPILE_DIR
      echo "Exiting." >&2
      exit 1
    fi
    HOME=$lasthome

    echo "====================="
    echo "Copying the resulting binary into the original directory '$CWD/../$MATLAB_BINARY_CALL'"
    cp -p $MATLAB_BINARY_CALL "$CWD/../$MATLAB_BINARY_CALL"
    echo "====================="
    echo "Removing the compilation temporary directory '$COMPILE_DIR'"
    cd "$CWD/.."
    rm -rf $COMPILE_DIR
    echo "====================="
    echo "We are now here:" `pwd`
    echo "Re-packing sources with the compiled binary in '$MATLAB_BINARY_CALL'"
    echo "   tar -cf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" $FILES_TO_DEPLOY $MATLAB_BINARY_CALL"
    tar -cf "$DEPLOY_DIR/$DEPLOY_ARCHIVE" $FILES_TO_DEPLOY $MATLAB_BINARY_CALL
    echo "====================="
  fi
fi
cd "$lastdir"
#
# End of packing of sources & compilation
#

export EXPID
export EXPPATH_SHORT
export ID

IDS_LENGTH=`echo $IDS | wc -w`
for ID in $IDS; do
  qsub -N "${EXPID}__${ID}" -l "walltime=$QUEUE" -v EXPID,ID,EXPPATH_SHORT $EXPPATH_SHORT/$EXPID/binary_task.sh
  if [ ! $? -eq 0 ] ; then
    echo "Nepodarilo se zadat ulohu segment ${ID}! Koncim."; exit 1
  else
    echo "Job ${EXPID}: ${ID} / ${IDS_LENGTH} submitted to the '$QUEUE' queue."
    touch "$EXPPATH_SHORT/$EXPID/queued_$ID"
  fi
done
