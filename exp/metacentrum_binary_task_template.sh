#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -l scratch=1gb

# it suppose the following variables set:
#
#   ID             -- integer idetificating this task's ID
#   EXPID          -- string with unique identifier of the whole experiment
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments
# 

# MATLAB Runtime environment
export LD_LIBRARY_PATH=/storage/plzen1/home/bajeluk/bin/mcr/v90/runtime/glnxa64:/storage/plzen1/home/bajeluk/bin/mcr/v90/bin/glnxa64:/storage/plzen1/home/bajeluk/bin/mcr/v90/sys/os/glnxa64:$LD_LIBRARY_PATH

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

cd "$EXPPATH_SHORT/.."
# ulimit -t unlimited

######### CALL #########
#
./metacentrum_task_matlab "$EXPID" "$EXPPATH_SHORT" $ID
#
########################

echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== FINISHED ===="
