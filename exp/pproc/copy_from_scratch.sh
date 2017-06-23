#!/bin/bash
#
# Copy results of an interrupted task from SCRATCH to the destination folder in STORAGE
#
# SYNTAX::
#
# copy_from_scratch.sh EXPID 1972093:minos17 2093483:tarkil8 ...

EXPID=$1
shift
JOBS="$*"

SCRATCH_SUBDIR_MASK="job_output/bbcomp_output/*.mat"
USER_SUBDIR="prg/surrogate-cmaes/exp/experiments/$EXPID/bbcomp_output/"

# Be aware: while/read does not work with SSH -- seems that SSH captures all the provided STDIN!
# while read TID ID SERVER; do

for JOB in $JOBS; do
  TID=${JOB%:*}
  SERVER=${JOB#*:}
  echo "TID: $TID SERVER: $SERVER"
  ssh $USER@$SERVER "cp -vp /scratch/$USER/job_${TID}.*/$SCRATCH_SUBDIR_MASK /storage/plzen1/home/$USER/$USER_SUBDIR"
  echo xxxxx result: $? xxxxxx
done
