#!/bin/bash
#
# Copy results of an interrupted task from SCRATCH to the destination folder in STORAGE
#
# SYNTAX::
#
# copy_from_scratch.sh EXPID 1972093:minos17 2093483:tarkil8 ...

if [ -n "$1" -a "$1" = "-n" ]; then
  DRY_RUN=1
  shift
else
  DRY_RUN=0
fi
EXPID=$1
shift
JOBS="$*"

SCRATCH_SUBDIR_MASK="job_output/bbcomp_output/*.mat"
USER_SUBDIR="prg/surrogate-cmaes/exp/experiments/$EXPID/bbcomp_output/"
MAX_FILES_TO_COPY=10

# Be aware: while/read does not work with SSH -- seems that SSH captures all the provided STDIN!
# while read TID ID SERVER; do

for JOB in $JOBS; do
  TID=${JOB%:*}
  SERVER=${JOB#*:}
  echo "TID: $TID SERVER: $SERVER"
  echo ssh $USER@$SERVER "cp -vp /scratch/$USER/job_${TID}.*/$SCRATCH_SUBDIR_MASK /storage/plzen1/home/$USER/$USER_SUBDIR"
  if [ $DRY_RUN -ne 1 ]; then
    ssh $USER@$SERVER "FILES=\`ls -tr /scratch/$USER/job_${TID}.*/$SCRATCH_SUBDIR_MASK | head -n$MAX_FILES_TO_COPY\`; test -n \"\$FILES\" && cp -vp \$FILES /storage/plzen1/home/$USER/$USER_SUBDIR"
  else
    ssh $USER@$SERVER "ls -ltr /scratch/$USER/job_${TID}.*/$SCRATCH_SUBDIR_MASK"
  fi
  echo xxxxx result: $? xxxxxx
done
