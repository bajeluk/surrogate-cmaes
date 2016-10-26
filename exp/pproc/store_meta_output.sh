#!/bin/bash
# move metacentrum output from surrogate-cmaes root to .tar file in folder exp/experiments/[exp_id]
# store_meta_output

# initialize

dir=`pwd`
# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $CWD/../..

EXPPATH="exp/experiments"
EXPID=$1

files=`ls ${EXPID}__*`
tarfile="$EXPPATH/$EXPID/${EXPID}__stdout.tgz"
if [ ! -z "$files" ]; then
  if [ -f "$tarfile" ]; then
    tar rzf "$tarfile" ${EXPID}__*
  else
    tar czf "$tarfile" ${EXPID}__*
  fi
  rm ${EXPID}__*
else
  echo "No files to process."
fi

cd "$dir"
