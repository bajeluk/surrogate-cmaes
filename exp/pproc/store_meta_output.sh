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
    unzipped="${tarfile/.tgz/.tar}"
    gunzip "$tarfile" # > "$unzipped"
    tar rf "$unzipped" ${EXPID}__* && rm ${EXPID}__*
    gzip -c "$unzipped" > "$tarfile" && rm "$unzipped"
  else
    tar czf "$tarfile" ${EXPID}__* && rm ${EXPID}__*
  fi
else
  echo "No files to process."
fi

cd "$dir"
