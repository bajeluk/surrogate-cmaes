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

# create list of files
FILE_LIST=/tmp/stdout_file_list
# find is used instead of ls due to possible large number of stdout files
find . -maxdepth 1 -type f -iname "${EXPID}__*" > $FILE_LIST

tarfile="$EXPPATH/$EXPID/${EXPID}__stdout.tgz"
if [ -s "${FILE_LIST}" ]; then
  if [ -f "$tarfile" ]; then
    unzipped="${tarfile/.tgz/.tar}"
    gunzip "$tarfile" # > "$unzipped"
    tar rf "$unzipped" -T ${FILE_LIST}
    # remove stdout files
    # xargs with "rm" is used due to possible large number of stdout files
    # -I{} rm {} is due to absence of newlines in FILE_LIST
    xargs -a ${FILE_LIST} -I{} rm {}
    gzip -c "$unzipped" > "$tarfile" && rm "$unzipped"
  else
    tar czf "$tarfile" -T ${FILE_LIST}
    # remove stdout files
    # xargs with "rm" is used due to possible large number of stdout files
    # -I{} rm {} is due to absence of newlines in FILE_LIST
    xargs -a ${FILE_LIST} -I{} rm {}
  fi
else
  echo "No files to process."
fi

# remove list of stdout files
rm $FILE_LIST
cd "$dir"
