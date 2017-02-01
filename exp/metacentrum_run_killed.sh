#!/bin/bash

if [ $1 == "-h" ]; then
  echo "usage:"
  echo "  $0 [EXPID] [NUMBER]"
  echo ""
  echo "where [NUMBER] is the highest possible experiment ID# (accessible from allids.txt)"
  exit 0
fi

EXPID="$1"
HIGHID="$2"
TMPFILE="/tmp/meta_run_killed_$$"

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd $CWD/experiments/$EXPID

ls -1 ${EXPID}_results_*[0-9].mat calculating_* > $TMPFILE

for i in `seq 1 $HIGHID`; do 
  RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
  # if ! grep -q "$RES_FILE"'$' $TMPFILE && ! grep -q "calculating_${i}"'$' $TMPFILE; then
  if ! grep -q "$RES_FILE"'$\|'"calculating_${i}"'$' $TMPFILE; then
    echo -n "$i "
  fi
done

echo ""
rm $TMPFILE
   
exit 0
