#!/bin/bash

TMPFILE=/tmp/strip_evals_$$

echo "Do you really want to overwrite the files in current dir?"
read prompt

if [ $prompt != "y" ]; then
  echo "OK, exiting."
  exit 1
fi

for f in `find . -name '*.*dat'`; do
	DIM=`echo $f | sed 's/.*DIM\([0-9]*\)\..*/\1/'`
	echo "============= $f ==============="
	awk '/^[0-9]/ { if ($1 <= '$DIM' * 100) print $0; }; /^[^0-9]/ { print $0 }' $f > $TMPFILE
	#  awk 'BEGIN { $first = 1; }; /^[0-9]/ { if ($1 <= '$DIM' * 100) print $0; else { if ($first == 1) print $0; $first = 0; } }; /^[^0-9]/ { print $0 }' $f
	mv $TMPFILE $f
done
