#!/bin/sh

i=$1
nMachines=$2

cd ~/prg/surrogate-cmaes
matlabcommand='/afs/ms/@sys/bin/matlab';

EXPPATH_SHORT=`pwd`"/exp/experiments"
EXP_ID="modeltrain_gp_01"
EXPPATH="$EXPPATH_SHORT/$EXP_ID"
LOGS=$EXPPATH"/"$EXP_ID"__log__u-pl"$i".txt"

testMatlabFinished () {
  if [ $? -eq 0 ]; then
    echo `date "+%Y-%m-%d %H:%M:%S"` " " **$EXP_ID** at [u-pl$1] $1 / $2 succeeded. >> ~/WWW/phd/cmaes.txt
  else
    echo `date "+%Y-%m-%d %H:%M:%S"` " " **$EXP_ID** at [u-pl$1] $1 / $2 !!!  FAILED !!! >> ~/WWW/phd/cmaes.txt
  fi
}

echo $matlabcommand -nodisplay -r "testModelParamsGP02('$EXPPATH_SHORT','$EXP_ID', $i, $nMachines); exit(0);"

echo "exit(1);" | nice -n 19 $matlabcommand -nodisplay -r "testModelParamsGP02('$EXPPATH_SHORT','$EXP_ID', $i, $nMachines); exit(0);" 2>&1  | tee -a $LOGS
testMatlabFinished $i $nMachines

sleep 10
