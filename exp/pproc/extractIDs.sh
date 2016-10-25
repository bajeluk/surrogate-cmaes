#!/bin/bash
#
# extractIDs.sh -- extract specified IDs from experiment into a new directory
# 
# The output directory is created in the exp/experiments/ and is named
# EXP_ID_[number-of-ids-specified]ids
#
# usage:
# exp/pproc/extractIDs.sh exp_doubleEC_01 [NUMBER OF ALGORITHMS] 3 4 5 6 7 10 11 12 24
#
# where [NUMBER OF ALGORITHMS] is the number of different tested settings
# 
# which will produce exp/experiments/exp_doubleEC_01_9ids
#

dir=`pwd`

EXP_ID=$1
N_ALGS=$2
shift 1

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$CWD/../experiments/"

outputDir=`pwd`"/${EXP_ID}_`expr ${#} - 1`ids"
if [ -d "$outputDir" ]; then
  echo "Output directory already exists. Exiting."
  exit 1
fi
mkdir "$outputDir"

maxID=`sed 's/ *$//;s/.* //' ${EXP_ID}/allids.txt`

allIDs=""
while shift; do
  if [ -z "$1" ]; then break; fi
  id=$1
  ids=`seq $id $N_ALGS $maxID | tr '\n' ',' | sed 's/,[^,]*$//'`
  
  if [ ! -z "$allIDs" ]; then
    allIDs="$allIDs,$ids"
  else
    allIDs="$ids"
  fi
done

echo "Following IDs will be processed:"
echo $allIDs | tr ',' ' '

cp -r ${EXP_ID}/scmaes_params.mat ${EXP_ID}/allids.txt ${EXP_ID}/cmaes_results $outputDir/
eval cp "${EXP_ID}/${EXP_ID}_"*"D_{$allIDs}."* $outputDir/
mkdir $outputDir/bbob_output
eval cp -r "${EXP_ID}/bbob_output/"*"D_{$allIDs}" $outputDir/bbob_output/
eval cp -r "${EXP_ID}/bbob_output/"*"D_{$allIDs}."* $outputDir/bbob_output/

cd $dir
