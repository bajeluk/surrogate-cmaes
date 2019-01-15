#!/bin/bash
#
# Job-submitting experiment for the first 10D GP model testing
# It is expected to be located in exp/experiments/

export EXPID='exp_DTSmodels_meta_02'
export DATASET="exp/experiments/dataset/DTS_meta_004.mat"

# Enable this option for using Matlab MCR compilated binaries:
export useMCR=1

# !!! THIS HAS TO BE HERE !!!
#     vvvvvvvvvvvvvvvvvvv
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
. $CWD/../metacentrum_testmodels_common.sh
#     ^^^^^^^^^^^^^^^^^^^
# !!! THIS HAS TO BE HERE !!!

# job submittion can be done via one of the two bash functions:
#
# (a) subtask [JOBNAME_SUFFIX]
#
#     submits a job with current-set $ID, $DIM, $FUNC, $INST and $OPTS
#
#
# (b) submit_sequence LOW_IDX STEP HIGH_IDX
#
#     call subtask() with 'modelOptionsIndices' set to a small
#     numbers of model-settings;
#     it submits a job for every $2 consecutive modelOptions' indices
#
#     e.g.: submit 11 4 22
#     submits jobs for the following modelOption indices:
#       - 11:14
#       - 15:18
#       - 19:22
#
#     Important: the HIGH_IDX should be the equal to (k*LOW_IDX) - 1
#                for some k


# critical characters has to be replaced in $OPTS:
# '|' with ',' and "%" with "'"
MODEL_IDS="1 2 3 4 5 6 8 9"

INST="[11 12 13 14 15]"
IDS="[1 2 3 4 5 6 8 9]" # 7 is ADD kernel
#FAILED_LINE=`exp/metacentrum_list.sh -e $EXPID {1..960} | grep -A1 "Errored out" | tail -n1`
#read -a FAILED <<< $FAILED_LINE
FAILED=(888 896 952)

if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "No failed ids found."
  exit
fi

echo -n "Submit ${#FAILED[@]} failed ids (${FAILED[@]})? [y/N] "
read ANS

if [[ $ANS = [yY] ]]; then
  echo "Submitting..."
  sleep 1
else
  exit
fi

ID=1

for DIM in 2 3 5 10 20; do
  case $DIM in
    [235])
      QUEUE="168:00:00"
      ;;
    10)
      QUEUE="168:00:00"
      ;;
    *)
      QUEUE="168:00:00"
      ;;
  esac

  for FUNC in `seq 1 24`; do
    for MODEL_ID in $MODEL_IDS; do
      OPTS="struct(%modelOptionsIndices%|$MODEL_ID)"

      echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
      if [[ $? -eq 0 ]]; then
        subtask
      else
        ID=$((ID+1))
      fi
    done
  done
done
